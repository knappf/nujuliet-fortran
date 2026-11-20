! Spherical Hartree-Fock solver (NuJuliet dev) ... counterpart of src/HF/HF_solver.jl (HF_solve), HF_allocate.jl,
    ! HF_density_operator.jl, HF_orbitals.jl, HF_energy.jl, HF_potential.jl & HF_export.jl ...
    ! Occupations are real numbers, re-determined in every HF iteration from the HF single-particle energies
    ! (HF_orbitals_make), the caller's Orb is updated in place (dev `Orb .= Orb_new`) ...
module nj_hf
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_o1b
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_v3b, only : V3b_no2b
    use nj_orbitals, only : orbitals_export
    use nj_functions, only : Psi_rad_LHO
    use nj_obdm, only : radial_grid
    implicit none
    private

    public :: HF_solve, HF_allocate, HF_density_operator, HF_orbital_ordering, HF_orbitals_make, HF_energy
    public :: HF_summary, HF_SPS_summary, HF_export, HF_radial_potential

contains

    ! Self-consistent iteration of the spherical HF equations ... Orb (occupations) is updated in place ...
    subroutine HF_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,h,C,Rho,Iteration)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(inout) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: T
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t), intent(out) :: h, C, Rho
        integer, intent(out) :: Iteration
        type(pnvector_t) :: SPE, SPE_old
        type(o1b_t) :: H_HF
        real(dp) :: epsilon, dE
        integer :: a_max, Iteration_max
        integer(i8) :: t0

        a_max = a_max_of_calc(Params)
        epsilon = Params%Calc%HF%Tol
        Iteration_max = Params%Calc%HF%Imax

        dE = 1.0_dp
        Iteration = 0
        allocate(SPE%p(a_max), SPE%n(a_max), SPE_old%p(a_max), SPE_old%n(a_max))
        SPE%p = 0.0_dp
        SPE%n = 0.0_dp
        SPE_old%p = 0.0_dp
        SPE_old%n = 0.0_dp

        ! Initial guess ... LHO orbitals with the orbitals_make occupations (possibly fractional) ...
        C%p = identity(a_max)
        C%n = identity(a_max)
        Rho = HF_density_operator(C,Orb)

        write(*,'(/,a)') 'Initializing the iteration of HF mean-field equations ...'
        call timer_start(t0)
        do while ((Iteration < Iteration_max) .and. (dE > epsilon))
            H_HF = HF_allocate(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

            call eig_sym(H_HF%p,SPE%p,C%p)
            call eig_sym(H_HF%n,SPE%n,C%n)

            ! Ordering into the LHO slots & new occupations (in place on Orb) ...
            call HF_orbitals_make(Params,SPE,C,Orb)

            Rho = HF_density_operator(C,Orb)

            dE = (sum(abs(SPE%p - SPE_old%p)) + sum(abs(SPE%n - SPE_old%n))) / real(2 * a_max,dp)
            SPE_old%p = SPE%p
            SPE_old%n = SPE%n
            Iteration = Iteration + 1

            write(*,'(a)') TAB // 'Iteration number: ' // cfmt_d(Iteration,4) // '   dE = ' // cfmt_f(dE,12,9) // ' MeV'
        end do
        call timer_report(t0)

        write(*,'(/,a)') 'Iteration of HF equations with NO2B NN+NNN interaction has finished ...'

        h%p = diagm(SPE%p)
        h%n = diagm(SPE%n)
    end subroutine HF_solve

    pure integer function a_max_of_calc(Params)
        type(params_t), intent(in) :: Params
        a_max_of_calc = (Params%Calc%Nmax + 1) * (Params%Calc%Nmax + 2) / 2
    end function a_max_of_calc

    ! Largest l & 2j of the (partially) occupied orbitals ... dev HF_allocate_indices / HF_energy ...
    subroutine HF_occupied_lj_max(Orb,l_max,j_max)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(out) :: l_max, j_max
        integer :: a
        l_max = -1
        j_max = -1
        do a = 1, size(Orb)
            if (Orb(a)%pO > 1.0e-8_dp .or. Orb(a)%nO > 1.0e-8_dp) then
                l_max = max(l_max,Orb(a)%l)
                j_max = max(j_max,Orb(a)%j)
            end if
        end do
    end subroutine HF_occupied_lj_max

    ! HF single-particle Hamiltonian for given density Rho ... NN + NO2B NNN + kinetic terms ...
        ! Port deviation 5: parallel over the (a,d) pairs (dev: serial (a,d), threads over (b,e)) ...
        ! (b,e) & c are pruned by l_max / j_max of the occupied orbitals (dev; Rho vanishes there, exact) ...
    function HF_allocate(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN) result(H)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t) :: H
        integer :: N_max, N_2max, N_3max, a_max, ad_count, be_count, ad_i, be_i, a, d, b, e, c, f, J, P_ab, P_abc
        integer :: n_a, l_a, j_a, n_d, l_d, j_d, n_b, l_b, j_b, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        integer :: l_max, j_max
        integer, allocatable :: ad(:,:), be(:,:)
        real(dp) :: A_mass, is_j_a_hat, J_j_a_hat, pHSum, nHSum, pRho_be, nRho_be, pRho_cf, nRho_cf
        real(dp) :: ME001, ME101, ME011, ME111, ME113
        character(len=32) :: CMS

        A_mass = real(Params%Calc%A,dp)
        N_max = Params%Calc%Nmax
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        CMS = Params%Calc%CMS
        a_max = (N_max + 1)*(N_max + 2)/2

        call HF_allocate_indices(Params,Orb,ad,ad_count,be,be_count,l_max,j_max)
        call o1b_alloc(H,a_max)

        !$omp parallel do schedule(dynamic) private(a,d,n_a,l_a,j_a,n_d,l_d,j_d,is_j_a_hat,pHSum,nHSum,be_i,b,e, &
        !$omp n_b,l_b,j_b,n_e,l_e,j_e,P_ab,pRho_be,nRho_be,J,J_j_a_hat,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f,pRho_cf,nRho_cf, &
        !$omp ME001,ME101,ME011,ME111,ME113)
        do ad_i = 1, ad_count
            a = ad(1,ad_i)
            d = ad(2,ad_i)
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            n_d = Orb(d)%n
            l_d = Orb(d)%l
            j_d = Orb(d)%j
            is_j_a_hat = 1.0_dp / real(j_a + 1,dp)
            pHSum = 0.0_dp
            nHSum = 0.0_dp
            do be_i = 1, be_count
                b = be(1,be_i)
                e = be(2,be_i)
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                n_e = Orb(e)%n
                l_e = Orb(e)%l
                j_e = Orb(e)%j
                P_ab = mod(l_a + l_b,2) + 1
                if ((2*(n_a + n_b) + l_a + l_b) > N_2max .or. l_b /= l_e .or. j_b /= j_e .or. &
                    (2*(n_d + n_e) + l_d + l_e) > N_2max .or. l_b > l_max .or. j_b > j_max) cycle
                pRho_be = Rho%p(b,e)
                nRho_be = Rho%n(b,e)
                do J = abs(j_d - j_e)/2, (j_d + j_e)/2
                    ! 2-body NN interaction part ...
                    if (P_ab == (mod(l_d + l_e, 2) + 1)) then
                        J_j_a_hat = real(2*J + 1,dp) * is_j_a_hat
                        pHSum = pHSum + J_j_a_hat * pRho_be * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                        pHSum = pHSum + J_j_a_hat * nRho_be * O2b_pn(a,b,d,e,J,P_ab,V_NN,Orb_NN)
                        nHSum = nHSum + J_j_a_hat * nRho_be * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                        nHSum = nHSum + J_j_a_hat * pRho_be * O2b_pn(b,a,e,d,J,P_ab,V_NN,Orb_NN)
                    end if

                    ! 3-body NNN interaction part ...
                    do c = 1, a_max
                        n_c = Orb(c)%n
                        l_c = Orb(c)%l
                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max .or. l_c > l_max) cycle
                        j_c = Orb(c)%j
                        if (j_c > j_max) cycle
                        do f = 1, a_max
                            n_f = Orb(f)%n
                            l_f = Orb(f)%l
                            if ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) > N_3max .or. l_c /= l_f) cycle
                            j_f = Orb(f)%j
                            if (j_c /= j_f) cycle
                            ! dev assigns P_abc = parity(l_d + l_e + l_f) here ... equal to parity(l_a + l_b + l_c),
                            ! since (a,d), (b,e) & (c,f) share l ...
                            P_abc = mod(l_d + l_e + l_f,2) + 1
                            pRho_cf = Rho%p(c,f)
                            nRho_cf = Rho%n(c,f)

                            ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                            ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                            ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                            ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                            ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P_abc,V_NNN,Orb,Orb_NNN)

                            pHSum = pHSum + is_j_a_hat * (0.5_dp*ME113*pRho_be*pRho_cf + 0.25_dp * (ME001 + &
                                    sqrt(1.0_dp/3.0_dp)*ME101 + sqrt(1.0_dp/3.0_dp)*ME011 + ME111/3.0_dp + &
                                    2.0_dp/3.0_dp*ME113)*nRho_be*nRho_cf + &
                                    1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*pRho_be*nRho_cf)

                            nHSum = nHSum + is_j_a_hat * (0.5_dp*ME113*nRho_be*nRho_cf + 0.25_dp * (ME001 + &
                                    sqrt(1.0_dp/3.0_dp)*ME101 + sqrt(1.0_dp/3.0_dp)*ME011 + ME111/3.0_dp + &
                                    2.0_dp/3.0_dp*ME113)*pRho_be*pRho_cf + &
                                    1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*nRho_be*pRho_cf)
                        end do
                    end do
                end do
            end do

            ! 1-body kinetic energy & Center-of-Mass motion correction ...
            if (CMS == 'CMS1+2B') then
                pHSum = pHSum + T%p(a,d) * (1.0_dp - 1.0_dp / A_mass)
                nHSum = nHSum + T%n(a,d) * (1.0_dp - 1.0_dp / A_mass)
            else if (CMS /= 'CMS2B') then
                pHSum = pHSum + T%p(a,d)
                nHSum = nHSum + T%n(a,d)
            end if

            H%p(a,d) = pHSum
            H%p(d,a) = pHSum
            H%n(a,d) = nHSum
            H%n(d,a) = nHSum
        end do
        !$omp end parallel do
    end function HF_allocate

    subroutine HF_allocate_indices(Params,Orb,ab,ab_count,de,de_count,l_max,j_max)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable, intent(out) :: ab(:,:), de(:,:)
        integer, intent(out) :: ab_count, de_count, l_max, j_max
        integer :: a_max, a, b, d, e, Pass

        a_max = a_max_of_calc(Params)
        do Pass = 1, 2
            ab_count = 0
            do a = 1, a_max
                do b = 1, a
                    if (Orb(a)%j == Orb(b)%j .and. Orb(a)%l == Orb(b)%l) then
                        ab_count = ab_count + 1
                        if (Pass == 2) then
                            ab(1,ab_count) = a
                            ab(2,ab_count) = b
                        end if
                    end if
                end do
            end do
            if (Pass == 1) allocate(ab(2,ab_count))
        end do

        de_count = a_max * a_max
        allocate(de(2,de_count))
        de_count = 0
        do d = 1, a_max
            do e = 1, a_max
                de_count = de_count + 1
                de(1,de_count) = d
                de(2,de_count) = e
            end do
        end do

        ! Maximal l & 2j of the occupied orbitals (current occupations, consistent with Rho) ...
        call HF_occupied_lj_max(Orb,l_max,j_max)
    end subroutine HF_allocate_indices

    ! Rho = sum_a O_a C(:,a) C(:,a)^T ... weights are the real occupations ... dev HF_density_operator(C,Orb) ...
    function HF_density_operator(C,Orb) result(Rho)
        type(o1b_t), intent(in) :: C
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t) :: Rho
        integer :: a_max, k, l, m, a
        real(dp) :: pO, nO

        a_max = size(Orb)
        call o1b_alloc(Rho,a_max)
        do m = 1, a_max
            a = Orb(m)%a
            pO = Orb(m)%pO
            nO = Orb(m)%nO
            do l = 1, a_max
                do k = 1, a_max
                    Rho%p(k,l) = Rho%p(k,l) + pO * (C%p(k,a) * C%p(l,a))
                    Rho%n(k,l) = Rho%n(k,l) + nO * (C%n(k,a) * C%n(l,a))
                end do
            end do
        end do
    end function HF_density_operator

    ! Assigns each HF orbital to an LHO orbital slot with matching l & 2j (in energy order) ... first part of the dev
        ! HF_orbitals_make (former HF_orbital_ordering, kept public for the qp solvers) ...
        ! Port deviation: an unassigned slot is a fatal error (Julia leaves Vector{Int}(undef) garbage) ...
    subroutine HF_orbital_ordering(Orb,a_max,C,SPE)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: a_max
        type(o1b_t), intent(inout) :: C
        type(pnvector_t), intent(inout) :: SPE
        integer :: pOrb_order(a_max), nOrb_order(a_max), a, cc, i
        logical :: pOrb_mask(a_max), nOrb_mask(a_max)
        real(dp) :: lp, jp, ln, jn

        pOrb_order = 0
        nOrb_order = 0
        pOrb_mask = .false.
        nOrb_mask = .false.
        do a = 1, a_max
            lp = 0.0_dp
            jp = 0.0_dp
            ln = 0.0_dp
            jn = 0.0_dp
            do i = 1, a_max
                lp = lp + C%p(i,a) * C%p(i,a) * real(Orb(i)%l,dp)
                jp = jp + C%p(i,a) * C%p(i,a) * real(Orb(i)%j,dp)
                ln = ln + C%n(i,a) * C%n(i,a) * real(Orb(i)%l,dp)
                jn = jn + C%n(i,a) * C%n(i,a) * real(Orb(i)%j,dp)
            end do
            do cc = 1, a_max
                if (.not. pOrb_mask(cc) .and. abs(Orb(cc)%l - lp) < 1.0e-8_dp .and. abs(Orb(cc)%j - jp) < 1.0e-8_dp) then
                    pOrb_order(cc) = a
                    pOrb_mask(cc) = .true.
                    exit
                end if
            end do
            do cc = 1, a_max
                if (.not. nOrb_mask(cc) .and. abs(Orb(cc)%l - ln) < 1.0e-8_dp .and. abs(Orb(cc)%j - jn) < 1.0e-8_dp) then
                    nOrb_order(cc) = a
                    nOrb_mask(cc) = .true.
                    exit
                end if
            end do
        end do
        if (any(pOrb_order == 0) .or. any(nOrb_order == 0)) &
            call fatal('HF orbital ordering failed ... HF orbitals with mixed l & j detected ...')

        SPE%p = SPE%p(pOrb_order)
        SPE%n = SPE%n(nOrb_order)
        C%p = C%p(:,pOrb_order)
        C%n = C%n(:,nOrb_order)
    end subroutine HF_orbital_ordering

    ! dev HF_orbitals_make ... ordering of the HF orbitals into the LHO slots, then new occupations by filling the
        ! slots in HF single-particle energy order (fractional last (sub)shell unless HF.ESH) ... Orb updated in place ...
    subroutine HF_orbitals_make(Params,SPE,C,Orb)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(inout) :: SPE
        type(o1b_t), intent(inout) :: C
        type(orb1b_t), intent(inout) :: Orb(:)
        real(dp), allocatable :: pO(:), nO(:)
        integer :: a_max, a

        a_max = a_max_of_calc(Params)
        call HF_orbital_ordering(Orb,a_max,C,SPE)

        allocate(pO(a_max), nO(a_max))
        call HF_occupations_fill(Params,SPE%p,Orb,real(Params%Calc%Z,dp),'proton','Z',pO)
        call HF_occupations_fill(Params,SPE%n,Orb,real(Params%Calc%A - Params%Calc%Z,dp),'neutron','N',nO)

        ! Orb .= Orb_new ... a, n, l, j unchanged, new occupations ...
        do a = 1, a_max
            Orb(a)%pO = pO(Orb(a)%a)
            Orb(a)%nO = nO(Orb(a)%a)
        end do
    end subroutine HF_orbitals_make

    ! Occupation rule of dev HF_orbitals_make (HF_orbitals.jl:66-128) for one isospin ... repeatedly the lowest-SPE
        ! unfilled slot (strict <, ties -> lowest index) gets min(1, remaining/(2j+1)) ... with ESH, slots with
        ! 2j+1 > remaining are skipped, and reaching a_max with such a slot is an error (dev throw, verbatim) ...
    subroutine HF_occupations_fill(Params,SPE,Orb,N_part,Kind,Symbol,O)
        type(params_t), intent(in) :: Params
        real(dp), intent(in) :: SPE(:)
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: N_part
        character(len=*), intent(in) :: Kind, Symbol
        real(dp), intent(out) :: O(:)
        logical, allocatable :: O_mask(:)
        logical :: ShellClosure, Found
        integer :: a_max, a, a_last
        real(dp) :: N_occ, E, E_last, j_hat, Occ

        a_max = a_max_of_calc(Params)
        ShellClosure = Params%Calc%HF%ESH
        allocate(O_mask(a_max))
        O = 0.0_dp
        O_mask = .false.

        N_occ = 0.0_dp
        do while (abs(N_part - N_occ) > 1.0e-8_dp)
            a_last = 1
            E_last = 1.0e5_dp
            Found = .false.
            do a = 1, a_max
                E = SPE(a)
                if ((E_last > E) .and. (.not. O_mask(a))) then
                    if (ShellClosure .and. ((N_part - N_occ) < real(Orb(a)%j + 1,dp))) then
                        if (a == a_max) call fatal('Expection: HF solver is not able to find closed-subshell solution for ' // &
                            Kind // ' sigle-particle orbitals with the given number of ' // Kind // 's ' // Symbol // ' = ' // &
                            jstr(N_part) // ' and the given number of single-particle orbitals a_max = ' // istr(a_max) // &
                            ' ... Please, try running the solver with the option ESH = .false. in &hf ...')
                        cycle
                    end if
                    E_last = E
                    a_last = a
                    Found = .true.
                end if
            end do
            ! Robustness guard (port only): dev falls back to slot 1 (already filled) and overwrites its occupation ...
            if (.not. Found) call fatal('HF_orbitals_make: no ' // Kind // ' orbital left to occupy (' // Symbol // ' = ' // &
                                        jstr(N_part) // ', occupied ' // jstr(N_occ) // ') ...')

            j_hat = real(Orb(a_last)%j + 1,dp)
            O_mask(a_last) = .true.
            Occ = min(1.0_dp, (N_part - N_occ) / j_hat)
            O(a_last) = Occ
            N_occ = N_occ + j_hat * Occ
        end do
    end subroutine HF_occupations_fill

    ! Total HF energy ... NN + NO2B NNN + kinetic terms ... a, b & c pruned by l_max / j_max of the occupied orbitals,
        ! J range max(|j_a-j_b|,|j_d-j_e|)/2 ... min(j_a+j_b,j_d+j_e)/2 (dev; the skipped terms vanish) ...
    function HF_energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN) result(E_HF)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        real(dp) :: E_HF
        integer :: A_mass, N_max, N_2max, N_3max, a_max, a, b, d, e, c, f, J, P, P3B, l_max, j_max
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_d, l_d, j_d, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        real(dp) :: Sum, Hat, pRho_ad, nRho_ad, pRho_be, nRho_be, pRho_cf, nRho_cf, me1, me2, me3, me4, me5
        character(len=32) :: CMS

        A_mass = Params%Calc%A
        N_max = Params%Calc%Nmax
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        CMS = Params%Calc%CMS
        a_max = (N_max + 1)*(N_max + 2)/2

        call HF_occupied_lj_max(Orb,l_max,j_max)

        write(*,'(/,a)') 'Calculating the total HF ground-state energy ...'
        E_HF = 0.0_dp

        !$omp parallel do schedule(dynamic) reduction(+:E_HF) private(Sum,n_a,l_a,j_a,b,n_b,l_b,j_b,P,d,l_d,j_d,n_d, &
        !$omp pRho_ad,nRho_ad,e,n_e,l_e,j_e,pRho_be,nRho_be,J,Hat,c,n_c,l_c,j_c,P3B,f,n_f,l_f,j_f,pRho_cf,nRho_cf, &
        !$omp me1,me2,me3,me4,me5)
        do a = 1, a_max
            Sum = 0.0_dp
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            if (l_a > l_max .or. j_a > j_max) cycle
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                if (l_b > l_max .or. j_b > j_max) cycle
                P = mod(l_a + l_b, 2) + 1
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                    do d = 1, a_max
                        l_d = Orb(d)%l
                        j_d = Orb(d)%j
                        if (l_a /= l_d .or. j_a /= j_d) cycle
                        n_d = Orb(d)%n
                        pRho_ad = Rho%p(a,d)
                        nRho_ad = Rho%n(a,d)
                        do e = 1, a_max
                            n_e = Orb(e)%n
                            l_e = Orb(e)%l
                            j_e = Orb(e)%j
                            if (l_b /= l_e .or. j_b /= j_e .or. (2*(n_d + n_e) + l_d + l_e) > N_2max) cycle
                            pRho_be = Rho%p(b,e)
                            nRho_be = Rho%n(b,e)
                            do J = max(abs(j_a - j_b)/2, abs(j_d - j_e)/2), min((j_a + j_b)/2, (j_d + j_e)/2)
                                Hat = real(2*J + 1,dp)
                                ! 2-body NN interaction ...
                                if (P == mod(l_d + l_e, 2) + 1) then
                                    Sum = Sum + 0.5_dp * Hat * pRho_ad * pRho_be * O2b_pp(a,b,d,e,J,P,V_NN,Orb,Orb_NN)
                                    Sum = Sum + 0.5_dp * Hat * nRho_ad * nRho_be * O2b_nn(a,b,d,e,J,P,V_NN,Orb,Orb_NN)
                                    Sum = Sum + Hat * pRho_ad * nRho_be * O2b_pn(a,b,d,e,J,P,V_NN,Orb_NN)
                                end if

                                ! 3-body NNN interaction ...
                                do c = 1, a_max
                                    n_c = Orb(c)%n
                                    l_c = Orb(c)%l
                                    if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max .or. l_c > l_max) cycle
                                    j_c = Orb(c)%j
                                    if (j_c > j_max) cycle
                                    P3B = mod(l_a + l_b + l_c, 2) + 1
                                    do f = 1, a_max
                                        n_f = Orb(f)%n
                                        l_f = Orb(f)%l
                                        j_f = Orb(f)%j
                                        if ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) > N_3max .or. l_c /= l_f .or. &
                                            j_c /= j_f .or. P3B /= (mod(l_d + l_e + l_f, 2) + 1)) cycle
                                        pRho_cf = Rho%p(c,f)
                                        nRho_cf = Rho%n(c,f)
                                        me1 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                        me2 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                        me3 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                        me4 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                        me5 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)
                                        Sum = Sum + 1.0_dp/6.0_dp * (me5 * pRho_ad * pRho_be * pRho_cf + &
                                              (2.0_dp * me4 + me5) * pRho_ad * pRho_be * nRho_cf + &
                                              (1.5_dp * me1 + sqrt(3.0_dp/4.0_dp) * me2 + sqrt(3.0_dp/4.0_dp) * me3 + &
                                              0.5_dp * me4 + me5) * pRho_ad * nRho_be * nRho_cf + &
                                              me5 * nRho_ad * nRho_be * nRho_cf)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
                if (CMS == 'CMS1+2B') then
                    Sum = Sum + (Rho%p(a,b) + Rho%n(a,b)) * T%p(a,b) * real(j_a + 1,dp) * (1.0_dp - 1.0_dp / real(A_mass,dp))
                else if (CMS /= 'CMS2B') then
                    Sum = Sum + (Rho%p(a,b) + Rho%n(a,b)) * T%p(a,b) * real(j_a + 1,dp)
                end if
            end do
            E_HF = E_HF + Sum
        end do
        !$omp end parallel do

        write(*,'(a)') TAB // 'HF energy    ...   E_HF = ' // jstr(E_HF) // ' MeV'
    end function HF_energy

    ! io/<Path>/HF/HF_summary.dat ... dev HF_summary ...
    subroutine HF_summary(Params,E_HF,T_HF,Iteration)
        type(params_t), intent(in) :: Params
        real(dp), intent(in) :: E_HF, T_HF
        integer, intent(in) :: Iteration
        integer :: u
        real(dp) :: A

        A = real(Params%Calc%A,dp)
        write(*,'(/,a)') 'Exporting HF calculation summary ...'
        open(newunit=u, file=io_path(Params,'/HF/HF_summary.dat'), position='append', action='write')
        call summary_header(u,Params,Params%Calc%A,Params%Calc%Z)
        write(u,'(/,a)') 'Spherical Hartree-Fock solution review:'
        write(u,'(/,a)') 'Number of Iterations = ' // cfmt_d(Iteration,5) // ', Precision = ' // &
                         cfmt_e(Params%Calc%HF%Tol,0,3) // ' ' // TAB // ' MeV'
        write(u,'(/,a)') 'E_HF = ' // cfmt_f(E_HF,15,8) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' Mean-field ground state energy'
        write(u,'(a)') 'T_HF = ' // cfmt_f(T_HF,15,8) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' Mean-field kinetic energy'
        write(u,'(/,a)') 'E_HF / A = ' // cfmt_f(E_HF / A,15,8) // TAB // ' MeV / Nucleon ' // TAB // TAB // ' ... ' // TAB // &
                         ' Mean-field ground state energy per nucleon'
        write(u,'(a)') 'T_HF / A = ' // cfmt_f(T_HF / A,15,8) // TAB // ' MeV / Nucleon ' // TAB // TAB // ' ... ' // TAB // &
                       ' Mean-field total kinetic energy per nucleon'
        close(u)
    end subroutine HF_summary

    subroutine summary_header(u,Params,A,Z)
        integer, intent(in) :: u, A, Z
        type(params_t), intent(in) :: Params
        integer :: N_max
        N_max = Params%Calc%Nmax
        write(u,'(a)') 'Nuclid data:    A = ' // istr(A) // ', Z = ' // istr(Z)
        write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                         ' , N_2max = ' // istr(Params%Calc%N2max) // ' , N_3max = ' // istr(Params%Calc%N3max) // &
                         ', J-basis size = ' // istr((N_max+1)*(N_max+2)/2) // ', M-basis size = ' // &
                         istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Combined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Only pure 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else
            write(u,'(/,a)') TAB // TAB // TAB // 'No Center-of-Mass System (CMS) motion correction is set ...'
        end if
    end subroutine summary_header

    ! Energy-ordered table of the HF single-particle states with the real HF occupations O_a ...
        ! Port deviation 4: the radial counters are sized by max(N2max+1, 2Nmax+1) in 2j ...
    subroutine HF_SPS_summary(Params,h,Orb)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: h
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), allocatable :: pSPE(:), nSPE(:)
        integer :: N_max, J_max, a_max, Ind_a, a, u, iso
        integer, allocatable :: Sort_Orbs(:), nr(:,:)
        real(dp) :: O_a, E_a

        N_max = Params%Calc%Nmax
        J_max = Params%Calc%N2max + 1
        a_max = (N_max + 1)*(N_max + 2)/2
        pSPE = diag(h%p)
        nSPE = diag(h%n)
        allocate(Sort_Orbs(a_max), nr(0:N_max,0:max(J_max,2*N_max+1)))

        write(*,'(/,a)') 'Exporting HF calculation summary ...'
        open(newunit=u, file=io_path(Params,'/HF/HF_summary.dat'), position='append', action='write')
        do iso = 1, 2
            nr = 0
            if (iso == 1) then
                write(u,'(/,a)') 'Proton single-particle states'
                call sortperm(pSPE,Sort_Orbs)
            else
                write(u,'(/,a)') 'Neutron single-particle states'
                call sortperm(nSPE,Sort_Orbs)
            end if
            write(u,'(a)') '____________________________________________________________________'
            write(u,'(a)') cfmt_s('n',4) // ' ' // cfmt_s('l',4) // ' ' // cfmt_s('2j',4) // ' ' // cfmt_s('O_a',8) // ' ' // &
                           cfmt_s('E',14)
            do Ind_a = 1, a_max
                a = Sort_Orbs(Ind_a)
                if (iso == 1) then
                    O_a = Orb(a)%pO
                    E_a = pSPE(a)
                else
                    O_a = Orb(a)%nO
                    E_a = nSPE(a)
                end if
                write(u,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // &
                               ' ' // cfmt_f(O_a,8,3) // ' ' // cfmt_f(E_a,14,5)
                nr(Orb(a)%l,Orb(a)%j) = nr(Orb(a)%l,Orb(a)%j) + 1
            end do
            write(u,'(a)') '____________________________________________________________________'
        end do
        close(u)
    end subroutine HF_SPS_summary

    ! io/<Path>/HF/bin/{orb,C,rho,H1b}.bin ... HF orbitals (final occupations), LHO -> HF coefficients, HF-basis
        ! density C^T Rho C & HF single-particle Hamiltonian diag(SPE) ...
    subroutine HF_export(Params,C,Rho,h,Orb)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: C, Rho, h
        type(orb1b_t), intent(in) :: Orb(:)
        call orbitals_export(Params,Orb,io_path(Params,'/HF/bin/orb.bin'))
        call O1b_export(Params,Orb,C,io_path(Params,'/HF/bin/C.bin'))
        call O1b_export(Params,Orb,Rho,io_path(Params,'/HF/bin/rho.bin'))
        call O1b_export(Params,Orb,h,io_path(Params,'/HF/bin/H1b.bin'))
    end subroutine HF_export

    ! Non-local HF radial potential V(r,r') = h_HF - T_1b per partial wave (l,j) & in total ... dev HF_potential.jl ...
        ! io/<Path>/HF/densities/HF_radial_potential.dat ... enabled by HF.V_HF ...
    subroutine HF_radial_potential(Params,Orb,C,T,h_N)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: C, T, h_N
        real(dp), parameter :: hc = 197.326980_dp, m_p = 938.272013_dp, m_n = 939.565346_dp
        real(dp), parameter :: Pi = 3.14159265358979323846264338327950288_dp
        character(len=*), parameter :: l_label = 'spdfghiklmnoqrtuvwxyz'
        real(dp), allocatable :: pV(:,:), nV(:,:), r_full(:), r_grid(:), pR(:,:), nR(:,:)
        real(dp), allocatable :: pV_lj(:,:,:), nV_lj(:,:,:), pV_tot(:,:), nV_tot(:,:)
        integer, allocatable :: lj_list(:,:), Block(:)
        character(len=16), allocatable :: lj_lab(:)
        character(len=:), allocatable :: Export_Path, Line
        real(dp) :: A_mass, hw, nu_p, nu_n
        integer :: a_max, N_r, N_lj, i, k, a, i_lj, l, j, u, N_B
        logical :: New
        type :: str_t
            character(len=:), allocatable :: s
        end type str_t
        type(str_t), allocatable :: Blocks(:)
        integer, parameter :: Chunk = 64
        integer :: i0, i1

        A_mass = real(Params%Calc%A,dp)
        hw = Params%Calc%hw
        a_max = a_max_of_calc(Params)

        write(*,'(/,a)') 'Preparing non-local radial HF mean-field potentials ...'

        ! HF potential in the LHO basis ... remove the 1-body kinetic term contained in h_N (leading block of T) ...
        pV = mmt(mm(C%p,h_N%p),C%p)
        nV = mmt(mm(C%n,h_N%n),C%n)
        if (Params%Calc%CMS == 'CMS1+2B') then
            pV = pV - T%p(1:a_max,1:a_max) * (1.0_dp - 1.0_dp / A_mass)
            nV = nV - T%n(1:a_max,1:a_max) * (1.0_dp - 1.0_dp / A_mass)
        else if (Params%Calc%CMS /= 'CMS2B') then
            pV = pV - T%p(1:a_max,1:a_max)
            nV = nV - T%n(1:a_max,1:a_max)
        end if

        ! Radial grid ... every 16th point of the rho_export grid ...
        r_full = radial_grid(0.0_dp + 1.0e-8_dp, 2.5_dp * 1.2_dp * A_mass**(1.0_dp/3.0_dp), 2**13 + 1)
        r_grid = r_full(1:size(r_full):16)
        N_r = size(r_grid)

        ! Radial LHO wave functions R_a(r) on the grid ...
        nu_p = 0.5_dp * m_p * hw / hc**2
        nu_n = 0.5_dp * m_n * hw / hc**2
        allocate(pR(a_max,N_r), nR(a_max,N_r))
        !$omp parallel do schedule(static) private(a)
        do i = 1, N_r
            do a = 1, a_max
                pR(a,i) = Psi_rad_LHO(r_grid(i),Orb(a)%n,Orb(a)%l,nu_p)
                nR(a,i) = Psi_rad_LHO(r_grid(i),Orb(a)%n,Orb(a)%l,nu_n)
            end do
        end do
        !$omp end parallel do

        ! Partial waves (l,j) of the basis, sorted by l, then 2j ...
        allocate(lj_list(2,a_max))
        N_lj = 0
        do a = 1, a_max
            New = .true.
            do i_lj = 1, N_lj
                if (lj_list(1,i_lj) == Orb(a)%l .and. lj_list(2,i_lj) == Orb(a)%j) New = .false.
            end do
            if (New) then
                N_lj = N_lj + 1
                lj_list(:,N_lj) = [Orb(a)%l, Orb(a)%j]
            end if
        end do
        call sort_lj(lj_list(:,1:N_lj))

        ! V_lj(r,r') = R_lj^T V_lj R_lj & V_tot = sum_lj (2j+1)/(4 pi) V_lj ...
        allocate(pV_lj(N_r,N_r,N_lj), nV_lj(N_r,N_r,N_lj), pV_tot(N_r,N_r), nV_tot(N_r,N_r), lj_lab(N_lj), Block(a_max))
        pV_tot = 0.0_dp
        nV_tot = 0.0_dp
        do i_lj = 1, N_lj
            l = lj_list(1,i_lj)
            j = lj_list(2,i_lj)
            N_B = 0
            do a = 1, a_max
                if (Orb(a)%l == l .and. Orb(a)%j == j) then
                    N_B = N_B + 1
                    Block(N_B) = a
                end if
            end do
            pV_lj(:,:,i_lj) = mm(mtm(pR(Block(1:N_B),:),pV(Block(1:N_B),Block(1:N_B))),pR(Block(1:N_B),:))
            nV_lj(:,:,i_lj) = mm(mtm(nR(Block(1:N_B),:),nV(Block(1:N_B),Block(1:N_B))),nR(Block(1:N_B),:))
            pV_tot = pV_tot + real(j + 1,dp) / (4.0_dp * Pi) * pV_lj(:,:,i_lj)
            nV_tot = nV_tot + real(j + 1,dp) / (4.0_dp * Pi) * nV_lj(:,:,i_lj)
            lj_lab(i_lj) = l_label(l+1:l+1) // istr(j) // '/2'
        end do

        ! Export ... columns r, r', pV_tot, nV_tot, then pV & nV of every partial wave, blank line after each r block ...
        Export_Path = io_path(Params,'/HF/densities/HF_radial_potential.dat')
        open(newunit=u, file=Export_Path, status='replace', action='write')
        write(u,'(a)') "# Non-local HF potential V(r,r') = h_HF - T_1b (CMS = " // trim(Params%Calc%CMS) // &
                       ") ... r, r' in fm, V in MeV fm^-3"
        write(u,'(a)') "# V_lj(r,r') = sum_ab R_a(r) V_ab R_b(r') per partial wave, V_tot = sum_lj (2j+1)/(4 pi) V_lj"
        Line = cfmt_s('r',14,.true.) // ' ' // cfmt_s("r'",14,.true.) // ' ' // cfmt_s('pV_tot',16,.true.) // ' ' // &
               cfmt_s('nV_tot',16,.true.)
        do i_lj = 1, N_lj
            Line = Line // ' ' // cfmt_s('pV_' // trim(lj_lab(i_lj)),16,.true.) // ' ' // cfmt_s('nV_' // trim(lj_lab(i_lj)),16,.true.)
        end do
        write(u,'(a)') Line
        close(u)
        ! Data rows ... the r blocks are formatted in parallel (chunks of r values) & appended as raw bytes in order ...
        open(newunit=u, file=Export_Path, access='stream', form='unformatted', position='append', action='write')
        allocate(Blocks(N_r))
        do i0 = 1, N_r, Chunk
            i1 = min(N_r, i0 + Chunk - 1)
            !$omp parallel do schedule(dynamic)
            do i = i0, i1
                Blocks(i)%s = format_block(i)
            end do
            !$omp end parallel do
            do i = i0, i1
                write(u) Blocks(i)%s
                deallocate(Blocks(i)%s)
            end do
        end do
        close(u)

        write(*,'(a)') TAB // 'Non-local radial HF mean-field potentials exported to ' // Export_Path // ' ...'

    contains

        ! Rows of the r = r_grid(i) block: "%-14.8f %-14.8f %-16.8e %-16.8e" + " %-16.8e %-16.8e" per wave, blank line after ...
        function format_block(i) result(B)
            integer, intent(in) :: i
            character(len=:), allocatable :: B
            character(len=:), allocatable :: Buf
            character(len=32*(2 + 2*N_lj)) :: Raw
            real(dp) :: X(2 + 2*N_lj)
            integer :: k, w, Pos, L, N_f, f, W_f, Off
            N_f = 2 + 2*N_lj
            L = 29 + 17*N_f + 1                             ! "%-14.8f %-14.8f" + N_f x " %-16.8e" + newline
            allocate(character(len=N_r*L + 1) :: Buf)
            Pos = 0
            do k = 1, N_r
                X(1) = pV_tot(i,k)
                X(2) = nV_tot(i,k)
                do w = 1, N_lj
                    X(1 + 2*w) = pV_lj(i,k,w)
                    X(2 + 2*w) = nV_lj(i,k,w)
                end do
                ! One internal write per row (fast), then C-like left justification & lower-case exponent per field ...
                write(Raw,'(2F16.8,*(ES17.8E2))') r_grid(i), r_grid(k), X(1:N_f)
                Buf(Pos+1:Pos+L) = ''
                do f = 1, 2 + N_f
                    if (f <= 2) then
                        W_f = 14
                        Off = (f - 1)*15
                        call put_field(Raw(16*(f-1)+1:16*f), W_f, Off, .false., r_grid(merge(i,k,f == 1)), Buf(Pos+1:Pos+L))
                    else
                        W_f = 16
                        Off = 30 + (f - 3)*17
                        call put_field(Raw(32 + 17*(f-3) + 1:32 + 17*(f-2)), W_f, Off, .true., X(f-2), Buf(Pos+1:Pos+L))
                    end if
                end do
                Buf(Pos+L:Pos+L) = achar(10)
                Pos = Pos + L
            end do
            Buf(Pos+1:Pos+1) = achar(10)
            B = Buf(1:Pos+1)
        end function format_block

        ! Left-justified field of width W_f at offset Off of Row ... E-format fields get a lower-case 'e'; a field the
            ! fast path cannot represent (3-digit exponent, overflow) falls back to cfmt_e ...
        subroutine put_field(Raw_f, W_f, Off, Is_E, x, Row)
            character(len=*), intent(in) :: Raw_f
            integer, intent(in) :: W_f, Off
            logical, intent(in) :: Is_E
            real(dp), intent(in) :: x
            character(len=*), intent(inout) :: Row
            character(len=:), allocatable :: Fld
            integer :: Pe
            Fld = trim(adjustl(Raw_f))
            if (Is_E) then
                Pe = index(Fld,'E')
                if (Pe == 0 .or. index(Fld,'*') > 0) then
                    Fld = trim(cfmt_e(x,W_f,8,.true.))
                else
                    Fld(Pe:Pe) = 'e'
                end if
            else
                if (Fld(1:1) == '.') Fld = '0' // Fld
                if (Fld(1:min(2,len(Fld))) == '-.') Fld = '-0' // Fld(2:)
                if (index(Fld,'*') > 0) Fld = trim(cfmt_f(x,W_f,8,.true.))
            end if
            Row(Off+1:Off+len(Fld)) = Fld
        end subroutine put_field

        ! Julia sort of (l,j) tuples ... lexicographic, insertion sort (few entries) ...
        subroutine sort_lj(X)
            integer, intent(inout) :: X(:,:)
            integer :: p, q, Tmp(2)
            do p = 2, size(X,2)
                Tmp = X(:,p)
                q = p - 1
                do while (q >= 1)
                    if (X(1,q) < Tmp(1) .or. (X(1,q) == Tmp(1) .and. X(2,q) <= Tmp(2))) exit
                    X(:,q+1) = X(:,q)
                    q = q - 1
                end do
                X(:,q+1) = Tmp
            end do
        end subroutine sort_lj
    end subroutine HF_radial_potential

end module nj_hf
