! LEGACY (removed in cleanup stage) ... the base-version spherical HF solver (old output layout io/<Path>/Bin/*,
    ! integer-like occupations), kept for the not-yet-ported solvers (BCS, HFB, HF_RPA, HF_RRPA, HF-EFA, EMPM, NCSM) ...
    ! counterpart of base Source/HF/HF_Solver.jl, HF_Allocate.jl, HF_Density_Operator.jl, HF_Ordering.jl, HF_Energy.jl &
    ! HF_Export.jl ... the dev HF is nj_hf ...
module nj_hf_legacy
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_o1b
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_v3b, only : V3b_no2b
    implicit none
    private

    public :: HF_solve, HF_allocate, HF_density_operator, HF_orbital_ordering, HF_energy
    public :: HF_summary, HF_SPS_summary, HF_export

contains

    ! Self-consistent iteration of the spherical HF equations ...
    subroutine HF_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,h,C,Rho,Iteration)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
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

        a_max = (Params%Calc%Nmax + 1) * (Params%Calc%Nmax + 2) / 2
        epsilon = Params%Calc%HF%Tol
        Iteration_max = Params%Calc%HF%Imax

        dE = 1.0_dp
        Iteration = 0
        allocate(SPE%p(a_max), SPE%n(a_max), SPE_old%p(a_max), SPE_old%n(a_max))
        SPE_old%p = 0.0_dp
        SPE_old%n = 0.0_dp

        ! Initial guess ... LHO orbitals ...
        C%p = identity(a_max)
        C%n = identity(a_max)
        Rho = HF_density_operator(a_max,C,Orb)

        write(*,'(/,a)') 'Initializing iteration of HF equations ...'
        call timer_start(t0)
        do while ((Iteration < Iteration_max) .and. (dE > epsilon))
            H_HF = HF_allocate(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)

            call eig_sym(H_HF%p,SPE%p,C%p)
            call eig_sym(H_HF%n,SPE%n,C%n)

            call HF_orbital_ordering(Orb,a_max,C,SPE)

            Rho = HF_density_operator(a_max,C,Orb)

            dE = (sum(abs(SPE%p - SPE_old%p)) + sum(abs(SPE%n - SPE_old%n))) / real(2 * a_max,dp)
            SPE_old%p = SPE%p
            SPE_old%n = SPE%n
            Iteration = Iteration + 1

            write(*,'(a)') TAB // 'HF iteration number: ' // cfmt_d(Iteration,4) // '   dE = ' // cfmt_f(dE,12,9) // ' MeV'
        end do
        call timer_report(t0)

        write(*,'(/,a)') 'Iteration of HF equations with NO2B NN+NNN interaction has finished ...'

        h%p = diagm(SPE%p)
        h%n = diagm(SPE%n)
    end subroutine HF_solve

    ! HF single-particle Hamiltonian for given density Rho ... NN + NO2B NNN + kinetic terms ...
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

        call HF_allocate_indices(Params,Orb,ad,ad_count,be,be_count)
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
                if ((2*(n_a + n_b) + l_a + l_b) > N_2max) cycle
                if (.not. (l_b == l_e .and. j_b == j_e .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max)) cycle
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
                        if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                        P_abc = mod(l_a + l_b + l_c, 2) + 1
                        j_c = Orb(c)%j
                        do f = 1, a_max
                            n_f = Orb(f)%n
                            l_f = Orb(f)%l
                            if ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                P_abc == (mod(l_d + l_e + l_f,2) + 1)) then
                                j_f = Orb(f)%j
                                if (j_c == j_f) then
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
                                end if
                            end if
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

    subroutine HF_allocate_indices(Params,Orb,ab,ab_count,de,de_count)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable, intent(out) :: ab(:,:), de(:,:)
        integer, intent(out) :: ab_count, de_count
        integer :: a_max, a, b, d, e, Pass

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
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
    end subroutine HF_allocate_indices

    function HF_density_operator(a_max,C,Orb) result(Rho)
        integer, intent(in) :: a_max
        type(o1b_t), intent(in) :: C
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t) :: Rho
        integer :: k, l, m
        real(dp) :: pSum, nSum

        call o1b_alloc(Rho,a_max)
        do k = 1, a_max
            do l = 1, a_max
                if (Orb(k)%l == Orb(l)%l .and. Orb(k)%j == Orb(l)%j) then
                    pSum = 0.0_dp
                    nSum = 0.0_dp
                    do m = 1, a_max
                        if ((Orb(m)%pO == 1) .and. Orb(m)%l == Orb(k)%l .and. Orb(m)%j == Orb(k)%j) &
                            pSum = pSum + C%p(k,m) * C%p(l,m)
                        if ((Orb(m)%nO == 1) .and. Orb(m)%l == Orb(k)%l .and. Orb(m)%j == Orb(k)%j) &
                            nSum = nSum + C%n(k,m) * C%n(l,m)
                    end do
                    Rho%p(k,l) = pSum
                    Rho%n(k,l) = nSum
                end if
            end do
        end do
    end function HF_density_operator

    ! Assigns each HF orbital to an LHO orbital slot with matching l & 2j (in energy order) ...
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

    ! Total HF energy ... NN + NO2B NNN + kinetic terms ...
    function HF_energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN) result(E_HF)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        real(dp) :: E_HF
        integer :: A_mass, N_max, N_2max, N_3max, a_max, a, b, d, e, c, f, J, P, P3B
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_d, l_d, j_d, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        real(dp) :: Sum, Hat, pRho_ad, nRho_ad, pRho_be, nRho_be, pRho_cf, nRho_cf, me1, me2, me3, me4, me5
        character(len=32) :: CMS

        A_mass = Params%Calc%A
        N_max = Params%Calc%Nmax
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        CMS = Params%Calc%CMS
        a_max = (N_max + 1)*(N_max + 2)/2

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
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                P = mod(l_a + l_b, 2) + 1
                j_b = Orb(b)%j
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                    do d = 1, a_max
                        l_d = Orb(d)%l
                        j_d = Orb(d)%j
                        if ((l_a == l_d) .and. (j_a == j_d)) then
                            n_d = Orb(d)%n
                            pRho_ad = Rho%p(a,d)
                            nRho_ad = Rho%n(a,d)
                            do e = 1, a_max
                                n_e = Orb(e)%n
                                l_e = Orb(e)%l
                                j_e = Orb(e)%j
                                if ((l_b == l_e .and. j_b == j_e) .and. ((2*(n_d + n_e) + l_d + l_e) <= N_2max)) then
                                    pRho_be = Rho%p(b,e)
                                    nRho_be = Rho%n(b,e)
                                    do J = abs(j_a - j_b)/2, (j_a + j_b)/2
                                        ! 2-body NN interaction ...
                                        if (mod(l_a + l_b, 2) == mod(l_d + l_e, 2)) then
                                            Hat = real(2*J + 1,dp)
                                            Sum = Sum + 0.5_dp * Hat * pRho_ad * pRho_be * O2b_pp(a,b,d,e,J,P,V_NN,Orb,Orb_NN)
                                            Sum = Sum + 0.5_dp * Hat * nRho_ad * nRho_be * O2b_nn(a,b,d,e,J,P,V_NN,Orb,Orb_NN)
                                            Sum = Sum + Hat * pRho_ad * nRho_be * O2b_pn(a,b,d,e,J,P,V_NN,Orb_NN)
                                        end if

                                        ! 3-body NNN interaction ...
                                        do c = 1, a_max
                                            n_c = Orb(c)%n
                                            l_c = Orb(c)%l
                                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max) then
                                                j_c = Orb(c)%j
                                                P3B = mod(l_a + l_b + l_c, 2) + 1
                                                do f = 1, a_max
                                                    n_f = Orb(f)%n
                                                    l_f = Orb(f)%l
                                                    j_f = Orb(f)%j
                                                    if ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                                        j_c == j_f .and. P3B == (mod(l_d + l_e + l_f, 2) + 1)) then
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
                                                    end if
                                                end do
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end if
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

    subroutine HF_summary(Params,E_HF,T_HF,Iteration)
        type(params_t), intent(in) :: Params
        real(dp), intent(in) :: E_HF, T_HF
        integer, intent(in) :: Iteration
        integer :: u, N_max
        real(dp) :: A

        A = real(Params%Calc%A,dp)
        N_max = Params%Calc%Nmax
        write(*,'(/,a)') 'Exporting HF calculation summary ...'
        open(newunit=u, file=io_path(Params,'/HF/HF_Summary.dat'), position='append', action='write')
        call summary_header(u,Params,Params%Calc%A,Params%Calc%Z)
        write(u,'(/,a)') 'Spherical Hartree-Fock solution review:'
        write(u,'(/,a)',advance='no') 'Number of Iterations = ' // cfmt_d(Iteration,5) // ', Precision = ' // &
                                      cfmt_e(Params%Calc%HF%Tol,0,3) // ' ' // TAB // ' MeV'
        write(u,'(/,a)') 'E_HF = ' // cfmt_f(E_HF,15,8) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' Mean-field ground state energy'
        write(u,'(a)') 'T_HF = ' // cfmt_f(T_HF,15,8) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' Mean-field kinetic energy'
        write(u,'(/,a)') 'E_HF / A = ' // cfmt_f(E_HF / A,15,8) // TAB // ' MeV / Nucleon ' // TAB // TAB // ' ... ' // TAB // &
                         ' Mean-field ground state energy per nucleon'
        write(u,'(a)') 'T_HF / A = ' // cfmt_f(T_HF / A,15,8) // TAB // ' MeV / Nucleon ' // TAB // TAB // ' ... ' // TAB // &
                       ' Mean-field total kinetic energy per nucleon'
        close(u)
    end subroutine HF_summary

    ! Common header of the HF & HF-EFA summary files ...
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

    subroutine HF_SPS_summary(Params,h,Orb)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: h
        type(orb1b_t), intent(in) :: Orb(:)
        call sps_table(Params,diag(h%p),diag(h%n),Orb,io_path(Params,'/HF/HF_Summary.dat'),'HF')
    end subroutine HF_SPS_summary

    subroutine HF_export(Params,C,Rho,h,Orb)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: C, Rho, h
        type(orb1b_t), intent(in) :: Orb(:)
        call O1b_export(Params,Orb,C,io_path(Params,'/Bin/C_HF.bin'))
        call O1b_export(Params,Orb,Rho,io_path(Params,'/Bin/Rho_HF.bin'))
        call O1b_export(Params,Orb,h,io_path(Params,'/Bin/h_HF.bin'))
    end subroutine HF_export

    ! Energy-ordered table of single-particle states ... occupations O_a from the naive filling ...
    subroutine sps_table(Params,pSPE,nSPE,Orb,File,Kind)
        type(params_t), intent(in) :: Params
        real(dp), intent(in) :: pSPE(:), nSPE(:)
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: File, Kind
        integer :: N_max, J_max, a_max, Ind_a, a, u, iso
        integer, allocatable :: Sort_Orbs(:), nr(:,:)
        real(dp) :: O_a

        N_max = Params%Calc%Nmax
        J_max = Params%Calc%N2max + 1
        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(Sort_Orbs(a_max), nr(0:N_max,0:max(J_max,2*N_max+1)))

        write(*,'(/,a)') 'Exporting ' // Kind // ' calculation summary ...'
        open(newunit=u, file=File, position='append', action='write')
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
                    O_a = merge(1.0_dp, 0.0_dp, Orb(a)%pO == 1)
                    write(u,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // &
                                   ' ' // cfmt_f(O_a,8,3) // ' ' // cfmt_f(pSPE(a),14,5)
                else
                    O_a = merge(1.0_dp, 0.0_dp, Orb(a)%nO == 1)
                    write(u,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // &
                                   ' ' // cfmt_f(O_a,8,3) // ' ' // cfmt_f(nSPE(a),14,5)
                end if
                nr(Orb(a)%l,Orb(a)%j) = nr(Orb(a)%l,Orb(a)%j) + 1
            end do
            write(u,'(a)') '____________________________________________________________________'
        end do
        close(u)
    end subroutine sps_table

end module nj_hf_legacy
