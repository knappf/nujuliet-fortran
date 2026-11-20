! Single-particle orbitals ... counterpart of src/common/orb.jl (NuJuliet dev; MatrixElements/Orb.jl before) ...
module nj_orbitals
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg, only : eig_sym, hungarian
    implicit none
    private

    public :: orbitals_make, orbitals_lj_make, orbitals_export, orbitals_import, orbitals_update
    public :: orbitals_1p1h_make, orb2b_1p1h_index, orbitals_2qp_make, orb2b_2qp_index, OBDM_natural_orbitals, a_max_of
    ! LEGACY (removed in cleanup stage) ... dead in dev ("Old remove soon"), still used by the old TDA/RPA/RRPA, HF-MBPT & EMPM ...
    public :: orbitals_ph_make, orbitals_ph_list, orbitals_one_phonon_make

    ! orbitals_export(Params,Orb,Export_Path) ... binary orbital table (dev) ... orbitals_export(Params,Orb) ... LEGACY
        ! human-readable tables of the HRF switch (no caller in dev, removed in cleanup stage) ...
    interface orbitals_export
        module procedure orbitals_export_bin, orbitals_export_hrf
    end interface orbitals_export

contains

    pure integer function a_max_of(N_max)
        integer, intent(in) :: N_max
        a_max_of = (N_max + 1)*(N_max + 2)/2
    end function a_max_of

    ! LHO single-particle orbitals in NuHamil ordering with the dev ground-state filling ... the orbitals of each isospin
        ! are filled in the order N, l = N, N-2, ..., j = l+1/2 first, each with min(1, remaining / (2j+1)), i.e. the first
        ! shell that does not fit gets the fractional (equal-filling) occupation ... Calc.Nmax orbitals, Int.Nmax ones with
        ! Use_Int = .true. (Julia keyword Int = true, used by the interaction readers) ...
        ! As in dev the orbitals after the last filled one get min(1, residual / (2j+1)) with a residual of ~1e-16 (not
        ! exactly 0), reproduced for bit parity (the arithmetic order is kept & j_hat*Occ is not contracted into an FMA) ...
        ! Deviation: dev loops forever if A-Z or Z exceed the capacity of the basis, the port stops with an error ...
    function orbitals_make(Params,Use_Int) result(Orb)
        type(params_t), intent(in) :: Params
        logical, intent(in), optional :: Use_Int
        type(orb1b_t), allocatable :: Orb(:)
        real(dp) :: A, Z, Z_occ, N_occ, Occ, j_hat
        real(dp), volatile :: Dj_occ
        integer :: N_max, a_max, a_i, N, l, n_r, j, b
        integer, allocatable :: Orb_arr(:,:)
        real(dp), allocatable :: Obc_occ(:,:)
        logical, allocatable :: Z_mask(:), N_mask(:)
        logical :: Z_fill, N_fill

        A = real(Params%Calc%A,dp)
        Z = real(Params%Calc%Z,dp)
        N_max = Params%Calc%Nmax
        if (present(Use_Int)) then
            if (Use_Int) N_max = Params%Int%Nmax
        end if
        a_max = a_max_of(N_max)

        write(*,'(/,a)') 'Constructing Linear-Harmonic Oscillator basis single-particle orbitals ...'
        allocate(Orb(a_max), Orb_arr(a_max,4), Obc_occ(a_max,2))
        Orb_arr = 0
        Obc_occ = 0.0_dp
        a_i = 0
        do N = 0, N_max
            do l = 0, N
                n_r = (N - l) / 2
                if (2*n_r + l == N) then
                    do j = abs(2*l - 1), 2*l + 1, 2
                        a_i = a_i + 1
                        Orb_arr(a_i,1) = a_i
                        Orb_arr(a_i,2) = n_r
                        Orb_arr(a_i,3) = l
                        Orb_arr(a_i,4) = j
                    end do
                end if
            end do
        end do

        if (abs(A) > 1.0e-8_dp .or. abs(Z) > 1.0e-8_dp) then
            Z_occ = 0.0_dp
            N_occ = 0.0_dp
            Z_fill = .true.
            N_fill = .true.
            allocate(Z_mask(a_max), N_mask(a_max))
            Z_mask = .false.
            N_mask = .false.

            ! Neutrons ... one pass visits every orbital once (mask), a second pass of the dev while loop would change
                ! nothing ...
            do N = 0, N_max
                do n_r = 0, N/2
                    l = N - 2*n_r
                    do j = 2*l + 1, abs(2*l - 1), -2
                        j_hat = real(j + 1,dp)
                        do b = 1, a_max
                            if (Orb_arr(b,2) == n_r .and. Orb_arr(b,3) == l .and. Orb_arr(b,4) == j .and. .not. N_mask(b)) then
                                Occ = min(1.0_dp, ((A - Z) - N_occ) / j_hat)
                                Obc_occ(b,2) = Occ
                                Dj_occ = j_hat * Occ
                                N_occ = N_occ + Dj_occ
                                N_mask(b) = .true.
                                if (abs((A - Z) - N_occ) < 1.0e-8_dp) N_fill = .false.
                            end if
                        end do
                    end do
                end do
            end do
            if (N_fill) call fatal('orbitals_make: A - Z = ' // istr(Params%Calc%A - Params%Calc%Z) // &
                                   ' neutrons do not fit into the basis Nmax = ' // istr(N_max) // ' ...')

            ! Protons ...
            do N = 0, N_max
                do n_r = 0, N/2
                    l = N - 2*n_r
                    do j = 2*l + 1, abs(2*l - 1), -2
                        j_hat = real(j + 1,dp)
                        do b = 1, a_max
                            if (Orb_arr(b,2) == n_r .and. Orb_arr(b,3) == l .and. Orb_arr(b,4) == j .and. .not. Z_mask(b)) then
                                Occ = min(1.0_dp, (Z - Z_occ) / j_hat)
                                Obc_occ(b,1) = Occ
                                Dj_occ = j_hat * Occ
                                Z_occ = Z_occ + Dj_occ
                                Z_mask(b) = .true.
                                if (abs(Z - Z_occ) < 1.0e-8_dp) Z_fill = .false.
                            end if
                        end do
                    end do
                end do
            end do
            if (Z_fill) call fatal('orbitals_make: Z = ' // istr(Params%Calc%Z) // ' protons do not fit into the basis Nmax = ' // &
                                   istr(N_max) // ' ...')
        end if

        do b = 1, a_max
            Orb(b) = orb1b_t(Orb_arr(b,1),Orb_arr(b,2),Orb_arr(b,3),Orb_arr(b,4),Obc_occ(b,1),Obc_occ(b,2))
        end do
        write(*,'(/,a)') 'Linear-Harmonic Oscillator basis single-particle orbitals initialized ...'
    end function orbitals_make

    ! Occupations from the diagonal of a density matrix ... pO_a = Rho%p(a,a), nO_a = Rho%n(a,a) ... counterpart of
        ! orbitals_update (which mutates its argument) ...
    subroutine orbitals_update(Orb,Rho)
        type(orb1b_t), intent(inout) :: Orb(:)
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), allocatable :: Orb_new(:)
        integer :: i, a
        allocate(Orb_new(size(Orb)))
        do i = 1, size(Orb)
            a = Orb(i)%a
            Orb_new(a) = orb1b_t(a,Orb(i)%n,Orb(i)%l,Orb(i)%j,Rho%p(a,a),Rho%n(a,a))
        end do
        Orb = Orb_new
    end subroutine orbitals_update

    ! Binary table of the single-particle orbitals ... a_max(Calc.Nmax) records Int64 a, n, l, j & Float64 pO, nO
        ! (48 bytes each, raw stream, no header) ... e.g. io/<Path>/HF/bin/orb.bin ...
    subroutine orbitals_export_bin(Params,Orb,Export_Path)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: Export_Path
        integer :: a_max, i, u
        a_max = a_max_of(Params%Calc%Nmax)
        if (size(Orb) < a_max) call fatal('orbitals_export: fewer orbitals than a_max(Calc.Nmax) ... ' // trim(Export_Path))
        open(newunit=u, file=Export_Path, access='stream', form='unformatted', status='replace', action='write')
        do i = 1, a_max
            write(u) int(Orb(i)%a,i8), int(Orb(i)%n,i8), int(Orb(i)%l,i8), int(Orb(i)%j,i8), Orb(i)%pO, Orb(i)%nO
        end do
        close(u)
    end subroutine orbitals_export_bin

    function orbitals_import(Params,Import_Path) result(Orb)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Import_Path
        type(orb1b_t), allocatable :: Orb(:)
        integer(i8) :: a, n, l, j
        real(dp) :: pO, nO
        integer :: a_max, i, u, Ios
        character(len=512) :: Msg
        a_max = a_max_of(Params%Calc%Nmax)
        allocate(Orb(a_max))
        open(newunit=u, file=Import_Path, access='stream', form='unformatted', status='old', action='read', &
             iostat=Ios, iomsg=Msg)
        if (Ios /= 0) call fatal('orbitals_import: cannot open ' // trim(Import_Path) // ' ... ' // trim(Msg))
        do i = 1, a_max
            read(u, iostat=Ios, iomsg=Msg) a, n, l, j, pO, nO
            if (Ios /= 0) call fatal('orbitals_import: cannot read ' // trim(Import_Path) // ' ... ' // trim(Msg))
            Orb(i) = orb1b_t(int(a),int(n),int(l),int(j),pO,nO)
        end do
        close(u)
    end function orbitals_import

    ! Orb_lj(l,j) lists orbitals of given l & 2j ... Julia Orb_lj[l+1,j] -> Orb_lj(l,j) here ...
    ! A subroutine (not a function) so that the (0:N_max, 1:2N_max+1) bounds reach the caller ...
    subroutine orbitals_lj_make(Params,Orb,Orb_lj)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(ivec_t), allocatable, intent(out) :: Orb_lj(:,:)
        integer :: N_max, a_max, l, j, a, lj_count
        N_max = Params%Calc%Nmax
        a_max = a_max_of(N_max)
        allocate(Orb_lj(0:N_max,1:2*N_max+1))
        do l = 0, N_max
            do j = 1, 2*N_max + 1
                lj_count = count(Orb(1:a_max)%l == l .and. Orb(1:a_max)%j == j)
                allocate(Orb_lj(l,j)%v(lj_count))
                lj_count = 0
                do a = 1, a_max
                    if (Orb(a)%l == l .and. Orb(a)%j == j) then
                        lj_count = lj_count + 1
                        Orb_lj(l,j)%v(lj_count) = a
                    end if
                end do
            end do
        end do
    end subroutine orbitals_lj_make

    ! LEGACY (removed in cleanup stage) ... HRF tables, no caller in dev ... occupations printed with jstr since they
        ! are Float64 now ...
    subroutine orbitals_export_hrf(Params,Orb)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer :: a_max, a, u
        real(dp), allocatable :: pSPE(:), nSPE(:)

        a_max = a_max_of(Params%Calc%Nmax)
        ! The original reads Bin/pE_HF.bin & Bin/nE_HF.bin, which no solver writes ... the HF single-particle
        ! energies are taken from the diagonal of Bin/h_HF.bin instead ...
        call read_h_HF_diag(Params,Orb,pSPE,nSPE)

        write(*,'(/,a)') 'Exporting LHO basis s.p. orbitals in human-readable format ...'
        open(newunit=u, file=io_path(Params,'/Orbitals_LHO.dat'), status='replace', action='write')
        write(u,'(a)') 'a' // TAB // 'n' // TAB // 'l' // TAB // '2j' // TAB // 'pOc' // TAB // 'nOc'
        do a = 1, a_max
            write(u,'(a)') istr(Orb(a)%a) // TAB // istr(Orb(a)%n) // TAB // istr(Orb(a)%l) // TAB // istr(Orb(a)%j) // &
                           TAB // ' ' // jstr(Orb(a)%pO) // TAB // ' ' // jstr(Orb(a)%nO)
        end do
        close(u)
        write(*,'(/,a)') 'Exported LHO basis s.p. orbitals in human-readable format ...'

        write(*,'(/,a)') 'Exporting HF basis s.p. orbitals in human-readable format ...'
        open(newunit=u, file=io_path(Params,'/Orbitals_HF.dat'), status='replace', action='write')
        write(u,'(a)') 'a' // TAB // 'l' // TAB // '2j' // TAB // 'pOc' // TAB // 'nOc' // TAB // 'E_p' // TAB // TAB // TAB // 'E_n'
        do a = 1, a_max
            write(u,'(a)') istr(Orb(a)%a) // TAB // istr(Orb(a)%l) // TAB // istr(Orb(a)%j) // TAB // ' ' // jstr(Orb(a)%pO) // &
                           TAB // ' ' // jstr(Orb(a)%nO) // TAB // jstr(pSPE(a)) // TAB // jstr(nSPE(a))
        end do
        close(u)
        write(*,'(/,a)') 'Exported HF basis s.p. orbitals in human-readable format ...'
    end subroutine orbitals_export_hrf

    ! Diagonal of the HF single-particle Hamiltonian stored in Bin/h_HF.bin ...
    subroutine read_h_HF_diag(Params,Orb,pSPE,nSPE)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), allocatable, intent(out) :: pSPE(:), nSPE(:)
        real(dp) :: ME
        integer :: a_max, a, b, u, iso

        a_max = a_max_of(Params%Calc%Nmax)
        allocate(pSPE(a_max), nSPE(a_max))
        open(newunit=u, file=io_path(Params,'/Bin/h_HF.bin'), access='stream', form='unformatted', status='old', action='read')
        do iso = 1, 2
            do a = 1, a_max
                do b = 1, a_max
                    if (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j) then
                        read(u) ME
                        if (a == b .and. iso == 1) pSPE(a) = ME
                        if (a == b .and. iso == 2) nSPE(a) = ME
                    end if
                end do
            end do
        end do
        close(u)
    end subroutine read_h_HF_diag

    ! LEGACY (removed in cleanup stage) ... particle & hole states with HF single-particle energies read from Bin/h_HF.bin ...
        ! holes are the orbitals with occupation 1, particles those with 0 (fractional occupations belong to neither) ...
    subroutine orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(pninteger_t), intent(out) :: N_Particle, N_Hole
        type(pnsvector_t), intent(out) :: Particle, Hole
        real(dp), allocatable :: pSPE(:), nSPE(:)
        integer :: a_max, a

        a_max = a_max_of(Params%Calc%Nmax)
        call read_h_HF_diag(Params,Orb,pSPE,nSPE)

        N_Hole%p = count(Orb(1:a_max)%pO == 1.0_dp)
        N_Particle%p = count(Orb(1:a_max)%pO == 0.0_dp)
        N_Hole%n = count(Orb(1:a_max)%nO == 1.0_dp)
        N_Particle%n = count(Orb(1:a_max)%nO == 0.0_dp)
        allocate(Particle%p(N_Particle%p), Particle%n(N_Particle%n), Hole%p(N_Hole%p), Hole%n(N_Hole%n))

        N_Particle = pninteger_t(0,0)
        N_Hole = pninteger_t(0,0)
        do a = 1, a_max
            if (Orb(a)%pO == 0.0_dp) then
                N_Particle%p = N_Particle%p + 1
                Particle%p(N_Particle%p) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,pSPE(a))
            end if
            if (Orb(a)%pO == 1.0_dp) then
                N_Hole%p = N_Hole%p + 1
                Hole%p(N_Hole%p) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,pSPE(a))
            end if
            if (Orb(a)%nO == 0.0_dp) then
                N_Particle%n = N_Particle%n + 1
                Particle%n(N_Particle%n) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,nSPE(a))
            end if
            if (Orb(a)%nO == 1.0_dp) then
                N_Hole%n = N_Hole%n + 1
                Hole%n(N_Hole%n) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,nSPE(a))
            end if
        end do
    end subroutine orbitals_ph_make

    ! LEGACY (removed in cleanup stage) ... particle-hole phonon indices ... ParticleHole%p(J,P,p,h) with J starting at 0 ...
    function orbitals_ph_list(J_max,N_Particle,Particle,N_Hole,Hole) result(ParticleHole)
        integer, intent(in) :: J_max
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(pnarray4_t) :: ParticleHole
        integer, allocatable :: N_Phonon(:,:)
        integer :: p, h, P_ph, J

        allocate(N_Phonon(0:J_max,2))
        N_Phonon = 0
        allocate(ParticleHole%p(0:J_max,2,N_Particle%p,N_Hole%p), ParticleHole%n(0:J_max,2,N_Particle%n,N_Hole%n))
        ParticleHole%p = 0
        ParticleHole%n = 0
        do p = 1, N_Particle%p
            do h = 1, N_Hole%p
                P_ph = mod(Particle%p(p)%l + Hole%p(h)%l,2) + 1
                do J = abs(Particle%p(p)%j - Hole%p(h)%j)/2, (Particle%p(p)%j + Hole%p(h)%j)/2
                    N_Phonon(J,P_ph) = N_Phonon(J,P_ph) + 1
                    ParticleHole%p(J,P_ph,p,h) = N_Phonon(J,P_ph)
                end do
            end do
        end do
        do p = 1, N_Particle%n
            do h = 1, N_Hole%n
                P_ph = mod(Particle%n(p)%l + Hole%n(h)%l,2) + 1
                do J = abs(Particle%n(p)%j - Hole%n(h)%j)/2, (Particle%n(p)%j + Hole%n(h)%j)/2
                    N_Phonon(J,P_ph) = N_Phonon(J,P_ph) + 1
                    ParticleHole%n(J,P_ph,p,h) = N_Phonon(J,P_ph)
                end do
            end do
        end do
    end function orbitals_ph_list

    ! LEGACY (removed in cleanup stage) ...
    subroutine orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole,N_Phonon,Phonon)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        integer, intent(out) :: N_Phonon
        type(phstate_t), allocatable, intent(out) :: Phonon(:)
        integer :: p, h, J, N_ph

        N_Phonon = 0
        do p = 1, N_Particle%p
            do h = 1, N_Hole%p
                N_Phonon = N_Phonon + (Particle%p(p)%j + Hole%p(h)%j)/2 - abs(Particle%p(p)%j - Hole%p(h)%j)/2 + 1
            end do
        end do
        do p = 1, N_Particle%n
            do h = 1, N_Hole%n
                N_Phonon = N_Phonon + (Particle%n(p)%j + Hole%n(h)%j)/2 - abs(Particle%n(p)%j - Hole%n(h)%j)/2 + 1
            end do
        end do

        allocate(Phonon(N_Phonon))
        N_ph = 0
        do p = 1, N_Particle%p
            do h = 1, N_Hole%p
                do J = abs(Particle%p(p)%j - Hole%p(h)%j)/2, (Particle%p(p)%j + Hole%p(h)%j)/2
                    N_ph = N_ph + 1
                    Phonon(N_ph) = phstate_t(p,h,-1,J,mod(Particle%p(p)%l + Hole%p(h)%l,2) + 1)
                end do
            end do
        end do
        do p = 1, N_Particle%n
            do h = 1, N_Hole%n
                do J = abs(Particle%n(p)%j - Hole%n(h)%j)/2, (Particle%n(p)%j + Hole%n(h)%j)/2
                    N_ph = N_ph + 1
                    Phonon(N_ph) = phstate_t(p,h,1,J,mod(Particle%n(p)%l + Hole%n(h)%l,2) + 1)
                end do
            end do
        end do
    end subroutine orbitals_one_phonon_make

    ! 1p-1h orbitals of TDA & RPA ... counterpart of orbitals_1p1h_make (orb.jl:348) ... a pair (a = particle, b = hole)
        ! of isospin T enters channel (J,P) if O_b > Tol_hol, O_a < Tol_par & O_b - O_a > 1e-8 ... order inside a channel:
        ! a outer, b inner (both over Orb), the proton pair (T = -1) before the neutron pair (T = +1) ... with fractional
        ! occupations an open shell with Tol_par <= O <= Tol_hol is neither particle nor hole (dropped, as in dev) ...
        ! The dev defaults of the keywords (Tol_par = 1e-4, Tol_hol = 0.999) are never used, all callers pass the
        ! tolerances of their solver ...
    function orbitals_1p1h_make(Params,Orb,Tol_par,Tol_hol) result(Orb_1p1h)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: Tol_par, Tol_hol
        type(orb2b_1p1h_t) :: Orb_1p1h
        integer :: J_max, a_max, J, P, ia, ib, a, b, Pass, N_ph_JP

        J_max = 2*Params%Calc%Nmax + 1
        a_max = maxval(Orb(:)%a)
        Orb_1p1h%a_max = a_max
        Orb_1p1h%J_max = J_max
        allocate(Orb_1p1h%N(2,0:J_max), Orb_1p1h%Ind(2,0:J_max), Orb_1p1h%Dic(a_max,a_max,0:J_max,2), &
                 Orb_1p1h%Par(a_max,a_max))
        Orb_1p1h%N = 0
        Orb_1p1h%Dic = 0
        do ia = 1, size(Orb)
            do ib = 1, size(Orb)
                Orb_1p1h%Par(Orb(ia)%a,Orb(ib)%a) = parity_of(Orb(ia)%l + Orb(ib)%l)
            end do
        end do

        ! Pass 1 counts, pass 2 fills (same order) ...
        do Pass = 1, 2
            do J = 0, J_max
                do P = 1, 2
                    N_ph_JP = 0
                    if (Pass == 2) allocate(Orb_1p1h%Ind(P,J)%v(Orb_1p1h%N(P,J)))
                    do ia = 1, size(Orb)
                        a = Orb(ia)%a
                        do ib = 1, size(Orb)
                            b = Orb(ib)%a
                            if (parity_of(Orb(ia)%l + Orb(ib)%l) /= P .or. .not. triangle(Orb(ia)%j,Orb(ib)%j,2*J)) cycle
                            if (Orb(ib)%pO > Tol_hol .and. Orb(ia)%pO < Tol_par .and. (Orb(ib)%pO - Orb(ia)%pO) > 1.0e-8_dp) then
                                N_ph_JP = N_ph_JP + 1
                                if (Pass == 2) then
                                    Orb_1p1h%Ind(P,J)%v(N_ph_JP) = ind2b_t(a,b,-1)
                                    Orb_1p1h%Dic(a,b,J,1) = N_ph_JP
                                end if
                            end if
                            if (Orb(ib)%nO > Tol_hol .and. Orb(ia)%nO < Tol_par .and. (Orb(ib)%nO - Orb(ia)%nO) > 1.0e-8_dp) then
                                N_ph_JP = N_ph_JP + 1
                                if (Pass == 2) then
                                    Orb_1p1h%Ind(P,J)%v(N_ph_JP) = ind2b_t(a,b,1)
                                    Orb_1p1h%Dic(a,b,J,2) = N_ph_JP
                                end if
                            end if
                        end do
                    end do
                    Orb_1p1h%N(P,J) = N_ph_JP
                end do
            end do
        end do
    end function orbitals_1p1h_make

    ! Index of the 1p-1h pair (a,b,J,P,T) in its channel ... 0 if absent (as get(Dic,key,0) in dev) ...
    pure integer function orb2b_1p1h_index(a,b,J,P,T,Orb_1p1h)
        integer, intent(in) :: a, b, J, P, T
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        orb2b_1p1h_index = 0
        if (a < 1 .or. a > Orb_1p1h%a_max .or. b < 1 .or. b > Orb_1p1h%a_max) return
        if (J < 0 .or. J > Orb_1p1h%J_max .or. (T /= -1 .and. T /= 1)) return
        if (Orb_1p1h%Par(a,b) /= P) return
        orb2b_1p1h_index = Orb_1p1h%Dic(a,b,J,(T + 3)/2)
    end function orb2b_1p1h_index

    ! 2-quasiparticle orbitals ... a >= b convention ... Orb_2qp%N(P,0:J_max) ... pp entries first, then nn entries ...
        ! counterpart of orbitals_2qp_make (orb.jl:465) incl. the dense lookup Dic (dev Dict of orb2b_2qp_key) ...
    function orbitals_2qp_make(Params,Orb) result(Orb_2qp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t) :: Orb_2qp
        integer :: N_max, a_max, J_max, a, b, J, P, T, Pass

        N_max = Params%Calc%Nmax
        a_max = a_max_of(N_max)
        J_max = 2*N_max + 1
        Orb_2qp%a_max = a_max
        Orb_2qp%J_max = J_max
        allocate(Orb_2qp%N(2,0:J_max), Orb_2qp%Ind(2,0:J_max), Orb_2qp%Dic(a_max,a_max,0:J_max,2))
        Orb_2qp%Dic = 0

        ! Pass 1 counts, pass 2 allocates the entries ... pp entries first, then nn entries ...
        do Pass = 1, 2
            Orb_2qp%N = 0
            do T = -1, 1, 2
                do a = 1, a_max
                    do b = 1, a
                        do J = abs(Orb(a)%j - Orb(b)%j)/2, (Orb(a)%j + Orb(b)%j)/2
                            if (a /= b .or. (a == b .and. mod(J,2) == 0)) then
                                P = parity_of(Orb(a)%l + Orb(b)%l)
                                Orb_2qp%N(P,J) = Orb_2qp%N(P,J) + 1
                                if (Pass == 2) then
                                    Orb_2qp%Ind(P,J)%v(Orb_2qp%N(P,J)) = ind2b_t(a,b,T)
                                    Orb_2qp%Dic(a,b,J,(T + 3)/2) = Orb_2qp%N(P,J)
                                end if
                            end if
                        end do
                    end do
                end do
            end do
            if (Pass == 1) then
                do P = 1, 2
                    do J = 0, J_max
                        allocate(Orb_2qp%Ind(P,J)%v(Orb_2qp%N(P,J)))
                    end do
                end do
            end if
        end do
    end function orbitals_2qp_make

    ! Index of the 2qp pair (a,b,J,P,T) ... 0 if absent ... (no caller in dev) ...
    pure integer function orb2b_2qp_index(a,b,J,P,T,Orb_2qp)
        integer, intent(in) :: a, b, J, P, T
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        integer :: k
        orb2b_2qp_index = 0
        if (a < 1 .or. a > Orb_2qp%a_max .or. b < 1 .or. b > Orb_2qp%a_max) return
        if (J < 0 .or. J > Orb_2qp%J_max .or. (T /= -1 .and. T /= 1) .or. P < 1 .or. P > 2) return
        k = Orb_2qp%Dic(a,b,J,(T + 3)/2)
        if (k == 0) return
        ! the parity is implied by l_a + l_b ... check the stored entry ...
        if (k > Orb_2qp%N(P,J)) return
        if (Orb_2qp%Ind(P,J)%v(k)%a /= a .or. Orb_2qp%Ind(P,J)%v(k)%b /= b .or. Orb_2qp%Ind(P,J)%v(k)%T /= T) return
        orb2b_2qp_index = k
    end function orb2b_2qp_index

    ! Natural orbitals of a 1-body density matrix ... counterpart of OBDM_natural_orbitals (orb.jl:551-650) ...
        ! Rho is symmetrized, cleaned outside the (l,j) blocks & split by +1e-10*a on the diagonal (the shift stays in the
        ! returned occupations), diagonalized (ascending), the eigenvectors are reordered by the Hungarian assignment on
        ! 1 - |C| with positive diagonal phases ... Orb_new(a) gets l & j from the weights |C(b,a)|^2 (n counted on the old
        ! Orb list), the occupations pO & nO from the eigenvalues, Rho_diag = diag(Occ) ...
        ! Fix F5: dev takes l & j of both isospins from the proton eigenvectors pC; the port derives them for the neutrons
        ! from nC as well & stops if the two disagree (one orb1b_t carries one l & j for both isospins) ...
    subroutine OBDM_natural_orbitals(Params,Orb,Rho,Orb_new,C,Rho_diag)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), allocatable, intent(out) :: Orb_new(:)
        type(o1b_t), intent(out) :: C, Rho_diag
        real(dp), allocatable :: pRho(:,:), nRho(:,:), pOcc(:), nOcc(:), pC(:,:), nC(:,:), Cost(:,:)
        integer, allocatable :: pReordering(:), nReordering(:)
        real(dp) :: lSum, jSum, nlSum, njSum
        integer :: N_max, a_max, a, b, l_a, j_a, n_a

        N_max = Params%Calc%Nmax
        a_max = a_max_of(N_max)
        allocate(pRho(a_max,a_max), nRho(a_max,a_max), pOcc(a_max), nOcc(a_max), pC(a_max,a_max), nC(a_max,a_max))
        allocate(Cost(a_max,a_max), pReordering(a_max), nReordering(a_max))

        ! Symmetrize, clean & regularize Rho ...
        pRho = 0.5_dp * (Rho%p(1:a_max,1:a_max) + transpose(Rho%p(1:a_max,1:a_max)))
        nRho = 0.5_dp * (Rho%n(1:a_max,1:a_max) + transpose(Rho%n(1:a_max,1:a_max)))
        do b = 1, a_max
            do a = 1, a_max
                if (Orb(a)%j /= Orb(b)%j .or. Orb(a)%l /= Orb(b)%l) then
                    pRho(a,b) = 0.0_dp
                    nRho(a,b) = 0.0_dp
                end if
            end do
        end do
        do a = 1, a_max
            pRho(a,a) = pRho(a,a) + 1.0e-10_dp * real(a,dp)
            nRho(a,a) = nRho(a,a) + 1.0e-10_dp * real(a,dp)
        end do

        ! Diagonalize Rho ...
        call eig_sym(pRho,pOcc,pC)
        call eig_sym(nRho,nOcc,nC)

        ! Hungarian reordering of the natural orbitals & positive diagonal phases ...
        Cost = 1.0_dp - abs(pC)
        call hungarian(Cost,pReordering)
        Cost = 1.0_dp - abs(nC)
        call hungarian(Cost,nReordering)
        pC = pC(:,pReordering)
        nC = nC(:,nReordering)
        pOcc = pOcc(pReordering)
        nOcc = nOcc(nReordering)
        do a = 1, a_max
            if (pC(a,a) < 0.0_dp) pC(:,a) = -pC(:,a)
            if (nC(a,a) < 0.0_dp) nC(:,a) = -nC(:,a)
        end do

        ! Density matrix in the natural-orbital basis ...
        C%p = pC
        C%n = nC
        allocate(Rho_diag%p(a_max,a_max), Rho_diag%n(a_max,a_max))
        Rho_diag%p = 0.0_dp
        Rho_diag%n = 0.0_dp
        do a = 1, a_max
            Rho_diag%p(a,a) = pOcc(a)
            Rho_diag%n(a,a) = nOcc(a)
        end do

        ! New single-particle orbitals ...
        allocate(Orb_new(a_max))
        do a = 1, a_max
            lSum = 0.0_dp
            jSum = 0.0_dp
            nlSum = 0.0_dp
            njSum = 0.0_dp
            do b = 1, a_max
                lSum = lSum + abs(pC(b,a))**2 * real(Orb(b)%l,dp)
                jSum = jSum + abs(pC(b,a))**2 * real(Orb(b)%j,dp)
                nlSum = nlSum + abs(nC(b,a))**2 * real(Orb(b)%l,dp)
                njSum = njSum + abs(nC(b,a))**2 * real(Orb(b)%j,dp)
            end do
            l_a = round_even(lSum)
            j_a = round_even(jSum)
            if (round_even(nlSum) /= l_a .or. round_even(njSum) /= j_a) &
                call fatal('OBDM_natural_orbitals: proton & neutron natural orbital ' // istr(a) // ' have different l/j')
            n_a = 0
            do b = 1, a - 1
                if (l_a == Orb(b)%l .and. j_a == Orb(b)%j) n_a = n_a + 1
            end do
            Orb_new(a) = orb1b_t(a,n_a,l_a,j_a,pOcc(a),nOcc(a))
        end do
    end subroutine OBDM_natural_orbitals

end module nj_orbitals
