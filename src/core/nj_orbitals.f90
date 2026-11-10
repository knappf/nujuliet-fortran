! Single-particle orbitals ... counterpart of MatrixElements/Orb.jl ...
module nj_orbitals
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    implicit none
    private

    public :: orbitals_make, orbitals_lj_make, orbitals_export, orbitals_ph_make, orbitals_ph_list
    public :: orbitals_one_phonon_make, orbitals_2qp_make, a_max_of

contains

    pure integer function a_max_of(N_max)
        integer, intent(in) :: N_max
        a_max_of = (N_max + 1)*(N_max + 2)/2
    end function a_max_of

    ! LHO single-particle orbitals in NuHamil ordering with naive ground-state filling ...
    function orbitals_make(Params) result(Orb)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        integer :: A, Z, N_max, a_max, a_i, N, l, n_r, j, b, N_proton, N_neutron
        integer, allocatable :: Orb_arr(:,:)
        logical :: Neutron_Filling, Proton_Filling

        A = Params%Calc%A
        Z = Params%Calc%Z
        N_max = Params%Int%Nmax
        a_max = a_max_of(N_max)

        write(*,'(/,a)') 'Preparing oscillator basis single-particle orbitals ...'
        allocate(Orb(a_max), Orb_arr(a_max,6))
        Orb_arr = 0
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

        if (A /= 0 .and. Z /= 0) then
            N_proton = 0
            N_neutron = 0
            Neutron_Filling = .true.
            Proton_Filling = .true.
            do while (Neutron_Filling)
                do N = 0, N_max
                    do n_r = 0, N/2
                        l = N - 2*n_r
                        do j = 2*l + 1, abs(2*l - 1), -2
                            if (N_neutron + (j + 1) <= (A - Z)) then
                                N_neutron = N_neutron + j + 1
                                do b = 1, a_max
                                    if (Orb_arr(b,2) == n_r .and. Orb_arr(b,3) == l .and. Orb_arr(b,4) == j) Orb_arr(b,6) = 1
                                end do
                            else
                                Neutron_Filling = .false.
                            end if
                        end do
                    end do
                end do
            end do
            do while (Proton_Filling)
                do N = 0, N_max
                    do n_r = 0, N/2
                        l = N - 2*n_r
                        do j = 2*l + 1, abs(2*l - 1), -2
                            if (N_proton + (j + 1) <= Z) then
                                N_proton = N_proton + j + 1
                                do b = 1, a_max
                                    if (Orb_arr(b,2) == n_r .and. Orb_arr(b,3) == l .and. Orb_arr(b,4) == j) Orb_arr(b,5) = 1
                                end do
                            else
                                Proton_Filling = .false.
                            end if
                        end do
                    end do
                end do
            end do
        end if

        do b = 1, a_max
            Orb(b) = orb1b_t(Orb_arr(b,1),Orb_arr(b,2),Orb_arr(b,3),Orb_arr(b,4),Orb_arr(b,5),Orb_arr(b,6))
        end do
        write(*,'(/,a)') 'Orbitals initialized ...'
    end function orbitals_make

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

    subroutine orbitals_export(Params,Orb)
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
                           TAB // ' ' // istr(Orb(a)%pO) // TAB // ' ' // istr(Orb(a)%nO)
        end do
        close(u)
        write(*,'(/,a)') 'Exported LHO basis s.p. orbitals in human-readable format ...'

        write(*,'(/,a)') 'Exporting HF basis s.p. orbitals in human-readable format ...'
        open(newunit=u, file=io_path(Params,'/Orbitals_HF.dat'), status='replace', action='write')
        write(u,'(a)') 'a' // TAB // 'l' // TAB // '2j' // TAB // 'pOc' // TAB // 'nOc' // TAB // 'E_p' // TAB // TAB // TAB // 'E_n'
        do a = 1, a_max
            write(u,'(a)') istr(Orb(a)%a) // TAB // istr(Orb(a)%l) // TAB // istr(Orb(a)%j) // TAB // ' ' // istr(Orb(a)%pO) // &
                           TAB // ' ' // istr(Orb(a)%nO) // TAB // jstr(pSPE(a)) // TAB // jstr(nSPE(a))
        end do
        close(u)
        write(*,'(/,a)') 'Exported HF basis s.p. orbitals in human-readable format ...'
    end subroutine orbitals_export

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

    ! Particle & hole states with HF single-particle energies read from Bin/h_HF.bin ...
    subroutine orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(pninteger_t), intent(out) :: N_Particle, N_Hole
        type(pnsvector_t), intent(out) :: Particle, Hole
        real(dp), allocatable :: pSPE(:), nSPE(:)
        integer :: a_max, a

        a_max = a_max_of(Params%Calc%Nmax)
        call read_h_HF_diag(Params,Orb,pSPE,nSPE)

        N_Hole%p = count(Orb(1:a_max)%pO == 1)
        N_Particle%p = count(Orb(1:a_max)%pO == 0)
        N_Hole%n = count(Orb(1:a_max)%nO == 1)
        N_Particle%n = count(Orb(1:a_max)%nO == 0)
        allocate(Particle%p(N_Particle%p), Particle%n(N_Particle%n), Hole%p(N_Hole%p), Hole%n(N_Hole%n))

        N_Particle = pninteger_t(0,0)
        N_Hole = pninteger_t(0,0)
        do a = 1, a_max
            if (Orb(a)%pO == 0) then
                N_Particle%p = N_Particle%p + 1
                Particle%p(N_Particle%p) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,pSPE(a))
            end if
            if (Orb(a)%pO == 1) then
                N_Hole%p = N_Hole%p + 1
                Hole%p(N_Hole%p) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,pSPE(a))
            end if
            if (Orb(a)%nO == 0) then
                N_Particle%n = N_Particle%n + 1
                Particle%n(N_Particle%n) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,nSPE(a))
            end if
            if (Orb(a)%nO == 1) then
                N_Hole%n = N_Hole%n + 1
                Hole%n(N_Hole%n) = spstate_t(Orb(a)%a,Orb(a)%l,Orb(a)%j,nSPE(a))
            end if
        end do
    end subroutine orbitals_ph_make

    ! Particle-hole phonon indices ... ParticleHole%p(J,P,p,h) with J starting at 0 ...
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

    ! 2-quasiparticle orbitals ... a >= b convention ... Orb_2qp%N(P,0:J_max) ...
    function orbitals_2qp_make(Params,Orb) result(Orb_2qp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(qporb2b_t) :: Orb_2qp
        integer :: N_max, a_max, J_max, a, b, J, P, T, Pass

        N_max = Params%Calc%Nmax
        a_max = a_max_of(N_max)
        J_max = 2*N_max + 1
        allocate(Orb_2qp%N(2,0:J_max), Orb_2qp%i(2,0:J_max))

        ! Pass 1 counts, pass 2 allocates the entries ... pp entries first, then nn entries ...
        do Pass = 1, 2
            Orb_2qp%N = 0
            do T = -1, 1, 2
                do a = 1, a_max
                    do b = 1, a
                        do J = abs(Orb(a)%j - Orb(b)%j)/2, (Orb(a)%j + Orb(b)%j)/2
                            if (a /= b .or. (a == b .and. mod(J,2) == 0)) then
                                P = mod(Orb(a)%l + Orb(b)%l,2) + 1
                                Orb_2qp%N(P,J) = Orb_2qp%N(P,J) + 1
                                if (Pass == 2) Orb_2qp%i(P,J)%v(Orb_2qp%N(P,J)) = ind2qp_t(a,b,T)
                            end if
                        end do
                    end do
                end do
            end do
            if (Pass == 1) then
                do P = 1, 2
                    do J = 0, J_max
                        allocate(Orb_2qp%i(P,J)%v(Orb_2qp%N(P,J)))
                    end do
                end do
            end if
        end do
    end function orbitals_2qp_make

end module nj_orbitals
