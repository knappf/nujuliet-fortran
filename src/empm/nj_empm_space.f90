! Restricted single-particle space of the EMPM chain ... HF levels selected by windows of the energy-ordered levels,
    ! counterpart of singpart_coup.dat & of the hole / particle ranges of input_tda_coup.dat of the original EMPM codes ...
module nj_empm_space
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_orbitals, only : orbitals_ph_make
    implicit none
    private

    ! HF level ... orbital a, isospin t (-1 proton, +1 neutron), energy rank k within its isospin (1 = deepest),
        ! LHO labels n, l, 2j of the orbital, occupation Occ (1 hole, 0 particle) & HF energy E ...
    type, public :: empm_level_t
        integer :: a = 0, t = 0, k = 0, n = 0, l = 0, j = 0, Occ = 0
        real(dp) :: E = 0.0_dp
    end type empm_level_t

    ! Levels of the windows Lev(1:N_lev) ... proton holes, proton particles, neutron holes, neutron particles,
        ! each block ordered by energy ... pHole, pPart, nHole & nPart hold the positions of the blocks in Lev ...
    type, public :: empm_space_t
        integer :: N_lev = 0
        type(empm_level_t), allocatable :: Lev(:)
        integer, allocatable :: pHole(:), pPart(:), nHole(:), nPart(:)
    end type empm_space_t

    public :: EMPM_space_make, EMPM_space_export, level_label, EMPM_npnh_dims

contains

    ! Number of npnh states of the window per J^P (J = 0 .. J_max, P = 1 +, 2 -), kp proton & n-kp neutron particle-hole
        ! pairs ... m-scheme counts N(2M, parity) from the products of the particle & hole subsets, d_J = N(2J) - N(2J+2)
        ! (the dimension of the complete n-phonon level, i.e. the rank of its overcomplete basis) ...
    function EMPM_npnh_dims(Space,n,J_max) result(d)
        type(empm_space_t), intent(in) :: Space
        integer, intent(in) :: n, J_max
        integer(i8) :: d(0:J_max,2)
        integer(i8), allocatable :: C(:,:,:,:), T(:,:)
        integer :: Mx, kp, k, M, P
        Mx = 0
        do k = 1, Space%N_lev
            Mx = max(Mx, Space%Lev(k)%j)
        end do
        Mx = 2*n*Mx
        ! C(k, 2M, parity, block) ... k m-states chosen in block 1 pPart, 2 pHole, 3 nPart, 4 nHole ...
        allocate(C(0:n,-Mx:Mx,0:1,4), T(-Mx:Mx,0:1))
        call subsets(Space%pPart,C(:,:,:,1))
        call subsets(Space%pHole,C(:,:,:,2))
        call subsets(Space%nPart,C(:,:,:,3))
        call subsets(Space%nHole,C(:,:,:,4))
        T = 0
        do kp = 0, n
            T = T + conv(conv(C(kp,:,:,1),C(kp,:,:,2)),conv(C(n-kp,:,:,3),C(n-kp,:,:,4)))
        end do
        d = 0
        do P = 1, 2
            do M = 0, J_max
                if (2*M + 2 <= Mx) then
                    d(M,P) = T(2*M,P-1) - T(2*M+2,P-1)
                else if (2*M <= Mx) then
                    d(M,P) = T(2*M,P-1)
                end if
            end do
        end do
    contains
        subroutine subsets(Pos,Cb)
            integer, intent(in) :: Pos(:)
            integer(i8), intent(out) :: Cb(0:,-Mx:,0:)
            integer :: i, m2, kk, q, j2, pl, mm, pp
            Cb = 0
            Cb(0,0,0) = 1
            do i = 1, size(Pos)
                j2 = Space%Lev(Pos(i))%j
                pl = mod(Space%Lev(Pos(i))%l,2)
                do m2 = -j2, j2, 2
                    do kk = n, 1, -1
                        do pp = 0, 1
                            do mm = -Mx, Mx
                                if (Cb(kk-1,mm,pp) == 0) cycle
                                q = mm + m2
                                if (abs(q) > Mx) cycle
                                Cb(kk,q,mod(pp+pl,2)) = Cb(kk,q,mod(pp+pl,2)) + Cb(kk-1,mm,pp)
                            end do
                        end do
                    end do
                end do
            end do
        end subroutine subsets
        function conv(X,Y) result(Z)
            integer(i8), intent(in) :: X(-Mx:,0:), Y(-Mx:,0:)
            integer(i8) :: Z(-Mx:Mx,0:1)
            integer :: M1, P1, M2, P2
            Z = 0
            do P1 = 0, 1
                do M1 = -Mx, Mx
                    if (X(M1,P1) == 0) cycle
                    do P2 = 0, 1
                        do M2 = max(-Mx,-Mx-M1), min(Mx,Mx-M1)
                            Z(M1+M2,mod(P1+P2,2)) = Z(M1+M2,mod(P1+P2,2)) + X(M1,P1)*Y(M2,P2)
                        end do
                    end do
                end do
            end do
        end function conv
    end function EMPM_npnh_dims

    subroutine EMPM_space_make(Params,Orb,Space)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(empm_space_t), intent(out) :: Space
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        type(empm_level_t), allocatable :: Lev_p(:), Lev_n(:)
        integer :: W_pHole(2), W_pPart(2), W_nHole(2), W_nPart(2), N_lev

        write(*,'(/,a)') 'Preparing the restricted EMPM single-particle space ...'
        call orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        Lev_p = ranked_levels(-1,Hole%p,Particle%p)
        Lev_n = ranked_levels(1,Hole%n,Particle%n)

        W_pHole = window(Params%Calc%EMPM%pHole,1,N_Hole%p)
        W_pPart = window(Params%Calc%EMPM%pPart,N_Hole%p + 1,size(Lev_p))
        W_nHole = window(Params%Calc%EMPM%nHole,1,N_Hole%n)
        W_nPart = window(Params%Calc%EMPM%nPart,N_Hole%n + 1,size(Lev_n))
        call check_window(Lev_p,W_pHole,1,'pHole')
        call check_window(Lev_p,W_pPart,0,'pPart')
        call check_window(Lev_n,W_nHole,1,'nHole')
        call check_window(Lev_n,W_nPart,0,'nPart')

        N_lev = 0
        allocate(Space%Lev(width(W_pHole) + width(W_pPart) + width(W_nHole) + width(W_nPart)))
        call add_block(Lev_p,W_pHole,Space%pHole)
        call add_block(Lev_p,W_pPart,Space%pPart)
        call add_block(Lev_n,W_nHole,Space%nHole)
        call add_block(Lev_n,W_nPart,Space%nPart)
        Space%N_lev = N_lev

        write(*,'(a)') TAB // 'Proton holes     ' // istr(W_pHole(1)) // ' .. ' // istr(W_pHole(2)) // TAB // &
                       'proton particles  ' // istr(W_pPart(1)) // ' .. ' // istr(W_pPart(2))
        write(*,'(a)') TAB // 'Neutron holes    ' // istr(W_nHole(1)) // ' .. ' // istr(W_nHole(2)) // TAB // &
                       'neutron particles ' // istr(W_nPart(1)) // ' .. ' // istr(W_nPart(2))
        write(*,'(a)') TAB // 'Number of levels in the restricted space = ' // istr(Space%N_lev)

    contains

        integer function width(W)
            integer, intent(in) :: W(2)
            width = max(W(2) - W(1) + 1, 0)
        end function width

        subroutine add_block(Lev,W,Pos)
            type(empm_level_t), intent(in) :: Lev(:)
            integer, intent(in) :: W(2)
            integer, allocatable, intent(out) :: Pos(:)
            integer :: k
            allocate(Pos(width(W)))
            do k = W(1), W(2)
                N_lev = N_lev + 1
                Space%Lev(N_lev) = Lev(k)
                Pos(k - W(1) + 1) = N_lev
            end do
        end subroutine add_block

    end subroutine EMPM_space_make

    ! All HF levels of isospin t ordered by energy ... stable sort, rank k ...
    function ranked_levels(t,H,Pt) result(Lev)
        integer, intent(in) :: t
        type(spstate_t), intent(in) :: H(:), Pt(:)
        type(empm_level_t), allocatable :: Lev(:), Tmp(:)
        integer, allocatable :: Perm(:)
        real(dp), allocatable :: E(:)
        integer :: N_h, k

        N_h = size(H)
        allocate(Tmp(N_h + size(Pt)), Perm(N_h + size(Pt)))
        do k = 1, N_h
            Tmp(k) = empm_level_t(H(k)%a,t,0,0,H(k)%l,H(k)%j,1,H(k)%E)
        end do
        do k = 1, size(Pt)
            Tmp(N_h + k) = empm_level_t(Pt(k)%a,t,0,0,Pt(k)%l,Pt(k)%j,0,Pt(k)%E)
        end do
        E = Tmp%E
        call sortperm(E,Perm)
        Lev = Tmp(Perm)
        do k = 1, size(Lev)
            Lev(k)%k = k
        end do
        ! Radial label n counted within each (l,j) block in the order of energy ...
        do k = 1, size(Lev)
            Lev(k)%n = count(Lev(1:k-1)%l == Lev(k)%l .and. Lev(1:k-1)%j == Lev(k)%j)
        end do
    end function ranked_levels

    ! Window (first, last) with zero bounds replaced by the defaults ...
    pure function window(Bounds,Lo,Hi) result(W)
        integer, intent(in) :: Bounds(2), Lo, Hi
        integer :: W(2)
        W = Bounds
        if (W(1) == 0) W(1) = Lo
        if (W(2) == 0) W(2) = Hi
    end function window

    subroutine check_window(Lev,W,Occ,Name)
        type(empm_level_t), intent(in) :: Lev(:)
        integer, intent(in) :: W(2), Occ
        character(len=*), intent(in) :: Name
        integer :: k
        if (W(1) > W(2)) return
        if (W(1) < 1 .or. W(2) > size(Lev)) call fatal('EMPM window ' // Name // ' = ' // istr(W(1)) // ', ' // istr(W(2)) // &
                                                       ' outside of the level range 1 .. ' // istr(size(Lev)))
        do k = W(1), W(2)
            if (Lev(k)%Occ /= Occ) call fatal('EMPM window ' // Name // ' contains the level ' // istr(k) // ' (' // &
                                              level_label(Lev(k)) // ', E = ' // cfmt_f(Lev(k)%E,0,5) // ' MeV) of wrong occupation')
        end do
    end subroutine check_window

    ! Spectroscopic label of the level ... e.g. 0h11/2 ...
    function level_label(Lev) result(s)
        type(empm_level_t), intent(in) :: Lev
        character(len=:), allocatable :: s
        character(len=*), parameter :: Orbit = 'spdfghijklmnoqrtuvwxyz'
        s = istr(Lev%n) // Orbit(Lev%l+1:Lev%l+1) // istr(Lev%j) // '/2'
    end function level_label

    subroutine EMPM_space_export(Params,Space)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        integer :: u, i
        character(len=1) :: tc

        open(newunit=u, file=io_path(Params,'/EMPM/Levels.dat'), status='replace', action='write')
        write(u,'(a)') 'Restricted EMPM single-particle space ... HF levels ordered by energy within each isospin (k = 1 deepest)'
        write(u,'(/,a)') cfmt_s('i',5) // cfmt_s('t',3) // cfmt_s('k',5) // cfmt_s('a',5) // cfmt_s('n',5) // cfmt_s('l',5) // &
                         cfmt_s('2j',5) // cfmt_s('Occ',5) // cfmt_s('E',15) // '   ' // 'Level'
        do i = 1, Space%N_lev
            tc = merge('p', 'n', Space%Lev(i)%t == -1)
            write(u,'(a)') cfmt_d(i,5) // cfmt_s(tc,3) // cfmt_d(Space%Lev(i)%k,5) // cfmt_d(Space%Lev(i)%a,5) // &
                           cfmt_d(Space%Lev(i)%n,5) // cfmt_d(Space%Lev(i)%l,5) // cfmt_d(Space%Lev(i)%j,5) // &
                           cfmt_d(Space%Lev(i)%Occ,5) // cfmt_f(Space%Lev(i)%E,15,5) // '   ' // level_label(Space%Lev(i))
        end do
        close(u)
    end subroutine EMPM_space_export

end module nj_empm_space
