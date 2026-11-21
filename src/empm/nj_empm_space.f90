! Restricted single-particle space of the EMPM chain ... HF levels selected by windows of the energy-ordered levels,
    ! counterpart of singpart_coup.dat & of the hole / particle ranges of input_tda_coup.dat of the original EMPM codes ...
module nj_empm_space
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_orbitals, only : orbitals_ph_make, orbitals_import
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
    public :: empm_path, empm_ref_path, EMPM_reference_check, EMPM_reference_orbitals

contains

    ! io/<Path>/EMPM // Suffix ... the EMPM output folder (dev-style "Name" folder: bin/, spectra/, transitions/, Levels.dat,
        ! TDA_summary.dat, 2Phonon/, ..., owned & wiped by EMPM_TDA) ...
    function empm_path(Params,Suffix) result(Path)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Suffix
        character(len=:), allocatable :: Path
        Path = io_path(Params,'/EMPM' // Suffix)
    end function empm_path

    ! io/<Path>/<&empm Reference> // Suffix ... the HF-like reference read by the EMPM (and the NCSM 'EMPM' Hamiltonian) ...
    function empm_ref_path(Params,Suffix) result(Path)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Suffix
        character(len=:), allocatable :: Path
        Path = ref_path(Params,Params%Calc%EMPM%Reference,Suffix)
    end function empm_ref_path

    ! .true. if the reference io/<Path>/<Reference>/bin exists (otherwise an error is printed, as the dev TDA does) ...
        ! H2b.bin (written by HF only with BMF = .true.) is checked when Need_H2b ...
    logical function EMPM_reference_check(Params,Need_H2b) result(Ok)
        type(params_t), intent(in) :: Params
        logical, intent(in) :: Need_H2b
        logical :: Ex
        Ok = is_dir(empm_ref_path(Params,'/bin'))
        if (.not. Ok) then
            write(*,'(/,a)') "Error! ... No precomputed reference ''" // empm_ref_path(Params,'') // "'' is available ... run " // &
                             trim(Params%Calc%EMPM%Reference) // ' first (&empm Reference) ...'
            return
        end if
        if (Need_H2b) then
            inquire(file=empm_ref_path(Params,'/bin/H2b.bin'), exist=Ex)
            if (.not. Ex) call fatal('EMPM ... the reference 2-body Hamiltonian ' // empm_ref_path(Params,'/bin/H2b.bin') // &
                                     ' is missing ... run HF with BMF = .true.')
        end if
    end function EMPM_reference_check

    ! Orbitals of the reference (orb.bin: LHO labels & the occupations of the reference, e.g. the Aufbau occupations of the
        ! dev HF) ...
    function EMPM_reference_orbitals(Params) result(Orb)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        if (.not. is_dir(empm_ref_path(Params,'/bin'))) call fatal('EMPM ... no reference ' // empm_ref_path(Params,'') // &
                                                                 ' ... run ' // trim(Params%Calc%EMPM%Reference) // ' first')
        Orb = orbitals_import(Params,empm_ref_path(Params,'/bin/orb.bin'))
    end function EMPM_reference_orbitals

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

    ! Orb ... the orbitals of the reference (EMPM_reference_orbitals), H1b_File (default the H1b.bin of &empm Reference) ...
        ! its diagonal gives the HF energies; it must be diagonal & the occupations integer (orbitals_ph_make) ...
    subroutine EMPM_space_make(Params,Orb,Space,H1b_File)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(empm_space_t), intent(out) :: Space
        character(len=*), intent(in), optional :: H1b_File
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        type(empm_level_t), allocatable :: Lev_p(:), Lev_n(:)
        integer :: W_pHole(2), W_pPart(2), W_nHole(2), W_nPart(2), N_lev

        write(*,'(/,a)') 'Preparing the restricted EMPM single-particle space ...'
        if (present(H1b_File)) then
            call orbitals_ph_make(Params,Orb,H1b_File,N_Particle,Particle,N_Hole,Hole)
        else
            call orbitals_ph_make(Params,Orb,empm_ref_path(Params,'/bin/H1b.bin'),N_Particle,Particle,N_Hole,Hole)
        end if
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

        open(newunit=u, file=empm_path(Params,'/Levels.dat'), status='replace', action='write')
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
