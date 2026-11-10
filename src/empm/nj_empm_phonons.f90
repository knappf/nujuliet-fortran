! TDA phonons of the restricted EMPM space as the 1-phonon building blocks of the EMPM ... counterpart of 1phonon/1f_*.dat
    ! & tda_r_overl.dat of the original EMPM codes, rebuilt from Bin/EMPM_TDA_*.bin & Bin/EMPM_F.bin ...
    ! Running phonon index over the (J,P) channels, J outer & P inner, as in EMPM/TDA_Structure.dat ... or with
        ! Legacy_Order the numbering of 1phonon/1f_states.dat of the original Tda: parity - before +, J inner, the spurious
        ! phonon the last one of the 1- channel ... phonons of one channel stay contiguous in both cases ...
module nj_empm_phonons
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_empm_space
    use nj_empm_fmat, only : fmat_t
    use nj_empm_tda, only : EMPM_TDA_phonons, EMPM_TDA_allocate
    implicit none
    private

    ! N phonons with J, parity P (1 = +, 2 = -) & energy E ... channel (J,P) holds the phonons First(J,P) .. First(J,P)+N_ch(J,P)-1 ...
        ! Cp(p,h,nu) & Cn(p,h,nu) amplitudes over the window particles x holes (positions in Space%pPart & Space%pHole etc.) ...
        ! i_cm spurious 1- phonon (0 = none), R(nu) = <nu|A|cm> for the 1- phonons ... E(i_cm) = R(i_cm) = <cm|A|cm>,
        ! the unprojected TDA energy of the spurious phonon used by the original EMPM codes ...
    type, public :: empm_phonons_t
        integer :: N = 0, J_max = 0, i_cm = 0
        integer, allocatable :: J(:), P(:)
        real(dp), allocatable :: E(:), R(:)
        integer, allocatable :: N_ch(:,:), First(:,:)
        real(dp), allocatable :: Cp(:,:,:), Cn(:,:,:)
        ! Cpt(h,p,nu) & Cnt(h,p,nu) ... the same amplitudes with the hole index first (unit stride over h for dgemm) ...
        real(dp), allocatable :: Cpt(:,:,:), Cnt(:,:,:)
    end type empm_phonons_t

    public :: EMPM_phonons_import

contains

    subroutine EMPM_phonons_import(Params,Space,F_ph,Ph)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        type(empm_phonons_t), intent(out) :: Ph
        type(phstate_t), allocatable :: Phonon(:)
        integer, allocatable :: N_nu(:,:), Loc(:)
        type(ivec_t), allocatable :: Orb_Phonon(:,:)
        type(rmat_t), allocatable :: A(:,:)
        real(dp), allocatable :: X(:,:), Rm(:,:)
        integer, allocatable :: Pos(:)
        integer :: N_Phonon, J_max, J, P, n, nu, k, iph, i, uE, uX, ip, Pc
        logical :: Spur_Ch

        write(*,'(/,a)') 'Importing the EMPM TDA phonons ...'
        J_max = 2*Params%Calc%Nmax + 1
        call EMPM_TDA_phonons(Params,Space,N_Phonon,Phonon,N_nu,Orb_Phonon)

        ! Local index of each level within its block ...
        allocate(Loc(Space%N_lev))
        Loc(Space%pHole) = [(k, k = 1, size(Space%pHole))]
        Loc(Space%pPart) = [(k, k = 1, size(Space%pPart))]
        Loc(Space%nHole) = [(k, k = 1, size(Space%nHole))]
        Loc(Space%nPart) = [(k, k = 1, size(Space%nPart))]

        Ph%N = N_Phonon
        Ph%J_max = 0
        allocate(Ph%J(Ph%N), Ph%P(Ph%N), Ph%E(Ph%N), Ph%R(Ph%N), Ph%N_ch(0:J_max,2), Ph%First(0:J_max,2))
        allocate(Ph%Cp(size(Space%pPart),size(Space%pHole),Ph%N), Ph%Cn(size(Space%nPart),size(Space%nHole),Ph%N))
        Ph%Cp = 0.0_dp
        Ph%Cn = 0.0_dp
        Ph%R = 0.0_dp
        Ph%N_ch = N_nu

        ! First phonon of each channel in the chosen numbering ...
        i = 0
        if (Params%Calc%EMPM%Legacy_Order) then
            do ip = 1, 2
                Pc = 3 - ip
                do J = 0, J_max
                    Ph%First(J,Pc) = i + 1
                    i = i + N_nu(J,Pc)
                end do
            end do
        else
            do J = 0, J_max
                do P = 1, 2
                    Ph%First(J,P) = i + 1
                    i = i + N_nu(J,P)
                end do
            end do
        end if

        open(newunit=uE, file=io_path(Params,'/Bin/EMPM_TDA_E.bin'), access='stream', form='unformatted', status='old', action='read')
        open(newunit=uX, file=io_path(Params,'/Bin/EMPM_TDA_X.bin'), access='stream', form='unformatted', status='old', action='read')
        ! Binary files hold the channels J outer & P inner ... solution nu of a channel goes to First + Pos(nu) - 1 ...
        do J = 0, J_max
            do P = 1, 2
                n = N_nu(J,P)
                if (n == 0) cycle
                Ph%J_max = J
                i = Ph%First(J,P) - 1
                Spur_Ch = Params%Calc%EMPM%Ortho .and. J == 1 .and. P == 2
                Pos = [(nu, nu = 1, n)]
                if (Spur_Ch .and. Params%Calc%EMPM%Legacy_Order) Pos = [n, (nu, nu = 1, n - 1)]
                allocate(X(n,n))
                read(uE) Ph%E(i+Pos)
                read(uX) X
                do nu = 1, n
                    Ph%J(i+Pos(nu)) = J
                    Ph%P(i+Pos(nu)) = P
                    do k = 1, n
                        iph = Orb_Phonon(J,P)%v(k)
                        if (Phonon(iph)%tz == -1) then
                            Ph%Cp(Loc(Phonon(iph)%p),Loc(Phonon(iph)%h),i+Pos(nu)) = X(k,nu)
                        else
                            Ph%Cn(Loc(Phonon(iph)%p),Loc(Phonon(iph)%h),i+Pos(nu)) = X(k,nu)
                        end if
                    end do
                end do
                ! Couplings of the 1- phonons to the spurious one (the first 1- solution) through the TDA matrix A ...
                if (Spur_Ch) then
                    if (.not. allocated(A)) call EMPM_TDA_allocate(Params,Space,F_ph,N_nu,Orb_Phonon,Phonon,A)
                    Rm = matmul(transpose(X),matmul(A(J,P)%v,X))
                    Ph%i_cm = i + Pos(1)
                    Ph%R(i+Pos) = Rm(:,1)
                end if
                deallocate(X)
            end do
        end do
        close(uE)
        close(uX)
        if (Ph%i_cm > 0) Ph%E(Ph%i_cm) = Ph%R(Ph%i_cm)
        Ph%Cpt = reshape(Ph%Cp,[size(Ph%Cp,2),size(Ph%Cp,1),Ph%N],order=[2,1,3])
        Ph%Cnt = reshape(Ph%Cn,[size(Ph%Cn,2),size(Ph%Cn,1),Ph%N],order=[2,1,3])
        write(*,'(a)') TAB // 'Number of phonons = ' // istr(Ph%N) // ' , J_max = ' // istr(Ph%J_max)
        if (Ph%i_cm > 0) write(*,'(a)') TAB // 'Spurious phonon i = ' // istr(Ph%i_cm) // ' , <cm|A|cm> = ' // &
                                        cfmt_f(Ph%E(Ph%i_cm),0,6) // ' MeV'
    end subroutine EMPM_phonons_import

end module nj_empm_phonons
