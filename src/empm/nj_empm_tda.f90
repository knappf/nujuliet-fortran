! Spherical HF-TDA in the restricted EMPM particle-hole space ... counterpart of tda_pn_coupl_ortog.f90 of the original EMPM codes ...
    ! A(ph,p'h') = (e_p - e_h) delta(ph,p'h') - (-1)^(j_p'+j_h'+J) F^J(p,h,h',p') ... identical to the A matrix of HF_RPA_allocate
        ! when the windows cover the full space ...
    ! Phonon channels are indexed (J,P) with J = 0 .. 2*Nmax+1 as in HF_RPA ... phstate_t%p & %h are positions in Space%Lev ...
module nj_empm_tda
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_tr1b, only : tr1b_t
    use nj_hf_rpa, only : HF_RPA_phonon_count, Spur_Ortho
    use nj_empm_space
    use nj_empm_fmat
    implicit none
    private

    public :: EMPM_TDA_phonons, EMPM_TDA_allocate, EMPM_TDA_diagonalize, EMPM_TDA_strengths

contains

    ! 1p-1h phonon basis ... protons first, particle outer & hole inner loop, as in phbase of the original ...
    subroutine EMPM_TDA_phonons(Params,Space,N_Phonon,Phonon,N_nu,Orb_Phonon)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        integer, intent(out) :: N_Phonon
        type(phstate_t), allocatable, intent(out) :: Phonon(:)
        integer, allocatable, intent(out) :: N_nu(:,:)
        type(ivec_t), allocatable, intent(out) :: Orb_Phonon(:,:)
        integer :: Pass

        do Pass = 1, 2
            N_Phonon = 0
            call add_isospin(Space%pPart,Space%pHole,-1)
            call add_isospin(Space%nPart,Space%nHole,1)
            if (Pass == 1) allocate(Phonon(N_Phonon))
        end do
        call HF_RPA_phonon_count(Params,N_Phonon,Phonon,N_nu,Orb_Phonon)

    contains

        subroutine add_isospin(Part,Hole,t)
            integer, intent(in) :: Part(:), Hole(:), t
            integer :: ip, ih, J
            type(empm_level_t) :: Lp, Lh
            do ip = 1, size(Part)
                Lp = Space%Lev(Part(ip))
                do ih = 1, size(Hole)
                    Lh = Space%Lev(Hole(ih))
                    do J = abs(Lp%j - Lh%j)/2, (Lp%j + Lh%j)/2
                        N_Phonon = N_Phonon + 1
                        if (Pass == 2) Phonon(N_Phonon) = phstate_t(Part(ip),Hole(ih),t,J,mod(Lp%l + Lh%l,2) + 1)
                    end do
                end do
            end do
        end subroutine add_isospin

    end subroutine EMPM_TDA_phonons

    subroutine EMPM_TDA_allocate(Params,Space,F_ph,N_nu,Orb_Phonon,Phonon,A)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(rmat_t), allocatable, intent(out) :: A(:,:)
        integer :: J_max, J, P, N_ph, ind_ph, ind_qg, ph, qg, i_p, i_h, i_q, i_g
        real(dp) :: ASum

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Allocating TDA matrices A from the F-matrix ...'
        allocate(A(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                allocate(A(J,P)%v(N_ph,N_ph))
                !$omp parallel do schedule(dynamic) private(ph,i_p,i_h,ind_qg,qg,i_q,i_g,ASum)
                do ind_ph = 1, N_ph
                    ph = Orb_Phonon(J,P)%v(ind_ph)
                    i_p = Phonon(ph)%p
                    i_h = Phonon(ph)%h
                    do ind_qg = 1, N_ph
                        qg = Orb_Phonon(J,P)%v(ind_qg)
                        i_q = Phonon(qg)%p
                        i_g = Phonon(qg)%h
                        ASum = 0.0_dp
                        if (ph == qg) ASum = Space%Lev(i_p)%E - Space%Lev(i_h)%E
                        ASum = ASum - phase((Space%Lev(i_q)%j + Space%Lev(i_g)%j)/2 + J) * Fmat_ME(F_ph,Space,i_p,i_h,i_g,i_q,J)
                        A(J,P)%v(ind_ph,ind_qg) = ASum
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        write(*,'(a)') TAB // 'Matrices A allocated ...'
    end subroutine EMPM_TDA_allocate

    ! TDA diagonalization with optional orthogonalization of the spurious 1- state ... as TDA_channel of HF_RPA_diagonalize,
        ! the spurious state is the first 1- solution with E = 0 ...
    subroutine EMPM_TDA_diagonalize(Params,Space,A,N_nu,Orb_Phonon,Phonon,TrOp,E_TDA,X_TDA)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        type(rmat_t), intent(in) :: A(0:,:)
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(tr1b_t), intent(in) :: TrOp
        type(rvec_t), allocatable, intent(out) :: E_TDA(:,:)
        type(rmat_t), allocatable, intent(out) :: X_TDA(:,:)
        integer :: J_max, J, P
        logical :: Spur

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Diagonalizing TDA matrices ...'
        allocate(E_TDA(0:J_max,2), X_TDA(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                Spur = Params%Calc%EMPM%Ortho .and. (J == 1) .and. (P == 2)
                call TDA_channel(A(J,P)%v,Spur,Orb_Phonon(J,P)%v,E_TDA(J,P)%v,X_TDA(J,P)%v)
            end do
        end do
        write(*,'(a)') TAB // 'Diagonalization done ...'

    contains

        subroutine TDA_channel(A_JP,Spur,Orb_Ph,E,X)
            real(dp), intent(in) :: A_JP(:,:)
            logical, intent(in) :: Spur
            integer, intent(in) :: Orb_Ph(:)
            real(dp), allocatable, intent(out) :: E(:), X(:,:)
            real(dp), allocatable :: Spur_State(:), V(:,:), U(:,:), A_ort(:,:)
            complex(dp), allocatable :: W(:), Z(:,:)
            integer :: N, Spur_Ind, k
            logical :: Is_Real

            N = size(A_JP,1)
            allocate(E(N), X(N,N))
            if (N == 0) return
            if (Spur .and. N > 1) then
                Spur_State = spurious_state(Orb_Ph)
                call Spur_Ortho(N,Spur_State,V,Spur_Ind)
                allocate(U(N,N-1))
                do k = 1, Spur_Ind - 1
                    U(:,k) = V(:,k)
                end do
                do k = Spur_Ind + 1, N
                    U(:,k-1) = V(:,k)
                end do
                A_ort = mtm(U,mm(A_JP,U))
                allocate(W(N-1), Z(N-1,N-1))
                call eig_general(A_ort,W,Z,Is_Real)
                E(1) = 0.0_dp
                E(2:N) = real(W)
                X(:,1) = Spur_State
                X(:,2:N) = mm(U,real(Z))
                write(*,'(a)') TAB // 'Spurious 1- state removed ... <spur|A|spur> = ' // &
                               cfmt_f(dot_product(Spur_State,matmul(A_JP,Spur_State)),0,6) // ' MeV'
            else
                allocate(W(N), Z(N,N))
                call eig_general(A_JP,W,Z,Is_Real)
                E = real(W)
                X = real(Z)
            end if
        end subroutine TDA_channel

        ! Normalized spurious (CM) state from the E1 operator restricted to the phonon basis ... as Spur_Ini ...
        function spurious_state(Orb_Ph) result(ph_CMS)
            integer, intent(in) :: Orb_Ph(:)
            real(dp), allocatable :: ph_CMS(:)
            integer :: k, a_p, a_h
            allocate(ph_CMS(size(Orb_Ph)))
            do k = 1, size(Orb_Ph)
                a_p = Space%Lev(Phonon(Orb_Ph(k))%p)%a
                a_h = Space%Lev(Phonon(Orb_Ph(k))%h)%a
                if (Phonon(Orb_Ph(k))%tz == 1) then
                    ph_CMS(k) = TrOp%E1%n(a_p,a_h)
                else
                    ph_CMS(k) = TrOp%E1%p(a_p,a_h)
                end if
            end do
            where (abs(ph_CMS) < 1.0e-6_dp) ph_CMS = 0.0_dp
            if (norm2(ph_CMS) == 0.0_dp) call fatal('No spurious 1- component in the restricted EMPM space ... set Ortho = .false.')
            ph_CMS = ph_CMS / norm2(ph_CMS)
        end function spurious_state

    end subroutine EMPM_TDA_diagonalize

    ! B(EL; 0 -> nu) of the TDA phonons in the window, L = 0..3 (E0 r^2, E1, E2, E3 of nj_tr1b) ... M_t = sum_ph
        ! <p||Q_L||h> X_ph (as HF_RPA_rM), B_ph = M_p^2, B_is = |M_p + M_n|^2/4, B_iv = |e_p M_p - e_n M_n|^2 (E1) or
        ! |M_p - M_n|^2/4 ... EMPM/Transitions/TDA_E<L>.dat (E, B_ph, B_is, B_iv as HF_RPA) ...
    subroutine EMPM_TDA_strengths(Params,Space,N_nu,Orb_Phonon,Phonon,TrOp,E_TDA,X_TDA)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(tr1b_t), intent(in) :: TrOp
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        real(dp) :: e_p, e_n, Mp, Mn, O, Biv
        integer :: L, P, nu, ph, i_ph, a_p, a_h, u
        e_p = real(Params%Calc%A - Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        do L = 0, 3
            P = mod(L,2) + 1
            if (.not. is_dir(io_path(Params,'/EMPM/Transitions'))) call make_dir(io_path(Params,'/EMPM/Transitions'))
            open(newunit=u, file=io_path(Params,'/EMPM/Transitions/TDA_E' // istr(L) // '.dat'), status='replace', action='write')
            write(u,'(4a21)') 'E                    ', 'B_ph                 ', 'B_is                 ', 'B_iv                 '
            do nu = 1, N_nu(L,P)
                Mp = 0.0_dp
                Mn = 0.0_dp
                do ph = 1, N_nu(L,P)
                    i_ph = Orb_Phonon(L,P)%v(ph)
                    a_p = Space%Lev(Phonon(i_ph)%p)%a
                    a_h = Space%Lev(Phonon(i_ph)%h)%a
                    select case (L)
                    case (0)
                        O = merge(TrOp%E0%p(a_p,a_h), TrOp%E0%n(a_p,a_h), Phonon(i_ph)%tz == -1)
                    case (1)
                        O = merge(TrOp%E1%p(a_p,a_h), TrOp%E1%n(a_p,a_h), Phonon(i_ph)%tz == -1)
                    case (2)
                        O = merge(TrOp%E2%p(a_p,a_h), TrOp%E2%n(a_p,a_h), Phonon(i_ph)%tz == -1)
                    case default
                        O = merge(TrOp%E3%p(a_p,a_h), TrOp%E3%n(a_p,a_h), Phonon(i_ph)%tz == -1)
                    end select
                    if (Phonon(i_ph)%tz == -1) then
                        Mp = Mp + O*X_TDA(L,P)%v(ph,nu)
                    else
                        Mn = Mn + O*X_TDA(L,P)%v(ph,nu)
                    end if
                end do
                if (L == 1) then
                    Biv = (e_p*Mp - e_n*Mn)**2
                else
                    Biv = 0.25_dp*(Mp - Mn)**2
                end if
                write(u,'(4f21.8)') E_TDA(L,P)%v(nu), Mp**2, 0.25_dp*(Mp + Mn)**2, Biv
            end do
            close(u)
        end do
        write(*,'(a)') TAB // 'TDA transition strengths B(E0..E3) written to EMPM/Transitions/TDA_E<L>.dat'
    end subroutine EMPM_TDA_strengths

end module nj_empm_tda
