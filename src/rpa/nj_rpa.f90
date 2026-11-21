! Random-phase approximation on a general (HF / correlated) reference (NuJuliet dev) ... counterpart of
    ! src/RPA/RPA_diagonalize.jl, RPA_rho.jl (RPA_OBDM_oY2 only), RPA_energy.jl, RPA_transitions.jl (RPA_rM, RPA_rM_calculate,
    ! RPA_rM_contract) & RPA_sigma.jl (RPA_sigma incl. fix F2) ... RPA_allocate / RPA_allocate_TBCF / RPA_spurious / RPA_rB are
    ! the shared routines of nj_ph_common ...
    ! Not ported (dead in dev): the Sigma variant of RPA_spurious_ini (RPA_spurious.jl:40-130), RPA_OBDM_sigma ...
module nj_rpa
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg, only : eig_sym, eig_general, mm, mtm, mmt, zmm
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn, O2b_initialize
    use nj_v2b, only : V2b_correlation_energy
    use nj_orbitals, only : orb2b_1p1h_index
    use nj_tr1b, only : tr1b_t
    use nj_ph_common
    use nj_tda, only : ph_rM_assemble
    use nj_rpa_se, only : RPA_allocate_SE
    implicit none
    private

    public :: rpa_energy_t, RPA_diagonalize, RPA_OBDM_oY2, RPA_energy, RPA_energy_QBA, RPA_energy_dress, RPA_energy_rho
    public :: RPA_rM, RPA_sigma, RPA_sigma_mix, RPA_allocate

    ! Julia NamedTuple (E_QBA, E_sigma, E_rho, E_dress, E_dens) ... NaN where not evaluated ...
    type :: rpa_energy_t
        real(dp) :: E_QBA = 0.0_dp, E_sigma = 0.0_dp, E_rho = 0.0_dp, E_dress = 0.0_dp, E_dens = 0.0_dp
    end type rpa_energy_t

contains

    ! RPA matrices A, B & N ... counterpart of RPA_allocate (RPA_allocate.jl:1-185): ph_allocate, for TBCF the Sigma box
        ! diagrams (RPA_allocate_TBCF) & with RPA.DEV the 1p-1h to 2p-2h self-energy coupling of A (RPA_allocate_SE,
        ! Mode = "OnShell") ... Rho is unused (as in dev) ... stage 4 (shared by RPA_solver & HF_ScRPA) ...
    subroutine RPA_allocate(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Rho,Sigma,TBCF,A,B,N)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: h_N, Rho
        type(o2b_t), intent(in) :: V_NN, Sigma
        logical, intent(in) :: TBCF
        type(rmat_t), allocatable, intent(out) :: A(:,:), B(:,:), N(:,:)

        write(*,'(/,a)') 'Allocating RPA matrices A, B & N ...'
        call ph_allocate(Orb,Orb_1p1h,Orb_NN,h_N,V_NN,.true.,A,B,N)
        if (TBCF) then
            call ph_allocate_tbcf(Params,Orb,Orb_1p1h,Orb_NN,V_NN,Sigma,.true.,Params%Calc%RPA%DEV,A,B,'RPA')
            if (Params%Calc%RPA%DEV) call RPA_allocate_SE(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,A)
        end if
        write(*,'(/,a)') TAB // 'RPA matrices A, B & N allocated ...'
    end subroutine RPA_allocate

    ! Linear mixing Sigma = gamma Sigma + (1 - gamma) Sigma_old of every allocated pp, pn & nn block (no-op for
        ! |1 - gamma| < 1e-12) ... counterpart of RPA_sigma_mix! (RPA_sigma.jl:420-434) ... stage 4 ...
    subroutine RPA_sigma_mix(Sigma,Sigma_old,gamma)
        type(o2b_t), intent(inout) :: Sigma
        type(o2b_t), intent(in) :: Sigma_old
        real(dp), intent(in) :: gamma
        real(dp) :: ca, cb
        integer :: i, k

        if (abs(1.0_dp - gamma) < 1.0e-12_dp) return
        ca = gamma
        cb = 1.0_dp - gamma
        do k = lbound(Sigma%pp,2), ubound(Sigma%pp,2)
            do i = lbound(Sigma%pp,1), ubound(Sigma%pp,1)
                if (allocated(Sigma%pp(i,k)%v)) Sigma%pp(i,k)%v = ca * Sigma%pp(i,k)%v + cb * Sigma_old%pp(i,k)%v
            end do
        end do
        do k = lbound(Sigma%pn,2), ubound(Sigma%pn,2)
            do i = lbound(Sigma%pn,1), ubound(Sigma%pn,1)
                if (allocated(Sigma%pn(i,k)%v)) Sigma%pn(i,k)%v = ca * Sigma%pn(i,k)%v + cb * Sigma_old%pn(i,k)%v
            end do
        end do
        do k = lbound(Sigma%nn,2), ubound(Sigma%nn,2)
            do i = lbound(Sigma%nn,1), ubound(Sigma%nn,1)
                if (allocated(Sigma%nn(i,k)%v)) Sigma%nn(i,k)%v = ca * Sigma%nn(i,k)%v + cb * Sigma_old%nn(i,k)%v
            end do
        end do
    end subroutine RPA_sigma_mix

    ! RPA eigenvalue problem per channel in the product form H = H_p H_m with the metric N ... counterpart of RPA_diagonalize
        ! (RPA_diagonalize.jl) ... the orthogonalized 1- channel solves the physical states in the complement of Spur & the
        ! spurious state from the 2x2 problem of Spur (stored as the 1st solution, X = Y = 0 for zero energy / norm) ...
    subroutine RPA_diagonalize(Params,Orb,Orb_1p1h,A,B,N,TrOp,E,X,Y,Stability)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rmat_t), intent(in) :: A(:,0:), B(:,0:), N(:,0:)
        type(tr1b_t), intent(in) :: TrOp
        type(cvec_t), allocatable, intent(out) :: E(:,:)
        type(cmat_t), allocatable, intent(out) :: X(:,:), Y(:,:)
        logical, intent(out) :: Stability
        real(dp), allocatable :: Nd(:), s(:), Km(:,:), H(:,:), Spur(:), U_sub(:,:), NU(:,:), N_ort(:,:), lam(:), V_N(:,:)
        real(dp), allocatable :: inv_sqrtN(:,:), U_invN(:,:), AB_U(:,:), H_p(:,:), H_m(:,:), Hr_re(:,:), Hr_im(:,:)
        complex(dp), allocatable :: W(:), R(:,:)
        integer :: J_max, J, P, N_ph, N_o, i, k
        logical :: Ortho, Is_Real
        real(dp) :: A_s, B_s, N_s, c_s, Hp_s, Hm_s, r_s, p_s, nrm_s
        complex(dp) :: E_s, x_s, y_s

        J_max = Orb_1p1h%J_max
        Ortho = Params%Calc%RPA%Ortho
        Stability = .true.
        allocate(E(2,0:J_max), X(2,0:J_max), Y(2,0:J_max))
        write(*,'(/,a)') 'Starting diagonalization of the RPA system of equations ...'

        do J = 0, J_max
            do P = 1, 2
                N_ph = Orb_1p1h%N(P,J)
                allocate(E(P,J)%v(N_ph), X(P,J)%v(N_ph,N_ph), Y(P,J)%v(N_ph,N_ph))
                write(*,'(a)') TAB // TAB // 'Diagonalizing:' // TAB // 'J = ' // istr(J) // ', P = ' // pm_str(P)
                if (N_ph == 0) cycle
                allocate(Nd(N_ph))
                do i = 1, N_ph
                    Nd(i) = N(P,J)%v(i,i)
                end do

                if ((.not. Ortho) .or. .not. (J == 1 .and. P == 2)) then
                    s = 1.0_dp / sqrt(Nd)
                    allocate(Km(N_ph,N_ph))
                    do k = 1, N_ph
                        do i = 1, N_ph
                            Km(i,k) = s(i)**2 * (A(P,J)%v(i,k) - B(P,J)%v(i,k)) * s(k)
                        end do
                    end do
                    H = mm(A(P,J)%v,Km)
                    H = H + mm(B(P,J)%v,Km)
                    do k = 1, N_ph
                        H(:,k) = H(:,k) * s
                    end do
                    allocate(W(N_ph), R(N_ph,N_ph))
                    call eig_general(H,W,R,Is_Real)
                    E(P,J)%v = sqrt(W)
                    ! (X - Y) = N^-1/2 R in X & E (X + Y) = K R in Y ...
                    do k = 1, N_ph
                        X(P,J)%v(:,k) = s * R(:,k)
                    end do
                    if (Is_Real) then
                        Y(P,J)%v = cmplx(mm(Km,real(R,dp)), 0.0_dp, kind=dp)
                    else
                        Y(P,J)%v = cmplx(mm(Km,real(R,dp)), mm(Km,aimag(R)), kind=dp)
                    end if
                    call amplitudes(1)
                    deallocate(Km, H, W, R, s)
                else
                    N_o = N_ph - 1
                    Spur = ph_spurious_ini(Orb,Orb_1p1h,TrOp)
                    U_sub = ph_spurious_basis(Spur)
                    if (N_o > 0) then
                        ! Symmetric inverse square root of N_ort = U_sub' N U_sub ...
                        allocate(NU(N_ph,N_o), lam(N_o), V_N(N_o,N_o))
                        do k = 1, N_o
                            NU(:,k) = Nd * U_sub(:,k)
                        end do
                        N_ort = mtm(U_sub,NU)
                        call eig_sym(N_ort,lam,V_N)
                        do k = 1, N_o
                            V_N(:,k) = V_N(:,k) * lam(k)**(-0.25_dp)
                        end do
                        inv_sqrtN = mmt(V_N,V_N)
                        U_invN = mm(U_sub,inv_sqrtN)
                        ! H_+- = U_invN' (A +- B) U_invN ...
                        AB_U = mm(A(P,J)%v,U_invN)
                        H_p = mtm(U_invN,AB_U)
                        AB_U = mm(B(P,J)%v,U_invN)
                        H_m = mtm(U_invN,AB_U)
                        H_p = H_p + H_m
                        H_m = H_p - 2.0_dp * H_m
                        H = mm(H_p,H_m)
                        allocate(W(N_o), R(N_o,N_o))
                        call eig_general(H,W,R,Is_Real)
                        E(P,J)%v(2:N_ph) = sqrt(W)
                        Hr_re = mm(H_m,real(R,dp))
                        if (Is_Real) then
                            X(P,J)%v(:,2:N_ph) = cmplx(mm(U_invN,real(R,dp)), 0.0_dp, kind=dp)
                            Y(P,J)%v(:,2:N_ph) = cmplx(mm(U_invN,Hr_re), 0.0_dp, kind=dp)
                        else
                            Hr_im = mm(H_m,aimag(R))
                            X(P,J)%v(:,2:N_ph) = cmplx(mm(U_invN,real(R,dp)), mm(U_invN,aimag(R)), kind=dp)
                            Y(P,J)%v(:,2:N_ph) = cmplx(mm(U_invN,Hr_re), mm(U_invN,Hr_im), kind=dp)
                        end if
                        call amplitudes(2)
                        deallocate(NU, lam, V_N, N_ort, inv_sqrtN, U_invN, AB_U, H_p, H_m, H, W, R, Hr_re)
                        if (allocated(Hr_im)) deallocate(Hr_im)
                    end if

                    ! 2x2 RPA problem in the spurious subspace ...
                    A_s = dot_product(Spur,matmul(A(P,J)%v,Spur))
                    B_s = dot_product(Spur,matmul(B(P,J)%v,Spur))
                    N_s = dot_product(Spur,Nd * Spur)
                    c_s = 1.0_dp / sqrt(N_s)
                    Hp_s = c_s * (A_s + B_s) * c_s
                    Hm_s = c_s * (A_s - B_s) * c_s
                    E_s = sqrt(cmplx(Hp_s * Hm_s, 0.0_dp, kind=dp))
                    r_s = c_s
                    p_s = c_s * Hm_s
                    x_s = 0.5_dp * p_s / E_s + 0.5_dp * r_s
                    y_s = 0.5_dp * p_s / E_s - 0.5_dp * r_s
                    nrm_s = N_s * (abs2(x_s) - abs2(y_s))
                    write(*,'(a)') TAB // TAB // TAB // 'Spurious 1- RPA Goldstone boson   ...   E = (' // &
                                   cstr(round_digits(real(E_s,dp),7),round_digits(aimag(E_s),7)) // ') MeV'
                    ! Robustness guard (dev leaves E(1) undefined for |E_s| or |nrm_s| == 1e-8 exactly): else-branch ...
                    if (abs(E_s) > 1.0e-8_dp .and. abs(nrm_s) > 1.0e-8_dp) then
                        E(P,J)%v(1) = E_s
                        X(P,J)%v(:,1) = ((1.0_dp / sqrt(abs(nrm_s))) * x_s) * Spur
                        Y(P,J)%v(:,1) = ((1.0_dp / sqrt(abs(nrm_s))) * y_s) * Spur
                    else
                        write(*,'(a)') TAB // TAB // TAB // 'Zero-norm spurious 1- level, setting E = 0 & X = Y = 0 ...'
                        E(P,J)%v(1) = (0.0_dp, 0.0_dp)
                        X(P,J)%v(:,1) = (0.0_dp, 0.0_dp)
                        Y(P,J)%v(:,1) = (0.0_dp, 0.0_dp)
                    end if
                    deallocate(Spur, U_sub)
                end if

                if (any(abs(aimag(E(P,J)%v)) > 1.0e-8_dp)) then
                    write(*,'(a)') TAB // TAB // TAB // 'RPA INSTABILITY DETECTED in channel:' // TAB // 'J = ' // istr(J) // &
                                   ', P = ' // pm_str(P)
                    Stability = .false.
                end if
                deallocate(Nd)
            end do
        end do
        write(*,'(a)') TAB // 'Finished diagonalization of the RPA system of equations ...'

    contains

        ! X = (p/E + r)/2, Y = (p/E - r)/2 from r = X - Y (stored in X) & p = E (X + Y) (stored in Y), then the metric
            ! normalization sum_ph N_ph (|X|^2 - |Y|^2) -> +-1 (X = Y = 0 for |norm| <= 1e-12) for the states nu >= nu_0 ...
        subroutine amplitudes(nu_0)
            integer, intent(in) :: nu_0
            integer :: nu, ph
            complex(dp) :: iE, rr, pp
            real(dp) :: Norm
            do nu = nu_0, N_ph
                iE = 1.0_dp / E(P,J)%v(nu)
                do ph = 1, N_ph
                    rr = X(P,J)%v(ph,nu)
                    pp = Y(P,J)%v(ph,nu)
                    X(P,J)%v(ph,nu) = (0.5_dp * iE) * pp + 0.5_dp * rr
                    Y(P,J)%v(ph,nu) = (0.5_dp * iE) * pp - 0.5_dp * rr
                end do
                Norm = 0.0_dp
                do ph = 1, N_ph
                    Norm = Norm + Nd(ph) * (abs2(X(P,J)%v(ph,nu)) - abs2(Y(P,J)%v(ph,nu)))
                end do
                if (abs(Norm) > 1.0e-12_dp) then
                    X(P,J)%v(:,nu) = X(P,J)%v(:,nu) * (1.0_dp / sqrt(abs(Norm)))
                    Y(P,J)%v(:,nu) = Y(P,J)%v(:,nu) * (1.0_dp / sqrt(abs(Norm)))
                else
                    X(P,J)%v(:,nu) = (0.0_dp, 0.0_dp)
                    Y(P,J)%v(:,nu) = (0.0_dp, 0.0_dp)
                    write(*,'(a)') TAB // TAB // TAB // 'Zero-norm RPA state detected, setting X = Y = 0 ...'
                end if
            end do
        end subroutine amplitudes

    end subroutine RPA_diagonalize

    ! Operator of isospin T (-1 proton part, +1 neutron part) ...
    pure real(dp) function op_t(Ops,k,a,b,T)
        type(o1b_t), intent(in) :: Ops(:)
        integer, intent(in) :: k, a, b, T
        if (T == -1) then
            op_t = Ops(k)%p(a,b)
        else
            op_t = Ops(k)%n(a,b)
        end if
    end function op_t

    pure real(dp) function abs2(z)
        complex(dp), intent(in) :: z
        abs2 = real(z,dp)**2 + aimag(z)**2
    end function abs2

    ! Julia string(::ComplexF64) ... "re + im im" / "re - im im" ...
    function cstr(re,im) result(s)
        real(dp), intent(in) :: re, im
        character(len=:), allocatable :: s
        if (im < 0.0_dp .or. (im == 0.0_dp .and. sign(1.0_dp,im) < 0.0_dp)) then
            s = jstr(re) // ' - ' // jstr(-im) // 'im'
        else
            s = jstr(re) // ' + ' // jstr(im) // 'im'
        end if
    end function cstr

    ! Diagonal RPA OBDM from the backward amplitudes, all states, weights (n_h - n_p)^2 ... counterpart of RPA_OBDM_oY2
        ! (RPA_rho.jl:1-73) ... per-JP buffers summed in the JP order (race-free & deterministic) ...
    function RPA_OBDM_oY2(Params,Orb,Orb_1p1h,Y) result(Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cmat_t), intent(in) :: Y(:,0:)
        type(o1b_t) :: Rho
        real(dp), allocatable :: Buf(:,:,:)
        integer :: a_max, N_JP, i_JP, J, P, ph, nu, a_p, a_h, j_p, j_h, a
        real(dp) :: S, n_p, n_h, Amp_p, Amp_h

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        N_JP = 2*(Orb_1p1h%J_max + 1)
        allocate(Rho%p(a_max,a_max), Rho%n(a_max,a_max), Buf(a_max,2,N_JP))
        Rho%p = 0.0_dp
        Rho%n = 0.0_dp
        do a = 1, a_max
            Rho%p(a,a) = Orb(a)%pO
            Rho%n(a,a) = Orb(a)%nO
        end do
        Buf = 0.0_dp

        !$omp parallel do schedule(dynamic) default(shared) private(J,P,ph,nu,a_p,a_h,j_p,j_h,S,n_p,n_h,Amp_p,Amp_h)
        do i_JP = 1, N_JP
            call JP_of(i_JP,J,P)
            do ph = 1, Orb_1p1h%N(P,J)
                a_p = Orb_1p1h%Ind(P,J)%v(ph)%a
                a_h = Orb_1p1h%Ind(P,J)%v(ph)%b
                j_p = Orb(a_p)%j
                j_h = Orb(a_h)%j
                S = 0.0_dp
                do nu = 1, Orb_1p1h%N(P,J)
                    S = S + abs2(Y(P,J)%v(ph,nu)) * hat_of(2*J)
                end do
                n_p = ph_occ(Orb,a_p,Orb_1p1h%Ind(P,J)%v(ph)%T)
                n_h = ph_occ(Orb,a_h,Orb_1p1h%Ind(P,J)%v(ph)%T)
                Amp_p = (n_h - n_p)**2 / hat_of(j_p)
                Amp_h = (n_h - n_p)**2 / hat_of(j_h)
                if (Orb_1p1h%Ind(P,J)%v(ph)%T == -1) then
                    Buf(a_h,1,i_JP) = Buf(a_h,1,i_JP) - Amp_h * S
                    Buf(a_p,1,i_JP) = Buf(a_p,1,i_JP) + Amp_p * S
                else
                    Buf(a_h,2,i_JP) = Buf(a_h,2,i_JP) - Amp_h * S
                    Buf(a_p,2,i_JP) = Buf(a_p,2,i_JP) + Amp_p * S
                end if
            end do
        end do
        !$omp end parallel do

        do i_JP = 1, N_JP
            do a = 1, a_max
                Rho%p(a,a) = Rho%p(a,a) + Buf(a,1,i_JP)
                Rho%n(a,a) = Rho%n(a,a) + Buf(a,2,i_JP)
            end do
        end do
    end function RPA_OBDM_oY2

    ! RPA ground-state correlation energies ... counterpart of RPA_energy (RPA_energy.jl) ... Density defaults to RPA.TBCF;
        ! Have_Sigma_AB: A & B contain the Sigma_AB box diagrams (their energy E_dress is evaluated explicitly) ...
    function RPA_energy(Params,Orb,Orb_1p1h,Orb_NN,E,X,Y,Rho,Rho_ref,h_N,Sigma,V_NN,Have_Sigma_AB,Sigma_AB,Density) result(E_corr)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(cvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(o1b_t), intent(in) :: Rho, Rho_ref, h_N
        type(o2b_t), intent(in) :: Sigma, V_NN, Sigma_AB
        logical, intent(in) :: Have_Sigma_AB
        logical, intent(in), optional :: Density
        type(rpa_energy_t) :: E_corr
        type(rmat_t), allocatable :: dA(:,:), dB(:,:)
        logical :: Dens
        real(dp) :: NaN
        integer :: J, P

        Dens = Params%Calc%RPA%TBCF
        if (present(Density)) Dens = Density
        NaN = ieee_value(1.0_dp, ieee_quiet_nan)

        write(*,'(/,a)') 'Calculating RPA correlation energy ...'
        E_corr%E_QBA = RPA_energy_QBA(Params,Orb,Orb_1p1h,E,X,Y)
        E_corr%E_sigma = NaN
        E_corr%E_rho = NaN
        E_corr%E_dress = NaN
        if (Dens) then
            call RPA_energy_rho(Params,Orb,Orb_NN,Rho,Rho_ref,h_N,Sigma,V_NN,E_corr%E_sigma,E_corr%E_rho)
            E_corr%E_dress = 0.0_dp
            if (Have_Sigma_AB) then
                allocate(dA(2,0:Orb_1p1h%J_max), dB(2,0:Orb_1p1h%J_max))
                do J = 0, Orb_1p1h%J_max
                    do P = 1, 2
                        allocate(dA(P,J)%v(Orb_1p1h%N(P,J),Orb_1p1h%N(P,J)), dB(P,J)%v(Orb_1p1h%N(P,J),Orb_1p1h%N(P,J)))
                        dA(P,J)%v = 0.0_dp
                        dB(P,J)%v = 0.0_dp
                    end do
                end do
                call ph_allocate_tbcf(Params,Orb,Orb_1p1h,Orb_NN,V_NN,Sigma_AB,.true.,Params%Calc%RPA%DEV,dA,dB,'RPA')
                ! DEV: 1p-1h to 2p-2h self-energy coupling of dA (RPA_allocate_SE, Mode = "OnShell") ...
                if (Params%Calc%RPA%DEV) call RPA_allocate_SE(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,dA)
                E_corr%E_dress = RPA_energy_dress(Params,Orb,Orb_1p1h,E,X,Y,dA,dB)
            end if
        end if
        E_corr%E_dens = E_corr%E_sigma + E_corr%E_rho + E_corr%E_dress

        write(*,'(a)') TAB // 'RPA correlation energy (Quasi-boson ... -sum w Y^+ N Y) ...                  E_QBA = ' // &
                       jstr(round_digits(E_corr%E_QBA,6)) // ' MeV'
        if (Dens) then
            write(*,'(a)') TAB // 'RPA correlation energy (Correlated 2-body correlation function) ...          E_sig = ' // &
                           jstr(round_digits(E_corr%E_sigma,6)) // ' MeV'
            write(*,'(a)') TAB // 'RPA correlation energy (Correlated 1-body density, h.drho) ...               E_rho = ' // &
                           jstr(round_digits(E_corr%E_rho,6)) // ' MeV'
            if (Have_Sigma_AB) &
                write(*,'(a)') TAB // 'RPA correlation energy (Sigma terms of A & B, Y^+ dA Y + X^+ dB Y) ...       D     = ' // &
                               jstr(round_digits(E_corr%E_dress,6)) // ' MeV'
            write(*,'(a)') TAB // 'RPA correlation energy (Density, E_sig + E_rho + D) ...                      E_RPA = ' // &
                           jstr(round_digits(E_corr%E_dens,6)) // ' MeV'
        end if
        write(*,'(a)') ''
    end function RPA_energy

    ! E_QBA = -sum_(J,P) (2J+1) sum_nu Re(E_nu) sum_ph D_ph |Y_ph,nu|^2 over the states with non-zero norm & real energy ...
    real(dp) function RPA_energy_QBA(Params,Orb,Orb_1p1h,E,X,Y) result(E_QBA)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        real(dp), allocatable :: E_JP(:), D(:)
        integer, allocatable :: N_drop(:)
        integer :: N_JP, i_JP, J, P, nu, ph, N_ph
        real(dp) :: SumE, Norm

        N_JP = 2*(Orb_1p1h%J_max + 1)
        allocate(E_JP(N_JP), N_drop(N_JP))
        E_JP = 0.0_dp
        N_drop = 0
        !$omp parallel do schedule(dynamic) default(shared) private(J,P,nu,ph,N_ph,SumE,Norm,D)
        do i_JP = 1, N_JP
            call JP_of(i_JP,J,P)
            N_ph = Orb_1p1h%N(P,J)
            if (N_ph == 0) cycle
            D = ph_metric(Orb,Orb_1p1h,P,J)
            SumE = 0.0_dp
            do nu = 1, N_ph
                Norm = 0.0_dp
                do ph = 1, N_ph
                    Norm = Norm + D(ph) * (abs2(X(P,J)%v(ph,nu)) - abs2(Y(P,J)%v(ph,nu)))
                end do
                if (abs(Norm) < 1.0e-8_dp) cycle
                if (abs(aimag(E(P,J)%v(nu))) > 1.0e-8_dp) then
                    N_drop(i_JP) = N_drop(i_JP) + 1
                    cycle
                end if
                do ph = 1, N_ph
                    SumE = SumE - hat_of(2*J) * real(E(P,J)%v(nu),dp) * D(ph) * abs2(Y(P,J)%v(ph,nu))
                end do
            end do
            E_JP(i_JP) = SumE
        end do
        !$omp end parallel do
        if (sum(N_drop) > 0) then
            write(*,'(a)') TAB // 'Warning! ' // istr(sum(N_drop)) // ' RPA solutions with complex energies (instability) are ' // &
                           'excluded from the quasi-boson energy ...'
            write(*,'(a)') TAB // '         (they are still contained in the 2-body correlation function Sigma ...)'
        end if
        E_QBA = sum(E_JP)
    end function RPA_energy_QBA

    ! E_dress = sum_(J,P) (2J+1) sum_nu Re(Y^+ dA Y + X^+ dB Y) over the same physical states as E_QBA ...
    real(dp) function RPA_energy_dress(Params,Orb,Orb_1p1h,E,X,Y,dA,dB) result(E_dress)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(rmat_t), intent(in) :: dA(:,0:), dB(:,0:)
        real(dp), allocatable :: E_JP(:), D(:)
        complex(dp), allocatable :: dAY(:), dBY(:)
        integer :: N_JP, i_JP, J, P, nu, ph, N_ph
        real(dp) :: SumD, Norm

        N_JP = 2*(Orb_1p1h%J_max + 1)
        allocate(E_JP(N_JP))
        E_JP = 0.0_dp
        !$omp parallel do schedule(dynamic) default(shared) private(J,P,nu,ph,N_ph,SumD,Norm,D,dAY,dBY)
        do i_JP = 1, N_JP
            call JP_of(i_JP,J,P)
            N_ph = Orb_1p1h%N(P,J)
            if (N_ph == 0) cycle
            D = ph_metric(Orb,Orb_1p1h,P,J)
            SumD = 0.0_dp
            do nu = 1, N_ph
                Norm = 0.0_dp
                do ph = 1, N_ph
                    Norm = Norm + D(ph) * (abs2(X(P,J)%v(ph,nu)) - abs2(Y(P,J)%v(ph,nu)))
                end do
                if (abs(Norm) < 1.0e-8_dp) cycle
                if (abs(aimag(E(P,J)%v(nu))) > 1.0e-8_dp) cycle
                dAY = matmul(dA(P,J)%v,Y(P,J)%v(:,nu))
                dBY = matmul(dB(P,J)%v,Y(P,J)%v(:,nu))
                SumD = SumD + hat_of(2*J) * real(dot_product(Y(P,J)%v(:,nu),dAY) + dot_product(X(P,J)%v(:,nu),dBY),dp)
            end do
            E_JP(i_JP) = SumD
        end do
        !$omp end parallel do
        E_dress = sum(E_JP)
    end function RPA_energy_dress

    ! E_sigma = V2b_correlation_energy(Sigma,V_NN), E_rho = sum_ab (2j_a+1) h_N_ab (Rho - Rho_ref)_ba (equal l & j) ...
    subroutine RPA_energy_rho(Params,Orb,Orb_NN,Rho,Rho_ref,h_N,Sigma,V_NN,E_sigma,E_rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: Rho, Rho_ref, h_N
        type(o2b_t), intent(in) :: Sigma, V_NN
        real(dp), intent(out) :: E_sigma, E_rho
        integer :: ia, ib, a, b

        E_sigma = V2b_correlation_energy(Params,Orb,Orb_NN,Sigma,V_NN)
        E_rho = 0.0_dp
        do ia = 1, size(Orb)
            a = Orb(ia)%a
            do ib = 1, size(Orb)
                b = Orb(ib)%a
                if (Orb(ib)%l /= Orb(ia)%l .or. Orb(ib)%j /= Orb(ia)%j) cycle
                E_rho = E_rho + hat_of(Orb(ia)%j) * (h_N%p(a,b) * (Rho%p(b,a) - Rho_ref%p(b,a)) + &
                                                     h_N%n(a,b) * (Rho%n(b,a) - Rho_ref%n(b,a)))
            end do
        end do
    end subroutine RPA_energy_rho

    ! Reduced transition matrix elements rM = M_X X + M_Y Y ... counterpart of RPA_rM, RPA_rM_calculate & RPA_rM_contract
        ! (RPA_transitions.jl:1-238) ... Transitions = 'Linear' (<0|[Q,M]|0> with the occupations of Orb) or 'Exact'
        ! (<0|Q M|0> with the OBDM Rho & for RPA.TBCF the RPA 2-body correlation function Sigma) ...
    function RPA_rM(Params,Orb,Orb_1p1h,Orb_NN,M,Rho,Sigma,X,Y) result(rM)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(tr1b_t), intent(in) :: M
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(in) :: Sigma
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(reduced_multipole_t) :: rM
        complex(dp), allocatable :: pE0(:,:), nE0(:,:), pE1(:,:), nE1(:,:), pE2(:,:), nE2(:,:), pE3(:,:), nE3(:,:)
        character(len=:), allocatable :: Transitions
        logical :: Exact, TBCF

        Transitions = trim(Params%Calc%RPA%Transitions)
        if (Transitions /= 'Exact' .and. Transitions /= 'Linear') &
            call fatal('Invalid value of "Params.Calc.RPA.Transitions" = "' // Transitions // '"; use "Exact" or "Linear".')
        Exact = (Transitions == 'Exact')
        TBCF = Params%Calc%RPA%TBCF

        write(*,'(/,a)') 'Calculating reduced transition matrix elements rM ... (' // Transitions // ')'
        call calc(0,1,ph_operators(M,0),pE0,nE0)
        call calc(1,2,ph_operators(M,1),pE1,nE1)
        call calc(2,1,ph_operators(M,2),pE2,nE2)
        call calc(3,2,ph_operators(M,3),pE3,nE3)
        call ph_rM_assemble(pE0,nE0,pE1,nE1,pE2,nE2,pE3,nE3,rM)
        write(*,'(a)') TAB // 'Reduced 1-body transition matrix elements rM calculated ...'

    contains

        subroutine calc(J,P,Ops,prM,nrM)
            integer, intent(in) :: J, P
            type(o1b_t), intent(in) :: Ops(:)
            complex(dp), allocatable, intent(out) :: prM(:,:), nrM(:,:)
            real(dp), allocatable :: pM_X(:,:), pM_Y(:,:), nM_X(:,:), nM_Y(:,:), Eye(:,:)
            type(o1b_t), allocatable :: G(:)
            real(dp), allocatable :: MX1(:), MY1(:), MX2(:), MY2(:)
            integer :: N_ph, N_o, i_o, i_ph, a_p, a_h, T_ph, l_p, j_p, l_h, j_h, ia, ib, a_a, l_a, j_a, a_b, l_b, j_b
            integer :: P_ha, I, a
            real(dp) :: Phase_hp, D_ph, Phase_X, Phase_Y, K_1, K_2, F, Phase_ab

            N_ph = Orb_1p1h%N(P,J)
            N_o = size(Ops)
            allocate(pM_X(N_o,N_ph), pM_Y(N_o,N_ph), nM_X(N_o,N_ph), nM_Y(N_o,N_ph), G(N_o))
            pM_X = 0.0_dp
            pM_Y = 0.0_dp
            nM_X = 0.0_dp
            nM_Y = 0.0_dp

            ! 1-body transition densities G = (1 - Rho) M Rho (Exact) ...
            if (Exact) then
                allocate(Eye(size(Rho%p,1),size(Rho%p,1)))
                Eye = 0.0_dp
                do a = 1, size(Eye,1)
                    Eye(a,a) = 1.0_dp
                end do
                do i_o = 1, N_o
                    G(i_o)%p = mm(mm(Eye - Rho%p,Ops(i_o)%p),Rho%p)
                    G(i_o)%n = mm(mm(Eye - Rho%n,Ops(i_o)%n),Rho%n)
                end do
            end if

            !$omp parallel do schedule(dynamic) default(shared) &
            !$omp private(a_p,a_h,T_ph,l_p,j_p,l_h,j_h,Phase_hp,MX1,MY1,MX2,MY2,D_ph,i_o,Phase_X,Phase_Y,ia,ib,a_a,l_a,j_a) &
            !$omp private(a_b,l_b,j_b,P_ha,K_1,K_2,I,F,Phase_ab)
            do i_ph = 1, N_ph
                a_p = Orb_1p1h%Ind(P,J)%v(i_ph)%a
                a_h = Orb_1p1h%Ind(P,J)%v(i_ph)%b
                T_ph = Orb_1p1h%Ind(P,J)%v(i_ph)%T
                l_p = Orb(a_p)%l
                j_p = Orb(a_p)%j
                l_h = Orb(a_h)%l
                j_h = Orb(a_h)%j
                Phase_hp = phase((j_p - j_h)/2)
                allocate(MX1(N_o), MY1(N_o), MX2(N_o), MY2(N_o))
                MX1 = 0.0_dp
                MY1 = 0.0_dp
                MX2 = 0.0_dp
                MY2 = 0.0_dp

                if (.not. Exact) then
                    D_ph = ph_occ(Orb,a_h,T_ph) - ph_occ(Orb,a_p,T_ph)
                    do i_o = 1, N_o
                        MX1(i_o) = op_t(Ops,i_o,a_p,a_h,T_ph) * D_ph * phase(J)
                        MY1(i_o) = op_t(Ops,i_o,a_h,a_p,T_ph) * D_ph * Phase_hp
                    end do
                else
                    do i_o = 1, N_o
                        if (T_ph == -1) then
                            MX1(i_o) = G(i_o)%p(a_p,a_h) * phase(J)
                            MY1(i_o) = -G(i_o)%p(a_h,a_p) * Phase_hp
                        else
                            MX1(i_o) = G(i_o)%n(a_p,a_h) * phase(J)
                            MY1(i_o) = -G(i_o)%n(a_h,a_p) * Phase_hp
                        end if
                    end do
                    if (TBCF) then
                        Phase_X = phase(1 + (j_p + j_h)/2)
                        Phase_Y = phase(J + (j_p + j_h)/2)
                        do ia = 1, size(Orb)
                            a_a = Orb(ia)%a
                            l_a = Orb(ia)%l
                            j_a = Orb(ia)%j
                            P_ha = parity_of(l_h + l_a)
                            do ib = 1, size(Orb)
                                a_b = Orb(ib)%a
                                l_b = Orb(ib)%l
                                j_b = Orb(ib)%j
                                if ((.not. triangle(j_a,j_b,2*J)) .or. parity_of(l_p + l_b) /= P_ha) cycle
                                K_1 = 0.0_dp
                                K_2 = 0.0_dp
                                do I = max(abs(j_h - j_a),abs(j_p - j_b))/2, min(j_h + j_a,j_p + j_b)/2
                                    F = phase(J + I + (j_p + j_a)/2) * hat_of(2*I) * f6j(j_p,j_h,2*J,j_a,j_b,2*I)
                                    if (T_ph == -1) then
                                        K_1 = K_1 + O2b_pp(a_h,a_a,a_p,a_b,I,P_ha,Sigma,Orb,Orb_NN) * F
                                        K_2 = K_2 + O2b_pn(a_h,a_a,a_p,a_b,I,P_ha,Sigma,Orb_NN) * F
                                    else
                                        K_1 = K_1 + O2b_nn(a_h,a_a,a_p,a_b,I,P_ha,Sigma,Orb,Orb_NN) * F
                                        K_2 = K_2 + O2b_pn(a_a,a_h,a_b,a_p,I,P_ha,Sigma,Orb_NN) * phase((j_p + j_a + j_h + j_b)/2) * F
                                    end if
                                end do
                                Phase_ab = phase((j_a - j_b)/2)
                                do i_o = 1, N_o
                                    MX1(i_o) = MX1(i_o) + Phase_X * op_t(Ops,i_o,a_a,a_b,T_ph) * K_1
                                    MX2(i_o) = MX2(i_o) + Phase_X * op_t(Ops,i_o,a_a,a_b,-T_ph) * K_2
                                    MY1(i_o) = MY1(i_o) + Phase_Y * Phase_ab * op_t(Ops,i_o,a_b,a_a,T_ph) * K_1
                                    MY2(i_o) = MY2(i_o) + Phase_Y * Phase_ab * op_t(Ops,i_o,a_b,a_a,-T_ph) * K_2
                                end do
                            end do
                        end do
                    end if
                end if

                ! Opposite-isospin contributions go to the opposite-isospin part of rM ...
                if (T_ph == -1) then
                    pM_X(:,i_ph) = MX1
                    pM_Y(:,i_ph) = MY1
                    nM_X(:,i_ph) = MX2
                    nM_Y(:,i_ph) = MY2
                else
                    nM_X(:,i_ph) = MX1
                    nM_Y(:,i_ph) = MY1
                    pM_X(:,i_ph) = MX2
                    pM_Y(:,i_ph) = MY2
                end if
                deallocate(MX1, MY1, MX2, MY2)
            end do
            !$omp end parallel do

            prM = zmm(cmplx(pM_X,0.0_dp,kind=dp),X(P,J)%v) + zmm(cmplx(pM_Y,0.0_dp,kind=dp),Y(P,J)%v)
            nrM = zmm(cmplx(nM_X,0.0_dp,kind=dp),X(P,J)%v) + zmm(cmplx(nM_Y,0.0_dp,kind=dp),Y(P,J)%v)

        end subroutine calc

    end function RPA_rM

    ! RPA 2-body correlation function Sigma from X & Y ... counterpart of RPA_sigma (RPA_sigma.jl) with fix F2 (the pn
        ! hole-particle-particle-hole term uses X, as the pp & nn branches; dev passes Y) ... particle O <= Tol_par, hole
        ! O >= Tol_hol (inclusive, unlike orbitals_1p1h_make; ported as is) ... parallel over the (Bra,Ket,T) list of each
        ! channel, every element written by one task ...
    function RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y,beta,Report) result(Sigma_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        real(dp), intent(in), optional :: beta
        logical, intent(in), optional :: Report
        type(o2b_t) :: Sigma_NN
        real(dp) :: pTol, hTol, bt
        logical :: Rep
        integer :: J, P, J_top, N_T0, N_T1, N_list, k, Bra, Ket, T, Ind, a, b, c, d, j_a, j_b, j_c, j_d, Q, Qx
        integer :: t_a, t_b, t_c, t_d
        integer, allocatable :: List(:,:)
        real(dp) :: O_a, O_b, O_c, O_d, Sum, Ph_ab, Ph_cd

        bt = 1.0_dp
        if (present(beta)) bt = beta
        Rep = .true.
        if (present(Report)) Rep = Report
        hTol = Params%Calc%RPA%Tol_hol
        pTol = Params%Calc%RPA%Tol_par

        call O2b_initialize(Params,Orb,Sigma_NN)
        if (Rep) write(*,'(/,a)') 'Allocating 2-body RPA correlation function matrix Sigma ...'

        J_top = min(2*Params%Calc%Nmax + 1, ubound(Orb_NN%N,3))
        do J = 0, J_top
            do P = 1, 2
                if (Rep) write(*,'(a)') TAB // 'Calculating   ...   J = ' // istr(J) // '/' // istr(2*Params%Calc%Nmax + 1) // &
                                        TAB // 'P = ' // pm_str(P)
                N_T0 = Orb_NN%N(1,P,J)
                N_T1 = Orb_NN%N(2,P,J)
                N_list = N_T0*(N_T0 + 1)/2 + N_T1*(N_T1 + 1)
                allocate(List(3,N_list))
                N_list = 0
                do Bra = 1, N_T1
                    do Ket = 1, Bra
                        N_list = N_list + 1
                        List(:,N_list) = [Bra, Ket, -1]
                    end do
                end do
                do Bra = 1, N_T0
                    do Ket = 1, Bra
                        N_list = N_list + 1
                        List(:,N_list) = [Bra, Ket, 0]
                    end do
                end do
                do Bra = 1, N_T1
                    do Ket = 1, Bra
                        N_list = N_list + 1
                        List(:,N_list) = [Bra, Ket, 1]
                    end do
                end do

                !$omp parallel do schedule(dynamic,64) default(shared) &
                !$omp private(Bra,Ket,T,Ind,a,b,c,d,j_a,j_b,j_c,j_d,Q,Qx,t_a,t_b,t_c,t_d,O_a,O_b,O_c,O_d,Sum,Ph_ab,Ph_cd)
                do k = 1, N_list
                    Bra = List(1,k)
                    Ket = List(2,k)
                    T = List(3,k)
                    Sum = 0.0_dp
                    if (T == 0) then
                        Ind = Bra + (Ket - 1) * N_T0 - (Ket * (Ket - 1))/2
                        a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                        b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                        c = Orb_NN%Ind(1,P,J)%v(1,Ket)
                        d = Orb_NN%Ind(1,P,J)%v(2,Ket)
                    else
                        Ind = Bra + (Ket - 1) * N_T1 - (Ket * (Ket - 1))/2
                        a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                        b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                        c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                        d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                    end if
                    j_a = Orb(a)%j
                    j_b = Orb(b)%j
                    j_c = Orb(c)%j
                    j_d = Orb(d)%j
                    Q = mod(Orb(a)%l + Orb(c)%l,2) + 1
                    Qx = mod(Orb(a)%l + Orb(d)%l,2) + 1
                    if (Q == mod(Orb(b)%l + Orb(d)%l,2) + 1) then
                        if (T == 0) then
                            O_a = Orb(a)%pO
                            O_b = Orb(b)%nO
                            O_c = Orb(c)%pO
                            O_d = Orb(d)%nO
                            t_a = -1
                            t_b = 1
                            t_c = -1
                            t_d = 1
                            ! 2p-2h terms ...
                            if (part(O_a) .and. part(O_b) .and. hole(O_c) .and. hole(O_d)) &
                                Sum = Sum + s_pphh(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q)
                            if (hole(O_a) .and. hole(O_b) .and. part(O_c) .and. part(O_d)) &
                                Sum = Sum + s_pphh(c,t_c,O_c,d,t_d,O_d,a,t_a,O_a,b,t_b,O_b,Q)
                            ! 1p-1h terms ... (hpph with Y in every branch, fix F2)
                            if (part(O_a) .and. hole(O_b) .and. hole(O_c) .and. part(O_d)) &
                                Sum = Sum + s_phhp(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q)
                            if (hole(O_a) .and. part(O_b) .and. part(O_c) .and. hole(O_d)) &
                                Sum = Sum + s_hpph(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q)
                        else
                            if (T == -1) then
                                O_a = Orb(a)%pO
                                O_b = Orb(b)%pO
                                O_c = Orb(c)%pO
                                O_d = Orb(d)%pO
                            else
                                O_a = Orb(a)%nO
                                O_b = Orb(b)%nO
                                O_c = Orb(c)%nO
                                O_d = Orb(d)%nO
                            end if
                            t_a = T
                            t_b = T
                            t_c = T
                            t_d = T
                            ! 2p-2h terms ... antisymmetrized projection 1/4 (1 - P_ab)(1 - P_cd) of the direct term ...
                            if (part(O_a) .and. part(O_b) .and. hole(O_c) .and. hole(O_d)) then
                                Ph_ab = phase(J + (j_a + j_b)/2)
                                Ph_cd = phase(J + (j_c + j_d)/2)
                                Sum = Sum + 0.25_dp * (s_pphh(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q) - Ph_cd * &
                                                       s_pphh(a,t_a,O_a,b,t_b,O_b,d,t_d,O_d,c,t_c,O_c,Qx) - Ph_ab * &
                                                       s_pphh(b,t_b,O_b,a,t_a,O_a,c,t_c,O_c,d,t_d,O_d,Qx) + Ph_ab * Ph_cd * &
                                                       s_pphh(b,t_b,O_b,a,t_a,O_a,d,t_d,O_d,c,t_c,O_c,Q))
                            end if
                            if (hole(O_a) .and. hole(O_b) .and. part(O_c) .and. part(O_d)) then
                                Ph_ab = phase(J + (j_a + j_b)/2)
                                Ph_cd = phase(J + (j_c + j_d)/2)
                                Sum = Sum + 0.25_dp * (s_pphh(c,t_c,O_c,d,t_d,O_d,a,t_a,O_a,b,t_b,O_b,Q) - Ph_cd * &
                                                       s_pphh(d,t_d,O_d,c,t_c,O_c,a,t_a,O_a,b,t_b,O_b,Qx) - Ph_ab * &
                                                       s_pphh(c,t_c,O_c,d,t_d,O_d,b,t_b,O_b,a,t_a,O_a,Qx) + Ph_ab * Ph_cd * &
                                                       s_pphh(d,t_d,O_d,c,t_c,O_c,b,t_b,O_b,a,t_a,O_a,Q))
                            end if
                            ! 1p-1h terms ...
                            if (part(O_a) .and. hole(O_b) .and. part(O_c) .and. hole(O_d)) &
                                Sum = Sum + (s_phhp(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q) - phase(J + (j_c + j_d)/2) * &
                                             s_phhp(a,t_a,O_a,b,t_b,O_b,d,t_d,O_d,c,t_c,O_c,Qx))
                            if (part(O_a) .and. hole(O_b) .and. hole(O_c) .and. part(O_d)) &
                                Sum = Sum + (s_phhp(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q) - phase(J + (j_c + j_d)/2) * &
                                             s_phhp(a,t_a,O_a,b,t_b,O_b,d,t_d,O_d,c,t_c,O_c,Qx))
                            if (hole(O_a) .and. part(O_b) .and. hole(O_c) .and. part(O_d)) &
                                Sum = Sum + (s_hpph(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q) - phase(J + (j_c + j_d)/2) * &
                                             s_hpph(a,t_a,O_a,b,t_b,O_b,d,t_d,O_d,c,t_c,O_c,Qx))
                            if (hole(O_a) .and. part(O_b) .and. part(O_c) .and. hole(O_d)) &
                                Sum = Sum + (s_hpph(a,t_a,O_a,b,t_b,O_b,c,t_c,O_c,d,t_d,O_d,Q) - phase(J + (j_c + j_d)/2) * &
                                             s_hpph(a,t_a,O_a,b,t_b,O_b,d,t_d,O_d,c,t_c,O_c,Qx))
                        end if
                    end if
                    select case (T)
                    case (-1)
                        Sigma_NN%pp(P,J)%v(Ind) = bt * Sum
                    case (0)
                        Sigma_NN%pn(P,J)%v(Ind) = bt * Sum
                    case default
                        Sigma_NN%nn(P,J)%v(Ind) = bt * Sum
                    end select
                end do
                !$omp end parallel do
                deallocate(List)
            end do
        end do
        if (Rep) write(*,'(a)') TAB // 'Succesfully allocated 2-body RPA correlation matrix Sigma ...'

    contains

        pure logical function part(O)
            real(dp), intent(in) :: O
            part = (O <= pTol)
        end function part

        pure logical function hole(O)
            real(dp), intent(in) :: O
            hole = (O >= hTol)
        end function hole

        ! Common reference (uncorrelated) terms of RPA_sigma_pphh / phhp / hpph ...
        pure real(dp) function s_ref(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,a_d,t_d)
            integer, intent(in) :: a_a, t_a, a_b, t_b, a_c, t_c, a_d, t_d
            real(dp), intent(in) :: O_a, O_b
            s_ref = 0.0_dp
            if (a_a == a_c .and. a_b == a_d) s_ref = s_ref - O_a * O_b
            if ((a_a == a_d .and. t_a == t_d) .and. (a_b == a_c .and. t_b == t_c)) &
                s_ref = s_ref - phase(J + (Orb(a_a)%j + Orb(a_b)%j)/2) * O_a
        end function s_ref

        ! RPA_sigma_pphh ... sum over the X-Y products (Y_ac X_bd*) ...
        real(dp) function s_pphh(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,O_c,a_d,t_d,O_d,PP) result(sg)
            integer, intent(in) :: a_a, t_a, a_b, t_b, a_c, t_c, a_d, t_d, PP
            real(dp), intent(in) :: O_a, O_b, O_c, O_d
            integer :: I, i_ac, i_bd, nu, ja, jb, jc, jd
            real(dp) :: Amp
            ja = Orb(a_a)%j
            jb = Orb(a_b)%j
            jc = Orb(a_c)%j
            jd = Orb(a_d)%j
            sg = 0.0_dp
            do I = max(abs(ja - jc)/2,abs(jb - jd)/2), min((ja + jc)/2,(jb + jd)/2)
                i_ac = orb2b_1p1h_index(a_a,a_c,I,PP,(t_a + t_c)/2,Orb_1p1h)
                i_bd = orb2b_1p1h_index(a_b,a_d,I,PP,(t_b + t_d)/2,Orb_1p1h)
                if (i_ac == 0 .or. i_bd == 0) cycle
                Amp = real(2*I + 1,dp) * phase(J + I + (jb + jc)/2) * f6j(ja,jb,2*J,jd,jc,2*I) * abs((O_c - O_a) * (O_d - O_b))
                do nu = 1, Orb_1p1h%N(PP,I)
                    sg = sg + Amp * real(Y(PP,I)%v(i_ac,nu) * conjg(X(PP,I)%v(i_bd,nu)),dp)
                end do
            end do
            sg = sg + s_ref(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,a_d,t_d)
        end function s_pphh

        ! RPA_sigma_phhp ... (Y_ac Y_db*) ...
        real(dp) function s_phhp(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,O_c,a_d,t_d,O_d,PP) result(sg)
            integer, intent(in) :: a_a, t_a, a_b, t_b, a_c, t_c, a_d, t_d, PP
            real(dp), intent(in) :: O_a, O_b, O_c, O_d
            integer :: I, i_ac, i_db, nu, ja, jb, jc, jd
            real(dp) :: Amp
            ja = Orb(a_a)%j
            jb = Orb(a_b)%j
            jc = Orb(a_c)%j
            jd = Orb(a_d)%j
            sg = 0.0_dp
            do I = max(abs(ja - jc)/2,abs(jb - jd)/2), min((ja + jc)/2,(jb + jd)/2)
                i_ac = orb2b_1p1h_index(a_a,a_c,I,PP,(t_a + t_c)/2,Orb_1p1h)
                i_db = orb2b_1p1h_index(a_d,a_b,I,PP,(t_b + t_d)/2,Orb_1p1h)
                if (i_ac == 0 .or. i_db == 0) cycle
                Amp = real(2*I + 1,dp) * phase(J + (jc + jd)/2) * f6j(ja,jb,2*J,jd,jc,2*I) * abs((O_c - O_a) * (O_b - O_d))
                do nu = 1, Orb_1p1h%N(PP,I)
                    sg = sg + Amp * real(Y(PP,I)%v(i_ac,nu) * conjg(Y(PP,I)%v(i_db,nu)),dp)
                end do
            end do
            sg = sg + s_ref(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,a_d,t_d)
        end function s_phhp

        ! RPA_sigma_hpph ... (Y_ca Y_bd*) ... fix F2: dev passes Y in the pn branch but X in the pp/nn branches; only Y keeps
            ! Sigma symmetric, <h p|Sigma|p' h'> must equal the transposed phhp element <p' h'|Sigma|h p> (Y_p'h Y_ph'*), and
            ! only the Ket <= Bra triangle is computed (docs/dev_port/upstream_reports/1_rpa_sigma_pp_nn_hpph.md) ...
        real(dp) function s_hpph(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,O_c,a_d,t_d,O_d,PP) result(sg)
            integer, intent(in) :: a_a, t_a, a_b, t_b, a_c, t_c, a_d, t_d, PP
            real(dp), intent(in) :: O_a, O_b, O_c, O_d
            integer :: I, i_ca, i_bd, nu, ja, jb, jc, jd
            real(dp) :: Amp
            ja = Orb(a_a)%j
            jb = Orb(a_b)%j
            jc = Orb(a_c)%j
            jd = Orb(a_d)%j
            sg = 0.0_dp
            do I = max(abs(ja - jc)/2,abs(jb - jd)/2), min((ja + jc)/2,(jb + jd)/2)
                i_ca = orb2b_1p1h_index(a_c,a_a,I,PP,(t_a + t_c)/2,Orb_1p1h)
                i_bd = orb2b_1p1h_index(a_b,a_d,I,PP,(t_b + t_d)/2,Orb_1p1h)
                if (i_ca == 0 .or. i_bd == 0) cycle
                Amp = real(2*I + 1,dp) * phase(J + (ja + jb)/2) * f6j(ja,jb,2*J,jd,jc,2*I) * abs((O_a - O_c) * (O_b - O_d))
                do nu = 1, Orb_1p1h%N(PP,I)
                    sg = sg + Amp * real(Y(PP,I)%v(i_ca,nu) * conjg(Y(PP,I)%v(i_bd,nu)),dp)
                end do
            end do
            sg = sg + s_ref(a_a,t_a,O_a,a_b,t_b,O_b,a_c,t_c,a_d,t_d)
        end function s_hpph

    end function RPA_sigma

end module nj_rpa
