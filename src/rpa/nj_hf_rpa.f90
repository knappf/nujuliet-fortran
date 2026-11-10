! Spherical HF-TDA & HF-RPA ... counterpart of HF_RPA/HF_RPA_Phonon_Count.jl, HF_RPA_Allocate.jl, HF_RPA_Spurious.jl,
    ! HF_RPA_Diagonalize.jl, HF_RPA_Energy.jl, HF_RPA_OBDM.jl & HF_RPA_Transitions.jl ...
    ! Phonon channels are indexed (J,P) with J = 0 .. 2*Nmax+1 ...
module nj_hf_rpa
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_tr1b, only : tr1b_t
    implicit none
    private

    public :: HF_RPA_phonon_count, HF_RPA_allocate, Spur_Ini, Spur_Ortho, HF_RPA_diagonalize, RPA_eigen_solve
    public :: RPA_normalize, RPA_eigen_solve_c, HF_RPA_energy, HF_RPA_OBDM, HF_RPA_rM, HF_RPA_rB, ph_orbitals

contains

    ! Orbital indices a_p & a_h of the particle-hole phonon ph ...
    pure subroutine ph_orbitals(ph,Phonon,Particle,Hole,a_p,a_h)
        integer, intent(in) :: ph
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        integer, intent(out) :: a_p, a_h
        if (Phonon(ph)%tz == -1) then
            a_p = Particle%p(Phonon(ph)%p)%a
            a_h = Hole%p(Phonon(ph)%h)%a
        else
            a_p = Particle%n(Phonon(ph)%p)%a
            a_h = Hole%n(Phonon(ph)%h)%a
        end if
    end subroutine ph_orbitals

    subroutine HF_RPA_phonon_count(Params,N_Phonon,Phonon,N_nu,Orb_Phonon)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_Phonon
        type(phstate_t), intent(in) :: Phonon(:)
        integer, allocatable, intent(out) :: N_nu(:,:)
        type(ivec_t), allocatable, intent(out) :: Orb_Phonon(:,:)
        integer :: J_max, J, P, ph, N_ph

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Counting all 1p-1h phonon states ...'
        allocate(N_nu(0:J_max,2), Orb_Phonon(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                N_ph = count(Phonon(1:N_Phonon)%J == J .and. Phonon(1:N_Phonon)%Par == P)
                N_nu(J,P) = N_ph
                allocate(Orb_Phonon(J,P)%v(N_ph))
                N_ph = 0
                do ph = 1, N_Phonon
                    if (Phonon(ph)%J == J .and. Phonon(ph)%Par == P) then
                        N_ph = N_ph + 1
                        Orb_Phonon(J,P)%v(N_ph) = ph
                    end if
                end do
            end do
        end do
        write(*,'(/,a)') 'Total number of phonons ... J_scheme phonon basis size = ' // istr(sum(N_nu))
    end subroutine HF_RPA_phonon_count

    ! TDA & RPA matrices A & B for each (J,P) channel ...
    subroutine HF_RPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,V_NN,A,B)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: V_NN
        type(rmat_t), allocatable, intent(out) :: A(:,:), B(:,:)
        integer :: J_max, J, P, N_ph, ind_ph, ind_qg, ph, qg, t_ph, t_qg, J_r, P_pg, P_pq
        type(spstate_t) :: sp, sh, sq, sg
        real(dp) :: ASum, BSum, ME

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Allocating TDA & RPA matrices A & B ...'
        allocate(A(0:J_max,2), B(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                allocate(A(J,P)%v(N_ph,N_ph), B(J,P)%v(N_ph,N_ph))
                !$omp parallel do schedule(dynamic) private(ph,t_ph,sp,sh,ind_qg,qg,t_qg,sq,sg,P_pg,P_pq,ASum,BSum,J_r,ME)
                do ind_ph = 1, N_ph
                    ph = Orb_Phonon(J,P)%v(ind_ph)
                    t_ph = Phonon(ph)%tz
                    call ph_states(ph,sp,sh)
                    do ind_qg = 1, N_ph
                        qg = Orb_Phonon(J,P)%v(ind_qg)
                        t_qg = Phonon(qg)%tz
                        call ph_states(qg,sq,sg)
                        P_pg = mod(sp%l + sg%l,2) + 1
                        P_pq = mod(sp%l + sq%l,2) + 1
                        ASum = 0.0_dp
                        BSum = 0.0_dp

                        if (t_ph == t_qg .and. ph == qg) ASum = ASum + (sp%E - sh%E)

                        if (mod(sp%l + sg%l, 2) == mod(sq%l + sh%l, 2)) then
                            do J_r = abs(sp%j - sg%j)/2, (sp%j + sg%j)/2
                                if ((abs(sq%j - sh%j) <= 2*J_r) .and. (2*J_r <= (sq%j + sh%j))) then
                                    if (t_ph == -1 .and. t_qg == -1) then
                                        ME = real(iphase((sh%j + sq%j)/2 + J_r) * (2*J_r + 1),dp) * &
                                             f6j(sp%j,sh%j,2*J,sq%j,sg%j,2*J_r) * O2b_pp(sp%a,sg%a,sh%a,sq%a,J_r,P_pg,V_NN,Orb,Orb_NN)
                                    else if (t_ph == -1 .and. t_qg == 1) then
                                        ME = real(iphase((sh%j + sq%j)/2 + J_r) * (2*J_r + 1),dp) * &
                                             f6j(sp%j,sh%j,2*J,sq%j,sg%j,2*J_r) * O2b_pn(sp%a,sg%a,sh%a,sq%a,J_r,P_pg,V_NN,Orb_NN)
                                    else if (t_ph == 1 .and. t_qg == -1) then
                                        ME = real(iphase((sg%j + sp%j)/2 + J_r) * (2*J_r + 1),dp) * &
                                             f6j(sq%j,sg%j,2*J,sp%j,sh%j,2*J_r) * O2b_pn(sq%a,sh%a,sg%a,sp%a,J_r,P_pg,V_NN,Orb_NN)
                                    else
                                        ME = real(iphase((sh%j + sq%j)/2 + J_r) * (2*J_r + 1),dp) * &
                                             f6j(sp%j,sh%j,2*J,sq%j,sg%j,2*J_r) * O2b_nn(sp%a,sg%a,sh%a,sq%a,J_r,P_pg,V_NN,Orb,Orb_NN)
                                    end if
                                    ASum = ASum + ME
                                end if
                            end do
                        end if

                        if (mod(sp%l + sq%l, 2) == mod(sh%l + sg%l, 2)) then
                            do J_r = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                                if ((abs(sh%j - sg%j) <= 2*J_r) .and. (2*J_r <= (sh%j + sg%j))) then
                                    if (t_ph == -1 .and. t_qg == -1) then
                                        ME = phase((sh%j + sq%j)/2 + J + J_r) * real(2*J_r + 1,dp) * &
                                             f6j(sp%j,sh%j,2*J,sg%j,sq%j,2*J_r) * O2b_pp(sp%a,sq%a,sh%a,sg%a,J_r,P_pq,V_NN,Orb,Orb_NN)
                                    else if (t_ph == -1 .and. t_qg == 1) then
                                        ME = phase((sh%j + sq%j)/2 + J + J_r) * real(2*J_r + 1,dp) * &
                                             f6j(sp%j,sh%j,2*J,sg%j,sq%j,2*J_r) * O2b_pn(sp%a,sq%a,sh%a,sg%a,J_r,P_pq,V_NN,Orb_NN)
                                    else if (t_ph == 1 .and. t_qg == -1) then
                                        ME = phase((sg%j + sp%j)/2 + J + J_r) * real(2*J_r + 1,dp) * &
                                             f6j(sq%j,sg%j,2*J,sh%j,sp%j,2*J_r) * O2b_pn(sq%a,sp%a,sg%a,sh%a,J_r,P_pq,V_NN,Orb_NN)
                                    else
                                        ME = phase((sh%j + sq%j)/2 + J + J_r) * real(2*J_r + 1,dp) * &
                                             f6j(sp%j,sh%j,2*J,sg%j,sq%j,2*J_r) * O2b_nn(sp%a,sq%a,sh%a,sg%a,J_r,P_pq,V_NN,Orb,Orb_NN)
                                    end if
                                    BSum = BSum + ME
                                end if
                            end do
                        end if

                        A(J,P)%v(ind_ph,ind_qg) = ASum
                        B(J,P)%v(ind_ph,ind_qg) = BSum
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        write(*,'(/,a)') 'Matrices A & B allocated ...'

    contains

        subroutine ph_states(ph,s_p,s_h)
            integer, intent(in) :: ph
            type(spstate_t), intent(out) :: s_p, s_h
            if (Phonon(ph)%tz == -1) then
                s_p = Particle%p(Phonon(ph)%p)
                s_h = Hole%p(Phonon(ph)%h)
            else
                s_p = Particle%n(Phonon(ph)%p)
                s_h = Hole%n(Phonon(ph)%h)
            end if
        end subroutine ph_states

    end subroutine HF_RPA_allocate

    ! Normalized spurious (CM) state in the 1- channel from the E1 operator ...
    function Spur_Ini(N_ph,Orb_Phonon,Phonon,Particle,Hole,TrOp) result(ph_CMS)
        integer, intent(in) :: N_ph
        integer, intent(in) :: Orb_Phonon(:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(tr1b_t), intent(in) :: TrOp
        real(dp) :: ph_CMS(N_ph)
        integer :: ph, Ind_ph, a_p, a_h
        do ph = 1, N_ph
            Ind_ph = Orb_Phonon(ph)
            call ph_orbitals(Ind_ph,Phonon,Particle,Hole,a_p,a_h)
            if (Phonon(Ind_ph)%tz == 1) then
                ph_CMS(ph) = TrOp%E1%n(a_p,a_h)
            else
                ph_CMS(ph) = TrOp%E1%p(a_p,a_h)
            end if
        end do
        where (abs(ph_CMS) < 1.0e-6_dp) ph_CMS = 0.0_dp
        ph_CMS = ph_CMS / norm2(ph_CMS)
    end function Spur_Ini

    ! Orthonormal basis with the spurious state as the column Spur_Ind ... projection & Gram-Schmidt ...
    subroutine Spur_Ortho(N_ph,Spur_State,U,Spur_Ind)
        integer, intent(in) :: N_ph
        real(dp), intent(in) :: Spur_State(:)
        real(dp), allocatable, intent(out) :: U(:,:)
        integer, intent(out) :: Spur_Ind
        real(dp), allocatable :: Id(:,:), u_v(:)
        real(dp) :: proj
        integer :: a, b, nu, mu

        Spur_Ind = 0
        Id = identity(N_ph)
        do a = 1, N_ph
            proj = 0.0_dp
            do b = 1, N_ph
                proj = proj + Spur_State(b) * Id(b,a)
            end do
            do b = 1, N_ph
                Id(b,a) = Id(b,a) - proj * Spur_State(b)
            end do
        end do
        do a = 1, N_ph
            Id(:,a) = Id(:,a) / norm2(Id(:,a))
        end do

        allocate(U(N_ph,N_ph))
        U = 0.0_dp
        U(:,1) = Id(:,1)
        do nu = 2, N_ph
            u_v = Id(:,nu)
            do mu = 1, nu - 1
                u_v = u_v - dot_product(U(:,mu),u_v) * U(:,mu)
            end do
            if (abs(norm2(u_v)) < 1.0e-4_dp) Spur_Ind = nu
            U(:,nu) = u_v / norm2(u_v)
        end do
        if (Spur_Ind == 0) call fatal('Spurious state orthogonalization failed ... no linearly dependent direction found ...')
        U(:,Spur_Ind) = Spur_State
    end subroutine Spur_Ortho

    ! Generic RPA eigenproblem for given A & B ... the three branches of HF_RPA_diagonalize ...
        ! (1) A - B positive definite, (2) A + B positive definite, (3) general complex case ...
    subroutine RPA_eigen_solve(A_JP,B_JP,E,X,Y)
        real(dp), intent(in) :: A_JP(:,:), B_JP(:,:)
        complex(dp), allocatable, intent(out) :: E(:), X(:,:), Y(:,:)
        real(dp), allocatable :: N_p(:,:), N_m(:,:), D_r(:), T_r(:,:), sqD_r(:), M_r(:,:)
        complex(dp), allocatable :: D(:), T(:,:), sqD(:), isqD(:), M(:,:), R(:,:), W(:), V(:,:)
        logical :: Is_Real
        integer :: N, ph

        N = size(A_JP,1)
        allocate(E(N), X(N,N), Y(N,N))
        if (N == 0) return
        N_p = A_JP + B_JP
        N_m = A_JP - B_JP
        allocate(W(N), R(N,N))

        if (isposdef(N_m) .or. isposdef(N_p)) then
            ! Branches (1) & (2) ... real symmetric decomposition of the positive-definite combination ...
            allocate(D_r(N), T_r(N,N))
            if (isposdef(N_m)) then
                call eig_sym(N_m,D_r,T_r)
                sqD_r = sqrt(D_r)
                M_r = mm(mm(mtm(scale_rows(T_r,sqD_r),N_p),T_r),diagm(sqD_r))
            else
                call eig_sym(N_p,D_r,T_r)
                sqD_r = sqrt(D_r)
                M_r = mm(mm(mtm(scale_rows(T_r,sqD_r),N_m),T_r),diagm(sqD_r))
            end if
            call eig_general(M_r,W,R,Is_Real)
            T = cmplx(T_r, 0.0_dp, kind=dp)
            sqD = cmplx(sqD_r, 0.0_dp, kind=dp)
            isqD = cmplx(1.0_dp / sqD_r, 0.0_dp, kind=dp)
        else
            ! Branch (3) ... general (possibly complex) decomposition of A - B ...
            allocate(D(N), T(N,N))
            call eig_general(N_m,D,T,Is_Real)
            sqD = sqrt(D)
            isqD = 1.0_dp / sqrt(D)
            M = zmm(zmm(zmhm(zscale_rows(T,sqD),cmplx(N_p,0.0_dp,kind=dp)),T),zdiagm(sqD))
            call eig_general_complex(M,W,R)
        end if

        ! Eigenvalues already sorted by (real, imag) ... E = sqrt(eigenvalue) ...
        E = sqrt(W)
        do ph = 1, N
            X(:,ph) = 0.5_dp * matmul(T, (sqD / E(ph) + isqD) * R(:,ph))
            Y(:,ph) = 0.5_dp * matmul(T, (sqD / E(ph) - isqD) * R(:,ph))
        end do

    contains

        ! diag(d) * T^T as (T * diag(d))^T operand ... i.e. returns T(:,k) * d(k) so that mtm gives diag(d) T^T ...
        function scale_rows(Tm,dv) result(S)
            real(dp), intent(in) :: Tm(:,:), dv(:)
            real(dp), allocatable :: S(:,:)
            integer :: k
            S = Tm
            do k = 1, size(dv)
                S(:,k) = Tm(:,k) * dv(k)
            end do
        end function scale_rows

        function zscale_rows(Tm,dv) result(S)
            complex(dp), intent(in) :: Tm(:,:), dv(:)
            complex(dp), allocatable :: S(:,:)
            integer :: k
            S = Tm
            do k = 1, size(dv)
                S(:,k) = Tm(:,k) * conjg(dv(k))
            end do
        end function zscale_rows

        function zdiagm(dv) result(Dm)
            complex(dp), intent(in) :: dv(:)
            complex(dp), allocatable :: Dm(:,:)
            integer :: k
            allocate(Dm(size(dv),size(dv)))
            Dm = (0.0_dp, 0.0_dp)
            do k = 1, size(dv)
                Dm(k,k) = dv(k)
            end do
        end function zdiagm

    end subroutine RPA_eigen_solve

    ! RPA eigenproblem for complex A & B (spurious-state projected RRPA channel) ... Julia's complex dispatch:
        ! Hermitian positive-definite A - B or A + B (zheevr), otherwise the general complex decomposition ...
    subroutine RPA_eigen_solve_c(A_JP,B_JP,E,X,Y)
        complex(dp), intent(in) :: A_JP(:,:), B_JP(:,:)
        complex(dp), allocatable, intent(out) :: E(:), X(:,:), Y(:,:)
        complex(dp), allocatable :: N_p(:,:), N_m(:,:), D(:), T(:,:), sqD(:), isqD(:), M(:,:), R(:,:), W(:)
        integer :: N, ph, k

        N = size(A_JP,1)
        allocate(E(N), X(N,N), Y(N,N))
        if (N == 0) return
        N_p = A_JP + B_JP
        N_m = A_JP - B_JP
        allocate(D(N), T(N,N), W(N), R(N,N))
        if (isposdef_c(N_m)) then
            call eig_general_complex(N_m,D,T)
            sqD = cmplx(sqrt(real(D)), 0.0_dp, kind=dp)
            M = N_p
        else if (isposdef_c(N_p)) then
            call eig_general_complex(N_p,D,T)
            sqD = cmplx(sqrt(real(D)), 0.0_dp, kind=dp)
            M = N_m
        else
            call eig_general_complex(N_m,D,T)
            sqD = sqrt(D)
            M = N_p
        end if
        isqD = 1.0_dp / sqD
        ! diag(sqD) T^H M T diag(sqD) ...
        M = zmm(zmhm(T,M),T)
        do k = 1, N
            M(k,:) = M(k,:) * sqD(k)
            M(:,k) = M(:,k) * sqD(k)
        end do
        call eig_general_complex(M,W,R)
        E = sqrt(W)
        do ph = 1, N
            X(:,ph) = 0.5_dp * matmul(T, (sqD / E(ph) + isqD) * R(:,ph))
            Y(:,ph) = 0.5_dp * matmul(T, (sqD / E(ph) - isqD) * R(:,ph))
        end do
    end subroutine RPA_eigen_solve_c

    ! Renormalization of RPA amplitudes (|X|^2 - |Y|^2 = +-1, or unit total norm for vanishing norm) ...
    ! Zero_Negative = .true. zeroes negative-norm modes (RRPA(0) & iteration variant) instead of multiplying by -i ...
    subroutine RPA_normalize(X,Y,Zero_Negative)
        complex(dp), intent(inout) :: X(:,:), Y(:,:)
        logical, intent(in), optional :: Zero_Negative
        integer :: nu, N
        real(dp) :: X_norm, Y_norm, Tot
        complex(dp) :: RPA_norm
        N = size(X,2)
        do nu = 1, N
            X_norm = sum(abs(X(:,nu))**2)
            Y_norm = sum(abs(Y(:,nu))**2)
            if ((X_norm - Y_norm) > 1.0e-8_dp) then
                RPA_norm = cmplx(1.0_dp / sqrt(abs(X_norm - Y_norm)), 0.0_dp, kind=dp)
            else if ((X_norm - Y_norm) < -1.0e-8_dp) then
                if (present(Zero_Negative)) then
                    if (Zero_Negative) then
                        X(:,nu) = (0.0_dp, 0.0_dp)
                        Y(:,nu) = (0.0_dp, 0.0_dp)
                        cycle
                    end if
                end if
                RPA_norm = (0.0_dp, -1.0_dp) * cmplx(1.0_dp / sqrt(abs(X_norm - Y_norm)), 0.0_dp, kind=dp)
            else
                Tot = sum(abs(X(:,nu))**2 + abs(Y(:,nu))**2)
                RPA_norm = cmplx(1.0_dp / sqrt(abs(Tot)), 0.0_dp, kind=dp)
            end if
            X(:,nu) = X(:,nu) * RPA_norm
            Y(:,nu) = Y(:,nu) * RPA_norm
        end do
    end subroutine RPA_normalize

    ! TDA & RPA diagonalization with optional orthogonalization of the spurious 1- state ...
    subroutine HF_RPA_diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,E_TDA,X_TDA,E_RPA,X_RPA,Y_RPA,Stability)
        type(params_t), intent(in) :: Params
        type(rmat_t), intent(in) :: A(0:,:), B(0:,:)
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(tr1b_t), intent(in) :: TrOp
        type(rvec_t), allocatable, intent(out) :: E_TDA(:,:)
        type(rmat_t), allocatable, intent(out) :: X_TDA(:,:)
        type(cvec_t), allocatable, intent(out) :: E_RPA(:,:)
        type(cmat_t), allocatable, intent(out) :: X_RPA(:,:), Y_RPA(:,:)
        logical, intent(out) :: Stability
        integer :: J_max, J, P, N_ph, nu
        logical :: Ortho, Spur

        J_max = 2*Params%Calc%Nmax + 1
        Ortho = Params%Calc%RPA%Ortho
        Stability = .true.
        write(*,'(/,a)') 'Diagonalizing RPA & TDA matrices ...'
        allocate(E_TDA(0:J_max,2), X_TDA(0:J_max,2), E_RPA(0:J_max,2), X_RPA(0:J_max,2), Y_RPA(0:J_max,2))

        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                Spur = Ortho .and. (J == 1) .and. (P == 2)
                call TDA_channel(A(J,P)%v,Spur,Orb_Phonon(J,P)%v,E_TDA(J,P)%v,X_TDA(J,P)%v)
                call RPA_channel(A(J,P)%v,B(J,P)%v,Spur,Orb_Phonon(J,P)%v,E_RPA(J,P)%v,X_RPA(J,P)%v,Y_RPA(J,P)%v)
                do nu = 1, N_ph
                    if (abs(aimag(E_RPA(J,P)%v(nu))) > 1.0e-8_dp) then
                        write(*,'(a)') TAB // TAB // TAB // 'RPA INSTABILITY DETECTED in channel:' // TAB // 'J = ' // istr(J) // ', P = +'
                        Stability = .false.
                    end if
                end do
            end do
        end do
        write(*,'(/,a)') 'Diagonalization done ...'

    contains

        subroutine TDA_channel(A_JP,Spur,Orb_Ph,E,X)
            real(dp), intent(in) :: A_JP(:,:)
            logical, intent(in) :: Spur
            integer, intent(in) :: Orb_Ph(:)
            real(dp), allocatable, intent(out) :: E(:), X(:,:)
            real(dp), allocatable :: Spur_State(:), V(:,:), U(:,:), A_ort(:,:)
            complex(dp), allocatable :: W(:), Z(:,:)
            integer :: N, Spur_Ind, a
            logical :: Is_Real

            N = size(A_JP,1)
            allocate(E(N), X(N,N))
            if (N == 0) return
            if (Spur) then
                Spur_State = Spur_Ini(N,Orb_Ph,Phonon,Particle,Hole,TrOp)
                call Spur_Ortho(N,Spur_State,V,Spur_Ind)
                allocate(U(N,N-1))
                do a = 1, Spur_Ind - 1
                    U(:,a) = V(:,a)
                end do
                do a = Spur_Ind + 1, N
                    U(:,a-1) = V(:,a)
                end do
                A_ort = mtm(U,mm(A_JP,U))
                allocate(W(N-1), Z(N-1,N-1))
                call eig_general(A_ort,W,Z,Is_Real)
                E(1) = 0.0_dp
                E(2:N) = real(W)
                X(:,1) = Spur_State
                X(:,2:N) = mm(U,real(Z))
            else
                allocate(W(N), Z(N,N))
                call eig_general(A_JP,W,Z,Is_Real)
                E = real(W)
                X = real(Z)
            end if
        end subroutine TDA_channel

        subroutine RPA_channel(A_JP,B_JP,Spur,Orb_Ph,E,X,Y)
            real(dp), intent(in) :: A_JP(:,:), B_JP(:,:)
            logical, intent(in) :: Spur
            integer, intent(in) :: Orb_Ph(:)
            complex(dp), allocatable, intent(out) :: E(:), X(:,:), Y(:,:)
            real(dp), allocatable :: Spur_State(:), V(:,:), U(:,:)
            complex(dp), allocatable :: E_o(:), X_o(:,:), Y_o(:,:), U_c(:,:)
            integer :: N, Spur_Ind, a

            N = size(A_JP,1)
            if (Spur .and. N > 0) then
                Spur_State = Spur_Ini(N,Orb_Ph,Phonon,Particle,Hole,TrOp)
                call Spur_Ortho(N,Spur_State,V,Spur_Ind)
                allocate(U(N,N-1))
                do a = 1, Spur_Ind - 1
                    U(:,a) = V(:,a)
                end do
                do a = Spur_Ind + 1, N
                    U(:,a-1) = V(:,a)
                end do
                call RPA_eigen_solve(mtm(U,mm(A_JP,U)),mtm(U,mm(B_JP,U)),E_o,X_o,Y_o)
                U_c = cmplx(U, 0.0_dp, kind=dp)
                allocate(E(N), X(N,N), Y(N,N))
                E(1) = (0.0_dp, 0.0_dp)
                E(2:N) = E_o
                X(:,1) = cmplx(Spur_State / sqrt(2.0_dp), 0.0_dp, kind=dp)
                Y(:,1) = cmplx(-1.0_dp * Spur_State / sqrt(2.0_dp), 0.0_dp, kind=dp)
                X(:,2:N) = zmm(U_c,X_o)
                Y(:,2:N) = zmm(U_c,Y_o)
            else
                call RPA_eigen_solve(A_JP,B_JP,E,X,Y)
            end if
            call RPA_normalize(X,Y)
        end subroutine RPA_channel

    end subroutine HF_RPA_diagonalize

    function HF_RPA_energy(Params,N_nu,E_RPA,Y_RPA) result(E_Corr_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: Y_RPA(0:,:)
        real(dp) :: E_Corr_RPA
        integer :: J_max, J, P, nu, ph, nu_min

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Calculating RPA correlation energy ...'
        E_Corr_RPA = 0.0_dp
        do J = 0, J_max
            do P = 1, 2
                nu_min = 1
                if (J == 1 .and. P == 2 .and. Params%Calc%RPA%Ortho) nu_min = 2
                do nu = nu_min, N_nu(J,P)
                    do ph = 1, N_nu(J,P)
                        E_Corr_RPA = E_Corr_RPA - real(E_RPA(J,P)%v(nu)) * abs(Y_RPA(J,P)%v(ph,nu))**2 * real(2*J + 1,dp)
                    end do
                end do
            end do
        end do
        write(*,'(/,a,/)') 'Correlation RPA energy is ...  E_RPA = ' // jstr(round_digits(E_Corr_RPA,6)) // ' MeV'
    end function HF_RPA_energy

    ! RPA one-body density matrix (quasi-boson approximation) in the HF basis ...
        ! The original accumulates the diagonal from several threads without synchronization; here serially ...
    function HF_RPA_OBDM(Params,Orb,Orb_Phonon,Phonon,Particle,Hole,N_nu,Y_RPA) result(Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        integer, intent(in) :: N_nu(0:,:)
        type(cmat_t), intent(in) :: Y_RPA(0:,:)
        type(o1b_t) :: Rho
        integer :: J_max, a_max, a, J, P, N_ph, ph, Ind_ph, a_p, a_h, j_p, j_h, nu
        real(dp) :: Sum, Y_norm

        J_max = 2*Params%Calc%Nmax + 1
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(Rho%p(a_max,a_max), Rho%n(a_max,a_max))
        Rho%p = 0.0_dp
        Rho%n = 0.0_dp
        do a = 1, a_max
            Rho%p(a,a) = real(Orb(a)%pO,dp)
            Rho%n(a,a) = real(Orb(a)%nO,dp)
        end do

        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                do ph = 1, N_ph
                    Ind_ph = Orb_Phonon(J,P)%v(ph)
                    call ph_orbitals(Ind_ph,Phonon,Particle,Hole,a_p,a_h)
                    j_p = Orb(a_p)%j
                    j_h = Orb(a_h)%j
                    Sum = 0.0_dp
                    do nu = 1, N_ph
                        Y_norm = norm2c(Y_RPA(J,P)%v(:,nu))**2
                        ! The spurious 1- state has |Y|^2 = 1/2 exactly (Y = -Spur_State/sqrt(2)); the original cut
                        ! |Y|^2 < 0.5 excludes it only up to rounding, so it is excluded explicitly here ...
                        if (abs(Y_norm) < 0.5_dp - 1.0e-8_dp) Sum = Sum + abs(Y_RPA(J,P)%v(ph,nu))**2 * real(2*J + 1,dp)
                    end do
                    if (Phonon(Ind_ph)%tz == -1) then
                        Rho%p(a_h,a_h) = Rho%p(a_h,a_h) - Sum / real(j_h + 1,dp)
                        Rho%p(a_p,a_p) = Rho%p(a_p,a_p) + Sum / real(j_p + 1,dp)
                    else
                        Rho%n(a_h,a_h) = Rho%n(a_h,a_h) - Sum / real(j_h + 1,dp)
                        Rho%n(a_p,a_p) = Rho%n(a_p,a_p) + Sum / real(j_p + 1,dp)
                    end if
                end do
            end do
        end do
    end function HF_RPA_OBDM

    pure real(dp) function norm2c(v)
        complex(dp), intent(in) :: v(:)
        norm2c = sqrt(sum(real(v)**2 + aimag(v)**2))
    end function norm2c

    ! Reduced multipole transition matrix elements rM for TDA & RPA phonons ...
    subroutine HF_RPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_TDA,X_RPA,Y_RPA,TrOp,rM_TDA,rM_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        type(cmat_t), intent(in) :: X_RPA(0:,:), Y_RPA(0:,:)
        type(tr1b_t), intent(in) :: TrOp
        type(reduced_multipole_t), intent(out) :: rM_TDA, rM_RPA

        write(*,'(/,a)') 'Calculating reduced transition matrix elements rM ...'
        call multipole(0,1,.false.,TrOp%E0,rM_TDA%E0,rM_RPA%E0)
        call multipole(1,2,.true.,TrOp%E1,rM_TDA%E1,rM_RPA%E1)
        call multipole(1,2,.true.,TrOp%E1_VC,rM_TDA%E1_VC,rM_RPA%E1_VC)
        call multipole(1,2,.true.,TrOp%E1_VS,rM_TDA%E1_VS,rM_RPA%E1_VS)
        call multipole(1,2,.true.,TrOp%E1_TC,rM_TDA%E1_TC,rM_RPA%E1_TC)
        call multipole(1,2,.true.,TrOp%E1_TS,rM_TDA%E1_TS,rM_RPA%E1_TS)
        call multipole(1,2,.true.,TrOp%E1_sTC,rM_TDA%E1_sTC,rM_RPA%E1_sTC)
        call multipole(1,2,.true.,TrOp%E1_sTS,rM_TDA%E1_sTS,rM_RPA%E1_sTS)
        call multipole(1,2,.true.,TrOp%E1_C,rM_TDA%E1_C,rM_RPA%E1_C)
        call multipole(1,2,.true.,TrOp%E1_sC,rM_TDA%E1_sC,rM_RPA%E1_sC)
        call multipole(2,1,.false.,TrOp%E2,rM_TDA%E2,rM_RPA%E2)
        call multipole(3,2,.false.,TrOp%E3,rM_TDA%E3,rM_RPA%E3)
        write(*,'(a)') TAB // 'Reduced 1-body transition matrix elements rM calculated ...'

    contains

        ! Odd = .true. for the E1 operators: TDA -O X, RPA O (-X + Y) ... otherwise TDA O X, RPA O (X + Y) ...
        subroutine multipole(J,P,Odd,O,M_TDA,M_RPA)
            integer, intent(in) :: J, P
            logical, intent(in) :: Odd
            type(o1b_t), intent(in) :: O
            type(pncvector_t), intent(out) :: M_TDA, M_RPA
            integer :: N_ph, nu, ph, i_ph, a_p, a_h
            complex(dp) :: pT, nT, pR, nR
            real(dp) :: Op

            N_ph = N_nu(J,P)
            allocate(M_TDA%p(N_ph), M_TDA%n(N_ph), M_RPA%p(N_ph), M_RPA%n(N_ph))
            !$omp parallel do private(pT,nT,pR,nR,ph,i_ph,a_p,a_h,Op)
            do nu = 1, N_ph
                pT = (0.0_dp, 0.0_dp)
                nT = (0.0_dp, 0.0_dp)
                pR = (0.0_dp, 0.0_dp)
                nR = (0.0_dp, 0.0_dp)
                do ph = 1, N_ph
                    i_ph = Orb_Phonon(J,P)%v(ph)
                    call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                    if (Phonon(i_ph)%tz == -1) then
                        Op = O%p(a_p,a_h)
                    else
                        Op = O%n(a_p,a_h)
                    end if
                    if (Phonon(i_ph)%tz == -1) then
                        if (Odd) then
                            pT = pT + cmplx(-1.0_dp * Op * X_TDA(J,P)%v(ph,nu), 0.0_dp, kind=dp)
                            pR = pR + cmplx(Op, 0.0_dp, kind=dp) * (-1.0_dp * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                        else
                            pT = pT + cmplx(Op * X_TDA(J,P)%v(ph,nu), 0.0_dp, kind=dp)
                            pR = pR + cmplx(Op, 0.0_dp, kind=dp) * (X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                        end if
                    else
                        if (Odd) then
                            nT = nT + cmplx(-1.0_dp * Op * X_TDA(J,P)%v(ph,nu), 0.0_dp, kind=dp)
                            nR = nR + cmplx(Op, 0.0_dp, kind=dp) * (-1.0_dp * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                        else
                            nT = nT + cmplx(Op * X_TDA(J,P)%v(ph,nu), 0.0_dp, kind=dp)
                            nR = nR + cmplx(Op, 0.0_dp, kind=dp) * (X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                        end if
                    end if
                end do
                M_TDA%p(nu) = pT
                M_TDA%n(nu) = nT
                M_RPA%p(nu) = pR
                M_RPA%n(nu) = nR
            end do
            !$omp end parallel do
        end subroutine multipole

    end subroutine HF_RPA_rM

    ! Reduced transition intensities B (physical, isoscalar & isovector) ...
    function HF_RPA_rB(Params,N_nu,rM,E_phonon) result(rB)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(reduced_multipole_t), intent(in) :: rM
        type(rvec_t), intent(in) :: E_phonon(0:,:)
        type(reduced_transition_t) :: rB
        real(dp), parameter :: g_p = 5.586_dp, g_n = -3.826_dp
        real(dp) :: e_p, e_n, E
        complex(dp) :: VCp, VCn, VSp, VSn, TCp, TCn, TSp, TSn, sTCp, sTCn, sTSp, sTSn, Cp, Cn, sCp, sCn, E1p, E1n
        integer :: N_ph, nu

        write(*,'(/,a)') 'Calculating the reduced transition intensities rB ...'
        call simple(rM%E0, N_nu(0,1), rB%E0)
        call simple(rM%E2, N_nu(2,1), rB%E2)
        call simple(rM%E3, N_nu(3,2), rB%E3)

        N_ph = N_nu(1,2)
        call alloc3(rB%E1,N_ph)
        call alloc3(rB%E1_V,N_ph)
        call alloc3(rB%E1_VC,N_ph)
        call alloc3(rB%E1_VS,N_ph)
        call alloc3(rB%E1_T,N_ph)
        call alloc3(rB%E1_TC,N_ph)
        call alloc3(rB%E1_TS,N_ph)
        call alloc3(rB%E1_C,N_ph)
        call alloc3(rB%E1_CC,N_ph)
        call alloc3(rB%E1_CS,N_ph)
        call alloc3(rB%E1_comp,N_ph)
        call alloc3(rB%E1_NLO_LWA,N_ph)

        e_p = real(Params%Calc%A - Params%Calc%Z,dp) / real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp) / real(Params%Calc%A,dp)

        do nu = 1, N_ph
            E = E_phonon(1,2)%v(nu)
            E1p = rM%E1%p(nu)
            E1n = rM%E1%n(nu)
            VCp = rM%E1_VC%p(nu)
            VCn = rM%E1_VC%n(nu)
            VSp = rM%E1_VS%p(nu)
            VSn = rM%E1_VS%n(nu)
            TCp = rM%E1_TC%p(nu)
            TCn = rM%E1_TC%n(nu)
            TSp = rM%E1_TS%p(nu)
            TSn = rM%E1_TS%n(nu)
            sTCp = rM%E1_sTC%p(nu)
            sTCn = rM%E1_sTC%n(nu)
            sTSp = rM%E1_sTS%p(nu)
            sTSn = rM%E1_sTS%n(nu)
            Cp = rM%E1_C%p(nu)
            Cn = rM%E1_C%n(nu)
            sCp = rM%E1_sC%p(nu)
            sCn = rM%E1_sC%n(nu)

            ! Physical components ...
            rB%E1%ph(nu) = abs(E1p)**2
            rB%E1_V%ph(nu) = abs(VCp + 0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn)**2
            rB%E1_VC%ph(nu) = abs(VCp)**2
            rB%E1_VS%ph(nu) = abs(0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn)**2
            rB%E1_T%ph(nu) = abs(TCp + 0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn)**2
            rB%E1_TC%ph(nu) = abs(TCp)**2
            rB%E1_TS%ph(nu) = abs(0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn)**2
            rB%E1_C%ph(nu) = abs(TCp + 0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn - (VCp + 0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn))**2
            rB%E1_CC%ph(nu) = abs(TCp - VCp)**2
            rB%E1_CS%ph(nu) = abs(0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn - (0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn))**2
            rB%E1_comp%ph(nu) = abs(Cp)**2
            rB%E1_NLO_LWA%ph(nu) = abs(E1p + E / HbarC * (TCp + 0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn))**2

            ! Isoscalar components ...
            rB%E1%is(nu) = 0.25_dp * abs(E1p + E1n)**2
            rB%E1_V%is(nu) = abs(0.5_dp * (VCp + VCn) + 0.125_dp * (g_p + g_n) * (VSp + VSn))**2
            rB%E1_VC%is(nu) = abs(0.5_dp * (VCp + VCn))**2
            rB%E1_VS%is(nu) = abs(0.125_dp * (g_p + g_n) * (VSp + VSn))**2
            rB%E1_T%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn) + 0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn))**2
            rB%E1_TC%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn))**2
            rB%E1_TS%is(nu) = abs(0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn))**2
            rB%E1_C%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn) + 0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn) - &
                             (0.5_dp * (VCp + VCn) + 0.125_dp * (g_p + g_n) * (VSp + VSn)))**2
            rB%E1_CC%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn) - (0.5_dp * (VCp + VCn)))**2
            rB%E1_CS%is(nu) = abs(0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn) - (0.125_dp * (g_p + g_n) * (VSp + VSn)))**2
            rB%E1_comp%is(nu) = 0.25_dp * abs(Cp + sCp + Cn + sCn)**2
            rB%E1_NLO_LWA%is(nu) = 0.25_dp * abs(E1p + E1n + E / HbarC * ((TCp + sTCp + TCn + sTCn) + 0.25_dp * (g_p + g_n) * &
                                   (TSp + sTSp + TSn + sTSn)))**2

            ! Isovector components ...
            rB%E1%iv(nu) = abs(e_p * E1p - e_n * E1n)**2
            rB%E1_V%iv(nu) = abs(0.5_dp * (VCp - VCn) + 0.125_dp * (g_p - g_n) * (VSp + VSn))**2
            rB%E1_VC%iv(nu) = abs(0.5_dp * (VCp - VCn))**2
            rB%E1_VS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (VSp + VSn))**2
            rB%E1_T%iv(nu) = abs(0.5_dp * (TCp - TCn) + 0.125_dp * (g_p - g_n) * (TSp + TSn))**2
            rB%E1_TC%iv(nu) = abs(0.5_dp * (TCp - TCn))**2
            rB%E1_TS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (TSp + TSn))**2
            rB%E1_C%iv(nu) = abs(0.5_dp * (TCp - TCn) + 0.125_dp * (g_p - g_n) * (TSp + TSn) - &
                             (0.5_dp * (VCp - VCn) + 0.125_dp * (g_p - g_n) * (VSp + VSn)))**2
            rB%E1_CC%iv(nu) = abs(0.5_dp * (TCp - TCn) - (0.5_dp * (VCp - VCn)))**2
            rB%E1_CS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (TSp + TSn) - (0.125_dp * (g_p - g_n) * (VSp + VSn)))**2
            rB%E1_comp%iv(nu) = 0.25_dp * abs(Cp - Cn)**2
            rB%E1_NLO_LWA%iv(nu) = abs(e_p * E1p - e_n * E1n + 0.5_dp * E / HbarC * ((TCp - TCn) + 0.125_dp * (g_p - g_n) * &
                                   (TSp + TSn)))**2
        end do

        write(*,'(a)') TAB // 'Reduced transition intensities rB evaluated ...'

    contains

        subroutine alloc3(T,N)
            type(transition_t), intent(out) :: T
            integer, intent(in) :: N
            allocate(T%ph(N), T%is(N), T%iv(N))
        end subroutine alloc3

        subroutine simple(M,N,T)
            type(pncvector_t), intent(in) :: M
            integer, intent(in) :: N
            type(transition_t), intent(out) :: T
            integer :: k
            allocate(T%ph(N), T%is(N), T%iv(N))
            do k = 1, N
                T%ph(k) = abs(M%p(k))**2
                T%is(k) = 0.25_dp * abs(M%p(k) + M%n(k))**2
                T%iv(k) = 0.25_dp * abs(M%p(k) - M%n(k))**2
            end do
        end subroutine simple

    end function HF_RPA_rB

end module nj_hf_rpa
