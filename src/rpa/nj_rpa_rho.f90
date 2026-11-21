! Renormalized-RPA one-body density matrix of HF-ScRPA (NuJuliet dev, stage 4) ... counterpart of src/RPA/RPA_rho.jl
    ! (RPA_OBDM_oY4, HF_ScRPA_OBDM_solve, RPA_OBDM_transform_amplitudes) & src/RPA/RPA_orbitals.jl (Orb1B_RPA_OBDM,
    ! RPA_OBDM_orbitals_make) ... RPA_OBDM_oY2 is the plain-RPA routine of nj_rpa ... replaces the
    ! HF_RRPA_OBDM_iteration of the base-version HF-RRPA (removed) ...
    ! Amplitude layout as RPA_diagonalize: X(P,J)%v(N_ph,N_ph) over the fixed Orb_1p1h list, metric normalization
    ! X^+ N X - Y^+ N Y = 1 with N = diag(n_h - n_p) ...
    ! Deviations (no change of results beyond rounding, docs/PORTING_NOTES.md "stage 4"):
    !   - RPA_OBDM_oY4: the M matrices are block-local (dev's caller-owned pM/nM storage is never read back); M & G of one
    !     isospin are built only from the 1p-1h rows of that isospin (M from the rows with W > 0, G = Y_T (1 + M) Y_T^+ on the
    !     rows of isospin T), only Re G is stored ... blocks whose X & Y are real (all imaginary parts exactly 0, i.e. no
    !     complex RPA solution in the block) use real BLAS, all other blocks complex BLAS as dev ...
    !   - RPA_OBDM_transform_amplitudes: U is block-diagonal in the isospin, X_new = U X is done per isospin block, with real
    !     BLAS on the real & (if non-zero) imaginary parts ...
    !   - BLAS inside the OpenMP block loops runs single-threaded per thread (mkl_set_num_threads_local, nj_blas_threads) ...
module nj_rpa_rho
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_blas_threads, only : blas_threads_local
    use nj_orbitals, only : orb2b_1p1h_index, OBDM_natural_orbitals, a_max_of
    use nj_o1b, only : O1b_import, O1b_transformation, O1b_diagonal, O1b_unit_matrix
    !$ use omp_lib
    implicit none
    private

    public :: RPA_OBDM_orbitals_make, RPA_OBDM_oY4, RPA_OBDM_transform_amplitudes, HF_ScRPA_OBDM_solve

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character(len=1), intent(in) :: Transa, Transb
            integer, intent(in) :: m, n, k, Lda, Ldb, Ldc
            real(dp), intent(in) :: Alpha, Beta, A(Lda,*), B(Ldb,*)
            real(dp), intent(inout) :: C(Ldc,*)
        end subroutine dgemm
        subroutine zgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character(len=1), intent(in) :: Transa, Transb
            integer, intent(in) :: m, n, k, Lda, Ldb, Ldc
            complex(dp), intent(in) :: Alpha, Beta, A(Lda,*), B(Ldb,*)
            complex(dp), intent(inout) :: C(Ldc,*)
        end subroutine zgemm
    end interface

contains

    ! Particle (|O| < Tol_par) & hole (|O| > Tol_hol) orbitals of the RPA OBDM, per isospin in Orb order ... counterpart of
        ! RPA_OBDM_orbitals_make (RPA_orbitals.jl:8-68) ... strict tests (unlike the inclusive ones of RPA_sigma; dev, B3) ...
    function RPA_OBDM_orbitals_make(Params,Orb) result(Orb_ph)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb1b_rpa_obdm_t) :: Orb_ph
        real(dp) :: Tol_par, Tol_hol
        integer :: ia

        Tol_par = Params%Calc%RPA%Tol_par
        Tol_hol = Params%Calc%RPA%Tol_hol
        Orb_ph%N_p%p = count(abs(Orb(:)%pO) < Tol_par)
        Orb_ph%N_h%p = count(abs(Orb(:)%pO) > Tol_hol)
        Orb_ph%N_p%n = count(abs(Orb(:)%nO) < Tol_par)
        Orb_ph%N_h%n = count(abs(Orb(:)%nO) > Tol_hol)
        allocate(Orb_ph%p%p(Orb_ph%N_p%p), Orb_ph%h%p(Orb_ph%N_h%p), Orb_ph%p%n(Orb_ph%N_p%n), Orb_ph%h%n(Orb_ph%N_h%n))
        Orb_ph%N_p = pninteger_t(0,0)
        Orb_ph%N_h = pninteger_t(0,0)
        do ia = 1, size(Orb)
            if (abs(Orb(ia)%pO) < Tol_par) then
                Orb_ph%N_p%p = Orb_ph%N_p%p + 1
                Orb_ph%p%p(Orb_ph%N_p%p) = Orb(ia)
            end if
            if (abs(Orb(ia)%pO) > Tol_hol) then
                Orb_ph%N_h%p = Orb_ph%N_h%p + 1
                Orb_ph%h%p(Orb_ph%N_h%p) = Orb(ia)
            end if
            if (abs(Orb(ia)%nO) < Tol_par) then
                Orb_ph%N_p%n = Orb_ph%N_p%n + 1
                Orb_ph%p%n(Orb_ph%N_p%n) = Orb(ia)
            end if
            if (abs(Orb(ia)%nO) > Tol_hol) then
                Orb_ph%N_h%n = Orb_ph%N_h%n + 1
                Orb_ph%h%n(Orb_ph%N_h%n) = Orb(ia)
            end if
        end do
    end function RPA_OBDM_orbitals_make

    ! (J,P) blocks of Orb_1p1h sorted by size, largest first (stable, as sortperm(...,rev = true) of the JP list J outer,
        ! P inner) ... only the load balance depends on this order ...
    subroutine jp_blocks_sorted(Orb_1p1h,JP)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        integer, allocatable, intent(out) :: JP(:,:)
        real(dp), allocatable :: Sz(:)
        integer, allocatable :: Perm(:)
        integer :: N_JP, i, J, P

        N_JP = 2*(Orb_1p1h%J_max + 1)
        allocate(Sz(N_JP), Perm(N_JP), JP(2,N_JP))
        i = 0
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                i = i + 1
                Sz(i) = real(Orb_1p1h%N(P,J),dp)
            end do
        end do
        call sortperm_rev(Sz,Perm)
        do i = 1, N_JP
            JP(1,i) = (Perm(i) - 1)/2
            JP(2,i) = mod(Perm(i) - 1,2) + 1
        end do
    end subroutine jp_blocks_sorted

    ! RPA OBDM at the O(Y^4) level ... counterpart of RPA_OBDM_oY4 (RPA_rho.jl:283-478) ...
        !   M_nu,mu = -1/2 sum_qg W_qg X_qg,nu X*_qg,mu   (W_qg = D_qg^2 for the particle-hole pairs of Orb_ph of isospin T with
        !                                                   D_qg = rho_g - rho_q >= 1e-8 from the diagonal of Rho, else 0)
        !   G = Y (1 + M) Y^+ ,  dRho_ab = +- 1/(2j_a+1) sum_(h|p,J) (2J+1) D D Re G   (particles a,b: +, holes a,b: -)
        ! Returns Rho + dRho (Rho is copied, dRho added at the (l,j)-equal particle & hole pairs of Orb_ph) ...
        ! dev's pM/nM work arrays are block-local here (never read by the callers) ...
    function RPA_OBDM_oY4(Params,Orb,Orb_ph,Orb_1p1h,Rho,X,Y) result(Rho_new)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb1b_rpa_obdm_t), intent(in) :: Orb_ph
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(o1b_t), intent(in) :: Rho
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(o1b_t) :: Rho_new
        logical, allocatable :: p_par(:), p_hol(:), n_par(:), n_hol(:)
        type(rmat_t), allocatable :: G(:,:,:)
        type(ivec_t), allocatable :: Loc(:,:)
        integer, allocatable :: JP(:,:), Tasks(:,:)
        real(dp), allocatable :: dRho(:)
        integer :: a_max, J_max, N_JP, i_JP, J, P, N_task, i_task, T, a, b, ia, ib, k, N_old, it
        integer :: i_par, i_hol, Cnt_par, Cnt_hol, Pass
        real(dp) :: Sum

        a_max = a_max_of(Params%Calc%Nmax)
        J_max = Orb_1p1h%J_max
        allocate(p_par(a_max), p_hol(a_max), n_par(a_max), n_hol(a_max))
        p_par = .false.
        p_hol = .false.
        n_par = .false.
        n_hol = .false.
        do k = 1, Orb_ph%N_p%p
            p_par(Orb_ph%p%p(k)%a) = .true.
        end do
        do k = 1, Orb_ph%N_h%p
            p_hol(Orb_ph%h%p(k)%a) = .true.
        end do
        do k = 1, Orb_ph%N_p%n
            n_par(Orb_ph%p%n(k)%a) = .true.
        end do
        do k = 1, Orb_ph%N_h%n
            n_hol(Orb_ph%h%n(k)%a) = .true.
        end do

        ! 1) Re G per (J,P) block & isospin (G(P,J,1) protons, G(P,J,2) neutrons) ... one block per thread, 1-thread BLAS ...
        allocate(G(2,0:J_max,2), Loc(2,0:J_max))
        call jp_blocks_sorted(Orb_1p1h,JP)
        N_JP = size(JP,2)
        !$omp parallel default(shared) private(i_JP,J,P,N_old)
        N_old = blas_threads_local(1)
        !$omp do schedule(dynamic,1)
        do i_JP = 1, N_JP
            J = JP(1,i_JP)
            P = JP(2,i_JP)
            call oY4_block(Orb_1p1h%N(P,J),Orb_1p1h%Ind(P,J)%v,X(P,J)%v,Y(P,J)%v,Rho,p_par,p_hol,n_par,n_hol, &
                           G(P,J,1)%v,G(P,J,2)%v,Loc(P,J)%v)
        end do
        !$omp end do
        N_old = blas_threads_local(N_old)
        !$omp end parallel

        ! 2) OBDM element list (T, Type, a, b) ... T = -1 then +1, particle pairs (Type 1) then hole pairs (Type 2), a outer ...
        do Pass = 1, 2
            N_task = 0
            do it = 1, 2
                T = 2*it - 3
                Cnt_par = merge(Orb_ph%N_p%p, Orb_ph%N_p%n, T == -1)
                Cnt_hol = merge(Orb_ph%N_h%p, Orb_ph%N_h%n, T == -1)
                do ia = 1, Cnt_par
                    do ib = 1, Cnt_par
                        a = ph_orb(Orb_ph,1,T,ia)
                        b = ph_orb(Orb_ph,1,T,ib)
                        if (Orb(a)%l /= Orb(b)%l .or. Orb(a)%j /= Orb(b)%j) cycle
                        N_task = N_task + 1
                        if (Pass == 2) Tasks(:,N_task) = [T, 1, a, b]
                    end do
                end do
                do ia = 1, Cnt_hol
                    do ib = 1, Cnt_hol
                        a = ph_orb(Orb_ph,2,T,ia)
                        b = ph_orb(Orb_ph,2,T,ib)
                        if (Orb(a)%l /= Orb(b)%l .or. Orb(a)%j /= Orb(b)%j) cycle
                        N_task = N_task + 1
                        if (Pass == 2) Tasks(:,N_task) = [T, 2, a, b]
                    end do
                end do
            end do
            if (Pass == 1) allocate(Tasks(4,N_task))
        end do
        allocate(dRho(N_task))

        ! 3) OBDM elements ... each task writes only dRho(i_task), sums in dev's order (h | p in Orb_ph order, J ascending) ...
        !$omp parallel do schedule(dynamic,16) default(shared) private(i_task,T,a,b,Sum,i_par,i_hol,Cnt_par,Cnt_hol)
        do i_task = 1, N_task
            T = Tasks(1,i_task)
            a = Tasks(3,i_task)
            b = Tasks(4,i_task)
            Sum = 0.0_dp
            if (Tasks(2,i_task) == 1) then
                Cnt_hol = merge(Orb_ph%N_h%p, Orb_ph%N_h%n, T == -1)
                do i_hol = 1, Cnt_hol
                    if (T == -1) then
                        Sum = Sum + elem_sum(a,b,Orb_ph%h%p(i_hol)%a,.true.,Rho%p,1)
                    else
                        Sum = Sum + elem_sum(a,b,Orb_ph%h%n(i_hol)%a,.true.,Rho%n,2)
                    end if
                end do
            else
                Cnt_par = merge(Orb_ph%N_p%p, Orb_ph%N_p%n, T == -1)
                do i_par = 1, Cnt_par
                    if (T == -1) then
                        Sum = Sum - elem_sum(a,b,Orb_ph%p%p(i_par)%a,.false.,Rho%p,1)
                    else
                        Sum = Sum - elem_sum(a,b,Orb_ph%p%n(i_par)%a,.false.,Rho%n,2)
                    end if
                end do
            end if
            dRho(i_task) = Sum / hat_of(Orb(a)%j)
        end do
        !$omp end parallel do

        ! 4) Rho_new = Rho + dRho ...
        Rho_new = Rho
        do i_task = 1, N_task
            a = Tasks(3,i_task)
            b = Tasks(4,i_task)
            if (Tasks(1,i_task) == -1) then
                Rho_new%p(a,b) = Rho_new%p(a,b) + dRho(i_task)
            else
                Rho_new%n(a,b) = Rho_new%n(a,b) + dRho(i_task)
            end if
        end do

    contains

        ! a-th particle (Kind 1) or hole (Kind 2) orbital index of isospin T ...
        pure integer function ph_orb(O,Kind,T_,i) result(a_)
            type(orb1b_rpa_obdm_t), intent(in) :: O
            integer, intent(in) :: Kind, T_, i
            if (Kind == 1) then
                if (T_ == -1) then
                    a_ = O%p%p(i)%a
                else
                    a_ = O%p%n(i)%a
                end if
            else
                if (T_ == -1) then
                    a_ = O%h%p(i)%a
                else
                    a_ = O%h%n(i)%a
                end if
            end if
        end function ph_orb

        ! Contribution of one summation orbital c (a hole for Par = .true., a particle otherwise) to the element (a,b):
            ! sum_J (2J+1) D_1 D_2 Re G[i_1,i_2] (0 if one of the D < 1e-8) ... dev's inner J loop ...
        real(dp) function elem_sum(a_,b_,c_,Par,Rho_T,iT) result(S)
            integer, intent(in) :: a_, b_, c_, iT
            logical, intent(in) :: Par
            real(dp), intent(in) :: Rho_T(:,:)
            real(dp) :: D_1, D_2
            integer :: P_, J_, i_1, i_2, T_
            S = 0.0_dp
            T_ = 2*iT - 3
            if (Par) then
                D_1 = Rho_T(c_,c_) - Rho_T(a_,a_)
                D_2 = Rho_T(c_,c_) - Rho_T(b_,b_)
            else
                D_1 = Rho_T(a_,a_) - Rho_T(c_,c_)
                D_2 = Rho_T(b_,b_) - Rho_T(c_,c_)
            end if
            if (D_1 < 1.0e-8_dp .or. D_2 < 1.0e-8_dp) return
            P_ = parity_of(Orb(a_)%l + Orb(c_)%l)
            do J_ = abs(Orb(a_)%j - Orb(c_)%j)/2, (Orb(a_)%j + Orb(c_)%j)/2
                if (Par) then
                    i_1 = orb2b_1p1h_index(a_,c_,J_,P_,T_,Orb_1p1h)
                    i_2 = orb2b_1p1h_index(b_,c_,J_,P_,T_,Orb_1p1h)
                else
                    i_1 = orb2b_1p1h_index(c_,a_,J_,P_,T_,Orb_1p1h)
                    i_2 = orb2b_1p1h_index(c_,b_,J_,P_,T_,Orb_1p1h)
                end if
                if (i_1 == 0 .or. i_2 == 0) cycle
                S = S + hat_of(2*J_) * D_1 * D_2 * G(P_,J_,iT)%v(Loc(P_,J_)%v(i_1),Loc(P_,J_)%v(i_2))
            end do
        end function elem_sum

    end function RPA_OBDM_oY4

    ! One (J,P) block of RPA_OBDM_oY4: Loc(i) = index of pair i among the pairs of its isospin; Gp/Gn = Re G restricted to the
        ! proton/neutron pairs ... M from the rows with W > 0 only (W = 0 rows add exact zeros in dev) ...
    subroutine oY4_block(N_ph,Ind,X,Y,Rho,p_par,p_hol,n_par,n_hol,Gp,Gn,Loc)
        integer, intent(in) :: N_ph
        type(ind2b_t), intent(in) :: Ind(:)
        complex(dp), intent(in) :: X(:,:), Y(:,:)
        type(o1b_t), intent(in) :: Rho
        logical, intent(in) :: p_par(:), p_hol(:), n_par(:), n_hol(:)
        real(dp), allocatable, intent(out) :: Gp(:,:), Gn(:,:)
        integer, allocatable, intent(out) :: Loc(:)
        integer, allocatable :: Rows(:), Rows_W(:)
        real(dp), allocatable :: W(:)
        logical :: Is_Real
        integer :: i, iT, T, N_T, N_W
        real(dp) :: D

        allocate(Loc(N_ph), W(N_ph))
        if (N_ph == 0) then
            allocate(Gp(0,0), Gn(0,0))
            return
        end if
        Is_Real = all(aimag(X) == 0.0_dp) .and. all(aimag(Y) == 0.0_dp)

        ! Weights D_qg^2 of the particle-hole pairs of Orb_ph with D_qg >= 1e-8 ...
        W = 0.0_dp
        do i = 1, N_ph
            if (Ind(i)%T == -1 .and. p_par(Ind(i)%a) .and. p_hol(Ind(i)%b)) then
                D = Rho%p(Ind(i)%b,Ind(i)%b) - Rho%p(Ind(i)%a,Ind(i)%a)
                if (D >= 1.0e-8_dp) W(i) = D**2
            else if (Ind(i)%T == 1 .and. n_par(Ind(i)%a) .and. n_hol(Ind(i)%b)) then
                D = Rho%n(Ind(i)%b,Ind(i)%b) - Rho%n(Ind(i)%a,Ind(i)%a)
                if (D >= 1.0e-8_dp) W(i) = D**2
            end if
        end do

        do iT = 1, 2
            T = 2*iT - 3
            N_T = count(Ind(1:N_ph)%T == T)
            N_W = count(Ind(1:N_ph)%T == T .and. W > 0.0_dp)
            allocate(Rows(N_T), Rows_W(N_W))
            N_T = 0
            N_W = 0
            do i = 1, N_ph
                if (Ind(i)%T /= T) cycle
                N_T = N_T + 1
                Rows(N_T) = i
                Loc(i) = N_T
                if (W(i) > 0.0_dp) then
                    N_W = N_W + 1
                    Rows_W(N_W) = i
                end if
            end do
            if (iT == 1) then
                call g_block(Gp)
            else
                call g_block(Gn)
            end if
            deallocate(Rows, Rows_W)
        end do

    contains

        ! Re G = Re[ Y_T (1 + M) Y_T^+ ] with M = -1/2 X_W^T diag(W) X_W^* (X_W = rows with W > 0) ...
        subroutine g_block(Gr)
            real(dp), allocatable, intent(out) :: Gr(:,:)
            real(dp), allocatable :: rM(:,:), rXW(:,:), rWX(:,:), rY(:,:), rK(:,:)
            complex(dp), allocatable :: zM(:,:), zXW(:,:), zWX(:,:), zY(:,:), zK(:,:), zG(:,:)
            integer :: k

            allocate(Gr(N_T,N_T))
            if (N_T == 0) return
            if (Is_Real) then
                allocate(rY(N_T,N_ph), rK(N_ph,N_T))
                rY = real(Y(Rows,1:N_ph),dp)
                rK = transpose(rY)
                if (N_W > 0) then
                    allocate(rM(N_ph,N_ph), rXW(N_W,N_ph), rWX(N_W,N_ph))
                    rXW = real(X(Rows_W,1:N_ph),dp)
                    do k = 1, N_ph
                        rWX(:,k) = W(Rows_W) * rXW(:,k)
                    end do
                    call dgemm('T','N',N_ph,N_ph,N_W,-0.5_dp,rXW,N_W,rWX,N_W,0.0_dp,rM,N_ph)
                    ! K = (1 + M) Y_T^T ...
                    call dgemm('N','T',N_ph,N_T,N_ph,1.0_dp,rM,N_ph,rY,N_T,1.0_dp,rK,N_ph)
                end if
                call dgemm('N','N',N_T,N_T,N_ph,1.0_dp,rY,N_T,rK,N_ph,0.0_dp,Gr,N_T)
            else
                allocate(zY(N_T,N_ph), zK(N_ph,N_T), zG(N_T,N_T))
                zY = Y(Rows,1:N_ph)
                zK = conjg(transpose(zY))
                if (N_W > 0) then
                    allocate(zM(N_ph,N_ph), zXW(N_W,N_ph), zWX(N_W,N_ph))
                    zXW = X(Rows_W,1:N_ph)
                    do k = 1, N_ph
                        zWX(:,k) = W(Rows_W) * conjg(zXW(:,k))
                    end do
                    call zgemm('T','N',N_ph,N_ph,N_W,(-0.5_dp,0.0_dp),zXW,N_W,zWX,N_W,(0.0_dp,0.0_dp),zM,N_ph)
                    call zgemm('N','C',N_ph,N_T,N_ph,(1.0_dp,0.0_dp),zM,N_ph,zY,N_T,(1.0_dp,0.0_dp),zK,N_ph)
                end if
                call zgemm('N','N',N_T,N_T,N_ph,(1.0_dp,0.0_dp),zY,N_T,zK,N_ph,(0.0_dp,0.0_dp),zG,N_T)
                Gr = real(zG,dp)
            end if
        end subroutine g_block

    end subroutine oY4_block

    ! Rotation of the RPA amplitudes into a new single-particle basis C (per isospin) ... counterpart of
        ! RPA_OBDM_transform_amplitudes (RPA_rho.jl:717-776): X_new[ph,nu] = sum_qg U[ph,qg] X[qg,nu],
        ! U[ph,qg] = C_T[q,p] C_T[g,h] (same isospin, else 0) over the fixed Orb_1p1h list ... returns new arrays ...
    subroutine RPA_OBDM_transform_amplitudes(Orb_1p1h,X,Y,C,X_new,Y_new)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(o1b_t), intent(in) :: C
        type(cmat_t), allocatable, intent(out) :: X_new(:,:), Y_new(:,:)
        integer, allocatable :: JP(:,:)
        integer :: N_JP, i_JP, J, P, N_old

        allocate(X_new(2,0:Orb_1p1h%J_max), Y_new(2,0:Orb_1p1h%J_max))
        call jp_blocks_sorted(Orb_1p1h,JP)
        N_JP = size(JP,2)
        !$omp parallel default(shared) private(i_JP,J,P,N_old)
        N_old = blas_threads_local(1)
        !$omp do schedule(dynamic,1)
        do i_JP = 1, N_JP
            J = JP(1,i_JP)
            P = JP(2,i_JP)
            call transform_block(Orb_1p1h%N(P,J),Orb_1p1h%Ind(P,J)%v,C,X(P,J)%v,X_new(P,J)%v)
            call transform_block(Orb_1p1h%N(P,J),Orb_1p1h%Ind(P,J)%v,C,Y(P,J)%v,Y_new(P,J)%v)
        end do
        !$omp end do
        N_old = blas_threads_local(N_old)
        !$omp end parallel
    end subroutine RPA_OBDM_transform_amplitudes

    ! A_new = U A for one block, per isospin (U is block diagonal in T), real BLAS on Re A & Im A ...
    subroutine transform_block(N_ph,Ind,C,A,A_new)
        integer, intent(in) :: N_ph
        type(ind2b_t), intent(in) :: Ind(:)
        type(o1b_t), intent(in) :: C
        complex(dp), intent(in) :: A(:,:)
        complex(dp), allocatable, intent(out) :: A_new(:,:)
        integer, allocatable :: Rows(:)
        real(dp), allocatable :: U(:,:), Ar(:,:), Ai(:,:), Br(:,:), Bi(:,:)
        integer :: iT, T, N_T, i, k, a_p, a_h, a_q, a_g
        logical :: Has_Im

        A_new = A
        if (N_ph == 0) return
        do iT = 1, 2
            T = 2*iT - 3
            N_T = count(Ind(1:N_ph)%T == T)
            if (N_T == 0) cycle
            allocate(Rows(N_T), U(N_T,N_T))
            N_T = 0
            do i = 1, N_ph
                if (Ind(i)%T == T) then
                    N_T = N_T + 1
                    Rows(N_T) = i
                end if
            end do
            do k = 1, N_T
                a_q = Ind(Rows(k))%a
                a_g = Ind(Rows(k))%b
                do i = 1, N_T
                    a_p = Ind(Rows(i))%a
                    a_h = Ind(Rows(i))%b
                    if (T == -1) then
                        U(i,k) = C%p(a_q,a_p) * C%p(a_g,a_h)
                    else
                        U(i,k) = C%n(a_q,a_p) * C%n(a_g,a_h)
                    end if
                end do
            end do
            Ar = real(A(Rows,1:N_ph),dp)
            allocate(Br(N_T,N_ph))
            call dgemm('N','N',N_T,N_ph,N_T,1.0_dp,U,N_T,Ar,N_T,0.0_dp,Br,N_T)
            Has_Im = any(aimag(A(Rows,1:N_ph)) /= 0.0_dp)
            if (Has_Im) then
                Ai = aimag(A(Rows,1:N_ph))
                allocate(Bi(N_T,N_ph))
                call dgemm('N','N',N_T,N_ph,N_T,1.0_dp,U,N_T,Ai,N_T,0.0_dp,Bi,N_T)
                A_new(Rows,1:N_ph) = cmplx(Br,Bi,kind=dp)
                deallocate(Ai, Bi)
            else
                A_new(Rows,1:N_ph) = cmplx(Br,0.0_dp,kind=dp)
            end if
            deallocate(Rows, U, Ar, Br)
        end do
    end subroutine transform_block

    ! Degeneracy-weighted Frobenius norm ||O||_w = sqrt( sum_a sum_b (2j_a+1) (pO_ab^2 + nO_ab^2) ) ... (a over Orb, b over a_max)
    real(dp) function obdm_norm(Orb,a_max,pO,nO) result(Nrm)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: a_max
        real(dp), intent(in) :: pO(:,:), nO(:,:)
        integer :: ia, a, b
        real(dp) :: Sum
        Sum = 0.0_dp
        do ia = 1, size(Orb)
            a = Orb(ia)%a
            do b = 1, a_max
                Sum = Sum + hat_of(Orb(ia)%j) * (pO(a,b)**2 + nO(a,b)**2)
            end do
        end do
        Nrm = sqrt(Sum)
    end function obdm_norm

    ! Self-consistent renormalized-RPA OBDM for fixed RPA amplitudes ... counterpart of HF_ScRPA_OBDM_solve
        ! (RPA_rho.jl:528-715) ... Bases LHO --C_LHO_HF--> HF --C_HF_NAT--> NAT (current) --C_tot--> NAT (new) ...
        ! Iteration (hard-coded as dev: tolerance 1e-8, 150 iterations, 2-cycle test after 100 iterations):
        !   Rho_new = Rho_ref + (oY4[Rho_old] - Rho_old) [diagonal if OBDM_diagonal], mixing alpha = OBDM_alpha, natural
        !   orbitals (OBDM_natural_orbitals, the unsplit Rho transformed & reduced to its diagonal), rotation of Rho_old,
        !   Rho_ref, X & Y, delta = ||Rho_new - Rho_old||_w / alpha / ||Rho_new - Rho_ref||_w ... divergence & 2-cycle guards,
        !   not converged: one unmixed explicit evaluation Rho_ref + dRho[Rho_ini] with the input X, Y ...
        ! Interface (Fortran form of `Rho, C_LHO_NAT, Orb = HF_ScRPA_OBDM_solve(...,X,Y,Rho,C_HF_NAT)`):
        !   Orb      in: current NAT orbitals; out: the same n,l,j with the new occupations (dev's new Orb_NAT vector),
        !   Rho      in: Rho_ini (current NAT basis); out: the new diagonal OBDM (new NAT basis),
        !   C_HF_NAT in: HF -> current NAT; out: HF -> new NAT (dev: in place),  C_LHO_NAT out: C_LHO_HF C_HF_NAT ...
        !   X, Y are not modified (rotated copies) ... rho.bin & C.bin are read from io/<Path>/<RPA.Reference>/bin on every call
        !   (as dev) ...
    subroutine HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT,C_LHO_NAT)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(inout) :: Orb(:)
        type(orb1b_rpa_obdm_t), intent(in) :: Orb_ph
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(o1b_t), intent(inout) :: Rho, C_HF_NAT
        type(o1b_t), intent(out) :: C_LHO_NAT
        real(dp), parameter :: epsilon = 1.0e-8_dp
        integer, parameter :: Iteration_max = 150
        type(cmat_t), allocatable :: X_it(:,:), Y_it(:,:), X_tmp(:,:), Y_tmp(:,:)
        type(o1b_t) :: Rho_ini, Rho_HF, C_LHO_HF, Rho_ref, Rho_old, Rho_ref_it, Rho_new, C_tot, C_step, Rho_expl
        character(len=:), allocatable :: IO_HF
        logical :: Diagonal_OBDM, Converged
        real(dp) :: alpha, Norm_ref, Norm_cor, delta, delta_abs, delta_odd, delta_even, delta_prev
        integer :: a_max, Iteration, ia, a

        a_max = a_max_of(Params%Calc%Nmax)
        Diagonal_OBDM = Params%Calc%RPA%OBDM_diagonal
        alpha = Params%Calc%RPA%OBDM_alpha
        Rho_ini = Rho

        ! Reference OBDM (HF basis) & LHO -> HF ... reference OBDM in the current NAT basis ...
        IO_HF = io_dir(Params,Params%Calc%RPA%Reference)
        Rho_HF = O1b_import(Params,Orb,IO_HF // '/bin/rho.bin')
        C_LHO_HF = O1b_import(Params,Orb,IO_HF // '/bin/C.bin')
        Rho_ref = O1b_transformation(Rho_HF,C_HF_NAT)

        X_it = X
        Y_it = Y
        Rho_old = Rho_ini
        Rho_ref_it = Rho_ref
        C_tot = O1b_unit_matrix(Params)
        Norm_ref = obdm_norm(Orb,a_max,Rho_ref%p,Rho_ref%n)

        Iteration = 0
        Converged = .false.
        delta = 1.0_dp
        delta_odd = 1.0_dp
        delta_even = 1.0_dp

        write(*,'(/,a)') TAB // 'Starting RPA O(Y^4) OBDM iteration ...'
        do while (Iteration < Iteration_max)
            ! 1) - 3) Rho_ref + dRho[Rho_old], diagonal approximation, linear mixing ...
            Rho_new = obdm_update(Rho_ref_it,Rho_old,X_it,Y_it)
            if (Diagonal_OBDM) Rho_new = O1b_diagonal(Rho_new)
            Rho_new%p = alpha * Rho_new%p + (1.0_dp - alpha) * Rho_old%p
            Rho_new%n = alpha * Rho_new%n + (1.0_dp - alpha) * Rho_old%n
            if (.not. (all(ieee_is_finite(Rho_new%p)) .and. all(ieee_is_finite(Rho_new%n)))) then
                write(*,'(/,a)') TAB // 'OBDM iteration stopped - diverging solution (Inf or NaN) ...'
                exit
            end if

            ! 4) Natural orbitals ... 5) rotation into the new NAT basis ...
            call obdm_natural(Rho_new,C_step)
            C_tot%p = matmul_dd(C_tot%p,C_step%p)
            C_tot%n = matmul_dd(C_tot%n,C_step%n)
            Rho_old = O1b_transformation(Rho_old,C_step)
            Rho_ref_it = O1b_transformation(Rho_ref_it,C_step)
            call RPA_OBDM_transform_amplitudes(Orb_1p1h,X_it,Y_it,C_step,X_tmp,Y_tmp)
            call move_alloc(X_tmp,X_it)
            call move_alloc(Y_tmp,Y_it)

            ! 6) Convergence measures ...
            delta_abs = obdm_norm(Orb,a_max,(Rho_new%p - Rho_old%p) / alpha,(Rho_new%n - Rho_old%n) / alpha)
            Norm_cor = obdm_norm(Orb,a_max,Rho_new%p - Rho_ref_it%p,Rho_new%n - Rho_ref_it%n)
            if (Norm_cor > 1.0e-12_dp) then
                delta = delta_abs / Norm_cor
            else
                delta = delta_abs
            end if
            Iteration = Iteration + 1
            write(*,'(a)') TAB // TAB // 'OBDM iteration:   Iteration number =  ' // cfmt_d(Iteration,4) // '   delta = ' // &
                           cfmt_f(delta,12,9)
            Rho_old = Rho_new

            if ((.not. ieee_is_finite(delta)) .or. delta_abs > Norm_ref) then
                write(*,'(/,a)') TAB // 'OBDM iteration stopped - diverging solution ...'
                exit
            end if
            if (delta < epsilon) then
                Converged = .true.
                exit
            end if
            if (Iteration > 100) then
                if (mod(Iteration,2) == 0) then
                    delta_prev = delta_even
                else
                    delta_prev = delta_odd
                end if
                if (abs(delta - delta_prev) < 1.0e-3_dp * delta) then
                    write(*,'(/,a)') TAB // 'OBDM iteration stopped - cycled degenerate solution found ...'
                    exit
                end if
                if (mod(Iteration,2) == 0) then
                    delta_even = delta
                else
                    delta_odd = delta
                end if
            end if
        end do

        if (Converged) then
            write(*,'(/,a)') TAB // 'RPA O(Y^4) OBDM iteration has converged ...'
            Rho = Rho_old
        else
            write(*,'(/,a)') TAB // 'RPA O(Y^4) OBDM iteration has not converged - using a single explicit evaluation of the OBDM ...'
            Rho_expl = obdm_update(Rho_ref,Rho_ini,X,Y)
            if (Diagonal_OBDM) Rho_expl = O1b_diagonal(Rho_expl)
            Rho = Rho_expl
            call obdm_natural(Rho,C_tot)
        end if

        ! New occupations (n,l,j of the input Orb) & transformation matrices ...
        do ia = 1, size(Orb)
            a = Orb(ia)%a
            Orb(ia)%pO = Rho%p(a,a)
            Orb(ia)%nO = Rho%n(a,a)
        end do
        C_HF_NAT%p = matmul_dd(C_HF_NAT%p,C_tot%p)
        C_HF_NAT%n = matmul_dd(C_HF_NAT%n,C_tot%n)
        C_LHO_NAT%p = matmul_dd(C_LHO_HF%p,C_HF_NAT%p)
        C_LHO_NAT%n = matmul_dd(C_LHO_HF%n,C_HF_NAT%n)

    contains

        ! Rho_ref + (oY4[Rho_o] - Rho_o) ...
        function obdm_update(Rho_r,Rho_o,X_u,Y_u) result(Rho_u)
            type(o1b_t), intent(in) :: Rho_r, Rho_o
            type(cmat_t), intent(in) :: X_u(:,0:), Y_u(:,0:)
            type(o1b_t) :: Rho_u, Rho_RPA
            Rho_RPA = RPA_OBDM_oY4(Params,Orb,Orb_ph,Orb_1p1h,Rho_o,X_u,Y_u)
            Rho_u%p = Rho_r%p + (Rho_RPA%p - Rho_o%p)
            Rho_u%n = Rho_r%n + (Rho_RPA%n - Rho_o%n)
        end function obdm_update

        ! C of OBDM_natural_orbitals (fixed Orb) & Rho <- diag(C' Rho C) of the unsplit Rho ...
        subroutine obdm_natural(Rho_io,C_o)
            type(o1b_t), intent(inout) :: Rho_io
            type(o1b_t), intent(out) :: C_o
            type(orb1b_t), allocatable :: Orb_tmp(:)
            type(o1b_t) :: Rho_d
            call OBDM_natural_orbitals(Params,Orb,Rho_io,Orb_tmp,C_o,Rho_d)
            Rho_io = O1b_diagonal(O1b_transformation(Rho_io,C_o))
        end subroutine obdm_natural

    end subroutine HF_ScRPA_OBDM_solve

    ! C = A B (dgemm) ...
    function matmul_dd(A,B) result(C)
        real(dp), intent(in) :: A(:,:), B(:,:)
        real(dp), allocatable :: C(:,:)
        allocate(C(size(A,1),size(B,2)))
        call dgemm('N','N',size(A,1),size(B,2),size(A,2),1.0_dp,A,size(A,1),B,size(B,1),0.0_dp,C,size(C,1))
    end function matmul_dd

end module nj_rpa_rho
