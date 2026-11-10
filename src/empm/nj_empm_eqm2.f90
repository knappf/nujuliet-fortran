! EMPM 2-phonon equations of motion for one J^P channel ... counterpart of eqm2_svd & eqm2_admat of the original EMPM codes,
    ! with the 1-phonon densities & redefined interaction computed on the fly (nj_empm_dens) instead of read from files ...
    ! Basis |alpha lambda; J> of all phonon pairs (alpha outer, lambda inner), the spurious states (alpha or lambda = cm) first ...
    ! Truncated rows ... all spurious states & the pairs with alpha >= lambda, E_alpha + E_lambda <= E2max & E_lambda <= E1max ...
    ! For a row (alpha,lambda) & a column (beta,lambda') of the full basis ...
        ! A = (E_a + E_l) d_ab d_ll' + d_ab R_ll' + d_ll' R_ab
        !   + (-1)^(J+Ja) sum_s {J Jl' Jb; s Ja Jl} sum_(i1,i2) rho_ab(s;i1,i2) V_ll'(s;i1,i2)
        ! D = d_ab d_ll' + d_al' d_bl {0 Jl Jl; J Ja Ja} sqrt((2Jl+1)(2Ja+1))
        !   - (-1)^(J+Jl'+Ja) sum_s {s Jl' Jl; J Ja Jb} sum_(i1,i2) (-1)^(j_i1-j_i2) rho_ab(s;i1,i2) rho_ll'(s;i2,i1)
        ! R_ll' = <l|A_TDA|l'> couplings to the spurious phonon (l /= l') ...
    ! Eigenproblem (A D^T) C = E D C on the linearly independent truncated rows (pivoted Cholesky of D) ... the spurious
        ! subspace S (SVD of the D elements of the spurious states) projected out exactly: physical states = eigenstates of
        ! P H P (P = 1 - S), spurious states = eigenstates of S H S, the couplings R = <physical|H|spurious> kept, so that
        ! physical & spurious states together are a complete set of the next phonon level (the shift of 1e8 MeV of
        ! ham_geev only in the limit) ...
module nj_empm_eqm2
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    use nj_empm_phonons, only : empm_phonons_t
    use nj_empm_dens
    use iso_fortran_env, only : int64
    !$ use omp_lib
    implicit none
    private

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
        subroutine dpstrf(Uplo,n,A,Lda,Piv,Rank,Tol,Work,Info)
            import :: dp
            character :: Uplo
            integer :: n, Lda, Piv(*), Rank, Info
            real(dp) :: A(Lda,*), Tol, Work(*)
        end subroutine dpstrf
        subroutine dgesvd(Jobu,Jobvt,m,n,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Info)
            import :: dp
            character :: Jobu, Jobvt
            integer :: m, n, Lda, Ldu, Ldvt, Lwork, Info
            real(dp) :: A(Lda,*), S(*), U(Ldu,*), Vt(Ldvt,*), Work(*)
        end subroutine dgesvd
        subroutine dgeqrf(m,n,A,Lda,Tau,Work,Lwork,Info)
            import :: dp
            integer :: m, n, Lda, Lwork, Info
            real(dp) :: A(Lda,*), Tau(*), Work(*)
        end subroutine dgeqrf
        subroutine dormqr(Side,Trans,m,n,k,A,Lda,Tau,C,Ldc,Work,Lwork,Info)
            import :: dp
            character :: Side, Trans
            integer :: m, n, k, Lda, Ldc, Lwork, Info
            real(dp) :: A(Lda,*), Tau(*), C(Ldc,*), Work(*)
        end subroutine dormqr
        subroutine dpotrs(Uplo,n,Nrhs,A,Lda,B,Ldb,Info)
            import :: dp
            character :: Uplo
            integer :: n, Nrhs, Lda, Ldb, Info
            real(dp) :: A(Lda,*), B(Ldb,*)
        end subroutine dpotrs
        subroutine dsyev(Jobz,Uplo,n,A,Lda,W,Work,Lwork,Info)
            import :: dp
            character :: Jobz, Uplo
            integer :: n, Lda, Lwork, Info
            real(dp) :: A(Lda,*), W(*), Work(*)
        end subroutine dsyev
        ! 64-bit integer LAPACK entry points of the LP64 MKL (the n^2 workspaces exceed 32-bit integers for n > ~32700) ...
        subroutine dsygvd_64(Itype,Jobz,Uplo,n,A,Lda,B,Ldb,W,Work,Lwork,Iwork,Liwork,Info)
            import :: dp, int64
            character :: Jobz, Uplo
            integer(int64) :: Itype, n, Lda, Ldb, Lwork, Liwork, Info, Iwork(*)
            real(dp) :: A(Lda,*), B(Ldb,*), W(*), Work(*)
        end subroutine dsygvd_64
        subroutine dgemm_64(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp, int64
            character :: Transa, Transb
            integer(int64) :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm_64
        subroutine dgemv_64(Trans,m,n,Alpha,A,Lda,x,Incx,Beta,y,Incy)
            import :: dp, int64
            character :: Trans
            integer(int64) :: m, n, Lda, Incx, Incy
            real(dp) :: Alpha, Beta, A(Lda,*), x(*), y(*)
        end subroutine dgemv_64
        subroutine dsyevd_64(Jobz,Uplo,n,A,Lda,W,Work,Lwork,Iwork,Liwork,Info)
            import :: dp, int64
            character :: Jobz, Uplo
            integer(int64) :: n, Lda, Lwork, Liwork, Info, Iwork(*)
            real(dp) :: A(Lda,*), W(*), Work(*)
        end subroutine dsyevd_64
    end interface

    real(dp), parameter :: Chol_Tol = 1.0e-6_dp, SVD_Tol = 1.0e-4_dp

    ! Wall-clock seconds spent in densities, interaction, row assembly & A D^T products of the current channel ...
    real(dp), save :: T_Cache = 0.0_dp, T_Row = 0.0_dp, T_H = 0.0_dp

    ! Full basis Al(k) (alpha), La(k) (lambda), Spur(k), K(alpha,lambda) full index ... Row(i) full index of truncated row i
        ! (2 phonons: the N_spur_tr spurious truncated rows first) ... Ks(m) full index of spurious state m = 1 .. N_spur ...
    type, public :: empm_basis2_t
        integer :: J = 0, P = 0, N_full = 0, N_spur = 0, N_tr = 0, N_spur_tr = 0
        integer, allocatable :: Al(:), La(:), K(:,:), Row(:), Ks(:)
        logical, allocatable :: Spur(:)
    end type empm_basis2_t

    ! Solution ... N_ind independent rows Ind(1:N_ind) (truncated indices), N_s spurious solutions (the last ones),
        ! energies E(1:N_ind) & amplitudes C(N_ind,N_ind) in the independent truncated basis (C^T D C = 1) ... the
        ! physical states (eigenstates of P H P) first, then the spurious ones (eigenstates of S H S), R(i,s) = <i|H|s> of
        ! physical i & spurious s ...
    type, public :: empm_solution2_t
        integer :: N_ind = 0, N_s = 0, N_spur_ind = 0
        integer, allocatable :: Ind(:)
        real(dp), allocatable :: E(:), C(:,:), R(:,:)
        real(dp) :: Asym = 0.0_dp
    end type empm_solution2_t

    public :: EMPM_basis2_make, EMPM_eqm2_build, EMPM_AD_rows, EMPM_eqm2_solve, EMPM_syevd, EMPM_chol_ind, EMPM_ind_solve, &
              EMPM_spur_directions, chol_piv

contains

    pure integer function par_sign(P)
        integer, intent(in) :: P
        par_sign = merge(1, -1, P == 1)
    end function par_sign

    ! Spur_Rows = 'all' keeps all spurious states as truncated rows, 'lambda' only those with lambda = spurious phonon
        ! (the variant of the original base.f90 used by its eqm2 runs), 'none' no spurious rows ... the pairs of the
        ! physical phonons span the whole physical subspace, the spurious states enter only as columns (the complete
        ! 1-phonon set of A D^T) & through their D elements Dsc (the constraints <cm lambda|psi> = 0) ...
    subroutine EMPM_basis2_make(Ph,J,P,E2max,E1max,Spur_Rows,Bs)
        type(empm_phonons_t), intent(in) :: Ph
        integer, intent(in) :: J, P
        real(dp), intent(in) :: E2max, E1max
        character(len=*), intent(in) :: Spur_Rows
        type(empm_basis2_t), intent(out) :: Bs
        integer :: a, l, k, Pass, n_s, n_o
        logical :: Sp

        Bs%J = J
        Bs%P = P
        allocate(Bs%K(Ph%N,Ph%N))
        Bs%K = 0
        ! Pass 1 counts, pass 2 fills ... spurious states first, both groups in the alpha-outer order ...
        do Pass = 1, 2
            n_s = 0
            n_o = 0
            do a = 1, Ph%N
                do l = 1, Ph%N
                    if (par_sign(Ph%P(a))*par_sign(Ph%P(l)) /= par_sign(P)) cycle
                    if (J < abs(Ph%J(a) - Ph%J(l)) .or. J > Ph%J(a) + Ph%J(l)) cycle
                    Sp = (a == Ph%i_cm .or. l == Ph%i_cm)
                    if (Sp) then
                        n_s = n_s + 1
                        k = n_s
                    else
                        n_o = n_o + 1
                        k = Bs%N_spur + n_o
                    end if
                    if (Pass == 2) then
                        Bs%Al(k) = a
                        Bs%La(k) = l
                        Bs%Spur(k) = Sp
                        Bs%K(a,l) = k
                    end if
                end do
            end do
            if (Pass == 1) then
                Bs%N_spur = n_s
                Bs%N_full = n_s + n_o
                allocate(Bs%Al(Bs%N_full), Bs%La(Bs%N_full), Bs%Spur(Bs%N_full))
            end if
        end do

        allocate(Bs%Row(Bs%N_full))
        Bs%N_tr = 0
        do k = 1, Bs%N_full
            a = Bs%Al(k)
            l = Bs%La(k)
            if (Bs%Spur(k)) then
                if (trim(Spur_Rows) == 'none') cycle
                if (trim(Spur_Rows) == 'lambda' .and. l /= Ph%i_cm) cycle
                Bs%N_tr = Bs%N_tr + 1
                Bs%N_spur_tr = Bs%N_spur_tr + 1
                Bs%Row(Bs%N_tr) = k
            else if (a >= l .and. Ph%E(a) + Ph%E(l) <= E2max .and. Ph%E(l) <= E1max) then
                Bs%N_tr = Bs%N_tr + 1
                Bs%Row(Bs%N_tr) = k
            end if
        end do
        if (all(trim(Spur_Rows) /= [character(len=6) :: 'all', 'lambda', 'none'])) call fatal('Invalid EMPM Spur_Rows = ' // trim(Spur_Rows))
        Bs%Row = Bs%Row(1:Bs%N_tr)
        Bs%Ks = [(k, k = 1, Bs%N_spur)]
    end subroutine EMPM_basis2_make

    ! H = A D^T, the truncated block of D & the D columns of the spurious states, with at most Mem_GB of A & D rows in memory ...
        ! Rows (sorted by lambda, so a block needs only the lambda cache of its own lambdas) are split into blocks of R rows ...
        ! block I computes A_I & D_I, the A_J & D_J of the later blocks are recomputed for each I ... H_IJ = A_I D_J^T &
        ! H_JI = A_J D_I^T are both computed, so that |AD - (AD)^T| checks every element of H ...
        ! Optional Vk(k,c) over the full basis gives Pk(i,c) = sum_k D_i(k) Vk(k,c) for the couplings to other subspaces ...
    subroutine EMPM_eqm2_build(Pairs,Ph,Bs,Mem_GB,Cache_GB,H,Dtr,Dsc,Vk,Pk)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_basis2_t), intent(in) :: Bs
        real(dp), intent(in) :: Mem_GB, Cache_GB
        real(dp), allocatable, intent(out) :: H(:,:), Dtr(:,:), Dsc(:,:)
        real(dp), intent(in), optional :: Vk(:,:)
        real(dp), allocatable, intent(out), optional :: Pk(:,:)
        real(dp), allocatable :: Ab(:,:), Db(:,:), Aj(:,:), Dj(:,:), G(:,:), Lam_key(:)
        integer, allocatable :: Order(:)
        integer :: N, R, n_blk, Ib, Jb, i0, i1, j0, j1, i
        integer(i8) :: t0

        N = Bs%N_tr
        allocate(Order(N), Lam_key(N))
        do i = 1, N
            Lam_key(i) = real(Bs%La(Bs%Row(i)),dp)
        end do
        call sortperm(Lam_key,Order)
        T_Cache = 0.0_dp
        T_Row = 0.0_dp
        T_H = 0.0_dp
        ! One block if all A & D rows fit into Mem_GB, otherwise blocks of R rows with A_I, D_I, A_J & D_J in memory ...
        if (16.0_dp*real(Bs%N_full,dp)*real(N,dp) <= Mem_GB*1.0e9_dp) then
            R = N
        else
            R = int(max(1.0_dp, Mem_GB*1.0e9_dp/(32.0_dp*real(Bs%N_full,dp))))
        end if
        n_blk = (N + R - 1)/R
        if (n_blk > 1) write(*,'(a)') TAB // 'Row blocks = ' // istr(n_blk) // ' of ' // istr(R) // ' rows (Mem_GB = ' // &
                                      cfmt_f(Mem_GB,0,1) // ') ... H_IJ & H_JI both computed, full symmetry check'
        allocate(H(N,N), Dtr(N,N), Dsc(Bs%N_spur,N))
        if (present(Vk)) allocate(Pk(N,size(Vk,2)))
        do Ib = 1, n_blk
            call timer_start(t0)
            i0 = (Ib - 1)*R + 1
            i1 = min(Ib*R, N)
            allocate(Ab(Bs%N_full,i1-i0+1), Db(Bs%N_full,i1-i0+1))
            call EMPM_AD_rows(Pairs,Ph,Bs,Order(i0:i1),Cache_GB,Db,Ab)
            call h_block(Ab,Db,Order(i0:i1),Order(i0:i1))
            call gather(Db,Order(i0:i1))
            ! Off-diagonal blocks H_IJ = A_I D_J^T & H_JI = A_J D_I^T, each element computed (no mirroring) ...
            do Jb = Ib + 1, n_blk
                j0 = (Jb - 1)*R + 1
                j1 = min(Jb*R, N)
                allocate(Aj(Bs%N_full,j1-j0+1), Dj(Bs%N_full,j1-j0+1))
                call EMPM_AD_rows(Pairs,Ph,Bs,Order(j0:j1),Cache_GB,Dj,Aj)
                call h_block(Ab,Dj,Order(i0:i1),Order(j0:j1))
                call h_block(Aj,Db,Order(j0:j1),Order(i0:i1))
                deallocate(Aj, Dj)
            end do
            deallocate(Ab, Db)
            if (n_blk > 1) call timer_report(t0,'row block ' // istr(Ib))
        end do
        write(*,'(a)') TAB // 'Seconds ... lambda cache (rho~ & V) ' // cfmt_f(T_Cache,0,1) // ' , alpha loop (rho_a & A, D rows) ' // &
                       cfmt_f(T_Row,0,1) // ' , A D^T ' // cfmt_f(T_H,0,1)

    contains

        subroutine h_block(Ab,Db,Ri,Rj)
            real(dp), intent(in) :: Ab(:,:), Db(:,:)
            integer, intent(in) :: Ri(:), Rj(:)
            call clock(T_H,.true.)
            allocate(G(size(Ri),size(Rj)))
            call dgemm('T','N',size(Ri),size(Rj),Bs%N_full,1.0_dp,Ab,Bs%N_full,Db,Bs%N_full,0.0_dp,G,size(Ri))
            H(Ri,Rj) = G
            deallocate(G)
            call clock(T_H,.false.)
        end subroutine h_block

        ! Dtr(i,j) = D_i(Row(j)) & Dsc(m,i) = D_i(Ks(m)) for the spurious states m = 1 .. N_spur ...
        subroutine gather(Db,Ri)
            real(dp), intent(in) :: Db(:,:)
            integer, intent(in) :: Ri(:)
            real(dp), allocatable :: Pq(:,:)
            integer :: q
            do q = 1, size(Ri)
                Dtr(Ri(q),:) = Db(Bs%Row,q)
                Dsc(:,Ri(q)) = Db(Bs%Ks,q)
            end do
            if (present(Vk)) then
                if (size(Vk,2) > 0) then
                    allocate(Pq(size(Ri),size(Vk,2)))
                    call dgemm('T','N',size(Ri),size(Vk,2),Bs%N_full,1.0_dp,Db,Bs%N_full,Vk,Bs%N_full,0.0_dp,Pq,size(Ri))
                    Pk(Ri,:) = Pq
                end if
            end if
        end subroutine gather

    end subroutine EMPM_eqm2_build

    ! Rows Rows(q) (truncated indices) of D & optionally A over the full basis ... Drow(:,q) & Arow(:,q) ...
        ! The lambdas of the rows are split into blocks whose V_l & rho~_l (cache) fit into Cache_GB ... for each block the
        ! cache is filled in parallel over lambda, then the rows are processed in parallel over alpha: rho_a once per alpha &
        ! block, then all rows (alpha, lambda of the block) ... everything inside the parallel loops is serial ...
    subroutine EMPM_AD_rows(Pairs,Ph,Bs,Rows,Cache_GB,Drow,Arow)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_basis2_t), intent(in) :: Bs
        integer, intent(in) :: Rows(:)
        real(dp), intent(in) :: Cache_GB
        real(dp), intent(out) :: Drow(:,:)
        real(dp), intent(out), optional :: Arow(:,:)
        type(phonon_block_t), allocatable :: Vc(:), Rc(:)
        type(phonon_block_t) :: Rho_l, Rho_a
        integer, allocatable :: Ch_bc(:,:), Lams(:), Slot(:), Qb(:), Order(:), G0(:)
        real(dp), allocatable :: Key(:)
        integer :: q, i, a, Jb, Pb, Jc, Pc, n_bc, nL, L0, L1, m, nq, ng, ig, iq
        real(dp) :: Bytes, Used
        integer(i8) :: t0
        logical :: Want_A
        integer, parameter :: Chunk = 4

        Want_A = present(Arow)
        ! Channel pairs (beta in (Jb,Pb), lambda' in (Jc,Pc)) of the basis states of J^P ...
        allocate(Ch_bc(4,4*(Ph%J_max + 1)**2))
        n_bc = 0
        do Jb = 0, Ph%J_max
            do Pb = 1, 2
                if (Ph%N_ch(Jb,Pb) == 0) cycle
                do Jc = abs(Bs%J - Jb), min(Bs%J + Jb, Ph%J_max)
                    do Pc = 1, 2
                        if (Ph%N_ch(Jc,Pc) == 0) cycle
                        if (par_sign(Pb)*par_sign(Pc) /= par_sign(Bs%P)) cycle
                        n_bc = n_bc + 1
                        Ch_bc(:,n_bc) = [Jb, Pb, Jc, Pc]
                    end do
                end do
            end do
        end do

        ! Distinct lambdas of the rows (ascending) ...
        allocate(Slot(Ph%N))
        Slot = 0
        do q = 1, size(Rows)
            Slot(Bs%La(Bs%Row(Rows(q)))) = 1
        end do
        Lams = pack([(i, i = 1, Ph%N)], Slot == 1)
        nL = size(Lams)

        L0 = 1
        do while (L0 <= nL)
            ! Block of lambdas within Cache_GB (at least one) ...
            Used = 0.0_dp
            L1 = L0 - 1
            do while (L1 < nL)
                Bytes = merge(2.0_dp, 1.0_dp, Want_A) * block_bytes(Lams(L1+1))
                if (L1 >= L0 .and. Used + Bytes > Cache_GB*1.0e9_dp) exit
                Used = Used + Bytes
                L1 = L1 + 1
            end do
            Slot = 0
            do m = L0, L1
                Slot(Lams(m)) = m - L0 + 1
            end do

            ! Cache V_l & rho~_l of the block ...
            call timer_start(t0)
            allocate(Rc(L1-L0+1))
            if (Want_A) allocate(Vc(L1-L0+1))
            !$omp parallel do schedule(dynamic) private(Rho_l)
            do m = L0, L1
                call EMPM_density(Pairs,Ph,Lams(m),Rho_l)
                call EMPM_density_transposed(Pairs,Rho_l,Rc(m-L0+1))
                if (Want_A) call EMPM_interaction(Pairs,Ph,Rho_l,Vc(m-L0+1))
            end do
            !$omp end parallel do
            T_Cache = T_Cache + elapsed(t0)

            ! Rows of the block grouped by alpha ...
            call timer_start(t0)
            nq = 0
            allocate(Qb(size(Rows)))
            do q = 1, size(Rows)
                if (Slot(Bs%La(Bs%Row(Rows(q)))) > 0) then
                    nq = nq + 1
                    Qb(nq) = q
                end if
            end do
            allocate(Key(nq), Order(nq))
            do iq = 1, nq
                Key(iq) = real(Bs%Al(Bs%Row(Rows(Qb(iq)))),dp)
            end do
            call sortperm(Key,Order)
            Qb(1:nq) = Qb(Order)
            ! Task starts G0(1:ng+1) ... rows of one alpha, at most Chunk of them per task (rho_a is cheap compared to the
                ! rows, large alpha groups would otherwise serialize the loop) ...
            allocate(G0(nq+1))
            ng = 0
            do iq = 1, nq
                if (iq == 1) then
                    ng = ng + 1
                    G0(ng) = iq
                else if (Key(Order(iq)) /= Key(Order(iq-1)) .or. iq - G0(ng) >= Chunk) then
                    ng = ng + 1
                    G0(ng) = iq
                end if
            end do
            G0(ng+1) = nq + 1
            !$omp parallel do schedule(dynamic) private(a,Rho_a,iq,q,m)
            do ig = 1, ng
                a = Bs%Al(Bs%Row(Rows(Qb(G0(ig)))))
                call EMPM_density(Pairs,Ph,a,Rho_a)
                do iq = G0(ig), G0(ig+1) - 1
                    q = Qb(iq)
                    m = Slot(Bs%La(Bs%Row(Rows(q))))
                    if (Want_A) then
                        call row_AD(Rho_a,Vc(m),Rc(m),q,a,Bs%La(Bs%Row(Rows(q))))
                    else
                        call row_AD(Rho_a,Rc(m),Rc(m),q,a,Bs%La(Bs%Row(Rows(q))))
                    end if
                end do
            end do
            !$omp end parallel do
            T_Row = T_Row + elapsed(t0)
            deallocate(Rc, Qb, Key, Order, G0)
            if (Want_A) deallocate(Vc)
            L0 = L1 + 1
        end do

    contains

        ! Bytes of a density-type block of phonon l ...
        real(dp) function block_bytes(l)
            integer, intent(in) :: l
            integer :: Jb, Pb, Sig
            block_bytes = 0.0_dp
            do Jb = 0, Ph%J_max
                do Pb = 1, 2
                    do Sig = abs(Ph%J(l) - Jb), min(Ph%J(l) + Jb, Pairs%L_max)
                        block_bytes = block_bytes + 8.0_dp*real(Pairs%Na(Sig,block_parity(Ph,l,Pb)),dp)*real(Ph%N_ch(Jb,Pb),dp)
                    end do
                end do
            end do
        end function block_bytes

        ! Row q = (al,la) ... V_l & Rt_l the interaction & transposed density of la (V_l unused without A) ...
        subroutine row_AD(Rho_a,V_l,Rt_l,q,al,la)
            type(phonon_block_t), intent(in) :: Rho_a, V_l, Rt_l
            integer, intent(in) :: q, al, la
            real(dp), allocatable :: G(:,:), GD(:,:)
            integer :: Ja, Jl, J, Jb, Pb, Jc, Pc, nb, nc, lo, hi, Sig, ib, ic, k, b, c, Ka, Pp, m
            real(dp) :: wgtA, wgtD

            Ja = Ph%J(al)
            Jl = Ph%J(la)
            J = Bs%J
            if (Want_A) Arow(:,q) = 0.0_dp
            Drow(:,q) = 0.0_dp
            ! G = sum_s w_s rho_a(:,s,b)^T X_l(:,s,c) ... one accumulating dgemm per sigma with the 6j weight as alpha ...
            do m = 1, n_bc
                Jb = Ch_bc(1,m)
                Pb = Ch_bc(2,m)
                Jc = Ch_bc(3,m)
                Pc = Ch_bc(4,m)
                lo = max(abs(Ja - Jb), abs(Jl - Jc))
                hi = min(Ja + Jb, Jl + Jc, Pairs%L_max)
                if (lo > hi) cycle
                nb = Ph%N_ch(Jb,Pb)
                nc = Ph%N_ch(Jc,Pc)
                Pp = Rho_a%Pi(Jb,Pb)
                allocate(GD(nb,nc))
                GD = 0.0_dp
                if (Want_A) then
                    allocate(G(nb,nc))
                    G = 0.0_dp
                end if
                do Sig = lo, hi
                    wgtD = -phase(J + Jc + Ja) * f6j(2*Sig,2*Jc,2*Jl,2*J,2*Ja,2*Jb)
                    Ka = Pairs%Na(Sig,Pp)
                    if (Ka == 0) cycle
                    if (wgtD /= 0.0_dp) call dgemm('T','N',nb,nc,Ka,wgtD,Rho_a%B(Jb,Pb)%M(Sig)%v,Ka, &
                                                   Rt_l%B(Jc,Pc)%M(Sig)%v,Ka,1.0_dp,GD,nb)
                    if (Want_A) then
                        wgtA = phase(J + Ja) * f6j(2*J,2*Jc,2*Jb,2*Sig,2*Ja,2*Jl)
                        if (wgtA /= 0.0_dp) call dgemm('T','N',nb,nc,Ka,wgtA,Rho_a%B(Jb,Pb)%M(Sig)%v,Ka, &
                                                       V_l%B(Jc,Pc)%M(Sig)%v,Ka,1.0_dp,G,nb)
                    end if
                end do
                ! Distinct channel pairs fill distinct columns k = K(beta,lambda') ...
                do ic = 1, nc
                    c = Ph%First(Jc,Pc) + ic - 1
                    do ib = 1, nb
                        b = Ph%First(Jb,Pb) + ib - 1
                        k = Bs%K(b,c)
                        Drow(k,q) = Drow(k,q) + GD(ib,ic)
                        if (Want_A) Arow(k,q) = Arow(k,q) + G(ib,ic)
                    end do
                end do
                deallocate(GD)
                if (Want_A) deallocate(G)
            end do

            ! Unperturbed energies, spurious couplings & the metric of the overcomplete basis ...
            k = Bs%K(al,la)
            Drow(k,q) = Drow(k,q) + 1.0_dp
            k = Bs%K(la,al)
            if (k > 0) Drow(k,q) = Drow(k,q) + f6j(0,2*Jl,2*Jl,2*J,2*Ja,2*Ja) * sqrt(real((2*Jl + 1)*(2*Ja + 1),dp))
            if (.not. Want_A) return
            k = Bs%K(al,la)
            Arow(k,q) = Arow(k,q) + Ph%E(al) + Ph%E(la)
            if (Ph%i_cm > 0) then
                do c = 1, Ph%N
                    if (c == la) cycle
                    if (la == Ph%i_cm .or. c == Ph%i_cm) then
                        k = Bs%K(al,c)
                        if (k > 0) Arow(k,q) = Arow(k,q) + merge(Ph%R(c), Ph%R(la), la == Ph%i_cm)
                    end if
                end do
                do b = 1, Ph%N
                    if (b == al) cycle
                    if (al == Ph%i_cm .or. b == Ph%i_cm) then
                        k = Bs%K(b,la)
                        if (k > 0) Arow(k,q) = Arow(k,q) + merge(Ph%R(b), Ph%R(al), al == Ph%i_cm)
                    end if
                end do
            end if
        end subroutine row_AD

    end subroutine EMPM_AD_rows

    ! Wall-clock seconds since t0 ...
    real(dp) function elapsed(t0)
        integer(i8), intent(in) :: t0
        integer(i8) :: c, Rate
        call system_clock(c,Rate)
        elapsed = real(c - t0,dp)/real(Rate,dp)
    end function elapsed

    ! Accumulates the wall-clock time between the calls with Start = .true. & .false. ...
    subroutine clock(T,Start)
        real(dp), intent(inout) :: T
        logical, intent(in) :: Start
        integer(i8) :: c, Rate
        call system_clock(c,Rate)
        if (Start) then
            T = T - real(c,dp)/real(Rate,dp)
        else
            T = T + real(c,dp)/real(Rate,dp)
        end if
    end subroutine clock

    ! Linear independence, spurious subspace & the generalized eigenproblem ... the steps of eqm.f90 & ham_geev ...
    ! Legacy_Spur reproduces the selection of the independent spurious states of the original eqm.f90: pivoted Cholesky of
        ! D(truncated rows 1 .. N_spur, full columns 1 .. N_spur), which mixes non-spurious rows in when not all spurious states
        ! are truncated rows (Spur_Rows = 'lambda') ... for one-to-one comparisons with the original runs only ...
        ! D_Spectrum ... eigenvalue spectrum of the whole truncated D in the log (a dense n^3 diagnostic) ...
    subroutine EMPM_eqm2_solve(Bs,H,Dtr,Dsc,Legacy_Spur,D_Spectrum,Sol)
        type(empm_basis2_t), intent(in) :: Bs
        real(dp), intent(in) :: H(:,:), Dtr(:,:), Dsc(:,:)
        logical, intent(in) :: Legacy_Spur, D_Spectrum
        type(empm_solution2_t), intent(out) :: Sol
        real(dp), allocatable :: Dd(:,:), Work(:), Dsp(:,:), Vs(:,:), Hr(:,:), Dr(:,:)
        integer, allocatable :: Piv(:), Nx_s(:), Ind_s(:), Sp_Rows(:), Mk(:)
        integer :: N, i, j, Rank, Info, n_s, n_ind
        real(dp) :: Asym

        N = Bs%N_tr
        ! Linearly independent spurious states ... pivoted Cholesky of their D block ... Nx_s(m) flags the spurious state
            ! m = 1 .. N_spur (full-basis index Ks(m)) ...
        allocate(Nx_s(Bs%N_spur))
        Nx_s = 0
        ! Spurious truncated rows Sp_Rows (anywhere among the rows) ...
        Sp_Rows = pack([(i, i = 1, N)], Bs%Spur(Bs%Row))
        if (Legacy_Spur .and. Bs%N_spur > 0) then
            if (N < Bs%N_spur) call fatal('Legacy_Spur needs at least N_spur truncated rows')
            Dd = transpose(Dsc(1:Bs%N_spur,1:Bs%N_spur))
            allocate(Piv(Bs%N_spur), Work(2*Bs%N_spur))
            call dpstrf('U',Bs%N_spur,Dd,Bs%N_spur,Piv,Rank,Chol_Tol,Work,Info)
            call chol_rank(Dd,Rank)
            Nx_s(Piv(1:Rank)) = 1
            deallocate(Piv, Work)
        else if (Bs%N_spur_tr > 0) then
            ! Mk(k) = m of the full-basis spurious state k = Ks(m) ...
            allocate(Mk(Bs%N_full))
            Mk = 0
            Mk(Bs%Ks) = [(i, i = 1, Bs%N_spur)]
            Dd = Dtr(Sp_Rows,Sp_Rows)
            allocate(Piv(Bs%N_spur_tr), Work(2*Bs%N_spur_tr))
            call dpstrf('U',Bs%N_spur_tr,Dd,Bs%N_spur_tr,Piv,Rank,Chol_Tol,Work,Info)
            call chol_rank(Dd,Rank)
            Nx_s(Mk(Bs%Row(Sp_Rows(Piv(1:Rank))))) = 1
            deallocate(Piv, Work)
        else
            ! No spurious rows (Spur_Rows = 'none') ... all spurious states, their linear dependence is left to the SVD ...
            Nx_s = 1
        end if
        Sol%N_spur_ind = sum(Nx_s)
        if (Bs%N_spur_tr > 0) call spectrum_report('spurious D block',Dtr(Sp_Rows,Sp_Rows))
        if (D_Spectrum) call spectrum_report('truncated D',Dtr)

        ! Linearly independent truncated states ...
        Dd = Dtr
        allocate(Piv(N), Work(2*N))
        call dpstrf('U',N,Dd,N,Piv,Rank,Chol_Tol,Work,Info)
        call chol_rank(Dd,Rank)
        n_ind = Rank
        Sol%N_ind = n_ind
        allocate(Sol%Ind(n_ind))
        Sol%Ind = Piv(1:Rank)
        call sort_int(Sol%Ind)
        deallocate(Piv, Work, Dd)

        Asym = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(max:Asym)
        do j = 1, N
            do i = 1, j - 1
                Asym = max(Asym, abs(H(i,j) - H(j,i)))
            end do
        end do
        !$omp end parallel do
        Sol%Asym = Asym

        ! Spurious directions ... SVD of the D elements between the independent rows & the independent spurious states
            ! (Dsc(m,:) of the spurious state m), projected out exactly by spur_gev ...
        allocate(Ind_s(Bs%N_spur))
        n_s = 0
        do i = 1, Bs%N_spur
            if (Nx_s(i) == 1) then
                n_s = n_s + 1
                Ind_s(n_s) = i
            end if
        end do
        Hr = H(Sol%Ind,Sol%Ind)
        Dr = Dtr(Sol%Ind,Sol%Ind)
        allocate(Dsp(n_s,n_ind))
        do j = 1, n_ind
            do i = 1, n_s
                Dsp(i,j) = Dsc(Ind_s(i),Sol%Ind(j))
            end do
        end do
        call EMPM_spur_directions(Dsp,Vs)
        call spur_gev(Hr,Dr,Vs,Sol)

    contains

        ! Rank of the pivoted Cholesky factor U ... only pivots U(k,k)^2 above Chol_Tol (dpstrf of MKL returns rank 1 for a
            ! numerically zero D, e.g. the Pauli-forbidden 3-phonon channels of 4He with J = 8, 9) ...
        subroutine chol_rank(U,Rank)
            real(dp), intent(in) :: U(:,:)
            integer, intent(inout) :: Rank
            do while (Rank > 0)
                if (U(Rank,Rank)**2 > Chol_Tol) exit
                Rank = Rank - 1
            end do
        end subroutine chol_rank

        ! Eigenvalue counts of a D block per decade ... the gap below Chol_Tol shows how sharp the linear dependence is ...
        subroutine spectrum_report(Name,M)
            character(len=*), intent(in) :: Name
            real(dp), intent(in) :: M(:,:)
            real(dp), allocatable :: Ev(:), Ms(:,:), Wk(:)
            real(dp) :: Wq1(1)
            integer :: d, n, Inf
            character(len=:), allocatable :: Line
            n = size(M,1)
            allocate(Ev(n))
            Ms = M
            call dsyev('N','U',n,Ms,n,Ev,Wq1,-1,Inf)
            allocate(Wk(int(Wq1(1))))
            call dsyev('N','U',n,Ms,n,Ev,Wk,size(Wk),Inf)
            Line = TAB // 'Eigenvalues of the ' // Name // ' ... min ' // cfmt_e(minval(Ev),0,2) // ' , max ' // &
                   cfmt_e(maxval(Ev),0,2) // ' , counts in [1e-k,1e-k+1) for k = 1 .. 10 :'
            do d = 1, 10
                Line = Line // ' ' // istr(count(Ev >= 10.0_dp**(-d) .and. Ev < 10.0_dp**(1-d)))
            end do
            write(*,'(a)') Line // ' , below 1e-10 : ' // istr(count(Ev < 1.0e-10_dp))
        end subroutine spectrum_report

        subroutine sort_int(v)
            integer, intent(inout) :: v(:)
            integer :: a, b, t
            do a = 2, size(v)
                t = v(a)
                b = a - 1
                do while (b >= 1)
                    if (v(b) <= t) exit
                    v(b+1) = v(b)
                    b = b - 1
                end do
                v(b+1) = t
            end do
        end subroutine sort_int

    end subroutine EMPM_eqm2_solve

    ! A C = B C diag(E), C^T B C = 1 ... eigenvectors C in A on return (dsygvd_64) ...
    subroutine sygvd(A,B,E)
        real(dp), intent(inout) :: A(:,:), B(:,:)
        real(dp), intent(out) :: E(:)
        real(dp), allocatable :: Work(:)
        integer(int64), allocatable :: Iwork(:)
        real(dp) :: Wq(1)
        integer(int64) :: n, Info, Iwq(1)
        n = size(A,1,kind=int64)
        if (n == 0) return
        call dsygvd_64(1_int64,'V','U',n,A,n,B,n,E,Wq,-1_int64,Iwq,-1_int64,Info)
        allocate(Work(int(Wq(1),int64)), Iwork(Iwq(1)))
        call dsygvd_64(1_int64,'V','U',n,A,n,B,n,E,Work,size(Work,kind=int64),Iwork,size(Iwork,kind=int64),Info)
        if (Info /= 0) call fatal('dsygvd failed in the EMPM 2-phonon eigenproblem ... info = ' // istr(int(Info)))
    end subroutine sygvd

    ! Independent rows of the truncated D (Ind_First) ... Dtr(:,i) = D row i (not modified), pivoted Cholesky chol_piv
        ! (pivots with remaining diagonal > Chol_Tol, as dpstrf & chol_rank of EMPM_eqm2_solve) ... Ind ascending (= sorted
        ! by the (n-1)-phonon state a), Dii = D(Ind,Ind) with D(i,j) = row i, column j as Dtr(Ind,Ind) of EMPM_eqm2_solve ...
    subroutine EMPM_chol_ind(Dtr,Ind,Dii)
        real(dp), intent(in) :: Dtr(:,:)
        integer, allocatable, intent(out) :: Ind(:)
        real(dp), allocatable, intent(out) :: Dii(:,:)
        integer, allocatable :: Piv(:)
        logical, allocatable :: Sel(:)
        integer :: N, Rank, i, j, r

        N = size(Dtr,1)
        allocate(Sel(N))
        call chol_piv(N,Dtr,Chol_Tol,Piv,Rank)
        Sel = .false.
        Sel(Piv(1:Rank)) = .true.
        allocate(Ind(Rank), Dii(Rank,Rank))
        r = 0
        do i = 1, N
            if (.not. Sel(i)) cycle
            r = r + 1
            Ind(r) = i
        end do
        !$omp parallel do schedule(dynamic)
        do j = 1, Rank
            do i = 1, j
                Dii(i,j) = Dtr(Ind(j),Ind(i))
                Dii(j,i) = Dii(i,j)
            end do
        end do
        !$omp end parallel do
    end subroutine EMPM_chol_ind

    ! Pivoted Cholesky of the symmetric positive semidefinite A (not modified, both triangles used as columns) ... greedy
        ! pivoting on the largest remaining diagonal d (as dpstrf), stop at max d <= Tol, only the pivots Piv(1:Rank)
        ! (original indices) are returned ... left-looking with a pool of the Np rows of largest d (the likely pivots), their
        ! Schur-complement columns A(:,pool) - L L(pool,:)^T kept up to date with one dgemm per panel of Nb pivots, within
        ! a panel a pivot column only needs the panel's own columns (dgemv) ... after each panel the pool is refilled with
        ! the current top rows (new columns by one dgemm), a pivot outside the pool (miss) gets its full column by dgemv ...
        ! cost ~ N Rank (Np + Rank/2) instead of N^2 Rank of the right-looking factorization ... L (live rows x pivots) is
        ! kept separately, dead rows (d <= Tol, d never grows) & pivots are compacted out once they are 1/8 of the rows ...
    subroutine chol_piv(n,A,Tol,Piv,Rank)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        real(dp), intent(in) :: Tol
        integer, allocatable, intent(out) :: Piv(:)
        integer, intent(out) :: Rank
        integer, parameter :: Nb = 64, Np_max = 2048
        real(dp), parameter :: Gone = -huge(1.0_dp)
        real(dp), allocatable :: L(:,:), d(:), Sc(:,:), S2(:,:), Lc(:,:), s(:), Tmp(:,:)
        integer, allocatable :: Live(:), Pr(:), P2(:), Tc(:), Keep(:), Pv(:), Ord(:), Inv(:)
        integer :: m, k, k0, t, p, i, j, Cap, nk, Miss, np, n2, nkeep, nt
        real(dp) :: x
        logical :: Done

        allocate(d(n), Live(n), Tc(n), s(n), Pv(n), Pr(0), Sc(n,0))
        do i = 1, n
            d(i) = A(i,i)
            Live(i) = i
        end do
        m = n
        k = 0
        Miss = 0
        np = 0
        Cap = min(n, max(Nb, n/4))
        allocate(L(n,Cap))
        Tc = 0
        Done = .false.
        do
            call refill()
            if (np == 0) exit
            k0 = k
            ! Panel ...
            do while (k - k0 < Nb)
                p = maxloc(d(1:m),1)
                if (d(p) <= Tol) then
                    Done = .true.
                    exit
                end if
                if (k + 1 > Cap) call grow()
                if (Tc(p) > 0) then
                    s(1:m) = Sc(1:m,Tc(p))
                    if (k > k0) call dgemv_64('N',int(m,int64),int(k-k0,int64),-1.0_dp,L(1,k0+1),int(n,int64),L(p,k0+1), &
                                              int(n,int64),1.0_dp,s,1_int64)
                else
                    Miss = Miss + 1
                    do i = 1, m
                        s(i) = A(Live(i),Live(p))
                    end do
                    if (k > 0) call dgemv_64('N',int(m,int64),int(k,int64),-1.0_dp,L,int(n,int64),L(p,1),int(n,int64), &
                                             1.0_dp,s,1_int64)
                end if
                k = k + 1
                x = 1.0_dp / sqrt(d(p))
                !$omp parallel do if (m > 20000)
                do i = 1, m
                    L(i,k) = s(i) * x
                    d(i) = d(i) - L(i,k)**2
                end do
                !$omp end parallel do
                Pv(k) = Live(p)
                d(p) = Gone
            end do
            if (Done .or. k == k0) exit
            ! Pool columns brought up to date with the panel ...
            if (np > 0) then
                Lc = L(Pr(1:np),k0+1:k)
                call dgemm_64('N','T',int(m,int64),int(np,int64),int(k-k0,int64),-1.0_dp,L(1,k0+1),int(n,int64),Lc, &
                              int(np,int64),1.0_dp,Sc,int(n,int64))
                deallocate(Lc)
            end if
            ! Compaction ... dead rows (d <= Tol) & pivots out of L, d, Live & the pool ...
            nk = count(d(1:m) > Tol)
            if (m - nk > max(Nb, m/8)) then
                Keep = pack([(i, i = 1, m)], d(1:m) > Tol)
                allocate(Inv(m))
                Inv = 0
                Inv(Keep) = [(i, i = 1, nk)]
                ! rows in ascending order ... Keep(i) >= i, column by column in place ...
                !$omp parallel do
                do j = 1, k
                    L(1:nk,j) = L(Keep,j)
                end do
                !$omp end parallel do
                !$omp parallel do
                do j = 1, np
                    Sc(1:nk,j) = Sc(Keep,j)
                end do
                !$omp end parallel do
                Tc(1:m) = 0
                nt = 0
                do t = 1, np
                    if (Inv(Pr(t)) == 0) cycle
                    nt = nt + 1
                    if (nt /= t) Sc(1:nk,nt) = Sc(1:nk,t)
                    Pr(nt) = Inv(Pr(t))
                    Tc(Pr(nt)) = nt
                end do
                np = nt
                d(1:nk) = d(Keep)
                Live(1:nk) = Live(Keep)
                m = nk
                deallocate(Inv)
            end if
        end do
        Rank = k
        Piv = Pv(1:k)
        if (Miss > 0) write(*,'(a)') TAB // 'Pivoted Cholesky ... ' // istr(Miss) // ' pivots outside the pool'

    contains

        ! Pool ... the live rows among the Np_max of largest d (above Tol) join, members stay while among the 2 Np_max
            ! largest (less recomputation near the edge) ... columns of the kept members copied, the new ones
            ! A(:,new) - L L(new,:)^T by one dgemm ...
        subroutine refill()
            integer :: nn, nl
            Tc(1:m) = 0
            do t = 1, np
                if (d(Pr(t)) > Tol) Tc(Pr(t)) = t
            end do
            nl = count(d(1:m) > Tol)
            allocate(Ord(m))
            call sortperm_rev(d(1:m),Ord)
            ! hysteresis ... members stay while among the top 2 Np_max rows, new rows enter from the top Np_max ...
            nkeep = count(Tc(Ord(1:min(2*Np_max,nl))) > 0)
            n2 = nkeep + count(Tc(Ord(1:min(Np_max,nl))) == 0)
            allocate(P2(n2), S2(n,n2))
            nkeep = 0
            do t = 1, min(2*Np_max,nl)
                if (Tc(Ord(t)) > 0) then
                    nkeep = nkeep + 1
                    P2(nkeep) = Ord(t)
                end if
            end do
            nn = nkeep
            do t = 1, min(Np_max,nl)
                if (Tc(Ord(t)) == 0) then
                    nn = nn + 1
                    P2(nn) = Ord(t)
                end if
            end do
            !$omp parallel do
            do t = 1, nkeep
                S2(1:m,t) = Sc(1:m,Tc(P2(t)))
            end do
            !$omp end parallel do
            !$omp parallel do private(i)
            do t = nkeep + 1, n2
                do i = 1, m
                    S2(i,t) = A(Live(i),Live(P2(t)))
                end do
            end do
            !$omp end parallel do
            if (k > 0 .and. n2 > nkeep) then
                Lc = L(P2(nkeep+1:n2),1:k)
                call dgemm_64('N','T',int(m,int64),int(n2-nkeep,int64),int(k,int64),-1.0_dp,L,int(n,int64),Lc, &
                              int(n2-nkeep,int64),1.0_dp,S2(1,nkeep+1),int(n,int64))
                deallocate(Lc)
            end if
            Tc(1:m) = 0
            do t = 1, n2
                Tc(P2(t)) = t
            end do
            call move_alloc(S2,Sc)
            call move_alloc(P2,Pr)
            np = n2
            deallocate(Ord)
        end subroutine refill

        ! L with 1.5 x the columns ...
        subroutine grow()
            Cap = min(n, int(1.5_dp*real(Cap,dp)) + Nb)
            allocate(Tmp(n,Cap))
            !$omp parallel do
            do j = 1, k
                Tmp(1:m,j) = L(1:m,j)
            end do
            !$omp end parallel do
            call move_alloc(Tmp,L)
        end subroutine grow

    end subroutine chol_piv

    ! Generalized eigenproblem (H, Dii) of the independent rows (Ind_First) ... |H - H^T|max as Asym, C^T Dii C = 1 ...
        ! Dsp(m,j) = D elements of the spurious states m & the independent rows j (optional), the spurious subspace projected
        ! out by spur_gev ... H & Dii destroyed ...
    subroutine EMPM_ind_solve(H,Dii,Sol,Dsp)
        real(dp), allocatable, intent(inout) :: H(:,:), Dii(:,:)
        type(empm_solution2_t), intent(inout) :: Sol
        real(dp), intent(in), optional :: Dsp(:,:)
        real(dp), allocatable :: Vs(:,:)
        real(dp) :: Asym
        integer :: i, j

        Asym = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(max:Asym)
        do j = 1, size(H,2)
            do i = 1, j - 1
                Asym = max(Asym, abs(H(i,j) - H(j,i)))
            end do
        end do
        !$omp end parallel do
        Sol%Asym = Asym
        Sol%N_spur_ind = 0
        if (present(Dsp)) then
            Sol%N_spur_ind = size(Dsp,1)
            call EMPM_spur_directions(Dsp,Vs)
        else
            allocate(Vs(0,size(H,1)))
        end if
        call spur_gev(H,Dii,Vs,Sol)
        deallocate(H, Dii)
    end subroutine EMPM_ind_solve

    ! Orthonormal spurious directions Vs(N_s,n) ... right singular vectors of Dsp(m,j) (D elements of the spurious states m
        ! & the independent rows j) with singular values above SVD_Tol ... the physical states x obey Vs x = 0 ...
    subroutine EMPM_spur_directions(Dsp,Vs)
        real(dp), intent(in) :: Dsp(:,:)
        real(dp), allocatable, intent(out) :: Vs(:,:)
        real(dp), allocatable :: A(:,:), Sv(:), U(:,:), Vt(:,:), Work(:)
        real(dp) :: Wq(1)
        integer :: m, n, k, Info

        m = size(Dsp,1)
        n = size(Dsp,2)
        k = min(m,n)
        if (k == 0) then
            allocate(Vs(0,n))
            return
        end if
        A = Dsp
        allocate(Sv(k), U(1,1), Vt(k,n))
        call dgesvd('N','S',m,n,A,m,Sv,U,1,Vt,k,Wq,-1,Info)
        allocate(Work(int(Wq(1))))
        call dgesvd('N','S',m,n,A,m,Sv,U,1,Vt,k,Work,size(Work),Info)
        if (Info /= 0) call fatal('dgesvd failed in the EMPM spurious subspace')
        Vs = Vt(1:count(abs(Sv) > SVD_Tol),:)
    end subroutine EMPM_spur_directions

    ! Generalized eigenproblem (H, D) on the independent rows with the spurious directions Vs (N_s x n, orthonormal rows)
        ! projected out exactly ... Vs^T = Q R (Householder), in y = Q^T x the physical subspace is y(1:N_s) = 0 ... physical
        ! states [0; Y_p] from (H_pp, D_pp), spurious states [I; Z] U_s with Z = -D_pp^-1 D_ps (D-orthogonal to all
        ! physical ones) & U_s from the N_s x N_s problem ([I;Z]^T H [I;Z], [I;Z]^T D [I;Z]) ... Sol%E = [E_p ; E_s],
        ! Sol%C = Q [Y_p , Y_s] (C^T D C = 1), Sol%R = Y_p^T H Y_s, Sol%N_s = N_s ... H & D destroyed ...
    subroutine spur_gev(H,D,Vs,Sol)
        real(dp), intent(inout) :: H(:,:), D(:,:)
        real(dp), intent(in) :: Vs(:,:)
        type(empm_solution2_t), intent(inout) :: Sol
        real(dp), allocatable :: Qf(:,:), Tau(:), Work(:), Hp(:,:), Dpp(:,:), Z(:,:), G(:,:), Hs(:,:), Ms(:,:), Ep(:), Es(:), &
                                 Gu(:,:), Y(:,:)
        real(dp) :: Wq(1)
        integer :: n, ns, np, Info

        n = size(H,1)
        ns = size(Vs,1)
        np = n - ns
        Sol%N_s = ns
        if (ns == 0) then
            allocate(Sol%E(n), Sol%R(n,0))
            call sygvd(H,D,Sol%E)
            Sol%C = H
            return
        end if
        Qf = transpose(Vs)
        allocate(Tau(ns))
        call dgeqrf(n,ns,Qf,n,Tau,Wq,-1,Info)
        allocate(Work(int(Wq(1))))
        call dgeqrf(n,ns,Qf,n,Tau,Work,size(Work),Info)
        if (Info /= 0) call fatal('dgeqrf failed in the EMPM spurious subspace')
        deallocate(Work)
        call apply_q('L','T',H)
        call apply_q('R','N',H)
        call apply_q('L','T',D)
        call apply_q('R','N',D)
        ! Physical states ... Hp -> Y_p, Dpp -> its Cholesky factor U (D_pp = U^T U) ...
        Hp = H(ns+1:n,ns+1:n)
        Dpp = D(ns+1:n,ns+1:n)
        allocate(Ep(np))
        call sygvd(Hp,Dpp,Ep)
        ! Spurious states ...
        Z = -D(ns+1:n,1:ns)
        if (np > 0) call dpotrs('U',np,ns,Dpp,np,Z,np,Info)
        if (Info /= 0) call fatal('dpotrs failed in the EMPM spurious subspace')
        G = H(:,1:ns)
        if (np > 0) call dgemm('N','N',n,ns,np,1.0_dp,H(:,ns+1:n),n,Z,np,1.0_dp,G,n)
        Hs = G(1:ns,:)
        Ms = D(1:ns,1:ns)
        if (np > 0) then
            call dgemm('T','N',ns,ns,np,1.0_dp,Z,np,G(ns+1,1),n,1.0_dp,Hs,ns)
            call dgemm('T','N',ns,ns,np,1.0_dp,D(ns+1:n,1:ns),np,Z,np,1.0_dp,Ms,ns)
        end if
        Hs = 0.5_dp*(Hs + transpose(Hs))
        Ms = 0.5_dp*(Ms + transpose(Ms))
        allocate(Es(ns))
        call sygvd(Hs,Ms,Es)
        ! Couplings R = Y_p^T (H [I;Z] U_s)_(physical rows) ...
        allocate(Gu(np,ns), Sol%R(np,ns))
        if (np > 0) then
            call dgemm('N','N',np,ns,ns,1.0_dp,G(ns+1,1),n,Hs,ns,0.0_dp,Gu,np)
            call dgemm('T','N',np,ns,np,1.0_dp,Hp,np,Gu,np,0.0_dp,Sol%R,np)
        end if
        allocate(Y(n,n))
        Y = 0.0_dp
        Y(ns+1:n,1:np) = Hp
        Y(1:ns,np+1:n) = Hs
        if (np > 0) call dgemm('N','N',np,ns,ns,1.0_dp,Z,np,Hs,ns,0.0_dp,Y(ns+1,np+1),n)
        call apply_q('L','N',Y)
        Sol%E = [Ep, Es]
        call move_alloc(Y,Sol%C)

    contains

        ! A <- Q^T A, Q A, A Q ... (Side, Trans of dormqr) ...
        subroutine apply_q(Side,Trans,A)
            character, intent(in) :: Side, Trans
            real(dp), intent(inout) :: A(:,:)
            real(dp), allocatable :: Wk(:)
            real(dp) :: Wq1(1)
            integer :: Inf
            call dormqr(Side,Trans,size(A,1),size(A,2),ns,Qf,n,Tau,A,size(A,1),Wq1,-1,Inf)
            allocate(Wk(int(Wq1(1))))
            call dormqr(Side,Trans,size(A,1),size(A,2),ns,Qf,n,Tau,A,size(A,1),Wk,size(Wk),Inf)
            if (Inf /= 0) call fatal('dormqr failed in the EMPM spurious subspace')
        end subroutine apply_q

    end subroutine spur_gev

    ! All eigenvalues E (ascending) & eigenvectors of the symmetric A (upper triangle), eigenvectors in A on return ...
        ! divide & conquer (dsyevd_64) ...
    subroutine EMPM_syevd(A,E)
        real(dp), intent(inout) :: A(:,:)
        real(dp), allocatable, intent(out) :: E(:)
        real(dp), allocatable :: Work(:)
        integer(int64), allocatable :: Iwork(:)
        real(dp) :: Wq(1)
        integer(int64) :: n, Info, Iwq(1)
        n = size(A,1,kind=int64)
        allocate(E(n))
        if (n == 0) return
        call dsyevd_64('V','U',n,A,n,E,Wq,-1_int64,Iwq,-1_int64,Info)
        allocate(Work(int(Wq(1),int64)), Iwork(Iwq(1)))
        call dsyevd_64('V','U',n,A,n,E,Work,size(Work,kind=int64),Iwork,size(Iwork,kind=int64),Info)
        if (Info /= 0) call fatal('dsyevd failed ... info = ' // istr(int(Info)))
    end subroutine EMPM_syevd

end module nj_empm_eqm2
