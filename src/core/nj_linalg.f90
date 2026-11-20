! LAPACK/BLAS wrappers reproducing the behaviour of Julia's LinearAlgebra calls used in NuJuliet ...
module nj_linalg
    use nj_kinds
    use nj_utils, only : sortperm, sortperm_rev, fatal
    implicit none
    private

    public :: eig_sym, eig_general, eig_general_complex, isposdef, mm, mtm, mmt, mtmm, identity, diagm, diag
    public :: zmm, zmhm, isposdef_c, hungarian
    public :: eig_sym_gen        ! stage 3 (TDA): generalized symmetric-definite problem ...

    interface
        subroutine dsygvd(Itype,Jobz,Uplo,n,A,Lda,B,Ldb,w,Work,Lwork,Iwork,Liwork,Info)
            import :: dp
            character :: Jobz, Uplo
            integer :: Itype, n, Lda, Ldb, Lwork, Liwork, Info
            real(dp) :: A(Lda,*), B(Ldb,*), w(*), Work(*)
            integer :: Iwork(*)
        end subroutine dsygvd
    end interface

    interface
        subroutine dsyevr(Jobz,Range,Uplo,n,A,Lda,Vl,Vu,Il,Iu,Abstol,m,w,z,Ldz,Isuppz,Work,Lwork,Iwork,Liwork,Info)
            import :: dp
            character :: Jobz, Range, Uplo
            integer :: n, Lda, Il, Iu, m, Ldz, Lwork, Liwork, Info
            integer :: Isuppz(*), Iwork(*)
            real(dp) :: A(Lda,*), Vl, Vu, Abstol, w(*), z(Ldz,*), Work(*)
        end subroutine dsyevr
        subroutine dgeevx(Balanc,Jobvl,Jobvr,Sense,n,A,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv, &
                          Work,Lwork,Iwork,Info)
            import :: dp
            character :: Balanc, Jobvl, Jobvr, Sense
            integer :: n, Lda, Ldvl, Ldvr, Ilo, Ihi, Lwork, Info
            integer :: Iwork(*)
            real(dp) :: A(Lda,*), Wr(*), Wi(*), Vl(Ldvl,*), Vr(Ldvr,*), Scale(*), Abnrm, Rconde(*), Rcondv(*), Work(*)
        end subroutine dgeevx
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
        subroutine zgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            complex(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine zgemm
        subroutine dpotrf(Uplo,n,A,Lda,Info)
            import :: dp
            character :: Uplo
            integer :: n, Lda, Info
            real(dp) :: A(Lda,*)
        end subroutine dpotrf
        subroutine zpotrf(Uplo,n,A,Lda,Info)
            import :: dp
            character :: Uplo
            integer :: n, Lda, Info
            complex(dp) :: A(Lda,*)
        end subroutine zpotrf
        subroutine zheevr(Jobz,Range,Uplo,n,A,Lda,Vl,Vu,Il,Iu,Abstol,m,w,z,Ldz,Isuppz,Work,Lwork,Rwork,Lrwork,Iwork,Liwork,Info)
            import :: dp
            character :: Jobz, Range, Uplo
            integer :: n, Lda, Il, Iu, m, Ldz, Lwork, Lrwork, Liwork, Info
            integer :: Isuppz(*), Iwork(*)
            real(dp) :: Vl, Vu, Abstol, w(*), Rwork(*)
            complex(dp) :: A(Lda,*), z(Ldz,*), Work(*)
        end subroutine zheevr
        subroutine zgeevx(Balanc,Jobvl,Jobvr,Sense,n,A,Lda,W,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv, &
                          Work,Lwork,Rwork,Info)
            import :: dp
            character :: Balanc, Jobvl, Jobvr, Sense
            integer :: n, Lda, Ldvl, Ldvr, Ilo, Ihi, Lwork, Info
            real(dp) :: Scale(*), Abnrm, Rconde(*), Rcondv(*), Rwork(*)
            complex(dp) :: A(Lda,*), W(*), Vl(Ldvl,*), Vr(Ldvr,*), Work(*)
        end subroutine zgeevx
    end interface

contains

    ! Julia eigen(Symmetric(A), sortby = +/-) ... syevr with the upper triangle ...
        ! Descending = .true. corresponds to sortby = - ... stable sorting of LAPACK's ascending output ...
    subroutine eig_sym(A,w,Z,Descending)
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(out) :: w(:)
        real(dp), intent(out) :: Z(:,:)
        logical, intent(in), optional :: Descending
        real(dp), allocatable :: B(:,:), Work(:), wt(:), Zt(:,:)
        integer, allocatable :: Iwork(:), Isuppz(:), Perm(:)
        real(dp) :: Work_q(1)
        integer :: n, m, Info, Lwork, Liwork, Iwork_q(1)

        n = size(A,1)
        if (n == 0) return
        allocate(B(n,n), Isuppz(2*n))
        B = A
        call dsyevr('V','A','U',n,B,n,0.0_dp,0.0_dp,0,0,-1.0_dp,m,w,Z,n,Isuppz,Work_q,-1,Iwork_q,-1,Info)
        Lwork = int(Work_q(1))
        Liwork = Iwork_q(1)
        allocate(Work(Lwork), Iwork(Liwork))
        call dsyevr('V','A','U',n,B,n,0.0_dp,0.0_dp,0,0,-1.0_dp,m,w,Z,n,Isuppz,Work,Lwork,Iwork,Liwork,Info)
        if (Info /= 0) call fatal('dsyevr failed')

        if (present(Descending)) then
            if (Descending) then
                allocate(Perm(n), wt(n), Zt(n,n))
                call sortperm_rev(w(1:n),Perm)
                wt = w(Perm)
                Zt = Z(:,Perm)
                w(1:n) = wt
                Z(:,1:n) = Zt
            end if
        end if
    end subroutine eig_sym

    ! Julia eigen!(Symmetric(A),Symmetric(B)) ... LAPACK sygvd(1,'V','U'), eigenvalues ascending (stable sort as Julia's
        ! sorteig!, a no-op on sygvd's ascending output), eigenvectors normalized to Z' B Z = 1 ... (stage 3, TDA_diagonalize)
    subroutine eig_sym_gen(A,B,w,Z)
        real(dp), intent(in) :: A(:,:), B(:,:)
        real(dp), intent(out) :: w(:)
        real(dp), intent(out) :: Z(:,:)
        real(dp), allocatable :: Bw(:,:), Work(:)
        integer, allocatable :: Iwork(:)
        real(dp) :: Work_q(1)
        integer :: n, Info, Lwork, Liwork, Iwork_q(1)

        n = size(A,1)
        if (n == 0) return
        Z(1:n,1:n) = A
        Bw = B
        call dsygvd(1,'V','U',n,Z,size(Z,1),Bw,n,w,Work_q,-1,Iwork_q,-1,Info)
        Lwork = int(Work_q(1))
        Liwork = Iwork_q(1)
        allocate(Work(Lwork), Iwork(Liwork))
        call dsygvd(1,'V','U',n,Z,size(Z,1),Bw,n,w,Work,Lwork,Iwork,Liwork,Info)
        if (Info /= 0) call fatal('dsygvd failed ... Info = ' // achar(48 + min(9,max(0,Info))))
    end subroutine eig_sym_gen

    ! Julia eigen(A::Matrix{Float64}) ... symmetric matrices dispatch to syevr, otherwise geevx('B','N','V','N') ...
        ! Eigenpairs are sorted by (real, imag) as with Julia's default eigsortby ...
        ! Is_Real = .true. when all eigenvalues are real ... W & V are then purely real ...
    subroutine eig_general(A,W,V,Is_Real)
        real(dp), intent(in) :: A(:,:)
        complex(dp), intent(out) :: W(:)
        complex(dp), intent(out) :: V(:,:)
        logical, intent(out) :: Is_Real
        real(dp), allocatable :: B(:,:), Wr(:), Wi(:), Vr(:,:), Vl(:,:), Scale(:), Rconde(:), Rcondv(:), Work(:), Zs(:,:)
        integer, allocatable :: Iwork(:), Perm(:)
        complex(dp), allocatable :: Wt(:), Vt(:,:)
        real(dp) :: Abnrm, Work_q(1)
        integer :: n, Ilo, Ihi, Info, Lwork, j

        n = size(A,1)
        if (n == 0) then
            Is_Real = .true.
            return
        end if

        if (all(A == transpose(A))) then
            allocate(Wr(n), Zs(n,n))
            call eig_sym(A,Wr,Zs)
            W = cmplx(Wr, 0.0_dp, kind=dp)
            V = cmplx(Zs, 0.0_dp, kind=dp)
            Is_Real = .true.
            return
        end if

        allocate(B(n,n), Wr(n), Wi(n), Vr(n,n), Vl(1,1), Scale(n), Rconde(n), Rcondv(n), Iwork(1))
        B = A
        call dgeevx('B','N','V','N',n,B,n,Wr,Wi,Vl,1,Vr,n,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work_q,-1,Iwork,Info)
        Lwork = int(Work_q(1))
        allocate(Work(Lwork))
        call dgeevx('B','N','V','N',n,B,n,Wr,Wi,Vl,1,Vr,n,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work,Lwork,Iwork,Info)
        if (Info /= 0) call fatal('dgeevx failed')

        Is_Real = all(Wi == 0.0_dp)
        if (Is_Real) then
            W = cmplx(Wr, 0.0_dp, kind=dp)
            V = cmplx(Vr, 0.0_dp, kind=dp)
        else
            W = cmplx(Wr, Wi, kind=dp)
            j = 1
            do while (j <= n)
                if (Wi(j) == 0.0_dp) then
                    V(:,j) = cmplx(Vr(:,j), 0.0_dp, kind=dp)
                else
                    V(:,j) = cmplx(Vr(:,j), Vr(:,j+1), kind=dp)
                    V(:,j+1) = cmplx(Vr(:,j), -Vr(:,j+1), kind=dp)
                    j = j + 1
                end if
                j = j + 1
            end do
        end if

        ! Stable sort by (real, imag) ...
        allocate(Perm(n), Wt(n), Vt(n,n))
        call sortperm_complex(W,Perm)
        Wt = W(Perm)
        Vt = V(:,Perm)
        W = Wt
        V = Vt
    end subroutine eig_general

    ! Julia isposdef(A::Matrix{Float64}) ... exact symmetry (ishermitian) & a successful Cholesky factorization ...
    logical function isposdef(A)
        real(dp), intent(in) :: A(:,:)
        real(dp), allocatable :: B(:,:)
        integer :: n, Info
        n = size(A,1)
        isposdef = .false.
        if (size(A,2) /= n) return
        if (.not. all(A == transpose(A))) return
        if (n == 0) then
            isposdef = .true.
            return
        end if
        B = A
        call dpotrf('U',n,B,n,Info)
        isposdef = (Info == 0)
    end function isposdef

    ! Julia isposdef(A::Matrix{ComplexF64}) ... exact Hermiticity & a successful Cholesky factorization ...
    logical function isposdef_c(A)
        complex(dp), intent(in) :: A(:,:)
        complex(dp), allocatable :: B(:,:)
        integer :: n, Info
        n = size(A,1)
        isposdef_c = .false.
        if (size(A,2) /= n) return
        if (.not. all(A == conjg(transpose(A)))) return
        if (n == 0) then
            isposdef_c = .true.
            return
        end if
        B = A
        call zpotrf('U',n,B,n,Info)
        isposdef_c = (Info == 0)
    end function isposdef_c

    ! Optimal assignment minimizing sum(Cost(i,Assign(i))) ... Hungarian algorithm (potentials, O(n^3)) ...
        ! Counterpart of hungarian() of Hungarian.jl ... Assign(i) is the column assigned to row i ...
    subroutine hungarian(Cost,Assign)
        real(dp), intent(in) :: Cost(:,:)
        integer, intent(out) :: Assign(:)
        real(dp), allocatable :: u(:), v(:), Minv(:)
        integer, allocatable :: Way(:), Pmatch(:)
        logical, allocatable :: Used(:)
        integer :: n, i, j, i0, j0, j1
        real(dp) :: Delta, Cur

        n = size(Cost,1)
        ! dev's in-house hungarian (functions.jl:711) also handles rectangular & non-finite costs, all callers pass square
            ! finite matrices ... the port supports the square case only ...
        if (size(Cost,2) /= n .or. size(Assign) < n) call fatal('hungarian: square cost matrix expected ...')
        allocate(u(0:n), v(0:n), Minv(0:n), Way(0:n), Pmatch(0:n), Used(0:n))
        u = 0.0_dp
        v = 0.0_dp
        Pmatch = 0
        Way = 0
        do i = 1, n
            Pmatch(0) = i
            j0 = 0
            Minv = huge(1.0_dp)
            Used = .false.
            do
                Used(j0) = .true.
                i0 = Pmatch(j0)
                Delta = huge(1.0_dp)
                j1 = 0
                do j = 1, n
                    if (.not. Used(j)) then
                        Cur = Cost(i0,j) - u(i0) - v(j)
                        if (Cur < Minv(j)) then
                            Minv(j) = Cur
                            Way(j) = j0
                        end if
                        if (Minv(j) < Delta) then
                            Delta = Minv(j)
                            j1 = j
                        end if
                    end if
                end do
                do j = 0, n
                    if (Used(j)) then
                        u(Pmatch(j)) = u(Pmatch(j)) + Delta
                        v(j) = v(j) - Delta
                    else
                        Minv(j) = Minv(j) - Delta
                    end if
                end do
                j0 = j1
                if (Pmatch(j0) == 0) exit
            end do
            do
                j1 = Way(j0)
                Pmatch(j0) = Pmatch(j1)
                j0 = j1
                if (j0 == 0) exit
            end do
        end do
        do j = 1, n
            Assign(Pmatch(j)) = j
        end do
    end subroutine hungarian

    ! Julia eigen(A::Matrix{ComplexF64}) ... Hermitian matrices dispatch to zheevr, otherwise zgeevx('B','N','V','N') ...
        ! Eigenpairs are sorted by (real, imag) ...
    subroutine eig_general_complex(A,W,V)
        complex(dp), intent(in) :: A(:,:)
        complex(dp), intent(out) :: W(:)
        complex(dp), intent(out) :: V(:,:)
        complex(dp), allocatable :: B(:,:), Work(:), Vl(:,:), Wt(:), Vt(:,:)
        real(dp), allocatable :: w_r(:), Rwork(:), Scale(:), Rconde(:), Rcondv(:)
        integer, allocatable :: Iwork(:), Isuppz(:), Perm(:)
        complex(dp) :: Work_q(1)
        real(dp) :: Rwork_q(1), Abnrm
        integer :: n, m, Info, Lwork, Lrwork, Liwork, Iwork_q(1), Ilo, Ihi

        n = size(A,1)
        if (n == 0) return
        B = A
        if (all(A == conjg(transpose(A)))) then
            allocate(w_r(n), Isuppz(2*n))
            call zheevr('V','A','U',n,B,n,0.0_dp,0.0_dp,0,0,-1.0_dp,m,w_r,V,n,Isuppz,Work_q,-1,Rwork_q,-1,Iwork_q,-1,Info)
            Lwork = int(real(Work_q(1)))
            Lrwork = int(Rwork_q(1))
            Liwork = Iwork_q(1)
            allocate(Work(Lwork), Rwork(Lrwork), Iwork(Liwork))
            call zheevr('V','A','U',n,B,n,0.0_dp,0.0_dp,0,0,-1.0_dp,m,w_r,V,n,Isuppz,Work,Lwork,Rwork,Lrwork,Iwork,Liwork,Info)
            if (Info /= 0) call fatal('zheevr failed')
            W = cmplx(w_r, 0.0_dp, kind=dp)
            return
        end if

        allocate(Vl(1,1), Scale(n), Rconde(n), Rcondv(n), Rwork(2*n))
        call zgeevx('B','N','V','N',n,B,n,W,Vl,1,V,n,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work_q,-1,Rwork,Info)
        Lwork = int(real(Work_q(1)))
        allocate(Work(Lwork))
        call zgeevx('B','N','V','N',n,B,n,W,Vl,1,V,n,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work,Lwork,Rwork,Info)
        if (Info /= 0) call fatal('zgeevx failed')

        allocate(Perm(n), Wt(n), Vt(n,n))
        call sortperm_complex(W,Perm)
        Wt = W(Perm)
        Vt = V(:,Perm)
        W = Wt
        V = Vt
    end subroutine eig_general_complex

    ! C = A * B for complex matrices ...
    function zmm(A,B) result(C)
        complex(dp), intent(in) :: A(:,:), B(:,:)
        complex(dp), allocatable :: C(:,:)
        allocate(C(size(A,1),size(B,2)))
        if (size(C) == 0) return
        if (size(A,2) == 0) then
            C = (0.0_dp, 0.0_dp)
            return
        end if
        call zgemm('N','N',size(A,1),size(B,2),size(A,2),(1.0_dp,0.0_dp),A,size(A,1),B,size(B,1),(0.0_dp,0.0_dp),C,size(C,1))
    end function zmm

    ! C = A^H * B for complex matrices ... Julia A' * B ...
    function zmhm(A,B) result(C)
        complex(dp), intent(in) :: A(:,:), B(:,:)
        complex(dp), allocatable :: C(:,:)
        allocate(C(size(A,2),size(B,2)))
        if (size(C) == 0) return
        if (size(A,1) == 0) then
            C = (0.0_dp, 0.0_dp)
            return
        end if
        call zgemm('C','N',size(A,2),size(B,2),size(A,1),(1.0_dp,0.0_dp),A,size(A,1),B,size(B,1),(0.0_dp,0.0_dp),C,size(C,1))
    end function zmhm

    ! Stable lexicographic (real, imag) sortperm ... insertion sort, sizes are moderate ...
    subroutine sortperm_complex(W,Perm)
        complex(dp), intent(in) :: W(:)
        integer, intent(out) :: Perm(:)
        integer :: i, j, k, n
        n = size(W)
        do i = 1, n
            Perm(i) = i
        end do
        do i = 2, n
            k = Perm(i)
            j = i - 1
            do while (j >= 1)
                if (.not. lex_less(W(k),W(Perm(j)))) exit
                Perm(j+1) = Perm(j)
                j = j - 1
            end do
            Perm(j+1) = k
        end do
    contains
        pure logical function lex_less(x,y)
            complex(dp), intent(in) :: x, y
            if (real(x) /= real(y)) then
                lex_less = real(x) < real(y)
            else
                lex_less = aimag(x) < aimag(y)
            end if
        end function lex_less
    end subroutine sortperm_complex

    ! C = A * B ...
    function mm(A,B) result(C)
        real(dp), intent(in) :: A(:,:), B(:,:)
        real(dp), allocatable :: C(:,:)
        allocate(C(size(A,1),size(B,2)))
        if (size(C) == 0) return
        if (size(A,2) == 0) then
            C = 0.0_dp
            return
        end if
        call dgemm('N','N',size(A,1),size(B,2),size(A,2),1.0_dp,A,size(A,1),B,size(B,1),0.0_dp,C,size(C,1))
    end function mm

    ! C = A' * B ...
    function mtm(A,B) result(C)
        real(dp), intent(in) :: A(:,:), B(:,:)
        real(dp), allocatable :: C(:,:)
        allocate(C(size(A,2),size(B,2)))
        if (size(C) == 0) return
        if (size(A,1) == 0) then
            C = 0.0_dp
            return
        end if
        call dgemm('T','N',size(A,2),size(B,2),size(A,1),1.0_dp,A,size(A,1),B,size(B,1),0.0_dp,C,size(C,1))
    end function mtm

    ! C = A * B' ...
    function mmt(A,B) result(C)
        real(dp), intent(in) :: A(:,:), B(:,:)
        real(dp), allocatable :: C(:,:)
        allocate(C(size(A,1),size(B,1)))
        if (size(C) == 0) return
        if (size(A,2) == 0) then
            C = 0.0_dp
            return
        end if
        call dgemm('N','T',size(A,1),size(B,1),size(A,2),1.0_dp,A,size(A,1),B,size(B,1),0.0_dp,C,size(C,1))
    end function mmt

    ! C = U' * A * U ... basis transformation ...
    function mtmm(U,A) result(C)
        real(dp), intent(in) :: U(:,:), A(:,:)
        real(dp), allocatable :: C(:,:)
        C = mm(mtm(U,A),U)
    end function mtmm

    function identity(n) result(E)
        integer, intent(in) :: n
        real(dp) :: E(n,n)
        integer :: i
        E = 0.0_dp
        do i = 1, n
            E(i,i) = 1.0_dp
        end do
    end function identity

    function diagm(d) result(E)
        real(dp), intent(in) :: d(:)
        real(dp) :: E(size(d),size(d))
        integer :: i
        E = 0.0_dp
        do i = 1, size(d)
            E(i,i) = d(i)
        end do
    end function diagm

    function diag(A) result(d)
        real(dp), intent(in) :: A(:,:)
        real(dp) :: d(min(size(A,1),size(A,2)))
        integer :: i
        do i = 1, size(d)
            d(i) = A(i,i)
        end do
    end function diag

end module nj_linalg
