! 64-bit integer LAPACK/BLAS entry points (dsygvd_64, dsyevd_64, dgemm_64, dgemv_64) for COMPILER=gnu ... the LP64 MKL
    ! provides them (used by nj_empm_eqm2 for n^2 workspaces beyond 32-bit integers), the reference LAPACK/BLAS does not:
    ! these wrappers call the default-integer routines & stop if an argument does not fit (sizes beyond the reach of the
    ! reference LAPACK anyway) ... not compiled for COMPILER=intel ...

module nj_lapack64_ref
    implicit none
    private

    public :: nj_lapack64_check

contains

    subroutine nj_lapack64_check(Name,Values)
        use nj_kinds, only : i8
        use nj_utils, only : fatal
        implicit none
        character(len=*), intent(in) :: Name
        integer(i8), intent(in) :: Values(:)
        if (any(abs(Values) > int(huge(1),i8))) call fatal(Name // ': argument beyond the default integer range of the' // &
            ' reference LAPACK/BLAS (COMPILER=gnu) ... use the MKL build (COMPILER=intel)')
    end subroutine nj_lapack64_check

end module nj_lapack64_ref

subroutine dsygvd_64(Itype,Jobz,Uplo,n,A,Lda,B,Ldb,W,Work,Lwork,Iwork,Liwork,Info)
    use nj_kinds, only : dp, i8
    use nj_lapack64_ref, only : nj_lapack64_check
    implicit none
    character :: Jobz, Uplo
    integer(i8) :: Itype, n, Lda, Ldb, Lwork, Liwork, Info, Iwork(*)
    real(dp) :: A(Lda,*), B(Ldb,*), W(*), Work(*)
    integer, allocatable :: Iw(:)
    integer :: Info4
    call nj_lapack64_check('dsygvd_64',[Itype,n,Lda,Ldb,Lwork,Liwork])
    allocate(Iw(max(1_i8,Liwork)))
    call dsygvd(int(Itype),Jobz,Uplo,int(n),A,int(Lda),B,int(Ldb),W,Work,int(Lwork),Iw,int(Liwork),Info4)
    Iwork(1:max(1_i8,Liwork)) = Iw
    Info = Info4
end subroutine dsygvd_64

subroutine dsyevd_64(Jobz,Uplo,n,A,Lda,W,Work,Lwork,Iwork,Liwork,Info)
    use nj_kinds, only : dp, i8
    use nj_lapack64_ref, only : nj_lapack64_check
    implicit none
    character :: Jobz, Uplo
    integer(i8) :: n, Lda, Lwork, Liwork, Info, Iwork(*)
    real(dp) :: A(Lda,*), W(*), Work(*)
    integer, allocatable :: Iw(:)
    integer :: Info4
    call nj_lapack64_check('dsyevd_64',[n,Lda,Lwork,Liwork])
    allocate(Iw(max(1_i8,Liwork)))
    call dsyevd(Jobz,Uplo,int(n),A,int(Lda),W,Work,int(Lwork),Iw,int(Liwork),Info4)
    Iwork(1:max(1_i8,Liwork)) = Iw
    Info = Info4
end subroutine dsyevd_64

subroutine dgemm_64(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
    use nj_kinds, only : dp, i8
    use nj_lapack64_ref, only : nj_lapack64_check
    implicit none
    character :: Transa, Transb
    integer(i8) :: m, n, k, Lda, Ldb, Ldc
    real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
    call nj_lapack64_check('dgemm_64',[m,n,k,Lda,Ldb,Ldc])
    call dgemm(Transa,Transb,int(m),int(n),int(k),Alpha,A,int(Lda),B,int(Ldb),Beta,C,int(Ldc))
end subroutine dgemm_64

subroutine dgemv_64(Trans,m,n,Alpha,A,Lda,x,Incx,Beta,y,Incy)
    use nj_kinds, only : dp, i8
    use nj_lapack64_ref, only : nj_lapack64_check
    implicit none
    character :: Trans
    integer(i8) :: m, n, Lda, Incx, Incy
    real(dp) :: Alpha, Beta, A(Lda,*), x(*), y(*)
    call nj_lapack64_check('dgemv_64',[m,n,Lda,Incx,Incy])
    call dgemv(Trans,int(m),int(n),Alpha,A,int(Lda),x,int(Incx),Beta,y,int(Incy))
end subroutine dgemv_64
