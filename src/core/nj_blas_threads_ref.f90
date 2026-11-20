! Thread-local BLAS thread count (stage 4, scrpa.md B7) ... reference-BLAS variant (COMPILER=gnu): single-threaded BLAS,
    ! nothing to control (see nj_blas_threads_mkl.f90) ...
module nj_blas_threads
    implicit none
    private

    public :: blas_threads_local

contains

    integer function blas_threads_local(N) result(N_old)
        integer, intent(in) :: N
        N_old = 0
        if (N < 0) N_old = 0
    end function blas_threads_local

end module nj_blas_threads
