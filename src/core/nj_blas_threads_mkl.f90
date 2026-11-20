! Thread-local BLAS thread count (stage 4, scrpa.md B7) ... MKL variant (COMPILER=intel): mkl_set_num_threads_local sets the
    ! number of MKL threads of the calling thread only (0 = back to the global setting), the thread-safe counterpart of the
    ! global BLAS.set_num_threads(1) of dev (C entry MKL_Set_Num_Threads_Local, by value; the lowercase symbol is the
    ! by-reference Fortran entry)'s RPA_OBDM_oY4 / RPA_OBDM_transform_amplitudes ... the reference-BLAS variant
    ! (nj_blas_threads_ref.f90, COMPILER=gnu) is a no-op ...
module nj_blas_threads
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    private

    public :: blas_threads_local

    interface
        integer(c_int) function mkl_set_num_threads_local(nt) bind(C, name='MKL_Set_Num_Threads_Local')
            import :: c_int
            integer(c_int), value :: nt
        end function mkl_set_num_threads_local
    end interface

contains

    ! Set the BLAS threads of the calling thread to N, return the previous thread-local value (pass it back to restore) ...
    integer function blas_threads_local(N) result(N_old)
        integer, intent(in) :: N
        N_old = int(mkl_set_num_threads_local(int(N,c_int)))
    end function blas_threads_local

end module nj_blas_threads
