! Kind parameters & physical constants shared by the whole NuJuliet port ...
module nj_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    integer, parameter, public :: sp = kind(1.0)
    integer, parameter, public :: qp = selected_real_kind(30)
    integer, parameter, public :: i8 = selected_int_kind(18)

    real(dp), parameter, public :: pi = 3.141592653589793238462643383279502884_dp

    ! Physical constants used throughout the original code ...
    real(dp), parameter, public :: HbarC = 197.326980_dp
    real(dp), parameter, public :: m_p = 938.272013_dp
    real(dp), parameter, public :: m_n = 939.565346_dp

    character(len=1), parameter, public :: TAB = achar(9)
end module nj_kinds
