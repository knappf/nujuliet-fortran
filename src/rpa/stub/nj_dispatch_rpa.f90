! Dispatcher of the &run actions handled by chain A (HF-RPA & HF-RRPA) ...
module nj_dispatch_rpa
    use nj_input, only : run_input_t
    implicit none
    private
    public :: dispatch_rpa
contains
    subroutine dispatch_rpa(Action,Inp,Handled)
        character(len=*), intent(in) :: Action
        type(run_input_t), intent(in) :: Inp
        logical, intent(out) :: Handled
        Handled = .true.
        select case (Action)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_rpa
end module nj_dispatch_rpa
