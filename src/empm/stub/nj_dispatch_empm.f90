! Dispatcher of the &run actions handled by chain C (EMPM) ...
module nj_dispatch_empm
    use nj_input, only : run_input_t
    implicit none
    private
    public :: dispatch_empm
contains
    subroutine dispatch_empm(Action,Inp,Handled)
        character(len=*), intent(in) :: Action
        type(run_input_t), intent(in) :: Inp
        logical, intent(out) :: Handled
        Handled = .true.
        select case (Action)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_empm
end module nj_dispatch_empm
