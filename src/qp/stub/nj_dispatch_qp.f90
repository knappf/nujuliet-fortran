! Dispatcher of the &run actions handled by chain B (BCS, HFB, QTDA, QRPA & pnQTDA) ... stub of a build without chain B ...
module nj_dispatch_qp
    use nj_input, only : run_input_t
    implicit none
    private
    public :: dispatch_qp
contains
    subroutine dispatch_qp(Action,Inp,Handled)
        character(len=*), intent(in) :: Action
        type(run_input_t), intent(in) :: Inp
        logical, intent(out) :: Handled
        Handled = .true.
        select case (Action)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_qp
end module nj_dispatch_qp
