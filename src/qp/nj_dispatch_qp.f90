! Dispatcher of the &run actions handled by chain B (BCS, HFB, QTDA & QRPA) ...
module nj_dispatch_qp
    use nj_input, only : run_input_t
    use nj_bcs, only : BCS
    use nj_hfb, only : HFB
    use nj_qtda, only : QTDA
    use nj_qrpa, only : QRPA
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
        case ('BCS')
            call BCS(Inp%Params)
        case ('HFB')
            call HFB(Inp%Params)
        case ('QTDA')
            call QTDA(Inp%Params)
        case ('QRPA')
            call QRPA(Inp%Params)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_qp
end module nj_dispatch_qp
