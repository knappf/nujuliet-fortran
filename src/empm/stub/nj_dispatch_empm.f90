! Dispatcher of the &run actions handled by chain C (EMPM: EMPM_TDA, EMPM_2PH .. EMPM_SF, EMPM) ... all read the dev
    ! reference io/<Path>/<&empm Reference>/bin/* & write io/<Path>/EMPM/ (stage 5) ... stub of a build without chain C ...
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
