! Dispatcher of the &run actions handled by chain C (EMPM) ...
module nj_dispatch_empm
    use nj_params, only : params_t
    use nj_input, only : run_input_t
    use nj_empm_driver, only : EMPM_TDA, EMPM_2PH, EMPM_3PH, EMPM_4PH, EMPM_FULL, EMPM_TRANS, EMPM_SF
    implicit none
    private
    public :: dispatch_empm
contains
    subroutine dispatch_empm(Action,Inp,Handled)
        character(len=*), intent(in) :: Action
        type(run_input_t), intent(in) :: Inp
        logical, intent(out) :: Handled
        type(params_t) :: P
        Handled = .true.
        select case (Action)
        case ('EMPM_TDA')
            call EMPM_TDA(Inp%Params)
        case ('EMPM_2PH')
            call EMPM_2PH(Inp%Params)
        case ('EMPM_3PH')
            call EMPM_3PH(Inp%Params)
        case ('EMPM_4PH')
            call EMPM_4PH(Inp%Params)
        case ('EMPM_FULL')
            call EMPM_FULL(Inp%Params)
        case ('EMPM_TRANS')
            call EMPM_TRANS(Inp%Params)
        case ('EMPM_SF')
            call EMPM_SF(Inp%Params)
        case ('EMPM')
            ! Complete EMPM ... TDA phonons & F-matrix, 2-phonon states, couplings & the full Hamiltonian ...
            call EMPM_TDA(Inp%Params)
            P = Inp%Params
            P%Calc%EMPM%Full_H = .true.
            call EMPM_2PH(P)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_empm
end module nj_dispatch_empm
