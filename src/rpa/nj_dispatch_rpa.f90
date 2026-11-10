! Dispatcher of the &run actions handled by chain A (HF-RPA & HF-RRPA) ...
module nj_dispatch_rpa
    use nj_input, only : run_input_t
    use nj_hf_rpa_driver, only : HF_RPA
    use nj_hf_rrpa_driver, only : HF_RRPA
    use nj_hf_rpa_external, only : phonon_transition_densities, phonon_transition_currents
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
        case ('HF_RPA')
            call HF_RPA(Inp%Params)
        case ('HF_RPA_transition_densities')
            call phonon_transition_densities(Inp%Params,'RPA',Inp%Densities%JP,Inp%Densities%nu_list(1:Inp%Densities%N_nu), &
                                             trim(Inp%Densities%File_Name))
        case ('HF_RPA_transition_currents')
            call phonon_transition_currents(Inp%Params,'RPA',Inp%Currents%JP,Inp%Currents%nu_list(1:Inp%Currents%N_nu), &
                                            trim(Inp%Currents%File_Name),Inp%Currents%M,Inp%Currents%xN_grid,Inp%Currents%zN_grid)
        case ('HF_RRPA')
            call HF_RRPA(Inp%Params)
        case ('HF_RRPA_transition_densities')
            call phonon_transition_densities(Inp%Params,'RRPA',Inp%Densities%JP,Inp%Densities%nu_list(1:Inp%Densities%N_nu), &
                                             trim(Inp%Densities%File_Name))
        case ('HF_RRPA_transition_currents')
            call phonon_transition_currents(Inp%Params,'RRPA',Inp%Currents%JP,Inp%Currents%nu_list(1:Inp%Currents%N_nu), &
                                            trim(Inp%Currents%File_Name),Inp%Currents%M,Inp%Currents%xN_grid,Inp%Currents%zN_grid)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_rpa
end module nj_dispatch_rpa
