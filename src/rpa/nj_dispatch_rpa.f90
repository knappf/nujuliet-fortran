! Dispatcher of the &run actions handled by chain A ... dev TDA, RPA, HF_ScRPA & the TDA/RPA externals ...
module nj_dispatch_rpa
    use nj_input, only : run_input_t
    use nj_tda_rpa_driver, only : TDA, RPA
    use nj_scrpa, only : HF_ScRPA
    use nj_ph_external
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
        case ('TDA')
            call TDA(Inp%Params)
        case ('RPA')
            call RPA(Inp%Params)
        case ('HF_ScRPA')
            call HF_ScRPA(Inp%Params)
        ! External / post-processing actions of stored TDA & RPA solutions (TDA_external.jl, RPA_external.jl) ...
        case ('TDA_transition_densities')
            call TDA_transition_densities(Inp%Params,Inp%Densities%JP,Inp%Densities%nu_list(1:Inp%Densities%N_nu), &
                                          trim(Inp%Densities%File_Name))
        case ('RPA_transition_densities')
            call RPA_transition_densities(Inp%Params,Inp%Densities%JP,Inp%Densities%nu_list(1:Inp%Densities%N_nu), &
                                          trim(Inp%Densities%File_Name))
        case ('TDA_transition_operator_densities')
            call TDA_transition_operator_densities(Inp%Params,Inp%Op_Densities%JP, &
                                                   Inp%Op_Densities%nu_list(1:Inp%Op_Densities%N_nu),trim(Inp%Op_Densities%File_Name))
        case ('RPA_transition_operator_densities')
            call RPA_transition_operator_densities(Inp%Params,Inp%Op_Densities%JP, &
                                                   Inp%Op_Densities%nu_list(1:Inp%Op_Densities%N_nu),trim(Inp%Op_Densities%File_Name))
        case ('TDA_transition_analysis')
            call TDA_transition_analysis(Inp%Params,Inp%Analysis%JP,Inp%Analysis%nu_list(1:Inp%Analysis%N_nu), &
                                         trim(Inp%Analysis%File_Name),Inp%Analysis%Cut)
        case ('RPA_transition_analysis')
            call RPA_transition_analysis(Inp%Params,Inp%Analysis%JP,Inp%Analysis%nu_list(1:Inp%Analysis%N_nu), &
                                         trim(Inp%Analysis%File_Name),Inp%Analysis%Cut)
        case ('TDA_transition_currents')
            call TDA_transition_currents(Inp%Params,Inp%Currents%JP,Inp%Currents%nu_list(1:Inp%Currents%N_nu), &
                                         trim(Inp%Currents%File_Name),Inp%Currents%M,Inp%Currents%xN_grid,Inp%Currents%zN_grid)
        case ('RPA_transition_currents')
            call RPA_transition_currents(Inp%Params,Inp%Currents%JP,Inp%Currents%nu_list(1:Inp%Currents%N_nu), &
                                         trim(Inp%Currents%File_Name),Inp%Currents%M,Inp%Currents%xN_grid,Inp%Currents%zN_grid)
        case ('TDA_transitions_recalculate')
            call TDA_transitions_recalculate(Inp%Params)
        case ('RPA_transitions_recalculate')
            call RPA_transitions_recalculate(Inp%Params)
        case default
            Handled = .false.
        end select
    end subroutine dispatch_rpa
end module nj_dispatch_rpa
