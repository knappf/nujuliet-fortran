! NuJuliet ... Fortran port of the Julia many-body nuclear structure solver ...
    ! Usage: nujuliet input1.nml [input2.nml ...]   (run from the directory containing io/)
program nujuliet
    use nj_kinds
    use nj_utils
    use nj_input
    use nj_hf_driver, only : HF
    use nj_dispatch_rpa, only : dispatch_rpa
    use nj_dispatch_qp, only : dispatch_qp
    use nj_dispatch_empm, only : dispatch_empm
    use nj_ncsm, only : NCSM
    use nj_v2b, only : V2b_parse
    use nj_v3b, only : V3b_no2b_parse
    implicit none
    type(run_input_t) :: Run
    type(run_input_t), allocatable :: Runs(:)
    character(len=1024) :: Arg
    integer :: i_arg, k
    logical :: Handled

    if (command_argument_count() < 1) then
        write(*,'(a)') 'Usage: nujuliet input1.nml [input2.nml ...]   (run from the directory containing io/; the input ' // &
                       'files run in sequence, outputs in io/<Path>/)'
        stop
    end if

    ! all input files are read & checked first: actions removed upstream are refused before anything runs (no partial
        ! run of an earlier file) ...
    allocate(Runs(command_argument_count()))
    do i_arg = 1, command_argument_count()
        call get_command_argument(i_arg, Arg)
        call read_input(trim(Arg), Runs(i_arg))
        do k = 1, Runs(i_arg)%N_Actions
            call removed_action_check(trim(Runs(i_arg)%Actions(k)))
        end do
    end do

    do i_arg = 1, command_argument_count()
        call get_command_argument(i_arg, Arg)
        Run = Runs(i_arg)
        write(*,'(/,a)') '=== NuJuliet input file: ' // trim(Arg) // ' ... output folder io/' // trim(Run%Params%Calc%Path)
        do k = 1, Run%N_Actions
            select case (trim(Run%Actions(k)))
            case ('HF')
                call HF(Run%Params)
            ! NCSM (Fortran-only) ... Hamiltonian 'EMPM' / Basis 'HF' read the dev reference of &empm Reference ...
            case ('NCSM')
                call NCSM(Run%Params)
            ! Interaction parse utilities (&parse) ... 'PARSE' = V2b_parse & V3b_no2b_parse (the latter only with an NNN_File) ...
            case ('V2b_parse')
                call parse_nn()
            case ('V3b_no2b_parse')
                call parse_nnn()
            case ('PARSE')
                call parse_nn()
                if (len_trim(Run%Params%Int%NNN_File) > 0) call parse_nnn()
            ! chain A (TDA, RPA, HF_ScRPA, TDA_/RPA_transition_*), chain B (BCS, HFB, QTDA, QRPA,
                ! pnQTDA), chain C (EMPM_TDA .. EMPM_SF, EMPM) ...
            case default
                call dispatch_rpa(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call dispatch_qp(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call dispatch_empm(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call fatal('Unknown solver in &run: ' // trim(Run%Actions(k)))
            end select
        end do
    end do

contains

    ! Actions of the base version that NuJuliet dev removed (and the Fortran-only HF_LEGACY of the port) ... fatal with the
        ! replacement ...
    subroutine removed_action_check(Action)
        character(len=*), intent(in) :: Action
        select case (Action)
        case ('HF_RPA')
            call fatal("'HF_RPA' was removed upstream, use 'TDA', 'RPA' (one &tda / &rpa per input file, see README)")
        case ('HF_RPA_transition_densities', 'HF_RPA_transition_currents')
            call fatal("'" // Action // "' was removed upstream, use 'TDA_" // Action(8:) // "' or 'RPA_" // Action(8:) // "'")
        case ('HF_RRPA')
            call fatal("'HF_RRPA' was removed upstream, use 'HF_ScRPA' (&rpa)")
        case ('HF_RRPA_transition_densities', 'HF_RRPA_transition_currents')
            call fatal("'" // Action // "' was removed upstream, use 'HF_ScRPA' and then 'RPA_" // Action(9:) // &
                       "' (&rpa Name of the HF_ScRPA run)")
        case ('HF_LEGACY')
            call fatal("'HF_LEGACY' (base-version HF writing io/<Path>/Bin) was removed, use 'HF' (io/<Path>/HF/bin)")
        end select
    end subroutine removed_action_check

    subroutine parse_nn()
        if (min(Run%Parse%N_max_Out,Run%Parse%N_2max_Out) < 0) &
            call fatal('V2b_parse needs N_max_Out & N_2max_Out in &parse ...')
        call V2b_parse(Run%Params,Run%Parse%N_max_Out,Run%Parse%N_2max_Out,trim(Run%Parse%NN_Out_File))
    end subroutine parse_nn

    subroutine parse_nnn()
        if (min(Run%Parse%N_max_Out,Run%Parse%N_2max_Out,Run%Parse%N_3max_Out) < 0) &
            call fatal('V3b_no2b_parse needs N_max_Out, N_2max_Out & N_3max_Out in &parse ...')
        call V3b_no2b_parse(Run%Params,Run%Parse%N_max_Out,Run%Parse%N_2max_Out,Run%Parse%N_3max_Out,trim(Run%Parse%NNN_Out_File))
    end subroutine parse_nnn

end program nujuliet
