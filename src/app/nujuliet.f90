! NuJuliet ... Fortran port of the Julia many-body nuclear structure solver ...
    ! Usage: nujuliet input1.nml [input2.nml ...]   (run from the directory containing io/)
program nujuliet
    use nj_kinds
    use nj_utils
    use nj_input
    use nj_hf_driver, only : HF
    use nj_hf_legacy_driver, only : HF_LEGACY      ! LEGACY (removed in cleanup stage)
    use nj_dispatch_rpa, only : dispatch_rpa
    use nj_dispatch_qp, only : dispatch_qp
    use nj_dispatch_empm, only : dispatch_empm
    use nj_ncsm, only : NCSM
    use nj_v2b, only : V2b_parse
    use nj_v3b, only : V3b_no2b_parse
    implicit none
    type(run_input_t) :: Run
    character(len=1024) :: Arg
    integer :: i_arg, k
    logical :: Handled, Legacy_HF

    if (command_argument_count() < 1) then
        write(*,'(a)') 'Usage: nujuliet input1.nml [input2.nml ...]   (run from the directory containing io/; the input ' // &
                       'files run in sequence, outputs in io/<Path>/)'
        stop
    end if

    do i_arg = 1, command_argument_count()
        call get_command_argument(i_arg, Arg)
        call read_input(trim(Arg), Run)
        write(*,'(/,a)') '=== NuJuliet input file: ' // trim(Arg) // ' ... output folder io/' // trim(Run%Params%Calc%Path)
        ! LEGACY (removed in cleanup stage) ... HF_RPA* & HF_RRPA* read the base-version HF output io/<Path>/Bin/*, so 'HF'
        ! runs the legacy HF when one of them is in the same &run (EMPM_*, EMPM & NCSM read the dev reference) ...
        Legacy_HF = .false.
        do k = 1, Run%N_Actions
            if (legacy_hf_consumer(trim(Run%Actions(k)))) Legacy_HF = .true.
        end do
        do k = 1, Run%N_Actions
            select case (trim(Run%Actions(k)))
            case ('HF')
                if (Legacy_HF) then
                    write(*,'(/,a)') "NOTICE: &run contains a legacy solver that reads io/<Path>/Bin/* ... 'HF' runs HF_LEGACY ..."
                    call HF_LEGACY(Run%Params)
                else
                    call HF(Run%Params)
                end if
            case ('HF_LEGACY')       ! LEGACY (removed in cleanup stage)
                call HF_LEGACY(Run%Params)
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
            ! chain A (TDA, RPA, HF_ScRPA, TDA_/RPA_transition_*, LEGACY HF_RPA* / HF_RRPA*), chain B (BCS, HFB, QTDA, QRPA,
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

    ! LEGACY (removed in cleanup stage) ... solvers that read the base-version HF output (io/<Path>/Bin/*) ... the EMPM & the
        ! NCSM read the dev reference io/<Path>/<&empm Reference>/bin/* since stage 5 (no longer routed) ...
    logical function legacy_hf_consumer(Action)
        character(len=*), intent(in) :: Action
        legacy_hf_consumer = (index(Action,'HF_RPA') == 1 .or. index(Action,'HF_RRPA') == 1)
    end function legacy_hf_consumer

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
