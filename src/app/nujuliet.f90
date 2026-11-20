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
    character(len=1024) :: Arg
    integer :: i_arg, k
    logical :: Handled

    if (command_argument_count() < 1) then
        write(*,'(a)') 'Usage: nujuliet input1.nml [input2.nml ...]'
        stop
    end if

    do i_arg = 1, command_argument_count()
        call get_command_argument(i_arg, Arg)
        call read_input(trim(Arg), Run)
        write(*,'(/,a)') '=== NuJuliet input file: ' // trim(Arg) // ' ... output folder io/' // trim(Run%Params%Calc%Path)
        do k = 1, Run%N_Actions
            select case (trim(Run%Actions(k)))
            case ('HF')
                call HF(Run%Params)
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
            case default
                call dispatch_rpa(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call dispatch_qp(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call dispatch_empm(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call fatal('Unknown solver in &run: ' // trim(Run%Actions(k)))
            end select
        end do
    end do

contains

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
