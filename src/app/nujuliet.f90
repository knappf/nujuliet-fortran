! NuJuliet ... Fortran port of the Julia many-body nuclear structure solver ...
    ! Usage: nujuliet input1.nml [input2.nml ...]   (run from the directory containing IO/)
program nujuliet
    use nj_kinds
    use nj_utils
    use nj_input
    use nj_hf_driver, only : HF
    use nj_dispatch_rpa, only : dispatch_rpa
    use nj_dispatch_qp, only : dispatch_qp
    use nj_dispatch_empm, only : dispatch_empm
    use nj_ncsm, only : NCSM
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
        write(*,'(/,a)') '=== NuJuliet input file: ' // trim(Arg) // ' ... output folder IO/' // trim(Run%Params%Calc%Path)
        do k = 1, Run%N_Actions
            select case (trim(Run%Actions(k)))
            case ('HF')
                call HF(Run%Params)
            case ('NCSM')
                call NCSM(Run%Params)
            case default
                call dispatch_rpa(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call dispatch_qp(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call dispatch_empm(trim(Run%Actions(k)),Run,Handled)
                if (.not. Handled) call fatal('Unknown solver in &run: ' // trim(Run%Actions(k)))
            end select
        end do
    end do
end program nujuliet
