! TDA & RPA entry points (NuJuliet dev) ... counterpart of src/TDA/TDA.jl, TDA_solver.jl (TDA_solver, TDA_solver_import),
    ! src/RPA/RPA.jl (RPA) & RPA_solver.jl (RPA_solver) ... both read the reference io/<Path>/<Reference>/bin/{orb,C,rho,H1b,
    ! H2b[,sigma]}.bin (HF or HF-MBPT) & write io/<Path>/<Name>/ ...
    ! Deviations (robustness guards, docs/PORTING_NOTES.md): TDA writes the marker <Name>/.TDA_output (dev checks for it but
    ! never writes it, so a TDA re-run with the same Name refused to overwrite); RPA stops when Reference_TBCF = .true. & the
    ! reference sigma.bin is missing (dev silently ignores it) ... Name == Reference is refused by the input reader ...
module nj_tda_rpa_driver
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_import, orbitals_1p1h_make
    use nj_o1b, only : O1b_import
    use nj_o2b, only : O2b_import
    use nj_obdm, only : rho_chR2
    use nj_tr1b
    use nj_ph_common
    use nj_tda
    use nj_rpa
    use nj_rpa_export, only : RPA_export
    use nj_rpa_se, only : RPA_allocate_SE
    !$ use omp_lib
    implicit none
    private

    public :: TDA, TDA_solver, RPA, RPA_solver
    ! Shared with HF_ScRPA (nj_scrpa, stage 4): identical header, thread report & output tree as RPA() ...
    public :: print_header, print_threads, make_output_dirs

contains

    subroutine print_header(Params,Title)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Title
        integer :: N_max
        N_max = Params%Calc%Nmax
        write(*,'(a)') Title
        write(*,'(/,a)') 'Calculation data:'
        write(*,'(a)') 'A = ' // istr(Params%Calc%A) // ' , Z = ' // istr(Params%Calc%Z) // ' , hw = ' // jstr(Params%Calc%hw) // &
                       ' MeV , N_max = ' // istr(N_max) // ' , J-scheme LHO basis size = ' // istr((N_max+1)*(N_max+2)/2) // &
                       ' , M-scheme LHO basis size = ' // istr((N_max+1)*(N_max+2)*(N_max+3)/6)
    end subroutine print_header

    subroutine print_threads()
        integer :: N_threads
        N_threads = 1
        !$ N_threads = omp_get_max_threads()
        if (N_threads > 1) then
            write(*,'(/,a)') '_________________________________________________________'
            write(*,'(a)') 'Multiple active threads detected ...'
            write(*,'(a)') 'Parallelization report:    Number of active threads = ' // istr(N_threads)
            write(*,'(a)') '_________________________________________________________'
        end if
    end subroutine print_threads

    ! Output directory tree io/<Path>/<Name>/{bin,densities,transitions/E0..E3,spectra,amplitudes} ...
    subroutine make_output_dirs(Out)
        character(len=*), intent(in) :: Out
        call make_dir(Out)
        call make_dir(Out // '/bin')
        call make_dir(Out // '/densities')
        call make_dir(Out // '/transitions/E0')
        call make_dir(Out // '/transitions/E1')
        call make_dir(Out // '/transitions/E2')
        call make_dir(Out // '/transitions/E3')
        call make_dir(Out // '/spectra')
        call make_dir(Out // '/amplitudes')
    end subroutine make_output_dirs

    logical function is_file(Path)
        character(len=*), intent(in) :: Path
        inquire(file=Path, exist=is_file)
    end function is_file

    ! Directory empty? ...
    logical function dir_is_empty(Path)
        character(len=*), intent(in) :: Path
        integer :: Stat
        call execute_command_line('test -z "$(ls -A ''' // Path // ''')"', exitstat=Stat)
        dir_is_empty = (Stat == 0)
    end function dir_is_empty

    ! TDA driver ... counterpart of TDA(Params) (TDA/TDA.jl) ...
    subroutine TDA(Params)
        type(params_t), intent(in) :: Params
        character(len=:), allocatable :: IO_ref, IO_TDA
        integer :: u
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        call print_header(Params,'Starting TDA calculations with residual NO2B NN+NNN interaction & for given reference solution ...')
        write(*,'(/,a)') "TDA parameters:    Name = ''" // trim(Params%Calc%TDA%Name) // "'' , Reference = ''" // &
                         trim(Params%Calc%TDA%Reference) // "'' , Reference_TBCF = " // lstr(Params%Calc%TDA%Reference_TBCF) // &
                         ' , Ortho = ' // lstr(Params%Calc%TDA%Ortho) // ' , Tol_hol = ' // jstr(Params%Calc%TDA%Tol_hol) // &
                         ' , Tol_par = ' // jstr(Params%Calc%TDA%Tol_par)
        call print_threads()

        IO_ref = io_dir(Params,Params%Calc%TDA%Reference)
        if (.not. is_dir(IO_ref)) then
            write(*,'(/,a)') "Error! ... No precomputed reference solutions available in IO directory ''" // io_path(Params,'') // &
                             "'' for given input path ''" // trim(Params%Calc%TDA%Reference) // "''... e.g. run HF solver first ..."
            return
        end if
        if (Params%Calc%TDA%Reference_TBCF .and. .not. is_file(IO_ref // '/bin/sigma.bin')) then
            write(*,'(/,a)') "Error! ... Reference_TBCF = true, but no reference 2-body correlation function ''" // IO_ref // &
                             "/bin/sigma.bin'' is available ..."
            return
        end if

        IO_TDA = io_dir(Params,Params%Calc%TDA%Name)
        if (trim(Params%Calc%TDA%Name) == trim(Params%Calc%TDA%Reference)) then
            write(*,'(/,a)') "Error! ... TDA output directory ''" // IO_TDA // "'' coincides with the reference directory ... " // &
                             'choose a different TDA_Parameters.Name ...'
            return
        end if
        if (is_dir(IO_TDA)) then
            if (is_file(IO_TDA // '/.TDA_output') .or. dir_is_empty(IO_TDA)) then
                call remove_dir(IO_TDA)
            else
                write(*,'(/,a)') "Error! ... Output directory ''" // IO_TDA // "'' exists & does not contain a previous TDA run ... " // &
                                 'it is NOT overwritten ... choose a different TDA_Parameters.Name ...'
                return
            end if
        end if
        call make_output_dirs(IO_TDA)
        ! Marker of a TDA output directory (robustness deviation: dev never writes it) ...
        open(newunit=u, file=IO_TDA // '/.TDA_output', status='replace', action='write')
        close(u)

        write(*,'(/,a)') 'Starting TDA calculation ...'
        call timer_start(t0)
        call TDA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'TDA calculation has finished ...'
    end subroutine TDA

    ! Reference import shared by TDA & RPA ... counterpart of TDA_solver_import & the import part of RPA_solver ...
    subroutine reference_import(Params,Reference,Read_Sigma,Tol_par,Tol_hol,Orb,Orb_1p1h,C,Rho,h_N,V_NN,Orb_NN,Sigma)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Reference
        logical, intent(in) :: Read_Sigma
        real(dp), intent(in) :: Tol_par, Tol_hol
        type(orb1b_t), allocatable, intent(out) :: Orb(:)
        type(orb2b_1p1h_t), intent(out) :: Orb_1p1h
        type(o1b_t), intent(out) :: C, Rho, h_N
        type(o2b_t), intent(out) :: V_NN, Sigma
        type(orb2b_t), intent(out) :: Orb_NN
        character(len=:), allocatable :: IO
        integer(i8) :: t0

        IO = io_dir(Params,Reference)
        Orb = orbitals_import(Params,IO // '/bin/orb.bin')
        C = O1b_import(Params,Orb,IO // '/bin/C.bin')
        Rho = O1b_import(Params,Orb,IO // '/bin/rho.bin')
        h_N = O1b_import(Params,Orb,IO // '/bin/H1b.bin')
        write(*,'(/,a)') 'Importing reference 2-body interaction ...'
        call timer_start(t0)
        call O2b_import(Params,Orb,IO // '/bin/H2b.bin',V_NN,Orb_NN)
        call timer_report(t0)
        if (Read_Sigma) then
            write(*,'(/,a)') 'Importing reference 2-body correlation function ...'
            call timer_start(t0)
            call O2b_import(Params,Orb,IO // '/bin/sigma.bin',Sigma)
            call timer_report(t0)
        end if
        Orb_1p1h = orbitals_1p1h_make(Params,Orb,Tol_par,Tol_hol)
    end subroutine reference_import

    ! TDA solver ... counterpart of TDA_solver (TDA_solver.jl) ...
    subroutine TDA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(orb2b_1p1h_t) :: Orb_1p1h
        type(o1b_t) :: C, Rho, h_N
        type(o2b_t) :: V_NN, Sigma
        type(orb2b_t) :: Orb_NN
        type(tr1b_t) :: TrOp
        type(rmat_t), allocatable :: A(:,:), B(:,:), N(:,:), X(:,:)
        type(rvec_t), allocatable :: E(:,:)
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        real(dp) :: chR2
        integer(i8) :: t0

        call reference_import(Params,Params%Calc%TDA%Reference,Params%Calc%TDA%Reference_TBCF,Params%Calc%TDA%Tol_par, &
                              Params%Calc%TDA%Tol_hol,Orb,Orb_1p1h,C,Rho,h_N,V_NN,Orb_NN,Sigma)

        ! Transition operators with the reference charge radius ...
        chR2 = rho_chR2(Params,Orb,C,Rho)
        TrOp = Tr1b_initialize(Params,Orb,chR2)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        write(*,'(/,a)') 'Allocating TDA matrices A & N ...'
        call timer_start(t0)
        call ph_allocate(Orb,Orb_1p1h,Orb_NN,h_N,V_NN,.false.,A,B,N)
        if (Params%Calc%TDA%Reference_TBCF) then
            allocate(B(2,0:Orb_1p1h%J_max))
            call ph_allocate_tbcf(Params,Orb,Orb_1p1h,Orb_NN,V_NN,Sigma,.false.,.false.,A,B,'TDA')
        end if
        write(*,'(/,a)') TAB // 'TDA matrices A & N allocated ...'
        call timer_report(t0)

        call timer_start(t0)
        call TDA_diagonalize(Params,Orb,Orb_1p1h,A,N,TrOp,E,X)
        call timer_report(t0)

        call timer_start(t0)
        rM = TDA_rM(Params,Orb,Orb_1p1h,Orb_NN,TrOp,Rho,Sigma,X)
        rB = ph_rB(Params,Orb_1p1h,rM,E(2,1)%v)
        call timer_report(t0)

        call TDA_export(Params,Orb,Orb_1p1h,E,X,rB,C,Rho)
    end subroutine TDA_solver

    ! RPA driver ... counterpart of RPA(Params) (RPA/RPA.jl:1-70) ... the output directory is replaced unconditionally (dev) ...
    subroutine RPA(Params)
        type(params_t), intent(in) :: Params
        character(len=:), allocatable :: IO_RPA
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        call print_header(Params,'Starting RPA calculations with residual NO2B NN+NNN interaction & for given reference solutions ...')
        call print_threads()

        if (.not. is_dir(io_dir(Params,Params%Calc%RPA%Reference))) then
            write(*,'(/,a)') "Error! ... No precomputed reference solutions available in IO directory ''" // io_path(Params,'') // &
                             "'' for given input path ''" // trim(Params%Calc%RPA%Reference) // "''... e.g. run HF solver first ..."
            return
        end if

        IO_RPA = io_dir(Params,Params%Calc%RPA%Name)
        if (is_dir(IO_RPA)) call remove_dir(IO_RPA)
        call make_output_dirs(IO_RPA)

        write(*,'(/,a)') 'Starting RPA calculation ...'
        call timer_start(t0)
        call RPA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'RPA calculation has finished ...'
    end subroutine RPA

    ! RPA solver ... counterpart of RPA_solver (RPA_solver.jl:1-92) ...
    subroutine RPA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(orb2b_1p1h_t) :: Orb_1p1h
        type(o1b_t) :: C, Rho, Rho_ref, h_N, Rho_LHO
        type(o2b_t) :: V_NN, Sigma_ref, Sigma
        type(orb2b_t) :: Orb_NN
        type(tr1b_t) :: TrOp
        type(rmat_t), allocatable :: A(:,:), B(:,:), N(:,:)
        type(cvec_t), allocatable :: E(:,:)
        type(cmat_t), allocatable :: X(:,:), Y(:,:)
        type(rpa_energy_t) :: E_corr
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        logical :: Stability, TBCF_allocate
        real(dp) :: chR2
        integer(i8) :: t0

        ! Robustness guard (dev silently ignores Reference_TBCF without the reference sigma.bin) ...
        TBCF_allocate = Params%Calc%RPA%Reference_TBCF
        if (TBCF_allocate .and. .not. is_file(io_dir(Params,Params%Calc%RPA%Reference) // '/bin/sigma.bin')) &
            call fatal('RPA: Reference_TBCF = .true., but no reference 2-body correlation function ' // &
                       io_dir(Params,Params%Calc%RPA%Reference) // '/bin/sigma.bin is available ...')

        call reference_import(Params,Params%Calc%RPA%Reference,TBCF_allocate,Params%Calc%RPA%Tol_par, &
                              Params%Calc%RPA%Tol_hol,Orb,Orb_1p1h,C,Rho,h_N,V_NN,Orb_NN,Sigma_ref)

        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        call timer_start(t0)
        call RPA_allocate(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Rho,Sigma_ref,TBCF_allocate,A,B,N)
        call timer_report(t0)

        call timer_start(t0)
        call RPA_diagonalize(Params,Orb,Orb_1p1h,A,B,N,TrOp,E,X,Y,Stability)
        call timer_report(t0)
        deallocate(A, B, N)

        ! RPA OBDM (reference basis, diagonal) & 2-body correlation function ...
        Rho_ref = Rho
        Rho = RPA_OBDM_oY2(Params,Orb,Orb_1p1h,Y)
        if (Params%Calc%RPA%TBCF) then
            call timer_start(t0)
            Sigma = RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y)
            call timer_report(t0)
        end if

        ! Transition operators with the RPA charge radius ...
        chR2 = rho_chR2(Params,Orb,C,Rho)
        TrOp = Tr1b_initialize(Params,Orb,chR2)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        call timer_start(t0)
        E_corr = RPA_energy(Params,Orb,Orb_1p1h,Orb_NN,E,X,Y,Rho,Rho_ref,h_N,Sigma,V_NN,TBCF_allocate,Sigma_ref)
        call timer_report(t0)

        call timer_start(t0)
        rM = RPA_rM(Params,Orb,Orb_1p1h,Orb_NN,TrOp,Rho,Sigma,X,Y)
        rB = ph_rB(Params,Orb_1p1h,rM,real(E(2,1)%v,dp))
        call timer_report(t0)

        ! RPA OBDM in the LHO basis ...
        Rho_LHO%p = matmul(matmul(C%p,Rho%p),transpose(C%p))
        Rho_LHO%n = matmul(matmul(C%n,Rho%n),transpose(C%n))

        call RPA_export(Params,Orb,Orb_1p1h,Orb_NN,Stability,E_corr,E,X,Y,rB,C,Rho_LHO,Sigma,h_N)
    end subroutine RPA_solver

end module nj_tda_rpa_driver
