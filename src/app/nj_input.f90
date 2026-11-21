! Namelist input ... replaces the keyword constructors Interaction_Parameters(...), Calculation_Parameters(...) etc.
    ! of the Julia scripts (counterpart of src/import_structures.jl of NuJuliet dev). Every group is optional; omitted
    ! values keep the defaults of the Julia structs ... one contained reader per group with its own local variables
    ! (the same item name may have different types in different groups, e.g. Transitions in &rpa & &ncsm) ...
module nj_input
    use nj_kinds
    use nj_params
    use nj_utils
    implicit none
    private

    integer, parameter, public :: Max_Actions = 64

    ! Extra arguments of the phonon density & current exports ...
    type, public :: export_request_t
        character(len=16) :: JP = '1-'
        integer :: M = 0
        integer :: N_nu = 0
        integer :: nu_list(256) = 0
        character(len=Str_Len) :: File_Name = ''
        integer :: xN_grid = 40
        integer :: zN_grid = 40
        real(dp) :: Cut = 0.1_dp      ! &transition_analysis (stage 3)
    end type export_request_t

    ! Output truncation & files of the interaction parse actions (V2b_parse / V3b_no2b_parse, &parse) ... -1 = not given ...
    type, public :: parse_request_t
        integer :: N_max_Out = -1, N_2max_Out = -1, N_3max_Out = -1
        character(len=Str_Len) :: NN_Out_File = 'io/NN_parsed.bin'
        character(len=Str_Len) :: NNN_Out_File = 'io/NNN_parsed.bin'
    end type parse_request_t

    type, public :: run_input_t
        type(params_t) :: Params
        integer :: N_Actions = 0
        character(len=64) :: Actions(Max_Actions) = ''
        type(export_request_t) :: Densities, Currents
        type(export_request_t) :: Op_Densities, Analysis   ! &transition_operator_densities, &transition_analysis (stage 3)
        type(parse_request_t) :: Parse
    end type run_input_t

    public :: read_input

contains

    subroutine read_input(File,Inp)
        character(len=*), intent(in) :: File
        type(run_input_t), intent(out) :: Inp
        integer :: u, Ios
        character(len=512) :: Msg

        open(newunit=u, file=File, status='old', action='read', iostat=Ios, iomsg=Msg)
        if (Ios /= 0) call fatal('Cannot open input file ' // trim(File) // ' ... ' // trim(Msg))

        call read_interaction(Inp%Params%Int)
        call read_calculation(Inp%Params%Calc)
        call read_hf(Inp%Params%Calc%HF)
        call read_hf_mbpt(Inp%Params%Calc%HF%MBPT)
        call read_hf_efa(Inp%Params%Calc%HF%EFA)
        call read_tda(Inp%Params%Calc%TDA)
        call read_rpa(Inp%Params%Calc%RPA)
        call read_rrpa(Inp%Params%Calc%RRPA)
        call read_bcs(Inp%Params%Calc%BCS)
        call read_hfb(Inp%Params%Calc%HFB)
        call read_qtda(Inp%Params%Calc%QTDA)
        call read_qrpa(Inp%Params%Calc%QRPA)
        call read_pnqtda(Inp%Params%Calc%pnQTDA)
        call read_empm(Inp%Params%Calc%EMPM)
        call read_ncsm(Inp%Params%Calc%NCSM)
        call read_run()
        call read_transition_densities(Inp%Densities)
        call read_transition_currents(Inp%Currents)
        call read_transition_operator_densities(Inp%Op_Densities)
        call read_transition_analysis(Inp%Analysis)
        call read_parse(Inp%Parse)

        close(u)

    contains

        ! A missing group (end of file) keeps the defaults, anything else is an input error ...
        subroutine check(Ios,Msg,Group)
            integer, intent(in) :: Ios
            character(len=*), intent(in) :: Msg, Group
            if (Ios == 0 .or. is_iostat_end(Ios)) return
            call fatal('Error in namelist &' // Group // ' of ' // trim(File) // ' ... ' // trim(Msg))
        end subroutine check

        ! Name & Reference are folder names below io/<Path>/ ... non-empty, no '/' ...
        subroutine check_folder_name(Name,Group,Item)
            character(len=*), intent(in) :: Name, Group, Item
            if (len_trim(Name) == 0) call fatal('&' // Group // ' ' // Item // ' must not be empty ...')
            if (index(Name,'/') > 0) call fatal('&' // Group // ' ' // Item // ' must not contain a / ... ' // trim(Name))
        end subroutine check_folder_name

        subroutine read_interaction(P)
            type(interaction_params_t), intent(inout) :: P
            character(len=Str_Len) :: NN_File, NNN_File
            real(dp) :: hw, cP2N, cP3N, cRes, pG2N, nG2N
            integer :: Nmax, N2max, N3max, Ios
            character(len=512) :: Msg
            namelist /interaction/ NN_File, NNN_File, hw, Nmax, N2max, N3max, cP2N, cP3N, cRes, pG2N, nG2N
            NN_File = P%NN_File
            NNN_File = P%NNN_File
            hw = P%hw
            Nmax = P%Nmax
            N2max = P%N2max
            N3max = P%N3max
            cP2N = P%cP2N
            cP3N = P%cP3N
            cRes = P%cRes
            pG2N = P%pG2N
            nG2N = P%nG2N
            rewind(u)
            read(u, nml=interaction, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'interaction')
            P%NN_File = NN_File
            P%NNN_File = NNN_File
            P%hw = hw
            P%Nmax = Nmax
            P%N2max = N2max
            P%N3max = N3max
            P%cP2N = cP2N
            P%cP3N = cP3N
            P%cRes = cRes
            P%pG2N = pG2N
            P%nG2N = nG2N
        end subroutine read_interaction

        subroutine read_calculation(P)
            type(calculation_params_t), intent(inout) :: P
            integer :: A, Z, Nmax, N2max, N3max, Ios
            real(dp) :: hw, CM_Lawson
            character(len=32) :: CMS
            character(len=Str_Len) :: Path
            character(len=512) :: Msg
            namelist /calculation/ A, Z, hw, Nmax, N2max, N3max, CMS, CM_Lawson, Path
            A = P%A
            Z = P%Z
            hw = P%hw
            Nmax = P%Nmax
            N2max = P%N2max
            N3max = P%N3max
            CMS = P%CMS
            CM_Lawson = P%CM_Lawson
            Path = ''
            rewind(u)
            read(u, nml=calculation, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'calculation')
            P%A = A
            P%Z = Z
            P%hw = hw
            P%Nmax = Nmax
            P%N2max = N2max
            P%N3max = N3max
            P%CMS = CMS
            P%CM_Lawson = CM_Lawson
            if (len_trim(Path) == 0) then
                P%Path = default_path(P)
            else
                P%Path = Path
            end if
        end subroutine read_calculation

        subroutine read_hf(P)
            type(hf_params_t), intent(inout) :: P
            real(dp) :: Tol
            integer :: Imax, Ios
            logical :: V_HF, BMF, HRF, ESH
            character(len=512) :: Msg
            namelist /hf/ Tol, Imax, V_HF, BMF, HRF, ESH
            Tol = P%Tol
            Imax = P%Imax
            V_HF = P%V_HF
            BMF = P%BMF
            HRF = P%HRF
            ESH = P%ESH
            rewind(u)
            read(u, nml=hf, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'hf')
            P%Tol = Tol
            P%Imax = Imax
            P%V_HF = V_HF
            P%BMF = BMF
            P%HRF = HRF
            P%ESH = ESH
        end subroutine read_hf

        ! Julia HF.MBPT = HF_MBPT_Parameters(...) ...
        subroutine read_hf_mbpt(P)
            type(hf_mbpt_params_t), intent(inout) :: P
            logical :: Export, TBCF, ScH
            integer :: Ios
            character(len=512) :: Msg
            namelist /hf_mbpt/ Export, TBCF, ScH
            Export = P%Export
            TBCF = P%TBCF
            ScH = P%ScH
            rewind(u)
            read(u, nml=hf_mbpt, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'hf_mbpt')
            P%Export = Export
            P%TBCF = TBCF
            P%ScH = ScH
        end subroutine read_hf_mbpt

        ! LEGACY (removed in cleanup stage) ... &hf_efa is gone in NuJuliet dev ...
        subroutine read_hf_efa(P)
            type(hf_efa_params_t), intent(inout) :: P
            integer :: A, Z, Ios
            logical :: BMF
            character(len=512) :: Msg
            namelist /hf_efa/ A, Z, BMF
            A = P%A
            Z = P%Z
            BMF = P%BMF
            rewind(u)
            read(u, nml=hf_efa, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'hf_efa')
            if (Ios == 0) write(*,'(a)') 'WARNING: &hf_efa is a legacy group (HF-EFA is removed in NuJuliet dev) ...'
            P%A = A
            P%Z = Z
            P%BMF = BMF
        end subroutine read_hf_efa

        subroutine read_tda(P)
            type(tda_params_t), intent(inout) :: P
            character(len=64) :: Name, Reference
            logical :: Reference_TBCF, Ortho
            real(dp) :: Tol_hol, Tol_par
            integer :: Ios
            character(len=512) :: Msg
            namelist /tda/ Name, Reference, Reference_TBCF, Ortho, Tol_hol, Tol_par
            Name = P%Name
            Reference = P%Reference
            Reference_TBCF = P%Reference_TBCF
            Ortho = P%Ortho
            Tol_hol = P%Tol_hol
            Tol_par = P%Tol_par
            rewind(u)
            read(u, nml=tda, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'tda')
            P%Name = Name
            P%Reference = Reference
            P%Reference_TBCF = Reference_TBCF
            P%Ortho = Ortho
            P%Tol_hol = Tol_hol
            P%Tol_par = Tol_par
            ! Dev refuses Name == Reference (the output folder would overwrite / delete the reference) ...
            if (trim(P%Name) == trim(P%Reference)) call fatal('&tda Name must differ from Reference ... ' // trim(P%Name))
            call check_folder_name(P%Name,'tda','Name')
            call check_folder_name(P%Reference,'tda','Reference')
        end subroutine read_tda

        subroutine read_rpa(P)
            type(rpa_params_t), intent(inout) :: P
            character(len=64) :: Name, Reference
            character(len=16) :: Transitions
            logical :: Reference_TBCF, RPA_Orb, Ortho, OBDM_diagonal, ScH, TBCF, DEV
            real(dp) :: Tol, Tol_hol, Tol_par, OBDM_alpha, TBCF_alpha, TBCF_beta
            integer :: Imax, Ios
            character(len=512) :: Msg
            namelist /rpa/ Name, Reference, Reference_TBCF, RPA_Orb, Ortho, Transitions, Tol, Tol_hol, Tol_par, Imax, &
                           OBDM_diagonal, OBDM_alpha, ScH, TBCF, TBCF_alpha, TBCF_beta, DEV
            Name = P%Name
            Reference = P%Reference
            Reference_TBCF = P%Reference_TBCF
            RPA_Orb = P%RPA_Orb
            Ortho = P%Ortho
            Transitions = P%Transitions
            Tol = P%Tol
            Tol_hol = P%Tol_hol
            Tol_par = P%Tol_par
            Imax = P%Imax
            OBDM_diagonal = P%OBDM_diagonal
            OBDM_alpha = P%OBDM_alpha
            ScH = P%ScH
            TBCF = P%TBCF
            TBCF_alpha = P%TBCF_alpha
            TBCF_beta = P%TBCF_beta
            DEV = P%DEV
            rewind(u)
            read(u, nml=rpa, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'rpa')
            P%Name = Name
            P%Reference = Reference
            P%Reference_TBCF = Reference_TBCF
            P%RPA_Orb = RPA_Orb
            P%Ortho = Ortho
            P%Transitions = Transitions
            P%Tol = Tol
            P%Tol_hol = Tol_hol
            P%Tol_par = Tol_par
            P%Imax = Imax
            P%OBDM_diagonal = OBDM_diagonal
            P%OBDM_alpha = OBDM_alpha
            P%ScH = ScH
            P%TBCF = TBCF
            P%TBCF_alpha = TBCF_alpha
            P%TBCF_beta = TBCF_beta
            P%DEV = DEV
            ! Dev RPA_external.jl:123 errors for other values (checked here already, before any solver runs) ...
            if (trim(P%Transitions) /= 'Exact' .and. trim(P%Transitions) /= 'Linear') &
                call fatal("&rpa Transitions must be 'Exact' or 'Linear' ... " // trim(P%Transitions))
            if (trim(P%Name) == trim(P%Reference)) call fatal('&rpa Name must differ from Reference ... ' // trim(P%Name))
            call check_folder_name(P%Name,'rpa','Name')
            call check_folder_name(P%Reference,'rpa','Reference')
        end subroutine read_rpa

        ! LEGACY (removed in cleanup stage) ... &rrpa is gone in NuJuliet dev ...
        subroutine read_rrpa(P)
            type(rrpa_params_t), intent(inout) :: P
            logical :: Ortho, ScV3N, ScOBH, dOBDM
            real(dp) :: Tol
            integer :: Imax, Ios
            character(len=512) :: Msg
            namelist /rrpa/ Ortho, Tol, Imax, ScV3N, ScOBH, dOBDM
            Ortho = P%Ortho
            Tol = P%Tol
            Imax = P%Imax
            ScV3N = P%ScV3N
            ScOBH = P%ScOBH
            dOBDM = P%dOBDM
            rewind(u)
            read(u, nml=rrpa, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'rrpa')
            if (Ios == 0) write(*,'(a)') 'WARNING: &rrpa is a legacy group (HF-RRPA is removed in NuJuliet dev) ...'
            P%Ortho = Ortho
            P%Tol = Tol
            P%Imax = Imax
            P%ScV3N = ScV3N
            P%ScOBH = ScOBH
            P%dOBDM = dOBDM
        end subroutine read_rrpa

        subroutine read_bcs(P)
            type(bcs_params_t), intent(inout) :: P
            logical :: ScBCS, BMF, HRF
            real(dp) :: Tol, q, pD0, nD0
            integer :: Imax, A0, Z0, Ios
            character(len=512) :: Msg
            namelist /bcs/ ScBCS, Tol, Imax, q, A0, Z0, pD0, nD0, BMF, HRF
            ScBCS = P%ScBCS
            Tol = P%Tol
            Imax = P%Imax
            q = P%q
            A0 = P%A0
            Z0 = P%Z0
            pD0 = P%pD0
            nD0 = P%nD0
            BMF = P%BMF
            HRF = P%HRF
            rewind(u)
            read(u, nml=bcs, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'bcs')
            P%ScBCS = ScBCS
            P%Tol = Tol
            P%Imax = Imax
            P%q = q
            P%A0 = A0
            P%Z0 = Z0
            P%pD0 = pD0
            P%nD0 = nD0
            P%BMF = BMF
            P%HRF = HRF
        end subroutine read_bcs

        ! Broy_Qmax & Broy_Bmax (removed in dev) are still accepted (LEGACY, removed in cleanup stage) & ignored ...
        subroutine read_hfb(P)
            type(hfb_params_t), intent(inout) :: P
            real(dp) :: Tol, dLmax, pL0, nL0, pa0, na0, Broy_Amax, Broy_Bmax, Broy_Tmax
            integer :: Imax, Broy_m, Broy_Qmax, Ios
            character(len=32) :: Pairing
            logical :: LNT, LNRes, BMF, HRF
            character(len=512) :: Msg
            namelist /hfb/ Tol, Imax, Pairing, dLmax, pL0, nL0, pa0, na0, Broy_m, Broy_Qmax, Broy_Amax, Broy_Bmax, &
                           Broy_Tmax, LNT, LNRes, BMF, HRF
            Tol = P%Tol
            Imax = P%Imax
            Pairing = P%Pairing
            dLmax = P%dLmax
            pL0 = P%pL0
            nL0 = P%nL0
            pa0 = P%pa0
            na0 = P%na0
            Broy_m = P%Broy_m
            Broy_Qmax = -huge(1)
            Broy_Amax = P%Broy_Amax
            Broy_Bmax = -huge(1.0_dp)
            Broy_Tmax = P%Broy_Tmax
            LNT = P%LNT
            LNRes = P%LNRes
            BMF = P%BMF
            HRF = P%HRF
            rewind(u)
            read(u, nml=hfb, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'hfb')
            if (Broy_Qmax /= -huge(1) .or. Broy_Bmax /= -huge(1.0_dp)) &
                write(*,'(a)') 'WARNING: &hfb Broy_Qmax & Broy_Bmax are removed in NuJuliet dev ... ignored ...'
            P%Tol = Tol
            P%Imax = Imax
            P%Pairing = Pairing
            P%dLmax = dLmax
            P%pL0 = pL0
            P%nL0 = nL0
            P%pa0 = pa0
            P%na0 = na0
            P%Broy_m = Broy_m
            P%Broy_Amax = Broy_Amax
            P%Broy_Tmax = Broy_Tmax
            P%LNT = LNT
            P%LNRes = LNRes
            P%BMF = BMF
            P%HRF = HRF
        end subroutine read_hfb

        subroutine read_qtda(P)
            type(qtda_params_t), intent(inout) :: P
            logical :: Ortho
            character(len=64) :: Reference
            integer :: Ios
            character(len=512) :: Msg
            namelist /qtda/ Ortho, Reference
            Ortho = P%Ortho
            Reference = P%Reference
            rewind(u)
            read(u, nml=qtda, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'qtda')
            P%Ortho = Ortho
            P%Reference = Reference
            call check_folder_name(P%Reference,'qtda','Reference')
        end subroutine read_qtda

        subroutine read_qrpa(P)
            type(qrpa_params_t), intent(inout) :: P
            logical :: Ortho
            character(len=64) :: Reference
            integer :: Ios
            character(len=512) :: Msg
            namelist /qrpa/ Ortho, Reference
            Ortho = P%Ortho
            Reference = P%Reference
            rewind(u)
            read(u, nml=qrpa, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'qrpa')
            P%Ortho = Ortho
            P%Reference = Reference
            call check_folder_name(P%Reference,'qrpa','Reference')
        end subroutine read_qrpa

        subroutine read_pnqtda(P)
            type(pnqtda_params_t), intent(inout) :: P
            character(len=64) :: Reference
            integer :: Ios
            character(len=512) :: Msg
            namelist /pnqtda/ Reference
            Reference = P%Reference
            rewind(u)
            read(u, nml=pnqtda, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'pnqtda')
            P%Reference = Reference
            call check_folder_name(P%Reference,'pnqtda','Reference')
        end subroutine read_pnqtda

        subroutine read_empm(P)
            type(empm_params_t), intent(inout) :: P
            character(len=64) :: Reference
            integer :: pHole(2), pPart(2), nHole(2), nPart(2)
            logical :: Ortho
            real(dp) :: E2max, Mem_GB, E1max, Cache_GB, E3max, E4max
            character(len=8) :: JP2(64), JP3(64), JP4(64), JPF(64)
            integer :: Full_N
            logical :: Save_X, Ind_First
            character(len=512) :: States_Dir
            logical :: Legacy_Order, Dry_Run, Legacy_Spur, Full_H, D_Spectrum
            character(len=8) :: Spur_Rows
            logical :: CM_Exact, CM_Dump, Ind_Greedy, CM_Variants, CM_Keep, CM_Ncm, Dens_Share
            character(len=8) :: Full_Diag
            logical :: Full_Svd
            real(dp) :: CM_Tol, CM_Beta, D_Disk_GB
            integer :: SF_Steps
            character(len=8) :: SF_Op, SF_Mode
            logical :: SF_Reorth, SF_Qrows
            integer :: Ios
            character(len=512) :: Msg
            namelist /empm/ Reference, pHole, pPart, nHole, nPart, Ortho, E2max, JP2, Mem_GB, E1max, Legacy_Order, Dry_Run, Spur_Rows, Legacy_Spur, Full_H, D_Spectrum, Cache_GB, &
                           E3max, JP3, E4max, JP4, Full_N, JPF, Save_X, Ind_First, States_Dir, CM_Exact, CM_Tol, CM_Dump, Ind_Greedy, CM_Variants, CM_Keep, Full_Diag, Full_Svd, CM_Ncm, CM_Beta, D_Disk_GB, Dens_Share, SF_Steps, SF_Op, SF_Reorth, SF_Mode, SF_Qrows
            Reference = P%Reference
            pHole = P%pHole
            pPart = P%pPart
            nHole = P%nHole
            nPart = P%nPart
            Ortho = P%Ortho
            E2max = P%E2max
            JP2 = P%JP2
            Mem_GB = P%Mem_GB
            E1max = P%E1max
            Legacy_Order = P%Legacy_Order
            Dry_Run = P%Dry_Run
            Spur_Rows = P%Spur_Rows
            Legacy_Spur = P%Legacy_Spur
            Full_H = P%Full_H
            D_Spectrum = P%D_Spectrum
            Cache_GB = P%Cache_GB
            E3max = P%E3max
            JP3 = P%JP3
            E4max = P%E4max
            JP4 = P%JP4
            Full_N = P%Full_N
            JPF = P%JPF
            Save_X = P%Save_X
            Ind_First = P%Ind_First
            States_Dir = P%States_Dir
            CM_Exact = P%CM_Exact
            CM_Tol = P%CM_Tol
            CM_Dump = P%CM_Dump
            Ind_Greedy = P%Ind_Greedy
            CM_Variants = P%CM_Variants
            CM_Keep = P%CM_Keep
            Full_Diag = P%Full_Diag
            Full_Svd = P%Full_Svd
            CM_Ncm = P%CM_Ncm
            CM_Beta = P%CM_Beta
            D_Disk_GB = P%D_Disk_GB
            Dens_Share = P%Dens_Share
            SF_Steps = P%SF_Steps
            SF_Op = P%SF_Op
            SF_Reorth = P%SF_Reorth
            SF_Mode = P%SF_Mode
            SF_Qrows = P%SF_Qrows
            rewind(u)
            read(u, nml=empm, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'empm')
            ! component-wise (a positional structure constructor breaks silently when a field is inserted) ...
            P%Reference = Reference
            P%pHole = pHole
            P%pPart = pPart
            P%nHole = nHole
            P%nPart = nPart
            P%Ortho = Ortho
            P%E2max = E2max
            P%JP2 = JP2
            P%Mem_GB = Mem_GB
            P%E1max = E1max
            P%Legacy_Order = Legacy_Order
            P%Dry_Run = Dry_Run
            P%Spur_Rows = Spur_Rows
            P%Legacy_Spur = Legacy_Spur
            P%Full_H = Full_H
            P%D_Spectrum = D_Spectrum
            P%Cache_GB = Cache_GB
            P%E3max = E3max
            P%JP3 = JP3
            P%E4max = E4max
            P%JP4 = JP4
            P%Full_N = Full_N
            P%JPF = JPF
            P%Save_X = Save_X
            P%Ind_First = Ind_First
            P%States_Dir = States_Dir
            P%CM_Exact = CM_Exact
            P%CM_Tol = CM_Tol
            P%CM_Dump = CM_Dump
            P%Ind_Greedy = Ind_Greedy
            P%CM_Variants = CM_Variants
            P%CM_Keep = CM_Keep
            P%Full_Diag = Full_Diag
            P%Full_Svd = Full_Svd
            P%CM_Ncm = CM_Ncm
            P%CM_Beta = CM_Beta
            P%D_Disk_GB = D_Disk_GB
            P%Dens_Share = Dens_Share
            P%SF_Steps = SF_Steps
            P%SF_Op = SF_Op
            P%SF_Reorth = SF_Reorth
            P%SF_Mode = SF_Mode
            P%SF_Qrows = SF_Qrows
            call check_folder_name(P%Reference,'empm','Reference')
        end subroutine read_empm

        subroutine read_ncsm(P)
            type(ncsm_params_t), intent(inout) :: P
            integer :: Nmax_List(16), M2, N_Eig, Lanczos_Max, Lanczos_Basis
            character(len=8) :: Hamiltonian, Matrix, Trunc, Basis
            real(dp) :: Lanczos_Tol, Beta_CM, Matrix_GB, Ecut_npnh(16)
            logical :: Single, Transitions, Window
            integer :: Dense_Max, Parity, Tr_Lanczos, Tr_Gs_Nmax, Tr_Ckpt
            character(len=8) :: Tr_From
            character(len=32) :: Tr_Ops
            character(len=8) :: Tr_Vec
            logical :: Tr_Reorth, Tr_Sym
            integer :: Ios
            character(len=512) :: Msg
            namelist /ncsm/ Nmax_List, Hamiltonian, M2, N_Eig, Lanczos_Max, Lanczos_Tol, Beta_CM, Matrix, Single, Lanczos_Basis, Matrix_GB, Trunc, Basis, Ecut_npnh, Transitions, Dense_Max, Window, Parity, Tr_From, Tr_Lanczos, Tr_Gs_Nmax, Tr_Ops, Tr_Vec, Tr_Reorth, Tr_Sym, Tr_Ckpt
            Nmax_List = P%Nmax_List
            Hamiltonian = P%Hamiltonian
            Trunc = P%Trunc
            Basis = P%Basis
            M2 = P%M2
            N_Eig = P%N_Eig
            Lanczos_Max = P%Lanczos_Max
            Lanczos_Tol = P%Lanczos_Tol
            Beta_CM = P%Beta_CM
            Matrix = P%Matrix
            Single = P%Single
            Lanczos_Basis = P%Lanczos_Basis
            Matrix_GB = P%Matrix_GB
            Ecut_npnh = P%Ecut_npnh
            Transitions = P%Transitions
            Dense_Max = P%Dense_Max
            Window = P%Window
            Parity = P%Parity
            Tr_From = P%Tr_From
            Tr_Lanczos = P%Tr_Lanczos
            Tr_Gs_Nmax = P%Tr_Gs_Nmax
            Tr_Ops = P%Tr_Ops
            Tr_Vec = P%Tr_Vec
            Tr_Reorth = P%Tr_Reorth
            Tr_Sym = P%Tr_Sym
            Tr_Ckpt = P%Tr_Ckpt
            rewind(u)
            read(u, nml=ncsm, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'ncsm')
            P = ncsm_params_t(Nmax_List,Hamiltonian,M2,N_Eig,Lanczos_Max,Lanczos_Tol,Beta_CM,Matrix,Single,Lanczos_Basis,Matrix_GB,Trunc,Basis,Ecut_npnh,Transitions,Dense_Max,Window,Parity,Tr_From,Tr_Lanczos,Tr_Gs_Nmax,Tr_Ops,Tr_Vec,Tr_Sym,Tr_Reorth,Tr_Ckpt)
        end subroutine read_ncsm

        subroutine read_run()
            character(len=64) :: Solvers(Max_Actions)
            integer :: k, Ios
            character(len=512) :: Msg
            namelist /run/ Solvers
            Solvers = ''
            rewind(u)
            read(u, nml=run, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'run')
            Inp%N_Actions = 0
            do k = 1, Max_Actions
                if (len_trim(Solvers(k)) > 0) then
                    Inp%N_Actions = Inp%N_Actions + 1
                    Inp%Actions(Inp%N_Actions) = adjustl(Solvers(k))
                end if
            end do
        end subroutine read_run

        subroutine read_transition_densities(X)
            type(export_request_t), intent(inout) :: X
            character(len=16) :: JP
            integer :: nu_list(256), Ios
            character(len=Str_Len) :: File_Name
            character(len=512) :: Msg
            namelist /transition_densities/ JP, nu_list, File_Name
            JP = X%JP
            nu_list = 0
            File_Name = 'Phonon_Densities'
            rewind(u)
            read(u, nml=transition_densities, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'transition_densities')
            X%JP = JP
            X%nu_list = nu_list
            X%N_nu = count(nu_list > 0)
            X%File_Name = File_Name
        end subroutine read_transition_densities

        ! Julia M = -1 (default) means M = 0 ...
        subroutine read_transition_currents(X)
            type(export_request_t), intent(inout) :: X
            character(len=16) :: JP
            integer :: M, nu_list(256), xN_grid, zN_grid, Ios
            character(len=Str_Len) :: File_Name
            character(len=512) :: Msg
            namelist /transition_currents/ JP, M, nu_list, File_Name, xN_grid, zN_grid
            JP = X%JP
            M = X%M
            nu_list = 0
            File_Name = 'Phonon_Currents'
            xN_grid = X%xN_grid
            zN_grid = X%zN_grid
            rewind(u)
            read(u, nml=transition_currents, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'transition_currents')
            X%JP = JP
            X%M = M
            X%nu_list = nu_list
            X%N_nu = count(nu_list > 0)
            X%File_Name = File_Name
            X%xN_grid = xN_grid
            X%zN_grid = zN_grid
        end subroutine read_transition_currents

        ! TDA_/RPA_transition_operator_densities (stage 3) ... dev kwargs JP, nu_list, File_Name ...
        subroutine read_transition_operator_densities(X)
            type(export_request_t), intent(inout) :: X
            character(len=16) :: JP
            integer :: nu_list(256), Ios
            character(len=Str_Len) :: File_Name
            character(len=512) :: Msg
            namelist /transition_operator_densities/ JP, nu_list, File_Name
            JP = X%JP
            nu_list = 0
            File_Name = 'Phonon_Transition_Operator_Densities'
            rewind(u)
            read(u, nml=transition_operator_densities, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'transition_operator_densities')
            X%JP = JP
            X%nu_list = nu_list
            X%N_nu = count(nu_list > 0)
            X%File_Name = File_Name
        end subroutine read_transition_operator_densities

        ! TDA_/RPA_transition_analysis (stage 3) ... dev kwargs JP, nu_list, File_Name (dev spelling "Transtion"), Cut ...
        subroutine read_transition_analysis(X)
            type(export_request_t), intent(inout) :: X
            character(len=16) :: JP
            integer :: nu_list(256), Ios
            character(len=Str_Len) :: File_Name
            real(dp) :: Cut
            character(len=512) :: Msg
            namelist /transition_analysis/ JP, nu_list, File_Name, Cut
            JP = X%JP
            nu_list = 0
            File_Name = 'Phonon_Transtion_Analysis'
            Cut = X%Cut
            rewind(u)
            read(u, nml=transition_analysis, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'transition_analysis')
            X%JP = JP
            X%nu_list = nu_list
            X%N_nu = count(nu_list > 0)
            X%File_Name = File_Name
            X%Cut = Cut
        end subroutine read_transition_analysis

        ! Output space & files of V2b_parse / V3b_no2b_parse (Julia function arguments N_max_Out, N_2max_Out, N_3max_Out,
            ! Out_File) ...
        subroutine read_parse(X)
            type(parse_request_t), intent(inout) :: X
            integer :: N_max_Out, N_2max_Out, N_3max_Out, Ios
            character(len=Str_Len) :: NN_Out_File, NNN_Out_File
            character(len=512) :: Msg
            namelist /parse/ N_max_Out, N_2max_Out, N_3max_Out, NN_Out_File, NNN_Out_File
            N_max_Out = X%N_max_Out
            N_2max_Out = X%N_2max_Out
            N_3max_Out = X%N_3max_Out
            NN_Out_File = X%NN_Out_File
            NNN_Out_File = X%NNN_Out_File
            rewind(u)
            read(u, nml=parse, iostat=Ios, iomsg=Msg)
            call check(Ios,Msg,'parse')
            X%N_max_Out = N_max_Out
            X%N_2max_Out = N_2max_Out
            X%N_3max_Out = N_3max_Out
            X%NN_Out_File = NN_Out_File
            X%NNN_Out_File = NNN_Out_File
        end subroutine read_parse

    end subroutine read_input

end module nj_input
