! Namelist input ... replaces the keyword constructors Interaction_Parameters(...), Calculation_Parameters(...) etc.
    ! of the Julia scripts. Every group is optional; omitted values keep the defaults of Import_Structures.jl ...
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
    end type export_request_t

    type, public :: run_input_t
        type(params_t) :: Params
        integer :: N_Actions = 0
        character(len=64) :: Actions(Max_Actions) = ''
        type(export_request_t) :: Densities, Currents
    end type run_input_t

    public :: read_input

contains

    subroutine read_input(File,Inp)
        character(len=*), intent(in) :: File
        type(run_input_t), intent(out) :: Inp
        type(interaction_params_t) :: I0
        type(calculation_params_t) :: C0
        type(export_request_t) :: X0
        integer :: u, Ios, k
        character(len=512) :: Msg

        ! &interaction ...
        character(len=Str_Len) :: NN_File, NNN_File
        real(dp) :: hw, cP2N, cP3N, cRes, pG2N, nG2N
        integer :: Nmax, N2max, N3max
        ! &calculation ...
        integer :: A, Z
        character(len=32) :: CMS
        real(dp) :: CM_Lawson
        character(len=Str_Len) :: Path
        ! &hf & &hf_efa ...
        real(dp) :: Tol
        integer :: Imax
        logical :: BMF, HRF
        ! &rpa, &qtda, &qrpa ...
        logical :: Ortho
        ! &rrpa ...
        logical :: ScV3N, ScOBH, dOBDM
        ! &bcs ...
        logical :: ScBCS
        real(dp) :: q, pD0, nD0
        integer :: A0, Z0
        ! &hfb ...
        character(len=32) :: Pairing
        real(dp) :: dLmax, pL0, nL0, pa0, na0, Broy_Amax, Broy_Bmax, Broy_Tmax
        integer :: Broy_m, Broy_Qmax
        logical :: LNT, LNRes
        ! &empm ...
        integer :: pHole(2), pPart(2), nHole(2), nPart(2)
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
        ! &ncsm ...
        integer :: Nmax_List(16), M2, N_Eig, Lanczos_Max, Lanczos_Basis
        character(len=8) :: Hamiltonian, Matrix, Trunc, Basis
        real(dp) :: Lanczos_Tol, Beta_CM, Matrix_GB, Ecut_npnh(16)
        logical :: Single, Transitions, Window
        integer :: Dense_Max, Parity, Tr_Lanczos, Tr_Gs_Nmax, Tr_Ckpt
        character(len=8) :: Tr_From
        character(len=32) :: Tr_Ops
        character(len=8) :: Tr_Vec
        logical :: Tr_Reorth, Tr_Sym
        ! &run ...
        character(len=64) :: Solvers(Max_Actions)
        ! &transition_densities & &transition_currents ...
        character(len=16) :: JP
        integer :: M, nu_list(256), xN_grid, zN_grid
        character(len=Str_Len) :: File_Name

        namelist /interaction/ NN_File, NNN_File, hw, Nmax, N2max, N3max, cP2N, cP3N, cRes, pG2N, nG2N
        namelist /calculation/ A, Z, hw, Nmax, N2max, N3max, CMS, CM_Lawson, Path
        namelist /hf/ Tol, Imax, BMF, HRF
        namelist /hf_efa/ A, Z, BMF
        namelist /rpa/ Ortho
        namelist /rrpa/ Ortho, Tol, Imax, ScV3N, ScOBH, dOBDM
        namelist /bcs/ ScBCS, Tol, Imax, q, A0, Z0, pD0, nD0, BMF, HRF
        namelist /hfb/ Tol, Imax, Pairing, dLmax, pL0, nL0, pa0, na0, Broy_m, Broy_Qmax, Broy_Amax, Broy_Bmax, &
                       Broy_Tmax, LNT, LNRes, BMF, HRF
        namelist /qtda/ Ortho
        namelist /qrpa/ Ortho
        namelist /empm/ pHole, pPart, nHole, nPart, Ortho, E2max, JP2, Mem_GB, E1max, Legacy_Order, Dry_Run, Spur_Rows, Legacy_Spur, Full_H, D_Spectrum, Cache_GB, &
                       E3max, JP3, E4max, JP4, Full_N, JPF, Save_X, Ind_First, States_Dir, CM_Exact, CM_Tol, CM_Dump, Ind_Greedy, CM_Variants, CM_Keep, Full_Diag, Full_Svd, CM_Ncm, CM_Beta, D_Disk_GB, Dens_Share, SF_Steps, SF_Op, SF_Reorth, SF_Mode, SF_Qrows
        namelist /ncsm/ Nmax_List, Hamiltonian, M2, N_Eig, Lanczos_Max, Lanczos_Tol, Beta_CM, Matrix, Single, Lanczos_Basis, Matrix_GB, Trunc, Basis, Ecut_npnh, Transitions, Dense_Max, Window, Parity, Tr_From, Tr_Lanczos, Tr_Gs_Nmax, Tr_Ops, Tr_Vec, Tr_Reorth, Tr_Sym, Tr_Ckpt
        namelist /run/ Solvers
        namelist /transition_densities/ JP, nu_list, File_Name
        namelist /transition_currents/ JP, M, nu_list, File_Name, xN_grid, zN_grid

        open(newunit=u, file=File, status='old', action='read', iostat=Ios, iomsg=Msg)
        if (Ios /= 0) call fatal('Cannot open input file ' // trim(File) // ' ... ' // trim(Msg))

        ! &interaction ...
        NN_File = I0%NN_File
        NNN_File = I0%NNN_File
        hw = I0%hw
        Nmax = I0%Nmax
        N2max = I0%N2max
        N3max = I0%N3max
        cP2N = I0%cP2N
        cP3N = I0%cP3N
        cRes = I0%cRes
        pG2N = I0%pG2N
        nG2N = I0%nG2N
        rewind(u)
        read(u, nml=interaction, iostat=Ios, iomsg=Msg)
        call check(Ios,'interaction')
        Inp%Params%Int = interaction_params_t(NN_File,NNN_File,hw,Nmax,N2max,N3max,cP2N,cP3N,cRes,pG2N,nG2N)

        ! &calculation ...
        A = C0%A
        Z = C0%Z
        hw = C0%hw
        Nmax = C0%Nmax
        N2max = C0%N2max
        N3max = C0%N3max
        CMS = C0%CMS
        CM_Lawson = C0%CM_Lawson
        Path = ''
        rewind(u)
        read(u, nml=calculation, iostat=Ios, iomsg=Msg)
        call check(Ios,'calculation')
        Inp%Params%Calc%A = A
        Inp%Params%Calc%Z = Z
        Inp%Params%Calc%hw = hw
        Inp%Params%Calc%Nmax = Nmax
        Inp%Params%Calc%N2max = N2max
        Inp%Params%Calc%N3max = N3max
        Inp%Params%Calc%CMS = CMS
        Inp%Params%Calc%CM_Lawson = CM_Lawson
        if (len_trim(Path) == 0) then
            Inp%Params%Calc%Path = default_path(Inp%Params%Calc)
        else
            Inp%Params%Calc%Path = Path
        end if

        ! &hf ...
        Tol = C0%HF%Tol
        Imax = C0%HF%Imax
        BMF = C0%HF%BMF
        HRF = C0%HF%HRF
        rewind(u)
        read(u, nml=hf, iostat=Ios, iomsg=Msg)
        call check(Ios,'hf')
        Inp%Params%Calc%HF%Tol = Tol
        Inp%Params%Calc%HF%Imax = Imax
        Inp%Params%Calc%HF%BMF = BMF
        Inp%Params%Calc%HF%HRF = HRF

        ! &hf_efa ...
        A = C0%HF%EFA%A
        Z = C0%HF%EFA%Z
        BMF = C0%HF%EFA%BMF
        rewind(u)
        read(u, nml=hf_efa, iostat=Ios, iomsg=Msg)
        call check(Ios,'hf_efa')
        Inp%Params%Calc%HF%EFA%A = A
        Inp%Params%Calc%HF%EFA%Z = Z
        Inp%Params%Calc%HF%EFA%BMF = BMF

        ! &rpa ...
        Ortho = C0%RPA%Ortho
        rewind(u)
        read(u, nml=rpa, iostat=Ios, iomsg=Msg)
        call check(Ios,'rpa')
        Inp%Params%Calc%RPA%Ortho = Ortho

        ! &rrpa ...
        Ortho = C0%RRPA%Ortho
        Tol = C0%RRPA%Tol
        Imax = C0%RRPA%Imax
        ScV3N = C0%RRPA%ScV3N
        ScOBH = C0%RRPA%ScOBH
        dOBDM = C0%RRPA%dOBDM
        rewind(u)
        read(u, nml=rrpa, iostat=Ios, iomsg=Msg)
        call check(Ios,'rrpa')
        Inp%Params%Calc%RRPA = rrpa_params_t(Ortho,Tol,Imax,ScV3N,ScOBH,dOBDM)

        ! &bcs ...
        ScBCS = C0%BCS%ScBCS
        Tol = C0%BCS%Tol
        Imax = C0%BCS%Imax
        q = C0%BCS%q
        A0 = C0%BCS%A0
        Z0 = C0%BCS%Z0
        pD0 = C0%BCS%pD0
        nD0 = C0%BCS%nD0
        BMF = C0%BCS%BMF
        HRF = C0%BCS%HRF
        rewind(u)
        read(u, nml=bcs, iostat=Ios, iomsg=Msg)
        call check(Ios,'bcs')
        Inp%Params%Calc%BCS = bcs_params_t(ScBCS,Tol,Imax,q,A0,Z0,pD0,nD0,BMF,HRF)

        ! &hfb ...
        Tol = C0%HFB%Tol
        Imax = C0%HFB%Imax
        Pairing = C0%HFB%Pairing
        dLmax = C0%HFB%dLmax
        pL0 = C0%HFB%pL0
        nL0 = C0%HFB%nL0
        pa0 = C0%HFB%pa0
        na0 = C0%HFB%na0
        Broy_m = C0%HFB%Broy_m
        Broy_Qmax = C0%HFB%Broy_Qmax
        Broy_Amax = C0%HFB%Broy_Amax
        Broy_Bmax = C0%HFB%Broy_Bmax
        Broy_Tmax = C0%HFB%Broy_Tmax
        LNT = C0%HFB%LNT
        LNRes = C0%HFB%LNRes
        BMF = C0%HFB%BMF
        HRF = C0%HFB%HRF
        rewind(u)
        read(u, nml=hfb, iostat=Ios, iomsg=Msg)
        call check(Ios,'hfb')
        Inp%Params%Calc%HFB = hfb_params_t(Tol,Imax,Pairing,dLmax,pL0,nL0,pa0,na0,Broy_m,Broy_Qmax,Broy_Amax, &
                                           Broy_Bmax,Broy_Tmax,LNT,LNRes,BMF,HRF)

        ! &qtda & &qrpa ...
        Ortho = C0%QTDA%Ortho
        rewind(u)
        read(u, nml=qtda, iostat=Ios, iomsg=Msg)
        call check(Ios,'qtda')
        Inp%Params%Calc%QTDA%Ortho = Ortho
        Ortho = C0%QRPA%Ortho
        rewind(u)
        read(u, nml=qrpa, iostat=Ios, iomsg=Msg)
        call check(Ios,'qrpa')
        Inp%Params%Calc%QRPA%Ortho = Ortho

        ! &empm ...
        pHole = C0%EMPM%pHole
        pPart = C0%EMPM%pPart
        nHole = C0%EMPM%nHole
        nPart = C0%EMPM%nPart
        Ortho = C0%EMPM%Ortho
        E2max = C0%EMPM%E2max
        JP2 = C0%EMPM%JP2
        Mem_GB = C0%EMPM%Mem_GB
        E1max = C0%EMPM%E1max
        Legacy_Order = C0%EMPM%Legacy_Order
        Dry_Run = C0%EMPM%Dry_Run
        Spur_Rows = C0%EMPM%Spur_Rows
        Legacy_Spur = C0%EMPM%Legacy_Spur
        Full_H = C0%EMPM%Full_H
        D_Spectrum = C0%EMPM%D_Spectrum
        Cache_GB = C0%EMPM%Cache_GB
        E3max = C0%EMPM%E3max
        JP3 = C0%EMPM%JP3
        E4max = C0%EMPM%E4max
        JP4 = C0%EMPM%JP4
        Full_N = C0%EMPM%Full_N
        JPF = C0%EMPM%JPF
        Save_X = C0%EMPM%Save_X
        Ind_First = C0%EMPM%Ind_First
        States_Dir = C0%EMPM%States_Dir
        CM_Exact = C0%EMPM%CM_Exact
        CM_Tol = C0%EMPM%CM_Tol
        CM_Dump = C0%EMPM%CM_Dump
        Ind_Greedy = C0%EMPM%Ind_Greedy
        CM_Variants = C0%EMPM%CM_Variants
        CM_Keep = C0%EMPM%CM_Keep
        Full_Diag = C0%EMPM%Full_Diag
        Full_Svd = C0%EMPM%Full_Svd
        CM_Ncm = C0%EMPM%CM_Ncm
        CM_Beta = C0%EMPM%CM_Beta
        D_Disk_GB = C0%EMPM%D_Disk_GB
        Dens_Share = C0%EMPM%Dens_Share
        SF_Steps = C0%EMPM%SF_Steps
        SF_Op = C0%EMPM%SF_Op
        SF_Reorth = C0%EMPM%SF_Reorth
        SF_Mode = C0%EMPM%SF_Mode
        SF_Qrows = C0%EMPM%SF_Qrows
        rewind(u)
        read(u, nml=empm, iostat=Ios, iomsg=Msg)
        call check(Ios,'empm')
        Inp%Params%Calc%EMPM = empm_params_t(pHole,pPart,nHole,nPart,Ortho,E2max,JP2,Mem_GB,E1max,Legacy_Order,Dry_Run,Spur_Rows,Legacy_Spur,Full_H,D_Spectrum,Cache_GB, &
                                                E3max,JP3,E4max,JP4,Full_N,JPF,Save_X,Ind_First,States_Dir,CM_Exact,CM_Tol,CM_Dump,Ind_Greedy,CM_Variants,CM_Keep,Full_Diag,Full_Svd,CM_Ncm,CM_Beta,D_Disk_GB,Dens_Share,SF_Steps,SF_Op,SF_Reorth,SF_Mode,SF_Qrows)

        ! &ncsm ...
        Nmax_List = C0%NCSM%Nmax_List
        Hamiltonian = C0%NCSM%Hamiltonian
        Trunc = C0%NCSM%Trunc
        Basis = C0%NCSM%Basis
        M2 = C0%NCSM%M2
        N_Eig = C0%NCSM%N_Eig
        Lanczos_Max = C0%NCSM%Lanczos_Max
        Lanczos_Tol = C0%NCSM%Lanczos_Tol
        Beta_CM = C0%NCSM%Beta_CM
        Matrix = C0%NCSM%Matrix
        Single = C0%NCSM%Single
        Lanczos_Basis = C0%NCSM%Lanczos_Basis
        Matrix_GB = C0%NCSM%Matrix_GB
        Ecut_npnh = C0%NCSM%Ecut_npnh
        Transitions = C0%NCSM%Transitions
        Dense_Max = C0%NCSM%Dense_Max
        Window = C0%NCSM%Window
        Parity = C0%NCSM%Parity
        Tr_From = C0%NCSM%Tr_From
        Tr_Lanczos = C0%NCSM%Tr_Lanczos
        Tr_Gs_Nmax = C0%NCSM%Tr_Gs_Nmax
        Tr_Ops = C0%NCSM%Tr_Ops
        Tr_Vec = C0%NCSM%Tr_Vec
        Tr_Reorth = C0%NCSM%Tr_Reorth
        Tr_Sym = C0%NCSM%Tr_Sym
        Tr_Ckpt = C0%NCSM%Tr_Ckpt
        rewind(u)
        read(u, nml=ncsm, iostat=Ios, iomsg=Msg)
        call check(Ios,'ncsm')
        Inp%Params%Calc%NCSM = ncsm_params_t(Nmax_List,Hamiltonian,M2,N_Eig,Lanczos_Max,Lanczos_Tol,Beta_CM,Matrix,Single,Lanczos_Basis,Matrix_GB,Trunc,Basis,Ecut_npnh,Transitions,Dense_Max,Window,Parity,Tr_From,Tr_Lanczos,Tr_Gs_Nmax,Tr_Ops,Tr_Vec,Tr_Sym,Tr_Reorth,Tr_Ckpt)

        ! &run ...
        Solvers = ''
        rewind(u)
        read(u, nml=run, iostat=Ios, iomsg=Msg)
        call check(Ios,'run')
        Inp%N_Actions = 0
        do k = 1, Max_Actions
            if (len_trim(Solvers(k)) > 0) then
                Inp%N_Actions = Inp%N_Actions + 1
                Inp%Actions(Inp%N_Actions) = adjustl(Solvers(k))
            end if
        end do

        ! &transition_densities ...
        JP = X0%JP
        nu_list = 0
        File_Name = 'Phonon_Densities'
        rewind(u)
        read(u, nml=transition_densities, iostat=Ios, iomsg=Msg)
        call check(Ios,'transition_densities')
        Inp%Densities%JP = JP
        Inp%Densities%nu_list = nu_list
        Inp%Densities%N_nu = count(nu_list > 0)
        Inp%Densities%File_Name = File_Name

        ! &transition_currents ...
        JP = X0%JP
        M = X0%M
        nu_list = 0
        File_Name = 'Phonon_Currents'
        xN_grid = X0%xN_grid
        zN_grid = X0%zN_grid
        rewind(u)
        read(u, nml=transition_currents, iostat=Ios, iomsg=Msg)
        call check(Ios,'transition_currents')
        Inp%Currents%JP = JP
        Inp%Currents%M = M
        Inp%Currents%nu_list = nu_list
        Inp%Currents%N_nu = count(nu_list > 0)
        Inp%Currents%File_Name = File_Name
        Inp%Currents%xN_grid = xN_grid
        Inp%Currents%zN_grid = zN_grid

        close(u)

    contains

        ! A missing group (end of file) keeps the defaults, anything else is an input error ...
        subroutine check(Ios,Group)
            integer, intent(in) :: Ios
            character(len=*), intent(in) :: Group
            if (Ios == 0 .or. is_iostat_end(Ios)) return
            call fatal('Error in namelist &' // Group // ' of ' // trim(File) // ' ... ' // trim(Msg))
        end subroutine check

    end subroutine read_input

end module nj_input
