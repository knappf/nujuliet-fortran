! Parameter structures ... counterparts of the Base.@kwdef structs of Import_Structures.jl with the same defaults ...
module nj_params
    use nj_kinds
    use nj_utils, only : jstr, istr
    implicit none
    private

    integer, parameter, public :: Str_Len = 1024

    type, public :: interaction_params_t
        character(len=Str_Len) :: NN_File = ''
        character(len=Str_Len) :: NNN_File = ''
        real(dp) :: hw = 0.0_dp
        integer :: Nmax = 0
        integer :: N2max = 0
        integer :: N3max = 0
        real(dp) :: cP2N = 1.0_dp
        real(dp) :: cP3N = 1.0_dp
        real(dp) :: cRes = 1.0_dp
        real(dp) :: pG2N = 0.0_dp
        real(dp) :: nG2N = 0.0_dp
    end type interaction_params_t

    ! LEGACY (removed in cleanup stage) ... HF-EFA is gone in NuJuliet dev (HF_orbitals_make & HF.ESH replace it) ...
    type, public :: hf_efa_params_t
        integer :: A = 0
        integer :: Z = 0
        logical :: BMF = .true.
    end type hf_efa_params_t

    ! Julia HF_MBPT_Parameters ... Export (correlated reference io/<Path>/HF-MBPT/bin/*), TBCF (2-body correlation
        ! function), ScH (polarised Hamiltonian) ...
    type, public :: hf_mbpt_params_t
        logical :: Export = .false.
        logical :: TBCF = .false.
        logical :: ScH = .false.
    end type hf_mbpt_params_t

    ! Julia HF_Parameters ... V_HF (export of the non-local HF radial potential) & ESH (enforced shell closure) are new,
        ! HRF has no reader left in dev (accepted, no effect) ...
    type, public :: hf_params_t
        real(dp) :: Tol = 1.0e-7_dp
        integer :: Imax = 100
        logical :: V_HF = .false.
        logical :: BMF = .true.
        logical :: HRF = .false.
        logical :: ESH = .false.
        type(hf_mbpt_params_t) :: MBPT
        type(hf_efa_params_t) :: EFA            ! LEGACY (removed in cleanup stage)
    end type hf_params_t

    ! Julia TDA_Parameters ... output folder io/<Path>/<Name>, reference io/<Path>/<Reference> ...
    type, public :: tda_params_t
        character(len=64) :: Name = 'TDA'
        character(len=64) :: Reference = 'HF'
        logical :: Reference_TBCF = .false.
        logical :: Ortho = .true.
        real(dp) :: Tol_hol = 0.7_dp
        real(dp) :: Tol_par = 0.1_dp
    end type tda_params_t

    ! Julia RPA_Parameters ... Transitions = 'Exact' or 'Linear', Tol & Imax of the ScRPA iteration, DEV developer switch ...
    type, public :: rpa_params_t
        character(len=64) :: Name = 'RPA'
        character(len=64) :: Reference = 'HF'
        logical :: Reference_TBCF = .false.
        logical :: RPA_Orb = .false.
        logical :: Ortho = .true.
        character(len=16) :: Transitions = 'Linear'
        real(dp) :: Tol = 1.0e-5_dp
        real(dp) :: Tol_hol = 0.7_dp
        real(dp) :: Tol_par = 0.1_dp
        integer :: Imax = 150
        logical :: OBDM_diagonal = .false.
        real(dp) :: OBDM_alpha = 0.6_dp
        logical :: ScH = .true.
        logical :: TBCF = .false.
        real(dp) :: TBCF_alpha = 0.8_dp
        real(dp) :: TBCF_beta = 5.0_dp
        logical :: DEV = .false.
    end type rpa_params_t

    ! LEGACY (removed in cleanup stage) ... HF-RRPA is gone in NuJuliet dev (HF_ScRPA replaces it) ...
    type, public :: rrpa_params_t
        logical :: Ortho = .true.
        real(dp) :: Tol = 1.0e-7_dp
        integer :: Imax = 150
        logical :: ScV3N = .true.
        logical :: ScOBH = .true.
        logical :: dOBDM = .false.
    end type rrpa_params_t

    type, public :: bcs_params_t
        logical :: ScBCS = .false.
        real(dp) :: Tol = 1.0e-7_dp
        integer :: Imax = 500
        real(dp) :: q = 0.05_dp
        integer :: A0 = 0
        integer :: Z0 = 0
        real(dp) :: pD0 = 0.5_dp
        real(dp) :: nD0 = 0.5_dp
        logical :: BMF = .true.
        logical :: HRF = .false.
    end type bcs_params_t

    type, public :: hfb_params_t
        real(dp) :: Tol = 1.0e-7_dp
        integer :: Imax = 150
        character(len=32) :: Pairing = 'Full'
        real(dp) :: dLmax = 0.5_dp
        real(dp) :: pL0 = -3.0_dp
        real(dp) :: nL0 = -3.0_dp
        real(dp) :: pa0 = 0.25_dp
        real(dp) :: na0 = 0.25_dp
        integer :: Broy_m = 8
        real(dp) :: Broy_Amax = 0.8_dp           ! dev default (0.9 before), Broy_Qmax & Broy_Bmax removed in dev
        real(dp) :: Broy_Tmax = 0.9_dp
        logical :: LNT = .false.
        logical :: LNRes = .true.
        logical :: BMF = .true.
        logical :: HRF = .false.
    end type hfb_params_t

    type, public :: qtda_params_t
        logical :: Ortho = .true.
        character(len=64) :: Reference = 'HFB'
    end type qtda_params_t

    type, public :: qrpa_params_t
        logical :: Ortho = .true.
        character(len=64) :: Reference = 'HFB'
    end type qrpa_params_t

    ! Julia pnQTDA_Parameters ... input folder io/<Path>/<Reference> of qpH1b.bin & qpH2b.bin ...
    type, public :: pnqtda_params_t
        character(len=64) :: Reference = 'HFB'
    end type pnqtda_params_t

    ! EMPM restricted particle-hole space ... windows (first, last) of the HF levels of each isospin ordered by
        ! energy (1 = deepest level), as in the input_tda_coup.dat of the original EMPM codes ...
        ! A zero bound means the full range ... holes 1 .. N_occupied, particles N_occupied+1 .. a_max ...
        ! 2-phonon step ... truncation E_alpha + E_lambda <= E2max, the J^P channels JP2 (e.g. '1-', '2+'; none = all) &
        ! the memory for the A & D rows Mem_GB (above it D rows are recomputed in blocks) ... E1max cut E_lambda <= E1max &
        ! Legacy_Order (phonon numbering of the original Tda) reproduce the truncated bases of the original eqm2 runs ...
        ! Dry_Run ... only the 2-phonon bases (dimensions, memory, Basis_*.dat), no A, D & eigenproblem ...
        ! Spur_Rows = 'all' or 'lambda' (spurious rows only with lambda = spurious phonon) & Legacy_Spur (spurious subspace
        ! selected as in the original eqm.f90) ... both for one-to-one comparisons with the original runs ... 'none' no
        ! spurious rows, the spurious subspace only from the D elements of the spurious states (smaller H & D) ...
        ! Full_H ... couplings to the 0- & 1-phonon subspaces & diagonalization of the full Hamiltonian of each channel ...
        ! D_Spectrum ... eigenvalue spectrum of the whole truncated D in the log (dense n^3 diagnostic, off by default) ...
        ! Cache_GB ... memory for the cached V_l & rho~_l of a block of lambdas (in addition to Mem_GB), more = fewer rho_a ...
        ! 3-phonon step (EMPM_3PH) ... truncated rows E_a + E_lambda <= E3max (a 2-phonon, lambda 1-phonon state), the J^P
        ! channels JP3 (none = all) ... 4-phonon step (EMPM_4PH) ... rows E_c + E_lambda <= E4max, channels JP4 ...
        ! full Hamiltonian with up to Full_N phonons (EMPM_FULL), channels JPF ... Save_X: EMPM_2PH writes the 2-phonon states
        ! (basis, E, C & X) to Bin/EMPM_S2_<JP>.bin for EMPM_3PH (large for big spaces, off by default) ...
        ! Ind_First ... 3- & 4-phonon steps: D of the truncated rows first, its pivoted Cholesky selects the independent rows,
        ! A & H = A D^T then only on them (same results, H & A D^T check on the independent block) ... .false. = H on all rows ...
        ! States_Dir ... folder of the n-phonon state files EMPM_S<n>_<JP>.bin (default Bin of the output folder), a leading
        ! environment variable is expanded (e.g. '$SCRATCHDIR/C12' ... local disk of a PBS job, gone after the job) ...
        ! Ind_Greedy ... with Ind_First: the independent rows selected block-wise (the rows of the (n-1)-phonon states a in a
        ! fixed random order, each block tested against the rows already selected by a pivoted Cholesky of its Schur
        ! complement in D) until their number reaches the dimension of the npnh space of the channel ... D of the visited rows
        ! only, never the whole N_tr x N_tr block ...
        ! CM_Exact ... EMPM_FULL keeps only the states with b_cm psi = 0 (exact CM condition, null space of the matrix of
        ! b_cm = HO lowering part of the CM coordinate between all EMPM states, requires Ortho = .false.): eigenvalues of
        ! N_cm = b+_cm b_cm below CM_Tol ... CM_Dump writes the b_cm blocks to Bin/EMPM_CMB_<JP>_<J'P'>.bin ...
        ! CM_Variants ... with CM_Exact also the weaker conditions, per truncation n: 'space' P_n b_cm psi = 0 (no component
        ! of b_cm psi in the <= n phonon space) & 'proj' P_(n-1) b_cm psi = 0 (physical = complement of the spurious space
        ! P_n b+_cm (<= n-1 phonon space), the generalization of the SVD treatment with the full operator) ... spectra &
        ! exact <N_cm> of their states in Bin/EMPM_FULL<n>_<JP>_E_<variant>.bin & _Ncm_<variant>.bin ...
    type, public :: empm_params_t
        integer :: pHole(2) = 0
        integer :: pPart(2) = 0
        integer :: nHole(2) = 0
        integer :: nPart(2) = 0
        logical :: Ortho = .true.
        real(dp) :: E2max = 1.0e8_dp
        character(len=8) :: JP2(64) = ''
        real(dp) :: Mem_GB = 100.0_dp
        real(dp) :: E1max = 1.0e8_dp
        logical :: Legacy_Order = .false.
        logical :: Dry_Run = .false.
        character(len=8) :: Spur_Rows = 'all'
        logical :: Legacy_Spur = .false.
        logical :: Full_H = .false.
        logical :: D_Spectrum = .false.
        real(dp) :: Cache_GB = 50.0_dp
        real(dp) :: E3max = 1.0e8_dp
        character(len=8) :: JP3(64) = ''
        real(dp) :: E4max = 1.0e8_dp
        character(len=8) :: JP4(64) = ''
        integer :: Full_N = 4
        character(len=8) :: JPF(64) = ''
        logical :: Save_X = .false.
        logical :: Ind_First = .true.
        character(len=512) :: States_Dir = ''
        logical :: CM_Exact = .false.
        real(dp) :: CM_Tol = 1.0e-4_dp
        logical :: CM_Dump = .false.
        logical :: Ind_Greedy = .false.
        logical :: CM_Variants = .false.
        ! CM_Keep (Ortho = .true. sets, EMPM_FULL) ... the cm phonon & the spurious n-phonon states stay in the basis with
            ! their couplings (R, Rc), so that CM_Exact / CM_Variants act on the complete sets in the CM-clean phonon basis ...
        logical :: CM_Keep = .false.
        ! Full_Diag (EMPM_FULL) ... which diagonalization of the full Hamiltonian per truncation n = 1 .. Full_N:
            ! 'unproj' (default) H on the complete <= n phonon basis, no CM projection -> EMPM/Full/<n>Phonon/Full_unproj_<JP>.dat
            !          (with Ortho and without CM_Keep the basis has no cm phonon & spurious states: this is the SVD result,
            !          written as Full_svd_<JP>.dat)
            ! 'proj'   only the CM projection P_(n-1) b_cm psi = 0 (orthogonal to the spurious space P_n b+_cm (<= n-1 phonons),
            !          formerly 'range'; complete sets: Ortho = .false. or Ortho + CM_Keep) -> Full_proj_<JP>.dat
            ! 'both'   both ... all variants have the same file structure: i, E (MeV, relative to E_HF), exact <N_cm> (-1 if not
            !          computed), phys = <N_cm> < 0.5, 0- .. 4-phonon weights; binaries Bin/EMPM_FULL<n>_<JP>_E_<lab>.bin,
            !          _Ncm_<lab>.bin (& _Z_<lab>.bin for unproj / proj), lowest states in Summary_<lab>.dat ...
        ! CM_Ncm ... exact <N_cm> of all states (complete N_cm matrix incl. the top-level completion, Full_N <= 3), computed
            ! whenever the sets are complete; .false. skips it (proj alone needs no b_cm within the top level, much cheaper);
            ! with CM_Keep also the SVD states (Full_svd_<JP>.dat) ... CM_Exact / CM_Variants ('exact' & 'space') are
            ! diagnostics only: in shell-truncated spaces they leave almost no states ...
        character(len=8) :: Full_Diag = 'unproj'
        ! Full_Svd (Ortho + CM_Keep) ... also the SVD result: H on the Ortho-physical basis states of each level (without the
            ! cm phonon & the spurious n-phonon states) -> Full_svd_<JP>.dat, with its exact <N_cm> if computed ...
        logical :: Full_Svd = .false.
        logical :: CM_Ncm = .true.
        ! CM_Beta (MeV, > 0) ... Lawson term in EMPM_FULL: H + CM_Beta N_cm diagonalized with the complete N_cm matrix
            ! ('lawson', and with Full_Diag 'proj'/'both' also on the proj space, 'projL'), E = <H> & <N_cm> to Full_<lab>_<JP>.dat
            ! (complete sets needed, the N_cm matrix is built as for CM_Ncm) ...
        real(dp) :: CM_Beta = 0.0_dp
        ! D_Disk_GB (Ind_Greedy) ... the D rows of the selected independent rows (N_full x N_ind, needed for H = A D^T & the
            ! export X = D^T C) kept in a scratch file next to the state files instead of memory when larger than D_Disk_GB
            ! (< 0 = automatic: larger than 0.3 Mem_GB) ... streamed back once per block of rows & in row blocks for X ...
        real(dp) :: D_Disk_GB = -1.0_dp
        ! Dens_Share (Ind_Greedy) ... the greedy selections of all channels of the step in one pass over the (n-1)-phonon
            ! states a, the densities of each a computed once for all channels (memory of the blocks split among them) ...
        logical :: Dens_Share = .false.
        ! EMPM_SF ... Lanczos strength function B(E1; 0+_gs -> 1-) with up to 3 phonons in the redundant 3-phonon basis:
            ! SF_Steps Lanczos steps, SF_Op = 'iv' / 'is' / 'ph', SF_Reorth full reorthogonalization ...
        integer :: SF_Steps = 100
        character(len=8) :: SF_Op = 'iv'
        logical :: SF_Reorth = .true.
        ! SF_Mode ... 'full' the whole redundant 3-phonon basis (A^T & D products, complete sets), 'ind' its independent
            ! rows (as EMPM_3PH, also with E2max / E3max) ...
        character(len=8) :: SF_Mode = 'full'
        ! SF_Qrows (SF_Mode = 'ind') ... keep all rows |0+ 2-phonon state x 1- phonon> regardless of E3max ...
        logical :: SF_Qrows = .false.
    end type empm_params_t

    ! No-core shell model in the m-scheme HO basis (action NCSM) ... Nmax_List the many-body truncations (HO quanta above
        ! the lowest configuration, -1 terminates the list; the parity is that of N0 + Nmax), Hamiltonian 'NN' (T_int + V_NN in
        ! the HO basis, no HF needed) or 'EMPM' (the normal-ordered EMPM Hamiltonian E_HF + e :a+a: + 1/4 V2B_HF :a+a+aa: of
        ! the HF run in io/<Path>, transformed back to the HO basis), M2 = 2M (-1 = A mod 2), N_Eig lowest states by Lanczos
        ! with full reorthogonalization (at most Lanczos_Max steps, residual < Lanczos_Tol MeV), Beta_CM the Lawson term
        ! H + Beta_CM N_cm (MeV per CM quantum, 0 = none) ...
    type, public :: ncsm_params_t
        integer :: Nmax_List(16) = -1
        character(len=8) :: Hamiltonian = 'EMPM'
        integer :: M2 = -1
        integer :: N_Eig = 10
        integer :: Lanczos_Max = 400
        real(dp) :: Lanczos_Tol = 1.0e-7_dp
        real(dp) :: Beta_CM = 0.0_dp
        character(len=8) :: Matrix = 'auto'      ! 'stored', 'otf' (on the fly) or 'auto' (stored if within Mem_GB)
        logical :: Single = .false.              ! stored values in single precision
        integer :: Lanczos_Basis = 0             ! thick-restart basis size (0: max(2 N_Eig + 50, 100))
        real(dp) :: Matrix_GB = 0.0_dp           ! limit of the stored matrix (0: 70% of MemAvailable)
        character(len=8) :: Trunc = 'Nmax'       ! 'Nmax' (HO quanta) or 'npnh' (all m-states of the model space e <= Nmax, at most
                                                 ! Nmax_List(i) particles outside the lowest A m-states, both parities)
        character(len=8) :: Basis = 'HO'         ! 'HO' or 'HF' (EMPM Hamiltonian kept in the HF orbitals, N_cm transformed; npnh
                                                 ! then counted w.r.t. the HF reference as in the EMPM)
        real(dp) :: Ecut_npnh(16) = -1.0_dp      ! 'npnh' & Basis 'HF' only: energy cuts (MeV) on the unperturbed HF excitation
                                                 ! energy sum(e_particles) - sum(e_holes) of the proton x neutron configurations,
                                                 ! Ecut_npnh(n) for the n-particle n-hole part (n = 1, 2, ...), e.g.
                                                 ! Ecut_npnh = 80, 70, 60, 50 (1p1h .. 4p4h); n beyond the last value given
                                                 ! repeat it (Ecut_npnh = 60: all n), a value < 0: no cut of that part
        logical :: Transitions = .false.         ! 'npnh', M2 = 0: B(EL; g.s. -> f), L = 0, 2 (parity +) & 1, 3 (parity -) from
                                                 ! the lowest state of parity + (the 0+ ground state)
        integer :: Dense_Max = 0                 ! full diagonalization (all states, LAPACK dsyevd) when the dimension <= Dense_Max
        logical :: Window = .false.              ! 'npnh', Basis 'HF': single-particle window of &empm (pHole, pPart, nHole, nPart, HF levels
                                                 ! by energy) ... holes below the window frozen (core), particles above dropped
        integer :: Parity = -1                   ! 'npnh': parities run, -1 both, 0 (+) or 1 (-) only
        character(len=8) :: Tr_From = 'gs'       ! Transitions from 'gs' (lowest state of parity +, both parities run) or 'HF'
                                                 ! (the reference determinant, E = E_HF)
        integer :: Tr_Lanczos = 0                ! > 0: strength-function Lanczos (iterations) from Q|g> for B_ph, B_is, B_iv
        integer :: Tr_Gs_Nmax = 0                ! Tr_From 'gs', > 0: the g.s. of the <= Tr_Gs_Nmax p-h run (Nmax_List) kept for the
                                                 ! larger runs (mapped by its SD occupations): g.s. <= n p n h, final states <= n' p n' h
        character(len=32) :: Tr_Ops = ''         ! Tr_Lanczos: strength functions computed, '' = all, or a list of E<L><ph|is|iv>,
                                                 ! e.g. 'E1iv' or 'E1iv,E1is' (the others skipped: fewer products of H)
        character(len=8) :: Tr_Vec = 'double'    ! Tr_Lanczos: precision of the stored Lanczos vectors, 'double' or 'single' (half the
                                                 ! memory; the current vectors w, products & coefficients stay double)
        logical :: Tr_Sym = .true.               ! Tr_Lanczos (2M = 0, Matrix 'otf', no eigenstates): rotation symmetry R_y(pi), half of the rows & of the species matrices
        logical :: Tr_Reorth = .true.            ! Tr_Lanczos: full reorthogonalization (all vectors kept, n x Tr_Lanczos) or none (3
                                                 ! vectors: plain Lanczos, moments & folded strength exact, ghost copies of converged poles)
        integer :: Tr_Ckpt = 0                   ! Tr_Lanczos: > 0 checkpoint every Tr_Ckpt steps (the Lanczos vectors appended to
                                                 ! <file>_vecs.bin every step, alpha, beta & the two current vectors to <file>_ckpt.bin)
                                                 ! & a run with a valid checkpoint resumes from its last step ... 0: none
    end type ncsm_params_t

    type, public :: calculation_params_t
        integer :: A = 0
        integer :: Z = 0
        real(dp) :: hw = 0.0_dp
        integer :: Nmax = 0
        integer :: N2max = 0
        integer :: N3max = 0
        character(len=32) :: CMS = 'CMS1+2B'
        ! CM_Lawson (MeV, 0 = none) ... beta N_cm added to V_NN when the interaction is read (V2b_read), N_cm = 1/A (sum_i N_i +
            ! sum_(i/=j) b+(i).b(j)) with the 1-body part folded into the 2-body one by 1/(A-1) ... HF, EMPM & NCSM 'EMPM' then
            ! work with H + beta N_cm (all energies include beta <N_cm>) ...
        real(dp) :: CM_Lawson = 0.0_dp
        character(len=Str_Len) :: Path = ''
        type(hf_params_t) :: HF
        type(tda_params_t) :: TDA
        type(rpa_params_t) :: RPA
        type(rrpa_params_t) :: RRPA              ! LEGACY (removed in cleanup stage)
        type(bcs_params_t) :: BCS
        type(hfb_params_t) :: HFB
        type(qtda_params_t) :: QTDA
        type(qrpa_params_t) :: QRPA
        type(pnqtda_params_t) :: pnQTDA
        type(empm_params_t) :: EMPM
        type(ncsm_params_t) :: NCSM
    end type calculation_params_t

    type, public :: params_t
        type(interaction_params_t) :: Int
        type(calculation_params_t) :: Calc
    end type params_t

    public :: default_path, io_path, io_dir

contains

    ! Default output folder name ... "A16_Z8_hw16.0_Nmax3_N2max6_N3max9_CMS1+2B" ...
    function default_path(Calc) result(Path)
        type(calculation_params_t), intent(in) :: Calc
        character(len=:), allocatable :: Path
        Path = 'A' // istr(Calc%A) // '_Z' // istr(Calc%Z) // '_hw' // jstr(Calc%hw) // &
               '_Nmax' // istr(Calc%Nmax) // '_N2max' // istr(Calc%N2max) // '_N3max' // istr(Calc%N3max) // &
               '_' // trim(Calc%CMS)
    end function default_path

    ! "io/" * Params.Calc.Path * Suffix ... the output root is lower-case io/ since NuJuliet dev (IO/ before) ...
    function io_path(Params,Suffix) result(Path)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Suffix
        character(len=:), allocatable :: Path
        Path = 'io/' // trim(Params%Calc%Path) // Suffix
    end function io_path

    ! "io/" * Params.Calc.Path * "/" * Sub ... per-solver folder, e.g. io_dir(Params,Params%Calc%RPA%Name) ...
    function io_dir(Params,Sub) result(Path)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Sub
        character(len=:), allocatable :: Path
        Path = 'io/' // trim(Params%Calc%Path) // '/' // trim(Sub)
    end function io_dir

end module nj_params
