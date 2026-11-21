! No-core shell model (NCSM) in the m-scheme HO basis ... a reference for the EMPM with the same matrix elements ...
    ! Basis: products |P>|N> of proton & neutron Slater determinants of HO m-states with total quanta <= N0 + Nmax &
        ! (-1)^(N0+Nmax) parity, M = M_P + M_N fixed ... the many-body vector is ordered by sectors (proton SD, neutron block
        ! (M_N, q_N)), neutron index fastest ...
    ! Operators in the plain form C0 + sum u a+a + 1/4 sum v_A a+a+aa (like particles, antisymmetrized m-scheme v_A) +
        ! sum v_pn a+_a a_c a+_b a_d (pn, direct), the m-scheme 2-body tables stored per pair block (2M, parity) of the pairs
        ! within the quanta cut (the largest pair quanta of the basis) ... v_m(ab,cd) = sum_J <ja ma jb mb|J M>
        ! <jc mc jd md|J M> V^J(ab,cd) from the NuJuliet J-coupled storage (the convention of tests/empm_fci_check.py) ...
    ! Hamiltonians ...
        ! 'NN'    T_int + V_NN in the HO basis: (1 - 1/A) T + V_NN incl. the 2-body CM correction (CMS1+2B, V2b_read) ...
        ! 'EMPM'  H = E_HF + sum e_a :a+a: + 1/4 sum H2b :a+a+aa: (H1b diagonal, residual 3-body dropped, as EMPM &
        !         tests/fci_sparse_cm.py) -> plain form u = e - V rho, C0 = -sum e rho + 1/2 V rho rho in the HF basis,
        !         transformed back to HO with C (u -> C u C^T, H2b by O2b_transformation with C^T) ... the HF-like reference
        !         io/<Path>/<&empm Reference>/bin/{orb,H1b,H2b,C,rho}.bin (stage 5; H1b diagonal, integer occupations) ...
    ! Observables: J^2 = sum j_i^2 + sum_(i/=j) j_i.j_j & N_cm = b+_cm.b_cm = 1/A (sum_i N_i + sum_(i/=j) b+(i).b(j)),
        ! b_mu = sqrt(2) (r_mu/b)_(N -> N-1) the HO lowering operator (radial wave functions positive at the origin, (l s) j
        ! coupling, as tests/empm_cm_check.py) ... Lawson term H + Beta_CM N_cm optional ...
    ! Engine (v2): one operator H + Beta_CM N_cm (2-body tables combined), its columns generated from the species matrices
        ! & the pn one-body jumps ... Matrix = 'stored' (lower half by columns, double or single precision values, built in
        ! two passes) or 'otf' (on the fly, no matrix), 'auto' = stored if it fits in Matrix_GB (default 70% of MemAvailable)
        ! ... thick-restart Lanczos (basis Lanczos_Basis vectors, full reorthogonalization) ... <J^2> & <N_cm> of the final
        ! states on the fly in blocks of 16 vectors, E = <H> = Ritz value of H + Beta N_cm - Beta <N_cm> ...
module nj_ncsm
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_hashmap
    use nj_wigner, only : wigner_init_float, fCG
    use nj_orbitals, only : orbitals_make, read_h1b_diag
    use nj_o1b, only : O1b_import, T1b
    use nj_o2b, only : O2b_import, O2b_transformation, O2b_pp, O2b_pn, O2b_nn
    use nj_v2b, only : V2b_read
    use nj_linalg, only : eig_sym
    use nj_tr1b, only : tr1b_t, Tr1b_initialize, Tr1b_transformation
    use nj_empm_space, only : empm_space_t, EMPM_space_make, empm_ref_path, EMPM_reference_orbitals
    use omp_lib
    implicit none
    private

    public :: NCSM

    integer, parameter :: i1 = selected_int_kind(2)
    integer, parameter :: f4 = kind(1.0)

    interface
        subroutine dgemv(Trans,m,n,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
            import :: dp
            character :: Trans
            integer :: m, n, Lda, Incx, Incy
            real(dp) :: Alpha, Beta, A(Lda,*), X(*), Y(*)
        end subroutine dgemv
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
    end interface

    abstract interface
        subroutine apply_i(x,y)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine apply_i
    end interface

    ! HO m-states k ... orbital a, 2m, quanta e, l, 2j, radial n ...
    ! One-body jump types (kc, ka) = a+_kc a_ka in buckets b = 2 (dm + dmax/2) + pt + 1, dm = (m2_kc - m2_ka)/2, pt the
        ! parity of l_kc + l_ka ... Bk(kc,ka) the bucket, Ix(kc,ka) the index within it, Tkc(b)%v, Tka(b)%v its types,
        ! Bp(b) the partner bucket of the pn coupling (dm -> -dm, same pt), b0 = dmax + 1 the diagonal bucket (dm = 0, pt 0) ...
    type :: sp_t
        integer :: ns = 0
        integer, allocatable :: a(:), m2(:), e(:), l(:), j2(:), n(:)
        ! binomials C(n,k), n = 0..ns, k = 0..KBIN, -1 where > huge(i8)/KBIN ... SD keys = combinatorial rank
        integer(i8), allocatable :: Bin(:,:)
        integer :: nbk = 0, dmax = 0, b0 = 0
        integer, allocatable :: Bk(:,:), Ix(:,:), Nty(:), Bp(:)
        type(ivec_t), allocatable :: Tkc(:), Tka(:)
    end type sp_t
    integer, parameter :: KBIN = 32
    integer, parameter :: i2 = selected_int_kind(4)
    integer, parameter :: IXOFF = 32768          ! Jx = Ix - IXOFF (2 bytes, up to 65535 types per bucket)

    ! Pairs of a kind (1 pp k1 < k2, 2 nn k1 < k2, 3 pn k1 proton, k2 neutron) within the quanta cut, by blocks (2M, P) ...
    type :: pairs_t
        integer :: np = 0, nb = 0
        integer, allocatable :: k1(:), k2(:), blk(:), loc(:), id(:,:)
        integer, allocatable :: bM2(:), bP(:)
        type(ivec_t), allocatable :: mem(:)
    end type pairs_t

    type :: tabs_t
        type(rmat_t), allocatable :: B(:)
    end type tabs_t

    ! Operator in the plain form ... u(k,k',t) & the 2-body tables per pair kind ... jump-type tables (pn_tables): W(b)%v
        ! (ixn, ixp) = v_pn(kc kc', ka ka') for the proton type ixp of bucket b & the neutron type ixn of Bp(b), uT(ix,t) the
        ! one-body u of the types of the diagonal bucket ...
    type :: op_t
        real(dp) :: C0 = 0.0_dp
        real(dp), allocatable :: u(:,:,:)
        type(tabs_t) :: T(3)
        type(rmat_t), allocatable :: W(:)
        real(dp), allocatable :: uT(:,:)
        integer, allocatable :: Bp(:)
    end type op_t

    ! Slater determinants of one species ... occupied m-states (ascending), quanta above N0, 2M, the one-body jumps
        ! a+_kc a_ka |s> = sg |s'> (incl. kc = ka) in the type buckets b of SPS: jumps of s with bucket b at
        ! Jo(b,s) .. Jo(b+1,s)-1, Jt = sg s' (the sign in the sign of the target), Jx = Ix(kc,ka) - IXOFF ...
    type :: sd_t
        integer :: Nt = 0, nsd = 0, N0 = 0, nd = 0
        integer, allocatable :: occ(:,:), q(:), m2(:)
        ! h(s) particles outside the reference (the lowest Nt m-states), L(s) the sector label: q (Trunc 'Nmax') or 2h + mod(q,2) ('npnh') ...
        integer, allocatable :: h(:), L(:)
        ! ex(s) unperturbed HF excitation energy sum(e_occ) - sum(e_ref) (MeV) ... only with the energy cut (Ecut_npnh) ...
        real(dp), allocatable :: ex(:)
        type(hashmap_t) :: Map
        integer(i8), allocatable :: Jo(:,:)   ! 64-bit offsets (> 2^31 jumps, e.g. 16O 4p4h in 5 shells)
        integer, allocatable :: Jt(:)
        integer(i2), allocatable :: Jx(:)
        ! Hole-particle form (Trunc 'npnh'): reference m-states RefL (ranks Rr), particle m-states ParL (ranks Rp), inactive
            ! m-states (Rr = Rp = 0, e.g. outside the window or of the other species) never occupied ... key = Off(h) +
            ! rank(holes) C(npar,h) + rank(particles), ranks combinatorial (Ck(n,k) = C(n,k), k <= Hk) ... Pt(k) particle
            ! m-state, Ac(k) active ... otherwise (HP false) the reference is 1..Nt & the key the rank of all occupations ...
        logical :: HP = .false.
        integer :: nref = 0, npar = 0, Hk = 0
        integer, allocatable :: Rr(:), Rp(:), RefL(:), ParL(:)
        integer(i8), allocatable :: Ck(:,:), Off(:)
        logical, allocatable :: Pt(:), Ac(:)
        ! He(h) the last SD with at most h particles outside the reference (hole-particle form: SDs in ascending h) ...
            ! the jumps of a bucket & the species-matrix columns are in ascending targets, i.e. ascending h': loops over
            ! them end at the first target outside the npnh truncation (col_gen) ...
        integer, allocatable :: He(:)
    end type sd_t

    ! Jumps of a chunk of SDs (sd_make): signed targets, types & buckets ...
    type :: jbuf_t
        integer, allocatable :: t(:)
        integer(i2), allocatable :: x(:), d(:)
    end type jbuf_t

    ! Species matrix column: partners r (lower half by the ranks) & values v (double) or vs (single) ...
    type :: col_t
        integer, allocatable :: r(:)
        real(dp), allocatable :: v(:)
        real(f4), allocatable :: vs(:)
    end type col_t

    ! Many-body basis ... state c = (proton SD sp(c), neutron SD sn(c)), St(sp,qn) the first index of the sector of sp &
        ! the neutron block (M - M_P, qn) (-1 if absent), Lc(sn) the position of sn in its block ... with the energy cut
        ! the blocks are in ascending neutron excitation energy & the sector (sp,qn) holds the first Nc(sp,qn) of them
        ! (E_P + E_N <= Ecut_npnh(h_P + h_N)) ... dim_all the number of products of the kept SDs before the pair cut ...
        ! Rn(sn) the rank of the neutron SD in the order of the sectors (blocks by L, then by position): for one sp the state
        ! index grows with Rn ... states are in ascending sp, so r >= c iff sp_r > sp_c or (sp_r = sp_c & Rn(sn_r) >= Rn(sn_c)) ...
    type :: mb_t
        integer :: dim = 0, M2 = 0, Nmax = 0
        integer(i8) :: dim_all = 0
        integer, allocatable :: sp(:), sn(:), St(:,:), Lc(:), Nc(:,:), Rn(:)
    end type mb_t

    ! Stored lower half ... per thread a list of segments (SEG entries each: rows r, values vs or vd), the columns of the
        ! thread (chunks of CHK columns, chunk k to thread mod(k,nt) + 1) in order, a column entirely within a segment
        ! (a new segment is started when it does not fit), Cnt(c) entries of column c ...
    ! Rows of a column ascending, as gaps: d = r - r_previous (r_previous = c - 1 at the start) - GOFF in 2 bytes for
        ! gaps 1 .. 65535, d = -GOFF marks a larger gap: the row is the next entry of the escape list e (ESEG per segment)
        ! ... Cnt(c) < 0: the column starts a new segment ...
    integer, parameter :: SEG = 2**22, ESEG = SEG/16, CHK = 1024, GOFF = 32768
    type :: seg_t
        integer :: n = 0
        integer(i2), allocatable :: d(:)
        integer, allocatable :: e(:)
        real(dp), allocatable :: vd(:)
        real(f4), allocatable :: vs(:)
    end type seg_t
    type :: tstore_t
        integer :: ns = 0
        type(seg_t), allocatable :: S(:)
    end type tstore_t

    ! Operator engine ... species matrices Mp, Mn (lower halves) & optionally the stored lower half with the thread
        ! buffers Yb of the symmetric products ...
    ! Truncation & basis of the current run ... Trunc 'Nmax' or 'npnh', Uhf(a,b,t) the HF coefficients (Basis 'HF') ...
    character(len=8), save :: Trunc = 'Nmax'
    real(dp), allocatable, save :: Uhf(:,:,:)
    ! Energy cut of the current run ('npnh', Basis 'HF') ... Ecn(n) the cut of the npnh configurations (huge = none) ...
    logical, save :: UseCut = .false.
    real(dp), allocatable, save :: Ecn(:)
    ! HF single-particle energies e(a,t) (diagonal of the reference H1b, set by hamiltonian_empm) ...
    real(dp), allocatable, save :: EpsHF(:,:)

    type :: eng_t
        type(col_t), allocatable :: Mp(:), Mn(:)
        logical :: Stored = .false., Single = .false.
        integer(i8) :: nnz = 0
        integer :: nt = 0
        integer, allocatable :: Cnt(:), ChT(:), ChS(:), ChP(:), ChE(:)   ! chunk k: owner thread, segment & positions (entries, escapes) at its start
        type(tstore_t), allocatable :: Ts(:)
        real(dp), allocatable :: Yb(:,:)
        ! Sector product (eng_mv_sec, no per-thread buffers) ... pass B tasks over tiles of the neutron blocks: group g =
            ! (proton 2M, neutron label Ln) with its proton SDs Gs(Go(g)+1 .. Go(g+1)) (sectors (sp,Ln) in the basis), task k:
            ! group Tg(k), positions Tl0(k)..Tl1(k) in the neutron block, proton SDs Gs(Tk0(k)..Tk1(k)), Ta(k) atomic updates
            ! (a tile split over several tasks) ...
        logical :: SecReady = .false.
        integer :: ntask = 0
        integer, allocatable :: GLn(:), Tg(:), Tl0(:), Tl1(:), Gs(:)
        integer(i8), allocatable :: Go(:), Tk0(:), Tk1(:)
        logical, allocatable :: Ta(:)
        ! pass A tasks: rows Ar0(k)..Ar1(k) of the proton SD As(k), the whole sector of sp (one task, plain updates) or a
            ! chunk of a big sector (> ABIG rows: several tasks, atomic neutron updates, Aa(k)) ...
        integer :: natask = 0
        integer, allocatable :: As(:), Ar0(:), Ar1(:)
        ! Pk(sn) = Lc(sn) 32 + L(sn): one random read for the row of a neutron partner ...
        integer, allocatable :: Pk(:)
        logical, allocatable :: Aa(:)
        ! Rotation symmetry R_y(pi) (2M = 0): R |s> = sg |s'>, TP(sp) = sg s' of the proton SDs, TN(sn) of the neutron SDs;
            ! a vector with R x = Eta x (Eta = (-1)^J): x(c') = Eta sg_P sg_N x(c) for c = (sp,sn), c' = (sp',sn') ... the
            ! products are computed for the rows with M_P > 0 or (M_P = 0 & sp <= sp') & copied to the others (sym_fill), the
            ! species matrices only for M_P >= 0 & M_N <= 0 ...
        ! Reduced blocks of the SDs with h = Hk = Nmax (paired only with the reference of the other species, RefP / RefN,
            ! a self-conjugate SD), RedP / RedN: columns only of the representatives a <= a' with K(a,b) = M(a,b) +
            ! eps sg_b M(a,b') (b a representative, eps = Eta sg_ref) & the columns of h < Hk with the representative partners
            ! of h = Hk only: the products add the mirror term y(s') += eps sg_s sum_a M(s,a) x(a) ... Eta fixed at the build ...
        logical :: RedP = .false., RedN = .false.
        integer :: RefP = 0, RefN = 0
        logical :: Sym = .false.
        integer :: Eta = 0
        integer, allocatable :: TP(:), TN(:)
    end type eng_t
    ! timings of the passes of eng_mv_sec (accumulated) ...
    real(dp), save :: TmA = 0.0_dp, TmB = 0.0_dp
    ! parts of pass A (NJ_MV = prof: timed separately once per basis) ...
    logical, save :: DoN = .true., DoPN = .true., Cnt = .false.

    ! Context of the H product passed to lanczos_tr (apply_h) ... module-level pointers instead of an internal procedure
        ! (an internal procedure as an actual argument needs an executable-stack trampoline: SIGSEGV with ifx) ...
    type(sd_t), pointer, save :: cP => null(), cNn => null()
    type(pairs_t), pointer, save :: cPr(:) => null()
    type(mb_t), pointer, save :: cB => null()
    type(op_t), pointer, save :: cH => null()
    type(eng_t), pointer, save :: cEH => null()

    ! Transitions (&ncsm Transitions) ... QL(a,b,t,L) = <a||Q_L||b> (Edmonds, t = 1 p, 2 n, L = 0..3: r^2 Y_0, r Y_1, r^2 Y_2,
        ! r^3 Y_3 of nj_tr1b) in the orbitals of the run (HF: transformed with C), QH the same in the HO basis (convention
        ! check), the ground state Xg (lowest state of parity +, M = 0) in its basis Bg & its energy Eg, kept for the parity - run
        ! (the SDs of 'npnh' do not depend on the parity: Bg indexes the same proton & neutron SDs) ...
    real(dp), allocatable, save :: QL(:,:,:,:), QH(:,:,:,:)
    type(mb_t), save :: Bg
    real(dp), allocatable, save :: Xg(:)
    real(dp), save :: Eg = 0.0_dp, Jg = 0.0_dp, Bosc = 0.0_dp
    ! Tr_Gs_Nmax > 0: the g.s. of the <= Tr_Gs_Nmax p-h run by its SD occupations (Gp, Gn: proton & neutron occ of the
        ! components, Gx the amplitudes) ... mapped into the SDs of the larger runs (gs_remap) ... CurNmax the current run ...
    integer, allocatable, save :: Gp(:,:), Gn(:,:)
    real(dp), allocatable, save :: Gx(:)
    integer, save :: CurNmax = 0
    logical, save :: HfBasis = .false.
    ! Reference & active orbitals per species (a, t = 1 p, 2 n): the HF occupations (Basis 'HF') or the &empm window
        ! (Window: frozen core below the hole window, particles above the particle window dropped) ... unallocated: the
        ! lowest Nt m-states ...
    logical, allocatable, save :: RefO(:,:), ActO(:,:)
    logical, save :: UseWin = .false.

contains

    subroutine NCSM(Params)
        type(params_t), intent(in) :: Params
        type(sp_t) :: SPS
        type(op_t) :: H, J2op, Nop
        real(dp), allocatable :: uH(:,:,:)
        type(o2b_t) :: V
        type(orb2b_t) :: Orb_NN
        type(orb1b_t), allocatable :: Orb(:)
        real(dp) :: C0, E_HF
        integer :: iN, Nmax, emax_sp, Z, N, A, Par
        integer(i8) :: t0
        character(len=8) :: Ham

        Ham = Params%Calc%NCSM%Hamiltonian
        Trunc = Params%Calc%NCSM%Trunc
        if (Trunc /= 'Nmax' .and. Trunc /= 'npnh') call fatal('NCSM ... unknown Trunc ' // trim(Trunc) // ' (Nmax or npnh)')
        if (Params%Calc%NCSM%Basis /= 'HO' .and. (Params%Calc%NCSM%Basis /= 'HF' .or. Ham /= 'EMPM')) &
            call fatal('NCSM ... Basis ' // trim(Params%Calc%NCSM%Basis) // ' (HO, or HF with the EMPM Hamiltonian)')
        A = Params%Calc%A
        Z = Params%Calc%Z
        N = A - Z
        if (Z < 1 .or. N < 1) call fatal('NCSM ... needs Z >= 1 & N >= 1')
        call wigner_init_float(2*Params%Calc%Nmax + 3,'Jmax',6)
        write(*,'(/,a)') 'Starting the NCSM in the m-scheme HO basis ... Hamiltonian ' // trim(Ham) // ' ...'
        if (.not. is_dir(io_path(Params,''))) call make_dir(io_path(Params,''))
        if (.not. is_dir(io_path(Params,'/NCSM'))) call make_dir(io_path(Params,'/NCSM'))

        ! Hamiltonian in the J-coupled HO storage, u(a,b,t) on the orbitals & C0 ...
        ! orbitals ... the HF-like reference io/<Path>/<&empm Reference>/bin/orb.bin for the EMPM Hamiltonian or the HF basis
            ! (stage 5; LHO labels in the LHO slots, Calc.Nmax size), else the LHO orbitals of orbitals_make (Calc.Nmax size,
            ! dev default; T1b & V2b_read depend on the labels only) ...
        if (Ham == 'EMPM' .or. Params%Calc%NCSM%Basis == 'HF') then
            Orb = EMPM_reference_orbitals(Params)
        else
            Orb = orbitals_make(Params)
        end if
        UseWin = Params%Calc%NCSM%Window
        if (allocated(RefO)) deallocate(RefO)
        if (allocated(ActO)) deallocate(ActO)
        if (UseWin) call window_make(Params,Orb)
        call timer_start(t0)
        E_HF = 0.0_dp
        select case (trim(Ham))
        case ('NN')
            call hamiltonian_nn(Params,Orb,V,Orb_NN,uH,C0)
        case ('EMPM')
            call hamiltonian_empm(Params,Orb,V,Orb_NN,uH,C0,E_HF)
        case default
            call fatal('NCSM ... unknown Hamiltonian ' // trim(Ham) // ' (NN or EMPM)')
        end select
        call timer_report(t0,'Hamiltonian in the J-coupled HO basis')
        if (Params%Calc%NCSM%Transitions) call trans_init(Params,Orb)

        do iN = 1, size(Params%Calc%NCSM%Nmax_List)
            Nmax = Params%Calc%NCSM%Nmax_List(iN)
            if (Nmax < 0) exit
            if (Trunc == 'npnh') then
                ! all m-states of the model space, at most Nmax particles outside the reference, both parities ...
                if (UseWin) then
                    call sp_make(Orb,Params%Calc%Nmax,SPS,ActO(:,1) .or. ActO(:,2))
                else
                    call sp_make(Orb,Params%Calc%Nmax,SPS)
                end if
                do Par = 0, 1
                    if (Params%Calc%NCSM%Parity >= 0 .and. Par /= Params%Calc%NCSM%Parity) cycle
                    ! Tr_Gs_Nmax: the g.s. run parity + only, the larger runs the other parities (from the kept g.s.) ...
                    if (Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_From == 'gs' .and. Params%Calc%NCSM%Tr_Gs_Nmax > 0) then
                        if (Nmax == Params%Calc%NCSM%Tr_Gs_Nmax .and. Par == 1) cycle
                        if (Nmax /= Params%Calc%NCSM%Tr_Gs_Nmax .and. Par == 0 .and. Params%Calc%NCSM%Parity < 0) cycle
                    end if
                    write(*,'(/,a)') '=== NCSM ' // istr(Nmax) // 'p' // istr(Nmax) // 'h, parity ' // merge('-','+',Par == 1) // &
                                     ' ... ' // trim(Params%Calc%NCSM%Basis) // ' m-states with e <= ' // istr(Params%Calc%Nmax)
                    call run_nmax(Params,Orb,SPS,V,Orb_NN,uH,C0,E_HF,Nmax,Par)
                end do
                cycle
            end if
            emax_sp = min(Params%Calc%Nmax, max(top_e(Z),top_e(N)) + Nmax)
            write(*,'(/,a)') '=== NCSM Nmax = ' // istr(Nmax) // ' ... HO m-states with e <= ' // istr(emax_sp)
            if (emax_sp < max(top_e(Z),top_e(N)) + Nmax) write(*,'(a)') TAB // 'Warning ... the single-particle space ' // &
                'e <= ' // istr(Params%Calc%Nmax) // ' cuts the Nmax space (not Nmax-complete) ...'
            call sp_make(Orb,emax_sp,SPS)
            call run_nmax(Params,Orb,SPS,V,Orb_NN,uH,C0,E_HF,Nmax,-1)
        end do
        write(*,'(/,a,/)') 'NCSM calculation has finished ...'

    contains

        ! Highest quanta of the lowest configuration of n nucleons of one species ...
        integer function top_e(n)
            integer, intent(in) :: n
            integer :: e, k
            k = 0
            e = 0
            do
                k = k + (e + 1)*(e + 2)
                if (k >= n) exit
                e = e + 1
            end do
            top_e = e
        end function top_e

    end subroutine NCSM

    ! One Nmax ... basis, operators, engine of H + Beta N_cm, Lanczos, observables & output ...
    subroutine run_nmax(Params,Orb,SPS,V,Orb_NN,uH,C0,E_HF,Nmax,ParIn)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(sp_t), intent(in) :: SPS
        type(o2b_t), intent(in) :: V
        type(orb2b_t), intent(in) :: Orb_NN
        real(dp), intent(in) :: uH(:,:,:), C0, E_HF
        integer, intent(in) :: Nmax, ParIn
        type(sd_t), target :: P, Nn
        integer :: Qcap, Hcap
        type(pairs_t), target :: Pr(3)
        type(op_t), target :: H
        type(op_t) :: Ncm
        type(mb_t), target :: B
        type(mb_t) :: B2
        type(eng_t), target :: EH
        type(eng_t) :: EO
        real(dp), allocatable :: E(:), Zv(:,:), Jv(:), Nv(:), Ev(:), Res(:), Xt(:,:)
        integer :: Z, Nneu, A, M2, Epair(3), t, k, Ne, Par, u, Nconv
        integer(i8) :: t0
        real(dp) :: Beta, MemGB
        logical :: Fits, Single
        character(len=8) :: Mode
        character(len=:), allocatable :: Fn, Tag, Hdr
        real(dp), allocatable :: EpsP(:), EpsN(:), ExmP(:), ExmN(:)
        integer :: Lg, n, NtP, NtN, KpM(2), KpP
        logical, allocatable :: RfS(:,:), AcS(:,:)
        logical :: SfOnly

        A = Params%Calc%A
        Z = Params%Calc%Z
        Nneu = A - Z
        M2 = Params%Calc%NCSM%M2
        if (M2 < 0) M2 = mod(A,2)
        Beta = Params%Calc%NCSM%Beta_CM
        CurNmax = Nmax
        ! reference & active m-states per species (HF occupations or the window), else the lowest Z, N m-states ...
        allocate(RfS(SPS%ns,2), AcS(SPS%ns,2))
        do k = 1, SPS%ns
            RfS(k,1) = k <= Z
            RfS(k,2) = k <= Nneu
        end do
        AcS = .true.
        if (allocated(RefO) .and. Trunc == 'npnh') then
            do t = 1, 2
                do k = 1, SPS%ns
                    RfS(k,t) = RefO(SPS%a(k),t)
                    AcS(k,t) = ActO(SPS%a(k),t) .or. RfS(k,t)
                end do
            end do
        end if
        NtP = count(RfS(:,1))
        NtN = count(RfS(:,2))
        ! SDs with h = Hcap pair only with the other reference (h = 0, M = 0, parity +): kept for 2M, 2M + 2 (<J^2>) & the
            ! parity of the run (not with transitions from the parity + g.s., which share the SDs of both parities) ...
        KpM = [M2, M2 + 2]
        ! strength functions only (no eigenstates, no <J^2>): the SDs of 2M + 2 & its basis are not needed ...
        SfOnly = Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_Gs_Nmax > 0 .and. Params%Calc%NCSM%Tr_Lanczos > 0 .and. &
                 Nmax /= Params%Calc%NCSM%Tr_Gs_Nmax
        if (SfOnly) KpM = M2
        KpP = -1
        if (Trunc == 'npnh') then
            KpP = mod(ParIn + sum(SPS%e,RfS(:,1)) + sum(SPS%e,RfS(:,2)),2)
            ! (the larger runs of Tr_Gs_Nmax keep the pruning: the kept g.s. has no SD with h = Hcap > Tr_Gs_Nmax) ...
            if (Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_From /= 'HF' .and. &
                .not. (Params%Calc%NCSM%Tr_Gs_Nmax > 0 .and. Nmax > Params%Calc%NCSM%Tr_Gs_Nmax)) KpP = -1
        end if
        if (UseWin) write(*,'(a)') TAB // 'window: ' // istr(NtP) // ' proton & ' // istr(NtN) // ' neutron active holes (m-states), ' // &
            istr(count(AcS(:,1) .and. .not. RfS(:,1))) // ' & ' // istr(count(AcS(:,2) .and. .not. RfS(:,2))) // &
            ' particle m-states (frozen core: ' // istr(Z - NtP) // ' p, ' // istr(Nneu - NtN) // ' n)'

        ! Slater determinants ...
        call timer_start(t0)
        Qcap = Nmax
        Hcap = huge(1)
        if (Trunc == 'npnh') then
            Qcap = 100000
            Hcap = Nmax
        end if
        ! Energy cut ('npnh', Basis 'HF') ... Ecn(n) = Ecut_npnh(min(n, last given)), Ecn(0) none; e = HF single-particle
            ! energies (diagonal of H1b); SDs pruned by their own excitation energy when every species has
            ! a gap (all reference e below all others: E_other >= 0) ...
        Tag = ''
        UseCut = Trunc == 'npnh' .and. any(Params%Calc%NCSM%Ecut_npnh >= 0.0_dp)
        if (UseCut) then
            if (Params%Calc%NCSM%Basis /= 'HF') call fatal('NCSM ... Ecut_npnh needs Trunc npnh & Basis HF')
            Lg = 0
            do n = 1, size(Params%Calc%NCSM%Ecut_npnh)
                if (Params%Calc%NCSM%Ecut_npnh(n) >= 0.0_dp) Lg = n
            end do
            if (allocated(Ecn)) deallocate(Ecn)
            allocate(Ecn(0:Nmax), ExmP(0:Nmax), ExmN(0:Nmax), EpsP(SPS%ns), EpsN(SPS%ns))
            Ecn(0) = huge(1.0_dp)
            do n = 1, Nmax
                Ecn(n) = Params%Calc%NCSM%Ecut_npnh(min(n,Lg))
                if (Ecn(n) < 0.0_dp) Ecn(n) = huge(1.0_dp)
            end do
            do k = 1, SPS%ns
                EpsP(k) = EpsHF(SPS%a(k),1)
                EpsN(k) = EpsHF(SPS%a(k),2)
            end do
            do k = 0, Nmax
                ExmP(k) = maxval(Ecn(k:Nmax))
            end do
            ExmN = ExmP
            if (.not. (gapm(EpsP,RfS(:,1),AcS(:,1)) .and. gapm(EpsN,RfS(:,2),AcS(:,2)))) then
                write(*,'(a)') TAB // 'Warning ... no HF gap (a reference e above a particle e): SDs not pruned by energy'
                ExmP = huge(1.0_dp)
                ExmN = huge(1.0_dp)
            end if
            Tag = '_Ec'
            do n = 1, Nmax
                if (n > 1) Tag = Tag // '-'
                if (Ecn(n) >= huge(1.0_dp)) then
                    Tag = Tag // 'x'
                else
                    Tag = Tag // istr(nint(Ecn(n)))
                end if
            end do
            if (all(Ecn(1:Nmax) == Ecn(1))) then     ! one cut for all n: _Ec<E>
                if (Ecn(1) >= huge(1.0_dp)) then
                    Tag = '_Ecx'
                else
                    Tag = '_Ec' // istr(nint(Ecn(1)))
                end if
            end if
            write(*,'(a)') TAB // 'energy cut ... E_P + E_N <= Ecut(n) for npnh, n = 1 .. ' // istr(Nmax) // ' (' // Tag(4:) // &
                           ' MeV), HF e: reference <= ' // cfmt_f(max(maxval(EpsP,RfS(:,1)),maxval(EpsN,RfS(:,2))),0,3) // &
                           ', particles >= ' // cfmt_f(min(minval(EpsP,AcS(:,1) .and. .not. RfS(:,1)), &
                                                           minval(EpsN,AcS(:,2) .and. .not. RfS(:,2))),0,3) // ' MeV'
            call sd_make(SPS,NtP,Qcap,Hcap,P,EpsP,ExmP,RfS(:,1),AcS(:,1),KpM,KpP)
            call sd_make(SPS,NtN,Qcap,Hcap,Nn,EpsN,ExmN,RfS(:,2),AcS(:,2),KpM,KpP)
            call timer_report(t0,'proton & neutron SDs & jumps (energy cut)')
        else
        call sd_make(SPS,NtP,Qcap,Hcap,P,Ref=RfS(:,1),Act=AcS(:,1),KeepM=KpM,KeepP=KpP)
        if (Nneu == Z .and. all(RfS(:,1) .eqv. RfS(:,2)) .and. all(AcS(:,1) .eqv. AcS(:,2))) then
            Nn = P            ! identical SDs & jumps (same s.p. space & truncation)
            call timer_report(t0,'proton SDs & jumps (neutrons copied)')
        else
            call sd_make(SPS,NtN,Qcap,Hcap,Nn,Ref=RfS(:,2),Act=AcS(:,2),KeepM=KpM,KeepP=KpP)
        end if
        end if
        Par = mod(P%N0 + Nn%N0 + Nmax,2)
        if (Trunc == 'npnh') Par = ParIn
        write(*,'(a)') TAB // istr(SPS%ns) // ' m-states per species, ' // istr(P%nsd) // ' proton & ' // istr(Nn%nsd) // &
                       ' neutron SDs, parity ' // merge('-','+',Par == 1) // ', 2M = ' // istr(M2)
        call basis_make(P,Nn,M2,Nmax,mod(Par + P%N0 + Nn%N0,2),B)
        write(*,'(a)') TAB // 'm-scheme dimension ' // istr(B%dim)
        if (UseCut) write(*,'(a)') TAB // '  (energy cut: ' // cfmt_f(100.0_dp*real(B%dim,dp)/real(max(B%dim_all,1_i8),dp),0,2) // &
                                   ' % of the ' // istr8(B%dim_all) // ' products of the SDs kept)'
        ! the basis of 2M + 2 (same sectors) for <J^2> = |J_+ z|^2 + M (M + 1) ...
        if (.not. SfOnly) call basis_make(P,Nn,M2 + 2,Nmax,mod(Par + P%N0 + Nn%N0,2),B2)
        call timer_report(t0,'basis & one-body jumps')
        call mem_report('SDs, jumps & basis')
        if (B%dim == 0) return

        ! Pair blocks within the quanta cut: the largest quanta of two nucleons in the basis ...
        if (Trunc == 'npnh') then
            Epair = 2*maxval(SPS%e)
        else
            Epair(1) = merge(SPS%e(max(Z,1)) + SPS%e(max(Z-1,1)), -1000, Z >= 2)
            Epair(2) = merge(SPS%e(max(Nneu,1)) + SPS%e(max(Nneu-1,1)), -1000, Nneu >= 2)
            Epair(3) = SPS%e(Z) + SPS%e(Nneu)
            Epair = min(Epair + Qcap, 2*maxval(SPS%e))
        end if
        do t = 1, 3
            call pairs_make(SPS,t,Epair(t),Pr(t))
        end do
        if (maxval(Epair) > Params%Calc%N2max) write(*,'(a)') TAB // 'Warning ... pair quanta up to ' // &
            istr(maxval(Epair)) // ' > N2max = ' // istr(Params%Calc%N2max) // ' ... missing 2-body matrix elements ...'

        ! Operators H + Beta N_cm & N_cm ...
        call timer_start(t0)
        call op_from_coupled(SPS,Orb,Pr,V,Orb_NN,uH,C0,H)
        call op_ncm(SPS,Pr,A,Ncm)
        if (Beta /= 0.0_dp) call op_axpy(H,Beta,Ncm)
        call pn_tables(SPS,Pr(3),H)
        call timer_report(t0,'m-scheme tables of H + Beta N_cm')
        call timer_start(t0)
        Single = Params%Calc%NCSM%Single
        if (UseCut) then
            call eng_make(SPS,P,Nn,Pr,H,B,EH,Single,Qcap,Hcap,EpsP,EpsN,ExmP,ExmN)
        else
            call eng_make(SPS,P,Nn,Pr,H,B,EH,Single,Qcap,Hcap,Sym=SfOnly .and. Params%Calc%NCSM%Tr_Sym .and. &
                          trim(Params%Calc%NCSM%Matrix) == 'otf',Eta=1 - 2*Par)
        end if
        call timer_report(t0,'species matrices of H + Beta N_cm')
        call mem_report('species matrices')

        ! Stored or on the fly ...
        Mode = Params%Calc%NCSM%Matrix
        MemGB = Params%Calc%NCSM%Matrix_GB
        if (MemGB <= 0.0_dp) MemGB = 0.7_dp*mem_available_gb()
        if (trim(Mode) /= 'otf') then
            call timer_start(t0)
            call eng_store(P,Nn,B,H,EH,Single,MemGB,Fits)
            call timer_report(t0,'stored matrix')
            if (.not. Fits) then
                ! Matrix = 'stored' insists on the stored form (on-the-fly products may take days in large spaces) ...
                if (trim(Mode) == 'stored') call fatal('NCSM ... stored matrix exceeds Matrix_GB (Matrix = ''stored'': no on-the-fly fallback)')
                write(*,'(a)') TAB // 'stored matrix exceeds Matrix_GB ... on-the-fly products instead'
                Mode = 'otf'
            else
                Mode = 'stored'
                deallocate(EH%Mp, EH%Mn)      ! not needed by the stored products
                call op_free(H)
            end if
        end if
        write(*,'(a)') TAB // 'products of H + Beta N_cm: ' // trim(Mode)
        call mem_report('before Lanczos')

        Fn = io_path(Params,'/NCSM/NCSM_' // trim(Params%Calc%NCSM%Hamiltonian) // '_Nmax' // istr(Nmax))
        ! Nmax in the HF basis (quanta e = 2n + l of the HF orbitals, i.e. HF major shells) ... NCSM_EMPM_HF_Nmax<N> ...
        if (Params%Calc%NCSM%Basis == 'HF') Fn = io_path(Params,'/NCSM/NCSM_' // trim(Params%Calc%NCSM%Hamiltonian) // &
                                                          '_HF_Nmax' // istr(Nmax))
        if (Trunc == 'npnh') Fn = io_path(Params,'/NCSM/NCSM_' // trim(Params%Calc%NCSM%Hamiltonian) // '_' // &
                                          trim(Params%Calc%NCSM%Basis) // '_' // istr(Nmax) // 'p' // istr(Nmax) // 'h' // merge('m','p',Par == 1))
        Fn = Fn // Tag
        Hdr = ''
        if (UseCut) Hdr = ' Ecut_npnh = ' // Tag(4:) // ' (n = 1..' // istr(Nmax) // ')'
        if (Beta /= 0.0_dp) Fn = Fn // '_B' // istr(nint(Beta))
        Fn = Fn // '.dat'

        ! Lanczos for H + Beta N_cm ...
        Ne = min(Params%Calc%NCSM%N_Eig, B%dim)
        ! the larger runs of Tr_Gs_Nmax with strength functions: no eigenstates (the converged Ritz states of the SF) ...
        if (Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_Gs_Nmax > 0 .and. Params%Calc%NCSM%Tr_Lanczos > 0 .and. &
            Nmax /= Params%Calc%NCSM%Tr_Gs_Nmax) Ne = 0
        call timer_start(t0)
        cP => P
        cNn => Nn
        cPr => Pr
        cB => B
        cH => H
        cEH => EH
        ! transitions from the HF reference (Tr_From 'HF'): Xg the reference determinant in the parity + basis of these
            ! SDs ... strength-function Lanczos from Q|g> (Tr_Lanczos) ...
        if (Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_From == 'HF') call hf_ref(P,Nn,M2,Nmax,E_HF)
        ! the g.s. of the smaller run (Tr_Gs_Nmax) in the parity + basis of these SDs ...
        if (Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_From == 'gs' .and. Params%Calc%NCSM%Tr_Gs_Nmax > 0 .and. &
            Nmax /= Params%Calc%NCSM%Tr_Gs_Nmax) call gs_remap(SPS,P,Nn,Nmax)
        if (Params%Calc%NCSM%Transitions .and. Params%Calc%NCSM%Tr_Lanczos > 0) then
            call trans_strength(Params,SPS,P,Nn,B,Par,Fn(1:len(Fn)-4))
            call timer_report(t0,'strength-function Lanczos')
        end if
        if (Ne <= 0) then
            nullify(cP, cNn, cPr, cB, cH, cEH)
            call eng_free(EH)
            call op_free(H)
            call op_free(Ncm)
            return
        end if
        call timer_start(t0)
        if (B%dim <= Params%Calc%NCSM%Dense_Max) then
            call dense_all(B%dim,E,Zv,Res,Nconv)
        else
            call lanczos_tr(B%dim,apply_h,Ne,Params%Calc%NCSM%Lanczos_Basis,Params%Calc%NCSM%Lanczos_Max, &
                            Params%Calc%NCSM%Lanczos_Tol,E,Zv,Res,Nconv)
        end if
        nullify(cP, cNn, cPr, cB, cH, cEH)
        call timer_report(t0,merge('dense   ','Lanczos ',B%dim <= Params%Calc%NCSM%Dense_Max))
        call mem_report('Lanczos')
        Ne = size(E)
        call eng_free(EH)
        call op_free(H)

        ! <N_cm> (2-body, lower-half expectation pass) & <J^2> (one-body: |J_+ z|^2 + M (M + 1)) ...
        call timer_start(t0)
        allocate(Ev(Ne), Jv(Ne), Nv(Ne), Xt(Ne,B%dim))
        Xt = transpose(Zv)
        if (UseWin) then
            Nv = 0.0_dp            ! N_cm of the window only: not computed
        else
        call pn_tables(SPS,Pr(3),Ncm)
        if (UseCut) then
            call eng_make(SPS,P,Nn,Pr,Ncm,B,EO,.false.,Qcap,Hcap,EpsP,EpsN,ExmP,ExmN)
        else
            call eng_make(SPS,P,Nn,Pr,Ncm,B,EO,.false.,Qcap,Hcap)
        end if
        call eng_expect(P,Nn,B,Ncm,EO,Ne,Xt,Nv)
        call eng_free(EO)
        end if
        call op_free(Ncm)
        call j2_expect(SPS,P,Nn,B,B2,Ne,Xt,Jv)
        deallocate(Xt)
        Ev = E - Beta*Nv
        call timer_report(t0,'<N_cm> & <J^2>')

        open(newunit=u, file=Fn, status='replace', action='write')
        write(u,'(a)') '# NCSM ' // trim(Params%Calc%NCSM%Hamiltonian) // ' A = ' // istr(A) // ' Z = ' // istr(Z) // &
                       ' hw = ' // jstr(Params%Calc%hw) // ' ' // merge('npnh','Nmax',Trunc == 'npnh') // ' = ' // istr(Nmax) // &
                       ' basis ' // trim(Params%Calc%NCSM%Basis) // ' parity ' // merge('-','+',Par == 1) // &
                       ' 2M = ' // istr(M2) // ' dim = ' // istr(B%dim) // ' Beta_CM = ' // jstr(Beta) // ' E_HF = ' // jstr(E_HF) // Hdr
        write(u,'(a)') '# k   E (MeV)   E - E_HF   J   <J^2>   <N_cm>   residual'
        write(*,'(a)') TAB // '  k        E (MeV)     E - E_HF       J      <J^2>     <N_cm>   residual'
        do k = 1, Ne
            write(u,'(i5,2f16.8,i5,f12.6,es14.5,es12.3e3)') k, Ev(k), Ev(k) - E_HF, nint(jval(Jv(k))), Jv(k), Nv(k), Res(k)
            if (k <= 60) write(*,'(a,i5,2f14.6,f8.3,f11.5,es11.3,es12.3e3)') TAB, k, Ev(k), Ev(k) - E_HF, jval(Jv(k)), &
                                                                            Jv(k), Nv(k), Res(k)
        end do
        close(u)
        if (Nconv < Ne) write(*,'(a)') TAB // 'Warning ... only ' // istr(Nconv) // ' of ' // istr(Ne) // ' states converged'
        write(*,'(a)') TAB // 'written ' // Fn
        if (Params%Calc%NCSM%Transitions) then
            call timer_start(t0)
            call trans_out(Params,SPS,P,Nn,B,Par,Ne,Zv,Ev,Jv,Fn(1:len(Fn)-4))
            call timer_report(t0,'transitions')
        end if

    contains

        real(dp) function jval(x)
            real(dp), intent(in) :: x
            jval = 0.5_dp*(-1.0_dp + sqrt(max(1.0_dp + 4.0_dp*x,0.0_dp)))
        end function jval

        ! all reference energies <= all active particle energies ...
        logical function gapm(e,rf,ac)
            real(dp), intent(in) :: e(:)
            logical, intent(in) :: rf(:), ac(:)
            gapm = .true.
            if (any(ac .and. .not. rf)) gapm = maxval(e,rf) <= minval(e,ac .and. .not. rf)
        end function gapm

        ! all reference energies e(1:nt) <= all others ...
        logical function gap(e,nt)
            real(dp), intent(in) :: e(:)
            integer, intent(in) :: nt
            gap = .true.
            if (nt < size(e)) gap = maxval(e(1:nt)) <= minval(e(nt+1:))
        end function gap

    end subroutine run_nmax

    ! Window of &empm (pHole, pPart, nHole, nPart: HF levels per isospin by energy, as EMPM_space_make) ... RefO the
        ! window holes, ActO the window levels (a, t) ...
    subroutine window_make(Params,Orb)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(empm_space_t) :: Space
        integer :: k, t
        if (Params%Calc%NCSM%Trunc /= 'npnh' .or. Params%Calc%NCSM%Basis /= 'HF' .or. &
            Params%Calc%NCSM%Hamiltonian /= 'EMPM') call fatal('NCSM ... Window needs Trunc npnh, Basis HF, Hamiltonian EMPM')
        if (Params%Calc%NCSM%Beta_CM /= 0.0_dp) call fatal('NCSM ... Window: no Lawson term (N_cm of the window only)')
        call EMPM_space_make(Params,Orb,Space)
        allocate(RefO(size(Orb),2), ActO(size(Orb),2))
        RefO = .false.
        ActO = .false.
        do k = 1, Space%N_lev
            t = merge(1,2,Space%Lev(k)%t == -1)
            ActO(Space%Lev(k)%a,t) = .true.
            RefO(Space%Lev(k)%a,t) = Space%Lev(k)%Occ == 1
        end do
        write(*,'(a)') TAB // 'NCSM window ... ' // istr(count(RefO(:,1))) // ' & ' // istr(count(RefO(:,2))) // &
                       ' hole orbitals, ' // istr(count(ActO(:,1) .and. .not. RefO(:,1))) // ' & ' // &
                       istr(count(ActO(:,2) .and. .not. RefO(:,2))) // ' particle orbitals (p & n) ... <N_cm> of the window only'
    end subroutine window_make

    ! ---------------------------------------------------------------------------------------------------------------
    ! Electric transitions from the ground state (&ncsm Transitions) ...

    ! Reduced one-body matrices of E0 (r^2), E1, E2, E3 (nj_tr1b, Edmonds convention) in the orbitals of the run ...
    subroutine trans_init(Params,Orb)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(tr1b_t) :: TrOp, TrH
        type(o1b_t) :: C
        integer :: a_max, M2
        M2 = Params%Calc%NCSM%M2
        if (M2 < 0) M2 = mod(Params%Calc%A,2)
        if (Params%Calc%NCSM%Trunc /= 'npnh' .or. M2 /= 0) call fatal('NCSM ... Transitions need Trunc = ''npnh'' & M2 = 0')
        if (Params%Calc%NCSM%Tr_From /= 'HF' .and. Params%Calc%NCSM%Tr_From /= 'gs') call fatal('NCSM ... Tr_From gs or HF')
        if (Params%Calc%NCSM%Tr_From == 'gs' .and. Params%Calc%NCSM%Parity == 1) &
            call fatal('NCSM ... Tr_From gs needs the parity + run (Parity = -1)')
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        TrH = Tr1b_initialize(Params,Orb,1.0_dp)
        HfBasis = Params%Calc%NCSM%Basis == 'HF'
        if (HfBasis) then
            C = O1b_import(Params,Orb,empm_ref_path(Params,'/bin/C.bin'))
            TrOp = Tr1b_transformation(Params,Orb,C,TrH)
        else
            TrOp = TrH
        end if
        if (allocated(QL)) deallocate(QL)
        if (allocated(QH)) deallocate(QH)
        allocate(QL(a_max,a_max,2,0:3), QH(a_max,a_max,2,0:3))
        call put(QL,TrOp)
        call put(QH,TrH)
        Bosc = sqrt(0.5_dp*(m_p + m_n)*Params%Calc%hw)/HbarC
        if (allocated(Xg)) deallocate(Xg)
        write(*,'(a)') TAB // 'transitions: E0 (r^2), E1, E2, E3 from the lowest state of parity + (' // &
                       merge('HF','HO',HfBasis) // ' orbitals)'
    contains
        subroutine put(Q,T)
            real(dp), intent(out) :: Q(:,:,:,0:)
            type(tr1b_t), intent(in) :: T
            Q(:,:,1,0) = T%E0%p
            Q(:,:,2,0) = T%E0%n
            Q(:,:,1,1) = T%E1%p
            Q(:,:,2,1) = T%E1%n
            Q(:,:,1,2) = T%E2%p
            Q(:,:,2,2) = T%E2%n
            Q(:,:,1,3) = T%E3%p
            Q(:,:,2,3) = T%E3%n
        end subroutine put
    end subroutine trans_init

    ! m-scheme <k|Q_L0|l> from the reduced matrix Q (Edmonds: <j m|T_L0|j' m> = <j' m L 0|j m> <j||T||j'> / sqrt(2j+1)) ...
    real(dp) function qm(SPS,Q,L,k,kl,t)
        type(sp_t), intent(in) :: SPS
        real(dp), intent(in) :: Q(:,:,:,0:)
        integer, intent(in) :: L, k, kl, t
        qm = 0.0_dp
        if (SPS%m2(k) /= SPS%m2(kl)) return
        if (mod(SPS%l(k) + SPS%l(kl) + L,2) /= 0) return
        if (abs(SPS%j2(k) - SPS%j2(kl)) > 2*L .or. SPS%j2(k) + SPS%j2(kl) < 2*L) return
        qm = fCG(SPS%j2(kl),2*L,SPS%j2(k),SPS%m2(kl),0,SPS%m2(k))*Q(SPS%a(k),SPS%a(kl),t,L)/sqrt(real(SPS%j2(k) + 1,dp))
    end function qm

    ! Checks of the m-scheme tables: hermiticity (all L) & E1 against the HO r_0 = b/sqrt(2) (b_0 + b_0^+) of b_mscheme
        ! (the (l s) j convention of the NCSM operators) ... printed once per run ...
    subroutine trans_check(SPS)
        type(sp_t), intent(in) :: SPS
        real(dp), allocatable :: Bm(:,:,:)
        real(dp) :: dh, dz, z, q, zmax
        integer :: k, kl, L_, t
        dh = 0.0_dp
        do L_ = 0, 3
            do t = 1, 2
                do k = 1, SPS%ns
                    do kl = 1, SPS%ns
                        dh = max(dh,abs(qm(SPS,QL,L_,k,kl,t) - qm(SPS,QL,L_,kl,k,t)))
                    end do
                end do
            end do
        end do
        call b_mscheme(SPS,Bm)
        dz = 0.0_dp
        zmax = 0.0_dp
        do k = 1, SPS%ns
            do kl = 1, SPS%ns
                z = sqrt(3.0_dp/(4.0_dp*pi))/Bosc/sqrt(2.0_dp)*(Bm(k,kl,0) + Bm(kl,k,0))
                q = qm(SPS,QH,1,k,kl,1)
                dz = max(dz,abs(q - z))
                zmax = max(zmax,abs(z))
            end do
        end do
        write(*,'(a,es10.2,a,es10.2,a,f8.4,a)') TAB // 'transition tables: max |q(k,l) - q(l,k)| = ', dh, &
            ', HO E1 vs sqrt(3/4pi) z of b_mscheme: max dev ', dz, ' (max |z| ', zmax, ' fm)'
        if (dz > 1.0e-8_dp*max(zmax,1.0_dp)) write(*,'(a)') TAB // 'Warning ... E1 table differs from the NCSM (l s) j convention'
    end subroutine trans_check

    ! yp, yn = Q_L0 (protons, neutrons) applied to the ground state Xg (basis Bg) in the target basis Bt ... the reverse jumps
        ! of j2_expect: y(cc) = sum over the jumps a+_kc a_ka |sd(cc)> = sg |s> of <ka|Q|kc> sg Xg(c(s)) ... the diagonal
        ! bucket dm = 0 with the parity type of L ...
    subroutine tr_apply(SPS,P,Nn,Bt,L,yp,yn)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: Bt
        integer, intent(in) :: L
        real(dp), intent(out) :: yp(:), yn(:)
        real(dp), allocatable :: vp(:), vn(:)
        integer :: br, ix, kc, ka, cc, spp, snp, s, st, c, L2
        integer(i8) :: j
        real(dp) :: v, ap, an
        br = SPS%b0 + mod(L,2)
        allocate(vp(SPS%Nty(br)), vn(SPS%Nty(br)))
        do ix = 1, SPS%Nty(br)
            kc = SPS%Tkc(br)%v(ix)
            ka = SPS%Tka(br)%v(ix)
            vp(ix) = qm(SPS,QL,L,ka,kc,1)
            vn(ix) = qm(SPS,QL,L,ka,kc,2)
        end do
        !$omp parallel do schedule(dynamic,256) private(spp,snp,s,st,c,L2,j,v,ap,an)
        do cc = 1, Bt%dim
            spp = Bt%sp(cc)
            snp = Bt%sn(cc)
            ap = 0.0_dp
            an = 0.0_dp
            L2 = Nn%L(snp)
            do j = P%Jo(br,spp), P%Jo(br+1,spp) - 1
                v = vp(int(P%Jx(j)) + IXOFF)
                if (v == 0.0_dp) cycle
                s = abs(P%Jt(j))
                st = Bg%St(s,L2)
                if (st < 0) cycle
                if (UseCut) then
                    if (Bg%Lc(snp) > Bg%Nc(s,L2)) cycle
                end if
                c = st + Bg%Lc(snp)
                ap = ap + v*real(sign(1,P%Jt(j)),dp)*Xg(c)
            end do
            do j = Nn%Jo(br,snp), Nn%Jo(br+1,snp) - 1
                v = vn(int(Nn%Jx(j)) + IXOFF)
                if (v == 0.0_dp) cycle
                s = abs(Nn%Jt(j))
                L2 = Nn%L(s)
                st = Bg%St(spp,L2)
                if (st < 0) cycle
                if (UseCut) then
                    if (Bg%Lc(s) > Bg%Nc(spp,L2)) cycle
                end if
                c = st + Bg%Lc(s)
                an = an + v*real(sign(1,Nn%Jt(j)),dp)*Xg(c)
            end do
            yp(cc) = ap
            yn(cc) = an
        end do
        !$omp end parallel do
    end subroutine tr_apply

    ! B(EL; g.s. -> k) for the states of this run (M = 0): M_t(k) = sqrt(2L+1) <k|Q_L0,t|g.s.> (= <J_k||Q_L,t||0>,
        ! J_k = L), B_ph = M_p^2, B_is = |M_p + M_n|^2/4, B_iv = |e_p M_p - e_n M_n|^2 (E1, e_p = N/A, e_n = Z/A) or
        ! |M_p - M_n|^2/4 (as HF_RPA_rB) ... sum rule (2L+1) |Q_L0,p g.s.|^2 = sum over all states of B_ph (complete only with
        ! all states, Dense_Max) ... files <Fb>_E<L>.dat ...
    subroutine trans_out(Params,SPS,P,Nn,B,Par,Ne,Zv,Ev,Jv,Fb)
        type(params_t), intent(in) :: Params
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        integer, intent(in) :: Par, Ne
        real(dp), intent(in) :: Zv(:,:), Ev(:), Jv(:)
        character(len=*), intent(in) :: Fb
        real(dp), allocatable :: yp(:), yn(:), Mp(:), Mn(:)
        real(dp) :: e_p, e_n, f, Sp, Sn, Bph, Bis, Biv, Tph, Tis, Tiv, Jk
        integer :: Ls(2), iL, L, k, u, A, Z, nJ
        logical, save :: Checked = .false.
        character(len=:), allocatable :: Fl
        A = Params%Calc%A
        Z = Params%Calc%Z
        e_p = real(A - Z,dp)/real(A,dp)
        e_n = real(Z,dp)/real(A,dp)
        if (.not. Checked) then
            call trans_check(SPS)
            Checked = .true.
        end if
        if (Params%Calc%NCSM%Tr_From == 'HF') then
            Ls = [0, 2]
            if (Par == 1) Ls = [1, 3]
        else if (Par == 0 .and. Params%Calc%NCSM%Tr_Gs_Nmax > 0 .and. CurNmax /= Params%Calc%NCSM%Tr_Gs_Nmax) then
            ! g.s. of the smaller run, mapped by gs_remap ...
            if (.not. allocated(Xg)) then
                write(*,'(a)') TAB // 'Warning ... transitions: no ground state of the Tr_Gs_Nmax run ... skipped'
                return
            end if
            Ls = [0, 2]
        else if (Par == 0) then
            Bg = B
            Xg = Zv(:,1)
            Eg = Ev(1)
            Jg = 0.5_dp*(-1.0_dp + sqrt(max(1.0_dp + 4.0_dp*Jv(1),0.0_dp)))
            write(*,'(a,f14.6,a,f8.4)') TAB // 'transitions: ground state E = ', Eg, ' MeV, J = ', Jg
            if (abs(Jg) > 1.0e-3_dp) write(*,'(a)') TAB // 'Warning ... the lowest state of parity + is not J = 0'
            if (Params%Calc%NCSM%Tr_Gs_Nmax > 0) then
                ! kept by the SD occupations for the larger runs ...
                if (allocated(Gp)) deallocate(Gp, Gn, Gx)
                allocate(Gp(size(P%occ,1),B%dim), Gn(size(Nn%occ,1),B%dim), Gx(B%dim))
                do k = 1, B%dim
                    Gp(:,k) = P%occ(:,B%sp(k))
                    Gn(:,k) = Nn%occ(:,B%sn(k))
                end do
                Gx = Xg
                write(*,'(a)') TAB // 'transitions: ground state of the ' // istr(CurNmax) // 'p' // istr(CurNmax) // &
                               'h run kept for the larger runs (' // istr(B%dim) // ' components)'
            end if
            Ls = [0, 2]
        else
            if (.not. allocated(Xg)) then
                write(*,'(a)') TAB // 'Warning ... transitions: no ground state (parity + not run first) ... skipped'
                return
            end if
            Ls = [1, 3]
        end if
        allocate(yp(B%dim), yn(B%dim), Mp(Ne), Mn(Ne))
        do iL = 1, 2
            L = Ls(iL)
            call tr_apply(SPS,P,Nn,B,L,yp,yn)
            f = sqrt(real(2*L + 1,dp))
            Sp = real(2*L + 1,dp)*dot_product(yp,yp)
            Sn = real(2*L + 1,dp)*dot_product(yn,yn)
            do k = 1, Ne
                Mp(k) = f*dot_product(Zv(:,k),yp)
                Mn(k) = f*dot_product(Zv(:,k),yn)
            end do
            Fl = Fb // '_E' // istr(L) // '.dat'
            open(newunit=u, file=Fl, status='replace', action='write')
            write(u,'(a)') '# NCSM B(E' // istr(L) // '; g.s. -> k), g.s. E = ' // jstr(Eg) // ' J = ' // jstr(Jg) // &
                           merge(' (HF reference)              ',' (lowest state of parity +)  ',Params%Calc%NCSM%Tr_From == 'HF') // &
                           ', states of parity ' // merge('-','+',Par == 1) // ', 2M = 0, dim = ' // &
                           istr(B%dim) // ', ' // istr(Ne) // ' states'
            write(u,'(a,2es16.8)') '# sum rule (2L+1)|Q_p g|^2, (2L+1)|Q_n g|^2 (e^2 fm^2L) = ', Sp, Sn
            write(u,'(a)') '# k   E (MeV)   Ex = E - E_gs   J   B_ph   B_is   B_iv   M_p   M_n'
            Tph = 0.0_dp
            Tis = 0.0_dp
            Tiv = 0.0_dp
            nJ = 0
            do k = 1, Ne
                Bph = Mp(k)**2
                Bis = 0.25_dp*(Mp(k) + Mn(k))**2
                if (L == 1) then
                    Biv = (e_p*Mp(k) - e_n*Mn(k))**2
                else
                    Biv = 0.25_dp*(Mp(k) - Mn(k))**2
                end if
                Jk = 0.5_dp*(-1.0_dp + sqrt(max(1.0_dp + 4.0_dp*Jv(k),0.0_dp)))
                write(u,'(i6,2f16.8,i5,3es16.8,2es16.8)') k, Ev(k), Ev(k) - Eg, nint(Jk), Bph, Bis, Biv, Mp(k), Mn(k)
                if (Par == 0 .and. k == 1 .and. Params%Calc%NCSM%Tr_From /= 'HF') cycle
                Tph = Tph + Bph
                Tis = Tis + Bis
                Tiv = Tiv + Biv
                if (nint(Jk) == L) nJ = nJ + 1
            end do
            close(u)
            write(*,'(a,i1,a,i7,a,3es12.4,a,es12.4,a)') TAB // 'E', L, ': ', nJ, ' states J = L, sum B_ph, B_is, B_iv =', &
                Tph, Tis, Tiv, ' (sum rule B_ph ', Sp, ')'
            write(*,'(a)') TAB // 'written ' // Fl
        end do
    end subroutine trans_out

    ! The g.s. of the smaller run (Tr_Gs_Nmax, kept by its SD occupations Gp, Gn, amplitudes Gx) in the parity + basis
        ! Bg of the current SDs: Xg(c) = Gx for the matching (proton SD, neutron SD) ... a component outside the current
        ! basis (other window, cut) is dropped with a warning (norm lost) ... Eg, Jg those of the smaller run ...
    subroutine gs_remap(SPS,P,Nn,Nmax)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        integer, intent(in) :: Nmax
        integer(i8) :: kp, kn
        integer :: i, sp, sn, c, L2, miss
        real(dp) :: lost
        if (allocated(Gx)) then
            if (size(Gp,1) /= size(P%occ,1) .or. size(Gn,1) /= size(Nn%occ,1)) call fatal('NCSM ... Tr_Gs_Nmax: other SD sizes')
        end if
        if (.not. allocated(Gx)) then
            write(*,'(a)') TAB // 'Warning ... Tr_Gs_Nmax: no ground state kept (its run must come first in Nmax_List)'
            if (allocated(Xg)) deallocate(Xg)
            return
        end if
        call basis_make(P,Nn,0,Nmax,mod(P%N0 + Nn%N0,2),Bg)
        if (allocated(Xg)) deallocate(Xg)
        allocate(Xg(Bg%dim))
        Xg = 0.0_dp
        miss = 0
        lost = 0.0_dp
        do i = 1, size(Gx)
            c = -1
            kp = sd_key(SPS,P,Gp(:,i))
            kn = sd_key(SPS,Nn,Gn(:,i))
            if (kp >= 0 .and. kn >= 0) then
                sp = hashmap_lookup(P%Map,kp)
                sn = hashmap_lookup(Nn%Map,kn)
                if (sp > 0 .and. sn > 0) then
                    L2 = Nn%L(sn)
                    c = Bg%St(sp,L2)
                    if (c >= 0) then
                        if (UseCut) then
                            if (Bg%Lc(sn) > Bg%Nc(sp,L2)) c = -1 - c
                        end if
                        if (c >= 0) c = c + Bg%Lc(sn)
                    end if
                end if
            end if
            if (c < 0) then
                miss = miss + 1
                lost = lost + Gx(i)**2
            else
                Xg(c) = Gx(i)
            end if
        end do
        write(*,'(a,f14.6,a,es10.2)') TAB // 'transitions: kept ground state (' // istr(size(Gx)) // ' components, E = ', Eg, &
            ' MeV) mapped into the ' // istr(Nmax) // 'p' // istr(Nmax) // 'h parity + basis, norm^2 lost ', lost
        if (miss > 0) write(*,'(a)') TAB // 'Warning ... ' // istr(miss) // ' g.s. components not in the current basis'
    end subroutine gs_remap

    ! The HF reference determinant as the initial state (Tr_From 'HF'): the parity + basis of the current SDs, Xg = unit
        ! vector of the reference (proton & neutron SDs with h = 0), Eg = E_HF (<HF|H|HF> of the EMPM Hamiltonian) ...
    subroutine hf_ref(P,Nn,M2,Nmax,E_HF)
        type(sd_t), intent(in) :: P, Nn
        integer, intent(in) :: M2, Nmax
        real(dp), intent(in) :: E_HF
        integer :: sp0, sn0, s, c
        if (M2 /= 0) call fatal('NCSM ... Tr_From HF needs M2 = 0')
        sp0 = 0
        sn0 = 0
        do s = 1, P%nsd
            if (P%h(s) == 0) sp0 = s
        end do
        do s = 1, Nn%nsd
            if (Nn%h(s) == 0) sn0 = s
        end do
        if (sp0 == 0 .or. sn0 == 0) call fatal('NCSM ... Tr_From HF: no reference SD')
        call basis_make(P,Nn,0,Nmax,mod(P%N0 + Nn%N0,2),Bg)
        c = Bg%St(sp0,Nn%L(sn0))
        if (c < 0) call fatal('NCSM ... Tr_From HF: reference not in the basis')
        c = c + Bg%Lc(sn0)
        if (allocated(Xg)) deallocate(Xg)
        allocate(Xg(Bg%dim))
        Xg = 0.0_dp
        Xg(c) = 1.0_dp
        Eg = E_HF
        Jg = 0.0_dp
        write(*,'(a,f14.6,a)') TAB // 'transitions from the HF reference, E_HF = ', Eg, ' MeV'
    end subroutine hf_ref

    ! Strength functions (Tr_Lanczos iterations, full reorthogonalization) of B_ph, B_is, B_iv from the pivots Q|g> in the
        ! target basis (current H context of apply_h) ... Ritz values th_k & B_k = S0 U(1,k)^2, S0 = (2L+1)|Q g|^2 (the total
        ! strength), residual beta_m |U(m,k)| ... the operator Q_L0|g> (g of J = 0) is pure J = L: only J = L states ...
        ! files <Fb>_E<L>_SF_<ph|is|iv>.dat ...
    subroutine trans_strength(Params,SPS,P,Nn,B,Par,Fb)
        type(params_t), intent(in) :: Params
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        integer, intent(in) :: Par
        character(len=*), intent(in) :: Fb
        real(dp), allocatable :: yp(:), yn(:), y(:), V(:,:), w(:), al(:), be(:), T(:,:), th(:), Uv(:,:), cf(:), qc(:), qp(:)
        real(f4), allocatable :: Vs(:,:)
        logical :: Reo, SgV
        real(dp) :: tA0, tB0
        real(dp) :: e_p, e_n, S0, f, m0, m1, eR
        integer :: Ls(2), iL, L, io, m, j, k, n, u, mk, ps, j0, nck
        integer(i8) :: ts, nb
        logical :: Sel(3), Ck
        character(len=2), parameter :: Nm(3) = ['ph', 'is', 'iv']
        character(len=:), allocatable :: Fl, Fc, Fv
        if (.not. allocated(Xg)) then
            write(*,'(a)') TAB // 'Warning ... strength functions: no initial state (Tr_From gs needs the parity + run first)'
            return
        end if
        n = B%dim
        e_p = real(Params%Calc%A - Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        Ls = [0, 2]
        if (Par == 1) Ls = [1, 3]
        m = min(Params%Calc%NCSM%Tr_Lanczos, n)
        ! Lanczos vectors: double & full reorthogonalization: V(n,m) (the vector of step j is V(:,j)) ... otherwise the
            ! current & previous vectors qc, qp in double & (reorthogonalization) the copies Vs(n,m) in single ...
        Reo = Params%Calc%NCSM%Tr_Reorth
        SgV = lcase(trim(Params%Calc%NCSM%Tr_Vec)) == 'single'
        ! Tr_Ckpt > 0: checkpoint & restart (bytes per stored Lanczos vector nb) ...
        nck = Params%Calc%NCSM%Tr_Ckpt
        Ck = nck > 0
        nb = int(n,i8)*merge(4_i8,8_i8,SgV)
        allocate(yp(n), yn(n), y(n), w(n), al(m), be(m), cf(m))
        if (Reo .and. .not. SgV) then
            allocate(V(n,m))
        else
            allocate(qc(n), qp(n))
            if (Reo) allocate(Vs(n,m))
        end if
        write(*,'(a)') TAB // 'strength-function Lanczos ... ' // istr(m) // ' steps, vectors ' // merge('single','double',SgV) // &
                       ', ' // merge('full reorthogonalization','no reorthogonalization  ',Reo) // ', stored vectors ' // &
                       cfmt_f(merge(real(m,dp),0.0_dp,Reo)*real(n,dp)*merge(4.0_dp,8.0_dp,SgV)/1.0e9_dp,0,1) // ' GB'
        do iL = 1, 2
            L = Ls(iL)
            ! Tr_Ops: only the listed strength functions (E<L><ph|is|iv>, case-insensitive), '' = all ...
            do io = 1, 3
                Sel(io) = len_trim(Params%Calc%NCSM%Tr_Ops) == 0 .or. &
                          index(lcase(Params%Calc%NCSM%Tr_Ops), 'e' // istr(L) // Nm(io)) > 0
            end do
            if (.not. any(Sel)) cycle
            call tr_apply(SPS,P,Nn,B,L,yp,yn)
            f = real(2*L + 1,dp)
            do io = 1, 3
                if (.not. Sel(io)) cycle
                select case (io)
                case (1)
                    y = yp
                case (2)
                    y = 0.5_dp*(yp + yn)
                case default
                    if (L == 1) then
                        y = e_p*yp - e_n*yn
                    else
                        y = 0.5_dp*(yp - yn)
                    end if
                end select
                S0 = f*dot_product(y,y)
                if (S0 <= 0.0_dp) cycle
                if (cEH%Sym) then
                    if (.not. cEH%SecReady) call eng_sec_make(P,Nn,B,cEH)
                    eR = sym_eta(P,B,cEH,y)
                    write(*,'(a)') TAB // 'rotation symmetry R_y(pi) of Q|g> (E' // istr(L) // ' ' // Nm(io) // '): <R> = ' // cfmt_f(eR,0,12)
                    if (abs(abs(eR) - 1.0_dp) > 1.0e-6_dp) call fatal('NCSM ... Q|g> is not symmetric under R_y(pi): set Tr_Sym = .false.')
                    if ((cEH%RedP .or. cEH%RedN) .and. nint(eR) /= cEH%Eta) call fatal('NCSM ... the species matrices were reduced for the other symmetry of Q|g>: set Tr_Sym = .false.')
                    cEH%Eta = nint(eR)
                    call sym_fill(P,B,cEH,y)
                end if
                if (allocated(V)) then
                    call pscale(n,1.0_dp/sqrt(pdot(n,y,y)),y,V(:,1))
                else
                    call pscale(n,1.0_dp/sqrt(pdot(n,y,y)),y,qc)
                end if
                Fl = Fb // '_E' // istr(L) // '_SF_' // Nm(io) // '.dat'
                Fc = Fb // '_E' // istr(L) // '_SF_' // Nm(io) // '_ckpt.bin'
                Fv = Fb // '_E' // istr(L) // '_SF_' // Nm(io) // '_vecs.bin'
                mk = 0
                j0 = 1
                if (Ck) call ck_restore()
                do j = j0, m
                    call timer_start(ts)
                    tA0 = TmA
                    tB0 = TmB
                    if (allocated(V)) then
                        call apply_h(V(:,j),w)
                        al(j) = pdot(n,V(:,j),w)
                        call paxpy(n,-al(j),V(:,j),w)
                        if (j > 1) call paxpy(n,-be(j-1),V(:,j-1),w)
                        do ps = 1, 2
                            call cgs_pass(n,j,V,w,cf)
                        end do
                    else
                        if (Reo) call pcopy_s(n,qc,Vs(:,j))
                        call apply_h(qc,w)
                        al(j) = pdot(n,qc,w)
                        call paxpy(n,-al(j),qc,w)
                        if (j > 1) call paxpy(n,-be(j-1),qp,w)
                        if (Reo) then
                            do ps = 1, 2
                                call cgs_pass_s(n,j,Vs,w,cf)
                            end do
                        end if
                    end if
                    be(j) = sqrt(pdot(n,w,w))
                    mk = j
                    ! the poles of the steps so far written after every step (a run stopped by the walltime keeps them) ...
                    call sf_write(.false.)
                    call sf_history()
                    call timer_report(ts,'  SF Lanczos step ' // istr(j) // ' of ' // istr(m) // ' (E' // istr(L) // ' ' // Nm(io) // &
                                      '), poles written, product passes ' // cfmt_f(TmA - tA0,0,1) // ' + ' // &
                                      cfmt_f(TmB - tB0,0,1) // ' s')
                    if (j == m .or. be(j) < 1.0e-9_dp*max(maxval(abs(al(1:j))),1.0_dp)) exit   ! Krylov space exhausted
                    if (allocated(V)) then
                        call pscale(n,1.0_dp/be(j),w,V(:,j+1))
                    else
                        call pcopy(n,qc,qp)
                        call pscale(n,1.0_dp/be(j),w,qc)
                    end if
                    ! checkpoint: vector j appended every step, the state after step j every nck steps ...
                    if (Ck) then
                        call timer_start(ts)
                        if (Reo) call ck_vec(j)
                        if (mod(j,nck) == 0) call ck_save(j)
                        call timer_report(ts,'  SF checkpoint, step ' // istr(j) // merge(' (state saved)','              ',mod(j,nck) == 0))
                    end if
                end do
                call sf_write(.true.)
            end do
        end do
    contains
        ! Restart from the checkpoint Fc (if present): the basis dimension, operator & S0 must match ... alpha, beta of the
            ! steps 1..j_, the current & previous vectors, the stored vectors 1..j_ from Fv; the _steps.dat file rebuilt;
            ! the loop continues at j0 = j_ + 1 (j_ < Tr_Lanczos: a finished run is extended by a larger Tr_Lanczos) ...
        subroutine ck_restore()
            integer :: uc, uv_, n_, L_, io_, j_, i_
            integer(i8) :: sz
            real(dp) :: S0_
            logical :: ex
            inquire(file=Fc, exist=ex)
            if (.not. ex) then
                write(*,'(a)') TAB // 'checkpoint ... none found (' // Fc // '), Lanczos from step 1'
                return
            end if
            open(newunit=uc, file=Fc, access='stream', form='unformatted', status='old', action='read')
            read(uc) n_, L_, io_, j_, S0_
            if (n_ /= n .or. L_ /= L .or. io_ /= io) call fatal('NCSM ... checkpoint ' // Fc // ' is for another basis or operator')
            if (abs(S0_ - S0) > 1.0e-8_dp*S0) call fatal('NCSM ... checkpoint ' // Fc // ': S0 differs (another g.s. or Hamiltonian?)')
            if (j_ >= m) call fatal('NCSM ... checkpoint ' // Fc // ' has ' // istr(j_) // ' steps: Tr_Lanczos must be larger')
            read(uc) al(1:j_), be(1:j_)
            if (allocated(V)) then
                read(uc) V(:,j_+1)
            else
                read(uc) qc
                read(uc) qp
            end if
            close(uc)
            if (Reo) then
                inquire(file=Fv, size=sz)
                if (sz < int(j_,i8)*nb) call fatal('NCSM ... checkpoint vectors ' // Fv // ' missing or shorter than ' // istr(j_) // ' vectors')
                open(newunit=uv_, file=Fv, access='stream', form='unformatted', status='old', action='read')
                do i_ = 1, j_
                    if (allocated(V)) then
                        read(uv_, pos=int(i_-1,i8)*nb + 1) V(:,i_)
                    else
                        read(uv_, pos=int(i_-1,i8)*nb + 1) Vs(:,i_)
                    end if
                end do
                close(uv_)
            end if
            do mk = 1, j_
                call sf_history()
            end do
            mk = j_
            call sf_write(.false.)
            j0 = j_ + 1
            write(*,'(a)') TAB // 'checkpoint ... restart after step ' // istr(j_) // ' (' // Fc // ')'
            flush(6)
        end subroutine ck_restore
        ! stored Lanczos vector j_ at its position in Fv (file new at j_ = 1) ...
        subroutine ck_vec(j_)
            integer, intent(in) :: j_
            integer :: uv_
            if (j_ == 1) then
                open(newunit=uv_, file=Fv, access='stream', form='unformatted', status='replace', action='write')
            else
                open(newunit=uv_, file=Fv, access='stream', form='unformatted', status='old', action='readwrite')
            end if
            if (allocated(V)) then
                write(uv_, pos=int(j_-1,i8)*nb + 1) V(:,j_)
            else
                write(uv_, pos=int(j_-1,i8)*nb + 1) Vs(:,j_)
            end if
            close(uv_)
        end subroutine ck_vec
        ! state after step j_ (alpha, beta, the next vector & for the 3-vector recurrence the current one) to Fc.tmp, then
            ! renamed to Fc (a run killed while writing keeps the previous checkpoint) ...
        subroutine ck_save(j_)
            integer, intent(in) :: j_
            integer :: uc, st
            open(newunit=uc, file=Fc // '.tmp', access='stream', form='unformatted', status='replace', action='write')
            write(uc) n, L, io, j_, S0
            write(uc) al(1:j_), be(1:j_)
            if (allocated(V)) then
                write(uc) V(:,j_+1)
            else
                write(uc) qc
                write(uc) qp
            end if
            close(uc)
            call execute_command_line('mv -f "' // Fc // '.tmp" "' // Fc // '"', exitstat=st)
            if (st /= 0) call fatal('NCSM ... checkpoint: cannot rename ' // Fc // '.tmp')
        end subroutine ck_save
        ! poles of the mk steps: eigenvalues of the tridiagonal T(1:mk), written to Fl (log line when Fin) ...
        subroutine sf_write(Fin)
            logical, intent(in) :: Fin
            integer :: i_, k_
            allocate(T(mk,mk), th(mk), Uv(mk,mk))
            T = 0.0_dp
            do i_ = 1, mk
                T(i_,i_) = al(i_)
                if (i_ < mk) then
                    T(i_,i_+1) = be(i_)
                    T(i_+1,i_) = be(i_)
                end if
            end do
            call eig_sym(T,th,Uv)
            open(newunit=u, file=Fl, status='replace', action='write')
            write(u,'(a)') '# NCSM strength function B_' // Nm(io) // '(E' // istr(L) // ') from ' // &
                           merge('the HF reference   ','the parity + g.s.  ',Params%Calc%NCSM%Tr_From == 'HF') // &
                           ' E_g = ' // jstr(Eg) // ', parity ' // merge('-','+',Par == 1) // ', 2M = 0, dim = ' // istr(n) // &
                           ', Lanczos iterations ' // istr(mk)
            write(u,'(a,es16.8)') '# total strength S0 = (2L+1)|Q g|^2 (e^2 fm^2L) = ', S0
            write(u,'(a)') '# k   E (MeV)   Ex = E - E_g   B   residual'
            m0 = 0.0_dp
            m1 = 0.0_dp
            do k_ = 1, mk
                write(u,'(i6,2f16.8,es18.8e3,es12.3e3)') k_, th(k_), th(k_) - Eg, S0*Uv(1,k_)**2, be(mk)*abs(Uv(mk,k_))
                m0 = m0 + S0*Uv(1,k_)**2
                m1 = m1 + S0*Uv(1,k_)**2*(th(k_) - Eg)
            end do
            close(u)
            if (Fin) write(*,'(a,i1,a,i5,a,es12.4,a,f10.4,a)') TAB // 'E', L, ' B_' // Nm(io) // ': ', mk, ' Lanczos steps, m0 = ', m0, &
                ', centroid ', m1/m0, ' MeV ... ' // Fl
            flush(6)
            deallocate(T, th, Uv)
        end subroutine sf_write
        ! history of the fragmentation: the poles of every step mk appended to <Fl>_steps.dat (block per step, blank line
            ! between blocks: step, k, E, Ex, B, residual) & the Lanczos coefficients alpha_j, beta_j (j <= mk) rewritten to
            ! <Fl>_coef.dat (the poles of any step can be rebuilt from them) ...
        subroutine sf_history()
            integer :: i_, k_, uh
            character(len=:), allocatable :: Fh
            Fh = Fl(1:len(Fl)-4)
            allocate(T(mk,mk), th(mk), Uv(mk,mk))
            T = 0.0_dp
            do i_ = 1, mk
                T(i_,i_) = al(i_)
                if (i_ < mk) then
                    T(i_,i_+1) = be(i_)
                    T(i_+1,i_) = be(i_)
                end if
            end do
            call eig_sym(T,th,Uv)
            if (mk == 1) then
                open(newunit=uh, file=Fh // '_steps.dat', status='replace', action='write')
                write(uh,'(a)') '# NCSM strength function B_' // Nm(io) // '(E' // istr(L) // ') poles after every Lanczos step,' // &
                                ' E_g = ' // jstr(Eg) // ', dim = ' // istr(n)
                write(uh,'(a,es16.8)') '# total strength S0 = (2L+1)|Q g|^2 (e^2 fm^2L) = ', S0
                write(uh,'(a)') '# step   k   E (MeV)   Ex = E - E_g   B   residual'
            else
                open(newunit=uh, file=Fh // '_steps.dat', status='old', position='append', action='write')
                write(uh,'(a)') ''
            end if
            do k_ = 1, mk
                write(uh,'(i6,i6,2f16.8,es18.8e3,es12.3e3)') mk, k_, th(k_), th(k_) - Eg, S0*Uv(1,k_)**2, be(mk)*abs(Uv(mk,k_))
            end do
            close(uh)
            open(newunit=uh, file=Fh // '_coef.dat', status='replace', action='write')
            write(uh,'(a)') '# Lanczos coefficients of B_' // Nm(io) // '(E' // istr(L) // '), E_g = ' // jstr(Eg) // ', dim = ' // istr(n)
            write(uh,'(a,es24.16)') '# S0 = ', S0
            write(uh,'(a)') '# j   alpha_j (MeV)   beta_j (MeV)'
            do i_ = 1, mk
                write(uh,'(i6,2es24.16)') i_, al(i_), be(i_)
            end do
            close(uh)
            deallocate(T, th, Uv)
        end subroutine sf_history
    end subroutine trans_strength

    pure function lcase(s) result(t)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: t
        integer :: i, c
        t = s
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= 65 .and. c <= 90) t(i:i) = achar(c + 32)
        end do
    end function lcase

    ! All eigenpairs of H + Beta N_cm by LAPACK (Dense_Max): the dense matrix from n products with unit vectors ...
    subroutine dense_all(n,E,Zv,Res,Nconv)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: E(:), Zv(:,:), Res(:)
        integer, intent(out) :: Nconv
        real(dp), allocatable :: Hd(:,:), x(:)
        real(dp) :: asy
        integer :: i
        allocate(Hd(n,n), x(n), E(n), Zv(n,n), Res(n))
        x = 0.0_dp
        do i = 1, n
            x(i) = 1.0_dp
            call apply_h(x,Hd(:,i))
            x(i) = 0.0_dp
        end do
        asy = maxval(abs(Hd - transpose(Hd)))
        write(*,'(a,i0,a,es10.2)') TAB // 'dense diagonalization of dimension ', n, ', max |H - H^T| = ', asy
        Hd = 0.5_dp*(Hd + transpose(Hd))
        call eig_sym(Hd,E,Zv)
        Res = 0.0_dp
        Nconv = n
    end subroutine dense_all

    ! y = (H + Beta N_cm) x ... the product of run_nmax for lanczos_tr (context in cP, cNn, cPr, cB, cH, cEH) ...
    subroutine apply_h(x,y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        character(len=16) :: Ev
        integer :: Ln_
        real(dp), allocatable :: y2(:)
        integer(i8) :: t1, t2, tr
        integer, save :: Checked = -1
        if (cEH%Stored) then
            call eng_mv_stored(cEH,cB%dim,x,y)
            return
        end if
        Ev = ''
        call get_environment_variable('NJ_MV',Ev,Ln_)
        if (trim(Ev) == 'old') then
            call eng_mv_otf(cP,cNn,cB,cH,cEH,x,y)
            return
        end if
        call eng_mv_sec(cP,cNn,cB,cH,cEH,x,y)
        if (trim(Ev) == 'cls' .and. Checked /= cB%dim) then
            Checked = cB%dim
            block
                real(dp) :: Sm(0:8,0:8)
                integer :: c_
                Sm = 0.0_dp
                do c_ = 1, cB%dim
                    Sm(cP%h(cB%sp(c_)),cNn%h(cB%sn(c_))) = Sm(cP%h(cB%sp(c_)),cNn%h(cB%sn(c_))) + x(c_)*y(c_)
                end do
                do c_ = 0, 3
                    write(*,'(a,i2,4es24.15)') 'CLS hp=', c_, Sm(c_,0:3)
                end do
            end block
        end if
        ! NJ_MV = check: the first product compared with the column product (eng_mv_otf, per-thread buffers) ...
        if (trim(Ev) == 'prof' .and. Checked /= cB%dim) then
            Checked = cB%dim
            allocate(y2(size(y)))
            DoPN = .false.
            call system_clock(t1,tr)
            call eng_mv_sec(cP,cNn,cB,cH,cEH,x,y2)
            call system_clock(t2)
            write(*,'(a)') TAB // 'NJ_MV prof ... neutron part only: ' // cfmt_f(real(t2 - t1,dp)/real(tr,dp),0,2) // ' s (+ pass B)'
            DoPN = .true.
            DoN = .false.
            Cnt = .true.
            call eng_mv_sec(cP,cNn,cB,cH,cEH,x,y2)
            Cnt = .false.
            call system_clock(t1,tr)
            call eng_mv_sec(cP,cNn,cB,cH,cEH,x,y2)
            call system_clock(t2)
            write(*,'(a)') TAB // 'NJ_MV prof ... pn part only: ' // cfmt_f(real(t2 - t1,dp)/real(tr,dp),0,2) // ' s (+ pass B)'
            DoN = .true.
            deallocate(y2)
        end if
        if (trim(Ev) == 'check' .and. Checked /= cB%dim) then
            Checked = cB%dim
            allocate(y2(size(y)))
            call system_clock(t1,tr)
            call eng_mv_otf(cP,cNn,cB,cH,cEH,x,y2)
            call system_clock(t2)
            write(*,'(a)') TAB // 'NJ_MV check ... max |y_sec - y_otf| = ' // cfmt_e(maxval(abs(y - y2)),0,3) // ', max |y| = ' // &
                           cfmt_e(maxval(abs(y2)),0,3) // ', otf product ' // cfmt_f(real(t2 - t1,dp)/real(tr,dp),0,2) // ' s'
            deallocate(y2)
            if (allocated(cEH%Yb)) deallocate(cEH%Yb)
        end if
    end subroutine apply_h

    ! One classical Gram-Schmidt pass, c = V(:,1:j)^T w, w <- w - V(:,1:j) c ... OpenMP over contiguous row blocks
        ! (sequential MKL dgemv per block: a threaded MKL dgemv on the tall n x j panel runs on one core), the block
        ! partials of c summed in thread order (reproducible) ...
    subroutine cgs_pass(n,j,V,w,c)
        integer, intent(in) :: n, j
        real(dp), intent(in) :: V(n,*)
        real(dp), intent(inout) :: w(n)
        real(dp), intent(out) :: c(j)
        real(dp), allocatable :: Cp(:,:)
        integer :: nt, tid, i0, i1, nr, t
        nt = omp_get_max_threads()
        allocate(Cp(j,nt))
        !$omp parallel num_threads(nt) private(tid,i0,i1,nr,t)
        tid = omp_get_thread_num()
        i0 = int(int(n,8)*tid/nt) + 1
        i1 = int(int(n,8)*(tid + 1)/nt)
        nr = i1 - i0 + 1
        Cp(:,tid+1) = 0.0_dp
        if (nr > 0) call dgemv('T',nr,j,1.0_dp,V(i0,1),n,w(i0),1,0.0_dp,Cp(1,tid+1),1)
        !$omp barrier
        !$omp single
        c = Cp(:,1)
        do t = 2, nt
            c = c + Cp(:,t)
        end do
        !$omp end single
        if (nr > 0) call dgemv('N',nr,j,-1.0_dp,V(i0,1),n,c,1,1.0_dp,w(i0),1)
        !$omp end parallel
        deallocate(Cp)
    end subroutine cgs_pass

    ! cgs_pass with the vectors V stored in single precision (sums in double) ... row blocks of the threads, sub-blocks
        ! of 4096 rows (w in cache over the j vectors) ...
    subroutine cgs_pass_s(n,j,V,w,c)
        integer, intent(in) :: n, j
        real(f4), intent(in) :: V(n,*)
        real(dp), intent(inout) :: w(n)
        real(dp), intent(out) :: c(j)
        real(dp), allocatable :: Cp(:,:)
        integer, parameter :: SB = 4096
        integer :: nt, tid, i0, i1, t, k, b0, b1, i
        real(dp) :: s
        nt = omp_get_max_threads()
        allocate(Cp(j,nt))
        !$omp parallel num_threads(nt) private(tid,i0,i1,t,k,b0,b1,i,s)
        tid = omp_get_thread_num()
        i0 = int(int(n,8)*tid/nt) + 1
        i1 = int(int(n,8)*(tid + 1)/nt)
        Cp(:,tid+1) = 0.0_dp
        do b0 = i0, i1, SB
            b1 = min(i1, b0 + SB - 1)
            do k = 1, j
                s = 0.0_dp
                !$omp simd reduction(+:s)
                do i = b0, b1
                    s = s + real(V(i,k),dp)*w(i)
                end do
                Cp(k,tid+1) = Cp(k,tid+1) + s
            end do
        end do
        !$omp barrier
        !$omp single
        c = Cp(:,1)
        do t = 2, nt
            c = c + Cp(:,t)
        end do
        !$omp end single
        do b0 = i0, i1, SB
            b1 = min(i1, b0 + SB - 1)
            do k = 1, j
                s = c(k)
                !$omp simd
                do i = b0, b1
                    w(i) = w(i) - s*real(V(i,k),dp)
                end do
            end do
        end do
        !$omp end parallel
        deallocate(Cp)
    end subroutine cgs_pass_s

    subroutine pcopy(n,x,y)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)
        integer :: i
        !$omp parallel do schedule(static)
        do i = 1, n
            y(i) = x(i)
        end do
        !$omp end parallel do
    end subroutine pcopy

    subroutine pcopy_s(n,x,y)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(f4), intent(out) :: y(n)
        integer :: i
        !$omp parallel do schedule(static)
        do i = 1, n
            y(i) = real(x(i),f4)
        end do
        !$omp end parallel do
    end subroutine pcopy_s

    ! Threaded vector kernels of the Lanczos loops (same static row blocks as cgs_pass: NUMA first touch) ...
    function pdot(n,x,y) result(s)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        real(dp) :: s
        integer :: i
        s = 0.0_dp
        !$omp parallel do schedule(static) reduction(+:s)
        do i = 1, n
            s = s + x(i)*y(i)
        end do
        !$omp end parallel do
    end function pdot

    subroutine paxpy(n,a,x,y)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, x(n)
        real(dp), intent(inout) :: y(n)
        integer :: i
        !$omp parallel do schedule(static)
        do i = 1, n
            y(i) = y(i) + a*x(i)
        end do
        !$omp end parallel do
    end subroutine paxpy

    subroutine pscale(n,a,x,y)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, x(n)
        real(dp), intent(out) :: y(n)
        integer :: i
        !$omp parallel do schedule(static)
        do i = 1, n
            y(i) = a*x(i)
        end do
        !$omp end parallel do
    end subroutine pscale

    subroutine op_free(Op)
        type(op_t), intent(inout) :: Op
        integer :: t
        if (allocated(Op%u)) deallocate(Op%u)
        if (allocated(Op%W)) deallocate(Op%W)
        if (allocated(Op%uT)) deallocate(Op%uT)
        do t = 1, 3
            if (allocated(Op%T(t)%B)) deallocate(Op%T(t)%B)
        end do
    end subroutine op_free

    subroutine eng_free(Eng)
        type(eng_t), intent(inout) :: Eng
        if (allocated(Eng%Mp)) deallocate(Eng%Mp)
        if (allocated(Eng%Mn)) deallocate(Eng%Mn)
        if (allocated(Eng%Cnt)) deallocate(Eng%Cnt)
        if (allocated(Eng%ChT)) deallocate(Eng%ChT, Eng%ChS, Eng%ChP, Eng%ChE)
        if (allocated(Eng%Ts)) deallocate(Eng%Ts)
        if (allocated(Eng%Yb)) deallocate(Eng%Yb)
        Eng%Stored = .false.
    end subroutine eng_free

    ! Current & peak resident memory (VmRSS, VmHWM of /proc/self/status) ...
    subroutine mem_report(Label)
        character(len=*), intent(in) :: Label
        integer :: u, Ios
        character(len=256) :: Line
        integer(i8) :: kb, Rss, Hwm
        Rss = -1
        Hwm = -1
        open(newunit=u, file='/proc/self/status', status='old', action='read', iostat=Ios)
        if (Ios /= 0) return
        do
            read(u,'(a)',iostat=Ios) Line
            if (Ios /= 0) exit
            if (index(Line,'VmRSS:') == 1) then
                read(Line(7:),*,iostat=Ios) kb
                if (Ios == 0) Rss = kb
            else if (index(Line,'VmHWM:') == 1) then
                read(Line(7:),*,iostat=Ios) kb
                if (Ios == 0) Hwm = kb
            end if
        end do
        close(u)
        if (Rss < 0) return
        write(*,'(a)') TAB // '  memory ... RSS ' // cfmt_f(real(Rss,dp)/976562.5_dp,0,1) // ' GB, peak ' // &
                       cfmt_f(real(Hwm,dp)/976562.5_dp,0,1) // ' GB (' // Label // ')'
    end subroutine mem_report

    function istr8(i) result(s)
        integer(i8), intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: Buf
        write(Buf,'(I0)') i
        s = trim(Buf)
    end function istr8

    ! ---------------------------------------------------------------------------------------------------------------
    ! Hamiltonians in the J-coupled HO storage ...

    subroutine hamiltonian_nn(Params,Orb,V,Orb_NN,uH,C0)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o2b_t), intent(out) :: V
        type(orb2b_t), intent(out) :: Orb_NN
        real(dp), allocatable, intent(out) :: uH(:,:,:)
        real(dp), intent(out) :: C0
        type(o1b_t) :: T
        integer :: a_max
        real(dp) :: f
        call V2b_read(Params,V,Orb_NN)
        T = T1b(Params,Orb)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        select case (trim(Params%Calc%CMS))
        case ('CMS1+2B')
            f = 1.0_dp - 1.0_dp/real(Params%Calc%A,dp)
        case ('CMS2B')
            f = 0.0_dp
        case default
            f = 1.0_dp
        end select
        allocate(uH(a_max,a_max,2))
        uH(:,:,1) = f * T%p(1:a_max,1:a_max)
        uH(:,:,2) = f * T%n(1:a_max,1:a_max)
        C0 = 0.0_dp
        write(*,'(a)') TAB // 'NN Hamiltonian ... one-body (' // cfmt_f(f,0,6) // ') T, V_NN with the CMS ' // &
                       trim(Params%Calc%CMS) // ' 2-body part'
    end subroutine hamiltonian_nn

    ! Plain form of the normal-ordered EMPM Hamiltonian in the HF basis, then back to HO ...
    subroutine hamiltonian_empm(Params,Orb,V,Orb_NN,uH,C0,E_HF)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o2b_t), intent(out) :: V
        type(orb2b_t), intent(out) :: Orb_NN
        real(dp), allocatable, intent(out) :: uH(:,:,:)
        real(dp), intent(out) :: C0, E_HF
        type(o1b_t) :: CHF, Ct, Rho
        real(dp), allocatable :: uF(:,:,:), occ(:,:), Cm(:,:), pSPE(:), nSPE(:)
        real(dp), parameter :: Tol_occ = 1.0e-8_dp
        integer :: a_max, a, b, c, J, t, Pab
        real(dp) :: s
        logical :: Ex

        ! the HF-like reference io/<Path>/<&empm Reference>/bin/{H2b,H1b,C,rho}.bin (stage 5) & E_HF of its HF_summary.dat ...
        if (.not. is_dir(empm_ref_path(Params,'/bin'))) call fatal('NCSM EMPM Hamiltonian ... no reference ' // &
            empm_ref_path(Params,'') // ' ... run ' // trim(Params%Calc%EMPM%Reference) // ' first (&empm Reference)')
        inquire(file=empm_ref_path(Params,'/bin/H2b.bin'), exist=Ex)
        if (.not. Ex) call fatal('NCSM EMPM Hamiltonian ... the reference 2-body Hamiltonian ' // &
                                 empm_ref_path(Params,'/bin/H2b.bin') // ' is missing ... run HF with BMF = .true.')
        ! H1b diagonal (HF-like, not HF-MBPT) first, then the integer occupations ...
        call read_h1b_diag(Params,Orb,empm_ref_path(Params,'/bin/H1b.bin'),pSPE,nSPE)
        Rho = O1b_import(Params,Orb,empm_ref_path(Params,'/bin/rho.bin'))
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(occ(a_max,2), uF(a_max,a_max,2))
        if (allocated(EpsHF)) deallocate(EpsHF)
        allocate(EpsHF(a_max,2))
        do a = 1, a_max
            occ(a,1) = Rho%p(a,a)
            occ(a,2) = Rho%n(a,a)
            EpsHF(a,1) = pSPE(a)
            EpsHF(a,2) = nSPE(a)
        end do
        ! the plain form & the SD reference need a determinant: occupations 0 or 1 within 1e-8 (used as read, as before) ...
        do t = 1, 2
            do a = 1, a_max
                if (occ(a,t) > Tol_occ .and. occ(a,t) < 1.0_dp - Tol_occ) &
                    call fatal('EMPM/NCSM need integer occupations: closed (sub)shell nucleus or &hf ESH=.true. ... ' // &
                               empm_ref_path(Params,'/bin/rho.bin') // ': orbital ' // istr(a) // ' (l = ' // istr(Orb(a)%l) // &
                               ', 2j = ' // istr(Orb(a)%j) // ') has the ' // trim(merge('proton ','neutron',t == 1)) // &
                               ' occupation ' // cfmt_e(occ(a,t),0,10))
            end do
        end do
        call O2b_import(Params,Orb,empm_ref_path(Params,'/bin/H2b.bin'),V,Orb_NN)
        CHF = O1b_import(Params,Orb,empm_ref_path(Params,'/bin/C.bin'))
        E_HF = read_ehf(empm_ref_path(Params,'/HF_summary.dat'))
        ! Window: the frozen core is in E_HF & e (normal order w.r.t. HF) ... H in the window = E_HF + e :a+a: + V :a+a+aa:
            ! normal-ordered w.r.t. the active reference -> plain form with the active holes only ... no window: the
            ! reference = the HF occupations ...
        if (UseWin) then
            do t = 1, 2
                do a = 1, a_max
                    occ(a,t) = merge(1.0_dp,0.0_dp,RefO(a,t))
                end do
            end do
        else if (Params%Calc%NCSM%Basis == 'HF') then
            allocate(RefO(a_max,2), ActO(a_max,2))
            RefO = occ(1:a_max,:) > 0.5_dp
            ActO = .true.
        end if

        ! u = e - sum_c n_c sum_J (2J+1)/(2j_a+1) V^J(ac,bc), C0 = -sum (2j+1) n e + 1/2 sum n n (2J+1) V^J ...
        uF = 0.0_dp
        C0 = 0.0_dp
        do t = 1, 2
            do a = 1, a_max
                uF(a,a,t) = EpsHF(a,t)
                C0 = C0 - real(Orb(a)%j + 1,dp) * occ(a,t) * uF(a,a,t)
            end do
        end do
        !$omp parallel do schedule(dynamic) private(b,c,J,s,Pab)
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%l /= Orb(b)%l .or. Orb(a)%j /= Orb(b)%j) cycle
                do c = 1, a_max
                    Pab = mod(Orb(a)%l + Orb(c)%l,2) + 1
                    if (occ(c,1) /= 0.0_dp) then
                        s = 0.0_dp
                        do J = abs(Orb(a)%j - Orb(c)%j)/2, (Orb(a)%j + Orb(c)%j)/2
                            s = s + real(2*J + 1,dp) * (vj(1,a,c,b,c,J,Pab) * 1.0_dp)
                        end do
                        uF(a,b,1) = uF(a,b,1) - occ(c,1) * s / real(Orb(a)%j + 1,dp)
                        s = 0.0_dp
                        do J = abs(Orb(a)%j - Orb(c)%j)/2, (Orb(a)%j + Orb(c)%j)/2
                            s = s + real(2*J + 1,dp) * vj(3,c,a,c,b,J,Pab)
                        end do
                        uF(a,b,2) = uF(a,b,2) - occ(c,1) * s / real(Orb(a)%j + 1,dp)
                    end if
                    if (occ(c,2) /= 0.0_dp) then
                        s = 0.0_dp
                        do J = abs(Orb(a)%j - Orb(c)%j)/2, (Orb(a)%j + Orb(c)%j)/2
                            s = s + real(2*J + 1,dp) * vj(2,a,c,b,c,J,Pab)
                        end do
                        uF(a,b,2) = uF(a,b,2) - occ(c,2) * s / real(Orb(a)%j + 1,dp)
                        s = 0.0_dp
                        do J = abs(Orb(a)%j - Orb(c)%j)/2, (Orb(a)%j + Orb(c)%j)/2
                            s = s + real(2*J + 1,dp) * vj(3,a,c,b,c,J,Pab)
                        end do
                        uF(a,b,1) = uF(a,b,1) - occ(c,2) * s / real(Orb(a)%j + 1,dp)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        do a = 1, a_max
            do c = 1, a_max
                Pab = mod(Orb(a)%l + Orb(c)%l,2) + 1
                do J = abs(Orb(a)%j - Orb(c)%j)/2, (Orb(a)%j + Orb(c)%j)/2
                    C0 = C0 + 0.5_dp * real(2*J + 1,dp) * (occ(a,1)*occ(c,1)*vj(1,a,c,a,c,J,Pab) + &
                                                           occ(a,2)*occ(c,2)*vj(2,a,c,a,c,J,Pab)) + &
                         real(2*J + 1,dp) * occ(a,1)*occ(c,2)*vj(3,a,c,a,c,J,Pab)
                end do
            end do
        end do

        ! Basis 'HF' ... H stays in the HF orbitals, C kept for N_cm (b & N transformed to HF) ...
        if (Params%Calc%NCSM%Basis == 'HF') then
            allocate(uH(a_max,a_max,2), Uhf(a_max,a_max,2))
            uH = uF
            Uhf(:,:,1) = CHF%p(1:a_max,1:a_max)
            Uhf(:,:,2) = CHF%n(1:a_max,1:a_max)
            write(*,'(a)') TAB // 'EMPM Hamiltonian ... E_HF = ' // cfmt_f(E_HF,0,8) // ' MeV, C0 = ' // cfmt_f(C0,0,8) // &
                           ' MeV (relative to E_HF), kept in the HF basis (N_cm transformed to HF)'
            C0 = C0 + E_HF
            return
        end if

        ! Back to HO ... u -> C u C^T, V -> transformation with C^T (O2b_transformation maps with C(i,a) on index i) ...
        allocate(uH(a_max,a_max,2), Cm(a_max,a_max))
        Cm = CHF%p(1:a_max,1:a_max)
        uH(:,:,1) = matmul(Cm, matmul(uF(:,:,1), transpose(Cm)))
        Cm = CHF%n(1:a_max,1:a_max)
        uH(:,:,2) = matmul(Cm, matmul(uF(:,:,2), transpose(Cm)))
        Ct%p = transpose(CHF%p)
        Ct%n = transpose(CHF%n)
        call O2b_transformation(Params,Orb,Orb_NN,V,Ct)
        write(*,'(a)') TAB // 'EMPM Hamiltonian ... E_HF = ' // cfmt_f(E_HF,0,8) // ' MeV, C0 = ' // cfmt_f(C0,0,8) // &
                       ' MeV (relative to E_HF), u & H2b transformed to HO'
        C0 = C0 + E_HF

    contains

        real(dp) function vj(kind,a1,a2,a3,a4,J,Pq)
            integer, intent(in) :: kind, a1, a2, a3, a4, J, Pq
            vj = coupled(V,Orb_NN,Orb,kind,a1,a2,a3,a4,J,Pq)
        end function vj

    end subroutine hamiltonian_empm

    real(dp) function read_ehf(File)
        character(len=*), intent(in) :: File
        character(len=512) :: Line
        integer :: u, Ios, k
        read_ehf = 0.0_dp
        open(newunit=u, file=File, status='old', action='read', iostat=Ios)
        if (Ios /= 0) then
            write(*,'(a)') TAB // 'Warning ... ' // File // ' not found, E_HF = 0'
            return
        end if
        do
            read(u,'(a)',iostat=Ios) Line
            if (Ios /= 0) exit
            k = index(Line,'E_HF =')
            if (k > 0 .and. index(Line,'/ A') == 0) read(Line(k+6:),*) read_ehf
        end do
        close(u)
    end function read_ehf

    ! J-coupled matrix element (kind 1 pp, 2 nn, 3 pn with a1, a3 protons) ... zero outside the stored pairs ...
    real(dp) function coupled(V,Orb_NN,Orb,kind,a1,a2,a3,a4,J,Pq)
        type(o2b_t), intent(in) :: V
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: kind, a1, a2, a3, a4, J, Pq
        coupled = 0.0_dp
        if (J > Orb_NN%J_max) return
        if (max(a1,a2,a3,a4) > Orb_NN%a_max) return
        if (kind == 3) then
            if (Orb_NN%Dic(a1,a2,J,0) == 0 .or. Orb_NN%Dic(a3,a4,J,0) == 0) return
            coupled = O2b_pn(a1,a2,a3,a4,J,Pq,V,Orb_NN)
        else
            if (Orb_NN%Dic(max(a1,a2),min(a1,a2),J,1) == 0 .or. Orb_NN%Dic(max(a3,a4),min(a3,a4),J,1) == 0) return
            if (kind == 1) then
                coupled = O2b_pp(a1,a2,a3,a4,J,Pq,V,Orb,Orb_NN)
            else
                coupled = O2b_nn(a1,a2,a3,a4,J,Pq,V,Orb,Orb_NN)
            end if
        end if
    end function coupled

    ! ---------------------------------------------------------------------------------------------------------------
    ! Single-particle m-states & pair blocks ...

    subroutine sp_make(Orb,emax,SPS,Keep)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: emax
        type(sp_t), intent(out) :: SPS
        logical, intent(in), optional :: Keep(:)
        integer :: a, m, k, Pass, i
        do Pass = 1, 2
            k = 0
            do a = 1, size(Orb)
                if (2*Orb(a)%n + Orb(a)%l > emax) cycle
                if (present(Keep)) then
                    if (.not. Keep(a)) cycle
                end if
                do m = -Orb(a)%j, Orb(a)%j, 2
                    k = k + 1
                    if (Pass == 2) then
                        SPS%a(k) = a
                        SPS%m2(k) = m
                        SPS%e(k) = 2*Orb(a)%n + Orb(a)%l
                        SPS%l(k) = Orb(a)%l
                        SPS%j2(k) = Orb(a)%j
                        SPS%n(k) = Orb(a)%n
                    end if
                end do
            end do
            if (Pass == 1) allocate(SPS%a(k), SPS%m2(k), SPS%e(k), SPS%l(k), SPS%j2(k), SPS%n(k))
        end do
        SPS%ns = k
        do k = 2, SPS%ns
            if (SPS%e(k) < SPS%e(k-1)) call fatal('NCSM ... orbitals not ordered by HO quanta')
        end do
        call jump_types(SPS)
        ! Pascal's triangle with an overflow marker (-1) ...
        allocate(SPS%Bin(0:SPS%ns,0:KBIN))
        SPS%Bin = 0
        SPS%Bin(:,0) = 1
        do k = 1, SPS%ns
            do i = 1, KBIN
                if (SPS%Bin(k-1,i-1) < 0 .or. SPS%Bin(k-1,i) < 0) then
                    SPS%Bin(k,i) = -1
                else if (SPS%Bin(k-1,i-1) > huge(1_i8)/KBIN - SPS%Bin(k-1,i)) then
                    SPS%Bin(k,i) = -1
                else
                    SPS%Bin(k,i) = SPS%Bin(k-1,i-1) + SPS%Bin(k-1,i)
                end if
            end do
        end do
    end subroutine sp_make

    ! Jump types (kc, ka) in buckets (2 dm, parity type) ...
    subroutine jump_types(SPS)
        type(sp_t), intent(inout) :: SPS
        integer :: kc, ka, b, d, pt, dn
        SPS%dmax = 2*maxval(SPS%j2)
        SPS%nbk = 2*(SPS%dmax + 1)
        SPS%b0 = SPS%dmax + 1
        allocate(SPS%Bk(SPS%ns,SPS%ns), SPS%Ix(SPS%ns,SPS%ns), SPS%Nty(SPS%nbk), SPS%Bp(SPS%nbk), SPS%Tkc(SPS%nbk), SPS%Tka(SPS%nbk))
        SPS%Nty = 0
        do ka = 1, SPS%ns
            do kc = 1, SPS%ns
                d = (SPS%m2(kc) - SPS%m2(ka) + SPS%dmax)/2 + 1
                pt = mod(SPS%l(kc) + SPS%l(ka),2)
                b = 2*(d - 1) + pt + 1
                SPS%Bk(kc,ka) = b
                SPS%Nty(b) = SPS%Nty(b) + 1
                SPS%Ix(kc,ka) = SPS%Nty(b)
            end do
        end do
        if (maxval(SPS%Nty) > 2*IXOFF - 1) call fatal('NCSM ... more than 65535 jump types in a bucket')
        do b = 1, SPS%nbk
            allocate(SPS%Tkc(b)%v(SPS%Nty(b)), SPS%Tka(b)%v(SPS%Nty(b)))
            d = (b - 1)/2 + 1
            pt = mod(b - 1,2)
            dn = SPS%dmax + 2 - d
            SPS%Bp(b) = 2*(dn - 1) + pt + 1
        end do
        do ka = 1, SPS%ns
            do kc = 1, SPS%ns
                b = SPS%Bk(kc,ka)
                SPS%Tkc(b)%v(SPS%Ix(kc,ka)) = kc
                SPS%Tka(b)%v(SPS%Ix(kc,ka)) = ka
            end do
        end do
    end subroutine jump_types

    subroutine pairs_make(SPS,kind,Ecut,Pr)
        type(sp_t), intent(in) :: SPS
        integer, intent(in) :: kind, Ecut
        type(pairs_t), intent(out) :: Pr
        integer :: k1, k2, Pass, n, b, Mmax, M2, Pq, i
        integer, allocatable :: Bid(:,:), Cnt(:)
        Mmax = 2*maxval(SPS%j2)
        allocate(Bid(-Mmax:Mmax,2), Pr%id(SPS%ns,SPS%ns))
        Bid = 0
        Pr%id = 0
        ! Blocks ...
        Pr%nb = 0
        do k1 = 1, SPS%ns
            do k2 = 1, SPS%ns
                if (.not. ok(k1,k2)) cycle
                M2 = SPS%m2(k1) + SPS%m2(k2)
                Pq = mod(SPS%l(k1) + SPS%l(k2),2) + 1
                if (Bid(M2,Pq) == 0) then
                    Pr%nb = Pr%nb + 1
                    Bid(M2,Pq) = Pr%nb
                end if
            end do
        end do
        allocate(Pr%bM2(Pr%nb), Pr%bP(Pr%nb), Pr%mem(Pr%nb), Cnt(Pr%nb))
        do M2 = -Mmax, Mmax
            do Pq = 1, 2
                if (Bid(M2,Pq) > 0) then
                    Pr%bM2(Bid(M2,Pq)) = M2
                    Pr%bP(Bid(M2,Pq)) = Pq
                end if
            end do
        end do
        do Pass = 1, 2
            n = 0
            Cnt = 0
            do k1 = 1, SPS%ns
                do k2 = 1, SPS%ns
                    if (.not. ok(k1,k2)) cycle
                    n = n + 1
                    b = Bid(SPS%m2(k1) + SPS%m2(k2), mod(SPS%l(k1) + SPS%l(k2),2) + 1)
                    Cnt(b) = Cnt(b) + 1
                    if (Pass == 2) then
                        Pr%k1(n) = k1
                        Pr%k2(n) = k2
                        Pr%blk(n) = b
                        Pr%loc(n) = Cnt(b)
                        Pr%id(k1,k2) = n
                        Pr%mem(b)%v(Cnt(b)) = n
                    end if
                end do
            end do
            if (Pass == 1) then
                allocate(Pr%k1(n), Pr%k2(n), Pr%blk(n), Pr%loc(n))
                do i = 1, Pr%nb
                    allocate(Pr%mem(i)%v(Cnt(i)))
                end do
            end if
        end do
        Pr%np = n
    contains
        logical function ok(x,y)
            integer, intent(in) :: x, y
            ok = (SPS%e(x) + SPS%e(y) <= Ecut)
            if (kind /= 3) ok = ok .and. (x < y)
        end function ok
    end subroutine pairs_make

    ! ---------------------------------------------------------------------------------------------------------------
    ! Operators: Hamiltonian from the J-coupled storage, J^2 & N_cm ...

    subroutine op_from_coupled(SPS,Orb,Pr,V,Orb_NN,uH,C0,Op)
        type(sp_t), intent(in) :: SPS
        type(orb1b_t), intent(in) :: Orb(:)
        type(pairs_t), intent(in) :: Pr(3)
        type(o2b_t), intent(in) :: V
        type(orb2b_t), intent(in) :: Orb_NN
        real(dp), intent(in) :: uH(:,:,:), C0
        type(op_t), intent(out) :: Op
        integer :: k1, k2, t, b, i, j, n, p1, p2, J2, Jmin, Jmax, M2, Pq, a, bb, c, d
        real(dp) :: s
        Op%C0 = C0
        allocate(Op%u(SPS%ns,SPS%ns,2))
        Op%u = 0.0_dp
        do k1 = 1, SPS%ns
            do k2 = 1, SPS%ns
                if (SPS%m2(k1) == SPS%m2(k2) .and. SPS%l(k1) == SPS%l(k2) .and. SPS%j2(k1) == SPS%j2(k2)) then
                    Op%u(k1,k2,:) = uH(SPS%a(k1),SPS%a(k2),:)
                end if
            end do
        end do
        do t = 1, 3
            allocate(Op%T(t)%B(Pr(t)%nb))
            do b = 1, Pr(t)%nb
                n = size(Pr(t)%mem(b)%v)
                allocate(Op%T(t)%B(b)%v(n,n))
                M2 = Pr(t)%bM2(b)
                Pq = Pr(t)%bP(b)
                !$omp parallel do schedule(dynamic) private(i,p1,p2,a,bb,c,d,Jmin,Jmax,s,J2)
                do j = 1, n
                    p2 = Pr(t)%mem(b)%v(j)
                    c = Pr(t)%k1(p2)
                    d = Pr(t)%k2(p2)
                    do i = 1, j
                        p1 = Pr(t)%mem(b)%v(i)
                        a = Pr(t)%k1(p1)
                        bb = Pr(t)%k2(p1)
                        Jmin = max(abs(SPS%j2(a) - SPS%j2(bb)), abs(SPS%j2(c) - SPS%j2(d)), abs(M2))
                        Jmax = min(SPS%j2(a) + SPS%j2(bb), SPS%j2(c) + SPS%j2(d))
                        s = 0.0_dp
                        do J2 = Jmin, Jmax, 2
                            if (mod(J2,2) /= 0) cycle
                            s = s + fCG(SPS%j2(a),SPS%j2(bb),J2,SPS%m2(a),SPS%m2(bb),M2) * fCG(SPS%j2(c),SPS%j2(d),J2,SPS%m2(c),SPS%m2(d),M2) * &
                                    coupled(V,Orb_NN,Orb,t,SPS%a(a),SPS%a(bb),SPS%a(c),SPS%a(d),J2/2,Pq)
                        end do
                        Op%T(t)%B(b)%v(i,j) = s
                        Op%T(t)%B(b)%v(j,i) = s
                    end do
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine op_from_coupled

    ! J^2 = sum j(j+1) + sum_(i/=j) j_i.j_j ... 2-body W = 2 j1.j2, like v_A = W - W_exch, pn W ...
    subroutine op_j2(SPS,Pr,Op)
        type(sp_t), intent(in) :: SPS
        type(pairs_t), intent(in) :: Pr(3)
        type(op_t), intent(out) :: Op
        integer :: k
        allocate(Op%u(SPS%ns,SPS%ns,2))
        Op%u = 0.0_dp
        do k = 1, SPS%ns
            Op%u(k,k,:) = 0.25_dp * real(SPS%j2(k)*(SPS%j2(k) + 2),dp)
        end do
        call op_two_body(SPS,Pr,Op,1)
    end subroutine op_j2

    ! N_cm = 1/A (sum_i N_i + sum_(i/=j) b+(i).b(j)) ... 2-body W = (b+(1).b(2) + b+(2).b(1))/A ...
    subroutine op_ncm(SPS,Pr,A,Op)
        type(sp_t), intent(in) :: SPS
        type(pairs_t), intent(in) :: Pr(3)
        integer, intent(in) :: A
        type(op_t), intent(out) :: Op
        integer :: k, kp, t, i
        real(dp) :: s
        allocate(Op%u(SPS%ns,SPS%ns,2))
        Op%u = 0.0_dp
        if (allocated(Uhf)) then
            ! N in the HF orbitals ... N(k,k') = sum_i C(i,a_k) e_i C(i,a_k') (same l, j, m) ...
            do t = 1, 2
                do k = 1, SPS%ns
                    do kp = 1, SPS%ns
                        if (SPS%m2(k) /= SPS%m2(kp) .or. SPS%l(k) /= SPS%l(kp) .or. SPS%j2(k) /= SPS%j2(kp)) cycle
                        s = 0.0_dp
                        do i = 1, SPS%ns
                            if (SPS%m2(i) /= SPS%m2(k) .or. SPS%l(i) /= SPS%l(k) .or. SPS%j2(i) /= SPS%j2(k)) cycle
                            s = s + Uhf(SPS%a(i),SPS%a(k),t) * real(SPS%e(i),dp) * Uhf(SPS%a(i),SPS%a(kp),t)
                        end do
                        Op%u(k,kp,t) = s / real(A,dp)
                    end do
                end do
            end do
        else
            do k = 1, SPS%ns
                Op%u(k,k,:) = real(SPS%e(k),dp) / real(A,dp)
            end do
        end if
        call op_two_body(SPS,Pr,Op,2,A)
    end subroutine op_ncm

    ! 2-body tables of J^2 (Which = 1) or N_cm (Which = 2) from single-particle matrix elements ...
    subroutine op_two_body(SPS,Pr,Op,Which,A)
        type(sp_t), intent(in) :: SPS
        type(pairs_t), intent(in) :: Pr(3)
        type(op_t), intent(inout) :: Op
        integer, intent(in) :: Which
        integer, intent(in), optional :: A
        real(dp), allocatable :: Bm(:,:,:), Jp(:,:), Jm(:,:), Bs(:,:,:,:), Cs(:,:)
        integer :: t, b, n, i, j, p1, p2, k1, k2, k3, k4, mu, s1, s2, sp
        real(dp) :: s, w
        if (Which == 1) then
            allocate(Jp(SPS%ns,SPS%ns), Jm(SPS%ns,SPS%ns))
            Jp = 0.0_dp
            do k1 = 1, SPS%ns
                do k2 = 1, SPS%ns
                    if (SPS%a(k1) == SPS%a(k2) .and. SPS%m2(k1) == SPS%m2(k2) + 2) &
                        Jp(k1,k2) = 0.5_dp*sqrt(real(SPS%j2(k2)*(SPS%j2(k2) + 2) - SPS%m2(k2)*(SPS%m2(k2) + 2),dp))
                end do
            end do
            Jm = transpose(Jp)
            w = 1.0_dp
        else
            call b_mscheme(SPS,Bm)
            w = 1.0_dp / real(A,dp)
            ! species copies, in the HF orbitals b -> C^T b C (C(i,k) = C(a_i,a_k) of the reference for the same l, j, m) ...
            allocate(Bs(SPS%ns,SPS%ns,-1:1,2))
            do sp = 1, 2
                Bs(:,:,:,sp) = Bm
                if (.not. allocated(Uhf)) cycle
                allocate(Cs(SPS%ns,SPS%ns))
                Cs = 0.0_dp
                do k1 = 1, SPS%ns
                    do k2 = 1, SPS%ns
                        if (SPS%m2(k1) == SPS%m2(k2) .and. SPS%l(k1) == SPS%l(k2) .and. SPS%j2(k1) == SPS%j2(k2)) &
                            Cs(k1,k2) = Uhf(SPS%a(k1),SPS%a(k2),sp)
                    end do
                end do
                do mu = -1, 1
                    Bs(:,:,mu,sp) = matmul(transpose(Cs), matmul(Bm(:,:,mu), Cs))
                end do
                deallocate(Cs)
            end do
        end if
        do t = 1, 3
            ! species of particles 1 & 2 (pp, nn, pn) ...
            s1 = merge(2,1,t == 2)
            s2 = merge(1,2,t == 1)
            allocate(Op%T(t)%B(Pr(t)%nb))
            do b = 1, Pr(t)%nb
                n = size(Pr(t)%mem(b)%v)
                allocate(Op%T(t)%B(b)%v(n,n))
                !$omp parallel do schedule(dynamic) private(i,p1,p2,k1,k2,k3,k4,s)
                do j = 1, n
                    p2 = Pr(t)%mem(b)%v(j)
                    k3 = Pr(t)%k1(p2)
                    k4 = Pr(t)%k2(p2)
                    do i = 1, n
                        p1 = Pr(t)%mem(b)%v(i)
                        k1 = Pr(t)%k1(p1)
                        k2 = Pr(t)%k2(p1)
                        s = W2(k1,k2,k3,k4)
                        if (t /= 3) s = s - W2(k1,k2,k4,k3)
                        Op%T(t)%B(b)%v(i,j) = w * s
                    end do
                end do
                !$omp end parallel do
            end do
        end do
    contains
        ! Direct 2-body matrix element <k1 k2|W|k3 k4> ...
        real(dp) function W2(x1,x2,x3,x4)
            integer, intent(in) :: x1, x2, x3, x4
            integer :: m
            W2 = 0.0_dp
            if (Which == 1) then
                if (x1 == x3 .and. x2 == x4) W2 = 0.25_dp * real(SPS%m2(x1)*SPS%m2(x2),dp)
                W2 = W2 + 0.5_dp*(Jp(x1,x3)*Jm(x2,x4) + Jm(x1,x3)*Jp(x2,x4))
                W2 = 2.0_dp * W2
            else
                do m = -1, 1
                    W2 = W2 + Bs(x3,x1,m,s1)*Bs(x2,x4,m,s2) + Bs(x1,x3,m,s1)*Bs(x4,x2,m,s2)
                end do
            end if
        end function W2
    end subroutine op_two_body

    ! <k|b_mu|k'> = sqrt(2) <k|r_mu/b|k'> for e_k = e_k' - 1 ... radial <n'l'|r|nl> analytic (positive at the origin),
        ! angular part by the (l s) j uncoupling & <l'm'|sqrt(4pi/3) Y_1mu|lm> = sqrt((2l+1)/(2l'+1)) <l0 10|l'0> <lm 1mu|l'm'> ...
    subroutine b_mscheme(SPS,Bm)
        type(sp_t), intent(in) :: SPS
        real(dp), allocatable, intent(out) :: Bm(:,:,:)
        integer :: k, kp, mu, ms, ml, mlp
        real(dp) :: rad, ang
        allocate(Bm(SPS%ns,SPS%ns,-1:1))
        Bm = 0.0_dp
        do kp = 1, SPS%ns
            do k = 1, SPS%ns
                if (SPS%e(kp) /= SPS%e(k) - 1 .or. abs(SPS%l(kp) - SPS%l(k)) /= 1) cycle
                rad = radial(SPS%n(kp),SPS%l(kp),SPS%n(k),SPS%l(k))
                if (rad == 0.0_dp) cycle
                do mu = -1, 1
                    if (SPS%m2(kp) /= SPS%m2(k) + 2*mu) cycle
                    ang = 0.0_dp
                    do ms = -1, 1, 2
                        mlp = SPS%m2(kp) - ms
                        ml = SPS%m2(k) - ms
                        if (abs(mlp) > 2*SPS%l(kp) .or. abs(ml) > 2*SPS%l(k)) cycle
                        ang = ang + fCG(2*SPS%l(kp),1,SPS%j2(kp),mlp,ms,SPS%m2(kp)) * fCG(2*SPS%l(k),1,SPS%j2(k),ml,ms,SPS%m2(k)) * &
                              sqrt(real(2*SPS%l(k) + 1,dp)/real(2*SPS%l(kp) + 1,dp)) * fCG(2*SPS%l(k),2,2*SPS%l(kp),0,0,0) * &
                              fCG(2*SPS%l(k),2,2*SPS%l(kp),ml,2*mu,mlp)
                    end do
                    Bm(kp,k,mu) = sqrt(2.0_dp) * rad * ang
                end do
            end do
        end do
    contains
        real(dp) function radial(n1,l1,n,l)
            integer, intent(in) :: n1, l1, n, l
            radial = 0.0_dp
            if (l1 == l + 1) then
                if (n1 == n) radial = sqrt(real(n + l,dp) + 1.5_dp)
                if (n1 == n - 1) radial = -sqrt(real(n,dp))
            else if (l1 == l - 1) then
                if (n1 == n) radial = sqrt(real(n + l,dp) + 0.5_dp)
                if (n1 == n + 1) radial = -sqrt(real(n + 1,dp))
            end if
        end function radial
    end subroutine b_mscheme

    ! ---------------------------------------------------------------------------------------------------------------
    ! Slater determinants, one-body jumps & the many-body basis ...

    ! Combinatorial rank of the ascending occupations oc (distinct, 1..ns) ... the SD hash key ...
    ! Hole-particle form (SDX%HP): Off(h) + rank(holes) C(npar,h) + rank(particles), -1 for an occupation outside the
        ! species list (inactive m-state, more than Hk holes) ... keys -1 must not be looked up (the empty marker) ...
    pure integer(i8) function sd_key(SPS,SDX,oc)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: SDX
        integer, intent(in) :: oc(:)
        integer :: ii, k, jr, nh, npp
        integer(i8) :: rh, rp
        sd_key = 0
        if (.not. SDX%HP) then
            do ii = 1, size(oc)
                sd_key = sd_key + SPS%Bin(oc(ii) - 1,ii)
            end do
            return
        end if
        sd_key = -1
        nh = 0
        npp = 0
        rh = 0
        rp = 0
        jr = 1
        do ii = 1, size(oc)
            k = oc(ii)
            if (SDX%Rr(k) > 0) then
                do while (SDX%RefL(jr) < k)
                    nh = nh + 1
                    if (nh > SDX%Hk) return
                    rh = rh + SDX%Ck(jr - 1,nh)
                    jr = jr + 1
                end do
                jr = jr + 1
            else if (SDX%Rp(k) > 0) then
                npp = npp + 1
                if (npp > SDX%Hk) return
                rp = rp + SDX%Ck(SDX%Rp(k) - 1,npp)
            else
                return
            end if
        end do
        do while (jr <= SDX%nref)
            nh = nh + 1
            if (nh > SDX%Hk) return
            rh = rh + SDX%Ck(jr - 1,nh)
            jr = jr + 1
        end do
        if (nh /= npp) return
        sd_key = SDX%Off(nh) + rh*SDX%Ck(SDX%npar,nh) + rp
    end function sd_key

    ! Target SD index of a+_kc a_(occ(ii)) |ss> & the sign xs ... -1 if outside the basis ...
    ! No array temporaries (called in threaded loops; a module procedure,
    ! as an internal one ifx rebuilt descriptors of the host's SDX/SPS components in static memory on every call): q from SDX%q, nb1 = ii - 1 (sorted), target list merged in tmp ...
    integer function sd_tgt(SPS,SDX,Nmax,ss,ii,kc_,xs,tmp)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: SDX
        integer, intent(in) :: Nmax, ss, ii, kc_
        integer, intent(out) :: xs
        integer, intent(inout) :: tmp(:)   ! work array (thread-private) ...
        integer :: nb2, ka_, a, b, kk
        integer(i8) :: kx
        logical :: ins
        ka_ = SDX%occ(ii,ss)
        if (kc_ == ka_) then
            sd_tgt = ss
            xs = 1
            return
        end if
        if (SDX%q(ss) - SPS%e(ka_) + SPS%e(kc_) > Nmax) then
            sd_tgt = -1
            return
        end if
        ! occupations without ka_, kc_ inserted in order; nb2 = number of the others below kc_ ...
        b = 0
        nb2 = 0
        ins = .false.
        do a = 1, SDX%Nt
            if (a == ii) cycle
            kk = SDX%occ(a,ss)
            if (.not. ins .and. kc_ < kk) then
                b = b + 1
                tmp(b) = kc_
                ins = .true.
            end if
            if (kk < kc_) nb2 = nb2 + 1
            b = b + 1
            tmp(b) = kk
        end do
        if (.not. ins) tmp(SDX%Nt) = kc_
        kx = sd_key(SPS,SDX,tmp)
        sd_tgt = -1
        if (kx >= 0) sd_tgt = hashmap_lookup(SDX%Map,kx)
        xs = 1 - 2*mod(ii - 1 + nb2,2)
    end function sd_tgt

    ! Eps (optional, energy cut): HF single-particle energies of the m-states of this species, Exmax(h) the largest
        ! excitation energy sum(Eps(occ)) - sum(Eps(1:Nt)) kept for SDs with h particles outside the reference ...
    subroutine sd_make(SPS,Nt,Nmax,Hcap,SDX,Eps,Exmax,Ref,Act,KeepM,KeepP)
        type(sp_t), intent(in) :: SPS
        integer, intent(in) :: Nt, Nmax, Hcap
        type(sd_t), intent(out) :: SDX
        real(dp), intent(in), optional :: Eps(:), Exmax(0:)
        logical, intent(in), optional :: Ref(:), Act(:)
        integer, intent(in), optional :: KeepM(:), KeepP
        logical :: Cut, isr, isa
        integer :: nr, npa, hcur, kk2
        integer, allocatable :: hsel(:), psel(:), Ps(:)
        real(dp) :: plim
        real(dp), allocatable :: ep(:)
        real(dp) :: Eref, exv, eminp
        real(dp), allocatable :: epl(:), exm(:), cum(:)
        integer, allocatable :: o(:), Cnt(:)
        integer :: Pass, n, s, i, kc, ka, d, sp, x
        integer(i8) :: j, Tot, tt
        integer, allocatable :: tmp(:)
        logical, allocatable :: Occ(:)
        ! jump pass: chunk buffers, thread buffers (target, type, bucket), chunk offsets ...
        type(jbuf_t), allocatable :: Bf(:)
        integer, allocatable :: Lt(:), t2(:)
        integer(i2), allocatable :: Lx(:), Ld(:), x2(:), d2(:)
        integer(i8), allocatable :: CO(:)
        integer(i8) :: k8, p8, nj
        integer :: ic, nch, m, hh, mt, qt, KP
        logical :: KM
        integer, parameter :: JCH = 1024
        call timer_start(tt)
        SDX%Nt = Nt
        SDX%HP = Trunc == 'npnh'
        allocate(SDX%Pt(SPS%ns), SDX%Ac(SPS%ns))
        if (SDX%HP) then
            ! hole-particle form: reference & particle lists, binomials up to Hk = min(Hcap, nref, npar), offsets per h ...
            allocate(SDX%Rr(SPS%ns), SDX%Rp(SPS%ns))
            SDX%Rr = 0
            SDX%Rp = 0
            nr = 0
            npa = 0
            do i = 1, SPS%ns
                isr = i <= Nt
                if (present(Ref)) isr = Ref(i)
                isa = .true.
                if (present(Act)) isa = Act(i)
                if (isr) then
                    nr = nr + 1
                    SDX%Rr(i) = nr
                else if (isa) then
                    npa = npa + 1
                    SDX%Rp(i) = npa
                end if
            end do
            if (nr /= Nt) call fatal('NCSM ... reference m-states ' // istr(nr) // ' /= particle number ' // istr(Nt))
            SDX%nref = nr
            SDX%npar = npa
            SDX%Hk = min(Hcap,nr,npa)
            allocate(SDX%RefL(nr), SDX%ParL(npa))
            do i = 1, SPS%ns
                if (SDX%Rr(i) > 0) SDX%RefL(SDX%Rr(i)) = i
                if (SDX%Rp(i) > 0) SDX%ParL(SDX%Rp(i)) = i
            end do
            allocate(SDX%Ck(0:max(nr,npa),0:SDX%Hk), SDX%Off(0:SDX%Hk+1))
            SDX%Ck = 0
            SDX%Ck(:,0) = 1
            do i = 1, max(nr,npa)
                do x = 1, SDX%Hk
                    if (real(SDX%Ck(i-1,x-1),dp) + real(SDX%Ck(i-1,x),dp) > 0.5_dp*real(huge(1_i8),dp)) &
                        call fatal('NCSM ... hole-particle keys overflow')
                    SDX%Ck(i,x) = SDX%Ck(i-1,x-1) + SDX%Ck(i-1,x)
                end do
            end do
            SDX%Off(0) = 0
            do x = 0, SDX%Hk
                if (real(SDX%Off(x),dp) + real(SDX%Ck(nr,x),dp)*real(SDX%Ck(npa,x),dp) > 0.5_dp*real(huge(1_i8),dp)) &
                    call fatal('NCSM ... hole-particle keys overflow')
                SDX%Off(x+1) = SDX%Off(x) + SDX%Ck(nr,x)*SDX%Ck(npa,x)
            end do
            SDX%Pt = SDX%Rp > 0
            SDX%Ac = SDX%Rr > 0 .or. SDX%Rp > 0
            SDX%N0 = sum(SPS%e(SDX%RefL))
        else
            do i = 1, SPS%ns
                SDX%Pt(i) = i > Nt
            end do
            SDX%Ac = .true.
            SDX%N0 = sum(SPS%e(1:Nt))
            ! keys = combinatorial rank of the sorted occupations, 0 .. C(ns,Nt) - 1 ...
            if (Nt > KBIN) call fatal('NCSM ... SD keys overflow (more than KBIN particles of a species)')
            if (SPS%Bin(SPS%ns,Nt) < 0) call fatal('NCSM ... SD keys overflow (C(ns,Nt) too large)')
        end if
        allocate(o(Nt))
        Cut = present(Eps)
        Eref = 0.0_dp
        if (Cut) then
            epl = Eps
            exm = Exmax
            if (SDX%HP) then
                Eref = sum(epl(SDX%RefL))
            else
                Eref = sum(epl(1:Nt))
            end if
            ! early pruning: cum(k) = sum of the reference energies e(1:k), eminp the lowest particle energy ...
            allocate(cum(0:Nt))
            cum(0) = 0.0_dp
            do i = 1, Nt
                cum(i) = cum(i-1) + epl(i)
            end do
            eminp = huge(1.0_dp)
            if (Nt < SPS%ns) eminp = minval(epl(Nt+1:))
        end if
        if (SDX%HP) then
            ! particles in ascending energy (with the cut, for pruning), else in list order ...
            allocate(hsel(max(SDX%Hk,1)), psel(max(SDX%Hk,1)), Ps(npa), ep(npa))
            do i = 1, npa
                Ps(i) = i
                ep(i) = 0.0_dp
                if (Cut) ep(i) = epl(SDX%ParL(i))
            end do
            if (Cut .and. npa > 0) call argsort_r(ep,Ps)
        end if
        do Pass = 1, 2
            n = 0
            if (SDX%HP) then
                do hcur = 0, SDX%Hk
                    call hrec(1,1)
                end do
            else
                call rec(1,1,0)
            end if
            if (Pass == 1) allocate(SDX%occ(Nt,n), SDX%q(n), SDX%m2(n), SDX%h(n), SDX%L(n))
            if (Pass == 1 .and. Cut) allocate(SDX%ex(n))
        end do
        SDX%nsd = n
        if (SDX%HP) then
            allocate(SDX%He(0:SDX%Hk))
            do x = 0, SDX%Hk
                SDX%He(x) = count(SDX%h <= x)
            end do
            do s = 2, n
                if (SDX%h(s) < SDX%h(s-1)) call fatal('NCSM ... SDs not in ascending h')
            end do
        end if
        call timer_report(tt,'  SDs: generation')
        call timer_start(tt)
        call SDX%Map%init(n)
        do s = 1, n
            call SDX%Map%set(sd_key(SPS,SDX,SDX%occ(:,s)), s)
        end do
        call timer_report(tt,'  SDs: hash map')
        call timer_start(tt)

        ! One-body jumps ... a+_kc a_ka |s> for occupied ka & kc empty or = ka ... one pass over chunks of JCH SDs
            ! (dynamic): the jumps of a chunk (signed target, type, bucket) into its buffer Bf(ic) & the counts per bucket
            ! into Jo(b+1,s), prefix sums (chunk totals, then per chunk), the buffers copied into Jt, Jx by bucket & each
            ! bucket sorted by target ... candidates outside the list are skipped before the lookup: h' > Hk, h' = Hcap
            ! with 2M or parity not kept (KeepM, KeepP, as emit), with the energy cut above Exmax(h') ...
        SDX%nd = SPS%nbk
        if (SDX%nd > huge(1_i2)) call fatal('NCSM ... too many jump buckets')
        allocate(SDX%Jo(SDX%nd + 1,n))
        SDX%Jo = 0
        KM = present(KeepM)
        KP = -1
        if (present(KeepP)) KP = KeepP
        nch = (n + JCH - 1)/JCH
        allocate(Bf(nch), CO(0:nch))
        !$omp parallel private(ic,s,i,ka,kc,d,x,sp,Cnt,tmp,Occ,exv,m,Lt,Lx,Ld,t2,x2,d2,hh,mt,qt)
        allocate(Cnt(SDX%nd), tmp(Nt), Occ(SPS%ns), Lt(64*JCH), Lx(64*JCH), Ld(64*JCH))
        Occ = .false.
        !$omp do schedule(dynamic,1)
        do ic = 1, nch
            m = 0
            do s = (ic - 1)*JCH + 1, min(n,ic*JCH)
                Cnt = 0
                do i = 1, Nt
                    Occ(SDX%occ(i,s)) = .true.
                end do
                do i = 1, Nt
                    ka = SDX%occ(i,s)
                    do kc = 1, SPS%ns
                        if (kc /= ka) then
                            if (Occ(kc)) cycle
                            if (.not. SDX%Ac(kc)) cycle
                            if (SDX%HP) then
                                hh = SDX%h(s) - merge(1,0,SDX%Pt(ka)) + merge(1,0,SDX%Pt(kc))
                                if (hh > SDX%Hk) cycle
                                if (hh == Hcap .and. hh > 0) then
                                    if (KM) then
                                        mt = SDX%m2(s) - SPS%m2(ka) + SPS%m2(kc)
                                        if (all(KeepM /= mt)) cycle
                                    end if
                                    if (KP >= 0) then
                                        qt = SDX%q(s) - SPS%e(ka) + SPS%e(kc)
                                        if (modulo(qt,2) /= KP) cycle
                                    end if
                                end if
                            end if
                            if (Cut) then
                                if (.not. inlist(s,ka,kc)) cycle
                            end if
                        end if
                        sp = sd_tgt(SPS,SDX,Nmax,s,i,kc,x,tmp)
                        if (sp < 0) cycle
                        d = SPS%Bk(kc,ka)
                        Cnt(d) = Cnt(d) + 1
                        m = m + 1
                        if (m > size(Lt)) then
                            ! thread buffers doubled (inline: an internal procedure would see the shared, not the private copies) ...
                            allocate(t2(2*size(Lt)), x2(2*size(Lt)), d2(2*size(Lt)))
                            t2(1:m-1) = Lt(1:m-1)
                            x2(1:m-1) = Lx(1:m-1)
                            d2(1:m-1) = Ld(1:m-1)
                            call move_alloc(t2,Lt)
                            call move_alloc(x2,Lx)
                            call move_alloc(d2,Ld)
                        end if
                        Lt(m) = x*sp
                        Lx(m) = int(SPS%Ix(kc,ka) - IXOFF,i2)
                        Ld(m) = int(d,i2)
                    end do
                end do
                do i = 1, Nt
                    Occ(SDX%occ(i,s)) = .false.
                end do
                do d = 1, SDX%nd
                    SDX%Jo(d+1,s) = int(Cnt(d),i8)
                end do
            end do
            allocate(Bf(ic)%t(m), Bf(ic)%x(m), Bf(ic)%d(m))
            Bf(ic)%t = Lt(1:m)
            Bf(ic)%x = Lx(1:m)
            Bf(ic)%d = Ld(1:m)
        end do
        !$omp end do
        deallocate(Cnt, tmp, Occ, Lt, Lx, Ld)
        !$omp end parallel
        CO(0) = 0
        do ic = 1, nch
            CO(ic) = CO(ic-1) + size(Bf(ic)%t,kind=i8)
        end do
        Tot = CO(nch)
        call timer_report(tt,'  SDs: jumps (one pass, ' // istr8(Tot) // ' jumps)')
        call timer_start(tt)
        allocate(SDX%Jt(Tot), SDX%Jx(Tot))
        !$omp parallel private(ic,s,d,j,k8,p8,nj,Cnt)
        allocate(Cnt(SDX%nd))
        !$omp do schedule(dynamic,1)
        do ic = 1, nch
            j = CO(ic-1)
            do s = (ic - 1)*JCH + 1, min(n,ic*JCH)
                SDX%Jo(1,s) = j + 1
                do d = 1, SDX%nd
                    SDX%Jo(d+1,s) = SDX%Jo(d,s) + SDX%Jo(d+1,s)
                end do
                j = SDX%Jo(SDX%nd+1,s) - 1
            end do
            p8 = 0
            do s = (ic - 1)*JCH + 1, min(n,ic*JCH)
                Cnt = 0
                nj = SDX%Jo(SDX%nd+1,s) - SDX%Jo(1,s)
                do k8 = p8 + 1, p8 + nj
                    d = int(Bf(ic)%d(k8))
                    j = SDX%Jo(d,s) + Cnt(d)
                    Cnt(d) = Cnt(d) + 1
                    SDX%Jt(j) = Bf(ic)%t(k8)
                    SDX%Jx(j) = Bf(ic)%x(k8)
                end do
                p8 = p8 + nj
                ! each bucket in ascending targets ...
                do d = 1, SDX%nd
                    if (SDX%Jo(d+1,s) - SDX%Jo(d,s) > 1) call jsort(SDX%Jt(SDX%Jo(d,s):SDX%Jo(d+1,s)-1),SDX%Jx(SDX%Jo(d,s):SDX%Jo(d+1,s)-1))
                end do
            end do
            deallocate(Bf(ic)%t, Bf(ic)%x, Bf(ic)%d)
        end do
        !$omp end do
        deallocate(Cnt)
        !$omp end parallel
        call timer_report(tt,'  SDs: jumps (prefix, fill & sort), ' // istr8(Tot) // ' jumps (' // &
                             cfmt_f(6.0_dp*real(Tot,dp)/1.0e9_dp,0,2) // ' GB)')
    contains
        ! Occupations ascending (reference m-states 1..Nt first, then particles) ... with the energy cut every partial
            ! occupation is bounded from below: holes known so far (reference states skipped), particles chosen so far &
            ! the lowest particle energy for each slot still to fill by a particle; pruned if above Exmax(holes) (Exmax
            ! is non-increasing in h & every further particle-hole pair costs >= 0 with an HF gap) ...
        recursive subroutine rec(pos,kmin,qs)
            integer, intent(in) :: pos, kmin, qs
            integer :: kk, nref, nh, ii
            real(dp) :: spart, sref, lb
            if (Cut .and. pos <= Nt) then
                nref = 0
                spart = 0.0_dp
                sref = 0.0_dp
                do ii = 1, pos - 1
                    if (o(ii) <= Nt) then
                        nref = nref + 1
                        sref = sref + epl(o(ii))
                    else
                        spart = spart + epl(o(ii))
                    end if
                end do
            end if
            if (pos > Nt) then
                if (qs - SDX%N0 <= Nmax .and. count(o > Nt) <= Hcap) then
                    if (Cut) then
                        exv = sum(epl(o)) - Eref
                        if (exv > exm(count(o > Nt)) + 1.0e-9_dp) return
                    end if
                    n = n + 1
                    if (Pass == 2) then
                        if (Cut) SDX%ex(n) = exv
                        SDX%occ(:,n) = o
                        SDX%q(n) = qs - SDX%N0
                        SDX%m2(n) = sum(SPS%m2(o))
                        SDX%h(n) = count(o > Nt)
                        SDX%L(n) = SDX%q(n)
                        if (Trunc == 'npnh') SDX%L(n) = 2*SDX%h(n) + mod(SDX%q(n),2)
                    end if
                end if
                return
            end if
            do kk = kmin, SPS%ns - (Nt - pos)
                if (count(o(1:pos-1) > Nt) + merge(1,0,kk > Nt) > Hcap) exit
                ! remaining particles at least at the quanta of kk ...
                if (qs + SPS%e(kk)*(Nt - pos + 1) - SDX%N0 > Nmax) exit
                if (Cut) then
                    if (kk > Nt) then       ! all reference decisions made: holes known, the other slots particles
                        nh = Nt - nref
                        lb = spart + epl(kk) + real(Nt - pos,dp)*eminp - (Eref - sref)
                    else                    ! holes so far: reference states below kk not taken
                        nh = kk - pos
                        lb = real(nh,dp)*eminp - (cum(kk-1) - sref)
                    end if
                    if (lb > exm(min(nh,ubound(exm,1))) + 1.0e-9_dp) cycle
                end if
                o(pos) = kk
                call rec(pos + 1,kk + 1,qs + SPS%e(kk))
            end do
        end subroutine rec



        ! Hole-particle generation: hcur holes (ranks hsel, lexicographic) x hcur particles (positions psel in Ps) ... with
            ! the cut sum(e_particles) <= Exmax(h) + sum(e_holes), Ps ascending in energy: pruned by the lower bound ...
        recursive subroutine hrec(d,st)
            integer, intent(in) :: d, st
            integer :: ii
            if (d > hcur) then
                plim = huge(1.0_dp)
                if (Cut) then
                    plim = exm(hcur)
                    do ii = 1, hcur
                        plim = plim + epl(SDX%RefL(hsel(ii)))
                    end do
                end if
                call prec(1,1,0.0_dp)
                return
            end if
            do ii = st, nr - (hcur - d)
                hsel(d) = ii
                call hrec(d + 1,ii + 1)
            end do
        end subroutine hrec

        recursive subroutine prec(d,st,es)
            integer, intent(in) :: d, st
            real(dp), intent(in) :: es
            integer :: ii
            if (d > hcur) then
                call emit()
                return
            end if
            do ii = st, npa - (hcur - d)
                if (Cut) then
                    if (es + real(hcur - d + 1,dp)*ep(Ps(ii)) > plim + 1.0e-9_dp) exit
                end if
                psel(d) = ii
                call prec(d + 1,ii + 1,es + ep(Ps(ii)))
            end do
        end subroutine prec

        subroutine emit()
            integer :: ii, jj, nn, qv
            nn = 0
            jj = 1
            do ii = 1, nr
                if (jj <= hcur) then
                    if (hsel(jj) == ii) then
                        jj = jj + 1
                        cycle
                    end if
                end if
                nn = nn + 1
                o(nn) = SDX%RefL(ii)
            end do
            do ii = 1, hcur
                nn = nn + 1
                o(nn) = SDX%ParL(Ps(psel(ii)))
            end do
            call isort(o)
            if (Cut) then
                exv = sum(epl(o)) - Eref
                if (exv > exm(hcur) + 1.0e-9_dp) return
            end if
            qv = sum(SPS%e(o)) - SDX%N0
            ! h = Hcap: partner = the other reference only ... 2M & parity of the run ...
            if (hcur == Hcap .and. hcur > 0) then
                if (present(KeepM)) then
                    if (all(KeepM /= sum(SPS%m2(o)))) return
                end if
                if (present(KeepP)) then
                    if (KeepP >= 0 .and. modulo(qv,2) /= KeepP) return
                end if
            end if
            n = n + 1
            if (Pass == 2) then
                if (Cut) SDX%ex(n) = exv
                SDX%occ(:,n) = o
                SDX%q(n) = qv
                SDX%m2(n) = sum(SPS%m2(o))
                SDX%h(n) = hcur
                SDX%L(n) = 2*hcur + modulo(qv,2)
            end if
        end subroutine emit

        ! a+_kc a_ka |s> within the energy cut of the list (h' <= Hcap & ex' <= Exmax(h')) ...
        logical function inlist(s_,ka_,kc_)
            integer, intent(in) :: s_, ka_, kc_
            integer :: hh
            hh = SDX%h(s_) - merge(1,0,SDX%Pt(ka_)) + merge(1,0,SDX%Pt(kc_))
            inlist = .false.
            if (hh > Hcap .or. hh > ubound(exm,1)) return
            inlist = SDX%ex(s_) - epl(ka_) + epl(kc_) <= exm(hh) + 1.0e-9_dp
        end function inlist
    end subroutine sd_make

    ! Jumps of a bucket (signed targets t, types x) into ascending |t| (Shell sort) ...
    pure subroutine jsort(t,x)
        integer, intent(inout) :: t(:)
        integer(i2), intent(inout) :: x(:)
        integer :: n, g, i, j, tv
        integer(i2) :: xv
        n = size(t)
        g = 1
        do while (g < n/3)
            g = 3*g + 1
        end do
        do while (g >= 1)
            do i = g + 1, n
                tv = t(i)
                xv = x(i)
                j = i - g
                do while (j >= 1)
                    if (abs(t(j)) <= abs(tv)) exit
                    t(j+g) = t(j)
                    x(j+g) = x(j)
                    j = j - g
                end do
                t(j+g) = tv
                x(j+g) = xv
            end do
            g = g/3
        end do
    end subroutine jsort

    pure subroutine isort(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, v
        do i = 2, size(x)
            v = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= v) exit
                x(j+1) = x(j)
                j = j - 1
            end do
            x(j+1) = v
        end do
    end subroutine isort

    ! Idx = permutation sorting x ascending (stable bottom-up merge sort) ...
    subroutine argsort_r(x,Idx)
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: Idx(:)
        integer, allocatable :: W(:)
        integer :: n, w_, lo, mi, hi, i, j, k
        n = size(x)
        do i = 1, n            ! explicit loops: no array-constructor temporaries on the stack (ifx, > 8 MB) ...
            Idx(i) = i
        end do
        allocate(W(n))
        w_ = 1
        do while (w_ < n)
            do lo = 1, n, 2*w_
                mi = min(lo + w_ - 1, n)
                hi = min(lo + 2*w_ - 1, n)
                if (mi >= hi) then
                    W(lo:hi) = Idx(lo:hi)
                    cycle
                end if
                i = lo
                j = mi + 1
                do k = lo, hi
                    if (j > hi) then
                        W(k) = Idx(i); i = i + 1
                    else if (i > mi) then
                        W(k) = Idx(j); j = j + 1
                    else if (x(Idx(j)) < x(Idx(i))) then
                        W(k) = Idx(j); j = j + 1
                    else
                        W(k) = Idx(i); i = i + 1
                    end if
                end do
            end do
            do i = 1, n
                Idx(i) = W(i)
            end do
            w_ = 2*w_
        end do
    end subroutine argsort_r

    subroutine basis_make(P,Nn,M2,Nmax,ParIn,B)
        type(sd_t), intent(in) :: P, Nn
        integer, intent(in) :: M2, Nmax, ParIn
        type(mb_t), intent(out) :: B
        integer :: Mmax, sp, sn, qn, Mn, c, Par, Lmax, i, f0, nc_, lo, hi, mid
        integer, allocatable :: Cnt(:,:), First(:,:), Lst(:), Ord(:)
        integer(i8) :: dim8
        real(dp) :: thr
        B%M2 = M2
        B%Nmax = Nmax
        Par = mod(Nmax,2)
        if (Trunc == 'npnh') Par = ParIn
        Lmax = max(maxval(Nn%L), maxval(P%L))
        Mmax = max(maxval(abs(Nn%m2)), maxval(abs(P%m2))) + abs(M2)
        ! neutron blocks (M_N, L_N), each in ascending excitation energy with the energy cut (else in SD order), Lc(sn) the
        ! position of sn in its block ...
        allocate(Cnt(-Mmax:Mmax,0:Lmax), B%Lc(Nn%nsd), First(-Mmax:Mmax,0:Lmax), Lst(Nn%nsd), Ord(Nn%nsd))
        Cnt = 0
        do sn = 1, Nn%nsd
            Cnt(Nn%m2(sn),Nn%L(sn)) = Cnt(Nn%m2(sn),Nn%L(sn)) + 1
        end do
        c = 0
        do qn = 0, Lmax
            do Mn = -Mmax, Mmax
                First(Mn,qn) = c
                c = c + Cnt(Mn,qn)
            end do
        end do
        if (UseCut) then
            call argsort_r(Nn%ex,Ord)
        else
            do sn = 1, Nn%nsd
                Ord(sn) = sn
            end do
        end if
        Cnt = 0
        do i = 1, Nn%nsd
            sn = Ord(i)
            Cnt(Nn%m2(sn),Nn%L(sn)) = Cnt(Nn%m2(sn),Nn%L(sn)) + 1
            B%Lc(sn) = Cnt(Nn%m2(sn),Nn%L(sn))
            Lst(First(Nn%m2(sn),Nn%L(sn)) + B%Lc(sn)) = sn
        end do
        allocate(B%Rn(Nn%nsd))
        do i = 1, Nn%nsd
            B%Rn(Lst(i)) = i
        end do
        ! sectors (sp, qn) ... with the energy cut the first Nc of the block (binary search, E_N <= Ecn(h_P + h_N) - E_P) ...
        allocate(B%St(P%nsd,0:Lmax), B%Nc(P%nsd,0:Lmax))
        B%St = -1
        B%Nc = 0
        dim8 = 0
        B%dim_all = 0
        do sp = 1, P%nsd
            Mn = M2 - P%m2(sp)
            if (abs(Mn) > Mmax) cycle
            do qn = 0, Lmax
                if (.not. accept(P%L(sp),qn)) cycle
                nc_ = Cnt(Mn,qn)
                if (nc_ == 0) cycle
                B%dim_all = B%dim_all + nc_
                if (UseCut) then
                    thr = Ecn(P%L(sp)/2 + qn/2) - P%ex(sp) + 1.0e-9_dp
                    f0 = First(Mn,qn)
                    lo = 0
                    hi = nc_
                    do while (lo < hi)
                        mid = (lo + hi + 1)/2
                        if (Nn%ex(Lst(f0 + mid)) <= thr) then
                            lo = mid
                        else
                            hi = mid - 1
                        end if
                    end do
                    nc_ = lo
                    if (nc_ == 0) cycle
                end if
                B%St(sp,qn) = int(dim8)
                B%Nc(sp,qn) = nc_
                dim8 = dim8 + nc_
            end do
        end do
        if (dim8 > huge(1)) call fatal('NCSM ... dimension exceeds 32-bit indices')
        B%dim = int(dim8)
        ! Inverse map ... states by sector ...
        allocate(B%sp(B%dim), B%sn(B%dim))
        do sp = 1, P%nsd
            Mn = M2 - P%m2(sp)
            if (abs(Mn) > Mmax) cycle
            do qn = 0, Lmax
                if (B%St(sp,qn) < 0) cycle
                do c = 1, B%Nc(sp,qn)
                    B%sp(B%St(sp,qn) + c) = sp
                    B%sn(B%St(sp,qn) + c) = Lst(First(Mn,qn) + c)
                end do
            end do
        end do
    contains
        ! Sector (proton label Lp, neutron label Ln) in the basis ... 'Nmax': q_P + q_N <= Nmax & the parity of Nmax,
            ! 'npnh': h_P + h_N <= Nmax (the npnh cut) & parity ParIn ...
        logical function accept(Lp,Ln)
            integer, intent(in) :: Lp, Ln
            if (Trunc == 'npnh') then
                accept = Lp/2 + Ln/2 <= Nmax .and. mod(Lp + Ln,2) == Par
            else
                accept = Lp + Ln <= Nmax .and. mod(Lp + Ln,2) == Par
            end if
        end function accept
    end subroutine basis_make

    ! ---------------------------------------------------------------------------------------------------------------
    ! Sparse matrices ...

    ! Like-particle part of Op on the SDs of species t, lower half: the partners s' with Rank(s') >= Rank(s) ... the
        ! candidates of a pair move (k1 k2) -> (k3 k4) are checked against the SD list before the hash lookup: k3, k4
        ! empty (occupation mask), quanta <= Nmaxq, h' <= Hcap & (energy cut) ex' <= Exmax(h'); with the cut the partner
        ! pairs of a block are taken in ascending energy & the loop ends at the first pair above Exmax(1) ...
        ! values single (vs) or double (v) ...
    subroutine species_matrix(SPS,SDX,Pr,Op,t,Single,Rank,Nmaxq,Hcap,Mt,Eps,Exmax,Skip,Tm,Eps_)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: SDX
        type(pairs_t), intent(in) :: Pr
        type(op_t), intent(in) :: Op
        integer, intent(in) :: t, Nmaxq, Hcap
        logical, intent(in) :: Single
        integer, intent(in) :: Rank(:)
        type(col_t), allocatable, intent(out) :: Mt(:)
        real(dp), intent(in), optional :: Eps(:), Exmax(0:)
        logical, intent(in), optional :: Skip(:)     ! SDs without a column (empty) ...
        ! reduced block of h = Hk (Eps_ = +-1, 0 = none): Tm the signed images ...
        integer, intent(in), optional :: Tm(:), Eps_
        integer :: epr, lim3
        logical :: Rd, S3
        integer, allocatable :: Pos(:), Rl(:), o(:)
        real(dp), allocatable :: Acc(:), Ep(:), Eb(:)
        logical, allocatable :: Occ(:)
        type(ivec_t), allocatable :: Ord(:), Lt(:,:)
        integer, allocatable :: Tq(:), Ty(:)
        logical :: Cut
        real(dp) :: exlim, base, v, sg
        integer :: ty_
        integer :: s, i1_, i2_, k1, k2, k3, k4, pid, b, q, qq, pq, sp, nl, Nt, nb1, nb2, ii, p3, p4, hh, h2, qn, h0, rs, b0, nsel
        integer(i8) :: j, kx
        Nt = SDX%Nt
        Cut = present(Eps)
        allocate(Mt(SDX%nsd))
        epr = 0
        if (present(Eps_)) epr = Eps_
        Rd = epr /= 0
        lim3 = huge(1)
        if (Rd) lim3 = SDX%He(SDX%Hk - 1)
        b0 = SPS%b0
        exlim = huge(1.0_dp)
        base = 0.0_dp
        if (Cut) then
            ! partner pairs of each block in ascending energy ...
            allocate(Ep(Pr%np), Ord(Pr%nb))
            do q = 1, Pr%np
                Ep(q) = Eps(Pr%k1(q)) + Eps(Pr%k2(q))
            end do
            do b = 1, Pr%nb
                nsel = size(Pr%mem(b)%v)
                allocate(Ord(b)%v(nsel), Eb(nsel))
                do q = 1, nsel
                    Eb(q) = Ep(Pr%mem(b)%v(q))
                end do
                call argsort_r(Eb,Ord(b)%v)
                deallocate(Eb)
            end do
            exlim = Exmax(min(1,ubound(Exmax,1)))
        end if
        ! block members split by the number of particle m-states (0, 1, 2) of the pair, inactive ones dropped, in the order
            ! of the scan (ascending energy with the cut) ... only the types within the hole cap are scanned ...
        allocate(Lt(0:2,Pr%nb))
        do b = 1, Pr%nb
            nsel = size(Pr%mem(b)%v)
            allocate(Tq(nsel), Ty(nsel))
            nl = 0
            do qq = 1, nsel
                q = qq
                if (Cut) q = Ord(b)%v(qq)
                pq = Pr%mem(b)%v(q)
                k3 = Pr%k1(pq)
                k4 = Pr%k2(pq)
                if (.not. (SDX%Ac(k3) .and. SDX%Ac(k4))) cycle
                nl = nl + 1
                Tq(nl) = q
                Ty(nl) = merge(1,0,SDX%Pt(k3)) + merge(1,0,SDX%Pt(k4))
            end do
            do ii = 0, 2
                Lt(ii,b)%v = pack(Tq(1:nl),Ty(1:nl) == ii)
            end do
            deallocate(Tq, Ty)
        end do
        !$omp parallel private(Pos,Rl,Acc,o,Occ,j,i1_,i2_,k1,k2,k3,k4,pid,b,q,qq,pq,sp,nl,v,sg,nb1,nb2,ii,p3,p4,hh,h2,qn,h0,rs,kx,ty_,S3) &
        !$omp firstprivate(base)
        allocate(Pos(SDX%nsd), Rl(4096), Acc(4096), o(Nt), Occ(SPS%ns))
        Pos = 0
        Occ = .false.
        !$omp do schedule(dynamic,16)
        do s = 1, SDX%nsd
            nl = 0
            rs = Rank(s)
            if (present(Skip)) then
                if (Skip(s)) then
                    allocate(Mt(s)%r(0))
                    if (Single) then
                        allocate(Mt(s)%vs(0))
                    else
                        allocate(Mt(s)%v(0))
                    end if
                    cycle
                end if
            end if
            S3 = .false.
            if (Rd) then
                S3 = s > lim3
                if (S3 .and. s > abs(Tm(s))) then       ! not a representative
                    allocate(Mt(s)%r(0))
                    if (Single) then
                        allocate(Mt(s)%vs(0))
                    else
                        allocate(Mt(s)%v(0))
                    end if
                    cycle
                end if
            end if
            do ii = 1, Nt
                Occ(SDX%occ(ii,s)) = .true.
            end do
            ! one-body (the diagonal bucket) ...
            do j = SDX%Jo(b0,s), SDX%Jo(b0+1,s) - 1
                sp = abs(SDX%Jt(j))
                v = Op%uT(int(SDX%Jx(j)) + IXOFF,t)*real(sign(1,SDX%Jt(j)),dp)
                if (Rd) then
                    if (.not. red_tgt(sp,v,S3)) cycle
                end if
                if (Rank(sp) < rs) cycle
                if (v /= 0.0_dp) call acc_add(Pos,Rl,Acc,nl,sp,v)
            end do
            ! two-body ...
            if (Nt >= 2) then
                h0 = SDX%h(s)
                do i1_ = 1, Nt - 1
                    k1 = SDX%occ(i1_,s)
                    do i2_ = i1_ + 1, Nt
                        k2 = SDX%occ(i2_,s)
                        pid = Pr%id(k1,k2)
                        if (pid == 0) cycle
                        b = Pr%blk(pid)
                        nb1 = i1_ + i2_ - 3
                        hh = h0 - merge(1,0,SDX%Pt(k1)) - merge(1,0,SDX%Pt(k2))
                        qn = SDX%q(s) - SPS%e(k1) - SPS%e(k2)
                        if (Cut) base = SDX%ex(s) - Eps(k1) - Eps(k2)
                        do ty_ = 0, min(2,Hcap - hh)
                        do qq = 1, size(Lt(ty_,b)%v)
                            q = Lt(ty_,b)%v(qq)
                            pq = Pr%mem(b)%v(q)
                            if (Cut) then
                                if (base + Ep(pq) > exlim + 1.0e-9_dp) exit
                            end if
                            k3 = Pr%k1(pq)
                            k4 = Pr%k2(pq)
                            if (.not. (SDX%Ac(k3) .and. SDX%Ac(k4))) cycle
                            if (k3 /= k1 .and. k3 /= k2) then
                                if (Occ(k3)) cycle
                            end if
                            if (k4 /= k1 .and. k4 /= k2) then
                                if (Occ(k4)) cycle
                            end if
                            if (qn + SPS%e(k3) + SPS%e(k4) > Nmaxq) cycle
                            h2 = hh + merge(1,0,SDX%Pt(k3)) + merge(1,0,SDX%Pt(k4))
                            if (h2 > Hcap) cycle
                            if (Cut) then
                                if (base + Ep(pq) > Exmax(min(h2,ubound(Exmax,1))) + 1.0e-9_dp) cycle
                            end if
                            v = Op%T(t)%B(b)%v(q,Pr%loc(pid))
                            if (v == 0.0_dp) cycle
                            o = SDX%occ(:,s)
                            o(i1_) = k3
                            o(i2_) = k4
                            call isort(o)
                            kx = sd_key(SPS,SDX,o)
                            if (kx < 0) cycle
                            sp = hashmap_lookup(SDX%Map,kx)
                            if (sp < 0) cycle
                            if (.not. Rd) then
                                if (Rank(sp) < rs) cycle
                            end if
                            p3 = 0
                            p4 = 0
                            do ii = 1, Nt
                                if (o(ii) == k3) p3 = ii
                                if (o(ii) == k4) p4 = ii
                            end do
                            nb2 = p3 + p4 - 3
                            sg = real(1 - 2*mod(nb1 + nb2,2),dp)
                            v = sg*v
                            if (Rd) then
                                if (.not. red_tgt(sp,v,S3)) cycle
                                if (Rank(sp) < rs) cycle
                            end if
                            call acc_add(Pos,Rl,Acc,nl,sp,v)
                        end do
                        end do
                    end do
                end do
            end if
            do ii = 1, Nt
                Occ(SDX%occ(ii,s)) = .false.
            end do
            if (nl > 1) call sort2(nl,Rl,Acc)       ! ascending partners (ascending h' in the hole-particle form)
            allocate(Mt(s)%r(nl))
            Mt(s)%r = Rl(1:nl)
            if (Single) then
                allocate(Mt(s)%vs(nl))
                Mt(s)%vs = real(Acc(1:nl),f4)
            else
                allocate(Mt(s)%v(nl))
                Mt(s)%v = Acc(1:nl)
            end if
            Pos(Rl(1:nl)) = 0
        end do
        !$omp end do
        deallocate(Pos, Rl, Acc, o, Occ)
        !$omp end parallel
    contains
        ! reduced block: the partner t of the column (S3: a representative of h = Hk) ... t of h = Hk & not a representative:
            ! t -> its image b with the value eps sg_b v (S3) or dropped (columns of h < Hk) ...
        logical function red_tgt(t_,v_,S3_)
            integer, intent(inout) :: t_
            real(dp), intent(inout) :: v_
            logical, intent(in) :: S3_
            red_tgt = .true.
            if (t_ <= lim3) return
            if (t_ <= abs(Tm(t_))) return
            if (.not. S3_) then
                red_tgt = .false.
                return
            end if
            v_ = v_*real(epr*sign(1,Tm(t_)),dp)
            t_ = abs(Tm(t_))
        end function red_tgt
    end subroutine species_matrix

    ! Jump-type tables of Op: W(b)%v(ixn,ixp) = v_pn(kc kc', ka ka') for the proton type (kc,ka) = ixp of bucket b & the
        ! neutron type (kc',ka') = ixn of the partner bucket Bp(b) (zero outside the pair quanta cut), uT(ix,t) the
        ! one-body u of the types of the diagonal bucket ...
    subroutine pn_tables(SPS,Pr3,Op)
        type(sp_t), intent(in) :: SPS
        type(pairs_t), intent(in) :: Pr3
        type(op_t), intent(inout) :: Op
        integer :: b, bn, ixp, ixn, kc, ka, kcn, kan, p1, p2, bb, t
        real(dp) :: v
        integer(i8) :: tot
        if (allocated(Op%W)) deallocate(Op%W)
        if (allocated(Op%uT)) deallocate(Op%uT)
        if (allocated(Op%Bp)) deallocate(Op%Bp)
        allocate(Op%W(SPS%nbk), Op%Bp(SPS%nbk))
        Op%Bp = SPS%Bp
        tot = 0
        do b = 1, SPS%nbk
            bn = SPS%Bp(b)
            allocate(Op%W(b)%v(SPS%Nty(bn),SPS%Nty(b)))
            tot = tot + int(SPS%Nty(bn),i8)*int(SPS%Nty(b),i8)
            !$omp parallel do schedule(dynamic,8) private(ixn,kc,ka,kcn,kan,p1,p2,bb,v)
            do ixp = 1, SPS%Nty(b)
                kc = SPS%Tkc(b)%v(ixp)
                ka = SPS%Tka(b)%v(ixp)
                do ixn = 1, SPS%Nty(bn)
                    kcn = SPS%Tkc(bn)%v(ixn)
                    kan = SPS%Tka(bn)%v(ixn)
                    v = 0.0_dp
                    p1 = Pr3%id(kc,kcn)
                    p2 = Pr3%id(ka,kan)
                    if (p1 > 0 .and. p2 > 0) then
                        bb = Pr3%blk(p1)
                        if (Pr3%blk(p2) == bb) v = Op%T(3)%B(bb)%v(Pr3%loc(p1),Pr3%loc(p2))
                    end if
                    Op%W(b)%v(ixn,ixp) = v
                end do
            end do
            !$omp end parallel do
        end do
        allocate(Op%uT(SPS%Nty(SPS%b0),2))
        do t = 1, 2
            do ixp = 1, SPS%Nty(SPS%b0)
                Op%uT(ixp,t) = Op%u(SPS%Tkc(SPS%b0)%v(ixp),SPS%Tka(SPS%b0)%v(ixp),t)
            end do
        end do
        write(*,'(a)') TAB // '  pn jump-type tables: ' // istr8(tot) // ' entries (' // cfmt_f(8.0e-9_dp*real(tot,dp),0,2) // ' GB)'
    end subroutine pn_tables

    ! ---------------------------------------------------------------------------------------------------------------
    ! Operator engine ... column c of a symmetric operator, lower half (rows r >= c), from the species matrices (lower
        ! halves by the ranks: partners sp' >= sp, Rn(sn') >= Rn(sn)) & the pn jumps (proton jump x neutron jump of the
        ! partner buckets, value W(b)(ixn,ixp), rows sp' > sp or sp' = sp & Rn(sn') >= Rn(sn)) ...
        ! Mode 1 (build): accumulated in Pos, Rl, Acc (nl entries) ... Mode 2 (expectation): S(1:nb) = sum_(r>=c) Op(r,c)
        ! X(1:nb,r) & Dg = Op(c,c) ... Mode 3 (product): S(1) = sum_(r>=c) Op(r,c) X(1,r) & Yb(r) += Op(r,c) X(1,c), r > c ...

    subroutine eng_make(SPS,P,Nn,Pr,Op,B,Eng,Single,Nmaxq,Hcap,EpsP,EpsN,ExmP,ExmN,Sym,Eta)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(pairs_t), intent(in) :: Pr(3)
        type(op_t), intent(in) :: Op
        type(mb_t), intent(in) :: B
        type(eng_t), intent(out) :: Eng
        logical, intent(in) :: Single
        integer, intent(in) :: Nmaxq, Hcap
        real(dp), intent(in), optional :: EpsP(:), EpsN(:), ExmP(0:), ExmN(0:)
        logical, intent(in), optional :: Sym
        integer, intent(in), optional :: Eta
        integer, allocatable :: Ip(:)
        integer :: epP, epN
        logical, allocatable :: SkP(:), SkN(:)
        integer :: i
        integer(i8) :: npp, nnn
        allocate(Ip(P%nsd))
        do i = 1, P%nsd
            Ip(i) = i
        end do
        Eng%Single = Single
        allocate(SkP(P%nsd), SkN(Nn%nsd))
        SkP = .false.
        SkN = .false.
        if (present(Sym)) then
            if (Sym) call sym_maps(SPS,P,Nn,B,Eng)
            if (Eng%Sym .and. present(Eta)) call red_setup(P,Nn,B,Eng,Eta)
            if (Eng%Sym) then
                SkP = P%m2 < 0
                SkN = Nn%m2 > 0
                if (Eng%RedP) write(*,'(a)') TAB // '  reduced proton block h = ' // istr(P%Hk) // ' (representatives, Eta = ' // istr(Eng%Eta) // ')'
                if (Eng%RedN) write(*,'(a)') TAB // '  reduced neutron block h = ' // istr(Nn%Hk) // ' (representatives, Eta = ' // istr(Eng%Eta) // ')'
                write(*,'(a)') TAB // '  rotation symmetry R_y(pi): rows with M_P >= 0 only, species matrices of ' // istr(count(.not. SkP)) // &
                               ' of ' // istr(P%nsd) // ' proton & ' // istr(count(.not. SkN)) // ' of ' // istr(Nn%nsd) // ' neutron SDs'
            else if (Sym) then
                write(*,'(a)') TAB // '  rotation symmetry R_y(pi) not used (SD lists not closed under m -> -m)'
            end if
        end if
        epP = 0
        epN = 0
        if (Eng%RedP) epP = Eng%Eta*sign(1,Eng%TN(Eng%RefN))
        if (Eng%RedN) epN = Eng%Eta*sign(1,Eng%TP(Eng%RefP))
        if (Eng%RedP .or. Eng%RedN) then
            call species_matrix(SPS,P,Pr(1),Op,1,Single,Ip,Nmaxq,Hcap,Eng%Mp,Skip=SkP,Tm=Eng%TP,Eps_=epP)
            call species_matrix(SPS,Nn,Pr(2),Op,2,Single,B%Rn,Nmaxq,Hcap,Eng%Mn,Skip=SkN,Tm=Eng%TN,Eps_=epN)
        else if (present(EpsP)) then
            call species_matrix(SPS,P,Pr(1),Op,1,Single,Ip,Nmaxq,Hcap,Eng%Mp,EpsP,ExmP,Skip=SkP)
            call species_matrix(SPS,Nn,Pr(2),Op,2,Single,B%Rn,Nmaxq,Hcap,Eng%Mn,EpsN,ExmN,Skip=SkN)
        else
            call species_matrix(SPS,P,Pr(1),Op,1,Single,Ip,Nmaxq,Hcap,Eng%Mp,Skip=SkP)
            call species_matrix(SPS,Nn,Pr(2),Op,2,Single,B%Rn,Nmaxq,Hcap,Eng%Mn,Skip=SkN)
        end if
        npp = 0
        do i = 1, P%nsd
            npp = npp + size(Eng%Mp(i)%r)
        end do
        nnn = 0
        do i = 1, Nn%nsd
            nnn = nnn + size(Eng%Mn(i)%r)
        end do
        write(*,'(a)') TAB // '  species matrices (lower halves): ' // istr8(npp) // ' proton & ' // istr8(nnn) // &
                       ' neutron entries, ' // cfmt_f(merge(8.0_dp,12.0_dp,Single)*real(npp + nnn,dp)/1.0e9_dp,0,2) // ' GB'
    end subroutine eng_make

    subroutine col_gen(c,P,Nn,B,Op,Eng,Mode,Pos,Rl,Acc,nl,nb,X,S,Dg,Yb)
        integer, intent(in) :: c, Mode, nb
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(op_t), intent(in) :: Op
        type(eng_t), intent(in) :: Eng
        integer, intent(inout) :: Pos(:), nl
        integer, allocatable, intent(inout) :: Rl(:)
        real(dp), allocatable, intent(inout) :: Acc(:)
        real(dp), intent(in) :: X(nb,*)
        real(dp), intent(inout) :: S(nb), Dg, Yb(*)
        integer :: sp, sn, j, bk, bn, st, sp2, sn2, r, Ln, L2, rn, ixp, limp, limn, lim2, h2, hx
        integer(i8) :: jp, jn
        real(dp) :: v, xc, sgp
        sp = B%sp(c)
        sn = B%sn(c)
        Ln = Nn%L(sn)
        rn = B%Rn(sn)
        ! npnh truncation (hole-particle form: SDs, jumps & species-matrix columns in ascending h): the last proton SD that
            ! pairs with sn & the last neutron SD that pairs with sp ... the loops end there ...
        limp = huge(1)
        limn = huge(1)
        hx = -1
        if (P%HP .and. Nn%HP) then
            hx = B%Nmax
            limp = P%He(min(max(hx - Nn%h(sn),0),P%Hk))
            limn = Nn%He(min(max(hx - P%h(sp),0),Nn%Hk))
        end if
        xc = 0.0_dp
        if (Mode == 3) xc = X(1,c)
        if (Op%C0 /= 0.0_dp) call put(c,Op%C0)
        ! protons (partners sp2 >= sp) ...
        do j = 1, size(Eng%Mp(sp)%r)
            sp2 = Eng%Mp(sp)%r(j)
            if (sp2 > limp) exit
            st = B%St(sp2,Ln)
            if (st < 0) cycle
            if (UseCut) then
                if (B%Lc(sn) > B%Nc(sp2,Ln)) cycle
            end if
            r = st + B%Lc(sn)
            if (Eng%Single) then
                v = real(Eng%Mp(sp)%vs(j),dp)
            else
                v = Eng%Mp(sp)%v(j)
            end if
            call put(r,v)
        end do
        ! neutrons (partners Rn(sn2) >= Rn(sn)) ...
        do j = 1, size(Eng%Mn(sn)%r)
            sn2 = Eng%Mn(sn)%r(j)
            if (sn2 > limn) exit
            L2 = Nn%L(sn2)
            st = B%St(sp,L2)
            if (st < 0) cycle
            if (UseCut) then
                if (B%Lc(sn2) > B%Nc(sp,L2)) cycle
            end if
            r = st + B%Lc(sn2)
            if (Eng%Single) then
                v = real(Eng%Mn(sn)%vs(j),dp)
            else
                v = Eng%Mn(sn)%v(j)
            end if
            call put(r,v)
        end do
        ! proton-neutron ...
        do bk = 1, P%nd
            if (P%Jo(bk+1,sp) == P%Jo(bk,sp)) cycle
            bn = Op%Bp(bk)
            if (Nn%Jo(bn+1,sn) == Nn%Jo(bn,sn)) cycle
            do jp = P%Jo(bk,sp), P%Jo(bk+1,sp) - 1
                sp2 = abs(P%Jt(jp))
                if (sp2 < sp) cycle
                lim2 = huge(1)
                if (hx >= 0) then
                    h2 = 0
                    do while (sp2 > P%He(h2))
                        h2 = h2 + 1
                    end do
                    lim2 = Nn%He(min(max(hx - h2,0),Nn%Hk))
                end if
                sgp = real(sign(1,P%Jt(jp)),dp)
                ixp = int(P%Jx(jp)) + IXOFF
                do jn = Nn%Jo(bn,sn), Nn%Jo(bn+1,sn) - 1
                    sn2 = abs(Nn%Jt(jn))
                    if (sn2 > lim2) exit
                    v = Op%W(bk)%v(int(Nn%Jx(jn)) + IXOFF,ixp)
                    if (v == 0.0_dp) cycle
                    L2 = Nn%L(sn2)
                    st = B%St(sp2,L2)
                    if (st < 0) cycle
                    if (UseCut) then
                        if (B%Lc(sn2) > B%Nc(sp2,L2)) cycle
                    end if
                    if (sp2 == sp) then
                        if (B%Rn(sn2) < rn) cycle
                    end if
                    r = st + B%Lc(sn2)
                    call put(r,v*sgp*real(sign(1,Nn%Jt(jn)),dp))
                end do
            end do
        end do
    contains
        subroutine put(r,v)
            integer, intent(in) :: r
            real(dp), intent(in) :: v
            if (Mode == 1) then
                call acc_add(Pos,Rl,Acc,nl,r,v)
            else if (Mode == 2) then
                S = S + v*X(:,r)
                if (r == c) Dg = Dg + v
            else
                S(1) = S(1) + v*X(1,r)
                if (r /= c) Yb(r) = Yb(r) + v*xc
            end if
        end subroutine put
    end subroutine col_gen

    ! Stored lower half ... one pass: the columns of each thread (chunks of CHK, cyclic) into its segments ... Fits =
        ! .false. (everything freed) when the segments would exceed MemGB (the product buffers Yb & y reserved) ...
    subroutine eng_store(P,Nn,B,Op,Eng,Single,MemGB,Fits)
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(op_t), intent(in) :: Op
        type(eng_t), intent(inout) :: Eng
        logical, intent(in) :: Single
        real(dp), intent(in) :: MemGB
        logical, intent(out) :: Fits
        integer, allocatable :: Pos(:), Rl(:)
        real(dp), allocatable :: Acc(:)
        real(dp) :: Xd(1,1), Sd(1), Yd(1), Dg, segGB, Gb
        integer :: nt, tid, nchunk, k, c, c0, c1, nl, m, j, cur, ps, maxseg, myseg, nseg, slots, ntr
        integer(i8) :: nnz, h1, h2, h3
        integer :: prev, d, pe, ne
        logical :: Abort, ab
        nt = omp_get_max_threads()
        Eng%nt = nt
        Eng%Single = Single
        nchunk = (B%dim + CHK - 1)/CHK
        allocate(Eng%Cnt(B%dim), Eng%Ts(nt), Eng%ChT(nchunk), Eng%ChS(nchunk), Eng%ChP(nchunk), Eng%ChE(nchunk))
        Eng%Cnt = 0
        Eng%ChT = 1
        Eng%ChS = 0
        Eng%ChP = SEG
        Eng%ChE = 0
        h1 = 0
        h2 = 0
        h3 = 0
        segGB = (real(SEG,dp)*merge(6.0_dp,10.0_dp,Single) + 4.0_dp*real(ESEG,dp))/1.0e9_dp
        maxseg = int(max(MemGB - 8.0e-9_dp*real(B%dim,dp)*real(nt + 1,dp), 0.0_dp)/segGB)
        slots = min(maxseg, 4*(maxseg/nt) + 16)
        nseg = 0
        Abort = .false.
        !$omp parallel private(tid,k,c,c0,c1,nl,m,j,cur,ps,myseg,Pos,Rl,Acc,Sd,Dg,ab,ntr,prev,d,pe,ne) reduction(+:h1,h2,h3)
        tid = omp_get_thread_num() + 1
        ntr = omp_get_num_threads()
        if (ntr /= nt) then
            !$omp atomic write
            Abort = .true.
        end if
        allocate(Eng%Ts(tid)%S(slots))
        allocate(Pos(B%dim), Rl(65536), Acc(65536))
        Pos = 0
        cur = 0
        ps = SEG
        pe = 0
        !$omp do schedule(dynamic,1)
        do k = 0, nchunk - 1
            !$omp atomic read
            ab = Abort
            if (ab) cycle
            Eng%ChT(k+1) = tid
            Eng%ChS(k+1) = cur
            Eng%ChP(k+1) = ps
            Eng%ChE(k+1) = pe
            c0 = k*CHK + 1
            c1 = min(B%dim, c0 + CHK - 1)
            do c = c0, c1
                nl = 0
                call col_gen(c,P,Nn,B,Op,Eng,1,Pos,Rl,Acc,nl,1,Xd,Sd,Dg,Yd)
                Pos(Rl(1:nl)) = 0
                m = 0
                do j = 1, nl
                    if (Acc(j) == 0.0_dp) cycle
                    m = m + 1
                    Rl(m) = Rl(j)
                    Acc(m) = Acc(j)
                end do
                if (m > SEG) call fatal('NCSM ... a column exceeds the segment size')
                if (m == 0) cycle
                call sort2(m,Rl,Acc)
                prev = c - 1
                ne = 0
                do j = 1, m
                    d = Rl(j) - prev
                    prev = Rl(j)
                    if (d <= 255) then
                        h1 = h1 + 1
                    else if (d <= 65535) then
                        h2 = h2 + 1
                    else
                        h3 = h3 + 1
                        ne = ne + 1
                    end if
                end do
                if (ne > ESEG) call fatal('NCSM ... a column exceeds the escape list of a segment')
                Eng%Cnt(c) = m
                if (SEG - ps < m .or. ESEG - pe < ne) then
                    Eng%Cnt(c) = -m
                    !$omp atomic capture
                    nseg = nseg + 1
                    myseg = nseg
                    !$omp end atomic
                    if (myseg > maxseg .or. cur + 1 > slots) then
                        !$omp atomic write
                        Abort = .true.
                        exit
                    end if
                    cur = cur + 1
                    ps = 0
                    pe = 0
                    allocate(Eng%Ts(tid)%S(cur)%d(SEG), Eng%Ts(tid)%S(cur)%e(ESEG))
                    if (Single) then
                        allocate(Eng%Ts(tid)%S(cur)%vs(SEG))
                    else
                        allocate(Eng%Ts(tid)%S(cur)%vd(SEG))
                    end if
                    Eng%Ts(tid)%ns = cur
                end if
                prev = c - 1
                do j = 1, m
                    d = Rl(j) - prev
                    prev = Rl(j)
                    if (d <= 65535) then
                        Eng%Ts(tid)%S(cur)%d(ps+j) = int(d - GOFF,i2)
                    else
                        Eng%Ts(tid)%S(cur)%d(ps+j) = int(-GOFF,i2)
                        pe = pe + 1
                        Eng%Ts(tid)%S(cur)%e(pe) = Rl(j)
                    end if
                end do
                if (Single) then
                    Eng%Ts(tid)%S(cur)%vs(ps+1:ps+m) = real(Acc(1:m),f4)
                else
                    Eng%Ts(tid)%S(cur)%vd(ps+1:ps+m) = Acc(1:m)
                end if
                ps = ps + m
                Eng%Ts(tid)%S(cur)%n = ps
            end do
        end do
        !$omp end do
        deallocate(Pos, Rl, Acc)
        !$omp end parallel
        if (Abort) then
            deallocate(Eng%Ts, Eng%Cnt)
            Fits = .false.
            write(*,'(a)') TAB // 'stored H: more than ' // istr(maxseg) // ' segments of ' // cfmt_f(segGB,0,3) // &
                           ' GB (limit ' // cfmt_f(MemGB,0,1) // ' GB)'
            return
        end if
        nnz = sum(int(abs(Eng%Cnt),i8))
        Gb = real(nseg,dp)*segGB + 8.0e-9_dp*real(B%dim,dp)*real(nt + 1,dp)
        write(*,'(a)') TAB // 'stored H: ' // istr8(nnz) // ' entries in the lower half (' // &
                       cfmt_f(real(2*nnz - B%dim,dp)/B%dim,0,1) // ' per column), ' // istr(nseg) // ' segments, ' // &
                       cfmt_f(Gb,0,2) // ' GB (' // merge('single','double',Single) // ' precision, incl. ' // istr(nt) // &
                       ' product buffers), limit ' // cfmt_f(MemGB,0,1) // ' GB'
        write(*,'(a)') TAB // '  row gaps (sorted columns): ' // cfmt_f(100.0_dp*real(h1,dp)/real(max(nnz,1_i8),dp),0,2) // ' % <= 255, ' // &
                       cfmt_f(100.0_dp*real(h2,dp)/real(max(nnz,1_i8),dp),0,2) // ' % in 2 bytes, ' // &
                       cfmt_f(100.0_dp*real(h3,dp)/real(max(nnz,1_i8),dp),0,2) // ' % larger'
        Eng%nnz = nnz
        Eng%Stored = .true.
        Fits = .true.
    end subroutine eng_store

    ! y = Op x with the stored lower half ... each thread its own columns (gather) & the transposed part into its buffer
        ! Yb(:,tid) (scatter), then the sum of the buffers ...
    subroutine eng_mv_stored(Eng,n,x,y)
        type(eng_t), intent(inout) :: Eng
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)
        integer :: nt, tid, nchunk, k, c, c0, c1, m, cur, ps, j, r, i0, i1, t, ow, pe, dl
        real(dp) :: s, xc, v
        nt = Eng%nt
        nchunk = (n + CHK - 1)/CHK
        if (.not. allocated(Eng%Yb)) allocate(Eng%Yb(n,nt))
        !$omp parallel private(tid,k,c,c0,c1,m,cur,ps,j,r,s,xc,v,i0,i1,t,ow,pe,dl)
        tid = omp_get_thread_num() + 1
        if (omp_get_num_threads() > nt) call fatal('NCSM ... more threads than at the build of the stored matrix')
        Eng%Yb(:,tid) = 0.0_dp
        !$omp do schedule(dynamic,4)
        do k = 0, nchunk - 1
            ow = Eng%ChT(k+1)
            cur = Eng%ChS(k+1)
            ps = Eng%ChP(k+1)
            pe = Eng%ChE(k+1)
            c0 = k*CHK + 1
            c1 = min(n, c0 + CHK - 1)
            do c = c0, c1
                m = Eng%Cnt(c)
                if (m == 0) cycle
                if (m < 0) then
                    m = -m
                    cur = cur + 1
                    ps = 0
                    pe = 0
                end if
                r = c - 1
                xc = x(c)
                s = 0.0_dp
                associate(Sg => Eng%Ts(ow)%S(cur))
                    if (Eng%Single) then
                        do j = ps + 1, ps + m
                            dl = int(Sg%d(j)) + GOFF
                            if (dl == 0) then
                                pe = pe + 1
                                r = Sg%e(pe)
                            else
                                r = r + dl
                            end if
                            v = real(Sg%vs(j),dp)
                            s = s + v*x(r)
                            if (r /= c) Eng%Yb(r,tid) = Eng%Yb(r,tid) + v*xc
                        end do
                    else
                        do j = ps + 1, ps + m
                            dl = int(Sg%d(j)) + GOFF
                            if (dl == 0) then
                                pe = pe + 1
                                r = Sg%e(pe)
                            else
                                r = r + dl
                            end if
                            v = Sg%vd(j)
                            s = s + v*x(r)
                            if (r /= c) Eng%Yb(r,tid) = Eng%Yb(r,tid) + v*xc
                        end do
                    end if
                end associate
                Eng%Yb(c,tid) = Eng%Yb(c,tid) + s
                ps = ps + m
            end do
        end do
        !$omp end do
        !$omp do schedule(static)
        do i0 = 1, n, 4096
            i1 = min(n, i0 + 4095)
            y(i0:i1) = Eng%Yb(i0:i1,1)
            do t = 2, nt
                y(i0:i1) = y(i0:i1) + Eng%Yb(i0:i1,t)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine eng_mv_stored

    ! y = Op x on the fly (symmetric: the lower half of every column, the transposed part into the thread buffers) ...
    subroutine eng_mv_otf(P,Nn,B,Op,Eng,x,y)
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(op_t), intent(in) :: Op
        type(eng_t), intent(inout) :: Eng
        real(dp), intent(in) :: x(B%dim)
        real(dp), intent(out) :: y(B%dim)
        integer :: nt, tid, c, nl, i0, i1, t, n
        integer :: Pos1(1)
        integer, allocatable :: Rl1(:)
        real(dp), allocatable :: Acc1(:)
        real(dp) :: S1(1), Dg
        n = B%dim
        nt = omp_get_max_threads()
        if (allocated(Eng%Yb)) then
            if (size(Eng%Yb,2) /= nt) deallocate(Eng%Yb)
        end if
        if (.not. allocated(Eng%Yb)) allocate(Eng%Yb(n,nt))
        !$omp parallel private(tid,c,nl,S1,Dg,Pos1,Rl1,Acc1,i0,i1,t)
        tid = omp_get_thread_num() + 1
        allocate(Rl1(1), Acc1(1))
        Eng%Yb(:,tid) = 0.0_dp
        !$omp do schedule(dynamic,256)
        do c = 1, n
            S1 = 0.0_dp
            Dg = 0.0_dp
            call col_gen(c,P,Nn,B,Op,Eng,3,Pos1,Rl1,Acc1,nl,1,x,S1,Dg,Eng%Yb(:,tid))
            Eng%Yb(c,tid) = Eng%Yb(c,tid) + S1(1)
        end do
        !$omp end do
        !$omp do schedule(static)
        do i0 = 1, n, 4096
            i1 = min(n, i0 + 4095)
            y(i0:i1) = Eng%Yb(i0:i1,1)
            do t = 2, nt
                y(i0:i1) = y(i0:i1) + Eng%Yb(i0:i1,t)
            end do
        end do
        !$omp end do
        deallocate(Rl1, Acc1)
        !$omp end parallel
    end subroutine eng_mv_otf

    ! ---------------------------------------------------------------------------------------------------------------
    ! Sector product y = Op x without per-thread buffers (eng_mv_otf needs n x threads doubles) ... every row written by
        ! one task only:
        ! pass A over the proton SDs sp (rows of sp contiguous): C0, the neutron species part (partners in the same sp:
        !   symmetric, both halves from the lower half) & the pn part gathered (all proton x neutron jumps, both directions:
        !   y(c) = sum_r Op(r,c) x(r)) ...
        ! pass B over tiles of the neutron blocks (group (2M_P, Ln), positions l0..l1): the proton species part, Mp x 1 on
        !   the tile, rows St(sp,Ln) + l contiguous in l ... a heavy tile (small neutron blocks, many proton SDs) split over
        !   several tasks with atomic updates ...
    ! ---------------------------------------------------------------------------------------------------------------
    ! Rotation R_y(pi) of the SDs: |a m> -> (-1)^(j-m) |a -m> (m-state k -> k - m2(k)), R |s> = sg |s'> with the sign of
        ! the phases & of the reordering ... TP, TN = sg s' ... Eng%Sym only if 2M = 0, no energy cut & all images are in
        ! the SD lists with the same label & the opposite M ...
    subroutine sym_maps(SPS,P,Nn,B,Eng)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(eng_t), intent(inout) :: Eng
        logical :: Ok
        Eng%Sym = .false.
        if (B%M2 /= 0 .or. UseCut) return
        call one(P,Eng%TP,Ok)
        if (Ok) call one(Nn,Eng%TN,Ok)
        if (.not. Ok) then
            if (allocated(Eng%TP)) deallocate(Eng%TP)
            if (allocated(Eng%TN)) deallocate(Eng%TN)
            return
        end if
        Eng%Sym = .true.
    contains
        subroutine one(SDX,T,Ok_)
            type(sd_t), intent(in) :: SDX
            integer, allocatable, intent(out) :: T(:)
            logical, intent(out) :: Ok_
            integer, allocatable :: o(:)
            integer :: s, ii, jj, k, ph, s2, tmp_, nbad
            integer(i8) :: kx
            allocate(T(SDX%nsd))
            nbad = 0
            !$omp parallel private(o,ii,jj,k,ph,s2,tmp_,kx) reduction(+:nbad)
            allocate(o(SDX%Nt))
            !$omp do schedule(static)
            do s = 1, SDX%nsd
                ph = 0
                do ii = 1, SDX%Nt
                    k = SDX%occ(ii,s)
                    ph = ph + (SPS%j2(k) - SPS%m2(k))/2
                    o(ii) = k - SPS%m2(k)
                end do
                do ii = 2, SDX%Nt           ! insertion sort, the swaps counted
                    jj = ii
                    do while (jj > 1)
                        if (o(jj-1) <= o(jj)) exit
                        tmp_ = o(jj)
                        o(jj) = o(jj-1)
                        o(jj-1) = tmp_
                        ph = ph + 1
                        jj = jj - 1
                    end do
                end do
                s2 = -1
                kx = sd_key(SPS,SDX,o)
                if (kx >= 0) s2 = hashmap_lookup(SDX%Map,kx)
                if (s2 >= 1) then
                    if (SDX%m2(s2) /= -SDX%m2(s) .or. SDX%L(s2) /= SDX%L(s)) s2 = -1
                end if
                if (s2 < 1) then
                    nbad = nbad + 1
                    T(s) = 0
                else
                    T(s) = s2*(1 - 2*mod(ph,2))
                end if
            end do
            !$omp end do
            deallocate(o)
            !$omp end parallel
            Ok_ = nbad == 0
            if (Ok_) then           ! an involution
                !$omp parallel do reduction(+:nbad)
                do s = 1, SDX%nsd
                    if (abs(T(abs(T(s)))) /= s) nbad = nbad + 1
                    if (sign(1,T(abs(T(s)))) /= sign(1,T(s))) nbad = nbad + 1
                end do
                !$omp end parallel do
                Ok_ = nbad == 0
            end if
        end subroutine one
    end subroutine sym_maps

    ! Reduced blocks of h = Hk (see eng_t) ... per species if Hk = Nmax (npnh, hole-particle form), the other species has
        ! one SD of h = 0 (the first, self-conjugate) & no SD of h = Hk is self-conjugate ...
    subroutine red_setup(P,Nn,B,Eng,Eta)
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(eng_t), intent(inout) :: Eng
        integer, intent(in) :: Eta
        character(len=8) :: Ev_
        if (.not. (P%HP .and. Nn%HP) .or. Trunc /= 'npnh') return
        if (P%Hk < 1 .or. Nn%Hk < 1) return
        if (P%He(0) /= 1 .or. Nn%He(0) /= 1) return
        if (abs(Eng%TP(1)) /= 1 .or. abs(Eng%TN(1)) /= 1) return
        Eng%RefP = 1
        Eng%RefN = 1
        Eng%RedP = P%Hk == B%Nmax .and. nosc(P,Eng%TP)
        Eng%RedN = Nn%Hk == B%Nmax .and. nosc(Nn,Eng%TN)
        call get_environment_variable('NJ_RED',Ev_)
        if (trim(Ev_) == 'p' .or. trim(Ev_) == 'none') Eng%RedN = .false.
        if (trim(Ev_) == 'n' .or. trim(Ev_) == 'none') Eng%RedP = .false.
        if (Eng%RedP .or. Eng%RedN) Eng%Eta = Eta
    contains
        logical function nosc(SDX,T)
            type(sd_t), intent(in) :: SDX
            integer, intent(in) :: T(:)
            integer :: s_
            nosc = .true.
            do s_ = SDX%He(SDX%Hk - 1) + 1, SDX%nsd
                if (abs(T(s_)) == s_) then
                    nosc = .false.
                    return
                end if
            end do
        end function nosc
    end subroutine red_setup

    ! Image of the state c = (sp,sn): R |c> = sg |c'>, c' in c2 & sg in the result ...
    integer function sym_img(B,Eng,sp,c,c2) result(sg)
        type(mb_t), intent(in) :: B
        type(eng_t), intent(in) :: Eng
        integer, intent(in) :: sp, c
        integer, intent(out) :: c2
        integer :: tn, pk
        tn = Eng%TN(B%sn(c))
        pk = Eng%Pk(abs(tn))
        c2 = B%St(abs(Eng%TP(sp)),iand(pk,31)) + ishft(pk,-5)
        sg = sign(1,Eng%TP(sp))*sign(1,tn)
    end function sym_img

    ! y of the rows not computed from their images: y(c) = Eta sg y(c') ... the rows of a self-conjugate proton SD
        ! (sp' = sp, M_P = 0) are all computed: the pairs made exactly symmetric (c' > c from c, zero if c' = c & Eta sg = -1) ...
    subroutine sym_fill(P,B,Eng,y)
        type(sd_t), intent(in) :: P
        type(mb_t), intent(in) :: B
        type(eng_t), intent(in) :: Eng
        real(dp), intent(inout) :: y(B%dim)
        integer :: sp, tp, Ln, c, c2, sg
        !$omp parallel do schedule(dynamic,64) private(tp,Ln,c,c2,sg)
        do sp = 1, P%nsd
            tp = abs(Eng%TP(sp))
            if (P%m2(sp) > 0 .or. (P%m2(sp) == 0 .and. tp > sp)) cycle
            do Ln = 0, ubound(B%St,2)
                if (B%St(sp,Ln) < 0) cycle
                if (tp /= sp) then
                    do c = B%St(sp,Ln) + 1, B%St(sp,Ln) + B%Nc(sp,Ln)
                        sg = sym_img(B,Eng,sp,c,c2)
                        y(c) = real(Eng%Eta*sg,dp)*y(c2)
                    end do
                else
                    do c = B%St(sp,Ln) + 1, B%St(sp,Ln) + B%Nc(sp,Ln)
                        sg = sym_img(B,Eng,sp,c,c2)
                        if (c2 > c) then
                            y(c2) = real(Eng%Eta*sg,dp)*y(c)
                        else if (c2 == c .and. Eng%Eta*sg < 0) then
                            y(c) = 0.0_dp
                        end if
                    end do
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine sym_fill

    ! <y|R|y>/<y|y> = sum_c sg(c) y(c) y(c') / sum_c y(c)^2 (Eta = +-1 for a vector of the symmetry) ...
    real(dp) function sym_eta(P,B,Eng,y) result(e)
        type(sd_t), intent(in) :: P
        type(mb_t), intent(in) :: B
        type(eng_t), intent(in) :: Eng
        real(dp), intent(in) :: y(B%dim)
        integer :: sp, Ln, c, c2, sg
        real(dp) :: a, d
        a = 0.0_dp
        d = 0.0_dp
        !$omp parallel do schedule(dynamic,64) private(Ln,c,c2,sg) reduction(+:a,d)
        do sp = 1, P%nsd
            do Ln = 0, ubound(B%St,2)
                if (B%St(sp,Ln) < 0) cycle
                do c = B%St(sp,Ln) + 1, B%St(sp,Ln) + B%Nc(sp,Ln)
                    sg = sym_img(B,Eng,sp,c,c2)
                    a = a + real(sg,dp)*y(c)*y(c2)
                    d = d + y(c)*y(c)
                end do
            end do
        end do
        !$omp end parallel do
        e = 0.0_dp
        if (d > 0.0_dp) e = a/d
    end function sym_eta

    subroutine eng_sec_make(P,Nn,B,Eng)
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(eng_t), intent(inout) :: Eng
        integer, parameter :: TILE = 64, ABIG = 4096, ACHK = 1024
        integer :: Lmax, Mm, sp, Ln, g, ng, k, nt, ntl, it, nch, ic, l0, l1, lm
        integer, allocatable :: Gi(:,:), Cn(:), Fill(:), Gl(:), Ord(:)
        integer(i8), allocatable :: Ew(:)
        integer(i8) :: wtot, wt, w, e0, ec, k0, kk, nsp
        real(dp), allocatable :: Wk(:)
        integer :: ntk
        integer, allocatable :: Tg2(:), Tl02(:), Tl12(:)
        integer(i8), allocatable :: Tk02(:), Tk12(:)
        logical, allocatable :: Ta2(:)
        Lmax = ubound(B%St,2)
        Mm = maxval(abs(P%m2))
        nt = omp_get_max_threads()
        ! groups (2M_P, Ln) with their sectors ...
        allocate(Gi(-Mm:Mm,0:Lmax))
        Gi = 0
        ng = 0
        do sp = 1, P%nsd
            if (Eng%Sym .and. P%m2(sp) < 0) cycle
            do Ln = 0, Lmax
                if (B%St(sp,Ln) < 0) cycle
                if (Gi(P%m2(sp),Ln) == 0) then
                    ng = ng + 1
                    Gi(P%m2(sp),Ln) = ng
                end if
            end do
        end do
        allocate(Cn(ng), Eng%GLn(ng), Eng%Go(ng+1), Ew(ng), Gl(ng))
        Cn = 0
        Ew = 0
        Gl = 0
        do sp = 1, P%nsd
            if (Eng%Sym .and. P%m2(sp) < 0) cycle
            do Ln = 0, Lmax
                if (B%St(sp,Ln) < 0) cycle
                g = Gi(P%m2(sp),Ln)
                Eng%GLn(g) = Ln
                Cn(g) = Cn(g) + 1
                Ew(g) = Ew(g) + size(Eng%Mp(sp)%r)
                Gl(g) = max(Gl(g),B%Nc(sp,Ln))
            end do
        end do
        Eng%Go(1) = 0
        do g = 1, ng
            Eng%Go(g+1) = Eng%Go(g) + Cn(g)
        end do
        allocate(Eng%Gs(Eng%Go(ng+1)), Fill(ng))
        Fill = 0
        do sp = 1, P%nsd
            if (Eng%Sym .and. P%m2(sp) < 0) cycle
            do Ln = 0, Lmax
                if (B%St(sp,Ln) < 0) cycle
                g = Gi(P%m2(sp),Ln)
                Fill(g) = Fill(g) + 1
                Eng%Gs(Eng%Go(g) + Fill(g)) = sp
            end do
        end do
        ! tasks: tiles of TILE positions, a tile above the target work wt split by entries over ceil(w/wt) tasks ...
        wtot = 0
        do g = 1, ng
            wtot = wtot + Ew(g)*int(Gl(g),i8)
        end do
        wt = max(wtot/(int(nt,i8)*32_i8),1_i8)
        ntk = 0
        do it = 1, 2          ! count, then fill
            k = 0
            do g = 1, ng
                lm = Gl(g)
                do l0 = 1, lm, TILE
                    l1 = min(lm, l0 + TILE - 1)
                    w = Ew(g)*int(l1 - l0 + 1,i8)
                    nch = 1
                    if (w > 4*wt) nch = int(min((w + wt - 1)/wt, Cn(g) + 0_i8, 4_i8*nt))
                    if (it == 1) then
                        k = k + nch
                        cycle
                    end if
                    ! chunks of the proton SDs of g by entries ...
                    e0 = 0
                    ec = 0
                    k0 = Eng%Go(g) + 1
                    ic = 1
                    do kk = Eng%Go(g) + 1, Eng%Go(g+1)
                        e0 = e0 + size(Eng%Mp(Eng%Gs(kk))%r)
                        ec = ec + size(Eng%Mp(Eng%Gs(kk))%r)
                        if ((ec*nch >= Ew(g)*ic .and. ic < nch) .or. kk == Eng%Go(g+1)) then
                            k = k + 1
                            Tg2(k) = g
                            Tl02(k) = l0
                            Tl12(k) = l1
                            Tk02(k) = k0
                            Tk12(k) = kk
                            Ta2(k) = nch > 1
                            Wk(k) = real(e0,dp)*real(l1 - l0 + 1,dp)
                            e0 = 0
                            k0 = kk + 1
                            ic = ic + 1
                        end if
                    end do
                end do
            end do
            if (it == 1) then
                ntk = k
                allocate(Tg2(ntk), Tl02(ntk), Tl12(ntk), Tk02(ntk), Tk12(ntk), Ta2(ntk), Wk(ntk))
            end if
        end do
        ntk = k
        ! heaviest tasks first (dynamic schedule) ...
        allocate(Ord(ntk))
        call argsort_r(-Wk(1:ntk),Ord)
        Eng%ntask = ntk
        allocate(Eng%Tg(ntk), Eng%Tl0(ntk), Eng%Tl1(ntk), Eng%Tk0(ntk), Eng%Tk1(ntk), Eng%Ta(ntk))
        do k = 1, ntk
            Eng%Tg(k) = Tg2(Ord(k))
            Eng%Tl0(k) = Tl02(Ord(k))
            Eng%Tl1(k) = Tl12(Ord(k))
            Eng%Tk0(k) = Tk02(Ord(k))
            Eng%Tk1(k) = Tk12(Ord(k))
            Eng%Ta(k) = Ta2(Ord(k))
        end do
        allocate(Eng%Pk(Nn%nsd))
        if (maxval(Nn%L) > 31 .or. maxval(B%Lc) >= 2**26) call fatal('NCSM ... sector product: L or Lc out of the packing range')
        do sp = 1, Nn%nsd
            Eng%Pk(sp) = B%Lc(sp)*32 + Nn%L(sp)
        end do
        ! pass A tasks ...
        do it = 1, 2
            k = 0
            do sp = 1, P%nsd
                if (Eng%Sym) then
                    if (P%m2(sp) < 0 .or. (P%m2(sp) == 0 .and. abs(Eng%TP(sp)) < sp)) cycle
                end if
                l0 = -1
                l1 = -2
                do Ln = 0, Lmax
                    if (B%St(sp,Ln) < 0) cycle
                    if (l0 < 0) l0 = B%St(sp,Ln) + 1
                    l1 = B%St(sp,Ln) + B%Nc(sp,Ln)
                end do
                if (l1 < l0) cycle
                if (l1 - l0 + 1 <= ABIG) then
                    k = k + 1
                    if (it == 2) then
                        Eng%As(k) = sp
                        Eng%Ar0(k) = l0
                        Eng%Ar1(k) = l1
                        Eng%Aa(k) = .false.
                    end if
                else
                    do lm = l0, l1, ACHK
                        k = k + 1
                        if (it == 2) then
                            Eng%As(k) = sp
                            Eng%Ar0(k) = lm
                            Eng%Ar1(k) = min(l1, lm + ACHK - 1)
                            Eng%Aa(k) = .true.
                        end if
                    end do
                end if
            end do
            if (it == 1) allocate(Eng%As(k), Eng%Ar0(k), Eng%Ar1(k), Eng%Aa(k))
        end do
        Eng%natask = k
        write(*,'(a)') TAB // 'sector product ... pass A ' // istr(k) // ' tasks, ' // istr(count(Eng%Aa)) // &
                       ' chunks of big sectors (> ' // istr(ABIG) // ' rows)'
        Eng%SecReady = .true.
        nsp = count(Eng%Ta)
        write(*,'(a)') TAB // 'sector product ... ' // istr(ng) // ' groups (2M_P, L_N), ' // istr(ntk) // ' tile tasks (' // &
                       istr8(nsp) // ' atomic), proton-part work ' // cfmt_f(real(wtot,dp)/1.0e9_dp,0,2) // ' G entries'
    end subroutine eng_sec_make

    subroutine eng_mv_sec(P,Nn,B,Op,Eng,x,y)
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(op_t), intent(in) :: Op
        type(eng_t), intent(inout) :: Eng
        real(dp), intent(in) :: x(B%dim)
        real(dp), intent(out) :: y(B%dim)
        integer :: pk, lc, i, HkP, jt, nq, q, nw, nwmax, ixq
        integer(i8) :: jg, jq
        logical :: Grp
        real(dp), allocatable :: Wc(:)
        integer, allocatable :: To(:), Tpk(:), Tw(:)
        integer(i8) :: jlo, jhi, jmid
        integer(i8), allocatable :: Jb(:), Je(:,:)
        integer(i8) :: n1, n2, n3, n4, n5
        integer, parameter :: ABIG_ = 4096, ACHK_ = 1024
        real(dp), allocatable :: Ac(:)
        integer, allocatable :: Stc(:)
        integer :: sp, qn, Lmax, c, c0, c1, sn, j, sn2, L2, st, r, limn, limp, hx, bk, bn, ixp, lim2, h2, sp2, k, g, Ln, &
                   a, b2, bb, l, sc, s2, hn
        integer(i8) :: jp, jn, kk
        real(dp) :: xc, acc, v, sgp, acc3
        integer(i8) :: t0, t1, t2, tr
        logical :: HP, At, Mir
        integer :: lim3n, lim3p, c2, sgm, LnRef
        if (.not. Eng%SecReady) call eng_sec_make(P,Nn,B,Eng)
        Lmax = ubound(B%St,2)
        HP = P%HP .and. Nn%HP
        hx = -1
        if (HP) hx = B%Nmax
        call system_clock(t0,tr)
        ! pass A ... y = 0, then the tasks: neutron part (symmetric, rows of the same sp) per row, then C0 & pn gathered
            ! with the proton jump outside & the rows of the task inside (x of the sector of sp2 & the neutron targets Pk of
            ! the rows reused across the proton jumps), sums in Ac(1:nr) ...
        call pscale(B%dim,0.0_dp,x,y)
        HkP = 0
        if (HP) HkP = P%Hk
        lim3n = huge(1)
        lim3p = huge(1)
        LnRef = -1
        if (Eng%RedN) lim3n = Nn%He(Nn%Hk - 1)
        if (Eng%RedP) then
            lim3p = P%He(P%Hk - 1)
            LnRef = Nn%L(Eng%RefN)
        end if
        nwmax = 1
        do k = 1, size(Op%W)
            nwmax = max(nwmax, size(Op%W(k)%v,1))
        end do
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        n5 = 0
        !$omp parallel private(sp,c,sn,j,sn2,L2,st,r,limn,bk,bn,ixp,lim2,h2,sp2,jp,jn,xc,acc,v,sgp,At,pk,lc,Ac,Stc,c0,c1,i,Jb,Je,jlo,jhi,jmid,jt,To,Tpk,Tw,nq,q,nw,jg,jq,ixq,Grp,Wc,Mir,acc3,c2,sgm) &
        !$omp reduction(+:n1,n2,n3,n4,n5)
        allocate(Ac(max(ABIG_,ACHK_)), Stc(0:Lmax), Jb(max(ABIG_,ACHK_)), Je(max(ABIG_,ACHK_),0:HkP), To(max(ABIG_,ACHK_)), &
                 Tpk(65536), Tw(65536), Wc(nwmax))
        !$omp do schedule(dynamic,8)
        do k = 1, Eng%natask
            sp = Eng%As(k)
            At = Eng%Aa(k)
            c0 = Eng%Ar0(k)
            c1 = Eng%Ar1(k)
            limn = huge(1)
            if (HP) limn = Nn%He(min(max(hx - P%h(sp),0),Nn%Hk))
            Mir = Eng%RedN .and. sp == Eng%RefP
            ! neutrons (partners Rn(sn2) >= Rn(sn), same sp) ...
            do c = c0, c1
                xc = x(c)
                acc = Op%C0*xc
                acc3 = 0.0_dp
                if (DoN) then
                sn = B%sn(c)
                do j = 1, size(Eng%Mn(sn)%r)
                    sn2 = Eng%Mn(sn)%r(j)
                    if (sn2 > limn) exit
                    pk = Eng%Pk(sn2)
                    L2 = iand(pk,31)
                    st = B%St(sp,L2)
                    if (st < 0) cycle
                    lc = ishft(pk,-5)
                    if (UseCut) then
                        if (lc > B%Nc(sp,L2)) cycle
                    end if
                    r = st + lc
                    if (Eng%Single) then
                        v = real(Eng%Mn(sn)%vs(j),dp)
                    else
                        v = Eng%Mn(sn)%v(j)
                    end if
                    if (r == c) then
                        acc = acc + v*xc
                    else
                        if (At) then
                            !$omp atomic
                            y(r) = y(r) + v*xc
                        else
                            y(r) = y(r) + v*xc
                        end if
                        acc = acc + v*x(r)
                        if (Mir) then
                            if (sn2 > lim3n .and. sn <= lim3n) acc3 = acc3 + v*x(r)
                        end if
                    end if
                end do
                ! mirror term of the reduced block: the partners of h = Hk are the representatives only ...
                if (acc3 /= 0.0_dp) then
                    sgm = Eng%Eta*sym_img(B,Eng,sp,c,c2)
                    if (c2 == c) then
                        acc = acc + real(sgm,dp)*acc3
                    else if (At) then
                        !$omp atomic
                        y(c2) = y(c2) + real(sgm,dp)*acc3
                    else
                        y(c2) = y(c2) + real(sgm,dp)*acc3
                    end if
                end if
                end if
                Ac(c - c0 + 1) = acc
            end do
            ! proton-neutron, gathered (all jumps) ... per row & bucket the neutron jumps Jb(i) .. Je(i,h2) within the npnh
                ! limit of a proton target with h2 particles (targets ascending): the inner loop without exits (simd) ...
            if (DoPN) then
            do bk = 1, P%nd
                if (P%Jo(bk+1,sp) == P%Jo(bk,sp)) cycle
                bn = Op%Bp(bk)
                ! the neutron jumps of the rows in the buffers Tpk (Pk of the target) & Tw (+-(type index)), row i at
                    ! To(i) + 0 .. Je(i,h2) - Jb(i) ...
                nq = 0
                do c = c0, c1
                    i = c - c0 + 1
                    sn = B%sn(c)
                    Jb(i) = Nn%Jo(bn,sn)
                    do h2 = 0, HkP
                        lim2 = huge(1)
                        if (HP) lim2 = Nn%He(min(max(hx - h2,0),Nn%Hk))
                        ! last jn with |Jt(jn)| <= lim2 (binary search) ...
                        jlo = Nn%Jo(bn,sn) - 1
                        jhi = Nn%Jo(bn+1,sn) - 1
                        do while (jlo < jhi)
                            jmid = jlo + (jhi - jlo + 1)/2
                            if (abs(Nn%Jt(jmid)) <= lim2) then
                                jlo = jmid
                            else
                                jhi = jmid - 1
                            end if
                        end do
                        Je(i,h2) = jlo
                    end do
                    if (Mir) then
                        if (sn > lim3n .and. sn > abs(Eng%TN(sn))) Je(i,0:HkP) = Jb(i) - 1
                    end if
                    To(i) = nq + 1
                    nq = nq + int(max(maxval(Je(i,0:HkP)) - Jb(i) + 1,0_i8))
                end do
                if (nq > size(Tpk)) then
                    deallocate(Tpk, Tw)
                    allocate(Tpk(2*nq), Tw(2*nq))
                end if
                do c = c0, c1
                    i = c - c0 + 1
                    q = To(i)
                    do jn = Jb(i), maxval(Je(i,0:HkP))
                        jt = Nn%Jt(jn)
                        Tpk(q) = Eng%Pk(abs(jt))
                        Tw(q) = sign(int(Nn%Jx(jn)) + IXOFF, jt)
                        q = q + 1
                    end do
                end do
                ! proton jumps to the same target (the diagonal ones a+_k a_k, k occupied) combined: Wc = sum sg W(:,ixp) ...
                nw = size(Op%W(bk)%v,1)
                jp = P%Jo(bk,sp)
                do while (jp <= P%Jo(bk+1,sp) - 1)
                    sp2 = abs(P%Jt(jp))
                    jg = jp
                    do while (jg + 1 <= P%Jo(bk+1,sp) - 1)
                        if (abs(P%Jt(jg+1)) /= sp2) exit
                        jg = jg + 1
                    end do
                    if (Cnt) n5 = n5 + 1
                    h2 = 0
                    if (HP) then
                        do while (sp2 > P%He(h2))
                            h2 = h2 + 1
                        end do
                    end if
                    do i = 0, Lmax
                        Stc(i) = B%St(sp2,i)
                    end do
                    ixp = int(P%Jx(jp)) + IXOFF
                    Grp = jg > jp
                    if (Grp) then
                        Wc(1:nw) = 0.0_dp
                        do jq = jp, jg
                            ixq = int(P%Jx(jq)) + IXOFF
                            if (P%Jt(jq) > 0) then
                                Wc(1:nw) = Wc(1:nw) + Op%W(bk)%v(1:nw,ixq)
                            else
                                Wc(1:nw) = Wc(1:nw) - Op%W(bk)%v(1:nw,ixq)
                            end if
                        end do
                        sgp = 1.0_dp
                    else
                        Wc(1:nw) = 0.0_dp
                        sgp = real(sign(1,P%Jt(jp)),dp)
                    end if
                    do c = c0, c1
                        i = c - c0 + 1
                        acc = 0.0_dp
                        if (Cnt) n2 = n2 + max(Je(i,h2) - Jb(i) + 1,0)
                        if (UseCut) then
                            do jn = Jb(i), Je(i,h2)
                                if (Grp) then
                                    v = Wc(int(Nn%Jx(jn)) + IXOFF)
                                else
                                    v = Op%W(bk)%v(int(Nn%Jx(jn)) + IXOFF,ixp)
                                end if
                                if (v == 0.0_dp) cycle
                                pk = Eng%Pk(abs(Nn%Jt(jn)))
                                L2 = iand(pk,31)
                                st = Stc(L2)
                                if (st < 0) cycle
                                lc = ishft(pk,-5)
                                if (lc > B%Nc(sp2,L2)) cycle
                                r = st + lc
                                if (Nn%Jt(jn) > 0) then
                                    acc = acc + v*x(r)
                                else
                                    acc = acc - v*x(r)
                                end if
                            end do
                        else if (Grp) then
                            !$omp simd reduction(+:acc) private(pk,st,v,jt)
                            do q = To(i), To(i) + int(Je(i,h2) - Jb(i))
                                jt = Tw(q)
                                pk = Tpk(q)
                                st = Stc(iand(pk,31))
                                v = Wc(abs(jt))
                                if (jt < 0) v = -v
                                if (st >= 0) acc = acc + v*x(st + ishft(pk,-5))
                            end do
                        else
                            !$omp simd reduction(+:acc) private(pk,st,v,jt)
                            do q = To(i), To(i) + int(Je(i,h2) - Jb(i))
                                jt = Tw(q)
                                pk = Tpk(q)
                                st = Stc(iand(pk,31))
                                v = Op%W(bk)%v(abs(jt),ixp)
                                if (jt < 0) v = -v
                                if (st >= 0) acc = acc + v*x(st + ishft(pk,-5))
                            end do
                        end if
                        Ac(i) = Ac(i) + sgp*acc
                    end do
                    jp = jg + 1
                end do
            end do
            end if
            if (At) then
                do c = c0, c1
                    !$omp atomic
                    y(c) = y(c) + Ac(c - c0 + 1)
                end do
            else
                do c = c0, c1
                    y(c) = y(c) + Ac(c - c0 + 1)
                end do
            end if
        end do
        !$omp end do
        deallocate(Ac, Stc, Jb, Je, To, Tpk, Tw, Wc)
        !$omp end parallel
        if (Cnt) write(*,'(a)') TAB // 'NJ_MV prof ... pn: proton jumps ' // istr8(n5) // ', neutron jumps within the limits ' // istr8(n2)
        call system_clock(t1)
        ! pass B ...
        !$omp parallel do schedule(dynamic,1) private(g,Ln,limp,hn,kk,sp,sc,a,bb,j,sp2,s2,b2,v,l,At,Mir,c2,sgm)
        do k = 1, Eng%ntask
            g = Eng%Tg(k)
            Ln = Eng%GLn(g)
            At = Eng%Ta(k)
            limp = huge(1)
            if (HP) then
                hn = Ln/2
                limp = P%He(min(max(hx - hn,0),P%Hk))
            end if
            do kk = Eng%Tk0(k), Eng%Tk1(k)
                sp = Eng%Gs(kk)
                sc = B%St(sp,Ln)
                a = Eng%Tl0(k)
                bb = min(Eng%Tl1(k),B%Nc(sp,Ln))
                if (bb < a) cycle
                Mir = Ln == LnRef .and. sp <= lim3p
                if (Mir) then
                    c2 = B%St(abs(Eng%TP(sp)),Ln)
                    sgm = Eng%Eta*sign(1,Eng%TP(sp))*sign(1,Eng%TN(Eng%RefN))
                end if
                do j = 1, size(Eng%Mp(sp)%r)
                    sp2 = Eng%Mp(sp)%r(j)
                    if (sp2 > limp) exit
                    s2 = B%St(sp2,Ln)
                    if (s2 < 0) cycle
                    b2 = min(bb,B%Nc(sp2,Ln))
                    if (b2 < a) cycle
                    if (Eng%Single) then
                        v = real(Eng%Mp(sp)%vs(j),dp)
                    else
                        v = Eng%Mp(sp)%v(j)
                    end if
                    if (Mir) then
                        if (sp2 > lim3p) then
                            do l = a, b2
                                !$omp atomic
                                y(c2+l) = y(c2+l) + real(sgm,dp)*v*x(s2+l)
                            end do
                        end if
                    end if
                    if (At) then
                        if (sp2 == sp) then
                            do l = a, b2
                                !$omp atomic
                                y(sc+l) = y(sc+l) + v*x(sc+l)
                            end do
                        else
                            do l = a, b2
                                !$omp atomic
                                y(s2+l) = y(s2+l) + v*x(sc+l)
                                !$omp atomic
                                y(sc+l) = y(sc+l) + v*x(s2+l)
                            end do
                        end if
                    else
                        if (sp2 == sp) then
                            y(sc+a:sc+b2) = y(sc+a:sc+b2) + v*x(sc+a:sc+b2)
                        else
                            y(s2+a:s2+b2) = y(s2+a:s2+b2) + v*x(sc+a:sc+b2)
                            y(sc+a:sc+b2) = y(sc+a:sc+b2) + v*x(s2+a:s2+b2)
                        end if
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        if (Eng%Sym) then
            if (Eng%Eta == 0) call fatal('NCSM ... symmetric product without the symmetry of the vector (Eta)')
            call sym_fill(P,B,Eng,y)
        end if
        call system_clock(t2)
        TmA = TmA + real(t1 - t0,dp)/real(tr,dp)
        TmB = TmB + real(t2 - t1,dp)/real(tr,dp)
    end subroutine eng_mv_sec

    ! Out(k) = <x_k|Op|x_k> for the nb vectors X(k,:) ... from the lower halves: sum_c x(c) (2 S(c) - Op(c,c) x(c)) ...
    subroutine eng_expect(P,Nn,B,Op,Eng,nb,X,Out)
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B
        type(op_t), intent(in) :: Op
        type(eng_t), intent(in) :: Eng
        integer, intent(in) :: nb
        real(dp), intent(in) :: X(nb,B%dim)
        real(dp), intent(out) :: Out(nb)
        integer :: c, nl
        integer :: Pos1(1)
        integer, allocatable :: Rl1(:)
        real(dp), allocatable :: Acc1(:), Ot(:), S(:)
        real(dp) :: Dg, Yd(1)
        Out = 0.0_dp
        !$omp parallel private(c,nl,S,Dg,Pos1,Rl1,Acc1,Ot)
        allocate(Rl1(1), Acc1(1), Ot(nb), S(nb))
        Ot = 0.0_dp
        !$omp do schedule(dynamic,256)
        do c = 1, B%dim
            S = 0.0_dp
            Dg = 0.0_dp
            call col_gen(c,P,Nn,B,Op,Eng,2,Pos1,Rl1,Acc1,nl,nb,X,S,Dg,Yd)
            Ot = Ot + X(:,c)*(2.0_dp*S - Dg*X(:,c))
        end do
        !$omp end do
        !$omp critical
        Out = Out + Ot
        !$omp end critical
        deallocate(Rl1, Acc1, Ot, S)
        !$omp end parallel
    end subroutine eng_expect

    ! Out(k) = <x_k|J^2|x_k> = |J_+ x_k|^2 + M (M + 1) ... J_+ x in the basis B2 (2M + 2) gathered from the reverse jumps
        ! (bucket dm = -1, parity type 0) of its states: <c'|a+_ka a_kc|c> = sign of the jump a+_kc a_ka of c' to c ...
    subroutine j2_expect(SPS,P,Nn,B,B2,nb,X,Out)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(mb_t), intent(in) :: B, B2
        integer, intent(in) :: nb
        real(dp), intent(in) :: X(nb,B%dim)
        real(dp), intent(out) :: Out(nb)
        real(dp), allocatable :: jT(:), y(:), Ot(:)
        integer :: br, ix, kc, ka, cc, spp, snp, s, st, c, L2
        integer(i8) :: j
        real(dp) :: v, Mr
        br = SPS%dmax - 1
        allocate(jT(SPS%Nty(br)))
        do ix = 1, SPS%Nty(br)
            kc = SPS%Tkc(br)%v(ix)
            ka = SPS%Tka(br)%v(ix)
            jT(ix) = 0.0_dp
            if (SPS%a(kc) == SPS%a(ka)) jT(ix) = 0.5_dp*sqrt(real(SPS%j2(kc)*(SPS%j2(kc) + 2) - SPS%m2(kc)*(SPS%m2(kc) + 2),dp))
        end do
        Out = 0.0_dp
        !$omp parallel private(cc,spp,snp,s,st,c,L2,j,v,y,Ot)
        allocate(y(nb), Ot(nb))
        Ot = 0.0_dp
        !$omp do schedule(dynamic,256)
        do cc = 1, B2%dim
            spp = B2%sp(cc)
            snp = B2%sn(cc)
            y = 0.0_dp
            L2 = Nn%L(snp)
            do j = P%Jo(br,spp), P%Jo(br+1,spp) - 1
                v = jT(int(P%Jx(j)) + IXOFF)
                if (v == 0.0_dp) cycle
                s = abs(P%Jt(j))
                st = B%St(s,L2)
                if (st < 0) cycle
                if (UseCut) then
                    if (B%Lc(snp) > B%Nc(s,L2)) cycle
                end if
                c = st + B%Lc(snp)
                y = y + v*real(sign(1,P%Jt(j)),dp)*X(:,c)
            end do
            do j = Nn%Jo(br,snp), Nn%Jo(br+1,snp) - 1
                v = jT(int(Nn%Jx(j)) + IXOFF)
                if (v == 0.0_dp) cycle
                s = abs(Nn%Jt(j))
                L2 = Nn%L(s)
                st = B%St(spp,L2)
                if (st < 0) cycle
                if (UseCut) then
                    if (B%Lc(s) > B%Nc(spp,L2)) cycle
                end if
                c = st + B%Lc(s)
                y = y + v*real(sign(1,Nn%Jt(j)),dp)*X(:,c)
            end do
            Ot = Ot + y*y
        end do
        !$omp end do
        !$omp critical
        Out = Out + Ot
        !$omp end critical
        deallocate(y, Ot)
        !$omp end parallel
        Mr = 0.5_dp*real(B%M2,dp)
        Out = Out + Mr*(Mr + 1.0_dp)
    end subroutine j2_expect

    ! Sorts k(1:n) ascending with the payload v (quicksort, insertion sort below 16) ...
    subroutine sort2(n,k,v)
        integer, intent(in) :: n
        integer, intent(inout) :: k(:)
        real(dp), intent(inout) :: v(:)
        integer :: Stk(2,64), top, lo, hi, i, j, pv, tk
        real(dp) :: tv
        top = 1
        Stk(1,1) = 1
        Stk(2,1) = n
        do while (top > 0)
            lo = Stk(1,top)
            hi = Stk(2,top)
            top = top - 1
            do while (hi - lo >= 16)
                pv = k((lo + hi)/2)
                i = lo
                j = hi
                do
                    do while (k(i) < pv)
                        i = i + 1
                    end do
                    do while (k(j) > pv)
                        j = j - 1
                    end do
                    if (i <= j) then
                        tk = k(i); k(i) = k(j); k(j) = tk
                        tv = v(i); v(i) = v(j); v(j) = tv
                        i = i + 1
                        j = j - 1
                    end if
                    if (i > j) exit
                end do
                ! the larger part on the stack, the smaller next ...
                if (j - lo < hi - i) then
                    top = top + 1
                    Stk(1,top) = i
                    Stk(2,top) = hi
                    hi = j
                else
                    top = top + 1
                    Stk(1,top) = lo
                    Stk(2,top) = j
                    lo = i
                end if
            end do
            do i = lo + 1, hi
                tk = k(i)
                tv = v(i)
                j = i - 1
                do while (j >= lo)
                    if (k(j) <= tk) exit
                    k(j+1) = k(j)
                    v(j+1) = v(j)
                    j = j - 1
                end do
                k(j+1) = tk
                v(j+1) = tv
            end do
        end do
    end subroutine sort2

    ! Accumulates x at row r of a column under construction (thread-local Pos, Rl, Acc, nl; Rl & Acc grown as needed) ...
    subroutine acc_add(Pos,Rl,Acc,nl,r,x)
        integer, intent(inout) :: Pos(:), nl
        integer, allocatable, intent(inout) :: Rl(:)
        real(dp), allocatable, intent(inout) :: Acc(:)
        integer, intent(in) :: r
        real(dp), intent(in) :: x
        integer, allocatable :: Ri(:)
        real(dp), allocatable :: Ai(:)
        if (Pos(r) == 0) then
            if (nl == size(Rl)) then
                allocate(Ri(2*nl), Ai(2*nl))
                Ri(1:nl) = Rl(1:nl)
                Ai(1:nl) = Acc(1:nl)
                call move_alloc(Ri,Rl)
                call move_alloc(Ai,Acc)
            end if
            nl = nl + 1
            Pos(r) = nl
            Rl(nl) = r
            Acc(nl) = 0.0_dp
        end if
        Acc(Pos(r)) = Acc(Pos(r)) + x
    end subroutine acc_add

    ! Op <- Op + a X (same pair blocks) ...
    subroutine op_axpy(Op,a,X)
        type(op_t), intent(inout) :: Op
        real(dp), intent(in) :: a
        type(op_t), intent(in) :: X
        integer :: t, b
        Op%C0 = Op%C0 + a*X%C0
        Op%u = Op%u + a*X%u
        do t = 1, 3
            do b = 1, size(Op%T(t)%B)
                Op%T(t)%B(b)%v = Op%T(t)%B(b)%v + a*X%T(t)%B(b)%v
            end do
        end do
    end subroutine op_axpy

    ! MemAvailable of /proc/meminfo in GB (0 if not readable) ...
    real(dp) function mem_available_gb()
        integer :: u, Ios
        character(len=256) :: Line
        integer(i8) :: kb
        mem_available_gb = 0.0_dp
        open(newunit=u, file='/proc/meminfo', status='old', action='read', iostat=Ios)
        if (Ios /= 0) return
        do
            read(u,'(a)',iostat=Ios) Line
            if (Ios /= 0) exit
            if (index(Line,'MemAvailable:') == 1) then
                read(Line(14:),*,iostat=Ios) kb
                if (Ios == 0) mem_available_gb = real(kb,dp)*1024.0_dp/1.0e9_dp
                exit
            end if
        end do
        close(u)
    end function mem_available_gb

    ! ---------------------------------------------------------------------------------------------------------------
    ! Thick-restart Lanczos (Wu & Simon) for the lowest Ne eigenpairs ... basis of at most Mb vectors with full
        ! reorthogonalization (twice), the projected matrix T = V^T A V from the Gram-Schmidt coefficients (exact Rayleigh-
        ! Ritz in the current basis, arrowhead after a restart) ... restart: the lowest k Ritz vectors kept (k = Ne + half
        ! of the free part) & the residual direction appended ... convergence |beta y_m,i| < Tol for the Ne lowest, checked at
        ! every restart & every 50 products ... at most Mmax products ...
    subroutine lanczos_tr(n,apply,Ne,Mb_in,Mmax,Tol,E,Zv,Res,Nconv)
        integer, intent(in) :: n, Ne, Mb_in, Mmax
        procedure(apply_i) :: apply
        real(dp), intent(in) :: Tol
        real(dp), allocatable, intent(out) :: E(:), Zv(:,:), Res(:)
        integer, intent(out) :: Nconv
        real(dp), allocatable :: Vk(:,:), w(:), h(:), g(:), T(:,:), th(:), Y(:,:), Rs(:)
        integer :: mb, k, j, jm, nmv, i, kk, nrs
        real(dp) :: tw0, tmv, tor, t1
        real(dp) :: bet
        logical :: Inv, Early
        mb = Mb_in
        if (mb <= 0) mb = max(2*Ne + 50, 100)
        mb = min(max(mb, Ne + 10), n)
        allocate(Vk(n,mb+1), w(n), h(mb), g(mb), T(mb,mb))
        write(*,'(a)') TAB // 'thick-restart Lanczos ... basis ' // istr(mb) // ' vectors (' // &
                       cfmt_f(8.0_dp*real(n,dp)*(mb + 1)/1.0e9_dp,0,2) // ' GB), at most ' // istr(Mmax) // ' products'
        do i = 1, n
            Vk(i,1) = sin(12.9898_dp*real(i,dp) + 0.37_dp) + 0.5_dp*cos(78.233_dp*real(i,dp))
        end do
        Vk(:,1) = Vk(:,1) / norm2(Vk(:,1))
        T = 0.0_dp
        tw0 = omp_get_wtime()
        tmv = 0.0_dp
        tor = 0.0_dp
        k = 0
        nmv = 0
        nrs = 0
        Nconv = 0
        bet = 0.0_dp
        jm = 0
        do
            Inv = .false.
            Early = .false.
            do j = k + 1, mb
                t1 = omp_get_wtime()
                call apply(Vk(:,j),w)
                tmv = tmv + omp_get_wtime() - t1
                t1 = omp_get_wtime()
                nmv = nmv + 1
                call cgs_pass(n,j,Vk,w,h)
                call cgs_pass(n,j,Vk,w,g)
                h(1:j) = h(1:j) + g(1:j)
                T(1:j,j) = h(1:j)
                T(j,1:j) = h(1:j)
                bet = sqrt(pdot(n,w,w))
                tor = tor + omp_get_wtime() - t1
                jm = j
                if (bet <= 1.0e-12_dp*max(1.0_dp,abs(h(j)))) then
                    Inv = .true.
                    exit
                end if
                call pscale(n,1.0_dp/bet,w,Vk(:,j+1))
                if (nmv >= Mmax) exit
                if (mod(nmv,50) == 0 .and. j < mb .and. j >= Ne) then
                    call ritz(j)
                    if (Nconv >= Ne) then
                        Early = .true.
                        exit
                    end if
                end if
            end do
            call ritz(jm)
            write(*,'(a)') TAB // 'Lanczos restart ' // istr(nrs) // ', products ' // istr(nmv) // ' ... E_1 = ' // &
                           cfmt_f(th(1),0,8) // ', E_' // istr(min(Ne,jm)) // ' = ' // cfmt_f(th(min(Ne,jm)),0,6) // &
                           ', converged ' // istr(Nconv) // ' of ' // istr(Ne) // ', ' // cfmt_f(omp_get_wtime() - tw0,0,1) // ' s'
            if (Nconv >= Ne .or. nmv >= Mmax .or. Inv .or. Early) exit
            ! thick restart ...
            kk = min(Ne + (mb - Ne)/2, mb - 5)
            call rotate(kk)
            Vk(:,kk+1) = Vk(:,jm+1)
            T = 0.0_dp
            do i = 1, kk
                T(i,i) = th(i)
            end do
            k = kk
            nrs = nrs + 1
        end do
        write(*,'(a)') TAB // '  Lanczos time: products ' // cfmt_f(tmv,0,1) // ' s (' // cfmt_f(tmv/max(nmv,1),0,3) // ' s each), orthogonalization ' // &
                       cfmt_f(tor,0,1) // ' s, total ' // cfmt_f(omp_get_wtime() - tw0,0,1) // ' s'
        k = min(Ne,jm)
        allocate(E(k), Res(k))
        E = th(1:k)
        Res = Rs(1:k)
        call rotate(k)
        allocate(Zv(n,k))
        Zv = Vk(:,1:k)
        deallocate(Vk)
        do i = 1, k
            Zv(:,i) = Zv(:,i) / norm2(Zv(:,i))
        end do
    contains
        subroutine ritz(mx)
            integer, intent(in) :: mx
            if (allocated(th)) deallocate(th, Y, Rs)
            allocate(th(mx), Y(mx,mx), Rs(mx))
            call eig_sym(T(1:mx,1:mx),th,Y)
            Rs = abs(bet*Y(mx,:))
            Nconv = count(Rs(1:min(Ne,mx)) < Tol)
        end subroutine ritz
        ! Vk(:,1:kx) <- Vk(:,1:jm) Y(1:jm,1:kx) in row chunks ...
        subroutine rotate(kx)
            integer, intent(in) :: kx
            real(dp), allocatable :: Tm(:,:)
            integer :: i0, i1, nr
            !$omp parallel private(Tm,i0,i1,nr)
            allocate(Tm(1024,kx))
            !$omp do schedule(static)
            do i0 = 1, n, 1024
                i1 = min(n,i0 + 1023)
                nr = i1 - i0 + 1
                call dgemm('N','N',nr,kx,jm,1.0_dp,Vk(i0,1),n,Y,jm,0.0_dp,Tm,1024)
                Vk(i0:i1,1:kx) = Tm(1:nr,1:kx)
            end do
            !$omp end do
            deallocate(Tm)
            !$omp end parallel
        end subroutine rotate
    end subroutine lanczos_tr

end module nj_ncsm
