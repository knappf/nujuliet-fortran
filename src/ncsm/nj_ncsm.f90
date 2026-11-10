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
        ! 'EMPM'  H = E_HF + sum e_a :a+a: + 1/4 sum V2B_HF :a+a+aa: (h_HF diagonal, residual 3-body dropped, as EMPM &
        !         tests/fci_sparse_cm.py) -> plain form u = e - V rho, C0 = -sum e rho + 1/2 V rho rho in the HF basis,
        !         transformed back to HO with C_HF (u -> C u C^T, V2B_HF by O2b_transformation with C^T) ...
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
    use nj_orbitals, only : orbitals_make
    use nj_o1b, only : O1b_import, T1b
    use nj_o2b, only : O2b_import, O2b_transformation, O2b_pp, O2b_pn, O2b_nn
    use nj_v2b, only : V2b_read
    use nj_linalg, only : eig_sym
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
    end type sd_t

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
    ! HF single-particle energies e(a,t) (diagonal of h_HF, set by hamiltonian_empm) ...
    real(dp), allocatable, save :: EpsHF(:,:)

    type :: eng_t
        type(col_t), allocatable :: Mp(:), Mn(:)
        logical :: Stored = .false., Single = .false.
        integer(i8) :: nnz = 0
        integer :: nt = 0
        integer, allocatable :: Cnt(:), ChT(:), ChS(:), ChP(:), ChE(:)   ! chunk k: owner thread, segment & positions (entries, escapes) at its start
        type(tstore_t), allocatable :: Ts(:)
        real(dp), allocatable :: Yb(:,:)
    end type eng_t

    ! Context of the H product passed to lanczos_tr (apply_h) ... module-level pointers instead of an internal procedure
        ! (an internal procedure as an actual argument needs an executable-stack trampoline: SIGSEGV with ifx) ...
    type(sd_t), pointer, save :: cP => null(), cNn => null()
    type(pairs_t), pointer, save :: cPr(:) => null()
    type(mb_t), pointer, save :: cB => null()
    type(op_t), pointer, save :: cH => null()
    type(eng_t), pointer, save :: cEH => null()

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
        Orb = orbitals_make(Params)
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

        do iN = 1, size(Params%Calc%NCSM%Nmax_List)
            Nmax = Params%Calc%NCSM%Nmax_List(iN)
            if (Nmax < 0) exit
            if (Trunc == 'npnh') then
                ! all m-states of the model space, at most Nmax particles outside the reference, both parities ...
                call sp_make(Orb,Params%Calc%Nmax,SPS)
                do Par = 0, 1
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
        integer :: Lg, n

        A = Params%Calc%A
        Z = Params%Calc%Z
        Nneu = A - Z
        M2 = Params%Calc%NCSM%M2
        if (M2 < 0) M2 = mod(A,2)
        Beta = Params%Calc%NCSM%Beta_CM

        ! Slater determinants ...
        call timer_start(t0)
        Qcap = Nmax
        Hcap = huge(1)
        if (Trunc == 'npnh') then
            Qcap = 100000
            Hcap = Nmax
        end if
        ! Energy cut ('npnh', Basis 'HF') ... Ecn(n) = Ecut_npnh(min(n, last given)), Ecn(0) none; e = HF single-particle
            ! energies (diagonal of h_HF); SDs pruned by their own excitation energy when every species has
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
            if (.not. (gap(EpsP,Z) .and. gap(EpsN,Nneu))) then
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
                           ' MeV), HF e: reference <= ' // cfmt_f(max(maxval(EpsP(1:Z)),maxval(EpsN(1:Nneu))),0,3) // &
                           ', particles >= ' // cfmt_f(min(minval(EpsP(Z+1:)),minval(EpsN(Nneu+1:))),0,3) // ' MeV'
            call sd_make(SPS,Z,Qcap,Hcap,P,EpsP,ExmP)
            call sd_make(SPS,Nneu,Qcap,Hcap,Nn,EpsN,ExmN)
            call timer_report(t0,'proton & neutron SDs & jumps (energy cut)')
        else
        call sd_make(SPS,Z,Qcap,Hcap,P)
        if (Nneu == Z) then
            Nn = P            ! identical SDs & jumps (same s.p. space & truncation)
            call timer_report(t0,'proton SDs & jumps (neutrons copied)')
        else
            call sd_make(SPS,Nneu,Qcap,Hcap,Nn)
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
        call basis_make(P,Nn,M2 + 2,Nmax,mod(Par + P%N0 + Nn%N0,2),B2)
        call timer_report(t0,'basis & one-body jumps')
        call mem_report('SDs, jumps & basis')
        if (B%dim == 0) return

        ! Pair blocks within the quanta cut: the largest quanta of two nucleons in the basis ...
        Epair(1) = merge(SPS%e(max(Z,1)) + SPS%e(max(Z-1,1)), -1000, Z >= 2)
        Epair(2) = merge(SPS%e(max(Nneu,1)) + SPS%e(max(Nneu-1,1)), -1000, Nneu >= 2)
        Epair(3) = SPS%e(Z) + SPS%e(Nneu)
        Epair = min(Epair + Qcap, 2*maxval(SPS%e))
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
            call eng_make(SPS,P,Nn,Pr,H,B,EH,Single,Qcap,Hcap)
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

        ! Lanczos for H + Beta N_cm ...
        Ne = min(Params%Calc%NCSM%N_Eig, B%dim)
        call timer_start(t0)
        cP => P
        cNn => Nn
        cPr => Pr
        cB => B
        cH => H
        cEH => EH
        call lanczos_tr(B%dim,apply_h,Ne,Params%Calc%NCSM%Lanczos_Basis,Params%Calc%NCSM%Lanczos_Max, &
                        Params%Calc%NCSM%Lanczos_Tol,E,Zv,Res,Nconv)
        nullify(cP, cNn, cPr, cB, cH, cEH)
        call timer_report(t0,'Lanczos')
        call mem_report('Lanczos')
        Ne = size(E)
        call eng_free(EH)
        call op_free(H)

        ! <N_cm> (2-body, lower-half expectation pass) & <J^2> (one-body: |J_+ z|^2 + M (M + 1)) ...
        call timer_start(t0)
        allocate(Ev(Ne), Jv(Ne), Nv(Ne), Xt(Ne,B%dim))
        Xt = transpose(Zv)
        call pn_tables(SPS,Pr(3),Ncm)
        if (UseCut) then
            call eng_make(SPS,P,Nn,Pr,Ncm,B,EO,.false.,Qcap,Hcap,EpsP,EpsN,ExmP,ExmN)
        else
            call eng_make(SPS,P,Nn,Pr,Ncm,B,EO,.false.,Qcap,Hcap)
        end if
        call eng_expect(P,Nn,B,Ncm,EO,Ne,Xt,Nv)
        call op_free(Ncm)
        call eng_free(EO)
        call j2_expect(SPS,P,Nn,B,B2,Ne,Xt,Jv)
        deallocate(Xt)
        Ev = E - Beta*Nv
        call timer_report(t0,'<N_cm> & <J^2>')

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

    contains

        real(dp) function jval(x)
            real(dp), intent(in) :: x
            jval = 0.5_dp*(-1.0_dp + sqrt(max(1.0_dp + 4.0_dp*x,0.0_dp)))
        end function jval

        ! all reference energies e(1:nt) <= all others ...
        logical function gap(e,nt)
            real(dp), intent(in) :: e(:)
            integer, intent(in) :: nt
            gap = .true.
            if (nt < size(e)) gap = maxval(e(1:nt)) <= minval(e(nt+1:))
        end function gap

    end subroutine run_nmax

    ! y = (H + Beta N_cm) x ... the product of run_nmax for lanczos_tr (context in cP, cNn, cPr, cB, cH, cEH) ...
    subroutine apply_h(x,y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        if (cEH%Stored) then
            call eng_mv_stored(cEH,cB%dim,x,y)
        else
            call eng_mv_otf(cP,cNn,cB,cH,cEH,x,y)
        end if
    end subroutine apply_h

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
        call V2b_read(Params,Orb,V,Orb_NN)
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
        type(o1b_t) :: h, CHF, Ct, Rho
        real(dp), allocatable :: uF(:,:,:), occ(:,:), Cm(:,:)
        integer :: a_max, a, b, c, J, t, Pab
        real(dp) :: s

        if (.not. is_dir(io_path(Params,'/Bin'))) call fatal('NCSM EMPM Hamiltonian ... no HF run in ' // io_path(Params,''))
        call O2b_import(Params,Orb,io_path(Params,'/Bin/V2B_HF.bin'),V,Orb_NN)
        h = O1b_import(Params,Orb,io_path(Params,'/Bin/h_HF.bin'))
        CHF = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        Rho = O1b_import(Params,Orb,io_path(Params,'/Bin/Rho_HF.bin'))
        E_HF = read_ehf(io_path(Params,'/HF/HF_Summary.dat'))
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(occ(a_max,2), uF(a_max,a_max,2))
        if (allocated(EpsHF)) deallocate(EpsHF)
        allocate(EpsHF(a_max,2))
        do a = 1, a_max
            occ(a,1) = Rho%p(a,a)
            occ(a,2) = Rho%n(a,a)
            EpsHF(a,1) = h%p(a,a)
            EpsHF(a,2) = h%n(a,a)
        end do

        ! u = e - sum_c n_c sum_J (2J+1)/(2j_a+1) V^J(ac,bc), C0 = -sum (2j+1) n e + 1/2 sum n n (2J+1) V^J ...
        uF = 0.0_dp
        C0 = 0.0_dp
        do t = 1, 2
            do a = 1, a_max
                if (t == 1) uF(a,a,t) = h%p(a,a)
                if (t == 2) uF(a,a,t) = h%n(a,a)
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

        ! Basis 'HF' ... H stays in the HF orbitals, C_HF kept for N_cm (b & N transformed to HF) ...
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
                       ' MeV (relative to E_HF), u & V2B_HF transformed to HO'
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

    subroutine sp_make(Orb,emax,SPS)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: emax
        type(sp_t), intent(out) :: SPS
        integer :: a, m, k, Pass, i
        do Pass = 1, 2
            k = 0
            do a = 1, size(Orb)
                if (2*Orb(a)%n + Orb(a)%l > emax) cycle
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
            ! species copies, in the HF orbitals b -> C^T b C (C(i,k) = C_HF(a_i,a_k) for the same l, j, m) ...
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
    pure integer(i8) function sd_key(SPS,oc)
        type(sp_t), intent(in) :: SPS
        integer, intent(in) :: oc(:)
        integer :: ii
        sd_key = 0
        do ii = 1, size(oc)
            sd_key = sd_key + SPS%Bin(oc(ii) - 1,ii)
        end do
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
        sd_tgt = hashmap_lookup(SDX%Map,sd_key(SPS,tmp))
        xs = 1 - 2*mod(ii - 1 + nb2,2)
    end function sd_tgt

    ! Eps (optional, energy cut): HF single-particle energies of the m-states of this species, Exmax(h) the largest
        ! excitation energy sum(Eps(occ)) - sum(Eps(1:Nt)) kept for SDs with h particles outside the reference ...
    subroutine sd_make(SPS,Nt,Nmax,Hcap,SDX,Eps,Exmax)
        type(sp_t), intent(in) :: SPS
        integer, intent(in) :: Nt, Nmax, Hcap
        type(sd_t), intent(out) :: SDX
        real(dp), intent(in), optional :: Eps(:), Exmax(0:)
        logical :: Cut
        real(dp) :: Eref, exv, eminp
        real(dp), allocatable :: epl(:), exm(:), cum(:)
        integer, allocatable :: o(:), Cnt(:)
        integer :: Pass, n, s, i, kc, ka, d, sp, x
        integer(i8) :: j, Tot, tt
        integer, allocatable :: tmp(:)
        logical, allocatable :: Occ(:)
        call timer_start(tt)
        SDX%Nt = Nt
        SDX%N0 = sum(SPS%e(1:Nt))
        ! keys = combinatorial rank of the sorted occupations, 0 .. C(ns,Nt) - 1 ...
        if (Nt > KBIN) call fatal('NCSM ... SD keys overflow (more than KBIN particles of a species)')
        if (SPS%Bin(SPS%ns,Nt) < 0) call fatal('NCSM ... SD keys overflow (C(ns,Nt) too large)')
        allocate(o(Nt))
        Cut = present(Eps)
        Eref = 0.0_dp
        if (Cut) then
            epl = Eps
            exm = Exmax
            Eref = sum(epl(1:Nt))
            ! early pruning: cum(k) = sum of the reference energies e(1:k), eminp the lowest particle energy ...
            allocate(cum(0:Nt))
            cum(0) = 0.0_dp
            do i = 1, Nt
                cum(i) = cum(i-1) + epl(i)
            end do
            eminp = huge(1.0_dp)
            if (Nt < SPS%ns) eminp = minval(epl(Nt+1:))
        end if
        do Pass = 1, 2
            n = 0
            call rec(1,1,0)
            if (Pass == 1) allocate(SDX%occ(Nt,n), SDX%q(n), SDX%m2(n), SDX%h(n), SDX%L(n))
            if (Pass == 1 .and. Cut) allocate(SDX%ex(n))
        end do
        SDX%nsd = n
        call timer_report(tt,'  SDs: generation')
        call timer_start(tt)
        call SDX%Map%init(n)
        do s = 1, n
            call SDX%Map%set(key(SDX%occ(:,s)), s)
        end do
        call timer_report(tt,'  SDs: hash map')
        call timer_start(tt)

        ! One-body jumps ... a+_kc a_ka |s> for occupied ka & kc empty or = ka ... threaded over s: counts per bucket
        ! (pass 1, into Jo(b+1,s)), prefix sum, fill (pass 2) ... Cnt & tmp private ... with the energy cut a target
        ! above Exmax(h') (or with h' > Hcap) is outside the list: skipped before the lookup ...
        SDX%nd = SPS%nbk
        allocate(SDX%Jo(SDX%nd + 1,n))
        SDX%Jo = 0
        !$omp parallel private(s,i,ka,kc,d,x,Cnt,tmp,Occ,exv)
        allocate(Cnt(SDX%nd), tmp(Nt), Occ(SPS%ns))
        Occ = .false.
        !$omp do schedule(dynamic,256)
        do s = 1, n
            Cnt = 0
            do i = 1, Nt
                Occ(SDX%occ(i,s)) = .true.
            end do
            do i = 1, Nt
                ka = SDX%occ(i,s)
                do kc = 1, SPS%ns
                    if (kc /= ka .and. Occ(kc)) cycle
                    if (Cut .and. kc /= ka) then
                        if (.not. inlist(s,ka,kc)) cycle
                    end if
                    if (sd_tgt(SPS,SDX,Nmax,s,i,kc,x,tmp) < 0) cycle
                    d = SPS%Bk(kc,ka)
                    Cnt(d) = Cnt(d) + 1
                end do
            end do
            do i = 1, Nt
                Occ(SDX%occ(i,s)) = .false.
            end do
            do d = 1, SDX%nd
                SDX%Jo(d+1,s) = int(Cnt(d),i8)
            end do
        end do
        !$omp end do
        deallocate(Cnt, tmp, Occ)
        !$omp end parallel
        Tot = 0
        do s = 1, n
            SDX%Jo(1,s) = Tot + 1
            do d = 1, SDX%nd
                SDX%Jo(d+1,s) = SDX%Jo(d,s) + SDX%Jo(d+1,s)
            end do
            Tot = SDX%Jo(SDX%nd+1,s) - 1
        end do
        call timer_report(tt,'  SDs: jumps pass 1 (count) & prefix')
        call timer_start(tt)
        allocate(SDX%Jt(Tot), SDX%Jx(Tot))
        !$omp parallel private(s,i,ka,kc,d,x,sp,j,Cnt,tmp,Occ,exv)
        allocate(Cnt(SDX%nd), tmp(Nt), Occ(SPS%ns))
        Occ = .false.
        !$omp do schedule(dynamic,256)
        do s = 1, n
            Cnt = 0
            do i = 1, Nt
                Occ(SDX%occ(i,s)) = .true.
            end do
            do i = 1, Nt
                ka = SDX%occ(i,s)
                do kc = 1, SPS%ns
                    if (kc /= ka .and. Occ(kc)) cycle
                    if (Cut .and. kc /= ka) then
                        if (.not. inlist(s,ka,kc)) cycle
                    end if
                    sp = sd_tgt(SPS,SDX,Nmax,s,i,kc,x,tmp)
                    if (sp < 0) cycle
                    d = SPS%Bk(kc,ka)
                    j = SDX%Jo(d,s) + Cnt(d)
                    Cnt(d) = Cnt(d) + 1
                    SDX%Jt(j) = x*sp
                    SDX%Jx(j) = int(SPS%Ix(kc,ka) - IXOFF,i2)
                end do
            end do
            do i = 1, Nt
                Occ(SDX%occ(i,s)) = .false.
            end do
        end do
        !$omp end do
        deallocate(Cnt, tmp, Occ)
        !$omp end parallel
        call timer_report(tt,'  SDs: jumps pass 2 (fill), ' // istr8(Tot) // ' jumps (' // &
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



        ! combinatorial rank of the ascending occupations oc (distinct, 1..ns) ...
        integer(i8) function key(oc)
            integer, intent(in) :: oc(:)
            integer :: ii
            key = 0
            do ii = 1, size(oc)
                key = key + SPS%Bin(oc(ii) - 1,ii)
            end do
        end function key

        ! a+_kc a_ka |s> within the energy cut of the list (h' <= Hcap & ex' <= Exmax(h')) ...
        logical function inlist(s_,ka_,kc_)
            integer, intent(in) :: s_, ka_, kc_
            integer :: hh
            hh = SDX%h(s_) - merge(1,0,ka_ > Nt) + merge(1,0,kc_ > Nt)
            inlist = .false.
            if (hh > Hcap .or. hh > ubound(exm,1)) return
            inlist = SDX%ex(s_) - epl(ka_) + epl(kc_) <= exm(hh) + 1.0e-9_dp
        end function inlist
    end subroutine sd_make

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
    subroutine species_matrix(SPS,SDX,Pr,Op,t,Single,Rank,Nmaxq,Hcap,Mt,Eps,Exmax)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: SDX
        type(pairs_t), intent(in) :: Pr
        type(op_t), intent(in) :: Op
        integer, intent(in) :: t, Nmaxq, Hcap
        logical, intent(in) :: Single
        integer, intent(in) :: Rank(:)
        type(col_t), allocatable, intent(out) :: Mt(:)
        real(dp), intent(in), optional :: Eps(:), Exmax(0:)
        integer, allocatable :: Pos(:), Rl(:), o(:)
        real(dp), allocatable :: Acc(:), Ep(:), Eb(:)
        logical, allocatable :: Occ(:)
        type(ivec_t), allocatable :: Ord(:)
        logical :: Cut
        real(dp) :: exlim, base, v, sg
        integer :: s, i1_, i2_, k1, k2, k3, k4, pid, b, q, qq, pq, sp, nl, Nt, nb1, nb2, ii, p3, p4, hh, h2, qn, h0, rs, b0, nsel
        integer(i8) :: j
        Nt = SDX%Nt
        Cut = present(Eps)
        allocate(Mt(SDX%nsd))
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
        !$omp parallel private(Pos,Rl,Acc,o,Occ,j,i1_,i2_,k1,k2,k3,k4,pid,b,q,qq,pq,sp,nl,v,sg,nb1,nb2,ii,p3,p4,hh,h2,qn,h0,rs) &
        !$omp firstprivate(base)
        allocate(Pos(SDX%nsd), Rl(4096), Acc(4096), o(Nt), Occ(SPS%ns))
        Pos = 0
        Occ = .false.
        !$omp do schedule(dynamic,16)
        do s = 1, SDX%nsd
            nl = 0
            rs = Rank(s)
            do ii = 1, Nt
                Occ(SDX%occ(ii,s)) = .true.
            end do
            ! one-body (the diagonal bucket) ...
            do j = SDX%Jo(b0,s), SDX%Jo(b0+1,s) - 1
                sp = abs(SDX%Jt(j))
                if (Rank(sp) < rs) cycle
                v = Op%uT(int(SDX%Jx(j)) + IXOFF,t)
                if (v /= 0.0_dp) call acc_add(Pos,Rl,Acc,nl,sp,v*real(sign(1,SDX%Jt(j)),dp))
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
                        hh = h0 - merge(1,0,k1 > Nt) - merge(1,0,k2 > Nt)
                        qn = SDX%q(s) - SPS%e(k1) - SPS%e(k2)
                        if (Cut) base = SDX%ex(s) - Eps(k1) - Eps(k2)
                        do qq = 1, size(Pr%mem(b)%v)
                            if (Cut) then
                                q = Ord(b)%v(qq)
                                pq = Pr%mem(b)%v(q)
                                if (base + Ep(pq) > exlim + 1.0e-9_dp) exit
                            else
                                q = qq
                                pq = Pr%mem(b)%v(q)
                            end if
                            k3 = Pr%k1(pq)
                            k4 = Pr%k2(pq)
                            if (k3 /= k1 .and. k3 /= k2) then
                                if (Occ(k3)) cycle
                            end if
                            if (k4 /= k1 .and. k4 /= k2) then
                                if (Occ(k4)) cycle
                            end if
                            if (qn + SPS%e(k3) + SPS%e(k4) > Nmaxq) cycle
                            h2 = hh + merge(1,0,k3 > Nt) + merge(1,0,k4 > Nt)
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
                            sp = hashmap_lookup(SDX%Map,keyo(o))
                            if (sp < 0) cycle
                            if (Rank(sp) < rs) cycle
                            p3 = 0
                            p4 = 0
                            do ii = 1, Nt
                                if (o(ii) == k3) p3 = ii
                                if (o(ii) == k4) p4 = ii
                            end do
                            nb2 = p3 + p4 - 3
                            sg = real(1 - 2*mod(nb1 + nb2,2),dp)
                            call acc_add(Pos,Rl,Acc,nl,sp,sg*v)
                        end do
                    end do
                end do
            end if
            do ii = 1, Nt
                Occ(SDX%occ(ii,s)) = .false.
            end do
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
        ! combinatorial rank of the ascending occupations oc (distinct, 1..ns) ...
        integer(i8) function keyo(oc)
            integer, intent(in) :: oc(:)
            integer :: ii
            keyo = 0
            do ii = 1, size(oc)
                keyo = keyo + SPS%Bin(oc(ii) - 1,ii)
            end do
        end function keyo
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

    subroutine eng_make(SPS,P,Nn,Pr,Op,B,Eng,Single,Nmaxq,Hcap,EpsP,EpsN,ExmP,ExmN)
        type(sp_t), intent(in) :: SPS
        type(sd_t), intent(in) :: P, Nn
        type(pairs_t), intent(in) :: Pr(3)
        type(op_t), intent(in) :: Op
        type(mb_t), intent(in) :: B
        type(eng_t), intent(out) :: Eng
        logical, intent(in) :: Single
        integer, intent(in) :: Nmaxq, Hcap
        real(dp), intent(in), optional :: EpsP(:), EpsN(:), ExmP(0:), ExmN(0:)
        integer, allocatable :: Ip(:)
        integer :: i
        integer(i8) :: npp, nnn
        allocate(Ip(P%nsd))
        do i = 1, P%nsd
            Ip(i) = i
        end do
        Eng%Single = Single
        if (present(EpsP)) then
            call species_matrix(SPS,P,Pr(1),Op,1,Single,Ip,Nmaxq,Hcap,Eng%Mp,EpsP,ExmP)
            call species_matrix(SPS,Nn,Pr(2),Op,2,Single,B%Rn,Nmaxq,Hcap,Eng%Mn,EpsN,ExmN)
        else
            call species_matrix(SPS,P,Pr(1),Op,1,Single,Ip,Nmaxq,Hcap,Eng%Mp)
            call species_matrix(SPS,Nn,Pr(2),Op,2,Single,B%Rn,Nmaxq,Hcap,Eng%Mn)
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
        integer :: sp, sn, j, bk, bn, st, sp2, sn2, r, Ln, L2, rn, ixp
        integer(i8) :: jp, jn
        real(dp) :: v, xc, sgp
        sp = B%sp(c)
        sn = B%sn(c)
        Ln = Nn%L(sn)
        rn = B%Rn(sn)
        xc = 0.0_dp
        if (Mode == 3) xc = X(1,c)
        if (Op%C0 /= 0.0_dp) call put(c,Op%C0)
        ! protons (partners sp2 >= sp) ...
        do j = 1, size(Eng%Mp(sp)%r)
            sp2 = Eng%Mp(sp)%r(j)
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
                sgp = real(sign(1,P%Jt(jp)),dp)
                ixp = int(P%Jx(jp)) + IXOFF
                do jn = Nn%Jo(bn,sn), Nn%Jo(bn+1,sn) - 1
                    v = Op%W(bk)%v(int(Nn%Jx(jn)) + IXOFF,ixp)
                    if (v == 0.0_dp) cycle
                    sn2 = abs(Nn%Jt(jn))
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
                call dgemv('T',n,j,1.0_dp,Vk,n,w,1,0.0_dp,h,1)
                call dgemv('N',n,j,-1.0_dp,Vk,n,h,1,1.0_dp,w,1)
                call dgemv('T',n,j,1.0_dp,Vk,n,w,1,0.0_dp,g,1)
                call dgemv('N',n,j,-1.0_dp,Vk,n,g,1,1.0_dp,w,1)
                h(1:j) = h(1:j) + g(1:j)
                T(1:j,j) = h(1:j)
                T(j,1:j) = h(1:j)
                bet = norm2(w)
                tor = tor + omp_get_wtime() - t1
                jm = j
                if (bet <= 1.0e-12_dp*max(1.0_dp,abs(h(j)))) then
                    Inv = .true.
                    exit
                end if
                Vk(:,j+1) = w / bet
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
