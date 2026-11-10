! n-phonon EMPM (n = 3, 4) ... basis |a lambda; J^P> of the physical (n-1)-phonon states a & the TDA phonons lambda,
    ! counterpart of eqm3_admat & phon_dens2 (ro7) of the original EMPM codes ... the (n-1)-phonon densities are computed on
    ! the fly from the (n-1)-phonon amplitudes C & X & the densities of the (n-2)-phonon states (the "lower level": the
    ! phonons themselves for n = 3, the 2-phonon states for n = 4), nothing is stored on disk ...
    ! (n-1)-phonon densities (ro7) ... row state a (amplitudes C^a), column state b (amplitudes X^b), both over the bases
        ! |mu nu> (mu the (n-2)-phonon alpha slot, nu the phonon lambda slot), Qa = sqrt(2Ja+1) ...
        ! rho_ab(s;i1,i2) = sum_(mu nu) X^b(mu nu) T^a_(mu nu)(s;i1,i2)
        ! T^a_(mu nu)(s) = sum_al (-1)^(j_mu+Ja+j_nu+s) Qa {s j_mu j_al; j_nu Ja Jb} C^a(al nu) rho_(al mu)(s)
        !                + sum_l1 (-1)^(j_mu+Jb+j_l1+s) Qa {s j_nu j_l1; j_mu Ja Jb} C^a(mu l1) rho_(l1 nu)(s)
        ! with the (n-2)-phonon densities rho_(al mu) of the lower level & the 1-phonon densities rho_(l1 nu) ...
    ! A & D rows (a lambda) as in the 2-phonon case with rho_(alpha beta) -> rho_ab & E_alpha -> E_a, the metric exchange
        ! term of eqm3_admat D += sum_g (-1)^(j_l'+j_l+Ja+Jb) {j_l j_g Jb; j_l' J Ja} X^a(g l') X^b(g l), g over the lower
        ! level ...
    ! Spurious states (Ortho) ... the (n-1)-phonon set is complete: the physical states (eigenstates of P H P) & the
        ! spurious ones (S H S) of each channel, with H^(n-1) = diag(E) + the couplings R(physical, spurious) ... the
        ! A rows get d_ll' H^(n-1)_ab + d_ab R_ll' (R_ll' = <l|A_TDA|cm> of the phonons) as in the 2-phonon case ... spurious
        ! n-phonon states |a lambda> with lambda = cm or a spurious, the n-phonon spurious subspace from their D elements ...
module nj_empm_eqm3
    use nj_kinds
    use nj_types
    use nj_utils
    use nj_wigner, only : f6j
    use nj_empm_phonons, only : empm_phonons_t
    use nj_empm_dens
    use nj_empm_eqm2
    implicit none
    private

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
        subroutine dpotrf(Uplo,n,A,Lda,Info)
            import :: dp
            character :: Uplo
            integer :: n, Lda, Info
            real(dp) :: A(Lda,*)
        end subroutine dpotrf
        subroutine dpotrs(Uplo,n,Nrhs,A,Lda,B,Ldb,Info)
            import :: dp
            character :: Uplo
            integer :: n, Nrhs, Lda, Ldb, Info
            real(dp) :: A(Lda,*), B(Ldb,*)
        end subroutine dpotrs
    end interface

    ! Physical n-phonon states of all channels (also the phonons themselves as a level without bases) ... N states, channel (J,P) holds First(J,P) .. First(J,P)+N_ch(J,P)-1 ...
        ! per channel the n-phonon basis Bs(J,P), the amplitudes C on the independent rows only (Kc(J,P)%v(r) full-basis index
        ! of independent row r, ascending, Cc(J,P)%v(r,i) its amplitude in state i, C = 0 on all other rows) &
        ! X(J,P)%v(k,i) = sqrt(2J+1) sum_(independent rows r) D_r(k) C(r,i) ...
        ! N_phys(J,P) physical states first in each channel, then the spurious ones (Spur), Rc(J,P)%v(i,s) = <i|H|s> of
        ! the physical state i & the spurious state s of the channel ...
    type, public :: empm_states2_t
        integer :: N = 0, J_max = 0
        integer, allocatable :: J(:), P(:), N_ch(:,:), First(:,:), N_phys(:,:)
        real(dp), allocatable :: E(:)
        logical, allocatable :: Spur(:)
        type(empm_basis2_t), allocatable :: Bs(:,:)
        type(ivec_t), allocatable :: Kc(:,:)
        type(rmat_t), allocatable :: Cc(:,:), X(:,:), Rc(:,:)
        ! Xt ... X transposed (n x N_full, rows X(k,:) contiguous for the exchange term of the rows; EMPM_states_read) ...
        type(rmat_t), allocatable :: Xt(:,:)
    end type empm_states2_t

    ! 1-phonon densities R, transposed densities Rt & redefined interactions V of all phonons ...
    type, public :: empm_cache1_t
        type(phonon_block_t), allocatable :: R(:), Rt(:), V(:)
    end type empm_cache1_t

    ! Marker of the format of EMPM_states_write (-2 the earlier one without spurious states) ...
    integer, parameter :: States_Format = -3

    ! Wall-clock seconds spent in the (n-1)-phonon densities, row assembly & A D^T products of the current channel ...
    real(dp), save :: T_Rho2 = 0.0_dp, T_Row3 = 0.0_dp, T_H3 = 0.0_dp
    ! profile of EMPM_densityN (thread-seconds): W6 table, T assembly, dgemm with X^b ...
    integer, save :: Dens_Alg = 1, Dens_Phase_Mode = 2, Dens_Nprint = 0, Dens_Nblk = 0
    logical, save :: Dens_Init = .false., Dens_Gemm = .false.
    real(dp), save :: Dens_Emax = 0.0_dp
    real(dp), save :: T_rG = 0.0_dp, T_rX = 0.0_dp, T_rS = 0.0_dp
    real(dp), save :: T_dW = 0.0_dp, T_dT = 0.0_dp, T_dG = 0.0_dp
    integer(i8), save :: F_dG = 0_i8

    type, public :: greedy_res_t
        integer, allocatable :: Ind(:), Ord(:)
        real(dp), allocatable :: Dii(:,:), H(:,:), Dr(:,:), Dsp(:,:)
        integer :: N_visit = 0
    end type greedy_res_t

    public :: EMPM_states1_make, EMPM_cache1_make, EMPM_basis3_make, EMPM_densityN, EMPM_eqm3_build, EMPM_D3_build, EMPM_H3_ind, EMPM_ind_greedy, &
              EMPM_ind_greedy_H, EMPM_ind_greedy_multi, EMPM_AD3_rows, EMPM_states_write, EMPM_states_write_D, EMPM_states_write_Df, EMPM_states_read

contains

    pure integer function par_sign(P)
        integer, intent(in) :: P
        par_sign = merge(1, -1, P == 1)
    end function par_sign

    ! Densities, transposed densities & redefined interactions of all phonons, in parallel ...
    subroutine EMPM_cache1_make(Pairs,Ph,C1)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(out) :: C1
        integer :: g
        allocate(C1%R(Ph%N), C1%Rt(Ph%N), C1%V(Ph%N))
        !$omp parallel do schedule(dynamic)
        do g = 1, Ph%N
            call EMPM_density(Pairs,Ph,g,C1%R(g))
            call EMPM_density_transposed(Pairs,C1%R(g),C1%Rt(g))
            call EMPM_interaction(Pairs,Ph,C1%R(g),C1%V(g))
        end do
        !$omp end parallel do
    end subroutine EMPM_cache1_make

    ! The phonons as the lowest level of states (N, J, P, E, channels; no bases) ...
    subroutine EMPM_states1_make(Ph,S1)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_states2_t), intent(out) :: S1
        integer :: g
        S1%N = Ph%N
        S1%J_max = Ph%J_max
        S1%J = Ph%J
        S1%P = Ph%P
        S1%E = Ph%E
        allocate(S1%N_ch(0:Ph%J_max,2), S1%First(0:Ph%J_max,2))
        S1%N_ch = Ph%N_ch(0:Ph%J_max,:)
        S1%First = Ph%First(0:Ph%J_max,:)
        S1%Spur = [(g == Ph%i_cm, g = 1, Ph%N)]
    end subroutine EMPM_states1_make

    ! n-phonon basis |a lambda; J^P> ... all couplings (Al = (n-1)-phonon state a, La = phonon lambda, a outer), truncated
        ! rows E_a + E_lambda <= E3max ... spurious states (lambda = cm or a spurious) all rows for Spur_Rows = 'all', no
        ! rows otherwise ...
    subroutine EMPM_basis3_make(Ph,S2,J,P,E3max,Spur_Rows,B3)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_states2_t), intent(in) :: S2
        integer, intent(in) :: J, P
        real(dp), intent(in) :: E3max
        character(len=*), intent(in) :: Spur_Rows
        type(empm_basis2_t), intent(out) :: B3
        integer :: a, l, k, Pass

        B3%J = J
        B3%P = P
        allocate(B3%K(S2%N,Ph%N))
        B3%K = 0
        do Pass = 1, 2
            k = 0
            do a = 1, S2%N
                do l = 1, Ph%N
                    if (par_sign(S2%P(a))*par_sign(Ph%P(l)) /= par_sign(P)) cycle
                    if (J < abs(S2%J(a) - Ph%J(l)) .or. J > S2%J(a) + Ph%J(l)) cycle
                    k = k + 1
                    if (Pass == 2) then
                        B3%Al(k) = a
                        B3%La(k) = l
                        B3%K(a,l) = k
                        B3%Spur(k) = (l == Ph%i_cm .or. S2%Spur(a))
                    end if
                end do
            end do
            if (Pass == 1) then
                B3%N_full = k
                allocate(B3%Al(k), B3%La(k), B3%Spur(k))
                B3%Spur = .false.
            end if
        end do
        B3%Ks = pack([(k, k = 1, B3%N_full)], B3%Spur)
        B3%N_spur = size(B3%Ks)
        allocate(B3%Row(B3%N_full))
        B3%N_tr = 0
        do k = 1, B3%N_full
            if (B3%Spur(k)) then
                if (trim(Spur_Rows) /= 'all') cycle
                B3%N_spur_tr = B3%N_spur_tr + 1
            else if (S2%E(B3%Al(k)) + Ph%E(B3%La(k)) > E3max) then
                cycle
            end if
            B3%N_tr = B3%N_tr + 1
            B3%Row(B3%N_tr) = k
        end do
        B3%Row = B3%Row(1:B3%N_tr)
    end subroutine EMPM_basis3_make

    ! Densities of a for all b ... NJ_DENS = old (EMPM_densityN_old: T^a X^b over the full basis of b), new (default,
        ! EMPM_densityN_T: U^a C^b over the independent rows of b), check (both, deviations reported, old returned) ...
    subroutine EMPM_densityN(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Sig_Only,Jb_Range)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        integer, intent(in) :: a
        type(phonon_block_t), intent(out) :: Rho2
        integer, intent(in), optional :: Sig_Only
        integer, intent(in), optional :: Jb_Range(2)
        type(phonon_block_t) :: Rn
        integer :: Jb, Pb, Sig
        real(dp) :: r, e, m
        character(len=32) :: Val
        integer :: L

        if (.not. Dens_Init) then
            !$omp critical (njd_init)
            if (.not. Dens_Init) then
                call get_environment_variable('NJ_DENS',Val,L)
                if (L > 0) then
                    select case (trim(Val))
                    case ('old')
                        Dens_Alg = 0
                    case ('check')
                        Dens_Alg = 2
                    case default
                        Dens_Alg = 1
                    end select
                end if
                call get_environment_variable('NJ_DENS_GEMM',Val,L)
                if (L > 0) Dens_Gemm = trim(Val) /= '0'
                call get_environment_variable('NJ_DENS_PHASE',Val,L)
                if (L > 0) read(Val,*) Dens_Phase_Mode
                write(*,'(a)') TAB // 'Densities ... ' // trim(merge('old  ', 'new  ', Dens_Alg == 0)) // &
                               merge(' (check against old)', '                    ', Dens_Alg == 2) // ' , phase mode ' // istr(Dens_Phase_Mode) // ' , U by dgemm ' // merge('yes','no ',Dens_Gemm)
                Dens_Init = .true.
            end if
            !$omp end critical (njd_init)
        end if
        if (Dens_Alg == 0 .or. present(Sig_Only)) then
            call EMPM_densityN_old(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Sig_Only,Jb_Range)
            return
        end if
        if (Dens_Alg == 1) then
            call EMPM_densityN_T(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Jb_Range)
            return
        end if
        call EMPM_densityN_old(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Jb_Range=Jb_Range)
        call EMPM_densityN_T(Pairs,Ph,C1,Rlow,Slow,S2,a,Rn,Jb_Range)
        do Jb = lbound(Rho2%B,1), ubound(Rho2%B,1)
            do Pb = 1, 2
                do Sig = Rho2%Lo(Jb,Pb), Rho2%Hi(Jb,Pb)
                    associate(O => Rho2%B(Jb,Pb)%M(Sig)%v, N => Rn%B(Jb,Pb)%M(Sig)%v)
                        if (size(O) == 0) cycle
                        m = maxval(abs(O))
                        r = sum(O*N)/max(sum(N*N),1.0e-300_dp)
                        e = maxval(abs(O - N))/max(m,1.0e-300_dp)
                        !$omp critical (njd_chk)
                        Dens_Nblk = Dens_Nblk + 1
                        Dens_Emax = max(Dens_Emax, merge(e, 0.0_dp, m > 1.0e-10_dp))
                        if (Dens_Nprint < 200 .and. m > 1.0e-10_dp .and. e > 1.0e-9_dp) then
                            Dens_Nprint = Dens_Nprint + 1
                            write(*,'(a,6i4,3es12.3)') ' DENSCHK a Ja Jb Pa Pb Sig  ratio  rel.err  max ', a, S2%J(a), Jb, S2%P(a), Pb, Sig, r, e, m
                        end if
                        if (mod(Dens_Nblk,20000) == 0) write(*,'(a,i10,es12.3)') ' DENSCHK blocks , max rel. deviation ', Dens_Nblk, Dens_Emax
                        !$omp end critical (njd_chk)
                    end associate
                end do
            end do
        end do
    end subroutine EMPM_densityN

    ! (n-1)-phonon densities rho_ab(s; pairs) of the row state a for all (n-1)-phonon states b (S2) ...
        ! Rho2%B(Jb,Pb)%M(s)%v(c,ib) on the allowed pairs c of (s, parity), as the 1-phonon blocks ... T^a per (s, channel of b)
        ! then one dgemm with X^b ... Rlow & Slow the densities & states of the lower level ... Serial (called inside the
        ! parallel loops of the callers) ... Sig_Only: only the rank s = Sig_Only (the recursion keeps s, so Rlow may hold
        ! only that rank too) ...
    subroutine EMPM_densityN_old(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Sig_Only,Jb_Range)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        integer, intent(in) :: a
        type(phonon_block_t), intent(out) :: Rho2
        integer, intent(in), optional :: Sig_Only
        ! Jb_Range ... only the channels Jb_Range(1) <= Jb <= Jb_Range(2) of b (the others empty, as channels without states) ...
        integer, intent(in), optional :: Jb_Range(2)
        real(dp), allocatable :: Ca(:), T(:,:), W6(:,:,:,:), Rt(:,:)
        integer, allocatable :: Pn(:), Kn(:), Pm(:), Km(:), Cn(:), Cm(:), Nz(:), Off(:)
        integer :: Slo, Shi, Kt, ip, jal
        real(dp) :: base
        integer :: Ja, Pa, ia, Jb, Pb, Sig, Ka, Pp, nb, k, kk, ic, mu, nu, al, l1, jmu, jnu, x, y, z, n1, nlo, nz_n, J1
        real(dp) :: Qa, w
        integer(i8) :: tp

        Ja = S2%J(a)
        Pa = S2%P(a)
        ia = a - S2%First(Ja,Pa) + 1
        Qa = sqrt(real(2*Ja + 1,dp))
        n1 = Ph%N
        nlo = Slow%N
        J1 = max(Ph%J_max, Slow%J_max)
        associate(Ba => S2%Bs(Ja,Pa))
            ! C^a on the full basis (zero off the independent rows Kc) ...
            allocate(Ca(Ba%N_full))
            Ca = 0.0_dp
            Ca(S2%Kc(Ja,Pa)%v) = S2%Cc(Ja,Pa)%v(:,ia)
            ! Non-zero components of C^a grouped by the lambda slot (Kn(Pn(g)+1 .. Pn(g+1)) with La = g) & by the alpha slot
                ! (Km(Pm(g)+1 .. Pm(g+1)) with Al = g) ...
            Nz = pack(S2%Kc(Ja,Pa)%v, S2%Cc(Ja,Pa)%v(:,ia) /= 0.0_dp)
            nz_n = size(Nz)
            allocate(Pn(n1+1), Pm(nlo+1), Kn(nz_n), Km(nz_n))
            Pn = 0
            Pm = 0
            do kk = 1, nz_n
                Pn(Ba%La(Nz(kk))+1) = Pn(Ba%La(Nz(kk))+1) + 1
                Pm(Ba%Al(Nz(kk))+1) = Pm(Ba%Al(Nz(kk))+1) + 1
            end do
            do k = 1, n1
                Pn(k+1) = Pn(k+1) + Pn(k)
            end do
            do k = 1, nlo
                Pm(k+1) = Pm(k+1) + Pm(k)
            end do
            Cn = Pn
            Cm = Pm
            do kk = 1, nz_n
                Cn(Ba%La(Nz(kk))) = Cn(Ba%La(Nz(kk))) + 1
                Kn(Cn(Ba%La(Nz(kk)))) = Nz(kk)
                Cm(Ba%Al(Nz(kk))) = Cm(Ba%Al(Nz(kk))) + 1
                Km(Cm(Ba%Al(Nz(kk)))) = Nz(kk)
            end do

            Rho2%g = a
            allocate(Rho2%Lo(0:S2%J_max,2), Rho2%Hi(0:S2%J_max,2), Rho2%Pi(0:S2%J_max,2), Rho2%B(0:S2%J_max,2))
            allocate(W6(0:J1,0:J1,0:J1,0:Pairs%L_max))
            do Jb = 0, S2%J_max
                do Pb = 1, 2
                    nb = S2%N_ch(Jb,Pb)
                    if (present(Jb_Range)) then
                        if (Jb < Jb_Range(1) .or. Jb > Jb_Range(2)) nb = 0
                    end if
                    Rho2%Lo(Jb,Pb) = abs(Ja - Jb)
                    Rho2%Hi(Jb,Pb) = min(Ja + Jb, Pairs%L_max)
                    if (present(Sig_Only)) then
                        Rho2%Lo(Jb,Pb) = max(Rho2%Lo(Jb,Pb), Sig_Only)
                        Rho2%Hi(Jb,Pb) = min(Rho2%Hi(Jb,Pb), Sig_Only)
                    end if
                    Rho2%Pi(Jb,Pb) = merge(1, 2, Pa == Pb)
                    if (nb == 0) Rho2%Hi(Jb,Pb) = Rho2%Lo(Jb,Pb) - 1
                    allocate(Rho2%B(Jb,Pb)%M(Rho2%Lo(Jb,Pb):Rho2%Hi(Jb,Pb)))
                    if (nb == 0 .or. Rho2%Lo(Jb,Pb) > Rho2%Hi(Jb,Pb)) cycle
                    ! W6(x,y,z,s) = {s x y; z Ja Jb} ...
                    call timer_start(tp)
                    do Sig = Rho2%Lo(Jb,Pb), Rho2%Hi(Jb,Pb)
                        do z = 0, J1
                            do y = 0, J1
                                do x = 0, J1
                                    W6(x,y,z,Sig) = f6j(2*Sig,2*x,2*y,2*z,2*Ja,2*Jb)
                                end do
                            end do
                        end do
                    end do
                    !$omp atomic
                    T_dW = T_dW + elapsed_s(tp)
                    Pp = Rho2%Pi(Jb,Pb)
                    ! all ranks s of the channel stacked (rows Off(s)+1 .. Off(s+1) of T), one sweep over the basis of b
                        ! & one dgemm with X^b ...
                    Slo = Rho2%Lo(Jb,Pb)
                    Shi = Rho2%Hi(Jb,Pb)
                    allocate(Off(Slo:Shi+1))
                    Off(Slo) = 0
                    do Sig = Slo, Shi
                        Off(Sig+1) = Off(Sig) + Pairs%Na(Sig,Pp)
                        allocate(Rho2%B(Jb,Pb)%M(Sig)%v(Pairs%Na(Sig,Pp),nb))
                    end do
                    Kt = Off(Shi+1)
                    associate(Bb => S2%Bs(Jb,Pb))
                        if (Kt > 0) then
                            call timer_start(tp)
                            allocate(T(Kt,Bb%N_full))
                            T = 0.0_dp
                            do k = 1, Bb%N_full
                                mu = Bb%Al(k)
                                nu = Bb%La(k)
                                jmu = Slow%J(mu)
                                jnu = Ph%J(nu)
                                ! C^a(al nu) rho_(al mu) ...
                                ip = mu - Slow%First(jmu,Slow%P(mu)) + 1
                                do ic = Pn(nu) + 1, Pn(nu+1)
                                    kk = Kn(ic)
                                    al = Ba%Al(kk)
                                    jal = Slow%J(al)
                                    base = Qa * Ca(kk)
                                    do Sig = Slo, Shi
                                        if (Off(Sig+1) == Off(Sig)) cycle
                                        w = W6(jmu,jal,jnu,Sig)
                                        if (w == 0.0_dp) cycle
                                        w = w * phase(jmu + Ja + jnu + Sig) * base
                                        T(Off(Sig)+1:Off(Sig+1),k) = T(Off(Sig)+1:Off(Sig+1),k) + &
                                            w * Rlow(al)%B(jmu,Slow%P(mu))%M(Sig)%v(:,ip)
                                    end do
                                end do
                                ! C^a(mu l1) rho_(l1 nu) ...
                                ip = nu - Ph%First(jnu,Ph%P(nu)) + 1
                                do ic = Pm(mu) + 1, Pm(mu+1)
                                    kk = Km(ic)
                                    l1 = Ba%La(kk)
                                    jal = Ph%J(l1)
                                    base = Qa * Ca(kk)
                                    do Sig = Slo, Shi
                                        if (Off(Sig+1) == Off(Sig)) cycle
                                        w = W6(jnu,jal,jmu,Sig)
                                        if (w == 0.0_dp) cycle
                                        w = w * phase(jmu + Jb + jal + Sig) * base
                                        T(Off(Sig)+1:Off(Sig+1),k) = T(Off(Sig)+1:Off(Sig+1),k) + &
                                            w * C1%R(l1)%B(jnu,Ph%P(nu))%M(Sig)%v(:,ip)
                                    end do
                                end do
                            end do
                            !$omp atomic
                            T_dT = T_dT + elapsed_s(tp)
                            call timer_start(tp)
                            allocate(Rt(Kt,nb))
                            call dgemm('N','N',Kt,nb,Bb%N_full,1.0_dp,T,Kt,S2%X(Jb,Pb)%v,Bb%N_full,0.0_dp,Rt,Kt)
                            do Sig = Slo, Shi
                                Rho2%B(Jb,Pb)%M(Sig)%v = Rt(Off(Sig)+1:Off(Sig+1),:)
                            end do
                            !$omp atomic
                            T_dG = T_dG + elapsed_s(tp)
                            !$omp atomic
                            F_dG = F_dG + 2_i8*Kt*nb*Bb%N_full
                            deallocate(T, Rt)
                        else
                            do Sig = Slo, Shi
                                Rho2%B(Jb,Pb)%M(Sig)%v = 0.0_dp
                            end do
                        end if
                    end associate
                    deallocate(Off)
                end do
            end do
        end associate
    end subroutine EMPM_densityN_old

    ! (n-1)-phonon densities rho_ab(s) of the row state a for all b, as EMPM_densityN, computed from the other side ...
        ! <a| ... |b> with the ket b expanded on its independent rows (C^b, N_ind of the channel) & the bra a by its
        ! overlaps X^a (full basis) ... U^a(s; pairs, j) on the independent rows j = (alpha nu) of the channel of b, then one
        ! dgemm U^a C^b (N_ind instead of N_full columns: X^b of EMPM_densityN) ... the result conjugated (transposed pairs,
        ! Sg & the phase Dens_Phase) to rho_ab ...
    subroutine EMPM_densityN_T(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Jb_Range)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        integer, intent(in) :: a
        type(phonon_block_t), intent(out) :: Rho2
        integer, intent(in), optional :: Jb_Range(2)
        real(dp), allocatable :: Xa(:), U(:,:), W6(:,:,:,:), Rt(:,:)
        integer, allocatable :: Pn(:), Kn(:), Pm(:), Km(:), Cn(:), Cm(:), Nz(:), Off(:)
        integer :: Ja, Pa, ia, Jb, Pb, Sig, nb, k, kk, ic, mu, nu, al, l1, jmu, jnu, jal, jl1, n1, nlo, nz_n, J1, Slo, Shi, Kt, &
                   ip, jj, nK, Pp, c, s, ct
        real(dp) :: Qb, w, base, phz
        integer(i8) :: tp

        Ja = S2%J(a)
        Pa = S2%P(a)
        ia = a - S2%First(Ja,Pa) + 1
        n1 = Ph%N
        nlo = Slow%N
        J1 = max(Ph%J_max, Slow%J_max)
        associate(Ba => S2%Bs(Ja,Pa))
            ! non-zero overlaps X^a of the bra grouped by the phonon (Kn) & by the lower state (Km) ...
            Xa = S2%X(Ja,Pa)%v(:,ia)
            Nz = pack([(k, k = 1, Ba%N_full)], Xa /= 0.0_dp)
            nz_n = size(Nz)
            allocate(Pn(n1+1), Pm(nlo+1), Kn(nz_n), Km(nz_n))
            Pn = 0
            Pm = 0
            do kk = 1, nz_n
                Pn(Ba%La(Nz(kk))+1) = Pn(Ba%La(Nz(kk))+1) + 1
                Pm(Ba%Al(Nz(kk))+1) = Pm(Ba%Al(Nz(kk))+1) + 1
            end do
            do k = 1, n1
                Pn(k+1) = Pn(k+1) + Pn(k)
            end do
            do k = 1, nlo
                Pm(k+1) = Pm(k+1) + Pm(k)
            end do
            Cn = Pn
            Cm = Pm
            do kk = 1, nz_n
                Cn(Ba%La(Nz(kk))) = Cn(Ba%La(Nz(kk))) + 1
                Kn(Cn(Ba%La(Nz(kk)))) = Nz(kk)
                Cm(Ba%Al(Nz(kk))) = Cm(Ba%Al(Nz(kk))) + 1
                Km(Cm(Ba%Al(Nz(kk)))) = Nz(kk)
            end do

            Rho2%g = a
            allocate(Rho2%Lo(0:S2%J_max,2), Rho2%Hi(0:S2%J_max,2), Rho2%Pi(0:S2%J_max,2), Rho2%B(0:S2%J_max,2))
            allocate(W6(0:J1,0:J1,0:J1,0:Pairs%L_max))
            do Jb = 0, S2%J_max
                do Pb = 1, 2
                    nb = S2%N_ch(Jb,Pb)
                    if (present(Jb_Range)) then
                        if (Jb < Jb_Range(1) .or. Jb > Jb_Range(2)) nb = 0
                    end if
                    Rho2%Lo(Jb,Pb) = abs(Ja - Jb)
                    Rho2%Hi(Jb,Pb) = min(Ja + Jb, Pairs%L_max)
                    Rho2%Pi(Jb,Pb) = merge(1, 2, Pa == Pb)
                    if (nb == 0) Rho2%Hi(Jb,Pb) = Rho2%Lo(Jb,Pb) - 1
                    allocate(Rho2%B(Jb,Pb)%M(Rho2%Lo(Jb,Pb):Rho2%Hi(Jb,Pb)))
                    if (nb == 0 .or. Rho2%Lo(Jb,Pb) > Rho2%Hi(Jb,Pb)) cycle
                    Slo = Rho2%Lo(Jb,Pb)
                    Shi = Rho2%Hi(Jb,Pb)
                    Pp = Rho2%Pi(Jb,Pb)
                    ! W6(x,y,z,s) = {s x y; z Jb Ja} (ket b, bra a) ...
                    call timer_start(tp)
                    do Sig = Slo, Shi
                        do jj = 0, J1
                            do k = 0, J1
                                do c = 0, J1
                                    W6(c,k,jj,Sig) = f6j(2*Sig,2*c,2*k,2*jj,2*Jb,2*Ja)
                                end do
                            end do
                        end do
                    end do
                    !$omp atomic
                    T_dW = T_dW + elapsed_s(tp)
                    allocate(Off(Slo:Shi+1))
                    Off(Slo) = 0
                    do Sig = Slo, Shi
                        Off(Sig+1) = Off(Sig) + Pairs%Na(Sig,Pp)
                        allocate(Rho2%B(Jb,Pb)%M(Sig)%v(Pairs%Na(Sig,Pp),nb))
                    end do
                    Kt = Off(Shi+1)
                    Qb = sqrt(real(2*Jb + 1,dp))
                    associate(Bb => S2%Bs(Jb,Pb), Kc => S2%Kc(Jb,Pb)%v)
                        nK = size(Kc)
                        if (Kt > 0 .and. nK > 0) then
                            call timer_start(tp)
                            allocate(U(Kt,nK))
                            U = 0.0_dp
                            if (Dens_Gemm) then
                                call U_gemm(Jb,Pb,Slo,Shi,Off,Kt,U)
                            else
                            do jj = 1, nK
                                k = Kc(jj)
                                al = Bb%Al(k)
                                l1 = Bb%La(k)
                                jal = Slow%J(al)
                                jl1 = Ph%J(l1)
                                ! ket (al l1), bra (mu l1) of a: rho^(n-2)_(al -> mu) ...
                                do ic = Pn(l1) + 1, Pn(l1+1)
                                    kk = Kn(ic)
                                    mu = Ba%Al(kk)
                                    jmu = Slow%J(mu)
                                    ip = mu - Slow%First(jmu,Slow%P(mu)) + 1
                                    base = Qb * Xa(kk)
                                    do Sig = Slo, Shi
                                        if (Off(Sig+1) == Off(Sig)) cycle
                                        w = W6(jmu,jal,jl1,Sig)
                                        if (w == 0.0_dp) cycle
                                        w = w * phase(jmu + Jb + jl1 + Sig) * base
                                        U(Off(Sig)+1:Off(Sig+1),jj) = U(Off(Sig)+1:Off(Sig+1),jj) + &
                                            w * Rlow(al)%B(jmu,Slow%P(mu))%M(Sig)%v(:,ip)
                                    end do
                                end do
                                ! ket (al l1), bra (al nu) of a: rho^1_(l1 -> nu) ...
                                do ic = Pm(al) + 1, Pm(al+1)
                                    kk = Km(ic)
                                    nu = Ba%La(kk)
                                    jnu = Ph%J(nu)
                                    ip = nu - Ph%First(jnu,Ph%P(nu)) + 1
                                    base = Qb * Xa(kk)
                                    do Sig = Slo, Shi
                                        if (Off(Sig+1) == Off(Sig)) cycle
                                        w = W6(jnu,jl1,jal,Sig)
                                        if (w == 0.0_dp) cycle
                                        w = w * phase(jal + Ja + jl1 + Sig) * base
                                        U(Off(Sig)+1:Off(Sig+1),jj) = U(Off(Sig)+1:Off(Sig+1),jj) + &
                                            w * C1%R(l1)%B(jnu,Ph%P(nu))%M(Sig)%v(:,ip)
                                    end do
                                end do
                            end do
                            end if
                            !$omp atomic
                            T_dT = T_dT + elapsed_s(tp)
                            call timer_start(tp)
                            allocate(Rt(Kt,nb))
                            call dgemm('N','N',Kt,nb,nK,1.0_dp,U,Kt,S2%Cc(Jb,Pb)%v,nK,0.0_dp,Rt,Kt)
                            ! <a| [pair s]^Sig |b> -> rho_ab ... transposed pair, Sg & the phase of the conjugation ...
                            do Sig = Slo, Shi
                                phz = Dens_Phase(Ja,Jb,Sig)
                                do ct = 1, Pairs%Na(Sig,Pp)
                                    s = Pairs%Al(Sig,Pp)%v(ct)
                                    c = Pairs%Ci(Pairs%Tr(s),Sig,Pp)
                                    if (c > 0) then
                                        Rho2%B(Jb,Pb)%M(Sig)%v(ct,:) = phz * Pairs%Sg(s) * Rt(Off(Sig)+c,:)
                                    else
                                        Rho2%B(Jb,Pb)%M(Sig)%v(ct,:) = 0.0_dp
                                    end if
                                end do
                            end do
                            !$omp atomic
                            T_dG = T_dG + elapsed_s(tp)
                            !$omp atomic
                            F_dG = F_dG + 2_i8*Kt*nb*nK
                            deallocate(U, Rt)
                        else
                            do Sig = Slo, Shi
                                Rho2%B(Jb,Pb)%M(Sig)%v = 0.0_dp
                            end do
                        end if
                    end associate
                    deallocate(Off)
                end do
            end do
        end associate

    contains

        ! U^a by dgemm ... term 1 per lower state al of the independent rows (consecutive columns) & channel of mu:
            ! rho^(n-2)_(al)(s; pairs, mu) X^a(mu nu), term 2 per phonon l1 & channel of nu: rho^1_(l1)(s; pairs, nu) X^a(al nu)
            ! (scattered columns) ... the basis of a holds all or none of the (mu nu) of a channel of mu (nu) ...
        subroutine U_gemm(Jb,Pb,Slo,Shi,Off,Ld,U)
            integer, intent(in) :: Jb, Pb, Slo, Shi, Off(Slo:Shi+1), Ld
            real(dp), intent(inout) :: U(Ld,*)
            real(dp), allocatable :: Y(:,:), Ys(:,:), Tq(:,:)
            integer, allocatable :: Lv(:), Pl(:), Byl(:), Cnt(:)
            integer :: nK, j0, j1, jj, al, jal, Jm, Pm_, nmu, m0, nv, iv, Sig, Ka, l1, jl1, i, Jn, Pn_, nnu, n0, m, q, k
            real(dp) :: cf, Qb
            logical :: Any_, Nz_
            Qb = sqrt(real(2*Jb + 1,dp))
            associate(Bb => S2%Bs(Jb,Pb), Kc => S2%Kc(Jb,Pb)%v, Ba => S2%Bs(Ja,Pa))
                nK = size(Kc)
                ! term 1 ... columns j0 .. j1 with the same lower state al ...
                j0 = 1
                do while (j0 <= nK)
                    al = Bb%Al(Kc(j0))
                    j1 = j0
                    do while (j1 < nK)
                        if (Bb%Al(Kc(j1+1)) /= al) exit
                        j1 = j1 + 1
                    end do
                    nv = j1 - j0 + 1
                    Lv = [(Bb%La(Kc(jj)), jj = j0, j1)]
                    jal = Slow%J(al)
                    do Jm = 0, Slow%J_max
                        do Pm_ = 1, 2
                            nmu = Slow%N_ch(Jm,Pm_)
                            if (nmu == 0) cycle
                            m0 = Slow%First(Jm,Pm_)
                            allocate(Y(nmu,nv))
                            Y = 0.0_dp
                            Any_ = .false.
                            do iv = 1, nv
                                if (Ba%K(m0,Lv(iv)) == 0) cycle
                                Any_ = .true.
                                do i = 1, nmu
                                    k = Ba%K(m0+i-1,Lv(iv))
                                    if (k > 0) Y(i,iv) = Xa(k)
                                end do
                            end do
                            if (Any_) then
                                allocate(Ys(nmu,nv))
                                do Sig = max(Slo,Rlow(al)%Lo(Jm,Pm_)), min(Shi,Rlow(al)%Hi(Jm,Pm_))
                                    Ka = Off(Sig+1) - Off(Sig)
                                    if (Ka == 0) cycle
                                    if (size(Rlow(al)%B(Jm,Pm_)%M(Sig)%v,1) /= Ka) cycle
                                    Nz_ = .false.
                                    do iv = 1, nv
                                        cf = W6(Jm,jal,Ph%J(Lv(iv)),Sig) * phase(Jm + Jb + Ph%J(Lv(iv)) + Sig) * Qb
                                        Ys(:,iv) = cf * Y(:,iv)
                                        if (cf /= 0.0_dp) Nz_ = .true.
                                    end do
                                    if (Nz_) call dgemm('N','N',Ka,nv,nmu,1.0_dp,Rlow(al)%B(Jm,Pm_)%M(Sig)%v,Ka,Ys,nmu,1.0_dp, &
                                                        U(Off(Sig)+1,j0),Ld)
                                end do
                                deallocate(Ys)
                            end if
                            deallocate(Y)
                        end do
                    end do
                    j0 = j1 + 1
                end do
                ! term 2 ... columns grouped by the phonon l1 ...
                allocate(Pl(Ph%N+1), Cnt(Ph%N), Byl(nK))
                Pl = 0
                do jj = 1, nK
                    Pl(Bb%La(Kc(jj))+1) = Pl(Bb%La(Kc(jj))+1) + 1
                end do
                Pl(1) = 0
                do l1 = 1, Ph%N
                    Pl(l1+1) = Pl(l1+1) + Pl(l1)
                end do
                Cnt = 0
                do jj = 1, nK
                    l1 = Bb%La(Kc(jj))
                    Cnt(l1) = Cnt(l1) + 1
                    Byl(Pl(l1)+Cnt(l1)) = jj
                end do
                do l1 = 1, Ph%N
                    m = Pl(l1+1) - Pl(l1)
                    if (m == 0) cycle
                    jl1 = Ph%J(l1)
                    do Jn = 0, Ph%J_max
                        do Pn_ = 1, 2
                            nnu = Ph%N_ch(Jn,Pn_)
                            if (nnu == 0) cycle
                            n0 = Ph%First(Jn,Pn_)
                            allocate(Y(nnu,m))
                            Y = 0.0_dp
                            Any_ = .false.
                            do q = 1, m
                                al = Bb%Al(Kc(Byl(Pl(l1)+q)))
                                if (Ba%K(al,n0) == 0) cycle
                                Any_ = .true.
                                do i = 1, nnu
                                    k = Ba%K(al,n0+i-1)
                                    if (k > 0) Y(i,q) = Xa(k)
                                end do
                            end do
                            if (Any_) then
                                allocate(Ys(nnu,m))
                                do Sig = max(Slo,C1%R(l1)%Lo(Jn,Pn_)), min(Shi,C1%R(l1)%Hi(Jn,Pn_))
                                    Ka = Off(Sig+1) - Off(Sig)
                                    if (Ka == 0) cycle
                                    if (size(C1%R(l1)%B(Jn,Pn_)%M(Sig)%v,1) /= Ka) cycle
                                    Nz_ = .false.
                                    do q = 1, m
                                        jal = Slow%J(Bb%Al(Kc(Byl(Pl(l1)+q))))
                                        cf = W6(Jn,jl1,jal,Sig) * phase(jal + Ja + jl1 + Sig) * Qb
                                        Ys(:,q) = cf * Y(:,q)
                                        if (cf /= 0.0_dp) Nz_ = .true.
                                    end do
                                    if (.not. Nz_) cycle
                                    allocate(Tq(Ka,m))
                                    call dgemm('N','N',Ka,m,nnu,1.0_dp,C1%R(l1)%B(Jn,Pn_)%M(Sig)%v,Ka,Ys,nnu,0.0_dp,Tq,Ka)
                                    do q = 1, m
                                        jj = Byl(Pl(l1)+q)
                                        U(Off(Sig)+1:Off(Sig+1),jj) = U(Off(Sig)+1:Off(Sig+1),jj) + Tq(:,q)
                                    end do
                                    deallocate(Tq)
                                end do
                                deallocate(Ys)
                            end if
                            deallocate(Y)
                        end do
                    end do
                end do
            end associate
        end subroutine U_gemm

    end subroutine EMPM_densityN_T

    ! Phase of the conjugation <a|O^s|b> -> <b|O^s|a> (with the transposed pairs & Sg) ... Dens_Phase_Mode selects the
        ! candidate during the check (NJ_DENS=check) ...
    real(dp) function Dens_Phase(Ja,Jb,Sig)
        integer, intent(in) :: Ja, Jb, Sig
        select case (Dens_Phase_Mode)
        case (1)
            Dens_Phase = 1.0_dp
        case (2)
            Dens_Phase = phase(Ja + Jb)
        case (3)
            Dens_Phase = phase(Ja + Jb + Sig)
        case (4)
            Dens_Phase = phase(Ja + Jb) * sqrt(real(2*Ja + 1,dp)/real(2*Jb + 1,dp))
        case (5)
            Dens_Phase = phase(Ja + Jb) * sqrt(real(2*Jb + 1,dp)/real(2*Ja + 1,dp))
        case default
            Dens_Phase = 1.0_dp
        end select
    end function Dens_Phase

    ! Independent rows without the whole truncated D (Ind_Greedy) ... the truncated rows grouped by their (n-1)-phonon state a,
        ! the groups in a fixed pseudo-random order, blocks of >= Nblk rows: full-basis D rows of the block (EMPM_AD3_rows, one
        ! density per a), Schur complement S = D_bb - D_bI D_II^-1 D_Ib against the rows I selected so far, pivoted Cholesky of
        ! S (absolute threshold Chol_Tol, as EMPM_chol_ind) adds its pivots to I ... stops when |I| = Target (the npnh
        ! dimension of the channel) or all rows are visited ... Ind ascending, Dii = D(Ind,Ind), N_visit rows visited ...
    subroutine EMPM_ind_greedy(Pairs,Ph,C1,Rlow,Slow,S2,B3,Target,Ind,Dii,N_visit)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B3
        integer, intent(in) :: Target
        integer, allocatable, intent(out) :: Ind(:)
        real(dp), allocatable, intent(out) :: Dii(:,:)
        integer, intent(out) :: N_visit
        real(dp), parameter :: Tol = 1.0e-6_dp
        integer, allocatable :: Cnt(:), Start(:), Rows(:), Perm(:), Blk(:), Sel(:), Piv(:), Ord(:), Tmp(:)
        real(dp), allocatable :: Db(:,:), Dbb(:,:), Cb(:,:), X(:,:), Lf(:,:), Dn(:,:), Key(:)
        integer :: N, Na, i, a, g, ng, nb, r, Rk, k, j, Info, Nblk
        integer(i8) :: Seed

        N = B3%N_tr
        Na = S2%N
        ! rows grouped by a ...
        allocate(Cnt(Na), Start(Na+1), Rows(N))
        Cnt = 0
        do i = 1, N
            a = B3%Al(B3%Row(i))
            Cnt(a) = Cnt(a) + 1
        end do
        Start(1) = 1
        do a = 1, Na
            Start(a+1) = Start(a) + Cnt(a)
        end do
        Cnt = 0
        do i = 1, N
            a = B3%Al(B3%Row(i))
            Rows(Start(a)+Cnt(a)) = i
            Cnt(a) = Cnt(a) + 1
        end do
        ! groups with rows in a fixed pseudo-random order (Fisher-Yates, LCG) ...
        Perm = pack([(a, a = 1, Na)], Start(2:Na+1) > Start(1:Na))
        ng = size(Perm)
        Seed = 20261103_i8
        do g = ng, 2, -1
            Seed = modulo(6364136223846793005_i8*Seed + 1442695040888963407_i8, huge(1_i8))
            j = 1 + int(modulo(Seed/65536_i8, int(g,i8)))
            k = Perm(g)
            Perm(g) = Perm(j)
            Perm(j) = k
        end do
        Nblk = max(256, Target/4)
        allocate(Sel(0), Dii(0,0))
        r = 0
        N_visit = 0
        g = 0
        do while (r < Target .and. g < ng)
            ! next block of whole groups ...
            allocate(Blk(0))
            do while (size(Blk) < Nblk .and. g < ng)
                g = g + 1
                a = Perm(g)
                Blk = [Blk, Rows(Start(a):Start(a+1)-1)]
            end do
            nb = size(Blk)
            N_visit = N_visit + nb
            allocate(Db(B3%N_full,nb))
            call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Blk,Db)
            Dbb = Db(B3%Row(Blk),:)
            if (r > 0) then
                Cb = Db(B3%Row(Sel),:)
                X = Cb
                Lf = Dii
                call dpotrf('L',r,Lf,r,Info)
                if (Info /= 0) call fatal('EMPM_ind_greedy ... D of the selected rows not positive definite (dpotrf ' // istr(Info) // ')')
                call dpotrs('L',r,nb,Lf,r,X,r,Info)
                call dgemm('T','N',nb,nb,r,-1.0_dp,Cb,r,X,r,1.0_dp,Dbb,nb)
                deallocate(X, Lf)
            end if
            call chol_piv(nb,Dbb,Tol,Piv,Rk)
            if (Rk > 0) then
                ! extended D(I,I) ...
                allocate(Dn(r+Rk,r+Rk))
                if (r > 0) then
                    Dn(1:r,1:r) = Dii
                    Dn(1:r,r+1:r+Rk) = Cb(:,Piv(1:Rk))
                    Dn(r+1:r+Rk,1:r) = transpose(Cb(:,Piv(1:Rk)))
                end if
                do j = 1, Rk
                    do i = 1, Rk
                        Dn(r+i,r+j) = Db(B3%Row(Blk(Piv(i))),Piv(j))
                    end do
                end do
                call move_alloc(Dn,Dii)
                Sel = [Sel, Blk(Piv(1:Rk))]
                r = r + Rk
            end if
            if (allocated(Cb)) deallocate(Cb)
            deallocate(Db, Dbb, Blk)
        end do
        ! ascending order ...
        allocate(Key(r))
        Key = real(Sel,dp)
        allocate(Ord(r))
        call sortperm(Key,Ord)
        Ind = Sel(Ord)
        Dii = Dii(Ord,Ord)
    end subroutine EMPM_ind_greedy

    ! Greedy independent rows & H = A D^T in one pass (Ind_Greedy) ... as EMPM_ind_greedy, but each block of rows gets its A
        ! rows too (one density per a, computed once): the D rows of the selected rows are kept (Dr, full basis, the X of the
        ! exported states), the A rows only per block: H(new, all selected) = A_new^T D_sel, H(old, new) by the symmetry of H
        ! (the asymmetry is checked on the diagonal blocks new x new) ... no second density pass, no stored A ... Dr is
        ! preallocated for Target rows, blocks of at most Nblk rows (Mem_GB / 4 for A & D of the block) ... Ind ascending,
        ! Dii = D(Ind,Ind), H & the columns of Dr in the same order ...
    ! Dr_File present: D of the selected rows not kept in memory (Dr not allocated) but written to Dr_File, column by column
        ! in the order of selection (Ord(i) = its column of the i-th row of the ascending Ind), H(new, old) streamed from it
        ! block by block ... Dsp = D(Ks, selected rows) (ascending Ind) of the spurious states returned instead ...
    subroutine EMPM_ind_greedy_H(Pairs,Ph,C1,Rlow,Slow,S2,B3,Target,Mem_GB,Ind,Dii,H,Dr,N_visit,Dr_File,Ord_File,Dsp)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B3
        integer, intent(in) :: Target
        real(dp), intent(in) :: Mem_GB
        integer, allocatable, intent(out) :: Ind(:)
        real(dp), allocatable, intent(out) :: Dii(:,:), H(:,:), Dr(:,:)
        integer, intent(out) :: N_visit
        character(len=*), intent(in), optional :: Dr_File
        integer, allocatable, intent(out), optional :: Ord_File(:)
        real(dp), allocatable, intent(out), optional :: Dsp(:,:)
        real(dp), parameter :: Tol = 1.0e-6_dp
        integer, allocatable :: Cnt(:), Start(:), Rows(:), Perm(:), Blk(:), Sel(:), Piv(:), Ord(:)
        real(dp), allocatable :: Db(:,:), Ab(:,:), Dbb(:,:), Cb(:,:), X(:,:), Lf(:,:), Key(:), Tmp(:,:), Col(:), Bf(:,:), Ds(:,:)
        integer :: N, Na, i, a, g, ng, nb, r, Rk, k, j, Info, Nblk, Nf, uD, q0, q1, Ns
        integer(i8) :: Seed, t1, t2
        logical, allocatable :: Done(:)
        logical :: Disk
        real(dp) :: T_io

        N = B3%N_tr
        Na = S2%N
        Nf = B3%N_full
        T_Rho2 = 0.0_dp
        T_Row3 = 0.0_dp
        T_H3 = 0.0_dp
        allocate(Cnt(Na), Start(Na+1), Rows(N))
        Cnt = 0
        do i = 1, N
            a = B3%Al(B3%Row(i))
            Cnt(a) = Cnt(a) + 1
        end do
        Start(1) = 1
        do a = 1, Na
            Start(a+1) = Start(a) + Cnt(a)
        end do
        Cnt = 0
        do i = 1, N
            a = B3%Al(B3%Row(i))
            Rows(Start(a)+Cnt(a)) = i
            Cnt(a) = Cnt(a) + 1
        end do
        ! the same pseudo-random order of the groups as EMPM_ind_greedy (same selected rows) ...
        Perm = pack([(a, a = 1, Na)], Start(2:Na+1) > Start(1:Na))
        ng = size(Perm)
        Seed = 20261103_i8
        do g = ng, 2, -1
            Seed = modulo(6364136223846793005_i8*Seed + 1442695040888963407_i8, huge(1_i8))
            j = 1 + int(modulo(Seed/65536_i8, int(g,i8)))
            k = Perm(g)
            Perm(g) = Perm(j)
            Perm(j) = k
        end do
        Nblk = max(256, Target/4)
        Nblk = min(Nblk, max(256, int(0.25_dp*Mem_GB*1.0e9_dp/(16.0_dp*real(Nf,dp)))))
        Disk = present(Dr_File)
        T_io = 0.0_dp
        if (Disk) then
            Ns = B3%N_spur
            allocate(Sel(Target), Dii(Target,Target), H(Target,Target), Ds(Ns,Target))
            open(newunit=uD, file=Dr_File, access='stream', form='unformatted', status='replace', action='readwrite')
        else
            allocate(Sel(Target), Dii(Target,Target), H(Target,Target), Dr(Nf,Target))
        end if
        r = 0
        N_visit = 0
        g = 0
        do while (r < Target .and. g < ng)
            allocate(Blk(0))
            do while (size(Blk) < Nblk .and. g < ng)
                g = g + 1
                a = Perm(g)
                Blk = [Blk, Rows(Start(a):Start(a+1)-1)]
            end do
            nb = size(Blk)
            N_visit = N_visit + nb
            allocate(Db(Nf,nb), Ab(Nf,nb))
            call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Blk,Db,Ab)
            Dbb = Db(B3%Row(Blk),:)
            if (r > 0) then
                Cb = Db(B3%Row(Sel(1:r)),:)
                X = Cb
                Lf = Dii(1:r,1:r)
                call dpotrf('L',r,Lf,r,Info)
                if (Info /= 0) call fatal('EMPM_ind_greedy_H ... D of the selected rows not positive definite (dpotrf ' // istr(Info) // ')')
                call dpotrs('L',r,nb,Lf,r,X,r,Info)
                call dgemm('T','N',nb,nb,r,-1.0_dp,Cb,r,X,r,1.0_dp,Dbb,nb)
                deallocate(X, Lf)
            end if
            call chol_piv(nb,Dbb,Tol,Piv,Rk)
            Rk = min(Rk, Target - r)
            if (Rk > 0) then
                if (r > 0) then
                    Dii(1:r,r+1:r+Rk) = Cb(:,Piv(1:Rk))
                    Dii(r+1:r+Rk,1:r) = transpose(Cb(:,Piv(1:Rk)))
                end if
                do j = 1, Rk
                    do i = 1, Rk
                        Dii(r+i,r+j) = Db(B3%Row(Blk(Piv(i))),Piv(j))
                    end do
                end do
                ! H(new, 1 .. r+Rk) = A_new^T D_sel, H(old, new) = H(new, old)^T ...
                call timer_start(t1)
                allocate(Tmp(Nf,Rk))
                !$omp parallel do
                do j = 1, Rk
                    Tmp(:,j) = Ab(:,Piv(j))
                end do
                !$omp end parallel do
                if (Disk) then
                    deallocate(Ab)
                    ! H(new, old) from the file in chunks of Nblk columns, H(new, new) from the block ...
                    if (r > 0) then
                        allocate(Bf(Nf,min(r,Nblk)))
                        do q0 = 1, r, Nblk
                            q1 = min(r, q0 + Nblk - 1)
                            call timer_start(t2)
                            read(uD, pos=1_i8 + 8_i8*int(Nf,i8)*int(q0-1,i8)) Bf(:,1:q1-q0+1)
                            T_io = T_io + elapsed_s(t2)
                            call dgemm('T','N',Rk,q1-q0+1,Nf,1.0_dp,Tmp,Nf,Bf,Nf,0.0_dp,H(r+1,q0),Target)
                        end do
                        deallocate(Bf)
                    end if
                    allocate(Bf(Nf,Rk))
                    !$omp parallel do
                    do j = 1, Rk
                        Bf(:,j) = Db(:,Piv(j))
                    end do
                    !$omp end parallel do
                    do j = 1, Rk
                        do i = 1, Ns
                            Ds(i,r+j) = Bf(B3%Ks(i),j)
                        end do
                    end do
                    call dgemm('T','N',Rk,Rk,Nf,1.0_dp,Tmp,Nf,Bf,Nf,0.0_dp,H(r+1,r+1),Target)
                    call timer_start(t2)
                    write(uD, pos=1_i8 + 8_i8*int(Nf,i8)*int(r,i8)) Bf
                    T_io = T_io + elapsed_s(t2)
                    deallocate(Bf)
                else
                    !$omp parallel do
                    do j = 1, Rk
                        Dr(:,r+j) = Db(:,Piv(j))
                    end do
                    !$omp end parallel do
                    call dgemm('T','N',Rk,r+Rk,Nf,1.0_dp,Tmp,Nf,Dr,Nf,0.0_dp,H(r+1,1),Target)
                end if
                deallocate(Tmp)
                if (r > 0) H(1:r,r+1:r+Rk) = transpose(H(r+1:r+Rk,1:r))
                T_H3 = T_H3 + elapsed_s(t1)
                Sel(r+1:r+Rk) = Blk(Piv(1:Rk))
                r = r + Rk
            end if
            if (allocated(Cb)) deallocate(Cb)
            if (allocated(Ab)) deallocate(Ab)
            deallocate(Db, Dbb, Blk)
        end do
        ! ascending order ... H, Dii & the columns of Dr permuted (Dr in place, column cycles) ...
        allocate(Key(r), Ord(r))
        Key = real(Sel(1:r),dp)
        call sortperm(Key,Ord)
        Ind = Sel(Ord)
        Dii = Dii(Ord,Ord)
        H = H(Ord,Ord)
        if (Disk) then
            close(uD)
            if (present(Ord_File)) Ord_File = Ord
            if (present(Dsp)) then
                allocate(Dsp(Ns,r))
                do j = 1, r
                    Dsp(:,j) = Ds(:,Ord(j))
                end do
            end if
            write(*,'(a)') TAB // 'D of the selected rows on disk ... ' // cfmt_f(8.0_dp*real(Nf,dp)*real(r,dp)/1.0e9_dp,0,2) // &
                           ' GB , I/O ' // cfmt_f(T_io,0,1) // ' s , ' // Dr_File
        else
        allocate(Done(r), Col(Nf))
        Done = .false.
        do i = 1, r
            if (Done(i) .or. Ord(i) == i) then
                Done(i) = .true.
                cycle
            end if
            ! cycle i <- Ord(i) <- Ord(Ord(i)) ... (new column j = old column Ord(j)) ...
            Col = Dr(:,i)
            j = i
            do
                Done(j) = .true.
                k = Ord(j)
                if (k == i) then
                    Dr(:,j) = Col
                    exit
                end if
                Dr(:,j) = Dr(:,k)
                j = k
            end do
        end do
        if (r < Target) then
            ! fewer independent rows than the npnh dimension (all rows visited) ...
            allocate(Tmp(Nf,r))
            Tmp = Dr(:,1:r)
            call move_alloc(Tmp,Dr)
        end if
        end if
        write(*,'(a)') TAB // 'Profile of the densities (thread-s) ... W6 ' // cfmt_f(T_dW,0,1) // ' , T ' // cfmt_f(T_dT,0,1) // &
                       ' , dgemm ' // cfmt_f(T_dG,0,1) // ' (' // cfmt_f(real(F_dG,dp)/max(T_dG,1.0e-9_dp)/1.0e9_dp,0,2) // ' GFlop/s per thread) ... rows: dgemm ' // cfmt_f(T_rG,0,1) // ' , exchange ' // cfmt_f(T_rX,0,1) // ' , scatter ' // cfmt_f(T_rS,0,1)
        write(*,'(a)') TAB // 'Seconds ... greedy A & D rows: densities ' // cfmt_f(T_Rho2,0,1) // ' (summed over threads) , rows ' // &
                       cfmt_f(T_Row3,0,1) // ' (summed over threads) , A D^T ' // cfmt_f(T_H3,0,1) // ' , block rows ' // istr(Nblk)
    end subroutine EMPM_ind_greedy_H

    ! Greedy independent rows & H = A D^T as EMPM_ind_greedy_H, for several channels B(1:nch) at once with shared densities
        ! (Dens_Share) ... one common pseudo-random order of the (n-1)-phonon states a, each block of a's gets its densities
        ! rho_ab once (the b channels of the channels that still need rows) & the A & D rows of all those channels from them ...
        ! the memory of the blocks split among the channels ... Disk(c): D of the selected rows of channel c in D_Files(c)
        ! (Res(c)%Ord, Res(c)%Dsp) instead of Res(c)%Dr ... with one channel the same rows as EMPM_ind_greedy_H ...
    subroutine EMPM_ind_greedy_multi(Pairs,Ph,C1,Rlow,Slow,S2,B,Target,Mem_GB,Disk,D_Files,Res)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B(:)
        integer, intent(in) :: Target(:)
        real(dp), intent(in) :: Mem_GB
        logical, intent(in) :: Disk(:)
        character(len=*), intent(in) :: D_Files(:)
        type(greedy_res_t), intent(out) :: Res(:)
        real(dp), parameter :: Tol = 1.0e-6_dp
        type :: gch_t
            integer, allocatable :: Start(:), Rows(:), Sel(:), Blk(:), Col0(:)
            real(dp), allocatable :: Ds(:,:), Db(:,:), Ab(:,:)
            integer :: r = 0, Nblk = 0, uD = -1, Nf = 0, Ns = 0, Cnt = 0
        end type gch_t
        type(gch_t), allocatable :: Gs(:)
        integer, allocatable :: Perm(:), Ablk(:)
        logical, allocatable :: Act(:), Need(:)
        type(phonon_block_t) :: Rho2
        integer :: nch, c, Na, i, a, k, j, g, ng, nab, ig, iq, q, Jlo, Jhi
        integer(i8) :: Seed, t1
        real(dp) :: tr, tw, T_io

        nch = size(B)
        Na = S2%N
        T_Rho2 = 0.0_dp
        T_Row3 = 0.0_dp
        T_H3 = 0.0_dp
        T_io = 0.0_dp
        allocate(Gs(nch), Act(nch), Need(nch))
        do c = 1, nch
            associate(Gc => Gs(c), B3 => B(c))
                Gc%Nf = B3%N_full
                Gc%Ns = B3%N_spur
                allocate(Gc%Start(Na+1), Gc%Rows(B3%N_tr))
                ! rows grouped by a (Start(a) .. Start(a+1)-1 in Rows) ...
                Gc%Start = 0
                do i = 1, B3%N_tr
                    a = B3%Al(B3%Row(i))
                    Gc%Start(a+1) = Gc%Start(a+1) + 1
                end do
                Gc%Start(1) = 1
                do a = 1, Na
                    Gc%Start(a+1) = Gc%Start(a+1) + Gc%Start(a)
                end do
                block
                    integer, allocatable :: Cn(:)
                    allocate(Cn(Na))
                    Cn = 0
                    do i = 1, B3%N_tr
                        a = B3%Al(B3%Row(i))
                        Gc%Rows(Gc%Start(a)+Cn(a)) = i
                        Cn(a) = Cn(a) + 1
                    end do
                end block
                Gc%Nblk = max(256, Target(c)/4)
                Gc%Nblk = min(Gc%Nblk, max(256, int(0.25_dp*Mem_GB*1.0e9_dp/(16.0_dp*real(Gc%Nf,dp)*real(nch,dp)))))
                allocate(Gc%Sel(Target(c)), Res(c)%Dii(Target(c),Target(c)), Res(c)%H(Target(c),Target(c)))
                if (Disk(c)) then
                    allocate(Gc%Ds(Gc%Ns,Target(c)))
                    open(newunit=Gc%uD, file=D_Files(c), access='stream', form='unformatted', status='replace', action='readwrite')
                else
                    allocate(Res(c)%Dr(Gc%Nf,Target(c)))
                end if
                Res(c)%N_visit = 0
            end associate
        end do
        ! common pseudo-random order of the a's with rows in any channel (the order of EMPM_ind_greedy for one channel) ...
        allocate(Perm(0))
        do a = 1, Na
            if (any([(Gs(c)%Start(a+1) > Gs(c)%Start(a), c = 1, nch)])) Perm = [Perm, a]
        end do
        ng = size(Perm)
        Seed = 20261103_i8
        do g = ng, 2, -1
            Seed = modulo(6364136223846793005_i8*Seed + 1442695040888963407_i8, huge(1_i8))
            j = 1 + int(modulo(Seed/65536_i8, int(g,i8)))
            k = Perm(g)
            Perm(g) = Perm(j)
            Perm(j) = k
        end do

        g = 0
        do
            Act = [(Gs(c)%r < Target(c), c = 1, nch)]
            if (.not. any(Act) .or. g >= ng) exit
            ! block of a's ... until an active channel has its block size ...
            allocate(Ablk(0))
            do c = 1, nch
                Gs(c)%Cnt = 0
            end do
            do while (g < ng)
                if (any([(Act(c) .and. Gs(c)%Cnt >= Gs(c)%Nblk, c = 1, nch)])) exit
                g = g + 1
                a = Perm(g)
                Ablk = [Ablk, a]
                do c = 1, nch
                    if (Act(c)) Gs(c)%Cnt = Gs(c)%Cnt + Gs(c)%Start(a+1) - Gs(c)%Start(a)
                end do
            end do
            nab = size(Ablk)
            do c = 1, nch
                associate(Gc => Gs(c))
                    if (.not. Act(c) .or. Gc%Cnt == 0) cycle
                    allocate(Gc%Blk(Gc%Cnt), Gc%Col0(nab+1), Gc%Db(Gc%Nf,Gc%Cnt), Gc%Ab(Gc%Nf,Gc%Cnt))
                    q = 0
                    do ig = 1, nab
                        a = Ablk(ig)
                        Gc%Col0(ig) = q
                        do i = Gc%Start(a), Gc%Start(a+1) - 1
                            q = q + 1
                            Gc%Blk(q) = Gc%Rows(i)
                        end do
                    end do
                    Gc%Col0(nab+1) = q
                    Res(c)%N_visit = Res(c)%N_visit + Gc%Cnt
                end associate
            end do
            ! densities of each a once, the rows of all active channels ...
            tr = 0.0_dp
            tw = 0.0_dp
            !$omp parallel do schedule(dynamic) private(a,Rho2,t1,Need,Jlo,Jhi,c,iq,q) reduction(+:tr,tw)
            do ig = 1, nab
                a = Ablk(ig)
                Need = [(Act(c) .and. Gs(c)%Start(a+1) > Gs(c)%Start(a), c = 1, nch)]
                if (.not. any(Need)) cycle
                Jlo = huge(1)
                Jhi = -1
                do c = 1, nch
                    if (.not. Need(c)) cycle
                    Jlo = min(Jlo, max(0, B(c)%J - Ph%J_max))
                    Jhi = max(Jhi, B(c)%J + Ph%J_max)
                end do
                call timer_start(t1)
                call EMPM_densityN(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Jb_Range=[Jlo, Jhi])
                tr = tr + elapsed_s(t1)
                call timer_start(t1)
                do c = 1, nch
                    if (.not. Need(c)) cycle
                    do iq = Gs(c)%Col0(ig) + 1, Gs(c)%Col0(ig+1)
                        q = Gs(c)%Blk(iq)
                        call row_AD3m(Pairs,Ph,C1,Slow,S2,B(c),Rho2,a,B(c)%La(B(c)%Row(q)),Gs(c)%Db(:,iq),Gs(c)%Ab(:,iq))
                    end do
                end do
                tw = tw + elapsed_s(t1)
            end do
            !$omp end parallel do
            T_Rho2 = T_Rho2 + tr
            T_Row3 = T_Row3 + tw
            do c = 1, nch
                if (.not. allocated(Gs(c)%Blk)) cycle
                call absorb(c)
                deallocate(Gs(c)%Blk, Gs(c)%Col0, Gs(c)%Db)
                if (allocated(Gs(c)%Ab)) deallocate(Gs(c)%Ab)
            end do
            deallocate(Ablk)
        end do
        do c = 1, nch
            call finish(c)
        end do
        write(*,'(a)') TAB // 'Profile of the densities (thread-s) ... W6 ' // cfmt_f(T_dW,0,1) // ' , T ' // cfmt_f(T_dT,0,1) // &
                       ' , dgemm ' // cfmt_f(T_dG,0,1) // ' (' // cfmt_f(real(F_dG,dp)/max(T_dG,1.0e-9_dp)/1.0e9_dp,0,2) // ' GFlop/s per thread) ... rows: dgemm ' // cfmt_f(T_rG,0,1) // ' , exchange ' // cfmt_f(T_rX,0,1) // ' , scatter ' // cfmt_f(T_rS,0,1)
        write(*,'(a)') TAB // 'Seconds ... shared greedy A & D rows of ' // istr(nch) // ' channels: densities ' // cfmt_f(T_Rho2,0,1) // &
                       ' (summed over threads) , rows ' // cfmt_f(T_Row3,0,1) // ' (summed over threads) , A D^T ' // cfmt_f(T_H3,0,1) // &
                       ' , I/O ' // cfmt_f(T_io,0,1) // ' , a visited ' // istr(g) // ' of ' // istr(ng)

    contains

        ! The block of channel c tested against its selected rows (Schur complement & pivoted Cholesky), new rows added to
            ! Dii, H & D (memory or file) ...
        subroutine absorb(c)
            integer, intent(in) :: c
            real(dp), allocatable :: Dbb(:,:), Cb(:,:), X(:,:), Lf(:,:), Tmp(:,:), Bf(:,:)
            integer, allocatable :: Piv(:)
            integer :: nb, r, Rk, i, j, Info, Tg, q0, q1
            integer(i8) :: t1, t2
            associate(Gc => Gs(c), B3 => B(c), Dii => Res(c)%Dii, H => Res(c)%H)
                nb = Gc%Cnt
                r = Gc%r
                Tg = Target(c)
                Dbb = Gc%Db(B3%Row(Gc%Blk),:)
                if (r > 0) then
                    Cb = Gc%Db(B3%Row(Gc%Sel(1:r)),:)
                    X = Cb
                    Lf = Dii(1:r,1:r)
                    call dpotrf('L',r,Lf,r,Info)
                    if (Info /= 0) call fatal('EMPM_ind_greedy_multi ... D of the selected rows not positive definite (dpotrf ' // istr(Info) // ')')
                    call dpotrs('L',r,nb,Lf,r,X,r,Info)
                    call dgemm('T','N',nb,nb,r,-1.0_dp,Cb,r,X,r,1.0_dp,Dbb,nb)
                    deallocate(X, Lf)
                end if
                call chol_piv(nb,Dbb,Tol,Piv,Rk)
                Rk = min(Rk, Tg - r)
                if (Rk <= 0) return
                if (r > 0) then
                    Dii(1:r,r+1:r+Rk) = Cb(:,Piv(1:Rk))
                    Dii(r+1:r+Rk,1:r) = transpose(Cb(:,Piv(1:Rk)))
                end if
                do j = 1, Rk
                    do i = 1, Rk
                        Dii(r+i,r+j) = Gc%Db(B3%Row(Gc%Blk(Piv(i))),Piv(j))
                    end do
                end do
                call timer_start(t1)
                allocate(Tmp(Gc%Nf,Rk))
                !$omp parallel do
                do j = 1, Rk
                    Tmp(:,j) = Gc%Ab(:,Piv(j))
                end do
                !$omp end parallel do
                deallocate(Gc%Ab)
                if (Disk(c)) then
                    if (r > 0) then
                        allocate(Bf(Gc%Nf,min(r,Gc%Nblk)))
                        do q0 = 1, r, Gc%Nblk
                            q1 = min(r, q0 + Gc%Nblk - 1)
                            call timer_start(t2)
                            read(Gc%uD, pos=1_i8 + 8_i8*int(Gc%Nf,i8)*int(q0-1,i8)) Bf(:,1:q1-q0+1)
                            T_io = T_io + elapsed_s(t2)
                            call dgemm('T','N',Rk,q1-q0+1,Gc%Nf,1.0_dp,Tmp,Gc%Nf,Bf,Gc%Nf,0.0_dp,H(r+1,q0),Tg)
                        end do
                        deallocate(Bf)
                    end if
                    allocate(Bf(Gc%Nf,Rk))
                    !$omp parallel do
                    do j = 1, Rk
                        Bf(:,j) = Gc%Db(:,Piv(j))
                    end do
                    !$omp end parallel do
                    do j = 1, Rk
                        do i = 1, Gc%Ns
                            Gc%Ds(i,r+j) = Bf(B3%Ks(i),j)
                        end do
                    end do
                    call dgemm('T','N',Rk,Rk,Gc%Nf,1.0_dp,Tmp,Gc%Nf,Bf,Gc%Nf,0.0_dp,H(r+1,r+1),Tg)
                    call timer_start(t2)
                    write(Gc%uD, pos=1_i8 + 8_i8*int(Gc%Nf,i8)*int(r,i8)) Bf
                    T_io = T_io + elapsed_s(t2)
                    deallocate(Bf)
                else
                    !$omp parallel do
                    do j = 1, Rk
                        Res(c)%Dr(:,r+j) = Gc%Db(:,Piv(j))
                    end do
                    !$omp end parallel do
                    call dgemm('T','N',Rk,r+Rk,Gc%Nf,1.0_dp,Tmp,Gc%Nf,Res(c)%Dr,Gc%Nf,0.0_dp,H(r+1,1),Tg)
                end if
                deallocate(Tmp)
                if (r > 0) H(1:r,r+1:r+Rk) = transpose(H(r+1:r+Rk,1:r))
                T_H3 = T_H3 + elapsed_s(t1)
                Gc%Sel(r+1:r+Rk) = Gc%Blk(Piv(1:Rk))
                Gc%r = r + Rk
            end associate
        end subroutine absorb

        ! Ascending order of the selected rows of channel c ... H, Dii permuted, D (memory: columns permuted in place;
            ! file: Ord & Dsp returned) ...
        subroutine finish(c)
            integer, intent(in) :: c
            real(dp), allocatable :: Key(:), Col(:), Tmp(:,:)
            integer, allocatable :: Ord(:)
            logical, allocatable :: Done(:)
            integer :: r, i, j, k
            r = Gs(c)%r
            allocate(Key(r), Ord(r))
            Key = real(Gs(c)%Sel(1:r),dp)
            call sortperm(Key,Ord)
            Res(c)%Ind = Gs(c)%Sel(Ord)
            Res(c)%Dii = Res(c)%Dii(Ord,Ord)
            Res(c)%H = Res(c)%H(Ord,Ord)
            if (Disk(c)) then
                close(Gs(c)%uD)
                Res(c)%Ord = Ord
                allocate(Res(c)%Dsp(Gs(c)%Ns,r))
                do j = 1, r
                    Res(c)%Dsp(:,j) = Gs(c)%Ds(:,Ord(j))
                end do
                deallocate(Gs(c)%Ds)
                return
            end if
            allocate(Done(r), Col(Gs(c)%Nf))
            Done = .false.
            associate(Dr => Res(c)%Dr)
                do i = 1, r
                    if (Done(i) .or. Ord(i) == i) then
                        Done(i) = .true.
                        cycle
                    end if
                    Col = Dr(:,i)
                    j = i
                    do
                        Done(j) = .true.
                        k = Ord(j)
                        if (k == i) then
                            Dr(:,j) = Col
                            exit
                        end if
                        Dr(:,j) = Dr(:,k)
                        j = k
                    end do
                end do
            end associate
            if (r < Target(c)) then
                allocate(Tmp(Gs(c)%Nf,r))
                Tmp = Res(c)%Dr(:,1:r)
                call move_alloc(Tmp,Res(c)%Dr)
            end if
        end subroutine finish

    end subroutine EMPM_ind_greedy_multi

    ! H = A D^T & the truncated block of D of an n-phonon channel, at most Mem_GB of A & D rows in memory ... rows sorted by
        ! the (n-1)-phonon state a (one density per a & block) ... H_IJ & H_JI both computed as in EMPM_eqm2_build ...
    subroutine EMPM_eqm3_build(Pairs,Ph,C1,Rlow,Slow,S2,B3,Mem_GB,H,Dtr,Dsc)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B3
        real(dp), intent(in) :: Mem_GB
        real(dp), allocatable, intent(out) :: H(:,:), Dtr(:,:), Dsc(:,:)
        real(dp), allocatable :: Ab(:,:), Db(:,:), Aj(:,:), Dj(:,:), G(:,:), Key(:)
        integer, allocatable :: Order(:)
        integer :: N, R, n_blk, Ib, Jb, i0, i1, j0, j1, i, q
        integer(i8) :: t0

        N = B3%N_tr
        allocate(Order(N), Key(N))
        do i = 1, N
            Key(i) = real(B3%Al(B3%Row(i)),dp)
        end do
        call sortperm(Key,Order)
        T_Rho2 = 0.0_dp
        T_Row3 = 0.0_dp
        T_H3 = 0.0_dp
        if (16.0_dp*real(B3%N_full,dp)*real(N,dp) <= Mem_GB*1.0e9_dp) then
            R = N
        else
            R = int(max(1.0_dp, Mem_GB*1.0e9_dp/(32.0_dp*real(B3%N_full,dp))))
        end if
        n_blk = (N + R - 1)/R
        if (n_blk > 1) write(*,'(a)') TAB // 'Row blocks = ' // istr(n_blk) // ' of ' // istr(R) // ' rows (Mem_GB = ' // &
                                      cfmt_f(Mem_GB,0,1) // ') ... H_IJ & H_JI both computed, full symmetry check'
        allocate(H(N,N), Dtr(N,N), Dsc(B3%N_spur,N))
        do Ib = 1, n_blk
            call timer_start(t0)
            i0 = (Ib - 1)*R + 1
            i1 = min(Ib*R, N)
            allocate(Ab(B3%N_full,i1-i0+1), Db(B3%N_full,i1-i0+1))
            call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Order(i0:i1),Db,Ab)
            call h_block(Ab,Db,Order(i0:i1),Order(i0:i1))
            do q = 1, i1 - i0 + 1
                Dtr(Order(i0+q-1),:) = Db(B3%Row,q)
                Dsc(:,Order(i0+q-1)) = Db(B3%Ks,q)
            end do
            do Jb = Ib + 1, n_blk
                j0 = (Jb - 1)*R + 1
                j1 = min(Jb*R, N)
                allocate(Aj(B3%N_full,j1-j0+1), Dj(B3%N_full,j1-j0+1))
                call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Order(j0:j1),Dj,Aj)
                call h_block(Ab,Dj,Order(i0:i1),Order(j0:j1))
                call h_block(Aj,Db,Order(j0:j1),Order(i0:i1))
                deallocate(Aj, Dj)
            end do
            deallocate(Ab, Db)
            if (n_blk > 1) call timer_report(t0,'row block ' // istr(Ib))
        end do
        write(*,'(a)') TAB // 'Seconds ... densities ' // cfmt_f(T_Rho2,0,1) // ' (summed over threads) , rows ' // &
                       cfmt_f(T_Row3,0,1) // ' (summed over threads) , A D^T ' // cfmt_f(T_H3,0,1)

    contains

        subroutine h_block(Ab,Db,Ri,Rj)
            real(dp), intent(in) :: Ab(:,:), Db(:,:)
            integer, intent(in) :: Ri(:), Rj(:)
            integer(i8) :: t1
            call timer_start(t1)
            allocate(G(size(Ri),size(Rj)))
            call dgemm('T','N',size(Ri),size(Rj),B3%N_full,1.0_dp,Ab,B3%N_full,Db,B3%N_full,0.0_dp,G,size(Ri))
            H(Ri,Rj) = G
            deallocate(G)
            T_H3 = T_H3 + elapsed_s(t1)
        end subroutine h_block

    end subroutine EMPM_eqm3_build

    ! Ind_First, pass 1 ... D of the truncated rows only (no A), Dtr(:,i) = D row i on the truncated columns, rows computed in
        ! blocks of at most Mem_GB (full-basis rows, consecutive rows = same (n-1)-phonon state a, as the basis is a outer) ...
    subroutine EMPM_D3_build(Pairs,Ph,C1,Rlow,Slow,S2,B3,Mem_GB,Dtr)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B3
        real(dp), intent(in) :: Mem_GB
        real(dp), allocatable, intent(out) :: Dtr(:,:)
        real(dp), allocatable :: Db(:,:)
        integer :: N, R, n_blk, Ib, i0, i1, i, q

        N = B3%N_tr
        T_Rho2 = 0.0_dp
        T_Row3 = 0.0_dp
        R = int(min(real(N,dp), max(1.0_dp, Mem_GB*1.0e9_dp/(8.0_dp*real(B3%N_full,dp)))))
        n_blk = (N + R - 1)/R
        write(*,'(a)') TAB // 'Independent rows first ... D of all truncated rows ' // &
                       cfmt_f(8.0_dp*real(N,dp)**2/1.0e9_dp,0,3) // ' GB , row blocks = ' // istr(n_blk) // ' of ' // istr(R)
        allocate(Dtr(N,N))
        do Ib = 1, n_blk
            i0 = (Ib - 1)*R + 1
            i1 = min(Ib*R, N)
            allocate(Db(B3%N_full,i1-i0+1))
            call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,[(i, i = i0, i1)],Db)
            !$omp parallel do
            do q = 1, i1 - i0 + 1
                Dtr(:,i0+q-1) = Db(B3%Row,q)
            end do
            !$omp end parallel do
            deallocate(Db)
        end do
        write(*,'(a)') TAB // 'Seconds ... D: densities ' // cfmt_f(T_Rho2,0,1) // ' (summed over threads) , rows ' // &
                       cfmt_f(T_Row3,0,1) // ' (summed over threads)'
    end subroutine EMPM_D3_build

    ! Ind_First, pass 2 ... A & D rows of the independent rows Ind (ascending = grouped by a) only, H = A D^T on them
        ! (r x r, H_IJ & H_JI both computed) & Dr = their D rows on the full basis (X of the exported states) ... all A rows
        ! at once if they fit in Mem_GB, else A in blocks after the D rows (D rows of the blocks computed again) ...
    subroutine EMPM_H3_ind(Pairs,Ph,C1,Rlow,Slow,S2,B3,Ind,Mem_GB,H,Dr)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B3
        integer, intent(in) :: Ind(:)
        real(dp), intent(in) :: Mem_GB
        real(dp), allocatable, intent(out) :: H(:,:), Dr(:,:)
        real(dp), allocatable :: Ab(:,:), Db(:,:)
        integer :: r, R_blk, n_blk, Ib, i0, i1
        integer(i8) :: t1

        r = size(Ind)
        T_Rho2 = 0.0_dp
        T_Row3 = 0.0_dp
        T_H3 = 0.0_dp
        allocate(H(r,r), Dr(B3%N_full,r))
        if (r == 0) return
        if (8.0_dp*real(B3%N_full,dp)*real(r,dp) <= Mem_GB*1.0e9_dp) then
            allocate(Ab(B3%N_full,r))
            call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Ind,Dr,Ab)
            call timer_start(t1)
            call dgemm('T','N',r,r,B3%N_full,1.0_dp,Ab,B3%N_full,Dr,B3%N_full,0.0_dp,H,r)
            T_H3 = T_H3 + elapsed_s(t1)
            deallocate(Ab)
        else
            call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Ind,Dr)
            R_blk = int(max(1.0_dp, Mem_GB*1.0e9_dp/(16.0_dp*real(B3%N_full,dp))))
            n_blk = (r + R_blk - 1)/R_blk
            write(*,'(a)') TAB // 'A row blocks = ' // istr(n_blk) // ' of ' // istr(R_blk) // ' rows (Mem_GB = ' // &
                           cfmt_f(Mem_GB,0,1) // ')'
            do Ib = 1, n_blk
                i0 = (Ib - 1)*R_blk + 1
                i1 = min(Ib*R_blk, r)
                allocate(Ab(B3%N_full,i1-i0+1), Db(B3%N_full,i1-i0+1))
                call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Ind(i0:i1),Db,Ab)
                call timer_start(t1)
                call dgemm('T','N',i1-i0+1,r,B3%N_full,1.0_dp,Ab,B3%N_full,Dr,B3%N_full,0.0_dp,H(i0,1),r)
                T_H3 = T_H3 + elapsed_s(t1)
                deallocate(Ab, Db)
            end do
        end if
        write(*,'(a)') TAB // 'Seconds ... A & D of the independent rows: densities ' // cfmt_f(T_Rho2,0,1) // &
                       ' (summed over threads) , rows ' // cfmt_f(T_Row3,0,1) // ' (summed over threads) , A D^T ' // &
                       cfmt_f(T_H3,0,1)
    end subroutine EMPM_H3_ind

    ! Rows Rows(q) (truncated indices) of D & optionally A over the full n-phonon basis, in parallel over the (n-1)-phonon
        ! states a (one density per group of consecutive rows with the same a) ...
    subroutine EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S2,B3,Rows,Drow,Arow)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S2
        type(empm_basis2_t), intent(in) :: B3
        integer, intent(in) :: Rows(:)
        real(dp), intent(out) :: Drow(:,:)
        real(dp), intent(out), optional :: Arow(:,:)
        logical :: Want_A
        type(phonon_block_t) :: Rho2
        integer, allocatable :: G0(:)
        integer :: ng, iq, ig, a
        integer(i8) :: t1
        real(dp) :: tr, tw

        Want_A = present(Arow)
        ! Groups of consecutive rows with the same a ...
        allocate(G0(size(Rows)+1))
        ng = 0
        do iq = 1, size(Rows)
            if (iq == 1) then
                ng = ng + 1
                G0(ng) = iq
            else if (B3%Al(B3%Row(Rows(iq))) /= B3%Al(B3%Row(Rows(iq-1)))) then
                ng = ng + 1
                G0(ng) = iq
            end if
        end do
        G0(ng+1) = size(Rows) + 1
        tr = 0.0_dp
        tw = 0.0_dp
        !$omp parallel do schedule(dynamic) private(a,Rho2,iq,t1) reduction(+:tr,tw)
        do ig = 1, ng
            a = B3%Al(B3%Row(Rows(G0(ig))))
            call timer_start(t1)
            ! only the channels b of the columns (b c; J): |J - J_c| <= J_b <= J + J_c, J_c <= the phonon J_max ...
            call EMPM_densityN(Pairs,Ph,C1,Rlow,Slow,S2,a,Rho2,Jb_Range=[max(0,B3%J-Ph%J_max), B3%J+Ph%J_max])
            tr = tr + elapsed_s(t1)
            call timer_start(t1)
            do iq = G0(ig), G0(ig+1) - 1
                call row_AD3(Rho2,iq,a,B3%La(B3%Row(Rows(iq))))
            end do
            tw = tw + elapsed_s(t1)
        end do
        !$omp end parallel do
        T_Rho2 = T_Rho2 + tr
        T_Row3 = T_Row3 + tw

    contains

        subroutine row_AD3(Rho2,q,a,la)
            type(phonon_block_t), intent(in) :: Rho2
            integer, intent(in) :: q, a, la
            if (Want_A) then
                call row_AD3m(Pairs,Ph,C1,Slow,S2,B3,Rho2,a,la,Drow(:,q),Arow(:,q))
            else
                call row_AD3m(Pairs,Ph,C1,Slow,S2,B3,Rho2,a,la,Drow(:,q))
            end if
        end subroutine row_AD3

    end subroutine EMPM_AD3_rows

    ! Row (a, la) of D (Dc) & optionally A (Ac) over the full n-phonon basis of B3, from the densities Rho2 of a ...
    subroutine row_AD3m(Pairs,Ph,C1,Slow,S2,B3,Rho2,a,la,Dc,Ac)
            type(empm_pairs_t), intent(in) :: Pairs
            type(empm_phonons_t), intent(in) :: Ph
            type(empm_cache1_t), intent(in) :: C1
            type(empm_states2_t), intent(in) :: Slow, S2
            type(empm_basis2_t), intent(in) :: B3
            real(dp), intent(out) :: Dc(:)
            real(dp), intent(out), optional :: Ac(:)
            logical :: Want_A
            type(phonon_block_t), intent(in) :: Rho2
            integer, intent(in) :: a, la
            real(dp), allocatable :: GA(:,:), GD(:,:), F6(:,:,:), Xg(:,:), Wx(:,:)
            integer, allocatable :: Gl(:), Kg(:)
            integer :: ng, ig
            logical :: Nzw
            integer :: Ja, Pa, ia, Jl, J, Jb, Pb, Jc, Pc, nb, nc, lo, hi, Sig, ib, ic, k, b, c, Ka, Pp, gm, kb, ka2, lp, Jg, &
                       np, ns, is
            real(dp) :: wgtA, wgtD, w, tg, tx, ts
            integer(i8) :: tq

            Want_A = present(Ac)
            tg = 0.0_dp
            tx = 0.0_dp
            ts = 0.0_dp
            Ja = S2%J(a)
            Pa = S2%P(a)
            ia = a - S2%First(Ja,Pa) + 1
            Jl = Ph%J(la)
            J = B3%J
            if (Want_A) Ac(:) = 0.0_dp
            Dc(:) = 0.0_dp
            ! F6(Jg,Jb,Jc) = {J_l Jg Jb; Jc J Ja} for the exchange term ...
            allocate(F6(0:Slow%J_max,0:S2%J_max,0:Ph%J_max))
            do Jc = 0, Ph%J_max
                do Jb = 0, S2%J_max
                    do Jg = 0, Slow%J_max
                        F6(Jg,Jb,Jc) = f6j(2*Jl,2*Jg,2*Jb,2*Jc,2*J,2*Ja)
                    end do
                end do
            end do
            do Jb = 0, S2%J_max
                do Pb = 1, 2
                    nb = S2%N_ch(Jb,Pb)
                    if (nb == 0) cycle
                    if (abs(J - Jb) > Ph%J_max) cycle
                    ! rows X^b(g la, :) of the basis states (g la) of the channel of b, gathered once (Xg(:,ig) contiguous) ...
                    call timer_start(tq)
                    associate(Bb => S2%Bs(Jb,Pb), Xb => S2%X(Jb,Pb)%v)
                        ng = count(Bb%K(:,la) /= 0)
                        allocate(Gl(ng), Kg(ng), Xg(nb,ng))
                        ng = 0
                        do gm = 1, Slow%N
                            if (Bb%K(gm,la) == 0) cycle
                            ng = ng + 1
                            Gl(ng) = gm
                            Kg(ng) = Bb%K(gm,la)
                        end do
                        if (allocated(S2%Xt)) then
                            do ig = 1, ng
                                Xg(:,ig) = S2%Xt(Jb,Pb)%v(:,Kg(ig))
                            end do
                        else
                            do ib = 1, nb
                                do ig = 1, ng
                                    Xg(ib,ig) = Xb(Kg(ig),ib)
                                end do
                            end do
                        end if
                    end associate
                    tx = tx + elapsed_s(tq)
                    do Jc = abs(J - Jb), min(J + Jb, Ph%J_max)
                        do Pc = 1, 2
                            nc = Ph%N_ch(Jc,Pc)
                            if (nc == 0) cycle
                            if (par_sign(Pb)*par_sign(Pc) /= par_sign(B3%P)) cycle
                            allocate(GA(nb,nc), GD(nb,nc))
                            GA = 0.0_dp
                            GD = 0.0_dp
                            ! sum_s w_s rho_ab(:,s)^T X_l(:,s,l') as in the 2-phonon rows ...
                            call timer_start(tq)
                            lo = max(abs(Ja - Jb), abs(Jl - Jc))
                            hi = min(Ja + Jb, Jl + Jc, Pairs%L_max)
                            Pp = Rho2%Pi(Jb,Pb)
                            do Sig = lo, hi
                                Ka = Pairs%Na(Sig,Pp)
                                if (Ka == 0) cycle
                                wgtD = -phase(J + Jc + Ja) * f6j(2*Sig,2*Jc,2*Jl,2*J,2*Ja,2*Jb)
                                if (wgtD /= 0.0_dp) call dgemm('T','N',nb,nc,Ka,wgtD,Rho2%B(Jb,Pb)%M(Sig)%v,Ka, &
                                                               C1%Rt(la)%B(Jc,Pc)%M(Sig)%v,Ka,1.0_dp,GD,nb)
                                if (.not. Want_A) cycle
                                wgtA = phase(J + Ja) * f6j(2*J,2*Jc,2*Jb,2*Sig,2*Ja,2*Jl)
                                if (wgtA /= 0.0_dp) call dgemm('T','N',nb,nc,Ka,wgtA,Rho2%B(Jb,Pb)%M(Sig)%v,Ka, &
                                                               C1%V(la)%B(Jc,Pc)%M(Sig)%v,Ka,1.0_dp,GA,nb)
                            end do
                            tg = tg + elapsed_s(tq)
                            call timer_start(tq)
                            ! Metric exchange sum_g (-1)^(j_l'+j_l+Ja+Jb) {j_l j_g Jb; j_l' J Ja} X^a(g l') X^b(g l) ...
                            ! W(ig,ic) = (-1)^(...) {..} X^a(g l', a) & GD += Xg W ...
                            if (ng > 0) then
                                allocate(Wx(ng,nc))
                                Wx = 0.0_dp
                                Nzw = .false.
                                associate(Ba => S2%Bs(Ja,Pa))
                                    do ic = 1, nc
                                        lp = Ph%First(Jc,Pc) + ic - 1
                                        do ig = 1, ng
                                            gm = Gl(ig)
                                            ka2 = Ba%K(gm,lp)
                                            if (ka2 == 0) cycle
                                            w = F6(Slow%J(gm),Jb,Jc)
                                            if (w == 0.0_dp) cycle
                                            Wx(ig,ic) = phase(Jc + Jl + Ja + Jb) * w * S2%X(Ja,Pa)%v(ka2,ia)
                                            Nzw = .true.
                                        end do
                                    end do
                                end associate
                                if (Nzw) call dgemm('N','N',nb,nc,ng,1.0_dp,Xg,nb,Wx,ng,1.0_dp,GD,nb)
                                deallocate(Wx)
                            end if
                            tx = tx + elapsed_s(tq)
                            call timer_start(tq)
                            do ic = 1, nc
                                c = Ph%First(Jc,Pc) + ic - 1
                                do ib = 1, nb
                                    b = S2%First(Jb,Pb) + ib - 1
                                    k = B3%K(b,c)
                                    Dc(k) = Dc(k) + GD(ib,ic)
                                    if (Want_A) Ac(k) = Ac(k) + GA(ib,ic)
                                end do
                            end do
                            deallocate(GA, GD)
                            ts = ts + elapsed_s(tq)
                        end do
                    end do
                    deallocate(Gl, Kg, Xg)
                end do
            end do
            !$omp atomic
            T_rG = T_rG + tg
            !$omp atomic
            T_rX = T_rX + tx
            !$omp atomic
            T_rS = T_rS + ts
            k = B3%K(a,la)
            Dc(k) = Dc(k) + 1.0_dp
            if (.not. Want_A) return
            Ac(k) = Ac(k) + S2%E(a) + Ph%E(la)
            ! d_ll' H^(n-1)_ab ... couplings of the physical & spurious (n-1)-phonon states of the channel of a ...
            np = S2%N_phys(Ja,Pa)
            ns = S2%N_ch(Ja,Pa) - np
            if (ns > 0) then
                if (ia <= np) then
                    do is = 1, ns
                        k = B3%K(S2%First(Ja,Pa)+np+is-1,la)
                        if (k > 0) Ac(k) = Ac(k) + S2%Rc(Ja,Pa)%v(ia,is)
                    end do
                else
                    do ib = 1, np
                        k = B3%K(S2%First(Ja,Pa)+ib-1,la)
                        if (k > 0) Ac(k) = Ac(k) + S2%Rc(Ja,Pa)%v(ib,ia-np)
                    end do
                end if
            end if
            ! d_ab R_ll' ... couplings of the 1- phonons to the spurious phonon ...
            if (Ph%i_cm > 0) then
                do c = 1, Ph%N
                    if (c == la) cycle
                    if (la == Ph%i_cm .or. c == Ph%i_cm) then
                        k = B3%K(a,c)
                        if (k > 0) Ac(k) = Ac(k) + merge(Ph%R(c), Ph%R(la), la == Ph%i_cm)
                    end if
                end do
            end if
    end subroutine row_AD3m

    ! One channel of a set of n-phonon states (raw stream) ... States_Format (< 0), N_full, n, N_ind, N_phys, Al(N_full)
        ! (the (n-1)-phonon state of each basis state), La(N_full) (the phonon), E(n), Kc(N_ind) (full-basis indices of the
        ! independent rows, ascending), C(N_ind,n) on them (zero on all other rows), X(N_full,n) = sqrt(2J+1)
        ! D(independent rows)^T C & R(N_phys,n-N_phys) ... the N_phys physical states first, then the spurious ones with
        ! their couplings R ... the input of the (n+1)-phonon step ... (format -2 without N_phys & R, all states physical,
        ! & the first format, without the marker, N_full, n, Al, La, E, C(N_full,n), X ... still read) ...
    subroutine EMPM_states_write(File,Bs,E,Kc,Cc,X,N_phys,R)
        character(len=*), intent(in) :: File
        type(empm_basis2_t), intent(in) :: Bs
        integer, intent(in) :: Kc(:), N_phys
        real(dp), intent(in) :: E(:), Cc(:,:), X(:,:), R(:,:)
        integer :: u
        open(newunit=u, file=File, access='stream', form='unformatted', status='replace', action='write')
        write(u) States_Format, Bs%N_full, size(E), size(Kc), N_phys
        write(u) Bs%Al, Bs%La
        write(u) E
        write(u) Kc
        write(u) Cc(:,1:size(E))
        write(u) X
        write(u) R(1:N_phys,1:size(E)-N_phys)
        close(u)
    end subroutine EMPM_states_write

    ! As EMPM_states_write (same file), with X = Scale * Dr C computed & written in blocks of columns (no full X in
        ! memory) ... Dr the full-basis D rows of the independent rows Kc ...
    subroutine EMPM_states_write_D(File,Bs,E,Kc,Cc,Dr,Scale,N_phys,R)
        character(len=*), intent(in) :: File
        type(empm_basis2_t), intent(in) :: Bs
        integer, intent(in) :: Kc(:), N_phys
        real(dp), intent(in) :: E(:), Cc(:,:), Dr(:,:), Scale, R(:,:)
        real(dp), allocatable :: Xb(:,:)
        integer :: u, n, nr, j0, j1, Nc
        n = size(E)
        nr = size(Dr,2)
        open(newunit=u, file=File, access='stream', form='unformatted', status='replace', action='write')
        write(u) States_Format, Bs%N_full, n, size(Kc), N_phys
        write(u) Bs%Al, Bs%La
        write(u) E
        write(u) Kc
        write(u) Cc(:,1:n)
        ! column blocks of about 2 GB ...
        Nc = max(1, min(n, int(2.0e9_dp/(8.0_dp*real(max(1,Bs%N_full),dp)))))
        allocate(Xb(Bs%N_full,Nc))
        do j0 = 1, n, Nc
            j1 = min(n, j0 + Nc - 1)
            if (nr > 0) call dgemm('N','N',Bs%N_full,j1-j0+1,nr,Scale,Dr,Bs%N_full,Cc(:,j0:j1),nr,0.0_dp,Xb,Bs%N_full)
            write(u) Xb(:,1:j1-j0+1)
        end do
        write(u) R(1:N_phys,1:n-N_phys)
        close(u)
    end subroutine EMPM_states_write_D

    ! As EMPM_states_write_D, the D rows read from Dr_File of EMPM_ind_greedy_H (Nf x nr, column Ord(i) = independent row i)
        ! ... X computed & written in blocks of basis rows (memory about Mem_GB / 4), Dr_File deleted at the end ...
    subroutine EMPM_states_write_Df(File,Bs,E,Kc,Cc,Dr_File,Ord,Scale,N_phys,R,Mem_GB)
        character(len=*), intent(in) :: File, Dr_File
        type(empm_basis2_t), intent(in) :: Bs
        integer, intent(in) :: Kc(:), N_phys, Ord(:)
        real(dp), intent(in) :: E(:), Cc(:,:), Scale, R(:,:), Mem_GB
        real(dp), allocatable :: Xb(:,:), Db(:,:), Cp(:,:)
        integer :: u, uD, n, nr, i0, i1, m, j, rb, Nf
        integer(i8) :: pX, t1
        real(dp) :: T_io
        n = size(E)
        nr = size(Ord)
        Nf = Bs%N_full
        open(newunit=u, file=File, access='stream', form='unformatted', status='replace', action='readwrite')
        write(u) States_Format, Nf, n, size(Kc), N_phys
        write(u) Bs%Al, Bs%La
        write(u) E
        write(u) Kc
        write(u) Cc(:,1:n)
        inquire(u, pos=pX)
        ! C with its rows in the order of the file columns ...
        allocate(Cp(nr,n))
        do j = 1, nr
            Cp(Ord(j),:) = Cc(j,1:n)
        end do
        rb = max(1, min(Nf, int(0.25_dp*Mem_GB*1.0e9_dp/(8.0_dp*real(nr + n,dp)))))
        allocate(Db(rb,nr), Xb(rb,n))
        open(newunit=uD, file=Dr_File, access='stream', form='unformatted', status='old', action='read')
        T_io = 0.0_dp
        do i0 = 1, Nf, rb
            i1 = min(Nf, i0 + rb - 1)
            m = i1 - i0 + 1
            call timer_start(t1)
            do j = 1, nr
                read(uD, pos=1_i8 + 8_i8*(int(Nf,i8)*int(j-1,i8) + int(i0-1,i8))) Db(1:m,j)
            end do
            T_io = T_io + elapsed_s(t1)
            if (nr > 0) call dgemm('N','N',m,n,nr,Scale,Db,rb,Cp,nr,0.0_dp,Xb,rb)
            call timer_start(t1)
            do j = 1, n
                write(u, pos=pX + 8_i8*(int(Nf,i8)*int(j-1,i8) + int(i0-1,i8))) Xb(1:m,j)
            end do
            T_io = T_io + elapsed_s(t1)
        end do
        close(uD, status='delete')
        write(u, pos=pX + 8_i8*int(Nf,i8)*int(n,i8)) R(1:N_phys,1:n-N_phys)
        close(u)
        write(*,'(a)') TAB // 'X from the D rows on disk ... row blocks of ' // istr(rb) // ' , I/O ' // cfmt_f(T_io,0,1) // ' s'
    end subroutine EMPM_states_write_Df

    ! Set of n-phonon states of all channels J = 0 .. J_max from the files Prefix // '<JP>.bin' of EMPM_states_write (J outer,
        ! P inner, as numbered when written) ... N_low (n-1)-phonon & N1 1-phonon states for the index K(a,lambda) of the bases ...
    subroutine EMPM_states_read(Prefix,J_max,N_low,N1,S)
        character(len=*), intent(in) :: Prefix
        integer, intent(in) :: J_max, N_low, N1
        type(empm_states2_t), intent(out) :: S
        type(rvec_t), allocatable :: Ech(:,:)
        integer :: J, P, u, N_full, n, i, k, Fmt, N_ind, Np
        real(dp), allocatable :: Cf(:,:)
        logical, allocatable :: Nz(:)
        logical :: Ex
        character(len=:), allocatable :: File

        S%J_max = 0
        allocate(S%N_ch(0:J_max,2), S%First(0:J_max,2), S%Bs(0:J_max,2), S%Kc(0:J_max,2), S%Cc(0:J_max,2), S%X(0:J_max,2), &
                 Ech(0:J_max,2), S%N_phys(0:J_max,2), S%Rc(0:J_max,2))
        S%N_ch = 0
        S%N_phys = 0
        do J = 0, J_max
            do P = 1, 2
                File = Prefix // istr(J) // merge('+','-',P == 1) // '.bin'
                inquire(file=File, exist=Ex)
                S%Bs(J,P)%J = J
                S%Bs(J,P)%P = P
                if (.not. Ex) then
                    allocate(S%Bs(J,P)%Al(0), S%Bs(J,P)%La(0), S%Kc(J,P)%v(0), S%Cc(J,P)%v(0,0), S%X(J,P)%v(0,0), Ech(J,P)%v(0), &
                             S%Rc(J,P)%v(0,0))
                    allocate(S%Bs(J,P)%K(N_low,N1))
                    S%Bs(J,P)%K = 0
                    cycle
                end if
                open(newunit=u, file=File, access='stream', form='unformatted', status='old', action='read')
                read(u) Fmt
                if (Fmt == States_Format) then
                    read(u) N_full, n, N_ind, Np
                else if (Fmt == -2) then
                    read(u) N_full, n, N_ind
                    Np = n
                else
                    N_full = Fmt
                    read(u) n
                    Np = n
                end if
                associate(Bs => S%Bs(J,P))
                    Bs%N_full = N_full
                    allocate(Bs%Al(N_full), Bs%La(N_full), Bs%K(N_low,N1), Ech(J,P)%v(n), S%X(J,P)%v(N_full,n))
                    read(u) Bs%Al, Bs%La
                    read(u) Ech(J,P)%v
                    if (Fmt < 0) then
                        allocate(S%Kc(J,P)%v(N_ind), S%Cc(J,P)%v(N_ind,n))
                        read(u) S%Kc(J,P)%v
                        read(u) S%Cc(J,P)%v
                    else
                        ! First format ... C on the full basis, kept on the rows with a non-zero amplitude ...
                        allocate(Cf(N_full,n))
                        read(u) Cf
                        Nz = any(Cf /= 0.0_dp, dim=2)
                        S%Kc(J,P)%v = pack([(k, k = 1, N_full)], Nz)
                        S%Cc(J,P)%v = Cf(S%Kc(J,P)%v,:)
                        deallocate(Cf)
                    end if
                    read(u) S%X(J,P)%v
                    allocate(S%Rc(J,P)%v(Np,n-Np))
                    if (Fmt == States_Format) read(u) S%Rc(J,P)%v
                    Bs%K = 0
                    do k = 1, N_full
                        Bs%K(Bs%Al(k),Bs%La(k)) = k
                    end do
                end associate
                close(u)
                S%N_ch(J,P) = n
                S%N_phys(J,P) = Np
                if (n > 0) S%J_max = J
            end do
        end do
        ! X transposed ...
        allocate(S%Xt(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                associate(Xv => S%X(J,P)%v)
                    allocate(S%Xt(J,P)%v(size(Xv,2),size(Xv,1)))
                    !$omp parallel do private(k)
                    do i = 1, size(Xv,2)
                        do k = 1, size(Xv,1)
                            S%Xt(J,P)%v(i,k) = Xv(k,i)
                        end do
                    end do
                    !$omp end parallel do
                end associate
            end do
        end do
        S%N = sum(S%N_ch)
        allocate(S%J(S%N), S%P(S%N), S%E(S%N), S%Spur(S%N))
        i = 0
        do J = 0, J_max
            do P = 1, 2
                S%First(J,P) = i + 1
                n = S%N_ch(J,P)
                S%J(i+1:i+n) = J
                S%P(i+1:i+n) = P
                S%E(i+1:i+n) = Ech(J,P)%v
                S%Spur(i+1:i+n) = [(k > S%N_phys(J,P), k = 1, n)]
                i = i + n
            end do
        end do
    end subroutine EMPM_states_read

    ! Wall-clock seconds since t0 ...
    real(dp) function elapsed_s(t0)
        integer(i8), intent(in) :: t0
        integer(i8) :: c, Rate
        call system_clock(c,Rate)
        elapsed_s = real(c - t0,dp)/real(Rate,dp)
    end function elapsed_s

end module nj_empm_eqm3
