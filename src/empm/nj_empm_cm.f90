! Exact centre-of-mass condition b_cm psi = 0 in the full EMPM Hamiltonian (EMPM_FULL with CM_Exact) ...
    ! b_cm,mu = A^(-1/2) sum_i b_mu(i), b_mu = sqrt(2) (r_mu/b)_(N -> N-1) the HO lowering part of r in the LHO basis,
    ! transformed to the HF basis (it has ph, hp, pp & hh parts there) ... L = sum_i b(i) is a rank-1 one-body operator, its
    ! reduced matrix elements (Edmonds, <J'M'|T_mu|JM> = <J M 1 mu|J' M'> <J'||T||J>/sqrt(2J'+1)) between the EMPM states
    ! of the complete sets of all levels n = 0 .. Full_N (Ortho = .false.) ...
        ! n -> n+1 (ph part)   <b||L||a>     = 1/sqrt3 sum_(lambda 1-) bL(lambda) X^b(a lambda)        (X of level n+1)
        ! n+1 -> n (hp part)   <a'||L||b>    = (-1)^(J'-J)/sqrt3 sum_(lambda 1-) bR(lambda) X^b(a' lambda)
        ! n -> n   (pp & hh)   <al||L||a>    = sum_(pairs s) wD(s) rho_(a al)(1; s) , wD(s) = <i1||L||i2>/sqrt3, by the
        !                      density recursion with the pair index contracted (EMPM_cm_dlevels, cm_dvec)
        ! 0 <-> 1              <lambda||L||0> = bL(lambda) , <0||L||lambda> = bR(lambda)
        ! bL(lambda) = sum_ph c^lambda_ph <p||L||h> & bR(lambda) = sum_ph c^lambda_ph <p||L+||h>, L+ the raising part (the
        ! hp part of L is the adjoint of the ph part of the raising operator) ...
    ! N_cm = b+_cm . b_cm = G2/A sum_J' |<J'||L||psi>|^2/(2J+1) on the channel J^P, the physical states its null space ...
    ! Phases & normalizations of the three block types (Sg_*) checked against the exact N_cm of the m-scheme FCI
    ! (tests/empm_cm_exact_check.py --fit) ...
module nj_empm_cm
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_tr1b, only : tr1b_t, Tr1b_initialize, Tr1b_transformation
    use nj_o1b, only : O1b_import
    use nj_empm_space
    use nj_empm_phonons, only : empm_phonons_t
    use nj_empm_dens, only : empm_pairs_t, phonon_block_t
    use nj_wigner, only : f6j
    use nj_empm_eqm3, only : empm_states2_t, empm_cache1_t
    implicit none
    private

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
    end interface

    ! bL, bR per phonon (non-zero only for the 1- phonons), wD per pair s, G2 = (b_mu / L_mu)^2 ...
    ! Top-level completion of N_cm (the rows of level n+1 that are not in the basis) ... only the ph part L_ph of L reaches
        ! level n+1 from level n, and sum_mu L_ph,mu+ L_ph,mu = sum_mu L_ph,mu L_ph,mu+ + C with the scalar 1-body commutator
        ! C = sum_mu [L_ph,mu+, L_ph,mu] (pp & hh parts only): c(h,h') = sum_p <p||L||h><p||L||h'>/(2j_h+1),
        ! c(p',p) = - sum_h <p'||L||h><p||L||h>/(2j_p+1) (Edmonds) ... wC per pair s = <i1||c||i2> = sqrt(2j+1) c(i1,i2) for the
        ! rank-0 density recursion, C0 = <HF|C|HF> = sum_h (2j_h+1) c(h,h) (a+_h a_h' = delta - a_h' a+_h) ...
    type, public :: empm_cm_t
        real(dp), allocatable :: bL(:), bR(:), wD(:), wC(:)
        real(dp) :: G2 = 1.0_dp, C0 = 0.0_dp
        integer :: A = 0
    end type empm_cm_t

    ! Signs of the ph, hp & pp/hh blocks by J' - J = -1, 0, 1 ... the hp block (adjoint of the ph block of the raising
        ! operator) carries (-1)^(J'-J), fixed (with the 1/sqrt3 factors & G2) by tests/empm_cm_exact_check.py --fit: N_cm of
        ! all FCI eigenstates of 4He (3 shells) reproduced to 1.3e-8 ...
    real(dp), parameter, public :: Sg_ph(-1:1) = [1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter, public :: Sg_hp(-1:1) = [-1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter, public :: Sg_d(-1:1) = [1.0_dp, 1.0_dp, 1.0_dp]

    public :: EMPM_cm_make, EMPM_cm_blocks, EMPM_cm_dims, EMPM_cm_dlevels, EMPM_cm_up, EMPM_cm_top

contains

    subroutine EMPM_cm_make(Params,Orb,Space,Pairs,Ph,Cm)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(empm_space_t), intent(in) :: Space
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cm_t), intent(out) :: Cm
        type(tr1b_t) :: TrOp, TrL, TrR
        type(o1b_t) :: C
        integer :: a, b, a_max, Na, Nb, l, p, h, ap, ah, s, k, r, i1, i2, a1, a2
        real(dp) :: b_osc

        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrL = TrOp
        TrR = TrOp
        a_max = size(TrOp%E1%p,1)
        do a = 1, a_max
            Na = 2*Orb(a)%n + Orb(a)%l
            do b = 1, a_max
                Nb = 2*Orb(b)%n + Orb(b)%l
                if (Na /= Nb - 1) then
                    TrL%E1%p(a,b) = 0.0_dp
                    TrL%E1%n(a,b) = 0.0_dp
                end if
                if (Na /= Nb + 1) then
                    TrR%E1%p(a,b) = 0.0_dp
                    TrR%E1%n(a,b) = 0.0_dp
                end if
            end do
        end do
        C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        TrL = Tr1b_transformation(Params,Orb,C,TrL)
        TrR = Tr1b_transformation(Params,Orb,C,TrR)

        ! b_mu = sqrt(2) b_osc r_mu (b_osc the inverse HO length), r_mu = sqrt(4 pi/3) r Y_1mu ...
        b_osc = sqrt(0.5_dp*(m_p + m_n)*Params%Calc%hw) / HbarC
        Cm%G2 = 2.0_dp * (4.0_dp*pi/3.0_dp) * b_osc**2
        Cm%A = Params%Calc%A

        allocate(Cm%bL(Ph%N), Cm%bR(Ph%N))
        Cm%bL = 0.0_dp
        Cm%bR = 0.0_dp
        if (Ph%J_max >= 1) then
            do l = Ph%First(1,2), Ph%First(1,2) + Ph%N_ch(1,2) - 1
                do h = 1, size(Space%pHole)
                    ah = Space%Lev(Space%pHole(h))%a
                    do p = 1, size(Space%pPart)
                        ap = Space%Lev(Space%pPart(p))%a
                        Cm%bL(l) = Cm%bL(l) + Ph%Cp(p,h,l) * TrL%E1%p(ap,ah)
                        Cm%bR(l) = Cm%bR(l) + Ph%Cp(p,h,l) * TrR%E1%p(ap,ah)
                    end do
                end do
                do h = 1, size(Space%nHole)
                    ah = Space%Lev(Space%nHole(h))%a
                    do p = 1, size(Space%nPart)
                        ap = Space%Lev(Space%nPart(p))%a
                        Cm%bL(l) = Cm%bL(l) + Ph%Cn(p,h,l) * TrL%E1%n(ap,ah)
                        Cm%bR(l) = Cm%bR(l) + Ph%Cn(p,h,l) * TrR%E1%n(ap,ah)
                    end do
                end do
            end do
        end if

        ! Pair weights ... s = Off(k) + i1 + (i2-1) Nk(k), blocks 1/2 proton/neutron particles, 3/4 holes ...
        allocate(Cm%wD(Pairs%S))
        Cm%wD = 0.0_dp
        do k = 1, 4
            do i2 = 1, Pairs%Nk(k)
                do i1 = 1, Pairs%Nk(k)
                    s = Pairs%Off(k) + i1 + (i2 - 1)*Pairs%Nk(k)
                    a1 = Space%Lev(Pairs%Pos(k)%v(i1))%a
                    a2 = Space%Lev(Pairs%Pos(k)%v(i2))%a
                    if (mod(k,2) == 1) then
                        Cm%wD(s) = TrL%E1%p(a1,a2) / sqrt(3.0_dp)
                    else
                        Cm%wD(s) = TrL%E1%n(a1,a2) / sqrt(3.0_dp)
                    end if
                end do
            end do
        end do
        ! Scalar commutator C of the ph part (top-level completion of N_cm) ...
        allocate(Cm%wC(Pairs%S))
        Cm%wC = 0.0_dp
        Cm%C0 = 0.0_dp
        do k = 1, 4
            do i2 = 1, Pairs%Nk(k)
                do i1 = 1, Pairs%Nk(k)
                    s = Pairs%Off(k) + i1 + (i2 - 1)*Pairs%Nk(k)
                    a1 = Space%Lev(Pairs%Pos(k)%v(i1))%a
                    a2 = Space%Lev(Pairs%Pos(k)%v(i2))%a
                    if (Orb(a1)%l /= Orb(a2)%l .or. Orb(a1)%j /= Orb(a2)%j) cycle
                    Cm%wC(s) = sqrt(real(Orb(a1)%j + 1,dp)) * cme(k,a1,a2)
                    if (k >= 3 .and. i1 == i2) Cm%C0 = Cm%C0 + real(Orb(a1)%j + 1,dp) * cme(k,a1,a2)
                end do
            end do
        end do
        r = count(Cm%bL /= 0.0_dp)
        write(*,'(a)') TAB // 'CM operator ... ' // istr(r) // ' 1- phonons with ph strength, |bL| = ' // &
                       cfmt_f(norm2(Cm%bL),0,6) // ' , |bR| = ' // cfmt_f(norm2(Cm%bR),0,6) // ' , |wD| = ' // &
                       cfmt_f(norm2(Cm%wD),0,6) // ' , G2 = ' // cfmt_f(Cm%G2,0,6) // ' , |wC| = ' // &
                       cfmt_f(norm2(Cm%wC),0,6) // ' , C0 = ' // cfmt_f(Cm%C0,0,6)
    contains
        ! c(a1,a2) of block k (1/2 proton/neutron particles, 3/4 holes) from the ph elements of L ...
        real(dp) function cme(k,a1,a2)
            integer, intent(in) :: k, a1, a2
            integer :: x, ax
            cme = 0.0_dp
            select case (k)
            case (1)
                do x = 1, size(Space%pHole)
                    ax = Space%Lev(Space%pHole(x))%a
                    cme = cme - TrL%E1%p(a1,ax) * TrL%E1%p(a2,ax)
                end do
            case (2)
                do x = 1, size(Space%nHole)
                    ax = Space%Lev(Space%nHole(x))%a
                    cme = cme - TrL%E1%n(a1,ax) * TrL%E1%n(a2,ax)
                end do
            case (3)
                do x = 1, size(Space%pPart)
                    ax = Space%Lev(Space%pPart(x))%a
                    cme = cme + TrL%E1%p(ax,a1) * TrL%E1%p(ax,a2)
                end do
            case (4)
                do x = 1, size(Space%nPart)
                    ax = Space%Lev(Space%nPart(x))%a
                    cme = cme + TrL%E1%n(ax,a1) * TrL%E1%n(ax,a2)
                end do
            end select
            cme = cme / real(Orb(a1)%j + 1,dp)
        end function cme
    end subroutine EMPM_cm_make

    ! Dimensions of the levels 0 .. 4 of channel (J,P) in EMPM_FULL with up to Nf phonons ...
    function EMPM_cm_dims(Ph,S2,S3,S4,Nf,J,P) result(Nb)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_states2_t), intent(in) :: S2, S3, S4
        integer, intent(in) :: Nf, J, P
        integer :: Nb(0:4)
        Nb = 0
        if (J == 0 .and. P == 1) Nb(0) = 1
        if (J <= Ph%J_max) Nb(1) = Ph%N_ch(J,P)
        if (Nf >= 2) then
            if (J <= ubound(S2%N_ch,1)) Nb(2) = S2%N_ch(J,P)
        end if
        if (Nf >= 3) then
            if (J <= ubound(S3%N_ch,1)) Nb(3) = S3%N_ch(J,P)
        end if
        if (Nf >= 4) then
            if (J <= ubound(S4%N_ch,1)) Nb(4) = S4%N_ch(J,P)
        end if
    end function EMPM_cm_dims

    ! Within-level matrices of the pp/hh part of L, contracted with the pair weights before the recursion (cheap: no
        ! rank-1 densities over the pairs are kept) ... Dl(n)%v(alpha,a) = <alpha||L_d||a> of all states of level n = 1 .. 3
        ! (level 1 from the phonon densities, level n from level n-1 by the recursion of EMPM_densityN at s = 1), the
        ! 4-phonon columns on the fly (cm_dvec) ...
    ! Rk ... rank of the operator (1: L_d with the weights wD, the default; 0: the scalar commutator C with the weights wC,
        ! positive parity, J' = J) ...
    subroutine EMPM_cm_dlevels(Pairs,Ph,C1,Cm,S1,S2,S3,Nf,Dl,Rk)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(empm_cm_t), intent(in) :: Cm
        type(empm_states2_t), intent(in) :: S1, S2, S3
        integer, intent(in) :: Nf
        type(rmat_t), intent(out) :: Dl(3)
        integer, intent(in), optional :: Rk
        integer :: g, Jb, Pb, Pp, cc, a, nb, b0, Jt, Pt, Rkk
        real(dp), allocatable :: Wt(:)
        Rkk = 1
        if (present(Rk)) Rkk = Rk
        if (Rkk == 1) then
            Wt = Cm%wD
        else
            Wt = Cm%wC
        end if
        allocate(Dl(1)%v(Ph%N,Ph%N))
        Dl(1)%v = 0.0_dp
        !$omp parallel do schedule(dynamic) private(Jb,Pb,Pp,cc,nb,b0)
        do g = 1, Ph%N
            associate(R => C1%R(g))
                do Jb = 0, Ph%J_max
                    do Pb = 1, 2
                        nb = Ph%N_ch(Jb,Pb)
                        if (nb == 0 .or. Rkk < R%Lo(Jb,Pb) .or. Rkk > R%Hi(Jb,Pb)) cycle
                        Pp = R%Pi(Jb,Pb)
                        b0 = Ph%First(Jb,Pb)
                        do cc = 1, Pairs%Na(Rkk,Pp)
                            Dl(1)%v(b0:b0+nb-1,g) = Dl(1)%v(b0:b0+nb-1,g) + Wt(Pairs%Al(Rkk,Pp)%v(cc)) * R%B(Jb,Pb)%M(Rkk)%v(cc,1:nb)
                        end do
                    end do
                end do
            end associate
        end do
        !$omp end parallel do
        if (Nf >= 2) call level(S1,S2,Dl(1)%v,Dl(2))
        if (Nf >= 3) call level(S2,S3,Dl(2)%v,Dl(3))
    contains
        subroutine level(Sl,S,Dlow,Dn)
            type(empm_states2_t), intent(in) :: Sl, S
            real(dp), intent(in) :: Dlow(:,:)
            type(rmat_t), intent(out) :: Dn
            real(dp), allocatable :: v(:), Tb(:,:)
            integer :: Ja, Pa, na, nb, nk, Nbat, i, i0, i1
            allocate(Dn%v(S%N,S%N))
            Dn%v = 0.0_dp
            ! per channel pair (Ja,Pa) -> (Jt,Pt): the t vectors of a batch of states a in parallel, then one dgemm with X ...
            do Pa = 1, 2
                do Ja = 0, S%J_max
                    na = S%N_ch(Ja,Pa)
                    if (na == 0) cycle
                    Pt = merge(3 - Pa, Pa, Rkk == 1)
                    do Jt = max(Ja - Rkk, 0), min(Ja + Rkk, S%J_max)
                        nb = S%N_ch(Jt,Pt)
                        if (nb == 0 .or. Ja + Jt < Rkk) cycle
                        nk = S%Bs(Jt,Pt)%N_full
                        Nbat = max(1, min(na, int(2.0e9_dp/(8.0_dp*real(nk,dp)))))
                        allocate(Tb(nk,Nbat))
                        do i0 = 1, na, Nbat
                            i1 = min(i0 + Nbat - 1, na)
                            !$omp parallel do schedule(dynamic) private(i,v)
                            do i = i0, i1
                                call cm_dvec(Ph,Sl,S,Dlow,Dl(1)%v,S%First(Ja,Pa)+i-1,Jt,Pt,v,Tb(:,i-i0+1),Rkk)
                            end do
                            !$omp end parallel do
                            call dgemm('T','N',nb,i1-i0+1,nk,1.0_dp,S%X(Jt,Pt)%v,nk,Tb,nk,0.0_dp, &
                                       Dn%v(S%First(Jt,Pt),S%First(Ja,Pa)+i0-1),S%N)
                        end do
                        deallocate(Tb)
                    end do
                end do
            end do
        end subroutine level
    end subroutine EMPM_cm_dlevels

    ! <beta||L_d||a> for the states beta of channel (Jt,Pt) of the level of a (states S, lower level Sl with the contracted
        ! matrix Dlow(mu,al) = <mu||L_d||al>, phonons D1(nu,l1)) ... the recursion of EMPM_densityN at s = 1 with the pair
        ! index contracted: t(k) = sum C^a(al nu) W {1 j_mu j_al; j_nu Ja Jt} (-1)^(j_mu+Ja+j_nu+1) Qa Dlow(mu,al)
        ! + sum C^a(mu l1) {1 j_nu j_l1; j_mu Ja Jt} (-1)^(j_mu+Jt+j_l1+1) Qa D1(nu,l1) over the basis k = (mu nu) of (Jt,Pt),
        ! v = X^T t ... serial ...
    subroutine cm_dvec(Ph,Sl,S,Dlow,D1,a,Jt,Pt,v,Tout,Rk)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_states2_t), intent(in) :: Sl, S
        real(dp), intent(in) :: Dlow(:,:), D1(:,:)
        integer, intent(in) :: a, Jt, Pt
        real(dp), allocatable, intent(out) :: v(:)
        real(dp), intent(out), optional :: Tout(:)
        integer, intent(in), optional :: Rk
        real(dp), allocatable :: Ca(:), t(:), W6(:,:,:)
        integer, allocatable :: Pn(:), Kn(:), Pm(:), Km(:), Cn(:), Cm2(:), Nz(:)
        integer :: Ja, Pa, ia, k, kk, ic, mu, nu, al, l1, jmu, jnu, x, y, z, n1, nlo, nz_n, J1, nb, r
        real(dp) :: Qa, w

        r = 1
        if (present(Rk)) r = Rk
        Ja = S%J(a)
        Pa = S%P(a)
        ia = a - S%First(Ja,Pa) + 1
        nb = S%N_ch(Jt,Pt)
        allocate(v(nb))
        v = 0.0_dp
        if (present(Tout)) Tout = 0.0_dp
        if (nb == 0 .or. abs(Ja - Jt) > r .or. Ja + Jt < r) return
        Qa = sqrt(real(2*Ja + 1,dp))
        n1 = Ph%N
        nlo = Sl%N
        J1 = max(Ph%J_max, Sl%J_max)
        associate(Ba => S%Bs(Ja,Pa), Bb => S%Bs(Jt,Pt))
            allocate(Ca(Ba%N_full))
            Ca = 0.0_dp
            Ca(S%Kc(Ja,Pa)%v) = S%Cc(Ja,Pa)%v(:,ia)
            Nz = pack(S%Kc(Ja,Pa)%v, S%Cc(Ja,Pa)%v(:,ia) /= 0.0_dp)
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
            Cm2 = Pm
            do kk = 1, nz_n
                Cn(Ba%La(Nz(kk))) = Cn(Ba%La(Nz(kk))) + 1
                Kn(Cn(Ba%La(Nz(kk)))) = Nz(kk)
                Cm2(Ba%Al(Nz(kk))) = Cm2(Ba%Al(Nz(kk))) + 1
                Km(Cm2(Ba%Al(Nz(kk)))) = Nz(kk)
            end do
            allocate(W6(0:J1,0:J1,0:J1))
            do z = 0, J1
                do y = 0, J1
                    do x = 0, J1
                        W6(x,y,z) = f6j(2*r,2*x,2*y,2*z,2*Ja,2*Jt)
                    end do
                end do
            end do
            allocate(t(Bb%N_full))
            t = 0.0_dp
            do k = 1, Bb%N_full
                mu = Bb%Al(k)
                nu = Bb%La(k)
                jmu = Sl%J(mu)
                jnu = Ph%J(nu)
                do ic = Pn(nu) + 1, Pn(nu+1)
                    kk = Kn(ic)
                    al = Ba%Al(kk)
                    w = W6(jmu,Sl%J(al),jnu)
                    if (w == 0.0_dp) cycle
                    t(k) = t(k) + w * phase(jmu + Ja + jnu + r) * Qa * Ca(kk) * Dlow(mu,al)
                end do
                do ic = Pm(mu) + 1, Pm(mu+1)
                    kk = Km(ic)
                    l1 = Ba%La(kk)
                    w = W6(jnu,Ph%J(l1),jmu)
                    if (w == 0.0_dp) cycle
                    t(k) = t(k) + w * phase(jmu + Jt + Ph%J(l1) + r) * Qa * Ca(kk) * D1(nu,l1)
                end do
            end do
            if (present(Tout)) then
                Tout = t
            else
                v = matmul(t,S%X(Jt,Pt)%v)
            end if
        end associate
    end subroutine cm_dvec

    ! The three block types of <(Jt,Pt) states|| L ||(J,P) states> over the bases [0 ; 1 ; .. ; Nf] of both channels ...
        ! Dl the contracted within-level matrices of levels 1 .. 3 (EMPM_cm_dlevels), 4-phonon columns by cm_dvec ...
    subroutine EMPM_cm_blocks(Ph,Cm,S1,S2,S3,S4,Dl,Nf,J,P,Jt,Pt,Bph,Bhp,Bd,Top)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cm_t), intent(in) :: Cm
        type(empm_states2_t), intent(in) :: S1, S2, S3, S4
        type(rmat_t), intent(in) :: Dl(3)
        integer, intent(in) :: Nf, J, P, Jt, Pt
        real(dp), allocatable, intent(out) :: Bph(:,:), Bhp(:,:), Bd(:,:)
        ! Top ... .false.: no pp/hh block within the top level Nf (proj only, Dl(Nf) not computed) ...
        logical, intent(in), optional :: Top
        logical :: Tp
        integer :: Nb(0:4), Mb(0:4), Oc(0:5), Ot(0:5), q, i, l0, nl, n, nk, Nbat, i0, i1
        real(dp), allocatable :: v(:), Tb(:,:)
        real(dp) :: r3

        r3 = 1.0_dp / sqrt(3.0_dp)
        Tp = .true.
        if (present(Top)) Tp = Top
        Nb = EMPM_cm_dims(Ph,S2,S3,S4,Nf,J,P)
        Mb = EMPM_cm_dims(Ph,S2,S3,S4,Nf,Jt,Pt)
        Oc(0) = 0
        Ot(0) = 0
        do q = 0, 4
            Oc(q+1) = Oc(q) + Nb(q)
            Ot(q+1) = Ot(q) + Mb(q)
        end do
        allocate(Bph(Ot(5),Oc(5)), Bhp(Ot(5),Oc(5)), Bd(Ot(5),Oc(5)))
        Bph = 0.0_dp
        Bhp = 0.0_dp
        Bd = 0.0_dp
        l0 = 0
        nl = 0
        if (Ph%J_max >= 1) then
            l0 = Ph%First(1,2)
            nl = Ph%N_ch(1,2)
        end if

        ! 0 -> 1 & 1 -> 0 ...
        if (Nb(0) == 1 .and. Jt == 1 .and. Pt == 2) Bph(Ot(1)+1:Ot(1)+Mb(1),1) = Cm%bL(Ph%First(1,2):Ph%First(1,2)+Mb(1)-1)
        if (Mb(0) == 1 .and. J == 1 .and. P == 2) Bhp(1,Oc(1)+1:Oc(1)+Nb(1)) = Cm%bR(Ph%First(1,2):Ph%First(1,2)+Nb(1)-1)

        ! Level 1 columns ... pp/hh within the phonons & ph to level 2 ...
        do i = 1, Nb(1)
            if (Mb(1) > 0 .and. (Tp .or. Nf > 1)) Bd(Ot(1)+1:Ot(1)+Mb(1),Oc(1)+i) = Dl(1)%v(Ph%First(Jt,Pt):Ph%First(Jt,Pt)+Mb(1)-1,Ph%First(J,P)+i-1)
            if (Nf >= 2) call phcol(S2,Ph%First(J,P)+i-1,Mb(2),Ot(2),Oc(1)+i)
        end do

        ! Levels n = 2 .. Nf ... hp to n-1, pp/hh within n, ph to n+1 ...
        do n = 2, Nf
            !$omp parallel do schedule(dynamic) private(i,v)
            do i = 1, Nb(n)
                select case (n)
                case (2)
                    call hpcol(S2,S1,i,Mb(1),Ot(1),Oc(2)+i)
                    if (Mb(2) > 0 .and. (Tp .or. Nf > 2)) Bd(Ot(2)+1:Ot(2)+Mb(2),Oc(2)+i) = Dl(2)%v(S2%First(Jt,Pt):S2%First(Jt,Pt)+Mb(2)-1, &
                                                                             S2%First(J,P)+i-1)
                    if (Nf >= 3) call phcol(S3,S2%First(J,P)+i-1,Mb(3),Ot(3),Oc(2)+i)
                case (3)
                    call hpcol(S3,S2,i,Mb(2),Ot(2),Oc(3)+i)
                    if (Mb(3) > 0 .and. (Tp .or. Nf > 3)) Bd(Ot(3)+1:Ot(3)+Mb(3),Oc(3)+i) = Dl(3)%v(S3%First(Jt,Pt):S3%First(Jt,Pt)+Mb(3)-1, &
                                                                             S3%First(J,P)+i-1)
                    if (Nf >= 4) call phcol(S4,S3%First(J,P)+i-1,Mb(4),Ot(4),Oc(3)+i)
                case (4)
                    call hpcol(S4,S3,i,Mb(3),Ot(3),Oc(4)+i)
                end select
            end do
            !$omp end parallel do
        end do
        ! 4-phonon pp/hh columns in batches ... t vectors (basis of (Jt,Pt)) in parallel, then one dgemm with X ...
        if (Nf >= 4 .and. Nb(4) > 0 .and. Mb(4) > 0 .and. Tp) then
            nk = S4%Bs(Jt,Pt)%N_full
            Nbat = max(1, min(Nb(4), int(2.0e9_dp/(8.0_dp*real(nk,dp)))))
            allocate(Tb(nk,Nbat))
            do i0 = 1, Nb(4), Nbat
                i1 = min(i0 + Nbat - 1, Nb(4))
                !$omp parallel do schedule(dynamic) private(i,v)
                do i = i0, i1
                    call cm_dvec(Ph,S3,S4,Dl(3)%v,Dl(1)%v,S4%First(J,P)+i-1,Jt,Pt,v,Tb(:,i-i0+1))
                end do
                !$omp end parallel do
                call dgemm('T','N',Mb(4),i1-i0+1,nk,1.0_dp,S4%X(Jt,Pt)%v,nk,Tb,nk,0.0_dp,Bd(Ot(4)+1,Oc(4)+i0),Ot(5))
            end do
            deallocate(Tb)
        end if

    contains

        ! ph ... column c = state a (global index of its level), rows the (Jt,Pt) states of the next level Su ...
        subroutine phcol(Su,a,m,Ot,c)
            type(empm_states2_t), intent(in) :: Su
            integer, intent(in) :: a, m, Ot, c
            integer :: ll, k
            if (m == 0) return
            do ll = l0, l0 + nl - 1
                if (Cm%bL(ll) == 0.0_dp) cycle
                k = Su%Bs(Jt,Pt)%K(a,ll)
                if (k == 0) cycle
                Bph(Ot+1:Ot+m,c) = Bph(Ot+1:Ot+m,c) + r3 * Cm%bL(ll) * Su%X(Jt,Pt)%v(k,1:m)
            end do
        end subroutine phcol

        ! hp ... column c = state i of channel (J,P) of level Sn, rows the (Jt,Pt) states a' of the level below Sl ...
        subroutine hpcol(Sn,Sl,i,m,Ot,c)
            type(empm_states2_t), intent(in) :: Sn, Sl
            integer, intent(in) :: i, m, Ot, c
            integer :: ll, k, ia, ap
            if (m == 0) return
            do ia = 1, m
                ap = Sl%First(Jt,Pt) + ia - 1
                do ll = l0, l0 + nl - 1
                    if (Cm%bR(ll) == 0.0_dp) cycle
                    k = Sn%Bs(J,P)%K(ap,ll)
                    if (k == 0) cycle
                    Bhp(Ot+ia,c) = Bhp(Ot+ia,c) + r3 * Cm%bR(ll) * Sn%X(J,P)%v(k,i)
                end do
            end do
        end subroutine hpcol

    end subroutine EMPM_cm_blocks

    ! G(alpha,gamma) = <alpha||L_ph||gamma> for the states alpha of channel (J,P) of level n (Su) & gamma of channel
        ! (Jt,Pt) of level n-1 (Sl; for n = 1 the ground state, for n = 2 the phonons) ... the ph block of EMPM_cm_blocks
        ! (phcol) with rows & columns of the other channel ...
    function EMPM_cm_up(Ph,Cm,Sl,Su,n,J,P,Jt,Pt) result(G)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cm_t), intent(in) :: Cm
        type(empm_states2_t), intent(in) :: Sl, Su
        integer, intent(in) :: n, J, P, Jt, Pt
        real(dp), allocatable :: G(:,:)
        integer :: m, ml, ia, a, ll, k
        if (n == 1) then
            m = 0
            if (J <= Ph%J_max) m = Ph%N_ch(J,P)
            ml = merge(1, 0, Jt == 0 .and. Pt == 1)
            allocate(G(m,ml))
            G = 0.0_dp
            if (ml == 1 .and. J == 1 .and. P == 2) G(:,1) = Cm%bL(Ph%First(1,2):Ph%First(1,2)+m-1)
            return
        end if
        m = 0
        if (J <= ubound(Su%N_ch,1)) m = Su%N_ch(J,P)
        ml = 0
        if (Jt <= ubound(Sl%N_ch,1)) ml = Sl%N_ch(Jt,Pt)
        allocate(G(m,ml))
        G = 0.0_dp
        if (m == 0 .or. ml == 0 .or. Ph%J_max < 1) return
        do ia = 1, ml
            a = Sl%First(Jt,Pt) + ia - 1
            do ll = Ph%First(1,2), Ph%First(1,2) + Ph%N_ch(1,2) - 1
                if (Cm%bL(ll) == 0.0_dp) cycle
                k = Su%Bs(J,P)%K(a,ll)
                if (k == 0) cycle
                G(:,ia) = G(:,ia) + Cm%bL(ll) / sqrt(3.0_dp) * Su%X(J,P)%v(k,1:m)
            end do
        end do
    end function EMPM_cm_up

    ! Missing part of N_cm on level n of channel (J,P) (the intermediate states of level n+1, reached only by L_ph):
        ! Mt = G2/(A(2J+1)) [sum_Jt G G^T + (2J+1) C0 1 + sqrt(2J+1) Dc] , G = EMPM_cm_up(level n), Dc = <||C_d||> the rank-0
        ! within-level matrix of level n (EMPM_cm_dlevels with Rk = 0) ... exact if level n-1 is complete (the identity
        ! L_ph+ L_ph = L_ph L_ph+ + C with L_ph+ mapping level n to level n-1) ...
    function EMPM_cm_top(Ph,Cm,Sl,Su,Dc,n,J,P) result(Mt)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cm_t), intent(in) :: Cm
        type(empm_states2_t), intent(in) :: Sl, Su
        real(dp), intent(in) :: Dc(:,:)
        integer, intent(in) :: n, J, P
        real(dp), allocatable :: Mt(:,:), G(:,:)
        integer :: m, Jt, Pt, i, f
        if (n == 1) then
            m = 0
            if (J <= Ph%J_max) m = Ph%N_ch(J,P)
            f = 0
            if (J <= Ph%J_max) f = Ph%First(J,P)
        else
            m = 0
            if (J <= ubound(Su%N_ch,1)) m = Su%N_ch(J,P)
            f = 0
            if (m > 0) f = Su%First(J,P)
        end if
        allocate(Mt(m,m))
        Mt = 0.0_dp
        if (m == 0) return
        Pt = 3 - P
        do Jt = max(J - 1, 0), J + 1
            G = EMPM_cm_up(Ph,Cm,Sl,Su,n,J,P,Jt,Pt)
            if (size(G,2) > 0) Mt = Mt + matmul(G,transpose(G))
        end do
        Mt = Mt + sqrt(real(2*J + 1,dp)) * Dc(f:f+m-1,f:f+m-1)
        do i = 1, m
            Mt(i,i) = Mt(i,i) + real(2*J + 1,dp) * Cm%C0
        end do
        Mt = Mt * Cm%G2 / real(Cm%A*(2*J + 1),dp)
    end function EMPM_cm_top

end module nj_empm_cm
