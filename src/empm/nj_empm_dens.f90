! 1-phonon densities & redefined phonon interaction of the EMPM computed on the fly ... counterparts of phon_dens1 (rop, roh)
    ! & phon_int (vintp, vinth) of the original EMPM codes, nothing is stored on disk ...
    ! Pair space s = (k,i1,i2) ... block k = 1 proton particles, 2 neutron particles, 3 proton holes, 4 neutron holes,
        ! s = Off(k) + i1 + (i2-1)*Nk(k) with i1, i2 local level indices of the block ...
    ! Phonon block of gamma ... B(Jb,Pb)%M(sigma)%v(c,beta) for the phonons beta of channel (Jb,Pb) & |J_gamma-Jb| <= sigma <=
        ! min(J_gamma+Jb, L_max) ... rho_{gamma beta}(sigma; i1,i2) or V_{gamma beta}(sigma; i1,i2) stored only for the pairs
        ! allowed by sigma: c = 1 .. Na(sigma,pi) indexes the pairs s = Al(sigma,pi)%v(c) with triangle(j_i1,j_i2,sigma) &
        ! parity (-1)^(l_i1+l_i2) = pi (1 even, 2 odd) = P_gamma P_beta, the only non-zero entries (16% of S for 208Pb) ...
    ! Densities (j's of the levels, J's of the phonons, c the TDA amplitudes) ...
        ! particles  rho(s;p1,p2) =   sum_h (-1)^(Jg+s+j_h+j_p1) Q {Jb s Jg; j_p2 j_h j_p1} c^g_{p2 h} c^b_{p1 h}
        ! holes      rho(s;h1,h2) = - sum_p (-1)^(Jb+j_p+j_h2)   Q {Jb s Jg; j_h1 j_p j_h2} c^g_{p h1} c^b_{p h2}
        ! Q = sqrt((2Jb+1)(2Jg+1)(2s+1)) ...
    ! Interaction ... V(s;k,l) = (-1)^(Jb+s) sum_(i1,i2) kappa rho(s;i1,i2) F^s(i1,i2,k,l), kappa = 1/2 for the density of
        ! the same block as (k,l), plus the pn exchange terms of vintp & vinth ...
module nj_empm_dens
    use nj_kinds
    use nj_types
    use nj_utils
    use nj_wigner, only : f6j
    use nj_empm_space
    use nj_empm_fmat, only : fmat_t, Fmat_ME
    use nj_empm_phonons, only : empm_phonons_t
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

    type, public :: r3_t
        real(dp), allocatable :: v(:,:,:)
    end type r3_t

    ! Pair space ... Pos(k)%v local level -> position in Space%Lev, jx(i,k) 2j, ix(i,k) index of 2j in the 6j table,
        ! lx(i,k) l, Tr(s) the transposed pair (i2,i1), Sg(s) = (-1)^(j_i1-j_i2), Na, Al & Ci(s,sigma,pi) (0 = not allowed)
        ! the allowed pairs & their compact index, Fkc(sigma,pi)%v(c,c') = kappa F^sigma(pair c', pair c) between allowed pairs,
        ! T6(a,b,c,Jb,s,Jg) = {Jb s Jg; j_a j_b j_c}, Ex(t)%v(:,:,Jb) exchange F^Jb blocks (see EMPM_interaction) ...
    type, public :: empm_pairs_t
        integer :: S = 0, L_max = 0, J_max = 0, nj = 0
        integer :: Off(4) = 0, Nk(4) = 0
        type(ivec_t) :: Pos(4)
        integer, allocatable :: jx(:,:), ix(:,:), lx(:,:), Tr(:), Na(:,:), Ci(:,:,:)
        type(ivec_t), allocatable :: Al(:,:)
        type(rmat_t), allocatable :: Fkc(:,:)
        real(dp), allocatable :: Sg(:), T6(:,:,:,:,:,:)
        type(r3_t) :: ExP(2), ExH(2)
    end type empm_pairs_t

    type, public :: sig_block_t
        type(rmat_t), allocatable :: M(:)
    end type sig_block_t

    ! Pi(Jb,Pb) parity index of the block's pairs ...
    type, public :: phonon_block_t
        integer :: g = 0
        integer, allocatable :: Lo(:,:), Hi(:,:), Pi(:,:)
        type(sig_block_t), allocatable :: B(:,:)
    end type phonon_block_t

    public :: EMPM_pairs_make, EMPM_density, EMPM_interaction, EMPM_density_transposed, block_parity

contains

    subroutine EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_pairs_t), intent(out) :: Pairs
        integer :: k, k2, i1, i2, s, s2, Lam, nmax, a, b, c, c2, Jb, Sig, Jg, t, tt, p, h, l, ph2, hp, n_p, n_h, m_p, m_h, Pp, n
        real(dp) :: Kappa

        Pairs%Pos(1)%v = Space%pPart
        Pairs%Pos(2)%v = Space%nPart
        Pairs%Pos(3)%v = Space%pHole
        Pairs%Pos(4)%v = Space%nHole
        nmax = 0
        do k = 1, 4
            Pairs%Nk(k) = size(Pairs%Pos(k)%v)
            nmax = max(nmax, Pairs%Nk(k))
            if (k > 1) Pairs%Off(k) = Pairs%Off(k-1) + Pairs%Nk(k-1)**2
        end do
        Pairs%S = Pairs%Off(4) + Pairs%Nk(4)**2
        Pairs%L_max = F_ph%L_max
        Pairs%J_max = Ph%J_max
        Pairs%nj = (maxval(Space%Lev(1:Space%N_lev)%j) + 1)/2

        allocate(Pairs%jx(nmax,4), Pairs%ix(nmax,4), Pairs%lx(nmax,4), Pairs%Tr(Pairs%S), Pairs%Sg(Pairs%S))
        do k = 1, 4
            do i1 = 1, Pairs%Nk(k)
                Pairs%jx(i1,k) = Space%Lev(Pairs%Pos(k)%v(i1))%j
                Pairs%lx(i1,k) = Space%Lev(Pairs%Pos(k)%v(i1))%l
                Pairs%ix(i1,k) = (Pairs%jx(i1,k) + 1)/2
            end do
            do i2 = 1, Pairs%Nk(k)
                do i1 = 1, Pairs%Nk(k)
                    s = pair(k,i1,i2)
                    Pairs%Tr(s) = pair(k,i2,i1)
                    Pairs%Sg(s) = phase((Pairs%jx(i1,k) - Pairs%jx(i2,k))/2)
                end do
            end do
        end do

        ! Allowed pairs per (sigma, parity) ...
        allocate(Pairs%Na(0:Pairs%L_max,2), Pairs%Al(0:Pairs%L_max,2), Pairs%Ci(Pairs%S,0:Pairs%L_max,2))
        Pairs%Ci = 0
        do Sig = 0, Pairs%L_max
            do Pp = 1, 2
                n = 0
                do k = 1, 4
                    do i2 = 1, Pairs%Nk(k)
                        do i1 = 1, Pairs%Nk(k)
                            if (allowed(k,i1,i2,Sig,Pp)) then
                                n = n + 1
                                Pairs%Ci(pair(k,i1,i2),Sig,Pp) = n
                            end if
                        end do
                    end do
                end do
                Pairs%Na(Sig,Pp) = n
                allocate(Pairs%Al(Sig,Pp)%v(n))
                do s = 1, Pairs%S
                    if (Pairs%Ci(s,Sig,Pp) > 0) Pairs%Al(Sig,Pp)%v(Pairs%Ci(s,Sig,Pp)) = s
                end do
            end do
        end do

        ! kappa F^L(pair c', pair c) between the allowed pairs ... c the pair of the interaction, c' the pair of the density ...
        allocate(Pairs%Fkc(0:Pairs%L_max,2))
        do Lam = 0, Pairs%L_max
            do Pp = 1, 2
                n = Pairs%Na(Lam,Pp)
                allocate(Pairs%Fkc(Lam,Pp)%v(n,n))
                !$omp parallel do schedule(dynamic) private(c2,s,s2,k,k2,i1,i2,a,b,Kappa)
                do c = 1, n
                    s = Pairs%Al(Lam,Pp)%v(c)
                    call unpair(s,k,i1,i2)
                    do c2 = 1, n
                        s2 = Pairs%Al(Lam,Pp)%v(c2)
                        call unpair(s2,k2,a,b)
                        Kappa = merge(0.5_dp, 1.0_dp, k2 == k)
                        Pairs%Fkc(Lam,Pp)%v(c,c2) = Kappa * Fmat_ME(F_ph,Space,Pairs%Pos(k2)%v(a),Pairs%Pos(k2)%v(b), &
                                                                    Pairs%Pos(k)%v(i1),Pairs%Pos(k)%v(i2),Lam)
                    end do
                end do
                !$omp end parallel do
            end do
        end do

        ! 6j table {Jb s Jg; j_a j_b j_c} ...
        allocate(Pairs%T6(Pairs%nj,Pairs%nj,Pairs%nj,0:Pairs%J_max,0:Pairs%L_max,0:Pairs%J_max))
        !$omp parallel do collapse(2) private(Sig,Jb,a,b,c)
        do Jg = 0, Pairs%J_max
            do Sig = 0, Pairs%L_max
                do Jb = 0, Pairs%J_max
                    do c = 1, Pairs%nj
                        do b = 1, Pairs%nj
                            do a = 1, Pairs%nj
                                Pairs%T6(a,b,c,Jb,Sig,Jg) = f6j(2*Jb,2*Sig,2*Jg,2*a-1,2*b-1,2*c-1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! pn exchange blocks ... for isospin t (1 p, 2 n) & t' the other one, columns (p',h') of t' amplitudes ...
            ! ExP(t)%v((h,l),(p',h'),Jb) = F^Jb(p',h',h,l) & ExH(t)%v((h1,p),(p',h'),Jb) = F^Jb(p',h',h1,p) ...
        do t = 1, 2
            tt = 3 - t
            n_p = Pairs%Nk(t)
            n_h = Pairs%Nk(t+2)
            m_p = Pairs%Nk(tt)
            m_h = Pairs%Nk(tt+2)
            allocate(Pairs%ExP(t)%v(n_h*n_p,m_p*m_h,0:Pairs%J_max), Pairs%ExH(t)%v(n_h*n_p,m_p*m_h,0:Pairs%J_max))
            do Jb = 0, Pairs%J_max
                do hp = 1, m_h
                    do ph2 = 1, m_p
                        do l = 1, n_p
                            do h = 1, n_h
                                Pairs%ExP(t)%v(h+(l-1)*n_h,ph2+(hp-1)*m_p,Jb) = &
                                    Fmat_ME(F_ph,Space,Pairs%Pos(tt)%v(ph2),Pairs%Pos(tt+2)%v(hp),Pairs%Pos(t+2)%v(h),Pairs%Pos(t)%v(l),Jb)
                            end do
                        end do
                        do p = 1, n_p
                            do h = 1, n_h
                                Pairs%ExH(t)%v(h+(p-1)*n_h,ph2+(hp-1)*m_p,Jb) = &
                                    Fmat_ME(F_ph,Space,Pairs%Pos(tt)%v(ph2),Pairs%Pos(tt+2)%v(hp),Pairs%Pos(t+2)%v(h),Pairs%Pos(t)%v(p),Jb)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        write(*,'(a)') TAB // 'EMPM pair space ... S = ' // istr(Pairs%S) // ' , sigma_max = ' // istr(Pairs%L_max) // &
                       ' , allowed pairs per (sigma,parity) on average ' // istr(nint(real(sum(Pairs%Na),dp)/real(size(Pairs%Na),dp)))

    contains

        pure integer function pair(k,i1,i2)
            integer, intent(in) :: k, i1, i2
            pair = Pairs%Off(k) + i1 + (i2 - 1)*Pairs%Nk(k)
        end function pair

        pure subroutine unpair(s,k,i1,i2)
            integer, intent(in) :: s
            integer, intent(out) :: k, i1, i2
            k = 4
            do while (s <= Pairs%Off(k))
                k = k - 1
            end do
            i1 = mod(s - Pairs%Off(k) - 1, Pairs%Nk(k)) + 1
            i2 = (s - Pairs%Off(k) - 1)/Pairs%Nk(k) + 1
        end subroutine unpair

        pure logical function allowed(k,i1,i2,Sig,Pp)
            integer, intent(in) :: k, i1, i2, Sig, Pp
            allowed = abs(Pairs%jx(i1,k) - Pairs%jx(i2,k))/2 <= Sig .and. Sig <= (Pairs%jx(i1,k) + Pairs%jx(i2,k))/2 .and. &
                      mod(Pairs%lx(i1,k) + Pairs%lx(i2,k),2) + 1 == Pp
        end function allowed

    end subroutine EMPM_pairs_make

    ! Parity index of the pairs of the block (g, channel Pb) ... 1 for P_g = P_b ...
    pure integer function block_parity(Ph,g,Pb)
        type(empm_phonons_t), intent(in) :: Ph
        integer, intent(in) :: g, Pb
        block_parity = merge(1, 2, Ph%P(g) == Pb)
    end function block_parity

    ! Empty phonon block of gamma with the sigma ranges of all channels ...
    subroutine block_allocate(Pairs,Ph,g,X)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        integer, intent(in) :: g
        type(phonon_block_t), intent(out) :: X
        integer :: Jb, Pb, Sig
        X%g = g
        allocate(X%Lo(0:Ph%J_max,2), X%Hi(0:Ph%J_max,2), X%Pi(0:Ph%J_max,2), X%B(0:Ph%J_max,2))
        do Jb = 0, Ph%J_max
            do Pb = 1, 2
                X%Lo(Jb,Pb) = abs(Ph%J(g) - Jb)
                X%Hi(Jb,Pb) = min(Ph%J(g) + Jb, Pairs%L_max)
                X%Pi(Jb,Pb) = block_parity(Ph,g,Pb)
                if (Ph%N_ch(Jb,Pb) == 0) X%Hi(Jb,Pb) = X%Lo(Jb,Pb) - 1
                allocate(X%B(Jb,Pb)%M(X%Lo(Jb,Pb):X%Hi(Jb,Pb)))
                do Sig = X%Lo(Jb,Pb), X%Hi(Jb,Pb)
                    allocate(X%B(Jb,Pb)%M(Sig)%v(Pairs%Na(Sig,X%Pi(Jb,Pb)),Ph%N_ch(Jb,Pb)))
                    X%B(Jb,Pb)%M(Sig)%v = 0.0_dp
                end do
            end do
        end do
    end subroutine block_allocate

    ! Densities of gamma = g for all beta as dgemm's ... for each channel (Jb,Pb), isospin t & the beta-side level ...
        ! particles  R(p2,s,beta) = sum_h M_p1(p2,s,h) c^beta(p1,h) , M_p1(p2,s,h) = (-1)^(Jg+s+j_h+j_p1) Q {Jb s Jg; j_p2 j_h j_p1} c^g(p2,h)
        ! holes      R(h1,s,beta) = sum_p M_h2(h1,s,p) c^beta(p,h2) , M_h2(h1,s,p) = -(-1)^(Jb+j_p+j_h2) Q {Jb s Jg; j_h1 j_p j_h2} c^g(p,h1)
        ! Serial (called inside the parallel loops of the callers) ...
    subroutine EMPM_density(Pairs,Ph,g,Rho)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        integer, intent(in) :: g
        type(phonon_block_t), intent(out) :: Rho
        real(dp), allocatable :: M(:,:,:), R(:,:,:)
        integer :: Jg, Jb, Pb, n, nS, b0, Sig, t, kp, kh, n_p, n_h, p1, p2, h, h1, h2, p, lo
        real(dp) :: Q

        call block_allocate(Pairs,Ph,g,Rho)
        Jg = Ph%J(g)
        do Jb = 0, Ph%J_max
            do Pb = 1, 2
                lo = Rho%Lo(Jb,Pb)
                nS = Rho%Hi(Jb,Pb) - lo + 1
                n = Ph%N_ch(Jb,Pb)
                if (nS <= 0 .or. n == 0) cycle
                b0 = Ph%First(Jb,Pb)
                do t = 1, 2
                    kp = t
                    kh = t + 2
                    n_p = Pairs%Nk(kp)
                    n_h = Pairs%Nk(kh)
                    if (n_p == 0 .or. n_h == 0) cycle
                    ! particles ...
                    allocate(M(n_p,lo:lo+nS-1,n_h), R(n_p,lo:lo+nS-1,n))
                    do p1 = 1, n_p
                        do h = 1, n_h
                            do Sig = lo, lo + nS - 1
                                Q = sqrt(real((2*Jb + 1)*(2*Jg + 1)*(2*Sig + 1),dp)) * &
                                    phase(Jg + Sig + (Pairs%jx(h,kh) + Pairs%jx(p1,kp))/2)
                                do p2 = 1, n_p
                                    M(p2,Sig,h) = Q * Pairs%T6(Pairs%ix(p2,kp),Pairs%ix(h,kh),Pairs%ix(p1,kp),Jb,Sig,Jg) * amp(t,p2,h,g)
                                end do
                            end do
                        end do
                        if (t == 1) then
                            call dgemm('N','N',n_p*nS,n,n_h,1.0_dp,M,n_p*nS,Ph%Cpt(1,p1,b0),n_h*n_p,0.0_dp,R,n_p*nS)
                        else
                            call dgemm('N','N',n_p*nS,n,n_h,1.0_dp,M,n_p*nS,Ph%Cnt(1,p1,b0),n_h*n_p,0.0_dp,R,n_p*nS)
                        end if
                        do p2 = 1, n_p
                            call put(Pairs%Off(kp)+p1+(p2-1)*n_p,R(p2,:,:))
                        end do
                    end do
                    deallocate(M, R)
                    ! holes ...
                    allocate(M(n_h,lo:lo+nS-1,n_p), R(n_h,lo:lo+nS-1,n))
                    do h2 = 1, n_h
                        do p = 1, n_p
                            do Sig = lo, lo + nS - 1
                                Q = -sqrt(real((2*Jb + 1)*(2*Jg + 1)*(2*Sig + 1),dp)) * phase(Jb + (Pairs%jx(p,kp) + Pairs%jx(h2,kh))/2)
                                do h1 = 1, n_h
                                    M(h1,Sig,p) = Q * Pairs%T6(Pairs%ix(h1,kh),Pairs%ix(p,kp),Pairs%ix(h2,kh),Jb,Sig,Jg) * amp(t,p,h1,g)
                                end do
                            end do
                        end do
                        if (t == 1) then
                            call dgemm('N','N',n_h*nS,n,n_p,1.0_dp,M,n_h*nS,Ph%Cp(1,h2,b0),n_p*n_h,0.0_dp,R,n_h*nS)
                        else
                            call dgemm('N','N',n_h*nS,n,n_p,1.0_dp,M,n_h*nS,Ph%Cn(1,h2,b0),n_p*n_h,0.0_dp,R,n_h*nS)
                        end if
                        do h1 = 1, n_h
                            call put(Pairs%Off(kh)+h1+(h2-1)*n_h,R(h1,:,:))
                        end do
                    end do
                    deallocate(M, R)
                end do
            end do
        end do

    contains

        ! Entries of pair s for all sigma of the block (Jb,Pb) ... only the allowed ones are kept ...
        subroutine put(s,Rs)
            integer, intent(in) :: s
            real(dp), intent(in) :: Rs(lo:,:)
            integer :: Sg, c
            do Sg = lo, lo + nS - 1
                c = Pairs%Ci(s,Sg,Rho%Pi(Jb,Pb))
                if (c > 0) Rho%B(Jb,Pb)%M(Sg)%v(c,:) = Rs(Sg,:)
            end do
        end subroutine put

        pure real(dp) function amp(t,p,h,i)
            integer, intent(in) :: t, p, h, i
            if (t == 1) then
                amp = Ph%Cp(p,h,i)
            else
                amp = Ph%Cn(p,h,i)
            end if
        end function amp

    end subroutine EMPM_density

    ! rho~(s,sigma,beta) = (-1)^(j_i1-j_i2) rho(s^T,sigma,beta) ... the density entering the D matrix (s^T allowed with s) ...
    subroutine EMPM_density_transposed(Pairs,Rho,Rt)
        type(empm_pairs_t), intent(in) :: Pairs
        type(phonon_block_t), intent(in) :: Rho
        type(phonon_block_t), intent(out) :: Rt
        integer :: Jb, Pb, Sig, c, s, Pp
        Rt%g = Rho%g
        Rt%Lo = Rho%Lo
        Rt%Hi = Rho%Hi
        Rt%Pi = Rho%Pi
        allocate(Rt%B(lbound(Rho%B,1):ubound(Rho%B,1),2))
        do Jb = lbound(Rho%B,1), ubound(Rho%B,1)
            do Pb = 1, 2
                Pp = Rho%Pi(Jb,Pb)
                allocate(Rt%B(Jb,Pb)%M(Rho%Lo(Jb,Pb):Rho%Hi(Jb,Pb)))
                do Sig = Rho%Lo(Jb,Pb), Rho%Hi(Jb,Pb)
                    allocate(Rt%B(Jb,Pb)%M(Sig)%v, mold=Rho%B(Jb,Pb)%M(Sig)%v)
                    do c = 1, Pairs%Na(Sig,Pp)
                        s = Pairs%Al(Sig,Pp)%v(c)
                        Rt%B(Jb,Pb)%M(Sig)%v(c,:) = Pairs%Sg(s) * Rho%B(Jb,Pb)%M(Sig)%v(Pairs%Ci(Pairs%Tr(s),Sig,Pp),:)
                    end do
                end do
            end do
        end do
    end subroutine EMPM_density_transposed

    ! Redefined interaction of gamma from its densities ... direct part by dgemm with Fk, pn exchange as in vintp & vinth ...
        ! particles of t   V(s;p,l)  -= sum_h  (-1)^(s+Jg+j_p+j_l)   Q {Jb s Jg; j_p j_h j_l}  c^g_{p h}  f_b(h,l)
        ! holes of t       V(s;h1,h) += sum_p  (-1)^(Jb+j_h1+j_h)    Q {Jb s Jg; j_h j_p j_h1} c^g_{p h}  g_b(h1,p)
        ! f_b(h,l) = sum_(p'h' of t') c^b_{p'h'} F^Jb(p',h',h,l) & g_b(h1,p) = sum_(p'h' of t') c^b_{p'h'} F^Jb(p',h',h1,p) ...
    subroutine EMPM_interaction(Pairs,Ph,Rho,V)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(phonon_block_t), intent(in) :: Rho
        type(phonon_block_t), intent(out) :: V
        real(dp), allocatable :: Fx(:,:), Gx(:,:)
        integer :: Jg, Jb, Pb, Sig, n, nS, g, t, kp, kh, n_p, n_h, nu, p, h, l, h1, b0, Pp, Ka, c
        real(dp) :: Q, Sum

        g = Rho%g
        Jg = Ph%J(g)
        call block_allocate(Pairs,Ph,g,V)
        do Jb = 0, Ph%J_max
            do Pb = 1, 2
                n = Ph%N_ch(Jb,Pb)
                nS = Rho%Hi(Jb,Pb) - Rho%Lo(Jb,Pb) + 1
                if (nS <= 0 .or. n == 0) cycle
                Pp = Rho%Pi(Jb,Pb)
                do Sig = Rho%Lo(Jb,Pb), Rho%Hi(Jb,Pb)
                    Ka = Pairs%Na(Sig,Pp)
                    if (Ka == 0) cycle
                    call dgemm('N','N',Ka,n,Ka,phase(Jb + Sig),Pairs%Fkc(Sig,Pp)%v,Ka, &
                               Rho%B(Jb,Pb)%M(Sig)%v,Ka,0.0_dp,V%B(Jb,Pb)%M(Sig)%v,Ka)
                end do
                b0 = Ph%First(Jb,Pb)
                do t = 1, 2
                    kp = t
                    kh = t + 2
                    n_p = Pairs%Nk(kp)
                    n_h = Pairs%Nk(kh)
                    if (n_p*n_h == 0 .or. size(Pairs%ExP(t)%v,2) == 0) cycle
                    allocate(Fx(n_h*n_p,n), Gx(n_h*n_p,n))
                    if (t == 1) then
                        Fx = matmul(Pairs%ExP(t)%v(:,:,Jb),reshape(Ph%Cn(:,:,b0:b0+n-1),[size(Pairs%ExP(t)%v,2),n]))
                        Gx = matmul(Pairs%ExH(t)%v(:,:,Jb),reshape(Ph%Cn(:,:,b0:b0+n-1),[size(Pairs%ExH(t)%v,2),n]))
                    else
                        Fx = matmul(Pairs%ExP(t)%v(:,:,Jb),reshape(Ph%Cp(:,:,b0:b0+n-1),[size(Pairs%ExP(t)%v,2),n]))
                        Gx = matmul(Pairs%ExH(t)%v(:,:,Jb),reshape(Ph%Cp(:,:,b0:b0+n-1),[size(Pairs%ExH(t)%v,2),n]))
                    end if
                    do nu = 1, n
                        do Sig = Rho%Lo(Jb,Pb), Rho%Hi(Jb,Pb)
                            Q = sqrt(real((2*Jb + 1)*(2*Jg + 1)*(2*Sig + 1),dp))
                            do l = 1, n_p
                                do p = 1, n_p
                                    Sum = 0.0_dp
                                    do h = 1, n_h
                                        Sum = Sum + phase(Sig + Jg + (Pairs%jx(p,kp) + Pairs%jx(l,kp))/2) * &
                                                    Pairs%T6(Pairs%ix(p,kp),Pairs%ix(h,kh),Pairs%ix(l,kp),Jb,Sig,Jg) * &
                                                    amp(t,p,h,g) * Fx(h+(l-1)*n_h,nu)
                                    end do
                                    c = Pairs%Ci(Pairs%Off(kp)+p+(l-1)*n_p,Sig,Pp)
                                    if (c > 0) V%B(Jb,Pb)%M(Sig)%v(c,nu) = V%B(Jb,Pb)%M(Sig)%v(c,nu) - Q * Sum
                                end do
                            end do
                            do h = 1, n_h
                                do h1 = 1, n_h
                                    Sum = 0.0_dp
                                    do p = 1, n_p
                                        Sum = Sum + phase(Jb + (Pairs%jx(h1,kh) + Pairs%jx(h,kh))/2) * &
                                                    Pairs%T6(Pairs%ix(h,kh),Pairs%ix(p,kp),Pairs%ix(h1,kh),Jb,Sig,Jg) * &
                                                    amp(t,p,h,g) * Gx(h1+(p-1)*n_h,nu)
                                    end do
                                    c = Pairs%Ci(Pairs%Off(kh)+h1+(h-1)*n_h,Sig,Pp)
                                    if (c > 0) V%B(Jb,Pb)%M(Sig)%v(c,nu) = V%B(Jb,Pb)%M(Sig)%v(c,nu) + Q * Sum
                                end do
                            end do
                        end do
                    end do
                    deallocate(Fx, Gx)
                end do
            end do
        end do

    contains

        pure real(dp) function amp(t,p,h,i)
            integer, intent(in) :: t, p, h, i
            if (t == 1) then
                amp = Ph%Cp(p,h,i)
            else
                amp = Ph%Cn(p,h,i)
            end if
        end function amp

    end subroutine EMPM_interaction

end module nj_empm_dens
