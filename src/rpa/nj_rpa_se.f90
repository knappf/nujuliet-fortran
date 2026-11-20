! Experimental DEV path of the RPA ... counterpart of RPA_allocate_SE (src/RPA/RPA_allocate.jl:1232-1868, NuJuliet dev) ...
    ! Energy-dependent 1p-1h <-> 2p-2h self-energy coupling of the matrix A, reached only with RPA.DEV = .true. & the
    ! TBCF allocation (Reference_TBCF with the reference sigma.bin), from RPA_allocate & from RPA_energy (dressing energy),
    ! both with Mode = "OnShell" & the defaults Regularization = "DSRG", Gamma = 2 MeV, dE_bin = Gamma/20, Chunk = 32768 ...
    !
    !   A(E) = A^(s) + Sigma(E) ,  Sigma(E)_ij = - sum_c A_12[i,c] A_12[j,c] g_c(E) - sum_c B_12[i,c] B_12[j,c] gb_c(E) ,
    !   g_c(E)  = w_lo,c f(Delta_b(c) - E) + w_hi,c f(Delta_b(c)+1 - E) ,  gb_c(E) = the same with Delta_b + E ,
    !   w_lo,c = (1 - t_c) / N_c , w_hi,c = t_c / N_c  (linear splitting of the pole Delta_c between the two neighbouring bins
    !   of the dev grid Delta_b = Delta_min + (b - 1) dE_bin, f(x) = (1 - exp(-(x/Gamma)^2)) / x) ,
    !   OnShell ... dA_ik = 1/2 (Sigma(w_i)_ik + Sigma(w_k)_ik) symmetrized, w_i = e_pp - e_hh of the 1p-1h state i ...
    !
    ! Deviation (algorithmic, same result up to round-off, docs/PORTING_NOTES.md): dev accumulates the spectral weights
    ! W[i,j,b] (N_ph^2 x N_bin per channel) & evaluates a full Sigma(w_i) for every row i (N_ph^3 N_bin); here the bin
    ! weights are contracted with f(Delta_b - w_i) per configuration first, dA = 1/2 (S + S^T) with
    ! S = (A_12 o G) A_12^T, G[i,c] = -g_c(w_i) ... one dgemm per block of configurations, no N_bin dimension ...
    ! The order of the degenerate 2p-2h configurations (Julia Dict iteration + stable sort) only changes the round-off.
    ! Not ported (unreachable in dev: both callers pass Mode = "OnShell" & the default regularization): the modes "QP"
    ! (SE_QP, root search) & "Fixed", the regularizations "Imaginary" & "None", RPA_allocate_SE_Old ...
module nj_rpa_se
    use, intrinsic :: iso_c_binding, only : c_double
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    implicit none
    private

    public :: RPA_allocate_SE

    interface
        pure function c_expm1(x) bind(C, name='expm1') result(y)
            import :: c_double
            real(c_double), value, intent(in) :: x
            real(c_double) :: y
        end function c_expm1
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
    end interface

    ! 2p-2h configuration ... particle & hole orbitals as indices into the Par & Hol lists, couplings, pole & 1/N_c ...
    type :: conf_t
        integer :: p_1 = 0, p_2 = 0, h_1 = 0, h_2 = 0, J12 = 0, J34 = 0
        real(dp) :: Delta = 0.0_dp, N_inv = 0.0_dp
    end type conf_t

    ! Orbital o = (a,T) ...
    type :: so_t
        integer :: a = 0, T = 0
    end type so_t

    ! Rows of a (J,P) channel keyed by an orbital (a,T) ... CSR lists of (row i, partner orbital) ...
    type :: rows_t
        integer, allocatable :: Start(:,:), Row(:), Partner_a(:)
    end type rows_t

contains

    subroutine RPA_allocate_SE(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,A)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: h_N
        type(o2b_t), intent(in) :: V_NN
        type(rmat_t), intent(inout) :: A(:,0:)
        real(dp), parameter :: Gamma = 2.0_dp, dE_bin = Gamma / 20.0_dp, Tol_ord = 1.0e-8_dp, Tol_N = 1.0e-8_dp
        integer, parameter :: Chunk = 32768
        type(so_t), allocatable :: Par(:), Hol(:)
        integer, allocatable :: PP_i1(:), PP_i2(:), PP_K(:), PP_Start(:,:,:), HH_i1(:), HH_i2(:), HH_K(:), HH_Start(:,:,:)
        type(conf_t), allocatable :: Conf(:), Conf_s(:)
        integer, allocatable :: Perm(:), Bin(:)
        real(dp), allocatable :: Delta_bin(:), W_lo(:), W_hi(:), w(:), A12(:,:), B12(:,:), G(:,:), S(:,:), Key(:)
        type(rows_t) :: Rows_H, Rows_P
        logical :: Overlap
        real(dp) :: Delta_min, Delta_max, e_min_p, e_max_p, e_min_h, e_max_h, x, t
        integer :: a_max, N_par, N_hol, N_bin, J, P, N_ph, N_conf, N_conf_max, K_max, i, k, c, c_0, nc, ia, b

        a_max = size(Orb)
        write(*,'(a)') TAB // TAB // 'Computing 1p-1h to 2p-2h self-energy coupling of 1p-1h A matrix ...'
        write(*,'(a)') TAB // TAB // TAB // 'Mode = OnShell, Regularization = DSRG, Gamma = ' // jstr(Gamma) // ' MeV, dE_bin = ' // &
                       jstr(dE_bin) // ' MeV'

        ! Particle & hole orbitals (same criteria as the 1p-1h space; the sets may overlap) ...
        allocate(Par(2*a_max), Hol(2*a_max))
        N_par = 0
        N_hol = 0
        do ia = 1, a_max
            if (Orb(ia)%pO < Params%Calc%RPA%Tol_par) call push(Par,N_par,Orb(ia)%a,-1)
            if (Orb(ia)%pO > Params%Calc%RPA%Tol_hol) call push(Hol,N_hol,Orb(ia)%a,-1)
            if (Orb(ia)%nO < Params%Calc%RPA%Tol_par) call push(Par,N_par,Orb(ia)%a,1)
            if (Orb(ia)%nO > Params%Calc%RPA%Tol_hol) call push(Hol,N_hol,Orb(ia)%a,1)
        end do
        if (N_par == 0 .or. N_hol == 0) call fatal('RPA_allocate_SE: empty particle or hole orbital set ...')

        Overlap = .false.
        do i = 1, N_par
            do k = 1, N_hol
                if (Par(i)%a == Hol(k)%a .and. Par(i)%T == Hol(k)%T) Overlap = .true.
            end do
        end do
        if (Overlap) write(*,'(a)') TAB // TAB // TAB // 'Overlapping particle & hole orbitals ... backward 2p-2h coupling ' // &
                                    'B_12 G(-E) B_21 is included in A ...'

        ! Grid of the 2p-2h pole energies ...
        e_min_p = huge(1.0_dp)
        e_max_p = -huge(1.0_dp)
        do i = 1, N_par
            e_min_p = min(e_min_p,e_so(Par(i),Par(i)))
            e_max_p = max(e_max_p,e_so(Par(i),Par(i)))
        end do
        e_min_h = huge(1.0_dp)
        e_max_h = -huge(1.0_dp)
        do i = 1, N_hol
            e_min_h = min(e_min_h,e_so(Hol(i),Hol(i)))
            e_max_h = max(e_max_h,e_so(Hol(i),Hol(i)))
        end do
        Delta_min = 2.0_dp * e_min_p - 2.0_dp * e_max_h - 2.0_dp * dE_bin
        Delta_max = 2.0_dp * e_max_p - 2.0_dp * e_min_h + 2.0_dp * dE_bin
        N_bin = ceiling((Delta_max - Delta_min) / dE_bin) + 2
        allocate(Delta_bin(N_bin))
        do b = 1, N_bin
            Delta_bin(b) = Delta_min + real(b - 1,dp) * dE_bin
        end do

        ! pp & hh pairs (o_1 <= o_2 in list order) bucketed by (T_1 + T_2, l_1 + l_2 mod 2, K) in insertion order ...
        K_max = maxval(Orb(:)%j)
        call pairs_make(Par,N_par,PP_i1,PP_i2,PP_K,PP_Start)
        call pairs_make(Hol,N_hol,HH_i1,HH_i2,HH_K,HH_Start)

        N_conf_max = 0
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                N_ph = Orb_1p1h%N(P,J)
                if (N_ph == 0) cycle

                call rows_make(Rows_H,Rows_P)
                call configurations_make()
                N_conf = size(Conf)
                N_conf_max = max(N_conf_max,N_conf)

                ! Sorted by the pole energy (stable) ... bins & weights ...
                allocate(Key(N_conf), Perm(N_conf))
                do c = 1, N_conf
                    Key(c) = Conf(c)%Delta
                end do
                call sortperm(Key,Perm)
                Conf_s = Conf(Perm)
                deallocate(Key, Perm)
                allocate(Bin(N_conf), W_lo(N_conf), W_hi(N_conf))
                do c = 1, N_conf
                    x = (Conf_s(c)%Delta - Delta_min) / dE_bin
                    Bin(c) = floor(x) + 1
                    t = x - real(Bin(c) - 1,dp)
                    W_lo(c) = (1.0_dp - t) * Conf_s(c)%N_inv
                    W_hi(c) = t * Conf_s(c)%N_inv
                end do

                ! On-shell energies of the 1p-1h states ...
                allocate(w(N_ph), S(N_ph,N_ph))
                do i = 1, N_ph
                    w(i) = e_ab(Orb_1p1h%Ind(P,J)%v(i)%a,Orb_1p1h%Ind(P,J)%v(i)%T,Orb_1p1h%Ind(P,J)%v(i)%a,Orb_1p1h%Ind(P,J)%v(i)%T) - &
                           e_ab(Orb_1p1h%Ind(P,J)%v(i)%b,Orb_1p1h%Ind(P,J)%v(i)%T,Orb_1p1h%Ind(P,J)%v(i)%b,Orb_1p1h%Ind(P,J)%v(i)%T)
                end do
                S = 0.0_dp

                ! Blocks of Chunk configurations ... S += (A_12 o G) A_12^T (+ the backward part) ...
                do c_0 = 1, N_conf, Chunk
                    nc = min(Chunk, N_conf - c_0 + 1)
                    allocate(A12(N_ph,nc), G(N_ph,nc))
                    if (Overlap) allocate(B12(N_ph,nc))
                    !$omp parallel do schedule(dynamic) default(shared) private(i)
                    do c = 1, nc
                        A12(:,c) = 0.0_dp
                        call A12_column(A12(:,c),Rows_H,Rows_P,Conf_s(c_0 + c - 1),J,.false.)
                        if (Overlap) then
                            B12(:,c) = 0.0_dp
                            call A12_column(B12(:,c),Rows_P,Rows_H,Conf_s(c_0 + c - 1),J,.true.)
                        end if
                        do i = 1, N_ph
                            G(i,c) = -A12(i,c) * (W_lo(c_0 + c - 1) * f_reg(Delta_bin(Bin(c_0 + c - 1)) - w(i)) + &
                                                  W_hi(c_0 + c - 1) * f_reg(Delta_bin(Bin(c_0 + c - 1) + 1) - w(i)))
                        end do
                    end do
                    !$omp end parallel do
                    call dgemm('N','T',N_ph,N_ph,nc,1.0_dp,G,N_ph,A12,N_ph,1.0_dp,S,N_ph)
                    if (Overlap) then
                        !$omp parallel do schedule(static) default(shared) private(i)
                        do c = 1, nc
                            do i = 1, N_ph
                                G(i,c) = -B12(i,c) * (W_lo(c_0 + c - 1) * f_reg(Delta_bin(Bin(c_0 + c - 1)) + w(i)) + &
                                                      W_hi(c_0 + c - 1) * f_reg(Delta_bin(Bin(c_0 + c - 1) + 1) + w(i)))
                            end do
                        end do
                        !$omp end parallel do
                        call dgemm('N','T',N_ph,N_ph,nc,1.0_dp,G,N_ph,B12,N_ph,1.0_dp,S,N_ph)
                        deallocate(B12)
                    end if
                    deallocate(A12, G)
                end do

                ! A_eff = A^(s) + 1/2 (Sigma_on + Sigma_on^T) , Sigma_on(i,k) = 1/2 (S(i,k) + S(k,i)) ...
                A(P,J)%v = A(P,J)%v + 0.5_dp * (S + transpose(S))
                deallocate(Bin, W_lo, W_hi, w, S, Conf, Conf_s)
            end do
        end do

        write(*,'(a)') TAB // TAB // TAB // 'Largest channel ... ' // istr(N_conf_max) // ' J-coupled 2p-2h configurations, ' // &
                       istr(N_bin) // ' pole bins'
        write(*,'(a)') TAB // TAB // TAB // 'Finished on evaluation of 1p-1h to 2p-2h self-energy contributions to the matrix ...'

    contains

        subroutine push(L,N,a_,T_)
            type(so_t), intent(inout) :: L(:)
            integer, intent(inout) :: N
            integer, intent(in) :: a_, T_
            N = N + 1
            L(N) = so_t(a_,T_)
        end subroutine push

        ! Bucket index of (T_1 + T_2) ... -2, 0, 2 -> 1, 2, 3 ...
        pure integer function t_ind(T_sum)
            integer, intent(in) :: T_sum
            t_ind = T_sum/2 + 2
        end function t_ind

        subroutine pairs_make(L,N,I1,I2,KK,Start)
            type(so_t), intent(in) :: L(:)
            integer, intent(in) :: N
            integer, allocatable, intent(out) :: I1(:), I2(:), KK(:), Start(:,:,:)
            integer, allocatable :: Cnt(:,:,:), Pos(:,:,:)
            integer :: i_1, i_2, K, N_pair, it, il, n_tot
            allocate(Cnt(3,0:1,0:K_max))
            Cnt = 0
            do i_1 = 1, N
                do i_2 = i_1, N
                    do K = abs(Orb(L(i_1)%a)%j - Orb(L(i_2)%a)%j)/2, (Orb(L(i_1)%a)%j + Orb(L(i_2)%a)%j)/2
                        if (i_1 == i_2 .and. mod(K,2) == 1) cycle
                        it = t_ind(L(i_1)%T + L(i_2)%T)
                        il = mod(Orb(L(i_1)%a)%l + Orb(L(i_2)%a)%l,2)
                        Cnt(it,il,K) = Cnt(it,il,K) + 1
                    end do
                end do
            end do
            ! Start(it,il,K) .. Start of the next bucket - 1 (Start has one extra slot at the end) ...
            allocate(Start(3,0:1,0:K_max+1), Pos(3,0:1,0:K_max))
            n_tot = 0
            do K = 0, K_max
                do il = 0, 1
                    do it = 1, 3
                        Start(it,il,K) = n_tot + 1
                        Pos(it,il,K) = n_tot
                        n_tot = n_tot + Cnt(it,il,K)
                    end do
                end do
            end do
            Start(:,:,K_max+1) = n_tot + 1
            N_pair = n_tot
            allocate(I1(N_pair), I2(N_pair), KK(N_pair))
            do i_1 = 1, N
                do i_2 = i_1, N
                    do K = abs(Orb(L(i_1)%a)%j - Orb(L(i_2)%a)%j)/2, (Orb(L(i_1)%a)%j + Orb(L(i_2)%a)%j)/2
                        if (i_1 == i_2 .and. mod(K,2) == 1) cycle
                        it = t_ind(L(i_1)%T + L(i_2)%T)
                        il = mod(Orb(L(i_1)%a)%l + Orb(L(i_2)%a)%l,2)
                        Pos(it,il,K) = Pos(it,il,K) + 1
                        I1(Pos(it,il,K)) = i_1
                        I2(Pos(it,il,K)) = i_2
                        KK(Pos(it,il,K)) = K
                    end do
                end do
            end do
        end subroutine pairs_make

        ! Last entry of bucket (it,il,K) ... the next bucket in the storage order (it fastest, then il, then K) starts after it ...
        pure integer function bucket_end(Start,it,il,K)
            integer, intent(in) :: Start(:,0:,0:), it, il, K
            if (it < 3) then
                bucket_end = Start(it+1,il,K) - 1
            else if (il == 0) then
                bucket_end = Start(1,1,K) - 1
            else
                bucket_end = Start(1,0,K+1) - 1
            end if
        end function bucket_end

        ! J-coupled 2p-2h configurations of the channel (J,P) ... counterpart of configurations_make ...
        subroutine configurations_make()
            integer :: it, il_12, il_34, J12, J34, k_pp, k_hh, N_c, p_1, p_2, h_1, h_2
            real(dp) :: Nc
            type(conf_t), allocatable :: Buf(:)
            allocate(Buf(1024))
            N_c = 0
            do it = 1, 3
                do il_12 = 0, 1
                    do J12 = 0, K_max
                        if (bucket_end(PP_Start,it,il_12,J12) < PP_Start(it,il_12,J12)) cycle
                        do il_34 = 0, 1
                            if (parity_of(il_12 + il_34) /= P) cycle
                            do J34 = abs(J12 - J), min(J12 + J,K_max)
                                if (bucket_end(HH_Start,it,il_34,J34) < HH_Start(it,il_34,J34)) cycle
                                do k_pp = PP_Start(it,il_12,J12), bucket_end(PP_Start,it,il_12,J12)
                                    p_1 = PP_i1(k_pp)
                                    p_2 = PP_i2(k_pp)
                                    do k_hh = HH_Start(it,il_34,J34), bucket_end(HH_Start,it,il_34,J34)
                                        h_1 = HH_i1(k_hh)
                                        h_2 = HH_i2(k_hh)
                                        if ((min(n_so(Hol(h_1)),n_so(Hol(h_2))) - max(n_so(Par(p_1)),n_so(Par(p_2)))) <= Tol_ord) cycle
                                        Nc = n_so(Hol(h_1)) * n_so(Hol(h_2)) * (1.0_dp - n_so(Par(p_1))) * (1.0_dp - n_so(Par(p_2))) - &
                                             (1.0_dp - n_so(Hol(h_1))) * (1.0_dp - n_so(Hol(h_2))) * n_so(Par(p_1)) * n_so(Par(p_2))
                                        if (abs(Nc) < Tol_N) cycle
                                        N_c = N_c + 1
                                        if (N_c > size(Buf)) call grow(Buf)
                                        Buf(N_c) = conf_t(p_1, p_2, h_1, h_2, J12, J34, &
                                                          e_so(Par(p_1),Par(p_1)) + e_so(Par(p_2),Par(p_2)) - &
                                                          e_so(Hol(h_1),Hol(h_1)) - e_so(Hol(h_2),Hol(h_2)), 1.0_dp / Nc)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            Conf = Buf(1:N_c)
        end subroutine configurations_make

        subroutine grow(Buf)
            type(conf_t), allocatable, intent(inout) :: Buf(:)
            type(conf_t), allocatable :: Tmp(:)
            allocate(Tmp(2*size(Buf)))
            Tmp(1:size(Buf)) = Buf
            call move_alloc(Tmp,Buf)
        end subroutine grow

        ! Rows of the 1p-1h states of the channel by their hole (Rows_H: partner = particle) & particle (Rows_P: partner =
            ! hole) orbital ... the reversed lists of dev (Rows_H_rev, Rows_P_rev) are Rows_P & Rows_H ...
        subroutine rows_make(RH,RP)
            type(rows_t), intent(out) :: RH, RP
            call rows_csr(RH,.true.)
            call rows_csr(RP,.false.)
        end subroutine rows_make

        subroutine rows_csr(R,By_hole)
            type(rows_t), intent(out) :: R
            logical, intent(in) :: By_hole
            integer, allocatable :: Cnt(:,:), Pos(:,:)
            integer :: i_, T_i, k_, x_, t_, n_
            type(ind2b_t) :: ph
            allocate(Cnt(a_max,2), R%Start(a_max+1,2), Pos(a_max,2))
            Cnt = 0
            do i_ = 1, N_ph
                ph = Orb_1p1h%Ind(P,J)%v(i_)
                T_i = (ph%T + 3)/2
                x_ = merge(ph%b, ph%a, By_hole)
                Cnt(x_,T_i) = Cnt(x_,T_i) + 1
            end do
            n_ = 0
            do t_ = 1, 2
                do k_ = 1, a_max
                    R%Start(k_,t_) = n_ + 1
                    Pos(k_,t_) = n_
                    n_ = n_ + Cnt(k_,t_)
                end do
                R%Start(a_max+1,t_) = n_ + 1
            end do
            allocate(R%Row(n_), R%Partner_a(n_))
            do i_ = 1, N_ph
                ph = Orb_1p1h%Ind(P,J)%v(i_)
                T_i = (ph%T + 3)/2
                x_ = merge(ph%b, ph%a, By_hole)
                Pos(x_,T_i) = Pos(x_,T_i) + 1
                R%Row(Pos(x_,T_i)) = i_
                R%Partner_a(Pos(x_,T_i)) = merge(ph%a, ph%b, By_hole)
            end do
        end subroutine rows_csr

        ! Column A_12[:,c] (B_12 for Reversed with the swapped row lists) ... counterpart of A12_column! ...
        subroutine A12_column(A_c,RH,RP,Cf,J_,Reversed)
            real(dp), intent(inout) :: A_c(:)
            type(rows_t), intent(in) :: RH, RP
            type(conf_t), intent(in) :: Cf
            integer, intent(in) :: J_
            logical, intent(in) :: Reversed
            type(so_t) :: p1, p2, h1, h2, o
            real(dp) :: Norm, Amp
            integer :: k_, i_
            p1 = Par(Cf%p_1)
            p2 = Par(Cf%p_2)
            h1 = Hol(Cf%h_1)
            h2 = Hol(Cf%h_2)
            Norm = 1.0_dp / sqrt((1.0_dp + delta_so(p1,p2)) * (1.0_dp + delta_so(h1,h2)))

            ! Hole h2 is the spectator ... alpha^(ph)_(h1;p1p2) ...
            do k_ = RH%Start(h2%a,(h2%T + 3)/2), RH%Start(h2%a + 1,(h2%T + 3)/2) - 1
                i_ = RH%Row(k_)
                o = so_t(RH%Partner_a(k_),h2%T)
                Amp = Norm * merge(phase((j_so(o) - j_so(h2))/2), 1.0_dp, Reversed)
                A_c(i_) = A_c(i_) + Amp * C_alpha(j_so(o),j_so(h2),j_so(h1),J_,Cf%J12,Cf%J34) * alpha(h2,h1,o,p1,p2,Cf%J12)
            end do
            ! Hole h1 is the spectator ... - alpha^(ph)_(h2;p1p2) ...
            do k_ = RH%Start(h1%a,(h1%T + 3)/2), RH%Start(h1%a + 1,(h1%T + 3)/2) - 1
                i_ = RH%Row(k_)
                o = so_t(RH%Partner_a(k_),h1%T)
                Amp = - Norm * merge(phase((j_so(o) - j_so(h1))/2), 1.0_dp, Reversed) * phase((j_so(h1) + j_so(h2))/2 - Cf%J34)
                A_c(i_) = A_c(i_) + Amp * C_alpha(j_so(o),j_so(h1),j_so(h2),J_,Cf%J12,Cf%J34) * alpha(h1,h2,o,p1,p2,Cf%J12)
            end do
            ! Particle p2 is the spectator ... beta^(ph)_(p1;h1h2) ...
            do k_ = RP%Start(p2%a,(p2%T + 3)/2), RP%Start(p2%a + 1,(p2%T + 3)/2) - 1
                i_ = RP%Row(k_)
                o = so_t(RP%Partner_a(k_),p2%T)
                Amp = Norm * merge(phase((j_so(p2) - j_so(o))/2), 1.0_dp, Reversed)
                A_c(i_) = A_c(i_) + Amp * C_beta(j_so(p2),j_so(o),j_so(p1),J_,Cf%J12,Cf%J34) * beta(p2,p1,o,h1,h2,Cf%J34)
            end do
            ! Particle p1 is the spectator ... - beta^(ph)_(p2;h1h2) ...
            do k_ = RP%Start(p1%a,(p1%T + 3)/2), RP%Start(p1%a + 1,(p1%T + 3)/2) - 1
                i_ = RP%Row(k_)
                o = so_t(RP%Partner_a(k_),p1%T)
                Amp = - Norm * merge(phase((j_so(p1) - j_so(o))/2), 1.0_dp, Reversed) * phase((j_so(p1) + j_so(p2))/2 - Cf%J12)
                A_c(i_) = A_c(i_) + Amp * C_beta(j_so(p1),j_so(o),j_so(p2),J_,Cf%J12,Cf%J34) * beta(p1,p2,o,h1,h2,Cf%J34)
            end do
        end subroutine A12_column

        pure real(dp) function n_so(o)
            type(so_t), intent(in) :: o
            if (o%T == -1) then
                n_so = Orb(o%a)%pO
            else
                n_so = Orb(o%a)%nO
            end if
        end function n_so

        pure real(dp) function nb_so(o)
            type(so_t), intent(in) :: o
            nb_so = 1.0_dp - n_so(o)
        end function nb_so

        pure integer function j_so(o)
            type(so_t), intent(in) :: o
            j_so = Orb(o%a)%j
        end function j_so

        pure real(dp) function delta_so(o_1,o_2)
            type(so_t), intent(in) :: o_1, o_2
            delta_so = merge(1.0_dp, 0.0_dp, o_1%a == o_2%a .and. o_1%T == o_2%T)
        end function delta_so

        ! Generalized mean field e_ab = h_N[a,b] (same isospin, l & j) ...
        pure real(dp) function e_ab(a_,T_a,b_,T_b)
            integer, intent(in) :: a_, T_a, b_, T_b
            e_ab = 0.0_dp
            if (T_a /= T_b .or. Orb(a_)%l /= Orb(b_)%l .or. Orb(a_)%j /= Orb(b_)%j) return
            if (T_a == -1) then
                e_ab = h_N%p(a_,b_)
            else
                e_ab = h_N%n(a_,b_)
            end if
        end function e_ab

        pure real(dp) function e_so(o_1,o_2)
            type(so_t), intent(in) :: o_1, o_2
            e_so = e_ab(o_1%a,o_1%T,o_2%a,o_2%T)
        end function e_so

        ! J-coupled antisymmetrized <o1 o2;K|V|o3 o4;K> for all isospin combinations ... counterpart of the local V ...
        pure real(dp) function V_so(o1,o2,o3,o4,K) result(V)
            type(so_t), intent(in) :: o1, o2, o3, o4
            integer, intent(in) :: K
            V = 0.0_dp
            if ((o1%T + o2%T) /= (o3%T + o4%T)) return
            if (o1%T == o2%T) then
                if (.not. allowed(o1,o2,o3,o4,K)) return
                if (Orb_NN%Dic(max(o1%a,o2%a),min(o1%a,o2%a),K,1) == 0 .or. &
                    Orb_NN%Dic(max(o3%a,o4%a),min(o3%a,o4%a),K,1) == 0) return
                if (o1%T == -1) then
                    V = O2b_pp(o1%a,o2%a,o3%a,o4%a,K,parity_of(Orb(o1%a)%l + Orb(o2%a)%l),V_NN,Orb,Orb_NN)
                else
                    V = O2b_nn(o1%a,o2%a,o3%a,o4%a,K,parity_of(Orb(o1%a)%l + Orb(o2%a)%l),V_NN,Orb,Orb_NN)
                end if
            else if (o1%T == -1) then
                if (o3%T == -1) then
                    V = V_pn(o1,o2,o3,o4,K)
                else
                    ! <pn|V|np> = -(-1)^((j_c + j_d)/2 - K) <pn|V|pn> ...
                    V = - phase((j_so(o3) + j_so(o4))/2 - K) * V_pn(o1,o2,o4,o3,K)
                end if
            else
                ! <np|V|..> = -(-1)^((j_a + j_b)/2 - K) <pn|V|..> ...
                if (o3%T == -1) then
                    V = - phase((j_so(o1) + j_so(o2))/2 - K) * V_pn(o2,o1,o3,o4,K)
                else
                    V = - phase((j_so(o1) + j_so(o2))/2 - K) * (- phase((j_so(o3) + j_so(o4))/2 - K)) * V_pn(o2,o1,o4,o3,K)
                end if
            end if
        end function V_so

        ! Charge-conserving parity, triangle & J-range tests of the local V (the order of the checks does not matter) ...
        pure logical function allowed(o1,o2,o3,o4,K)
            type(so_t), intent(in) :: o1, o2, o3, o4
            integer, intent(in) :: K
            allowed = .false.
            if (parity_of(Orb(o1%a)%l + Orb(o2%a)%l) /= parity_of(Orb(o3%a)%l + Orb(o4%a)%l)) return
            if (.not. triangle(j_so(o1),j_so(o2),2*K) .or. .not. triangle(j_so(o3),j_so(o4),2*K)) return
            if (K > ubound(Orb_NN%Dic,3)) return
            allowed = .true.
        end function allowed

        ! <p n|V|p n> (o1, o3 protons) ...
        pure real(dp) function V_pn(o1,o2,o3,o4,K)
            type(so_t), intent(in) :: o1, o2, o3, o4
            integer, intent(in) :: K
            V_pn = 0.0_dp
            if (.not. allowed(o1,o2,o3,o4,K)) return
            if (Orb_NN%Dic(o1%a,o2%a,K,0) == 0 .or. Orb_NN%Dic(o3%a,o4%a,K,0) == 0) return
            V_pn = O2b_pn(o1%a,o2%a,o3%a,o4%a,K,parity_of(Orb(o1%a)%l + Orb(o2%a)%l),V_NN,Orb_NN)
        end function V_pn

        pure real(dp) function Fh(h,h_1,p,p_1,p_2)
            type(so_t), intent(in) :: h, h_1, p, p_1, p_2
            Fh = 0.5_dp * (n_so(h) - n_so(p)) * (n_so(h_1) * nb_so(p_1) * nb_so(p_2) + nb_so(h_1) * n_so(p_1) * n_so(p_2))
            Fh = Fh + 0.5_dp * (n_so(h) * n_so(h_1) * nb_so(p_1) * nb_so(p_2) - nb_so(h) * nb_so(h_1) * n_so(p_1) * n_so(p_2))
        end function Fh

        pure real(dp) function Fp(p,p_1,h,h_1,h_2)
            type(so_t), intent(in) :: p, p_1, h, h_1, h_2
            Fp = 0.5_dp * (n_so(h) - n_so(p)) * (nb_so(p_1) * n_so(h_1) * n_so(h_2) + n_so(p_1) * nb_so(h_1) * nb_so(h_2))
            Fp = Fp + 0.5_dp * (nb_so(p) * nb_so(p_1) * n_so(h_1) * n_so(h_2) - n_so(p) * n_so(p_1) * nb_so(h_1) * nb_so(h_2))
        end function Fp

        pure real(dp) function G_occ(h,h_1,p,p_1)
            type(so_t), intent(in) :: h, h_1, p, p_1
            G_occ = (n_so(h) - n_so(p)) * (n_so(h_1) - n_so(p_1))
        end function G_occ

        ! Vertex alpha^(ph)_(h1;p1p2) = <h1 p;K|alpha|p1 p2;K> (spectator hole h) ...
        pure real(dp) function alpha(h,h_1,p,p_1,p_2,K) result(ME)
            type(so_t), intent(in) :: h, h_1, p, p_1, p_2
            integer, intent(in) :: K
            ME = Fh(h,h_1,p,p_1,p_2) * V_so(h_1,p,p_1,p_2,K)
            ME = ME + 0.5_dp * delta_so(p,p_2) * e_so(h_1,p_1) * G_occ(h,h_1,p,p_1)
            ME = ME - 0.5_dp * phase((j_so(p) + j_so(p_2))/2 - K) * delta_so(p,p_1) * e_so(h_1,p_2) * G_occ(h,h_1,p,p_2)
        end function alpha

        ! Vertex beta^(ph)_(p1;h1h2) = <h p1;K|beta|h1 h2;K> (spectator particle p) ...
        pure real(dp) function beta(p,p_1,h,h_1,h_2,K)
            type(so_t), intent(in) :: p, p_1, h, h_1, h_2
            integer, intent(in) :: K
            beta = Fp(p,p_1,h,h_1,h_2) * V_so(h,p_1,h_1,h_2,K)
        end function beta

        ! Recoupling of the vertices to [(p1 p2)_J12 (h1 h2)_J34]_J ...
        real(dp) function C_alpha(j_p,j_h,j_x,J_,J12,J34)
            integer, intent(in) :: j_p, j_h, j_x, J_, J12, J34
            C_alpha = phase((j_p + j_h)/2 + (j_p + j_x)/2 + J_ + J12) * sqrt(hat_of(2*J12) * hat_of(2*J34)) * &
                      f6j(j_p,j_h,2*J_,2*J34,2*J12,j_x)
        end function C_alpha

        real(dp) function C_beta(j_p,j_h,j_x,J_,J12,J34)
            integer, intent(in) :: j_p, j_h, j_x, J_, J12, J34
            C_beta = phase((j_p + j_h)/2 + (j_p + j_x)/2 + (j_h + j_x)/2 + J_ + J12 + J34) * &
                     sqrt(hat_of(2*J12) * hat_of(2*J34)) * f6j(j_p,j_h,2*J_,2*J34,2*J12,j_x)
        end function C_beta

        ! DSRG-regularized f(x) ~ 1/x ...
        pure real(dp) function f_reg(x_)
            real(dp), intent(in) :: x_
            if (x_ == 0.0_dp) then
                f_reg = 0.0_dp
            else
                f_reg = -c_expm1(-(x_ / Gamma)**2) / x_
            end if
        end function f_reg

    end subroutine RPA_allocate_SE

end module nj_rpa_se
