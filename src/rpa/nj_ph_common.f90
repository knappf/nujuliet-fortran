! Common 1p-1h building blocks of TDA & RPA (NuJuliet dev) ... counterpart of src/RPA/RPA_spurious.jl (== TDA_spurious.jl),
    ! RPA_allocate.jl (RPA_allocate, RPA_allocate_TBCF; == TDA_allocate.jl with the B parts switched off), the shared
    ! RPA_rB / TDA_rB of RPA_transitions.jl / TDA_transitions.jl & the shared *_transitions_export of RPA_export.jl /
    ! TDA_export.jl ...
    ! Storage conventions: channel arrays X(P,0:J_max) (P = 1 positive, 2 negative parity) as Orb_1p1h%N(P,J); the
    ! dev iteration order JP_initialize(J_max) is J outer, P inner ...
module nj_ph_common
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg, only : mm, mtm, spurious_ortho
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn, O2b_temp, O2b_temp_initialize
    use nj_tr1b, only : tr1b_t
    use nj_obdm, only : rho_rN
    implicit none
    private

    public :: ph_occ, ph_metric, ph_spurious_ini, ph_spurious_basis
    public :: ph_allocate, ph_allocate_tbcf, ph_rB, ph_transitions_export
    public :: ph_operators, ph_summary_header, ph_dimensions_block, ph_plot_spectrum, lstr, pm_str, JP_of

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
    end interface

contains

    ! Julia string(::Bool) ...
    function lstr(L) result(s)
        logical, intent(in) :: L
        character(len=:), allocatable :: s
        s = merge('true ', 'false', L)
        s = trim(s)
    end function lstr

    pure function pm_str(P) result(s)
        integer, intent(in) :: P
        character(len=1) :: s
        s = merge('+', '-', P == 1)
    end function pm_str

    ! i_JP -> (J,P) in the dev order JP_initialize (J outer, P inner) ...
    pure subroutine JP_of(i_JP,J,P)
        integer, intent(in) :: i_JP
        integer, intent(out) :: J, P
        J = (i_JP - 1) / 2
        P = mod(i_JP - 1, 2) + 1
    end subroutine JP_of

    ! Operator lists of the multipoles ... E0 [E0], E1 [E1,E1_VC,E1_VS,E1_TC,E1_TS,E1_sTC,E1_sTS,E1_C,E1_sC], E2 [E2], E3 [E3] ...
        ! (built element by element) ...
    function ph_operators(M,Multipole) result(Ops)
        type(tr1b_t), intent(in) :: M
        integer, intent(in) :: Multipole
        type(o1b_t), allocatable :: Ops(:)
        select case (Multipole)
        case (0)
            allocate(Ops(1))
            Ops(1) = M%E0
        case (1)
            allocate(Ops(9))
            Ops(1) = M%E1
            Ops(2) = M%E1_VC
            Ops(3) = M%E1_VS
            Ops(4) = M%E1_TC
            Ops(5) = M%E1_TS
            Ops(6) = M%E1_sTC
            Ops(7) = M%E1_sTS
            Ops(8) = M%E1_C
            Ops(9) = M%E1_sC
        case (2)
            allocate(Ops(1))
            Ops(1) = M%E2
        case default
            allocate(Ops(1))
            Ops(1) = M%E3
        end select
    end function ph_operators

    ! Occupation of orbital a with isospin T (-1 proton, +1 neutron) ...
    pure real(dp) function ph_occ(Orb,a,T)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: a, T
        if (T == -1) then
            ph_occ = Orb(a)%pO
        else
            ph_occ = Orb(a)%nO
        end if
    end function ph_occ

    ! Metric N = diag(D), D_ph = n_h - n_p of the 1p-1h channel (J,P) ...
    function ph_metric(Orb,Orb_1p1h,P,J) result(D)
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        integer, intent(in) :: P, J
        real(dp), allocatable :: D(:)
        integer :: i_ph
        type(ind2b_t) :: ph
        allocate(D(Orb_1p1h%N(P,J)))
        do i_ph = 1, Orb_1p1h%N(P,J)
            ph = Orb_1p1h%Ind(P,J)%v(i_ph)
            D(i_ph) = ph_occ(Orb,ph%b,ph%T) - ph_occ(Orb,ph%a,ph%T)
        end do
    end function ph_metric

    ! Spurious center-of-mass 1- vector ... counterpart of RPA_spurious_ini (3-argument method) == TDA_spurious_ini ...
        ! Spur_ph = (n_h - n_p) <p||E1||h>, |Spur_ph| < 1e-14 -> 0, normalized to 1 ...
    function ph_spurious_ini(Orb,Orb_1p1h,TrOp) result(Spur)
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(tr1b_t), intent(in) :: TrOp
        real(dp), allocatable :: Spur(:)
        integer :: N_ph, i_ph, a_p, a_h
        real(dp) :: D_ph

        N_ph = Orb_1p1h%N(2,1)
        allocate(Spur(N_ph))
        do i_ph = 1, N_ph
            a_p = Orb_1p1h%Ind(2,1)%v(i_ph)%a
            a_h = Orb_1p1h%Ind(2,1)%v(i_ph)%b
            if (Orb_1p1h%Ind(2,1)%v(i_ph)%T == 1) then
                D_ph = Orb(a_h)%nO - Orb(a_p)%nO
                Spur(i_ph) = D_ph * TrOp%E1%n(a_p,a_h)
            else
                D_ph = Orb(a_h)%pO - Orb(a_p)%pO
                Spur(i_ph) = D_ph * TrOp%E1%p(a_p,a_h)
            end if
        end do
        do i_ph = 1, N_ph
            if (abs(Spur(i_ph)) < 1.0e-14_dp) Spur(i_ph) = 0.0_dp
        end do
        Spur = Spur / norm2(Spur)
    end function ph_spurious_ini

    ! The closed-form Gram-Schmidt TDA_spurious_ortho! == RPA_spurious_ortho! is nj_linalg::spurious_ortho (core, stage 5:
        ! shared with the EMPM) ...

    ! U from spurious_ortho with the spurious column removed ... U_sub = U(:,1:N_ph-1) (dev shifts the columns
        ! Spur_ind+1..N_ph one to the left) ...
    function ph_spurious_basis(Spur) result(U_sub)
        real(dp), intent(in) :: Spur(:)
        real(dp), allocatable :: U_sub(:,:), U(:,:)
        integer :: N_ph, Spur_ind, a
        N_ph = size(Spur)
        allocate(U(N_ph,N_ph))
        call spurious_ortho(U,Spur,Spur_ind)
        do a = Spur_ind + 1, N_ph
            U(:,a-1) = U(:,a)
        end do
        U_sub = U(:,1:N_ph-1)
    end function ph_spurious_basis

    ! RPA matrices A, B & the metric N ... counterpart of RPA_allocate (RPA_allocate.jl:1-185; Want_B = .false. gives
        ! TDA_allocate, TDA_allocate.jl:1-137) without the TBCF part (ph_allocate_tbcf) ... lower triangle i_qg <= i_ph,
        ! mirrored ... parallel over i_ph (row i_ph & column i_ph of the lower/upper triangle are written by one task) ...
    subroutine ph_allocate(Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Want_B,A,B,N)
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: h_N
        type(o2b_t), intent(in) :: V_NN
        logical, intent(in) :: Want_B
        type(rmat_t), allocatable, intent(out) :: A(:,:), B(:,:), N(:,:)
        integer :: J_max, J, P, N_ph, i_ph, i_qg, a_p, a_h, T_ph, a_q, a_g, T_qg, l_p, j_p, l_h, j_h, l_q, j_q, l_g, j_g
        integer :: P_pg, P_pq, I
        real(dp) :: n_p, n_h, n_q, n_g, dn_ph, dn_qg, d_pq, d_hg, d_pg, d_hq, ASum, BSum, NSum, I_hat, Amp
        real(dp) :: h_pq, h_hg, h_pg, h_hq

        J_max = Orb_1p1h%J_max
        allocate(A(2,0:J_max), N(2,0:J_max))
        if (Want_B) allocate(B(2,0:J_max))

        do J = 0, J_max
            do P = 1, 2
                N_ph = Orb_1p1h%N(P,J)
                allocate(A(P,J)%v(N_ph,N_ph), N(P,J)%v(N_ph,N_ph))
                A(P,J)%v = 0.0_dp
                N(P,J)%v = 0.0_dp
                if (Want_B) then
                    allocate(B(P,J)%v(N_ph,N_ph))
                    B(P,J)%v = 0.0_dp
                end if

                !$omp parallel do schedule(dynamic) default(shared) &
                !$omp private(i_qg,a_p,a_h,T_ph,a_q,a_g,T_qg,l_p,j_p,l_h,j_h,l_q,j_q,l_g,j_g,P_pg,P_pq,I) &
                !$omp private(n_p,n_h,n_q,n_g,dn_ph,dn_qg,d_pq,d_hg,d_pg,d_hq,ASum,BSum,NSum,I_hat,Amp,h_pq,h_hg,h_pg,h_hq)
                do i_ph = 1, N_ph
                    a_p = Orb_1p1h%Ind(P,J)%v(i_ph)%a
                    a_h = Orb_1p1h%Ind(P,J)%v(i_ph)%b
                    T_ph = Orb_1p1h%Ind(P,J)%v(i_ph)%T
                    l_p = Orb(a_p)%l
                    j_p = Orb(a_p)%j
                    l_h = Orb(a_h)%l
                    j_h = Orb(a_h)%j
                    n_p = ph_occ(Orb,a_p,T_ph)
                    n_h = ph_occ(Orb,a_h,T_ph)
                    dn_ph = n_h - n_p

                    do i_qg = 1, i_ph
                        a_q = Orb_1p1h%Ind(P,J)%v(i_qg)%a
                        a_g = Orb_1p1h%Ind(P,J)%v(i_qg)%b
                        T_qg = Orb_1p1h%Ind(P,J)%v(i_qg)%T
                        l_q = Orb(a_q)%l
                        j_q = Orb(a_q)%j
                        l_g = Orb(a_g)%l
                        j_g = Orb(a_g)%j
                        n_q = ph_occ(Orb,a_q,T_qg)
                        n_g = ph_occ(Orb,a_g,T_qg)
                        dn_qg = n_g - n_q

                        d_pq = real(kronecker_delta(a_p,a_q),dp)
                        d_hg = real(kronecker_delta(a_h,a_g),dp)
                        d_pg = real(kronecker_delta(a_p,a_g),dp)
                        d_hq = real(kronecker_delta(a_h,a_q),dp)

                        P_pg = mod(l_p + l_g,2) + 1
                        P_pq = mod(l_p + l_q,2) + 1

                        ASum = 0.0_dp
                        BSum = 0.0_dp
                        NSum = 0.0_dp

                        ! 1-body mean-field term h & the overlap matrix N ...
                        if (T_ph == T_qg) then
                            if (T_ph == -1) then
                                h_pq = h_N%p(a_p,a_q)
                                h_hg = h_N%p(a_h,a_g)
                                h_pg = h_N%p(a_p,a_g)
                                h_hq = h_N%p(a_h,a_q)
                            else
                                h_pq = h_N%n(a_p,a_q)
                                h_hg = h_N%n(a_h,a_g)
                                h_pg = h_N%n(a_p,a_g)
                                h_hq = h_N%n(a_h,a_q)
                            end if
                            ASum = ASum + 0.5_dp * (dn_ph + dn_qg) * (h_pq * d_hg - h_hg * d_pq)
                            BSum = BSum + 0.5_dp * phase(J + (j_p + j_h)/2) * (dn_ph - dn_qg) * (h_pg * d_hq - h_hq * d_pg)
                            NSum = NSum + d_pq * d_hg * dn_ph
                        end if

                        ! 2-body scattering TDA matrix A ...
                        if (mod(l_p + l_g,2) == mod(l_q + l_h,2)) then
                            do I = max(abs(j_p - j_g),abs(j_q - j_h))/2, min(j_p + j_g,j_q + j_h)/2
                                I_hat = real(2*I + 1,dp)
                                if (T_ph == -1 .and. T_qg == -1) then
                                    Amp = dn_ph * dn_qg * phase((j_h + j_q)/2 + I) * I_hat * f6j(j_p,j_h,2*J,j_q,j_g,2*I)
                                    ASum = ASum + Amp * O2b_pp(a_p,a_g,a_h,a_q,I,P_pg,V_NN,Orb,Orb_NN)
                                else if (T_ph == -1 .and. T_qg == 1) then
                                    Amp = dn_ph * dn_qg * phase((j_h + j_q)/2 + I) * I_hat * f6j(j_p,j_h,2*J,j_q,j_g,2*I)
                                    ASum = ASum + Amp * O2b_pn(a_p,a_g,a_h,a_q,I,P_pg,V_NN,Orb_NN)
                                else if (T_ph == 1 .and. T_qg == -1) then
                                    Amp = dn_ph * dn_qg * phase((j_g + j_p)/2 + I) * I_hat * f6j(j_q,j_g,2*J,j_p,j_h,2*I)
                                    ASum = ASum + Amp * O2b_pn(a_q,a_h,a_g,a_p,I,P_pg,V_NN,Orb_NN)
                                else
                                    Amp = dn_ph * dn_qg * phase((j_h + j_q)/2 + I) * I_hat * f6j(j_p,j_h,2*J,j_q,j_g,2*I)
                                    ASum = ASum + Amp * O2b_nn(a_p,a_g,a_h,a_q,I,P_pg,V_NN,Orb,Orb_NN)
                                end if
                            end do
                        end if

                        ! 2-body correlation matrix B ...
                        if (Want_B) then
                            if (mod(l_p + l_q,2) == mod(l_h + l_g,2)) then
                                do I = max(abs(j_p - j_q),abs(j_h - j_g))/2, min(j_p + j_q,j_h + j_g)/2
                                    I_hat = real(2*I + 1,dp)
                                    if (T_ph == -1 .and. T_qg == -1) then
                                        Amp = dn_ph * dn_qg * phase((j_h + j_q)/2 + J + I) * I_hat * f6j(j_p,j_h,2*J,j_g,j_q,2*I)
                                        BSum = BSum + Amp * O2b_pp(a_p,a_q,a_h,a_g,I,P_pq,V_NN,Orb,Orb_NN)
                                    else if (T_ph == -1 .and. T_qg == 1) then
                                        Amp = dn_ph * dn_qg * phase((j_h + j_q)/2 + J + I) * I_hat * f6j(j_p,j_h,2*J,j_g,j_q,2*I)
                                        BSum = BSum + Amp * O2b_pn(a_p,a_q,a_h,a_g,I,P_pq,V_NN,Orb_NN)
                                    else if (T_ph == 1 .and. T_qg == -1) then
                                        Amp = dn_ph * dn_qg * phase((j_g + j_p)/2 + J + I) * I_hat * f6j(j_q,j_g,2*J,j_h,j_p,2*I)
                                        BSum = BSum + Amp * O2b_pn(a_q,a_p,a_g,a_h,I,P_pq,V_NN,Orb_NN)
                                    else
                                        Amp = dn_ph * dn_qg * phase((j_h + j_q)/2 + J + I) * I_hat * f6j(j_p,j_h,2*J,j_g,j_q,2*I)
                                        BSum = BSum + Amp * O2b_nn(a_p,a_q,a_h,a_g,I,P_pq,V_NN,Orb,Orb_NN)
                                    end if
                                end do
                            end if
                        end if

                        A(P,J)%v(i_ph,i_qg) = ASum
                        N(P,J)%v(i_ph,i_qg) = NSum
                        if (Want_B) B(P,J)%v(i_ph,i_qg) = BSum
                        if (i_ph /= i_qg) then
                            A(P,J)%v(i_qg,i_ph) = ASum
                            N(P,J)%v(i_qg,i_ph) = NSum
                            if (Want_B) B(P,J)%v(i_qg,i_ph) = BSum
                        end if
                    end do
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine ph_allocate

    ! C = A1 B1' + A2 B2' (+ A3 B3' + A4 B4') ... the dgemm sequence of RPA_allocate_TBCF (mul! with beta = 1) ...
    subroutine gemm_nt_sum(C,A1,B1,A2,B2,A3,B3,A4,B4)
        real(dp), intent(inout) :: C(:,:)
        real(dp), intent(in) :: A1(:,:), B1(:,:), A2(:,:), B2(:,:)
        real(dp), intent(in), optional :: A3(:,:), B3(:,:), A4(:,:), B4(:,:)
        integer :: n
        n = size(C,1)
        if (n == 0) return
        call dgemm('N','T',n,n,n,1.0_dp,A1,n,B1,n,0.0_dp,C,n)
        call dgemm('N','T',n,n,n,1.0_dp,A2,n,B2,n,1.0_dp,C,n)
        if (present(A3)) call dgemm('N','T',n,n,n,1.0_dp,A3,n,B3,n,1.0_dp,C,n)
        if (present(A4)) call dgemm('N','T',n,n,n,1.0_dp,A4,n,B4,n,1.0_dp,C,n)
    end subroutine gemm_nt_sum

    ! 2-body correlation function (box diagram) contributions to A (& B for Want_B) ... counterpart of RPA_allocate_TBCF
        ! (RPA_allocate.jl:187-566; Want_B = .false. & DEV = .false. gives TDA_allocate_TBCF, TDA_allocate.jl:139-446) ...
        ! DEV adds the 1-body polarization contact terms pf/nf (RPA.DEV only) ...
    subroutine ph_allocate_tbcf(Params,Orb,Orb_1p1h,Orb_NN,V_NN,Sigma,Want_B,DEV,A,B,Kind)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: V_NN, Sigma
        logical, intent(in) :: Want_B, DEV
        type(rmat_t), intent(inout) :: A(:,0:)
        type(rmat_t), intent(inout) :: B(:,0:)
        character(len=*), intent(in) :: Kind
        type(o2b_temp_t) :: W1_NN, W2_NN
        type(orb2b_temp_t) :: Orb_NN_t
        real(dp), allocatable :: ppXV(:,:), ppXS(:,:), pnXV(:,:), pnXS(:,:), nnXV(:,:), nnXS(:,:), npXV(:,:), npXS(:,:)
        real(dp), allocatable :: mpnXV(:,:), mpnXS(:,:), A_JP(:,:), B_JP(:,:), pf(:,:), nf(:,:)
        integer :: J_t_max, J, P, Nt, Bra, Col, a_k, a_l, a_b, a_d, l_k, j_k, l_l, j_l, l_b, j_b, l_d, j_d, P_kb, J_1, a_max
        integer :: N_ph, i_ph, i_qg, a_p, a_h, T_ph, a_q, a_g, T_qg, l_p, j_p, l_h, j_h, l_q, j_q, l_g, j_g, P_pg, P_pq, I
        integer :: k_ord, l_ord
        real(dp) :: Phase_kl, Amp_J1, Amp, ppV, ppS, pnV, pnS, nnV, nnS, npV, npS, mpnV, mpnS, ASum, BSum, Amp_A, Amp_B
        real(dp) :: pfSum, nfSum

        a_max = size(Orb)
        if (Want_B) then
            write(*,'(a)') TAB // TAB // 'Computing 2-body correlation function contributions to the ' // Kind // ' matrices A & B ...'
        else
            write(*,'(a)') TAB // TAB // 'Computing 2-body correlation function contributions to the ' // Kind // ' matrix A ...'
        end if

        ! Precompute box diagram 1 ... cross-coupled X[Bra,Col] = sum_J_1 (2J_1 + 1) {j_k j_l J; j_d j_b J_1} <kb;J_1|O|ld;J_1> ...
        call O2b_temp_initialize(Params,Orb,W1_NN,Orb_NN_t)
        J_t_max = ubound(Orb_NN_t%N,2)
        do J = 0, J_t_max
            do P = 1, 2
                Nt = Orb_NN_t%N(P,J)
                allocate(ppXV(Nt,Nt), ppXS(Nt,Nt), pnXV(Nt,Nt), pnXS(Nt,Nt), nnXV(Nt,Nt), nnXS(Nt,Nt), npXV(Nt,Nt), npXS(Nt,Nt), &
                         mpnXV(Nt,Nt), mpnXS(Nt,Nt))
                !$omp parallel do schedule(dynamic) default(shared) &
                !$omp private(Col,a_k,a_l,a_b,a_d,l_k,j_k,l_l,j_l,l_b,j_b,l_d,j_d,P_kb,J_1,Phase_kl,Amp_J1,Amp) &
                !$omp private(ppV,ppS,pnV,pnS,nnV,nnS,npV,npS,mpnV,mpnS)
                do Bra = 1, Nt
                    a_k = Orb_NN_t%Ind(P,J)%v(1,Bra)
                    a_l = Orb_NN_t%Ind(P,J)%v(2,Bra)
                    l_k = Orb(a_k)%l
                    j_k = Orb(a_k)%j
                    l_l = Orb(a_l)%l
                    j_l = Orb(a_l)%j
                    Phase_kl = phase((j_k + j_l)/2)
                    do Col = 1, Nt
                        a_b = Orb_NN_t%Ind(P,J)%v(1,Col)
                        a_d = Orb_NN_t%Ind(P,J)%v(2,Col)
                        l_b = Orb(a_b)%l
                        j_b = Orb(a_b)%j
                        l_d = Orb(a_d)%l
                        j_d = Orb(a_d)%j
                        P_kb = parity_of(l_k + l_b)
                        ppV = 0.0_dp
                        ppS = 0.0_dp
                        pnV = 0.0_dp
                        pnS = 0.0_dp
                        nnV = 0.0_dp
                        nnS = 0.0_dp
                        npV = 0.0_dp
                        npS = 0.0_dp
                        mpnV = 0.0_dp
                        mpnS = 0.0_dp
                        do J_1 = max(abs(j_k - j_b),abs(j_l - j_d))/2, min(j_k + j_b,j_l + j_d)/2
                            Amp_J1 = hat_of(2*J_1) * f6j(j_k,j_l,2*J,j_d,j_b,2*J_1)
                            Amp = phase(J_1) * Amp_J1
                            ppV = ppV + Amp * O2b_pp(a_k,a_b,a_l,a_d,J_1,P_kb,V_NN,Orb,Orb_NN)
                            ppS = ppS + Amp * O2b_pp(a_k,a_b,a_l,a_d,J_1,P_kb,Sigma,Orb,Orb_NN)
                            pnV = pnV + Amp * O2b_pn(a_k,a_b,a_l,a_d,J_1,P_kb,V_NN,Orb_NN)
                            pnS = pnS + Amp * O2b_pn(a_k,a_b,a_l,a_d,J_1,P_kb,Sigma,Orb_NN)
                            nnV = nnV + Amp * O2b_nn(a_k,a_b,a_l,a_d,J_1,P_kb,V_NN,Orb,Orb_NN)
                            nnS = nnS + Amp * O2b_nn(a_k,a_b,a_l,a_d,J_1,P_kb,Sigma,Orb,Orb_NN)
                            npV = npV + Phase_kl * Amp * O2b_pn(a_b,a_k,a_d,a_l,J_1,P_kb,V_NN,Orb_NN)
                            npS = npS + Phase_kl * Amp * O2b_pn(a_b,a_k,a_d,a_l,J_1,P_kb,Sigma,Orb_NN)
                            mpnV = mpnV + Amp_J1 * O2b_pn(a_k,a_b,a_d,a_l,J_1,P_kb,V_NN,Orb_NN)
                            mpnS = mpnS + Amp_J1 * O2b_pn(a_k,a_b,a_d,a_l,J_1,P_kb,Sigma,Orb_NN)
                        end do
                        ppXV(Bra,Col) = ppV
                        ppXS(Bra,Col) = ppS
                        pnXV(Bra,Col) = pnV
                        pnXS(Bra,Col) = pnS
                        nnXV(Bra,Col) = nnV
                        nnXS(Bra,Col) = nnS
                        npXV(Bra,Col) = npV
                        npXS(Bra,Col) = npS
                        mpnXV(Bra,Col) = mpnV
                        mpnXS(Bra,Col) = mpnS
                    end do
                end do
                !$omp end parallel do
                ! W1[Bra,Ket] = sum_Col X_V[Bra,Col] X_S[Ket,Col] + X_S[Bra,Col] X_V[Ket,Col] ...
                call gemm_nt_sum(W1_NN%pp(P,J)%v,ppXV,ppXS,ppXS,ppXV,pnXV,pnXS,pnXS,pnXV)
                call gemm_nt_sum(W1_NN%pn(P,J)%v,mpnXV,mpnXS,mpnXS,mpnXV)
                call gemm_nt_sum(W1_NN%nn(P,J)%v,nnXV,nnXS,nnXS,nnXV,npXV,npXS,npXS,npXV)
                deallocate(ppXV, ppXS, pnXV, pnXS, nnXV, nnXS, npXV, npXS, mpnXV, mpnXS)
            end do
        end do

        ! Precompute box diagram 2 ... direct <kl;J|O|mn;J>, W2 = V S' + S V' ...
        call O2b_temp_initialize(Params,Orb,W2_NN)
        do J = 0, J_t_max
            do P = 1, 2
                Nt = Orb_NN_t%N(P,J)
                allocate(ppXV(Nt,Nt), ppXS(Nt,Nt), pnXV(Nt,Nt), pnXS(Nt,Nt), nnXV(Nt,Nt), nnXS(Nt,Nt))
                !$omp parallel do schedule(dynamic) default(shared) private(Col,a_k,a_l,a_b,a_d)
                do Bra = 1, Nt
                    a_k = Orb_NN_t%Ind(P,J)%v(1,Bra)
                    a_l = Orb_NN_t%Ind(P,J)%v(2,Bra)
                    do Col = 1, Nt
                        a_b = Orb_NN_t%Ind(P,J)%v(1,Col)
                        a_d = Orb_NN_t%Ind(P,J)%v(2,Col)
                        ppXV(Bra,Col) = O2b_pp(a_k,a_l,a_b,a_d,J,P,V_NN,Orb,Orb_NN)
                        ppXS(Bra,Col) = O2b_pp(a_k,a_l,a_b,a_d,J,P,Sigma,Orb,Orb_NN)
                        pnXV(Bra,Col) = O2b_pn(a_k,a_l,a_b,a_d,J,P,V_NN,Orb_NN)
                        pnXS(Bra,Col) = O2b_pn(a_k,a_l,a_b,a_d,J,P,Sigma,Orb_NN)
                        nnXV(Bra,Col) = O2b_nn(a_k,a_l,a_b,a_d,J,P,V_NN,Orb,Orb_NN)
                        nnXS(Bra,Col) = O2b_nn(a_k,a_l,a_b,a_d,J,P,Sigma,Orb,Orb_NN)
                    end do
                end do
                !$omp end parallel do
                call gemm_nt_sum(W2_NN%pp(P,J)%v,ppXV,ppXS,ppXS,ppXV)
                call gemm_nt_sum(W2_NN%pn(P,J)%v,pnXV,pnXS,pnXS,pnXV)
                call gemm_nt_sum(W2_NN%nn(P,J)%v,nnXV,nnXS,nnXS,nnXV)
                deallocate(ppXV, ppXS, pnXV, pnXS, nnXV, nnXS)
            end do
        end do

        ! DEV only: 1-body polarization potential pf, nf (contact terms O(V sigma) of A, cf. RPA_allocate_SE) ...
        allocate(pf(a_max,a_max), nf(a_max,a_max))
        pf = 0.0_dp
        nf = 0.0_dp
        if (DEV) then
            !$omp parallel do schedule(dynamic) default(shared) &
            !$omp private(l_ord,a_k,a_l,l_k,j_k,l_l,j_l,a_b,l_b,j_b,P_kb,I,Amp,pfSum,nfSum)
            do k_ord = 1, a_max
                do l_ord = k_ord, a_max
                    a_k = Orb(k_ord)%a
                    a_l = Orb(l_ord)%a
                    l_k = Orb(k_ord)%l
                    j_k = Orb(k_ord)%j
                    l_l = Orb(l_ord)%l
                    j_l = Orb(l_ord)%j
                    if (j_k /= j_l) cycle
                    pfSum = 0.0_dp
                    nfSum = 0.0_dp
                    do a_b = 1, a_max
                        l_b = Orb(a_b)%l
                        j_b = Orb(a_b)%j
                        P_kb = parity_of(l_k + l_b)
                        if (P_kb /= parity_of(l_l + l_b)) cycle
                        do I = abs(j_k - j_b)/2, (j_k + j_b)/2
                            Amp = 0.5_dp * hat_of(2*I) / hat_of(j_k)
                            pfSum = pfSum + (0.5_dp * Amp * O2b_temp(a_k,a_b,a_l,a_b,I,P_kb,-1,W2_NN,Orb_NN_t) + &
                                             Amp * O2b_temp(a_k,a_b,a_l,a_b,I,P_kb,0,W2_NN,Orb_NN_t))
                            nfSum = nfSum + (0.5_dp * Amp * O2b_temp(a_k,a_b,a_l,a_b,I,P_kb,1,W2_NN,Orb_NN_t) + &
                                             phase(1 + (j_k + j_l)/2) * Amp * O2b_temp(a_b,a_k,a_b,a_l,I,P_kb,0,W2_NN,Orb_NN_t))
                        end do
                    end do
                    pf(a_k,a_l) = pfSum
                    nf(a_k,a_l) = nfSum
                    if (a_k /= a_l) then
                        pf(a_l,a_k) = pfSum
                        nf(a_l,a_k) = nfSum
                    end if
                end do
            end do
            !$omp end parallel do
        end if

        ! Box-diagram contributions per 1p-1h channel ... lower triangle, mirrored, parallel over i_ph ...
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                N_ph = Orb_1p1h%N(P,J)
                allocate(A_JP(N_ph,N_ph), B_JP(N_ph,N_ph))
                A_JP = 0.0_dp
                B_JP = 0.0_dp
                !$omp parallel do schedule(dynamic) default(shared) &
                !$omp private(i_qg,a_p,a_h,T_ph,a_q,a_g,T_qg,l_p,j_p,l_h,j_h,l_q,j_q,l_g,j_g,P_pg,P_pq,I,Amp,ASum,BSum) &
                !$omp private(Amp_A,Amp_B)
                do i_ph = 1, N_ph
                    a_p = Orb_1p1h%Ind(P,J)%v(i_ph)%a
                    a_h = Orb_1p1h%Ind(P,J)%v(i_ph)%b
                    T_ph = Orb_1p1h%Ind(P,J)%v(i_ph)%T
                    l_p = Orb(a_p)%l
                    j_p = Orb(a_p)%j
                    l_h = Orb(a_h)%l
                    j_h = Orb(a_h)%j
                    do i_qg = 1, i_ph
                        a_q = Orb_1p1h%Ind(P,J)%v(i_qg)%a
                        a_g = Orb_1p1h%Ind(P,J)%v(i_qg)%b
                        T_qg = Orb_1p1h%Ind(P,J)%v(i_qg)%T
                        l_q = Orb(a_q)%l
                        j_q = Orb(a_q)%j
                        l_g = Orb(a_g)%l
                        j_g = Orb(a_g)%j
                        P_pg = parity_of(l_p + l_g)
                        P_pq = parity_of(l_p + l_q)
                        ASum = 0.0_dp
                        BSum = 0.0_dp

                        if (DEV .and. T_ph == T_qg) then
                            ! 1-body polarization diagrams of A ...
                            Amp_A = phase(1)
                            if (T_ph == -1) then
                                if (a_p == a_q) ASum = ASum + Amp_A * pf(a_h,a_g)
                                if (a_h == a_g) ASum = ASum + Amp_A * pf(a_p,a_q)
                            else
                                if (a_p == a_q) ASum = ASum + Amp_A * nf(a_h,a_g)
                                if (a_h == a_g) ASum = ASum + Amp_A * nf(a_p,a_q)
                            end if
                            ! ... & of B ...
                            Amp_B = phase(1 - J + (j_p + j_h)/2)
                            if (T_ph == -1) then
                                if (a_p == a_g) BSum = BSum + Amp_B * pf(a_h,a_q)
                                if (a_h == a_q) BSum = BSum + Amp_B * pf(a_p,a_g)
                            else
                                if (a_p == a_g) BSum = BSum + Amp_B * nf(a_h,a_q)
                                if (a_h == a_q) BSum = BSum + Amp_B * nf(a_p,a_g)
                            end if
                        end if

                        ! A ... box diagram 1 ...
                        do I = max(abs(j_p - j_q),abs(j_h - j_g))/2, min(j_p + j_q,j_h + j_g)/2
                            if (T_ph == -1 .and. T_qg == -1) then
                                Amp = phase(2 + J + I + (j_h + j_g)/2) * hat_of(2*I) * f6j(j_p,j_q,2*I,j_g,j_h,2*J)
                                ASum = ASum + Amp * O2b_temp(a_p,a_q,a_h,a_g,I,P_pq,-1,W1_NN,Orb_NN_t)
                            else if (T_ph == -1 .and. T_qg == 1) then
                                Amp = phase(2 + J + I + (j_h + j_q)/2) * hat_of(2*I) * f6j(j_p,j_q,2*I,j_g,j_h,2*J)
                                ASum = ASum + Amp * O2b_temp(a_p,a_q,a_h,a_g,I,P_pq,0,W1_NN,Orb_NN_t)
                            else if (T_ph == 1 .and. T_qg == -1) then
                                Amp = phase(2 + J + I + (j_g + j_p)/2) * hat_of(2*I) * f6j(j_q,j_p,2*I,j_h,j_g,2*J)
                                ASum = ASum + Amp * O2b_temp(a_q,a_p,a_g,a_h,I,P_pq,0,W1_NN,Orb_NN_t)
                            else
                                Amp = phase(2 + J + I + (j_h + j_g)/2) * hat_of(2*I) * f6j(j_p,j_q,2*I,j_g,j_h,2*J)
                                ASum = ASum + Amp * O2b_temp(a_p,a_q,a_h,a_g,I,P_pq,1,W1_NN,Orb_NN_t)
                            end if
                        end do
                        ! A ... box diagram 2 ...
                        do I = max(abs(j_p - j_g),abs(j_q - j_h))/2, min(j_p + j_g,j_q + j_h)/2
                            Amp = phase(1 + I + (j_h + j_q)/2) * hat_of(2*I) * f6j(j_p,j_h,2*J,j_q,j_g,2*I)
                            if (T_ph == -1 .and. T_qg == -1) then
                                ASum = ASum + 0.5_dp * Amp * O2b_temp(a_p,a_g,a_h,a_q,I,P_pg,-1,W2_NN,Orb_NN_t)
                            else if (T_ph == -1 .and. T_qg == 1) then
                                ASum = ASum + Amp * O2b_temp(a_p,a_g,a_h,a_q,I,P_pg,0,W2_NN,Orb_NN_t)
                            else if (T_ph == 1 .and. T_qg == -1) then
                                ASum = ASum + phase((j_p + j_h + j_q + j_g)/2) * Amp * O2b_temp(a_g,a_p,a_q,a_h,I,P_pg,0,W2_NN,Orb_NN_t)
                            else
                                ASum = ASum + 0.5_dp * Amp * O2b_temp(a_p,a_g,a_h,a_q,I,P_pg,1,W2_NN,Orb_NN_t)
                            end if
                        end do

                        if (Want_B) then
                            ! B ... box diagram 1 ...
                            do I = max(abs(j_p - j_g),abs(j_h - j_q))/2, min(j_p + j_g,j_h + j_q)/2
                                if (T_ph == -1 .and. T_qg == -1) then
                                    Amp = phase(1 + I + (j_h + j_g)/2) * hat_of(2*I) * f6j(j_p,j_g,2*I,j_q,j_h,2*J)
                                    BSum = BSum + Amp * O2b_temp(a_p,a_g,a_h,a_q,I,P_pg,-1,W1_NN,Orb_NN_t)
                                else if (T_ph == -1 .and. T_qg == 1) then
                                    Amp = phase(1 + I + (j_h + j_q)/2) * hat_of(2*I) * f6j(j_p,j_g,2*I,j_q,j_h,2*J)
                                    BSum = BSum + Amp * O2b_temp(a_p,a_g,a_h,a_q,I,P_pg,0,W1_NN,Orb_NN_t)
                                else if (T_ph == 1 .and. T_qg == -1) then
                                    Amp = phase(1 + I + (j_g + j_p)/2) * hat_of(2*I) * f6j(j_q,j_h,2*I,j_p,j_g,2*J)
                                    BSum = BSum + Amp * O2b_temp(a_q,a_h,a_g,a_p,I,P_pg,0,W1_NN,Orb_NN_t)
                                else
                                    Amp = phase(1 + I + (j_h + j_g)/2) * hat_of(2*I) * f6j(j_p,j_g,2*I,j_q,j_h,2*J)
                                    BSum = BSum + Amp * O2b_temp(a_p,a_g,a_h,a_q,I,P_pg,1,W1_NN,Orb_NN_t)
                                end if
                            end do
                            ! B ... box diagram 2 ...
                            do I = max(abs(j_p - j_q),abs(j_h - j_g))/2, min(j_p + j_q,j_h + j_g)/2
                                Amp = phase(J + I + (j_h + j_q)/2) * hat_of(2*I) * f6j(j_p,j_h,2*J,j_g,j_q,2*I)
                                if (T_ph == -1 .and. T_qg == -1) then
                                    BSum = BSum + 0.5_dp * Amp * O2b_temp(a_p,a_q,a_h,a_g,I,P_pq,-1,W2_NN,Orb_NN_t)
                                else if (T_ph == -1 .and. T_qg == 1) then
                                    BSum = BSum + Amp * O2b_temp(a_p,a_q,a_h,a_g,I,P_pq,0,W2_NN,Orb_NN_t)
                                else if (T_ph == 1 .and. T_qg == -1) then
                                    BSum = BSum + phase((j_p + j_h + j_q + j_g)/2) * Amp * O2b_temp(a_q,a_p,a_g,a_h,I,P_pq,0,W2_NN,Orb_NN_t)
                                else
                                    BSum = BSum + 0.5_dp * Amp * O2b_temp(a_p,a_q,a_h,a_g,I,P_pq,1,W2_NN,Orb_NN_t)
                                end if
                            end do
                        end if

                        A_JP(i_ph,i_qg) = ASum
                        B_JP(i_ph,i_qg) = BSum
                        if (i_ph /= i_qg) then
                            A_JP(i_qg,i_ph) = ASum
                            B_JP(i_qg,i_ph) = BSum
                        end if
                    end do
                end do
                !$omp end parallel do
                A(P,J)%v = A(P,J)%v + A_JP
                if (Want_B) B(P,J)%v = B(P,J)%v + B_JP
                deallocate(A_JP, B_JP)
            end do
        end do

        if (Want_B) then
            write(*,'(a)') TAB // TAB // TAB // 'Finished on evaluation of 2-body correlation function contributions to the ' // &
                           Kind // ' matrices A & B ...'
        else
            write(*,'(a)') TAB // TAB // TAB // 'Finished on evaluation of 2-body correlation function contributions to the ' // &
                           Kind // ' matrix A ...'
        end if
    end subroutine ph_allocate_tbcf

    ! Reduced transition intensities B (physical, isoscalar & isovector) ... counterpart of RPA_rB == TDA_rB
        ! (RPA_transitions.jl:240-430) incl. the isovector E1 fixes of f3f7a33 ... E1 = energies of the (J,P) = (1,-) channel ...
    function ph_rB(Params,Orb_1p1h,rM,E1) result(rB)
        type(params_t), intent(in) :: Params
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(reduced_multipole_t), intent(in) :: rM
        real(dp), intent(in) :: E1(:)
        type(reduced_transition_t) :: rB
        real(dp), parameter :: hc = 197.326980_dp, g_p = 5.586_dp, g_n = -3.826_dp
        real(dp) :: e_p, e_n, E
        complex(dp) :: E1p, E1n, VCp, VCn, VSp, VSn, TCp, TCn, TSp, TSn, sTCp, sTCn, sTSp, sTSn, Cp, Cn, sCp, sCn
        integer :: N_nu, nu

        write(*,'(/,a)') 'Calculating the reduced transition intensities rB ...'
        call simple(rM%E0, Orb_1p1h%N(1,0), rB%E0)
        call simple(rM%E2, Orb_1p1h%N(1,2), rB%E2)
        call simple(rM%E3, Orb_1p1h%N(2,3), rB%E3)

        N_nu = Orb_1p1h%N(2,1)
        call alloc3(rB%E1,N_nu)
        call alloc3(rB%E1_V,N_nu)
        call alloc3(rB%E1_VC,N_nu)
        call alloc3(rB%E1_VS,N_nu)
        call alloc3(rB%E1_T,N_nu)
        call alloc3(rB%E1_TC,N_nu)
        call alloc3(rB%E1_TS,N_nu)
        call alloc3(rB%E1_C,N_nu)
        call alloc3(rB%E1_CC,N_nu)
        call alloc3(rB%E1_CS,N_nu)
        call alloc3(rB%E1_comp,N_nu)
        call alloc3(rB%E1_NLO_LWA,N_nu)

        e_p = real(Params%Calc%A - Params%Calc%Z,dp) / real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp) / real(Params%Calc%A,dp)

        do nu = 1, N_nu
            E = E1(nu)
            E1p = rM%E1%p(nu)
            E1n = rM%E1%n(nu)
            VCp = rM%E1_VC%p(nu)
            VCn = rM%E1_VC%n(nu)
            VSp = rM%E1_VS%p(nu)
            VSn = rM%E1_VS%n(nu)
            TCp = rM%E1_TC%p(nu)
            TCn = rM%E1_TC%n(nu)
            TSp = rM%E1_TS%p(nu)
            TSn = rM%E1_TS%n(nu)
            sTCp = rM%E1_sTC%p(nu)
            sTCn = rM%E1_sTC%n(nu)
            sTSp = rM%E1_sTS%p(nu)
            sTSn = rM%E1_sTS%n(nu)
            Cp = rM%E1_C%p(nu)
            Cn = rM%E1_C%n(nu)
            sCp = rM%E1_sC%p(nu)
            sCn = rM%E1_sC%n(nu)

            ! Physical components ...
            rB%E1%ph(nu) = abs(E1p)**2
            rB%E1_V%ph(nu) = abs(VCp + 0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn)**2
            rB%E1_VC%ph(nu) = abs(VCp)**2
            rB%E1_VS%ph(nu) = abs(0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn)**2
            rB%E1_T%ph(nu) = abs(TCp + 0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn)**2
            rB%E1_TC%ph(nu) = abs(TCp)**2
            rB%E1_TS%ph(nu) = abs(0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn)**2
            rB%E1_C%ph(nu) = abs(TCp + 0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn - (VCp + 0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn))**2
            rB%E1_CC%ph(nu) = abs(TCp - VCp)**2
            rB%E1_CS%ph(nu) = abs(0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn - (0.5_dp * g_p * VSp + 0.5_dp * g_n * VSn))**2
            rB%E1_comp%ph(nu) = abs(Cp)**2
            rB%E1_NLO_LWA%ph(nu) = abs(E1p + E / hc * (TCp + 0.5_dp * g_p * TSp + 0.5_dp * g_n * TSn))**2

            ! Isoscalar components ...
            rB%E1%is(nu) = 0.25_dp * abs(E1p + E1n)**2
            rB%E1_V%is(nu) = abs(0.5_dp * (VCp + VCn) + 0.125_dp * (g_p + g_n) * (VSp + VSn))**2
            rB%E1_VC%is(nu) = abs(0.5_dp * (VCp + VCn))**2
            rB%E1_VS%is(nu) = abs(0.125_dp * (g_p + g_n) * (VSp + VSn))**2
            rB%E1_T%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn) + 0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn))**2
            rB%E1_TC%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn))**2
            rB%E1_TS%is(nu) = abs(0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn))**2
            rB%E1_C%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn) + 0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn) - &
                             (0.5_dp * (VCp + VCn) + 0.125_dp * (g_p + g_n) * (VSp + VSn)))**2
            rB%E1_CC%is(nu) = abs(0.5_dp * (TCp + sTCp + TCn + sTCn) - (0.5_dp * (VCp + VCn)))**2
            rB%E1_CS%is(nu) = abs(0.125_dp * (g_p + g_n) * (TSp + sTSp + TSn + sTSn) - (0.125_dp * (g_p + g_n) * (VSp + VSn)))**2
            rB%E1_comp%is(nu) = 0.25_dp * abs(Cp + sCp + Cn + sCn)**2
            rB%E1_NLO_LWA%is(nu) = 0.25_dp * abs(E1p + E1n + E / hc * ((TCp + sTCp + TCn + sTCn) + 0.25_dp * (g_p + g_n) * &
                                   (TSp + sTSp + TSn + sTSn)))**2

            ! Isovector components ... (f3f7a33: VS.p - VS.n, TS.p - TS.n & 0.25 (g_p - g_n) in the NLO term)
            rB%E1%iv(nu) = abs(e_p * E1p - e_n * E1n)**2
            rB%E1_V%iv(nu) = abs(0.5_dp * (VCp - VCn) + 0.125_dp * (g_p - g_n) * (VSp - VSn))**2
            rB%E1_VC%iv(nu) = abs(0.5_dp * (VCp - VCn))**2
            rB%E1_VS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (VSp - VSn))**2
            rB%E1_T%iv(nu) = abs(0.5_dp * (TCp - TCn) + 0.125_dp * (g_p - g_n) * (TSp - TSn))**2
            rB%E1_TC%iv(nu) = abs(0.5_dp * (TCp - TCn))**2
            rB%E1_TS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (TSp - TSn))**2
            rB%E1_C%iv(nu) = abs(0.5_dp * (TCp - TCn) + 0.125_dp * (g_p - g_n) * (TSp - TSn) - &
                             (0.5_dp * (VCp - VCn) + 0.125_dp * (g_p - g_n) * (VSp - VSn)))**2
            rB%E1_CC%iv(nu) = abs(0.5_dp * (TCp - TCn) - (0.5_dp * (VCp - VCn)))**2
            rB%E1_CS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (TSp - TSn) - (0.125_dp * (g_p - g_n) * (VSp - VSn)))**2
            rB%E1_comp%iv(nu) = 0.25_dp * abs(Cp - Cn)**2
            rB%E1_NLO_LWA%iv(nu) = abs(e_p * E1p - e_n * E1n + 0.5_dp * E / hc * ((TCp - TCn) + 0.25_dp * (g_p - g_n) * &
                                   (TSp - TSn)))**2
        end do

        write(*,'(a)') TAB // 'Reduced transition intensities rB evaluated ...'

    contains

        subroutine alloc3(T,N)
            type(transition_t), intent(out) :: T
            integer, intent(in) :: N
            allocate(T%ph(N), T%is(N), T%iv(N))
        end subroutine alloc3

        subroutine simple(M,N,T)
            type(pncvector_t), intent(in) :: M
            integer, intent(in) :: N
            type(transition_t), intent(out) :: T
            integer :: k
            allocate(T%ph(N), T%is(N), T%iv(N))
            do k = 1, N
                T%ph(k) = abs(M%p(k))**2
                T%is(k) = 0.25_dp * abs(M%p(k) + M%n(k))**2
                T%iv(k) = 0.25_dp * abs(M%p(k) - M%n(k))**2
            end do
        end subroutine simple

    end function ph_rB

    ! Title, nuclide, calculation & CMS lines of the TDA/RPA summary files ...
    subroutine ph_summary_header(u,Params,Title)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Title
        integer :: N_max
        N_max = Params%Calc%Nmax
        write(u,'(a)') Title
        write(u,'(/,a)') 'Nuclide data:    A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)
        write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                         ', s.p. J-basis size = ' // istr((N_max+1)*(N_max+2)/2) // ', s.p. M-basis size = ' // &
                         istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Combined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Only pure 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else
            write(u,'(/,a)') TAB // TAB // TAB // 'No Center-of-Mass System (CMS) motion correction is set ...'
        end if
    end subroutine ph_summary_header

    ! Dimension tables of the 1p-1h subspaces ...
    subroutine ph_dimensions_block(u,Orb_1p1h,Kind)
        integer, intent(in) :: u
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        character(len=*), intent(in) :: Kind
        integer :: J, N_M_tot
        write(u,'(/,a)') 'Information on dimensions of 1-phonon 1p-1h ' // Kind // ' subspaces ...'
        N_M_tot = 0
        do J = 0, Orb_1p1h%J_max
            N_M_tot = N_M_tot + (Orb_1p1h%N(1,J) + Orb_1p1h%N(2,J)) * (2*J + 1)
        end do
        write(u,'(/,a)') 'Total number of 1 phonon states in J-scheme = ' // istr(sum(Orb_1p1h%N))
        write(u,'(a)') 'Total number of 1 phonon states in M-scheme = ' // istr(N_M_tot)
        write(u,'(/,a)') 'J' // TAB // 'N_phonon'
        do J = 0, Orb_1p1h%J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_1p1h%N(1,J) + Orb_1p1h%N(2,J))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = +'
        do J = 0, Orb_1p1h%J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_1p1h%N(1,J))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = -'
        do J = 0, Orb_1p1h%J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_1p1h%N(2,J))
        end do
    end subroutine ph_dimensions_block

    ! Plot-ready spectrum (*_plot.dat) ... ground state + all N levels sorted by energy (all N+1 rows, dev) ...
    subroutine ph_plot_spectrum(Path,Orb_1p1h,E)
        character(len=*), intent(in) :: Path
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E(:,0:)
        real(dp), allocatable :: Sol(:,:), E_full(:)
        integer, allocatable :: Sort(:)
        real(dp), parameter :: Delta = 1.0_dp
        integer :: N_ph, nu_count, J, P, nu, u

        N_ph = sum(Orb_1p1h%N)
        allocate(Sol(4,N_ph+1), Sort(N_ph+1))
        Sol(:,1) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
        nu_count = 1
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                do nu = 1, Orb_1p1h%N(P,J)
                    nu_count = nu_count + 1
                    Sol(:,nu_count) = [real(J,dp), real(P,dp), E(P,J)%v(nu), E(P,J)%v(nu)]
                end do
            end do
        end do
        E_full = Sol(3,:)
        call sortperm(E_full,Sort)
        Sol = Sol(:,Sort)
        do nu = 2, N_ph + 1
            if (Sol(4,nu) < Sol(4,nu-1)) then
                Sol(4,nu) = Sol(4,nu-1) + Delta
            else if (abs(Sol(4,nu) - Sol(4,nu-1)) < Delta) then
                Sol(4,nu) = Sol(4,nu-1) + Delta
            end if
        end do

        open(newunit=u, file=Path, status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // &
                       cfmt_s('label',12,.true.)
        do nu = 1, N_ph + 1
            write(u,'(a)') cfmt_d(round_even(Sol(1,nu)),5,.true.) // ' ' // &
                           cfmt_s(merge('+', '-', abs(Sol(2,nu) - 1.0_dp) < 1.0e-3_dp),5,.true.) // ' ' // &
                           cfmt_f(Sol(3,nu),12,6,.true.) // ' ' // cfmt_f(Sol(4,nu),12,6,.true.)
        end do
        close(u)
    end subroutine ph_plot_spectrum

    ! Electromagnetic moments m_k, TRK sum rules & the transition tables ... counterpart of RPA_transitions_export
        ! (RPA_export.jl:525-993) & TDA_transitions_export (TDA_export.jl:360-819) ... Rho_ref_basis is the OBDM in the
        ! TDA/RPA (reference) basis, E the real parts of the energies, Out_Dir = io/<Path>/<Name> ...
    subroutine ph_transitions_export(Params,Kind,Out_Dir,Ortho,Orb,Orb_1p1h,C,Rho_ref_basis,E,rB)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind, Out_Dir
        logical, intent(in) :: Ortho
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(o1b_t), intent(in) :: C, Rho_ref_basis
        type(rvec_t), intent(in) :: E(:,0:)
        type(reduced_transition_t), intent(in) :: rB
        real(dp), parameter :: hc = 197.326980_dp, pmc2 = 938.272013_dp
        real(dp) :: pR2, pR4, Z, TRK(0:3), m(5,3,0:3)
        character(len=6), parameter :: Units(0:3) = ['fm^4  ', 'fm^2  ', 'fm^4  ', 'fm^6  ']
        character(len=10), parameter :: Kinds(3) = ['physical  ', 'isoscalar ', 'isovector ']
        integer :: u, iso, X
        character(len=:), allocatable :: D

        Z = real(Params%Calc%Z,dp)
        pR2 = rho_rN(Params,2,Rho_ref_basis,C,Orb)
        pR4 = rho_rN(Params,4,Rho_ref_basis,C,Orb)

        write(*,'(/,a)') 'Preparing export of information on ' // Kind // ' 1-phonon electromagnetic transitions ...'

        TRK(0) = hc**2 / pmc2 / (2.0_dp * pi) * Z * pR2
        TRK(1) = hc**2 / pmc2 * 9.0_dp / (8.0_dp * pi) * Z
        TRK(2) = hc**2 / pmc2 * 50.0_dp / (8.0_dp * pi) * Z * pR2
        TRK(3) = hc**2 / pmc2 * 147.0_dp / (8.0_dp * pi) * Z * pR4

        call moments(rB%E0,E(1,0)%v,.false.,m(:,:,0))
        call moments(rB%E1,E(2,1)%v,Ortho,m(:,:,1))
        call moments(rB%E2,E(1,2)%v,.false.,m(:,:,2))
        call moments(rB%E3,E(2,3)%v,.false.,m(:,:,3))

        open(newunit=u, file=Out_Dir // '/' // Kind // '_summary.dat', position='append', action='write')
        write(u,'(/,a)') 'Review of ' // Kind // ' electromagnetic transition moments m_k ...'
        do iso = 1, 3
            do X = 0, 3
                write(u,'(/,a)') TAB // 'E' // istr(X) // ' ' // trim(Kinds(iso)) // ':'
                write(u,'(a)') 'm_-1  = ' // cfmt_f(m(1,iso,X),20,8) // TAB // 'e^2 ' // trim(Units(X)) // ' / MeV'
                write(u,'(a)') 'm_0   = ' // cfmt_f(m(2,iso,X),20,8) // TAB // 'e^2 ' // trim(Units(X))
                write(u,'(a)') 'm_1   = ' // cfmt_f(m(3,iso,X),20,8) // TAB // 'e^2 ' // trim(Units(X)) // ' MeV'
                write(u,'(a)') 'm_2   = ' // cfmt_f(m(4,iso,X),20,8) // TAB // 'e^2 ' // trim(Units(X)) // ' MeV^2'
                write(u,'(a)') 'm_3   = ' // cfmt_f(m(5,iso,X),20,8) // TAB // 'e^2 ' // trim(Units(X)) // ' MeV^3'
            end do
        end do
        if (Kind == 'TDA') then
            write(u,'(/,a,/)') 'Physical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values (reference values, not ' // &
                               'fulfilled exactly by TDA) ...'
        else
            write(u,'(/,a,/)') 'Review of ' // Kind // ' physical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values ...'
        end if
        do X = 0, 3
            write(u,'(a)') 'S_E' // istr(X) // ' = ' // cfmt_f(TRK(X),20,8) // TAB // 'e^2 ' // trim(Units(X)) // ' MeV'
        end do
        close(u)

        write(*,'(/,a)') 'Preparing export of ' // Kind // ' transition intensities ...'
        D = Out_Dir // '/transitions'
        call table(D // '/E0/' // Kind // '_E0.dat',E(1,0)%v,rB%E0)
        call table(D // '/E1/' // Kind // '_E1.dat',E(2,1)%v,rB%E1)
        call table(D // '/E1/' // Kind // '_E1_V.dat',E(2,1)%v,rB%E1_V)
        call table(D // '/E1/' // Kind // '_E1_V_conv.dat',E(2,1)%v,rB%E1_VC)
        call table(D // '/E1/' // Kind // '_E1_V_spin.dat',E(2,1)%v,rB%E1_VS)
        call table(D // '/E1/' // Kind // '_E1_T.dat',E(2,1)%v,rB%E1_T)
        call table(D // '/E1/' // Kind // '_E1_T_conv.dat',E(2,1)%v,rB%E1_TC)
        call table(D // '/E1/' // Kind // '_E1_T_spin.dat',E(2,1)%v,rB%E1_TS)
        call table(D // '/E1/' // Kind // '_E1_C.dat',E(2,1)%v,rB%E1_C)
        call table(D // '/E1/' // Kind // '_E1_C_conv.dat',E(2,1)%v,rB%E1_CC)
        call table(D // '/E1/' // Kind // '_E1_C_spin.dat',E(2,1)%v,rB%E1_CS)
        call table(D // '/E1/' // Kind // '_E1_c.dat',E(2,1)%v,rB%E1_comp)
        call table(D // '/E1/' // Kind // '_E1_NLO_LWA.dat',E(2,1)%v,rB%E1_NLO_LWA)
        call table(D // '/E2/' // Kind // '_E2.dat',E(1,2)%v,rB%E2)
        call table(D // '/E3/' // Kind // '_E3.dat',E(2,3)%v,rB%E3)
        write(*,'(a)') TAB // 'Export of information on ' // Kind // ' 1-phonon electromagnetic transitions completed ...'

    contains

        ! m(1:5) = (m_-1, m_0, m_1, m_2, m_3) for ph, is & iv ... Skip_First skips the spurious 1- level (Ortho) ...
        subroutine moments(T,Ev,Skip_First,mm_out)
            type(transition_t), intent(in) :: T
            real(dp), intent(in) :: Ev(:)
            logical, intent(in) :: Skip_First
            real(dp), intent(out) :: mm_out(5,3)
            real(dp) :: B(3), En
            integer :: nu
            mm_out = 0.0_dp
            do nu = 1, size(Ev)
                if (Skip_First .and. nu == 1) cycle
                En = Ev(nu)
                B = [T%ph(nu), T%is(nu), T%iv(nu)]
                if (abs(En) > 1.0e-3_dp) mm_out(1,:) = mm_out(1,:) + B / En
                mm_out(2,:) = mm_out(2,:) + B
                mm_out(3,:) = mm_out(3,:) + En * B
                mm_out(4,:) = mm_out(4,:) + En**2 * B
                mm_out(5,:) = mm_out(5,:) + En**3 * B
            end do
        end subroutine moments

        subroutine table(Path,Ev,T)
            character(len=*), intent(in) :: Path
            real(dp), intent(in) :: Ev(:)
            type(transition_t), intent(in) :: T
            integer :: uu, nu
            open(newunit=uu, file=Path, status='replace', action='write')
            write(uu,'(a)') cfmt_s('E',20,.true.) // ' ' // cfmt_s('B_ph',20,.true.) // ' ' // cfmt_s('B_is',20,.true.) // ' ' // &
                            cfmt_s('B_iv',20,.true.)
            do nu = 1, size(Ev)
                write(uu,'(a)') cfmt_f(Ev(nu),20,8,.true.) // ' ' // cfmt_f(T%ph(nu),20,8,.true.) // ' ' // &
                                cfmt_f(T%is(nu),20,8,.true.) // ' ' // cfmt_f(T%iv(nu),20,8,.true.)
            end do
            close(uu)
        end subroutine table

    end subroutine ph_transitions_export

end module nj_ph_common
