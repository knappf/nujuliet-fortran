! Common parts of the QTDA & QRPA solvers ... QTDA_Allocate.jl / QRPA_Allocate.jl, *_Spurious.jl, *_Transitions.jl and
    ! the shared parts of *_Export.jl (which are identical in the original up to the labels) ...
module nj_qex
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_tr1b, only : qptr1b_t
    use nj_obdm, only : rho_rN
    use nj_o1b, only : O1b_import
    use nj_qpo2b
    implicit none
    private

    public :: qex_allocate, spurious_PN_initialize, spurious_CM_initialize, spurious_orthogonalize
    public :: qex_rM, qex_rB, qex_transitions_export, qex_spectrum_plot_export, qex_dimensions_summary
    public :: qex_file_exists, qex_reference_check

contains

    ! QTDA/QRPA matrix A (& optionally B) for all (J,P) channels up to J = 2*Nmax + 1 ...
    subroutine qex_allocate(Params,Orb,Orb_NN,Orb_2qp,H_N,H_NN,Amat,Bmat)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(qpo1b_t), intent(in) :: H_N
        type(qpo2b_t), intent(in) :: H_NN
        type(rmat_t), allocatable, intent(out) :: Amat(:,:)
        type(rmat_t), allocatable, intent(out), optional :: Bmat(:,:)
        integer :: J_max, J, P, N_qp, i_qp, j_qp, a, b, c, d, T_ab, T_cd, JP
        real(dp) :: ASum, BSum, Norm, E_a, E_b
        logical :: With_B

        With_B = present(Bmat)
        J_max = 2*Params%Calc%Nmax + 1
        allocate(Amat(2,0:J_max))
        if (With_B) allocate(Bmat(2,0:J_max))
        do J = 0, J_max
            do P = 1, 2
                allocate(Amat(P,J)%v(Orb_2qp%N(P,J),Orb_2qp%N(P,J)))
                if (With_B) allocate(Bmat(P,J)%v(Orb_2qp%N(P,J),Orb_2qp%N(P,J)))
            end do
        end do

        !$omp parallel do schedule(dynamic) private(J,P,N_qp,i_qp,j_qp,a,b,c,d,T_ab,T_cd,ASum,BSum,Norm,E_a,E_b)
        do JP = 1, 2*(J_max + 1)
            J = (JP - 1)/2
            P = mod(JP - 1,2) + 1
            N_qp = Orb_2qp%N(P,J)
            do i_qp = 1, N_qp
                a = Orb_2qp%Ind(P,J)%v(i_qp)%a
                b = Orb_2qp%Ind(P,J)%v(i_qp)%b
                T_ab = Orb_2qp%Ind(P,J)%v(i_qp)%T
                do j_qp = 1, i_qp
                    c = Orb_2qp%Ind(P,J)%v(j_qp)%a
                    d = Orb_2qp%Ind(P,J)%v(j_qp)%b
                    T_cd = Orb_2qp%Ind(P,J)%v(j_qp)%T
                    ASum = 0.0_dp
                    BSum = 0.0_dp
                    Norm = 1.0_dp / sqrt(real((1 + kronecker_delta(a,b)) * (1 + kronecker_delta(c,d)),dp))
                    if (a == c .and. b == d .and. T_ab == T_cd) then
                        if (T_ab == -1) then
                            E_a = H_N%qp11%p(a,a)
                            E_b = H_N%qp11%p(b,b)
                        else
                            E_a = H_N%qp11%n(a,a)
                            E_b = H_N%qp11%n(b,b)
                        end if
                        ASum = ASum + (E_a + E_b)
                    end if
                    if (T_ab == -1 .and. T_cd == -1) then
                        ASum = ASum + Norm * qpO2b_22_pp(a,b,c,d,J,P,H_NN,Orb,Orb_NN)
                        if (With_B) BSum = BSum + Norm * qpO2b_40_pp(a,b,c,d,J,P,H_NN,Orb,Orb_NN)
                    else if (T_ab == -1 .and. T_cd == 1) then
                        ASum = ASum + Norm * qpO2b_2002_pn(a,b,c,d,J,P,H_NN,Orb,Orb_NN)
                        if (With_B) BSum = BSum + Norm * qpO2b_40_pn(a,b,c,d,J,P,H_NN,Orb,Orb_NN)
                    else if (T_ab == 1 .and. T_cd == -1) then
                        ASum = ASum + Norm * qpO2b_2002_pn(c,d,a,b,J,P,H_NN,Orb,Orb_NN)
                        if (With_B) BSum = BSum + Norm * qpO2b_40_pn(c,d,a,b,J,P,H_NN,Orb,Orb_NN)
                    else
                        ASum = ASum + Norm * qpO2b_22_nn(a,b,c,d,J,P,H_NN,Orb,Orb_NN)
                        if (With_B) BSum = BSum + Norm * qpO2b_40_nn(a,b,c,d,J,P,H_NN,Orb,Orb_NN)
                    end if
                    Amat(P,J)%v(i_qp,j_qp) = ASum
                    Amat(P,J)%v(j_qp,i_qp) = ASum
                    if (With_B) then
                        Bmat(P,J)%v(i_qp,j_qp) = BSum
                        Bmat(P,J)%v(j_qp,i_qp) = BSum
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine qex_allocate

    ! 0+ spurious (particle-number) vector ...
    function spurious_PN_initialize(Orb,Orb_2qp,qpN) result(Spur)
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(qpo1b_t), intent(in) :: qpN
        real(dp), allocatable :: Spur(:)
        integer :: qp, a, b
        allocate(Spur(Orb_2qp%N(1,0)))
        do qp = 1, size(Spur)
            a = Orb_2qp%Ind(1,0)%v(qp)%a
            b = Orb_2qp%Ind(1,0)%v(qp)%b
            if (Orb_2qp%Ind(1,0)%v(qp)%T == -1) then
                Spur(qp) = sqrt(real(Orb(a)%j + 1,dp)) * qpN%qp20%p(a,b)
            else
                Spur(qp) = sqrt(real(Orb(a)%j + 1,dp)) * qpN%qp20%n(a,b)
            end if
        end do
        Spur = Spur / norm2(Spur)
    end function spurious_PN_initialize

    ! 1- spurious (centre-of-mass) vector ...
    function spurious_CM_initialize(Orb_2qp,qpTrOp) result(Spur)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(qptr1b_t), intent(in) :: qpTrOp
        real(dp), allocatable :: Spur(:)
        integer :: qp, a, b
        allocate(Spur(Orb_2qp%N(2,1)))
        do qp = 1, size(Spur)
            a = Orb_2qp%Ind(2,1)%v(qp)%a
            b = Orb_2qp%Ind(2,1)%v(qp)%b
            if (Orb_2qp%Ind(2,1)%v(qp)%T == -1) then
                Spur(qp) = qpTrOp%E1%qp20%p(a,b)
            else
                Spur(qp) = qpTrOp%E1%qp20%n(a,b)
            end if
        end do
        Spur = Spur / norm2(Spur)
    end function spurious_CM_initialize

    ! Basis of the 2qp space with the spurious vector projected out ... returns V & the index of the dropped vector ...
    subroutine spurious_orthogonalize(N_qp,Spur,V,Spur_ind)
        integer, intent(in) :: N_qp
        real(dp), intent(in) :: Spur(:)
        real(dp), intent(out) :: V(:,:)
        integer, intent(out) :: Spur_ind
        real(dp), allocatable :: Iden(:,:), w(:)
        integer :: a, nu, mu

        allocate(Iden(N_qp,N_qp))
        Iden = 0.0_dp
        do a = 1, N_qp
            Iden(a,a) = 1.0_dp
        end do
        do a = 1, N_qp
            Iden(:,a) = Iden(:,a) - dot_product(Spur,Iden(:,a)) * Spur
        end do
        do a = 1, N_qp
            Iden(:,a) = Iden(:,a) / norm2(Iden(:,a))
        end do

        Spur_ind = 0
        V(:,1) = Iden(:,1)
        do nu = 2, N_qp
            w = Iden(:,nu)
            do mu = 1, nu - 1
                w = w - dot_product(V(:,mu),w) * V(:,mu)
            end do
            if (abs(norm2(w)) < 1.0e-4_dp) Spur_ind = nu
            w = w / norm2(w)
            V(:,nu) = w
        end do
        if (Spur_ind == 0) call fatal('Spurious state orthogonalization failed ... no spurious direction found ...')
        V(:,Spur_ind) = Spur
    end subroutine spurious_orthogonalize

    ! Reduced transition matrix elements ... Amp = (X - Y)/sqrt(1 + delta_ab) (Y = 0 for QTDA) ...
    function qex_rM(Orb_2qp,X,Y,qpTrOp) result(rM)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(qptr1b_t), intent(in) :: qpTrOp
        type(reduced_multipole_t) :: rM

        call channel(0,1,qpTrOp%E0,rM%E0)
        call channel(1,2,qpTrOp%E1,rM%E1)
        call channel(1,2,qpTrOp%E1_VC,rM%E1_VC)
        call channel(1,2,qpTrOp%E1_VS,rM%E1_VS)
        call channel(1,2,qpTrOp%E1_TC,rM%E1_TC)
        call channel(1,2,qpTrOp%E1_TS,rM%E1_TS)
        call channel(1,2,qpTrOp%E1_C,rM%E1_C)
        call channel(1,2,qpTrOp%E1_sTC,rM%E1_sTC)
        call channel(1,2,qpTrOp%E1_sTS,rM%E1_sTS)
        call channel(1,2,qpTrOp%E1_sC,rM%E1_sC)
        call channel(2,1,qpTrOp%E2,rM%E2)
        call channel(3,2,qpTrOp%E3,rM%E3)
    contains
        subroutine channel(J,P,Op,M)
            integer, intent(in) :: J, P
            type(qpo1b_t), intent(in) :: Op
            type(pncvector_t), intent(out) :: M
            integer :: N_qp, nu, i_qp, a, b
            complex(dp) :: Amp, pS, nS
            N_qp = Orb_2qp%N(P,J)
            allocate(M%p(N_qp), M%n(N_qp))
            !$omp parallel do private(pS,nS,i_qp,a,b,Amp)
            do nu = 1, N_qp
                pS = (0.0_dp,0.0_dp)
                nS = (0.0_dp,0.0_dp)
                do i_qp = 1, N_qp
                    a = Orb_2qp%Ind(P,J)%v(i_qp)%a
                    b = Orb_2qp%Ind(P,J)%v(i_qp)%b
                    Amp = (X(P,J)%v(i_qp,nu) - Y(P,J)%v(i_qp,nu)) * (1.0_dp / sqrt(1.0_dp + kronecker_delta(a,b)))
                    if (Orb_2qp%Ind(P,J)%v(i_qp)%T == -1) then
                        pS = pS + Amp * Op%qp20%p(a,b)
                    else
                        nS = nS + Amp * Op%qp20%n(a,b)
                    end if
                end do
                M%p(nu) = pS
                M%n(nu) = nS
            end do
            !$omp end parallel do
        end subroutine channel
    end function qex_rM

    ! Reduced transition intensities B ... E is the (real part of the) phonon energy per channel ...
    function qex_rB(Params,Orb_2qp,rM,E) result(rB)
        type(params_t), intent(in) :: Params
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(reduced_multipole_t), intent(in) :: rM
        type(rvec_t), intent(in) :: E(:,0:)
        type(reduced_transition_t) :: rB
        real(dp), parameter :: hc = HbarC, g_p = 5.586_dp, g_n = -3.826_dp
        real(dp) :: e_p, e_n, En
        integer :: nu, N_qp

        ! E0 ...
        N_qp = Orb_2qp%N(1,0)
        call alloc3(rB%E0,N_qp)
        do nu = 1, N_qp
            rB%E0%ph(nu) = abs(rM%E0%p(nu))**2
            rB%E0%is(nu) = 0.25_dp * abs(rM%E0%p(nu) + rM%E0%n(nu))**2
            rB%E0%iv(nu) = 0.25_dp * abs(rM%E0%p(nu) - rM%E0%n(nu))**2
        end do

        ! E1 ...
        N_qp = Orb_2qp%N(2,1)
        call alloc3(rB%E1,N_qp)
        call alloc3(rB%E1_V,N_qp)
        call alloc3(rB%E1_VC,N_qp)
        call alloc3(rB%E1_VS,N_qp)
        call alloc3(rB%E1_T,N_qp)
        call alloc3(rB%E1_TC,N_qp)
        call alloc3(rB%E1_TS,N_qp)
        call alloc3(rB%E1_C,N_qp)
        call alloc3(rB%E1_CC,N_qp)
        call alloc3(rB%E1_CS,N_qp)
        call alloc3(rB%E1_comp,N_qp)
        call alloc3(rB%E1_NLO_LWA,N_qp)
        e_p = real(Params%Calc%A - Params%Calc%Z,dp) / real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp) / real(Params%Calc%A,dp)
        associate (M => rM)
        do nu = 1, N_qp
            En = E(2,1)%v(nu)
            ! Physical components ...
            rB%E1%ph(nu) = abs(M%E1%p(nu))**2
            rB%E1_V%ph(nu) = abs(M%E1_VC%p(nu) + 0.5_dp * g_p * M%E1_VS%p(nu) + 0.5_dp * g_n * M%E1_VS%n(nu))**2
            rB%E1_VC%ph(nu) = abs(M%E1_VC%p(nu))**2
            rB%E1_VS%ph(nu) = abs(0.5_dp * g_p * M%E1_VS%p(nu) + 0.5_dp * g_n * M%E1_VS%n(nu))**2
            rB%E1_T%ph(nu) = abs(M%E1_TC%p(nu) + 0.5_dp * g_p * M%E1_TS%p(nu) + 0.5_dp * g_n * M%E1_TS%n(nu))**2
            rB%E1_TC%ph(nu) = abs(M%E1_TC%p(nu))**2
            rB%E1_TS%ph(nu) = abs(0.5_dp * g_p * M%E1_TS%p(nu) + 0.5_dp * g_n * M%E1_TS%n(nu))**2
            rB%E1_C%ph(nu) = abs(M%E1_TC%p(nu) + 0.5_dp * g_p * M%E1_TS%p(nu) + 0.5_dp * g_n * M%E1_TS%n(nu) - &
                             (M%E1_VC%p(nu) + 0.5_dp * g_p * M%E1_VS%p(nu) + 0.5_dp * g_n * M%E1_VS%n(nu)))**2
            rB%E1_CC%ph(nu) = abs(M%E1_TC%p(nu) - (M%E1_VC%p(nu)))**2
            rB%E1_CS%ph(nu) = abs(0.5_dp * g_p * M%E1_TS%p(nu) + 0.5_dp * g_n * M%E1_TS%n(nu) - &
                              (0.5_dp * g_p * M%E1_VS%p(nu) + 0.5_dp * g_n * M%E1_VS%n(nu)))**2
            rB%E1_comp%ph(nu) = abs(M%E1_C%p(nu))**2
            rB%E1_NLO_LWA%ph(nu) = abs(M%E1%p(nu) + En / hc * (M%E1_TC%p(nu) + 0.5_dp * g_p * M%E1_TS%p(nu) + &
                                   0.5_dp * g_n * M%E1_TS%n(nu)))**2
            ! Isoscalar components ...
            rB%E1%is(nu) = 0.25_dp * abs(M%E1%p(nu) + M%E1%n(nu))**2
            rB%E1_V%is(nu) = abs(0.5_dp * (M%E1_VC%p(nu) + M%E1_VC%n(nu)) + 0.125_dp * (g_p + g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu)))**2
            rB%E1_VC%is(nu) = abs(0.5_dp * (M%E1_VC%p(nu) + M%E1_VC%n(nu)))**2
            rB%E1_VS%is(nu) = abs(0.125_dp * (g_p + g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu)))**2
            rB%E1_T%is(nu) = abs(0.5_dp * (M%E1_TC%p(nu) + M%E1_sTC%p(nu) + M%E1_TC%n(nu) + M%E1_sTC%n(nu)) + &
                             0.125_dp * (g_p + g_n) * (M%E1_TS%p(nu) + M%E1_sTS%p(nu) + M%E1_TS%n(nu) + M%E1_sTS%n(nu)))**2
            rB%E1_TC%is(nu) = abs(0.5_dp * (M%E1_TC%p(nu) + M%E1_sTC%p(nu) + M%E1_TC%n(nu) + M%E1_sTC%n(nu)))**2
            rB%E1_TS%is(nu) = abs(0.125_dp * (g_p + g_n) * (M%E1_TS%p(nu) + M%E1_sTS%p(nu) + M%E1_TS%n(nu) + M%E1_sTS%n(nu)))**2
            rB%E1_C%is(nu) = abs(0.5_dp * (M%E1_TC%p(nu) + M%E1_sTC%p(nu) + M%E1_TC%n(nu) + M%E1_sTC%n(nu)) + &
                             0.125_dp * (g_p + g_n) * (M%E1_TS%p(nu) + M%E1_sTS%p(nu) + M%E1_TS%n(nu) + M%E1_sTS%n(nu)) - &
                             (0.5_dp * (M%E1_VC%p(nu) + M%E1_VC%n(nu)) + 0.125_dp * (g_p + g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu))))**2
            rB%E1_CC%is(nu) = abs(0.5_dp * (M%E1_TC%p(nu) + M%E1_sTC%p(nu) + M%E1_TC%n(nu) + M%E1_sTC%n(nu)) - &
                              (0.5_dp * (M%E1_VC%p(nu) + M%E1_VC%n(nu))))**2
            rB%E1_CS%is(nu) = abs(0.125_dp * (g_p + g_n) * (M%E1_TS%p(nu) + M%E1_sTS%p(nu) + M%E1_TS%n(nu) + M%E1_sTS%n(nu)) - &
                              (0.125_dp * (g_p + g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu))))**2
            rB%E1_comp%is(nu) = 0.25_dp * abs(M%E1_C%p(nu) + M%E1_sC%p(nu) + M%E1_C%n(nu) + M%E1_sC%n(nu))**2
            rB%E1_NLO_LWA%is(nu) = 0.25_dp * abs(M%E1%p(nu) + M%E1%n(nu) + En / hc * ((M%E1_TC%p(nu) + M%E1_sTC%p(nu) + &
                                   M%E1_TC%n(nu) + M%E1_sTC%n(nu)) + 0.25_dp * (g_p + g_n) * (M%E1_TS%p(nu) + M%E1_sTS%p(nu) + &
                                   M%E1_TS%n(nu) + M%E1_sTS%n(nu))))**2
            ! Isovector components ...
            rB%E1%iv(nu) = abs(e_p * M%E1%p(nu) - e_n * M%E1%n(nu))**2
            rB%E1_V%iv(nu) = abs(0.5_dp * (M%E1_VC%p(nu) - M%E1_VC%n(nu)) + 0.125_dp * (g_p - g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu)))**2
            rB%E1_VC%iv(nu) = abs(0.5_dp * (M%E1_VC%p(nu) - M%E1_VC%n(nu)))**2
            rB%E1_VS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu)))**2
            rB%E1_T%iv(nu) = abs(0.5_dp * (M%E1_TC%p(nu) - M%E1_TC%n(nu)) + 0.125_dp * (g_p - g_n) * (M%E1_TS%p(nu) + M%E1_TS%n(nu)))**2
            rB%E1_TC%iv(nu) = abs(0.5_dp * (M%E1_TC%p(nu) - M%E1_TC%n(nu)))**2
            rB%E1_TS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (M%E1_TS%p(nu) + M%E1_TS%n(nu)))**2
            rB%E1_C%iv(nu) = abs(0.5_dp * (M%E1_TC%p(nu) - M%E1_TC%n(nu)) + 0.125_dp * (g_p - g_n) * (M%E1_TS%p(nu) + M%E1_TS%n(nu)) - &
                             (0.5_dp * (M%E1_VC%p(nu) - M%E1_VC%n(nu)) + 0.125_dp * (g_p - g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu))))**2
            rB%E1_CC%iv(nu) = abs(0.5_dp * (M%E1_TC%p(nu) - M%E1_TC%n(nu)) - (0.5_dp * (M%E1_VC%p(nu) - M%E1_VC%n(nu))))**2
            rB%E1_CS%iv(nu) = abs(0.125_dp * (g_p - g_n) * (M%E1_TS%p(nu) + M%E1_TS%n(nu)) - &
                              (0.125_dp * (g_p - g_n) * (M%E1_VS%p(nu) + M%E1_VS%n(nu))))**2
            rB%E1_comp%iv(nu) = 0.25_dp * abs(M%E1_C%p(nu) - M%E1_C%n(nu))**2
            rB%E1_NLO_LWA%iv(nu) = abs(e_p * M%E1%p(nu) - e_n * M%E1%n(nu) + 0.5_dp * En / hc * ((M%E1_TC%p(nu) - M%E1_TC%n(nu)) + &
                                   0.125_dp * (g_p - g_n) * (M%E1_TS%p(nu) + M%E1_TS%n(nu))))**2
        end do
        end associate

        ! E2 & E3 ...
        N_qp = Orb_2qp%N(1,2)
        call alloc3(rB%E2,N_qp)
        do nu = 1, N_qp
            rB%E2%ph(nu) = abs(rM%E2%p(nu))**2
            rB%E2%is(nu) = 0.25_dp * abs(rM%E2%p(nu) + rM%E2%n(nu))**2
            rB%E2%iv(nu) = 0.25_dp * abs(rM%E2%p(nu) - rM%E2%n(nu))**2
        end do
        N_qp = Orb_2qp%N(2,3)
        call alloc3(rB%E3,N_qp)
        do nu = 1, N_qp
            rB%E3%ph(nu) = abs(rM%E3%p(nu))**2
            rB%E3%is(nu) = 0.25_dp * abs(rM%E3%p(nu) + rM%E3%n(nu))**2
            rB%E3%iv(nu) = 0.25_dp * abs(rM%E3%p(nu) - rM%E3%n(nu))**2
        end do
    contains
        subroutine alloc3(Tr,N)
            type(transition_t), intent(out) :: Tr
            integer, intent(in) :: N
            allocate(Tr%ph(N), Tr%is(N), Tr%iv(N))
        end subroutine alloc3
    end function qex_rB

    ! Dimension tables of the summary file ...
    subroutine qex_dimensions_summary(u,Orb_2qp,J_max,Label)
        integer, intent(in) :: u, J_max
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        character(len=*), intent(in) :: Label
        integer :: J, N_qp
        N_qp = 0
        do J = 0, J_max
            N_qp = N_qp + (Orb_2qp%N(1,J) + Orb_2qp%N(2,J)) * (2*J + 1)
        end do
        write(u,'(/,a)') 'Information on dimensions of ' // Label
        write(u,'(/,a)') 'Total number of 1-phonon quasiparticle states in J-scheme = ' // istr(sum(Orb_2qp%N))
        write(u,'(a)') 'Total number of 1-phonon quasiparticle states in M-scheme = ' // istr(N_qp)
        write(u,'(/,a)') 'J' // TAB // 'N_phonon'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_2qp%N(1,J) + Orb_2qp%N(2,J))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = +'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_2qp%N(1,J))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = -'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_2qp%N(2,J))
        end do
    end subroutine qex_dimensions_summary

    ! inquire(exist) ... Julia isfile ...
    logical function qex_file_exists(Path)
        character(len=*), intent(in) :: Path
        inquire(file=trim(Path), exist=qex_file_exists)
    end function qex_file_exists

    ! Precondition of dev QTDA.jl / QRPA.jl / pnQTDA.jl: the given reference files exist in io/<Path>/<Reference>/bin ...
        ! otherwise the dev error text is printed & .false. returned ...
    logical function qex_reference_check(Params,Ref,Files) result(Ok)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Ref
        character(len=*), intent(in) :: Files(:)
        integer :: k
        Ok = .true.
        do k = 1, size(Files)
            if (.not. qex_file_exists(io_path(Params,'/' // trim(Ref) // '/bin/' // trim(Files(k))))) Ok = .false.
        end do
        if (.not. Ok) write(*,'(/,a)') 'Error! ... No precomputed reference quasiparticle mean-field (HFB, BCS) solution is ' // &
                                       'available in given Input_File path ... run HFB/BCS solver first ...'
    end function qex_reference_check

    ! Plot-ready spectrum ... dev writes N_qp of the N_qp + 1 sorted rows (including the 0+ ground state) ... file
        ! spectra/QTDA_plot.dat (QTDA) or spectra/QRPA_Plot.dat (QRPA), as dev ...
    subroutine qex_spectrum_plot_export(Params,Orb_2qp,E,Label)
        type(params_t), intent(in) :: Params
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(rvec_t), intent(in) :: E(:,0:)
        character(len=*), intent(in) :: Label
        real(dp), allocatable :: Ex(:,:)
        integer, allocatable :: Sort(:)
        integer :: J_max, N_qp, qp_count, J, P, nu, u
        character(len=1) :: Ps

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Preparing plot-ready export of ' // Label // ' spectra ...'
        N_qp = sum(Orb_2qp%N)
        allocate(Ex(4,N_qp+1), Sort(N_qp+1))
        Ex(:,1) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
        qp_count = 1
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, Orb_2qp%N(P,J)
                    qp_count = qp_count + 1
                    Ex(:,qp_count) = [real(J,dp), real(P,dp), E(P,J)%v(nu), E(P,J)%v(nu)]
                end do
            end do
        end do
        call sortperm(Ex(3,:),Sort)
        Ex = Ex(:,Sort)
        do nu = 2, N_qp + 1
            if (Ex(4,nu) < Ex(4,nu-1)) then
                Ex(4,nu) = Ex(4,nu-1) + 1.0_dp
            else if (abs(Ex(4,nu) - Ex(4,nu-1)) < 1.0_dp) then
                Ex(4,nu) = Ex(4,nu-1) + 1.0_dp
            end if
        end do

        write(*,'(/,a)') 'Performing the export of plot-ready ' // Label // ' spectra ...'
        if (Label == 'QTDA') then
            open(newunit=u, file=io_path(Params,'/QTDA/spectra/QTDA_plot.dat'), status='replace', action='write')
        else
            open(newunit=u, file=io_path(Params,'/' // Label // '/spectra/' // Label // '_Plot.dat'), status='replace', action='write')
        end if
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // cfmt_s('label',12,.true.)
        do nu = 1, N_qp
            Ps = merge('+', '-', abs(Ex(2,nu) - 1.0_dp) < 1.0e-3_dp)
            write(u,'(a)') cfmt_d(nint(Ex(1,nu)),5,.true.) // ' ' // cfmt_s(Ps,5,.true.) // ' ' // cfmt_f(Ex(3,nu),12,6,.true.) // &
                           ' ' // cfmt_f(Ex(4,nu),12,6,.true.)
        end do
        close(u)
    end subroutine qex_spectrum_plot_export

    ! Sum-rule moments m_k, TRK values & B(EX) tables ... QTDA_transitions_export / QRPA_transitions_export ...
    subroutine qex_transitions_export(Params,Orb,Orb_2qp,Rho,E,rB,Label,Ref)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(o1b_t), intent(in) :: Rho
        type(rvec_t), intent(in) :: E(:,0:)
        type(reduced_transition_t), intent(in) :: rB
        character(len=*), intent(in) :: Label, Ref       ! Ref ... reference folder (Calc%QTDA/QRPA%Reference)
        type(o1b_t) :: C
        real(dp) :: pR2, pR4, Zr, m(5,3,4), TRK(4)
        character(len=3) :: Unit_MeV
        character(len=:), allocatable :: Dir
        integer :: u, k, iso, X
        character(len=*), parameter :: Kind_Name(3) = [character(len=10) :: 'physical', 'isoscalar', 'isovector']
        character(len=*), parameter :: Fm(4) = [character(len=6) :: 'fm^4', 'fm^2', 'fm^4', 'fm^6']

        Zr = real(Params%Calc%Z,dp)
        if (Label == 'QTDA') write(*,'(/,a)') 'Preparing export of information on QTDA 1-phonon electromagnetic transitions ...'
        C = O1b_import(Params,Orb,io_path(Params,'/' // trim(Ref) // '/bin/C.bin'))
        pR2 = rho_rN(Params,2,Rho,C,Orb)
        pR4 = rho_rN(Params,4,Rho,C,Orb)
        if (Label == 'QRPA') write(*,'(/,a)') 'Preparing export of information on QRPA 1-phonon electromagnetic transitions ...'

        TRK(1) = HbarC**2 / m_p / (2.0_dp * pi) * Zr * pR2
        TRK(2) = HbarC**2 / m_p * 9.0_dp / (8.0_dp * pi) * Zr
        TRK(3) = HbarC**2 / m_p * 50.0_dp / (8.0_dp * pi) * Zr * pR2
        TRK(4) = HbarC**2 / m_p * 147.0_dp / (8.0_dp * pi) * Zr * pR4

        ! m(k, iso, X) ... k = 1..5 <-> m_-1, m_0, m_1, m_2, m_3 ... X = E0..E3 ...
        m = 0.0_dp
        call moments(E(1,0)%v, rB%E0, 1, .true.)
        call moments(E(2,1)%v, rB%E1, 2, .true.)
        call moments(E(1,2)%v, rB%E2, 3, .false.)
        call moments(E(2,3)%v, rB%E3, 4, .false.)

        open(newunit=u, file=io_path(Params,'/' // Label // '/' // Label // '_summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Review of ' // Label // ' electrogmagnetic transition moments m_k ...'
        do iso = 1, 3
            do X = 1, 4
                if (X == 1 .and. iso > 1) write(u,'(/,/,/)',advance='no')
                write(u,'(/,a)') TAB // 'E' // istr(X-1) // ' ' // trim(Kind_Name(iso)) // ':'
                write(u,'(a)') 'm _-1 = ' // cfmt_f(m(1,iso,X),20,8) // TAB // 'e^2 ' // trim(Fm(X)) // ' / MeV'
                write(u,'(a)') 'm _0  = ' // cfmt_f(m(2,iso,X),20,8) // TAB // 'e^2 ' // trim(Fm(X))
                write(u,'(a)') 'm _1  = ' // cfmt_f(m(3,iso,X),20,8) // TAB // 'e^2 ' // trim(Fm(X)) // ' MeV'
                write(u,'(a)') 'm _2  = ' // cfmt_f(m(4,iso,X),20,8) // TAB // 'e^2 ' // trim(Fm(X)) // ' MeV^2'
                write(u,'(a)') 'm _3  = ' // cfmt_f(m(5,iso,X),20,8) // TAB // 'e^2 ' // trim(Fm(X)) // ' MeV^3'
            end do
        end do
        ! The QRPA summary of the original spells the TRK unit "Mev" ...
        Unit_MeV = merge('Mev', 'MeV', Label == 'QRPA')
        write(u,'(/,/,/,/,a,/)') 'Review of ' // Label // ' physical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values ...'
        do X = 1, 4
            write(u,'(a)') 'S_E' // istr(X-1) // ' = ' // cfmt_f(TRK(X),20,8) // TAB // 'e^2 ' // trim(Fm(X)) // ' ' // Unit_MeV
        end do
        close(u)

        write(*,'(a)') TAB // 'Preparing export of ' // Label // ' 1-phonon electromagnetic transition intensities ...'
        Dir = io_path(Params,'/' // Label // '/transitions/')
        call table(Dir // 'E0/' // Label // '_E0.dat', E(1,0)%v, rB%E0)
        call table(Dir // 'E1/' // Label // '_E1.dat', E(2,1)%v, rB%E1)
        call table(Dir // 'E1/' // Label // '_E1_V.dat', E(2,1)%v, rB%E1_V)
        call table(Dir // 'E1/' // Label // '_E1_V_conv.dat', E(2,1)%v, rB%E1_VC)
        call table(Dir // 'E1/' // Label // '_E1_V_spin.dat', E(2,1)%v, rB%E1_VS)
        call table(Dir // 'E1/' // Label // '_E1_T.dat', E(2,1)%v, rB%E1_T)
        call table(Dir // 'E1/' // Label // '_E1_T_conv.dat', E(2,1)%v, rB%E1_TC)
        call table(Dir // 'E1/' // Label // '_E1_T_spin.dat', E(2,1)%v, rB%E1_TS)
        call table(Dir // 'E1/' // Label // '_E1_C.dat', E(2,1)%v, rB%E1_C)
        call table(Dir // 'E1/' // Label // '_E1_C_conv.dat', E(2,1)%v, rB%E1_CC)
        call table(Dir // 'E1/' // Label // '_E1_C_spin.dat', E(2,1)%v, rB%E1_CS)
        call table(Dir // 'E1/' // Label // '_E1_c.dat', E(2,1)%v, rB%E1_comp)
        call table(Dir // 'E1/' // Label // '_E1_NLO_LWA.dat', E(2,1)%v, rB%E1_NLO_LWA)
        call table(Dir // 'E2/' // Label // '_E2.dat', E(1,2)%v, rB%E2)
        ! The QRPA E3 table of the original is named "QRPAE3.dat" (no underscore) ...
        if (Label == 'QRPA') then
            call table(Dir // 'E3/QRPAE3.dat', E(2,3)%v, rB%E3)
        else
            call table(Dir // 'E3/' // Label // '_E3.dat', E(2,3)%v, rB%E3)
        end if
        write(*,'(a)') TAB // TAB // 'Export of information on ' // Label // ' 1-phonon electromagnetic transitions completed ...'

    contains

        ! E0 & E1 accumulate m_-1 & m_0 only for |E| > 1e-3, E2 & E3 accumulate m_0 always ...
        subroutine moments(En,Tr,X,M0_cut)
            real(dp), intent(in) :: En(:)
            type(transition_t), intent(in) :: Tr
            integer, intent(in) :: X
            logical, intent(in) :: M0_cut
            real(dp) :: B(3), Ev
            integer :: nu, i
            do nu = 1, size(En)
                Ev = En(nu)
                B = [Tr%ph(nu), Tr%is(nu), Tr%iv(nu)]
                do i = 1, 3
                    if (abs(Ev) > 1.0e-3_dp) then
                        m(1,i,X) = m(1,i,X) + B(i) / Ev
                        if (M0_cut) m(2,i,X) = m(2,i,X) + B(i)
                    end if
                    if (.not. M0_cut) m(2,i,X) = m(2,i,X) + B(i)
                    m(3,i,X) = m(3,i,X) + Ev * B(i)
                    m(4,i,X) = m(4,i,X) + Ev**2 * B(i)
                    m(5,i,X) = m(5,i,X) + Ev**3 * B(i)
                end do
            end do
        end subroutine moments

        subroutine table(File,En,Tr)
            character(len=*), intent(in) :: File
            real(dp), intent(in) :: En(:)
            type(transition_t), intent(in) :: Tr
            integer :: uf, nu
            open(newunit=uf, file=File, status='replace', action='write')
            write(uf,'(a)') cfmt_s('E',20,.true.) // ' ' // cfmt_s('B_ph',20,.true.) // ' ' // cfmt_s('B_is',20,.true.) // ' ' // &
                            cfmt_s('B_iv',20,.true.)
            do nu = 1, size(En)
                write(uf,'(a)') cfmt_f(En(nu),20,8,.true.) // ' ' // cfmt_f(Tr%ph(nu),20,8,.true.) // ' ' // &
                                cfmt_f(Tr%is(nu),20,8,.true.) // ' ' // cfmt_f(Tr%iv(nu),20,8,.true.)
            end do
            close(uf)
        end subroutine table

    end subroutine qex_transitions_export

end module nj_qex
