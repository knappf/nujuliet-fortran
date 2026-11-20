! Electromagnetic & weak 1-body transition operators ... counterpart of src/common/Tr1b.jl (NuJuliet dev; MatrixElements/Tr1B.jl
    ! before) ... dev additions: Fermi & Gamow-Teller operators F & GT (proton-neutron mixing in the basis transformations),
    ! the convective E1 vortical/toroidal fix (rME_nabla_dot_rN_YJL, s = sigma/2 in rME_nabla_cross_S_dot_rN_YJL, 0.5 of the
    ! convective VC/TC/sTC) & the quasiparticle 11/20 forms with O_lk & the time-odd phase ...
module nj_tr1b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : fCG, f6j
    use nj_functions, only : radial_moment_LHO
    use nj_o1b, only : o1b_alloc
    implicit none
    private

    ! Julia struct Tr1B ... F & GT: %p = (p||O||n), %n = (n||O||p) ...
    type, public :: tr1b_t
        type(o1b_t) :: E0, E1, E2, E3, E1_VC, E1_VS, E1_TC, E1_TS, E1_sTC, E1_sTS, E1_C, E1_sC, M1, M2, M3, F, GT
    end type tr1b_t

    ! Julia struct qpTr1B ...
    type, public :: qptr1b_t
        type(qpo1b_t) :: E0, E1, E2, E3, E1_VC, E1_VS, E1_TC, E1_TS, E1_sTC, E1_sTS, E1_C, E1_sC, M1, M2, M3, F, GT
    end type qptr1b_t

    public :: Tr1b_initialize, Tr1b_transformation, qpTr1b_initialize

contains

    function Tr1b_initialize(Params,Orb,chR2) result(TrOp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: chR2
        type(tr1b_t) :: TrOp
        write(*,'(/,a)') 'Constructing electromagnetic & weak transition operators E^lambda, M^lambda, M^F & M^GT in the LHO basis ...'
        call Tr1b_EX_initialize(Params,Orb,TrOp)
        call Tr1b_E1_NLO_initialize(Params,Orb,chR2,TrOp)
        call Tr1b_MX_initialize(Params,Orb,TrOp)
        call Tr1b_BX_initialize(Params,Orb,TrOp)
        write(*,'(a)') TAB // 'Transition 1-body electromagnetic & weak operators in the LHO basis succesfully constructed ...'
    end function Tr1b_initialize

    pure real(dp) function b_osc_of(hw)
        real(dp), intent(in) :: hw
        b_osc_of = sqrt(0.5_dp * (m_p + m_n) * hw) / HbarC
    end function b_osc_of

    ! Reduced EX matrix element ...
    pure real(dp) function rEX(n_a,l_a,j_a,n_b,l_b,j_b,X,b_osc)
        integer, intent(in) :: n_a, l_a, j_a, n_b, l_b, j_b, X
        real(dp), intent(in) :: b_osc
        real(dp) :: Amp
        Amp = phase(X + (j_a - 1)/2) * sqrt(real((j_a + 1) * (j_b + 1),dp)) * fCG(j_a,j_b,2*X,1,-1,0) / sqrt(4.0_dp * pi)
        if (X == 0) then
            rEX = Amp * radial_moment_LHO(X+2,n_a,l_a,n_b,l_b,b_osc)
        else
            rEX = Amp * radial_moment_LHO(X,n_a,l_a,n_b,l_b,b_osc)
        end if
    end function rEX

    subroutine Tr1b_EX_initialize(Params,Orb,TrOp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(tr1b_t), intent(inout) :: TrOp
        integer :: a_max, a, b, n_a, l_a, j_a, n_b, l_b, j_b
        real(dp) :: b_osc, rME

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        b_osc = b_osc_of(Params%Calc%hw)
        write(*,'(a)') TAB // 'Constructing 1-body EX transition operators in the LHO basis ...'
        call o1b_alloc(TrOp%E0,a_max)
        call o1b_alloc(TrOp%E1,a_max)
        call o1b_alloc(TrOp%E2,a_max)
        call o1b_alloc(TrOp%E3,a_max)

        !$omp parallel do private(n_a,l_a,j_a,b,n_b,l_b,j_b,rME)
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                if (mod(l_a + l_b,2) == 0) then
                    rME = rEX(n_a,l_a,j_a,n_b,l_b,j_b,0,b_osc)
                    TrOp%E0%p(a,b) = rME
                    TrOp%E0%n(a,b) = rME
                end if
                if (mod(l_a + l_b + 1,2) == 0) then
                    rME = rEX(n_a,l_a,j_a,n_b,l_b,j_b,1,b_osc)
                    TrOp%E1%p(a,b) = rME
                    TrOp%E1%n(a,b) = rME
                end if
                if (mod(l_a + l_b + 2,2) == 0) then
                    rME = rEX(n_a,l_a,j_a,n_b,l_b,j_b,2,b_osc)
                    TrOp%E2%p(a,b) = rME
                    TrOp%E2%n(a,b) = rME
                end if
                if (mod(l_a + l_b + 3,2) == 0) then
                    rME = rEX(n_a,l_a,j_a,n_b,l_b,j_b,3,b_osc)
                    TrOp%E3%p(a,b) = rME
                    TrOp%E3%n(a,b) = rME
                end if
            end do
        end do
        !$omp end parallel do
        write(*,'(a)') TAB // 'EX 1-body transition operators succesfully constructed in the LHO basis ...'
    end subroutine Tr1b_EX_initialize

    ! Radial moment <a| r^N |b> in the LHO basis ...
    pure real(dp) function rME_rN(n_a,l_a,n_b,l_b,N,hw)
        integer, intent(in) :: n_a, l_a, n_b, l_b, N
        real(dp), intent(in) :: hw
        rME_rN = radial_moment_LHO(N,n_a,l_a,n_b,l_b,b_osc_of(hw))
    end function rME_rN

    ! Radial matrix element <a| r^N d/dr |b> in the LHO basis ...
    pure real(dp) function rME_rN_dr(n_a,l_a,n_b,l_b,N,hw)
        integer, intent(in) :: n_a, l_a, n_b, l_b, N
        real(dp), intent(in) :: hw
        real(dp) :: b_osc
        b_osc = b_osc_of(hw)
        rME_rN_dr = real(l_b,dp) * radial_moment_LHO(N-1,n_a,l_a,n_b,l_b,b_osc) - b_osc**2 * radial_moment_LHO(N+1,n_a,l_a,n_b,l_b,b_osc)
        if (1 <= n_b) rME_rN_dr = rME_rN_dr - 2.0_dp * b_osc * sqrt(real(n_b,dp)) * radial_moment_LHO(N,n_a,l_a,n_b-1,l_b+1,b_osc)
    end function rME_rN_dr

    ! Convective current combination <a|| A . (psi_a^+ grad psi_b - (grad psi_a)^+ psi_b) ||b>, i.e. the (Ket - Bra) term ... (dev
        ! fix of the Bra terms: sign & sqrt((l_a+1)/(2l_a+1)) of the l_a+1 branch, without the l_a >= 1 condition) ...
    pure real(dp) function rME_nabla_dot_rN_YJL(a,b,N,L,J,hw,Orb) result(rME)
        integer, intent(in) :: a, b, N, L, J
        real(dp), intent(in) :: hw
        type(orb1b_t), intent(in) :: Orb(:)
        integer :: n_a, l_a, j_a, n_b, l_b, j_b
        real(dp) :: Amp_Bra, Amp_Ket, ME
        n_a = Orb(a)%n
        l_a = Orb(a)%l
        j_a = Orb(a)%j
        n_b = Orb(b)%n
        l_b = Orb(b)%l
        j_b = Orb(b)%j
        rME = 0.0_dp
        Amp_Bra = phase(l_b + (j_b + 1)/2) * f6j(2*l_b,1,j_b,j_a,2*J,2*l_a) * &
                  sqrt(real((2*J + 1) * (2*L + 1) * (j_a + 1) * (j_b + 1) * (2*l_a + 1) * (2*l_b + 1),dp) / (4.0_dp * pi))
        Amp_Ket = phase(l_a + (j_b + 1)/2) * f6j(2*l_a,1,j_a,j_b,2*J,2*l_b) * &
                  sqrt(real((2*J + 1) * (2*L + 1) * (j_a + 1) * (j_b + 1) * (2*l_a + 1) * (2*l_b + 1),dp) / (4.0_dp * pi))
        if (abs(L - l_b) <= (l_a - 1) .and. (l_a - 1) <= (L + l_b)) then
            ME = Amp_Bra * fCG(2*l_b,2*L,2*l_a-2,0,0,0) * f6j(2*L,2,2*J,2*l_a,2*l_b,2*l_a-2) * sqrt(real(l_a,dp) / real(2*l_a + 1,dp)) * &
                 (rME_rN_dr(n_b,l_b,n_a,l_a,N,hw) + real(l_a + 1,dp) * rME_rN(n_a,l_a,n_b,l_b,N-1,hw))
            rME = rME - ME
        end if
        if (abs(L - l_b) <= (l_a + 1) .and. (l_a + 1) <= (L + l_b)) then
            ME = Amp_Bra * fCG(2*l_b,2*L,2*l_a+2,0,0,0) * f6j(2*L,2,2*J,2*l_a,2*l_b,2*l_a+2) * sqrt(real(l_a + 1,dp) / real(2*l_a + 1,dp)) * &
                 (rME_rN_dr(n_b,l_b,n_a,l_a,N,hw) - real(l_a,dp) * rME_rN(n_a,l_a,n_b,l_b,N-1,hw))
            rME = rME + ME
        end if
        if (abs(L - l_a) <= (l_b - 1) .and. (l_b - 1) <= (L + l_a)) then
            ME = Amp_Ket * fCG(2*l_a,2*L,2*l_b-2,0,0,0) * f6j(2*L,2,2*J,2*l_b,2*l_a,2*l_b-2) * sqrt(real(l_b,dp) / real(2*l_b + 1,dp)) * &
                 (rME_rN_dr(n_a,l_a,n_b,l_b,N,hw) + real(l_b + 1,dp) * rME_rN(n_a,l_a,n_b,l_b,N-1,hw))
            rME = rME + ME
        end if
        if (abs(L - l_a) <= (l_b + 1) .and. (l_b + 1) <= (L + l_a)) then
            ME = Amp_Ket * fCG(2*l_a,2*L,2*l_b+2,0,0,0) * f6j(2*L,2,2*J,2*l_b,2*l_a,2*l_b+2) * sqrt(real(l_b + 1,dp) / real(2*l_b + 1,dp)) * &
                 (rME_rN_dr(n_a,l_a,n_b,l_b,N,hw) - real(l_b,dp) * rME_rN(n_a,l_a,n_b,l_b,N-1,hw))
            rME = rME - ME
        end if
    end function rME_nabla_dot_rN_YJL

    ! Normalized to the spin s = sigma/2, i.e. <a|| A . curl(psi_a^+ s psi_b) ||b> (dev: 16 pi instead of 8 pi) ...
    pure real(dp) function rME_nabla_cross_S_dot_rN_YJL(a,b,N,L,J,hw,Orb) result(rME)
        integer, intent(in) :: a, b, N, L, J
        real(dp), intent(in) :: hw
        type(orb1b_t), intent(in) :: Orb(:)
        integer :: n_a, l_a, j_a, n_b, l_b, j_b
        real(dp) :: Amp
        n_a = Orb(a)%n
        l_a = Orb(a)%l
        j_a = Orb(a)%j
        n_b = Orb(b)%n
        l_b = Orb(b)%l
        j_b = Orb(b)%j
        rME = 0.0_dp
        if (abs(l_a - l_b) > J .or. J > (l_a + l_b)) return
        Amp = phase(l_a + l_b + J + (j_b + 3)/2) * sqrt(real((j_a + 1) * (j_b + 1) * (2*l_a + 1) * (2*l_b + 1),dp) / &
              (16.0_dp * pi * real(J * (J + 1),dp))) * &
              real((l_a - l_b) * (l_a + l_b + 1) - ((j_a - j_b) * (j_a + j_b + 2))/4,dp) * fCG(2*l_a,2*l_b,2*J,0,0,0) * &
              f6j(2*l_a,j_a,1,j_b,2*l_b,2*J)
        if (abs(Amp) < 1.0e-12_dp) return
        if (L == J + 1) rME = rME + Amp * sqrt(real(J,dp) / real(2*J + 1,dp)) * (rME_rN_dr(n_a,l_a,n_b,l_b,N,hw) + &
                              rME_rN_dr(n_b,l_b,n_a,l_a,N,hw) - real(J,dp) * rME_rN(n_a,l_a,n_b,l_b,N-1,hw))
        if (L == J - 1) rME = rME + Amp * sqrt(real(J + 1,dp) / real(2*J + 1,dp)) * (rME_rN_dr(n_a,l_a,n_b,l_b,N,hw) + &
                              rME_rN_dr(n_b,l_b,n_a,l_a,N,hw) + real(J + 1,dp) * rME_rN(n_a,l_a,n_b,l_b,N-1,hw))
    end function rME_nabla_cross_S_dot_rN_YJL

    pure real(dp) function rME_rN_YL(a,b,N,L,hw,Orb) result(rME)
        integer, intent(in) :: a, b, N, L
        real(dp), intent(in) :: hw
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp) :: rY, rN
        rME = 0.0_dp
        if (mod(Orb(a)%l + Orb(b)%l + L,2) == 0) then
            rY = phase(L + (Orb(a)%j - 1)/2) * sqrt(real((Orb(a)%j + 1) * (Orb(b)%j + 1),dp) / (4.0_dp * pi)) * &
                 fCG(Orb(a)%j,Orb(b)%j,2*L,1,-1,0)
            rN = radial_moment_LHO(N,Orb(a)%n,Orb(a)%l,Orb(b)%n,Orb(b)%l,b_osc_of(hw))
            rME = rN * rY
        end if
    end function rME_rN_YL

    subroutine Tr1b_E1_NLO_initialize(Params,Orb,chR2,TrOp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: chR2
        type(tr1b_t), intent(inout) :: TrOp
        integer :: a_max, a, b, l_a, j_a, l_b, j_b, iso
        real(dp) :: hw, hc_mc2, rG_dot_r2_Y10, rG_dot_Y10, rG_dot_r2_Y12, rGxS_dot_Y10, rGxS_dot_r2_Y10, rGxS_dot_r2_Y12
        real(dp) :: rr1_Y1, rr3_Y1, VC, VS, TC, TS, sTC, sTS, CC, sC
        real(dp), parameter :: hc_pmc2 = HbarC / m_p, hc_nmc2 = HbarC / m_n

        hw = Params%Calc%hw
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(a)') TAB // 'Constructing 1-body E1 vortical transition operators in the LHO basis ...'
        call o1b_alloc(TrOp%E1_VC,a_max)
        call o1b_alloc(TrOp%E1_VS,a_max)
        call o1b_alloc(TrOp%E1_TC,a_max)
        call o1b_alloc(TrOp%E1_TS,a_max)
        call o1b_alloc(TrOp%E1_sTC,a_max)
        call o1b_alloc(TrOp%E1_sTS,a_max)
        call o1b_alloc(TrOp%E1_C,a_max)
        call o1b_alloc(TrOp%E1_sC,a_max)

        !$omp parallel do schedule(dynamic) private(l_a,j_a,b,l_b,j_b,rG_dot_r2_Y10,rG_dot_Y10,rG_dot_r2_Y12,rGxS_dot_Y10, &
        !$omp rGxS_dot_r2_Y10,rGxS_dot_r2_Y12,rr1_Y1,rr3_Y1,iso,hc_mc2,VC,VS,TC,TS,sTC,sTS,CC,sC)
        do a = 1, a_max
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do b = 1, a_max
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                ! E1 selection rule ...
                if (.not. (mod(l_a + l_b + 1,2) == 0 .and. abs(j_b - 2) <= j_a .and. j_a <= (j_b + 2))) cycle

                rG_dot_r2_Y10 = rME_nabla_dot_rN_YJL(a,b,2,0,1,hw,Orb)
                rG_dot_Y10 = rME_nabla_dot_rN_YJL(a,b,0,0,1,hw,Orb)
                rG_dot_r2_Y12 = rME_nabla_dot_rN_YJL(a,b,2,2,1,hw,Orb)
                rGxS_dot_Y10 = rME_nabla_cross_S_dot_rN_YJL(a,b,0,0,1,hw,Orb)
                rGxS_dot_r2_Y10 = rME_nabla_cross_S_dot_rN_YJL(a,b,2,0,1,hw,Orb)
                rGxS_dot_r2_Y12 = rME_nabla_cross_S_dot_rN_YJL(a,b,2,2,1,hw,Orb)
                rr1_Y1 = rME_rN_YL(a,b,1,1,hw,Orb)
                rr3_Y1 = rME_rN_YL(a,b,3,1,hw,Orb)

                ! Units & normalization (Kvasil et al., PRC 84, 034303 (2011)) ... convective parts (VC, TC, sTC) 0.5 hbar/mc
                    ! <a||A . (psi_a^+ grad psi_b - (grad psi_a)^+ psi_b)||b>, magnetization parts (VS, TS, sTS) hbar/mc
                    ! <a||A . curl(psi_a^+ s psi_b)||b>, combined as VC + 0.5 g_s VS downstream ...
                do iso = 1, 2
                    hc_mc2 = merge(hc_pmc2, hc_nmc2, iso == 1)
                    VC = 0.5_dp * hc_mc2 * sqrt(3.0_dp / 2.0_dp) / 5.0_dp * rG_dot_r2_Y12
                    VS = hc_mc2 * sqrt(3.0_dp / 2.0_dp) / 5.0_dp * rGxS_dot_r2_Y12
                    TC = 0.5_dp * hc_mc2 * (sqrt(3.0_dp / 2.0_dp) / 15.0_dp * rG_dot_r2_Y12 + sqrt(3.0_dp) / 6.0_dp * rG_dot_r2_Y10)
                    TS = hc_mc2 * (sqrt(3.0_dp / 2.0_dp) / 15.0_dp * rGxS_dot_r2_Y12 + sqrt(3.0_dp) / 6.0_dp * rGxS_dot_r2_Y10)
                    sTC = -0.5_dp * hc_mc2 * sqrt(3.0_dp) / 6.0_dp * chR2 * rG_dot_Y10
                    sTS = -hc_mc2 * sqrt(3.0_dp) / 6.0_dp * chR2 * rGxS_dot_Y10
                    CC = 1.0_dp / 10.0_dp * rr3_Y1
                    sC = -1.0_dp / 6.0_dp * chR2 * rr1_Y1
                    if (iso == 1) then
                        TrOp%E1_VC%p(a,b) = VC
                        TrOp%E1_VS%p(a,b) = VS
                        TrOp%E1_TC%p(a,b) = TC
                        TrOp%E1_TS%p(a,b) = TS
                        TrOp%E1_sTC%p(a,b) = sTC
                        TrOp%E1_sTS%p(a,b) = sTS
                        TrOp%E1_C%p(a,b) = CC
                        TrOp%E1_sC%p(a,b) = sC
                    else
                        TrOp%E1_VC%n(a,b) = VC
                        TrOp%E1_VS%n(a,b) = VS
                        TrOp%E1_TC%n(a,b) = TC
                        TrOp%E1_TS%n(a,b) = TS
                        TrOp%E1_sTC%n(a,b) = sTC
                        TrOp%E1_sTS%n(a,b) = sTS
                        TrOp%E1_C%n(a,b) = CC
                        TrOp%E1_sC%n(a,b) = sC
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        write(*,'(a)') TAB // 'E1 vortical 1-body transition operators succesfully constructed in the LHO basis ...'
    end subroutine Tr1b_E1_NLO_initialize

    ! Reduced MX matrix element ... T = -1 proton, T = 1 neutron ...
    pure real(dp) function rMX(n_a,l_a,j_a,n_b,l_b,j_b,X,T,b_osc)
        integer, intent(in) :: n_a, l_a, j_a, n_b, l_b, j_b, X, T
        real(dp), intent(in) :: b_osc
        real(dp), parameter :: mu_N = 0.10515_dp, g_p = 5.586_dp, g_n = -3.826_dp
        real(dp) :: kappa, Amp
        kappa = 0.5_dp * real(iphase((j_a+1)/2 + l_a) * (j_a + 1) + iphase((j_b+1)/2 + l_b) * (j_b + 1),dp)
        if (T == 1) then
            Amp = mu_N * phase(X + (j_a - 1)/2) * sqrt(real((j_a + 1) * (j_b + 1),dp)) * fCG(j_a,j_b,2*X,1,-1,0) * &
                  (real(X,dp) - kappa) * (-0.5_dp * g_n) / sqrt(4.0_dp * pi)
        else
            Amp = mu_N * phase(X + (j_a - 1)/2) * sqrt(real((j_a + 1) * (j_b + 1),dp)) * fCG(j_a,j_b,2*X,1,-1,0) * &
                  (real(X,dp) - kappa) * (1.0_dp + kappa / (real(X,dp) + 1.0_dp) - 0.5_dp * g_p) / sqrt(4.0_dp * pi)
        end if
        rMX = Amp * radial_moment_LHO(X,n_a,l_a,n_b,l_b,b_osc)
    end function rMX

    subroutine Tr1b_MX_initialize(Params,Orb,TrOp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(tr1b_t), intent(inout) :: TrOp
        integer :: a_max, a, b, n_a, l_a, j_a, n_b, l_b, j_b
        real(dp) :: b_osc

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        b_osc = b_osc_of(Params%Calc%hw)
        write(*,'(a)') TAB // 'Constructing 1-body MX transition operators in the LHO basis ...'
        call o1b_alloc(TrOp%M1,a_max)
        call o1b_alloc(TrOp%M2,a_max)
        call o1b_alloc(TrOp%M3,a_max)

        !$omp parallel do private(n_a,l_a,j_a,b,n_b,l_b,j_b)
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                if (mod(l_a + l_b + 2,2) == 0) then
                    TrOp%M1%p(a,b) = rMX(n_a,l_a,j_a,n_b,l_b,j_b,1,-1,b_osc)
                    TrOp%M1%n(a,b) = rMX(n_a,l_a,j_a,n_b,l_b,j_b,1,1,b_osc)
                end if
                if (mod(l_a + l_b + 3,2) == 0) then
                    TrOp%M2%p(a,b) = rMX(n_a,l_a,j_a,n_b,l_b,j_b,2,-1,b_osc)
                    TrOp%M2%n(a,b) = rMX(n_a,l_a,j_a,n_b,l_b,j_b,2,1,b_osc)
                end if
                if (mod(l_a + l_b + 4,2) == 0) then
                    TrOp%M3%p(a,b) = rMX(n_a,l_a,j_a,n_b,l_b,j_b,3,-1,b_osc)
                    TrOp%M3%n(a,b) = rMX(n_a,l_a,j_a,n_b,l_b,j_b,3,1,b_osc)
                end if
            end do
        end do
        !$omp end parallel do
        write(*,'(a)') TAB // 'MX 1-body transition operators succesfully constructed in the LHO basis ...'
    end subroutine Tr1b_MX_initialize

    ! Fermi & Gamow-Teller transition operators ... counterpart of Tr1b_BX_initialize ...
        ! (a||F||b) = sqrt(j_a+1) for n,l,j equal, (a||GT||b) = sqrt(2 (j_a+1)(j_b+1)) (-1)^(l_a+(j_a+3)/2) {1/2 1/2 1; j_b j_a l_a}
        ! for n_a = n_b, l_a = l_b, |j_a-j_b| <= 2 ...
    subroutine Tr1b_BX_initialize(Params,Orb,TrOp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(tr1b_t), intent(inout) :: TrOp
        integer :: a_max, a, b, n_a, l_a, j_a, n_b, l_b, j_b
        real(dp) :: rME

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(a)') TAB // 'Constructing 1-body Fermi & Gamow-Teller BX transition operators in the LHO basis ...'
        call o1b_alloc(TrOp%F,a_max)
        call o1b_alloc(TrOp%GT,a_max)

        !$omp parallel do private(n_a,l_a,j_a,b,n_b,l_b,j_b,rME)
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                if (l_a == l_b .and. j_a == j_b .and. n_a == n_b) then
                    rME = sqrt(real(j_a + 1,dp))
                    TrOp%F%p(a,b) = rME
                    TrOp%F%n(a,b) = rME
                end if
                if (l_a == l_b .and. n_a == n_b .and. abs(j_a - j_b) <= 2) then
                    rME = sqrt(2.0_dp * real(j_a + 1,dp) * real(j_b + 1,dp)) * phase(l_a + (j_a + 3)/2) * f6j(1,1,2,j_b,j_a,2*l_a)
                    TrOp%GT%p(a,b) = rME
                    TrOp%GT%n(a,b) = rME
                end if
            end do
        end do
        !$omp end parallel do
        write(*,'(a)') TAB // 'Fermi & Gamow-Teller BX 1-body transition operators succesfully constructed in the LHO basis ...'
    end subroutine Tr1b_BX_initialize

    ! Transformation of all transition operators to the target basis C ... F & GT mix the bases, F'%p = C%p^T F%p C%n,
        ! F'%n = C%n^T F%n C%p ... O' = C^T O C restricted to (l,j) blocks ...
    function Tr1b_transformation(Params,Orb,C,TrOp) result(TrOp_new)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: C
        type(tr1b_t), intent(in) :: TrOp
        type(tr1b_t) :: TrOp_new
        integer :: a_max

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(/,a)') 'Transforming 1-body transition operators from the reference basis to the given target basis ...'
        TrOp_new%E0 = transform(TrOp%E0)
        TrOp_new%E1 = transform(TrOp%E1)
        TrOp_new%E2 = transform(TrOp%E2)
        TrOp_new%E3 = transform(TrOp%E3)
        TrOp_new%E1_VC = transform(TrOp%E1_VC)
        TrOp_new%E1_VS = transform(TrOp%E1_VS)
        TrOp_new%E1_TC = transform(TrOp%E1_TC)
        TrOp_new%E1_TS = transform(TrOp%E1_TS)
        TrOp_new%E1_sTC = transform(TrOp%E1_sTC)
        TrOp_new%E1_sTS = transform(TrOp%E1_sTS)
        TrOp_new%E1_C = transform(TrOp%E1_C)
        TrOp_new%E1_sC = transform(TrOp%E1_sC)
        TrOp_new%M1 = transform(TrOp%M1)
        TrOp_new%M2 = transform(TrOp%M2)
        TrOp_new%M3 = transform(TrOp%M3)
        TrOp_new%F = transform_pn(TrOp%F)
        TrOp_new%GT = transform_pn(TrOp%GT)
        write(*,'(a)') TAB // 'Transition 1-body operators transformed to the given target basis ...'
    contains
        function transform(O) result(O_new)
            type(o1b_t), intent(in) :: O
            type(o1b_t) :: O_new
            integer :: a, b, k, l
            real(dp) :: pSum, nSum
            call o1b_alloc(O_new,a_max)
            !$omp parallel do private(b,pSum,nSum,k,l)
            do a = 1, a_max
                do b = 1, a_max
                    pSum = 0.0_dp
                    nSum = 0.0_dp
                    do k = 1, a_max
                        if (Orb(a)%l == Orb(k)%l .and. Orb(a)%j == Orb(k)%j) then
                            do l = 1, a_max
                                if (Orb(b)%l == Orb(l)%l .and. Orb(b)%j == Orb(l)%j) then
                                    pSum = pSum + O%p(k,l) * C%p(k,a) * C%p(l,b)
                                    nSum = nSum + O%n(k,l) * C%n(k,a) * C%n(l,b)
                                end if
                            end do
                        end if
                    end do
                    O_new%p(a,b) = pSum
                    O_new%n(a,b) = nSum
                end do
            end do
            !$omp end parallel do
        end function transform

        function transform_pn(O) result(O_new)
            type(o1b_t), intent(in) :: O
            type(o1b_t) :: O_new
            integer :: a, b, k, l
            real(dp) :: pSum, nSum
            call o1b_alloc(O_new,a_max)
            !$omp parallel do private(b,pSum,nSum,k,l)
            do a = 1, a_max
                do b = 1, a_max
                    pSum = 0.0_dp
                    nSum = 0.0_dp
                    do k = 1, a_max
                        if (Orb(a)%l == Orb(k)%l .and. Orb(a)%j == Orb(k)%j) then
                            do l = 1, a_max
                                if (Orb(b)%l == Orb(l)%l .and. Orb(b)%j == Orb(l)%j) then
                                    pSum = pSum + O%p(k,l) * C%p(k,a) * C%n(l,b)
                                    nSum = nSum + O%n(k,l) * C%n(k,a) * C%p(l,b)
                                end if
                            end do
                        end if
                    end do
                    O_new%p(a,b) = pSum
                    O_new%n(a,b) = nSum
                end do
            end do
            !$omp end parallel do
        end function transform_pn
    end function Tr1b_transformation

    ! Transition operators in the quasiparticle representation (11 & 20 parts) ...
    function qpTr1b_initialize(Params,Orb,chR2,C,U,V) result(qpTrOp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: chR2
        type(o1b_t), intent(in) :: C, U, V
        type(qptr1b_t) :: qpTrOp
        type(tr1b_t) :: TrOp

        write(*,'(/,a)') 'Preparing 1-body electromagnetic & weak transition operators in quasiparticle representation ...'
        TrOp = Tr1b_initialize(Params,Orb,chR2)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        ! The current operators VC, VS, TC, TS, sTC & sTS are time-odd (T_odd = .true.) ...
        qpTrOp%E0 = qp_transform(TrOp%E0,0,.true.,.false.)
        qpTrOp%E1 = qp_transform(TrOp%E1,1,.true.,.false.)
        qpTrOp%E2 = qp_transform(TrOp%E2,2,.true.,.false.)
        qpTrOp%E3 = qp_transform(TrOp%E3,3,.true.,.false.)
        qpTrOp%E1_VC = qp_transform(TrOp%E1_VC,1,.true.,.true.)
        qpTrOp%E1_VS = qp_transform(TrOp%E1_VS,1,.true.,.true.)
        qpTrOp%E1_TC = qp_transform(TrOp%E1_TC,1,.true.,.true.)
        qpTrOp%E1_TS = qp_transform(TrOp%E1_TS,1,.true.,.true.)
        qpTrOp%E1_C = qp_transform(TrOp%E1_C,1,.true.,.false.)
        qpTrOp%E1_sTC = qp_transform(TrOp%E1_sTC,1,.true.,.true.)
        qpTrOp%E1_sTS = qp_transform(TrOp%E1_sTS,1,.true.,.true.)
        qpTrOp%E1_sC = qp_transform(TrOp%E1_sC,1,.true.,.false.)
        qpTrOp%M1 = qp_transform(TrOp%M1,1,.false.,.false.)
        qpTrOp%M2 = qp_transform(TrOp%M2,2,.false.,.false.)
        qpTrOp%M3 = qp_transform(TrOp%M3,3,.false.,.false.)
        qpTrOp%F = qp_transform_bx(TrOp%F,0)
        qpTrOp%GT = qp_transform_bx(TrOp%GT,1)
        write(*,'(a)') TAB // 'Electromagnetic & weak 1-body transition operators succesfully evaluated in the quasiparticle representation ...'

    contains

        ! qpTr1b_EX_initialize (Electric = .true.) & qpTr1b_MX_initialize (Electric = .false.) ... dev forms:
            !   O^11_ab = sum_kl [O_kl U_ka U_lb + (-1)^((j_a+j_b)/2+X) O_lk V_ka V_lb] (no c-number term for X = 0),
            !   O^20_ab = sum_kl O_kl [Phase_T V_ka U_lb + (-1)^X U_ka V_lb] (EX, X > 0; Phase_T = -1 for time-odd operators),
            !   O^20_ab = sum_kl O_kl [V_ka U_lb + U_ka V_lb] (E0), O^20_ab = sum_kl O_kl [V_ka U_lb - (-1)^X U_ka V_lb] (MX) ...
        function qp_transform(O,X,Electric,T_odd) result(qpO)
            type(o1b_t), intent(in) :: O
            integer, intent(in) :: X
            logical, intent(in) :: Electric, T_odd
            type(qpo1b_t) :: qpO
            integer :: a_max, a, b, k, l, l_a, j_a, l_b, j_b, Par_req
            real(dp) :: Phase_11, Phase_20, Phase_T, p11, n11, p20, n20, pO_kl, nO_kl, pO_lk, nO_lk
            real(dp) :: pU_ka, pV_ka, pU_lb, pV_lb, nU_ka, nV_ka, nU_lb, nV_lb

            a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
            call o1b_alloc(qpO%qp11,a_max)
            call o1b_alloc(qpO%qp20,a_max)
            Par_req = merge(0, 1, Electric)
            Phase_T = merge(-1.0_dp, 1.0_dp, T_odd)

            !$omp parallel do private(l_a,j_a,b,l_b,j_b,p11,n11,p20,n20,k,l,Phase_11,Phase_20, &
            !$omp pU_ka,pV_ka,pU_lb,pV_lb,nU_ka,nV_ka,nU_lb,nV_lb,pO_kl,nO_kl,pO_lk,nO_lk)
            do a = 1, a_max
                l_a = Orb(a)%l
                j_a = Orb(a)%j
                do b = 1, a_max
                    l_b = Orb(b)%l
                    j_b = Orb(b)%j
                    if (.not. ((abs(j_a - j_b)/2 <= X) .and. (X <= (j_a + j_b)/2) .and. (mod(l_a + l_b + X,2) == Par_req))) cycle
                    p11 = 0.0_dp
                    n11 = 0.0_dp
                    p20 = 0.0_dp
                    n20 = 0.0_dp
                    do k = 1, a_max
                        if (.not. (l_a == Orb(k)%l .and. j_a == Orb(k)%j)) cycle
                        do l = 1, a_max
                            if (.not. (Orb(l)%l == l_b .and. Orb(l)%j == j_b)) cycle
                            Phase_11 = phase((j_a + j_b)/2 + X)
                            Phase_20 = phase(X)
                            pU_ka = U%p(k,a)
                            pV_ka = V%p(k,a)
                            pU_lb = U%p(l,b)
                            pV_lb = V%p(l,b)
                            nU_ka = U%n(k,a)
                            nV_ka = V%n(k,a)
                            nU_lb = U%n(l,b)
                            nV_lb = V%n(l,b)
                            pO_kl = O%p(k,l)
                            nO_kl = O%n(k,l)
                            pO_lk = O%p(l,k)
                            nO_lk = O%n(l,k)
                            p11 = p11 + (pO_kl * pU_ka * pU_lb + Phase_11 * pO_lk * pV_ka * pV_lb)
                            n11 = n11 + (nO_kl * nU_ka * nU_lb + Phase_11 * nO_lk * nV_ka * nV_lb)
                            if (.not. Electric) then
                                p20 = p20 + pO_kl * (pV_ka * pU_lb - Phase_20 * pU_ka * pV_lb)
                                n20 = n20 + nO_kl * (nV_ka * nU_lb - Phase_20 * nU_ka * nV_lb)
                            else if (X == 0) then
                                p20 = p20 + pO_kl * (pV_ka * pU_lb + pU_ka * pV_lb)
                                n20 = n20 + nO_kl * (nV_ka * nU_lb + nU_ka * nV_lb)
                            else
                                p20 = p20 + pO_kl * (Phase_T * pV_ka * pU_lb + Phase_20 * pU_ka * pV_lb)
                                n20 = n20 + nO_kl * (Phase_T * nV_ka * nU_lb + Phase_20 * nU_ka * nV_lb)
                            end if
                        end do
                    end do
                    qpO%qp11%p(a,b) = p11
                    qpO%qp11%n(a,b) = n11
                    qpO%qp20%p(a,b) = p20
                    qpO%qp20%n(a,b) = n20
                end do
            end do
            !$omp end parallel do
        end function qp_transform

        ! Fermi (X = 0) & Gamow-Teller (X = 1) in the quasiparticle representation ... counterpart of qpTr1b_BX_initialize ...
            ! UU from the own block, VV from the opposite (transposed) block ...
        function qp_transform_bx(O,X) result(qpO)
            type(o1b_t), intent(in) :: O
            integer, intent(in) :: X
            type(qpo1b_t) :: qpO
            integer :: a_max, a, b, k, l, l_a, j_a, l_b, j_b
            real(dp) :: Phase_11, p11, n11, p20, n20
            real(dp) :: pU_ka, pV_ka, pU_lb, pV_lb, nU_ka, nV_ka, nU_lb, nV_lb

            a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
            call o1b_alloc(qpO%qp11,a_max)
            call o1b_alloc(qpO%qp20,a_max)

            !$omp parallel do private(l_a,j_a,b,l_b,j_b,p11,n11,p20,n20,k,l,Phase_11, &
            !$omp pU_ka,pV_ka,pU_lb,pV_lb,nU_ka,nV_ka,nU_lb,nV_lb)
            do a = 1, a_max
                l_a = Orb(a)%l
                j_a = Orb(a)%j
                do b = 1, a_max
                    l_b = Orb(b)%l
                    j_b = Orb(b)%j
                    if (.not. ((abs(j_a - j_b)/2 <= X) .and. (X <= (j_a + j_b)/2) .and. (mod(l_a + l_b,2) == 0))) cycle
                    Phase_11 = phase((j_a + j_b)/2 + X)
                    p11 = 0.0_dp
                    n11 = 0.0_dp
                    p20 = 0.0_dp
                    n20 = 0.0_dp
                    do k = 1, a_max
                        if (.not. (l_a == Orb(k)%l .and. j_a == Orb(k)%j)) cycle
                        pU_ka = U%p(k,a)
                        pV_ka = V%p(k,a)
                        nU_ka = U%n(k,a)
                        nV_ka = V%n(k,a)
                        do l = 1, a_max
                            if (.not. (Orb(l)%l == l_b .and. Orb(l)%j == j_b)) cycle
                            pU_lb = U%p(l,b)
                            pV_lb = V%p(l,b)
                            nU_lb = U%n(l,b)
                            nV_lb = V%n(l,b)
                            p11 = p11 + (O%p(k,l) * pU_ka * nU_lb + Phase_11 * O%n(l,k) * pV_ka * nV_lb)
                            n11 = n11 + (O%n(k,l) * nU_ka * pU_lb + Phase_11 * O%p(l,k) * nV_ka * pV_lb)
                            p20 = p20 + O%p(k,l) * pU_ka * nV_lb
                            n20 = n20 + O%n(k,l) * nU_ka * pV_lb
                        end do
                    end do
                    qpO%qp11%p(a,b) = p11
                    qpO%qp11%n(a,b) = n11
                    qpO%qp20%p(a,b) = p20
                    qpO%qp20%n(a,b) = n20
                end do
            end do
            !$omp end parallel do
        end function qp_transform_bx

    end function qpTr1b_initialize

end module nj_tr1b
