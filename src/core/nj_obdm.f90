! Radial densities & radii from the 1-body density matrix ... counterpart of src/common/rho.jl (NuJuliet dev;
    ! MatrixElements/OBDM.jl before, OBDM_* renamed rho_*) ...
    ! Rho is given in the LHO basis, C is the transformation from the LHO to the reference basis ...
module nj_obdm
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : fCG
    use nj_functions, only : Psi_rad_LHO, integrate_trap, radial_moment_matrix_LHO
    use nj_linalg, only : mtmm
    implicit none
    private

    public :: rho_export, rho_radial_cms, rho_radial_density_grid, rho_rN, rho_chR2
    public :: rho_radial_chdensity_grid, rho_radial_radii, rho_radial_density_export, rho_radial_summary
    ! LEGACY (removed in cleanup stage) ... the former OBDM_* names, still used by the not-yet-ported solvers ...
    public :: OBDM_export, OBDM_radial_cms, OBDM_radial_density_grid, OBDM_rN, OBDM_chR2

    interface OBDM_export
        module procedure rho_export
    end interface OBDM_export
    interface OBDM_radial_cms
        module procedure rho_radial_cms
    end interface OBDM_radial_cms
    interface OBDM_radial_density_grid
        module procedure rho_radial_density_grid
    end interface OBDM_radial_density_grid
    interface OBDM_rN
        module procedure rho_rN
    end interface OBDM_rN
    interface OBDM_chR2
        module procedure rho_chR2
    end interface OBDM_chR2

    integer, parameter :: N_Sampling = 2**13 + 1

contains

    subroutine rho_export(Params,Orb,Summary_File,Densities_File,Rho_LHO,C)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: Summary_File, Densities_File
        type(o1b_t), intent(in) :: Rho_LHO, C
        type(o1b_t) :: Rho
        real(dp) :: R_CMS(4), kR2, pR2, nR2, chR2
        real(dp), allocatable :: r_grid(:), pRho_grid(:), nRho_grid(:), chRho_grid(:)

        ! Transform Rho from LHO to the reference basis ...
        Rho%p = mtmm(C%p,Rho_LHO%p)
        Rho%n = mtmm(C%n,Rho_LHO%n)

        if (Params%Calc%CMS == 'CMS1+2B' .or. Params%Calc%CMS == 'CMS2B') then
            R_CMS = rho_radial_cms(Params,Rho,C,Orb)
        else
            R_CMS = 0.0_dp
        end if

        call rho_radial_density_grid(Params,Rho,C,Orb,r_grid,pRho_grid,nRho_grid)
        call rho_radial_chdensity_grid(Params,Rho,C,Orb,R_CMS,r_grid,pRho_grid,nRho_grid,chRho_grid,kR2)
        call rho_radial_radii(r_grid,pRho_grid,nRho_grid,chRho_grid,pR2,nR2,chR2)
        call rho_radial_density_export(Densities_File,r_grid,pRho_grid,nRho_grid,chRho_grid)
        call rho_radial_summary(Params,Summary_File,pR2,nR2,chR2,R_CMS,kR2)

        write(*,'(/,a)') 'Radial densities and corresponding radii sucessfully exported ...'
    end subroutine rho_export

    function rho_radial_cms(Params,Rho,C,Orb) result(R_CMS)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, C
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp) :: R_CMS(4)
        real(dp) :: hw, A_mass, Z, Amp, pME, nME, pR2_CMS1, pR2_CMS2, nR2_CMS1, nR2_CMS2
        real(dp) :: pR1Sum, pR2Sum, nR1Sum, nR2Sum
        real(dp), allocatable :: Radial_r1_LHO(:,:), Radial_r2_LHO(:,:), pRadial_r1(:,:), pRadial_r2(:,:), &
                                 nRadial_r1(:,:), nRadial_r2(:,:)
        integer :: a_max, a, b, k, l, l_a, j_a, l_b, j_b

        hw = Params%Calc%hw
        A_mass = real(Params%Calc%A,dp)
        Z = real(Params%Calc%Z,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2

        Radial_r1_LHO = radial_moment_matrix_LHO(1,hw,Orb(1:a_max))
        Radial_r2_LHO = radial_moment_matrix_LHO(2,hw,Orb(1:a_max))
        allocate(pRadial_r1(a_max,a_max), pRadial_r2(a_max,a_max), nRadial_r1(a_max,a_max), nRadial_r2(a_max,a_max))

        !$omp parallel do private(l_a,j_a,b,l_b,j_b,pR1Sum,pR2Sum,nR1Sum,nR2Sum,k,l)
        do a = 1, a_max
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do b = 1, a_max
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                pR1Sum = 0.0_dp
                pR2Sum = 0.0_dp
                nR1Sum = 0.0_dp
                nR2Sum = 0.0_dp
                do k = 1, a_max
                    if ((Orb(k)%l == l_a) .and. (Orb(k)%j == j_a)) then
                        do l = 1, a_max
                            if ((Orb(l)%l == l_b) .and. (Orb(l)%j == j_b)) then
                                pR1Sum = pR1Sum + C%p(k,a) * C%p(l,b) * Radial_r1_LHO(k,l)
                                pR2Sum = pR2Sum + C%p(k,a) * C%p(l,b) * Radial_r2_LHO(k,l)
                                nR1Sum = nR1Sum + C%n(k,a) * C%n(l,b) * Radial_r1_LHO(k,l)
                                nR2Sum = nR2Sum + C%n(k,a) * C%n(l,b) * Radial_r2_LHO(k,l)
                            end if
                        end do
                    end if
                end do
                pRadial_r1(a,b) = pR1Sum
                pRadial_r2(a,b) = pR2Sum
                nRadial_r1(a,b) = nR1Sum
                nRadial_r2(a,b) = nR2Sum
            end do
        end do
        !$omp end parallel do

        pR2_CMS1 = 0.0_dp
        pR2_CMS2 = 0.0_dp
        nR2_CMS1 = 0.0_dp
        nR2_CMS2 = 0.0_dp
        do a = 1, a_max
            j_a = Orb(a)%j
            do b = 1, a_max
                j_b = Orb(b)%j
                Amp = fCG(j_a,j_b,2,1,-1,0) * fCG(j_b,j_a,2,1,-1,0)
                pME = Amp * ((Rho%p(a,a) * Rho%p(b,b) * abs(pRadial_r1(a,b))**2) * (2.0_dp / (Z*A_mass) - 1.0_dp / A_mass**2) - &
                      1.0_dp / A_mass**2 * Rho%n(a,a) * Rho%n(b,b) * abs(nRadial_r1(a,b))**2)
                nME = Amp * ((Rho%n(a,a) * Rho%n(b,b) * abs(nRadial_r1(a,b))**2) * (2.0_dp / (Z*A_mass) - 1.0_dp / A_mass**2) - &
                      1.0_dp / A_mass**2 * Rho%p(a,a) * Rho%p(b,b) * abs(pRadial_r1(a,b))**2)
                pR2_CMS2 = pR2_CMS2 + pME
                nR2_CMS2 = nR2_CMS2 + nME
            end do
            pME = ((1.0_dp / A_mass**2 - 2.0_dp / (Z*A_mass)) * Rho%p(a,a) * pRadial_r2(a,a) + Rho%n(a,a) * nRadial_r2(a,a) / A_mass**2) * &
                  real(j_a + 1,dp)
            nME = ((1.0_dp / A_mass**2 - 2.0_dp / ((A_mass-Z)*A_mass)) * Rho%n(a,a) * nRadial_r2(a,a) + Rho%p(a,a) * pRadial_r2(a,a) / A_mass**2) * &
                  real(j_a + 1,dp)
            pR2_CMS1 = pR2_CMS1 + pME
            nR2_CMS1 = nR2_CMS1 + nME
        end do

        R_CMS = [pR2_CMS1, pR2_CMS2, nR2_CMS1, nR2_CMS2]
    end function rho_radial_cms

    ! Julia collect(range(r1, stop = r2, length = N)) ...
    function radial_grid(r1,r2,N) result(r_grid)
        real(dp), intent(in) :: r1, r2
        integer, intent(in) :: N
        real(dp), allocatable :: r_grid(:)
        integer :: i
        allocate(r_grid(N))
        do i = 1, N
            ! Julia's LinRange-like StepRangeLen evaluation ... (r1*(N-i) + r2*(i-1)) / (N-1) ...
            r_grid(i) = (r1 * real(N - i,dp) + r2 * real(i - 1,dp)) / real(N - 1,dp)
        end do
    end function radial_grid

    subroutine rho_radial_density_grid(Params,Rho,C,Orb,r_grid,pRho_grid,nRho_grid)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, C
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), allocatable, intent(out) :: r_grid(:), pRho_grid(:), nRho_grid(:)
        real(dp) :: hw, A_mass, nu_proton, nu_neutron, r, pSum, nSum, pRad, nRad
        integer :: a_max, i, a, k

        hw = Params%Calc%hw
        A_mass = real(Params%Calc%A,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        nu_proton = 0.5_dp * m_p * hw / HbarC**2
        nu_neutron = 0.5_dp * m_n * hw / HbarC**2

        r_grid = radial_grid(0.0_dp + 1.0e-8_dp, 2.5_dp * 1.2_dp * A_mass**(1.0_dp/3.0_dp), N_Sampling)
        allocate(pRho_grid(N_Sampling), nRho_grid(N_Sampling))

        !$omp parallel do private(r,pSum,nSum,a,pRad,nRad,k)
        do i = 1, N_Sampling
            r = r_grid(i)
            pSum = 0.0_dp
            nSum = 0.0_dp
            do a = 1, a_max
                pRad = 0.0_dp
                nRad = 0.0_dp
                do k = 1, a_max
                    if ((Orb(a)%j == Orb(k)%j) .and. (Orb(a)%l == Orb(k)%l)) then
                        pRad = pRad + Psi_rad_LHO(r,Orb(k)%n,Orb(k)%l,nu_proton) * C%p(k,a)
                        nRad = nRad + Psi_rad_LHO(r,Orb(k)%n,Orb(k)%l,nu_neutron) * C%n(k,a)
                    end if
                end do
                pSum = pSum + Rho%p(a,a) * pRad**2 * real(Orb(a)%j + 1,dp)
                nSum = nSum + Rho%n(a,a) * nRad**2 * real(Orb(a)%j + 1,dp)
            end do
            pRho_grid(i) = pSum / (4.0_dp * pi)
            nRho_grid(i) = nSum / (4.0_dp * pi)
        end do
        !$omp end parallel do
    end subroutine rho_radial_density_grid

    subroutine rho_radial_chdensity_grid(Params,Rho,C,Orb,R_CMS,r_grid,pRho_grid,nRho_grid,chRho_grid,kR2)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, C
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: R_CMS(4), r_grid(:), pRho_grid(:), nRho_grid(:)
        real(dp), allocatable, intent(out) :: chRho_grid(:)
        real(dp), intent(out) :: kR2
        real(dp), parameter :: k_p = 1.793_dp, k_n = 1.913_dp
        real(dp) :: hw, Z, nu_proton, nu_neutron, pR_CMS, nR_CMS, r, chSum, pRad, nRad, SO, chME, kME
        integer :: a_max, i, a, k, l_a, j_a

        hw = Params%Calc%hw
        Z = real(Params%Calc%Z,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        nu_proton = 0.5_dp * m_p * hw / HbarC**2
        nu_neutron = 0.5_dp * m_n * hw / HbarC**2
        pR_CMS = R_CMS(1) + R_CMS(2)
        nR_CMS = R_CMS(3) + R_CMS(4)

        allocate(chRho_grid(N_Sampling))
        chRho_grid = 0.0_dp
        kR2 = 0.0_dp

        !$omp parallel do schedule(dynamic,16) private(r,chSum,a,l_a,j_a,pRad,nRad,k,SO,chME)
        do i = 1, N_Sampling
            r = r_grid(i)
            chSum = 0.0_dp
            do a = 1, a_max
                l_a = Orb(a)%l
                j_a = Orb(a)%j
                pRad = 0.0_dp
                nRad = 0.0_dp
                do k = 1, a_max
                    if ((j_a == Orb(k)%j) .and. (l_a == Orb(k)%l)) then
                        pRad = pRad + Psi_rad_LHO(r,Orb(k)%n,Orb(k)%l,nu_proton) * C%p(k,a)
                        nRad = nRad + Psi_rad_LHO(r,Orb(k)%n,Orb(k)%l,nu_neutron) * C%n(k,a)
                    end if
                end do
                ! Spin-orbit term ...
                SO = real(l_a * kronecker_delta(j_a,2*l_a+1) - (l_a+1) * kronecker_delta(j_a,2*l_a-1),dp)
                chME = (k_n / m_n**2 * Rho%n(a,a) * nRad**2 + k_p / m_p**2 * Rho%p(a,a) * pRad**2) * &
                       HbarC**2 * SO * real(j_a + 1,dp) / Z
                chSum = chSum + chME
            end do
            chSum = chSum + fold_pRho(r,r_grid,pRho_grid,pR_CMS)
            chSum = chSum + fold_nRho(r,r_grid,nRho_grid,nR_CMS)
            chRho_grid(i) = chRho_grid(i) + chSum
        end do
        !$omp end parallel do

        ! Anomalous magnetic moment contribution to the charge radius ...
        do a = 1, a_max
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            SO = real(l_a * kronecker_delta(j_a,2*l_a+1) - (l_a+1) * kronecker_delta(j_a,2*l_a-1),dp)
            kME = HbarC**2 * (k_p / m_p**2 * Rho%p(a,a) + k_n / m_n**2 * Rho%n(a,a)) * real(j_a + 1,dp)**2 * SO / (4.0_dp * pi * Z)
            kR2 = kR2 + kME
        end do
    end subroutine rho_radial_chdensity_grid

    subroutine rho_radial_radii(r_grid,pRho_grid,nRho_grid,chRho_grid,pR2,nR2,chR2)
        real(dp), intent(in) :: r_grid(:), pRho_grid(:), nRho_grid(:), chRho_grid(:)
        real(dp), intent(out) :: pR2, nR2, chR2
        pR2 = integrate_trap(r_grid, r_grid**4 * pRho_grid) / integrate_trap(r_grid, r_grid**2 * pRho_grid)
        nR2 = integrate_trap(r_grid, r_grid**4 * nRho_grid) / integrate_trap(r_grid, r_grid**2 * nRho_grid)
        chR2 = integrate_trap(r_grid, r_grid**4 * chRho_grid) / integrate_trap(r_grid, r_grid**2 * chRho_grid)
    end subroutine rho_radial_radii

    subroutine rho_radial_density_export(Densities_File,r_grid,pRho,nRho,chRho)
        character(len=*), intent(in) :: Densities_File
        real(dp), intent(in) :: r_grid(:), pRho(:), nRho(:), chRho(:)
        integer :: u, i, N
        N = min(size(r_grid), size(pRho), size(nRho), size(chRho))
        open(newunit=u, file=Densities_File, status='replace', action='write')
        write(u,'(a)') cfmt_s('r',14,.true.) // ' ' // cfmt_s('pRho',14,.true.) // ' ' // cfmt_s('nRho',14,.true.) // ' ' // &
                       cfmt_s('chRho',14,.true.)
        do i = 1, N
            write(u,'(a)') cfmt_f(r_grid(i),14,8,.true.) // ' ' // cfmt_f(pRho(i),14,8,.true.) // ' ' // &
                           cfmt_f(nRho(i),14,8,.true.) // ' ' // cfmt_f(chRho(i),14,8,.true.)
        end do
        close(u)
    end subroutine rho_radial_density_export

    subroutine rho_radial_summary(Params,Summary_File,pR2,nR2,chR2,R_CMS,kR2)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Summary_File
        real(dp), intent(in) :: pR2, nR2, chR2, R_CMS(4), kR2
        real(dp) :: A_mass, Z, N, rpDF, rnDF, A13
        integer :: u, Pass

        A_mass = real(Params%Calc%A,dp)
        Z = real(Params%Calc%Z,dp)
        N = real(Params%Calc%A - Params%Calc%Z,dp)
        rpDF = 0.5_dp * (HbarC / m_p)**2
        rnDF = 0.5_dp * (HbarC / m_n)**2

        write(*,'(/,a)') 'Resulting point-proton, charged & point-neutron radii ...'
        write(*,'(a)') TAB // 'r_p  = ' // jstr(round_sig(sqrt(pR2 + R_CMS(1) + R_CMS(2) + rpDF),6)) // ' fm ' // TAB // TAB // &
                       ' ... ' // TAB // ' Point-proton radius'
        write(*,'(a)') TAB // 'r_n  = ' // jstr(round_sig(sqrt(nR2 + R_CMS(3) + R_CMS(4) + rnDF),6)) // ' fm ' // TAB // TAB // &
                       ' ... ' // TAB // ' Point-neutron radius'
        write(*,'(a)') TAB // 'r_ch = ' // jstr(round_sig(sqrt(pR2 + R_CMS(1) + R_CMS(2) - 0.106_dp * N / Z + 0.8414_dp - 0.143_dp),6)) // &
                       ' fm ' // TAB // TAB // ' ... ' // TAB // ' Standard charge radius'
        write(*,'(a)') TAB // 'r_ch = ' // jstr(round_sig(sqrt(chR2 + R_CMS(1) + R_CMS(2) + rpDF),6)) // ' fm ' // TAB // TAB // &
                       ' ... ' // TAB // ' Folded charge radius'
        write(*,'(a)') "For more details see the summary file ... ''" // Summary_File // "'' ..."

        write(*,'(/,a)') 'Exporting information on radii ...'
        open(newunit=u, file=Summary_File, position='append', action='write')
        do Pass = 1, 2
            if (Pass == 1) then
                write(u,'(/,/,a)') 'Resulting nuclear radii ...'
                A13 = 1.0_dp
            else
                ! dev (rho.jl:276-282) normalizes the radius, sqrt(R2) / A^1/3 (sqrt(R2 / A^1/3) before) ...
                write(u,'(/,/,a)') 'Resulting nuclear radii normalized by the scale factor 1/A^1/3 ...'
                A13 = A_mass**(1.0_dp/3.0_dp)
            end if
            write(u,'(/,a)') 'r_p   = ' // cfmt_f(sqrt(abs(pR2 + R_CMS(1) + R_CMS(2) + rpDF)) / A13,12,6) // ' fm ' // TAB // TAB // &
                             ' ... ' // TAB // ' Point-proton radius'
            write(u,'(a)') 'r_n   = ' // cfmt_f(sqrt(abs(nR2 + R_CMS(3) + R_CMS(4) + rnDF)) / A13,12,6) // ' fm ' // TAB // TAB // &
                           ' ... ' // TAB // ' Point-neutron radius'
            write(u,'(a)') 'r_ch  = ' // cfmt_f(sqrt(abs(pR2 + R_CMS(1) + R_CMS(2) - 0.106_dp * N / Z + 0.8414_dp - 0.143_dp)) / A13,12,6) // &
                           ' fm ' // TAB // TAB // ' ... ' // TAB // ' Standard charge radius'
            write(u,'(a)') 'r_ch  = ' // cfmt_f(sqrt(abs(chR2 + R_CMS(1) + R_CMS(2) + rpDF)) / A13,12,6) // ' fm ' // TAB // TAB // &
                           ' ... ' // TAB // ' Folded charge radius'
            write(u,'(/,a)') 'r_p  = ' // cfmt_f(sqrt(abs(pR2)) / A13,12,6) // ' fm ' // TAB // TAB // ' ... ' // TAB // &
                             ' Uncorrected point-proton radius'
            write(u,'(a)') 'r_ch = ' // cfmt_f(sqrt(abs(chR2)) / A13,12,6) // ' fm ' // TAB // TAB // ' ... ' // TAB // &
                           ' Uncorrected folded charge radius'
            write(u,'(a)') 'r_n  = ' // cfmt_f(sqrt(abs(nR2)) / A13,12,6) // ' fm ' // TAB // TAB // ' ... ' // TAB // &
                           ' Uncorrected point-neutron radius'
            write(u,'(/,a,/)') 'Details on radii corrections ...'
            write(u,'(a)') 'r_p2_cm1^2 = ' // cfmt_f(R_CMS(1),12,6) // ' fm^2'
            write(u,'(a)') 'r_p2_cm2^2 = ' // cfmt_f(R_CMS(2),12,6) // ' fm^2'
            write(u,'(a)') 'r_n2_cm1^2 = ' // cfmt_f(R_CMS(3),12,6) // ' fm^2'
            write(u,'(a)') 'r_n2_cm2^2 = ' // cfmt_f(R_CMS(4),12,6) // ' fm^2'
            write(u,'(a)') 'r_ch2_k^2  = ' // cfmt_f(kR2,12,6) // ' fm^2'
            write(u,'(a)') 'r_pDF2^2   = ' // cfmt_f(rpDF,12,6) // ' fm^2'
            write(u,'(a)') 'r_nDF2^2   = ' // cfmt_f(rnDF,12,6) // ' fm^2'
        end do
        close(u)
    end subroutine rho_radial_summary

    ! Folding with the finite proton & neutron charge distributions ...
    real(dp) function fold_pRho(r,r_grid,pRho_grid,pR_CMS)
        real(dp), intent(in) :: r, r_grid(:), pRho_grid(:), pR_CMS
        real(dp) :: Rho(size(r_grid))
        integer :: i
        do i = 1, size(r_grid)
            Rho(i) = convolution_p(r_grid(i),r,pR_CMS) * pRho_grid(i)
        end do
        fold_pRho = integrate_trap(r_grid,Rho)
    end function fold_pRho

    real(dp) function fold_nRho(r,r_grid,nRho_grid,nR_CMS)
        real(dp), intent(in) :: r, r_grid(:), nRho_grid(:), nR_CMS
        real(dp) :: Rho(size(r_grid))
        integer :: i
        do i = 1, size(r_grid)
            Rho(i) = convolution_n(r_grid(i),r,nR_CMS) * nRho_grid(i)
        end do
        fold_nRho = integrate_trap(r_grid,Rho)
    end function fold_nRho

    pure real(dp) function convolution_p(x,r,pR_CMS)
        real(dp), intent(in) :: x, r, pR_CMS
        real(dp), parameter :: a_1 = 0.506373_dp, a_2 = 0.327922_dp, a_3 = 0.165705_dp
        real(dp) :: r_1, r_2, r_3
        r_1 = sqrt(abs(0.431566_dp + pR_CMS + 0.5_dp * (HbarC / m_p)**2))
        r_2 = sqrt(abs(0.139140_dp + pR_CMS + 0.5_dp * (HbarC / m_p)**2))
        r_3 = sqrt(abs(1.525540_dp + pR_CMS + 0.5_dp * (HbarC / m_p)**2))
        convolution_p = x / (r * sqrt(pi)) * (a_1 / r_1 * (exp(-((r-x)/r_1)**2) - exp(-((r+x)/r_1)**2)) + &
                                              a_2 / r_2 * (exp(-((r-x)/r_2)**2) - exp(-((r+x)/r_2)**2)) + &
                                              a_3 / r_3 * (exp(-((r-x)/r_3)**2) - exp(-((r+x)/r_3)**2)))
    end function convolution_p

    pure real(dp) function convolution_n(x,r,nR_CMS)
        real(dp), intent(in) :: x, r, nR_CMS
        real(dp) :: r_p, r_m
        r_p = sqrt(abs(0.4828_dp - 0.038664_dp + nR_CMS + 0.5_dp * (HbarC / m_n)**2))
        r_m = sqrt(abs(0.4828_dp + 0.038664_dp + nR_CMS + 0.5_dp * (HbarC / m_n)**2))
        convolution_n = x / (r * sqrt(pi)) * ((exp(-((r-x)/r_p)**2) - exp(-((r+x)/r_p)**2)) / r_p - &
                                              (exp(-((r-x)/r_m)**2) - exp(-((r+x)/r_m)**2)) / r_m)
    end function convolution_n

    real(dp) function rho_rN(Params,N,Rho,C,Orb)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N
        type(o1b_t), intent(in) :: Rho, C
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), allocatable :: r_grid(:), pRho_grid(:)
        real(dp) :: hw, A_mass, nu_proton, r, pSum, pRad
        integer :: a_max, i, a, k

        hw = Params%Calc%hw
        A_mass = real(Params%Calc%A,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        nu_proton = 0.5_dp * m_p * hw / HbarC**2
        r_grid = radial_grid(0.0_dp + 1.0e-8_dp, 2.5_dp * 1.2_dp * A_mass**(1.0_dp/3.0_dp), N_Sampling)
        allocate(pRho_grid(N_Sampling))

        !$omp parallel do private(r,pSum,a,pRad,k)
        do i = 1, N_Sampling
            r = r_grid(i)
            pSum = 0.0_dp
            do a = 1, a_max
                pRad = 0.0_dp
                do k = 1, a_max
                    if ((Orb(a)%j == Orb(k)%j) .and. (Orb(a)%l == Orb(k)%l)) then
                        pRad = pRad + Psi_rad_LHO(r,Orb(k)%n,Orb(k)%l,nu_proton) * C%p(k,a)
                    end if
                end do
                pSum = pSum + real(Orb(a)%j + 1,dp) * Rho%p(a,a) * pRad**2
            end do
            pRho_grid(i) = pSum / (4.0_dp * pi)
        end do
        !$omp end parallel do

        rho_rN = integrate_trap(r_grid, r_grid**(N+2) * pRho_grid) / integrate_trap(r_grid, r_grid**2 * pRho_grid)
    end function rho_rN

    real(dp) function rho_chR2(Params,Orb,C,Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: C, Rho
        real(dp) :: R_CMS(4), Z, N, pR2, nR2, chR2
        real(dp), allocatable :: r_grid(:), pRho_grid(:), nRho_grid(:)

        Z = real(Params%Calc%Z,dp)
        N = real(Params%Calc%A - Params%Calc%Z,dp)
        if (Params%Calc%CMS == 'CMS1+2B' .or. Params%Calc%CMS == 'CMS2B') then
            R_CMS = rho_radial_cms(Params,Rho,C,Orb)
        else
            R_CMS = 0.0_dp
        end if
        write(*,'(/,a)') 'Calculating the charge-radius squared chR2 ...'
        call rho_radial_density_grid(Params,Rho,C,Orb,r_grid,pRho_grid,nRho_grid)
        call rho_radial_radii(r_grid,pRho_grid,nRho_grid,pRho_grid,pR2,nR2,chR2)
        rho_chR2 = pR2 + R_CMS(1) + R_CMS(2) - 0.106_dp * N / Z + 0.8414_dp - 0.143_dp
        write(*,'(a)') TAB // 'The charge-radius squared chR2 computed ...'
    end function rho_chR2

end module nj_obdm
