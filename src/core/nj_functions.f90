! General functions ... counterpart of Functions/Functions.jl ...
module nj_functions
    use nj_kinds
    use nj_types
    use nj_utils, only : phase, kronecker_delta, fatal
    use nj_wigner, only : fCG, f6j, f9j
    implicit none
    private

    public :: generalized_laguerre, doublefactorial, factorial_i8, Psi_rad_LHO, integrate_trap
    public :: radial_moment_matrix_LHO, radial_moment_LHO, differentiate_Oh6
    public :: associated_legendre_polynomial, normalize_spherical_angles
    public :: scalar_spherical_harmonic_function, vector_spherical_harmonic_function
    public :: J1b_convective, J1b_spin

contains

    pure real(dp) function generalized_laguerre(n,k,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: k, x
        real(dp) :: GL, GL1, GL2
        integer :: j
        if (n == 0) then
            generalized_laguerre = 1.0_dp
        else if (n == 1) then
            generalized_laguerre = 1.0_dp + k - x
        else
            GL1 = 1.0_dp
            GL2 = 1.0_dp + k - x
            GL = GL2
            do j = 2, n
                GL = ((k + real(2*j - 1,dp) - x) * GL2 + (-k - real(j,dp) + 1.0_dp) * GL1) / real(j,dp)
                GL1 = GL2
                GL2 = GL
            end do
            generalized_laguerre = GL
        end if
    end function generalized_laguerre

    pure integer(i8) function doublefactorial(Number)
        integer, intent(in) :: Number
        integer :: m
        doublefactorial = 1_i8
        do m = merge(2, 1, modulo(Number,2) == 0), Number, 2
            doublefactorial = doublefactorial * int(m,i8)
        end do
    end function doublefactorial

    pure integer(i8) function factorial_i8(n)
        integer, intent(in) :: n
        integer :: m
        factorial_i8 = 1_i8
        do m = 2, n
            factorial_i8 = factorial_i8 * int(m,i8)
        end do
    end function factorial_i8

    ! Radial LHO wave function ...
    pure real(dp) function Psi_rad_LHO(r,n,l,nu)
        real(dp), intent(in) :: r, nu
        integer, intent(in) :: n, l
        real(dp) :: N_norm
        N_norm = sqrt(sqrt(2.0_dp * nu**3 / pi) * 2.0_dp**(n + 2*l + 3) * real(factorial_i8(n),dp) * nu**l / &
                      real(doublefactorial(2*n + 2*l + 1),dp))
        Psi_rad_LHO = N_norm * exp(-nu * r**2) * r**l * generalized_laguerre(n, real(l,dp) + 0.5_dp, 2.0_dp * nu * r**2)
    end function Psi_rad_LHO

    pure real(dp) function integrate_trap(x,y)
        real(dp), intent(in) :: x(:), y(:)
        integer :: i
        real(dp) :: dx
        integrate_trap = 0.0_dp
        do i = 1, size(x) - 1
            dx = x(i+1) - x(i)
            integrate_trap = integrate_trap + 0.5_dp * (y(i) + y(i+1)) * dx
        end do
    end function integrate_trap

    function radial_moment_matrix_LHO(Lambda,hw,Orb) result(Radial_lambda)
        integer, intent(in) :: Lambda
        real(dp), intent(in) :: hw
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), allocatable :: Radial_lambda(:,:)
        real(dp) :: b_osc
        integer :: a, b, a_max
        b_osc = 1.0_dp / (HbarC / sqrt(0.5_dp * (m_p + m_n) * hw))
        a_max = size(Orb)
        allocate(Radial_lambda(a_max,a_max))
        !$omp parallel do private(b) schedule(dynamic)
        do a = 1, a_max
            do b = 1, a_max
                Radial_lambda(a,b) = radial_moment_LHO(Lambda,Orb(a)%n,Orb(a)%l,Orb(b)%n,Orb(b)%l,b_osc)
            end do
        end do
        !$omp end parallel do
    end function radial_moment_matrix_LHO

    pure real(dp) function radial_moment_LHO(k,n_a,l_a,n_b,l_b,b_osc)
        integer, intent(in) :: k, n_a, l_a, n_b, l_b
        real(dp), intent(in) :: b_osc
        real(dp) :: gamma(171), dfactorial(171), gi(171), Amp, Sum
        integer :: n, m_a, m_b

        gamma(1) = 1.0_dp
        gamma(2) = sqrt(pi)
        gamma(3) = 1.0_dp
        do n = 3, 170
            gamma(n+1) = (real(n,dp)/2.0_dp - 1.0_dp) * gamma(n-1)
        end do

        dfactorial(1) = 1.0_dp
        dfactorial(2) = 1.0_dp
        do n = 2, 170
            dfactorial(n+1) = real(n,dp) * dfactorial(n-1)
        end do

        gi(1) = sqrt(pi) / 2.0_dp
        gi(2) = 0.5_dp
        do n = 2, 170
            if (modulo(n,2) == 0) then
                gi(n+1) = sqrt(pi) * dfactorial(n) / 2.0_dp**(1 + n/2)
            else
                gi(n+1) = dfactorial(n) / 2.0_dp**(1 + (n-1)/2)
            end if
        end do

        Amp = 4.0_dp / (b_osc**(k+0) * sqrt(pi)) * sqrt(real(factorial_i8(n_a) * factorial_i8(n_b) * &
              2_i8**(n_a + n_b + l_a + l_b),dp) / dfactorial(2*n_a + 2*l_a + 2) / dfactorial(2*n_b + 2*l_b + 2))

        Sum = 0.0_dp
        do m_a = 0, n_a
            do m_b = 0, n_b
                Sum = Sum + phase(m_a + m_b) * gamma(2*n_a + 2*l_a + 3 + 1) * gamma(2*n_b + 2*l_b + 3 + 1) * &
                      gi(k + 0 + 2 + l_a + l_b + 2*m_a + 2*m_b + 1) / real(factorial_i8(m_a),dp) / &
                      real(factorial_i8(m_b),dp) / real(factorial_i8(n_a - m_a),dp) / &
                      real(factorial_i8(n_b - m_b),dp) / gamma(2*l_a + 2*m_a + 3 + 1) / gamma(2*l_b + 2*m_b + 3 + 1)
            end do
        end do

        radial_moment_LHO = Sum * Amp
    end function radial_moment_LHO

    ! Central, forward & backward 6th-order differences with 7-point stencils ...
    function differentiate_Oh6(x,f) result(df)
        real(dp), intent(in) :: x(:), f(:)
        real(dp), allocatable :: df(:)
        integer :: n, j
        real(dp) :: dx
        n = size(x)
        if (n < 7) call fatal('The numerical 6th order O(h^6) numerical differentiation requires at least 7 grid points ...')
        if (size(f) /= n) call fatal('Grid dimensions do not matchh ... numerical differentation is not possible ...')
        dx = x(2) - x(1)
        allocate(df(n))
        do j = 4, n - 3
            df(j) = (f(j+3) - 9.0_dp*f(j+2) + 45.0_dp*f(j+1) - 45.0_dp*f(j-1) + 9.0_dp*f(j-2) - f(j-3)) / (60.0_dp*dx)
        end do
        df(1) = (-147.0_dp*f(1) + 360.0_dp*f(2) - 450.0_dp*f(3) + 400.0_dp*f(4) - 225.0_dp*f(5) + 72.0_dp*f(6) &
                 - 10.0_dp*f(7)) / (60.0_dp*dx)
        df(2) = (-147.0_dp*f(2) + 360.0_dp*f(3) - 450.0_dp*f(4) + 400.0_dp*f(5) - 225.0_dp*f(6) + 72.0_dp*f(7) &
                 - 10.0_dp*f(8)) / (60.0_dp*dx)
        df(3) = (-147.0_dp*f(3) + 360.0_dp*f(4) - 450.0_dp*f(5) + 400.0_dp*f(6) - 225.0_dp*f(7) + 72.0_dp*f(8) &
                 - 10.0_dp*f(9)) / (60.0_dp*dx)
        df(n) = (147.0_dp*f(n) - 360.0_dp*f(n-1) + 450.0_dp*f(n-2) - 400.0_dp*f(n-3) + 225.0_dp*f(n-4) &
                 - 72.0_dp*f(n-5) + 10.0_dp*f(n-6)) / (60.0_dp*dx)
        df(n-1) = (147.0_dp*f(n-1) - 360.0_dp*f(n-2) + 450.0_dp*f(n-3) - 400.0_dp*f(n-4) + 225.0_dp*f(n-5) &
                   - 72.0_dp*f(n-6) + 10.0_dp*f(n-7)) / (60.0_dp*dx)
        df(n-2) = (147.0_dp*f(n-2) - 360.0_dp*f(n-3) + 450.0_dp*f(n-4) - 400.0_dp*f(n-5) + 225.0_dp*f(n-6) &
                   - 72.0_dp*f(n-7) + 10.0_dp*f(n-8)) / (60.0_dp*dx)
    end function differentiate_Oh6

    pure real(dp) function associated_legendre_polynomial(n,m,x)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x
        real(dp) :: Lmm, Lmmp1, Ln, somx2, f
        integer :: i
        Lmm = 1.0_dp
        if (m > 0) then
            somx2 = sqrt(max(0.0_dp, 1.0_dp - x**2))
            f = 1.0_dp
            do i = 1, m
                Lmm = Lmm * (-f * somx2)
                f = f + 2.0_dp
            end do
        end if
        if (n == m) then
            associated_legendre_polynomial = Lmm
            return
        end if
        Lmmp1 = x * real(2*m + 1,dp) * Lmm
        if (n == m + 1) then
            associated_legendre_polynomial = Lmmp1
            return
        end if
        Ln = 0.0_dp
        do i = m + 2, n
            Ln = (real(2*i - 1,dp) * x * Lmmp1 - real(i + m - 1,dp) * Lmm) / real(i - m,dp)
            Lmm = Lmmp1
            Lmmp1 = Ln
        end do
        associated_legendre_polynomial = Ln
    end function associated_legendre_polynomial

    pure subroutine normalize_spherical_angles(theta,phi,theta_norm,phi_norm)
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: theta_norm, phi_norm
        theta_norm = modulo(theta, 2.0_dp * pi)
        phi_norm = modulo(phi, 2.0_dp * pi)
        if (theta_norm > pi) then
            theta_norm = 2.0_dp * pi - theta_norm
            phi_norm = modulo(phi_norm + pi, 2.0_dp * pi)
        end if
    end subroutine normalize_spherical_angles

    pure complex(dp) function scalar_spherical_harmonic_function(l,m,theta_in,phi_in) result(Ylm)
        integer, intent(in) :: l, m
        real(dp), intent(in) :: theta_in, phi_in
        real(dp) :: theta, phi, x, Llm, r, N_norm
        complex(dp) :: Ph, Ylpm
        integer :: m_abs, i
        m_abs = abs(m)
        if (m_abs > l) then
            Ylm = (0.0_dp, 0.0_dp)
            return
        end if
        call normalize_spherical_angles(theta_in,phi_in,theta,phi)
        x = min(max(cos(theta), -1.0_dp), 1.0_dp)
        Llm = associated_legendre_polynomial(l,m_abs,x)
        r = 1.0_dp
        do i = l - m_abs + 1, l + m_abs
            r = r / real(i,dp)
        end do
        N_norm = sqrt(real(2*l + 1,dp) * r / (4.0_dp * pi))
        if (m >= 0) then
            Ph = cmplx(cos(real(m,dp) * phi), sin(real(m,dp) * phi), kind=dp)
            Ylm = Ph * N_norm * Llm
        else
            Ph = cmplx(cos(real(m_abs,dp) * phi), sin(real(m_abs,dp) * phi), kind=dp)
            Ylpm = Ph * N_norm * Llm
            if (modulo(m_abs,2) == 0) then
                Ylm = conjg(Ylpm)
            else
                Ylm = -conjg(Ylpm)
            end if
        end if
    end function scalar_spherical_harmonic_function

    ! Cartesian components of the vector spherical harmonic Y_JLM ...
    pure subroutine vector_spherical_harmonic_function(J,L,M,theta_in,phi_in,Yx,Yy,Yz)
        integer, intent(in) :: J, L, M
        real(dp), intent(in) :: theta_in, phi_in
        complex(dp), intent(out) :: Yx, Yy, Yz
        real(dp) :: theta, phi, isqrt2, CG
        complex(dp) :: CG_Ylm
        complex(dp), parameter :: im = (0.0_dp, 1.0_dp)
        integer :: mu, mm
        call normalize_spherical_angles(theta_in,phi_in,theta,phi)
        Yx = (0.0_dp, 0.0_dp)
        Yy = (0.0_dp, 0.0_dp)
        Yz = (0.0_dp, 0.0_dp)
        isqrt2 = 1.0_dp / sqrt(2.0_dp)
        do mu = -1, 1
            mm = M - mu
            if (abs(mm) <= L) then
                CG = fCG(2*L,2,2*J,2*mm,2*mu,2*M)
                if (abs(CG) > 1.0e-12_dp) then
                    CG_Ylm = CG * scalar_spherical_harmonic_function(L,mm,theta,phi)
                    if (mu == 1) then
                        Yx = Yx - CG_Ylm * isqrt2
                        Yy = Yy - CG_Ylm * isqrt2 * im
                    else if (mu == 0) then
                        Yz = Yz + CG_Ylm
                    else
                        Yx = Yx + CG_Ylm * isqrt2
                        Yy = Yy - CG_Ylm * isqrt2 * im
                    end if
                end if
            end if
        end do
    end subroutine vector_spherical_harmonic_function

    ! 1-body convective current ... returns complex Cartesian components ...
    subroutine J1b_convective(r,theta,phi,l_a,j_a,m_a,R_a,dR_a,l_b,j_b,m_b,R_b,dR_b,J_x,J_y,J_z)
        real(dp), intent(in) :: r, theta, phi, R_a, dR_a, R_b, dR_b
        integer, intent(in) :: l_a, j_a, m_a, l_b, j_b, m_b
        complex(dp), intent(out) :: J_x, J_y, J_z
        real(dp) :: Amp_1, Amp_2
        J_x = (0.0_dp, 0.0_dp)
        J_y = (0.0_dp, 0.0_dp)
        J_z = (0.0_dp, 0.0_dp)

        ! Bra term ...
        Amp_1 = fCG(2*l_a,1,j_a,m_a-1,1,m_a) * fCG(2*l_b,1,j_b,m_b-1,1,m_b)
        if (abs(Amp_1) > 1.0e-12_dp) call add_terms(Amp_1, julia_div(m_a-1,2), julia_div(m_b-1,2))

        ! Ket term ...
        Amp_2 = fCG(2*l_a,1,j_a,m_a+1,-1,m_a) * fCG(2*l_b,1,j_b,m_b+1,-1,m_b)
        if (abs(Amp_2) > 1.0e-12_dp) call add_terms(Amp_2, julia_div(m_a+1,2), julia_div(m_b+1,2))

    contains

        subroutine add_terms(Amp,ma2,mb2)
            real(dp), intent(in) :: Amp
            integer, intent(in) :: ma2, mb2
            complex(dp) :: Y, Y_x, Y_y, Y_z, Jc
            ! First bracket ...
            if (l_a > 0) then
                Y = scalar_spherical_harmonic_function(l_b,mb2,theta,phi)
                call vector_spherical_harmonic_function(l_a,l_a-1,-ma2,theta,phi,Y_x,Y_y,Y_z)
                Jc = sqrt(real(l_a,dp) / real(2*l_a + 1,dp)) * (dR_a + real(l_a + 1,dp) / r * R_a) * R_b * Y * Amp * phase(ma2)
                J_x = J_x + Jc * Y_x
                J_y = J_y + Jc * Y_y
                J_z = J_z + Jc * Y_z
                call vector_spherical_harmonic_function(l_a,l_a+1,-ma2,theta,phi,Y_x,Y_y,Y_z)
                Jc = -sqrt(real(l_a - 1,dp) / real(2*l_a + 1,dp)) * (dR_a - real(l_a,dp) / r * R_a) * R_b * Y * Amp * phase(ma2)
                J_x = J_x + Jc * Y_x
                J_y = J_y + Jc * Y_y
                J_z = J_z + Jc * Y_z
            end if
            ! Second bracket ...
            if (l_b > 0) then
                Y = phase(ma2) * scalar_spherical_harmonic_function(l_a,-ma2,theta,phi)
                call vector_spherical_harmonic_function(l_b,l_b-1,mb2,theta,phi,Y_x,Y_y,Y_z)
                Jc = -sqrt(real(l_b,dp) / real(2*l_b + 1,dp)) * (dR_b + real(l_b + 1,dp) / r * R_b) * R_a * Y * Amp
                J_x = J_x + Jc * Y_x
                J_y = J_y + Jc * Y_y
                J_z = J_z + Jc * Y_z
                call vector_spherical_harmonic_function(l_b,l_b+1,mb2,theta,phi,Y_x,Y_y,Y_z)
                Jc = sqrt(real(l_b - 1,dp) / real(2*l_b + 1,dp)) * (dR_b - real(l_b,dp) / r * R_b) * R_a * Y * Amp
                J_x = J_x + Jc * Y_x
                J_y = J_y + Jc * Y_y
                J_z = J_z + Jc * Y_z
            end if
        end subroutine add_terms

    end subroutine J1b_convective

    ! 1-body spin current ...
    subroutine J1b_spin(r,theta,phi,l_a,j_a,m_a,R_a,dR_a,l_b,j_b,m_b,R_b,dR_b,J_x,J_y,J_z)
        real(dp), intent(in) :: r, theta, phi, R_a, dR_a, R_b, dR_b
        integer, intent(in) :: l_a, j_a, m_a, l_b, j_b, m_b
        complex(dp), intent(out) :: J_x, J_y, J_z
        complex(dp) :: Y_x, Y_y, Y_z
        real(dp) :: Amp, C_la_lb_L, Amp_JL, ME
        integer :: M, L, JJ

        J_x = (0.0_dp, 0.0_dp)
        J_y = (0.0_dp, 0.0_dp)
        J_z = (0.0_dp, 0.0_dp)
        M = julia_div(m_b - m_a,2)
        Amp = phase(l_a + julia_div(j_b + m_a,2)) * sqrt(real((j_a + 1) * (j_b + 1) * (2*l_a + 1) * (2*l_b + 1),dp) * &
              3.0_dp / (8.0_dp * pi))

        do L = abs(l_a - l_b), l_a + l_b
            C_la_lb_L = fCG(2*l_a,2*l_b,2*L,0,0,0)
            do JJ = abs(j_a - j_b)/2, (j_a + j_b)/2
                if (abs(M) > JJ) cycle
                Amp_JL = Amp * C_la_lb_L * phase(JJ) * fCG(j_a,j_b,2*JJ,-m_a,m_b,2*M) * f9j(2*l_a,j_a,1,2*l_b,j_b,1,2*L,2*JJ,2)
                if (abs(Amp_JL) < 1.0e-12_dp) cycle
                if (JJ == (L - 1) .and. JJ > 0) then
                    call vector_spherical_harmonic_function(JJ,JJ,M,theta,phi,Y_x,Y_y,Y_z)
                    ME = Amp_JL * sqrt(real(JJ,dp) / real(2*JJ + 1,dp)) * (dR_a * R_b + R_a * dR_b + real(JJ + 2,dp) / r * R_a * R_b)
                    call acc(ME)
                end if
                if (JJ == L .and. JJ > 0) then
                    call vector_spherical_harmonic_function(JJ,JJ+1,M,theta,phi,Y_x,Y_y,Y_z)
                    ME = Amp_JL * sqrt(real(JJ,dp) / real(2*JJ + 1,dp)) * (dR_a * R_b + R_a * dR_b - real(JJ,dp) / r * R_a * R_b)
                    call acc(ME)
                end if
                if (JJ == L) then
                    call vector_spherical_harmonic_function(JJ,JJ-1,M,theta,phi,Y_x,Y_y,Y_z)
                    ME = Amp_JL * sqrt(real(JJ + 1,dp) / real(2*JJ + 1,dp)) * (dR_a * R_b + R_a * dR_b + &
                         real(JJ + 1,dp) / r * R_a * R_b)
                    call acc(ME)
                end if
                if (JJ == (L + 1)) then
                    call vector_spherical_harmonic_function(JJ,JJ,M,theta,phi,Y_x,Y_y,Y_z)
                    ME = Amp_JL * sqrt(real(JJ + 1,dp) / real(2*JJ + 1,dp)) * (dR_a * R_b + R_a * dR_b - &
                         real(JJ - 1,dp) / r * R_a * R_b)
                    call acc(ME)
                end if
            end do
        end do

    contains

        subroutine acc(ME)
            real(dp), intent(in) :: ME
            J_x = J_x + ME * Y_x
            J_y = J_y + ME * Y_y
            J_z = J_z + ME * Y_z
        end subroutine acc

    end subroutine J1b_spin

    ! Julia div(a,b) truncates toward zero ... same as Fortran integer division ...
    pure integer function julia_div(a,b)
        integer, intent(in) :: a, b
        julia_div = a / b
    end function julia_div

end module nj_functions
