! Spherical quasiparticle random-phase approximation ... counterpart of dev QRPA/*.jl ... input io/<Path>/<Reference>/bin,
    ! output io/<Path>/QRPA/{bin,densities,transitions,spectra,amplitudes} ...
module nj_qrpa
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make, orbitals_2qp_make
    use nj_o1b
    use nj_obdm, only : rho_chR2, rho_export
    use nj_tr1b, only : qptr1b_t, qpTr1b_initialize
    use nj_qpo2b, only : qpO2b_import
    use nj_hfb, only : HFB_particle_number
    use nj_qex
    use nj_qtda, only : qex_header
    !$ use omp_lib
    implicit none
    private

    public :: QRPA

contains

    subroutine QRPA(Params)
        type(params_t), intent(in) :: Params
        integer :: N_max, N_threads
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        N_max = Params%Calc%Nmax
        write(*,'(a)') 'Starting spherical QRPA calculation with residual NO2B NN+NNN interaction based on available ' // &
                       'reference quasiparticle mean-field calculation'
        write(*,'(/,a)') 'Calculation data:'
        write(*,'(a)') 'A = ' // istr(Params%Calc%A) // ' , Z = ' // istr(Params%Calc%Z) // ' , hw = ' // jstr(Params%Calc%hw) // &
                       ' MeV , N_max = ' // istr(N_max) // ' , J-scheme LHO basis size = ' // istr((N_max+1)*(N_max+2)/2) // &
                       ' , M-scheme LHO basis size = ' // istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        N_threads = 1
        !$ N_threads = omp_get_max_threads()
        if (N_threads > 1) then
            write(*,'(/,a)') '_________________________________________________________'
            write(*,'(a)') 'Multiple active threads detected ...'
            write(*,'(a)') 'Parallelization report:    Number of active threads = ' // istr(N_threads)
            write(*,'(a)') '_________________________________________________________'
        end if
        write(*,'(a)') trim(Params%Calc%Path)
        if (.not. qex_reference_check(Params,Params%Calc%QRPA%Reference,[character(len=5) :: 'C.bin', 'U.bin', 'V.bin'])) return
        if (is_dir(io_path(Params,'/QRPA'))) call remove_dir(io_path(Params,'/QRPA'))
        call make_dir(io_path(Params,'/QRPA/bin'))
        call make_dir(io_path(Params,'/QRPA/densities'))
        call make_dir(io_path(Params,'/QRPA/transitions/E0'))
        call make_dir(io_path(Params,'/QRPA/transitions/E1'))
        call make_dir(io_path(Params,'/QRPA/transitions/E2'))
        call make_dir(io_path(Params,'/QRPA/transitions/E3'))
        call make_dir(io_path(Params,'/QRPA/spectra'))
        call make_dir(io_path(Params,'/QRPA/amplitudes'))

        write(*,'(/,a)') 'Starting QRPA calculation ...'
        call timer_start(t0)
        call QRPA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'Spherical QRPA solver run completed ...'
    end subroutine QRPA

    subroutine QRPA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(orb2b_2qp_t) :: Orb_2qp
        type(o1b_t) :: C, U, V, Rho
        type(qpo1b_t) :: H_N, qpN
        type(qpo2b_t) :: H_NN
        type(orb2b_t) :: Orb_NN
        type(qptr1b_t) :: qpTrOp
        type(rmat_t), allocatable :: A(:,:), B(:,:)
        type(cvec_t), allocatable :: E(:,:)
        type(rvec_t), allocatable :: E_re(:,:)
        type(cmat_t), allocatable :: X(:,:), Y(:,:)
        type(ivec_t), allocatable :: A_phonon(:,:)
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        real(dp) :: chR2, E_corr
        logical :: Stability
        integer :: J, P
        character(len=:), allocatable :: Ref

        Ref = '/' // trim(Params%Calc%QRPA%Reference) // '/bin/'
        Orb = orbitals_make(Params)
        Orb_2qp = orbitals_2qp_make(Params,Orb)
        C = O1b_import(Params,Orb,io_path(Params,Ref // 'C.bin'))
        U = O1b_import(Params,Orb,io_path(Params,Ref // 'U.bin'))
        V = O1b_import(Params,Orb,io_path(Params,Ref // 'V.bin'))
        H_N = qpO1B_import(Params,Orb,io_path(Params,Ref // 'qpH1b.bin'))
        call qpO2b_import(Params,Orb,io_path(Params,Ref // 'qpH2b.bin'),H_NN,Orb_NN)
        qpTrOp = qpTr1b_initialize(Params,Orb,1.0_dp,C,U,V)
        qpN = qpN1B_make(Params,U,V)

        write(*,'(/,a)') 'Allocating the QRPA matrices A & B ...'
        call qex_allocate(Params,Orb,Orb_NN,Orb_2qp,H_N,H_NN,A,B)
        write(*,'(/,a)') 'The QRPA matrices A & B succesfully allocated ...'

        call QRPA_diagonalize(Params,Orb,Orb_2qp,A,B,qpN,qpTrOp,E,X,Y,Stability)
        A_phonon = QRPA_A_phonon_particle_number(Params,Orb,Orb_2qp,U,V,X,Y)
        Rho = QRPA_OBDM(Params,Orb,Orb_2qp,U,V,Y,A_phonon)
        chR2 = rho_chR2(Params,Orb,C,Rho)
        qpTrOp = qpTr1b_initialize(Params,Orb,chR2,C,U,V)

        allocate(E_re(2,0:ubound(E,2)))
        do J = 0, ubound(E,2)
            do P = 1, 2
                E_re(P,J)%v = real(E(P,J)%v)
            end do
        end do
        write(*,'(/,a)') 'Calculating the QRPA 1-phonon reduced transition matrix elements rM^lambda ...'
        rM = qex_rM(Orb_2qp,X,Y,qpTrOp)
        write(*,'(/,a)') 'Calculating the QRPA 1-phonon reduced transition intensities rB^lambda ...'
        rB = qex_rB(Params,Orb_2qp,rM,E_re)
        E_corr = QRPA_energy(Params,Orb_2qp,E,Y,A_phonon)

        call QRPA_export(Params,Orb,Orb_2qp,Stability,E_corr,Rho,C,E,E_re,X,Y,A_phonon,rB)
        write(*,'(/,a,/)') 'QRPA solver executed properly ...'
    end subroutine QRPA_solver

    ! Solution of one (A - B)(A + B) eigenvalue problem ... E = sqrt(E^2), X & Y from R & (A + B) R ...
    subroutine solve_block(Am,Bm,E,X,Y)
        real(dp), intent(in) :: Am(:,:), Bm(:,:)
        complex(dp), intent(out) :: E(:), X(:,:), Y(:,:)
        real(dp), allocatable :: N_p(:,:), N_m(:,:), H(:,:), R(:,:), N_p_R(:,:)
        complex(dp), allocatable :: W(:), Vc(:,:)
        complex(dp) :: sE, iE, sqrtE
        logical :: Is_Real
        integer :: n, nu

        n = size(Am,1)
        N_p = Am + Bm
        N_m = Am - Bm
        H = mm(N_m,N_p)
        allocate(W(n), Vc(n,n))
        call eig_general(H,W,Vc,Is_Real)
        ! The original stores the eigenvectors in a real buffer ... complex ones would raise an InexactError there ...
        if (any(aimag(Vc) /= 0.0_dp)) call fatal('Complex eigenvectors of (A-B)(A+B) in QRPA ... not supported by the original either ...')
        R = real(Vc)
        do nu = 1, n
            E(nu) = sqrt(W(nu))
        end do
        N_p_R = mm(N_p,R)
        do nu = 1, n
            sqrtE = sqrt(E(nu))
            sE = 0.5_dp * sqrtE
            iE = 0.5_dp / sqrtE
            X(:,nu) = sE * R(:,nu) + iE * N_p_R(:,nu)
            Y(:,nu) = sE * R(:,nu) - iE * N_p_R(:,nu)
        end do
    end subroutine solve_block

    subroutine QRPA_diagonalize(Params,Orb,Orb_2qp,Amat,Bmat,qpN,qpTrOp,E,X,Y,Stability)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(rmat_t), intent(in) :: Amat(:,0:), Bmat(:,0:)
        type(qpo1b_t), intent(in) :: qpN
        type(qptr1b_t), intent(in) :: qpTrOp
        type(cvec_t), allocatable, intent(out) :: E(:,:)
        type(cmat_t), allocatable, intent(out) :: X(:,:), Y(:,:)
        logical, intent(out) :: Stability
        real(dp), allocatable :: Spur(:), Vb(:,:), Us(:,:), A_ort(:,:), B_ort(:,:)
        complex(dp), allocatable :: E_ort(:), X_ort(:,:), Y_ort(:,:), Usc(:,:)
        complex(dp) :: Norm
        real(dp) :: X_norm, Y_norm, T_norm
        integer :: J_max, J, P, N_qp, Spur_ind, a, nu

        J_max = 2*Params%Calc%Nmax + 1
        Stability = .true.
        allocate(E(2,0:J_max), X(2,0:J_max), Y(2,0:J_max))
        write(*,'(/,a)') 'Diagonalizing the QRPA system ...'
        do J = 0, J_max
            do P = 1, 2
                N_qp = Orb_2qp%N(P,J)
                write(*,'(a)') TAB // TAB // 'Diagonalizing the block:' // TAB // 'J = ' // istr(J) // ', P = ' // merge('+','-',P == 1)
                allocate(E(P,J)%v(N_qp), X(P,J)%v(N_qp,N_qp), Y(P,J)%v(N_qp,N_qp))
                if (N_qp == 0) cycle

                if (.not. Params%Calc%QRPA%Ortho .or. .not. ((J == 0 .and. P == 1) .or. (J == 1 .and. P == 2))) then
                    call solve_block(Amat(P,J)%v,Bmat(P,J)%v,E(P,J)%v,X(P,J)%v,Y(P,J)%v)
                else
                    if (J == 0) then
                        Spur = spurious_PN_initialize(Orb,Orb_2qp,qpN)
                    else
                        Spur = spurious_CM_initialize(Orb_2qp,qpTrOp)
                    end if
                    allocate(Vb(N_qp,N_qp), Us(N_qp,N_qp-1))
                    call spurious_orthogonalize(N_qp,Spur,Vb,Spur_ind)
                    do a = 1, Spur_ind - 1
                        Us(:,a) = Vb(:,a)
                    end do
                    do a = Spur_ind + 1, N_qp
                        Us(:,a-1) = Vb(:,a)
                    end do
                    A_ort = mtm(Us,mm(Amat(P,J)%v,Us))
                    B_ort = mtm(Us,mm(Bmat(P,J)%v,Us))
                    allocate(E_ort(N_qp-1), X_ort(N_qp-1,N_qp-1), Y_ort(N_qp-1,N_qp-1))
                    call solve_block(A_ort,B_ort,E_ort,X_ort,Y_ort)
                    Usc = cmplx(Us, 0.0_dp, kind=dp)
                    E(P,J)%v(1) = (0.0_dp,0.0_dp)
                    E(P,J)%v(2:N_qp) = E_ort
                    X(P,J)%v(:,1) = cmplx(Spur / sqrt(2.0_dp), 0.0_dp, kind=dp)
                    Y(P,J)%v(:,1) = cmplx(-Spur / sqrt(2.0_dp), 0.0_dp, kind=dp)
                    X(P,J)%v(:,2:N_qp) = matmul(Usc,X_ort)
                    Y(P,J)%v(:,2:N_qp) = matmul(Usc,Y_ort)
                    deallocate(Vb, Us, E_ort, X_ort, Y_ort)
                end if

                ! Renormalization of X & Y ...
                do nu = 1, N_qp
                    X_norm = sum(abs(X(P,J)%v(:,nu))**2)
                    Y_norm = sum(abs(Y(P,J)%v(:,nu))**2)
                    if ((X_norm - Y_norm) > 1.0e-8_dp) then
                        Norm = cmplx(1.0_dp / sqrt(abs(X_norm - Y_norm)), 0.0_dp, kind=dp)
                    else if ((X_norm - Y_norm) < -1.0e-8_dp) then
                        Norm = (0.0_dp,-1.0_dp) * cmplx(1.0_dp / sqrt(abs(X_norm - Y_norm)), 0.0_dp, kind=dp)
                    else
                        T_norm = sum(abs(X(P,J)%v(:,nu))**2 + abs(Y(P,J)%v(:,nu))**2)
                        Norm = cmplx(1.0_dp / sqrt(abs(T_norm)), 0.0_dp, kind=dp)
                    end if
                    X(P,J)%v(:,nu) = X(P,J)%v(:,nu) * Norm
                    Y(P,J)%v(:,nu) = Y(P,J)%v(:,nu) * Norm
                end do

                do nu = 1, N_qp
                    if (abs(aimag(E(P,J)%v(nu))) > 1.0e-8_dp) then
                        write(*,'(a)') TAB // TAB // TAB // 'QRPA instability detected in channel:' // TAB // 'J = ' // istr(J) // &
                                       ', P = ' // merge('+','-',P == 1)
                        Stability = .false.
                    end if
                end do
            end do
        end do
        write(*,'(a)') TAB // 'The system of QRPA equations has been succesfully solved ...'
    end subroutine QRPA_diagonalize

    ! Mask of phonons with non-vanishing particle-number (20) overlap ... "A-body" levels ...
    function QRPA_A_phonon_particle_number(Params,Orb,Orb_2qp,U,V,X,Y) result(A_phonon)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(o1b_t), intent(in) :: U, V
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(ivec_t), allocatable :: A_phonon(:,:)
        real(dp), allocatable :: pE(:,:), pO(:,:), nE(:,:), nO(:,:)
        integer :: a_max, J_max, a, b, k, l, J, P, nu, qp, T_ab
        complex(dp) :: Sum, Amp

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        J_max = 2*Params%Calc%Nmax + 1
        allocate(pE(a_max,a_max), pO(a_max,a_max), nE(a_max,a_max), nO(a_max,a_max))
        pE = 0.0_dp
        pO = 0.0_dp
        nE = 0.0_dp
        nO = 0.0_dp
        do a = 1, a_max
            do b = 1, a_max
                do k = 1, a_max
                    if (.not. (Orb(a)%l == Orb(k)%l .and. Orb(a)%j == Orb(k)%j)) cycle
                    do l = 1, a_max
                        if (.not. (Orb(l)%l == Orb(b)%l .and. Orb(l)%j == Orb(b)%j)) cycle
                        pE(a,b) = pE(a,b) + (V%p(k,a) * U%p(l,b) + U%p(k,a) * V%p(l,b))
                        nE(a,b) = nE(a,b) + (V%n(k,a) * U%n(l,b) + U%n(k,a) * V%n(l,b))
                        pO(a,b) = pO(a,b) + (V%p(k,a) * U%p(l,b) - U%p(k,a) * V%p(l,b))
                        nO(a,b) = nO(a,b) + (V%n(k,a) * U%n(l,b) - U%n(k,a) * V%n(l,b))
                    end do
                end do
            end do
        end do

        allocate(A_phonon(2,0:J_max))
        do J = 0, J_max
            do P = 1, 2
                allocate(A_phonon(P,J)%v(Orb_2qp%N(P,J)))
                A_phonon(P,J)%v = 0
                do nu = 1, Orb_2qp%N(P,J)
                    Sum = (0.0_dp,0.0_dp)
                    do qp = 1, Orb_2qp%N(P,J)
                        a = Orb_2qp%Ind(P,J)%v(qp)%a
                        b = Orb_2qp%Ind(P,J)%v(qp)%b
                        T_ab = Orb_2qp%Ind(P,J)%v(qp)%T
                        Amp = cmplx(1.0_dp / sqrt(real(1 + kronecker_delta(a,b),dp)), 0.0_dp, kind=dp) * &
                              (X(P,J)%v(qp,nu) - Y(P,J)%v(qp,nu))
                        if (T_ab == -1) then
                            if (mod(J,2) == 0) then
                                Sum = Sum + Amp * pE(a,b)
                            else
                                Sum = Sum + Amp * pO(a,b)
                            end if
                        else
                            if (mod(J,2) == 0) then
                                Sum = Sum + Amp * nE(a,b)
                            else
                                Sum = Sum + Amp * nO(a,b)
                            end if
                        end if
                    end do
                    if (abs(Sum) > 1.0e-8_dp) A_phonon(P,J)%v(nu) = 1
                end do
            end do
        end do
    end function QRPA_A_phonon_particle_number

    function QRPA_energy(Params,Orb_2qp,E,Y,A_phonon) result(E_corr)
        type(params_t), intent(in) :: Params
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(cvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: Y(:,0:)
        type(ivec_t), intent(in) :: A_phonon(:,0:)
        real(dp) :: E_corr, Sum, Nf
        integer :: J, P, qp, nu, a, b

        write(*,'(/,a)') 'Calculating the total QRPA correlation energy ...'
        E_corr = 0.0_dp
        do J = 0, 2*Params%Calc%Nmax + 1
            do P = 1, 2
                Sum = 0.0_dp
                do qp = 1, Orb_2qp%N(P,J)
                    a = Orb_2qp%Ind(P,J)%v(qp)%a
                    b = Orb_2qp%Ind(P,J)%v(qp)%b
                    Nf = -1.0_dp / real(1 + kronecker_delta(a,b),dp)
                    do nu = 1, Orb_2qp%N(P,J)
                        if (A_phonon(P,J)%v(nu) == 0) cycle
                        Sum = Sum + Nf * real(E(P,J)%v(nu)) * abs(Y(P,J)%v(qp,nu))**2
                    end do
                end do
                E_corr = E_corr + real(2*J + 1,dp) * Sum
            end do
        end do
        write(*,'(/,a)') TAB // 'Total QRPA correlation energy ...  E_QRPA = ' // cfmt_f(E_corr,12,6) // ' MeV'
        write(*,'(a)') TAB // 'Calculation of the total QRPA correlation energy completed ...'
    end function QRPA_energy

    ! Diagonal OBDM corrected by the QRPA ground-state quasiparticle occupations ...
    function QRPA_OBDM(Params,Orb,Orb_2qp,U,V,Y,A_phonon) result(Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(o1b_t), intent(in) :: U, V
        type(cmat_t), intent(in) :: Y(:,0:)
        type(ivec_t), intent(in) :: A_phonon(:,0:)
        type(o1b_t) :: Rho
        real(dp), allocatable :: pN(:), nN(:)
        real(dp) :: J_hat, Amp_a, Amp_b, aSum, bSum, pS, nS
        integer :: a_max, J, P, qp, nu, a, b, mu

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(pN(a_max), nN(a_max))
        pN = 0.0_dp
        nN = 0.0_dp
        do J = 0, 2*Params%Calc%Nmax + 1
            do P = 1, 2
                J_hat = real(2*J + 1,dp)
                do qp = 1, Orb_2qp%N(P,J)
                    a = Orb_2qp%Ind(P,J)%v(qp)%a
                    b = Orb_2qp%Ind(P,J)%v(qp)%b
                    Amp_a = J_hat / real(Orb(a)%j + 1,dp) / sqrt(real(1 + kronecker_delta(a,b),dp))
                    Amp_b = J_hat / real(Orb(b)%j + 1,dp) / sqrt(real(1 + kronecker_delta(a,b),dp))
                    aSum = 0.0_dp
                    bSum = 0.0_dp
                    do nu = 1, Orb_2qp%N(P,J)
                        if (A_phonon(P,J)%v(nu) == 0) cycle
                        if (b <= a) aSum = aSum + Amp_a * abs(Y(P,J)%v(qp,nu))**2
                        if (a <= b) bSum = bSum + Amp_b * abs(Y(P,J)%v(qp,nu))**2
                    end do
                    if (Orb_2qp%Ind(P,J)%v(qp)%T == -1) then
                        pN(a) = pN(a) + aSum
                        pN(b) = pN(b) + bSum
                    else
                        nN(a) = nN(a) + aSum
                        nN(b) = nN(b) + bSum
                    end if
                end do
            end do
        end do

        call o1b_alloc(Rho,a_max)
        do a = 1, a_max
            pS = 0.0_dp
            nS = 0.0_dp
            do mu = 1, a_max
                pS = pS + V%p(a,mu) * V%p(a,mu) * (1.0_dp - pN(mu)) + U%p(a,mu) * U%p(a,mu) * pN(mu)
                nS = nS + V%n(a,mu) * V%n(a,mu) * (1.0_dp - nN(mu)) + U%n(a,mu) * U%n(a,mu) * nN(mu)
            end do
            Rho%p(a,a) = pS
            Rho%n(a,a) = nS
        end do
    end function QRPA_OBDM

    subroutine QRPA_export(Params,Orb,Orb_2qp,Stability,E_corr,Rho,C,E,E_re,X,Y,A_phonon,rB)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        logical, intent(in) :: Stability
        real(dp), intent(in) :: E_corr
        type(o1b_t), intent(in) :: Rho, C
        type(cvec_t), intent(in) :: E(:,0:)
        type(rvec_t), intent(in) :: E_re(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(ivec_t), intent(in) :: A_phonon(:,0:)
        type(reduced_transition_t), intent(in) :: rB
        real(dp), allocatable :: E_ord(:), x_ord(:), y_ord(:)
        integer, allocatable :: J_ord(:), P_ord(:), Sort(:)
        integer :: J_max, u, J, P, nu, N_tot, k
        real(dp) :: Z, N

        J_max = 2*Params%Calc%Nmax + 1

        ! Summary ...
        write(*,'(/,a)') 'Preparing QRPA summary file ...'
        open(newunit=u, file=io_path(Params,'/QRPA/QRPA_summary.dat'), status='replace', action='write')
        write(u,'(a)') 'Spherical Quasiparticle Random-Phase Approximation review:'
        call qex_header(u,Params,Params%Calc%QRPA%Ortho)
        write(u,'(/,a)') 'Stability of the QRPA system:'
        if (Stability) then
            write(u,'(a)') TAB // 'QRPA system is STABLE ... no complex eigenvalues detected ...'
        else
            write(u,'(a)') TAB // 'QRPA system is UNSTABLE ... complex eigenvalues detected ...'
        end if
        call qex_dimensions_summary(u,Orb_2qp,J_max,'QRPA 1-phonon 2qp subspaces ...')
        write(u,'(/,a)') 'E_corr     = ' // cfmt_f(E_corr,12,6) // ' ' // TAB // ' MeV ... QRPA ground-state correlation energy'
        write(u,'(a)') 'E_corr / A = ' // cfmt_f(E_corr / real(Params%Calc%A,dp),12,6) // ' ' // TAB // &
                       ' MeV ... QRPA ground-state correlation energy per nuclon'
        close(u)
        write(*,'(/,a)') 'QRPA summary exported ...'

        ! Binary export ... complex X, Y & E ...
        write(*,'(/,a)') 'Performing binary export of QRPA energies E & amplitudes X & Y ...'
        open(newunit=u, file=io_path(Params,'/QRPA/bin/QRPA_X.bin'), access='stream', form='unformatted', status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(X(P,J)%v) > 0) write(u) X(P,J)%v
            end do
        end do
        close(u)
        open(newunit=u, file=io_path(Params,'/QRPA/bin/QRPA_Y.bin'), access='stream', form='unformatted', status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(Y(P,J)%v) > 0) write(u) Y(P,J)%v
            end do
        end do
        close(u)
        open(newunit=u, file=io_path(Params,'/QRPA/bin/QRPA_E.bin'), access='stream', form='unformatted', status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(E(P,J)%v) > 0) write(u) E(P,J)%v
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Binary export of QRPA energies E amplitudes X & Y completed ...'

        ! Spectrum ...
        write(*,'(/,a)') 'Preparing export of QRPA spectra ...'
        open(newunit=u, file=io_path(Params,'/QRPA/spectra/QRPA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('Re E',20,.true.) // ' ' // &
                       cfmt_s('Im E',20,.true.) // ' ' // cfmt_s('A-body',12,.true.)
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, Orb_2qp%N(P,J)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(merge('+','-',P == 1),5,.true.) // ' ' // &
                                   cfmt_f(real(E(P,J)%v(nu)),20,6,.true.) // ' ' // cfmt_f(aimag(E(P,J)%v(nu)),20,6,.true.) // ' ' // &
                                   cfmt_s(trim(merge('true ','false',A_phonon(P,J)%v(nu) == 1)),12,.true.)
                end do
                write(u,'(/)')
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of QRPA spectra succesfully finished ...'

        call qex_spectrum_plot_export(Params,Orb_2qp,E_re,'QRPA')

        ! Norms of the amplitudes ... per channel & energy ordered ...
        write(*,'(/,a)') 'Preparing export of the norms of the QRPA amplitudes X & Y ...'
        open(newunit=u, file=io_path(Params,'/QRPA/amplitudes/QRPA_amplitudes_ordered.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('E',12,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('|Y|^2',12,.true.)
        do J = 0, J_max
            do P = 1, 2
                write(u,'(/,a)') 'J = ' // istr(J) // ',' // TAB // 'P = ' // merge('+','-',P == 1)
                do nu = 1, Orb_2qp%N(P,J)
                    write(u,'(a)') cfmt_f(real(E(P,J)%v(nu)),12,8,.true.) // ' ' // &
                                   cfmt_f(sqrt(sum(abs(X(P,J)%v(:,nu))**2))**2,12,8,.true.) // ' ' // &
                                   cfmt_f(sqrt(sum(abs(Y(P,J)%v(:,nu))**2))**2,12,8,.true.)
                end do
            end do
        end do
        close(u)

        N_tot = sum(Orb_2qp%N)
        allocate(E_ord(N_tot), x_ord(N_tot), y_ord(N_tot), J_ord(N_tot), P_ord(N_tot), Sort(N_tot))
        k = 0
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, Orb_2qp%N(P,J)
                    k = k + 1
                    E_ord(k) = real(E(P,J)%v(nu))
                    x_ord(k) = sum(abs(X(P,J)%v(:,nu))**2)
                    y_ord(k) = sum(abs(Y(P,J)%v(:,nu))**2)
                    J_ord(k) = J
                    P_ord(k) = P
                end do
            end do
        end do
        call sortperm(E_ord,Sort)
        open(newunit=u, file=io_path(Params,'/QRPA/amplitudes/QRPA_amplitudes.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('Ind',12,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // cfmt_s('J',12,.true.) // ' ' // &
                       cfmt_s('P',12,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('|Y|^2',12,.true.)
        do nu = 1, N_tot
            k = Sort(nu)
            write(u,'(a)') cfmt_d(nu,12,.true.) // ' ' // cfmt_f(E_ord(k),12,6,.true.) // ' ' // cfmt_d(J_ord(k),12,.true.) // ' ' // &
                           cfmt_s(merge('+','-',P_ord(k) == 1),12,.true.) // ' ' // cfmt_f(x_ord(k),12,6,.true.) // ' ' // &
                           cfmt_f(y_ord(k),12,6,.true.)
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of QRPA norms of X & Y completed ...'

        ! Particle numbers ...
        write(*,'(/,a)') 'Calculating the QRPA corrections to particle numbers ...'
        call HFB_particle_number(Params,Rho,Orb,Z,N)
        write(*,'(a)') TAB // 'Values of QRPA corrected particle numbers are ...'
        write(*,'(a)') TAB // TAB // 'dZ = ' // cfmt_f(Z,12,6)
        write(*,'(a)') TAB // TAB // 'dN = ' // cfmt_f(N,12,6)
        write(*,'(/,a)') TAB // 'Exporting the corrected QRPA particle numbers ...'
        open(newunit=u, file=io_path(Params,'/QRPA/QRPA_summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'QRPA particle number values Z & N:'
        write(u,'(/,a)') 'Z = ' // cfmt_f(Z,12,6)
        write(u,'(a)') 'N = ' // cfmt_f(N,12,6)
        close(u)

        call rho_export(Params,Orb,io_path(Params,'/QRPA/QRPA_summary.dat'), &
                        io_path(Params,'/QRPA/densities/QRPA_radial_densities.dat'),Rho,C)
        call qex_transitions_export(Params,Orb,Orb_2qp,Rho,E_re,rB,'QRPA',Params%Calc%QRPA%Reference)
    end subroutine QRPA_export

end module nj_qrpa
