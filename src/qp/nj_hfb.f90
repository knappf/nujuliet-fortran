! Spherical HFB solver ... counterpart of dev HFB/*.jl (HFB.jl, HFB_solver.jl, HFB_Broyden.jl, HFB_chemical_potential.jl,
    ! HFB_ordering.jl, HFB_canonical_basis.jl, HFB_export.jl, ...) ... output io/<Path>/HFB/{bin,densities} ...
module nj_hfb
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make
    use nj_o1b
    use nj_o2b
    use nj_v2b, only : V2b_read, V2b_residual_no2b
    use nj_v3b, only : V3b_no2b_read, V3b_no2b
    use nj_obdm, only : rho_export
    use nj_hf, only : HF_orbitals_make
    use nj_hf_driver, only : print_run_header
    use nj_qpo2b
    use nj_qp_common, only : pairing_energy
    implicit none
    private

    public :: HFB, HFB_density_operator, HFB_particle_number

    ! Modified Broyden mixing state ... counterpart of dev HFB_Broyden & HFB_Broyden_Vector (HFB_Broyden.jl) ... the 4
        ! components (last index) are 1 pRho, 2 nRho, 3 pKappa, 4 nKappa; M linearized (l,j)-block elements (a <= b),
        ! m_hist history vectors ...
    type :: broyden_t
        integer :: M = 0, m_hist = 0
        integer, allocatable :: Key(:,:)                    ! (2,M) ... a, b with a <= b
        real(dp), allocatable :: w(:)                       ! (M) weights
        real(dp), allocatable :: x(:,:), y(:,:), z(:,:), r(:,:), s(:,:)       ! (M,4)
        real(dp), allocatable :: Xh(:,:,:), Rh(:,:,:)                       ! (M,m_hist,4)
    end type broyden_t

    interface
        subroutine dsysv(Uplo,n,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
            import :: dp
            character(len=1), intent(in) :: Uplo
            integer, intent(in) :: n, Nrhs, Lda, Ldb, Lwork
            real(dp), intent(inout) :: A(Lda,*), B(Ldb,*), Work(*)
            integer, intent(out) :: Ipiv(*), Info
        end subroutine dsysv
    end interface

contains

    subroutine HFB(Params)
        type(params_t), intent(in) :: Params
        integer(i8) :: t0

        call wigner_init_float(Params%Int%N2max + 2, 'Jmax', 6)
        write(*,'(a)') 'Starting spherical HFB calculation with NO2B NN+NNN interaction ...'
        call print_run_header(Params)

        if (is_dir(io_path(Params,''))) call remove_dir(io_path(Params,''))
        call make_dir(io_path(Params,''))
        call make_dir(io_path(Params,'/HFB'))
        call make_dir(io_path(Params,'/HFB/bin'))
        call make_dir(io_path(Params,'/HFB/densities'))

        write(*,'(/,a)') 'Starting spherical Hartree-Fock-Bogoliubov calculation ...'
        call timer_start(t0)
        call HFB_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'All calculations have finished ...'
    end subroutine HFB

    subroutine HFB_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: T, U, V, C, Rho, Kappa, H, Delta, Rho_c, Kappa_c
        type(o2b_t) :: V_NN, V_NN_res
        type(orb2b_t) :: Orb_NN
        type(v3b_t) :: V_NNN
        type(orb3b_t) :: Orb_NNN
        type(pnvector_t) :: SQE
        type(pnfloat_t) :: Lambda, dA
        type(qpo1b_t) :: H_N
        type(qpo2b_t) :: H_NN
        real(dp) :: E_HFB(4), T_HFB
        logical :: Convergence
        integer :: Iteration
        integer(i8) :: t0

        Orb = orbitals_make(Params)
        T = T1b(Params,Orb)
        call timer_start(t0)
        call V2b_read(Params,V_NN,Orb_NN)
        call timer_report(t0)
        call timer_start(t0)
        call V3b_no2b_read(Params,V_NNN,Orb_NNN)
        call timer_report(t0)

        call timer_start(t0)
        call HFB_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Lambda,SQE,U,V,C,Rho,Kappa,H,Delta,Convergence,Iteration)
        call timer_report(t0)

        E_HFB = HFB_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,.false.)
        T_HFB = T1b_energy(Params,Rho,Orb,T)
        dA = HFB_particle_number_dispersion(Params,Rho,Orb)
        call HFB_summary(Params,E_HFB,T_HFB,Lambda,dA,Convergence,Iteration)

        call rho_export(Params,Orb,io_path(Params,'/HFB/HFB_summary.dat'), &
                        io_path(Params,'/HFB/densities/HFB_radial_densities.dat'),Rho,C)
        Rho_c%p = mtmm(C%p,Rho%p)
        Rho_c%n = mtmm(C%n,Rho%n)
        call HFB_summary_SQS(Params,SQE,Rho_c,Orb)

        if (Params%Calc%HFB%BMF) then
            H_N%qp11%p = diagm(SQE%p)
            H_N%qp11%n = diagm(SQE%n)
            call o1b_alloc(H_N%qp20,size(SQE%p))
            call timer_start(t0)
            call V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)
            call O2b_move(V_NN_res,V_NN)
            call timer_report(t0)
            if (Params%Calc%HFB%LNT .and. Params%Calc%HFB%LNRes) call HFB_Lipkin_Nogami_V2b_residual_no2b(Params,Orb,Orb_NN,V_NN)

            U%p = mtm(C%p,U%p)
            U%n = mtm(C%n,U%n)
            V%p = mtm(C%p,V%p)
            V%n = mtm(C%n,V%n)

            call timer_start(t0)
            call qpO2b(Params,Orb,Orb_NN,V_NN,C,U,V,H_NN)
            call timer_report(t0)

            Kappa_c%p = mtmm(C%p,Kappa%p)
            Kappa_c%n = mtmm(C%n,Kappa%n)
            call O1b_export(Params,Orb,C,io_path(Params,'/HFB/bin/C.bin'))
            call O1b_export(Params,Orb,U,io_path(Params,'/HFB/bin/U.bin'))
            call O1b_export(Params,Orb,V,io_path(Params,'/HFB/bin/V.bin'))
            call O1b_export(Params,Orb,Rho_c,io_path(Params,'/HFB/bin/Rho.bin'))
            call O1b_export(Params,Orb,Kappa_c,io_path(Params,'/HFB/bin/Kappa.bin'))
            call qpO1B_export(Params,Orb,H_N,io_path(Params,'/HFB/bin/qpH1b.bin'))
            call qpO2b_export(Params,Orb_NN,H_NN,io_path(Params,'/HFB/bin/qpH2b.bin'))

            call qp_BMBPT_energy(Params,Orb,Orb_NN,H_N,H_NN,'HFB',io_path(Params,'/HFB/HFB_summary.dat'))
        end if
    end subroutine HFB_solver

    ! Diagonalization of the HFB equations with the pairing treatment selected by Params%Calc%HFB%Pairing ...
    subroutine diagonalize(Params,H,Delta,Orb,SQE,U,V)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: H, Delta
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnvector_t), intent(out) :: SQE
        type(o1b_t), intent(out) :: U, V
        select case (trim(Params%Calc%HFB%Pairing))
        case ('Full')
            call HFB_diagonalize(Params,H,Delta,Orb,SQE,U,V)
        case ('MCA')
            call HFB_diagonalize_MCA(Params,H,Delta,Orb,SQE,U,V)
        case ('BCS')
            call HFB_diagonalize_BCS(Params,H,Delta,Orb,SQE,U,V)
        case default
            call fatal('Unknown HFB pairing option ' // trim(Params%Calc%HFB%Pairing))
        end select
    end subroutine diagonalize

    subroutine HFB_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Lambda,SQE,U,V,C,Rho,Kappa,H,Delta,Convergence,Iteration)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: T
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(pnfloat_t), intent(out) :: Lambda
        type(pnvector_t), intent(out) :: SQE
        type(o1b_t), intent(out) :: U, V, C, Rho, Kappa, H, Delta
        logical, intent(out) :: Convergence
        integer, intent(out) :: Iteration
        type(o1b_t) :: U_1, U_2, V_1, V_2, Rho_1, Rho_2, Kappa_1, Kappa_2, H_1, H_2, Delta_1, Delta_2, Rho_Temp, Kappa_Temp
        type(pnvector_t) :: SQE_1, SQE_2, SQE_old
        type(pnfloat_t) :: Lambda_1, Lambda_2, Lambda_H
        type(broyden_t) :: Broyden
        real(dp) :: Z_target, N_target, epsilon, dE, d2E, dZ, dN, Z, N, Z_1, Z_2, N_1, N_2, E_HFB_1(4), E_HFB_2(4), Eye_shift
        ! Bracket Lambda_lo < Lambda < Lambda_hi (index 1 p, 2 n) & the shell-gap Lambda ... dev NaN = unknown is
            ! represented by the logical flags Has_* ...
        real(dp) :: Lambda_lo(2), Lambda_hi(2), Lambda_gap(2)
        logical :: Has_lo(2), Has_hi(2), Has_gap(2)
        integer :: a_max, Iteration_max, L, a
        logical :: Degeneracy

        Z_target = real(Params%Calc%Z,dp)
        N_target = real(Params%Calc%A - Params%Calc%Z,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        epsilon = Params%Calc%HFB%Tol
        Iteration = 0
        Iteration_max = Params%Calc%HFB%Imax
        Convergence = .false.
        Degeneracy = .false.
        dE = 0.1_dp
        d2E = 0.1_dp
        dZ = 0.1_dp
        dN = 0.1_dp
        E_HFB_1 = 0.0_dp
        E_HFB_2 = 0.0_dp
        allocate(SQE%p(a_max), SQE%n(a_max))
        SQE%p = 0.0_dp
        SQE%n = 0.0_dp
        call o1b_alloc(Rho_1,a_max)
        call o1b_alloc(Kappa_1,a_max)

        Lambda = pnfloat_t(Params%Calc%HFB%pL0, Params%Calc%HFB%nL0)
        call HFB_density_operator_initialize(Params,Orb,Rho,Kappa)
        call HFB_Broyden_initialize(Params,Rho,Kappa,Orb,Broyden)
        SQE_old = SQE

        Z_1 = 0.0_dp
        Z_2 = 0.0_dp
        N_1 = 0.0_dp
        N_2 = 0.0_dp
        Lambda_1 = pnfloat_t(0.75_dp * Lambda%p, 0.75_dp * Lambda%n)
        Lambda_2 = Lambda
        Lambda_H = Lambda               ! chemical potential currently contained in the field H
        Lambda_lo = 0.0_dp
        Lambda_hi = 0.0_dp
        Lambda_gap = 0.0_dp
        Has_lo = .false.
        Has_hi = .false.
        Has_gap = .false.

        write(*,'(/,a)') 'Starting iteration of HFB equations ...'
        do while (((dE > epsilon) .or. (dZ > epsilon) .or. (dN > epsilon)) .and. (Iteration < Iteration_max))
            ! Safeguarded secant iterations for the chemical potentials ...
            do L = 1, Iteration_max
                ! Reset of the bracket ... fields H & Delta have changed ...
                if (L == 1) then
                    Has_lo = .false.
                    Has_hi = .false.
                end if
                if ((L == 1) .and. (Iteration == 0)) then
                    call HFB_allocate(Params,Lambda_1,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,H_1,Delta_1)
                    call HFB_allocate(Params,Lambda_2,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,H_2,Delta_2)
                    call diagonalize(Params,H_1,Delta_1,Orb,SQE_1,U_1,V_1)
                    call diagonalize(Params,H_2,Delta_2,Orb,SQE_2,U_2,V_2)
                    call HFB_density_operator(Params,U_1,V_1,Orb,Rho_1,Kappa_1)
                    call HFB_density_operator(Params,U_2,V_2,Orb,Rho_2,Kappa_2)
                    call HFB_particle_number(Params,Rho_1,Orb,Z_1,N_1)
                    call HFB_particle_number(Params,Rho_2,Orb,Z_2,N_2)
                    Z = Z_2
                    N = N_2
                    Lambda = HFB_chemical_potential_secant(Params,Lambda_1,pnfloat_t(Z_target - Z_1, N_target - N_1), &
                                                           Lambda_2,pnfloat_t(Z_target - Z_2, N_target - N_2), &
                                                           Lambda_lo,Lambda_hi,Has_lo,Has_hi)
                    Lambda_1 = Lambda_2
                    Lambda_2 = Lambda
                    Lambda_H = Lambda_1           ! H = H_2 was built with the old Lambda_2 (now Lambda_1)
                    SQE = SQE_2
                    U = U_2
                    V = V_2
                    H = H_2
                    Delta = Delta_2
                end if

                ! Closed-shell systems (collapsed pairing) ... Lambda in the middle of the shell gap ...
                if (L == 1) call HFB_chemical_potential_gap(Params,H,Lambda_H,Delta,Orb,Lambda_gap,Has_gap)
                if (Has_gap(1)) Lambda%p = Lambda_gap(1)
                if (Has_gap(2)) Lambda%n = Lambda_gap(2)
                Lambda_2 = Lambda

                ! Shift of H for the updated chemical potential ... dev: by Lambda_H - Lambda (base: Lambda_1 - Lambda) ...
                do a = 1, a_max
                    Eye_shift = Lambda_H%p - Lambda%p
                    H%p(a,a) = H%p(a,a) + Eye_shift
                    Eye_shift = Lambda_H%n - Lambda%n
                    H%n(a,a) = H%n(a,a) + Eye_shift
                end do
                Lambda_H = Lambda

                call diagonalize(Params,H,Delta,Orb,SQE,U,V)
                call HFB_density_operator(Params,U,V,Orb,Rho_Temp,Kappa_Temp)
                call HFB_particle_number(Params,Rho_Temp,Orb,Z,N)
                Z_1 = Z_2
                N_1 = N_2
                Z_2 = Z
                N_2 = N
                Lambda = HFB_chemical_potential_secant(Params,Lambda_1,pnfloat_t(Z_target - Z_1, N_target - N_1), &
                                                       Lambda_2,pnfloat_t(Z_target - Z_2, N_target - N_2), &
                                                       Lambda_lo,Lambda_hi,Has_lo,Has_hi)
                Lambda_1 = Lambda_2
                Lambda_2 = Lambda
                dZ = abs(Z_target - Z)
                dN = abs(N_target - N)
                if (dZ < epsilon .and. dN < epsilon) exit
            end do

            ! New densities & their modified Broyden mixing ...
            call HFB_density_operator(Params,U,V,Orb,Rho,Kappa)
            call HFB_Broyden_update(Params,Iteration,Rho,Kappa,Broyden)

            call HFB_allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,H,Delta)
            Lambda_H = Lambda
            if (Params%Calc%HFB%LNT .and. Iteration > 15) call HFB_Lipkin_Nogami(Params,Orb,Orb_NN,Orb_NNN,Rho,Kappa,H,V_NN,V_NNN)
            call diagonalize(Params,H,Delta,Orb,SQE,U,V)
            call HFB_particle_number(Params,Rho,Orb,Z,N)

            d2E = abs(dE - (sum(abs(SQE%p - SQE_old%p)) + sum(abs(SQE%n - SQE_old%n))) / real(2 * a_max,dp))
            dE = (sum(abs(SQE%p - SQE_old%p)) + sum(abs(SQE%n - SQE_old%n))) / real(2 * a_max,dp)
            dZ = abs(Z_target - Z)
            dN = abs(N_target - N)

            if (.not. Degeneracy) then
                Iteration = Iteration + 1
                write(*,'(a)') TAB // 'HFB iteration number: ' // cfmt_d(Iteration,4) // '   dE = ' // cfmt_f(dE,12,9) // &
                               ' MeV,   dZ = ' // cfmt_f(dZ,12,9) // ',   dN = ' // cfmt_f(dN,12,9)
            end if
            SQE_old = SQE

            ! Resolve degenerate solutions ...
            if (Degeneracy) then
                H_2 = H
                Delta_2 = Delta
                Lambda_2 = Lambda
                call diagonalize(Params,H_2,Delta_2,Orb,SQE_2,U_2,V_2)
                call HFB_density_operator(Params,U_2,V_2,Orb,Rho_2,Kappa_2)
                E_HFB_2 = HFB_energy(Params,Rho_2,Kappa_2,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,.true.)
                if (sum(E_HFB_1) < sum(E_HFB_2)) then
                    U = U_1
                    V = V_1
                    Lambda = Lambda_1
                    H = H_1
                    Delta = Delta_1
                    SQE = SQE_1
                    Rho = Rho_1
                    Kappa = Kappa_1
                else
                    U = U_2
                    V = V_2
                    Lambda = Lambda_2
                    H = H_2
                    Delta = Delta_2
                    SQE = SQE_2
                    Rho = Rho_2
                    Kappa = Kappa_2
                end if
                write(*,'(/,a,/)') 'The degenerate solution has been chosen ... HFB iteration terminates'
                Convergence = .true.
                exit
            end if

            if (dE < epsilon .and. dZ < epsilon .and. dN < epsilon) then
                write(*,'(/,a)') 'HFB iteration with residual NN interaction has converged ...'
                Convergence = .true.
                exit
            end if

            ! Detection of a degenerate solution ... Rho_1 & Kappa_1 are those of the secant initialization as in the original ...
            if (d2E < 1.0e-3_dp * dE .and. dZ < 1.0e-3_dp .and. dN < 1.0e-3_dp .and. Iteration > 100) then
                write(*,'(/,a)') 'HFB iteration stucked at a degenerate solution ... Degenerate solutions will be analyzed ...'
                Degeneracy = .true.
                H_1 = H
                Delta_1 = Delta
                Lambda_1 = Lambda
                call diagonalize(Params,H_1,Delta_1,Orb,SQE_1,U_1,V_1)
                E_HFB_1 = HFB_energy(Params,Rho_1,Kappa_1,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,.true.)
            end if

            if (Iteration == Iteration_max) then
                write(*,'(/,a,/)') 'Maximum number of HFB iterations reached ... Iteration terminates & convergence was not reached ...'
                Convergence = .false.
                exit
            end if
        end do

        if (Params%Calc%HFB%Pairing == 'BCS' .or. Params%Calc%HFB%Pairing == 'MCA') then
            C = HFB_canonical_basis_BCS(Params,H,Orb)
        else
            C = HFB_canonical_basis(Params,Rho,Orb)
        end if

        call HFB_orbital_ordering(Params,SQE,U,V,Orb,.true.)
    end subroutine HFB_solve

    ! Builds & diagonalizes the 2a_max x 2a_max HFB matrix of one isospin ... positive energy half ...
    subroutine bogoliubov_eig(Hm,Dm,a_max,Uo,Vo,E)
        real(dp), intent(in) :: Hm(:,:), Dm(:,:)
        integer, intent(in) :: a_max
        real(dp), intent(out) :: Uo(:,:), Vo(:,:), E(:)
        real(dp), allocatable :: M(:,:), W(:,:), Ev(:)
        integer :: a
        allocate(M(2*a_max,2*a_max), W(2*a_max,2*a_max), Ev(2*a_max))
        M(1:a_max,1:a_max) = Hm
        M(1:a_max,a_max+1:2*a_max) = Dm
        M(a_max+1:2*a_max,1:a_max) = Dm
        M(a_max+1:2*a_max,a_max+1:2*a_max) = -Hm
        call eig_sym(M,Ev,W)
        do a = 1, a_max
            E(a) = Ev(a+a_max)
            Uo(:,a) = W(1:a_max,a+a_max)
            Vo(:,a) = W(a_max+1:2*a_max,a+a_max)
        end do
    end subroutine bogoliubov_eig

    subroutine clean(M,Tol)
        real(dp), intent(inout) :: M(:,:)
        real(dp), intent(in) :: Tol
        where (abs(M) < Tol) M = 0.0_dp
    end subroutine clean

    ! Rho = V V^T restricted to the (l,j) blocks ...
    function block_rho(Vm,Orb) result(R)
        real(dp), intent(in) :: Vm(:,:)
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp) :: R(size(Vm,1),size(Vm,1))
        integer :: a, b, c, n
        real(dp) :: Sum
        n = size(Vm,1)
        R = 0.0_dp
        do a = 1, n
            do b = 1, n
                if ((Orb(a)%j == Orb(b)%j) .and. (Orb(a)%l == Orb(b)%l)) then
                    Sum = 0.0_dp
                    do c = 1, n
                        if ((Orb(a)%j == Orb(c)%j) .and. (Orb(a)%l == Orb(c)%l)) Sum = Sum + Vm(a,c) * Vm(b,c)
                    end do
                    R(a,b) = Sum
                end if
            end do
        end do
    end function block_rho

    ! U & V from the BMZ decomposition U = C U_BMZ D, V = C V_BMZ D ... Full = .true. for HFB_diagonalize, .false. for MCA ...
    subroutine bmz_amplitudes(Orb,Cm,Occ,Um,Vm,Full)
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: Cm(:,:), Occ(:)
        real(dp), intent(inout) :: Um(:,:), Vm(:,:)
        logical, intent(in) :: Full
        real(dp), allocatable :: U_BMZ(:,:), V_BMZ(:,:), CU(:,:), CV(:,:), D(:,:)
        real(dp) :: pu, pv, Ph
        integer :: n, a, b
        n = size(Occ)
        allocate(U_BMZ(n,n), V_BMZ(n,n))
        U_BMZ = 0.0_dp
        V_BMZ = 0.0_dp
        do a = 1, n
            Ph = phase(Orb(a)%l)
            if (Full) then
                pv = sqrt(min(1.0_dp,Occ(a)))
                pu = Ph * sqrt(abs(1.0_dp - pv**2))
            else
                pu = Ph * sqrt(max(0.0_dp, 1.0_dp - Occ(a)))
                pv = sqrt(min(1.0_dp,Occ(a)))
            end if
            U_BMZ(a,a) = pu
            V_BMZ(a,a) = pv
        end do
        CU = mtm(Cm,Um)
        CV = mtm(Cm,Vm)
        D = identity(n)
        do a = 1, n
            do b = 1, n
                if (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j) then
                    if (abs(U_BMZ(a,a)) > abs(V_BMZ(a,a))) then
                        D(a,b) = CU(a,b) / U_BMZ(a,a)
                    else
                        D(a,b) = CV(a,b) / V_BMZ(a,a)
                    end if
                end if
            end do
        end do
        Um = mm(mm(Cm,U_BMZ),D)
        Vm = mm(mm(Cm,V_BMZ),D)
    end subroutine bmz_amplitudes

    subroutine clean_occ(Occ,Tol)
        real(dp), intent(inout) :: Occ(:)
        real(dp), intent(in) :: Tol
        integer :: a
        do a = 1, size(Occ)
            if (Occ(a) < Tol) then
                Occ(a) = 0.0_dp
            else if ((1.0_dp - Occ(a)) < Tol) then
                Occ(a) = 1.0_dp
            end if
        end do
    end subroutine clean_occ

    subroutine HFB_diagonalize(Params,H,Delta,Orb,SQE,U,V)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: H, Delta
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnvector_t), intent(out) :: SQE
        type(o1b_t), intent(out) :: U, V
        type(o1b_t) :: Cc
        type(pnvector_t) :: Occ
        real(dp), allocatable :: pRho(:,:), nRho(:,:)
        real(dp), parameter :: Tol = 1.0e-11_dp, nTol = 1.0e-13_dp
        real(dp) :: R
        integer :: N_max, a_max, a, b

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(U%p(a_max,a_max), U%n(a_max,a_max), V%p(a_max,a_max), V%n(a_max,a_max), SQE%p(a_max), SQE%n(a_max))
        call bogoliubov_eig(H%p,Delta%p,a_max,U%p,V%p,SQE%p)
        call bogoliubov_eig(H%n,Delta%n,a_max,U%n,V%n,SQE%n)
        call clean(U%p,nTol)
        call clean(V%p,nTol)
        call clean(U%n,nTol)
        call clean(V%n,nTol)
        call HFB_orbital_ordering(Params,SQE,U,V,Orb,.false.)

        pRho = block_rho(V%p,Orb(1:a_max))
        nRho = block_rho(V%n,Orb(1:a_max))
        pRho = 0.5_dp * (pRho + transpose(pRho))
        nRho = 0.5_dp * (nRho + transpose(nRho))
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%l /= Orb(b)%l .or. Orb(a)%j /= Orb(b)%j) then
                    pRho(a,b) = 0.0_dp
                    nRho(a,b) = 0.0_dp
                end if
            end do
            R = 10.0_dp * real(Orb(a)%l*(N_max+1) + (Orb(a)%j+1)/2,dp)
            pRho(a,a) = pRho(a,a) + R
            nRho(a,a) = nRho(a,a) + R
        end do

        allocate(Occ%p(a_max), Occ%n(a_max), Cc%p(a_max,a_max), Cc%n(a_max,a_max))
        call eig_sym(pRho,Occ%p,Cc%p,Descending=.true.)
        call eig_sym(nRho,Occ%n,Cc%n,Descending=.true.)
        call HFB_canonical_orbital_ordering(Params,Occ,Cc,Orb)

        do a = 1, a_max
            R = 10.0_dp * real(Orb(a)%l*(N_max+1) + (Orb(a)%j+1)/2,dp)
            Occ%p(a) = Occ%p(a) - R
            Occ%n(a) = Occ%n(a) - R
        end do
        call clean_occ(Occ%p,Tol)
        call clean_occ(Occ%n,Tol)

        call bmz_amplitudes(Orb(1:a_max),Cc%p,Occ%p,U%p,V%p,.true.)
        call bmz_amplitudes(Orb(1:a_max),Cc%n,Occ%n,U%n,V%n,.true.)
        call clean(U%p,nTol)
        call clean(V%p,nTol)
        call clean(U%n,nTol)
        call clean(V%n,nTol)
    end subroutine HFB_diagonalize

    subroutine HFB_diagonalize_MCA(Params,H,Delta,Orb,SQE,U,V)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: H, Delta
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnvector_t), intent(out) :: SQE
        type(o1b_t), intent(out) :: U, V
        type(o1b_t) :: Cc
        type(pnvector_t) :: Occ
        real(dp) :: Tol
        real(dp), parameter :: nTol = 1.0e-14_dp
        integer :: a_max, Pass

        Tol = Params%Calc%HFB%Tol
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(U%p(a_max,a_max), U%n(a_max,a_max), V%p(a_max,a_max), V%n(a_max,a_max), SQE%p(a_max), SQE%n(a_max))
        allocate(Occ%p(a_max), Occ%n(a_max), Cc%p(a_max,a_max), Cc%n(a_max,a_max))

        ! The original builds the canonical-basis-diagonal Delta after the first pass but then
        ! diagonalizes again with the unmodified Delta ... reproduced as is ...
        do Pass = 1, 2
            call bogoliubov_eig(H%p,Delta%p,a_max,U%p,V%p,SQE%p)
            call bogoliubov_eig(H%n,Delta%n,a_max,U%n,V%n,SQE%n)
            if (Pass == 2) then
                call clean(U%p,nTol)
                call clean(V%p,nTol)
                call clean(U%n,nTol)
                call clean(V%n,nTol)
            end if
            call HFB_orbital_ordering(Params,SQE,U,V,Orb,.false.)
            call mca_canonical(V%p,Occ%p,Cc%p)
            call mca_canonical(V%n,Occ%n,Cc%n)
            call HFB_canonical_orbital_ordering(Params,Occ,Cc,Orb)
            call clean_occ(Occ%p,Tol)
            call clean_occ(Occ%n,Tol)
        end do

        call bmz_amplitudes(Orb(1:a_max),Cc%p,Occ%p,U%p,V%p,.false.)
        call bmz_amplitudes(Orb(1:a_max),Cc%n,Occ%n,U%n,V%n,.false.)
        call clean(U%p,nTol)
        call clean(V%p,nTol)
        call clean(U%n,nTol)
        call clean(V%n,nTol)
    contains
        subroutine mca_canonical(Vm,Oc,Cm)
            real(dp), intent(in) :: Vm(:,:)
            real(dp), intent(out) :: Oc(:), Cm(:,:)
            real(dp), allocatable :: R(:,:)
            integer :: a
            R = block_rho(Vm,Orb(1:a_max))
            R = 0.5_dp * (R + transpose(R))
            do a = 1, a_max
                R(a,a) = R(a,a) + 1.0e-11_dp * real(a_max - a + 1,dp)
            end do
            call eig_sym(R,Oc,Cm,Descending=.true.)
        end subroutine mca_canonical
    end subroutine HFB_diagonalize_MCA

    subroutine HFB_diagonalize_BCS(Params,H,Delta,Orb,SQE,U,V)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: H, Delta
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnvector_t), intent(out) :: SQE
        type(o1b_t), intent(out) :: U, V
        type(o1b_t) :: Cc
        type(pnvector_t) :: SPE
        type(orb1b_t), allocatable :: Orb_tmp(:)
        real(dp), allocatable :: pDelta(:,:), nDelta(:,:)
        real(dp) :: Ph, pE_a, pD_a, nE_a, nD_a
        integer :: a_max, a

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(U,a_max)
        call o1b_alloc(V,a_max)
        allocate(SQE%p(a_max), SQE%n(a_max), SPE%p(a_max), SPE%n(a_max), Cc%p(a_max,a_max), Cc%n(a_max,a_max))
        call eig_sym(H%p,SPE%p,Cc%p)
        call eig_sym(H%n,SPE%n,Cc%n)
        ! dev HF_orbitals_make: same ordering as before & target occupations written into Orb (unused by HFB, a local
            ! copy here) ... fatal for HF.ESH with an open-shell target as dev ...
        Orb_tmp = Orb
        call HF_orbitals_make(Params,SPE,Cc,Orb_tmp)
        pDelta = mtmm(Cc%p,Delta%p)
        nDelta = mtmm(Cc%n,Delta%n)
        do a = 1, a_max
            Ph = phase(Orb(a)%l)
            pE_a = SPE%p(a)
            pD_a = pDelta(a,a)
            U%p(a,a) = Ph * sqrt(0.5_dp * (1.0_dp + pE_a / sqrt(pE_a**2 + pD_a**2)))
            V%p(a,a) = sqrt(0.5_dp * (1.0_dp - pE_a / sqrt(pE_a**2 + pD_a**2)))
            SQE%p(a) = sqrt(pE_a**2 + pD_a**2)
            nE_a = SPE%n(a)
            nD_a = nDelta(a,a)
            U%n(a,a) = Ph * sqrt(0.5_dp * (1.0_dp + nE_a / sqrt(nE_a**2 + nD_a**2)))
            V%n(a,a) = sqrt(0.5_dp * (1.0_dp - nE_a / sqrt(nE_a**2 + nD_a**2)))
            SQE%n(a) = sqrt(nE_a**2 + nD_a**2)
        end do
        U%p = mm(Cc%p,U%p)
        V%p = mm(Cc%p,V%p)
        U%n = mm(Cc%n,U%n)
        V%n = mm(Cc%n,V%n)
    end subroutine HFB_diagonalize_BCS

    ! Matches the quasiparticle vectors (U;V) to the (l,j)-blocks of the orbitals ... counterpart of HFB_orbital_ordering
        ! (dev HFB_ordering.jl) ... optional final sort by |U|^2 ascending (not for Pairing = 'BCS') ...
    subroutine HFB_orbital_ordering(Params,SQE,U,V,Orb,Final_Ordering)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(inout) :: SQE
        type(o1b_t), intent(inout) :: U, V
        type(orb1b_t), intent(in) :: Orb(:)
        logical, intent(in) :: Final_Ordering
        logical :: Sort_U
        Sort_U = Final_Ordering .and. Params%Calc%HFB%Pairing /= 'BCS'
        call order_one(SQE%p,U%p,V%p)
        call order_one(SQE%n,U%n,V%n)
    contains
        subroutine order_one(E,Um,Vm)
            real(dp), intent(inout) :: E(:), Um(:,:), Vm(:,:)
            real(dp), allocatable :: Nrm(:), UV(:,:), UV_new(:,:)
            integer, allocatable :: Perm(:), Ord(:)
            integer :: n, a
            n = size(E)
            if (Sort_U) then
                allocate(Nrm(n), Perm(n))
                do a = 1, n
                    Nrm(a) = sqrt(sum(Um(:,a)**2))**2
                end do
                call sortperm(Nrm,Perm)
                E = E(Perm)
                Um = Um(:,Perm)
                Vm = Vm(:,Perm)
            end if
            allocate(UV(2*n,n), Ord(n))
            UV(1:n,:) = Um
            UV(n+1:2*n,:) = Vm
            call HFB_orbital_matching(UV,Orb(1:n),Ord,UV_new)
            E = E(Ord)
            Um = UV_new(1:n,Ord)
            Vm = UV_new(n+1:2*n,Ord)
        end subroutine order_one
    end subroutine HFB_orbital_ordering

    ! Reorders canonical orbitals (columns of C) & occupations to the (l,j)-blocks ... HFB_canonical_orbital_ordering ...
    subroutine HFB_canonical_orbital_ordering(Params,Occ,C,Orb)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(inout) :: Occ
        type(o1b_t), intent(inout) :: C
        type(orb1b_t), intent(in) :: Orb(:)
        call canonical_one(Occ%p,C%p)
        call canonical_one(Occ%n,C%n)
    contains
        subroutine canonical_one(O,Cm)
            real(dp), intent(inout) :: O(:), Cm(:,:)
            real(dp), allocatable :: C_new(:,:)
            integer, allocatable :: Ord(:)
            allocate(Ord(size(O)))
            call HFB_orbital_matching(Cm,Orb(1:size(O)),Ord,C_new)
            O = O(Ord)
            Cm = C_new(:,Ord)
        end subroutine canonical_one
    end subroutine HFB_canonical_orbital_ordering

    ! Ordering of a canonical basis C ... HFB_canonical_basis_orbital_ordering ...
    subroutine HFB_canonical_basis_orbital_ordering(Cm,Orb)
        real(dp), intent(inout) :: Cm(:,:)
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), allocatable :: C_new(:,:)
        integer, allocatable :: Ord(:)
        allocate(Ord(size(Cm,2)))
        call HFB_orbital_matching(Cm,Orb(1:size(Cm,2)),Ord,C_new)
        Cm = C_new(:,Ord)
    end subroutine HFB_canonical_basis_orbital_ordering

    ! Matching of the vectors X(:,a) (rows: LHO components, a_max or 2*a_max rows, row r belongs to orbital
        ! mod(r-1,a_max)+1) to the (l,j)-blocks ... counterpart of HFB_orbital_matching (dev HFB_ordering.jl:79-185) ...
        ! weights W(B,a), argmax block, Hungarian assignment if the block counts do not fit, Order(b) = a (k-th vector of
        ! a block -> k-th orbital of the block, input order kept) & Lowdin projection of the vectors of mixed blocks ...
    subroutine HFB_orbital_matching(X,Orb,Order,X_new)
        real(dp), intent(in) :: X(:,:)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(out) :: Order(:)
        real(dp), allocatable, intent(out) :: X_new(:,:)
        real(dp), parameter :: Tol_pure = 1.0e-10_dp
        integer, allocatable :: Blk(:), Vector_block(:), Block_size(:), Block_count(:), Assignment(:)
        integer, allocatable :: Orbs(:), Vecs(:), Rows(:)
        real(dp), allocatable :: W(:,:), Cost(:,:), P(:,:), S(:,:), Z(:,:), Ev(:), Sm(:,:)
        integer :: a_max, Row_max, B_max, a, b, r, iB, k, nV, nR

        a_max = size(Orb)
        Row_max = size(X,1)

        ! (l,j)-blocks numbered by first appearance ...
        allocate(Blk(a_max))
        Blk = 0
        B_max = 0
        do a = 1, a_max
            do b = 1, a - 1
                if ((Orb(b)%l == Orb(a)%l) .and. (Orb(b)%j == Orb(a)%j)) then
                    Blk(a) = Blk(b)
                    exit
                end if
            end do
            if (Blk(a) == 0) then
                B_max = B_max + 1
                Blk(a) = B_max
            end if
        end do

        ! Normalized block weights of the vectors ...
        allocate(W(B_max,a_max))
        W = 0.0_dp
        do a = 1, a_max
            do r = 1, Row_max
                b = mod(r - 1,a_max) + 1
                W(Blk(b),a) = W(Blk(b),a) + X(r,a)**2
            end do
            W(:,a) = W(:,a) / sum(W(:,a))
        end do

        ! Block of the largest weight (first maximum) ...
        allocate(Vector_block(a_max))
        do a = 1, a_max
            Vector_block(a) = 1
            do iB = 2, B_max
                if (W(iB,a) > W(Vector_block(a),a)) Vector_block(a) = iB
            end do
        end do

        ! Wrong block populations ... Hungarian assignment of vectors a to orbitals b with the cost 1 - W(Block(b),a) ...
        allocate(Block_size(B_max), Block_count(B_max))
        do iB = 1, B_max
            Block_size(iB) = count(Blk == iB)
            Block_count(iB) = count(Vector_block == iB)
        end do
        if (any(Block_size /= Block_count)) then
            allocate(Cost(a_max,a_max), Assignment(a_max))
            do b = 1, a_max
                do a = 1, a_max
                    Cost(b,a) = 1.0_dp - W(Blk(b),a)
                end do
            end do
            call hungarian(Cost,Assignment)
            do b = 1, a_max
                Vector_block(Assignment(b)) = Blk(b)
            end do
        end if

        ! Order ... vectors of block B (input order) go to the orbitals of block B (LHO order) ...
        Order = 0
        do iB = 1, B_max
            Orbs = pack([(a, a = 1, a_max)], Blk == iB)
            Vecs = pack([(a, a = 1, a_max)], Vector_block == iB)
            Order(Orbs) = Vecs
        end do

        ! Lowdin projection of the vectors of mixed blocks P_new = P (P^T P)^(-1/2), zero outside the block ...
        X_new = X
        do iB = 1, B_max
            Vecs = pack([(a, a = 1, a_max)], Vector_block == iB)
            nV = size(Vecs)
            if (all(W(iB,Vecs) > (1.0_dp - Tol_pure))) cycle
            Rows = pack([(r, r = 1, Row_max)], [(Blk(mod(r - 1,a_max) + 1) == iB, r = 1, Row_max)])
            nR = size(Rows)
            P = X(Rows,Vecs)
            S = matmul(transpose(P),P)
            allocate(Ev(nV), Z(nV,nV))
            call eig_sym(S,Ev,Z)
            allocate(Sm(nV,nV))
            Sm = 0.0_dp
            do k = 1, nV
                Sm(:,k) = Z(:,k) / sqrt(Ev(k))
            end do
            Sm = matmul(Sm,transpose(Z))
            do k = 1, nV
                X_new(:,Vecs(k)) = 0.0_dp
            end do
            X_new(Rows,Vecs) = matmul(P,Sm)
            deallocate(Ev, Z, Sm)
        end do
    end subroutine HFB_orbital_matching

    function HFB_canonical_basis(Params,Rho,Orb) result(C)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t) :: C
        real(dp), allocatable :: pRho(:,:), nRho(:,:), w(:)
        real(dp) :: R
        integer :: N_max, a_max, a, b

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        pRho = 0.5_dp * (Rho%p + transpose(Rho%p))
        nRho = 0.5_dp * (Rho%n + transpose(Rho%n))
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%l /= Orb(b)%l .or. Orb(a)%j /= Orb(b)%j) then
                    pRho(a,b) = 0.0_dp
                    nRho(a,b) = 0.0_dp
                end if
            end do
            R = 10.0_dp * real(Orb(a)%l*(N_max+1) + (Orb(a)%j+1)/2,dp)
            pRho(a,a) = pRho(a,a) + R
            nRho(a,a) = nRho(a,a) + R
        end do
        allocate(w(a_max), C%p(a_max,a_max), C%n(a_max,a_max))
        call eig_sym(pRho,w,C%p,Descending=.true.)
        call eig_sym(nRho,w,C%n,Descending=.true.)
        call clean(C%p,1.0e-13_dp)
        call clean(C%n,1.0e-13_dp)
        call HFB_canonical_basis_orbital_ordering(C%p,Orb)
        call HFB_canonical_basis_orbital_ordering(C%n,Orb)
    end function HFB_canonical_basis

    function HFB_canonical_basis_BCS(Params,H,Orb) result(C)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: H
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t) :: C
        real(dp), allocatable :: pH(:,:), nH(:,:), w(:)
        integer :: a_max, a, b

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        pH = 0.5_dp * (H%p + transpose(H%p))
        nH = 0.5_dp * (H%n + transpose(H%n))
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%j /= Orb(b)%j .or. Orb(a)%l /= Orb(b)%l) then
                    pH(a,b) = 0.0_dp
                    nH(a,b) = 0.0_dp
                end if
            end do
        end do
        allocate(w(a_max), C%p(a_max,a_max), C%n(a_max,a_max))
        call eig_sym(pH,w,C%p)
        call eig_sym(nH,w,C%n)
        call clean(C%p,1.0e-13_dp)
        call clean(C%n,1.0e-13_dp)
        call HFB_canonical_basis_orbital_ordering(C%p,Orb)
        call HFB_canonical_basis_orbital_ordering(C%n,Orb)
    end function HFB_canonical_basis_BCS

    ! ---------------------------------------------------------------------------------------------------------------
    ! Modified Broyden mixing of (Rho, Kappa) ... counterpart of dev HFB_Broyden.jl (D. D. Johnson, PRB 38, 12807) ...
        ! one joint Jacobian for pRho, nRho, pKappa, nKappa ...

    ! Mapping (a <= b in the same (l,j)-block), weights & the initial input vector ... HFB_Broyden_initialize ...
    subroutine HFB_Broyden_initialize(Params,Rho,Kappa,Orb,B)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, Kappa
        type(orb1b_t), intent(in) :: Orb(:)
        type(broyden_t), intent(out) :: B
        integer :: a_max, a, b2, i, Pass

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        B%m_hist = Params%Calc%HFB%Broy_m
        do Pass = 1, 2
            i = 0
            do a = 1, a_max
                do b2 = a, a_max
                    if ((Orb(a)%l == Orb(b2)%l) .and. (Orb(a)%j == Orb(b2)%j)) then
                        i = i + 1
                        if (Pass == 2) B%Key(:,i) = [a, b2]
                    end if
                end do
            end do
            if (Pass == 1) then
                B%M = i
                allocate(B%Key(2,i))
            end if
        end do
        allocate(B%w(B%M))
        do i = 1, B%M
            B%w(i) = sqrt(real(Orb(B%Key(1,i))%j + 1,dp))
            if (B%Key(1,i) /= B%Key(2,i)) B%w(i) = B%w(i) * sqrt(2.0_dp)
        end do
        allocate(B%x(B%M,4), B%y(B%M,4), B%z(B%M,4), B%r(B%M,4), B%s(B%M,4))
        allocate(B%Xh(B%M,B%m_hist,4), B%Rh(B%M,B%m_hist,4))
        call HFB_Broyden_linearize(B,Rho,Kappa,B%x)
        B%y = 0.0_dp
        B%z = 0.0_dp
        B%r = 0.0_dp
        B%s = 0.0_dp
        B%Xh = 0.0_dp
        B%Rh = 0.0_dp
    end subroutine HFB_Broyden_initialize

    ! Weighted, (a,b)-(b,a) symmetrized elements of the (l,j)-blocks ... HFB_Broyden_linearize ...
    subroutine HFB_Broyden_linearize(B,Rho,Kappa,F)
        type(broyden_t), intent(in) :: B
        type(o1b_t), intent(in) :: Rho, Kappa
        real(dp), intent(out) :: F(:,:)
        integer :: i, a, b2
        do i = 1, B%M
            a = B%Key(1,i)
            b2 = B%Key(2,i)
            F(i,1) = 0.5_dp * B%w(i) * (Rho%p(a,b2) + Rho%p(b2,a))
            F(i,2) = 0.5_dp * B%w(i) * (Rho%n(a,b2) + Rho%n(b2,a))
            F(i,3) = 0.5_dp * B%w(i) * (Kappa%p(a,b2) + Kappa%p(b2,a))
            F(i,4) = 0.5_dp * B%w(i) * (Kappa%n(a,b2) + Kappa%n(b2,a))
        end do
    end subroutine HFB_Broyden_linearize

    ! Symmetric block matrices from a linearized vector (weights removed, zero elsewhere) ... HFB_Broyden_reconstruct ...
    subroutine HFB_Broyden_reconstruct(B,a_max,x,Rho,Kappa)
        type(broyden_t), intent(in) :: B
        integer, intent(in) :: a_max
        real(dp), intent(in) :: x(:,:)
        type(o1b_t), intent(out) :: Rho, Kappa
        integer :: i, a, b2
        call o1b_alloc(Rho,a_max)
        call o1b_alloc(Kappa,a_max)
        do i = 1, B%M
            a = B%Key(1,i)
            b2 = B%Key(2,i)
            Rho%p(a,b2) = x(i,1) / B%w(i)
            Rho%n(a,b2) = x(i,2) / B%w(i)
            Kappa%p(a,b2) = x(i,3) / B%w(i)
            Kappa%n(a,b2) = x(i,4) / B%w(i)
            Rho%p(b2,a) = Rho%p(a,b2)
            Rho%n(b2,a) = Rho%n(a,b2)
            Kappa%p(b2,a) = Kappa%p(a,b2)
            Kappa%n(b2,a) = Kappa%n(a,b2)
        end do
    end subroutine HFB_Broyden_reconstruct

    ! Joint norm of (pRho, nRho, pKappa, nKappa) ... HFB_Broyden_norm ...
    pure real(dp) function HFB_Broyden_norm(x) result(Nrm)
        real(dp), intent(in) :: x(:,:)
        Nrm = sqrt(sum(x(:,1)**2) + sum(x(:,2)**2) + sum(x(:,3)**2) + sum(x(:,4)**2))
    end function HFB_Broyden_norm

    ! Newest normalized history vectors (cyclic slot k = mod(Iteration,m) + 1) ... HFB_Broyden_history_update ...
    subroutine HFB_Broyden_history_update(Iteration,B)
        integer, intent(in) :: Iteration
        type(broyden_t), intent(inout) :: B
        real(dp), allocatable :: dR(:,:)
        real(dp) :: dR_norm
        integer :: k
        k = mod(Iteration,B%m_hist) + 1
        dR = B%r - B%s
        dR_norm = HFB_Broyden_norm(dR)
        if (dR_norm < 1.0e-14_dp) return
        B%Xh(:,k,:) = (B%x - B%y) / dR_norm
        B%Rh(:,k,:) = dR / dR_norm
    end subroutine HFB_Broyden_history_update

    ! Coefficients gamma = (w0^2 I + dR^T dR)^(-1) dR^T r ... HFB_Broyden_gamma (Julia Symmetric(A) \ c = dsysv 'U') ...
    function HFB_Broyden_gamma(B) result(gamma)
        type(broyden_t), intent(in) :: B
        real(dp) :: gamma(B%m_hist)
        real(dp), parameter :: w0 = 0.01_dp
        real(dp), allocatable :: A(:,:), c(:,:), Work(:)
        integer, allocatable :: Ipiv(:)
        integer :: k, m, Info, Lwork

        m = B%m_hist
        allocate(A(m,m), c(m,1), Ipiv(m))
        A = 0.0_dp
        c = 0.0_dp
        do k = 1, 4
            A = A + matmul(transpose(B%Rh(:,:,k)),B%Rh(:,:,k))
            c(:,1) = c(:,1) + matmul(transpose(B%Rh(:,:,k)),B%r(:,k))
        end do
        do k = 1, m
            A(k,k) = A(k,k) + w0**2
        end do
        Lwork = max(1,64*m)
        allocate(Work(Lwork))
        call dsysv('U',m,1,A,m,Ipiv,c,m,Work,Lwork,Info)
        gamma = c(:,1)
        if (Info /= 0 .or. .not. all(ieee_is_finite(gamma))) gamma = 0.0_dp
    end function HFB_Broyden_gamma

    ! Broyden update of the output densities (Rho, Kappa) of iteration Iteration ... HFB_Broyden_update ... plain
        ! iteration before Iteration_Broyden (1, or 17 with LNT where the history is restarted) ...
    subroutine HFB_Broyden_update(Params,Iteration,Rho,Kappa,B)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: Iteration
        type(o1b_t), intent(inout) :: Rho, Kappa
        type(broyden_t), intent(inout) :: B
        real(dp), allocatable :: F(:,:)
        real(dp) :: alpha, T_max, r_norm, s_norm
        real(dp), allocatable :: gamma(:)
        integer :: Iteration_Broyden, k, a_max

        a_max = size(Rho%p,1)
        alpha = Params%Calc%HFB%Broy_Amax
        T_max = Params%Calc%HFB%Broy_Tmax
        Iteration_Broyden = 1
        if (Params%Calc%HFB%LNT) Iteration_Broyden = 17

        allocate(F(B%M,4))
        call HFB_Broyden_linearize(B,Rho,Kappa,F)
        B%s = B%r
        B%r = F - B%x

        if (Iteration < Iteration_Broyden) then
            B%z = F
        else
            r_norm = HFB_Broyden_norm(B%r)
            s_norm = HFB_Broyden_norm(B%s)
            if ((Iteration == Iteration_Broyden) .and. (Iteration_Broyden > 1)) then
                B%Xh = 0.0_dp
                B%Rh = 0.0_dp
            else if (r_norm > s_norm / T_max) then
                B%Xh = 0.0_dp
                B%Rh = 0.0_dp
                call HFB_Broyden_history_update(Iteration,B)
            else
                call HFB_Broyden_history_update(Iteration,B)
            end if
            gamma = HFB_Broyden_gamma(B)
            do k = 1, 4
                B%z(:,k) = B%x(:,k) + alpha * B%r(:,k) - matmul(B%Xh(:,:,k) + alpha * B%Rh(:,:,k), gamma)
            end do
        end if

        ! Collapsed pairing ... keep the mixed Kappa exactly zero ...
        if (sqrt(sum(F(:,3)**2)) < 1.0e-12_dp) B%z(:,3) = 0.0_dp
        if (sqrt(sum(F(:,4)**2)) < 1.0e-12_dp) B%z(:,4) = 0.0_dp

        B%y = B%x
        B%x = B%z
        call HFB_Broyden_reconstruct(B,a_max,B%z,Rho,Kappa)
    end subroutine HFB_Broyden_update

    subroutine HFB_allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,H,Delta)
        type(params_t), intent(in) :: Params
        type(pnfloat_t), intent(in) :: Lambda
        type(o1b_t), intent(in) :: Rho, Kappa, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t), intent(out) :: H, Delta
        integer :: a_max, N_2max, N_3max, ad_i, ad_count, Pass, a, d, b, e, c, f, J, P_ab, P_abc, P_abe, P_adc
        integer :: n_a, l_a, j_a, n_d, l_d, j_d, n_b, l_b, j_b, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        integer, allocatable :: ad(:,:)
        real(dp) :: A_mass, p2, p3, j_a_hat, is_j_a_hat, J_j_a_hat, pHSum, nHSum, pRho_be, nRho_be, pRho_cf, nRho_cf
        real(dp) :: pKappa_be, nKappa_be, pKappa_cf, nKappa_cf, A_NNN_Amp, Amp, ME001, ME101, ME011, ME111, ME113
        real(dp) :: pD2, nD2, pD3, nD3, pDSum, nDSum
        character(len=32) :: CMS

        A_mass = real(Params%Calc%A,dp)
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        p2 = Params%Int%cP2N
        p3 = Params%Int%cP3N
        CMS = Params%Calc%CMS
        call o1b_alloc(H,a_max)
        call o1b_alloc(Delta,a_max)

        do Pass = 1, 2
            ad_count = 0
            do a = 1, a_max
                do d = 1, a
                    if (Orb(a)%j == Orb(d)%j .and. Orb(a)%l == Orb(d)%l) then
                        ad_count = ad_count + 1
                        if (Pass == 2) ad(:,ad_count) = [a,d]
                    end if
                end do
            end do
            if (Pass == 1) allocate(ad(2,ad_count))
        end do

        !$omp parallel do schedule(dynamic) private(a,d,n_a,l_a,j_a,n_d,l_d,j_d,j_a_hat,is_j_a_hat,pHSum,nHSum,b,e, &
        !$omp n_b,l_b,j_b,n_e,l_e,j_e,P_ab,pRho_be,nRho_be,J,J_j_a_hat,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f,pRho_cf,nRho_cf, &
        !$omp ME001,ME101,ME011,ME111,ME113,pKappa_be,nKappa_be,P_abe,A_NNN_Amp,pKappa_cf,nKappa_cf,Amp,pD2,nD2,pD3,nD3, &
        !$omp pDSum,nDSum,P_adc)
        do ad_i = 1, ad_count
            a = ad(1,ad_i)
            d = ad(2,ad_i)
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            n_d = Orb(d)%n
            l_d = Orb(d)%l
            j_d = Orb(d)%j
            j_a_hat = sqrt(real(j_a + 1,dp))
            is_j_a_hat = 1.0_dp / (real(j_a,dp) + 1.0_dp)

            ! Single-particle field H ...
            pHSum = 0.0_dp
            nHSum = 0.0_dp
            do b = 1, a_max
                do e = 1, a_max
                    n_b = Orb(b)%n
                    l_b = Orb(b)%l
                    j_b = Orb(b)%j
                    n_e = Orb(e)%n
                    l_e = Orb(e)%l
                    j_e = Orb(e)%j
                    P_ab = mod(l_a + l_b,2) + 1
                    if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                        if (l_b == l_e .and. j_b == j_e .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max) then
                            pRho_be = Rho%p(b,e)
                            nRho_be = Rho%n(b,e)
                            do J = abs(j_d - j_e)/2, (j_d + j_e)/2
                                if (mod(l_a + l_b, 2) == mod(l_d + l_e, 2)) then
                                    J_j_a_hat = real(2*J + 1,dp) * is_j_a_hat
                                    pHSum = pHSum + J_j_a_hat * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN) * pRho_be
                                    pHSum = pHSum + J_j_a_hat * O2b_pn(a,b,d,e,J,P_ab,V_NN,Orb_NN) * nRho_be
                                    nHSum = nHSum + J_j_a_hat * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN) * nRho_be
                                    nHSum = nHSum + J_j_a_hat * O2b_pn(b,a,e,d,J,P_ab,V_NN,Orb_NN) * pRho_be
                                end if
                                do c = 1, a_max
                                    n_c = Orb(c)%n
                                    l_c = Orb(c)%l
                                    if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                                    P_abc = mod(l_a + l_b + l_c, 2) + 1
                                    j_c = Orb(c)%j
                                    do f = 1, a_max
                                        n_f = Orb(f)%n
                                        l_f = Orb(f)%l
                                        if (.not. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                                   P_abc == (mod(l_d + l_e + l_f,2) + 1))) cycle
                                        j_f = Orb(f)%j
                                        if (j_c /= j_f) cycle
                                        pRho_cf = Rho%p(c,f)
                                        nRho_cf = Rho%n(c,f)
                                        ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P_abc,V_NNN,Orb,Orb_NNN)
                                        pHSum = pHSum + is_j_a_hat * (0.5_dp*ME113*pRho_be*pRho_cf + 0.25_dp * (ME001 + &
                                                sqrt(1.0_dp/3.0_dp)*ME101 + sqrt(1.0_dp/3.0_dp)*ME011 + ME111/3.0_dp + &
                                                2.0_dp/3.0_dp*ME113)*nRho_be*nRho_cf + 1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*pRho_be*nRho_cf)
                                        nHSum = nHSum + is_j_a_hat * (0.5_dp*ME113*nRho_be*nRho_cf + 0.25_dp * (ME001 + &
                                                sqrt(1.0_dp/3.0_dp)*ME101 + sqrt(1.0_dp/3.0_dp)*ME011 + ME111/3.0_dp + &
                                                2.0_dp/3.0_dp*ME113)*pRho_be*pRho_cf + 1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*nRho_be*pRho_cf)
                                    end do
                                end do
                            end do
                        end if
                    end if

                    ! Anomalous density-dependent part ... NNN only ...
                    if ((2*(n_b + n_e) + l_b + l_e) <= N_2max .and. (2*(n_a + n_b + n_e) + l_a + l_b + l_e) <= N_3max .and. &
                        l_b == l_e .and. j_b == j_e) then
                        pKappa_be = Kappa%p(b,e)
                        nKappa_be = Kappa%n(b,e)
                        P_abe = mod(l_a + l_b + l_e, 2) + 1
                        do c = 1, a_max
                            n_c = Orb(c)%n
                            l_c = Orb(c)%l
                            j_c = Orb(c)%j
                            A_NNN_Amp = 0.25_dp * sqrt(real((j_b + 1) * (j_c + 1),dp))
                            do f = 1, a_max
                                n_f = Orb(f)%n
                                l_f = Orb(f)%l
                                if (.not. ((2*(n_c + n_d + n_f) + l_c + l_d + l_f) <= N_3max .and. &
                                           P_abe == (mod(l_c + l_d + l_f,2) + 1) .and. l_c == l_f)) cycle
                                j_f = Orb(f)%j
                                if (j_c /= j_f) cycle
                                pKappa_cf = Kappa%p(c,f)
                                nKappa_cf = Kappa%n(c,f)
                                ME111 = V3b_no2b(b,e,a,1,c,f,d,1,0,1,P_abe,V_NNN,Orb,Orb_NNN)
                                ME113 = V3b_no2b(b,e,a,1,c,f,d,1,0,3,P_abe,V_NNN,Orb,Orb_NNN)
                                pHSum = pHSum + p3 * A_NNN_Amp * (ME113 * pKappa_be * pKappa_cf + &
                                        (2.0_dp * ME111 + ME113) / 3.0_dp * nKappa_be * nKappa_cf)
                                nHSum = nHSum + p3 * A_NNN_Amp * (ME113 * nKappa_be * nKappa_cf + &
                                        (2.0_dp * ME111 + ME113) / 3.0_dp * pKappa_be * pKappa_cf)
                            end do
                        end do
                    end if
                end do
            end do

            if (CMS == 'CMS1+2B') then
                pHSum = pHSum + T%p(a,d) * (1.0_dp - 1.0_dp / A_mass)
                nHSum = nHSum + T%n(a,d) * (1.0_dp - 1.0_dp / A_mass)
            else if (CMS /= 'CMS2B') then
                pHSum = pHSum + T%p(a,d)
                nHSum = nHSum + T%n(a,d)
            end if
            if (a == d) then
                pHSum = pHSum - Lambda%p
                nHSum = nHSum - Lambda%n
            end if
            H%p(a,d) = pHSum
            H%p(d,a) = pHSum
            H%n(a,d) = nHSum
            H%n(d,a) = nHSum

            ! Pairing field Delta ...
            if ((2*(n_a + n_d) + l_a + l_d) <= N_2max) then
                pDSum = 0.0_dp
                nDSum = 0.0_dp
                do b = 1, a_max
                    do e = 1, a_max
                        n_b = Orb(b)%n
                        l_b = Orb(b)%l
                        j_b = Orb(b)%j
                        n_e = Orb(e)%n
                        l_e = Orb(e)%l
                        j_e = Orb(e)%j
                        Amp = 0.5_dp * sqrt(real(j_b + 1,dp)) / j_a_hat
                        pD2 = 0.0_dp
                        nD2 = 0.0_dp
                        pD3 = 0.0_dp
                        nD3 = 0.0_dp
                        if ((2*(n_b + n_e) + l_b + l_e) <= N_2max .and. (l_b == l_e) .and. (j_b == j_e)) then
                            pKappa_be = Kappa%p(b,e)
                            nKappa_be = Kappa%n(b,e)
                            pD2 = pD2 + Amp * pKappa_be * O2b_pp(a,d,b,e,0,1,V_NN,Orb,Orb_NN)
                            nD2 = nD2 + Amp * nKappa_be * O2b_nn(a,d,b,e,0,1,V_NN,Orb,Orb_NN)
                            do c = 1, a_max
                                n_c = Orb(c)%n
                                l_c = Orb(c)%l
                                if ((2*(n_a + n_d + n_c) + l_a + l_d + l_c) > N_3max) cycle
                                P_adc = mod(l_a + l_d + l_c, 2) + 1
                                j_c = Orb(c)%j
                                do f = 1, a_max
                                    n_f = Orb(f)%n
                                    l_f = Orb(f)%l
                                    if (.not. ((2*(n_b + n_e + n_f) + l_b + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                               P_adc == (mod(l_b + l_e + l_f,2) + 1))) cycle
                                    j_f = Orb(f)%j
                                    if (j_c /= j_f) cycle
                                    pRho_cf = Rho%p(c,f)
                                    nRho_cf = Rho%n(c,f)
                                    ME111 = V3b_no2b(a,d,c,1,b,e,f,1,0,1,P_adc,V_NNN,Orb,Orb_NNN)
                                    ME113 = V3b_no2b(a,d,c,1,b,e,f,1,0,3,P_adc,V_NNN,Orb,Orb_NNN)
                                    pD3 = pD3 + Amp * (ME113 * pRho_cf + (2.0_dp * ME111 + ME113) / 3.0_dp * nRho_cf) * pKappa_be
                                    nD3 = nD3 + Amp * (ME113 * nRho_cf + (2.0_dp * ME111 + ME113) / 3.0_dp * pRho_cf) * nKappa_be
                                end do
                            end do
                        end if
                        pDSum = pDSum + (p2 * pD2 + p3 * pD3)
                        nDSum = nDSum + (p2 * nD2 + p3 * nD3)
                    end do
                end do
                Delta%p(a,d) = pDSum
                Delta%p(d,a) = pDSum
                Delta%n(a,d) = nDSum
                Delta%n(d,a) = nDSum
            end if
        end do
        !$omp end parallel do
    end subroutine HFB_allocate

    ! Initial Fermi-like occupations with particle-number adjustment & a Kappa guess ...
    subroutine HFB_density_operator_initialize(Params,Orb,Rho,Kappa)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(out) :: Rho, Kappa
        real(dp) :: Z_target, N_target, Z, N, pa0, na0, dZ, dN
        integer :: N_max, a_max, pN_val, nN_val, NN_a, NN_i, N_tot, a, b, i

        Z_target = real(Params%Calc%Z,dp)
        N_target = real(Params%Calc%A - Params%Calc%Z,dp)
        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2

        pN_val = 0
        nN_val = 0
        do NN_a = 0, N_max
            N_tot = 0
            do NN_i = 0, NN_a
                N_tot = N_tot + (NN_i + 1)*(NN_i + 2)
            end do
            if (real(N_tot,dp) >= Z_target) then
                pN_val = NN_a
                exit
            end if
        end do
        do NN_a = 0, N_max
            N_tot = 0
            do NN_i = 0, NN_a
                N_tot = N_tot + (NN_i + 1)*(NN_i + 2)
            end do
            if (real(N_tot,dp) >= N_target) then
                nN_val = NN_a
                exit
            end if
        end do

        call o1b_alloc(Rho,a_max)
        call o1b_alloc(Kappa,a_max)
        pa0 = Params%Calc%HFB%pa0
        na0 = Params%Calc%HFB%na0
        do a = 1, a_max
            NN_a = 2*Orb(a)%n + Orb(a)%l
            Rho%p(a,a) = exp(-real(NN_a - pN_val,dp) / pa0) / (1.0_dp + exp(-real(NN_a - pN_val,dp) / pa0))
            Rho%n(a,a) = exp(-real(NN_a - nN_val,dp) / na0) / (1.0_dp + exp(-real(NN_a - nN_val,dp) / na0))
        end do

        call HFB_particle_number(Params,Rho,Orb,Z,N)

        if (Z > Z_target) then
            do i = 1, 500
                if (.not. (Z > Z_target)) exit
                do a = 1, a_max
                    if ((Rho%p(a,a) - 0.001_dp) >= 1.0e-7_dp) then
                        dZ = 0.001_dp * real(Orb(a)%j + 1,dp)
                        if ((Z - dZ) > Z_target) then
                            Z = Z - dZ
                            Rho%p(a,a) = Rho%p(a,a) - 0.001_dp
                        end if
                    end if
                end do
            end do
        else if (Z < Z_target) then
            do i = 1, 500
                if (.not. (Z < Z_target)) exit
                do a = 1, a_max
                    NN_a = 2*Orb(a)%n + Orb(a)%l
                    if ((NN_a >= pN_val) .and. (NN_a < (pN_val + 3)) .and. ((Rho%p(a,a) + 0.001_dp) < 1.0_dp)) then
                        dZ = 0.001_dp * real(Orb(a)%j + 1,dp)
                        if ((Z + dZ) < Z_target) then
                            Z = Z + dZ
                            Rho%p(a,a) = Rho%p(a,a) + 0.001_dp
                        end if
                    end if
                end do
            end do
        end if

        if (N > N_target) then
            do i = 1, 500
                if (.not. (N > N_target)) exit
                do a = 1, a_max
                    if ((Rho%n(a,a) - 0.001_dp) >= 1.0e-7_dp) then
                        dN = 0.001_dp * real(Orb(a)%j + 1,dp)
                        if ((N - dN) > N_target) then
                            N = N - dN
                            Rho%n(a,a) = Rho%n(a,a) - 0.001_dp
                        end if
                    end if
                end do
            end do
        else if (N < N_target) then
            do i = 1, 500
                if (.not. (N < N_target)) exit
                do a = 1, a_max
                    NN_a = 2*Orb(a)%n + Orb(a)%l
                    if ((NN_a >= nN_val) .and. (NN_a < (nN_val + 3)) .and. ((Rho%n(a,a) + 0.001_dp) < 1.0_dp)) then
                        dN = 0.001_dp * real(Orb(a)%j + 1,dp)
                        if ((N + dN) < N_target) then
                            N = N + dN
                            Rho%n(a,a) = Rho%n(a,a) + 0.001_dp
                        end if
                    end if
                end do
            end do
        end if

        do a = 1, a_max
            Kappa%p(a,a) = phase(Orb(a)%l) * sqrt(Rho%p(a,a) - Rho%p(a,a)**2)
            Kappa%n(a,a) = phase(Orb(a)%l) * sqrt(Rho%n(a,a) - Rho%n(a,a)**2)
        end do
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%j == Orb(b)%j .and. Orb(a)%l == Orb(b)%l .and. a /= b) then
                    Kappa%p(a,b) = Kappa%p(a,b) + phase(Orb(a)%l) * sqrt(Kappa%p(a,a) * Kappa%p(b,b))
                    Kappa%n(a,b) = Kappa%n(a,b) + phase(Orb(a)%l) * sqrt(Kappa%n(a,a) * Kappa%n(b,b))
                end if
            end do
        end do
    end subroutine HFB_density_operator_initialize

    ! Rho = V V^T & Kappa = V U^T restricted to the (l,j) blocks ...
    subroutine HFB_density_operator(Params,U,V,Orb,Rho,Kappa)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: U, V
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(out) :: Rho, Kappa
        integer :: a_max, a, b, c
        real(dp) :: pR, nR, pK, nK

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(Rho,a_max)
        call o1b_alloc(Kappa,a_max)
        do a = 1, a_max
            do b = 1, a_max
                if (.not. ((Orb(a)%j == Orb(b)%j) .and. (Orb(a)%l == Orb(b)%l))) cycle
                pR = 0.0_dp
                nR = 0.0_dp
                pK = 0.0_dp
                nK = 0.0_dp
                do c = 1, a_max
                    if ((Orb(a)%j == Orb(c)%j) .and. (Orb(a)%l == Orb(c)%l)) then
                        pR = pR + V%p(a,c) * V%p(b,c)
                        pK = pK + V%p(a,c) * U%p(b,c)
                        nR = nR + V%n(a,c) * V%n(b,c)
                        nK = nK + V%n(a,c) * U%n(b,c)
                    end if
                end do
                Rho%p(a,b) = pR
                Kappa%p(a,b) = pK
                Rho%n(a,b) = nR
                Kappa%n(a,b) = nK
            end do
        end do
    end subroutine HFB_density_operator

    ! Safeguarded secant step for protons & neutrons ... counterpart of HFB_chemical_potential_secant (dev
        ! HFB_chemical_potential.jl) ... the brackets (index 1 p, 2 n) are updated, Has_* = .false. is dev's NaN ...
    function HFB_chemical_potential_secant(Params,Lambda_1,dA_1,Lambda_2,dA_2,Lambda_lo,Lambda_hi,Has_lo,Has_hi) result(L)
        type(params_t), intent(in) :: Params
        type(pnfloat_t), intent(in) :: Lambda_1, dA_1, Lambda_2, dA_2
        real(dp), intent(inout) :: Lambda_lo(2), Lambda_hi(2)
        logical, intent(inout) :: Has_lo(2), Has_hi(2)
        type(pnfloat_t) :: L
        L%p = HFB_chemical_potential_secant_step(Lambda_1%p,dA_1%p,Lambda_2%p,dA_2%p,Lambda_lo(1),Lambda_hi(1), &
                                                 Has_lo(1),Has_hi(1),Params%Calc%HFB%dLmax,Params%Calc%HFB%Tol)
        L%n = HFB_chemical_potential_secant_step(Lambda_1%n,dA_1%n,Lambda_2%n,dA_2%n,Lambda_lo(2),Lambda_hi(2), &
                                                 Has_lo(2),Has_hi(2),Params%Calc%HFB%dLmax,Params%Calc%HFB%Tol)
    end function HFB_chemical_potential_secant

    ! dA = A_target - A(Lambda), A increasing in Lambda ... (L2,dA2) is the latest point ...
    function HFB_chemical_potential_secant_step(L1,dA1,L2,dA2,Lo,Hi,Has_lo,Has_hi,dLmax,Tol) result(L)
        real(dp), intent(in) :: L1, dA1, L2, dA2, dLmax, Tol
        real(dp), intent(inout) :: Lo, Hi
        logical, intent(inout) :: Has_lo, Has_hi
        real(dp) :: L, Denom, dL
        if (abs(dA2) < Tol) then
            L = L2
            return
        end if
        if (dA2 > 0.0_dp) then
            Lo = L2
            Has_lo = .true.
        else
            Hi = L2
            Has_hi = .true.
        end if
        Denom = dA2 - dA1
        dL = 0.0_dp
        if (Denom /= 0.0_dp .and. L1 /= L2 .and. ieee_is_finite(Denom)) dL = -dA2 * (L2 - L1) / Denom
        if (.not. (dL * dA2 > 0.0_dp)) dL = julia_sign(dA2) * dLmax
        if (abs(dL) > dLmax) dL = julia_sign(dL) * dLmax
        L = L2 + dL
        if (Has_lo .and. Has_hi) then
            if (L <= Lo .or. L >= Hi) L = 0.5_dp * (Lo + Hi)
        end if
    contains
        pure real(dp) function julia_sign(x)
            real(dp), intent(in) :: x
            julia_sign = 0.0_dp
            if (x > 0.0_dp) julia_sign = 1.0_dp
            if (x < 0.0_dp) julia_sign = -1.0_dp
        end function julia_sign
    end function HFB_chemical_potential_secant_step

    ! Chemical potential in the middle of the shell gap for a collapsed pairing field ... HFB_chemical_potential_gap ...
        ! Has_gap = .false. is dev's NaN (not a closed shell, Lambda from the secant method) ...
    subroutine HFB_chemical_potential_gap(Params,H,Lambda_H,Delta,Orb,Lambda_gap,Has_gap)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: H, Delta
        type(pnfloat_t), intent(in) :: Lambda_H
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(out) :: Lambda_gap(2)
        logical, intent(out) :: Has_gap(2)
        real(dp) :: Tol
        Tol = Params%Calc%HFB%Tol
        Lambda_gap = 0.0_dp
        Has_gap = .false.
        if (maxval(abs(Delta%p)) < Tol) call HFB_chemical_potential_gap_species(H%p,Lambda_H%p,real(Params%Calc%Z,dp), &
                                                                                 Orb,Lambda_gap(1),Has_gap(1))
        if (maxval(abs(Delta%n)) < Tol) call HFB_chemical_potential_gap_species(H%n,Lambda_H%n, &
                                                 real(Params%Calc%A - Params%Calc%Z,dp),Orb,Lambda_gap(2),Has_gap(2))
    end subroutine HFB_chemical_potential_gap

    ! Levels of h = H + Lambda_H per (l,j) block (first appearance order, ascending eigenvalues), stable sort, filling ...
    subroutine HFB_chemical_potential_gap_species(Hm,Lambda_H,A_target,Orb,Lambda,Found)
        real(dp), intent(in) :: Hm(:,:), Lambda_H, A_target
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(out) :: Lambda
        logical, intent(out) :: Found
        real(dp), allocatable :: SPE(:), Deg(:), Hb(:,:), Ev(:), Zv(:,:)
        integer, allocatable :: Blk(:), Order(:)
        integer :: a_max, a, b, n, nS, i
        logical :: First
        real(dp) :: A_filled

        a_max = size(Orb)
        allocate(SPE(a_max), Deg(a_max))
        nS = 0
        do a = 1, a_max
            First = .true.
            do b = 1, a - 1
                if (Orb(b)%l == Orb(a)%l .and. Orb(b)%j == Orb(a)%j) then
                    First = .false.
                    exit
                end if
            end do
            if (.not. First) cycle
            Blk = pack([(b, b = 1, a_max)], [((Orb(b)%l == Orb(a)%l .and. Orb(b)%j == Orb(a)%j), b = 1, a_max)])
            n = size(Blk)
            Hb = Hm(Blk,Blk)
            if (allocated(Ev)) deallocate(Ev, Zv)
            allocate(Ev(n), Zv(n,n))
            call eig_sym(Hb,Ev,Zv)
            do i = 1, n
                nS = nS + 1
                SPE(nS) = Ev(i) + Lambda_H
                Deg(nS) = real(Orb(a)%j + 1,dp)
            end do
        end do

        Lambda = 0.0_dp
        Found = .false.
        allocate(Order(nS))
        call sortperm(SPE(1:nS),Order)
        A_filled = 0.0_dp
        do i = 1, nS - 1
            A_filled = A_filled + Deg(Order(i))
            if (abs(A_filled - A_target) < 1.0e-8_dp) then
                Lambda = 0.5_dp * (SPE(Order(i)) + SPE(Order(i+1)))
                Found = .true.
                return
            else if (A_filled > A_target) then
                exit
            end if
        end do
    end subroutine HFB_chemical_potential_gap_species

    subroutine HFB_particle_number(Params,Rho,Orb,Z,N)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(out) :: Z, N
        integer :: a
        Z = 0.0_dp
        N = 0.0_dp
        do a = 1, (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
            Z = Z + Rho%p(a,a) * real(Orb(a)%j + 1,dp)
            N = N + Rho%n(a,a) * real(Orb(a)%j + 1,dp)
        end do
    end subroutine HFB_particle_number

    function HFB_particle_number_dispersion(Params,Rho,Orb) result(dA)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnfloat_t) :: dA
        real(dp) :: Tol, dZ, dN, jh
        integer :: a_max, a, b

        Tol = min(1.0e-6_dp, Params%Calc%HFB%Tol)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(/,a)') 'Calculating the HFB dispersion of proton & neutron particle numbers ...'
        dZ = 0.0_dp
        dN = 0.0_dp
        do a = 1, a_max
            jh = real(Orb(a)%j + 1,dp)
            do b = 1, a_max
                if ((Orb(a)%j == Orb(b)%j) .and. (Orb(a)%l == Orb(b)%l)) then
                    dZ = dZ - 2.0_dp * jh * Rho%p(a,b) * Rho%p(b,a)
                    dN = dN - 2.0_dp * jh * Rho%n(a,b) * Rho%n(b,a)
                end if
            end do
            dZ = dZ + 2.0_dp * jh * Rho%p(a,a)
            dN = dN + 2.0_dp * jh * Rho%n(a,a)
        end do
        if (dZ < 0.0_dp .and. abs(dZ) < Tol) dZ = 0.0_dp
        if (dN < 0.0_dp .and. abs(dN) < Tol) dN = 0.0_dp
        dZ = sqrt(abs(dZ))
        dN = sqrt(abs(dN))
        write(*,'(/,a)') 'HFB dispersion of nucleon numbers are ...'
        write(*,'(a)') TAB // 'dZ = ' // jstr(round_digits(dZ,5))
        write(*,'(a)') TAB // 'dN = ' // jstr(round_digits(dN,5))
        dA = pnfloat_t(dZ,dN)
    end function HFB_particle_number_dispersion

    ! Lipkin-Nogami correction of the single-particle field H (Delta is left untouched as in the original) ...
    subroutine HFB_Lipkin_Nogami(Params,Orb,Orb_NN,Orb_NNN,Rho,Kappa,H,V_NN,V_NNN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: Rho, Kappa
        type(o1b_t), intent(inout) :: H
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(pnfloat_t) :: Lambda_2
        integer :: a_max, a, u

        Lambda_2 = HFB_Lipkin_Nogami_Lambda_2(Params,Orb,Orb_NN,Orb_NNN,Rho,Kappa,V_NN,V_NNN)
        open(newunit=u, file=io_path(Params,'/HFB/bin/lambda_2.bin'), access='stream', form='unformatted', status='replace', action='write')
        write(u) Lambda_2%p, Lambda_2%n
        close(u)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        H%p = H%p + 4.0_dp * Lambda_2%p * Rho%p
        H%n = H%n + 4.0_dp * Lambda_2%n * Rho%n
        do a = 1, a_max
            H%p(a,a) = H%p(a,a) - 2.0_dp * Lambda_2%p
            H%n(a,a) = H%n(a,a) - 2.0_dp * Lambda_2%n
        end do
    end subroutine HFB_Lipkin_Nogami

    function HFB_Lipkin_Nogami_Lambda_2_import(Params) result(Lambda_2)
        type(params_t), intent(in) :: Params
        type(pnfloat_t) :: Lambda_2
        integer :: u
        open(newunit=u, file=io_path(Params,'/HFB/bin/lambda_2.bin'), access='stream', form='unformatted', status='old', action='read')
        read(u) Lambda_2%p, Lambda_2%n
        close(u)
    end function HFB_Lipkin_Nogami_Lambda_2_import

    function HFB_Lipkin_Nogami_Lambda_2(Params,Orb,Orb_NN,Orb_NNN,Rho,Kappa,V_NN,V_NNN) result(Lambda_2)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: Rho, Kappa
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(pnfloat_t) :: Lambda_2
        real(dp), allocatable :: pRb(:,:), nRb(:,:), pKb(:,:), nKb(:,:), pM2(:,:), nM2(:,:), Tr_a(:,:)
        real(dp) :: p2, p3, Tol, pTr1, pTr2, nTr1, nTr2, pRho_ad, nRho_ad, pRho_be, nRho_be, pRho_cf, nRho_cf, J_hat
        real(dp) :: ME111, ME113, pKappa_ad, nKappa_ad, pKappa_be, nKappa_be, Amp, j_a_hat, pDen, nDen, pTS, nTS
        integer :: a_max, N_2max, N_3max, a, b, c, d, e, f, J, P_ab, P_abc, P_adc
        integer :: n_a, l_a, j_a, n_d, l_d, j_d, n_b, l_b, j_b, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f

        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        p2 = Params%Int%cP2N
        p3 = Params%Int%cP3N
        Tol = Params%Calc%HFB%Tol

        allocate(pRb(a_max,a_max), nRb(a_max,a_max), pKb(a_max,a_max), nKb(a_max,a_max))
        pRb = 0.0_dp
        nRb = 0.0_dp
        pKb = 0.0_dp
        nKb = 0.0_dp
        do a = 1, a_max
            do b = 1, a_max
                if (.not. (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j)) cycle
                do c = 1, a_max
                    if (Orb(a)%l == Orb(c)%l .and. Orb(a)%j == Orb(c)%j) then
                        pRb(a,b) = pRb(a,b) - Rho%p(a,c) * Rho%p(c,b)
                        nRb(a,b) = nRb(a,b) - Rho%n(a,c) * Rho%n(c,b)
                        pKb(a,b) = pKb(a,b) + Rho%p(a,c) * Kappa%p(c,b)
                        nKb(a,b) = nKb(a,b) + Rho%n(a,c) * Kappa%n(c,b)
                    end if
                end do
                pRb(a,b) = pRb(a,b) + Rho%p(a,b)
                nRb(a,b) = nRb(a,b) + Rho%n(a,b)
            end do
        end do

        ! Only diagonal (a == d, b == e) terms contribute to the traces ...
        pTr1 = 0.0_dp
        pTr2 = 0.0_dp
        nTr1 = 0.0_dp
        nTr2 = 0.0_dp
        ! Deterministic (thread-count independent) traces: per-a partial sums, added in a order afterwards (the LN field
            ! feeds back into the iteration, a reduction would make 1 vs 8 threads differ) ...
        allocate(Tr_a(4,a_max))
        !$omp parallel do schedule(dynamic) private(pTr1,pTr2,nTr1,nTr2,d,n_a,l_a,j_a,n_d,l_d,j_d,j_a_hat, &
        !$omp pRho_ad,nRho_ad,b,e,n_b,l_b,j_b,n_e,l_e,j_e,P_ab,pRho_be,nRho_be,J,J_hat,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f, &
        !$omp pRho_cf,nRho_cf,ME111,ME113,pKappa_ad,nKappa_ad,Amp,pKappa_be,nKappa_be,P_adc)
        do a = 1, a_max
            d = a
            pTr1 = 0.0_dp
            pTr2 = 0.0_dp
            nTr1 = 0.0_dp
            nTr2 = 0.0_dp
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            n_d = n_a
            l_d = l_a
            j_d = j_a
            j_a_hat = sqrt(real(j_a + 1,dp))
            pRho_ad = pRb(a,d)
            nRho_ad = nRb(a,d)
            do b = 1, a_max
                e = b
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                n_e = n_b
                l_e = l_b
                j_e = j_b
                P_ab = mod(l_a + l_b,2) + 1
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max) then
                    pRho_be = pRb(b,e)
                    nRho_be = nRb(b,e)
                    do J = abs(j_d - j_e)/2, (j_d + j_e)/2
                        J_hat = real(2*J + 1,dp)
                        pTr1 = pTr1 + J_hat * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN) * pRho_ad * pRho_be
                        nTr1 = nTr1 + J_hat * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN) * nRho_ad * nRho_be
                        do c = 1, a_max
                            n_c = Orb(c)%n
                            l_c = Orb(c)%l
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                            P_abc = mod(l_a + l_b + l_c, 2) + 1
                            j_c = Orb(c)%j
                            do f = 1, a_max
                                n_f = Orb(f)%n
                                l_f = Orb(f)%l
                                if (.not. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                           P_abc == (mod(l_d + l_e + l_f,2) + 1))) cycle
                                j_f = Orb(f)%j
                                if (j_c /= j_f) cycle
                                pRho_cf = Rho%p(c,f)
                                nRho_cf = Rho%n(c,f)
                                ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P_abc,V_NNN,Orb,Orb_NNN)
                                pTr1 = pTr1 + (ME113 * pRho_cf + (2.0_dp * ME111 + ME113) / 3.0_dp * nRho_cf) * pRho_ad * pRho_be
                                nTr1 = nTr1 + (ME113 * nRho_cf + (2.0_dp * ME111 + ME113) / 3.0_dp * pRho_cf) * nRho_ad * nRho_be
                            end do
                        end do
                    end do
                end if
            end do

            if ((2*(n_a + n_d) + l_a + l_d) <= N_2max) then
                pKappa_ad = Kappa%p(a,d) - pKb(a,d)
                nKappa_ad = Kappa%n(a,d) - nKb(a,d)
                do b = 1, a_max
                    e = b
                    n_b = Orb(b)%n
                    l_b = Orb(b)%l
                    j_b = Orb(b)%j
                    n_e = n_b
                    l_e = l_b
                    Amp = sqrt(real(j_b + 1,dp)) * j_a_hat
                    if ((2*(n_b + n_e) + l_b + l_e) > N_2max) cycle
                    pKappa_be = pKb(b,e)
                    nKappa_be = nKb(b,e)
                    pTr2 = pTr2 - p2 * Amp * pKappa_ad * pKappa_be * O2b_pp(a,d,b,e,0,1,V_NN,Orb,Orb_NN)
                    nTr2 = nTr2 - p2 * Amp * nKappa_ad * nKappa_be * O2b_nn(a,d,b,e,0,1,V_NN,Orb,Orb_NN)
                    do c = 1, a_max
                        n_c = Orb(c)%n
                        l_c = Orb(c)%l
                        if ((2*(n_a + n_d + n_c) + l_a + l_d + l_c) > N_3max) cycle
                        P_adc = mod(l_a + l_d + l_c, 2) + 1
                        j_c = Orb(c)%j
                        do f = 1, a_max
                            n_f = Orb(f)%n
                            l_f = Orb(f)%l
                            if (.not. ((2*(n_b + n_e + n_f) + l_b + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                       P_adc == (mod(l_b + l_e + l_f,2) + 1))) cycle
                            j_f = Orb(f)%j
                            if (j_c /= j_f) cycle
                            pRho_cf = Rho%p(c,f)
                            nRho_cf = Rho%n(c,f)
                            ME111 = V3b_no2b(a,d,c,1,b,e,f,1,0,1,P_adc,V_NNN,Orb,Orb_NNN)
                            ME113 = V3b_no2b(a,d,c,1,b,e,f,1,0,3,P_adc,V_NNN,Orb,Orb_NNN)
                            pTr2 = pTr2 - p3 * Amp * (ME113 * pRho_cf + (2.0_dp * ME111 + ME113) / 3.0_dp * nRho_cf) * pKappa_ad * pKappa_be
                            nTr2 = nTr2 - p3 * Amp * (ME113 * nRho_cf + (2.0_dp * ME111 + ME113) / 3.0_dp * pRho_cf) * nKappa_ad * nKappa_be
                        end do
                    end do
                end do
            end if
            Tr_a(:,a) = [pTr1, pTr2, nTr1, nTr2]
        end do
        !$omp end parallel do
        pTr1 = 0.0_dp
        pTr2 = 0.0_dp
        nTr1 = 0.0_dp
        nTr2 = 0.0_dp
        do a = 1, a_max
            pTr1 = pTr1 + Tr_a(1,a)
            pTr2 = pTr2 + Tr_a(2,a)
            nTr1 = nTr1 + Tr_a(3,a)
            nTr2 = nTr2 + Tr_a(4,a)
        end do

        pM2 = mm(pRb,pRb)
        nM2 = mm(nRb,nRb)
        pDen = 0.0_dp
        nDen = 0.0_dp
        pTS = 0.0_dp
        nTS = 0.0_dp
        do a = 1, a_max
            pDen = pDen - 4.0_dp * real(Orb(a)%j + 1,dp) * pM2(a,a)
            pTS = pTS + real(Orb(a)%j + 1,dp) * pRb(a,a)
            nDen = nDen - 4.0_dp * real(Orb(a)%j + 1,dp) * nM2(a,a)
            nTS = nTS + real(Orb(a)%j + 1,dp) * nRb(a,a)
        end do
        pDen = pDen + 2.0_dp * pTS**2
        nDen = nDen + 2.0_dp * nTS**2
        Lambda_2 = pnfloat_t(0.0_dp, 0.0_dp)
        if (abs(pDen) > max(Tol,1.0e-4_dp)) Lambda_2%p = (pTr1 + pTr2) / pDen
        if (abs(nDen) > max(Tol,1.0e-4_dp)) Lambda_2%n = (nTr1 + nTr2) / nDen
    end function HFB_Lipkin_Nogami_Lambda_2

    real(dp) function HFB_Lipkin_Nogami_energy(Params,Orb,Lambda_2,Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnfloat_t), intent(in) :: Lambda_2
        type(o1b_t), intent(in) :: Rho
        real(dp), allocatable :: pR2(:,:), nR2(:,:)
        integer :: a
        pR2 = mm(Rho%p,Rho%p)
        nR2 = mm(Rho%n,Rho%n)
        HFB_Lipkin_Nogami_energy = 0.0_dp
        do a = 1, (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
            HFB_Lipkin_Nogami_energy = HFB_Lipkin_Nogami_energy - 2.0_dp * real(Orb(a)%j + 1,dp) * Lambda_2%p * (Rho%p(a,a) - pR2(a,a))
            HFB_Lipkin_Nogami_energy = HFB_Lipkin_Nogami_energy - 2.0_dp * real(Orb(a)%j + 1,dp) * Lambda_2%n * (Rho%n(a,a) - nR2(a,a))
        end do
    end function HFB_Lipkin_Nogami_energy

    subroutine HFB_Lipkin_Nogami_V2b_residual_no2b(Params,Orb,Orb_NN,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(inout) :: V_NN
        type(pnfloat_t) :: Lambda_2
        integer :: J, P, N_T1, Bra, Ket, Ind, a, b, c, d
        real(dp) :: Ph, ME

        Lambda_2 = HFB_Lipkin_Nogami_Lambda_2_import(Params)
        write(*,'(/,a,/)') 'Adding the LN N^2 correction to the residual NN interaction ...'
        do J = 0, Params%Calc%N2max + 1
            do P = 1, 2
                N_T1 = Orb_NN%N(2,P,J)
                do Bra = 1, N_T1
                    a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                    b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                    do Ket = 1, Bra
                        Ind = packed_index(Bra,Ket,N_T1)
                        c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                        d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                        Ph = phase((Orb(c)%j + Orb(d)%j)/2 - J)
                        ME = -2.0_dp * Params%Int%cRes * real(kronecker_delta(a,c) * kronecker_delta(b,d) - &
                             nint(Ph) * kronecker_delta(a,d) * kronecker_delta(b,c),dp) / real(1 + kronecker_delta(a,b),dp)
                        V_NN%pp(P,J)%v(Ind) = V_NN%pp(P,J)%v(Ind) + Lambda_2%p * ME
                        V_NN%nn(P,J)%v(Ind) = V_NN%nn(P,J)%v(Ind) + Lambda_2%n * ME
                    end do
                end do
            end do
        end do
    end subroutine HFB_Lipkin_Nogami_V2b_residual_no2b

    ! Total HFB energy ... [mean field, NN pairing, NNN pairing, Lipkin-Nogami] ...
    function HFB_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Quiet) result(E_HFB)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, Kappa, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        logical, intent(in) :: Quiet
        real(dp) :: E_HFB(4)

        if (.not. Quiet) write(*,'(/,a)') 'Calculating total HFB ground-state energy ...'
        E_HFB = 0.0_dp
        call pairing_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,E_HFB(1),E_HFB(2),E_HFB(3),.true.)
        if (Params%Calc%HFB%LNT) E_HFB(4) = HFB_Lipkin_Nogami_energy(Params,Orb,HFB_Lipkin_Nogami_Lambda_2_import(Params),Rho)
        if (Quiet) return
        write(*,'(/,a)') TAB // 'HFB energy                 ...   E_HFB = ' // jstr(sum(E_HFB)) // ' MeV'
        write(*,'(a)') TAB // 'Mean-field energy          ...   E_MF  = ' // jstr(E_HFB(1)) // ' MeV'
        write(*,'(a)') TAB // 'Pairing energy             ...   E_Par = ' // jstr(E_HFB(2) + E_HFB(3)) // ' MeV'
        write(*,'(a)') TAB // 'NN pairing energy        ...   E_Par2N = ' // jstr(E_HFB(2)) // ' MeV'
        write(*,'(a)') TAB // 'NNN pairing energy       ...   E_Par3N = ' // jstr(E_HFB(3)) // ' MeV'
        if (Params%Calc%HFB%LNT) write(*,'(a)') TAB // 'Lipkin-Nogami energy     ...   E_LN  = ' // jstr(E_HFB(4)) // ' MeV'
    end function HFB_energy

    subroutine HFB_summary(Params,E_HFB,T_HFB,Lambda,dA,Convergence,Iteration)
        type(params_t), intent(in) :: Params
        real(dp), intent(in) :: E_HFB(4), T_HFB
        type(pnfloat_t), intent(in) :: Lambda, dA
        logical, intent(in) :: Convergence
        integer, intent(in) :: Iteration
        type(pnfloat_t) :: Lambda_2
        real(dp) :: A
        integer :: u, N_max
        character(len=*), parameter :: M = TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' '
        character(len=*), parameter :: MN = TAB // ' MeV / Nucleon ' // TAB // TAB // ' ... ' // TAB // ' '

        A = real(Params%Calc%A,dp)
        N_max = Params%Calc%Nmax
        write(*,'(/,a)') 'Exporting HFB calculation summary ...'
        open(newunit=u, file=io_path(Params,'/HFB/HFB_summary.dat'), position='append', action='write')
        write(u,'(a)') 'Spherical Hartree-Fock-Bogoliubov mean-field calculation summary'
        write(u,'(/,a)') 'Target nuclid:       A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)
        write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                         ' , N_2max = ' // istr(Params%Calc%N2max) // ' , N_3max = ' // istr(Params%Calc%N3max) // &
                         ', J-basis size = ' // istr((N_max+1)*(N_max+2)/2) // ', M-basis size = ' // istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Combined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Only pure 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else
            write(u,'(/,a)') TAB // TAB // TAB // 'No Center-of-Mass System (CMS) motion correction is set ...'
        end if
        write(u,'(/,a)') 'Pairing interaction parameters:'
        write(u,'(a)') 'cP2N     = ' // jstr(Params%Int%cP2N)
        write(u,'(a)') 'cP3N     = ' // jstr(Params%Int%cP3N)
        if (Convergence) then
            write(u,'(/,a)') 'HFB equations converged in ' // istr(Iteration) // ' iterations with precision epsilon = ' // &
                             jstr(Params%Calc%HFB%Tol) // ' ...'
        else
            write(u,'(/,/,a)') 'HFB equations did NOT converge in ' // istr(Iteration) // ' iterations with precision epsilon = ' // &
                               jstr(Params%Calc%HFB%Tol) // ' ...'
        end if
        write(u,'(/,a)') 'Spherical HFB solution review:'
        write(u,'(/,a)') 'Number of HFB Iterations = ' // cfmt_d(Iteration,5) // ', Convergence precision =  ' // &
                         cfmt_e(Params%Calc%HFB%Tol,0,3) // ' (dE <-> MeV, dN <-> Nucleons)'
        write(u,'(/,a)') 'E_HFB      = ' // cfmt_f(sum(E_HFB),15,8) // M // 'HFB mean-field ground-state energy'
        write(u,'(a)') 'E_MF       = ' // cfmt_f(E_HFB(1),15,8) // M // 'mean-field ground-state energy'
        write(u,'(a)') 'E_Par      = ' // cfmt_f(E_HFB(2) + E_HFB(3),15,8) // M // 'pairing ground-state energy'
        write(u,'(a)') 'E_Par (2N) = ' // cfmt_f(E_HFB(2),15,8) // M // '2-body pairing ground-state energy'
        write(u,'(a)') 'E_Par (3N) = ' // cfmt_f(E_HFB(3),15,8) // M // '3-body pairing ground-state energy'
        if (Params%Calc%HFB%LNT) write(u,'(a)') 'E_LN       = ' // cfmt_f(E_HFB(4),15,8) // M // &
                                                'HFB mean-field Lipkin-Nogami ground-state energy correction'
        write(u,'(a)') 'T_HFB      = ' // cfmt_f(T_HFB,15,8) // TAB // ' MeV' // TAB // TAB // ' ... ' // TAB // &
                       ' HFB mean-field ground-state kinetic energy'
        write(u,'(/,a)') 'E_HFB / A      = ' // cfmt_f(sum(E_HFB) / A,15,8) // MN // 'HFB mean-field ground-state energy per nucleon'
        write(u,'(a)') 'E_MF  / A      = ' // cfmt_f(E_HFB(1) / A,15,8) // MN // 'mean-field ground-state energy per nucleon'
        write(u,'(a)') 'E_Par / A      = ' // cfmt_f((E_HFB(2) + E_HFB(3)) / A,15,8) // MN // 'pairing ground-state energy per nucleon'
        write(u,'(a)') 'E_Par / A (2N) = ' // cfmt_f(E_HFB(2) / A,15,8) // MN // '2-body pairing ground-state energy per nucleon'
        write(u,'(a)') 'E_Par / A (3N) = ' // cfmt_f(E_HFB(3) / A,15,8) // MN // '3-body pairing ground-state energy per nucleon'
        write(u,'(a)') 'T_HFB / A      = ' // cfmt_f(T_HFB / A,15,8) // TAB // ' MeV / Nucleon' // TAB // TAB // ' ... ' // TAB // &
                       ' HFB mean-field ground-state kinetic energy per nucleon'
        write(u,'(/,/,a)') 'Values of HFB chemical potentials:'
        write(u,'(/,a)') 'pLambda =  ' // cfmt_f(Lambda%p,15,8) // ' MeV'
        write(u,'(a)') 'nLambda =  ' // cfmt_f(Lambda%n,15,8) // ' MeV'
        if (Params%Calc%HFB%LNT) then
            Lambda_2 = HFB_Lipkin_Nogami_Lambda_2_import(Params)
            write(u,'(/,/,a)') 'Values of LN-HFB coefficients:'
            write(u,'(/,a)') 'pLambda_2 =  ' // cfmt_f(Lambda_2%p,15,8) // ' MeV'
            write(u,'(a)') 'nLambda_2 =  ' // cfmt_f(Lambda_2%n,15,8) // ' MeV'
        end if
        write(u,'(/,a)') 'Dispersion of HFB particle numbers dZ & dN:'
        write(u,'(/,a)') 'dZ = ' // cfmt_f(dA%p,12,6)
        write(u,'(a)') 'dN = ' // cfmt_f(dA%n,12,6)
        write(u,'(/,a)') 'dZ / Z = ' // cfmt_f(dA%p / real(Params%Calc%Z,dp),12,6)
        write(u,'(a)') 'dN / N = ' // cfmt_f(dA%n / real(Params%Calc%A - Params%Calc%Z,dp),12,6)
        close(u)
    end subroutine HFB_summary

    subroutine HFB_summary_SQS(Params,SQE,Rho,Orb)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SQE
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable :: nr(:,:)
        integer :: N_max, a_max, a, u, iso
        character(len=*), parameter :: Line = '_______________________________________________________________________________________________________'

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(nr(0:N_max,0:max(Params%Calc%N2max+1,2*N_max+1)))
        write(*,'(/,a)') 'Exporting information on single-quasiparticle orbits ...'
        open(newunit=u, file=io_path(Params,'/HFB/HFB_summary.dat'), position='append', action='write')
        do iso = 1, 2
            nr = 0
            if (iso == 1) then
                write(u,'(/,a)') 'Proton single-quasiparticle states'
            else
                write(u,'(/,a)') 'Neutron single-quasiparticle states'
            end if
            write(u,'(a)') Line
            write(u,'(a)') cfmt_s('n',4) // ' ' // cfmt_s('l',4) // ' ' // cfmt_s('2j',4) // ' ' // cfmt_s('O_a',10) // ' ' // cfmt_s('E',14)
            do a = 1, a_max
                if (iso == 1) then
                    write(u,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // ' ' // &
                                   cfmt_f(Rho%p(a,a),10,3) // ' ' // cfmt_f(SQE%p(a),14,5)
                else
                    write(u,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // ' ' // &
                                   cfmt_f(Rho%n(a,a),10,3) // ' ' // cfmt_f(SQE%n(a),14,5)
                end if
                nr(Orb(a)%l,Orb(a)%j) = nr(Orb(a)%l,Orb(a)%j) + 1
            end do
            write(u,'(a)') Line
        end do
        close(u)
    end subroutine HFB_summary_SQS

end module nj_hfb
