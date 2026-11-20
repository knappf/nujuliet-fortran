! Spherical HFB solver ... counterpart of HFB/*.jl ...
module nj_hfb
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
    use nj_obdm, only : OBDM_export
    use nj_hf, only : HF_orbital_ordering
    use nj_hf_driver, only : print_run_header
    use nj_qpo2b
    use nj_qp_common, only : pairing_energy
    implicit none
    private

    public :: HFB, HFB_density_operator, HFB_particle_number

contains

    subroutine HFB(Params)
        type(params_t), intent(in) :: Params
        integer(i8) :: t0

        call wigner_init_float(Params%Int%N2max + 2, 'Jmax', 6)
        write(*,'(a)') 'Starting spherical HFB calculation with NO2B NN+NNN interaction ...'
        call print_run_header(Params)

        if (is_dir(io_path(Params,''))) call remove_dir(io_path(Params,''))
        call make_dir(io_path(Params,''))
        call make_dir(io_path(Params,'/Bin'))
        call make_dir(io_path(Params,'/HFB'))
        call make_dir(io_path(Params,'/HFB/Densities'))

        write(*,'(/,a)') 'Starting spherical HFB calculation of designated nuclid ...'
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

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
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

        call OBDM_export(Params,Orb,io_path(Params,'/HFB/HFB_Summary.dat'), &
                         io_path(Params,'/HFB/Densities/HFB_Radial_Densities.dat'),Rho,C)
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
            call O1b_export(Params,Orb,C,io_path(Params,'/Bin/C.bin'))
            call O1b_export(Params,Orb,U,io_path(Params,'/Bin/U.bin'))
            call O1b_export(Params,Orb,V,io_path(Params,'/Bin/V.bin'))
            call O1b_export(Params,Orb,Rho_c,io_path(Params,'/Bin/Rho.bin'))
            call O1b_export(Params,Orb,Kappa_c,io_path(Params,'/Bin/Kappa.bin'))
            call qpO1B_export(Params,Orb,H_N,io_path(Params,'/Bin/qpH1B.bin'))
            call qpO2b_export(Params,Orb_NN,H_NN,io_path(Params,'/Bin/qpH2B.bin'))

            call qp_BMBPT_energy(Params,Orb,Orb_NN,H_N,H_NN,'HFB',io_path(Params,'/HFB/HFB_Summary.dat'))
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
        type(pnfloat_t) :: Lambda_1, Lambda_2
        real(dp) :: Z_target, N_target, epsilon, dE, d2E, dZ, dN, Z, N, Z_1, Z_2, N_1, N_2, E_HFB_1(4), E_HFB_2(4), Eye_shift
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
        SQE_old = SQE

        Z_1 = 0.0_dp
        Z_2 = 0.0_dp
        N_1 = 0.0_dp
        N_2 = 0.0_dp
        Lambda_1 = pnfloat_t(0.75_dp * Lambda%p, 0.75_dp * Lambda%n)
        Lambda_2 = Lambda

        write(*,'(/,a)') 'Starting iteration of HFB equations ...'
        do while (((dE > epsilon) .or. (dZ > epsilon) .or. (dN > epsilon)) .and. (Iteration < Iteration_max))
            ! Secant iterations for the chemical potentials ...
            do L = 1, Iteration_max
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
                    Lambda = HFB_Lambda_secant(Params,Lambda_1,pnfloat_t(Z_target - Z_1, N_target - N_1), &
                                               Lambda_2,pnfloat_t(Z_target - Z_2, N_target - N_2))
                    Lambda_1 = Lambda_2
                    Lambda_2 = Lambda
                    SQE = SQE_2
                    U = U_2
                    V = V_2
                    H = H_2
                    Delta = Delta_2
                end if

                ! Shift of H for the updated chemical potential ...
                do a = 1, a_max
                    Eye_shift = Lambda_1%p - Lambda%p
                    H%p(a,a) = H%p(a,a) + Eye_shift
                    Eye_shift = Lambda_1%n - Lambda%n
                    H%n(a,a) = H%n(a,a) + Eye_shift
                end do

                call diagonalize(Params,H,Delta,Orb,SQE,U,V)
                call HFB_density_operator(Params,U,V,Orb,Rho_Temp,Kappa_Temp)
                call HFB_particle_number(Params,Rho_Temp,Orb,Z,N)
                Z_1 = Z_2
                N_1 = N_2
                Z_2 = Z
                N_2 = N
                Lambda = HFB_Lambda_secant(Params,Lambda_1,pnfloat_t(Z_target - Z_1, N_target - N_1), &
                                           Lambda_2,pnfloat_t(Z_target - Z_2, N_target - N_2))
                Lambda_1 = Lambda_2
                Lambda_2 = Lambda
                dZ = abs(Z_target - Z)
                dN = abs(N_target - N)
                if (dZ < epsilon .and. dN < epsilon) exit
            end do

            ! New densities ... the Broyden update of the original reduces to a block symmetrization (see notes) ...
            call HFB_density_operator(Params,U,V,Orb,Rho,Kappa)
            call HFB_Broyden_update(Orb,Rho,Kappa)

            call HFB_allocate(Params,Lambda,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,H,Delta)
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
        real(dp), allocatable :: pDelta(:,:), nDelta(:,:)
        real(dp) :: Ph, pE_a, pD_a, nE_a, nD_a
        integer :: a_max, a

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(U,a_max)
        call o1b_alloc(V,a_max)
        allocate(SQE%p(a_max), SQE%n(a_max), SPE%p(a_max), SPE%n(a_max), Cc%p(a_max,a_max), Cc%n(a_max,a_max))
        call eig_sym(H%p,SPE%p,Cc%p)
        call eig_sym(H%n,SPE%n,Cc%n)
        call HF_orbital_ordering(Orb,a_max,Cc,SPE)
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

    ! Assigns quasiparticle states to orbital slots with matching <l> & <2j> ... optional final sort by |U|^2 ...
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
            real(dp), allocatable :: Nrm(:), lv(:), jv(:)
            integer, allocatable :: Perm(:), Ord(:)
            logical, allocatable :: Mask(:)
            integer :: n, a, b
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
            allocate(lv(n), jv(n), Ord(n), Mask(n))
            Ord = 0
            Mask = .false.
            do a = 1, n
                lv(a) = sum(real(Orb(1:n)%l,dp) * (Um(:,a)**2 + Vm(:,a)**2))
                jv(a) = sum(real(Orb(1:n)%j,dp) * (Um(:,a)**2 + Vm(:,a)**2))
            end do
            do a = 1, n
                do b = 1, n
                    if (.not. Mask(b) .and. abs(real(Orb(b)%l,dp) - lv(a)) < 1.0e-3_dp .and. &
                        abs(real(Orb(b)%j,dp) - jv(a)) < 1.0e-3_dp) then
                        Ord(b) = a
                        Mask(b) = .true.
                        exit
                    end if
                end do
            end do
            if (any(Ord == 0)) call fatal('HFB quasiparticle ordering failed ... mixed l & j detected ...')
            E = E(Ord)
            Um = Um(:,Ord)
            Vm = Vm(:,Ord)
        end subroutine order_one
    end subroutine HFB_orbital_ordering

    ! Reorders canonical orbitals (columns of C) & occupations to the orbital slots with matching l & 2j ...
    subroutine HFB_canonical_orbital_ordering(Params,Occ,C,Orb)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(inout) :: Occ
        type(o1b_t), intent(inout) :: C
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable :: Ord(:)
        Ord = canonical_order(C%p,Orb)
        Occ%p = Occ%p(Ord)
        C%p = C%p(:,Ord)
        Ord = canonical_order(C%n,Orb)
        Occ%n = Occ%n(Ord)
        C%n = C%n(:,Ord)
    end subroutine HFB_canonical_orbital_ordering

    function canonical_order(Cm,Orb) result(Ord)
        real(dp), intent(in) :: Cm(:,:)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable :: Ord(:)
        logical, allocatable :: Mask(:)
        real(dp) :: lv, jv
        integer :: n, a, b
        n = size(Cm,1)
        allocate(Ord(n), Mask(n))
        Ord = 0
        Mask = .false.
        do a = 1, n
            lv = sum(real(Orb(1:n)%l,dp) * Cm(:,a)**2)
            jv = sum(real(Orb(1:n)%j,dp) * Cm(:,a)**2)
            do b = 1, n
                if (.not. Mask(b) .and. abs(real(Orb(b)%l,dp) - lv) < 1.0e-3_dp .and. abs(real(Orb(b)%j,dp) - jv) < 1.0e-3_dp) then
                    Ord(b) = a
                    Mask(b) = .true.
                    exit
                end if
            end do
        end do
        if (any(Ord == 0)) call fatal('HFB canonical orbital ordering failed ... mixed l & j detected ...')
    end function canonical_order

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
        C%p = C%p(:,canonical_order(C%p,Orb))
        C%n = C%n(:,canonical_order(C%n,Orb))
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
        C%p = C%p(:,canonical_order(C%p,Orb))
        C%n = C%n(:,canonical_order(C%n,Orb))
    end function HFB_canonical_basis_BCS

    ! The Broyden update of the original sets the trial densities equal to the new densities (the mixing terms are
    ! commented out there), so its net effect is the block symmetrization & noise cleaning of HFB_Broyden_reconstruct ...
    subroutine HFB_Broyden_update(Orb,Rho,Kappa)
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(inout) :: Rho, Kappa
        call sym_block(Rho%p)
        call sym_block(Rho%n)
        call sym_block(Kappa%p)
        call sym_block(Kappa%n)
    contains
        subroutine sym_block(M)
            real(dp), intent(inout) :: M(:,:)
            real(dp), allocatable :: R(:,:)
            integer :: n, a, b
            n = size(M,1)
            allocate(R(n,n))
            R = 0.0_dp
            do a = 1, n
                do b = 1, n
                    if (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j) R(a,b) = 0.5_dp * (M(a,b) + M(b,a))
                end do
            end do
            where (abs(R) < 1.0e-14_dp) R = 0.0_dp
            M = R
        end subroutine sym_block
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

    function HFB_Lambda_secant(Params,Lambda_1,dA_1,Lambda_2,dA_2) result(L)
        type(params_t), intent(in) :: Params
        type(pnfloat_t), intent(in) :: Lambda_1, dA_1, Lambda_2, dA_2
        type(pnfloat_t) :: L
        real(dp), parameter :: Denom_eps = 1.0e-7_dp
        real(dp) :: dLmax, dp_L, dn_L
        dLmax = Params%Calc%HFB%dLmax
        L = Lambda_2
        if (abs(dA_2%p - dA_1%p) > Denom_eps) then
            dp_L = (Lambda_1%p * dA_2%p - Lambda_2%p * dA_1%p) / (dA_2%p - dA_1%p) - L%p
        else if (abs(dA_2%p) <= huge(1.0_dp)) then
            dp_L = 0.1_dp * dA_2%p
        else
            dp_L = 0.05_dp * L%p
        end if
        if (abs(dA_2%n - dA_1%n) > Denom_eps) then
            dn_L = (Lambda_1%n * dA_2%n - Lambda_2%n * dA_1%n) / (dA_2%n - dA_1%n) - L%n
        else if (abs(dA_2%n) <= huge(1.0_dp)) then
            dn_L = 0.1_dp * dA_2%n
        else
            dn_L = 0.05_dp * L%n
        end if
        if (abs(dp_L) > dLmax) dp_L = julia_sign(dp_L) * dLmax
        if (abs(dn_L) > dLmax) dn_L = julia_sign(dn_L) * dLmax
        L%p = L%p + dp_L
        L%n = L%n + dn_L
    contains
        pure real(dp) function julia_sign(x)
            real(dp), intent(in) :: x
            julia_sign = 0.0_dp
            if (x > 0.0_dp) julia_sign = 1.0_dp
            if (x < 0.0_dp) julia_sign = -1.0_dp
        end function julia_sign
    end function HFB_Lambda_secant

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
        open(newunit=u, file=io_path(Params,'/Bin/Lambda_2.bin'), access='stream', form='unformatted', status='replace', action='write')
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
        open(newunit=u, file=io_path(Params,'/Bin/Lambda_2.bin'), access='stream', form='unformatted', status='old', action='read')
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
        real(dp), allocatable :: pRb(:,:), nRb(:,:), pKb(:,:), nKb(:,:), pM2(:,:), nM2(:,:)
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
        !$omp parallel do schedule(dynamic) reduction(+:pTr1,pTr2,nTr1,nTr2) private(d,n_a,l_a,j_a,n_d,l_d,j_d,j_a_hat, &
        !$omp pRho_ad,nRho_ad,b,e,n_b,l_b,j_b,n_e,l_e,j_e,P_ab,pRho_be,nRho_be,J,J_hat,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f, &
        !$omp pRho_cf,nRho_cf,ME111,ME113,pKappa_ad,nKappa_ad,Amp,pKappa_be,nKappa_be,P_adc)
        do a = 1, a_max
            d = a
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
        end do
        !$omp end parallel do

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
        open(newunit=u, file=io_path(Params,'/HFB/HFB_Summary.dat'), position='append', action='write')
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
        open(newunit=u, file=io_path(Params,'/HFB/HFB_Summary.dat'), position='append', action='write')
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
