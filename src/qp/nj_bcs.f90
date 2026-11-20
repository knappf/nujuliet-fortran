! Spherical HF-BCS solver ... counterpart of BCS/*.jl ...
module nj_bcs
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
    use nj_hf, only : HF_solve, HF_energy, HF_orbital_ordering
    use nj_hf_driver, only : print_run_header
    use nj_qpo2b
    use nj_qp_common, only : pairing_energy
    implicit none
    private

    public :: BCS, BCS_reference

contains

    subroutine BCS(Params)
        type(params_t), intent(in) :: Params
        type(params_t) :: Params_Ref
        integer(i8) :: t0

        call wigner_init_float(Params%Int%N2max + 2, 'Jmax', 6)
        write(*,'(a)') 'Starting BCS calculation with NO2B NN+NNN interaction ...'
        call print_run_header(Params)

        if (is_dir(io_path(Params,''))) call remove_dir(io_path(Params,''))
        call make_dir(io_path(Params,''))
        call make_dir(io_path(Params,'/Bin'))
        call make_dir(io_path(Params,'/BCS'))
        call make_dir(io_path(Params,'/BCS/Densities'))

        Params_Ref = BCS_reference(Params)
        write(*,'(/,a)') 'Reference closed-shell nucleus selected for BCS calculation is:'
        write(*,'(a)') '     A = ' // istr(Params_Ref%Calc%A) // ',     Z = ' // istr(Params_Ref%Calc%Z)

        write(*,'(/,a)') 'Starting BCS calculation ...'
        call timer_start(t0)
        call BCS_solver(Params,Params_Ref)
        call timer_report(t0)
        write(*,'(/,a,/)') 'All calculations have finished ...'
    end subroutine BCS

    ! Reference closed-(sub)shell nucleus ... user given (A0, Z0) or the closest lower magic numbers ...
    function BCS_reference(Params) result(Params_Ref)
        type(params_t), intent(in) :: Params
        type(params_t) :: Params_Ref
        integer, parameter :: pMagic(6) = [2,8,20,28,50,82], nMagic(7) = [2,8,20,28,50,82,126]
        integer :: N_target, Z_target, N_ref, Z_ref, dN, dZ, k

        Params_Ref = Params
        ! QTDA & QRPA parameters are not forwarded by the original ... defaults ...
        Params_Ref%Calc%QTDA = qtda_params_t()
        Params_Ref%Calc%QRPA = qrpa_params_t()
        if (Params%Calc%BCS%A0 > 0 .and. Params%Calc%BCS%Z0 > 0) then
            Params_Ref%Calc%A = Params%Calc%BCS%A0
            Params_Ref%Calc%Z = Params%Calc%BCS%Z0
        else
            N_target = Params%Calc%A - Params%Calc%Z
            Z_target = Params%Calc%Z
            N_ref = 0
            Z_ref = 0
            dN = 1000
            dZ = 1000
            do k = 1, size(pMagic)
                if (pMagic(k) <= Z_target .and. (Z_target - pMagic(k)) < dZ) then
                    Z_ref = pMagic(k)
                    dZ = Z_target - pMagic(k)
                end if
            end do
            do k = 1, size(nMagic)
                if (nMagic(k) <= N_target .and. (N_target - nMagic(k)) < dN) then
                    N_ref = nMagic(k)
                    dN = N_target - nMagic(k)
                end if
            end do
            Params_Ref%Calc%A = N_ref + Z_ref
            Params_Ref%Calc%Z = Z_ref
        end if
    end function BCS_reference

    subroutine BCS_solver(Params,Params_Ref)
        type(params_t), intent(in) :: Params, Params_Ref
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: T, C, Rho, Kappa, U_m, V_m
        type(o2b_t) :: V_NN, V_NN_res
        type(orb2b_t) :: Orb_NN
        type(v3b_t) :: V_NNN
        type(orb3b_t) :: Orb_NNN
        type(pnvector_t) :: SPE, SQE, U, V, Delta
        type(pnfloat_t) :: Lambda, dA
        type(qpo1b_t) :: H_N
        type(qpo2b_t) :: H_NN
        real(dp) :: E_MF, E_BCS(2), T_BCS
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
        call HF_BCS_solve(Params,Params_Ref,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Lambda,SPE,SQE,C,U,V,Rho,Kappa,Delta, &
                          Convergence,Iteration)
        call timer_report(t0)

        call BCS_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,E_MF,E_BCS)
        T_BCS = T1b_energy(Params,Rho,Orb,T)
        dA = BCS_particle_number_dispersion(Params,U,V,Orb)
        call BCS_summary(Params,Params_Ref,E_MF,E_BCS,T_BCS,Lambda,dA,Convergence,Iteration)

        call OBDM_export(Params_Ref,Orb,io_path(Params,'/BCS/BCS_Summary.dat'), &
                         io_path(Params,'/BCS/Densities/BCS_Radial_Densities.dat'),Rho,C)
        call BCS_summary_SQS(Params_Ref,SPE,SQE,U,V)

        if (Params%Calc%BCS%BMF) then
            H_N%qp11%p = diagm(SQE%p)
            H_N%qp11%n = diagm(SQE%n)
            call o1b_alloc(H_N%qp20,size(SQE%p))
            call timer_start(t0)
            call V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)
            call O2b_move(V_NN_res,V_NN)
            call timer_report(t0)
            U_m%p = diagm(U%p)
            U_m%n = diagm(U%n)
            V_m%p = diagm(V%p)
            V_m%n = diagm(V%n)
            call timer_start(t0)
            call qpO2b(Params,Orb,Orb_NN,V_NN,C,U_m,V_m,H_NN)
            call timer_report(t0)
            call BCS_export(Params,Orb,Orb_NN,C,U_m,V_m,H_N,H_NN)
            call qp_BMBPT_energy(Params,Orb,Orb_NN,H_N,H_NN,'HF-BCS',io_path(Params,'/BCS/BCS_Summary.dat'))
        end if
    end subroutine BCS_solver

    subroutine HF_BCS_solve(Params,Params_Ref,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Lambda,SPE,SQE,C,U,V,Rho,Kappa,Delta, &
                            Convergence,Iteration)
        type(params_t), intent(in) :: Params, Params_Ref
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: T
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(pnfloat_t), intent(out) :: Lambda
        type(pnvector_t), intent(out) :: SPE, SQE, U, V, Delta
        type(o1b_t), intent(out) :: C, Rho, Kappa
        logical, intent(out) :: Convergence
        integer, intent(out) :: Iteration
        type(orb1b_t), allocatable :: Orb_HF(:)
        type(o1b_t) :: h
        type(o2b_t) :: V_NN_Res
        type(pnvector_t) :: SPE_Old
        real(dp) :: Z_target, N_target, epsilon, Z, N, dZ, dN, dE, E_HF
        integer :: a_max, Iteration_max, Iteration_HF, Iteration_BCS
        integer(i8) :: t0

        Z_target = real(Params%Calc%Z,dp)
        N_target = real(Params%Calc%A - Params%Calc%Z,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        epsilon = Params%Calc%BCS%Tol
        Iteration = 0
        Iteration_max = Params%Calc%BCS%Imax
        Convergence = .false.

        allocate(Delta%p(a_max), Delta%n(a_max))
        Delta%p = Params%Calc%BCS%pD0
        Delta%n = Params%Calc%BCS%nD0

        write(*,'(/,a,/)') 'Starting iteration of HF-BCS with NO2B NN+NNN interaction ...'

        ! Reference HF calculation ... occupations of the reference closed-shell nucleus ...
        call timer_start(t0)
        Orb_HF = orbitals_make(Params_Ref,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        ! The original also reads the NN interaction once more here (V_NN_Ref) without ever using it ... skipped ...
        write(*,'(/,a)') 'Solving the HF equations for reference closed-shell system ...'
        write(*,'(a,/)') 'Reference nucleus:     A = ' // istr(Params_Ref%Calc%A) // ',     Z = ' // istr(Params_Ref%Calc%Z)
        call HF_solve(Params,Orb_HF,Orb_NN,Orb_NNN,T,V_NN,V_NNN,h,C,Rho,Iteration_HF)
        SPE%p = diag(h%p)
        SPE%n = diag(h%n)
        write(*,'(/,a)') 'Calculating the HF mean-field ground-state energy ... sanity check ...'
        E_HF = HF_energy(Params,Rho,Orb_HF,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
        write(*,'(/,a)') 'Making density-depenent residual NN interaction ... s-wave channel (J = 0) ...'
        call BCS_V2b_res(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,C,Rho,V_NN_Res)
        Lambda = HF_BCS_initialize_chemical_potential(Params_Ref,SPE,Orb_HF)
        call timer_report(t0)

        write(*,'(/,a)') 'Initializing the HF-BCS approximation ...'
        Z = Z_target
        N = N_target
        dZ = abs(Z_target - Z) + 0.1_dp
        dN = abs(N_target - N) + 0.1_dp

        write(*,'(/,a,/)') 'Starting iteration of BCS equations ...'
        do while (((abs(dZ) > epsilon) .or. (abs(dN) > epsilon)) .and. (Iteration < Iteration_max))
            if (Iteration == 0) then
                SQE = BCS_allocate_SQE(Params,SPE,Lambda,Delta)
                call BCS_allocate_amplitudes(Params,Lambda,SPE,Delta,Orb,U,V)
                call BCS_particle_number(Params,V,Orb,Z,N)
                dZ = Z_target - Z
                dN = N_target - N
                Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnfloat_t(dZ,dN),Orb)
            else
                Delta = BCS_allocate_Delta(Params,U,V,Orb,Orb_NN,V_NN_Res)
                call BCS_allocate_amplitudes(Params,Lambda,SPE,Delta,Orb,U,V)
                call BCS_particle_number(Params,V,Orb,Z,N)
                Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnfloat_t(dZ,dN),Orb)
            end if

            Iteration = Iteration + 1
            dZ = Z_target - Z
            dN = N_target - N
            write(*,'(a)') TAB // 'BCS iteration number:   ' // istr(Iteration) // '   dZ = ' // jstr(round_sig(dZ,8)) // &
                           '   &   dN = ' // jstr(round_sig(dN,8))
            write(*,'(a)') TAB // TAB // 'Current particle numbers     ...     Z = ' // jstr(Z) // ', N = ' // jstr(N)
            write(*,'(a)') TAB // TAB // 'Current chemical potentials        ...     pLambda = ' // jstr(round_digits(Lambda%p,6)) // &
                           ' MeV, nLambda = ' // jstr(round_digits(Lambda%n,6)) // ' MeV'

            if ((abs(dZ) < epsilon) .and. (abs(dN) < epsilon)) then
                if (.not. Params%Calc%BCS%ScBCS) then
                    Convergence = .true.
                    write(*,'(/,a)') 'BCS iteration with residual NN interaction has converged ...'
                else
                    write(*,'(/,a)') 'Initial ScBCS iteration with residual NN interaction has converged ...'
                end if
                exit
            end if
            if (Iteration == Iteration_max) then
                Convergence = .false.
                write(*,'(/,a)') 'BCS iteration with residual NN interaction has terminated ...'
                write(*,'(/,a,/)') TAB // ' (!!!) CONVERGENCE WAS NOT REACHED (!!!) '
                exit
            end if
        end do

        ! Self-consistent iteration of the BCS equations ...
        if (Params%Calc%BCS%ScBCS) then
            write(*,'(/,a,/)') 'Starting self-consistent iteration of BCS equations ...'
            Iteration = 0
            dE = 1.0_dp
            allocate(SPE_Old%p(a_max), SPE_Old%n(a_max))
            SPE_Old%p = 0.0_dp
            SPE_Old%n = 0.0_dp
            do while ((dE > epsilon) .and. (Iteration < Iteration_max))
                call BCS_density_operator(Params,U,V,Rho,Kappa)
                Rho%p = mmt(mm(C%p,Rho%p),C%p)
                Rho%n = mmt(mm(C%n,Rho%n),C%n)
                Kappa%p = mmt(mm(C%p,Kappa%p),C%p)
                Kappa%n = mmt(mm(C%n,Kappa%n),C%n)

                h = HF_BCS_allocate(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
                call eig_sym(h%p,SPE%p,C%p)
                call eig_sym(h%n,SPE%n,C%n)
                call HF_orbital_ordering(Orb,a_max,C,SPE)

                ! The original silences the terminal output of this call ...
                call BCS_V2b_res(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,C,Rho,V_NN_Res)

                Iteration_BCS = 0
                dZ = 0.1_dp
                dN = 0.1_dp
                write(*,'(a,/)') TAB // 'Starting iteration of BCS equations ...'
                do while (((dZ > epsilon) .or. (dN > epsilon)) .and. (Iteration_BCS < Iteration_max))
                    if (Iteration_BCS == 0) then
                        call BCS_allocate_amplitudes(Params,Lambda,SPE,Delta,Orb,U,V)
                        call BCS_particle_number(Params,V,Orb,Z,N)
                        dZ = abs(Z_target - Z)
                        dN = abs(N_target - N)
                        Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnfloat_t(dZ,dN),Orb)
                    else
                        Delta = BCS_allocate_Delta(Params,U,V,Orb,Orb_NN,V_NN_Res)
                        call BCS_allocate_amplitudes(Params,Lambda,SPE,Delta,Orb,U,V)
                        call BCS_particle_number(Params,V,Orb,Z,N)
                        Lambda = BCS_Lambda(Params,SPE,Delta,Lambda,pnfloat_t(dZ,dN),Orb)
                    end if
                    Iteration_BCS = Iteration_BCS + 1
                    dZ = abs(real(Params%Calc%Z,dp) - Z)
                    dN = abs(real(Params%Calc%A - Params%Calc%Z,dp) - N)
                end do

                Iteration = Iteration + 1
                dE = (sum(abs(SPE%p - SPE_Old%p)) + sum(abs(SPE%n - SPE_Old%n))) / real(2 * a_max,dp)
                write(*,'(a)') TAB // 'ScBCS teration number:   ' // istr(Iteration) // '   dE = ' // jstr(round_sig(dE,8)) // ' MeV' // &
                               ',   dZ = ' // jstr(round_sig(dZ,8)) // '   &   dN = ' // jstr(round_sig(dN,8))
                write(*,'(a)') TAB // TAB // 'Current particle numbers     ...     Z = ' // jstr(Z) // ', N = ' // jstr(N)
                write(*,'(a)') TAB // TAB // 'Current chemical potentials        ...     pLambda = ' // jstr(round_digits(Lambda%p,6)) // &
                               ' MeV, nLambda = ' // jstr(round_digits(Lambda%n,6)) // ' MeV'
                SPE_Old%p = SPE%p
                SPE_Old%n = SPE%n
                if ((abs(dZ) < epsilon) .and. (abs(dN) < epsilon) .and. dE < epsilon) then
                    Convergence = .true.
                    write(*,'(/,a)') 'ScBCS iteration with residual NN interaction has converged ...'
                    exit
                end if
                if (Iteration == Iteration_max) then
                    Convergence = .false.
                    write(*,'(/,a)') 'ScBCS iteration with residual NN interaction has terminated ...'
                    write(*,'(/,a,/)') TAB // ' (!!!) CONVERGENCE WAS NOT REACHED (!!!) '
                    exit
                end if
            end do
        end if

        SQE = BCS_allocate_SQE(Params,SPE,Lambda,Delta)
        call BCS_allocate_amplitudes(Params,Lambda,SPE,Delta,Orb,U,V)
        call BCS_density_operator(Params,U,V,Rho,Kappa)
        Rho%p = mmt(mm(C%p,Rho%p),C%p)
        Rho%n = mmt(mm(C%n,Rho%n),C%n)
        Kappa%p = mmt(mm(C%p,Kappa%p),C%p)
        Kappa%n = mmt(mm(C%n,Kappa%n),C%n)
    end subroutine HF_BCS_solve


    ! J = 0 residual NN interaction with the NO2B NNN contraction, transformed to the HF basis C ...
    subroutine BCS_V2b_res(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,C,Rho,V_NN_Res)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t), intent(in) :: C, Rho
        type(o2b_t), intent(out) :: V_NN_Res
        integer :: a_max, N_3max, J, P, N_T0, N_T1, Bra, Ket, Ind, a, b, d, e, cc, f, P3B
        integer :: n_a, l_a, n_b, l_b, n_d, l_d, n_e, l_e, n_c, l_c, j_c, n_f, l_f, j_f
        real(dp) :: p2, p3, Hat, pnSum, ppSum, nnSum, ME001, ME101, ME011, ME111, ME113
        integer(i8) :: t0

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        N_3max = Params%Calc%N3max
        p2 = Params%Int%cP2N
        p3 = Params%Int%cP3N
        V_NN_Res = V_NN

        write(*,'(/,a,/)') 'Starting calculation of residual NN interaction...'
        write(*,'(/,a)') 'Making density dependent residual 2-body interaction...'
        call timer_start(t0)
        J = 0
        do P = 1, 2
            if (P == 1) then
                write(*,'(a)') 'Calculating   ...   J = 0/' // istr(Params%Calc%N2max+1) // TAB // 'P = +'
            else
                write(*,'(a)') 'Calculating   ...   J = 0/' // istr(Params%Calc%N2max+1) // TAB // 'P = -'
            end if
            N_T0 = Orb_NN%N(1,P,J)
            N_T1 = Orb_NN%N(2,P,J)
            Hat = 1.0_dp / real(2*J + 1,dp)
            !$omp parallel do schedule(dynamic) private(a,b,n_a,l_a,n_b,l_b,Ket,Ind,d,e,n_d,l_d,n_e,l_e,pnSum,ppSum,nnSum, &
            !$omp cc,n_c,l_c,j_c,P3B,f,n_f,l_f,j_f,ME001,ME101,ME011,ME111,ME113)
            do Bra = 1, max(N_T0,N_T1)
                if (Bra <= N_T0) then
                    a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                    b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                    l_a = Orb(a)%l
                    n_a = Orb(a)%n
                    l_b = Orb(b)%l
                    n_b = Orb(b)%n
                    do Ket = 1, Bra
                        Ind = packed_index(Bra,Ket,N_T0)
                        d = Orb_NN%Ind(1,P,J)%v(1,Ket)
                        e = Orb_NN%Ind(1,P,J)%v(2,Ket)
                        l_d = Orb(d)%l
                        n_d = Orb(d)%n
                        l_e = Orb(e)%l
                        n_e = Orb(e)%n
                        pnSum = 0.0_dp
                        do cc = 1, a_max
                            n_c = Orb(cc)%n
                            l_c = Orb(cc)%l
                            j_c = Orb(cc)%j
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                            P3B = mod(l_a + l_b + l_c, 2) + 1
                            do f = 1, a_max
                                n_f = Orb(f)%n
                                l_f = Orb(f)%l
                                j_f = Orb(f)%j
                                if ((l_c == l_f) .and. (j_c == j_f) .and. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)) then
                                    ME001 = V3b_no2b(a,b,cc,0,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME101 = V3b_no2b(a,b,cc,1,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME011 = V3b_no2b(a,b,cc,0,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME111 = V3b_no2b(a,b,cc,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME113 = V3b_no2b(a,b,cc,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)
                                    pnSum = pnSum + Hat * ((0.5_dp * ME001 - 1.0_dp / sqrt(12.0_dp) * ME101 - 1.0_dp / sqrt(12.0_dp) * ME011 + &
                                            1.0_dp / 6.0_dp * ME111 + 1.0_dp / 3.0_dp * ME113) * Rho%p(cc,f) + (0.5_dp * ME001 + &
                                            1.0_dp / sqrt(12.0_dp) * ME101 + 1.0_dp / sqrt(12.0_dp) * ME011 + 1.0_dp / 6.0_dp * ME111 + &
                                            1.0_dp / 3.0_dp * ME113) * Rho%n(cc,f))
                                end if
                            end do
                        end do
                        V_NN_Res%pn(P,J)%v(Ind) = p2 * V_NN_Res%pn(P,J)%v(Ind) + p3 * pnSum
                    end do
                end if

                if (Bra <= N_T1) then
                    a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                    b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                    l_a = Orb(a)%l
                    n_a = Orb(a)%n
                    l_b = Orb(b)%l
                    n_b = Orb(b)%n
                    do Ket = 1, Bra
                        Ind = packed_index(Bra,Ket,N_T1)
                        d = Orb_NN%Ind(2,P,J)%v(1,Ket)
                        e = Orb_NN%Ind(2,P,J)%v(2,Ket)
                        l_d = Orb(d)%l
                        n_d = Orb(d)%n
                        l_e = Orb(e)%l
                        n_e = Orb(e)%n
                        ppSum = 0.0_dp
                        nnSum = 0.0_dp
                        do cc = 1, a_max
                            n_c = Orb(cc)%n
                            l_c = Orb(cc)%l
                            j_c = Orb(cc)%j
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                            P3B = mod(l_a + l_b + l_c, 2) + 1
                            do f = 1, a_max
                                n_f = Orb(f)%n
                                l_f = Orb(f)%l
                                j_f = Orb(f)%j
                                if ((l_c == l_f) .and. (j_c == j_f) .and. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)) then
                                    ME111 = V3b_no2b(a,b,cc,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                    ME113 = V3b_no2b(a,b,cc,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)
                                    ppSum = ppSum + Hat * (ME113 * Rho%p(cc,f) + (2.0_dp/3.0_dp * ME111 + 1.0_dp/3.0_dp * ME113) * Rho%n(cc,f))
                                    nnSum = nnSum + Hat * (ME113 * Rho%n(cc,f) + (2.0_dp/3.0_dp * ME111 + 1.0_dp/3.0_dp * ME113) * Rho%p(cc,f))
                                end if
                            end do
                        end do
                        V_NN_Res%pp(P,J)%v(Ind) = p2 * V_NN_Res%pp(P,J)%v(Ind) + p3 * ppSum
                        V_NN_Res%nn(P,J)%v(Ind) = p2 * V_NN_Res%nn(P,J)%v(Ind) + p3 * nnSum
                    end do
                end if
            end do
            !$omp end parallel do
        end do
        call timer_report(t0)

        write(*,'(/,a)') 'Transforming residual 2-body interaction from the LHO to the target basis...'
        call O2b_transformation(Params,Orb,Orb_NN,V_NN_Res,C)
        write(*,'(/,a,/)') 'Residual 2-body interaction ready...'
    end subroutine BCS_V2b_res

    subroutine BCS_density_operator(Params,U,V,Rho,Kappa)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: U, V
        type(o1b_t), intent(out) :: Rho, Kappa
        integer :: a_max, a
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(Rho,a_max)
        call o1b_alloc(Kappa,a_max)
        do a = 1, a_max
            Rho%p(a,a) = V%p(a)**2
            Kappa%p(a,a) = V%p(a) * U%p(a)
            Rho%n(a,a) = V%n(a)**2
            Kappa%n(a,a) = V%n(a) * U%n(a)
        end do
    end subroutine BCS_density_operator

    ! Chemical potentials for fixed Delta ... simple fixed-step iteration of the original ...
    function BCS_Lambda(Params,SPE,Delta,Lambda,dA,Orb) result(L)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SPE, Delta
        type(pnfloat_t), intent(in) :: Lambda, dA
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnfloat_t) :: L
        real(dp) :: Z_target, N_target, epsilon, dZ, dN, Z, N, V_a
        integer :: a_max, i, a

        Z_target = real(Params%Calc%Z,dp)
        N_target = real(Params%Calc%A - Params%Calc%Z,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        epsilon = Params%Calc%BCS%Tol
        L = Lambda

        if (abs(dA%p) > epsilon) then
            dZ = dA%p
            do i = 1, 500
                L%p = L%p + Params%Calc%BCS%q * dZ
                Z = 0.0_dp
                do a = 1, a_max
                    V_a = sqrt(0.5_dp * (1.0_dp - (SPE%p(a) - L%p) / sqrt((SPE%p(a) - L%p)**2 + Delta%p(a)**2)))
                    Z = Z + real(Orb(a)%j + 1,dp) * V_a**2
                end do
                dZ = Z_target - Z
                if (abs(dZ) < epsilon) exit
            end do
        end if

        if (abs(dA%n) > epsilon) then
            dN = dA%n
            do i = 1, 500
                L%n = L%n + Params%Calc%BCS%q * dN
                N = 0.0_dp
                do a = 1, a_max
                    V_a = sqrt(0.5_dp * (1.0_dp - (SPE%n(a) - L%n) / sqrt((SPE%n(a) - L%n)**2 + Delta%n(a)**2)))
                    N = N + real(Orb(a)%j + 1,dp) * V_a**2
                end do
                dN = N_target - N
                if (abs(dN) < epsilon) exit
            end do
        end if
    end function BCS_Lambda

    function HF_BCS_initialize_chemical_potential(Params,SPE,Orb) result(Lambda)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SPE
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnfloat_t) :: Lambda
        real(dp) :: pEps_particle, pEps_hole, nEps_particle, nEps_hole
        integer :: a_max, a

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        pEps_particle = 1000.0_dp
        pEps_hole = -1000.0_dp
        nEps_particle = 1000.0_dp
        nEps_hole = -1000.0_dp
        do a = 1, a_max
            if (Orb(a)%pO == 1) then
                if (pEps_hole < SPE%p(a)) pEps_hole = SPE%p(a)
            else if (Orb(a)%pO == 0) then
                if (pEps_particle > SPE%p(a)) pEps_particle = SPE%p(a)
            end if
            if (Orb(a)%nO == 1) then
                if (nEps_hole < SPE%n(a)) nEps_hole = SPE%n(a)
            else if (Orb(a)%nO == 0) then
                if (nEps_particle > SPE%n(a)) nEps_particle = SPE%n(a)
            end if
        end do
        Lambda = pnfloat_t(0.5_dp * (pEps_hole - pEps_particle), 0.5_dp * (nEps_hole - nEps_particle))
    end function HF_BCS_initialize_chemical_potential

    subroutine BCS_particle_number(Params,V,Orb,Z,N)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: V
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(out) :: Z, N
        integer :: a
        Z = 0.0_dp
        N = 0.0_dp
        do a = 1, (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
            Z = Z + real(Orb(a)%j + 1,dp) * V%p(a)**2
            N = N + real(Orb(a)%j + 1,dp) * V%n(a)**2
        end do
    end subroutine BCS_particle_number

    function BCS_particle_number_dispersion(Params,U,V,Orb) result(dA)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: U, V
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnfloat_t) :: dA
        real(dp) :: dZ, dN
        integer :: a
        write(*,'(/,a)') 'Calculating the BCS dispersion of proton & neutron particle numbers ...'
        dZ = 0.0_dp
        dN = 0.0_dp
        do a = 1, (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
            dZ = dZ + 2.0_dp * real(Orb(a)%j + 1,dp) * U%p(a)**2 * V%p(a)**2
            dN = dN + 2.0_dp * real(Orb(a)%j + 1,dp) * U%n(a)**2 * V%n(a)**2
        end do
        dZ = sqrt(dZ)
        dN = sqrt(dN)
        write(*,'(/,a)') 'BCS dispersion of nucleon numbers are ...'
        write(*,'(a)') TAB // 'dZ = ' // jstr(round_digits(dZ,5))
        write(*,'(a)') TAB // 'dN = ' // jstr(round_digits(dN,5))
        dA = pnfloat_t(dZ,dN)
    end function BCS_particle_number_dispersion

    ! Self-consistent HF field including the anomalous NNN contribution ... HF_BCS_allocate ...
    function HF_BCS_allocate(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN) result(H)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, Kappa, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t) :: H
        integer :: a_max, N_2max, N_3max, ad_i, a, d, b, e, c, f, J, P_ab, P_abc, P_abe, Pass, ad_count
        integer :: n_a, l_a, j_a, n_d, l_d, j_d, n_b, l_b, j_b, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        integer, allocatable :: ad(:,:)
        real(dp) :: A_mass, p3, is_j_a_hat, J_j_a_hat, pHSum, nHSum, pRho_be, nRho_be, pRho_cf, nRho_cf
        real(dp) :: pKappa_be, nKappa_be, pKappa_cf, nKappa_cf, A_NNN_Amp, ME001, ME101, ME011, ME111, ME113
        character(len=32) :: CMS

        A_mass = real(Params%Calc%A,dp)
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        p3 = Params%Int%cP3N
        CMS = Params%Calc%CMS
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(H,a_max)

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

        !$omp parallel do schedule(dynamic) private(a,d,n_a,l_a,j_a,n_d,l_d,j_d,is_j_a_hat,pHSum,nHSum,b,e,n_b,l_b,j_b, &
        !$omp n_e,l_e,j_e,pRho_be,nRho_be,P_ab,J,J_j_a_hat,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f,pRho_cf,nRho_cf, &
        !$omp ME001,ME101,ME011,ME111,ME113,pKappa_be,nKappa_be,P_abe,A_NNN_Amp,pKappa_cf,nKappa_cf)
        do ad_i = 1, ad_count
            a = ad(1,ad_i)
            d = ad(2,ad_i)
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            n_d = Orb(d)%n
            l_d = Orb(d)%l
            j_d = Orb(d)%j
            is_j_a_hat = 1.0_dp / (real(j_a,dp) + 1.0_dp)
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
                    ! Normal density-dependent part ...
                    if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                        if (l_b == l_e .and. j_b == j_e .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max) then
                            pRho_be = Rho%p(b,e)
                            nRho_be = Rho%n(b,e)
                            P_ab = mod(l_a + l_b,2) + 1
                            do J = abs(j_d - j_e)/2, (j_d + j_e)/2
                                if (mod(l_a + l_b, 2) == mod(l_d + l_e, 2)) then
                                    J_j_a_hat = real(2*J + 1,dp) * is_j_a_hat
                                    pHSum = pHSum + J_j_a_hat * pRho_be * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                                    pHSum = pHSum + J_j_a_hat * nRho_be * O2b_pn(a,b,d,e,J,P_ab,V_NN,Orb_NN)
                                    nHSum = nHSum + J_j_a_hat * nRho_be * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                                    nHSum = nHSum + J_j_a_hat * pRho_be * O2b_pn(b,a,e,d,J,P_ab,V_NN,Orb_NN)
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
                    if (l_b == l_e .and. j_b == j_e .and. (2*(n_b + n_e) + l_b + l_e) <= N_2max .and. &
                        (2*(n_a + n_b + n_e) + l_a + l_b + l_e) <= N_3max) then
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
                                        1.0_dp / 3.0_dp * (2.0_dp * ME111 + ME113) * nKappa_be * nKappa_cf)
                                nHSum = nHSum + p3 * A_NNN_Amp * (ME113 * nKappa_be * nKappa_cf + &
                                        1.0_dp / 3.0_dp * (2.0_dp * ME111 + ME113) * pKappa_be * pKappa_cf)
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
            H%p(a,d) = pHSum
            H%p(d,a) = pHSum
            H%n(a,d) = nHSum
            H%n(d,a) = nHSum
        end do
        !$omp end parallel do
    end function HF_BCS_allocate

    function BCS_allocate_SQE(Params,SPE,Lambda,Delta) result(SQE)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SPE, Delta
        type(pnfloat_t), intent(in) :: Lambda
        type(pnvector_t) :: SQE
        write(*,'(/,a)') 'Evalutiang BCS single-quasiparticle energies ...'
        SQE%p = sqrt((SPE%p - Lambda%p)**2 + Delta%p**2)
        SQE%n = sqrt((SPE%n - Lambda%n)**2 + Delta%n**2)
    end function BCS_allocate_SQE

    subroutine BCS_allocate_amplitudes(Params,Lambda,SPE,Delta,Orb,U,V)
        type(params_t), intent(in) :: Params
        type(pnfloat_t), intent(in) :: Lambda
        type(pnvector_t), intent(in) :: SPE, Delta
        type(orb1b_t), intent(in) :: Orb(:)
        type(pnvector_t), intent(inout) :: U, V
        integer :: a_max, a
        real(dp) :: Amp, pME, nME

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        if (allocated(U%p)) deallocate(U%p, U%n)
        if (allocated(V%p)) deallocate(V%p, V%n)
        allocate(U%p(a_max), U%n(a_max), V%p(a_max), V%n(a_max))
        do a = 1, a_max
            Amp = phase(Orb(a)%l)
            pME = (SPE%p(a) - Lambda%p) / sqrt((SPE%p(a) - Lambda%p)**2 + Delta%p(a)**2)
            V%p(a) = sqrt(0.5_dp * (1.0_dp - pME))
            U%p(a) = Amp * sqrt(0.5_dp * (1.0_dp + pME))
            nME = (SPE%n(a) - Lambda%n) / sqrt((SPE%n(a) - Lambda%n)**2 + Delta%n(a)**2)
            V%n(a) = sqrt(0.5_dp * (1.0_dp - nME))
            U%n(a) = Amp * sqrt(0.5_dp * (1.0_dp + nME))
        end do
    end subroutine BCS_allocate_amplitudes

    function BCS_allocate_Delta(Params,U,V,Orb,Orb_NN_Res,V_NN_Res) result(Delta)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: U, V
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN_Res
        type(o2b_t), intent(in) :: V_NN_Res
        type(pnvector_t) :: Delta
        integer :: a_max, a, b
        real(dp) :: pSum, nSum, ja_jb_hat

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(Delta%p(a_max), Delta%n(a_max))
        do a = 1, a_max
            pSum = 0.0_dp
            nSum = 0.0_dp
            do b = 1, a_max
                ja_jb_hat = sqrt(real(Orb(b)%j + 1,dp) / real(Orb(a)%j + 1,dp))
                pSum = pSum - 0.5_dp * ja_jb_hat * O2b_pp(a,a,b,b,0,1,V_NN_Res,Orb,Orb_NN_Res) * U%p(b) * V%p(b)
                nSum = nSum - 0.5_dp * ja_jb_hat * O2b_nn(a,a,b,b,0,1,V_NN_Res,Orb,Orb_NN_Res) * U%n(b) * V%n(b)
            end do
            Delta%p(a) = pSum
            Delta%n(a) = nSum
        end do
    end function BCS_allocate_Delta

    subroutine BCS_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,E_MF,E_BCS)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, Kappa, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        real(dp), intent(out) :: E_MF, E_BCS(2)
        real(dp) :: Sum2N, Sum3N
        write(*,'(/,a)') 'Calculating the total mean-field + BCS pairing ground-state energy ...'
        call pairing_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,E_MF,Sum2N,Sum3N,.false.)
        E_BCS = [Sum2N, Sum3N]
        write(*,'(a)') TAB // 'MF energy      ...   E_MF  = ' // jstr(E_MF) // ' MeV'
        write(*,'(a)') TAB // 'BCS energy     ...   E_BCS = ' // jstr(sum(E_BCS)) // ' MeV'
        write(*,'(a)') TAB // 'BCS 2N energy ...   E_BCS_2N = ' // jstr(E_BCS(1)) // ' MeV'
        write(*,'(a)') TAB // 'BCS 3N energy ...   E_BCS_3N = ' // jstr(E_BCS(2)) // ' MeV'
    end subroutine BCS_energy

    subroutine BCS_summary(Params,Params_Ref,E_MF,E_BCS,T_BCS,Lambda,dA,Convergence,Iteration)
        type(params_t), intent(in) :: Params, Params_Ref
        real(dp), intent(in) :: E_MF, E_BCS(2), T_BCS
        type(pnfloat_t), intent(in) :: Lambda, dA
        logical, intent(in) :: Convergence
        integer, intent(in) :: Iteration
        integer :: u, N_max
        real(dp) :: A
        character(len=*), parameter :: M = TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' '
        character(len=*), parameter :: MN = TAB // ' MeV / Nucleon ' // TAB // TAB // ' ... ' // TAB // ' '

        A = real(Params%Calc%A,dp)
        N_max = Params%Calc%Nmax
        write(*,'(/,a)') 'Exporting BCS calculation summary ...'
        open(newunit=u, file=io_path(Params,'/BCS/BCS_Summary.dat'), position='append', action='write')
        write(u,'(a)') 'Target nuclid:       A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)
        write(u,'(a)') 'Reference nuclid:    A = ' // istr(Params_Ref%Calc%A) // ', Z = ' // istr(Params_Ref%Calc%Z)
        write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                         ' , N_2max = ' // istr(Params%Calc%N2max) // ' , N_3max = ' // istr(Params%Calc%N3max) // &
                         ', J-basis size = ' // istr((N_max+1)*(N_max+2)/2) // ', M-basis size = ' // istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        call cms_lines(u,Params)
        write(u,'(a)') TAB // TAB // TAB // 'CMS correction is calculated with respect to the reference (closed-subshell) nucleus! ...'
        if (Convergence) then
            write(u,'(/,/,a)') 'BCS equations converged in ' // istr(Iteration) // ' iterations with precision epsilon = ' // &
                               jstr(Params%Calc%BCS%Tol) // ' ...'
        else
            write(u,'(/,/,a)') 'BCS equations did NOT converge in ' // istr(Iteration) // ' iterations with precision epsilon = ' // &
                               jstr(Params%Calc%BCS%Tol) // ' ...'
        end if
        write(u,'(/,/,a)') 'BCS solution review:'
        write(u,'(/,a)') 'E_tot      = ' // cfmt_f(E_MF + sum(E_BCS),15,8) // M // 'Total BCS ground-state energy'
        write(u,'(a)') 'E_MF       = ' // cfmt_f(E_MF,15,8) // M // 'Mean-field ground-state energy'
        write(u,'(a)') 'E_BCS      = ' // cfmt_f(sum(E_BCS),15,8) // M // 'BCS pairing ground-state energy'
        write(u,'(a)') 'E_BCS (2N) = ' // cfmt_f(E_BCS(1),15,8) // M // 'BCS pairing NN ground-state energy'
        write(u,'(a)') 'E_BCS (3N) = ' // cfmt_f(E_BCS(2),15,8) // M // 'BCS pairing NNN ground-state energy'
        write(u,'(a)') 'T_BCS      = ' // cfmt_f(T_BCS,15,8) // M // 'BCS mean-field ground-state kinetic energy'
        write(u,'(/,a)') 'E_tot / A      = ' // cfmt_f((E_MF + sum(E_BCS)) / A,15,8) // MN // 'Total BCS ground-state energy per nucleon'
        write(u,'(a)') 'E_MF / A       = ' // cfmt_f(E_MF / A,15,8) // MN // 'Mean-field ground-state energy per nucleon'
        write(u,'(a)') 'E_BCS / A      = ' // cfmt_f(sum(E_BCS) / A,15,8) // MN // 'BCS pairing ground-state energy per nucleon'
        write(u,'(a)') 'E_BCS / A (2N) = ' // cfmt_f(E_BCS(1) / A,15,8) // MN // 'BCS pairing NN ground-state energy per nucleon'
        write(u,'(a)') 'E_BCS / A (3N) = ' // cfmt_f(E_BCS(2) / A,15,8) // MN // 'BCS pairing NNN ground-state energy per nucleon'
        write(u,'(a)') 'T_BCS / A      = ' // cfmt_f(T_BCS / A,15,8) // MN // 'BCS mean-field ground-state kinetic energy per nucleon'
        write(u,'(/,/,a)') 'Values of BCS chemical potentials:'
        write(u,'(/,a)') 'pLambda =  ' // cfmt_f(Lambda%p,15,8) // ' MeV'
        write(u,'(a)') 'nLambda =  ' // cfmt_f(Lambda%n,15,8) // ' MeV'
        write(u,'(/,a)') 'BCS particle number dispersions dZ & dN:'
        write(u,'(/,a)') 'dZ = ' // jstr(round_sig(dA%p,6))
        write(u,'(a)') 'dN = ' // jstr(round_sig(dA%n,6))
        write(u,'(/,a)') 'dZ / Z = ' // jstr(round_sig(dA%p / real(Params%Calc%Z,dp),6))
        write(u,'(a)') 'dN / N = ' // jstr(round_sig(dA%n / real(Params%Calc%A - Params%Calc%Z,dp),6))
        close(u)
    end subroutine BCS_summary

    subroutine BCS_summary_SQS(Params,SPE,SQE,U,V)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SPE, SQE, U, V
        type(orb1b_t), allocatable :: Orb(:)
        integer, allocatable :: nr(:,:)
        integer :: N_max, a_max, a, u_f, iso
        real(dp) :: O_a

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        allocate(nr(0:N_max,0:max(Params%Calc%N2max+1,2*N_max+1)))
        write(*,'(/,a)') 'Exporting information on single-(quasi)particle orbits ...'
        open(newunit=u_f, file=io_path(Params,'/BCS/BCS_Summary.dat'), position='append', action='write')
        do iso = 1, 2
            nr = 0
            if (iso == 1) then
                write(u_f,'(/,a)') 'Proton single-(quasi)particle states'
            else
                write(u_f,'(/,a)') 'Neutron single-(quasi)particle states'
            end if
            write(u_f,'(a)') '____________________________________________________________________'
            write(u_f,'(a)') cfmt_s('n',4) // ' ' // cfmt_s('l',4) // ' ' // cfmt_s('2j',4) // ' ' // cfmt_s('O_a',8) // ' ' // &
                             cfmt_s('e',14) // ' ' // cfmt_s('U2',14) // ' ' // cfmt_s('V2',14) // ' ' // cfmt_s('E',14)
            do a = 1, a_max
                if (iso == 1) then
                    O_a = merge(1.0_dp, 0.0_dp, Orb(a)%pO == 1)
                    write(u_f,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // ' ' // &
                                     cfmt_f(O_a,8,3) // ' ' // cfmt_f(SPE%p(a),14,5) // ' ' // cfmt_f(U%p(a)**2,14,3) // ' ' // &
                                     cfmt_f(V%p(a)**2,14,3) // ' ' // cfmt_f(SQE%p(a),14,5)
                else
                    O_a = merge(1.0_dp, 0.0_dp, Orb(a)%nO == 1)
                    write(u_f,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // ' ' // &
                                     cfmt_f(O_a,8,3) // ' ' // cfmt_f(SPE%n(a),14,5) // ' ' // cfmt_f(U%n(a)**2,14,3) // ' ' // &
                                     cfmt_f(V%n(a)**2,14,3) // ' ' // cfmt_f(SQE%n(a),14,5)
                end if
                nr(Orb(a)%l,Orb(a)%j) = nr(Orb(a)%l,Orb(a)%j) + 1
            end do
            write(u_f,'(a)') '____________________________________________________________________'
        end do
        close(u_f)
    end subroutine BCS_summary_SQS

    subroutine BCS_export(Params,Orb,Orb_NN,C,U,V,H_N,H_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: C, U, V
        type(qpo1b_t), intent(in) :: H_N
        type(qpo2b_t), intent(in) :: H_NN
        call O1b_export(Params,Orb,C,io_path(Params,'/Bin/C.bin'))
        call O1b_export(Params,Orb,U,io_path(Params,'/Bin/U.bin'))
        call O1b_export(Params,Orb,V,io_path(Params,'/Bin/V.bin'))
        call qpO1B_export(Params,Orb,H_N,io_path(Params,'/Bin/qpH1B.bin'))
        call qpO2b_export(Params,Orb_NN,H_NN,io_path(Params,'/Bin/qpH2B.bin'))
    end subroutine BCS_export

    subroutine cms_lines(u,Params)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Combined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Only pure 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else
            write(u,'(/,a)') TAB // TAB // TAB // 'No Center-of-Mass System (CMS) motion correction is set ...'
        end if
    end subroutine cms_lines

end module nj_bcs
