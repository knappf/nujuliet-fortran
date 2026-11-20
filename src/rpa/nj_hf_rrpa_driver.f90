! HF-RRPA entry point ... counterpart of HF_RRPA/HF_RRPA.jl & HF_RRPA_Solver.jl ...
module nj_hf_rrpa_driver
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make, orbitals_ph_make, orbitals_one_phonon_make
    use nj_o1b, only : O1b_import, T1b
    use nj_o2b, only : O2b_import
    use nj_v2b, only : V2b_read
    use nj_v3b, only : V3b_no2b_read
    use nj_obdm, only : OBDM_chR2
    use nj_tr1b
    use nj_hf_rpa, only : HF_RPA_phonon_count, HF_RPA_allocate, HF_RPA_rB
    use nj_hf_rrpa
    use nj_hf_rrpa_export, only : HF_RRPA_export
    !$ use omp_lib
    implicit none
    private

    public :: HF_RRPA

contains

    subroutine HF_RRPA(Params)
        type(params_t), intent(in) :: Params
        integer :: N_max, N_threads
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        N_max = Params%Calc%Nmax
        write(*,'(a)') 'Starting Renormalized-RPA calculations with residual NO2B NN+NNN interaction'
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
        write(*,'(/,a)') "Input/Output path to the calculation reads ... ''io/" // trim(Params%Calc%Path) // "''"
        if (.not. is_dir(io_path(Params,''))) then
            write(*,'(/,a)') 'Error! ... No precomputed HF solution is available in given Input_File path ... run HF solver first ...'
            return
        end if
        call remove_dir(io_path(Params,'/RRPA'))
        call make_dir(io_path(Params,'/RRPA/Densities'))
        call make_dir(io_path(Params,'/RRPA/Transitions/E0'))
        call make_dir(io_path(Params,'/RRPA/Transitions/E1'))
        call make_dir(io_path(Params,'/RRPA/Transitions/E2'))
        call make_dir(io_path(Params,'/RRPA/Transitions/E3'))
        call make_dir(io_path(Params,'/RRPA/Spectra'))
        call make_dir(io_path(Params,'/RRPA/Amplitudes'))

        write(*,'(/,a)') 'Starting RRPA calculation ...'
        call timer_start(t0)
        call HF_RRPA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'All calculations have finished ...'
    end subroutine HF_RRPA

    subroutine HF_RRPA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: T, Rho, C, h
        type(o2b_t) :: V_NN_Bare
        type(orb2b_t) :: Orb_NN_Bare
        type(v3b_t) :: V_NNN_Bare
        type(orb3b_t) :: Orb_NNN_Bare
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        type(phstate_t), allocatable :: Phonon(:)
        integer :: N_Phonon, J, P, J_max
        integer, allocatable :: N_nu(:,:)
        type(ivec_t), allocatable :: Orb_Phonon(:,:)
        type(cvec_t), allocatable :: E_RPA(:,:)
        type(cmat_t), allocatable :: X_RPA(:,:), Y_RPA(:,:)
        type(rvec_t), allocatable :: E_re(:,:)
        type(tr1b_t) :: TrOp
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        real(dp) :: chR2, E_corr
        integer(i8) :: t0

        J_max = 2*Params%Calc%Nmax + 1
        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call timer_start(t0)
        T = T1b(Params,Orb)
        call timer_report(t0)
        call timer_start(t0)
        call V2b_read(Params,V_NN_Bare,Orb_NN_Bare)
        call timer_report(t0)
        call timer_start(t0)
        call V3b_no2b_read(Params,V_NNN_Bare,Orb_NNN_Bare)
        call timer_report(t0)
        call orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        call orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole,N_Phonon,Phonon)
        call HF_RPA_phonon_count(Params,N_Phonon,Phonon,N_nu,Orb_Phonon)

        call timer_start(t0)
        call HF_RRPA_solve(Params,N_nu,Orb_Phonon,Phonon,N_Particle,Particle,N_Hole,Hole,Orb,Orb_NN_Bare,Orb_NNN_Bare,T, &
                           V_NN_Bare,V_NNN_Bare,E_RPA,X_RPA,Y_RPA,Rho,C,h)
        call timer_report(t0)

        chR2 = OBDM_chR2(Params,Orb,C,Rho)
        TrOp = Tr1b_initialize(Params,Orb,chR2)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        call HF_RRPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Rho,TrOp,X_RPA,Y_RPA,rM)
        allocate(E_re(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                E_re(J,P)%v = real(E_RPA(J,P)%v)
            end do
        end do
        rB = HF_RPA_rB(Params,N_nu,rM,E_re)
        E_corr = HF_RRPA_energy(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)

        call timer_start(t0)
        call HF_RRPA_export(Params,Orb,N_nu,E_corr,C,Rho,h,E_RPA,X_RPA,Y_RPA,rB)
        call timer_report(t0)
    end subroutine HF_RRPA_solver

    ! Self-consistent RRPA iteration ... OBDM, effective 1-body Hamiltonian & residual interaction ...
    subroutine HF_RRPA_solve(Params,N_nu,Orb_Phonon,Phonon,N_Particle,Particle,N_Hole,Hole,Orb,Orb_NN_Bare,Orb_NNN_Bare,T, &
                             V_NN_Bare,V_NNN_Bare,E_RPA,X_RPA,Y_RPA,Rho,C_tot,H)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN_Bare
        type(orb3b_t), intent(in) :: Orb_NNN_Bare
        type(o1b_t), intent(in) :: T
        type(o2b_t), intent(in) :: V_NN_Bare
        type(v3b_t), intent(in) :: V_NNN_Bare
        type(cvec_t), allocatable, intent(out) :: E_RPA(:,:)
        type(cmat_t), allocatable, intent(out) :: X_RPA(:,:), Y_RPA(:,:)
        type(o1b_t), intent(out) :: Rho, C_tot, H
        type(o2b_t) :: V_NN_Res
        type(orb2b_t) :: Orb_NN_Res
        type(o1b_t) :: C_LHO_HF, C_HF_RRPA, Rho_LHO, Rho_HF, Rho_0, C, C_LHO_RRPA
        type(rmat_t), allocatable :: A(:,:), B(:,:)
        type(cvec_t), allocatable :: E_new_v(:,:)
        type(cmat_t), allocatable :: X_new(:,:), Y_new(:,:)
        type(tr1b_t) :: TrOp
        real(dp) :: eta, epsilon, E_old, E_new, Sum
        integer :: a_max, ia, ib, Iteration, Iteration_max
        integer(i8) :: t0

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        eta = 1.0_dp
        epsilon = Params%Calc%RRPA%Tol
        Iteration = 0
        Iteration_max = Params%Calc%RRPA%Imax
        E_old = 2.0_dp
        E_new = 1.0_dp

        write(*,'(/,a)') 'Importing Residual 2-body interaction ...'
        call timer_start(t0)
        call O2b_import(Params,Orb,io_path(Params,'/Bin/V2B_HF.bin'),V_NN_Res,Orb_NN_Res)
        call timer_report(t0)
        C_LHO_HF = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        C_HF_RRPA%p = identity(a_max)
        C_HF_RRPA%n = identity(a_max)
        allocate(Rho%p(a_max,a_max), Rho%n(a_max,a_max))
        Rho%p = 0.0_dp
        Rho%n = 0.0_dp
        allocate(Rho_HF%p(a_max,a_max), Rho_HF%n(a_max,a_max))
        Rho_HF%p = 0.0_dp
        Rho_HF%n = 0.0_dp
        do ia = 1, a_max
            Rho_HF%p(ia,ia) = real(Orb(ia)%pO,dp)
            Rho_HF%n(ia,ia) = real(Orb(ia)%nO,dp)
        end do
        Rho_0 = Rho_HF
        allocate(H%p(a_max,a_max), H%n(a_max,a_max))
        H%p = 0.0_dp
        H%n = 0.0_dp

        ! Standard RPA ... RRPA(0) ...
        write(*,'(/,a)') 'Allocating RPA matrices A & B ...'
        call timer_start(t0)
        call HF_RPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN_Res,V_NN_Res,A,B)
        call timer_report(t0)
        write(*,'(/,a)') 'Solving the RPA Generalized-Eigenvalue Problem ...'
        call timer_start(t0)
        call HF_RRPA0_diagonalize(Params,A,B,N_nu,E_RPA,X_RPA,Y_RPA)
        call timer_report(t0)

        do while ((eta > epsilon) .and. (Iteration < Iteration_max))
            call timer_start(t0)
            call HF_RRPA_OBDM_iteration(Params,N_nu,N_Particle,Particle,N_Hole,Hole,X_RPA,Y_RPA,Rho_HF,Rho_0,Rho,C)
            call timer_report(t0)

            C_HF_RRPA%p = mm(C_HF_RRPA%p,C%p)
            C_HF_RRPA%n = mm(C_HF_RRPA%n,C%n)
            Rho_0%p = mtmm(C%p,Rho_0%p)
            Rho_0%n = mtmm(C%n,Rho_0%n)

            C_LHO_RRPA%p = mm(C_LHO_HF%p,C_HF_RRPA%p)
            C_LHO_RRPA%n = mm(C_LHO_HF%n,C_HF_RRPA%n)
            Rho_LHO%p = mmt(mm(C_LHO_RRPA%p,Rho%p),C_LHO_RRPA%p)
            Rho_LHO%n = mmt(mm(C_LHO_RRPA%n,Rho%n),C_LHO_RRPA%n)
            Rho_LHO%p = 0.5_dp * (Rho_LHO%p + transpose(Rho_LHO%p))
            Rho_LHO%n = 0.5_dp * (Rho_LHO%n + transpose(Rho_LHO%n))

            ! Effective 1-body Hamiltonian in the LHO basis ... then in the HF-RRPA basis ...
            if (Params%Calc%RRPA%ScOBH) then
                H = HF_RRPA_h1b(Params,Rho_LHO,Orb,Orb_NN_Bare,Orb_NNN_Bare,T,V_NN_Bare,V_NNN_Bare)
            else
                H = HF_RRPA_h1b(Params,Rho_HF,Orb,Orb_NN_Bare,Orb_NNN_Bare,T,V_NN_Bare,V_NNN_Bare)
            end if
            H%p = mtmm(C_LHO_RRPA%p,H%p)
            H%n = mtmm(C_LHO_RRPA%n,H%n)
            H%p = 0.5_dp * (H%p + transpose(H%p))
            H%n = 0.5_dp * (H%n + transpose(H%n))

            write(*,'(/,a)') 'Updating the residual interaction ...'
            call timer_start(t0)
            if (Params%Calc%RRPA%ScV3N) then
                call HF_RRPA_V2b(Params,Orb,Orb_NN_Bare,Orb_NNN_Bare,V_NN_Bare,V_NNN_Bare,Rho_LHO,C_LHO_RRPA,V_NN_Res)
            else
                ! The original tests a non-existent RRPA_Parameters field here ... the basis update is applied ...
                call HF_RRPA_V2b_transform(Params,Orb,Orb_NN_Res,V_NN_Res,C)
            end if
            call timer_report(t0)

            write(*,'(/,a)') 'Allocating RPA matrices A & B ...'
            call timer_start(t0)
            call HF_RRPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN_Res,V_NN_Res,Rho,H,A,B)
            call timer_report(t0)

            write(*,'(/,a)') 'Solving the RPA Generalized-Eigenvalue Problem ...'
            call timer_start(t0)
            call HF_RRPA0_diagonalize(Params,A,B,N_nu,E_RPA,X_RPA,Y_RPA)
            call timer_report(t0)

            Sum = 0.0_dp
            do ia = 1, a_max
                do ib = 1, a_max
                    Sum = Sum + (abs(Rho%p(ia,ib) - Rho_0%p(ia,ib)) + abs(Rho%n(ia,ib) - Rho_0%n(ia,ib))) / real(2*a_max**2,dp)
                end do
            end do

            E_new = HF_RRPA_energy(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA,Quiet=.true.)

            ! E_old is never updated in the original ... kept as is ...
            if (abs(Sum - eta) < epsilon .and. (Iteration > 50) .and. (E_new < E_old)) then
                write(*,'(/,a)') 'HF-RRPA iteration finished - cycled degenerate solution found ...'
                write(*,'(a)') TAB // 'HThe solution with lower total energy was picked ...'
                exit
            else if ((abs(Sum) < 1.0e2_dp * epsilon) .and. (Iteration > 75) .and. (E_new < E_old)) then
                write(*,'(/,a)') 'RRPA iteration has finished ... Solution converged with lesser precisions 100.0 * Tol ...'
                write(*,'(a)') TAB // 'HThe solution with lower total energy was picked ...'
                exit
            else if ((abs(Sum) < 1.0e3_dp * epsilon) .and. (Iteration > 100) .and. (E_new < E_old)) then
                write(*,'(/,a)') 'RRPA iteration has finished ... Solution converged with lesser precisions 1000.0 * Tol ...'
                write(*,'(a)') TAB // 'HThe solution with lower total energy was picked ...'
                exit
            end if

            eta = Sum
            Iteration = Iteration + 1
            write(*,'(a)') TAB // 'RRPA iteration:   Iteration Number = ' // istr(Iteration) // ',' // TAB // 'eta = ' // jstr(eta)
            Rho_0 = Rho
        end do

        write(*,'(/,a)') 'HF-RRPA iteration has terminated ...'

        ! Final calculation with/without orthogonalization of the spurious state ...
        write(*,'(/,a)') 'Allocating RPA matrices A & B ...'
        call HF_RRPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN_Res,V_NN_Res,Rho,H,A,B)
        C_tot%p = mm(C_LHO_HF%p,C_HF_RRPA%p)
        C_tot%n = mm(C_LHO_HF%n,C_HF_RRPA%n)
        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb,C_tot,TrOp)
        write(*,'(/,a)') 'Solving the RPA Generalized-Eigenvalue Problem ...'
        call HF_RRPA_diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,Rho,X_RPA,Y_RPA,E_new_v,X_new,Y_new)
        call move_alloc(E_new_v,E_RPA)
        call move_alloc(X_new,X_RPA)
        call move_alloc(Y_new,Y_RPA)

        ! Intermediate (printed) correlation energies ... as in the original ...
        call HF_RRPA_energy_bosonic(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)
        call HF_RRPA_energy_density(Params,N_nu,N_Particle,Particle,N_Hole,Hole,A,E_RPA,Y_RPA)
    end subroutine HF_RRPA_solve

end module nj_hf_rrpa_driver
