! HF-based self-consistent RPA (HF-ScRPA, NuJuliet dev, stage 4) ... counterpart of src/RPA/RPA.jl (HF_ScRPA, lines 72-141)
    ! & src/RPA/RPA_solver.jl (HF_ScRPA_solver with the local functions HF_ScRPA_solver_beta & ScRPA_terminate_summary,
    ! lines 94-418) ... replaces the LEGACY HF_RRPA chain (nj_hf_rrpa*, still reachable until the cleanup stage) ...
    ! Reads io/<Path>/<RPA.Reference>/bin/{orb,C,rho,H1b,H2b}.bin (HF or HF-MBPT), writes io/<Path>/<RPA.Name>/ in the layout of
    ! plain RPA (RPA_export with C = C_LHO_NAT, Rho = Rho_LHO, h = h_N in the NAT basis) ... a non-converged outer iteration
    ! writes only RPA_summary.dat (termination summary) ...
    ! Building blocks: nj_rpa (RPA_allocate incl. TBCF & DEV, RPA_diagonalize, RPA_energy[_QBA], RPA_sigma[_mix], RPA_rM),
    ! nj_rpa_rho (RPA_OBDM_orbitals_make, HF_ScRPA_OBDM_solve), nj_rpa_hamiltonian (RPA_Hamiltonian) ...
    ! Deviations (docs/PORTING_NOTES.md "stage 4 part 2", none changes results):
    !   - Name == Reference is refused by the input reader (&rpa, B2; dev deletes the reference before reading it) ...
    !   - ScH = false: the bare V_NN / V_NNN are not read (decision 6, B5) ... Orb_NN then comes from the import of H2b.bin,
    !     which builds the identical 2-body orbital table (O2b_initialize over the first a_max(Calc.Nmax) orbitals) ...
    !   - the redundant final RPA_Hamiltonian of dev for TBCF (B5: same Rho, Sigma & C as step 4 of the last iteration) is
    !     kept (cheap compared to the loop, keeps the printed log identical) ...
    !   - dev's unused Delta_E / Delta_E_old bookkeeping is not ported ...
module nj_scrpa
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_import, orbitals_1p1h_make, orbitals_update, a_max_of
    use nj_o1b, only : O1b_import, T1b, O1b_transformation, O1b_transpose, O1b_unit_matrix
    use nj_o2b, only : O2b_import, O2b_transformation
    use nj_v2b, only : V2b_read
    use nj_v3b, only : V3b_no2b_read
    use nj_obdm, only : rho_chR2
    use nj_tr1b
    use nj_ph_common, only : ph_rB, lstr
    use nj_rpa, only : rpa_energy_t, RPA_allocate, RPA_diagonalize, RPA_energy, RPA_energy_QBA, RPA_rM, RPA_sigma, &
                       RPA_sigma_mix
    use nj_rpa_rho, only : RPA_OBDM_orbitals_make, HF_ScRPA_OBDM_solve
    use nj_rpa_hamiltonian, only : RPA_Hamiltonian
    use nj_rpa_export, only : RPA_export
    use nj_tda_rpa_driver, only : print_header, print_threads, make_output_dirs
    implicit none
    private

    public :: HF_ScRPA, HF_ScRPA_solver

contains

    ! HF-ScRPA driver ... counterpart of HF_ScRPA(Params) (RPA.jl:72-141) ... the output directory is replaced unconditionally
        ! (dev; Name /= Reference is guaranteed by the &rpa reader) ...
    subroutine HF_ScRPA(Params)
        type(params_t), intent(in) :: Params
        character(len=:), allocatable :: IO_RPA
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        call print_header(Params,'Starting HF-based Self-consistent RPA (HF-ScRPA) calculations with residual NO2B NN+NNN ' // &
                          'interaction & for given reference solution ...')
        call print_threads()

        if (.not. is_dir(io_dir(Params,Params%Calc%RPA%Reference))) then
            write(*,'(/,a)') "Error! ... No precomputed reference solutions available in IO directory ''" // io_path(Params,'') // &
                             "'' for given input path ''" // trim(Params%Calc%RPA%Reference) // "''... e.g. run HF solver first ..."
            return
        end if
        ! B2 guard (also in the &rpa reader): never delete the reference ...
        if (trim(Params%Calc%RPA%Name) == trim(Params%Calc%RPA%Reference)) &
            call fatal('HF_ScRPA: RPA Name must differ from Reference ... ' // trim(Params%Calc%RPA%Name))

        IO_RPA = io_dir(Params,Params%Calc%RPA%Name)
        if (is_dir(IO_RPA)) call remove_dir(IO_RPA)
        call make_output_dirs(IO_RPA)

        write(*,'(/,a)') 'Starting HF-ScRPA calculation ...'
        call timer_start(t0)
        call HF_ScRPA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'HF-ScRPA calculation has finished ...'
    end subroutine HF_ScRPA

    ! Ramp (damping) of the 2-body correlation function ... beta_i = tanh((i/TBCF_beta)^2), 1 above 0.95 ... counterpart of the
        ! local HF_ScRPA_solver_beta ...
    real(dp) function scrpa_beta(Params,i) result(beta_i)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: i
        beta_i = tanh((real(i,dp) / Params%Calc%RPA%TBCF_beta)**2)
        if (beta_i > 0.95_dp) beta_i = 1.0_dp
    end function scrpa_beta

    ! HF-ScRPA solver ... counterpart of HF_ScRPA_solver (RPA_solver.jl:94-418) ...
        !   setup: reference orbitals, bare H (only for ScH), reference C, Rho, H1b, H2b, frozen 1p-1h space Orb_1p1h & OBDM
        !          particle / hole orbitals Orb_ph, initial RPA on the reference ...
        !   loop:  1) TBCF: Sigma from X, Y (ramp beta, mixing TBCF_alpha), 2) self-consistent OBDM & NAT basis
        !          (HF_ScRPA_OBDM_solve, TBCF: Sigma rotated to the new NAT basis), 4) h_N, V_NN (RPA_Hamiltonian), 5) RPA in the
        !          NAT basis, 6) relative change of E_QBA below Tol in 2 consecutive iterations (with beta = 1 for TBCF) ...
        !   final: RPA, Sigma (always), energies, transitions & export ... not converged: termination summary only ...
    subroutine HF_ScRPA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb_ref(:), Orb(:)
        type(orb2b_1p1h_t) :: Orb_1p1h
        type(orb1b_rpa_obdm_t) :: Orb_ph
        type(orb2b_t) :: Orb_NN
        type(orb3b_t) :: Orb_NNN
        type(o1b_t) :: T, C_LHO_HF, Rho_ref, h_N, Rho, C_HF_NAT, C_LHO_NAT, C_HF_NAT_old, C_step, Rho_LHO
        type(o2b_t) :: V_NN_bare, V_NN, Sigma, Sigma_old, Sigma_fin, No_Sigma
        type(v3b_t) :: V_NNN_bare
        type(tr1b_t) :: TrOp_LHO, TrOp
        type(rmat_t), allocatable :: A(:,:), B(:,:), N(:,:)
        type(cvec_t), allocatable :: E(:,:)
        type(cmat_t), allocatable :: X(:,:), Y(:,:)
        type(rpa_energy_t) :: E_corr
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        real(dp), allocatable :: E_QBA_list(:)
        logical, allocatable :: Stability_list(:)
        character(len=:), allocatable :: IO
        logical :: Stability, TBCF, ScH, Converged
        real(dp) :: epsilon, alpha_TBCF, E_QBA_old, E_QBA, dE_QBA, beta, chR2
        integer :: Iteration_max, Iteration, N_converged, N_iterations
        integer(i8) :: t0

        IO = io_dir(Params,Params%Calc%RPA%Reference)
        epsilon = Params%Calc%RPA%Tol
        Iteration_max = Params%Calc%RPA%Imax
        alpha_TBCF = Params%Calc%RPA%TBCF_alpha
        TBCF = Params%Calc%RPA%TBCF
        ScH = Params%Calc%RPA%ScH

        ! Reference single-particle orbitals ...
        Orb_ref = orbitals_import(Params,IO // '/bin/orb.bin')

        ! Bare nuclear Hamiltonian (LHO basis) ... only used for ScH = true (decision 6: not read for ScH = false) ...
        call timer_start(t0)
        T = T1b(Params,Orb_ref)
        call timer_report(t0)
        if (ScH) then
            call timer_start(t0)
            call V2b_read(Params,V_NN_bare,Orb_NN)
            call timer_report(t0)
            call timer_start(t0)
            call V3b_no2b_read(Params,V_NNN_bare,Orb_NNN)
            call timer_report(t0)
        end if

        ! Reference basis LHO -> HF, OBDM, 1-body Hamiltonian & 2-body interaction ...
        C_LHO_HF = O1b_import(Params,Orb_ref,IO // '/bin/C.bin')
        Rho_ref = O1b_import(Params,Orb_ref,IO // '/bin/rho.bin')
        h_N = O1b_import(Params,Orb_ref,IO // '/bin/H1b.bin')
        write(*,'(/,a)') 'Importing reference 2-body interaction ...'
        call timer_start(t0)
        if (ScH) then
            call O2b_import(Params,Orb_ref,IO // '/bin/H2b.bin',V_NN)
        else
            call O2b_import(Params,Orb_ref,IO // '/bin/H2b.bin',V_NN,Orb_NN)
        end if
        call timer_report(t0)

        ! Frozen 1p-1h space & particle / hole orbitals of the OBDM (reference occupations) ...
        Orb_1p1h = orbitals_1p1h_make(Params,Orb_ref,Params%Calc%RPA%Tol_par,Params%Calc%RPA%Tol_hol)
        Orb_ph = RPA_OBDM_orbitals_make(Params,Orb_ref)

        ! 1-body transition operators (LHO basis, unit charge radius) & in the reference basis ...
        TrOp_LHO = Tr1b_initialize(Params,Orb_ref,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb_ref,C_LHO_HF,TrOp_LHO)

        ! Initial RPA on top of the reference solution ...
        write(*,'(/,a)') 'Solving the initial RPA equations ...'
        call timer_start(t0)
        call RPA_allocate(Params,Orb_ref,Orb_1p1h,Orb_NN,h_N,V_NN,Rho_ref,No_Sigma,.false.,A,B,N)
        call timer_report(t0)
        call timer_start(t0)
        call RPA_diagonalize(Params,Orb_ref,Orb_1p1h,A,B,N,TrOp,E,X,Y,Stability)
        call timer_report(t0)
        deallocate(A, B, N)
        if (.not. Stability) write(*,'(/,a)') TAB // 'Warning! The initial RPA solution is unstable ...'

        E_QBA_old = RPA_energy_QBA(Params,Orb_ref,Orb_1p1h,E,X,Y)
        allocate(E_QBA_list(0:max(Iteration_max,0)), Stability_list(max(Iteration_max,1)))
        E_QBA_list(0) = E_QBA_old

        ! Self-consistent iteration ... start from the reference solution ...
        Orb = Orb_ref
        Rho = Rho_ref
        C_HF_NAT = O1b_unit_matrix(Params)
        C_LHO_NAT = C_LHO_HF

        Iteration = 1
        Converged = .false.
        E_QBA = E_QBA_old
        dE_QBA = 1.0_dp
        N_converged = 0
        beta = scrpa_beta(Params,1)

        write(*,'(/,a)') 'Starting HF-ScRPA iteration ...'
        do while (Iteration <= Iteration_max)
            C_HF_NAT_old = C_HF_NAT

            ! 1) 2-body correlation function from the current amplitudes (current NAT basis), damped by beta & mixed ...
            if (TBCF) then
                call timer_start(t0)
                Sigma = RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y,beta=beta,Report=.false.)
                if (Iteration > 1) call RPA_sigma_mix(Sigma,Sigma_old,alpha_TBCF)
                call timer_report(t0)
            end if

            ! 2) Self-consistent OBDM ... Rho, Orb (occupations) & C_HF_NAT updated, C_LHO_NAT = C_LHO_HF C_HF_NAT ...
            call timer_start(t0)
            call HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT,C_LHO_NAT)
            call timer_report(t0)

            ! Sigma into the new NAT basis (C_step = old NAT -> new NAT) ... Sigma_old = Sigma ...
            if (TBCF) then
                C_step%p = matmul(transpose(C_HF_NAT_old%p),C_HF_NAT%p)
                C_step%n = matmul(transpose(C_HF_NAT_old%n),C_HF_NAT%n)
                call timer_start(t0)
                call O2b_transformation(Params,Orb,Orb_NN,Sigma,C_step,Report=.false.)
                call timer_report(t0)
                Sigma_old = Sigma
            end if

            ! 4) Hamiltonian & transition operators in the new NAT basis ...
            call timer_start(t0)
            call RPA_Hamiltonian(Params,Orb,Orb_NN,Orb_NNN,C_LHO_NAT,Rho,Sigma,T,V_NN_bare,V_NNN_bare,h_N,V_NN)
            call timer_report(t0)
            TrOp = Tr1b_transformation(Params,Orb,C_LHO_NAT,TrOp_LHO)

            ! 5) RPA in the new NAT basis with the new occupations ...
            call timer_start(t0)
            call RPA_allocate(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Rho,Sigma,TBCF,A,B,N)
            call timer_report(t0)
            call timer_start(t0)
            call RPA_diagonalize(Params,Orb,Orb_1p1h,A,B,N,TrOp,E,X,Y,Stability)
            call timer_report(t0)
            deallocate(A, B, N)

            ! 6) Relative change of the quasi-boson correlation energy ...
            E_QBA = RPA_energy_QBA(Params,Orb,Orb_1p1h,E,X,Y)
            if (abs(E_QBA) > 1.0e-12_dp) then
                dE_QBA = abs(E_QBA - E_QBA_old) / abs(E_QBA)
            else
                dE_QBA = abs(E_QBA - E_QBA_old)
            end if
            E_QBA_old = E_QBA
            E_QBA_list(Iteration) = E_QBA
            Stability_list(Iteration) = Stability

            write(*,'(a)') TAB // 'HF-ScRPA iteration:   Iteration number = ' // cfmt_d(Iteration,4) // '   E0 = ' // &
                           cfmt_f(E_QBA,12,9) // ' MeV   dE0_rel = ' // cfmt_f(dE_QBA,12,9)

            if (dE_QBA < epsilon .and. (.not. TBCF .or. beta == 1.0_dp)) then
                N_converged = N_converged + 1
            else
                N_converged = 0
            end if
            if (N_converged >= 2 .and. (.not. TBCF .or. beta == 1.0_dp)) then
                Converged = .true.
                exit
            end if

            beta = scrpa_beta(Params,Iteration + 1)
            Iteration = Iteration + 1
        end do

        if (Converged) then
            N_iterations = Iteration
        else
            N_iterations = Iteration - 1
        end if
        write(*,'(/,a)') 'HF-ScRPA iteration finished after ' // cfmt_d(N_iterations,4) // ' iterations ...   E0 = ' // &
                         cfmt_f(E_QBA,12,9) // ' MeV   dE0_rel = ' // cfmt_f(dE_QBA,12,9)
        if (Converged) then
            write(*,'(/,a)') 'HF-ScRPA iteration has converged ...'
        else
            write(*,'(/,a)') 'HF-ScRPA iteration has NOT converged ... Run terminated ...'
            call ScRPA_terminate_summary(Params,N_iterations,E_QBA_list(0:N_iterations),Stability_list(1:N_iterations), &
                                         scrpa_beta(Params,N_iterations))
            return
        end if

        ! Final RPA for the final OBDM (& correlation function), final NAT basis ...
        write(*,'(/,a)') 'Solving the final HF-ScRPA equations ...'
        if (TBCF) then
            call timer_start(t0)
            call RPA_Hamiltonian(Params,Orb,Orb_NN,Orb_NNN,C_LHO_NAT,Rho,Sigma,T,V_NN_bare,V_NNN_bare,h_N,V_NN)
            call timer_report(t0)
        end if
        if (allocated(Sigma_old%pp)) call o2b_free(Sigma_old)

        call timer_start(t0)
        call RPA_allocate(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Rho,Sigma,TBCF,A,B,N)
        call timer_report(t0)
        call timer_start(t0)
        call RPA_diagonalize(Params,Orb,Orb_1p1h,A,B,N,TrOp,E,X,Y,Stability)
        call timer_report(t0)
        deallocate(A, B, N)
        if (.not. Stability) write(*,'(/,a)') 'Warning! The final HF-ScRPA solution is unstable ...'

        ! 2-body correlation function of the final amplitudes (always, beta = 1; Sigma stays the Sigma_AB of A & B) ...
        call timer_start(t0)
        Sigma_fin = RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y)
        call timer_report(t0)

        ! Transition operators with the HF-ScRPA charge radius ...
        chR2 = rho_chR2(Params,Orb,C_LHO_NAT,Rho)
        TrOp = Tr1b_initialize(Params,Orb,chR2)
        TrOp = Tr1b_transformation(Params,Orb,C_LHO_NAT,TrOp)

        ! Occupations from the final OBDM (no-op, Orb already carries diag(Rho); as dev) ...
        call orbitals_update(Orb,Rho)

        ! Correlation energies (Rho - Rho_ref in the final NAT basis) ...
        call timer_start(t0)
        E_corr = RPA_energy(Params,Orb,Orb_1p1h,Orb_NN,E,X,Y,Rho,O1b_transformation(Rho_ref,C_HF_NAT),h_N,Sigma_fin,V_NN, &
                            TBCF,Sigma,Density=.true.)
        call timer_report(t0)

        call timer_start(t0)
        rM = RPA_rM(Params,Orb,Orb_1p1h,Orb_NN,TrOp,Rho,Sigma_fin,X,Y)
        rB = ph_rB(Params,Orb_1p1h,rM,real(E(2,1)%v,dp))
        call timer_report(t0)

        ! OBDM in the LHO basis ... C_LHO_NAT Rho C_LHO_NAT^T ...
        Rho_LHO = O1b_transformation(Rho,O1b_transpose(C_LHO_NAT))

        call RPA_export(Params,Orb,Orb_1p1h,Orb_NN,Stability,E_corr,E,X,Y,rB,C_LHO_NAT,Rho_LHO,Sigma_fin,h_N)
    end subroutine HF_ScRPA_solver

    subroutine o2b_free(O)
        type(o2b_t), intent(inout) :: O
        if (allocated(O%pp)) deallocate(O%pp)
        if (allocated(O%pn)) deallocate(O%pn)
        if (allocated(O%nn)) deallocate(O%nn)
    end subroutine o2b_free

    ! Summary of a non-converged HF-ScRPA iteration io/<Path>/<Name>/RPA_summary.dat ... counterpart of the local
        ! ScRPA_terminate_summary (RPA_solver.jl:117-174) ... E_QBA_list(0:N) (0 = initial RPA), Stability_list(1:N) ...
    subroutine ScRPA_terminate_summary(Params,N_iterations,E_QBA_list,Stability_list,beta)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_iterations
        real(dp), intent(in) :: E_QBA_list(0:)
        logical, intent(in) :: Stability_list(:)
        real(dp), intent(in) :: beta
        character(len=:), allocatable :: Output_Path
        integer :: u, N_max, a_max, N_E, i
        real(dp) :: Delta_E, dE_rel

        N_max = Params%Calc%Nmax
        a_max = a_max_of(N_max)
        Output_Path = io_dir(Params,Params%Calc%RPA%Name) // '/RPA_summary.dat'

        write(*,'(/,a)') 'Writing HF-ScRPA summary file (non-converged iteration) ...'
        associate (R => Params%Calc%RPA)
            open(newunit=u, file=Output_Path, status='replace', action='write')
            write(u,'(a)') 'Spherical self-consistent Random-Phase Approximation (HF-ScRPA) review:'
            write(u,'(/,a)') '!!! The HF-ScRPA iteration has NOT converged ... no RPA results were exported !!!'
            write(u,'(/,a)') 'Nuclide data:    A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)
            write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                             ', s.p. J-basis size = ' // istr(a_max) // ', s.p. M-basis size = ' // &
                             istr((N_max+1)*(N_max+2)*(N_max+3)/6)
            write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
            if (R%Ortho) then
                write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (1-) is included ...'
            else
                write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (1-) is NOT included ...'
            end if

            write(u,'(/,a)') 'HF-ScRPA parameters:'
            write(u,'(a)') TAB // "Reference = ''" // trim(R%Reference) // "'',   Reference_TBCF = " // lstr(R%Reference_TBCF) // &
                           ',   ScH = ' // lstr(R%ScH)
            write(u,'(a)') TAB // 'Tol = ' // jstr(R%Tol) // ',   Imax = ' // istr(R%Imax) // ',   Tol_hol = ' // jstr(R%Tol_hol) // &
                           ',   Tol_par = ' // jstr(R%Tol_par)
            write(u,'(a)') TAB // 'OBDM_alpha = ' // jstr(R%OBDM_alpha) // ',   OBDM_diagonal = ' // lstr(R%OBDM_diagonal)
            write(u,'(a)') TAB // 'TBCF = ' // lstr(R%TBCF) // ',   TBCF_alpha = ' // jstr(R%TBCF_alpha) // ',   TBCF_beta = ' // &
                           jstr(R%TBCF_beta)

            write(u,'(/,a)') 'Convergence information:'
            write(u,'(a)') TAB // 'Iterations performed = ' // istr(N_iterations) // ' (Imax = ' // istr(R%Imax) // ')'
            if (R%TBCF) then
                write(u,'(a)') TAB // 'Criterion ... relative change of E_QBA below Tol in two consecutive iterations, with beta = 1 ...'
                write(u,'(a)') TAB // 'Last beta (ramp of Sigma) = ' // jstr(beta)
            else
                write(u,'(a)') TAB // 'Criterion ... relative change of E_QBA below Tol in two consecutive iterations ...'
            end if
            write(u,'(a)') TAB // 'Unstable iterations (complex RPA eigenvalues) = ' // istr(count(.not. Stability_list)) // ' of ' // &
                           istr(size(Stability_list))

            N_E = size(E_QBA_list)
            if (N_E >= 2) then
                Delta_E = E_QBA_list(N_E-1) - E_QBA_list(N_E-2)
                write(u,'(a)') TAB // 'Last E0 = ' // jstr(E_QBA_list(N_E-1)) // ' MeV,   last change = ' // jstr(abs(Delta_E)) // &
                               ' MeV (relative ' // jstr(abs(Delta_E) / max(abs(E_QBA_list(N_E-1)),1.0e-12_dp)) // ')'
            end if

            write(u,'(/,a)') 'Last iterations ... (Iteration, E_QBA [MeV], relative change, stable)'
            ! Julia i = max(2,N_E-9):N_E (1-based list) ... here k = i - 1 = iteration number, list index k ...
            do i = max(2,N_E-9), N_E
                dE_rel = abs(E_QBA_list(i-1) - E_QBA_list(i-2)) / max(abs(E_QBA_list(i-1)),1.0e-12_dp)
                write(u,'(a)') istr(i-1) // TAB // jstr(E_QBA_list(i-1)) // TAB // jstr(dE_rel) // TAB // lstr(Stability_list(i-1))
            end do
            close(u)
        end associate
        write(*,'(a)') TAB // "HF-ScRPA summary file exported ... ''" // Output_Path // "''"
    end subroutine ScRPA_terminate_summary

end module nj_scrpa
