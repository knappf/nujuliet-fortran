! LEGACY (removed in cleanup stage) ... the base-version HF entry point, &run action 'HF_LEGACY' (old output layout
    ! io/<Path>/Bin/{C_HF,Rho_HF,h_HF,V2B_HF}.bin, io/<Path>/HF/HF_Summary.dat, HF-EFA, HRF export, HF-MBPT(2)) for the
    ! not-yet-ported solvers that read io/<Path>/Bin/* (HF_RPA, HF_RRPA, EMPM, NCSM) ... counterpart of base
    ! Source/HF/HF.jl & HF_solver() of Source/HF/HF_Solver.jl ... the dev HF is nj_hf_driver ...
module nj_hf_legacy_driver
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make, orbitals_export
    use nj_o1b, only : T1b, T1b_energy
    use nj_o2b, only : O2b_transformation, O2b_export, O2b_export_HRF, O2b_move
    use nj_v2b, only : V2b_read, V2b_residual_no2b
    use nj_v3b, only : V3b_no2b_read
    use nj_obdm, only : OBDM_export
    use nj_hf_legacy
    use nj_hf_mbpt_legacy, only : HF_MBPT
    use nj_hf_driver, only : print_run_header
    use nj_hf_efa, only : HF_EFA
    !$ use omp_lib
    implicit none
    private

    public :: HF_LEGACY

contains

    subroutine HF_LEGACY(Params)
        type(params_t), intent(in) :: Params
        integer(i8) :: t0

        call wigner_init_float(Params%Int%N2max + 2, 'Jmax', 6)

        write(*,'(a)') 'Starting calculations with NO2B NN+NNN interaction'
        write(*,'(a)') 'LEGACY HF (HF_LEGACY): base-version HF with the io/<Path>/Bin/* output layout ...'
        call print_run_header(Params)

        ! Fresh output directories ...
        if (is_dir(io_path(Params,''))) call remove_dir(io_path(Params,''))
        call make_dir(io_path(Params,''))
        call make_dir(io_path(Params,'/Bin'))
        call make_dir(io_path(Params,'/HF'))
        call make_dir(io_path(Params,'/HF/Densities'))
        if (Params%Calc%HF%EFA%A /= 0 .and. Params%Calc%HF%EFA%Z /= 0) then
            call make_dir(io_path(Params,'/HF-EFA'))
            call make_dir(io_path(Params,'/HF-EFA/Densities'))
        end if

        write(*,'(/,a)') 'Starting spherical Hartree-Fock calculation ...'
        call timer_start(t0)
        call HF_solver_legacy(Params)
        call timer_report(t0)

        write(*,'(/,a,/)') 'All calculations have finished ...'
    end subroutine HF_LEGACY

    subroutine HF_solver_legacy(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: T, h, C, Rho, Rho_HF
        type(o2b_t) :: V_NN, V_NN_res
        type(orb2b_t) :: Orb_NN
        type(v3b_t) :: V_NNN
        type(orb3b_t) :: Orb_NNN
        type(pnvector_t) :: SPE
        real(dp) :: E_HF, T_HF
        integer :: Iteration
        integer(i8) :: t0

        ! Single-particle orbitals ... NuHamil ordering ...
        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)

        call timer_start(t0)
        T = T1b(Params,Orb)
        call timer_report(t0)

        call timer_start(t0)
        call V2b_read(Params,V_NN,Orb_NN)
        call timer_report(t0)

        call timer_start(t0)
        call V3b_no2b_read(Params,V_NNN,Orb_NNN)
        call timer_report(t0)

        call timer_start(t0)
        call HF_solve(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,h,C,Rho,Iteration)
        call timer_report(t0)

        call timer_start(t0)
        E_HF = HF_energy(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
        call timer_report(t0)

        call timer_start(t0)
        T_HF = T1b_energy(Params,Rho,Orb,T)
        call timer_report(t0)

        call HF_summary(Params,E_HF,T_HF,Iteration)

        ! Radial HF densities & radii ...
        call timer_start(t0)
        call OBDM_export(Params,Orb,io_path(Params,'/HF/HF_Summary.dat'), &
                         io_path(Params,'/HF/Densities/HF_Radial_Densities.dat'),Rho,C)
        call timer_report(t0)

        call HF_SPS_summary(Params,h,Orb)

        ! Export of HF orbitals C, HF-basis density & HF single-particle Hamiltonian ...
        Rho_HF%p = mtmm(C%p,Rho%p)
        Rho_HF%n = mtmm(C%n,Rho%n)
        call HF_export(Params,C,Rho_HF,h,Orb)

        ! HF-EFA calculation for a neighbouring nucleus ...
        if (Params%Calc%HF%EFA%A /= 0 .and. Params%Calc%HF%EFA%Z /= 0) then
            SPE%p = diag(h%p)
            SPE%n = diag(h%n)
            call timer_start(t0)
            call HF_EFA(Params,Orb,Orb_NNN,SPE,C,V_NNN)
            call timer_report(t0)
        end if

        ! Beyond-mean-field HF-MBPT & export of the residual NN interaction ...
        if (Params%Calc%HF%BMF) then
            call timer_start(t0)
            ! LEGACY: the residual replaces V_NN here (the dev HF keeps the bare V_NN for HF_MBPT_export, HF group) ...
            call V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)
            call O2b_move(V_NN_res,V_NN)
            call timer_report(t0)

            call timer_start(t0)
            call O2b_transformation(Params,Orb,Orb_NN,V_NN,C)
            call timer_report(t0)

            call timer_start(t0)
            call O2b_export(Params,Orb_NN,V_NN,io_path(Params,'/Bin/V2B_HF.bin'))
            call timer_report(t0)

            call timer_start(t0)
            call HF_MBPT(Params,Orb,Orb_NN,V_NN)
            call timer_report(t0)

            if (Params%Calc%HF%HRF) then
                write(*,'(a)') 'Exporting the HF single-particle orbitals & residual 2-body NN interaction in the ' // &
                               'Human-Readable-Format (HRF) ...'
                call orbitals_export(Params,Orb)
                call O2b_export_HRF(Params,Orb_NN,V_NN,io_path(Params,'/Bin/V2B_HF.dat'))
            end if
        end if
    end subroutine HF_solver_legacy

end module nj_hf_legacy_driver
