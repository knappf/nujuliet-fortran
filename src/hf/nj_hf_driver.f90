! HF entry point (NuJuliet dev) ... counterpart of src/HF/HF.jl & HF_solver() of src/HF/HF_solver.jl ...
    ! Outputs in io/<Path>/HF/{HF_summary.dat, bin/{orb,C,rho,H1b,H2b}.bin, densities/*} ... the base-version HF output
    ! (io/<Path>/Bin/*, HF-EFA, HRF tables) is gone, as in dev ...
module nj_hf_driver
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make
    use nj_o1b, only : T1b, T1b_energy
    use nj_o2b, only : O2b_transformation, O2b_export
    use nj_v2b, only : V2b_read, V2b_residual_no2b
    use nj_v3b, only : V3b_no2b_read
    use nj_obdm, only : rho_export
    use nj_hf
    use nj_hf_mbpt, only : HF_MBPT, HF_MBPT_occupation_check, HF_MBPT_export
    !$ use omp_lib
    implicit none
    private

    public :: HF, HF_solver, print_run_header

contains

    ! Common header printed by the solvers ... interaction & calculation data, CMS option & threads ...
    subroutine print_run_header(Params)
        type(params_t), intent(in) :: Params
        integer :: N_int, N_calc, N_threads

        N_int = Params%Int%Nmax
        N_calc = Params%Calc%Nmax
        write(*,'(/,a)') 'Interaction data:'
        write(*,'(a)') 'hw = ' // jstr(Params%Int%hw) // ' , N_max = ' // istr(N_int) // ' , N_2max = ' // istr(Params%Int%N2max) // &
                       ' , N_3max = ' // istr(Params%Int%N3max) // ', J-scheme basis size = ' // istr((N_int+1)*(N_int+2)/2) // &
                       ', M-scheme basis size = ' // istr((N_int+1)*(N_int+2)*(N_int+3)/6)
        write(*,'(/,a)') 'Calculation data:'
        write(*,'(a)') 'A = ' // istr(Params%Calc%A) // ' , Z = ' // istr(Params%Calc%Z) // ' , N_max = ' // istr(N_calc) // &
                       ' , N_2max = ' // istr(Params%Calc%N2max) // ' , N_3max = ' // istr(Params%Calc%N3max) // &
                       ', J-scheme basis size = ' // istr((N_calc+1)*(N_calc+2)/2) // ', M-scheme basis size = ' // &
                       istr((N_calc+1)*(N_calc+2)*(N_calc+3)/6)
        write(*,'(a)') 'Center of mass correction option is set to:    ' // trim(Params%Calc%CMS)
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(*,'(/,a)') 'Combined 1-body + 2-body center of mass motion correction is included ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(*,'(/,a)') 'Only pure 2-body center of mass motion correction is included ...'
        else
            write(*,'(/,a)') 'No center of mass motion correction is included ...'
        end if

        N_threads = 1
        !$ N_threads = omp_get_max_threads()
        if (N_threads > 1) then
            write(*,'(/,a)') '_________________________________________________________'
            write(*,'(a)') 'Multiple active threads detected ...'
            write(*,'(a)') 'Parallelization report:    Number of active threads = ' // istr(N_threads)
            write(*,'(a)') '_________________________________________________________'
        end if
    end subroutine print_run_header

    subroutine HF(Params)
        type(params_t), intent(in) :: Params
        integer(i8) :: t0

        call wigner_init_float(Params%Int%N2max + 2, 'Jmax', 6)

        write(*,'(a)') 'Starting calculations with NO2B NN+NNN interaction'
        call print_run_header(Params)

        ! Fresh output directories ... the whole io/<Path> is wiped (as in dev, including other solvers' output) ...
        if (is_dir(io_path(Params,''))) call remove_dir(io_path(Params,''))
        call make_dir(io_path(Params,''))
        call make_dir(io_path(Params,'/HF'))
        call make_dir(io_path(Params,'/HF/bin'))
        call make_dir(io_path(Params,'/HF/densities'))

        write(*,'(/,a)') 'Starting spherical Hartree-Fock calculation ...'
        call timer_start(t0)
        call HF_solver(Params)
        call timer_report(t0)

        write(*,'(/,a,/)') 'All calculations have finished ...'
    end subroutine HF

    subroutine HF_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: T, h, C, Rho, Rho_HF
        type(o2b_t) :: V_NN, V_NN_res
        type(orb2b_t) :: Orb_NN
        type(v3b_t) :: V_NNN
        type(orb3b_t) :: Orb_NNN
        real(dp) :: E_HF, T_HF
        integer :: Iteration
        integer(i8) :: t0

        ! Single-particle LHO orbitals (Calc.Nmax) ... NuHamil ordering, naive (possibly fractional) filling ...
        call timer_start(t0)
        Orb = orbitals_make(Params)
        call timer_report(t0)

        call timer_start(t0)
        T = T1b(Params,Orb)
        call timer_report(t0)

        call timer_start(t0)
        call V2b_read(Params,V_NN,Orb_NN)
        call timer_report(t0)

        call timer_start(t0)
        call V3b_no2b_read(Params,V_NNN,Orb_NNN)
        call timer_report(t0)

        ! Self-consistent HF ... Orb carries the final HF occupations afterwards (dev in-place update) ...
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
        call rho_export(Params,Orb,io_path(Params,'/HF/HF_summary.dat'), &
                        io_path(Params,'/HF/densities/HF_radial_densities.dat'),Rho,C)
        call timer_report(t0)

        call HF_SPS_summary(Params,h,Orb)

        ! Export of the HF orbitals (orb.bin), C, HF-basis density C^T Rho C & HF single-particle Hamiltonian ...
        Rho_HF%p = mtmm(C%p,Rho%p)
        Rho_HF%n = mtmm(C%n,Rho%n)
        call HF_export(Params,C,Rho_HF,h,Orb)

        ! Non-local HF radial potential ...
        if (Params%Calc%HF%V_HF) then
            call timer_start(t0)
            call HF_radial_potential(Params,Orb,C,T,h)
            call timer_report(t0)
        end if

        ! Beyond-mean-field HF-MBPT & export of the residual NN interaction ...
        if (Params%Calc%HF%BMF) then
            ! Density-dependent residual NN interaction (NO2B) ... V_NN stays the bare interaction (HF_MBPT_export) ...
            call timer_start(t0)
            call V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)
            call timer_report(t0)
            ! Port: the bare V_NN is only needed by HF_MBPT_export ... freed early otherwise (memory) ...
            if (.not. Params%Calc%HF%MBPT%Export) deallocate(V_NN%pp, V_NN%pn, V_NN%nn)

            call timer_start(t0)
            call O2b_transformation(Params,Orb,Orb_NN,V_NN_res,C)
            call timer_report(t0)

            call timer_start(t0)
            call O2b_export(Params,Orb_NN,V_NN_res,io_path(Params,'/HF/bin/H2b.bin'))
            call timer_report(t0)

            if (HF_MBPT_occupation_check(Params,Rho_HF)) then
                call timer_start(t0)
                call HF_MBPT(Params,Orb,Orb_NN,h,V_NN_res)
                call timer_report(t0)

                if (Params%Calc%HF%MBPT%Export) then
                    call timer_start(t0)
                    call HF_MBPT_export(Params,Orb,Orb_NN,Orb_NNN,C,T,h,V_NN,V_NN_res,V_NNN)
                    call timer_report(t0)
                end if
            end if
        end if
    end subroutine HF_solver

end module nj_hf_driver
