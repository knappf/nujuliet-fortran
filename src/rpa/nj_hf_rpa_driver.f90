! HF-TDA & HF-RPA entry point ... counterpart of HF_RPA/HF_RPA.jl & HF_RPA_Solver.jl ...
module nj_hf_rpa_driver
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make, orbitals_ph_make, orbitals_one_phonon_make
    use nj_o1b, only : O1b_import
    use nj_o2b, only : O2b_import
    use nj_obdm, only : OBDM_chR2
    use nj_tr1b
    use nj_hf_rpa
    use nj_hf_rpa_export, only : HF_RPA_export
    !$ use omp_lib
    implicit none
    private

    public :: HF_RPA, HF_RPA_solver

contains

    subroutine HF_RPA(Params)
        type(params_t), intent(in) :: Params
        integer :: N_max, N_threads
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        N_max = Params%Calc%Nmax
        write(*,'(a)') 'Starting RPA (TDA) calculations with residual NO2B NN+NNN interaction'
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
        if (.not. is_dir(io_path(Params,''))) then
            write(*,'(/,a)') 'Error! ... No precomputed HF solution is available in given Input_File path ... run HF solver first ...'
            return
        end if
        call remove_dir(io_path(Params,'/RPA'))
        call make_dir(io_path(Params,'/RPA/Densities'))
        call make_dir(io_path(Params,'/RPA/Transitions/E0'))
        call make_dir(io_path(Params,'/RPA/Transitions/E1'))
        call make_dir(io_path(Params,'/RPA/Transitions/E2'))
        call make_dir(io_path(Params,'/RPA/Transitions/E3'))
        call make_dir(io_path(Params,'/RPA/Spectra'))
        call make_dir(io_path(Params,'/RPA/Amplitudes'))

        write(*,'(/,a)') 'Starting RPA & TDA calculation ...'
        call timer_start(t0)
        call HF_RPA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'All calculations have finished ...'
    end subroutine HF_RPA

    subroutine HF_RPA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: C, Rho_RPA, Rho_LHO
        type(o2b_t) :: V_NN
        type(orb2b_t) :: Orb_NN
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        type(phstate_t), allocatable :: Phonon(:)
        integer :: N_Phonon, J, P, J_max
        integer, allocatable :: N_nu(:,:)
        type(ivec_t), allocatable :: Orb_Phonon(:,:)
        type(tr1b_t) :: TrOp
        type(rmat_t), allocatable :: A(:,:), B(:,:), X_TDA(:,:)
        type(rvec_t), allocatable :: E_TDA(:,:), E_RPA_re(:,:)
        type(cvec_t), allocatable :: E_RPA(:,:)
        type(cmat_t), allocatable :: X_RPA(:,:), Y_RPA(:,:)
        type(reduced_multipole_t) :: rM_TDA, rM_RPA
        type(reduced_transition_t) :: rB_TDA, rB_RPA
        logical :: Stability
        real(dp) :: chR2, E_corr
        integer(i8) :: t0

        J_max = 2*Params%Calc%Nmax + 1
        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        write(*,'(/,a)') 'Importing Residual 2-body interaction ...'
        call timer_start(t0)
        call O2b_import(Params,Orb,io_path(Params,'/Bin/V2B_HF.bin'),V_NN,Orb_NN)
        call timer_report(t0)
        call orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        call orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole,N_Phonon,Phonon)

        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        call HF_RPA_phonon_count(Params,N_Phonon,Phonon,N_nu,Orb_Phonon)

        call timer_start(t0)
        call HF_RPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,V_NN,A,B)
        call timer_report(t0)

        call timer_start(t0)
        call HF_RPA_diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,E_TDA,X_TDA,E_RPA,X_RPA,Y_RPA,Stability)
        call timer_report(t0)

        Rho_RPA = HF_RPA_OBDM(Params,Orb,Orb_Phonon,Phonon,Particle,Hole,N_nu,Y_RPA)

        ! RPA charge radius & transition operators with the corresponding finite-size corrections ...
        chR2 = OBDM_chR2(Params,Orb,C,Rho_RPA)
        TrOp = Tr1b_initialize(Params,Orb,chR2)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        ! RPA OBDM in the LHO basis ...
        Rho_LHO%p = mmt(mm(C%p,Rho_RPA%p),C%p)
        Rho_LHO%n = mmt(mm(C%n,Rho_RPA%n),C%n)

        E_corr = HF_RPA_energy(Params,N_nu,E_RPA,Y_RPA)

        call HF_RPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,X_TDA,X_RPA,Y_RPA,TrOp,rM_TDA,rM_RPA)
        rB_TDA = HF_RPA_rB(Params,N_nu,rM_TDA,E_TDA)
        allocate(E_RPA_re(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                E_RPA_re(J,P)%v = real(E_RPA(J,P)%v)
            end do
        end do
        rB_RPA = HF_RPA_rB(Params,N_nu,rM_RPA,E_RPA_re)

        call HF_RPA_export(Params,Orb,N_nu,E_corr,E_TDA,E_RPA,X_TDA,X_RPA,Y_RPA,rB_TDA,rB_RPA,C,Rho_LHO,Stability)
    end subroutine HF_RPA_solver

end module nj_hf_rpa_driver
