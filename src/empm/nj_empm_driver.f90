! EMPM entry points ... first step of the Equation of Motion Phonon Method: F-matrix & TDA phonons in a restricted
    ! particle-hole space on top of a precomputed HF solution, counterpart of the Fmat & Tda programs of the original EMPM codes ...
module nj_empm_driver
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make
    use nj_o1b, only : O1b_import
    use nj_o2b, only : O2b_import
    use nj_tr1b
    use nj_empm_space
    use nj_empm_fmat
    use nj_empm_tda
    use nj_empm_export, only : EMPM_TDA_export
    use nj_empm_phonons
    use nj_empm_dens, only : empm_pairs_t, EMPM_pairs_make, phonon_block_t
    use nj_empm_eqm2
    use nj_empm_full
    use nj_empm_eqm3
    use nj_empm_cm
    use nj_linalg, only : mm, mtm, eig_sym
    !$ use omp_lib
    implicit none
    private

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
    end interface

    public :: EMPM_TDA, EMPM_2PH, EMPM_3PH, EMPM_4PH, EMPM_FULL, EMPM_TRANS, EMPM_SF

contains

    subroutine EMPM_TDA(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o1b_t) :: C
        type(o2b_t) :: V_NN
        type(orb2b_t) :: Orb_NN
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(tr1b_t) :: TrOp
        type(phstate_t), allocatable :: Phonon(:)
        integer :: N_Phonon
        integer, allocatable :: N_nu(:,:)
        type(ivec_t), allocatable :: Orb_Phonon(:,:)
        type(rmat_t), allocatable :: A(:,:), X_TDA(:,:)
        type(rvec_t), allocatable :: E_TDA(:,:)
        integer(i8) :: t0, t1

        call wigner_init_float(75,'Jmax',9)
        write(*,'(a)') 'Starting EMPM TDA calculation in the restricted particle-hole space ...'
        write(*,'(a)') trim(Params%Calc%Path)
        if (.not. is_dir(io_path(Params,''))) then
            write(*,'(/,a)') 'Error! ... No precomputed HF solution is available in given Input_File path ... run HF solver first ...'
            return
        end if
        call remove_dir(io_path(Params,'/EMPM'))
        call make_dir(io_path(Params,'/EMPM/Spectra'))
        call timer_start(t0)

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        call EMPM_space_make(Params,Orb,Space)
        call EMPM_space_export(Params,Space)

        write(*,'(/,a)') 'Importing Residual 2-body interaction ...'
        call timer_start(t1)
        call O2b_import(Params,Orb,io_path(Params,'/Bin/V2B_HF.bin'),V_NN,Orb_NN)
        call timer_report(t1)

        call timer_start(t1)
        call Fmat_make(Params,Orb,Orb_NN,V_NN,Space,F_ph)
        call timer_report(t1)
        deallocate(V_NN%pp, V_NN%pn, V_NN%nn)
        call Fmat_export(F_ph,io_path(Params,'/Bin/EMPM_F.bin'))

        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        call EMPM_TDA_phonons(Params,Space,N_Phonon,Phonon,N_nu,Orb_Phonon)
        call timer_start(t1)
        call EMPM_TDA_allocate(Params,Space,F_ph,N_nu,Orb_Phonon,Phonon,A)
        call timer_report(t1)
        call timer_start(t1)
        call EMPM_TDA_diagonalize(Params,Space,A,N_nu,Orb_Phonon,Phonon,TrOp,E_TDA,X_TDA)
        call timer_report(t1)

        call EMPM_TDA_export(Params,Space,F_ph,N_nu,Orb_Phonon,Phonon,E_TDA,X_TDA)
        call EMPM_TDA_strengths(Params,Space,N_nu,Orb_Phonon,Phonon,TrOp,E_TDA,X_TDA)
        call timer_report(t0)
        write(*,'(/,a,/)') 'EMPM TDA calculation has finished ...'
    end subroutine EMPM_TDA

    ! 2-phonon EMPM on top of EMPM_TDA ... phonons & F-matrix from Bin/EMPM_TDA_*.bin & Bin/EMPM_F.bin, densities &
        ! redefined interaction on the fly ...
    subroutine EMPM_2PH(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(empm_phonons_t) :: Ph
        type(empm_pairs_t) :: Pairs
        type(empm_basis2_t) :: Bs
        type(empm_solution2_t) :: Sol
        type(empm_full_t) :: Full
        real(dp), allocatable :: H(:,:), Dtr(:,:), Dsc(:,:), Vk(:,:), Pk(:,:)
        integer, allocatable :: Ch(:,:), A1(:)
        integer :: N0
        integer :: k, uS, uY, uZ, N2s
        integer(i8) :: t0, t1, t2

        call wigner_init_float(75,'Jmax',9)
        !$ call omp_set_max_active_levels(1)
        write(*,'(a)') 'Starting EMPM 2-phonon calculation ...'
        write(*,'(a)') trim(Params%Calc%Path)
        if (.not. is_dir(io_path(Params,'/EMPM'))) then
            write(*,'(/,a)') 'Error! ... No EMPM TDA solution is available in given Input_File path ... run EMPM_TDA first ...'
            return
        end if
        call remove_dir(io_path(Params,'/EMPM/2Phonon'))
        call make_dir(io_path(Params,'/EMPM/2Phonon'))
        call timer_start(t0)

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call EMPM_space_make(Params,Orb,Space)
        call Fmat_import(Space,io_path(Params,'/Bin/EMPM_F.bin'),F_ph)
        call EMPM_phonons_import(Params,Space,F_ph,Ph)
        write(*,'(/,a)') 'Preparing the EMPM pair space, 6j table & F blocks ...'
        call timer_start(t1)
        call EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        call timer_report(t1)
        Ch = channels(Params%Calc%EMPM%JP2,2*Ph%J_max)

        open(newunit=uS, file=io_path(Params,'/EMPM/2Phonon/Spectrum.dat'), status='replace', action='write')
        write(uS,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',20,.true.) // ' ' // 'Spurious'
        open(newunit=uY, file=io_path(Params,'/EMPM/2Phonon/Summary.dat'), status='replace', action='write')
        write(uY,'(a)') 'EMPM 2-phonon calculation ... truncation E_alpha + E_lambda <= ' // cfmt_f(Params%Calc%EMPM%E2max,0,4) // &
                        ' MeV , number of 1-phonon states = ' // istr(Ph%N)
        write(uY,'(/,a)') cfmt_s('J^P',6) // cfmt_s('N_full',12) // cfmt_s('N_tr',10) // cfmt_s('N_spur',8) // &
                          cfmt_s('N_ind',10) // cfmt_s('N_sp_ind',10) // cfmt_s('N_sp_sol',10) // cfmt_s('|AD-DA|max',14) // &
                          cfmt_s('seconds',12)

        N2s = 0
        if (Params%Calc%EMPM%Save_X) then
            open(newunit=uZ, file=io_path(Params,'/EMPM/2Phonon/States2.dat'), status='replace', action='write')
            write(uZ,'(a)') '2-phonon states ... index (channels in the order of the run), J, P, E, spurious (1) ... basis, C, X ' // &
                            '& the couplings of the spurious states in ' // states_dir(Params) // '/EMPM_S2_<JP>.bin'
            call make_dir(states_dir(Params))
            write(*,'(a)') TAB // 'State files ... ' // states_dir(Params)
        end if
        do k = 1, size(Ch,2)
            call timer_start(t1)
            call EMPM_basis2_make(Ph,Ch(1,k),Ch(2,k),Params%Calc%EMPM%E2max,Params%Calc%EMPM%E1max,Params%Calc%EMPM%Spur_Rows,Bs)
            write(*,'(/,a)') 'EMPM 2-phonon channel J^P = ' // jp_str(Bs%J,Bs%P) // ' ... full basis = ' // istr(Bs%N_full) // &
                             ' , truncated = ' // istr(Bs%N_tr) // ' , spurious = ' // istr(Bs%N_spur) // ' (rows ' // istr(Bs%N_spur_tr) // &
                             ') , all A & D rows ' // &
                             cfmt_f(2.0_dp*8.0_dp*real(Bs%N_full,dp)*real(Bs%N_tr,dp)/1.0e9_dp,0,3) // ' GB'
            if (Bs%N_tr == 0) cycle
            if (Params%Calc%EMPM%Dry_Run) then
                call export_basis()
                cycle
            end if
            if (Params%Calc%EMPM%Full_H) then
                call timer_start(t2)
                call EMPM_coupling_columns(Ph,Bs%J,Bs%P,A1,N0)
                call EMPM_coupling_V(Space,F_ph,Pairs,Ph,Bs,A1,N0,Vk)
                call timer_report(t2,'couplings V to ' // istr(size(A1)) // ' 1-phonon states' // merge(' & ground state ', &
                                  '                ', N0 == 1))
            end if
            call timer_start(t2)
            if (Params%Calc%EMPM%Full_H) then
                call EMPM_eqm2_build(Pairs,Ph,Bs,Params%Calc%EMPM%Mem_GB,Params%Calc%EMPM%Cache_GB,H,Dtr,Dsc,Vk,Pk)
                deallocate(Vk)
            else
                call EMPM_eqm2_build(Pairs,Ph,Bs,Params%Calc%EMPM%Mem_GB,Params%Calc%EMPM%Cache_GB,H,Dtr,Dsc)
            end if
            call timer_report(t2,'A D^T & D')
            call timer_start(t2)
            call EMPM_eqm2_solve(Bs,H,Dtr,Dsc,Params%Calc%EMPM%Legacy_Spur,Params%Calc%EMPM%D_Spectrum,Sol)
            call timer_report(t2,'eigenproblem')
            deallocate(H, Dtr, Dsc)
            write(*,'(a)') TAB // 'Independent states = ' // istr(Sol%N_ind) // ' , independent spurious states = ' // &
                           istr(Sol%N_spur_ind) // ' , spurious solutions = ' // istr(Sol%N_s) // ' , |AD - (AD)^T|max = ' // &
                           cfmt_e(Sol%Asym,0,3)
            call export_channel()
            if (Params%Calc%EMPM%Save_X) call export_states()
            if (Params%Calc%EMPM%Full_H) then
                call timer_start(t2)
                call EMPM_full_solve(Ph,Bs,Sol,A1,N0,Pk,Full)
                call timer_report(t2,'full Hamiltonian (dimension ' // istr(Full%N) // ')')
                call export_full()
                deallocate(Pk)
            end if
        end do
        close(uS)
        close(uY)
        if (Params%Calc%EMPM%Save_X) close(uZ)
        call timer_report(t0)
        write(*,'(/,a,/)') 'EMPM 2-phonon calculation has finished ...'

    contains

        ! 2-phonon states of the channel for the 3-phonon step ... basis, E, C on the independent rows (Kc = their full-basis
            ! indices) &
            ! X = sqrt(2J+1) D(independent rows)^T C in Bin/EMPM_S2_<JP>.bin, the list in EMPM/2Phonon/States2.dat ... the
            ! physical states & (Ortho) the spurious ones with their couplings R, a complete set of the 2-phonon space ...
        subroutine export_states()
            real(dp), allocatable :: Dr(:,:), X2(:,:)
            integer :: n2, i
            integer(i8) :: t3
            call timer_start(t3)
            n2 = Sol%N_ind
            allocate(Dr(Bs%N_full,Sol%N_ind), X2(Bs%N_full,n2))
            call EMPM_AD_rows(Pairs,Ph,Bs,Sol%Ind,Params%Calc%EMPM%Cache_GB,Dr)
            if (n2 > 0) call dgemm('N','N',Bs%N_full,n2,Sol%N_ind,sqrt(real(2*Bs%J + 1,dp)),Dr,Bs%N_full,Sol%C,Sol%N_ind, &
                                   0.0_dp,X2,Bs%N_full)
            call EMPM_states_write(states_file(Params,2,jp_str(Bs%J,Bs%P)),Bs,Sol%E(1:n2),Bs%Row(Sol%Ind),Sol%C,X2, &
                                   Sol%N_ind-Sol%N_s,Sol%R)
            do i = 1, n2
                write(uZ,'(a)') cfmt_d(N2s+i,8) // cfmt_d(Bs%J,5) // cfmt_s(merge('+','-',Bs%P == 1),4) // cfmt_f(Sol%E(i),16,6) // &
                                cfmt_d(merge(1, 0, i > Sol%N_ind - Sol%N_s),4)
            end do
            N2s = N2s + n2
            call timer_report(t3,'X & C of the 2-phonon states exported (' // istr(n2) // ' states, ' // &
                              cfmt_f(16.0_dp*real(Bs%N_full,dp)*real(n2,dp)/1.0e9_dp,0,3) // ' GB)')
        end subroutine export_states

        subroutine export_channel()
            integer :: u, i, n_phys
            integer(i8) :: t3, Rate
            character(len=:), allocatable :: Tag

            Tag = jp_str(Bs%J,Bs%P)
            n_phys = Sol%N_ind - Sol%N_s
            do i = 1, Sol%N_ind
                write(uS,'(a)') cfmt_d(Bs%J,5,.true.) // ' ' // cfmt_s(merge('+','-',Bs%P == 1),5,.true.) // ' ' // &
                                cfmt_f(Sol%E(i),20,6,.true.) // ' ' // istr(merge(1, 0, i > n_phys))
            end do
            write(uS,'(a)') ''
            write(uS,'(a)') ''
            call system_clock(t3,Rate)
            write(uY,'(a)') cfmt_s(Tag,6) // cfmt_d(Bs%N_full,12) // cfmt_d(Bs%N_tr,10) // cfmt_d(Bs%N_spur,8) // &
                            cfmt_d(Sol%N_ind,10) // cfmt_d(Sol%N_spur_ind,10) // cfmt_d(Sol%N_s,10) // cfmt_e(Sol%Asym,14,3) // &
                            cfmt_f(real(t3 - t1,dp)/real(Rate,dp),12,2)

            call export_basis(Sol%Ind)

            ! Energies & amplitudes C(independent row, solution) ... the independent rows are the Ind = 1 lines of Basis_*.dat ...
            open(newunit=u, file=io_path(Params,'/Bin/EMPM_2PH_' // Tag // '_E.bin'), access='stream', form='unformatted', &
                 status='replace', action='write')
            write(u) Sol%E
            close(u)
            open(newunit=u, file=io_path(Params,'/Bin/EMPM_2PH_' // Tag // '_C.bin'), access='stream', form='unformatted', &
                 status='replace', action='write')
            write(u) Sol%C
            close(u)
        end subroutine export_channel

        ! Full Hamiltonian ... spectrum with the 0-, 1- & 2-phonon weights & the eigenvectors ...
        subroutine export_full()
            integer :: u, i
            character(len=:), allocatable :: Tag
            Tag = jp_str(Bs%J,Bs%P)
            open(newunit=u, file=io_path(Params,'/EMPM/2Phonon/Full_' // Tag // '.dat'), status='replace', action='write')
            write(u,'(a)') 'EMPM full Hamiltonian J^P = ' // Tag // ' ... basis: ' // istr(Full%N0) // ' ground state, ' // &
                           istr(Full%N1) // ' 1-phonon states (phonon indices in Full_' // Tag // '_1phon.dat), ' // &
                           istr(Full%N2) // ' physical 2-phonon states'
            write(u,'(/,a)') cfmt_s('i',8) // cfmt_s('E',16) // cfmt_s('w0',12) // cfmt_s('w1',12) // cfmt_s('w2',12)
            do i = 1, Full%N
                write(u,'(a)') cfmt_d(i,8) // cfmt_f(Full%E(i),16,6) // cfmt_f(sum(Full%Z(1:Full%N0,i)**2),12,6) // &
                               cfmt_f(sum(Full%Z(Full%N0+1:Full%N0+Full%N1,i)**2),12,6) // &
                               cfmt_f(sum(Full%Z(Full%N0+Full%N1+1:Full%N,i)**2),12,6)
            end do
            close(u)
            open(newunit=u, file=io_path(Params,'/EMPM/2Phonon/Full_' // Tag // '_1phon.dat'), status='replace', action='write')
            write(u,'(a)') 'Position in the full basis & 1-phonon index (as in Spectrum numbering of the phonons) & energy'
            do i = 1, Full%N1
                write(u,'(a)') cfmt_d(Full%N0+i,8) // cfmt_d(Full%A1(i),8) // cfmt_f(Ph%E(Full%A1(i)),16,6)
            end do
            close(u)
            open(newunit=u, file=io_path(Params,'/Bin/EMPM_FULL_' // Tag // '_E.bin'), access='stream', form='unformatted', &
                 status='replace', action='write')
            write(u) Full%E
            close(u)
            open(newunit=u, file=io_path(Params,'/Bin/EMPM_FULL_' // Tag // '_Z.bin'), access='stream', form='unformatted', &
                 status='replace', action='write')
            write(u) Full%Z
            close(u)
        end subroutine export_full

        ! Truncated basis with spurious & (when known) independence flags ...
        subroutine export_basis(Ind)
            integer, intent(in), optional :: Ind(:)
            integer :: u, i
            character(len=:), allocatable :: Tag
            logical, allocatable :: Is_Ind(:)

            Tag = jp_str(Bs%J,Bs%P)
            allocate(Is_Ind(Bs%N_tr))
            Is_Ind = .false.
            if (present(Ind)) Is_Ind(Ind) = .true.
            open(newunit=u, file=io_path(Params,'/EMPM/2Phonon/Basis_' // Tag // '.dat'), status='replace', action='write')
            write(u,'(a)') 'EMPM 2-phonon truncated basis |alpha lambda; J^P = ' // Tag // '> ... i truncated index, ' // &
                           'k full-basis index, phonon indices ' // merge('of the original Tda (Legacy_Order)  ', &
                                                                         'as in EMPM/TDA_Structure.dat        ', &
                                                                         Params%Calc%EMPM%Legacy_Order)
            write(u,'(/,a)') cfmt_s('i',8) // cfmt_s('k',10) // cfmt_s('Spur',6) // cfmt_s('Ind',5) // cfmt_s('lambda',8) // &
                             cfmt_s('E_lambda',12) // cfmt_s('alpha',8) // cfmt_s('E_alpha',12)
            do i = 1, Bs%N_tr
                write(u,'(a)') cfmt_d(i,8) // cfmt_d(Bs%Row(i),10) // cfmt_d(merge(1,0,Bs%Spur(Bs%Row(i))),6) // &
                               cfmt_d(merge(1,0,Is_Ind(i)),5) // cfmt_d(Bs%La(Bs%Row(i)),8) // &
                               cfmt_f(Ph%E(Bs%La(Bs%Row(i))),12,5) // cfmt_d(Bs%Al(Bs%Row(i)),8) // cfmt_f(Ph%E(Bs%Al(Bs%Row(i))),12,5)
            end do
            close(u)
        end subroutine export_basis

    end subroutine EMPM_2PH

    ! 3-phonon EMPM on top of EMPM_2PH (Save_X) ... the 2-phonon states of all channels (basis, E, C & X; with Ortho the
        ! physical & spurious states with their couplings Rc, a complete set) from EMPM_S2_*.bin, the 2-phonon densities
        ! computed on the fly ... with Ortho the spurious 3-phonon subspace (cm phonon x 2-phonon states & phonons x
        ! spurious 2-phonon states) is projected out as in the 2-phonon step (spur_gev) ...
    subroutine EMPM_3PH(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(empm_phonons_t) :: Ph
        type(empm_pairs_t) :: Pairs
        type(empm_states2_t) :: S1, S2
        type(empm_cache1_t) :: C1
        integer :: uY
        integer(i8) :: t0, t1

        call wigner_init_float(75,'Jmax',9)
        !$ call omp_set_max_active_levels(1)
        write(*,'(a)') 'Starting EMPM 3-phonon calculation ...'
        write(*,'(a)') trim(Params%Calc%Path)
        if (.not. is_dir(io_path(Params,'/EMPM'))) then
            write(*,'(/,a)') 'Error! ... No EMPM TDA solution is available in given Input_File path ... run EMPM_TDA first ...'
            return
        end if
        call remove_dir(io_path(Params,'/EMPM/3Phonon'))
        call make_dir(io_path(Params,'/EMPM/3Phonon'))
        call timer_start(t0)

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call EMPM_space_make(Params,Orb,Space)
        call Fmat_import(Space,io_path(Params,'/Bin/EMPM_F.bin'),F_ph)
        call EMPM_phonons_import(Params,Space,F_ph,Ph)
        write(*,'(/,a)') 'Preparing the EMPM pair space, 6j table & F blocks ...'
        call timer_start(t1)
        call EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        call timer_report(t1)

        write(*,'(/,a)') 'Importing the 2-phonon states of EMPM_2PH (Bin/EMPM_S2_*.bin) ...'
        call EMPM_states_read(states_prefix(Params,2),2*Ph%J_max,Ph%N,Ph%N,S2)
        if (S2%N == 0) call fatal('EMPM_3PH ... no 2-phonon states in Bin/EMPM_S2_*.bin, run EMPM_2PH with Save_X = .true. first')
        write(*,'(a)') TAB // 'Number of 2-phonon states = ' // istr(S2%N) // ' , physical ' // istr(count(.not. S2%Spur)) // &
                       ' , spurious ' // istr(count(S2%Spur))
        open(newunit=uY, file=io_path(Params,'/EMPM/3Phonon/Summary.dat'), status='replace', action='write')
        write(uY,'(a)') 'EMPM 3-phonon calculation ... ' // istr(S2%N) // ' 2-phonon states from EMPM_2PH , 3-phonon rows ' // &
                        'E_a + E_lambda <= ' // cfmt_f(Params%Calc%EMPM%E3max,0,4) // ' MeV , number of 1-phonon states = ' // &
                        istr(Ph%N)

        write(*,'(/,a)') 'EMPM 1-phonon densities & redefined interactions of all phonons ...'
        call timer_start(t1)
        call EMPM_cache1_make(Pairs,Ph,C1)
        call timer_report(t1)

        call EMPM_states1_make(Ph,S1)
        call nph_channels(Params,3,Pairs,Ph,C1,C1%R,S1,S2,Params%Calc%EMPM%E3max,Params%Calc%EMPM%JP3,uY,Space)
        close(uY)
        call timer_report(t0)
        write(*,'(/,a,/)') 'EMPM 3-phonon calculation has finished ...'
    end subroutine EMPM_3PH

    ! Channels of the n-phonon step (n = 3, 4) ... basis |a lambda> of the (n-1)-phonon states S & the phonons, A & D with
        ! the densities of S computed from the lower level (Rlow, Slow), eigenproblem, spectrum & summary (unit uY), and the
        ! n-phonon states (basis, E, C & X) in Bin/EMPM_S<n>_<JP>.bin for the (n+1)-phonon step ...
    subroutine nph_channels(Params,Nph,Pairs,Ph,C1,Rlow,Slow,S,Emax,JP,uY,Space)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        integer, intent(in) :: Nph, uY
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_cache1_t), intent(in) :: C1
        type(phonon_block_t), intent(in) :: Rlow(:)
        type(empm_states2_t), intent(in) :: Slow, S
        real(dp), intent(in) :: Emax
        character(len=*), intent(in) :: JP(:)
        type(empm_basis2_t) :: B
        type(empm_solution2_t) :: Sol
        real(dp), allocatable :: H(:,:), Dtr(:,:), Dsc(:,:), Dr(:,:), Xn(:,:), Dii(:,:), Dsp(:,:)
        integer, allocatable :: Ch(:,:)
        integer :: k, i, uS, u, nn, Ns
        integer(i8) :: t1, t2, t3, Rate
        character(len=:), allocatable :: Tag, Dir, Sn
        integer(i8), allocatable :: Dnp(:,:)
        integer :: N_visit, Tgt
        integer, allocatable :: Ord_D(:)
        character(len=:), allocatable :: D_File
        real(dp) :: D_GB, D_Lim
        type(greedy_res_t), allocatable :: GR(:)
        type(empm_basis2_t), allocatable :: Bsh(:)
        integer, allocatable :: Kmap(:), Tsh(:)
        logical, allocatable :: Dsk(:)
        character(len=1024), allocatable :: Fsh(:)
        integer :: nsh

        Sn = istr(Nph)
        allocate(Dnp(0:Nph*Ph%J_max+1,2))
        Dnp(:,:) = EMPM_npnh_dims(Space,Nph,Nph*Ph%J_max + 1)
        call make_dir(states_dir(Params))
        write(*,'(a)') TAB // 'State files ... ' // states_dir(Params)
        Dir = '/EMPM/' // Sn // 'Phonon'
        Ch = channels(JP,S%J_max + Ph%J_max)
        open(newunit=uS, file=io_path(Params,Dir // '/Spectrum.dat'), status='replace', action='write')
        write(uS,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',20,.true.) // ' ' // 'Spurious'
        write(uY,'(/,a)') Sn // '-phonon channels ...'
        write(uY,'(a)') cfmt_s('J^P',6) // cfmt_s('N_full',12) // cfmt_s('N_tr',10) // cfmt_s('N_spur',8) // cfmt_s('N_ind',10) // &
                        cfmt_s('N_sp_sol',10) // cfmt_s('|AD-DA|max',14) // cfmt_s('seconds',12)
        open(newunit=u, file=io_path(Params,Dir // '/States' // Sn // '.dat'), status='replace', action='write')
        write(u,'(a)') Sn // '-phonon states ... index (channels in the order of the run), J, P, E, spurious (1) ... C, X & ' // &
                       'the couplings of the spurious states in ' // states_dir(Params) // '/EMPM_S' // Sn // '_<JP>.bin'
        Ns = 0

        ! Dens_Share ... the greedy selections of all channels first, in one pass with shared densities ...
        allocate(Kmap(size(Ch,2)))
        Kmap = 0
        if (Params%Calc%EMPM%Dens_Share .and. Params%Calc%EMPM%Ind_First .and. Params%Calc%EMPM%Ind_Greedy .and. &
            .not. Params%Calc%EMPM%Dry_Run) then
            allocate(Bsh(size(Ch,2)), Tsh(size(Ch,2)), Dsk(size(Ch,2)), Fsh(size(Ch,2)))
            nsh = 0
            do k = 1, size(Ch,2)
                call EMPM_basis3_make(Ph,S,Ch(1,k),Ch(2,k),Emax,Params%Calc%EMPM%Spur_Rows,B)
                if (B%N_tr == 0) cycle
                nsh = nsh + 1
                Kmap(k) = nsh
                Tgt = B%N_tr
                if (B%J <= ubound(Dnp,1)) Tgt = int(min(int(B%N_tr,i8), Dnp(B%J,B%P)))
                Tsh(nsh) = Tgt
                D_GB = 8.0_dp*real(B%N_full,dp)*real(Tgt,dp)/1.0e9_dp
                D_Lim = Params%Calc%EMPM%D_Disk_GB
                if (D_Lim < 0.0_dp) D_Lim = 0.3_dp*Params%Calc%EMPM%Mem_GB/real(size(Ch,2),dp)
                Dsk(nsh) = D_GB > D_Lim
                Fsh(nsh) = states_file(Params,Nph,jp_str(B%J,B%P)) // '.D.tmp'
                Bsh(nsh) = B
            end do
            write(*,'(/,a)') 'EMPM ' // Sn // '-phonon greedy selection of ' // istr(nsh) // ' channels with shared densities ...'
            allocate(GR(nsh))
            call timer_start(t2)
            call EMPM_ind_greedy_multi(Pairs,Ph,C1,Rlow,Slow,S,Bsh(1:nsh),Tsh(1:nsh),Params%Calc%EMPM%Mem_GB,Dsk(1:nsh), &
                                       Fsh(1:nsh),GR)
            call timer_report(t2,'shared greedy selection & A D^T of ' // istr(nsh) // ' channels')
            deallocate(Bsh)
        end if

        do k = 1, size(Ch,2)
            call timer_start(t1)
            call EMPM_basis3_make(Ph,S,Ch(1,k),Ch(2,k),Emax,Params%Calc%EMPM%Spur_Rows,B)
            Tag = jp_str(B%J,B%P)
            write(*,'(/,a)') 'EMPM ' // Sn // '-phonon channel J^P = ' // Tag // ' ... full basis = ' // istr(B%N_full) // &
                             ' , truncated = ' // istr(B%N_tr) // ' , spurious = ' // istr(B%N_spur) // ' (rows ' // &
                             istr(B%N_spur_tr) // ') , all A & D rows ' // &
                             cfmt_f(2.0_dp*8.0_dp*real(B%N_full,dp)*real(B%N_tr,dp)/1.0e9_dp,0,3) // ' GB'
            if (B%N_tr == 0 .or. Params%Calc%EMPM%Dry_Run) cycle
            if (Params%Calc%EMPM%Ind_First) then
                ! D of all truncated rows -> independent rows (pivoted Cholesky) -> A, D & H = A D^T on them only ...
                Sol = empm_solution2_t()
                Tgt = B%N_tr
                if (B%J <= ubound(Dnp,1)) Tgt = int(min(int(B%N_tr,i8), Dnp(B%J,B%P)))
                if (Kmap(k) > 0) then
                    ! from the shared pass ...
                    i = Kmap(k)
                    call move_alloc(GR(i)%Ind,Sol%Ind)
                    call move_alloc(GR(i)%Dii,Dii)
                    call move_alloc(GR(i)%H,H)
                    if (allocated(GR(i)%Dr)) call move_alloc(GR(i)%Dr,Dr)
                    if (allocated(GR(i)%Ord)) then
                        call move_alloc(GR(i)%Ord,Ord_D)
                        call move_alloc(GR(i)%Dsp,Dsp)
                        D_File = trim(Fsh(i))
                    end if
                    Sol%N_ind = size(Sol%Ind)
                    write(*,'(a)') TAB // 'independent rows (shared greedy) ... ' // istr(Sol%N_ind) // ' of ' // istr(B%N_tr) // &
                                   ' , npnh dimension ' // istr(Tsh(i)) // ' , rows visited ' // istr(GR(i)%N_visit)
                    if (Sol%N_ind /= Tsh(i)) write(*,'(a)') TAB // 'Warning ... independent rows /= npnh dimension'
                else if (Params%Calc%EMPM%Ind_Greedy) then
                    ! block-wise selection up to the npnh dimension, no N_tr x N_tr D ...
                    ! (A rows in the same pass, H = A D^T accumulated, only D of the selected rows kept) ...
                    call timer_start(t2)
                    ! D of the selected rows on disk when larger than D_Disk_GB (automatic: 0.3 Mem_GB) ...
                    D_GB = 8.0_dp*real(B%N_full,dp)*real(Tgt,dp)/1.0e9_dp
                    D_Lim = Params%Calc%EMPM%D_Disk_GB
                    if (D_Lim < 0.0_dp) D_Lim = 0.3_dp*Params%Calc%EMPM%Mem_GB
                    if (D_GB > D_Lim) then
                        D_File = states_file(Params,Nph,Tag) // '.D.tmp'
                        call EMPM_ind_greedy_H(Pairs,Ph,C1,Rlow,Slow,S,B,Tgt,Params%Calc%EMPM%Mem_GB,Sol%Ind,Dii,H,Dr,N_visit, &
                                               D_File,Ord_D,Dsp)
                    else
                        call EMPM_ind_greedy_H(Pairs,Ph,C1,Rlow,Slow,S,B,Tgt,Params%Calc%EMPM%Mem_GB,Sol%Ind,Dii,H,Dr,N_visit)
                    end if
                    Sol%N_ind = size(Sol%Ind)
                    call timer_report(t2,'independent rows (greedy) & A D^T ... ' // istr(Sol%N_ind) // ' of ' // istr(B%N_tr) // &
                                      ' , npnh dimension ' // istr(Tgt) // ' , rows visited ' // istr(N_visit))
                    if (Sol%N_ind /= Tgt) write(*,'(a)') TAB // 'Warning ... independent rows /= npnh dimension'
                else
                call timer_start(t2)
                call EMPM_D3_build(Pairs,Ph,C1,Rlow,Slow,S,B,Params%Calc%EMPM%Mem_GB,Dtr)
                call timer_report(t2,'D of the truncated rows')
                call timer_start(t2)
                call EMPM_chol_ind(Dtr,Sol%Ind,Dii)
                deallocate(Dtr)
                Sol%N_ind = size(Sol%Ind)
                call timer_report(t2,'independent rows (pivoted Cholesky) ... ' // istr(Sol%N_ind) // ' of ' // istr(B%N_tr) // &
                                  ' , npnh dimension ' // istr(Tgt))
                end if
                if (.not. Params%Calc%EMPM%Ind_Greedy) then
                    call timer_start(t2)
                    call EMPM_H3_ind(Pairs,Ph,C1,Rlow,Slow,S,B,Sol%Ind,Params%Calc%EMPM%Mem_GB,H,Dr)
                    call timer_report(t2,'A D^T of the independent rows')
                end if
                call timer_start(t2)
                if (B%N_spur > 0) then
                    ! D elements of the spurious states & the independent rows (full-basis D rows Dr) ...
                    if (.not. allocated(Ord_D)) Dsp = Dr(B%Ks,:)
                    call EMPM_ind_solve(H,Dii,Sol,Dsp)
                    deallocate(Dsp)
                else
                    call EMPM_ind_solve(H,Dii,Sol)
                end if
                call timer_report(t2,'eigenproblem')
            else
                call timer_start(t2)
                call EMPM_eqm3_build(Pairs,Ph,C1,Rlow,Slow,S,B,Params%Calc%EMPM%Mem_GB,H,Dtr,Dsc)
                call timer_report(t2,'A D^T & D')
                call timer_start(t2)
                call EMPM_eqm2_solve(B,H,Dtr,Dsc,.false.,Params%Calc%EMPM%D_Spectrum,Sol)
                call timer_report(t2,'eigenproblem')
                deallocate(H, Dtr, Dsc)
            end if
            write(*,'(a)') TAB // 'Independent states = ' // istr(Sol%N_ind) // ' , spurious solutions = ' // istr(Sol%N_s) // &
                           ' , |AD - (AD)^T|max = ' // cfmt_e(Sol%Asym,0,3)
            do i = 1, Sol%N_ind
                write(uS,'(a)') cfmt_d(B%J,5,.true.) // ' ' // cfmt_s(merge('+','-',B%P == 1),5,.true.) // ' ' // &
                                cfmt_f(Sol%E(i),20,6,.true.) // ' ' // istr(merge(1, 0, i > Sol%N_ind - Sol%N_s))
            end do
            write(uS,'(a)') ''
            write(uS,'(a)') ''
            call system_clock(t3,Rate)
            write(uY,'(a)') cfmt_s(Tag,6) // cfmt_d(B%N_full,12) // cfmt_d(B%N_tr,10) // cfmt_d(B%N_spur,8) // cfmt_d(Sol%N_ind,10) // &
                            cfmt_d(Sol%N_s,10) // cfmt_e(Sol%Asym,14,3) // cfmt_f(real(t3 - t1,dp)/real(Rate,dp),12,2)
            flush(uS)
            flush(uY)
            ! n-phonon states for the (n+1)-phonon step ... X = sqrt(2J+1) D(independent rows)^T C & C on the independent rows
                ! (Kc = their full-basis indices) ...
                ! the physical & the spurious states with their couplings R (a complete set) ...
            call timer_start(t2)
            nn = Sol%N_ind
            if (allocated(Ord_D)) then
                ! X = sqrt(2J+1) Dr C from the D rows on disk, in blocks of basis rows ...
                call EMPM_states_write_Df(states_file(Params,Nph,Tag),B,Sol%E(1:nn),B%Row(Sol%Ind),Sol%C,D_File,Ord_D, &
                                          sqrt(real(2*B%J + 1,dp)),Sol%N_ind-Sol%N_s,Sol%R,Params%Calc%EMPM%Mem_GB)
                deallocate(Ord_D)
            else
            if (.not. allocated(Dr)) then
                allocate(Dr(B%N_full,Sol%N_ind))
                call EMPM_AD3_rows(Pairs,Ph,C1,Rlow,Slow,S,B,Sol%Ind,Dr)
            end if
            ! X = sqrt(2J+1) Dr C written in column blocks (no full X in memory) ...
            call EMPM_states_write_D(states_file(Params,Nph,Tag),B,Sol%E(1:nn),B%Row(Sol%Ind),Sol%C,Dr, &
                                     sqrt(real(2*B%J + 1,dp)),Sol%N_ind-Sol%N_s,Sol%R)
            deallocate(Dr)
            end if
            do i = 1, nn
                write(u,'(a)') cfmt_d(Ns+i,8) // cfmt_d(B%J,5) // cfmt_s(merge('+','-',B%P == 1),4) // cfmt_f(Sol%E(i),16,6) // &
                               cfmt_d(merge(1, 0, i > Sol%N_ind - Sol%N_s),4)
            end do
            Ns = Ns + nn
            call timer_report(t2,'X & C of the ' // Sn // '-phonon states exported (' // istr(nn) // ' states, ' // &
                              cfmt_f(8.0_dp*(real(B%N_full,dp) + real(Sol%N_ind,dp))*real(nn,dp)/1.0e9_dp,0,3) // ' GB)')
        end do
        close(u)
        close(uS)
    end subroutine nph_channels

    ! 4-phonon EMPM on top of a complete EMPM_3PH run ... the 2- & 3-phonon states (basis, E, C & X) from Bin/EMPM_S2_*.bin &
        ! Bin/EMPM_S3_*.bin, the 2-phonon densities of all 2-phonon states cached, the 3-phonon densities on the fly ...
        ! basis |c lambda> of the 3-phonon states c & the phonons, rows E_c + E_lambda <= E4max, channels JP4 ...
    subroutine EMPM_4PH(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(empm_phonons_t) :: Ph
        type(empm_pairs_t) :: Pairs
        type(empm_states2_t) :: S1, S2, S3
        type(empm_cache1_t) :: C1
        type(phonon_block_t), allocatable :: R2(:)
        integer :: a, uY
        integer(i8) :: t0, t1

        call wigner_init_float(75,'Jmax',9)
        !$ call omp_set_max_active_levels(1)
        write(*,'(a)') 'Starting EMPM 4-phonon calculation ...'
        write(*,'(a)') trim(Params%Calc%Path)
        if (.not. is_dir(io_path(Params,'/EMPM'))) then
            write(*,'(/,a)') 'Error! ... No EMPM TDA solution is available in given Input_File path ... run EMPM_TDA first ...'
            return
        end if
        call remove_dir(io_path(Params,'/EMPM/4Phonon'))
        call make_dir(io_path(Params,'/EMPM/4Phonon'))
        call timer_start(t0)

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call EMPM_space_make(Params,Orb,Space)
        call Fmat_import(Space,io_path(Params,'/Bin/EMPM_F.bin'),F_ph)
        call EMPM_phonons_import(Params,Space,F_ph,Ph)
        write(*,'(/,a)') 'Preparing the EMPM pair space, 6j table & F blocks ...'
        call timer_start(t1)
        call EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        call timer_report(t1)

        write(*,'(/,a)') 'Importing the 2- & 3-phonon states of EMPM_3PH (Bin/EMPM_S2_*.bin & Bin/EMPM_S3_*.bin) ...'
        call EMPM_states_read(states_prefix(Params,2),2*Ph%J_max,Ph%N,Ph%N,S2)
        if (S2%N == 0) call fatal('EMPM_4PH ... no 2-phonon states in Bin/EMPM_S2_*.bin, run EMPM_3PH first')
        call EMPM_states_read(states_prefix(Params,3),3*Ph%J_max,S2%N,Ph%N,S3)
        if (S3%N == 0) call fatal('EMPM_4PH ... no 3-phonon states in Bin/EMPM_S3_*.bin, run EMPM_3PH first')
        write(*,'(a)') TAB // '2-phonon states = ' // istr(S2%N) // ' (spurious ' // istr(count(S2%Spur)) // ') , 3-phonon ' // &
                       'states = ' // istr(S3%N) // ' (spurious ' // istr(count(S3%Spur)) // ')'

        open(newunit=uY, file=io_path(Params,'/EMPM/4Phonon/Summary.dat'), status='replace', action='write')
        write(uY,'(a)') 'EMPM 4-phonon calculation ... ' // istr(S2%N) // ' 2-phonon & ' // istr(S3%N) // &
                        ' 3-phonon states from EMPM_3PH , 4-phonon rows E_c + E_lambda <= ' // &
                        cfmt_f(Params%Calc%EMPM%E4max,0,4) // ' MeV , number of 1-phonon states = ' // istr(Ph%N)

        write(*,'(/,a)') 'EMPM 1-phonon densities & redefined interactions of all phonons, 2-phonon densities of all ' // &
                         '2-phonon states ...'
        call timer_start(t1)
        call EMPM_cache1_make(Pairs,Ph,C1)
        call EMPM_states1_make(Ph,S1)
        allocate(R2(S2%N))
        !$omp parallel do schedule(dynamic)
        do a = 1, S2%N
            call EMPM_densityN(Pairs,Ph,C1,C1%R,S1,S2,a,R2(a))
        end do
        !$omp end parallel do
        call timer_report(t1)

        call nph_channels(Params,4,Pairs,Ph,C1,R2,S2,S3,Params%Calc%EMPM%E4max,Params%Calc%EMPM%JP4,uY,Space)
        close(uY)
        call timer_report(t0)
        write(*,'(/,a,/)') 'EMPM 4-phonon calculation has finished ...'
    end subroutine EMPM_4PH

    ! Full EMPM Hamiltonian per J^P channel, diagonalized sequentially with up to 1, 2, .. Full_N (<= 4) phonons, on top of
        ! complete EMPM_2PH (Save_X), EMPM_3PH & EMPM_4PH runs ... basis [ground (0+) ; 1-phonon ; 2- ; 3- ; 4-phonon states],
        ! H built once up to Full_N, its leading [0 ; .. ; n] block diagonalized for n = 1 .. Full_N, each with its own outputs
        ! in EMPM/Full/<n>Phonon & Bin/EMPM_FULL<n>_<JP>_*.bin ... diagonal the energies of each subspace, couplings
        ! n+1 <- n: H = X^T V with V(alpha lambda; a) = (-1)^(Ja+J_alpha+J_lambda)/(2Ja+1) sum_s rho_(a alpha)(J_lambda; s) W(s,lambda)
        ! (n-phonon densities of a, as vintn12 & vintn23) ... n+2 <- n: only the pair-creation part V20 of H, which commutes with
        ! the phonon operators, G_0 = h20 = <gamma|H|0> of the 0+ 2-phonon states, G_n(beta; alpha) = 1/sqrt(2J+1)
        ! sum_(k = mu nu) C^alpha_k sum_d G_(n-1)(d; mu) X^beta(d nu) (n = 1 the h31 of nondiag13) ...
    subroutine EMPM_FULL(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(empm_phonons_t) :: Ph
        type(empm_pairs_t) :: Pairs
        type(empm_states2_t) :: S1, S2, S3, S4
        type(empm_cache1_t) :: C1
        type(phonon_block_t), allocatable :: R2(:)
        type(phonon_block_t) :: Rho
        type(rmat_t), allocatable :: G1(:,:)
        real(dp), allocatable :: W(:,:), Vk(:,:), H20(:), Hf(:,:), Hs(:,:), Ef(:), Vc(:,:), Blk(:,:), Uv(:)
        real(dp), allocatable :: Mcm(:,:), Bph(:,:), Bhp(:,:), Bd(:,:), Bt(:,:), Ncm(:), Ucm(:,:), Ncm0(:), Hp(:,:), Ml(:,:,:)
        type(rmat_t) :: Dl(3), Dc(3)
        type(empm_cm_t) :: Cm
        real(dp), allocatable :: Mt(:,:)
        character(len=8) :: Chk
        logical :: CMx, CMv, CMp, CMn, CMa, CMl, Keep, DoU
        character(len=8) :: Fd
        character(len=:), allocatable :: LabU
        real(dp), allocatable :: Ur(:,:)
        integer :: Jt, Pt, Nbc(0:4), Mbc(0:4), Np, uC, q, iv, kr, nph, ii
        integer, allocatable :: Ch(:,:), A1(:), I1(:)
        integer :: Nf, J, P, k, i, a, n0, n1, n2, n3, n4, o1, o2, o3, o4, N, u, J_all, kk, mu, nu, d, jm, pm, ia, nd, nn, Ns, ic
        integer :: uS(4), uP(4), Nb(0:4)
        integer(i8) :: t0, t1
        character(len=:), allocatable :: Tag
        real(dp) :: wt(0:4)

        call wigner_init_float(75,'Jmax',9)
        !$ call omp_set_max_active_levels(1)
        Nf = Params%Calc%EMPM%Full_N
        if (Nf < 1 .or. Nf > 4) call fatal('EMPM_FULL ... Full_N must be 1, 2, 3 or 4')
        write(*,'(a)') 'Starting the full EMPM Hamiltonian, diagonalized with up to 1 .. ' // istr(Nf) // ' phonons ...'
        write(*,'(a)') trim(Params%Calc%Path)
        call remove_dir(io_path(Params,'/EMPM/Full'))
        do nn = 1, Nf
            call make_dir(io_path(Params,'/EMPM/Full/' // istr(nn) // 'Phonon'))
        end do
        call timer_start(t0)

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call EMPM_space_make(Params,Orb,Space)
        call Fmat_import(Space,io_path(Params,'/Bin/EMPM_F.bin'),F_ph)
        call EMPM_phonons_import(Params,Space,F_ph,Ph)
        ! Ortho ... all spurious states are left out of the basis: the spurious 1- phonon i_cm (the normalized E1 vector of
            ! EMPM_TDA) & the spurious n-phonon states (the last ones of each channel, n_phys = N_phys(J,P) physical) ... the
            ! couplings sum over the complete sets (physical & spurious) of the intermediate states ...
        ! CM_Keep ... with the Ortho = .true. sets all states stay in the basis (the CM-clean phonons, the cm phonon & the
            ! spurious n-phonon states), the within-level couplings R(nu) = <nu|A|cm> & Rc = <physical|H|spurious> filled in,
            ! so that the CM conditions (CM_Exact, CM_Variants) act on the complete sets in this basis ...
        Keep = Params%Calc%EMPM%CM_Keep .and. Ph%i_cm > 0
        if (Ph%i_cm > 0 .and. .not. Keep) write(*,'(a)') TAB // 'Spurious 1- phonon i = ' // istr(Ph%i_cm) // &
                                        ' & the spurious n-phonon states left out of the basis'
        if (Keep) write(*,'(a)') TAB // 'CM_Keep ... the spurious 1- phonon i = ' // istr(Ph%i_cm) // ' & the spurious ' // &
                                 'n-phonon states kept in the basis with their couplings'
        call EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        ! Diagonalizations (Full_Diag): DoU the full H without CM projection ('unproj' / 'both'), CMp the CM projection
            ! 'proj' ('proj' / 'both', or the diagnostics CM_Exact & CM_Variants) ... CMx exact spectrum (CM_Exact), CMv 'space'
            ! (CM_Exact & CM_Variants), CMn the complete N_cm matrix (exact <N_cm> of all states): CM_Ncm with complete sets,
            ! always with CM_Exact & CM_Beta ... without CMn the within-level b_cm of the top level Full_N is not needed ...
        Fd = Params%Calc%EMPM%Full_Diag
        if (Fd /= 'unproj' .and. Fd /= 'proj' .and. Fd /= 'both') call fatal('EMPM_FULL ... Full_Diag must be ''unproj'', ' // &
                                                                              '''proj'' or ''both''')
        DoU = Fd == 'unproj' .or. Fd == 'both'
        CMx = Params%Calc%EMPM%CM_Exact
        CMv = CMx .and. Params%Calc%EMPM%CM_Variants
        CMp = Fd == 'proj' .or. Fd == 'both' .or. CMv
        ! unprojected spectrum of an Ortho basis without CM_Keep (no cm phonon & spurious states) = the SVD result ...
        LabU = trim(merge('svd   ','unproj',Ph%i_cm > 0 .and. .not. Keep))
        ! CMl ... Lawson term (CM_Beta > 0): H + beta N_cm on the complete N_cm matrix ('lawson') & with CMp also on the
            ! proj space ('projL') ...
        CMl = Params%Calc%EMPM%CM_Beta > 0.0_dp
        if ((CMx .or. CMp .or. CMl) .and. Ph%i_cm > 0 .and. .not. Keep) &
            call fatal('EMPM_FULL ... Full_Diag = ''proj''/''both'', CM_Exact & CM_Beta need complete sets (Ortho = .false., ' // &
                       'or Ortho = .true. with CM_Keep)')
        CMn = CMx .or. CMl .or. (Params%Calc%EMPM%CM_Ncm .and. (Ph%i_cm == 0 .or. Keep))
        CMa = CMx .or. CMp .or. CMl .or. CMn

        write(*,'(/,a)') 'Importing the n-phonon states (Bin/EMPM_S<n>_*.bin) ...'
        J_all = Ph%J_max
        if (Nf >= 2) then
            ! Full_N = 2 ... X of the 2-phonon states per channel only when needed (Lazy, load_x2 / free_x2), as for the 4-phonon
                ! states with Full_N = 4 (208Pb: X of all 2-phonon channels 311 GB) ...
            call EMPM_states_read(states_prefix(Params,2),2*Ph%J_max,Ph%N,Ph%N,S2,Lazy=(Nf == 2),No_Xt=.true.)
            if (S2%N == 0) call fatal('EMPM_FULL ... no 2-phonon states in Bin/EMPM_S2_*.bin, run EMPM_2PH with Save_X first')
            J_all = max(J_all, S2%J_max)
        end if
        if (Nf >= 3) then
            call EMPM_states_read(states_prefix(Params,3),3*Ph%J_max,S2%N,Ph%N,S3,No_Xt=.true.)
            if (S3%N == 0) call fatal('EMPM_FULL ... no 3-phonon states in Bin/EMPM_S3_*.bin, run EMPM_3PH first')
            J_all = max(J_all, S3%J_max)
        end if
        if (Nf >= 4) then
            ! Lazy ... X of the 4-phonon states per channel only when needed (EMPM_states_load_X): X(J,P) while J^P is built
                ! & coupled to N_cm, X(Jt,Pt) of one intermediate channel at a time in EMPM_cm_blocks ...
            call EMPM_states_read(states_prefix(Params,4),4*Ph%J_max,S3%N,Ph%N,S4,Lazy=.true.,No_Xt=.true.)
            if (S4%N == 0) call fatal('EMPM_FULL ... no 4-phonon states in Bin/EMPM_S4_*.bin, run EMPM_4PH first')
            J_all = max(J_all, S4%J_max)
        end if
        write(*,'(a)') TAB // 'States ... 1-phonon ' // istr(Ph%N) // &
                       merge(' , 2-phonon ' // istr(S2%N), repeat(' ',len(' , 2-phonon ' // istr(S2%N))), Nf >= 2) // &
                       merge(' , 3-phonon ' // istr(S3%N), repeat(' ',len(' , 3-phonon ' // istr(S3%N))), Nf >= 3) // &
                       merge(' , 4-phonon ' // istr(S4%N), repeat(' ',len(' , 4-phonon ' // istr(S4%N))), Nf >= 4)
        if (S2%Lazy) write(*,'(a)') TAB // 'Lazy ... X of the 2-phonon states loaded per channel'
        call full_mem_report('n-phonon states imported')

        write(*,'(/,a)') 'Densities, redefined interactions & the couplings to the ground state ...'
        call timer_start(t1)
        ! phonon cache C1 (densities, transposed densities & redefined interactions of all phonons, all ranks) only for
            ! Full_N >= 3 (the 3 <- 2 & 4 <- 3 couplings, the 2-phonon densities) ... with Full_N <= 2 the CM matrices are built
            ! from one phonon density at a time (EMPM_cm_dlevels_nc) ...
        if (Nf >= 3) call EMPM_cache1_make(Pairs,Ph,C1)
        call EMPM_states1_make(Ph,S1)
        call EMPM_phonon_W(Space,F_ph,Pairs,Ph,W)
        if (Nf >= 4) then
            allocate(R2(S2%N))
            !$omp parallel do schedule(dynamic)
            do a = 1, S2%N
                call EMPM_densityN(Pairs,Ph,C1,C1%R,S1,S2,a,R2(a))
            end do
            !$omp end parallel do
        end if
        ! Exact CM condition ... b_cm operator & the rank-1 densities of all 2- & 3-phonon states ...
        if (CMa) then
            write(*,'(a)') TAB // 'CM conditions ... ' // trim(merge('exact ','      ',CMx) // merge(' space','      ',CMv) // &
                           merge(' proj ','      ',CMp)) // ' (null spaces, eigenvalues < ' // &
                           cfmt_e(Params%Calc%EMPM%CM_Tol,0,1) // ')' // trim(merge(' , exact <N_cm> of the states ', ' , no <N_cm> (proj only)      ', CMn))
            call EMPM_cm_make(Params,Orb,Space,Pairs,Ph,Cm)
            call timer_start(t1)
            if (Nf <= 2) then
                ! both ranks in one pass (rank 0 for the top-level completion with CMn) ...
                call EMPM_cm_dlevels_nc(Pairs,Ph,Cm,S1,S2,max(1,merge(Nf,Nf-1,CMn)),merge(Nf,0,CMn),Dl,Dc)
                call timer_report(t1,'within-level b_cm matrices of the 1- .. ' // istr(max(1,merge(Nf,Nf-1,CMn))) // &
                                  '-phonon states' // trim(merge(' & of the commutator C (top-level completion of N_cm)', &
                                  '                                                     ', CMn)) // ', no phonon cache')
                call full_mem_report('within-level CM matrices')
            else
            call EMPM_cm_dlevels(Pairs,Ph,C1,Cm,S1,S2,S3,merge(Nf,Nf-1,CMn),Dl)
            call timer_report(t1,'within-level b_cm matrices of the 1- .. ' // istr(min(3,merge(Nf,Nf-1,CMn))) // &
                              '-phonon states')
            end if
            ! top-level completion of N_cm (level Full_N + 1 not in the basis): rank-0 within-level matrices of the commutator
                ! C = sum_mu [L_ph+, L_ph] up to level Full_N (EMPM_cm_top) ...
            if (CMn .and. Nf >= 3) then
                if (Nf > 3) then
                    ! level Full_N + 1 exists only with more than Full_N hole m-states (4He, Full_N = 4: none, <N_cm> exact)
                    if (sum(Space%Lev%j + 1, mask = Space%Lev%Occ == 1) > Nf) &
                        write(*,'(a)') TAB // 'Warning ... no top-level completion of N_cm for Full_N = 4 (<N_cm> of the 4-phonon ' // &
                                       'components a lower bound)'
                else
                    call timer_start(t1)
                    call EMPM_cm_dlevels(Pairs,Ph,C1,Cm,S1,S2,S3,Nf,Dc,0)
                    call timer_report(t1,'within-level matrices of the commutator C (top-level completion of N_cm)')
                end if
            end if
        end if
        ! h20(gamma) = <gamma|H|0> of the 0+ 2-phonon states ...
        if (Nf >= 2) then
            allocate(H20(S2%N_ch(0,1)))
        else
            allocate(H20(0))
        end if
        H20 = 0.0_dp
        if (size(H20) > 0) then
            allocate(A1(0))
            call EMPM_coupling_V(Space,F_ph,Pairs,Ph,S2%Bs(0,1),A1,1,Vk)
            call load_x2(0,1)
            H20 = matmul(transpose(S2%X(0,1)%v),Vk(:,1))
            call EMPM_states_free_X(S2,0,1)
            deallocate(A1)
        end if
        ! G1(J,P)%v(beta,lambda) = <beta|H|lambda> of the 3-phonon states beta & the phonons lambda of (J,P) ...
        if (Nf >= 3) then
            allocate(G1(0:Ph%J_max,2))
            do J = 0, Ph%J_max
                do P = 1, 2
                    n1 = Ph%N_ch(J,P)
                    n3 = 0
                    if (J <= S3%J_max) n3 = S3%N_ch(J,P)
                    allocate(G1(J,P)%v(n3,n1))
                    G1(J,P)%v = 0.0_dp
                    do i = 1, n1
                        do k = 1, S2%N_ch(0,1)
                            kk = S3%Bs(J,P)%K(S2%First(0,1)+k-1,Ph%First(J,P)+i-1)
                            if (kk > 0 .and. n3 > 0) G1(J,P)%v(:,i) = G1(J,P)%v(:,i) + H20(k)*S3%X(J,P)%v(kk,:)/sqrt(real(2*J + 1,dp))
                        end do
                    end do
                end do
            end do
        end if
        call timer_report(t1)

        ! Summary_<lab>.dat (lab = unproj or svd, proj) ... per channel the dimensions & the lowest state (the lowest physical
            ! one, <N_cm> < 0.5, when <N_cm> is known) ...
        do nn = 1, Nf
            if (DoU) call summary_open(uS(nn),LabU)
            if (CMp) call summary_open(uP(nn),'proj')
        end do
        Ch = channels(Params%Calc%EMPM%JPF,J_all)
        do k = 1, size(Ch,2)
            J = Ch(1,k)
            P = Ch(2,k)
            Tag = jp_str(J,P)
            n0 = merge(1, 0, J == 0 .and. P == 1)
            ! 1-phonon basis I1 ... the phonons of (J,P) without the spurious one ...
            if (J <= Ph%J_max) then
                I1 = pack([(Ph%First(J,P)+i-1, i = 1, Ph%N_ch(J,P))], [(Ph%First(J,P)+i-1, i = 1, Ph%N_ch(J,P))] /= Ph%i_cm &
                          .or. Keep)
            else
                allocate(I1(0))
            end if
            n1 = size(I1)
            n2 = 0
            if (Nf >= 2) then
                if (J <= ubound(S2%N_ch,1)) n2 = merge(S2%N_ch(J,P), S2%N_phys(J,P), Keep)
            end if
            n3 = 0
            if (Nf >= 3) then
                if (J <= ubound(S3%N_ch,1)) n3 = merge(S3%N_ch(J,P), S3%N_phys(J,P), Keep)
            end if
            n4 = 0
            if (Nf >= 4) then
                if (J <= ubound(S4%N_ch,1)) n4 = merge(S4%N_ch(J,P), S4%N_phys(J,P), Keep)
            end if
            o1 = n0
            o2 = o1 + n1
            o3 = o2 + n2
            o4 = o3 + n3
            N = o4 + n4
            if (N == 0) then
                deallocate(I1)
                cycle
            end if
            call timer_start(t1)
            allocate(Hf(N,N))
            Hf = 0.0_dp
            ! Diagonal ... the energies of the subspaces (0 for the HF ground state) ...
            do i = 1, n1
                Hf(o1+i,o1+i) = Ph%E(I1(i))
            end do
            do i = 1, n2
                Hf(o2+i,o2+i) = S2%E(S2%First(J,P)+i-1)
            end do
            do i = 1, n3
                Hf(o3+i,o3+i) = S3%E(S3%First(J,P)+i-1)
            end do
            do i = 1, n4
                Hf(o4+i,o4+i) = S4%E(S4%First(J,P)+i-1)
            end do
            ! CM_Keep ... within-level couplings: <nu|A|cm> of the 1- phonons, <physical|H|spurious> of the n-phonon sets
                ! (lower triangle, the spurious states after the physical ones) ...
            if (Keep) then
                do i = 1, n1
                    if (I1(i) /= Ph%i_cm) cycle
                    do ii = 1, n1
                        if (ii /= i) Hf(o1+max(i,ii),o1+min(i,ii)) = Ph%R(I1(ii))
                    end do
                end do
                call keep_rc(S2,n2,o2,Nf >= 2)
                call keep_rc(S3,n3,o3,Nf >= 3)
                call keep_rc(S4,n4,o4,Nf >= 4)
            end if
            ! 2 <- 0 ...
            if (n0 == 1 .and. n2 > 0) Hf(o2+1:o2+n2,1) = H20(1:n2)
            ! X of the 2-phonon states of J^P (Lazy S2, Full_N = 2) ...
            if (n2 > 0) call load_x2(J,P)
            ! 2 <- 1 ...
            if (n1 > 0 .and. n2 > 0) then
                A1 = I1
                call EMPM_coupling_V(Space,F_ph,Pairs,Ph,S2%Bs(J,P),A1,0,Vk)
                Hf(o2+1:o2+n2,o1+1:o1+n1) = mtm(S2%X(J,P)%v(:,1:n2),Vk)
            end if
            ! 3 <- 1 ...
            if (n1 > 0 .and. n3 > 0) Hf(o3+1:o3+n3,o1+1:o1+n1) = G1(J,P)%v(1:n3,I1-Ph%First(J,P)+1)
            ! 3 <- 2 ...
            if (n2 > 0 .and. n3 > 0) then
                allocate(Vc(S3%Bs(J,P)%N_full,n2))
                !$omp parallel do schedule(dynamic) private(Rho)
                do ia = 1, n2
                    call EMPM_densityN(Pairs,Ph,C1,C1%R,S1,S2,S2%First(J,P)+ia-1,Rho)
                    call EMPM_coupling_Vn(Pairs,Ph,W,Rho,J,S2%J_max,S2%N_ch,S2%First,S3%Bs(J,P),Vc(:,ia))
                end do
                !$omp end parallel do
                Hf(o3+1:o3+n3,o2+1:o2+n2) = mtm(S3%X(J,P)%v(:,1:n3),Vc)
                deallocate(Vc)
            end if
            ! X of the 4-phonon states of J^P (Lazy) ...
            if (Nf >= 4 .and. n4 > 0) call load_x4(J,P)
            ! 4 <- 3 ...
            if (n3 > 0 .and. n4 > 0) then
                allocate(Vc(S4%Bs(J,P)%N_full,n3))
                !$omp parallel do schedule(dynamic) private(Rho)
                do ia = 1, n3
                    call EMPM_densityN(Pairs,Ph,C1,R2,S2,S3,S3%First(J,P)+ia-1,Rho)
                    call EMPM_coupling_Vn(Pairs,Ph,W,Rho,J,S3%J_max,S3%N_ch,S3%First,S4%Bs(J,P),Vc(:,ia))
                end do
                !$omp end parallel do
                Hf(o4+1:o4+n4,o3+1:o3+n3) = mtm(S4%X(J,P)%v(:,1:n4),Vc)
                deallocate(Vc)
            end if
            ! 4 <- 2 ... G2(beta,alpha) = 1/sqrt(2J+1) sum_(k = mu nu) C^alpha_k sum_d G1(d,mu) X^beta(d nu) ...
            if (n2 > 0 .and. n4 > 0) then
                allocate(Blk(n4,n2), Uv(n4))
                Blk = 0.0_dp
                associate(B2 => S2%Bs(J,P), B4 => S4%Bs(J,P))
                    do ic = 1, size(S2%Kc(J,P)%v)
                        kk = S2%Kc(J,P)%v(ic)
                        mu = B2%Al(kk)
                        nu = B2%La(kk)
                        jm = Ph%J(mu)
                        pm = Ph%P(mu)
                        if (jm > S3%J_max) cycle
                        nd = S3%N_ch(jm,pm)
                        Uv = 0.0_dp
                        do d = 1, nd
                            i = B4%K(S3%First(jm,pm)+d-1,nu)
                            if (i == 0) cycle
                            Uv = Uv + G1(jm,pm)%v(d,mu-Ph%First(jm,pm)+1) * S4%X(J,P)%v(i,1:n4)
                        end do
                        do ia = 1, n2
                            Blk(:,ia) = Blk(:,ia) + S2%Cc(J,P)%v(ic,ia) * Uv
                        end do
                    end do
                end associate
                Hf(o4+1:o4+n4,o2+1:o2+n2) = Blk / sqrt(real(2*J + 1,dp))
                deallocate(Blk, Uv)
            end if
            ! Symmetric completion, then the sequential diagonalizations with up to nn = 1 .. Nf phonons ...
            do i = 1, N
                Hf(i,i+1:N) = Hf(i+1:N,i)
            end do
            call timer_report(t1,'full Hamiltonian J^P = ' // Tag // ' built (dimension ' // istr(N) // ' = ' // istr(n0) // &
                              ' + ' // istr(n1) // ' + ' // istr(n2) // ' + ' // istr(n3) // ' + ' // istr(n4) // ')')
            ! completeness of the intermediate channels: N_cm (& proj) of J^P need the states of J' = J-1 .. J+1 (opposite
                ! parity) on every level q = 2 .. Full_N; a channel of level q whose file is missing although it can be built
                ! from the loaded level q-1 & the phonons was not computed (JP<q>) -> warning ...
            if (Nf >= 4 .and. .not. CMa) call EMPM_states_free_X(S4,J,P)
            if (CMa) call check_channels()
            ! N_cm on the channel basis, Mcm = G2/(A (2J+1)) sum_J' B^T B with B = <J' states||L||J states> ...
            if (CMa) then
                call timer_start(t1)
                Nbc = EMPM_cm_dims(Ph,S2,S3,S4,Nf,J,P)
                if (sum(Nbc) /= N) call fatal('EMPM_FULL ... CM basis mismatch in J^P = ' // Tag)
                if (CMn) then
                    allocate(Mcm(N,N))
                    Mcm = 0.0_dp
                end if
                ! space & proj ... Ml(:,:,q) = the part of Mcm from the intermediate states of level q ...
                if (CMp) then
                    allocate(Ml(N,N,0:Nf))
                    Ml = 0.0_dp
                end if
                Pt = 3 - P
                do Jt = max(J - 1, 0), J + 1
                    if (J == 0 .and. Jt == 0) cycle
                    if (Nf >= 4) call load_x4(Jt,Pt)
                    if (Nf == 2) call load_x2(Jt,Pt)
                    call EMPM_cm_blocks(Ph,Cm,S1,S2,S3,S4,Dl,Nf,J,P,Jt,Pt,Bph,Bhp,Bd,CMn)
                    if (Nf >= 4) call EMPM_states_free_X(S4,Jt,Pt)
                    if (Nf == 2) call EMPM_states_free_X(S2,Jt,Pt)
                    if (size(Bph,1) == 0) cycle
                    Bt = Sg_ph(Jt-J)*Bph + Sg_hp(Jt-J)*Bhp + Sg_d(Jt-J)*Bd
                    if (CMp) then
                        Mbc = EMPM_cm_dims(Ph,S2,S3,S4,Nf,Jt,Pt)
                        kr = 0
                        do q = 0, Nf
                            if (Mbc(q) > 0) Ml(:,:,q) = Ml(:,:,q) + mtm(Bt(kr+1:kr+Mbc(q),:),Bt(kr+1:kr+Mbc(q),:))
                            kr = kr + Mbc(q)
                        end do
                    end if
                    if (CMn) Mcm = Mcm + mtm(Bt,Bt)
                    if (Params%Calc%EMPM%CM_Dump) then
                        open(newunit=uC, file=io_path(Params,'/Bin/EMPM_CMB_' // Tag // '_' // jp_str(Jt,Pt) // '.bin'), &
                             access='stream', form='unformatted', status='replace', action='write')
                        write(uC) size(Bph,1), size(Bph,2)
                        write(uC) Bph, Bhp, Bd
                        close(uC)
                    end if
                end do
                if (Nf >= 4) call EMPM_states_free_X(S4,J,P)
                if (CMn) Mcm = Mcm * Cm%G2 / real(Cm%A*(2*J + 1),dp)
                if (CMp) Ml = Ml * Cm%G2 / real(Cm%A*(2*J + 1),dp)
                ! top-level completion: the part of N_cm from the intermediate states of level Full_N + 1 ...
                if (CMn .and. Nf <= 3) then
                    if (Nbc(Nf) > 0) then
                        Mt = cm_top(Nf)
                        q = sum(Nbc(0:Nf-1))
                        Mcm(q+1:q+Nbc(Nf),q+1:q+Nbc(Nf)) = Mcm(q+1:q+Nbc(Nf),q+1:q+Nbc(Nf)) + Mt
                    end if
                    ! check (NJ_CMTOP_CHECK set, with CM_Proj): the completion of level Full_N - 1 against its exact part
                        ! Ml(level Full_N) from the rows in the basis ...
                    call get_environment_variable('NJ_CMTOP_CHECK',Chk)
                    if (len_trim(Chk) > 0 .and. CMp .and. Nf >= 2) then
                        if (Nbc(Nf-1) > 0) then
                            Mt = cm_top(Nf-1)
                            q = sum(Nbc(0:Nf-2))
                            write(*,'(a)') TAB // 'CMTOP check J^P = ' // Tag // ' level ' // istr(Nf-1) // ' ... max |completion - exact| = ' // &
                                cfmt_e(maxval(abs(Mt - Ml(q+1:q+Nbc(Nf-1),q+1:q+Nbc(Nf-1),Nf))),0,3) // ' , max |exact| = ' // &
                                cfmt_e(maxval(abs(Ml(q+1:q+Nbc(Nf-1),q+1:q+Nbc(Nf-1),Nf))),0,3) // ' , trace completion / exact = ' // &
                                cfmt_e(sum([(Mt(i,i), i = 1, Nbc(Nf-1))]),0,6) // ' / ' // &
                                cfmt_e(sum([(Ml(q+i,q+i,Nf), i = 1, Nbc(Nf-1))]),0,6)
                        end if
                    end if
                end if
                if (Params%Calc%EMPM%CM_Dump) then
                    open(newunit=uC, file=io_path(Params,'/Bin/EMPM_CMH_' // Tag // '.bin'), access='stream', &
                         form='unformatted', status='replace', action='write')
                    write(uC) N, [n0, n1, n2, n3, n4]
                    write(uC) Hf
                    close(uC)
                end if
                call timer_report(t1,trim(merge('N_cm & b_cm parts of J^P =  ','b_cm parts (proj) of J^P =  ',CMn)) // ' ' // Tag // ' built')
            end if
            do nn = 1, Nf
                Nb = [n0, n1, n2, n3, n4]
                Nb(nn+1:4) = 0
                Ns = sum(Nb)
                if (Ns == 0) cycle
                call timer_start(t1)
                ! H without CM projection (Full_Diag 'unproj' / 'both') ...
                if (DoU) then
                    Hs = Hf(1:Ns,1:Ns)
                    call EMPM_syevd(Hs,Ef)
                    call state_ncm(Hs,Ncm0)
                    call write_full(LabU,Ef,Hs,Ncm0,.true.)
                    deallocate(Hs, Ef, Ncm0)
                end if
                ! Full_Svd (with CM_Keep) ... the SVD (Ortho) result: H on the physical states only, with its exact <N_cm> ...
                if (Keep .and. Params%Calc%EMPM%Full_Svd) call svd_ncm(nn)
                ! H on null spaces: 'space' ... intermediate levels <= nn, 'proj' ... <= nn-1, 'exact' ... complete N_cm ...
                if (CMv) then
                    Ucm = sum(Ml(1:Ns,1:Ns,0:nn),dim=3)
                    call cm_null('space')
                end if
                if (CMp) then
                    Ucm = sum(Ml(1:Ns,1:Ns,0:nn-1),dim=3)
                    call cm_null('proj')
                end if
                ! Lawson term H + beta N_cm: on the <= nn phonon space & on the proj space ...
                if (CMl) then
                    call lawson('lawson')
                    if (CMp) call lawson('projL',Ur)
                end if
                if (allocated(Ur)) deallocate(Ur)
                if (CMx) then
                    Ucm = Mcm(1:Ns,1:Ns)
                    call cm_null('exact')
                end if
                call timer_report(t1,'J^P = ' // Tag // ' with up to ' // istr(nn) // ' phonons (dimension ' // istr(Ns) // ')')
            end do
            if (allocated(Mcm)) deallocate(Mcm)
            if (allocated(Ml)) deallocate(Ml)
            deallocate(Hf, I1)
            if (Nf >= 4) call EMPM_states_free_X(S4,J,P)
            if (Nf == 2) then
                call full_mem_report('channel ' // Tag)
                call EMPM_states_free_X(S2,J,P)
            end if
        end do
        do nn = 1, Nf
            if (DoU) close(uS(nn))
            if (CMp) close(uP(nn))
        end do
        call timer_report(t0)
        write(*,'(/,a,/)') 'Full EMPM Hamiltonian has finished ...'

    contains

        ! X of the 2-phonon states of channel (Jx,Px) from the state file (Lazy S2), with time & size in the log ...
        subroutine load_x2(Jx,Px)
            integer, intent(in) :: Jx, Px
            integer(i8) :: tx
            if (.not. S2%Lazy) return
            if (Jx > ubound(S2%N_ch,1)) return
            if (S2%N_ch(Jx,Px) == 0 .or. size(S2%X(Jx,Px)%v) > 0) return
            call timer_start(tx)
            call EMPM_states_load_X(S2,Jx,Px)
            call timer_report(tx,'X of the 2-phonon states J^P = ' // jp_str(Jx,Px) // ' loaded (' // &
                              cfmt_f(8.0_dp*real(size(S2%X(Jx,Px)%v,kind=i8),dp)/1.0e9_dp,0,3) // ' GB)')
        end subroutine load_x2

        ! Current & peak resident memory (VmRSS, VmHWM of /proc/self/status) in the log, flushed ...
        subroutine full_mem_report(Label)
            character(len=*), intent(in) :: Label
            integer :: um, Ios
            character(len=256) :: Line
            integer(i8) :: kb, Rss, Hwm
            Rss = -1
            Hwm = -1
            open(newunit=um, file='/proc/self/status', status='old', action='read', iostat=Ios)
            if (Ios /= 0) return
            do
                read(um,'(a)',iostat=Ios) Line
                if (Ios /= 0) exit
                if (index(Line,'VmRSS:') == 1) then
                    read(Line(7:),*,iostat=Ios) kb
                    if (Ios == 0) Rss = kb
                else if (index(Line,'VmHWM:') == 1) then
                    read(Line(7:),*,iostat=Ios) kb
                    if (Ios == 0) Hwm = kb
                end if
            end do
            close(um)
            if (Rss < 0) return
            write(*,'(a)') TAB // 'memory ... RSS ' // cfmt_f(real(Rss,dp)/976562.5_dp,0,1) // ' GB, peak ' // &
                           cfmt_f(real(Hwm,dp)/976562.5_dp,0,1) // ' GB (' // Label // ')'
            flush(6)
        end subroutine full_mem_report

        ! X of the 4-phonon states of channel (Jx,Px) from the state file (Lazy S4), with time & size in the log ...
        subroutine load_x4(Jx,Px)
            integer, intent(in) :: Jx, Px
            integer(i8) :: tx
            if (Jx > ubound(S4%N_ch,1)) return
            if (S4%N_ch(Jx,Px) == 0 .or. size(S4%X(Jx,Px)%v) > 0) return
            call timer_start(tx)
            call EMPM_states_load_X(S4,Jx,Px)
            call timer_report(tx,'X of the 4-phonon states J^P = ' // jp_str(Jx,Px) // ' loaded (' // &
                              cfmt_f(8.0_dp*real(size(S4%X(Jx,Px)%v,kind=i8),dp)/1.0e9_dp,0,3) // ' GB)')
        end subroutine load_x4

        ! warning for intermediate channels (J' = J-1 .. J+1, opposite parity) of the levels 2 .. Full_N that were not computed ...
        subroutine check_channels()
            integer :: q, Jx, Px
            logical :: Ex
            Px = 3 - P
            do q = 2, Nf
                do Jx = max(J - 1, 0), J + 1
                    if (J == 0 .and. Jx == 0) cycle
                    inquire(file=states_file(Params,q,jp_str(Jx,Px)), exist=Ex)
                    if (Ex) cycle
                    if (.not. can_build(q,Jx,Px)) cycle
                    write(*,'(a)') TAB // 'Warning ... J^P = ' // Tag // ' : the ' // istr(q) // '-phonon channel ' // &
                                   jp_str(Jx,Px) // ' was not computed (add it to JP' // istr(q) // '), <N_cm> of ' // Tag // &
                                   ' & the CM projection are incomplete'
                end do
            end do
        end subroutine check_channels

        ! whether a channel (Jx,Px) of level q can be built from the loaded states of level q-1 & the phonons ...
        logical function can_build(q,Jx,Px)
            integer, intent(in) :: q, Jx, Px
            integer :: Ja, Pa, Jl, Pl, Jm
            can_build = .false.
            Jm = Ph%J_max
            if (q == 3) Jm = S2%J_max
            if (q == 4) Jm = S3%J_max
            do Ja = 0, Jm
                do Pa = 1, 2
                    if (q == 2) then
                        if (Ph%N_ch(Ja,Pa) == 0) cycle
                    else if (q == 3) then
                        if (Ja > ubound(S2%N_ch,1)) cycle
                        if (S2%N_ch(Ja,Pa) == 0) cycle
                    else
                        if (Ja > ubound(S3%N_ch,1)) cycle
                        if (S3%N_ch(Ja,Pa) == 0) cycle
                    end if
                    do Jl = 0, Ph%J_max
                        do Pl = 1, 2
                            if (Ph%N_ch(Jl,Pl) == 0) cycle
                            if (merge(1,2,Pa == Pl) /= Px) cycle
                            if (Jx < abs(Ja - Jl) .or. Jx > Ja + Jl) cycle
                            can_build = .true.
                            return
                        end do
                    end do
                end do
            end do
        end function can_build

        ! top-level completion of N_cm on level n of the current channel ...
        function cm_top(n) result(Mo)
            integer, intent(in) :: n
            real(dp), allocatable :: Mo(:,:)
            select case (n)
            case (1)
                Mo = EMPM_cm_top(Ph,Cm,S1,S1,Dc(1)%v,1,J,P)
            case (2)
                Mo = EMPM_cm_top(Ph,Cm,S1,S2,Dc(2)%v,2,J,P)
            case (3)
                Mo = EMPM_cm_top(Ph,Cm,S2,S3,Dc(3)%v,3,J,P)
            end select
        end function cm_top

        ! <physical|H|spurious> block Rc of the n-phonon set S in the basis (offset o, n states of channel J^P) ...
        subroutine keep_rc(S,n,o,Use)
            type(empm_states2_t), intent(in) :: S
            integer, intent(in) :: n, o
            logical, intent(in) :: Use
            integer :: i2, s2
            if (.not. Use .or. n == 0) return
            nph = S%N_phys(J,P)
            do s2 = 1, n - nph
                do i2 = 1, nph
                    Hf(o+nph+s2,o+i2) = S%Rc(J,P)%v(i2,s2)
                end do
            end do
        end subroutine keep_rc

        ! CM condition on the <= nn phonon space (dimension Ns): Ucm = matrix of the CM condition (in: Mcm = exact N_cm, or
            ! its part from the intermediate levels <= nn / <= nn-1), H diagonalized on its null space (eigenvalues < CM_Tol)
            ! ... the states (eigenvalues, eigenvectors on the <= nn basis) & their exact <N_cm> to Full_<Lab>_<JP>.dat
            ! (write_full), the dimensions & <N_cm> range in the log ...
        subroutine cm_null(Lab)
            character(len=*), intent(in) :: Lab
            real(dp), allocatable :: Nc(:), Tq(:,:)
            character(len=:), allocatable :: Ln
            integer :: i2, uq
            call EMPM_syevd(Ucm,Ncm)
            Np = count(Ncm < Params%Calc%EMPM%CM_Tol)
            if (Np > 0) then
                Hp = mtm(Ucm(:,1:Np),mm(Hf(1:Ns,1:Ns),Ucm(:,1:Np)))
                call EMPM_syevd(Hp,Ef)
                Hs = mm(Ucm(:,1:Np),Hp)
                deallocate(Hp)
                call state_ncm(Hs,Nc)
            else
                allocate(Hs(Ns,0), Ef(0), Nc(0))
            end if
            Ln = ''
            if (size(Nc) > 0) Ln = ' , exact <N_cm> of the lowest / all states ' // cfmt_e(Nc(1),0,2) // ' / ' // &
                             cfmt_e(minval(Nc),0,2) // ' .. ' // cfmt_e(maxval(Nc),0,2)
            write(*,'(a)') TAB // 'J^P = ' // Tag // ' , up to ' // istr(nn) // ' phonons ... CM ' // Lab // ' ... N_cm = 0 ' // &
                           'states ' // istr(Np) // ' of ' // istr(Ns) // ' , eigenvalues ' // cfmt_e(Ncm(1),0,2) // ' .. ' // &
                           cfmt_e(Ncm(Ns),0,2) // ' , gap ' // cfmt_e(Ncm(max(Np,1)),0,2) // ' | ' // &
                           cfmt_e(Ncm(min(Np+1,Ns)),0,2) // Ln
            if (Lab == 'proj' .and. CMl) Ur = Ucm(:,1:Np)
            call write_full(Lab,Ef,Hs,Nc,Lab == 'proj')
            deallocate(Hs, Ef)
            deallocate(Ucm, Ncm, Nc)
        end subroutine cm_null

        ! exact <N_cm> of the states (columns of Z on the <= nn basis), empty without CMn ...
        subroutine state_ncm(Z,Nc)
            real(dp), intent(in) :: Z(:,:)
            real(dp), allocatable, intent(out) :: Nc(:)
            real(dp), allocatable :: Tq(:,:)
            integer :: i2
            if (.not. CMn .or. size(Z,2) == 0) then
                allocate(Nc(0))
                return
            end if
            Tq = mm(Mcm(1:Ns,1:Ns),Z)
            Nc = [(dot_product(Z(:,i2),Tq(:,i2)), i2 = 1, size(Z,2))]
        end subroutine state_ncm

        ! Output of one diagonalization of the full H with up to nn phonons, the same for all variants (Lab = unproj, proj,
            ! svd, lawson, projL, exact, space): EMPM/Full/<nn>Phonon/Full_<Lab>_<JP>.dat ... i, E (MeV, relative to E_HF),
            ! exact <N_cm> (-1 without CMn), phys = <N_cm> < 0.5 (1 for all states without <N_cm>), 0- .. 4-phonon weights;
            ! Bin/EMPM_FULL<nn>_<JP>_E_<Lab>.bin & _Ncm_<Lab>.bin, with Zb also _Z_<Lab>.bin (Nb, then the eigenvectors on
            ! the <= nn basis); the lowest states in the log & Summary_<Lab>.dat (unproj / svd, proj) ...
        subroutine write_full(Lab,Ev,Z,Nc,Zb)
            character(len=*), intent(in) :: Lab
            real(dp), intent(in) :: Ev(:), Z(:,:), Nc(:)
            logical, intent(in) :: Zb
            integer :: i2, q, uq, Off(0:5), Nph, i1, us1
            real(dp) :: w(0:4), Nv
            logical :: Sm
            character(len=:), allocatable :: Ds, Dir, Bf
            select case (Lab)
            case ('unproj')
                Ds = 'no CM projection'
            case ('proj')
                Ds = 'CM projection P_(n-1) b_cm psi = 0'
            case ('svd')
                Ds = 'SVD (Ortho) physical states'
            case ('exact')
                Ds = 'exact CM condition b_cm psi = 0'
            case ('space')
                Ds = 'CM condition P_n b_cm psi = 0'
            case ('lawson')
                Ds = 'Lawson H + beta N_cm (beta ' // cfmt_f(Params%Calc%EMPM%CM_Beta,0,1) // ' MeV), E = <H>'
            case ('projL')
                Ds = 'Lawson H + beta N_cm on the proj space (beta ' // cfmt_f(Params%Calc%EMPM%CM_Beta,0,1) // ' MeV), E = <H>'
            case default
                Ds = Lab
            end select
            Off(0) = 0
            do q = 0, 4
                Off(q+1) = Off(q) + Nb(q)
            end do
            Nph = size(Ev)
            if (size(Nc) > 0) Nph = count(Nc < 0.5_dp)
            Dir = '/EMPM/Full/' // istr(nn) // 'Phonon'
            open(newunit=uq, file=io_path(Params,Dir // '/Full_' // Lab // '_' // Tag // '.dat'), status='replace', action='write')
            write(uq,'(a)') 'EMPM_FULL J^P = ' // Tag // ' with up to ' // istr(nn) // ' phonons, ' // Ds // ' ... basis ' // &
                            istr(Nb(0)) // ' + ' // istr(Nb(1)) // ' + ' // istr(Nb(2)) // ' + ' // istr(Nb(3)) // ' + ' // &
                            istr(Nb(4)) // ' (0- .. 4-phonon), ' // istr(size(Ev)) // ' states, ' // istr(Nph) // ' physical' // &
                            trim(merge(' (<N_cm> < 0.5)       ',' (<N_cm> not computed)',size(Nc) > 0))
            write(uq,'(/,a)') cfmt_s('i',8) // cfmt_s('E',16) // cfmt_s('N_cm',14) // cfmt_s('phys',6) // cfmt_s('w0',12) // &
                              cfmt_s('w1',12) // cfmt_s('w2',12) // cfmt_s('w3',12) // cfmt_s('w4',12)
            ! the main spectra (unproj or svd without CM_Keep, proj) go to the summaries & the log ...
            Sm = .false.
            us1 = 0
            if (Lab == 'proj' .and. CMp) then
                Sm = .true.
                us1 = uP(nn)
            else if (Lab == LabU .and. DoU) then
                Sm = .true.
                us1 = uS(nn)
            end if
            if (Sm) write(*,'(a)') TAB // 'J^P = ' // Tag // ' , up to ' // istr(nn) // ' phonons, ' // Lab // &
                                        ' ... lowest states:      E        N_cm        w0        w1        w2        w3        w4'
            i1 = 0
            do i2 = 1, size(Ev)
                do q = 0, 4
                    w(q) = sum(Z(Off(q)+1:Off(q+1),i2)**2)
                end do
                Nv = -1.0_dp
                if (size(Nc) > 0) Nv = Nc(i2)
                write(uq,'(a)') cfmt_d(i2,8) // cfmt_f(Ev(i2),16,6) // cfmt_f(Nv,14,8) // cfmt_d(merge(1,0,Nv < 0.5_dp),6) // &
                                cfmt_f(w(0),12,6) // cfmt_f(w(1),12,6) // cfmt_f(w(2),12,6) // cfmt_f(w(3),12,6) // cfmt_f(w(4),12,6)
                if (Sm .and. i2 <= 5) write(*,'(a)') TAB // TAB // cfmt_d(i2,4) // cfmt_f(Ev(i2),14,5) // cfmt_f(Nv,12,5) // &
                    cfmt_f(w(0),10,5) // cfmt_f(w(1),10,5) // cfmt_f(w(2),10,5) // cfmt_f(w(3),10,5) // cfmt_f(w(4),10,5)
                if (Sm .and. i1 == 0 .and. Nv < 0.5_dp) then
                    i1 = i2
                    write(us1,'(a)') cfmt_s(Tag,6) // cfmt_d(Nb(0),4) // cfmt_d(Nb(1),6) // cfmt_d(Nb(2),8) // cfmt_d(Nb(3),8) // &
                                     cfmt_d(Nb(4),8) // cfmt_d(size(Ev),8) // cfmt_d(Nph,8) // cfmt_d(i2,6) // cfmt_f(Ev(i2),14,5) // &
                                     cfmt_f(Nv,12,6) // cfmt_f(w(0),10,5) // cfmt_f(w(1),10,5) // cfmt_f(w(2),10,5) // &
                                     cfmt_f(w(3),10,5) // cfmt_f(w(4),10,5)
                    flush(us1)
                end if
            end do
            close(uq)
            Bf = '/Bin/EMPM_FULL' // istr(nn) // '_' // Tag // '_'
            open(newunit=uq, file=io_path(Params,Bf // 'E_' // Lab // '.bin'), access='stream', form='unformatted', &
                 status='replace', action='write')
            write(uq) Ev
            close(uq)
            if (size(Nc) > 0) then
                open(newunit=uq, file=io_path(Params,Bf // 'Ncm_' // Lab // '.bin'), access='stream', form='unformatted', &
                     status='replace', action='write')
                write(uq) Nc
                close(uq)
            end if
            if (Zb) then
                open(newunit=uq, file=io_path(Params,Bf // 'Z_' // Lab // '.bin'), access='stream', form='unformatted', &
                     status='replace', action='write')
                write(uq) Nb
                write(uq) Z
                close(uq)
            end if
        end subroutine write_full

        ! Summary_<Lab>.dat of truncation nn ...
        subroutine summary_open(us1,Lab)
            integer, intent(out) :: us1
            character(len=*), intent(in) :: Lab
            open(newunit=us1, file=io_path(Params,'/EMPM/Full/' // istr(nn) // 'Phonon/Summary_' // Lab // '.dat'), &
                 status='replace', action='write')
            write(us1,'(a)') 'EMPM_FULL with up to ' // istr(nn) // ' phonons, ' // Lab // ' ... per J^P the dimensions n0 .. n4, ' // &
                             'the number of states & of physical states (<N_cm> < 0.5), the lowest physical state: index, E ' // &
                             '(relative to E_HF), <N_cm> (-1 not computed) & its 0- .. 4-phonon weights'
            write(us1,'(/,a)') cfmt_s('J^P',6) // cfmt_s('n0',4) // cfmt_s('n1',6) // cfmt_s('n2',8) // cfmt_s('n3',8) // &
                               cfmt_s('n4',8) // cfmt_s('N',8) // cfmt_s('N_phys',8) // cfmt_s('i',6) // cfmt_s('E',14) // &
                               cfmt_s('N_cm',12) // cfmt_s('w0',10) // cfmt_s('w1',10) // cfmt_s('w2',10) // cfmt_s('w3',10) // &
                               cfmt_s('w4',10)
        end subroutine summary_open

        ! CM_Keep (with <N_cm>) ... the SVD (Ortho = .true.) states: H diagonalized on the physical states of each level (the
            ! ground state, the phonons without the cm one, the first N_phys states of the n-phonon channels), their exact
            ! <N_cm> (Mcm of the complete basis, intermediate levels <= Full_N) ... Full_svd_<JP>.dat (write_full), the range
            ! of <N_cm> in the log ...
        subroutine svd_ncm(nn)
            integer, intent(in) :: nn
            integer, allocatable :: Ix(:)
            real(dp), allocatable :: Hx(:,:), Ex(:), Zx(:,:), Tq(:,:), Nc(:)
            integer :: i2, uq
            allocate(Ix(0))
            if (n0 == 1) Ix = [Ix, 1]
            if (nn >= 1) Ix = [Ix, pack([(o1+i2, i2 = 1, n1)], I1 /= Ph%i_cm)]
            if (nn >= 2 .and. n2 > 0) Ix = [Ix, [(o2+i2, i2 = 1, S2%N_phys(J,P))]]
            if (nn >= 3 .and. n3 > 0) Ix = [Ix, [(o3+i2, i2 = 1, S3%N_phys(J,P))]]
            if (nn >= 4 .and. n4 > 0) Ix = [Ix, [(o4+i2, i2 = 1, S4%N_phys(J,P))]]
            if (size(Ix) == 0) return
            Hx = Hf(Ix,Ix)
            call EMPM_syevd(Hx,Ex)
            allocate(Zx(Ns,size(Ix)))
            Zx = 0.0_dp
            Zx(Ix,:) = Hx
            call state_ncm(Zx,Nc)
            write(*,'(a)') TAB // 'J^P = ' // Tag // ' , up to ' // istr(nn) // ' phonons ... SVD (physical) states ' // &
                           istr(size(Ix)) // ' of ' // istr(Ns) // ' , lowest E ' // cfmt_f(Ex(1),0,5)
            if (size(Nc) > 0) write(*,'(a)') TAB // TAB // 'exact <N_cm> of the lowest / all SVD states ' // cfmt_e(Nc(1),0,2) // &
                                            ' / ' // cfmt_e(minval(Nc),0,2) // ' .. ' // cfmt_e(maxval(Nc),0,2)
            call write_full('svd',Ex,Zx,Nc,.false.)
        end subroutine svd_ncm

        ! Lawson term (CM_Beta) on the <= nn phonon space (or on the subspace spanned by the columns of Ub, e.g. the proj
            ! space): H + beta N_cm diagonalized, E = <H> = eigenvalue - beta <N_cm> & <N_cm> in the order of the eigenvalues
            ! (as the NCSM) ... Full_<Lab>_<JP>.dat (write_full); N_cm exact for nn < Full_N and, with the
            ! top-level completion (EMPM_cm_top), also at nn = Full_N <= 3 (lower bound for Full_N = 4, exact if no level 5 exists, e.g. 4He) ...
        subroutine lawson(Lab,Ub)
            character(len=*), intent(in) :: Lab
            real(dp), intent(in), optional :: Ub(:,:)
            real(dp), allocatable :: Hl(:,:), El(:), Zl(:,:), Tq(:,:), Nc(:)
            real(dp) :: Bt0
            integer :: i2, uq, nl
            Bt0 = Params%Calc%EMPM%CM_Beta
            Hl = Hf(1:Ns,1:Ns) + Bt0*Mcm(1:Ns,1:Ns)
            if (present(Ub)) then
                if (size(Ub,2) == 0) return
                Hl = mtm(Ub,mm(Hl,Ub))
            end if
            call EMPM_syevd(Hl,El)
            if (present(Ub)) then
                Zl = mm(Ub,Hl)
            else
                call move_alloc(Hl,Zl)
            end if
            nl = size(El)
            Tq = mm(Mcm(1:Ns,1:Ns),Zl)
            Nc = [(dot_product(Zl(:,i2),Tq(:,i2)), i2 = 1, nl)]
            El = El - Bt0*Nc
            write(*,'(a)') TAB // 'J^P = ' // Tag // ' , up to ' // istr(nn) // ' phonons ... CM ' // Lab // ' (beta ' // &
                           cfmt_f(Bt0,0,1) // ' MeV) ... ' // istr(nl) // ' states, lowest E = <H> ' // cfmt_f(El(1),0,5) // &
                           ' , <N_cm> ' // cfmt_e(Nc(1),0,2)
            call write_full(Lab,El,Zl,Nc,.false.)
        end subroutine lawson


    end subroutine EMPM_FULL

    ! B(E1; 0+_gs -> 1-_k) of the full EMPM states (after EMPM_FULL with Full_N <= 2, JPF with 0+ & 1-, CM_Keep or no CM
        ! phonon) ... reduced matrix elements <1- basis||Q_t||0+ basis> of the E1 operator of protons & neutrons (t = p, n)
        ! between the EMPM bases [0 ; 1-phonon ; 2-phonon] by the machinery of the CM operator (EMPM_cm_blocks: ph part
        ! 1/sqrt3 sum_lambda bL X, hp part, pp/hh part by the density recursion; Q Hermitian: bR = bL, EMPM_op_make), the
        ! within-level matrices only for the block 1- x 0+ (EMPM_cm_dlevels_nc with Sel) ... M_t(k) = Z_1-(:,k)^T B_t
        ! Z_0+(:,gs) with the eigenvectors of Bin/EMPM_FULL<n>_<JP>_Z_<lab>.bin (lab = unproj or svd, proj), n = 1 .. Full_N,
        ! gs = the lowest 0+ state with <N_cm> < 0.5 (the lowest one without <N_cm>) ... B_ph = M_p^2, B_is = (M_p + M_n)^2/4,
        ! B_iv = (e_p M_p - e_n M_n)^2 (recoil charges N/A, Z/A), as EMPM_TDA_strengths & the NCSM transitions ...
        ! EMPM/Transitions/Full<n>_<lab>_E1.dat (k, E, Ex = E - E_gs, N_cm, phys, B_ph, B_is, B_iv) ...
    subroutine EMPM_TRANS(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(empm_phonons_t) :: Ph
        type(empm_pairs_t) :: Pairs
        type(empm_states2_t) :: S1, S2, S3, S4
        type(empm_cm_t) :: Op(2)
        type(rmat_t) :: Dl(3), Dc(3), Bt(2)
        real(dp), allocatable :: Bph(:,:), Bhp(:,:), Bd(:,:), Z0(:,:), Z1(:,:), E0(:), E1(:), N0(:), N1(:), v(:), Mt(:,:)
        integer :: Nf, nn, it, Nb0(0:4), Nb1(0:4), Mb(0:4), Mb1(0:4), r1, c0, g, k, u, il, n0s, n1s, ig
        real(dp), allocatable :: Psi(:)
        real(dp) :: Eg, Ng
        character(len=:), allocatable :: Gdesc, Fsuf
        integer(i8) :: t0, t1
        logical :: Keep, Ex0, Ex1
        real(dp) :: e_p, e_n, B(3), S(3), Sp(3), M1
        character(len=8) :: Labs(2)
        character(len=:), allocatable :: Lab

        call wigner_init_float(75,'Jmax',9)
        !$ call omp_set_max_active_levels(1)
        Nf = Params%Calc%EMPM%Full_N
        if (Nf < 1 .or. Nf > 2) call fatal('EMPM_TRANS ... Full_N must be 1 or 2')
        write(*,'(a)') 'Starting B(E1; 0+ -> 1-) of the full EMPM states with up to 1 .. ' // istr(Nf) // ' phonons ...'
        write(*,'(a)') trim(Params%Calc%Path)
        call timer_start(t0)
        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call EMPM_space_make(Params,Orb,Space)
        call Fmat_import(Space,io_path(Params,'/Bin/EMPM_F.bin'),F_ph)
        call EMPM_phonons_import(Params,Space,F_ph,Ph)
        Keep = Params%Calc%EMPM%CM_Keep .and. Ph%i_cm > 0
        if (Ph%i_cm > 0 .and. .not. Keep) call fatal('EMPM_TRANS ... needs complete sets (Ortho = .false., or Ortho with CM_Keep)')
        if (Ph%J_max < 1) call fatal('EMPM_TRANS ... no 1- phonons')
        call EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        if (Nf >= 2) then
            call EMPM_states_read(states_prefix(Params,2),2*Ph%J_max,Ph%N,Ph%N,S2,Lazy=.true.,No_Xt=.true.)
            if (S2%N == 0) call fatal('EMPM_TRANS ... no 2-phonon states in Bin/EMPM_S2_*.bin')
        end if
        call EMPM_states1_make(Ph,S1)
        e_p = real(Params%Calc%A - Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp)/real(Params%Calc%A,dp)

        ! <1- basis||Q_t||0+ basis> on [0 ; 1 ; 2] of both channels (Full_N levels) ...
        Mb = EMPM_cm_dims(Ph,S2,S3,S4,Nf,1,2)
        do it = 1, 2
            call timer_start(t1)
            call EMPM_op_make(Params,Orb,Space,Pairs,Ph,it,Op(it))
            call EMPM_cm_dlevels_nc(Pairs,Ph,Op(it),S1,S2,Nf,0,Dl,Dc,Sel=[0,1,1,2])
            if (Nf >= 2) then
                call load_x(0,1)
                call load_x(1,2)
            end if
            call EMPM_cm_blocks(Ph,Op(it),S1,S2,S3,S4,Dl,Nf,0,1,1,2,Bph,Bhp,Bd,.true.)
            Bt(it)%v = Sg_ph(1)*Bph + Sg_hp(1)*Bhp + Sg_d(1)*Bd
            deallocate(Bph, Bhp, Bd)
            if (Nf >= 2) then
                call EMPM_states_free_X(S2,0,1)
                call EMPM_states_free_X(S2,1,2)
            end if
            call timer_report(t1,'<1-||Q||0+> of the ' // trim(merge('protons ','neutrons',it == 1)) // ' (' // &
                              istr(size(Bt(it)%v,1)) // ' x ' // istr(size(Bt(it)%v,2)) // ')')
            call trans_mem('E1 blocks')
        end do

        if (.not. is_dir(io_path(Params,'/EMPM/Transitions'))) call make_dir(io_path(Params,'/EMPM/Transitions'))
        Labs = [character(len=8) :: 'unproj', 'proj']
        if (Ph%i_cm > 0 .and. .not. Keep) Labs(1) = 'svd'
        do nn = 1, Nf
            do il = 1, 2
                Lab = trim(Labs(il))
                inquire(file=zfile(nn,'0+',Lab), exist=Ex0)
                inquire(file=zfile(nn,'1-',Lab), exist=Ex1)
                if (.not. (Ex0 .and. Ex1)) cycle
                call read_z(nn,'0+',Lab,Nb0,Z0,E0,N0)
                call read_z(nn,'1-',Lab,Nb1,Z1,E1,N1)
                Mb1 = Mb
                Mb1(nn+1:4) = 0
                if (any(Nb1 /= Mb1)) call fatal('EMPM_TRANS ... basis of 1- in ' // zfile(nn,'1-',Lab) // ' does not match')
                r1 = sum(Nb1)
                c0 = sum(Nb0)
                ! ground state: the lowest physical 0+ ...
                g = 1
                if (size(N0) > 0) then
                    do g = 1, size(E0)
                        if (N0(g) < 0.5_dp) exit
                    end do
                    if (g > size(E0)) call fatal('EMPM_TRANS ... no physical 0+ state (' // Lab // ')')
                end if
                n1s = size(E1)
                ! initial state: ig = 1 the correlated g.s. of the same diagonalization, ig = 2 |HF> (0+ basis state 1, n >= 2) ...
                do ig = 1, merge(2,1,nn >= 2 .and. Nb0(0) == 1)
                if (ig == 1) then
                    Psi = Z0(:,g)
                    Eg = E0(g)
                    Ng = merge(N0(min(g,max(size(N0),1))),-1.0_dp,size(N0) > 0)
                    Gdesc = 'g.s. = 0+ state ' // istr(g)
                    Fsuf = ''
                else
                    Psi = [1.0_dp, (0.0_dp, k = 2, c0)]
                    Eg = 0.0_dp
                    Ng = -1.0_dp
                    Gdesc = 'initial state |HF> (0+ basis state 1)'
                    Fsuf = '_fromHF'
                end if
                allocate(Mt(n1s,2))
                do it = 1, 2
                    v = matmul(Bt(it)%v(1:r1,1:c0),Psi)
                    Mt(:,it) = matmul(v,Z1)
                end do
                open(newunit=u, file=io_path(Params,'/EMPM/Transitions/Full' // istr(nn) // '_' // Lab // '_E1' // Fsuf // '.dat'), &
                     status='replace', action='write')
                S = 0.0_dp
                Sp = 0.0_dp
                M1 = 0.0_dp
                do k = 1, n1s
                    B = [Mt(k,1)**2, 0.25_dp*(Mt(k,1) + Mt(k,2))**2, (e_p*Mt(k,1) - e_n*Mt(k,2))**2]
                    S = S + B
                    if (size(N1) == 0) then
                        Sp = Sp + B
                    else if (N1(k) < 0.5_dp) then
                        Sp = Sp + B
                        M1 = M1 + (E1(k) - Eg)*B(3)
                    end if
                end do
                write(u,'(a)') '# EMPM_FULL B(E1; 0+_gs -> 1-_k), up to ' // istr(nn) // ' phonons, ' // Lab // ' ... ' // Gdesc // &
                               ' , E_gs - E_HF = ' // cfmt_f(Eg,0,6) // ' , <N_cm>_gs = ' // &
                               cfmt_e(Ng,0,3) // ' , bases 0+ ' // istr(c0) // ' , 1- ' // istr(r1)
                write(u,'(a)') '# sums B_ph, B_is, B_iv (e^2 fm^2): all 1- states ' // cfmt_f(S(1),0,5) // ' ' // cfmt_f(S(2),0,5) // &
                               ' ' // cfmt_f(S(3),0,5) // ' , physical (<N_cm> < 0.5) ' // cfmt_f(Sp(1),0,5) // ' ' // &
                               cfmt_f(Sp(2),0,5) // ' ' // cfmt_f(Sp(3),0,5) // ' , centroid B_iv (physical) ' // &
                               cfmt_f(M1/max(Sp(3),1.0e-30_dp),0,4) // ' MeV'
                write(u,'(a)') '#      k        E - E_HF          Ex = E - E_gs      N_cm  phys          B_ph            B_is            B_iv'
                do k = 1, n1s
                    B = [Mt(k,1)**2, 0.25_dp*(Mt(k,1) + Mt(k,2))**2, (e_p*Mt(k,1) - e_n*Mt(k,2))**2]
                    write(u,'(i8,2f18.8,es12.3,i6,3es16.7)') k, E1(k), E1(k) - Eg, &
                        merge(N1(min(k,max(size(N1),1))),-1.0_dp,size(N1) > 0), merge(1, 0, size(N1) == 0 .or. N1(min(k,max(size(N1),1))) < 0.5_dp), B
                end do
                close(u)
                write(*,'(a)') TAB // 'up to ' // istr(nn) // ' phonons, ' // Lab // ' ... ' // Gdesc // ' (E - E_HF = ' // &
                               cfmt_f(Eg,0,4) // ') , sum B_ph / B_iv (physical 1-) = ' // cfmt_f(Sp(1),0,4) // ' / ' // &
                               cfmt_f(Sp(3),0,4) // ' , centroid B_iv ' // cfmt_f(M1/max(Sp(3),1.0e-30_dp),0,3) // ' MeV'
                deallocate(Mt)
                end do
            end do
        end do
        call timer_report(t0)
        write(*,'(/,a,/)') 'EMPM transitions have finished ...'

    contains

        function zfile(n,Tg,Lb) result(F)
            integer, intent(in) :: n
            character(len=*), intent(in) :: Tg, Lb
            character(len=:), allocatable :: F
            F = io_path(Params,'/Bin/EMPM_FULL' // istr(n) // '_' // Tg // '_Z_' // Lb // '.bin')
        end function zfile

        ! eigenvectors (Nb, Z), energies & <N_cm> (empty if not written) of one EMPM_FULL diagonalization ...
        subroutine read_z(n,Tg,Lb,Nb,Z,E,Nc)
            integer, intent(in) :: n
            character(len=*), intent(in) :: Tg, Lb
            integer, intent(out) :: Nb(0:4)
            real(dp), allocatable, intent(out) :: Z(:,:), E(:), Nc(:)
            integer :: uz, Ns, Nst
            integer(i8) :: Sz
            character(len=:), allocatable :: Fb
            logical :: Exn
            Fb = io_path(Params,'/Bin/EMPM_FULL' // istr(n) // '_' // Tg // '_')
            inquire(file=Fb // 'E_' // Lb // '.bin', size=Sz)
            Nst = int(Sz/8_i8)
            allocate(E(Nst))
            open(newunit=uz, file=Fb // 'E_' // Lb // '.bin', access='stream', form='unformatted', status='old', action='read')
            read(uz) E
            close(uz)
            open(newunit=uz, file=Fb // 'Z_' // Lb // '.bin', access='stream', form='unformatted', status='old', action='read')
            read(uz) Nb
            Ns = sum(Nb)
            allocate(Z(Ns,Nst))
            read(uz) Z
            close(uz)
            inquire(file=Fb // 'Ncm_' // Lb // '.bin', exist=Exn)
            if (Exn) then
                allocate(Nc(Nst))
                open(newunit=uz, file=Fb // 'Ncm_' // Lb // '.bin', access='stream', form='unformatted', status='old', action='read')
                read(uz) Nc
                close(uz)
            else
                allocate(Nc(0))
            end if
        end subroutine read_z

        subroutine load_x(Jx,Px)
            integer, intent(in) :: Jx, Px
            integer(i8) :: tx
            if (Jx > ubound(S2%N_ch,1)) return
            if (S2%N_ch(Jx,Px) == 0 .or. size(S2%X(Jx,Px)%v) > 0) return
            call timer_start(tx)
            call EMPM_states_load_X(S2,Jx,Px)
            call timer_report(tx,'X of the 2-phonon states J^P = ' // jp_str(Jx,Px) // ' loaded (' // &
                              cfmt_f(8.0_dp*real(size(S2%X(Jx,Px)%v,kind=i8),dp)/1.0e9_dp,0,3) // ' GB)')
        end subroutine load_x

        subroutine trans_mem(Label)
            character(len=*), intent(in) :: Label
            integer :: um, Ios
            character(len=256) :: Line
            integer(i8) :: kb, Rss, Hwm
            Rss = -1
            Hwm = -1
            open(newunit=um, file='/proc/self/status', status='old', action='read', iostat=Ios)
            if (Ios /= 0) return
            do
                read(um,'(a)',iostat=Ios) Line
                if (Ios /= 0) exit
                if (index(Line,'VmRSS:') == 1) then
                    read(Line(7:),*,iostat=Ios) kb
                    if (Ios == 0) Rss = kb
                else if (index(Line,'VmHWM:') == 1) then
                    read(Line(7:),*,iostat=Ios) kb
                    if (Ios == 0) Hwm = kb
                end if
            end do
            close(um)
            if (Rss < 0) return
            write(*,'(a)') TAB // 'memory ... RSS ' // cfmt_f(real(Rss,dp)/976562.5_dp,0,1) // ' GB, peak ' // &
                           cfmt_f(real(Hwm,dp)/976562.5_dp,0,1) // ' GB (' // Label // ')'
            flush(6)
        end subroutine trans_mem
    end subroutine EMPM_TRANS

    ! Strength function B(E1; 0+_gs -> 1-) of the EMPM with up to 3 phonons by Lanczos, without the 3-phonon eigenproblem ...
        ! g.s. = the lowest 0+ of EMPM_FULL with Full_N = 2 (Bin/EMPM_FULL2_0+_Z_unproj.bin), 1- space = [1-phonon ;
        ! 2-phonon states (EMPM_S2) ; redundant 3-phonon basis |a lambda; 1-> of all 2-phonon states a & phonons lambda] ...
        ! a vector of the 3-phonon sector is kept by its coefficients c on the redundant basis (not unique): the equations of
        ! motion give H|i> = sum_k A(i,k)|k>, so H acts as c -> A^T c, and the scalar product is c^T D d with the metric D ...
        ! neither the inverse / Cholesky of D nor the 3-phonon states are needed: one product with A^T & one with D per step ...
        ! couplings in the same (contravariant) form, from the formulas of EMPM_FULL without the projection X^T:
        ! H|a_2> -> sqrt(2J+1) V(:,a) (EMPM_coupling_Vn), H|lambda> -> h20(gamma) on |gamma lambda> (gamma the 0+ 2-phonon
        ! states), Q_ph|a_2; 0+> -> bL(lambda) on |a lambda>; the rest of Q|g> on [1 ; 2] as EMPM_TRANS ...
        ! SF_Mode = 'full': complete sets only (Ortho = .false., all rows: no E3max) ... dense A & D (16 N3^2 bytes) ...
        ! SF_Mode = 'ind': the independent rows of the 3-phonon basis as in EMPM_3PH (also truncated: E2max, E3max),
        ! orthonormalized by the Cholesky factor of D, still without the 3-phonon eigenproblem & the full diagonalization ...
        ! outputs EMPM/SF/SF_E1_<op>.dat (poles), _steps.dat (poles of every step) & _coef.dat (alpha, beta; energies
        ! relative to E_HF), as the NCSM strength function ... SF_Steps, SF_Op = 'iv' / 'is' / 'ph', SF_Reorth ...
    subroutine EMPM_SF(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(empm_space_t) :: Space
        type(fmat_t) :: F_ph
        type(empm_phonons_t) :: Ph
        type(empm_pairs_t) :: Pairs
        type(empm_states2_t) :: S1, S2, S3, S4
        type(empm_cache1_t) :: C1
        type(empm_cm_t) :: Op(2)
        type(rmat_t) :: Dl(3), Dc(3)
        type(empm_basis2_t) :: B3
        type(phonon_block_t) :: Rho
        real(dp), allocatable :: Bph(:,:), Bhp(:,:), Bd(:,:), Bt(:,:), Z0(:,:), E0(:), Psi(:), W(:,:), Vk(:,:), H20(:), H21(:,:)
        real(dp), allocatable :: Am(:,:), Dm(:,:), G1(:,:), Gc(:,:), E1(:), E2(:)
        real(dp), allocatable :: Q12(:,:), Q3(:,:), P3(:,:), w12(:), w3(:), dw3(:), x12(:), x3(:), y12(:), y3(:), px(:), py(:), t3(:)
        real(dp), allocatable :: al(:), be(:), T(:,:), th(:), Uv(:,:)
        integer, allocatable :: A1(:), Rows(:)
        integer :: Nb0(0:4), n1, n2, n12, N3, n20, m, mk, j, k, i, ia, a, ll, it, g, ps, u, l0
        real(dp) :: e_p, e_n, cf(2), Eg, S0, s, b2, hxy, hyx, m0, m1, eu
        integer(i8) :: t0, t1
        character(len=:), allocatable :: Nm, Fl
        logical :: Reo, Prj, Full
        character(len=:), allocatable :: DimS
        real(dp), allocatable :: Hi(:,:), Dii(:,:), Dr(:,:), Dtr(:,:), Cv(:,:), Yv(:,:), wy(:)
        integer, allocatable :: Ind(:)
        integer(i8), allocatable :: Dnp(:,:)
        integer :: r, Tgt, N_visit, Info, nc
        real(dp) :: Asym
        interface
            subroutine dgemv(Trans,m,n,Alpha,A,Lda,x,Incx,Beta,y,Incy)
                import :: dp
                character :: Trans
                integer :: m, n, Lda, Incx, Incy
                real(dp) :: Alpha, Beta, A(Lda,*), x(*), y(*)
            end subroutine dgemv
            subroutine dsymv(Uplo,n,Alpha,A,Lda,x,Incx,Beta,y,Incy)
                import :: dp
                character :: Uplo
                integer :: n, Lda, Incx, Incy
                real(dp) :: Alpha, Beta, A(Lda,*), x(*), y(*)
            end subroutine dsymv
            subroutine dpotrf(Uplo,n,A,Lda,Info)
                import :: dp
                character :: Uplo
                integer :: n, Lda, Info
                real(dp) :: A(Lda,*)
            end subroutine dpotrf
            subroutine dsygst(Itype,Uplo,n,A,Lda,B,Ldb,Info)
                import :: dp
                character :: Uplo
                integer :: Itype, n, Lda, Ldb, Info
                real(dp) :: A(Lda,*), B(Ldb,*)
            end subroutine dsygst
            subroutine dtrsm(Side,Uplo,Transa,Diag,m,n,Alpha,A,Lda,B,Ldb)
                import :: dp
                character :: Side, Uplo, Transa, Diag
                integer :: m, n, Lda, Ldb
                real(dp) :: Alpha, A(Lda,*), B(Ldb,*)
            end subroutine dtrsm
        end interface

        call wigner_init_float(75,'Jmax',9)
        !$ call omp_set_max_active_levels(1)
        m = Params%Calc%EMPM%SF_Steps
        Nm = trim(adjustl(Params%Calc%EMPM%SF_Op))
        Reo = Params%Calc%EMPM%SF_Reorth
        select case (trim(adjustl(Params%Calc%EMPM%SF_Mode)))
        case ('full')
            Full = .true.
        case ('ind')
            Full = .false.
        case default
            call fatal('EMPM_SF ... SF_Mode must be ''full'' or ''ind''')
        end select
        write(*,'(a)') 'Starting the EMPM strength function B_' // Nm // '(E1; 0+_gs -> 1-) with up to 3 phonons by Lanczos (' // &
                       istr(m) // ' steps), ' // trim(merge('redundant 3-phonon basis              ','independent rows of the 3-phonon basis',Full)) // ' ...'
        write(*,'(a)') trim(Params%Calc%Path)
        if (m < 1) call fatal('EMPM_SF ... SF_Steps must be positive')
        e_p = real(Params%Calc%A - Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        e_n = real(Params%Calc%Z,dp)/real(Params%Calc%A,dp)
        select case (Nm)
        case ('iv')
            cf = [e_p, -e_n]
        case ('is')
            cf = [0.5_dp, 0.5_dp]
        case ('ph')
            cf = [1.0_dp, 0.0_dp]
        case default
            call fatal('EMPM_SF ... SF_Op must be ''iv'', ''is'' or ''ph''')
        end select
        call timer_start(t0)
        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call EMPM_space_make(Params,Orb,Space)
        call Fmat_import(Space,io_path(Params,'/Bin/EMPM_F.bin'),F_ph)
        call EMPM_phonons_import(Params,Space,F_ph,Ph)
        if (Ph%i_cm > 0) call fatal('EMPM_SF ... needs the complete sets of Ortho = .false.')
        if (Ph%J_max < 1) call fatal('EMPM_SF ... no 1- phonons')
        call EMPM_pairs_make(Space,F_ph,Ph,Pairs)
        call EMPM_states_read(states_prefix(Params,2),2*Ph%J_max,Ph%N,Ph%N,S2)
        if (S2%N == 0) call fatal('EMPM_SF ... no 2-phonon states in Bin/EMPM_S2_*.bin, run EMPM_2PH with Save_X first')
        call EMPM_states1_make(Ph,S1)
        n1 = Ph%N_ch(1,2)
        l0 = Ph%First(1,2)
        n2 = 0
        if (ubound(S2%N_ch,1) >= 1) n2 = S2%N_ch(1,2)
        n20 = S2%N_ch(0,1)
        n12 = n1 + n2
        call EMPM_basis3_make(Ph,S2,1,2,merge(huge(1.0_dp),Params%Calc%EMPM%E3max,Full),'all',B3)
        N3 = B3%N_full
        ! SF_Qrows (SF_Mode = 'ind') ... the rows |a lambda> of the 0+ 2-phonon states a & the 1- phonons lambda are kept
            ! whatever their energy: they carry Q_ph acting on the 2-phonon part of the g.s., i.e. the g.s. correlations
            ! of the 1-phonon states; without them the 1- states miss the correlation energy the g.s. has ...
        if (Params%Calc%EMPM%SF_Qrows .and. .not. Full) then
            block
                logical, allocatable :: Kp(:)
                integer :: kq
                allocate(Kp(N3))
                Kp = .false.
                Kp(B3%Row) = .true.
                kq = count(Kp)
                do k = 1, N3
                    if (S2%J(B3%Al(k)) == 0 .and. S2%P(B3%Al(k)) == 1 .and. Ph%J(B3%La(k)) == 1 .and. Ph%P(B3%La(k)) == 2) Kp(k) = .true.
                end do
                B3%Row = pack([(k, k = 1, N3)], Kp)
                B3%N_tr = size(B3%Row)
                write(*,'(a)') TAB // 'SF_Qrows ... ' // istr(B3%N_tr - kq) // ' rows |0+ 2-phonon state x 1- phonon> added to the ' // &
                               istr(kq) // ' rows below E3max'
            end block
        end if
        write(*,'(a)') TAB // '1- space ... 1-phonon ' // istr(n1) // ' , 2-phonon ' // istr(n2) // ' , redundant 3-phonon basis ' // &
                       istr(N3) // ' (' // istr(S2%N) // ' 2-phonon states x ' // istr(Ph%N) // ' phonons)' // &
                       trim(merge(' , A & D (GB) ' // cfmt_f(16.0_dp*real(N3,dp)**2/1.0e9_dp,12,3), &
                                  ' , rows       ' // cfmt_d(B3%N_tr,12), Full))

        if (.not. Full) write(*,'(a)') TAB // 'memory ... D of the rows ' // cfmt_f(8.0_dp*real(B3%N_tr,dp)**2/1.0e9_dp,0,1) // &
                                       ' GB , D rows on the full basis up to ' // cfmt_f(8.0_dp*real(B3%N_tr,dp)*real(N3,dp)/1.0e9_dp,0,1) // &
                                       ' GB (if all rows are independent) , distinct 2-phonon states of the rows ' // &
                                       istr(nrow_states())
        if (Params%Calc%EMPM%Dry_Run) return
        ! ground state: the lowest 0+ of EMPM_FULL with up to 2 phonons, basis [HF ; 0+ phonons ; 0+ 2-phonon states] ...
        call read_gs()
        if (Nb0(0) /= 1 .or. Nb0(1) /= Ph%N_ch(0,1) .or. Nb0(2) /= n20 .or. any(Nb0(3:4) /= 0)) &
            call fatal('EMPM_SF ... the 0+ basis of Bin/EMPM_FULL2_0+_Z_unproj.bin does not match (EMPM_FULL with Full_N = 2)')
        write(*,'(a)') TAB // 'g.s. = 0+ state ' // istr(g) // ' of EMPM_FULL (<= 2 phonons) , E - E_HF = ' // cfmt_f(Eg,0,6) // &
                       ' MeV , weights 0 / 1 / 2 phonons = ' // cfmt_f(Psi(1)**2,0,4) // ' / ' // &
                       cfmt_f(sum(Psi(2:1+Nb0(1))**2),0,4) // ' / ' // cfmt_f(sum(Psi(2+Nb0(1):)**2),0,4)

        ! Q|g> ... [1 ; 2] by the blocks of EMPM_TRANS, the 3-phonon part = Q_ph on the 2-phonon components of g ...
        call timer_start(t1)
        allocate(x12(n12), x3(N3))
        x12 = 0.0_dp
        x3 = 0.0_dp
        do it = 1, 2
            call EMPM_op_make(Params,Orb,Space,Pairs,Ph,it,Op(it))
            if (cf(it) == 0.0_dp) cycle
            call EMPM_cm_dlevels_nc(Pairs,Ph,Op(it),S1,S2,2,0,Dl,Dc,Sel=[0,1,1,2])
            call EMPM_cm_blocks(Ph,Op(it),S1,S2,S3,S4,Dl,2,0,1,1,2,Bph,Bhp,Bd,.true.)
            Bt = Sg_ph(1)*Bph + Sg_hp(1)*Bhp + Sg_d(1)*Bd
            if (size(Bt,1) /= n12 .or. size(Bt,2) /= size(Psi)) call fatal('EMPM_SF ... operator blocks do not match the bases')
            x12 = x12 + cf(it)*matmul(Bt,Psi)
            do ia = 1, n20
                a = S2%First(0,1) + ia - 1
                do ll = l0, l0 + n1 - 1
                    k = B3%K(a,ll)
                    if (k > 0) x3(k) = x3(k) + cf(it)*Sg_ph(1)*Op(it)%bL(ll)*Psi(1+Nb0(1)+ia)
                end do
            end do
            deallocate(Bph, Bhp, Bd, Bt)
        end do
        call timer_report(t1,'Q|g> on the 1- space')

        ! Hamiltonian blocks ...
        call timer_start(t1)
        call EMPM_cache1_make(Pairs,Ph,C1)
        call EMPM_phonon_W(Space,F_ph,Pairs,Ph,W)
        E1 = Ph%E(l0:l0+n1-1)
        allocate(E2(n2), H21(n2,n1), H20(n20))
        if (n2 > 0) E2 = S2%E(S2%First(1,2):S2%First(1,2)+n2-1)
        ! 2 <- 0 (0+) & 2 <- 1 ...
        allocate(A1(0))
        call EMPM_coupling_V(Space,F_ph,Pairs,Ph,S2%Bs(0,1),A1,1,Vk)
        H20 = matmul(transpose(S2%X(0,1)%v),Vk(:,1))
        deallocate(A1, Vk)
        H21 = 0.0_dp
        if (n2 > 0) then
            A1 = [(l0 + i - 1, i = 1, n1)]
            call EMPM_coupling_V(Space,F_ph,Pairs,Ph,S2%Bs(1,2),A1,0,Vk)
            H21 = mtm(S2%X(1,2)%v(:,1:n2),Vk)
            deallocate(Vk)
        end if
        ! 3 <- 1: H|lambda> = sum_gamma h20(gamma) |gamma lambda> ...
        allocate(G1(N3,n1))
        G1 = 0.0_dp
        do i = 1, n1
            do k = 1, n20
                j = B3%K(S2%First(0,1)+k-1,l0+i-1)
                if (j > 0) G1(j,i) = H20(k)
            end do
        end do
        ! 3 <- 2: H|a> = sqrt(2J+1) sum_k V(k,a) |k> ...
        allocate(Gc(N3,n2))
        !$omp parallel do schedule(dynamic) private(Rho)
        do ia = 1, n2
            call EMPM_densityN(Pairs,Ph,C1,C1%R,S1,S2,S2%First(1,2)+ia-1,Rho)
            call EMPM_coupling_Vn(Pairs,Ph,W,Rho,1,S2%J_max,S2%N_ch,S2%First,B3,Gc(:,ia))
        end do
        !$omp end parallel do
        Gc = sqrt(3.0_dp)*Gc
        call timer_report(t1,'couplings 2 <- 1, 3 <- 1 & 3 <- 2')
        if (.not. Full) then
            ! Independent rows I of the 3-phonon basis (pivoted Cholesky of D, or the greedy selection with Ind_Greedy),
                ! H_II = A_I D_I^T & D_II on them as in EMPM_3PH, but no eigenproblem: D_II = U^T U, the orthonormal states
                ! e = sum_i |i> (U^-1)_i, H~ = U^-T H_II U^-1, and a vector with the coefficients g on the full basis has the
                ! components U^-T D_I g ... Lanczos on [1 ; 2 ; e] with the Euclidean metric ... also for truncated bases
                ! (E2max, E3max), with the approximations of the EMPM for incomplete sets ...
            call timer_start(t1)
            allocate(Dnp(0:3*Ph%J_max+1,2))
            Dnp(:,:) = EMPM_npnh_dims(Space,3,3*Ph%J_max + 1)
            Tgt = int(min(int(B3%N_tr,i8), Dnp(1,2)))
            if (Params%Calc%EMPM%Ind_Greedy) then
                call EMPM_ind_greedy_H(Pairs,Ph,C1,C1%R,S1,S2,B3,Tgt,Params%Calc%EMPM%Mem_GB,Ind,Dii,Hi,Dr,N_visit)
                write(*,'(a)') TAB // 'independent rows (greedy) ... ' // istr(size(Ind)) // ' of ' // istr(B3%N_tr) // &
                               ' , 3p3h dimension of 1- ' // istr(int(Dnp(1,2))) // ' , rows visited ' // istr(N_visit)
            else
                call EMPM_D3_build(Pairs,Ph,C1,C1%R,S1,S2,B3,Params%Calc%EMPM%Mem_GB,Dtr)
                call EMPM_chol_ind(Dtr,Ind,Dii)
                deallocate(Dtr)
                call EMPM_H3_ind(Pairs,Ph,C1,C1%R,S1,S2,B3,Ind,Params%Calc%EMPM%Mem_GB,Hi,Dr)
                write(*,'(a)') TAB // 'independent rows (pivoted Cholesky) ... ' // istr(size(Ind)) // ' of ' // istr(B3%N_tr) // &
                               ' , 3p3h dimension of 1- ' // istr(int(Dnp(1,2)))
            end if
            r = size(Ind)
            call timer_report(t1,'independent rows, A D^T & D on them')
            call timer_start(t1)
            Asym = 0.0_dp
            !$omp parallel do schedule(dynamic) reduction(max:Asym) private(i)
            do j = 1, r
                do i = 1, j - 1
                    Asym = max(Asym, abs(Hi(i,j) - Hi(j,i)))
                    Hi(i,j) = 0.5_dp*(Hi(i,j) + Hi(j,i))
                end do
            end do
            !$omp end parallel do
            write(*,'(a)') TAB // '|AD - (AD)^T|max = ' // cfmt_e(Asym,0,3)
            call dpotrf('U',r,Dii,r,Info)
            if (Info /= 0) call fatal('EMPM_SF ... D of the independent rows is not positive definite, info = ' // istr(Info))
            call dsygst(1,'U',r,Hi,r,Dii,r,Info)
            if (Info /= 0) call fatal('EMPM_SF ... dsygst failed, info = ' // istr(Info))
            ! couplings & start vector: [G1 | Gc | x3] -> U^-T Dr^T [..] ...
            nc = n12 + 1
            allocate(Cv(N3,nc), Yv(r,nc))
            Cv(:,1:n1) = G1
            if (n2 > 0) Cv(:,n1+1:n12) = Gc
            Cv(:,nc) = x3
            call dgemm('T','N',r,nc,N3,1.0_dp,Dr,N3,Cv,N3,0.0_dp,Yv,r)
            call dtrsm('L','U','T','N',r,nc,1.0_dp,Dii,r,Yv,r)
            deallocate(Cv, Dr, Dii, G1, Gc)
            call timer_report(t1,'Cholesky of D, H & the couplings on the orthonormal states (' // istr(r) // ')')
            S0 = dot_product(x12,x12) + dot_product(Yv(:,nc),Yv(:,nc))
            write(*,'(a)') TAB // 'total strength S0 = |Q g|^2 = ' // cfmt_e(S0,0,10) // ' e^2 fm^2 (1- / 2- / 3-phonon parts ' // &
                           cfmt_f(sum(x12(1:n1)**2)/S0,0,6) // ' / ' // cfmt_f(sum(x12(n1+1:)**2)/S0,0,6) // ' / ' // &
                           cfmt_f(dot_product(Yv(:,nc),Yv(:,nc))/S0,0,6) // ')'
            call make_dir(io_path(Params,'/EMPM/SF'))
            Fl = io_path(Params,'/EMPM/SF/SF_E1_' // Nm // '.dat')
            DimS = istr(n1) // ' + ' // istr(n2) // ' + ' // istr(r) // ' (independent rows of ' // istr(B3%N_tr) // ')'
            allocate(Q12(n12,m+1), Q3(r,m+1), al(m), be(m), w12(n12), wy(r))
            Q12(:,1) = x12/sqrt(S0)
            Q3(:,1) = Yv(:,nc)/sqrt(S0)
            mk = 0
            do j = 1, m
                call timer_start(t1)
                w12(1:n1) = E1*Q12(1:n1,j)
                if (n2 > 0) then
                    w12(1:n1) = w12(1:n1) + matmul(Q12(n1+1:n12,j),H21)
                    w12(n1+1:n12) = E2*Q12(n1+1:n12,j) + matmul(H21,Q12(1:n1,j))
                end if
                call dgemv('T',r,n12,1.0_dp,Yv,r,Q3(:,j),1,1.0_dp,w12,1)
                call dsymv('U',r,1.0_dp,Hi,r,Q3(:,j),1,0.0_dp,wy,1)
                call dgemv('N',r,n12,1.0_dp,Yv,r,Q12(:,j),1,1.0_dp,wy,1)
                al(j) = dot_product(Q12(:,j),w12) + dot_product(Q3(:,j),wy)
                w12 = w12 - al(j)*Q12(:,j)
                wy = wy - al(j)*Q3(:,j)
                if (j > 1) then
                    w12 = w12 - be(j-1)*Q12(:,j-1)
                    wy = wy - be(j-1)*Q3(:,j-1)
                end if
                if (Reo) then
                    do ps = 1, 2
                        do k = 1, j
                            s = dot_product(Q12(:,k),w12) + dot_product(Q3(:,k),wy)
                            w12 = w12 - s*Q12(:,k)
                            wy = wy - s*Q3(:,k)
                        end do
                    end do
                end if
                be(j) = sqrt(dot_product(w12,w12) + dot_product(wy,wy))
                mk = j
                call sf_write()
                call timer_report(t1,'  SF Lanczos step ' // istr(j) // ' of ' // istr(m) // ' , alpha = ' // cfmt_f(al(j),0,6) // &
                                  ' , beta = ' // cfmt_f(be(j),0,6))
                if (j == m .or. be(j) < 1.0e-9_dp*max(maxval(abs(al(1:j))),1.0_dp)) exit
                Q12(:,j+1) = w12/be(j)
                Q3(:,j+1) = wy/be(j)
            end do
        else
            if (B3%N_tr /= N3) call fatal('EMPM_SF ... incomplete rows of the 3-phonon basis')
            ! A & D of the redundant basis: Am(:,i) = row i of A, Dm(:,i) = row i of D ...
            call timer_start(t1)
            allocate(Am(N3,N3), Dm(N3,N3), Rows(N3))
            Rows = [(i, i = 1, N3)]
            call EMPM_AD3_rows(Pairs,Ph,C1,C1%R,S1,S2,B3,Rows,Dm,Am)
            call timer_report(t1,'A & D of the redundant 3-phonon basis (' // istr(N3) // ' rows)')

            ! Check: <x|H y> = <H x|y> for two random vectors ...
            allocate(y12(n12), y3(N3), w12(n12), w3(N3), dw3(N3), px(N3), py(N3), t3(N3))
            block
                real(dp), allocatable :: a12(:), a3(:), b12(:), b3(:), r12(:)
                allocate(a12(n12), a3(N3), b12(n12), b3(N3), r12(n12))
                call random_number(y12)
                call random_number(y3)
                call random_number(r12)
                call random_number(t3)
                y12 = y12 - 0.5_dp
                y3 = y3 - 0.5_dp
                r12 = r12 - 0.5_dp
                t3 = t3 - 0.5_dp
                call dmul(t3,px)
                call dmul(y3,py)
                call hmul(y12,y3,py,a12,a3)
                call hmul(r12,t3,px,b12,b3)
                hxy = dot_product(r12,a12) + dot_product(px,a3)
                hyx = dot_product(y12,b12) + dot_product(py,b3)
                write(*,'(a)') TAB // 'Hermiticity check ... <x|H y> = ' // cfmt_e(hxy,0,12) // ' , <H x|y> = ' // cfmt_e(hyx,0,12) // &
                               ' , relative difference ' // cfmt_e(abs(hxy - hyx)/max(abs(hxy),1.0e-300_dp),0,3)
            end block

            ! Lanczos ... vectors q_j = [q12 ; c3] with p_j = D c3 ...
            ! complete sets: |a lambda> = O+_lambda |a> with charge-conserving phonons, so sum_i |i><i| = N_p(pi) N_h(pi) +
                ! N_p(nu) N_h(nu) on the 3p3h space: the metric D has only the eigenvalues 0 (redundancy), 5 (2 + 1 proton /
                ! neutron ph pairs) & 9 (3 + 0), and P = D (14 - D) / 45 is the projector on its range ... the null-space part
                ! of c (which the metric does not see, but which grows from step to step and ruins the precision) is removed by
                ! c -> P c in every step, for one more product with D ...
            call dmul(t3,px)
            call dmul(px,py)
            call dmul(py,y3)
            eu = sqrt(sum((y3 - 14.0_dp*py + 45.0_dp*px)**2)/max(sum(y3**2),1.0e-300_dp))
            Prj = eu < 1.0e-8_dp
            write(*,'(a)') TAB // 'metric ... D (D - 5) (D - 9) = 0 on a random vector to ' // cfmt_e(eu,0,3) // ' (relative)'
            if (Prj) then
                write(*,'(a)') TAB // 'metric ... the null space of D is projected out in every step, P = D (14 - D) / 45'
            else
                write(*,'(a)') TAB // 'Warning ... D does not have the eigenvalues 0, 5 & 9 only (incomplete sets?), no projection ' // &
                               'of the null space: the precision is lost after some tens of steps'
            end if
            call dmul(x3,px)
            if (Prj) then
                call dmul(px,py)
                x3 = (14.0_dp*px - py)/45.0_dp
            end if
            S0 = dot_product(x12,x12) + dot_product(x3,px)
            write(*,'(a)') TAB // 'total strength S0 = |Q g|^2 = ' // cfmt_e(S0,0,10) // ' e^2 fm^2 (1- / 2- / 3-phonon parts ' // &
                           cfmt_f(sum(x12(1:n1)**2)/S0,0,6) // ' / ' // cfmt_f(sum(x12(n1+1:)**2)/S0,0,6) // ' / ' // &
                           cfmt_f(dot_product(x3,px)/S0,0,6) // ')'
            call make_dir(io_path(Params,'/EMPM/SF'))
            Fl = io_path(Params,'/EMPM/SF/SF_E1_' // Nm // '.dat')
            DimS = istr(n1) // ' + ' // istr(n2) // ' + ' // istr(N3) // ' (redundant)'
            allocate(Q12(n12,m+1), Q3(N3,m+1), P3(N3,m+1), al(m), be(m))
            Q12(:,1) = x12/sqrt(S0)
            Q3(:,1) = x3/sqrt(S0)
            P3(:,1) = px/sqrt(S0)
            mk = 0
            do j = 1, m
                call timer_start(t1)
                call hmul(Q12(:,j),Q3(:,j),P3(:,j),w12,w3)
                call dmul(w3,dw3)
                if (Prj) then
                    call dmul(dw3,t3)
                    w3 = (14.0_dp*dw3 - t3)/45.0_dp
                end if
                al(j) = dot_product(Q12(:,j),w12) + dot_product(P3(:,j),w3)
                call axpy3(-al(j),j)
                if (j > 1) call axpy3(-be(j-1),j-1)
                if (Reo) then
                    do ps = 1, 2
                        do k = 1, j
                            s = dot_product(Q12(:,k),w12) + dot_product(P3(:,k),w3)
                            call axpy3(-s,k)
                        end do
                    end do
                end if
                b2 = dot_product(w12,w12) + dot_product(w3,dw3)
                be(j) = sqrt(max(b2,0.0_dp))
                mk = j
                call sf_write()
                eu = sqrt(dot_product(Q3(:,j),Q3(:,j)))
                call timer_report(t1,'  SF Lanczos step ' // istr(j) // ' of ' // istr(m) // ' , alpha = ' // cfmt_f(al(j),0,6) // &
                                  ' , beta = ' // cfmt_f(be(j),0,6) // ' , |c3| = ' // cfmt_e(eu,0,3) // ' , <c3|D|c3> = ' // &
                                  cfmt_e(dot_product(Q3(:,j),P3(:,j)),0,3))
                if (j == m .or. be(j) < 1.0e-9_dp*max(maxval(abs(al(1:j))),1.0_dp)) exit
                Q12(:,j+1) = w12/be(j)
                Q3(:,j+1) = w3/be(j)
                P3(:,j+1) = dw3/be(j)
            end do
        end if
        write(*,'(a)') TAB // 'B_' // Nm // '(E1): ' // istr(mk) // ' Lanczos steps, m0 = ' // cfmt_e(m0,0,6) // ' , centroid ' // &
                       cfmt_f(m1/m0,0,4) // ' MeV ... ' // Fl
        call timer_report(t0)
        write(*,'(/,a,/)') 'EMPM strength function has finished ...'

    contains

        ! number of distinct 2-phonon states a among the rows ...
        integer function nrow_states()
            logical, allocatable :: Us(:)
            allocate(Us(S2%N))
            Us = .false.
            Us(B3%Al(B3%Row)) = .true.
            nrow_states = count(Us)
        end function nrow_states

        ! p = D c ...
        subroutine dmul(c,p)
            real(dp), intent(in) :: c(:)
            real(dp), intent(out) :: p(:)
            call dgemv('T',N3,N3,1.0_dp,Dm,N3,c,1,0.0_dp,p,1)
        end subroutine dmul

        ! [o12 ; o3] = H [v12 ; c3] with p = D c3 ...
        subroutine hmul(v12,c3,p,o12,o3)
            real(dp), intent(in) :: v12(:), c3(:), p(:)
            real(dp), intent(out) :: o12(:), o3(:)
            o12(1:n1) = E1*v12(1:n1)
            call dgemv('T',N3,n1,1.0_dp,G1,N3,p,1,1.0_dp,o12,1)
            call dgemv('N',N3,N3,1.0_dp,Am,N3,c3,1,0.0_dp,o3,1)
            call dgemv('N',N3,n1,1.0_dp,G1,N3,v12,1,1.0_dp,o3,1)
            if (n2 > 0) then
                o12(1:n1) = o12(1:n1) + matmul(v12(n1+1:n12),H21)
                o12(n1+1:n12) = E2*v12(n1+1:n12) + matmul(H21,v12(1:n1))
                call dgemv('T',N3,n2,1.0_dp,Gc,N3,p,1,1.0_dp,o12(n1+1:n12),1)
                call dgemv('N',N3,n2,1.0_dp,Gc,N3,v12(n1+1:n12),1,1.0_dp,o3,1)
            end if
        end subroutine hmul

        ! w += s q_k (with D w) ...
        subroutine axpy3(s_,k_)
            real(dp), intent(in) :: s_
            integer, intent(in) :: k_
            integer :: i_
            w12 = w12 + s_*Q12(:,k_)
            !$omp parallel do
            do i_ = 1, N3
                w3(i_) = w3(i_) + s_*Q3(i_,k_)
                dw3(i_) = dw3(i_) + s_*P3(i_,k_)
            end do
            !$omp end parallel do
        end subroutine axpy3

        ! g.s. from Bin/EMPM_FULL2_0+_{E,Z,Ncm}_unproj.bin: the lowest 0+ (with <N_cm> < 0.5 when <N_cm> is there) ...
        subroutine read_gs()
            integer :: uz, Nst
            integer(i8) :: Sz
            real(dp), allocatable :: Nc(:)
            character(len=:), allocatable :: Fb
            logical :: Ex
            Fb = io_path(Params,'/Bin/EMPM_FULL2_0+_')
            inquire(file=Fb // 'Z_unproj.bin', exist=Ex)
            if (.not. Ex) call fatal('EMPM_SF ... no ' // Fb // 'Z_unproj.bin, run EMPM_FULL with Full_N = 2 first')
            inquire(file=Fb // 'E_unproj.bin', size=Sz)
            Nst = int(Sz/8_i8)
            allocate(E0(Nst))
            open(newunit=uz, file=Fb // 'E_unproj.bin', access='stream', form='unformatted', status='old', action='read')
            read(uz) E0
            close(uz)
            open(newunit=uz, file=Fb // 'Z_unproj.bin', access='stream', form='unformatted', status='old', action='read')
            read(uz) Nb0
            allocate(Z0(sum(Nb0),Nst))
            read(uz) Z0
            close(uz)
            g = 1
            inquire(file=Fb // 'Ncm_unproj.bin', exist=Ex)
            if (Ex) then
                allocate(Nc(Nst))
                open(newunit=uz, file=Fb // 'Ncm_unproj.bin', access='stream', form='unformatted', status='old', action='read')
                read(uz) Nc
                close(uz)
                do g = 1, Nst
                    if (Nc(g) < 0.5_dp) exit
                end do
                if (g > Nst) g = 1
            end if
            Psi = Z0(:,g)
            Eg = E0(g)
        end subroutine read_gs

        ! poles of the mk steps (file Fl), their history (_steps.dat) & the Lanczos coefficients (_coef.dat) ...
        subroutine sf_write()
            integer :: i_, k_, uh
            character(len=:), allocatable :: Fh
            Fh = Fl(1:len(Fl)-4)
            allocate(T(mk,mk), th(mk), Uv(mk,mk))
            T = 0.0_dp
            do i_ = 1, mk
                T(i_,i_) = al(i_)
                if (i_ < mk) then
                    T(i_,i_+1) = be(i_)
                    T(i_+1,i_) = be(i_)
                end if
            end do
            call eig_sym(T,th,Uv)
            open(newunit=u, file=Fl, status='replace', action='write')
            write(u,'(a)') '# EMPM strength function B_' // Nm // '(E1) from the 0+ g.s. (<= 2 phonons), 1- with <= 3 phonons, E_g - E_HF = ' // &
                           cfmt_f(Eg,0,10) // ', dim = ' // DimS // ', Lanczos iterations ' // istr(mk)
            write(u,'(a,es16.8)') '# total strength S0 = |<1-||Q||g>|^2 (e^2 fm^2) = ', S0
            write(u,'(a)') '# k   E - E_HF (MeV)   Ex = E - E_g   B   residual'
            m0 = 0.0_dp
            m1 = 0.0_dp
            do k_ = 1, mk
                write(u,'(i6,2f16.8,es18.8e3,es12.3e3)') k_, th(k_), th(k_) - Eg, S0*Uv(1,k_)**2, be(mk)*abs(Uv(mk,k_))
                m0 = m0 + S0*Uv(1,k_)**2
                m1 = m1 + S0*Uv(1,k_)**2*(th(k_) - Eg)
            end do
            close(u)
            if (mk == 1) then
                open(newunit=uh, file=Fh // '_steps.dat', status='replace', action='write')
                write(uh,'(a)') '# EMPM strength function B_' // Nm // '(E1) poles after every Lanczos step, E_g - E_HF = ' // &
                                cfmt_f(Eg,0,10) // ', dim = ' // DimS
                write(uh,'(a,es16.8)') '# total strength S0 = ', S0
                write(uh,'(a)') '# step   k   E - E_HF (MeV)   Ex = E - E_g   B   residual'
            else
                open(newunit=uh, file=Fh // '_steps.dat', status='old', position='append', action='write')
                write(uh,'(a)') ''
            end if
            do k_ = 1, mk
                write(uh,'(i6,i6,2f16.8,es18.8e3,es12.3e3)') mk, k_, th(k_), th(k_) - Eg, S0*Uv(1,k_)**2, be(mk)*abs(Uv(mk,k_))
            end do
            close(uh)
            open(newunit=uh, file=Fh // '_coef.dat', status='replace', action='write')
            write(uh,'(a,es24.16,a)') '# Lanczos coefficients of B_' // Nm // '(E1), E_g - E_HF = ', Eg, ', dim = ' // DimS
            write(uh,'(a,es24.16)') '# S0 = ', S0
            write(uh,'(a)') '# j   alpha_j - E_HF (MeV)   beta_j (MeV)'
            do i_ = 1, mk
                write(uh,'(i6,2es24.16)') i_, al(i_), be(i_)
            end do
            close(uh)
            deallocate(T, th, Uv)
        end subroutine sf_write
    end subroutine EMPM_SF

    ! Folder of the n-phonon state files EMPM_S<n>_<JP>.bin ... States_Dir (a leading $VARIABLE expanded) or Bin of the
        ! output folder ...
    function states_dir(Params) result(Dir)
        type(params_t), intent(in) :: Params
        character(len=:), allocatable :: Dir
        character(len=1024) :: Val
        integer :: k, L

        Dir = trim(adjustl(Params%Calc%EMPM%States_Dir))
        if (len(Dir) == 0) then
            Dir = io_path(Params,'/Bin')
            return
        end if
        if (Dir(1:1) == '$') then
            k = scan(Dir,'/')
            if (k == 0) k = len(Dir) + 1
            call get_environment_variable(Dir(2:k-1),Val,L)
            if (L == 0) call fatal('States_Dir ... environment variable ' // Dir(2:k-1) // ' is not set')
            Dir = trim(Val) // Dir(k:)
        end if
    end function states_dir

    function states_prefix(Params,n) result(Prefix)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: n
        character(len=:), allocatable :: Prefix
        Prefix = states_dir(Params) // '/EMPM_S' // istr(n) // '_'
    end function states_prefix

    function states_file(Params,n,Tag) result(File)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: n
        character(len=*), intent(in) :: Tag
        character(len=:), allocatable :: File
        File = states_prefix(Params,n) // Tag // '.bin'
    end function states_file

    function jp_str(J,P) result(s)
        integer, intent(in) :: J, P
        character(len=:), allocatable :: s
        s = istr(J) // merge('+', '-', P == 1)
    end function jp_str

    ! Channels (J,P) of the JP2 list ... all J = 0 .. J_max of both parities for an empty list ...
    function channels(JP2,J_max) result(Ch)
        character(len=*), intent(in) :: JP2(:)
        integer, intent(in) :: J_max
        integer, allocatable :: Ch(:,:)
        integer :: k, n, J, Ios
        character(len=:), allocatable :: s

        n = count(len_trim(JP2) > 0)
        if (n == 0) then
            allocate(Ch(2,2*(J_max + 1)))
            n = 0
            do J = 0, J_max
                Ch(:,n+1) = [J, 1]
                Ch(:,n+2) = [J, 2]
                n = n + 2
            end do
            return
        end if
        allocate(Ch(2,n))
        n = 0
        do k = 1, size(JP2)
            s = trim(adjustl(JP2(k)))
            if (len(s) == 0) cycle
            n = n + 1
            read(s(1:len(s)-1),*,iostat=Ios) J
            if (Ios /= 0 .or. index('+-',s(len(s):len(s))) == 0) call fatal('Invalid EMPM channel JP2 = ' // s)
            Ch(:,n) = [J, merge(1, 2, s(len(s):len(s)) == '+')]
        end do
    end function channels

end module nj_empm_driver
