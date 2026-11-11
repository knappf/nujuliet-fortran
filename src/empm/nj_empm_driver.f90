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
    use nj_linalg, only : mm, mtm
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

    public :: EMPM_TDA, EMPM_2PH, EMPM_3PH, EMPM_4PH, EMPM_FULL

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

        Orb = orbitals_make(Params)
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

        Orb = orbitals_make(Params)
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

        Orb = orbitals_make(Params)
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

        Orb = orbitals_make(Params)
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

        Orb = orbitals_make(Params)
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
            call EMPM_states_read(states_prefix(Params,2),2*Ph%J_max,Ph%N,Ph%N,S2,No_Xt=.true.)
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

        write(*,'(/,a)') 'Densities, redefined interactions & the couplings to the ground state ...'
        call timer_start(t1)
        call EMPM_cache1_make(Pairs,Ph,C1)
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
            call EMPM_cm_dlevels(Pairs,Ph,C1,Cm,S1,S2,S3,merge(Nf,Nf-1,CMn),Dl)
            call timer_report(t1,'within-level b_cm matrices of the 1- .. ' // istr(min(3,merge(Nf,Nf-1,CMn))) // &
                              '-phonon states')
            ! top-level completion of N_cm (level Full_N + 1 not in the basis): rank-0 within-level matrices of the commutator
                ! C = sum_mu [L_ph+, L_ph] up to level Full_N (EMPM_cm_top) ...
            if (CMn) then
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
            H20 = matmul(transpose(S2%X(0,1)%v),Vk(:,1))
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
                    call EMPM_cm_blocks(Ph,Cm,S1,S2,S3,S4,Dl,Nf,J,P,Jt,Pt,Bph,Bhp,Bd,CMn)
                    if (Nf >= 4) call EMPM_states_free_X(S4,Jt,Pt)
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
        end do
        do nn = 1, Nf
            if (DoU) close(uS(nn))
            if (CMp) close(uP(nn))
        end do
        call timer_report(t0)
        write(*,'(/,a,/)') 'Full EMPM Hamiltonian has finished ...'

    contains

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
