! Spherical quasiparticle Tamm-Dancoff approximation ... counterpart of QTDA/*.jl ...
module nj_qtda
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make, orbitals_2qp_make
    use nj_o1b
    use nj_obdm, only : OBDM_chR2
    use nj_tr1b, only : qptr1b_t, qpTr1b_initialize
    use nj_qpo2b, only : qpO2b_import
    use nj_hfb, only : HFB_density_operator
    use nj_qex
    !$ use omp_lib
    implicit none
    private

    public :: QTDA, qex_header

contains

    subroutine QTDA(Params)
        type(params_t), intent(in) :: Params
        integer :: N_max, N_threads
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        N_max = Params%Calc%Nmax
        write(*,'(a)') 'Starting spherical QTDA calculation with residual NO2B NN+NNN interaction based on available ' // &
                       'reference quasiparticle mean-field calculation'
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
            write(*,'(/,a)') 'Error! ... No precomputed reference quasiparticle mean-field (HFB, BCS) solution is available in ' // &
                             'given Input_File path ... run HFB/BCS solver first ...'
            return
        end if
        if (is_dir(io_path(Params,'/QTDA'))) call remove_dir(io_path(Params,'/QTDA'))
        call make_dir(io_path(Params,'/QTDA/Densities'))
        call make_dir(io_path(Params,'/QTDA/Transitions/E0'))
        call make_dir(io_path(Params,'/QTDA/Transitions/E1'))
        call make_dir(io_path(Params,'/QTDA/Transitions/E2'))
        call make_dir(io_path(Params,'/QTDA/Transitions/E3'))
        call make_dir(io_path(Params,'/QTDA/Spectra'))

        write(*,'(/,a)') 'Starting QTDA calculation ...'
        call timer_start(t0)
        call QTDA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'Spherical QTDA solver run completed ...'
    end subroutine QTDA

    subroutine QTDA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(orb2b_2qp_t) :: Orb_2qp
        type(o1b_t) :: C, U, V, Rho, Kappa
        type(qpo1b_t) :: H_N, qpN
        type(qpo2b_t) :: H_NN
        type(orb2b_t) :: Orb_NN
        type(qptr1b_t) :: qpTrOp
        type(rmat_t), allocatable :: A(:,:)
        type(rvec_t), allocatable :: E_QTDA(:,:)
        type(cmat_t), allocatable :: X(:,:), Y(:,:)
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        real(dp) :: chR2
        integer :: J, P

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        Orb_2qp = orbitals_2qp_make(Params,Orb)
        C = O1b_import(Params,Orb,io_path(Params,'/Bin/C.bin'))
        U = O1b_import(Params,Orb,io_path(Params,'/Bin/U.bin'))
        V = O1b_import(Params,Orb,io_path(Params,'/Bin/V.bin'))
        call HFB_density_operator(Params,U,V,Orb,Rho,Kappa)
        H_N = qpO1B_import(Params,Orb,io_path(Params,'/Bin/qpH1B.bin'))
        call qpO2b_import(Params,Orb,io_path(Params,'/Bin/qpH2B.bin'),H_NN,Orb_NN)
        chR2 = OBDM_chR2(Params,Orb,C,Rho)
        qpTrOp = qpTr1b_initialize(Params,Orb,chR2,C,U,V)
        qpN = qpN1B_make(Params,U,V)

        write(*,'(/,a)') 'Allocating the QTDA matrix A ...'
        call qex_allocate(Params,Orb,Orb_NN,Orb_2qp,H_N,H_NN,A)
        write(*,'(/,a)') 'The QTDA matrix A succesfully allocated ...'

        call QTDA_diagonalize(Params,Orb,Orb_2qp,A,qpN,qpTrOp,E_QTDA,X)

        ! Amplitudes as complex X & zero Y for the shared transition routines ...
        allocate(Y(2,0:ubound(X,2)))
        do J = 0, ubound(X,2)
            do P = 1, 2
                allocate(Y(P,J)%v(size(X(P,J)%v,1),size(X(P,J)%v,2)))
                Y(P,J)%v = (0.0_dp,0.0_dp)
            end do
        end do
        write(*,'(/,a)') 'Calculating the QTDA 1-phonon reduced transition matrix elements rM^lambda ...'
        rM = qex_rM(Orb_2qp,X,Y,qpTrOp)
        write(*,'(a)') TAB // 'QTDA 1-phonon reduced transition matrix elements rM^lambda succesfully calculated ...'
        write(*,'(/,a)') 'Calculating the QTDA 1-phonon reduced transition intensities rB^lambda ...'
        rB = qex_rB(Params,Orb_2qp,rM,E_QTDA)
        write(*,'(a)') TAB // 'QTDA reduced transition intensities rB evaluated ...'

        call QTDA_export(Params,Orb,Orb_2qp,Rho,E_QTDA,X,rB)
        write(*,'(a,/)') TAB // 'QTDA solver executed properly ...'
    end subroutine QTDA_solver

    subroutine QTDA_diagonalize(Params,Orb,Orb_2qp,Amat,qpN,qpTrOp,E,X)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(rmat_t), intent(in) :: Amat(:,0:)
        type(qpo1b_t), intent(in) :: qpN
        type(qptr1b_t), intent(in) :: qpTrOp
        type(rvec_t), allocatable, intent(out) :: E(:,:)
        type(cmat_t), allocatable, intent(out) :: X(:,:)
        real(dp), allocatable :: Spur(:), Vb(:,:), Us(:,:), A_ort(:,:), Z(:,:)
        complex(dp), allocatable :: Wc(:), Vc(:,:)
        logical :: Is_Real
        integer :: J_max, J, P, N_qp, Spur_ind, a

        J_max = 2*Params%Calc%Nmax + 1
        allocate(E(2,0:J_max), X(2,0:J_max))
        write(*,'(/,a)') 'Diagonalizing the QTDA matrix Amat ...'
        do J = 0, J_max
            do P = 1, 2
                N_qp = Orb_2qp%N(P,J)
                allocate(E(P,J)%v(N_qp), X(P,J)%v(N_qp,N_qp))
                if (N_qp == 0) cycle
                if (.not. Params%Calc%QTDA%Ortho .or. .not. ((J == 0 .and. P == 1) .or. (J == 1 .and. P == 2))) then
                    ! Amat is exactly symmetric ... eigen! dispatches to the symmetric solver ...
                    allocate(Z(N_qp,N_qp))
                    call eig_sym(Amat(P,J)%v,E(P,J)%v,Z)
                    X(P,J)%v = cmplx(Z, 0.0_dp, kind=dp)
                    deallocate(Z)
                else
                    if (J == 0) then
                        Spur = spurious_PN_initialize(Orb,Orb_2qp,qpN)
                    else
                        Spur = spurious_CM_initialize(Orb_2qp,qpTrOp)
                    end if
                    allocate(Vb(N_qp,N_qp), Us(N_qp,N_qp-1))
                    call spurious_orthogonalize(N_qp,Spur,Vb,Spur_ind)
                    do a = 1, Spur_ind - 1
                        Us(:,a) = Vb(:,a)
                    end do
                    do a = Spur_ind + 1, N_qp
                        Us(:,a-1) = Vb(:,a)
                    end do
                    A_ort = mtm(Us,mm(Amat(P,J)%v,Us))
                    ! Plain Matrix in the original ... symmetric solver only if exactly symmetric, otherwise geevx ...
                    allocate(Wc(N_qp-1), Vc(N_qp-1,N_qp-1))
                    call eig_general(A_ort,Wc,Vc,Is_Real)
                    if (.not. Is_Real) call fatal('Complex QTDA eigenvalues in the orthogonalized subspace ...')
                    E(P,J)%v(1) = 0.0_dp
                    E(P,J)%v(2:N_qp) = real(Wc)
                    X(P,J)%v(:,1) = cmplx(Spur, 0.0_dp, kind=dp)
                    X(P,J)%v(:,2:N_qp) = cmplx(mm(Us,real(Vc)), 0.0_dp, kind=dp)
                    deallocate(Vb, Us, Wc, Vc)
                end if
            end do
        end do
        write(*,'(/,a)') 'The QTDA matrix Amat succesfully diagonalized ...'
    end subroutine QTDA_diagonalize

    subroutine QTDA_export(Params,Orb,Orb_2qp,Rho,E,X,rB)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t), intent(in) :: Orb_2qp
        type(o1b_t), intent(in) :: Rho
        type(rvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: X(:,0:)
        type(reduced_transition_t), intent(in) :: rB
        integer :: N_max, J_max, u, J, P, nu

        N_max = Params%Calc%Nmax
        J_max = 2*N_max + 1

        ! Summary ...
        write(*,'(/,a)') 'Preparing QTDA summary file ...'
        open(newunit=u, file=io_path(Params,'/QTDA/QTDA_Summary.dat'), status='replace', action='write')
        write(u,'(a)') 'Spherical Quasiparticle Tamm-Dancoff Approximation review:'
        call qex_header(u,Params,Params%Calc%QTDA%Ortho)
        call qex_dimensions_summary(u,Orb_2qp,J_max,'QTDA 1-phonon subspaces ...')
        close(u)
        write(*,'(/,a)') 'QTDA summary exported ...'

        ! Binary export ... X column by column, then E ...
        write(*,'(/,a)') 'Performing binary export of QTDA energies E & amplitudes X ...'
        open(newunit=u, file=io_path(Params,'/Bin/QTDA_X.bin'), access='stream', form='unformatted', status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(X(P,J)%v) > 0) write(u) real(X(P,J)%v)
            end do
        end do
        close(u)
        open(newunit=u, file=io_path(Params,'/Bin/QTDA_E.bin'), access='stream', form='unformatted', status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(E(P,J)%v) > 0) write(u) E(P,J)%v
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Binary export of QTDA energies E amplitudes X completed ...'

        ! Spectrum ...
        write(*,'(/,a)') 'Preparing export of QTDA spectra ...'
        open(newunit=u, file=io_path(Params,'/QTDA/Spectra/QTDA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',10,.true.)
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, Orb_2qp%N(P,J)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(merge('+','-',P == 1),5,.true.) // ' ' // &
                                   cfmt_f(E(P,J)%v(nu),10,6,.true.)
                end do
                write(u,'(/)')
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of QTDA spectra succesfully finished ...'

        call qex_spectrum_plot_export(Params,Orb_2qp,E,'QTDA')
        call qex_transitions_export(Params,Orb,Orb_2qp,Rho,E,rB,'QTDA')
    end subroutine QTDA_export

    ! Common header of the QTDA & QRPA summaries ...
    subroutine qex_header(u,Params,Ortho)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        logical, intent(in) :: Ortho
        integer :: N_max
        N_max = Params%Calc%Nmax
        write(u,'(/,a)') 'Nuclid data:    A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)
        write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                         ', s.p. J-basis size = ' // istr((N_max+1)*(N_max+2)/2) // ', s.p. M-basis size = ' // &
                         istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Combined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Only pure 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else
            write(u,'(/,a)') TAB // TAB // TAB // 'No Center-of-Mass System (CMS) motion correction is set ...'
        end if
        if (Ortho) then
            write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (0+, 1-) is included ...'
        else
            write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (0+, 1-) is NOT included ...'
        end if
    end subroutine qex_header

end module nj_qtda
