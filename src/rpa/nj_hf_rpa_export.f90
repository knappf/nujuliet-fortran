! Export of HF-TDA & HF-RPA results ... counterpart of HF_RPA/HF_RPA_Export.jl ...
module nj_hf_rpa_export
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg, only : mtmm
    use nj_obdm, only : OBDM_export, OBDM_rN
    use nj_hf_legacy, only : HF_density_operator   ! LEGACY (removed in cleanup stage)
    implicit none
    private

    public :: HF_RPA_export, write_transition_file, em_moments, write_moments, plot_spectrum_file
    public :: spectrum_block_end, rpa_summary_header, dimensions_block, amplitudes_export

contains

    subroutine HF_RPA_export(Params,Orb,N_nu,E_corr,E_TDA,E_RPA,X_TDA,X_RPA,Y_RPA,rB_TDA,rB_RPA,C,Rho_RPA,Stability)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: N_nu(0:,:)
        real(dp), intent(in) :: E_corr
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        type(cmat_t), intent(in) :: X_RPA(0:,:), Y_RPA(0:,:)
        type(reduced_transition_t), intent(in) :: rB_TDA, rB_RPA
        type(o1b_t), intent(in) :: C, Rho_RPA
        logical, intent(in) :: Stability
        integer(i8) :: t0

        call timer_start(t0)
        call HF_RPA_summary(Params,N_nu,E_corr,Stability)
        call timer_report(t0)
        call timer_start(t0)
        call HF_RPA_binary_export(Params,N_nu,E_TDA,E_RPA,X_TDA,X_RPA,Y_RPA)
        call timer_report(t0)
        call timer_start(t0)
        call HF_RPA_spectrum_export(Params,N_nu,E_TDA,E_RPA)
        call timer_report(t0)
        call timer_start(t0)
        call HF_RPA_plot_spectrum_export(Params,N_nu,E_TDA,E_RPA)
        call timer_report(t0)
        call timer_start(t0)
        write(*,'(/,a)') 'Preparing export norms of the RPA amplitudes X & Y ...'
        call amplitudes_export(Params,N_nu,E_RPA,X_RPA,Y_RPA,'RPA')
        write(*,'(a)') TAB // 'Export of RPA norms of X & Y completed ...'
        call timer_report(t0)
        call timer_start(t0)
        call OBDM_export(Params,Orb,io_path(Params,'/RPA/RPA_Summary.dat'), &
                         io_path(Params,'/RPA/Densities/RPA_Radial_Densities.dat'),Rho_RPA,C)
        call timer_report(t0)
        call timer_start(t0)
        call HF_RPA_transitions_export(Params,Orb,N_nu,C,Rho_RPA,E_TDA,E_RPA,rB_TDA,rB_RPA)
        call timer_report(t0)
    end subroutine HF_RPA_export

    ! Header lines shared by the TDA/RPA/RRPA summary files ...
    subroutine rpa_summary_header(u,Params,Title,Ortho)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Title
        logical, intent(in) :: Ortho
        integer :: N_max
        N_max = Params%Calc%Nmax
        write(u,'(a)') Title
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
            write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (1-) is included ...'
        else
            write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (1-) is NOT included ...'
        end if
    end subroutine rpa_summary_header

    subroutine dimensions_block(u,N_nu,J_max,Kind)
        integer, intent(in) :: u, J_max
        integer, intent(in) :: N_nu(0:,:)
        character(len=*), intent(in) :: Kind
        integer :: J, N_M_tot
        write(u,'(/,a)') 'Information on dimensions of 1-phonon 1p-1h ' // Kind // ' subspaces ...'
        N_M_tot = 0
        do J = 0, J_max
            N_M_tot = N_M_tot + (N_nu(J,1) + N_nu(J,2)) * (2*J + 1)
        end do
        write(u,'(/,a)') 'Total number of 1 phonon states in J-scheme = ' // istr(sum(N_nu))
        write(u,'(a)') 'Total number of 1 phonon states in M-scheme = ' // istr(N_M_tot)
        write(u,'(/,a)') 'J' // TAB // 'N_phonon'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(N_nu(J,1) + N_nu(J,2))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = +'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(N_nu(J,1))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = -'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(N_nu(J,2))
        end do
    end subroutine dimensions_block

    subroutine HF_RPA_summary(Params,N_nu,E_corr,Stability)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        real(dp), intent(in) :: E_corr
        logical, intent(in) :: Stability
        integer :: u, J_max

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Preparing TDA & HF-RPA summary files ...'

        open(newunit=u, file=io_path(Params,'/RPA/TDA_Summary.dat'), status='replace', action='write')
        call rpa_summary_header(u,Params,'Spherical Tamm-Dancoff Approximation review:',Params%Calc%RPA%Ortho)
        call dimensions_block(u,N_nu,J_max,'TDA')
        close(u)

        open(newunit=u, file=io_path(Params,'/RPA/RPA_Summary.dat'), status='replace', action='write')
        call rpa_summary_header(u,Params,'Spherical Random-Phase Approximation review:',Params%Calc%RPA%Ortho)
        if (Stability) then
            write(u,'(a)') TAB // 'RPA system is STABLE ... no complex eigenvalues detected ...'
        else
            write(u,'(a)') TAB // 'RPA system is UNSTABLE ... complex eigenvalues detected ...'
        end if
        call dimensions_block(u,N_nu,J_max,'RPA')
        write(u,'(/,a)',advance='no') 'E_corr     = ' // cfmt_f(E_corr,12,6) // ' ' // TAB // ' MeV ... RPA ground-state correlation energy'
        write(u,'(/,a)') 'E_corr / A = ' // cfmt_f(E_corr / real(Params%Calc%A,dp),12,6) // ' ' // TAB // &
                         ' MeV ... RPA ground-state correlation energy per nuclon'
        close(u)
        write(*,'(a)') TAB // 'RPA & TDA summary files exported ...'
    end subroutine HF_RPA_summary

    ! Raw Float64 / ComplexF64 streams of energies & amplitudes, channel by channel ...
    subroutine HF_RPA_binary_export(Params,N_nu,E_TDA,E_RPA,X_TDA,X_RPA,Y_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        type(cmat_t), intent(in) :: X_RPA(0:,:), Y_RPA(0:,:)
        integer :: u1, u2, u3, u4, u5, J, P, J_max

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Performing binary export of TDA & RPA energies E & amplitudes X & Y ...'
        open(newunit=u1, file=io_path(Params,'/Bin/TDA_X.bin'), access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u2, file=io_path(Params,'/Bin/TDA_E.bin'), access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u3, file=io_path(Params,'/Bin/RPA_X.bin'), access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u4, file=io_path(Params,'/Bin/RPA_Y.bin'), access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u5, file=io_path(Params,'/Bin/RPA_E.bin'), access='stream', form='unformatted', status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (N_nu(J,P) == 0) cycle
                write(u1) X_TDA(J,P)%v
                write(u2) E_TDA(J,P)%v
                write(u3) X_RPA(J,P)%v
                write(u4) Y_RPA(J,P)%v
                write(u5) E_RPA(J,P)%v
            end do
        end do
        close(u1)
        close(u2)
        close(u3)
        close(u4)
        close(u5)
        write(*,'(a)') TAB // 'Binary export of TDA & RPA energies E amplitudes X & Y completed ...'
    end subroutine HF_RPA_binary_export

    ! Julia println(io, "\n") after each (J,P) block ... two line feeds ...
    subroutine spectrum_block_end(u)
        integer, intent(in) :: u
        write(u,'(a)') ''
        write(u,'(a)') ''
    end subroutine spectrum_block_end

    subroutine HF_RPA_spectrum_export(Params,N_nu,E_TDA,E_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        integer :: u, J, P, nu, J_max
        character(len=1) :: Pc

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Preparing export of TDA & RPA spectra ...'
        open(newunit=u, file=io_path(Params,'/RPA/Spectra/TDA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',20,.true.)
        do J = 0, J_max
            do P = 1, 2
                Pc = merge('+', '-', P == 1)
                do nu = 1, N_nu(J,P)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(Pc,5,.true.) // ' ' // cfmt_f(E_TDA(J,P)%v(nu),20,6,.true.)
                end do
                call spectrum_block_end(u)
            end do
        end do
        close(u)

        open(newunit=u, file=io_path(Params,'/RPA/Spectra/RPA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('Re E',20,.true.) // ' ' // &
                       cfmt_s('Im E',20,.true.)
        do J = 0, J_max
            do P = 1, 2
                Pc = merge('+', '-', P == 1)
                do nu = 1, N_nu(J,P)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(Pc,5,.true.) // ' ' // &
                                   cfmt_f(real(E_RPA(J,P)%v(nu)),20,6,.true.) // ' ' // cfmt_f(aimag(E_RPA(J,P)%v(nu)),20,6,.true.)
                end do
                call spectrum_block_end(u)
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of TDA & RPA spectra succesfully finished ...'
    end subroutine HF_RPA_spectrum_export

    ! Plot-ready spectrum ... ground state + all levels sorted by energy with shifted label positions ...
        ! As in the original only the first N of the N+1 sorted entries are written ...
    subroutine plot_spectrum_file(Path,N_nu,J_max,E)
        character(len=*), intent(in) :: Path
        integer, intent(in) :: N_nu(0:,:), J_max
        type(rvec_t), intent(in) :: E(0:,:)
        real(dp), allocatable :: Sol(:,:), E_full(:)
        integer, allocatable :: Sort(:)
        real(dp), parameter :: Delta = 1.0_dp
        integer :: N_ph, nu_count, J, P, nu, u
        character(len=1) :: Pc

        N_ph = sum(N_nu)
        allocate(Sol(4,N_ph+1), Sort(N_ph+1))
        Sol(:,1) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
        nu_count = 1
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, N_nu(J,P)
                    nu_count = nu_count + 1
                    Sol(:,nu_count) = [real(J,dp), real(P,dp), E(J,P)%v(nu), E(J,P)%v(nu)]
                end do
            end do
        end do
        E_full = Sol(3,:)
        call sortperm(E_full,Sort)
        Sol = Sol(:,Sort)
        do nu = 2, N_ph + 1
            if (Sol(4,nu) < Sol(4,nu-1)) then
                Sol(4,nu) = Sol(4,nu-1) + 1.0_dp
            else if (abs(Sol(4,nu) - Sol(4,nu-1)) < Delta) then
                Sol(4,nu) = Sol(4,nu-1) + 1.0_dp
            end if
        end do

        open(newunit=u, file=Path, status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // &
                       cfmt_s('label',12,.true.)
        do nu = 1, N_ph
            Pc = merge('+', '-', abs(Sol(2,nu) - 1.0_dp) < 1.0e-3_dp)
            write(u,'(a)') cfmt_d(round_even(Sol(1,nu)),5,.true.) // ' ' // cfmt_s(Pc,5,.true.) // ' ' // &
                           cfmt_f(Sol(3,nu),12,6,.true.) // ' ' // cfmt_f(Sol(4,nu),12,6,.true.)
        end do
        close(u)
    end subroutine plot_spectrum_file

    subroutine HF_RPA_plot_spectrum_export(Params,N_nu,E_TDA,E_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(rvec_t), allocatable :: E_re(:,:)
        integer :: J_max, J, P

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Preparing plot-ready export of RPA & TDA solutions ...'
        write(*,'(a)') TAB // 'Performing the export of plot-ready TDA spectra ...'
        call plot_spectrum_file(io_path(Params,'/RPA/Spectra/TDA_Plot.dat'),N_nu,J_max,E_TDA)
        write(*,'(a)') TAB // TAB // 'TDA plot-ready spectra successfully exported ...'
        allocate(E_re(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                E_re(J,P)%v = real(E_RPA(J,P)%v)
            end do
        end do
        write(*,'(a)') TAB // 'Performing the export of plot-ready RPA spectra ...'
        call plot_spectrum_file(io_path(Params,'/RPA/Spectra/RPA_Plot.dat'),N_nu,J_max,E_re)
        write(*,'(a)') TAB // TAB // 'RPA plot-ready spectra successfully exported ...'
    end subroutine HF_RPA_plot_spectrum_export

    subroutine amplitudes_export(Params,N_nu,E_RPA,X_RPA,Y_RPA,Kind)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: X_RPA(0:,:), Y_RPA(0:,:)
        character(len=*), intent(in) :: Kind
        real(dp), allocatable :: E_ord(:), x_ord(:), y_ord(:)
        integer, allocatable :: J_ord(:), P_ord(:), Sort(:)
        integer :: u, J, P, nu, J_max, N_tot, cnt

        J_max = 2*Params%Calc%Nmax + 1
        open(newunit=u, file=io_path(Params,'/' // Kind // '/Amplitudes/' // Kind // '_Amplitudes_Ordered.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('E',12,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('|Y|^2',12,.true.)
        do J = 0, J_max
            do P = 1, 2
                write(u,'(/,a)') 'J = ' // istr(J) // ',' // TAB // 'P = ' // merge('+', '-', P == 1)
                do nu = 1, N_nu(J,P)
                    write(u,'(a)') cfmt_f(real(E_RPA(J,P)%v(nu)),12,8,.true.) // ' ' // &
                                   cfmt_f(sum(abs(X_RPA(J,P)%v(:,nu))**2),12,8,.true.) // ' ' // &
                                   cfmt_f(sum(abs(Y_RPA(J,P)%v(:,nu))**2),12,8,.true.)
                end do
            end do
        end do
        close(u)

        N_tot = sum(N_nu)
        allocate(E_ord(N_tot), x_ord(N_tot), y_ord(N_tot), J_ord(N_tot), P_ord(N_tot), Sort(N_tot))
        cnt = 0
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, N_nu(J,P)
                    cnt = cnt + 1
                    E_ord(cnt) = real(E_RPA(J,P)%v(nu))
                    x_ord(cnt) = sum(abs(X_RPA(J,P)%v(:,nu))**2)
                    y_ord(cnt) = sum(abs(Y_RPA(J,P)%v(:,nu))**2)
                    J_ord(cnt) = J
                    P_ord(cnt) = P
                end do
            end do
        end do
        call sortperm(E_ord,Sort)
        open(newunit=u, file=io_path(Params,'/' // Kind // '/Amplitudes/' // Kind // '_Amplitudes.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('Ind',12,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // cfmt_s('J',12,.true.) // ' ' // &
                       cfmt_s('P',12,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('|Y|^2',12,.true.)
        do nu = 1, N_tot
            cnt = Sort(nu)
            write(u,'(a)') cfmt_d(nu,12,.true.) // ' ' // cfmt_f(E_ord(cnt),12,6,.true.) // ' ' // cfmt_d(J_ord(cnt),12,.true.) // &
                           ' ' // cfmt_s(merge('+', '-', P_ord(cnt) == 1),12,.true.) // ' ' // cfmt_f(x_ord(cnt),12,6,.true.) // &
                           ' ' // cfmt_f(y_ord(cnt),12,6,.true.)
        end do
        close(u)
    end subroutine amplitudes_export

    ! Energy-weighted moments m_-1 ... m_3 of a strength distribution ... m(1:5) = (m_-1, m_0, m_1, m_2, m_3) ...
        ! M0_Guarded = .true. includes m_0 only for |E| > 1e-3 (as done for the TDA E1 moments in the original) ...
    subroutine em_moments(E,B,M0_Guarded,m)
        real(dp), intent(in) :: E(:), B(:)
        logical, intent(in) :: M0_Guarded
        real(dp), intent(out) :: m(5)
        integer :: nu
        m = 0.0_dp
        do nu = 1, size(E)
            if (abs(E(nu)) > 1.0e-3_dp) then
                m(1) = m(1) + B(nu) / E(nu)
                if (M0_Guarded) m(2) = m(2) + B(nu)
            end if
            if (.not. M0_Guarded) m(2) = m(2) + B(nu)
            m(3) = m(3) + E(nu) * B(nu)
            m(4) = m(4) + E(nu)**2 * B(nu)
            m(5) = m(5) + E(nu)**3 * B(nu)
        end do
    end subroutine em_moments

    subroutine write_moments(u,Label,Unit,m)
        integer, intent(in) :: u
        character(len=*), intent(in) :: Label, Unit
        real(dp), intent(in) :: m(5)
        write(u,'(/,a)') TAB // Label
        write(u,'(a)') 'm _-1 = ' // cfmt_f(m(1),20,8) // TAB // 'e^2 ' // Unit // ' / MeV'
        write(u,'(a)') 'm _0  = ' // cfmt_f(m(2),20,8) // TAB // 'e^2 ' // Unit
        write(u,'(a)') 'm _1  = ' // cfmt_f(m(3),20,8) // TAB // 'e^2 ' // Unit // ' MeV'
        write(u,'(a)') 'm _2  = ' // cfmt_f(m(4),20,8) // TAB // 'e^2 ' // Unit // ' MeV^2'
        write(u,'(a)') 'm _3  = ' // cfmt_f(m(5),20,8) // TAB // 'e^2 ' // Unit // ' MeV^3'
    end subroutine write_moments

    ! "%-20s ..." header & "%-20.8f ..." rows of E, B_ph, B_is, B_iv ...
    subroutine write_transition_file(Path,E,T)
        character(len=*), intent(in) :: Path
        real(dp), intent(in) :: E(:)
        type(transition_t), intent(in) :: T
        integer :: u, nu
        open(newunit=u, file=Path, status='replace', action='write')
        write(u,'(a)') cfmt_s('E',20,.true.) // ' ' // cfmt_s('B_ph',20,.true.) // ' ' // cfmt_s('B_is',20,.true.) // ' ' // &
                       cfmt_s('B_iv',20,.true.)
        do nu = 1, size(E)
            write(u,'(a)') cfmt_f(E(nu),20,8,.true.) // ' ' // cfmt_f(T%ph(nu),20,8,.true.) // ' ' // &
                           cfmt_f(T%is(nu),20,8,.true.) // ' ' // cfmt_f(T%iv(nu),20,8,.true.)
        end do
        close(u)
    end subroutine write_transition_file

    subroutine HF_RPA_transitions_export(Params,Orb,N_nu,C,Rho_RPA,E_TDA,E_RPA,rB_TDA,rB_RPA)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: N_nu(0:,:)
        type(o1b_t), intent(in) :: C, Rho_RPA
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(reduced_transition_t), intent(in) :: rB_TDA, rB_RPA
        type(o1b_t) :: Rho_HF
        real(dp) :: pR2_TDA, pR4_TDA, pR2_RPA, pR4_RPA, Z, TRK_TDA(4), TRK_RPA(4)
        real(dp) :: mT(5,3,0:3), mR(5,3,0:3)
        real(dp), allocatable :: E0(:), E1(:), E2(:), E3(:)
        character(len=6), parameter :: Units(0:3) = ['fm^4  ', 'fm^2  ', 'fm^4  ', 'fm^6  ']
        character(len=11), parameter :: Kinds(3) = ['physical:  ', 'isoscalar: ', 'isovector: ']
        integer :: a_max, u, iso, X
        character(len=:), allocatable :: D

        Z = real(Params%Calc%Z,dp)
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2

        ! HF OBDM in the reference basis ...
        Rho_HF = HF_density_operator(a_max,C,Orb)
        Rho_HF%p = mtmm(C%p,Rho_HF%p)
        Rho_HF%n = mtmm(C%n,Rho_HF%n)

        pR2_TDA = OBDM_rN(Params,2,Rho_HF,C,Orb)
        pR4_TDA = OBDM_rN(Params,4,Rho_HF,C,Orb)
        pR2_RPA = OBDM_rN(Params,2,Rho_RPA,C,Orb)
        pR4_RPA = OBDM_rN(Params,4,Rho_RPA,C,Orb)

        write(*,'(/,a)') 'Preparing export of information on TDA & RPA 1-phonon electromagnetic transitions ...'

        TRK_TDA = trk(pR2_TDA,pR4_TDA)
        TRK_RPA = trk(pR2_RPA,pR4_RPA)

        ! Electromagnetic moments m_k for the physical, isoscalar & isovector strengths ...
        call strengths(rB_TDA%E0,E_TDA(0,1)%v,.false.,mT(:,:,0))
        call strengths(rB_TDA%E1,E_TDA(1,2)%v,.true.,mT(:,:,1))
        call strengths(rB_TDA%E2,E_TDA(2,1)%v,.false.,mT(:,:,2))
        call strengths(rB_TDA%E3,E_TDA(3,2)%v,.false.,mT(:,:,3))
        E0 = real(E_RPA(0,1)%v)
        E1 = real(E_RPA(1,2)%v)
        E2 = real(E_RPA(2,1)%v)
        E3 = real(E_RPA(3,2)%v)
        call strengths(rB_RPA%E0,E0,.false.,mR(:,:,0))
        call strengths(rB_RPA%E1,E1,.false.,mR(:,:,1))
        call strengths(rB_RPA%E2,E2,.false.,mR(:,:,2))
        call strengths(rB_RPA%E3,E3,.false.,mR(:,:,3))

        open(newunit=u, file=io_path(Params,'/RPA/TDA_Summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Review of TDA electromagnetic transition moments m_k ...'
        do iso = 1, 3
            do X = 0, 3
                call write_moments(u,'E' // istr(X) // ' ' // trim(Kinds(iso)),trim(Units(X)),mT(:,iso,X))
            end do
        end do
        write(u,'(/,a,/)') 'Review of TDA physical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values ...'
        do X = 0, 3
            write(u,'(a)') 'S_E' // istr(X) // ' = ' // cfmt_f(TRK_TDA(X+1),20,8) // TAB // 'e^2 ' // trim(Units(X)) // ' Mev'
        end do
        close(u)

        open(newunit=u, file=io_path(Params,'/RPA/RPA_Summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Review of RPA electromagnetic transition moments m_k ...'
        do iso = 1, 3
            do X = 0, 3
                call write_moments(u,'E' // istr(X) // ' ' // trim(Kinds(iso)),trim(Units(X)),mR(:,iso,X))
            end do
        end do
        write(u,'(/,a,/)') 'Review of RPA physical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values ...'
        do X = 0, 3
            write(u,'(a)') 'S_E' // istr(X) // ' = ' // cfmt_f(TRK_RPA(X+1),20,8) // TAB // 'e^2 ' // trim(Units(X))
        end do
        close(u)

        write(*,'(/,a)') 'Preparing export of RPA & TDA transition intensities ...'
        D = io_path(Params,'/RPA/Transitions')
        call write_transition_file(D // '/E0/TDA_E0.dat',E_TDA(0,1)%v,rB_TDA%E0)
        call write_transition_file(D // '/E1/TDA_E1.dat',E_TDA(1,2)%v,rB_TDA%E1)
        call write_transition_file(D // '/E1/TDA_E1_V.dat',E_TDA(1,2)%v,rB_TDA%E1_V)
        call write_transition_file(D // '/E1/TDA_E1_V_conv.dat',E_TDA(1,2)%v,rB_TDA%E1_VC)
        call write_transition_file(D // '/E1/TDA_E1_V_spin.dat',E_TDA(1,2)%v,rB_TDA%E1_VS)
        call write_transition_file(D // '/E1/TDA_E1_T.dat',E_TDA(1,2)%v,rB_TDA%E1_T)
        call write_transition_file(D // '/E1/TDA_E1_T_conv.dat',E_TDA(1,2)%v,rB_TDA%E1_TC)
        call write_transition_file(D // '/E1/TDA_E1_T_spin.dat',E_TDA(1,2)%v,rB_TDA%E1_TS)
        call write_transition_file(D // '/E1/TDA_E1_C.dat',E_TDA(1,2)%v,rB_TDA%E1_C)
        call write_transition_file(D // '/E1/TDA_E1_NLO_LWA.dat',E_TDA(1,2)%v,rB_TDA%E1_NLO_LWA)
        call write_transition_file(D // '/E2/TDA_E2.dat',E_TDA(2,1)%v,rB_TDA%E2)
        call write_transition_file(D // '/E3/TDA_E3.dat',E_TDA(3,2)%v,rB_TDA%E3)

        call write_transition_file(D // '/E0/RPA_E0.dat',E0,rB_RPA%E0)
        call write_transition_file(D // '/E1/RPA_E1.dat',E1,rB_RPA%E1)
        call write_transition_file(D // '/E1/RPA_E1_V.dat',E1,rB_RPA%E1_V)
        call write_transition_file(D // '/E1/RPA_E1_V_conv.dat',E1,rB_RPA%E1_VC)
        call write_transition_file(D // '/E1/RPA_E1_V_spin.dat',E1,rB_RPA%E1_VS)
        call write_transition_file(D // '/E1/RPA_E1_T.dat',E1,rB_RPA%E1_T)
        call write_transition_file(D // '/E1/RPA_E1_T_conv.dat',E1,rB_RPA%E1_TC)
        call write_transition_file(D // '/E1/RPA_E1_T_spin.dat',E1,rB_RPA%E1_TS)
        call write_transition_file(D // '/E1/RPA_E1_C.dat',E1,rB_RPA%E1_C)
        call write_transition_file(D // '/E1/RPA_E1_C_conv.dat',E1,rB_RPA%E1_CC)
        call write_transition_file(D // '/E1/RPA_E1_C_spin.dat',E1,rB_RPA%E1_CS)
        call write_transition_file(D // '/E1/RPA_E1_c.dat',E1,rB_RPA%E1_comp)
        call write_transition_file(D // '/E1/RPA_E1_NLO_LWA.dat',E1,rB_RPA%E1_NLO_LWA)
        call write_transition_file(D // '/E2/RPA_E2.dat',E2,rB_RPA%E2)
        call write_transition_file(D // '/E3/RPA_E3.dat',E3,rB_RPA%E3)

        write(*,'(a)') TAB // 'Export of information on TDA & RPA 1-phonon electromagnetic transitions completed ...'

    contains

        ! Thomas-Reiche-Kuhn sum rules for E0 ... E3 ...
        function trk(pR2,pR4) result(S)
            real(dp), intent(in) :: pR2, pR4
            real(dp) :: S(4)
            S(1) = HbarC**2 / m_p / (2.0_dp * pi) * Z * pR2
            S(2) = HbarC**2 / m_p * 9.0_dp / (8.0_dp * pi) * Z
            S(3) = HbarC**2 / m_p * 50.0_dp / (8.0_dp * pi) * Z * pR2
            S(4) = HbarC**2 / m_p * 147.0_dp / (8.0_dp * pi) * Z * pR4
        end function trk

        subroutine strengths(T,E,Guard,m)
            type(transition_t), intent(in) :: T
            real(dp), intent(in) :: E(:)
            logical, intent(in) :: Guard
            real(dp), intent(out) :: m(5,3)
            call em_moments(E,T%ph,Guard,m(:,1))
            call em_moments(E,T%is,Guard,m(:,2))
            call em_moments(E,T%iv,Guard,m(:,3))
        end subroutine strengths

    end subroutine HF_RPA_transitions_export

end module nj_hf_rpa_export
