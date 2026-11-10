! Export of HF-RRPA results ... counterpart of HF_RRPA/HF_RRPA_Export.jl ...
module nj_hf_rrpa_export
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_o1b, only : O1b_export
    use nj_obdm, only : OBDM_export, OBDM_rN
    use nj_hf_rpa_export, only : rpa_summary_header, dimensions_block, amplitudes_export, spectrum_block_end, &
                                 plot_spectrum_file, write_transition_file, em_moments, write_moments
    implicit none
    private

    public :: HF_RRPA_export

contains

    subroutine HF_RRPA_export(Params,Orb,N_nu,E_corr,C,Rho,h,E_RPA,X_RPA,Y_RPA,rB)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: N_nu(0:,:)
        real(dp), intent(in) :: E_corr
        type(o1b_t), intent(in) :: C, Rho, h
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: X_RPA(0:,:), Y_RPA(0:,:)
        type(reduced_transition_t), intent(in) :: rB
        type(o1b_t) :: Rho_LHO
        type(rvec_t), allocatable :: E_re(:,:)
        integer :: J_max, J, P, u
        integer(i8) :: t0

        J_max = 2*Params%Calc%Nmax + 1
        call HF_RRPA_export_summary(Params,N_nu,E_corr)
        call O1b_export(Params,Orb,C,io_path(Params,'/Bin/C_RRPA.bin'))
        call O1b_export(Params,Orb,Rho,io_path(Params,'/Bin/Rho_RRPA.bin'))
        call O1b_export(Params,Orb,h,io_path(Params,'/Bin/h_RRPA.bin'))

        call timer_start(t0)
        write(*,'(/,a)') 'Performing binary export of RRPA energies E & amplitudes X & Y ...'
        call write_stream(io_path(Params,'/Bin/RRPA_X.bin'),1)
        call write_stream(io_path(Params,'/Bin/RRPA_Y.bin'),2)
        call write_stream(io_path(Params,'/Bin/RRPA_E.bin'),3)
        write(*,'(a)') TAB // 'Binary export of RRPA energies E amplitudes X & Y completed ...'
        call timer_report(t0)

        write(*,'(/,a)') 'Preparing export norms of the RRPA amplitudes X & Y ...'
        call amplitudes_export(Params,N_nu,E_RPA,X_RPA,Y_RPA,'RRPA')
        write(*,'(a)') TAB // 'Export of RRPA norms of X & Y completed ...'

        ! Spectrum ...
        write(*,'(/,a)') 'Preparing export of RRPA spectra ...'
        open(newunit=u, file=io_path(Params,'/RRPA/Spectra/RRPA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('Re E',20,.true.) // ' ' // &
                       cfmt_s('Im E',20,.true.)
        do J = 0, J_max
            do P = 1, 2
                call spectrum_rows(J,P)
                call spectrum_block_end(u)
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of RRPA spectra succesfully finished ...'

        write(*,'(/,a)') 'Preparing plot-ready export of RRPA solutions ...'
        allocate(E_re(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                E_re(J,P)%v = real(E_RPA(J,P)%v)
            end do
        end do
        write(*,'(a)') TAB // 'Performing the export of plot-ready RPA spectra ...'
        call plot_spectrum_file(io_path(Params,'/RRPA/Spectra/RRPA_Plot.dat'),N_nu,J_max,E_re)
        write(*,'(a)') TAB // TAB // 'RRPA plot-ready spectra successfully exported ...'

        call HF_RRPA_export_occupation(Params,Orb,Rho)

        call timer_start(t0)
        Rho_LHO%p = mmt(mm(C%p,Rho%p),C%p)
        Rho_LHO%n = mmt(mm(C%n,Rho%n),C%n)
        call OBDM_export(Params,Orb,io_path(Params,'/RRPA/RRPA_Summary.dat'), &
                         io_path(Params,'/RRPA/Densities/RRPA_Radial_Densities.dat'),Rho_LHO,C)
        call timer_report(t0)

        call timer_start(t0)
        call HF_RRPA_export_single_particle_spectrum(Params,Orb,h)
        call timer_report(t0)

        call HF_RRPA_export_transitions(Params,Orb,N_nu,C,Rho,E_re,rB)

    contains

        subroutine write_stream(Path,What)
            character(len=*), intent(in) :: Path
            integer, intent(in) :: What
            integer :: us, JJ, PP
            open(newunit=us, file=Path, access='stream', form='unformatted', status='replace', action='write')
            do JJ = 0, J_max
                do PP = 1, 2
                    if (N_nu(JJ,PP) == 0) cycle
                    select case (What)
                    case (1)
                        write(us) X_RPA(JJ,PP)%v
                    case (2)
                        write(us) Y_RPA(JJ,PP)%v
                    case default
                        write(us) E_RPA(JJ,PP)%v
                    end select
                end do
            end do
            close(us)
        end subroutine write_stream

        subroutine spectrum_rows(JJ,PP)
            integer, intent(in) :: JJ, PP
            integer :: nu
            do nu = 1, N_nu(JJ,PP)
                write(u,'(a)') cfmt_d(JJ,5,.true.) // ' ' // cfmt_s(merge('+','-',PP == 1),5,.true.) // ' ' // &
                               cfmt_f(real(E_RPA(JJ,PP)%v(nu)),20,6,.true.) // ' ' // cfmt_f(aimag(E_RPA(JJ,PP)%v(nu)),20,6,.true.)
            end do
        end subroutine spectrum_rows

    end subroutine HF_RRPA_export

    subroutine HF_RRPA_export_summary(Params,N_nu,E_corr)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        real(dp), intent(in) :: E_corr
        integer :: u
        write(*,'(/,a)') 'Preparing RRPA summary ...'
        open(newunit=u, file=io_path(Params,'/RRPA/RRPA_Summary.dat'), status='replace', action='write')
        call rpa_summary_header(u,Params,'Spherical Renormalized-Random-Phase Approximation review:',Params%Calc%RRPA%Ortho)
        call dimensions_block(u,N_nu,2*Params%Calc%Nmax + 1,'RRPA')
        write(u,'(/,a)',advance='no') 'E_corr     = ' // cfmt_f(E_corr,12,6) // ' ' // TAB // ' MeV ... RPA ground-state correlation energy'
        write(u,'(/,a)') 'E_corr / A = ' // cfmt_f(E_corr / real(Params%Calc%A,dp),12,6) // ' ' // TAB // &
                         ' MeV ... RPA ground-state correlation energy per nuclon'
        write(u,'(/,a)') TAB // 'To get the total correct ground-state energy add the mean-field energy ...'
        close(u)
        write(*,'(/,a)') 'RRPA summary exported ...'
    end subroutine HF_RRPA_export_summary

    subroutine HF_RRPA_export_occupation(Params,Orb,Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: Rho
        real(dp), allocatable :: pO_HF(:), nO_HF(:), pO(:), nO(:)
        integer, allocatable :: pS(:), nS(:)
        integer :: a_max, a, u

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(pO_HF(a_max), nO_HF(a_max), pO(a_max), nO(a_max), pS(a_max), nS(a_max))
        do a = 1, a_max
            pO_HF(a) = real(Orb(a)%pO,dp)
            nO_HF(a) = real(Orb(a)%nO,dp)
            pO(a) = round_digits(Rho%p(a,a),4)
            nO(a) = round_digits(Rho%n(a,a),4)
        end do
        call sortperm_rev(pO,pS)
        call sortperm_rev(nO,nS)
        open(newunit=u, file=io_path(Params,'/RRPA/Densities/RRPA_Occupation.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('a',6,.true.) // ' ' // cfmt_s('a_p',6,.true.) // ' ' // cfmt_s('pN_HF',12,.true.) // ' ' // &
                       cfmt_s('pN_RRPA',12,.true.) // ' ' // cfmt_s('pN_dep',12,.true.) // ' ' // cfmt_s('a_n',6,.true.) // ' ' // &
                       cfmt_s('nN_HF',12,.true.) // ' ' // cfmt_s('nN_RRPA',12,.true.) // ' ' // cfmt_s('nN_dep',12,.true.)
        do a = 1, a_max
            write(u,'(a)') cfmt_d(a,6,.true.) // ' ' // cfmt_d(pS(a),6,.true.) // ' ' // cfmt_f(pO_HF(pS(a)),12,6,.true.) // ' ' // &
                           cfmt_f(pO(pS(a)),12,6,.true.) // ' ' // cfmt_f(pO(pS(a)) - pO_HF(pS(a)),12,6,.true.) // ' ' // &
                           cfmt_d(nS(a),6,.true.) // ' ' // cfmt_f(nO_HF(nS(a)),12,6,.true.) // ' ' // &
                           cfmt_f(nO(nS(a)),12,6,.true.) // ' ' // cfmt_f(nO(nS(a)) - nO_HF(nS(a)),12,6,.true.)
        end do
        close(u)
    end subroutine HF_RRPA_export_occupation

    ! Effective RRPA single-particle energies with (l, 2j) assigned from the eigenvector content ...
    subroutine HF_RRPA_export_single_particle_spectrum(Params,Orb,h)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: h
        real(dp), allocatable :: E(:), Cm(:,:)
        integer, allocatable :: nr(:,:), l_q(:), j_q(:)
        integer :: N_max, a_max, a, b, u, iso
        real(dp) :: lS, jS

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(E(a_max), Cm(a_max,a_max), nr(0:N_max,1:2*N_max+1), l_q(a_max), j_q(a_max))
        write(*,'(/,a)') 'Preparing export of effective RRPA single-particle energies ...'
        open(newunit=u, file=io_path(Params,'/RRPA/RRPA_Summary.dat'), position='append', action='write')
        do iso = 1, 2
            if (iso == 1) then
                call eig_sym(h%p,E,Cm)
            else
                call eig_sym(h%n,E,Cm)
            end if
            do a = 1, a_max
                lS = 0.0_dp
                jS = 0.0_dp
                do b = 1, a_max
                    lS = lS + Orb(b)%l * Cm(b,a)**2
                    jS = jS + Orb(b)%j * Cm(b,a)**2
                end do
                l_q(a) = round_even(lS)
                j_q(a) = round_even(jS)
            end do
            nr = 0
            if (iso == 1) then
                write(u,'(/,a)') 'Proton RRPA single-particle states ...'
            else
                write(u,'(/,a)') 'Neutron RRPA single-particle states ...'
            end if
            write(u,'(a)') '____________________________________________________________________'
            write(u,'(a)') cfmt_s('n',4) // ' ' // cfmt_s('l',4) // ' ' // cfmt_s('2j',4) // ' ' // cfmt_s('E',14)
            do a = 1, a_max
                write(u,'(a)') cfmt_d(nr(l_q(a),j_q(a)),4) // ' ' // cfmt_d(l_q(a),4) // ' ' // cfmt_d(j_q(a),4) // ' ' // &
                               cfmt_f(E(a),14,6)
                nr(l_q(a),j_q(a)) = nr(l_q(a),j_q(a)) + 1
            end do
            write(u,'(a)') '____________________________________________________________________'
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of RRPA effective single-particle energies finished ...'
    end subroutine HF_RRPA_export_single_particle_spectrum

    subroutine HF_RRPA_export_transitions(Params,Orb,N_nu,C,Rho,E,rB)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: N_nu(0:,:)
        type(o1b_t), intent(in) :: C, Rho
        type(rvec_t), intent(in) :: E(0:,:)
        type(reduced_transition_t), intent(in) :: rB
        character(len=6), parameter :: Units(0:3) = ['fm^4  ', 'fm^2  ', 'fm^4  ', 'fm^6  ']
        character(len=11), parameter :: Kinds(3) = ['physical:  ', 'isoscalar: ', 'isovector: ']
        real(dp) :: pR2, pR4, Z, TRK(4), m(5,3,0:3)
        integer :: u, iso, X
        character(len=:), allocatable :: D

        Z = real(Params%Calc%Z,dp)
        pR2 = OBDM_rN(Params,2,Rho,C,Orb)
        pR4 = OBDM_rN(Params,4,Rho,C,Orb)
        write(*,'(/,a)') 'Preparing export of information on RRPA 1-phonon electromagnetic transitions ...'
        TRK(1) = HbarC**2 / m_p / (2.0_dp * pi) * Z * pR2
        TRK(2) = HbarC**2 / m_p * 9.0_dp / (8.0_dp * pi) * Z
        TRK(3) = HbarC**2 / m_p * 50.0_dp / (8.0_dp * pi) * Z * pR2
        TRK(4) = HbarC**2 / m_p * 147.0_dp / (8.0_dp * pi) * Z * pR4

        ! m_0 of E0 is accumulated only for |E| > 1e-3, as in the original ...
        call strengths(rB%E0,E(0,1)%v,.true.,m(:,:,0))
        call strengths(rB%E1,E(1,2)%v,.false.,m(:,:,1))
        call strengths(rB%E2,E(2,1)%v,.false.,m(:,:,2))
        call strengths(rB%E3,E(3,2)%v,.false.,m(:,:,3))

        open(newunit=u, file=io_path(Params,'/RRPA/RRPA_Summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Review of RRPA electrogmagnetic transition moments m_k ...'
        do iso = 1, 3
            do X = 0, 3
                call write_moments(u,'E' // istr(X) // ' ' // trim(Kinds(iso)),trim(Units(X)),m(:,iso,X))
            end do
        end do
        write(u,'(/,a,/)') 'Review of RRPA physical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values ...'
        do X = 0, 3
            write(u,'(a)') 'S_E' // istr(X) // ' = ' // cfmt_f(TRK(X+1),20,8) // TAB // 'e^2 ' // trim(Units(X))
        end do
        close(u)

        write(*,'(/,a)') 'Preparing export of RRPA transition intensities ...'
        D = io_path(Params,'/RRPA/Transitions')
        call write_transition_file(D // '/E0/RRPA_E0.dat',E(0,1)%v,rB%E0)
        call write_transition_file(D // '/E1/RRPA_E1.dat',E(1,2)%v,rB%E1)
        call write_transition_file(D // '/E1/RRPA_E1_V.dat',E(1,2)%v,rB%E1_V)
        call write_transition_file(D // '/E1/RRPA_E1_V_conv.dat',E(1,2)%v,rB%E1_VC)
        call write_transition_file(D // '/E1/RRPA_E1_V_spin.dat',E(1,2)%v,rB%E1_VS)
        call write_transition_file(D // '/E1/RRPA_E1_T.dat',E(1,2)%v,rB%E1_T)
        call write_transition_file(D // '/E1/RRPA_E1_T_conv.dat',E(1,2)%v,rB%E1_TC)
        call write_transition_file(D // '/E1/RRPA_E1_T_spin.dat',E(1,2)%v,rB%E1_TS)
        call write_transition_file(D // '/E1/RRPA_E1_C.dat',E(1,2)%v,rB%E1_C)
        call write_transition_file(D // '/E1/RRPA_E1_C_conv.dat',E(1,2)%v,rB%E1_CC)
        call write_transition_file(D // '/E1/RRPA_E1_C_spin.dat',E(1,2)%v,rB%E1_CS)
        call write_transition_file(D // '/E1/RRPA_E1_c.dat',E(1,2)%v,rB%E1_comp)
        call write_transition_file(D // '/E1/RRPA_E1_NLO_LWA.dat',E(1,2)%v,rB%E1_NLO_LWA)
        call write_transition_file(D // '/E2/RRPA_E2.dat',E(2,1)%v,rB%E2)
        call write_transition_file(D // '/E3/RRPA_E3.dat',E(3,2)%v,rB%E3)
        write(*,'(a)') TAB // 'Export of information on RRPA 1-phonon electromagnetic transitions completed ...'

    contains

        subroutine strengths(T,Ev,Guard,mm5)
            type(transition_t), intent(in) :: T
            real(dp), intent(in) :: Ev(:)
            logical, intent(in) :: Guard
            real(dp), intent(out) :: mm5(5,3)
            call em_moments(Ev,T%ph,Guard,mm5(:,1))
            call em_moments(Ev,T%is,Guard,mm5(:,2))
            call em_moments(Ev,T%iv,Guard,mm5(:,3))
        end subroutine strengths

    end subroutine HF_RRPA_export_transitions

end module nj_hf_rrpa_export
