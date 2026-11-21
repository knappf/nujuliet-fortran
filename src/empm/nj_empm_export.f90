! Export of the EMPM restricted-space TDA ... summary, binary energies & amplitudes, spectrum & phonon structure,
    ! counterparts of tda_pn_cp.log, 1phonon/1f_*.dat & phon_struct.dat of the original EMPM codes ...
module nj_empm_export
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_empm_space
    use nj_empm_fmat, only : fmat_t
    implicit none
    private

    public :: EMPM_TDA_export

contains

    subroutine EMPM_TDA_export(Params,Space,F_ph,N_nu,Orb_Phonon,Phonon,E_TDA,X_TDA)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        integer(i8) :: t0

        call timer_start(t0)
        call EMPM_TDA_summary(Params,Space,F_ph,N_nu,E_TDA)
        call EMPM_TDA_binary_export(Params,N_nu,E_TDA,X_TDA)
        call EMPM_TDA_spectrum_export(Params,N_nu,E_TDA)
        call EMPM_TDA_structure_export(Params,Space,N_nu,Orb_Phonon,Phonon,E_TDA,X_TDA)
        call timer_report(t0)
    end subroutine EMPM_TDA_export

    ! EMPM/TDA_summary.dat ... the text of the dev TDA_summary (nj_tda; header, Ortho, stability & the lowest excitation,
        ! dimensions; local copies so that the EMPM depends on core modules only) with the EMPM items (reference, windows,
        ! F-matrix) instead of the TDA parameters ...
    subroutine EMPM_TDA_summary(Params,Space,F_ph,N_nu,E_TDA)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        integer, intent(in) :: N_nu(0:,:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        integer :: u, J, P, nu, N_neg, J_min, P_min, J_max
        real(dp) :: E_min
        logical :: Ortho

        write(*,'(/,a)') 'Preparing the EMPM TDA summary file ...'
        Ortho = Params%Calc%EMPM%Ortho
        J_max = 2*Params%Calc%Nmax + 1
        N_neg = 0
        E_min = huge(1.0_dp)
        J_min = -1
        P_min = 0
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, N_nu(J,P)
                    if (Ortho .and. J == 1 .and. P == 2 .and. nu == 1 .and. N_nu(J,P) > 1) cycle
                    if (E_TDA(J,P)%v(nu) < -1.0e-8_dp) N_neg = N_neg + 1
                    if (E_TDA(J,P)%v(nu) < E_min) then
                        E_min = E_TDA(J,P)%v(nu)
                        J_min = J
                        P_min = P
                    end if
                end do
            end do
        end do

        open(newunit=u, file=empm_path(Params,'/TDA_summary.dat'), status='replace', action='write')
        call summary_header(u,Params,'Spherical Tamm-Dancoff Approximation in the restricted EMPM space review:')
        write(u,'(/,a)') 'EMPM parameters:'
        write(u,'(a)') TAB // "Reference = ''" // trim(Params%Calc%EMPM%Reference) // "'',   Ortho = " // trim(merge('true ','false',Ortho))
        write(u,'(/,a)') 'Restricted particle-hole space ... windows of the HF levels ordered by energy (1 = deepest level):'
        call window_line('proton holes     ',Space%pHole)
        call window_line('proton particles ',Space%pPart)
        call window_line('neutron holes    ',Space%nHole)
        call window_line('neutron particles',Space%nPart)
        write(u,'(/,a)') 'F-matrix:    L_max = ' // istr(F_ph%L_max) // ' , number of level pairs = ' // istr(sum(F_ph%N))
        if (Ortho) then
            write(u,'(/,a)') 'Orthogonalization of spurious center-of-mass levels (1-) is included ...'
        else
            write(u,'(/,a)') 'Orthogonalization of spurious center-of-mass levels (1-) is NOT included ...'
        end if
        if (N_neg == 0) then
            write(u,'(a)') TAB // 'TDA spectrum is STABLE ... no negative eigenvalues detected ...'
        else
            write(u,'(a)') TAB // 'TDA spectrum is UNSTABLE ... ' // istr(N_neg) // ' negative eigenvalue(s) detected ...'
        end if
        if (J_min >= 0) write(u,'(a)') TAB // 'Lowest TDA excitation energy:   E = ' // cfmt_f(E_min,12,6) // ' MeV   (J = ' // &
                                       istr(J_min) // ', P = ' // merge('+','-',P_min == 1) // ')'
        call dimensions_block(u,N_nu,J_max,'TDA')
        close(u)
        write(*,'(a)') TAB // 'EMPM TDA summary file exported ...'

    contains

        subroutine window_line(Name,Pos)
            character(len=*), intent(in) :: Name
            integer, intent(in) :: Pos(:)
            if (size(Pos) == 0) then
                write(u,'(a)') TAB // Name // '   none'
            else
                write(u,'(a)') TAB // Name // '   k = ' // cfmt_d(Space%Lev(Pos(1))%k,3) // ' .. ' // &
                               cfmt_d(Space%Lev(Pos(size(Pos)))%k,3) // '   E = ' // cfmt_f(Space%Lev(Pos(1))%E,10,4) // ' .. ' // &
                               cfmt_f(Space%Lev(Pos(size(Pos)))%E,10,4) // ' MeV'
            end if
        end subroutine window_line

    end subroutine EMPM_TDA_summary

    ! Raw Float64 streams of energies & amplitudes, channel by channel (J outer, P inner), empty channels skipped ...
    subroutine EMPM_TDA_binary_export(Params,N_nu,E_TDA,X_TDA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        integer :: u1, u2, J, P

        write(*,'(/,a)') 'Performing binary export of EMPM TDA energies E & amplitudes X ...'
        open(newunit=u1, file=empm_path(Params,'/bin/EMPM_TDA_X.bin'), access='stream', form='unformatted', status='replace', &
             action='write')
        open(newunit=u2, file=empm_path(Params,'/bin/EMPM_TDA_E.bin'), access='stream', form='unformatted', status='replace', &
             action='write')
        do J = 0, 2*Params%Calc%Nmax + 1
            do P = 1, 2
                if (N_nu(J,P) == 0) cycle
                write(u1) X_TDA(J,P)%v
                write(u2) E_TDA(J,P)%v
            end do
        end do
        close(u1)
        close(u2)
        write(*,'(a)') TAB // 'Binary export of EMPM TDA energies E & amplitudes X completed ...'
    end subroutine EMPM_TDA_binary_export

    subroutine EMPM_TDA_spectrum_export(Params,N_nu,E_TDA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        integer :: u, J, P, nu
        character(len=1) :: Pc

        open(newunit=u, file=empm_path(Params,'/spectra/TDA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',20,.true.)
        do J = 0, 2*Params%Calc%Nmax + 1
            do P = 1, 2
                Pc = merge('+', '-', P == 1)
                do nu = 1, N_nu(J,P)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(Pc,5,.true.) // ' ' // cfmt_f(E_TDA(J,P)%v(nu),20,6,.true.)
                end do
                write(u,'(a)') ''
                write(u,'(a)') ''
            end do
        end do
        close(u)
    end subroutine EMPM_TDA_spectrum_export

    ! Phonon structure ... running phonon index (J outer, P inner, nu), proton & neutron norms & the ph components
        ! with |X| >= 0.05 ordered by |X| ...
    subroutine EMPM_TDA_structure_export(Params,Space,N_nu,Orb_Phonon,Phonon,E_TDA,X_TDA)
        type(params_t), intent(in) :: Params
        type(empm_space_t), intent(in) :: Space
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(rvec_t), intent(in) :: E_TDA(0:,:)
        type(rmat_t), intent(in) :: X_TDA(0:,:)
        real(dp), parameter :: X_min = 0.05_dp
        integer, allocatable :: Perm(:)
        integer :: u, J, P, nu, k, ph, N_ph, Ind
        real(dp) :: Norm_p
        character(len=1) :: Pc, tc

        open(newunit=u, file=empm_path(Params,'/TDA_Structure.dat'), status='replace', action='write')
        write(u,'(a)') 'EMPM TDA phonon structure ... components with |X| >= ' // cfmt_f(X_min,0,2)
        Ind = 0
        do J = 0, 2*Params%Calc%Nmax + 1
            do P = 1, 2
                Pc = merge('+', '-', P == 1)
                N_ph = N_nu(J,P)
                if (N_ph == 0) cycle
                allocate(Perm(N_ph))
                do nu = 1, N_ph
                    Ind = Ind + 1
                    Norm_p = 0.0_dp
                    do k = 1, N_ph
                        if (Phonon(Orb_Phonon(J,P)%v(k))%tz == -1) Norm_p = Norm_p + X_TDA(J,P)%v(k,nu)**2
                    end do
                    write(u,'(/,a)') '____________________________________________________________'
                    write(u,'(a)') 'i = ' // istr(Ind) // TAB // 'J = ' // istr(J) // Pc // TAB // 'nu = ' // istr(nu) // TAB // &
                                   'E = ' // cfmt_f(E_TDA(J,P)%v(nu),0,5) // ' MeV' // TAB // 'proton norm = ' // &
                                   cfmt_f(Norm_p,0,5) // TAB // 'neutron norm = ' // cfmt_f(1.0_dp - Norm_p,0,5)
                    call sortperm(-abs(X_TDA(J,P)%v(:,nu)),Perm)
                    do k = 1, N_ph
                        if (abs(X_TDA(J,P)%v(Perm(k),nu)) < X_min) exit
                        ph = Orb_Phonon(J,P)%v(Perm(k))
                        tc = merge('p', 'n', Phonon(ph)%tz == -1)
                        write(u,'(a)') '   ' // tc // '   ' // cfmt_s(level_label(Space%Lev(Phonon(ph)%p)),8,.true.) // &
                                       cfmt_s(level_label(Space%Lev(Phonon(ph)%h)),8,.true.) // '^-1   ' // &
                                       cfmt_f(X_TDA(J,P)%v(Perm(k),nu),10,5)
                    end do
                end do
                deallocate(Perm)
            end do
        end do
        close(u)
    end subroutine EMPM_TDA_structure_export

    ! Header of the summary ... copy of ph_summary_header (dev TDA/RPA summaries) ...
    subroutine summary_header(u,Params,Title)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Title
        integer :: N_max
        N_max = Params%Calc%Nmax
        write(u,'(a)') Title
        write(u,'(/,a)') 'Nuclide data:    A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)
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
    end subroutine summary_header

    ! Dimension tables of the 1p-1h subspaces ... copy of ph_dimensions_block on the EMPM channel counts N_nu(J,P) ...
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

end module nj_empm_export
