! Export of the RPA results (NuJuliet dev) ... counterpart of src/RPA/RPA_export.jl ... outputs in io/<Path>/<Name>/:
    ! RPA_summary.dat, bin/{RPA_X,RPA_Y,RPA_E,orb,C,rho[,H1b][,sigma]}.bin, spectra/{RPA,RPA_plot}.dat,
    ! amplitudes/{RPA_amplitudes_ordered,RPA_amplitudes}.dat, densities/{<Name>_radial_densities,RPA_occupation}.dat,
    ! transitions/E*/RPA_*.dat ...
module nj_rpa_export
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg, only : mm
    use nj_orbitals, only : orbitals_import, orbitals_export, orbitals_update, OBDM_natural_orbitals
    use nj_o1b, only : O1b_export, O1b_transformation, O1b_diagonal
    use nj_o2b, only : O2b_export
    use nj_obdm, only : rho_export
    use nj_ph_common
    use nj_rpa, only : rpa_energy_t
    implicit none
    private

    public :: RPA_export

contains

    ! Rho_LHO = RPA OBDM in the LHO basis, h = reference 1-body Hamiltonian (reference basis), Sigma = RPA Sigma (TBCF) ...
    subroutine RPA_export(Params,Orb,Orb_1p1h,Orb_NN,Stability,E_corr,E,X,Y,rB,C,Rho_LHO,Sigma,h)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        logical, intent(in) :: Stability
        type(rpa_energy_t), intent(in) :: E_corr
        type(cvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        type(reduced_transition_t), intent(in) :: rB
        type(o1b_t), intent(in) :: C, Rho_LHO, h
        type(o2b_t), intent(in) :: Sigma
        type(rvec_t), allocatable :: E_re(:,:)
        character(len=:), allocatable :: Out_Dir
        integer :: J, P
        integer(i8) :: t0

        Out_Dir = io_dir(Params,Params%Calc%RPA%Name)
        allocate(E_re(2,0:Orb_1p1h%J_max))
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                E_re(P,J)%v = real(E(P,J)%v,dp)
            end do
        end do

        call timer_start(t0)
        call RPA_summary(Params,Orb_1p1h,E_corr,Stability,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        call RPA_binary_export(Orb_1p1h,E,X,Y,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        call RPA_spectrum_export(Orb_1p1h,E,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        write(*,'(/,a)') 'Preparing plot-ready export of RPA solutions ...'
        write(*,'(a)') TAB // 'Performing the export of plot-ready RPA spectra ...'
        call ph_plot_spectrum(Out_Dir // '/spectra/RPA_plot.dat',Orb_1p1h,E_re)
        write(*,'(a)') TAB // TAB // 'RPA plot-ready spectra successfully exported ...'
        call timer_report(t0)
        call timer_start(t0)
        call RPA_amplitudes_export(Orb,Orb_1p1h,E_re,X,Y,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        call rho_export(Params,Orb,Out_Dir // '/RPA_summary.dat', &
                        Out_Dir // '/densities/' // trim(Params%Calc%RPA%Name) // '_radial_densities.dat',Rho_LHO,C)
        call timer_report(t0)
        call timer_start(t0)
        call RPA_orbitals_export(Params,Orb,C,Rho_LHO,h,Out_Dir)
        call timer_report(t0)
        if (Params%Calc%RPA%TBCF) then
            call timer_start(t0)
            call O2b_export(Params,Orb_NN,Sigma,Out_Dir // '/bin/sigma.bin')
            call timer_report(t0)
        end if
        call timer_start(t0)
        call ph_transitions_export(Params,'RPA',Out_Dir,Params%Calc%RPA%Ortho,Orb,Orb_1p1h,C,O1b_transformation(Rho_LHO,C), &
                                   E_re,rB)
        call timer_report(t0)
    end subroutine RPA_export

    subroutine RPA_summary(Params,Orb_1p1h,E_corr,Stability,Out_Dir)
        type(params_t), intent(in) :: Params
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rpa_energy_t), intent(in) :: E_corr
        logical, intent(in) :: Stability
        character(len=*), intent(in) :: Out_Dir
        integer :: u

        write(*,'(/,a)') 'Writing RPA summary files ...'
        open(newunit=u, file=Out_Dir // '/RPA_summary.dat', status='replace', action='write')
        call ph_summary_header(u,Params,'Spherical Random-Phase Approximation review:')
        if (Params%Calc%RPA%Ortho) then
            write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (1-) is included ...'
        else
            write(u,'(/,a)') 'Orthogonalization of spurious Goldstone levels (1-) is NOT included ...'
        end if
        if (Stability) then
            write(u,'(a)') TAB // 'RPA system is STABLE ... no complex eigenvalues detected ...'
        else
            write(u,'(a)') TAB // 'RPA system is UNSTABLE ... complex eigenvalues detected ...'
        end if
        call ph_dimensions_block(u,Orb_1p1h,'RPA')
        write(u,'(/,/,a)') 'RPA ground-state correlation energy ...'
        write(u,'(/,a)',advance='no') 'E_corr     = ' // cfmt_f(E_corr%E_QBA,12,6) // ' ' // TAB // &
                                      ' MeV ... RPA ground-state correlation energy (quasi-boson)'
        if (.not. ieee_is_nan(E_corr%E_dens)) &
            write(u,'(/,a)',advance='no') 'E_corr     = ' // cfmt_f(E_corr%E_dens,12,6) // ' ' // TAB // &
                                          ' MeV ... RPA ground-state correlation energy (density)'
        write(u,'(a)') ''
        close(u)
        write(*,'(a)') TAB // 'RPA summary file exported ...'
    end subroutine RPA_summary

    ! bin/RPA_X.bin, RPA_Y.bin (row-major ComplexF64, permutedims per channel) & RPA_E.bin ...
    subroutine RPA_binary_export(Orb_1p1h,E,X,Y,Out_Dir)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cvec_t), intent(in) :: E(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        character(len=*), intent(in) :: Out_Dir
        integer :: u1, u2, u3, J, P

        write(*,'(/,a)') 'Performing binary export of RPA energies E & amplitudes X & Y ...'
        open(newunit=u1, file=Out_Dir // '/bin/RPA_X.bin', access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u2, file=Out_Dir // '/bin/RPA_Y.bin', access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u3, file=Out_Dir // '/bin/RPA_E.bin', access='stream', form='unformatted', status='replace', action='write')
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                if (Orb_1p1h%N(P,J) == 0) cycle
                write(u1) transpose(X(P,J)%v)
                write(u2) transpose(Y(P,J)%v)
                write(u3) E(P,J)%v
            end do
        end do
        close(u1)
        close(u2)
        close(u3)
        write(*,'(a)') TAB // 'Binary export of RPA energies E, amplitudes X & Y completed ...'
    end subroutine RPA_binary_export

    subroutine RPA_spectrum_export(Orb_1p1h,E,Out_Dir)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(cvec_t), intent(in) :: E(:,0:)
        character(len=*), intent(in) :: Out_Dir
        integer :: u, J, P, nu

        write(*,'(/,a)') 'Preparing export of RPA spectra ...'
        open(newunit=u, file=Out_Dir // '/spectra/RPA.dat', status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('Re E',20,.true.) // ' ' // &
                       cfmt_s('Im E',20,.true.)
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                do nu = 1, Orb_1p1h%N(P,J)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(pm_str(P),5,.true.) // ' ' // &
                                   cfmt_f(real(E(P,J)%v(nu),dp),20,6,.true.) // ' ' // cfmt_f(aimag(E(P,J)%v(nu)),20,6,.true.)
                end do
                write(u,'(a)') ''
                write(u,'(a)') ''
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of RPA spectra successfully finished ...'
    end subroutine RPA_spectrum_export

    ! amplitudes/RPA_amplitudes_ordered.dat & RPA_amplitudes.dat ... metric norms sum (n_h - n_p)|X|^2 & sum (n_h - n_p)|Y|^2 ...
    subroutine RPA_amplitudes_export(Orb,Orb_1p1h,E_re,X,Y,Out_Dir)
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E_re(:,0:)
        type(cmat_t), intent(in) :: X(:,0:), Y(:,0:)
        character(len=*), intent(in) :: Out_Dir
        real(dp), allocatable :: dn(:), L_E(:), L_x(:), L_y(:)
        integer, allocatable :: L_J(:), L_P(:), Sort(:)
        integer :: u, J, P, nu, ph, N_lev, cnt, k
        real(dp) :: x2, y2

        write(*,'(/,a)') 'Preparing export norms of the RPA amplitudes X & Y ...'
        N_lev = sum(Orb_1p1h%N)
        allocate(L_E(N_lev), L_x(N_lev), L_y(N_lev), L_J(N_lev), L_P(N_lev), Sort(N_lev))

        open(newunit=u, file=Out_Dir // '/amplitudes/RPA_amplitudes_ordered.dat', status='replace', action='write')
        write(u,'(a)') cfmt_s('E',12,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('|Y|^2',12,.true.)
        cnt = 0
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                write(u,'(/,a)') 'J = ' // istr(J) // ',' // TAB // 'P = ' // pm_str(P)
                dn = ph_metric(Orb,Orb_1p1h,P,J)
                do nu = 1, Orb_1p1h%N(P,J)
                    x2 = 0.0_dp
                    y2 = 0.0_dp
                    do ph = 1, Orb_1p1h%N(P,J)
                        x2 = x2 + dn(ph) * (real(X(P,J)%v(ph,nu),dp)**2 + aimag(X(P,J)%v(ph,nu))**2)
                        y2 = y2 + dn(ph) * (real(Y(P,J)%v(ph,nu),dp)**2 + aimag(Y(P,J)%v(ph,nu))**2)
                    end do
                    write(u,'(a)') cfmt_f(E_re(P,J)%v(nu),12,8,.true.) // ' ' // cfmt_f(x2,12,8,.true.) // ' ' // &
                                   cfmt_f(y2,12,8,.true.)
                    cnt = cnt + 1
                    L_E(cnt) = E_re(P,J)%v(nu)
                    L_x(cnt) = x2
                    L_y(cnt) = y2
                    L_J(cnt) = J
                    L_P(cnt) = P
                end do
            end do
        end do
        close(u)

        call sortperm(L_E,Sort)
        open(newunit=u, file=Out_Dir // '/amplitudes/RPA_amplitudes.dat', status='replace', action='write')
        write(u,'(a)') cfmt_s('Ind',12,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // cfmt_s('J',12,.true.) // ' ' // &
                       cfmt_s('P',12,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('|Y|^2',12,.true.)
        do nu = 1, N_lev
            k = Sort(nu)
            write(u,'(a)') cfmt_d(nu,12,.true.) // ' ' // cfmt_f(L_E(k),12,6,.true.) // ' ' // cfmt_d(L_J(k),12,.true.) // ' ' // &
                           cfmt_s(pm_str(L_P(k)),12,.true.) // ' ' // cfmt_f(L_x(k),12,6,.true.) // ' ' // cfmt_f(L_y(k),12,6,.true.)
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of RPA norms of X & Y completed ...'
    end subroutine RPA_amplitudes_export

    ! RPA natural orbitals (summary table, occupation file) & the binary orb/C/rho(/H1b) export ... counterpart of
        ! RPA_orbitals_export & RPA_natural_orbitals_export (RPA_export.jl:54-153) ...
    subroutine RPA_orbitals_export(Params,Orb,C,Rho_LHO,h,Out_Dir)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: C, Rho_LHO, h
        character(len=*), intent(in) :: Out_Dir
        type(o1b_t) :: Rho_RPA, C_NAT, Rho_diag, Rho_NAT, h_NAT, C_LHO_NAT
        type(orb1b_t), allocatable :: Orb_tmp(:), Orb_ref(:), Orb_NAT(:)
        real(dp), allocatable :: pO(:), nO(:), pE(:), nE(:)
        integer :: N_max, a_max, a, u

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        write(*,'(/,a)') 'Preparing export of RPA single-particle orbitals ...'

        Rho_RPA = O1b_transformation(Rho_LHO,C)
        call OBDM_natural_orbitals(Params,Orb,Rho_RPA,Orb_tmp,C_NAT,Rho_diag)
        Rho_NAT = O1b_transformation(Rho_RPA,C_NAT)
        h_NAT = O1b_transformation(h,C_NAT)
        allocate(pO(a_max), nO(a_max), pE(a_max), nE(a_max))
        do a = 1, a_max
            pO(a) = Rho_NAT%p(a,a)
            nO(a) = Rho_NAT%n(a,a)
            pE(a) = h_NAT%p(a,a)
            nE(a) = h_NAT%n(a,a)
        end do
        Orb_ref = orbitals_import(Params,io_dir(Params,Params%Calc%RPA%Reference) // '/bin/orb.bin')

        open(newunit=u, file=Out_Dir // '/RPA_summary.dat', position='append', action='write')
        write(u,'(/,a)') 'RPA single-particle orbitals ... natural orbitals of the RPA OBDM, E = <a|h|a> centroids of the ' // &
                         'mean-field Hamiltonian h ...'
        call table('Proton',pO,pE)
        call table('Neutron',nO,nE)
        close(u)

        open(newunit=u, file=Out_Dir // '/densities/RPA_occupation.dat', status='replace', action='write')
        write(u,'(a)') cfmt_s('a',6,.true.) // ' ' // cfmt_s('n',6,.true.) // ' ' // cfmt_s('l',6,.true.) // ' ' // &
                       cfmt_s('2j',6,.true.) // ' ' // cfmt_s('pO_ref',12,.true.) // ' ' // cfmt_s('pO_RPA',12,.true.) // ' ' // &
                       cfmt_s('pO_dep',12,.true.) // ' ' // cfmt_s('nO_ref',12,.true.) // ' ' // cfmt_s('nO_RPA',12,.true.) // ' ' // &
                       cfmt_s('nO_dep',12,.true.)
        do a = 1, a_max
            write(u,'(a)') cfmt_d(a,6,.true.) // ' ' // cfmt_d(Orb(a)%n,6,.true.) // ' ' // cfmt_d(Orb(a)%l,6,.true.) // ' ' // &
                           cfmt_d(Orb(a)%j,6,.true.) // ' ' // cfmt_f(Orb_ref(a)%pO,12,6,.true.) // ' ' // &
                           cfmt_f(pO(a),12,6,.true.) // ' ' // cfmt_f(pO(a) - Orb_ref(a)%pO,12,6,.true.) // ' ' // &
                           cfmt_f(Orb_ref(a)%nO,12,6,.true.) // ' ' // cfmt_f(nO(a),12,6,.true.) // ' ' // &
                           cfmt_f(nO(a) - Orb_ref(a)%nO,12,6,.true.)
        end do
        close(u)

        if (Params%Calc%RPA%RPA_Orb) then
            ! RPA natural orbitals as a new reference basis ...
            Rho_NAT = O1b_diagonal(O1b_transformation(Rho_RPA,C_NAT))
            Orb_NAT = Orb
            call orbitals_update(Orb_NAT,Rho_NAT)
            C_LHO_NAT%p = mm(C%p,C_NAT%p)
            C_LHO_NAT%n = mm(C%n,C_NAT%n)
            call orbitals_export(Params,Orb_NAT,Out_Dir // '/bin/orb.bin')
            call O1b_export(Params,Orb_NAT,C_LHO_NAT,Out_Dir // '/bin/C.bin')
            call O1b_export(Params,Orb_NAT,Rho_NAT,Out_Dir // '/bin/rho.bin')
            call O1b_export(Params,Orb_NAT,h_NAT,Out_Dir // '/bin/H1b.bin')
        else
            call orbitals_export(Params,Orb,Out_Dir // '/bin/orb.bin')
            call O1b_export(Params,Orb,C,Out_Dir // '/bin/C.bin')
            call O1b_export(Params,Orb,Rho_RPA,Out_Dir // '/bin/rho.bin')
        end if
        write(*,'(a)') TAB // 'Export of RPA single-particle orbitals finished ...'

    contains

        subroutine table(Label,O_t,E_t)
            character(len=*), intent(in) :: Label
            real(dp), intent(in) :: O_t(:), E_t(:)
            integer, allocatable :: n_count(:,:), Sort(:)
            integer :: i, b
            real(dp) :: Sum
            allocate(n_count(0:N_max,0:2*N_max+1), Sort(a_max))
            n_count = 0
            write(u,'(/,a)') Label // ' single-particle states'
            write(u,'(a)') '____________________________________________________________________'
            write(u,'(a)') cfmt_s('n',4) // ' ' // cfmt_s('l',4) // ' ' // cfmt_s('2j',4) // ' ' // cfmt_s('O_a',10) // ' ' // &
                           cfmt_s('E',14)
            call sortperm(E_t,Sort)
            do i = 1, a_max
                b = Sort(i)
                write(u,'(a)') cfmt_d(n_count(Orb(b)%l,Orb(b)%j),4) // ' ' // cfmt_d(Orb(b)%l,4) // ' ' // cfmt_d(Orb(b)%j,4) // &
                               ' ' // cfmt_f(O_t(b),10,6) // ' ' // cfmt_f(E_t(b),14,6)
                n_count(Orb(b)%l,Orb(b)%j) = n_count(Orb(b)%l,Orb(b)%j) + 1
            end do
            write(u,'(a)') '____________________________________________________________________'
            Sum = 0.0_dp
            do b = 1, a_max
                Sum = Sum + hat_of(Orb(b)%j) * O_t(b)
            end do
            write(u,'(a)') cfmt_s('Particle number ...',14,.true.) // ' ' // cfmt_f(Sum,15,6)
        end subroutine table

    end subroutine RPA_orbitals_export

end module nj_rpa_export
