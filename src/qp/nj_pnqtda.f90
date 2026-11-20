! Spherical proton-neutron quasiparticle Tamm-Dancoff approximation (charge-exchange 2qp states |a_p b_n; J P>) on top of
    ! an HFB/BCS reference ... counterpart of dev pnQTDA/*.jl (pnQTDA.jl, pnQTDA_solver.jl, pnQTDA_orbitals.jl,
    ! pnQTDA_allocate.jl, pnQTDA_diagonalize.jl, pnQTDA_export.jl) ... input io/<Path>/<Reference>/bin/{qpH1b,qpH2b}.bin,
    ! output io/<Path>/pnQTDA/{bin,spectra} ... no B matrix, no spurious states, no transitions, no densities ...
    ! pnQTDA_import_binary (pnQTDA_external.jl) has no caller in dev and is not ported ...
module nj_pnqtda
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_make
    use nj_o1b
    use nj_qpo2b, only : qpO2b_import, qpO2b_1111_pn
    use nj_qex, only : qex_reference_check
    !$ use omp_lib
    implicit none
    private

    public :: pnQTDA, orbitals_pn2qp_make

contains

    ! Driver ... counterpart of pnQTDA (pnQTDA.jl) ...
    subroutine pnQTDA(Params)
        type(params_t), intent(in) :: Params
        integer :: N_max, N_threads
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        N_max = Params%Calc%Nmax
        write(*,'(a)') 'Starting spherical pnQTDA calculation with residual NO2B NN+NNN interaction based on available ' // &
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
        if (.not. qex_reference_check(Params,Params%Calc%pnQTDA%Reference,[character(len=9) :: 'qpH1b.bin', 'qpH2b.bin'])) return
        if (is_dir(io_path(Params,'/pnQTDA'))) call remove_dir(io_path(Params,'/pnQTDA'))
        call make_dir(io_path(Params,'/pnQTDA'))
        call make_dir(io_path(Params,'/pnQTDA/bin'))
        call make_dir(io_path(Params,'/pnQTDA/spectra'))

        write(*,'(/,a)') 'Starting pnQTDA calculation ...'
        call timer_start(t0)
        call pnQTDA_solver(Params)
        call timer_report(t0)
        write(*,'(/,a,/)') 'Spherical pnQTDA solver run completed ...'
    end subroutine pnQTDA

    ! Counterpart of pnQTDA_solver (pnQTDA_solver.jl) ... C, U & V are not needed ...
    subroutine pnQTDA_solver(Params)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(orb2b_2qp_t) :: Orb_pn2qp
        type(qpo1b_t) :: H_N
        type(qpo2b_t) :: H_NN
        type(orb2b_t) :: Orb_NN
        type(rmat_t), allocatable :: Amat(:,:), X(:,:)
        type(rvec_t), allocatable :: E(:,:)
        character(len=:), allocatable :: Ref

        ! The pn pair set is not truncated by N2max, while H^22 pn(11,11) exists only for pairs within N2max (dev: KeyError
            ! of O2b_lookup in qpO2b_1111_pn) ... fatal with a clear message instead ...
        if (Params%Calc%N2max < 2*Params%Calc%Nmax) call fatal('pnQTDA: N2max = ' // istr(Params%Calc%N2max) // &
            ' < 2*Nmax = ' // istr(2*Params%Calc%Nmax) // ' is not supported (the pn 2qp pairs need H^22 pn(11,11) ' // &
            'up to 2*Nmax; dev throws a KeyError here)')

        Ref = '/' // trim(Params%Calc%pnQTDA%Reference) // '/bin/'
        Orb = orbitals_make(Params)
        Orb_pn2qp = orbitals_pn2qp_make(Params,Orb)
        H_N = qpO1B_import(Params,Orb,io_path(Params,Ref // 'qpH1b.bin'))
        call qpO2b_import(Params,Orb,io_path(Params,Ref // 'qpH2b.bin'),H_NN,Orb_NN)

        call pnQTDA_allocate(Params,Orb_NN,Orb_pn2qp,H_N,H_NN,Amat)
        call pnQTDA_diagonalize(Params,Amat,E,X)
        call pnQTDA_export(Params,Orb_pn2qp,E,X)
        write(*,'(a,/)') TAB // 'pnQTDA solver executed properly ...'
    end subroutine pnQTDA_solver

    ! pn 2qp orbitals ... counterpart of orbitals_pn2qp_make (pnQTDA_orbitals.jl) ... a proton, b neutron, all (a,b),
        ! J = |j_a - j_b|/2 .. (j_a + j_b)/2, T = 0, no N2max truncation ... the dense Dic is not built (dev Dict unused) ...
    function orbitals_pn2qp_make(Params,Orb) result(Orb_pn2qp)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_2qp_t) :: Orb_pn2qp
        integer :: N_max, a_max, J_max, a, b, J, P, Pass

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        J_max = 2*N_max + 1
        Orb_pn2qp%a_max = a_max
        Orb_pn2qp%J_max = J_max
        allocate(Orb_pn2qp%N(2,0:J_max), Orb_pn2qp%Ind(2,0:J_max))
        do Pass = 1, 2
            Orb_pn2qp%N = 0
            do a = 1, a_max
                do b = 1, a_max
                    P = parity_of(Orb(a)%l + Orb(b)%l)
                    do J = abs(Orb(a)%j - Orb(b)%j)/2, (Orb(a)%j + Orb(b)%j)/2
                        Orb_pn2qp%N(P,J) = Orb_pn2qp%N(P,J) + 1
                        if (Pass == 2) Orb_pn2qp%Ind(P,J)%v(Orb_pn2qp%N(P,J)) = ind2b_t(a,b,0)
                    end do
                end do
            end do
            if (Pass == 1) then
                do P = 1, 2
                    do J = 0, J_max
                        allocate(Orb_pn2qp%Ind(P,J)%v(Orb_pn2qp%N(P,J)))
                    end do
                end do
            end if
        end do
    end function orbitals_pn2qp_make

    ! A(i,j) = delta_ac delta_bd (E_a^p + E_b^n) + H^22 pn(11,11)(a,b,c,d;J) ... counterpart of pnQTDA_allocate ... one
        ! (J,P) block per task (each block written by one thread only) ...
    subroutine pnQTDA_allocate(Params,Orb_NN,Orb_pn2qp,H_N,H_NN,Amat)
        type(params_t), intent(in) :: Params
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb2b_2qp_t), intent(in) :: Orb_pn2qp
        type(qpo1b_t), intent(in) :: H_N
        type(qpo2b_t), intent(in) :: H_NN
        type(rmat_t), allocatable, intent(out) :: Amat(:,:)
        integer :: J_max, JP, J, P, N_qp, i_qp, j_qp, a, b, c, d
        real(dp) :: ASum

        J_max = 2*Params%Calc%Nmax + 1
        allocate(Amat(2,0:J_max))
        do J = 0, J_max
            do P = 1, 2
                allocate(Amat(P,J)%v(Orb_pn2qp%N(P,J),Orb_pn2qp%N(P,J)))
            end do
        end do
        write(*,'(/,a)') 'Allocating the pnQTDA matrix A ...'
        !$omp parallel do schedule(dynamic) private(J,P,N_qp,i_qp,j_qp,a,b,c,d,ASum)
        do JP = 1, 2*(J_max + 1)
            J = (JP - 1)/2
            P = mod(JP - 1,2) + 1
            N_qp = Orb_pn2qp%N(P,J)
            do i_qp = 1, N_qp
                a = Orb_pn2qp%Ind(P,J)%v(i_qp)%a
                b = Orb_pn2qp%Ind(P,J)%v(i_qp)%b
                do j_qp = 1, i_qp
                    c = Orb_pn2qp%Ind(P,J)%v(j_qp)%a
                    d = Orb_pn2qp%Ind(P,J)%v(j_qp)%b
                    ASum = 0.0_dp
                    if (a == c .and. b == d) ASum = ASum + (H_N%qp11%p(a,a) + H_N%qp11%n(b,b))
                    ASum = ASum + qpO2b_1111_pn(a,b,c,d,J,P,H_NN,Orb_NN)
                    Amat(P,J)%v(i_qp,j_qp) = ASum
                    if (i_qp /= j_qp) Amat(P,J)%v(j_qp,i_qp) = ASum
                end do
            end do
        end do
        !$omp end parallel do
        write(*,'(/,a)') 'The pnQTDA matrix Amat succesfully allocated ...'
    end subroutine pnQTDA_allocate

    ! eigen!(A, sortby = +) of the exactly symmetric blocks ... counterpart of pnQTDA_diagonalize ...
    subroutine pnQTDA_diagonalize(Params,Amat,E,X)
        type(params_t), intent(in) :: Params
        type(rmat_t), intent(in) :: Amat(:,0:)
        type(rvec_t), allocatable, intent(out) :: E(:,:)
        type(rmat_t), allocatable, intent(out) :: X(:,:)
        integer :: J_max, J, P, N

        J_max = 2*Params%Calc%Nmax + 1
        allocate(E(2,0:J_max), X(2,0:J_max))
        write(*,'(/,a)') 'Diagonalizing the pnQTDA matrix Amat ...'
        do J = 0, J_max
            do P = 1, 2
                N = size(Amat(P,J)%v,1)
                allocate(E(P,J)%v(N), X(P,J)%v(N,N))
                if (N > 0) call eig_sym(Amat(P,J)%v,E(P,J)%v,X(P,J)%v)
            end do
        end do
        write(*,'(/,a)') 'The pnQTDA matrix Amat succesfully diagonalized ...'
    end subroutine pnQTDA_diagonalize

    ! Summary, binaries (X column-major, then E; J outer, P inner), spectrum & plot-ready spectrum ... pnQTDA_export.jl ...
    subroutine pnQTDA_export(Params,Orb_pn2qp,E,X)
        type(params_t), intent(in) :: Params
        type(orb2b_2qp_t), intent(in) :: Orb_pn2qp
        type(rvec_t), intent(in) :: E(:,0:)
        type(rmat_t), intent(in) :: X(:,0:)
        integer :: N_max, J_max, u, J, P, nu, N_qp, qp_count
        real(dp), allocatable :: Ex(:,:)
        integer, allocatable :: Sort(:)
        character(len=1) :: Ps

        N_max = Params%Calc%Nmax
        J_max = 2*N_max + 1

        ! Summary ...
        write(*,'(/,a)') 'Preparing pnQTDA summary file ...'
        open(newunit=u, file=io_path(Params,'/pnQTDA/pnQTDA_summary.dat'), status='replace', action='write')
        write(u,'(a)') 'Spherical proton-neutron Quasiparticle Tamm-Dancoff Approximation review:'
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
        write(u,'(/,a)') "Reference quasiparticle mean-field solution:   Reference = ''" // trim(Params%Calc%pnQTDA%Reference) // "''"
        N_qp = 0
        do J = 0, J_max
            do P = 1, 2
                N_qp = N_qp + Orb_pn2qp%N(P,J) * (2*J + 1)
            end do
        end do
        write(u,'(/,a)') 'Information on dimensions of pnQTDA 1-phonon subspaces ...'
        write(u,'(/,a)') 'Total number of 1-phonon pn quasiparticle states in J-scheme = ' // istr(sum(Orb_pn2qp%N))
        write(u,'(a)') 'Total number of 1-phonon pn quasiparticle states in M-scheme = ' // istr(N_qp)
        write(u,'(/,a)') 'J' // TAB // 'N_phonon'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_pn2qp%N(1,J) + Orb_pn2qp%N(2,J))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = +'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_pn2qp%N(1,J))
        end do
        write(u,'(/,a)') 'Subspace' // TAB // 'P = -'
        do J = 0, J_max
            write(u,'(a)') istr(J) // TAB // istr(Orb_pn2qp%N(2,J))
        end do
        close(u)
        write(*,'(/,a)') 'pnQTDA summary exported ...'

        ! Binary export ...
        write(*,'(/,a)') 'Performing binary export of pnQTDA energies E & amplitudes X ...'
        open(newunit=u, file=io_path(Params,'/pnQTDA/bin/pnQTDA_X.bin'), access='stream', form='unformatted', &
             status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(X(P,J)%v) > 0) write(u) X(P,J)%v
            end do
        end do
        close(u)
        open(newunit=u, file=io_path(Params,'/pnQTDA/bin/pnQTDA_E.bin'), access='stream', form='unformatted', &
             status='replace', action='write')
        do J = 0, J_max
            do P = 1, 2
                if (size(E(P,J)%v) > 0) write(u) E(P,J)%v
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Binary export of pnQTDA energies E amplitudes X completed ...'

        ! Spectrum ...
        write(*,'(/,a)') 'Preparing export of pnQTDA spectra ...'
        open(newunit=u, file=io_path(Params,'/pnQTDA/spectra/pnQTDA.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',10,.true.)
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, Orb_pn2qp%N(P,J)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(merge('+','-',P == 1),5,.true.) // ' ' // &
                                   cfmt_f(E(P,J)%v(nu),10,6,.true.)
                end do
                write(u,'(/)')
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of pnQTDA spectra succesfully finished ...'

        ! Plot-ready spectrum ... no ground-state row, all N rows are written (unlike QTDA/QRPA) ...
        write(*,'(/,a)') 'Preparing plot-ready export of pnQTDA spectra ...'
        N_qp = sum(Orb_pn2qp%N)
        allocate(Ex(4,N_qp), Sort(N_qp))
        qp_count = 0
        do J = 0, J_max
            do P = 1, 2
                do nu = 1, Orb_pn2qp%N(P,J)
                    qp_count = qp_count + 1
                    Ex(:,qp_count) = [real(J,dp), real(P,dp), E(P,J)%v(nu), E(P,J)%v(nu)]
                end do
            end do
        end do
        if (N_qp > 0) then
            call sortperm(Ex(3,:),Sort)
            Ex = Ex(:,Sort)
        end if
        do nu = 2, N_qp
            if (Ex(4,nu) < Ex(4,nu-1)) then
                Ex(4,nu) = Ex(4,nu-1) + 1.0_dp
            else if (abs(Ex(4,nu) - Ex(4,nu-1)) < 1.0_dp) then
                Ex(4,nu) = Ex(4,nu-1) + 1.0_dp
            end if
        end do
        write(*,'(/,a)') 'Performing the export of plot-ready pnQTDA spectra ...'
        open(newunit=u, file=io_path(Params,'/pnQTDA/spectra/pnQTDA_plot.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',12,.true.) // ' ' // &
                       cfmt_s('label',12,.true.)
        do nu = 1, N_qp
            Ps = merge('+', '-', abs(Ex(2,nu) - 1.0_dp) < 1.0e-3_dp)
            write(u,'(a)') cfmt_d(nint(Ex(1,nu)),5,.true.) // ' ' // cfmt_s(Ps,5,.true.) // ' ' // cfmt_f(Ex(3,nu),12,6,.true.) // &
                           ' ' // cfmt_f(Ex(4,nu),12,6,.true.)
        end do
        close(u)
        write(*,'(a)') TAB // 'pnQTDA plot-ready spectra successfully exported ...'
    end subroutine pnQTDA_export

end module nj_pnqtda
