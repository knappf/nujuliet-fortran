! External / post-processing routines of stored TDA & RPA solutions (NuJuliet dev) ... counterpart of src/TDA/TDA_external.jl
    ! & src/RPA/RPA_external.jl: transition densities, transition operator densities, analysis of the phonon levels,
    ! transition current maps & recalculation of the transitions ... (&run actions TDA_/RPA_transition_densities,
    ! _transition_operator_densities, _transition_analysis, _transition_currents, _transitions_recalculate) ...
    !
    ! Inputs: TDA ... io/<Path>/<Reference>/bin/{orb,C,rho[,sigma]}.bin & io/<Path>/<Name>/bin/TDA_{X,E}.bin ,
    !         RPA ... io/<Path>/<Name>/bin/{orb,C,rho[,sigma],RPA_X,RPA_Y,RPA_E}.bin (RPA_Orb = .false. only) ...
    ! Both share one transition density matrix Gamma[(a,b,t),nu] = D_X X + D_Y Y (TDA: the exact kernel of TDA_rM, no Y;
    ! RPA: RPA_transition_kernel, Linear or Exact incl. the RPA Sigma) ... the TDA & RPA variants differ only where dev
    ! differs (gauge of real/complex amplitudes, the phase & normalization of the current pair maps, the analysis columns).
    ! TDA_transition_operator_densities evaluates <a|r^X Y_J|b> for a < b as Amp(a,b) directly instead of dev's mirrored
    ! (-1)^((j_a-j_b)/2) Amp(b,a) (identical up to round-off, no a_max^2 x N_grid buffer) ...
    ! Robustness guards (docs/PORTING_NOTES.md): fatal for a JP / nu outside the stored solutions & for missing input files
    ! (dev: BoundsError / SystemError), levels/, densities/ & transitions/E0..E3 are created with mkdir -p ...
module nj_ph_external
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : wigner_init_float, fCG, f6j
    use nj_functions, only : Psi_rad_LHO, differentiate_Oh6, J1b_convective, J1b_spin, phonon_weights
    use nj_orbitals, only : orbitals_import, orbitals_1p1h_make
    use nj_o1b, only : O1b_import
    use nj_o2b, only : O2b_import, O2b_pp, O2b_pn, O2b_nn
    use nj_obdm, only : rho_chR2
    use nj_tr1b
    use nj_ph_common
    use nj_tda, only : TDA_transition_kernel, TDA_rM
    use nj_rpa, only : RPA_rM
    !$ use omp_lib
    implicit none
    private

    public :: TDA_transition_densities, TDA_transition_operator_densities, TDA_transition_analysis, TDA_transition_currents
    public :: TDA_transitions_recalculate
    public :: RPA_transition_densities, RPA_transition_operator_densities, RPA_transition_analysis, RPA_transition_currents
    public :: RPA_transitions_recalculate

    real(dp), parameter :: hc = 197.326980_dp, m_n = 939.565346_dp, m_p = 938.272013_dp
    real(dp), parameter :: g_n = -3.82608545_dp, g_p = 5.585694713_dp

    ! Imported TDA or RPA data ... (Kind = 'TDA' or 'RPA') ...
    type :: ext_data_t
        character(len=3) :: Kind = ''
        character(len=:), allocatable :: IO
        type(orb1b_t), allocatable :: Orb(:)
        type(orb2b_1p1h_t) :: Orb_1p1h
        type(o1b_t) :: C, Rho
        type(o2b_t) :: Sigma
        type(orb2b_t) :: Orb_NN
        type(tr1b_t) :: TrOp
        type(rmat_t), allocatable :: X_T(:,:)
        type(rvec_t), allocatable :: E_T(:,:)
        type(cmat_t), allocatable :: X(:,:), Y(:,:)
        type(cvec_t), allocatable :: E(:,:)
    end type ext_data_t

    ! Transition density matrix of the selected levels of one multipole ...
    type :: ext_gamma_t
        integer :: J = 0, P = 0
        type(ind2b_t), allocatable :: Pairs(:)
        real(dp), allocatable :: D_X(:,:), D_Y(:,:), f_N(:)
        complex(dp), allocatable :: Gamma(:,:)
        complex(dp), allocatable :: pM(:), nM(:), isM(:), ivM(:)
    end type ext_gamma_t

contains

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Public entry points ...
    ! ---------------------------------------------------------------------------------------------------------------------

    subroutine TDA_transition_densities(Params,JP,nu_list,File_Name)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:)
        call ext_densities(Params,'TDA',JP,nu_list,File_Name,.false.)
    end subroutine TDA_transition_densities

    subroutine RPA_transition_densities(Params,JP,nu_list,File_Name)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:)
        call ext_densities(Params,'RPA',JP,nu_list,File_Name,.false.)
    end subroutine RPA_transition_densities

    subroutine TDA_transition_operator_densities(Params,JP,nu_list,File_Name)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:)
        call ext_densities(Params,'TDA',JP,nu_list,File_Name,.true.)
    end subroutine TDA_transition_operator_densities

    subroutine RPA_transition_operator_densities(Params,JP,nu_list,File_Name)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:)
        call ext_densities(Params,'RPA',JP,nu_list,File_Name,.true.)
    end subroutine RPA_transition_operator_densities

    subroutine TDA_transition_analysis(Params,JP,nu_list,File_Name,Cut)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:)
        real(dp), intent(in) :: Cut
        call ext_analysis(Params,'TDA',JP,nu_list,File_Name,Cut)
    end subroutine TDA_transition_analysis

    subroutine RPA_transition_analysis(Params,JP,nu_list,File_Name,Cut)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:)
        real(dp), intent(in) :: Cut
        call ext_analysis(Params,'RPA',JP,nu_list,File_Name,Cut)
    end subroutine RPA_transition_analysis

    subroutine TDA_transition_currents(Params,JP,nu_list,File_Name,M,xN_grid,zN_grid)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:), M, xN_grid, zN_grid
        call ext_currents(Params,'TDA',JP,nu_list,File_Name,M,xN_grid,zN_grid)
    end subroutine TDA_transition_currents

    subroutine RPA_transition_currents(Params,JP,nu_list,File_Name,M,xN_grid,zN_grid)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: JP, File_Name
        integer, intent(in) :: nu_list(:), M, xN_grid, zN_grid
        call ext_currents(Params,'RPA',JP,nu_list,File_Name,M,xN_grid,zN_grid)
    end subroutine RPA_transition_currents

    ! Recalculation of the TDA transitions from the stored solution ... counterpart of TDA_transitions_recalculate ...
    subroutine TDA_transitions_recalculate(Params)
        type(params_t), intent(in) :: Params
        type(ext_data_t) :: Dt
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        call ext_import(Params,'TDA',Dt)
        call transitions_dirs(Dt%IO)
        call timer_start(t0)
        rM = TDA_rM(Params,Dt%Orb,Dt%Orb_1p1h,Dt%Orb_NN,Dt%TrOp,Dt%Rho,Dt%Sigma,Dt%X_T)
        call timer_report(t0)
        call timer_start(t0)
        rB = ph_rB(Params,Dt%Orb_1p1h,rM,Dt%E_T(2,1)%v)
        call timer_report(t0)
        call summary_truncate(Dt%IO // '/TDA_summary.dat',NEW_LINE('a') // 'Review of TDA electromagnetic transition moments')
        call timer_start(t0)
        call ph_transitions_export(Params,'TDA',Dt%IO,Params%Calc%TDA%Ortho,Dt%Orb,Dt%Orb_1p1h,Dt%C,Dt%Rho,Dt%E_T,rB)
        call timer_report(t0)
    end subroutine TDA_transitions_recalculate

    ! Recalculation of the RPA transitions from the stored solution ... counterpart of RPA_transitions_recalculate ...
    subroutine RPA_transitions_recalculate(Params)
        type(params_t), intent(in) :: Params
        type(ext_data_t) :: Dt
        type(reduced_multipole_t) :: rM
        type(reduced_transition_t) :: rB
        type(rvec_t), allocatable :: E_re(:,:)
        integer :: J, P
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        call ext_import(Params,'RPA',Dt)
        call transitions_dirs(Dt%IO)
        call timer_start(t0)
        rM = RPA_rM(Params,Dt%Orb,Dt%Orb_1p1h,Dt%Orb_NN,Dt%TrOp,Dt%Rho,Dt%Sigma,Dt%X,Dt%Y)
        call timer_report(t0)
        allocate(E_re(2,0:Dt%Orb_1p1h%J_max))
        do J = 0, Dt%Orb_1p1h%J_max
            do P = 1, 2
                E_re(P,J)%v = real(Dt%E(P,J)%v,dp)
            end do
        end do
        call timer_start(t0)
        rB = ph_rB(Params,Dt%Orb_1p1h,rM,E_re(2,1)%v)
        call timer_report(t0)
        call summary_truncate(Dt%IO // '/RPA_summary.dat',NEW_LINE('a') // 'Review of RPA electromagnetic transition moments')
        ! (dev passes Rho_LHO = C Rho C^T, which RPA_transitions_export transforms back to the RPA basis) ...
        call timer_start(t0)
        call ph_transitions_export(Params,'RPA',Dt%IO,Params%Calc%RPA%Ortho,Dt%Orb,Dt%Orb_1p1h,Dt%C,Dt%Rho,E_re,rB)
        call timer_report(t0)
    end subroutine RPA_transitions_recalculate

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Import ...
    ! ---------------------------------------------------------------------------------------------------------------------

    ! Counterpart of TDA_import_binary + the import blocks of TDA_external.jl & of RPA_external_import / RPA_import_binary ...
    subroutine ext_import(Params,Kind,Dt)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind
        type(ext_data_t), intent(out) :: Dt
        character(len=:), allocatable :: IO_in
        logical :: Read_Sigma
        real(dp) :: chR2
        integer :: J, P, N, u1, u2, u3
        real(dp), allocatable :: Buf_r(:,:)
        complex(dp), allocatable :: Buf_c(:,:)

        Dt%Kind = Kind
        if (Kind == 'TDA') then
            Dt%IO = io_dir(Params,Params%Calc%TDA%Name)
            IO_in = io_dir(Params,Params%Calc%TDA%Reference)
            Read_Sigma = Params%Calc%TDA%Reference_TBCF
        else
            if (trim(Params%Calc%RPA%Transitions) /= 'Exact' .and. trim(Params%Calc%RPA%Transitions) /= 'Linear') &
                call fatal('Invalid Params.Calc.RPA.Transitions = "' // trim(Params%Calc%RPA%Transitions) // &
                           '"; use "Exact" or "Linear".')
            if (Params%Calc%RPA%RPA_Orb) &
                call fatal('RPA external routines: Params.Calc.RPA.RPA_Orb = true ... orb.bin, C.bin & rho.bin of the RPA ' // &
                           'directory are given in the RPA natural-orbital basis, which does not match the basis of the RPA ' // &
                           'amplitudes X & Y ... use the RPA solution exported with RPA_Orb = false.')
            Dt%IO = io_dir(Params,Params%Calc%RPA%Name)
            IO_in = Dt%IO
            Read_Sigma = (trim(Params%Calc%RPA%Transitions) == 'Exact') .and. Params%Calc%RPA%TBCF
        end if

        Dt%Orb = orbitals_import(Params,need(IO_in // '/bin/orb.bin'))
        if (Kind == 'TDA') then
            Dt%Orb_1p1h = orbitals_1p1h_make(Params,Dt%Orb,Params%Calc%TDA%Tol_par,Params%Calc%TDA%Tol_hol)
        else
            Dt%Orb_1p1h = orbitals_1p1h_make(Params,Dt%Orb,Params%Calc%RPA%Tol_par,Params%Calc%RPA%Tol_hol)
        end if
        Dt%C = O1b_import(Params,Dt%Orb,need(IO_in // '/bin/C.bin'))
        Dt%Rho = O1b_import(Params,Dt%Orb,need(IO_in // '/bin/rho.bin'))
        if (Read_Sigma) then
            if (Kind == 'TDA') then
                write(*,'(/,a)') 'Importing reference 2-body correlation function ...'
            else
                write(*,'(/,a)') 'Importing RPA 2-body correlation function ...'
            end if
            call O2b_import(Params,Dt%Orb,need(IO_in // '/bin/sigma.bin'),Dt%Sigma,Dt%Orb_NN)
        end if

        ! Solutions (row-major binaries, one block per non-empty channel in the dev JP order) ...
        write(*,'(/,a)') 'Importing ' // Kind // ' solutions from binary format ...'
        if (Kind == 'TDA') then
            allocate(Dt%X_T(2,0:Dt%Orb_1p1h%J_max), Dt%E_T(2,0:Dt%Orb_1p1h%J_max))
            open(newunit=u1, file=need(Dt%IO // '/bin/TDA_X.bin'), access='stream', form='unformatted', status='old', action='read')
            open(newunit=u2, file=need(Dt%IO // '/bin/TDA_E.bin'), access='stream', form='unformatted', status='old', action='read')
            do J = 0, Dt%Orb_1p1h%J_max
                do P = 1, 2
                    N = Dt%Orb_1p1h%N(P,J)
                    allocate(Dt%X_T(P,J)%v(N,N), Dt%E_T(P,J)%v(N))
                    if (N == 0) cycle
                    allocate(Buf_r(N,N))
                    read(u1) Buf_r
                    Dt%X_T(P,J)%v = transpose(Buf_r)
                    deallocate(Buf_r)
                    read(u2) Dt%E_T(P,J)%v
                end do
            end do
            close(u1)
            close(u2)
        else
            allocate(Dt%X(2,0:Dt%Orb_1p1h%J_max), Dt%Y(2,0:Dt%Orb_1p1h%J_max), Dt%E(2,0:Dt%Orb_1p1h%J_max))
            open(newunit=u1, file=need(Dt%IO // '/bin/RPA_X.bin'), access='stream', form='unformatted', status='old', action='read')
            open(newunit=u2, file=need(Dt%IO // '/bin/RPA_Y.bin'), access='stream', form='unformatted', status='old', action='read')
            open(newunit=u3, file=need(Dt%IO // '/bin/RPA_E.bin'), access='stream', form='unformatted', status='old', action='read')
            do J = 0, Dt%Orb_1p1h%J_max
                do P = 1, 2
                    N = Dt%Orb_1p1h%N(P,J)
                    allocate(Dt%X(P,J)%v(N,N), Dt%Y(P,J)%v(N,N), Dt%E(P,J)%v(N))
                    if (N == 0) cycle
                    allocate(Buf_c(N,N))
                    read(u1) Buf_c
                    Dt%X(P,J)%v = transpose(Buf_c)
                    read(u2) Buf_c
                    Dt%Y(P,J)%v = transpose(Buf_c)
                    deallocate(Buf_c)
                    read(u3) Dt%E(P,J)%v
                end do
            end do
            close(u1)
            close(u2)
            close(u3)
        end if
        write(*,'(a)') TAB // Kind // ' solution binaries have been succesufully imported ...'

        ! Electromagnetic operators with the charge radius of the imported OBDM (as in TDA_solver / RPA_solver) ...
        chR2 = rho_chR2(Params,Dt%Orb,Dt%C,Dt%Rho)
        Dt%TrOp = Tr1b_initialize(Params,Dt%Orb,chR2)
        Dt%TrOp = Tr1b_transformation(Params,Dt%Orb,Dt%C,Dt%TrOp)
    end subroutine ext_import

    ! Input file check (robustness guard) ...
    function need(Path) result(P)
        character(len=*), intent(in) :: Path
        character(len=:), allocatable :: P
        logical :: Ex
        inquire(file=Path, exist=Ex)
        if (.not. Ex) call fatal('External ' // 'TDA/RPA routine: input file ' // Path // ' does not exist ...')
        P = Path
    end function need

    ! "J+" / "J-" -> (J, P) ... counterpart of RPA_external_JP ... + range guard of the stored channels ...
    subroutine ext_JP(JP,Dt,J,P)
        character(len=*), intent(in) :: JP
        type(ext_data_t), intent(in) :: Dt
        integer, intent(out) :: J, P
        character(len=:), allocatable :: S
        integer :: Ios
        S = trim(adjustl(JP))
        if (len(S) < 2) call fatal('Invalid JP = "' // S // '"; use J+ or J- ...')
        if (S(len(S):len(S)) == '+') then
            P = 1
        else if (S(len(S):len(S)) == '-') then
            P = 2
        else
            call fatal('Invalid JP = "' // S // '"; use J+ or J- ...')
        end if
        read(S(1:len(S)-1),*,iostat=Ios) J
        if (Ios /= 0) call fatal('Invalid JP = "' // S // '"; use J+ or J- ...')
        if (J < 0 .or. J > Dt%Orb_1p1h%J_max) call fatal('JP = "' // S // '" is outside the stored ' // Dt%Kind // ' channels ...')
    end subroutine ext_JP

    subroutine check_nu(Dt,J,P,nu_list)
        type(ext_data_t), intent(in) :: Dt
        integer, intent(in) :: J, P, nu_list(:)
        integer :: k
        if (size(nu_list) == 0) call fatal('Empty nu_list ...')
        do k = 1, size(nu_list)
            if (nu_list(k) < 1 .or. nu_list(k) > Dt%Orb_1p1h%N(P,J)) &
                call fatal('nu = ' // istr(nu_list(k)) // ' is outside the ' // istr(Dt%Orb_1p1h%N(P,J)) // ' ' // Dt%Kind // &
                           ' levels of JP = ' // istr(J) // pm_str(P) // ' ...')
        end do
    end subroutine check_nu

    ! sE = string(round(E,digits=3)) of level nu ...
    function level_E(Dt,J,P,nu) result(s)
        type(ext_data_t), intent(in) :: Dt
        integer, intent(in) :: J, P, nu
        character(len=:), allocatable :: s
        if (Dt%Kind == 'TDA') then
            s = jstr(round_digits(Dt%E_T(P,J)%v(nu),3))
        else
            s = jstr(round_digits(real(Dt%E(P,J)%v(nu),dp),3))
        end if
    end function level_E

    pure logical function is_EJ(J,P)
        integer, intent(in) :: J, P
        is_EJ = (J == 0 .and. P == 1) .or. (J == 1 .and. P == 2) .or. (J == 2 .and. P == 1) .or. (J == 3 .and. P == 2)
    end function is_EJ

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Kernels & multipoles ...
    ! ---------------------------------------------------------------------------------------------------------------------

    ! RPA transition kernels D_X & D_Y on all (a,b,t) of the multipole, only non-vanishing pairs kept ... counterpart of
        ! RPA_transition_kernel (RPA_external.jl:173-321) ... (threads over the ph columns, each writes its own column)
    subroutine RPA_transition_kernel(Params,J,P,Orb,Orb_1p1h,Orb_NN,Rho,Sigma,Pairs,D_X,D_Y)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: J, P
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(in) :: Sigma
        type(ind2b_t), allocatable, intent(out) :: Pairs(:)
        real(dp), allocatable, intent(out) :: D_X(:,:), D_Y(:,:)
        type(ind2b_t), allocatable :: Pairs_all(:)
        integer, allocatable :: Pair_Map(:,:,:)
        real(dp), allocatable :: D_X_all(:,:), D_Y_all(:,:)
        logical, allocatable :: Keep(:)
        logical :: Exact, TBCF
        integer :: a_max, N_ph, N_all, t, ia, ib, i_ph, a_p, a_h, T_ph, l_p, j_p, l_h, j_h, t_1, t_2, c, d, k
        integer :: a_a, l_a, j_a, a_b, l_b, j_b, P_ha, I
        real(dp) :: Phase_hp, D_ph, Phase_X, Phase_Y, K_1, K_2, F, Phase_ab

        Exact = (trim(Params%Calc%RPA%Transitions) == 'Exact')
        TBCF = Exact .and. Params%Calc%RPA%TBCF
        N_ph = Orb_1p1h%N(P,J)
        a_max = size(Orb)

        allocate(Pairs_all(2*a_max*a_max), Pair_Map(a_max,a_max,2))
        Pair_Map = 0
        N_all = 0
        do t = -1, 1, 2
            do ia = 1, a_max
                do ib = 1, a_max
                    if (parity_of(Orb(ia)%l + Orb(ib)%l) == P .and. triangle(Orb(ia)%j,Orb(ib)%j,2*J)) then
                        N_all = N_all + 1
                        Pairs_all(N_all) = ind2b_t(Orb(ia)%a,Orb(ib)%a,t)
                        Pair_Map(Orb(ia)%a,Orb(ib)%a,(t + 3)/2) = N_all
                    end if
                end do
            end do
        end do

        allocate(D_X_all(N_all,N_ph), D_Y_all(N_all,N_ph))
        D_X_all = 0.0_dp
        D_Y_all = 0.0_dp

        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(a_p,a_h,T_ph,l_p,j_p,l_h,j_h,t_1,t_2,Phase_hp,D_ph,c,d,Phase_X,Phase_Y,ia,ib,a_a,l_a,j_a,a_b,l_b,j_b) &
        !$omp private(P_ha,K_1,K_2,I,F,Phase_ab)
        do i_ph = 1, N_ph
            a_p = Orb_1p1h%Ind(P,J)%v(i_ph)%a
            a_h = Orb_1p1h%Ind(P,J)%v(i_ph)%b
            T_ph = Orb_1p1h%Ind(P,J)%v(i_ph)%T
            if (T_ph /= -1 .and. T_ph /= 1) cycle
            l_p = Orb(a_p)%l
            j_p = Orb(a_p)%j
            l_h = Orb(a_h)%l
            j_h = Orb(a_h)%j
            t_1 = (T_ph + 3)/2
            t_2 = (-T_ph + 3)/2
            Phase_hp = phase((j_p - j_h)/2)

            if (.not. Exact) then
                ! Linearized transition densities ...
                D_ph = ph_occ(Orb,a_h,T_ph) - ph_occ(Orb,a_p,T_ph)
                D_X_all(Pair_Map(a_p,a_h,t_1),i_ph) = D_X_all(Pair_Map(a_p,a_h,t_1),i_ph) + phase(J) * D_ph
                D_Y_all(Pair_Map(a_h,a_p,t_1),i_ph) = D_Y_all(Pair_Map(a_h,a_p,t_1),i_ph) + Phase_hp * D_ph
            else
                ! 1-body part ... (1 - Rho) M Rho ...
                do c = 1, a_max
                    do d = 1, a_max
                        if (Orb(c)%l == l_p .and. Orb(c)%j == j_p .and. Orb(d)%l == l_h .and. Orb(d)%j == j_h) &
                            D_X_all(Pair_Map(c,d,t_1),i_ph) = D_X_all(Pair_Map(c,d,t_1),i_ph) + &
                                phase(J) * (real(kronecker_delta(c,a_p),dp) - rho_t(a_p,c,T_ph)) * rho_t(d,a_h,T_ph)
                        if (Orb(c)%l == l_h .and. Orb(c)%j == j_h .and. Orb(d)%l == l_p .and. Orb(d)%j == j_p) &
                            D_Y_all(Pair_Map(c,d,t_1),i_ph) = D_Y_all(Pair_Map(c,d,t_1),i_ph) - &
                                Phase_hp * (real(kronecker_delta(c,a_h),dp) - rho_t(a_h,c,T_ph)) * rho_t(d,a_p,T_ph)
                    end do
                end do

                ! 2-body part (RPA 2-body correlation function) ... as in RPA_rM ...
                if (TBCF) then
                    Phase_X = phase(1 + (j_p + j_h)/2)
                    Phase_Y = phase(J + (j_p + j_h)/2)
                    do ia = 1, a_max
                        a_a = Orb(ia)%a
                        l_a = Orb(ia)%l
                        j_a = Orb(ia)%j
                        P_ha = parity_of(l_h + l_a)
                        do ib = 1, a_max
                            a_b = Orb(ib)%a
                            l_b = Orb(ib)%l
                            j_b = Orb(ib)%j
                            if ((.not. triangle(j_a,j_b,2*J)) .or. parity_of(l_p + l_b) /= P_ha) cycle
                            K_1 = 0.0_dp
                            K_2 = 0.0_dp
                            do I = max(abs(j_h - j_a),abs(j_p - j_b))/2, min(j_h + j_a,j_p + j_b)/2
                                F = phase(J + I + (j_p + j_a)/2) * hat_of(2*I) * f6j(j_p,j_h,2*J,j_a,j_b,2*I)
                                if (T_ph == -1) then
                                    K_1 = K_1 + O2b_pp(a_h,a_a,a_p,a_b,I,P_ha,Sigma,Orb,Orb_NN) * F
                                    K_2 = K_2 + O2b_pn(a_h,a_a,a_p,a_b,I,P_ha,Sigma,Orb_NN) * F
                                else
                                    K_1 = K_1 + O2b_nn(a_h,a_a,a_p,a_b,I,P_ha,Sigma,Orb,Orb_NN) * F
                                    K_2 = K_2 + O2b_pn(a_a,a_h,a_b,a_p,I,P_ha,Sigma,Orb_NN) * phase((j_p + j_a + j_h + j_b)/2) * F
                                end if
                            end do
                            Phase_ab = phase((j_a - j_b)/2)
                            D_X_all(Pair_Map(a_a,a_b,t_1),i_ph) = D_X_all(Pair_Map(a_a,a_b,t_1),i_ph) + Phase_X * K_1
                            D_X_all(Pair_Map(a_a,a_b,t_2),i_ph) = D_X_all(Pair_Map(a_a,a_b,t_2),i_ph) + Phase_X * K_2
                            D_Y_all(Pair_Map(a_b,a_a,t_1),i_ph) = D_Y_all(Pair_Map(a_b,a_a,t_1),i_ph) + Phase_Y * Phase_ab * K_1
                            D_Y_all(Pair_Map(a_b,a_a,t_2),i_ph) = D_Y_all(Pair_Map(a_b,a_a,t_2),i_ph) + Phase_Y * Phase_ab * K_2
                        end do
                    end do
                end if
            end if
        end do
        !$omp end parallel do

        ! Keep only the pairs with a non-vanishing kernel ...
        allocate(Keep(N_all))
        do k = 1, N_all
            Keep(k) = any(D_X_all(k,:) /= 0.0_dp) .or. any(D_Y_all(k,:) /= 0.0_dp)
        end do
        Pairs = pack(Pairs_all(1:N_all),Keep)
        allocate(D_X(size(Pairs),N_ph), D_Y(size(Pairs),N_ph))
        i_ph = 0
        do k = 1, N_all
            if (.not. Keep(k)) cycle
            i_ph = i_ph + 1
            D_X(i_ph,:) = D_X_all(k,:)
            D_Y(i_ph,:) = D_Y_all(k,:)
        end do

    contains

        pure real(dp) function rho_t(a,b,T_)
            integer, intent(in) :: a, b, T_
            if (T_ == -1) then
                rho_t = Rho%p(a,b)
            else
                rho_t = Rho%n(a,b)
            end if
        end function rho_t

    end subroutine RPA_transition_kernel

    ! Kernel, metric & transition density matrix Gamma of the selected levels ... TDA: TDA_occupation_factors & D X ,
        ! RPA: RPA_transition_kernel & RPA_transition_density_matrix ...
    subroutine ext_gamma(Params,Dt,J,P,nu_list,G)
        type(params_t), intent(in) :: Params
        type(ext_data_t), intent(in) :: Dt
        integer, intent(in) :: J, P, nu_list(:)
        type(ext_gamma_t), intent(out) :: G
        integer :: N_ph, k
        real(dp), allocatable :: Xs(:,:)
        complex(dp), allocatable :: Xc(:,:), Yc(:,:)

        G%J = J
        G%P = P
        N_ph = Dt%Orb_1p1h%N(P,J)
        if (Dt%Kind == 'TDA') then
            call TDA_transition_kernel(Params,J,P,Dt%Orb,Dt%Orb_1p1h,Dt%Orb_NN,Dt%Rho,Dt%Sigma,G%Pairs,G%D_X)
            allocate(G%D_Y(size(G%Pairs),N_ph))
            G%D_Y = 0.0_dp
            allocate(Xs(N_ph,size(nu_list)))
            do k = 1, size(nu_list)
                Xs(:,k) = Dt%X_T(P,J)%v(:,nu_list(k))
            end do
            G%Gamma = cmplx(matmul(G%D_X,Xs),0.0_dp,kind=dp)
        else
            call RPA_transition_kernel(Params,J,P,Dt%Orb,Dt%Orb_1p1h,Dt%Orb_NN,Dt%Rho,Dt%Sigma,G%Pairs,G%D_X,G%D_Y)
            allocate(Xc(N_ph,size(nu_list)), Yc(N_ph,size(nu_list)))
            do k = 1, size(nu_list)
                Xc(:,k) = Dt%X(P,J)%v(:,nu_list(k))
                Yc(:,k) = Dt%Y(P,J)%v(:,nu_list(k))
            end do
            G%Gamma = matmul(cmplx(G%D_X,0.0_dp,kind=dp),Xc) + matmul(cmplx(G%D_Y,0.0_dp,kind=dp),Yc)
        end if
        G%f_N = ph_metric(Dt%Orb,Dt%Orb_1p1h,P,J)
    end subroutine ext_gamma

    ! Reduced matrix elements <nu||M||0>/hat(J) of the p, n, is & iv EJ operators ... counterpart of RPA_external_multipoles
        ! (& the identical inline blocks of TDA_external.jl) ...
    subroutine ext_multipoles(Params,Dt,G)
        type(params_t), intent(in) :: Params
        type(ext_data_t), intent(in) :: Dt
        type(ext_gamma_t), intent(inout) :: G
        type(o1b_t) :: Op
        real(dp) :: A, Z, e_p, e_n
        complex(dp) :: Amp, MME, pS, nS, isS, ivS
        integer :: N_phonon, nu, ab, a_a, a_b

        A = real(Params%Calc%A,dp)
        Z = real(Params%Calc%Z,dp)
        e_p = (A - Z) / A
        e_n = Z / A
        N_phonon = size(G%Gamma,2)
        allocate(G%pM(N_phonon), G%nM(N_phonon), G%isM(N_phonon), G%ivM(N_phonon))
        G%pM = (0.0_dp,0.0_dp)
        G%nM = (0.0_dp,0.0_dp)
        G%isM = (0.0_dp,0.0_dp)
        G%ivM = (0.0_dp,0.0_dp)
        if (.not. is_EJ(G%J,G%P)) return
        select case (G%J)
        case (0)
            Op = Dt%TrOp%E0
        case (1)
            Op = Dt%TrOp%E1
        case (2)
            Op = Dt%TrOp%E2
        case default
            Op = Dt%TrOp%E3
        end select

        do nu = 1, N_phonon
            pS = (0.0_dp,0.0_dp)
            nS = (0.0_dp,0.0_dp)
            isS = (0.0_dp,0.0_dp)
            ivS = (0.0_dp,0.0_dp)
            do ab = 1, size(G%Pairs)
                a_a = G%Pairs(ab)%a
                a_b = G%Pairs(ab)%b
                Amp = G%Gamma(ab,nu) / sqrt(real(2*G%J + 1,dp))
                if (G%Pairs(ab)%T == -1) then
                    MME = Op%p(a_a,a_b) * Amp
                    pS = pS + MME
                    if (G%J == 1) then
                        isS = isS + 0.5_dp * (Dt%TrOp%E1_C%p(a_a,a_b) + Dt%TrOp%E1_sC%p(a_a,a_b)) * Amp
                        ivS = ivS + e_p * MME
                    else
                        isS = isS + 0.5_dp * MME
                        ivS = ivS + 0.5_dp * MME
                    end if
                else if (G%Pairs(ab)%T == 1) then
                    MME = Op%n(a_a,a_b) * Amp
                    nS = nS + MME
                    if (G%J == 1) then
                        isS = isS + 0.5_dp * (Dt%TrOp%E1_C%n(a_a,a_b) + Dt%TrOp%E1_sC%n(a_a,a_b)) * Amp
                        ivS = ivS - e_n * MME
                    else
                        isS = isS + 0.5_dp * MME
                        ivS = ivS - 0.5_dp * MME
                    end if
                end if
            end do
            G%pM(nu) = pS
            G%nM(nu) = nS
            G%isM(nu) = isS
            G%ivM(nu) = ivS
        end do
    end subroutine ext_multipoles

    ! Gauge of the exported densities ... conj(pM)/|pM| (1 if pM = 0 or for a non-EJ multipole) ... for the real TDA
        ! amplitudes this is dev's pM/|pM| ...
    function ext_gauge(G) result(Gauge)
        type(ext_gamma_t), intent(in) :: G
        complex(dp), allocatable :: Gauge(:)
        integer :: nu
        allocate(Gauge(size(G%pM)))
        Gauge = (1.0_dp,0.0_dp)
        if (.not. is_EJ(G%J,G%P)) return
        do nu = 1, size(G%pM)
            if (abs(G%pM(nu)) > 0.0_dp) Gauge(nu) = conjg(G%pM(nu)) / abs(G%pM(nu))
        end do
    end function ext_gauge

    ! Radial reference-basis orbitals R_a(r) = sum_k C[k,a] R_k^LHO(r) ... counterpart of RPA_external_orbitals_grid ...
    subroutine ext_orbitals_grid(Orb,C,r_grid,nu_proton,nu_neutron,pR,nR)
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: C
        real(dp), intent(in) :: r_grid(:), nu_proton, nu_neutron
        real(dp), allocatable, intent(out) :: pR(:,:), nR(:,:)
        integer :: a_max, N_grid, i, a, k
        real(dp) :: pS, nS
        a_max = size(Orb)
        N_grid = size(r_grid)
        allocate(pR(a_max,N_grid), nR(a_max,N_grid))
        !$omp parallel do schedule(static) default(shared) private(a,k,pS,nS)
        do i = 1, N_grid
            do a = 1, a_max
                pS = 0.0_dp
                nS = 0.0_dp
                do k = 1, a_max
                    if (Orb(a)%l /= Orb(k)%l .or. Orb(a)%j /= Orb(k)%j) cycle
                    pS = pS + C%p(k,a) * Psi_rad_LHO(r_grid(i),Orb(k)%n,Orb(k)%l,nu_proton)
                    nS = nS + C%n(k,a) * Psi_rad_LHO(r_grid(i),Orb(k)%n,Orb(k)%l,nu_neutron)
                end do
                pR(a,i) = pS
                nR(a,i) = nS
            end do
        end do
        !$omp end parallel do
    end subroutine ext_orbitals_grid

    ! Julia collect(range(a, stop = b, length = N)) ...
    function lin_range(a,b,N) result(g)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: N
        real(dp), allocatable :: g(:)
        integer :: i
        allocate(g(N))
        do i = 1, N
            g(i) = (a * real(N - i,dp) + b * real(i - 1,dp)) / real(N - 1,dp)
        end do
    end function lin_range

    ! <a||Y_J||b> ... zero for odd l_a + l_b + J ...
    real(dp) function YJ_reduced(Orb,a,b,J)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: a, b, J
        YJ_reduced = 0.0_dp
        if (mod(Orb(a)%l + Orb(b)%l + J,2) == 1) return
        YJ_reduced = phase(J + (Orb(a)%j - 1)/2) * sqrt(real((Orb(a)%j + 1) * (Orb(b)%j + 1),dp)) * &
                     fCG(Orb(a)%j,Orb(b)%j,2*J,1,-1,0) / sqrt(4.0_dp * pi)
    end function YJ_reduced

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Transition densities & transition operator densities ...
    ! ---------------------------------------------------------------------------------------------------------------------

    ! Counterpart of {TDA,RPA}_transition_densities (Operator = .false.) & {TDA,RPA}_transition_operator_densities ...
    subroutine ext_densities(Params,Kind,JP,nu_list,File_Name,Operator)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind, JP, File_Name
        integer, intent(in) :: nu_list(:)
        logical, intent(in) :: Operator
        type(ext_data_t) :: Dt
        type(ext_gamma_t) :: G
        real(dp), allocatable :: r_grid(:), pR(:,:), nR(:,:), pRho(:,:), nRho(:,:), Amp_ab(:)
        complex(dp), allocatable :: Gauge(:)
        real(dp) :: A, hw, nu_proton, nu_neutron, r_min, r_max, Amp, pS, nS, Rho_ab
        integer :: J, P, Power, N_grid, N_phonon, i, nu, ab
        character(len=8) :: H(4)
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        A = real(Params%Calc%A,dp)
        hw = Params%Calc%hw
        nu_proton = 0.5_dp * m_p * hw / hc**2
        nu_neutron = 0.5_dp * m_n * hw / hc**2
        r_min = 0.0_dp + 1.0e-8_dp
        r_max = 2.5_dp * 1.2_dp * A**(1.0_dp/3.0_dp)
        N_grid = 2**13 + 1

        call ext_import(Params,Kind,Dt)
        call ext_JP(JP,Dt,J,P)
        Power = 0
        if (Operator) then
            if (.not. is_EJ(J,P)) &
                call fatal('The given values of J & P are not supported for radial transition operator density calculations ...')
            Power = merge(2, J, J == 0)
        end if
        call check_nu(Dt,J,P,nu_list)
        N_phonon = size(nu_list)

        call timer_start(t0)
        call ext_gamma(Params,Dt,J,P,nu_list,G)
        call ext_multipoles(Params,Dt,G)

        r_grid = lin_range(r_min,r_max,N_grid)
        call ext_orbitals_grid(Dt%Orb,Dt%C,r_grid,nu_proton,nu_neutron,pR,nR)

        ! r^(2+Power)/(2J+1) sum_ab Re(g Gamma_ab) <a||Y_J||b> R_a R_b ...
        allocate(Amp_ab(size(G%Pairs)))
        do ab = 1, size(G%Pairs)
            Amp_ab(ab) = YJ_reduced(Dt%Orb,G%Pairs(ab)%a,G%Pairs(ab)%b,J)
        end do
        Gauge = ext_gauge(G)
        allocate(pRho(N_phonon,N_grid), nRho(N_phonon,N_grid))
        !$omp parallel do schedule(static) default(shared) private(Amp,nu,pS,nS,ab,Rho_ab)
        do i = 1, N_grid
            Amp = r_grid(i)**(2 + Power) / real(2*J + 1,dp)
            do nu = 1, N_phonon
                pS = 0.0_dp
                nS = 0.0_dp
                do ab = 1, size(G%Pairs)
                    if (Amp_ab(ab) == 0.0_dp) cycle
                    Rho_ab = real(Gauge(nu) * G%Gamma(ab,nu),dp) * Amp_ab(ab)
                    if (G%Pairs(ab)%T == -1) then
                        pS = pS + Rho_ab * pR(G%Pairs(ab)%a,i) * pR(G%Pairs(ab)%b,i)
                    else if (G%Pairs(ab)%T == 1) then
                        nS = nS + Rho_ab * nR(G%Pairs(ab)%a,i) * nR(G%Pairs(ab)%b,i)
                    end if
                end do
                pRho(nu,i) = Amp * pS
                nRho(nu,i) = Amp * nS
            end do
        end do
        !$omp end parallel do

        if (Operator) then
            H = [character(len=8) :: 'pM', 'nM', 'isM', 'ivM']
        else
            H = [character(len=8) :: 'rho_p', 'rho_n', 'rho_is', 'rho_iv']
        end if
        call make_dir(Dt%IO // '/densities')
        call export_densities()
        call timer_report(t0)

    contains

        subroutine export_densities()
            integer :: u, nu_, i_
            real(dp) :: pS_, nS_, isS_, ivS_
            do nu_ = 1, N_phonon
                open(newunit=u, file=Dt%IO // '/densities/' // Kind // '_' // File_Name // '_JP' // trim(adjustl(JP)) // '_nu' // &
                     istr(nu_list(nu_)) // '_E' // level_E(Dt,J,P,nu_list(nu_)) // '.dat', status='replace', action='write')
                write(u,'(a)') row_s(['r       ', H])
                do i_ = 1, N_grid
                    write(u,'(a)') row_f([r_grid(i_), pRho(nu_,i_), nRho(nu_,i_), 0.5_dp * (pRho(nu_,i_) + nRho(nu_,i_)), &
                                          0.5_dp * (pRho(nu_,i_) - nRho(nu_,i_))])
                end do
                close(u)
            end do
            write(*,'(a)') '"Export of individual phonon radial transition densities successfuly completed ..."'

            if (.not. is_EJ(J,P)) return
            call phonon_weights(G%pM)
            call phonon_weights(G%nM)
            call phonon_weights(G%isM)
            call phonon_weights(G%ivM)
            open(newunit=u, file=Dt%IO // '/densities/' // Kind // '_' // File_Name // '_JP' // trim(adjustl(JP)) // '_Averaged.dat', &
                 status='replace', action='write')
            write(u,'(a)') row_s([character(len=8) :: 'r', 'rho_p', 'rho_n', 'rho_is', 'rho_iv'])
            do i_ = 1, N_grid
                pS_ = 0.0_dp
                nS_ = 0.0_dp
                isS_ = 0.0_dp
                ivS_ = 0.0_dp
                do nu_ = 1, N_phonon
                    pS_ = pS_ + real(G%pM(nu_),dp) * pRho(nu_,i_)
                    nS_ = nS_ + real(G%nM(nu_),dp) * nRho(nu_,i_)
                    isS_ = isS_ + 0.5_dp * real(G%isM(nu_),dp) * (pRho(nu_,i_) + nRho(nu_,i_))
                    ivS_ = ivS_ + 0.5_dp * real(G%ivM(nu_),dp) * (pRho(nu_,i_) - nRho(nu_,i_))
                end do
                write(u,'(a)') row_f([r_grid(i_), pS_, nS_, isS_, ivS_])
            end do
            close(u)
            write(*,'(a)') '"Export of averaged phonon radial transitions densities successfuly completed ..."'
        end subroutine export_densities

    end subroutine ext_densities

    ! "%-20s %-20s ..." ...
    function row_s(H) result(s)
        character(len=*), intent(in) :: H(:)
        character(len=:), allocatable :: s
        integer :: k
        s = cfmt_s(trim(H(1)),20,.true.)
        do k = 2, size(H)
            s = s // ' ' // cfmt_s(trim(H(k)),20,.true.)
        end do
    end function row_s

    ! "%-20.8f %-20.8f ..." ...
    function row_f(v) result(s)
        real(dp), intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: k
        s = cfmt_f(v(1),20,8,.true.)
        do k = 2, size(v)
            s = s // ' ' // cfmt_f(v(k),20,8,.true.)
        end do
    end function row_f

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Analysis of the phonon levels ...
    ! ---------------------------------------------------------------------------------------------------------------------

    ! Counterpart of TDA_transition_analysis (TDA_external.jl:759-907) & RPA_transition_analysis (RPA_external.jl:647-741) ...
    subroutine ext_analysis(Params,Kind,JP,nu_list,File_Name,Cut)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind, JP, File_Name
        integer, intent(in) :: nu_list(:)
        real(dp), intent(in) :: Cut
        type(ext_data_t) :: Dt
        type(ext_gamma_t) :: G
        type(o1b_t) :: Op
        real(dp), allocatable :: pW_X(:), pW_Y(:), nW_X(:), nW_Y(:)
        real(dp) :: M_1, M_2, Amp_M, x_T, pM_ph, nM_ph
        complex(dp) :: x_R, y_R, Amp, X_ph, Y_ph
        integer :: J, P, N_ph, nu_ind, nu, ph, ab, u, a_p, a_h, t_ph
        logical :: Have_Op

        call wigner_init_float(75,'Jmax',9)
        call ext_import(Params,Kind,Dt)
        call ext_JP(JP,Dt,J,P)
        call check_nu(Dt,J,P,nu_list)
        call make_dir(Dt%IO // '/levels')

        N_ph = Dt%Orb_1p1h%N(P,J)
        call ext_gamma(Params,Dt,J,P,nu_list,G)

        Have_Op = is_EJ(J,P)
        if (Have_Op) then
            select case (J)
            case (0)
                Op = Dt%TrOp%E0
            case (1)
                Op = Dt%TrOp%E1
            case (2)
                Op = Dt%TrOp%E2
            case default
                Op = Dt%TrOp%E3
            end select
        end if

        if (Kind == 'RPA') then
            ! W^t_X[ph] = sum_(a,b,t) M^t_ab D_X[(a,b,t),ph] (& W_Y) ...
            allocate(pW_X(N_ph), pW_Y(N_ph), nW_X(N_ph), nW_Y(N_ph))
            pW_X = 0.0_dp
            pW_Y = 0.0_dp
            nW_X = 0.0_dp
            nW_Y = 0.0_dp
            if (Have_Op) then
                do ab = 1, size(G%Pairs)
                    do ph = 1, N_ph
                        if (G%Pairs(ab)%T == -1) then
                            pW_X(ph) = pW_X(ph) + Op%p(G%Pairs(ab)%a,G%Pairs(ab)%b) * G%D_X(ab,ph)
                            pW_Y(ph) = pW_Y(ph) + Op%p(G%Pairs(ab)%a,G%Pairs(ab)%b) * G%D_Y(ab,ph)
                        else
                            nW_X(ph) = nW_X(ph) + Op%n(G%Pairs(ab)%a,G%Pairs(ab)%b) * G%D_X(ab,ph)
                            nW_Y(ph) = nW_Y(ph) + Op%n(G%Pairs(ab)%a,G%Pairs(ab)%b) * G%D_Y(ab,ph)
                        end if
                    end do
                end do
            end if
        end if

        do nu_ind = 1, size(nu_list)
            nu = nu_list(nu_ind)
            open(newunit=u, file=Dt%IO // '/levels/' // Kind // '_' // File_Name // '_JP' // trim(adjustl(JP)) // '_nu' // &
                 istr(nu) // '_E' // level_E(Dt,J,P,nu) // '.dat', status='replace', action='write')
            if (Kind == 'TDA') then
                write(u,'(a)') i8s('a_p') // ' ' // i8s('a_h') // ' ' // i8s('T_ph') // ' ' // cfmt_s('X_ph',20,.true.) // ' ' // &
                               cfmt_s('M_ph',20,.true.) // ' ' // cfmt_s('M_ph(-t)',20,.true.)
            else
                write(u,'(a)') i8s('a_p') // ' ' // i8s('a_h') // ' ' // i8s('T_ph') // ' ' // cfmt_s('Re X_ph',20,.true.) // ' ' // &
                               cfmt_s('Re Y_ph',20,.true.) // ' ' // cfmt_s('Re M_ph',20,.true.) // ' ' // &
                               cfmt_s('Re M_ph(-t)',20,.true.)
            end if
            do ph = 1, N_ph
                a_p = Dt%Orb_1p1h%Ind(P,J)%v(ph)%a
                a_h = Dt%Orb_1p1h%Ind(P,J)%v(ph)%b
                t_ph = Dt%Orb_1p1h%Ind(P,J)%v(ph)%T
                if (Kind == 'TDA') then
                    x_T = sqrt(G%f_N(ph)) * Dt%X_T(P,J)%v(ph,nu)
                    if (x_T**2 <= Cut) cycle
                    M_1 = 0.0_dp
                    M_2 = 0.0_dp
                    if (Have_Op) then
                        do ab = 1, size(G%Pairs)
                            if (G%D_X(ab,ph) == 0.0_dp) cycle
                            if (G%Pairs(ab)%T == -1) then
                                Amp_M = G%D_X(ab,ph) * Dt%X_T(P,J)%v(ph,nu) * Op%p(G%Pairs(ab)%a,G%Pairs(ab)%b)
                            else
                                Amp_M = G%D_X(ab,ph) * Dt%X_T(P,J)%v(ph,nu) * Op%n(G%Pairs(ab)%a,G%Pairs(ab)%b)
                            end if
                            if (G%Pairs(ab)%T == t_ph) then
                                M_1 = M_1 + Amp_M
                            else
                                M_2 = M_2 + Amp_M
                            end if
                        end do
                    end if
                    write(u,'(a)') i8d(a_p) // ' ' // i8d(a_h) // ' ' // i8d(t_ph) // ' ' // row_f([x_T, M_1, M_2])
                else
                    x_R = sqrt(G%f_N(ph)) * Dt%X(P,J)%v(ph,nu)
                    y_R = sqrt(G%f_N(ph)) * Dt%Y(P,J)%v(ph,nu)
                    Amp = phase(J) * x_R + y_R
                    if (real(Amp,dp)**2 + aimag(Amp)**2 <= Cut) cycle
                    X_ph = Dt%X(P,J)%v(ph,nu)
                    Y_ph = Dt%Y(P,J)%v(ph,nu)
                    pM_ph = real(pW_X(ph) * X_ph + pW_Y(ph) * Y_ph,dp)
                    nM_ph = real(nW_X(ph) * X_ph + nW_Y(ph) * Y_ph,dp)
                    if (t_ph == -1) then
                        M_1 = pM_ph
                        M_2 = nM_ph
                    else
                        M_1 = nM_ph
                        M_2 = pM_ph
                    end if
                    write(u,'(a)') i8d(a_p) // ' ' // i8d(a_h) // ' ' // i8d(t_ph) // ' ' // &
                                   row_f([real(x_R,dp), real(y_R,dp), M_1, M_2])
                end if
            end do
            close(u)
        end do
        write(*,'(a)') '"Export of information on individual phonon levels successfuly completed ..."'

    contains

        function i8s(S) result(r)
            character(len=*), intent(in) :: S
            character(len=:), allocatable :: r
            r = cfmt_s(S,8,.true.)
        end function i8s

        function i8d(k) result(r)
            integer, intent(in) :: k
            character(len=:), allocatable :: r
            r = cfmt_s(istr(k),8,.true.)
        end function i8d

    end subroutine ext_analysis

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Transition current maps ...
    ! ---------------------------------------------------------------------------------------------------------------------

    ! Counterpart of TDA_transition_currents (TDA_external.jl:909-1462) & RPA_transition_currents (RPA_external.jl:743-1052) ...
        ! TDA vs RPA (as in dev, UPSTREAM_BUGS): pair maps with (-1)^((j_b - m_b)/2) resp. (-1)^((j_b + m_b)/2) (and the
        ! |cg| < 1e-14 skip of RPA), level amplitudes g (-1)^J Gamma resp. -Re(g Gamma)/hat(J) ...
    subroutine ext_currents(Params,Kind,JP,nu_list,File_Name,M,xN_grid,zN_grid)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind, JP, File_Name
        integer, intent(in) :: nu_list(:), M, xN_grid, zN_grid
        type(ext_data_t) :: Dt
        type(ext_gamma_t) :: G
        integer, parameter :: rN_grid = 2**13 + 1, tN_grid = 256
        real(dp), allocatable :: r_grid(:), theta_grid(:), pR(:,:), nR(:,:), pdR(:,:), ndR(:,:), x_map(:,:)
        real(dp), allocatable :: Jc_ab(:,:,:), Js_ab(:,:,:), pJc(:,:,:), nJc(:,:,:), pJs(:,:,:), nJs(:,:,:), Lv(:,:)
        integer, allocatable :: i_map(:,:)
        complex(dp), allocatable :: Gauge(:)
        real(dp) :: A, hw, nu_proton, nu_neutron, r_min, r_max, x_min, x_max, z_min, z_max, dx, dz, dr, dtheta
        real(dp) :: x, z, r, theta, Norm_J
        integer :: J, P, M_proj, N_phonon, N_ab, N_map, ia_, ind_x, ind_z, n_r, m_t, i_N, k, nu, ab
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        A = real(Params%Calc%A,dp)
        hw = Params%Calc%hw
        nu_proton = 0.5_dp * m_p * hw / hc**2
        nu_neutron = 0.5_dp * m_n * hw / hc**2
        r_min = 0.01_dp
        r_max = 3.0_dp * A**(1.0_dp/3.0_dp)
        M_proj = merge(0, M, M == -1)
        if (xN_grid < 1 .or. zN_grid < 1) call fatal('xN_grid & zN_grid must be positive ...')

        call ext_import(Params,Kind,Dt)
        call ext_JP(JP,Dt,J,P)
        call check_nu(Dt,J,P,nu_list)
        N_phonon = size(nu_list)
        call ext_gamma(Params,Dt,J,P,nu_list,G)
        call ext_multipoles(Params,Dt,G)
        N_ab = size(G%Pairs)

        r_grid = lin_range(r_min,r_max,rN_grid)
        theta_grid = lin_range(0.0_dp,2.0_dp*pi,tN_grid)
        call ext_orbitals_grid(Dt%Orb,Dt%C,r_grid,nu_proton,nu_neutron,pR,nR)
        allocate(pdR(size(Dt%Orb),rN_grid), ndR(size(Dt%Orb),rN_grid))
        !$omp parallel do schedule(static) default(shared)
        do ia_ = 1, size(Dt%Orb)
            pdR(ia_,:) = differentiate_Oh6(r_grid,pR(ia_,:))
            ndR(ia_,:) = differentiate_Oh6(r_grid,nR(ia_,:))
        end do
        !$omp end parallel do

        ! Cartesian map points snapped to the (r, theta) grid ...
        x_min = -0.8_dp * r_max
        x_max = 0.8_dp * r_max
        z_min = -0.8_dp * r_max
        z_max = 0.8_dp * r_max
        dx = (x_max - x_min) / real(xN_grid,dp)
        dz = (z_max - z_min) / real(zN_grid,dp)
        N_map = xN_grid * zN_grid
        allocate(i_map(2,N_map), x_map(2,N_map))
        dr = (r_max - r_min) / real(rN_grid - 1,dp)
        dtheta = 2.0_dp * pi / real(tN_grid - 1,dp)
        call map_grid()

        ! Orbital pair current maps (each map point owned by one thread) ...
        allocate(Jc_ab(3,N_ab,N_map), Js_ab(3,N_ab,N_map))
        call timer_start(t0)
        !$omp parallel do schedule(dynamic) default(shared)
        do i_N = 1, N_map
            call pair_maps(i_N)
        end do
        !$omp end parallel do
        call timer_report(t0)

        ! Level current maps ...
        Gauge = ext_gauge(G)
        Norm_J = 1.0_dp / sqrt(real(2*J + 1,dp))
        allocate(pJc(3,N_phonon,N_map), nJc(3,N_phonon,N_map), pJs(3,N_phonon,N_map), nJs(3,N_phonon,N_map))
        allocate(Lv(N_ab,N_phonon))
        do nu = 1, N_phonon
            do ab = 1, N_ab
                if (Kind == 'TDA') then
                    Lv(ab,nu) = real(Gauge(nu),dp) * phase(J) * real(G%Gamma(ab,nu),dp)
                else
                    Lv(ab,nu) = -Norm_J * real(Gauge(nu) * G%Gamma(ab,nu),dp)
                end if
            end do
        end do
        !$omp parallel do schedule(static) default(shared) private(nu,ab,k)
        do i_N = 1, N_map
            do nu = 1, N_phonon
                pJc(:,nu,i_N) = 0.0_dp
                nJc(:,nu,i_N) = 0.0_dp
                pJs(:,nu,i_N) = 0.0_dp
                nJs(:,nu,i_N) = 0.0_dp
                do ab = 1, N_ab
                    if (G%Pairs(ab)%T == -1) then
                        do k = 1, 3
                            pJc(k,nu,i_N) = pJc(k,nu,i_N) + Lv(ab,nu) * Jc_ab(k,ab,i_N)
                            pJs(k,nu,i_N) = pJs(k,nu,i_N) + Lv(ab,nu) * Js_ab(k,ab,i_N)
                        end do
                    else if (G%Pairs(ab)%T == 1) then
                        do k = 1, 3
                            nJc(k,nu,i_N) = nJc(k,nu,i_N) + Lv(ab,nu) * Jc_ab(k,ab,i_N)
                            nJs(k,nu,i_N) = nJs(k,nu,i_N) + Lv(ab,nu) * Js_ab(k,ab,i_N)
                        end do
                    end if
                end do
            end do
        end do
        !$omp end parallel do

        ! Units [c] ...
        pJc = pJc * 0.5_dp * hc / m_p
        nJc = nJc * 0.5_dp * hc / m_n
        pJs = pJs * hc / m_p
        nJs = nJs * hc / m_n

        call make_dir(Dt%IO // '/densities')
        do nu = 1, N_phonon
            call write_level('Convective',pJc,nJc,nu)
            call write_level('Spin',pJs,nJs,nu)
        end do
        write(*,'(a)') '"Export of individual phonon level transition current maps successfuly completed ..."'

        if (.not. is_EJ(J,P)) return
        call phonon_weights(G%pM)
        call phonon_weights(G%nM)
        call phonon_weights(G%isM)
        call phonon_weights(G%ivM)
        call write_averaged('Convective',pJc,nJc)
        call write_averaged('Spin',pJs,nJs)
        write(*,'(a)') '"Export of averaged phonon level transition current maps successfuly completed ..."'

    contains

        ! Map points snapped to the nearest (r, theta) node ... the products are rounded separately (volatile, no FMA
            ! contraction) & the loop is not vectorized (scalar atan2) so that the ties of round(theta/dtheta) (e.g. x = 0,
            ! theta = pi, pi/dtheta = 127.5) resolve exactly as in Julia ...
        subroutine map_grid()
            real(dp), volatile :: t_x, t_z, x2, z2, q_r, q_t
            integer :: k_
            k_ = 0
            !DIR$ NOVECTOR
            do ind_x = 1, xN_grid
                t_x = real(ind_x,dp) * dx
                x = x_min + t_x
                !DIR$ NOVECTOR
                do ind_z = 1, zN_grid
                    t_z = real(ind_z,dp) * dz
                    z = z_min + t_z
                    x2 = x * x
                    z2 = z * z
                    r = sqrt(x2 + z2)
                    theta = atan2(x,z)
                    if (theta < 0.0_dp) theta = theta + 2.0_dp*pi
                    q_r = (r - r_min) / dr
                    q_t = (theta - 0.0_dp) / dtheta
                    n_r = min(max(round_even(q_r) + 1, 1), rN_grid)
                    m_t = min(max(round_even(q_t) + 1, 1), tN_grid)
                    k_ = k_ + 1
                    i_map(1,k_) = n_r
                    i_map(2,k_) = m_t
                    x_map(1,k_) = r_grid(n_r)
                    x_map(2,k_) = theta_grid(m_t)
                end do
            end do
        end subroutine map_grid

        ! J_ab = sum_(m_a, m_b) cg(j_a j_b J; m_a m_b M) phase * J1b(a m_a, b -m_b) at map point N_ ...
        subroutine pair_maps(N_)
            integer, intent(in) :: N_
            integer :: nr_, ab_, a_a, a_b, l_a, j_a, l_b, j_b, m_a, m_b
            real(dp) :: r_, th_, gs, R_a, R_b, dR_a, dR_b, Amp, Jc(3), Js(3)
            complex(dp) :: jx, jy, jz
            nr_ = i_map(1,N_)
            r_ = r_grid(nr_)
            th_ = theta_grid(i_map(2,N_))
            do ab_ = 1, N_ab
                a_a = G%Pairs(ab_)%a
                a_b = G%Pairs(ab_)%b
                l_a = Dt%Orb(a_a)%l
                j_a = Dt%Orb(a_a)%j
                l_b = Dt%Orb(a_b)%l
                j_b = Dt%Orb(a_b)%j
                if (G%Pairs(ab_)%T == -1) then
                    R_a = pR(a_a,nr_)
                    R_b = pR(a_b,nr_)
                    dR_a = pdR(a_a,nr_)
                    dR_b = pdR(a_b,nr_)
                    gs = 0.5_dp * g_p
                else
                    R_a = nR(a_a,nr_)
                    R_b = nR(a_b,nr_)
                    dR_a = ndR(a_a,nr_)
                    dR_b = ndR(a_b,nr_)
                    gs = 0.5_dp * g_n
                end if
                Jc = 0.0_dp
                Js = 0.0_dp
                do m_a = -j_a, j_a, 2
                    m_b = 2*M_proj - m_a
                    if (abs(m_b) > j_b) cycle
                    if (Kind == 'TDA') then
                        Amp = fCG(j_a,j_b,2*J,m_a,m_b,2*M_proj) * phase((j_b - m_b)/2)
                    else
                        Amp = fCG(j_a,j_b,2*J,m_a,m_b,2*M_proj) * phase((j_b + m_b)/2)
                        if (abs(Amp) < 1.0e-14_dp) cycle
                    end if
                    call J1b_convective(r_,th_,0.0_dp,l_a,j_a,m_a,R_a,dR_a,l_b,j_b,-m_b,R_b,dR_b,jx,jy,jz)
                    Jc(1) = Jc(1) + Amp * real(jx,dp)
                    Jc(2) = Jc(2) + Amp * real(jy,dp)
                    Jc(3) = Jc(3) + Amp * real(jz,dp)
                    call J1b_spin(r_,th_,0.0_dp,l_a,j_a,m_a,R_a,dR_a,l_b,j_b,-m_b,R_b,dR_b,jx,jy,jz)
                    Js(1) = Js(1) + Amp * real(gs * jx,dp)
                    Js(2) = Js(2) + Amp * real(gs * jy,dp)
                    Js(3) = Js(3) + Amp * real(gs * jz,dp)
                end do
                Jc_ab(:,ab_,N_) = Jc
                Js_ab(:,ab_,N_) = Js
            end do
        end subroutine pair_maps

        subroutine write_level(Label,pJ,nJ,nu_)
            character(len=*), intent(in) :: Label
            real(dp), intent(in) :: pJ(:,:,:), nJ(:,:,:)
            integer, intent(in) :: nu_
            integer :: u, N_
            real(dp) :: xx, zz
            open(newunit=u, file=Dt%IO // '/densities/' // Kind // '_' // Label // '_' // File_Name // '_JP' // &
                 trim(adjustl(JP)) // '_nu' // istr(nu_list(nu_)) // '_E' // level_E(Dt,J,P,nu_list(nu_)) // '.dat', &
                 status='replace', action='write')
            write(u,'(a)') row_s([character(len=8) :: 'x', 'z', 'pJ_x', 'pJ_z', 'nJ_x', 'nJ_z', 'isJ_x', 'isJ_z', 'ivJ_x', 'ivJ_z'])
            do N_ = 1, N_map
                xx = x_map(1,N_) * sin(x_map(2,N_))
                zz = x_map(1,N_) * cos(x_map(2,N_))
                write(u,'(a)') row_f([xx, zz, pJ(1,nu_,N_), pJ(3,nu_,N_), nJ(1,nu_,N_), nJ(3,nu_,N_), &
                                      0.5_dp * (pJ(1,nu_,N_) + nJ(1,nu_,N_)), 0.5_dp * (pJ(3,nu_,N_) + nJ(3,nu_,N_)), &
                                      0.5_dp * (pJ(1,nu_,N_) - nJ(1,nu_,N_)), 0.5_dp * (pJ(3,nu_,N_) - nJ(3,nu_,N_))])
            end do
            close(u)
        end subroutine write_level

        subroutine write_averaged(Label,pJ,nJ)
            character(len=*), intent(in) :: Label
            real(dp), intent(in) :: pJ(:,:,:), nJ(:,:,:)
            integer :: u, N_, nu_
            real(dp) :: xx, zz, S(12), pm, nm, ism, ivm, pJ_x, pJ_z, nJ_x, nJ_z
            open(newunit=u, file=Dt%IO // '/densities/' // Kind // '_' // Label // '_' // File_Name // '_JP' // &
                 trim(adjustl(JP)) // '_Averaged.dat', status='replace', action='write')
            write(u,'(a)') row_s([character(len=8) :: 'x', 'z', 'ppJ_x', 'ppJ_z', 'pnJ_x', 'pnJ_z', 'npJ_x', 'npJ_z', 'nnJ_x', &
                                  'nnJ_z', 'isJ_x', 'isJ_z', 'ivJ_x', 'ivJ_z'])
            do N_ = 1, N_map
                xx = x_map(1,N_) * sin(x_map(2,N_))
                zz = x_map(1,N_) * cos(x_map(2,N_))
                S = 0.0_dp
                do nu_ = 1, N_phonon
                    pJ_x = pJ(1,nu_,N_)
                    pJ_z = pJ(3,nu_,N_)
                    nJ_x = nJ(1,nu_,N_)
                    nJ_z = nJ(3,nu_,N_)
                    pm = real(G%pM(nu_),dp)
                    nm = real(G%nM(nu_),dp)
                    ism = real(G%isM(nu_),dp)
                    ivm = real(G%ivM(nu_),dp)
                    S(1) = S(1) + pm * pJ_x
                    S(2) = S(2) + pm * pJ_z
                    S(3) = S(3) + pm * nJ_x
                    S(4) = S(4) + pm * nJ_z
                    S(5) = S(5) + nm * pJ_x
                    S(6) = S(6) + nm * pJ_z
                    S(7) = S(7) + nm * nJ_x
                    S(8) = S(8) + nm * nJ_z
                    S(9) = S(9) + 0.5_dp * ism * (pJ_x + nJ_x)
                    S(10) = S(10) + 0.5_dp * ism * (pJ_z + nJ_z)
                    S(11) = S(11) + 0.5_dp * ivm * (pJ_x - nJ_x)
                    S(12) = S(12) + 0.5_dp * ivm * (pJ_z - nJ_z)
                end do
                write(u,'(a)') row_f([xx, zz, S])
            end do
            close(u)
        end subroutine write_averaged

    end subroutine ext_currents

    ! ---------------------------------------------------------------------------------------------------------------------
    ! Recalculation helpers ...
    ! ---------------------------------------------------------------------------------------------------------------------

    subroutine transitions_dirs(IO)
        character(len=*), intent(in) :: IO
        call make_dir(IO // '/transitions/E0')
        call make_dir(IO // '/transitions/E1')
        call make_dir(IO // '/transitions/E2')
        call make_dir(IO // '/transitions/E3')
    end subroutine transitions_dirs

    ! Cuts <Name>_summary.dat before the first Marker (dev: findfirst + write of Summary[1:first-1]) ...
    subroutine summary_truncate(File,Marker)
        character(len=*), intent(in) :: File, Marker
        character(len=:), allocatable :: Buf
        integer :: u, Sz, k
        logical :: Ex
        inquire(file=File, exist=Ex, size=Sz)
        if (.not. Ex) return
        allocate(character(len=Sz) :: Buf)
        if (Sz > 0) then
            open(newunit=u, file=File, access='stream', form='unformatted', status='old', action='read')
            read(u) Buf
            close(u)
        end if
        k = index(Buf,Marker)
        if (k == 0) return
        open(newunit=u, file=File, access='stream', form='unformatted', status='replace', action='write')
        if (k > 1) write(u) Buf(1:k-1)
        close(u)
    end subroutine summary_truncate

end module nj_ph_external
