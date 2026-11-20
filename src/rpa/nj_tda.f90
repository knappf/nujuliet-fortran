! Tamm-Dancoff approximation on a general (HF / correlated) reference (NuJuliet dev) ... counterpart of
    ! src/TDA/TDA_diagonalize.jl, TDA_transitions.jl (TDA_rM, TDA_transition_kernel, TDA_pair_operators) & TDA_export.jl ...
    ! TDA_allocate / TDA_allocate_TBCF / TDA_spurious / TDA_rB are the shared routines of nj_ph_common ...
module nj_tda
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg, only : eig_sym, eig_sym_gen, mm, mtm
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_tr1b, only : tr1b_t
    use nj_ph_common
    implicit none
    private

    public :: TDA_diagonalize, TDA_rM, TDA_transition_kernel, TDA_export, ph_rM_assemble

contains

    ! Solves A X = E N X per channel ... counterpart of TDA_diagonalize (TDA_diagonalize.jl) ... non-orthogonalized
        ! channels: H = N^-1/2 A N^-1/2 (dsyevr), X = N^-1/2 R; spurious 1- channel (Ortho, N_ph >= 2): generalized problem in
        ! the complement of Spur (dsygvd), the spurious level E = <Spur|A|Spur>/<Spur|N|Spur> stored as the 1st solution ...
    subroutine TDA_diagonalize(Params,Orb,Orb_1p1h,A,N,TrOp,E,X)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rmat_t), intent(in) :: A(:,0:), N(:,0:)
        type(tr1b_t), intent(in) :: TrOp
        type(rvec_t), allocatable, intent(out) :: E(:,:)
        type(rmat_t), allocatable, intent(out) :: X(:,:)
        real(dp), allocatable :: Nd(:), s(:), H(:,:), Spur(:), U_sub(:,:), AU(:,:), A_ort(:,:), N_ort(:,:), w(:), Z(:,:)
        integer :: J_max, J, P, N_ph, i, k, N_o
        real(dp) :: A_spur, N_spur, E_spur
        logical :: Ortho

        J_max = Orb_1p1h%J_max
        Ortho = Params%Calc%TDA%Ortho
        allocate(E(2,0:J_max), X(2,0:J_max))
        write(*,'(/,a)') 'Starting diagonalization of the TDA system of equations ...'

        do J = 0, J_max
            do P = 1, 2
                N_ph = Orb_1p1h%N(P,J)
                write(*,'(a)') TAB // TAB // 'Diagonalizing:' // TAB // 'J = ' // istr(J) // ', P = ' // pm_str(P)
                allocate(E(P,J)%v(N_ph), X(P,J)%v(N_ph,N_ph))
                if (N_ph == 0) cycle
                allocate(Nd(N_ph))
                do i = 1, N_ph
                    Nd(i) = N(P,J)%v(i,i)
                end do

                if ((.not. Ortho) .or. .not. (J == 1 .and. P == 2) .or. N_ph < 2) then
                    s = 1.0_dp / sqrt(Nd)
                    allocate(H(N_ph,N_ph), w(N_ph), Z(N_ph,N_ph))
                    do k = 1, N_ph
                        do i = 1, N_ph
                            H(i,k) = s(i) * A(P,J)%v(i,k) * s(k)
                        end do
                    end do
                    call eig_sym(H,w,Z)
                    E(P,J)%v = w
                    do k = 1, N_ph
                        X(P,J)%v(:,k) = s * Z(:,k)
                    end do
                    deallocate(H, w, Z)
                else
                    N_o = N_ph - 1
                    Spur = ph_spurious_ini(Orb,Orb_1p1h,TrOp)
                    U_sub = ph_spurious_basis(Spur)
                    AU = mm(A(P,J)%v,U_sub)
                    A_ort = mtm(U_sub,AU)
                    do k = 1, N_o
                        AU(:,k) = Nd * U_sub(:,k)
                    end do
                    N_ort = mtm(U_sub,AU)
                    allocate(w(N_o), Z(N_o,N_o))
                    call eig_sym_gen(A_ort,N_ort,w,Z)
                    E(P,J)%v(2:N_ph) = w
                    X(P,J)%v(:,2:N_ph) = mm(U_sub,Z)

                    ! 1x1 problem in the spurious subspace ...
                    A_spur = dot_product(Spur,matmul(A(P,J)%v,Spur))
                    N_spur = dot_product(Spur,Nd * Spur)
                    E_spur = A_spur / N_spur
                    write(*,'(a)') TAB // TAB // TAB // 'Spurious 1- TDA level   ...   E = (' // jstr(round_digits(E_spur,7)) // ') MeV'
                    E(P,J)%v(1) = E_spur
                    X(P,J)%v(:,1) = Spur / sqrt(N_spur)
                    deallocate(w, Z, AU, A_ort, N_ort, U_sub, Spur)
                end if

                if (any(E(P,J)%v < -1.0e-8_dp)) &
                    write(*,'(a)') TAB // TAB // TAB // 'Negative TDA eigenvalue detected in channel:' // TAB // 'J = ' // istr(J) // &
                                   ', P = ' // pm_str(P)
                deallocate(Nd)
            end do
        end do
        write(*,'(a)') TAB // 'Finished diagonalization of the TDA system of equations ...'
    end subroutine TDA_diagonalize

    ! Exact TDA transition kernel D of the multipole (J,P) ... Gamma^t_ab(nu) = sum_ph D[(a,b,t),ph] X[ph,nu] ... counterpart
        ! of TDA_transition_kernel (TDA_transitions.jl) ... pairs: the 1p-1h pairs + the pairs coupled by a non-diagonal Rho
        ! (Reference_TBCF = false), or all (a,b,t) of parity P coupled to J (Reference_TBCF = true) ...
    subroutine TDA_transition_kernel(Params,J,P,Orb,Orb_1p1h,Orb_NN,Rho,Sigma,Pairs,D)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: J, P
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(in) :: Sigma
        type(ind2b_t), allocatable, intent(out) :: Pairs(:)
        real(dp), allocatable, intent(out) :: D(:,:)
        real(dp), parameter :: Tol_rho = 1.0e-12_dp
        type(ind2b_t), allocatable :: Buf(:)
        integer, allocatable :: Pair_Map(:,:,:)
        logical :: TBCF
        integer :: N_ph, a_max, N_pair, i_ph, c, d_, a_p, a_h, T_ph, t, ia, ib, t_1, t_2, l_p, j_p, l_h, j_h, k
        integer :: a_a, a_b, l_a, j_a, j_b, P_ha, I, i_pair
        real(dp) :: R_pc, Phase_X, K_1, K_2, F

        TBCF = Params%Calc%TDA%Reference_TBCF
        N_ph = Orb_1p1h%N(P,J)
        a_max = size(Orb)
        allocate(Pair_Map(a_max,a_max,2), Buf(N_ph + 2*a_max*a_max))
        Pair_Map = 0
        N_pair = 0

        if (.not. TBCF) then
            do i_ph = 1, N_ph
                N_pair = N_pair + 1
                Buf(N_pair) = Orb_1p1h%Ind(P,J)%v(i_ph)
                Pair_Map(Buf(N_pair)%a,Buf(N_pair)%b,(Buf(N_pair)%T + 3)/2) = N_pair
            end do
            do i_ph = 1, N_ph
                a_p = Orb_1p1h%Ind(P,J)%v(i_ph)%a
                a_h = Orb_1p1h%Ind(P,J)%v(i_ph)%b
                T_ph = Orb_1p1h%Ind(P,J)%v(i_ph)%T
                do c = 1, a_max
                    do d_ = 1, a_max
                        if (abs(real(kronecker_delta(c,a_p),dp) - rho_t(a_p,c,T_ph)) * abs(rho_t(d_,a_h,T_ph)) > Tol_rho) then
                            if (Pair_Map(c,d_,(T_ph + 3)/2) == 0) then
                                N_pair = N_pair + 1
                                Buf(N_pair) = ind2b_t(c,d_,T_ph)
                                Pair_Map(c,d_,(T_ph + 3)/2) = N_pair
                            end if
                        end if
                    end do
                end do
            end do
        else
            do t = -1, 1, 2
                do ia = 1, a_max
                    do ib = 1, a_max
                        if (parity_of(Orb(ia)%l + Orb(ib)%l) == P .and. triangle(Orb(ia)%j,Orb(ib)%j,2*J)) then
                            N_pair = N_pair + 1
                            Buf(N_pair) = ind2b_t(Orb(ia)%a,Orb(ib)%a,t)
                            Pair_Map(Orb(ia)%a,Orb(ib)%a,(t + 3)/2) = N_pair
                        end if
                    end do
                end do
            end do
        end if
        Pairs = Buf(1:N_pair)

        allocate(D(N_pair,N_ph))
        D = 0.0_dp

        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(a_p,a_h,T_ph,l_p,j_p,l_h,j_h,t_1,t_2,c,d_,R_pc,i_pair,Phase_X,k,a_a,a_b,l_a,j_a,j_b,P_ha,K_1,K_2,I,F)
        do i_ph = 1, N_ph
            a_p = Orb_1p1h%Ind(P,J)%v(i_ph)%a
            a_h = Orb_1p1h%Ind(P,J)%v(i_ph)%b
            T_ph = Orb_1p1h%Ind(P,J)%v(i_ph)%T
            l_p = Orb(a_p)%l
            j_p = Orb(a_p)%j
            l_h = Orb(a_h)%l
            j_h = Orb(a_h)%j
            t_1 = (T_ph + 3)/2
            t_2 = (-T_ph + 3)/2

            ! 1-body part ... (1 - Rho)_pc Rho_dh ...
            do c = 1, a_max
                if (Orb(c)%l /= l_p .or. Orb(c)%j /= j_p) cycle
                R_pc = real(kronecker_delta(c,a_p),dp) - rho_t(a_p,c,T_ph)
                do d_ = 1, a_max
                    if (Orb(d_)%l /= l_h .or. Orb(d_)%j /= j_h) cycle
                    i_pair = Pair_Map(c,d_,t_1)
                    if (i_pair == 0) cycle
                    D(i_pair,i_ph) = D(i_pair,i_ph) + R_pc * rho_t(d_,a_h,T_ph) * phase(J)
                end do
            end do

            ! 2-body part (reference 2-body correlation function Sigma) ...
            if (TBCF) then
                Phase_X = phase(1 + (j_p + j_h)/2)
                do k = 1, N_pair
                    if (Pairs(k)%T /= T_ph) cycle
                    a_a = Pairs(k)%a
                    a_b = Pairs(k)%b
                    l_a = Orb(a_a)%l
                    j_a = Orb(a_a)%j
                    j_b = Orb(a_b)%j
                    P_ha = parity_of(l_h + l_a)
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
                    D(Pair_Map(a_a,a_b,t_1),i_ph) = D(Pair_Map(a_a,a_b,t_1),i_ph) + Phase_X * K_1
                    D(Pair_Map(a_a,a_b,t_2),i_ph) = D(Pair_Map(a_a,a_b,t_2),i_ph) + Phase_X * K_2
                end do
            end if
        end do
        !$omp end parallel do

    contains

        pure real(dp) function rho_t(a,b,T)
            integer, intent(in) :: a, b, T
            if (T == -1) then
                rho_t = Rho%p(a,b)
            else
                rho_t = Rho%n(a,b)
            end if
        end function rho_t

    end subroutine TDA_transition_kernel

    ! Exact reduced transition matrix elements w.r.t. the reference ground state ... counterpart of TDA_rM & TDA_rM_calculate ...
    function TDA_rM(Params,Orb,Orb_1p1h,Orb_NN,M,Rho,Sigma,X) result(rM)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(orb2b_t), intent(in) :: Orb_NN
        type(tr1b_t), intent(in) :: M
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(in) :: Sigma
        type(rmat_t), intent(in) :: X(:,0:)
        type(reduced_multipole_t) :: rM
        complex(dp), allocatable :: pE0(:,:), nE0(:,:), pE1(:,:), nE1(:,:), pE2(:,:), nE2(:,:), pE3(:,:), nE3(:,:)

        write(*,'(/,a)') 'Calculating reduced transition matrix elements rM ... (Exact)'
        call calc(0,1,ph_operators(M,0),pE0,nE0)
        call calc(1,2,ph_operators(M,1),pE1,nE1)
        call calc(2,1,ph_operators(M,2),pE2,nE2)
        call calc(3,2,ph_operators(M,3),pE3,nE3)
        call ph_rM_assemble(pE0,nE0,pE1,nE1,pE2,nE2,pE3,nE3,rM)
        write(*,'(a)') TAB // 'Reduced 1-body transition matrix elements rM calculated ...'

    contains

        subroutine calc(J,P,Ops,prM,nrM)
            integer, intent(in) :: J, P
            type(o1b_t), intent(in) :: Ops(:)
            complex(dp), allocatable, intent(out) :: prM(:,:), nrM(:,:)
            type(ind2b_t), allocatable :: Pairs(:)
            real(dp), allocatable :: D(:,:), Gamma(:,:), pM_pair(:,:), nM_pair(:,:)
            integer :: k, i_o, N_o
            call TDA_transition_kernel(Params,J,P,Orb,Orb_1p1h,Orb_NN,Rho,Sigma,Pairs,D)
            N_o = size(Ops)
            allocate(pM_pair(N_o,size(Pairs)), nM_pair(N_o,size(Pairs)))
            pM_pair = 0.0_dp
            nM_pair = 0.0_dp
            do k = 1, size(Pairs)
                do i_o = 1, N_o
                    if (Pairs(k)%T == -1) then
                        pM_pair(i_o,k) = Ops(i_o)%p(Pairs(k)%a,Pairs(k)%b)
                    else
                        nM_pair(i_o,k) = Ops(i_o)%n(Pairs(k)%a,Pairs(k)%b)
                    end if
                end do
            end do
            Gamma = mm(D,X(P,J)%v)
            prM = cmplx(mm(pM_pair,Gamma), 0.0_dp, kind=dp)
            nrM = cmplx(mm(nM_pair,Gamma), 0.0_dp, kind=dp)
        end subroutine calc

    end function TDA_rM

    ! Reduced multipole structure from the per-multipole (operator, state) arrays ... (shared by TDA_rM & RPA_rM) ...
    subroutine ph_rM_assemble(pE0,nE0,pE1,nE1,pE2,nE2,pE3,nE3,rM)
        complex(dp), intent(in) :: pE0(:,:), nE0(:,:), pE1(:,:), nE1(:,:), pE2(:,:), nE2(:,:), pE3(:,:), nE3(:,:)
        type(reduced_multipole_t), intent(out) :: rM
        call pn_row(rM%E0,pE0,nE0,1)
        call pn_row(rM%E1,pE1,nE1,1)
        call pn_row(rM%E2,pE2,nE2,1)
        call pn_row(rM%E3,pE3,nE3,1)
        call pn_row(rM%E1_VC,pE1,nE1,2)
        call pn_row(rM%E1_VS,pE1,nE1,3)
        call pn_row(rM%E1_TC,pE1,nE1,4)
        call pn_row(rM%E1_TS,pE1,nE1,5)
        call pn_row(rM%E1_sTC,pE1,nE1,6)
        call pn_row(rM%E1_sTS,pE1,nE1,7)
        call pn_row(rM%E1_C,pE1,nE1,8)
        call pn_row(rM%E1_sC,pE1,nE1,9)
    contains

        ! Row k of the (operator, state) arrays ...
        subroutine pn_row(V,pM,nM,k)
            type(pncvector_t), intent(out) :: V
            complex(dp), intent(in) :: pM(:,:), nM(:,:)
            integer, intent(in) :: k
            integer :: nu
            allocate(V%p(size(pM,2)), V%n(size(nM,2)))
            do nu = 1, size(pM,2)
                V%p(nu) = pM(k,nu)
                V%n(nu) = nM(k,nu)
            end do
        end subroutine pn_row

    end subroutine ph_rM_assemble

    ! Export of the TDA results ... counterpart of TDA_export (TDA_export.jl) ... Out_Dir = io/<Path>/<Name>, Rho = reference
        ! OBDM in the reference basis ...
    subroutine TDA_export(Params,Orb,Orb_1p1h,E,X,rB,C,Rho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E(:,0:)
        type(rmat_t), intent(in) :: X(:,0:)
        type(reduced_transition_t), intent(in) :: rB
        type(o1b_t), intent(in) :: C, Rho
        character(len=:), allocatable :: Out_Dir
        integer(i8) :: t0

        Out_Dir = io_dir(Params,Params%Calc%TDA%Name)
        call timer_start(t0)
        call TDA_summary(Params,Orb_1p1h,E,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        call TDA_binary_export(Orb_1p1h,E,X,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        call TDA_spectrum_export(Orb_1p1h,E,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        write(*,'(/,a)') 'Preparing plot-ready export of TDA solutions ...'
        write(*,'(a)') TAB // 'Performing the export of plot-ready TDA spectra ...'
        call ph_plot_spectrum(Out_Dir // '/spectra/TDA_plot.dat',Orb_1p1h,E)
        write(*,'(a)') TAB // TAB // 'TDA plot-ready spectra successfully exported ...'
        call timer_report(t0)
        call timer_start(t0)
        call TDA_amplitudes_export(Orb,Orb_1p1h,E,X,Out_Dir)
        call timer_report(t0)
        call timer_start(t0)
        call ph_transitions_export(Params,'TDA',Out_Dir,Params%Calc%TDA%Ortho,Orb,Orb_1p1h,C,Rho,E,rB)
        call timer_report(t0)
    end subroutine TDA_export

    subroutine TDA_summary(Params,Orb_1p1h,E,Out_Dir)
        type(params_t), intent(in) :: Params
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E(:,0:)
        character(len=*), intent(in) :: Out_Dir
        integer :: u, J, P, nu, N_neg, J_min, P_min
        real(dp) :: E_min
        logical :: Ortho

        Ortho = Params%Calc%TDA%Ortho
        write(*,'(/,a)') 'Writing TDA summary files ...'
        N_neg = 0
        E_min = huge(1.0_dp)
        J_min = -1
        P_min = 0
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                do nu = 1, size(E(P,J)%v)
                    if (Ortho .and. J == 1 .and. P == 2 .and. nu == 1) cycle
                    if (E(P,J)%v(nu) < -1.0e-8_dp) N_neg = N_neg + 1
                    if (E(P,J)%v(nu) < E_min) then
                        E_min = E(P,J)%v(nu)
                        J_min = J
                        P_min = P
                    end if
                end do
            end do
        end do

        open(newunit=u, file=Out_Dir // '/TDA_summary.dat', status='replace', action='write')
        call ph_summary_header(u,Params,'Spherical Tamm-Dancoff Approximation review:')
        write(u,'(/,a)') 'TDA parameters:'
        write(u,'(a)') TAB // "Name = ''" // trim(Params%Calc%TDA%Name) // "'',   Reference = ''" // trim(Params%Calc%TDA%Reference) // &
                       "'',   Reference_TBCF = " // lstr(Params%Calc%TDA%Reference_TBCF) // ',   Ortho = ' // lstr(Ortho)
        write(u,'(a)') TAB // '1p-1h space tolerances:   Tol_hol = ' // jstr(Params%Calc%TDA%Tol_hol) // ',   Tol_par = ' // &
                       jstr(Params%Calc%TDA%Tol_par)
        if (Params%Calc%TDA%Reference_TBCF) then
            write(u,'(a)') TAB // 'Reference 2-body correlation function is INCLUDED in the matrix A & in the (exact) transition ' // &
                           'matrix elements ...'
        else
            write(u,'(a)') TAB // 'Reference 2-body correlation function is NOT included ...'
        end if
        write(u,'(a)') TAB // "Transition matrix elements ... ''Exact'' with respect to the reference ground state ..."
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
                                       istr(J_min) // ', P = ' // pm_str(P_min) // ')'
        call ph_dimensions_block(u,Orb_1p1h,'TDA')
        close(u)
        write(*,'(a)') TAB // 'TDA summary file exported ...'
    end subroutine TDA_summary

    ! bin/TDA_X.bin (row-major Float64, permutedims per channel) & bin/TDA_E.bin ...
    subroutine TDA_binary_export(Orb_1p1h,E,X,Out_Dir)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E(:,0:)
        type(rmat_t), intent(in) :: X(:,0:)
        character(len=*), intent(in) :: Out_Dir
        integer :: u1, u2, J, P

        write(*,'(/,a)') 'Performing binary export of TDA energies E & amplitudes X ...'
        open(newunit=u1, file=Out_Dir // '/bin/TDA_X.bin', access='stream', form='unformatted', status='replace', action='write')
        open(newunit=u2, file=Out_Dir // '/bin/TDA_E.bin', access='stream', form='unformatted', status='replace', action='write')
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                if (Orb_1p1h%N(P,J) == 0) cycle
                write(u1) transpose(X(P,J)%v)
                write(u2) E(P,J)%v
            end do
        end do
        close(u1)
        close(u2)
        write(*,'(a)') TAB // 'Binary export of TDA energies E & amplitudes X completed ...'
    end subroutine TDA_binary_export

    subroutine TDA_spectrum_export(Orb_1p1h,E,Out_Dir)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E(:,0:)
        character(len=*), intent(in) :: Out_Dir
        integer :: u, J, P, nu

        write(*,'(/,a)') 'Preparing export of TDA spectra ...'
        open(newunit=u, file=Out_Dir // '/spectra/TDA.dat', status='replace', action='write')
        write(u,'(a)') cfmt_s('J',5,.true.) // ' ' // cfmt_s('P',5,.true.) // ' ' // cfmt_s('E',20,.true.)
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                do nu = 1, Orb_1p1h%N(P,J)
                    write(u,'(a)') cfmt_d(J,5,.true.) // ' ' // cfmt_s(pm_str(P),5,.true.) // ' ' // cfmt_f(E(P,J)%v(nu),20,6,.true.)
                end do
                write(u,'(a)') ''
                write(u,'(a)') ''
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'Export of TDA spectra successfully finished ...'
    end subroutine TDA_spectrum_export

    ! amplitudes/TDA_amplitudes_ordered.dat & TDA_amplitudes.dat (energy ordered) ... metric norm sum (n_h - n_p)|X|^2 &
        ! the dominant 1p-1h component ...
    subroutine TDA_amplitudes_export(Orb,Orb_1p1h,E,X,Out_Dir)
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_1p1h_t), intent(in) :: Orb_1p1h
        type(rvec_t), intent(in) :: E(:,0:)
        type(rmat_t), intent(in) :: X(:,0:)
        character(len=*), intent(in) :: Out_Dir
        real(dp), allocatable :: L_E(:), L_x(:), L_w(:), dn(:)
        integer, allocatable :: L_J(:), L_P(:), L_a(:), L_b(:), L_T(:), Sort(:), Ident(:)
        integer :: N_lev, cnt, J, P, nu, ph, i_max
        real(dp) :: w, w_max, x2

        write(*,'(/,a)') 'Preparing export of the TDA amplitudes X ...'
        N_lev = sum(Orb_1p1h%N)
        allocate(L_E(N_lev), L_x(N_lev), L_w(N_lev), L_J(N_lev), L_P(N_lev), L_a(N_lev), L_b(N_lev), L_T(N_lev), Sort(N_lev), &
                 Ident(N_lev))
        cnt = 0
        do J = 0, Orb_1p1h%J_max
            do P = 1, 2
                dn = ph_metric(Orb,Orb_1p1h,P,J)
                do nu = 1, Orb_1p1h%N(P,J)
                    cnt = cnt + 1
                    x2 = 0.0_dp
                    w_max = -1.0_dp
                    i_max = 0
                    do ph = 1, Orb_1p1h%N(P,J)
                        w = dn(ph) * X(P,J)%v(ph,nu)**2
                        x2 = x2 + w
                        if (w > w_max) then
                            w_max = w
                            i_max = ph
                        end if
                    end do
                    L_E(cnt) = E(P,J)%v(nu)
                    L_x(cnt) = x2
                    L_J(cnt) = J
                    L_P(cnt) = P
                    L_a(cnt) = Orb_1p1h%Ind(P,J)%v(i_max)%a
                    L_b(cnt) = Orb_1p1h%Ind(P,J)%v(i_max)%b
                    L_T(cnt) = Orb_1p1h%Ind(P,J)%v(i_max)%T
                    L_w(cnt) = w_max
                end do
            end do
        end do
        do cnt = 1, N_lev
            Ident(cnt) = cnt
        end do
        call write_file(Out_Dir // '/amplitudes/TDA_amplitudes_ordered.dat',Ident)
        call sortperm(L_E,Sort)
        call write_file(Out_Dir // '/amplitudes/TDA_amplitudes.dat',Sort)
        write(*,'(a)') TAB // 'Export of TDA amplitudes X completed ...'

    contains

        subroutine write_file(Path,Order)
            character(len=*), intent(in) :: Path
            integer, intent(in) :: Order(:)
            integer :: u, i, k
            character(len=:), allocatable :: p_str, h_str
            open(newunit=u, file=Path, status='replace', action='write')
            write(u,'(a)') cfmt_s('Ind',8,.true.) // ' ' // cfmt_s('E',14,.true.) // ' ' // cfmt_s('J',4,.true.) // ' ' // &
                           cfmt_s('P',4,.true.) // ' ' // cfmt_s('|X|^2',12,.true.) // ' ' // cfmt_s('t',6,.true.) // ' ' // &
                           cfmt_s('p (n,l,2j)',14,.true.) // ' ' // cfmt_s('h (n,l,2j)',14,.true.) // ' ' // cfmt_s('max weight',12,.true.)
            do i = 1, size(Order)
                k = Order(i)
                p_str = '(' // istr(Orb(L_a(k))%n) // ',' // istr(Orb(L_a(k))%l) // ',' // istr(Orb(L_a(k))%j) // ')'
                h_str = '(' // istr(Orb(L_b(k))%n) // ',' // istr(Orb(L_b(k))%l) // ',' // istr(Orb(L_b(k))%j) // ')'
                write(u,'(a)') cfmt_d(i,8,.true.) // ' ' // cfmt_f(L_E(k),14,6,.true.) // ' ' // cfmt_d(L_J(k),4,.true.) // ' ' // &
                               cfmt_s(pm_str(L_P(k)),4,.true.) // ' ' // cfmt_f(L_x(k),12,6,.true.) // ' ' // &
                               cfmt_s(merge('p', 'n', L_T(k) == -1),6,.true.) // ' ' // cfmt_s(p_str,14,.true.) // ' ' // &
                               cfmt_s(h_str,14,.true.) // ' ' // cfmt_f(L_w(k),12,6,.true.)
            end do
            close(u)
        end subroutine write_file

    end subroutine TDA_amplitudes_export

end module nj_tda
