! Spherical HF Renormalized-RPA ... counterpart of HF_RRPA/HF_RRPA_Hamiltonian.jl, HF_RRPA_Allocate.jl,
    ! HF_RRPA_Diagonalize.jl, HF_RRPA_Spurious.jl, HF_RRPA_OBDM.jl, HF_RRPA_Energy.jl & HF_RRPA_Transitions.jl ...
module nj_hf_rrpa
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_wigner, only : f6j
    use nj_orbitals, only : orbitals_ph_list
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn, O2b_transformation
    use nj_v2b, only : V2b_residual_no2b
    use nj_v3b, only : V3b_no2b
    use nj_tr1b, only : tr1b_t
    use nj_hf_rpa, only : RPA_eigen_solve, RPA_eigen_solve_c, RPA_normalize, ph_orbitals
    implicit none
    private

    public :: HF_RRPA_h1b, HF_RRPA_V2b, HF_RRPA_V2b_transform, HF_RRPA_allocate, HF_RRPA0_diagonalize
    public :: HF_RRPA_diagonalize, HF_RRPA_OBDM_iteration, HF_RRPA_energy, HF_RRPA_energy_density
    public :: HF_RRPA_energy_bosonic, HF_RRPA_rM

contains

    ! Mean-field 1-body Hamiltonian for density Rho (LHO basis) ... NN + NO2B NNN + kinetic terms ...
    function HF_RRPA_h1b(Params,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN) result(H)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t) :: H
        integer :: A_mass, N_2max, N_3max, a_max, a, d, b, e, J, c, f, P_ab, P_abc
        integer :: n_a, l_a, j_a, n_d, l_d, j_d, n_b, l_b, j_b, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        real(dp) :: Hat, Hat_2b, pSum, nSum, pRho_be, nRho_be, pRho_cf, nRho_cf, ME001, ME101, ME011, ME111, ME113
        character(len=32) :: CMS

        A_mass = Params%Calc%A
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        CMS = Params%Calc%CMS
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(H%p(a_max,a_max), H%n(a_max,a_max))
        H%p = 0.0_dp
        H%n = 0.0_dp

        !$omp parallel do schedule(dynamic) private(n_a,l_a,j_a,Hat,d,l_d,j_d,n_d,pSum,nSum,b,n_b,l_b,P_ab,j_b,e,n_e,l_e,j_e, &
        !$omp pRho_be,nRho_be,J,Hat_2b,c,n_c,l_c,P_abc,j_c,f,n_f,l_f,j_f,pRho_cf,nRho_cf,ME001,ME101,ME011,ME111,ME113)
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            Hat = 1.0_dp / (real(j_a,dp) + 1.0_dp)
            do d = 1, a
                l_d = Orb(d)%l
                j_d = Orb(d)%j
                if (.not. (l_a == l_d .and. j_a == j_d)) cycle
                n_d = Orb(d)%n
                pSum = 0.0_dp
                nSum = 0.0_dp
                do b = 1, a_max
                    n_b = Orb(b)%n
                    l_b = Orb(b)%l
                    if ((2*(n_a + n_b) + l_a + l_b) > N_2max) cycle
                    P_ab = mod(l_a + l_b, 2) + 1
                    j_b = Orb(b)%j
                    do e = 1, a_max
                        n_e = Orb(e)%n
                        l_e = Orb(e)%l
                        j_e = Orb(e)%j
                        if (.not. (l_b == l_e .and. j_b == j_e .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max)) cycle
                        pRho_be = Rho%p(b,e)
                        nRho_be = Rho%n(b,e)
                        do J = abs(j_d - j_e)/2, (j_d + j_e)/2
                            if (mod(l_a + l_b, 2) == mod(l_d + l_e, 2)) then
                                Hat_2b = (real(2*J,dp) + 1.0_dp) * Hat
                                pSum = pSum + Hat_2b * pRho_be * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                                pSum = pSum + Hat_2b * nRho_be * O2b_pn(a,b,d,e,J,P_ab,V_NN,Orb_NN)
                                nSum = nSum + Hat_2b * nRho_be * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                                nSum = nSum + Hat_2b * pRho_be * O2b_pn(b,a,e,d,J,P_ab,V_NN,Orb_NN)
                            end if
                            do c = 1, a_max
                                n_c = Orb(c)%n
                                l_c = Orb(c)%l
                                if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                                P_abc = mod(l_a + l_b + l_c, 2) + 1
                                j_c = Orb(c)%j
                                do f = 1, a_max
                                    n_f = Orb(f)%n
                                    l_f = Orb(f)%l
                                    if (.not. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                               P_abc == (mod(l_d + l_e + l_f,2) + 1))) cycle
                                    j_f = Orb(f)%j
                                    if (j_c /= j_f) cycle
                                    pRho_cf = Rho%p(c,f)
                                    nRho_cf = Rho%n(c,f)
                                    ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                    ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                    ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                    ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                    ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P_abc,V_NNN,Orb,Orb_NNN)
                                    pSum = pSum + Hat * (0.5_dp*ME113*pRho_be*pRho_cf + 0.25_dp * (ME001 + &
                                           sqrt(1.0_dp/3.0_dp)*ME101 + sqrt(1.0_dp/3.0_dp)*ME011 + 1.0_dp/3.0_dp*ME111 + &
                                           2.0_dp/3.0_dp*ME113)*nRho_be*nRho_cf + 1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*pRho_be*nRho_cf)
                                    nSum = nSum + Hat * (0.5_dp*ME113*nRho_be*nRho_cf + 0.25_dp * (ME001 + &
                                           sqrt(1.0_dp/3.0_dp)*ME101 + sqrt(1.0_dp/3.0_dp)*ME011 + 1.0_dp/3.0_dp*ME111 + &
                                           2.0_dp/3.0_dp*ME113)*pRho_be*pRho_cf + 1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*nRho_be*pRho_cf)
                                end do
                            end do
                        end do
                    end do
                end do
                if (CMS == 'CMS1+2B') then
                    H%p(a,d) = pSum + T%p(a,d) * (1.0_dp - 1.0_dp/real(A_mass,dp))
                    H%n(a,d) = nSum + T%n(a,d) * (1.0_dp - 1.0_dp/real(A_mass,dp))
                else if (CMS == 'CMS2B') then
                    H%p(a,d) = pSum
                    H%n(a,d) = nSum
                else
                    H%p(a,d) = pSum + T%p(a,d)
                    H%n(a,d) = nSum + T%n(a,d)
                end if
                H%p(d,a) = H%p(a,d)
                H%n(d,a) = H%n(a,d)
            end do
        end do
        !$omp end parallel do
    end function HF_RRPA_h1b

    ! Density-dependent residual NN interaction for Rho (LHO) transformed to the basis C ... output silenced ...
    subroutine HF_RRPA_V2b(Params,Orb,Orb_NN_Bare,Orb_NNN_Bare,V_NN_Bare,V_NNN_Bare,Rho,C,V_NN_Res)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN_Bare
        type(orb3b_t), intent(in) :: Orb_NNN_Bare
        type(o2b_t), intent(in) :: V_NN_Bare
        type(v3b_t), intent(in) :: V_NNN_Bare
        type(o1b_t), intent(in) :: Rho, C
        type(o2b_t), intent(out) :: V_NN_Res
        logical :: Verbose_Old
        Verbose_Old = Verbose
        Verbose = .false.
        call V2b_residual_no2b(Params,Orb,Orb_NN_Bare,Orb_NNN_Bare,Rho,V_NN_Bare,V_NNN_Bare,V_NN_Res)
        call O2b_transformation(Params,Orb,Orb_NN_Bare,V_NN_Res,C)
        Verbose = Verbose_Old
    end subroutine HF_RRPA_V2b

    subroutine HF_RRPA_V2b_transform(Params,Orb,Orb_NN_Res,V_NN_Res,C)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN_Res
        type(o2b_t), intent(inout) :: V_NN_Res
        type(o1b_t), intent(in) :: C
        logical :: Verbose_Old
        Verbose_Old = Verbose
        Verbose = .false.
        call O2b_transformation(Params,Orb,Orb_NN_Res,V_NN_Res,C)
        Verbose = Verbose_Old
    end subroutine HF_RRPA_V2b_transform

    ! Renormalized RPA matrices A & B ... occupation-weighted with the effective 1-body Hamiltonian H ...
    subroutine HF_RRPA_allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,V_NN,Rho,H,A,B)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: V_NN
        type(o1b_t), intent(in) :: Rho, H
        type(rmat_t), allocatable, intent(out) :: A(:,:), B(:,:)
        integer :: J_max, J, P, N_ph, ind_ph, ind_qg, ph, qg, t_ph, t_qg, J_r, P_pg, P_pq
        type(spstate_t) :: sp, sh, sq, sg
        real(dp) :: ASum, BSum, ME, Amp

        J_max = 2*Params%Calc%Nmax + 1
        allocate(A(0:J_max,2), B(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                allocate(A(J,P)%v(N_ph,N_ph), B(J,P)%v(N_ph,N_ph))
                !$omp parallel do schedule(dynamic) private(ph,t_ph,sp,sh,ind_qg,qg,t_qg,sq,sg,P_pg,P_pq,ASum,BSum,J_r,ME,Amp)
                do ind_ph = 1, N_ph
                    ph = Orb_Phonon(J,P)%v(ind_ph)
                    t_ph = Phonon(ph)%tz
                    call states(ph,sp,sh)
                    do ind_qg = 1, N_ph
                        qg = Orb_Phonon(J,P)%v(ind_qg)
                        t_qg = Phonon(qg)%tz
                        call states(qg,sq,sg)
                        ASum = 0.0_dp
                        BSum = 0.0_dp
                        P_pg = mod(sp%l + sg%l,2) + 1
                        P_pq = mod(sp%l + sq%l,2) + 1

                        if (t_ph == t_qg) then
                            if (t_ph == -1) then
                                Amp = 0.5_dp * (sqrt((Rho%p(sh%a,sh%a) - Rho%p(sp%a,sp%a)) / (Rho%p(sg%a,sg%a) - Rho%p(sq%a,sq%a))) + &
                                      sqrt((Rho%p(sg%a,sg%a) - Rho%p(sq%a,sq%a)) / (Rho%p(sh%a,sh%a) - Rho%p(sp%a,sp%a)))) * &
                                      (H%p(sp%a,sq%a) * kronecker_delta(sh%a,sg%a) - H%p(sh%a,sg%a) * kronecker_delta(sp%a,sq%a))
                            else
                                Amp = 0.5_dp * (sqrt((Rho%n(sh%a,sh%a) - Rho%n(sp%a,sp%a)) / (Rho%n(sg%a,sg%a) - Rho%n(sq%a,sq%a))) + &
                                      sqrt((Rho%n(sg%a,sg%a) - Rho%n(sq%a,sq%a)) / (Rho%n(sh%a,sh%a) - Rho%n(sp%a,sp%a)))) * &
                                      (H%n(sp%a,sq%a) * kronecker_delta(sh%a,sg%a) - H%n(sh%a,sg%a) * kronecker_delta(sp%a,sq%a))
                            end if
                            ASum = ASum + Amp
                        end if

                        if (mod(sp%l + sg%l, 2) == mod(sq%l + sh%l, 2)) then
                            do J_r = abs(sp%j - sg%j)/2, (sp%j + sg%j)/2
                                if (.not. ((abs(sq%j - sh%j) <= 2*J_r) .and. (2*J_r <= (sq%j + sh%j)))) cycle
                                Amp = weight(t_ph,sp%a,sh%a,t_qg,sq%a,sg%a)
                                if (t_ph == -1 .and. t_qg == -1) then
                                    ME = Amp * iphase((sh%j + sq%j)/2 + J_r) * real(2*J_r + 1,dp) * f6j(sp%j,sh%j,2*J,sq%j,sg%j,2*J_r) * &
                                         O2b_pp(sp%a,sg%a,sh%a,sq%a,J_r,P_pg,V_NN,Orb,Orb_NN)
                                else if (t_ph == -1 .and. t_qg == 1) then
                                    ME = Amp * iphase((sh%j + sq%j)/2 + J_r) * real(2*J_r + 1,dp) * f6j(sp%j,sh%j,2*J,sq%j,sg%j,2*J_r) * &
                                         O2b_pn(sp%a,sg%a,sh%a,sq%a,J_r,P_pg,V_NN,Orb_NN)
                                else if (t_ph == 1 .and. t_qg == -1) then
                                    ME = Amp * iphase((sg%j + sp%j)/2 + J_r) * real(2*J_r + 1,dp) * f6j(sq%j,sg%j,2*J,sp%j,sh%j,2*J_r) * &
                                         O2b_pn(sq%a,sh%a,sg%a,sp%a,J_r,P_pg,V_NN,Orb_NN)
                                else
                                    ME = Amp * iphase((sh%j + sq%j)/2 + J_r) * real(2*J_r + 1,dp) * f6j(sp%j,sh%j,2*J,sq%j,sg%j,2*J_r) * &
                                         O2b_nn(sp%a,sg%a,sh%a,sq%a,J_r,P_pg,V_NN,Orb,Orb_NN)
                                end if
                                ASum = ASum + ME
                            end do
                        end if

                        if (mod(sp%l + sq%l, 2) == mod(sh%l + sg%l, 2)) then
                            do J_r = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                                if (.not. ((abs(sh%j - sg%j) <= 2*J_r) .and. (2*J_r <= (sh%j + sg%j)))) cycle
                                Amp = weight(t_ph,sp%a,sh%a,t_qg,sq%a,sg%a)
                                if (t_ph == -1 .and. t_qg == -1) then
                                    ME = Amp * phase((sh%j + sq%j)/2 + J + J_r) * real(2*J_r + 1,dp) * f6j(sp%j,sh%j,2*J,sg%j,sq%j,2*J_r) * &
                                         O2b_pp(sp%a,sq%a,sh%a,sg%a,J_r,P_pq,V_NN,Orb,Orb_NN)
                                else if (t_ph == -1 .and. t_qg == 1) then
                                    ME = Amp * phase((sh%j + sq%j)/2 + J + J_r) * real(2*J_r + 1,dp) * f6j(sp%j,sh%j,2*J,sg%j,sq%j,2*J_r) * &
                                         O2b_pn(sp%a,sq%a,sh%a,sg%a,J_r,P_pq,V_NN,Orb_NN)
                                else if (t_ph == 1 .and. t_qg == -1) then
                                    ME = Amp * phase((sg%j + sp%j)/2 + J + J_r) * real(2*J_r + 1,dp) * f6j(sq%j,sg%j,2*J,sh%j,sp%j,2*J_r) * &
                                         O2b_pn(sq%a,sp%a,sg%a,sh%a,J_r,P_pq,V_NN,Orb_NN)
                                else
                                    ME = Amp * phase((sh%j + sq%j)/2 + J + J_r) * real(2*J_r + 1,dp) * f6j(sp%j,sh%j,2*J,sg%j,sq%j,2*J_r) * &
                                         O2b_nn(sp%a,sq%a,sh%a,sg%a,J_r,P_pq,V_NN,Orb,Orb_NN)
                                end if
                                BSum = BSum + ME
                            end do
                        end if
                        A(J,P)%v(ind_ph,ind_qg) = ASum
                        B(J,P)%v(ind_ph,ind_qg) = BSum
                    end do
                end do
                !$omp end parallel do
            end do
        end do

    contains

        subroutine states(ph,s_p,s_h)
            integer, intent(in) :: ph
            type(spstate_t), intent(out) :: s_p, s_h
            if (Phonon(ph)%tz == -1) then
                s_p = Particle%p(Phonon(ph)%p)
                s_h = Hole%p(Phonon(ph)%h)
            else
                s_p = Particle%n(Phonon(ph)%p)
                s_h = Hole%n(Phonon(ph)%h)
            end if
        end subroutine states

        ! sqrt((rho_hh - rho_pp) (rho_gg - rho_qq)) with the isospin of each pair ...
        real(dp) function weight(t1,a_p,a_h,t2,a_q,a_g)
            integer, intent(in) :: t1, a_p, a_h, t2, a_q, a_g
            real(dp) :: D1, D2
            if (t1 == -1) then
                D1 = Rho%p(a_h,a_h) - Rho%p(a_p,a_p)
            else
                D1 = Rho%n(a_h,a_h) - Rho%n(a_p,a_p)
            end if
            if (t2 == -1) then
                D2 = Rho%p(a_g,a_g) - Rho%p(a_q,a_q)
            else
                D2 = Rho%n(a_g,a_g) - Rho%n(a_q,a_q)
            end if
            weight = sqrt(D1 * D2)
        end function weight

    end subroutine HF_RRPA_allocate

    ! RRPA(0) & iteration diagonalization ... negative-norm, |Y|^2 > 0.501 & complex-energy modes are zeroed ...
    subroutine HF_RRPA0_diagonalize(Params,A,B,N_nu,E_RPA,X_RPA,Y_RPA)
        type(params_t), intent(in) :: Params
        type(rmat_t), intent(in) :: A(0:,:), B(0:,:)
        integer, intent(in) :: N_nu(0:,:)
        type(cvec_t), allocatable, intent(out) :: E_RPA(:,:)
        type(cmat_t), allocatable, intent(out) :: X_RPA(:,:), Y_RPA(:,:)
        integer :: J_max, J, P, nu
        J_max = 2*Params%Calc%Nmax + 1
        allocate(E_RPA(0:J_max,2), X_RPA(0:J_max,2), Y_RPA(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                call RPA_eigen_solve(A(J,P)%v,B(J,P)%v,E_RPA(J,P)%v,X_RPA(J,P)%v,Y_RPA(J,P)%v)
                call RPA_normalize(X_RPA(J,P)%v,Y_RPA(J,P)%v,Zero_Negative=.true.)
                do nu = 1, N_nu(J,P)
                    if (sum(abs(Y_RPA(J,P)%v(:,nu))**2) > 0.501_dp) then
                        X_RPA(J,P)%v(:,nu) = (0.0_dp, 0.0_dp)
                        Y_RPA(J,P)%v(:,nu) = (0.0_dp, 0.0_dp)
                    end if
                end do
                do nu = 1, N_nu(J,P)
                    if (abs(aimag(E_RPA(J,P)%v(nu))) > 1.0e-3_dp) then
                        X_RPA(J,P)%v(:,nu) = (0.0_dp, 0.0_dp)
                        Y_RPA(J,P)%v(:,nu) = (0.0_dp, 0.0_dp)
                    end if
                end do
            end do
        end do
    end subroutine HF_RRPA0_diagonalize

    ! Final RRPA diagonalization with optional orthogonalization of the spurious 1- state ...
    subroutine HF_RRPA_diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,Rho,X_0,Y_0,E_RPA,X_RPA,Y_RPA)
        type(params_t), intent(in) :: Params
        type(rmat_t), intent(in) :: A(0:,:), B(0:,:)
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(tr1b_t), intent(in) :: TrOp
        type(o1b_t), intent(in) :: Rho
        type(cmat_t), intent(in) :: X_0(0:,:), Y_0(0:,:)
        type(cvec_t), allocatable, intent(out) :: E_RPA(:,:)
        type(cmat_t), allocatable, intent(out) :: X_RPA(:,:), Y_RPA(:,:)
        complex(dp), allocatable :: Spur(:), V(:,:), U(:,:), E_o(:), X_o(:,:), Y_o(:,:)
        integer :: J_max, J, P, N_ph, Spur_Ind, k

        J_max = 2*Params%Calc%Nmax + 1
        allocate(E_RPA(0:J_max,2), X_RPA(0:J_max,2), Y_RPA(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                if (Params%Calc%RRPA%Ortho .and. J == 1 .and. P == 2 .and. N_ph > 0) then
                    Spur = spurious_initialize(N_ph,Orb_Phonon(J,P)%v,X_0(J,P)%v,Y_0(J,P)%v)
                    call spurious_orthogonalize(N_ph,Spur,V,Spur_Ind)
                    allocate(U(N_ph,N_ph-1))
                    do k = 1, Spur_Ind - 1
                        U(:,k) = V(:,k)
                    end do
                    do k = Spur_Ind + 1, N_ph
                        U(:,k-1) = V(:,k)
                    end do
                    ! Plain transpose U^T (not the adjoint) as in the original ...
                    call RPA_eigen_solve_c(matmul(transpose(U),matmul(cmplx(A(J,P)%v,0.0_dp,kind=dp),U)), &
                                           matmul(transpose(U),matmul(cmplx(B(J,P)%v,0.0_dp,kind=dp),U)),E_o,X_o,Y_o)
                    allocate(E_RPA(J,P)%v(N_ph), X_RPA(J,P)%v(N_ph,N_ph), Y_RPA(J,P)%v(N_ph,N_ph))
                    E_RPA(J,P)%v(1) = (0.0_dp, 0.0_dp)
                    E_RPA(J,P)%v(2:N_ph) = E_o
                    X_RPA(J,P)%v(:,1) = Spur
                    Y_RPA(J,P)%v(:,1) = (0.0_dp, 0.0_dp)
                    X_RPA(J,P)%v(:,2:N_ph) = zmm(U,X_o)
                    Y_RPA(J,P)%v(:,2:N_ph) = zmm(U,Y_o)
                    deallocate(U)
                else
                    call RPA_eigen_solve(A(J,P)%v,B(J,P)%v,E_RPA(J,P)%v,X_RPA(J,P)%v,Y_RPA(J,P)%v)
                end if
                call RPA_normalize(X_RPA(J,P)%v,Y_RPA(J,P)%v)
            end do
        end do

    contains

        ! Spurious state from the E1 operator weighted by the previous RRPA amplitudes ...
        function spurious_initialize(N_ph,Orb_Ph,X0,Y0) result(S)
            integer, intent(in) :: N_ph, Orb_Ph(:)
            complex(dp), intent(in) :: X0(:,:), Y0(:,:)
            complex(dp) :: S(N_ph)
            integer :: ph, qg, nu, Ind_ph, Ind_qg, a_q, a_g
            complex(dp) :: Acc
            real(dp) :: Op
            do ph = 1, N_ph
                Ind_ph = Orb_Ph(ph)
                Acc = (0.0_dp, 0.0_dp)
                do nu = 1, N_ph
                    do qg = 1, N_ph
                        Ind_qg = Orb_Ph(qg)
                        if (Phonon(Ind_qg)%tz /= Phonon(Ind_ph)%tz) cycle
                        call ph_orbitals(Ind_qg,Phonon,Particle,Hole,a_q,a_g)
                        if (Phonon(Ind_ph)%tz == 1) then
                            Op = TrOp%E1%n(a_q,a_g) * sqrt(Rho%n(a_g,a_g) - Rho%n(a_q,a_q))
                        else
                            Op = TrOp%E1%p(a_q,a_g) * sqrt(Rho%p(a_g,a_g) - Rho%p(a_q,a_q))
                        end if
                        Acc = Acc + Op * (conjg(Y0(qg,nu)) - conjg(X0(qg,nu))) * (X0(ph,nu) - Y0(ph,nu))
                    end do
                end do
                S(ph) = Acc
            end do
            where (abs(S) < 1.0e-6_dp) S = (0.0_dp, 0.0_dp)
            S = S / sqrt(sum(abs(S)**2))
        end function spurious_initialize

        subroutine spurious_orthogonalize(N_ph,S,Uo,Ind)
            integer, intent(in) :: N_ph
            complex(dp), intent(in) :: S(:)
            complex(dp), allocatable, intent(out) :: Uo(:,:)
            integer, intent(out) :: Ind
            complex(dp), allocatable :: Id(:,:), u_v(:)
            complex(dp) :: proj
            integer :: a, b, nu, mu
            Ind = 0
            allocate(Id(N_ph,N_ph))
            Id = (0.0_dp, 0.0_dp)
            do a = 1, N_ph
                Id(a,a) = (1.0_dp, 0.0_dp)
            end do
            do a = 1, N_ph
                proj = (0.0_dp, 0.0_dp)
                do b = 1, N_ph
                    proj = proj + S(b) * Id(b,a)
                end do
                do b = 1, N_ph
                    Id(b,a) = Id(b,a) - proj * S(b)
                end do
            end do
            do a = 1, N_ph
                Id(:,a) = Id(:,a) / sqrt(sum(abs(Id(:,a))**2))
            end do
            allocate(Uo(N_ph,N_ph))
            Uo = (0.0_dp, 0.0_dp)
            Uo(:,1) = Id(:,1)
            do nu = 2, N_ph
                u_v = Id(:,nu)
                do mu = 1, nu - 1
                    u_v = u_v - dot_product(Uo(:,mu),u_v) * Uo(:,mu)
                end do
                if (sqrt(sum(abs(u_v)**2)) < 1.0e-4_dp) Ind = nu
                Uo(:,nu) = u_v / sqrt(sum(abs(u_v)**2))
            end do
            if (Ind == 0) call fatal('RRPA spurious state orthogonalization failed ...')
            Uo(:,Ind) = S
        end subroutine spurious_orthogonalize

    end subroutine HF_RRPA_diagonalize

    ! Self-consistent RRPA OBDM in its natural-orbital basis ... X & Y are transformed in place (as in the original) ...
    subroutine HF_RRPA_OBDM_iteration(Params,N_nu,N_Particle,Particle,N_Hole,Hole,X_RPA,Y_RPA,Rho_HF,Rho_Ini,Rho_out,C_out)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(cmat_t), intent(inout) :: X_RPA(0:,:), Y_RPA(0:,:)
        type(o1b_t), intent(in) :: Rho_HF, Rho_Ini
        type(o1b_t), intent(out) :: Rho_out, C_out
        real(dp), parameter :: Tol_Float32 = 1.0e-8_dp
        real(dp), allocatable :: pRho_old(:,:), nRho_old(:,:), pRho(:,:), nRho(:,:), pC(:,:), nC(:,:), pN(:), nN(:)
        real(dp), allocatable :: pC_RRPA(:,:), nC_RRPA(:,:)
        type(pnarray4_t) :: ParticleHole
        type(obdm_iteration_pairs_t) :: Pairs
        type(cmat_t), allocatable :: pM(:,:), nM(:,:)
        real(dp) :: delta, delta_odd, delta_even, epsilon, dSum
        integer :: J_max, a_max, Iteration, Iteration_max, J, P, a, b

        J_max = 2*Params%Calc%Nmax + 1
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        pC_RRPA = identity(a_max)
        nC_RRPA = identity(a_max)
        pRho_old = Rho_Ini%p
        nRho_old = Rho_Ini%n
        allocate(pRho(a_max,a_max), nRho(a_max,a_max), pC(a_max,a_max), nC(a_max,a_max), pN(a_max), nN(a_max))
        pRho = 0.0_dp
        nRho = 0.0_dp
        delta = 1.0_dp
        delta_odd = 1.0_dp
        delta_even = 1.0_dp
        epsilon = Params%Calc%RRPA%Tol
        Iteration = 0
        Iteration_max = Params%Calc%RRPA%Imax

        ParticleHole = orbitals_ph_list(J_max,N_Particle,Particle,N_Hole,Hole)
        Pairs = iteration_list()
        allocate(pM(0:J_max,2), nM(0:J_max,2))
        do J = 0, J_max
            do P = 1, 2
                allocate(pM(J,P)%v(N_nu(J,P),N_nu(J,P)), nM(J,P)%v(N_nu(J,P),N_nu(J,P)))
            end do
        end do

        write(*,'(/,a)') 'Starting HF-RRPA OBDM calculation ...'
        do while ((delta > epsilon) .and. (Iteration_max > Iteration))
            call obdm(pRho,nRho)

            pRho = 0.5_dp * (pRho + transpose(pRho))
            nRho = 0.5_dp * (nRho + transpose(nRho))
            if (Params%Calc%RRPA%dOBDM) then
                do a = 1, a_max
                    do b = 1, a_max
                        if (a /= b) then
                            pRho(a,b) = 0.0_dp
                            nRho(a,b) = 0.0_dp
                        end if
                    end do
                end do
            end if
            where (abs(pRho) < Tol_Float32) pRho = 0.0_dp
            where (abs(nRho) < Tol_Float32) nRho = 0.0_dp

            call eig_sym(pRho,pN,pC)
            call eig_sym(nRho,nN,nC)
            where (abs(pC) < Tol_Float32) pC = 0.0_dp
            where (abs(nC) < Tol_Float32) nC = 0.0_dp

            call reorder(pN,pC)
            call reorder(nN,nC)

            if ((.not. occupation_check()) .or. Iteration == Iteration_max - 1) then
                write(*,'(/,a)') 'OBDM occupation finished - Singular solution - HF-RRPA iteration is about to finish with the ' // &
                                 'last stable OBDM solution ...'
                Rho_out = Rho_Ini
                C_out%p = identity(a_max)
                C_out%n = identity(a_max)
                return
            end if

            pC_RRPA = mm(pC_RRPA,pC)
            nC_RRPA = mm(nC_RRPA,nC)
            where (abs(pC_RRPA) < Tol_Float32) pC_RRPA = 0.0_dp
            where (abs(nC_RRPA) < Tol_Float32) nC_RRPA = 0.0_dp

            pRho = diagm(pN)
            nRho = diagm(nN)
            pRho_old = mtmm(pC,pRho_old)
            nRho_old = mtmm(nC,nRho_old)

            dSum = 0.0_dp
            do a = 1, a_max
                do b = 1, a_max
                    dSum = dSum + (abs(pRho(a,b) - pRho_old(a,b)) + abs(nRho(a,b) - nRho_old(a,b))) / real(2*a_max**2,dp)
                end do
            end do

            call transform_amplitudes(pC,nC)

            delta = dSum
            Iteration = Iteration + 1

            if (mod(Iteration,2) == 0 .and. Iteration > 100) then
                if (abs(delta - delta_even) < 10.0_dp * epsilon) then
                    write(*,'(/,a)') 'OBDM iteration finished - cycled degenerate solution found ...'
                    exit
                end if
                delta_even = delta
            else if (mod(Iteration,2) == 1 .and. Iteration > 100) then
                if (abs(delta - delta_odd) < 10.0_dp * epsilon) then
                    write(*,'(/,a)') 'OBDM iteration finished - cycled degenerate solution found ...'
                    exit
                end if
                delta_odd = delta
            end if
            if (Iteration > 150 .and. delta < 100.0_dp * epsilon) then
                write(*,'(/,a)') 'OBDM iteration finished with lower precisions 10 x epsilon ...'
                exit
            end if
            if (Iteration > 200 .and. delta < 1000.0_dp * epsilon) then
                write(*,'(/,a)') 'OBDM iteration finished with lower precisions 100 x epsilon ...'
                exit
            end if
            if (Iteration > 250 .and. delta < 10000.0_dp * epsilon) then
                write(*,'(/,a)') 'OBDM iteration finished with lower precisions 1000 x epsilon ...'
                exit
            end if

            write(*,'(a)') TAB // 'OBDM iteration:   Iteration Number = ' // istr(Iteration) // ',' // TAB // '   delta = ' // jstr(delta)
            pRho_old = pRho
            nRho_old = nRho
        end do

        write(*,'(/,a)') 'HF-RRPA OBDM iteration has finished ...'
        Rho_out%p = pRho
        Rho_out%n = nRho
        C_out%p = pC_RRPA
        C_out%n = nC_RRPA

    contains

        ! Particle-particle & hole-hole pairs (a >= b, same l & j) ... types 1..4 = p particles, p holes, n particles, n holes ...
        function iteration_list() result(Pr)
            type(obdm_iteration_pairs_t) :: Pr
            integer :: Pass, k
            do Pass = 1, 2
                k = 0
                call add_pairs(Particle%p,1,Pass,k,Pr)
                call add_pairs(Hole%p,2,Pass,k,Pr)
                call add_pairs(Particle%n,3,Pass,k,Pr)
                call add_pairs(Hole%n,4,Pass,k,Pr)
                if (Pass == 1) then
                    Pr%N = k
                    allocate(Pr%T(k), Pr%a(k), Pr%b(k))
                end if
            end do
        end function iteration_list

        subroutine add_pairs(S,Typ,Pass,k,Pr)
            type(spstate_t), intent(in) :: S(:)
            integer, intent(in) :: Typ, Pass
            integer, intent(inout) :: k
            type(obdm_iteration_pairs_t), intent(inout) :: Pr
            integer :: x, y
            do x = 1, size(S)
                do y = 1, size(S)
                    if (S(x)%l == S(y)%l .and. S(x)%j == S(y)%j .and. S(x)%a >= S(y)%a) then
                        k = k + 1
                        if (Pass == 2) then
                            Pr%T(k) = Typ
                            Pr%a(k) = x
                            Pr%b(k) = y
                        end if
                    end if
                end do
            end do
        end subroutine add_pairs

        ! New OBDM from the RRPA amplitudes ... HF_RRPA_OBDM ...
        subroutine obdm(pR,nR)
            real(dp), intent(out) :: pR(:,:), nR(:,:)
            integer :: JJ, PP, N_ph, nu, mu, q, g, qg, Index, Typ, x, y, a_x, a_y, j_x, l_x, h, a_h, l_h, j_h, J_l, P_l
            integer :: ph, yh
            complex(dp) :: pS, nS
            real(dp) :: Sum, Amp_x, Amp_y

            pR = Rho_HF%p
            nR = Rho_HF%n

            ! Inner loop matrices M ...
            !$omp parallel do collapse(2) schedule(dynamic) private(N_ph,nu,mu,pS,nS,q,g,qg)
            do JJ = 0, J_max
                do PP = 1, 2
                    N_ph = N_nu(JJ,PP)
                    do nu = 1, N_ph
                        do mu = 1, N_ph
                            pS = (0.0_dp, 0.0_dp)
                            nS = (0.0_dp, 0.0_dp)
                            do q = 1, N_Particle%p
                                do g = 1, N_Hole%p
                                    if ((PP == (mod(Particle%p(q)%l + Hole%p(g)%l,2) + 1)) .and. &
                                        (abs(Particle%p(q)%j - Hole%p(g)%j) <= 2*JJ) .and. (2*JJ <= (Particle%p(q)%j + Hole%p(g)%j))) then
                                        qg = ParticleHole%p(JJ,PP,q,g)
                                        pS = pS + real(-0.5_dp * (pRho_old(Hole%p(g)%a,Hole%p(g)%a) - pRho_old(Particle%p(q)%a, &
                                             Particle%p(q)%a)),dp) * X_RPA(JJ,PP)%v(qg,mu) * conjg(X_RPA(JJ,PP)%v(qg,nu))
                                    end if
                                end do
                            end do
                            do q = 1, N_Particle%n
                                do g = 1, N_Hole%n
                                    if ((PP == (mod(Particle%n(q)%l + Hole%n(g)%l,2) + 1)) .and. &
                                        (abs(Particle%n(q)%j - Hole%n(g)%j) <= 2*JJ) .and. (2*JJ <= (Particle%n(q)%j + Hole%n(g)%j))) then
                                        qg = ParticleHole%n(JJ,PP,q,g)
                                        nS = nS + real(-0.5_dp * (nRho_old(Hole%n(g)%a,Hole%n(g)%a) - nRho_old(Particle%n(q)%a, &
                                             Particle%n(q)%a)),dp) * X_RPA(JJ,PP)%v(qg,mu) * conjg(X_RPA(JJ,PP)%v(qg,nu))
                                    end if
                                end do
                            end do
                            pM(JJ,PP)%v(nu,mu) = pS
                            nM(JJ,PP)%v(nu,mu) = nS
                        end do
                    end do
                end do
            end do
            !$omp end parallel do

            ! Pair loop ... each pair writes its own (a_x,a_y) & (a_y,a_x) entries ...
            !$omp parallel do schedule(dynamic) private(Typ,x,y,a_x,a_y,l_x,j_x,Sum,h,a_h,l_h,j_h,Amp_x,Amp_y,P_l,J_l,N_ph,ph,yh,nu,mu)
            do Index = 1, Pairs%N
                Typ = Pairs%T(Index)
                x = Pairs%a(Index)
                y = Pairs%b(Index)
                Sum = 0.0_dp
                select case (Typ)
                case (1, 3)
                    ! Particles ... sum over holes ...
                    if (Typ == 1) then
                        a_x = Particle%p(x)%a; l_x = Particle%p(x)%l; j_x = Particle%p(x)%j
                        a_y = Particle%p(y)%a
                    else
                        a_x = Particle%n(x)%a; l_x = Particle%n(x)%l; j_x = Particle%n(x)%j
                        a_y = Particle%n(y)%a
                    end if
                    do h = 1, merge(N_Hole%p, N_Hole%n, Typ == 1)
                        if (Typ == 1) then
                            a_h = Hole%p(h)%a; l_h = Hole%p(h)%l; j_h = Hole%p(h)%j
                            Amp_x = pRho_old(a_h,a_h) - pRho_old(a_x,a_x)
                            Amp_y = pRho_old(a_h,a_h) - pRho_old(a_y,a_y)
                        else
                            a_h = Hole%n(h)%a; l_h = Hole%n(h)%l; j_h = Hole%n(h)%j
                            Amp_x = nRho_old(a_h,a_h) - nRho_old(a_x,a_x)
                            Amp_y = nRho_old(a_h,a_h) - nRho_old(a_y,a_y)
                        end if
                        P_l = mod(l_x + l_h,2) + 1
                        do J_l = abs(j_x - j_h)/2, (j_x + j_h)/2
                            N_ph = N_nu(J_l,P_l)
                            if (Typ == 1) then
                                ph = ParticleHole%p(J_l,P_l,x,h)
                                yh = ParticleHole%p(J_l,P_l,y,h)
                            else
                                ph = ParticleHole%n(J_l,P_l,x,h)
                                yh = ParticleHole%n(J_l,P_l,y,h)
                            end if
                            call pair_sum(merge(1,2,Typ == 1),J_l,P_l,N_ph,ph,yh,Amp_x,Amp_y,a_x > a_y,.false.,Sum)
                        end do
                    end do
                    Sum = Sum / real(j_x + 1,dp)
                    if (Typ == 1) then
                        pR(a_x,a_y) = pR(a_x,a_y) + Sum
                        if (a_x > a_y) pR(a_y,a_x) = pR(a_y,a_x) + Sum
                    else
                        nR(a_x,a_y) = nR(a_x,a_y) + Sum
                        if (a_x > a_y) nR(a_y,a_x) = nR(a_y,a_x) + Sum
                    end if
                case default
                    ! Holes ... sum over particles ...
                    if (Typ == 2) then
                        a_x = Hole%p(x)%a; l_x = Hole%p(x)%l; j_x = Hole%p(x)%j
                        a_y = Hole%p(y)%a
                    else
                        a_x = Hole%n(x)%a; l_x = Hole%n(x)%l; j_x = Hole%n(x)%j
                        a_y = Hole%n(y)%a
                    end if
                    do h = 1, merge(N_Particle%p, N_Particle%n, Typ == 2)
                        if (Typ == 2) then
                            a_h = Particle%p(h)%a; l_h = Particle%p(h)%l; j_h = Particle%p(h)%j
                            Amp_x = pRho_old(a_x,a_x) - pRho_old(a_h,a_h)
                            Amp_y = pRho_old(a_y,a_y) - pRho_old(a_h,a_h)
                        else
                            a_h = Particle%n(h)%a; l_h = Particle%n(h)%l; j_h = Particle%n(h)%j
                            Amp_x = nRho_old(a_x,a_x) - nRho_old(a_h,a_h)
                            Amp_y = nRho_old(a_y,a_y) - nRho_old(a_h,a_h)
                        end if
                        P_l = mod(l_h + l_x,2) + 1
                        do J_l = abs(j_h - j_x)/2, (j_h + j_x)/2
                            N_ph = N_nu(J_l,P_l)
                            if (Typ == 2) then
                                ph = ParticleHole%p(J_l,P_l,h,x)
                                yh = ParticleHole%p(J_l,P_l,h,y)
                            else
                                ph = ParticleHole%n(J_l,P_l,h,x)
                                yh = ParticleHole%n(J_l,P_l,h,y)
                            end if
                            ! The neutron-hole branch of the original uses Y(pg,mu) in the mirrored term ... kept as is ...
                            call pair_sum(merge(1,2,Typ == 2),J_l,P_l,N_ph,ph,yh,Amp_x,Amp_y,a_x > a_y,Typ == 4,Sum)
                        end do
                    end do
                    Sum = Sum / real(j_x + 1,dp)
                    if (Typ == 2) then
                        pR(a_x,a_y) = pR(a_x,a_y) - Sum
                        if (a_x > a_y) pR(a_y,a_x) = pR(a_y,a_x) - Sum
                    else
                        nR(a_x,a_y) = nR(a_x,a_y) - Sum
                        if (a_x > a_y) nR(a_y,a_x) = nR(a_y,a_x) - Sum
                    end if
                end select
            end do
            !$omp end parallel do
        end subroutine obdm

        ! Inner (nu,mu) sums of a pair ... Iso = 1 protons (pM), 2 neutrons (nM) ...
        subroutine pair_sum(Iso,J_l,P_l,N_ph,ph,yh,Amp_x,Amp_y,Mirror,Quirk,Sum)
            integer, intent(in) :: Iso, J_l, P_l, N_ph, ph, yh
            real(dp), intent(in) :: Amp_x, Amp_y
            logical, intent(in) :: Mirror, Quirk
            real(dp), intent(inout) :: Sum
            integer :: nu, mu
            complex(dp) :: Y_ph_nu, Y_yh_nu, Y_yh_mu, Y_ph_mu, Mnm, Mmn
            real(dp) :: sq
            sq = sqrt(Amp_x * Amp_y)
            do nu = 1, N_ph
                Y_ph_nu = Y_RPA(J_l,P_l)%v(ph,nu)
                Y_yh_nu = Y_RPA(J_l,P_l)%v(yh,nu)
                do mu = 1, N_ph
                    Y_yh_mu = Y_RPA(J_l,P_l)%v(yh,mu)
                    if (Iso == 1) then
                        Mnm = pM(J_l,P_l)%v(nu,mu)
                        Mmn = pM(J_l,P_l)%v(mu,nu)
                    else
                        Mnm = nM(J_l,P_l)%v(nu,mu)
                        Mmn = nM(J_l,P_l)%v(mu,nu)
                    end if
                    Sum = Sum + real(2*J_l + 1,dp) * (sq * real(Mnm * Y_ph_nu * conjg(Y_yh_mu)))
                    if (Mirror) then
                        if (Quirk) then
                            Y_ph_mu = Y_RPA(J_l,P_l)%v(yh,mu)
                        else
                            Y_ph_mu = Y_RPA(J_l,P_l)%v(ph,mu)
                        end if
                        Sum = Sum + real(2*J_l + 1,dp) * (sq * real(Mmn * Y_ph_mu * conjg(Y_yh_nu)))
                    end if
                end do
                Sum = Sum + real(2*J_l + 1,dp) * (sq * real(Y_ph_nu * conjg(Y_yh_nu)))
            end do
        end subroutine pair_sum

        ! Hungarian reordering of natural orbitals onto the previous basis & sign fixing ...
        subroutine reorder(Nv,Cm)
            real(dp), intent(inout) :: Nv(:), Cm(:,:)
            integer :: Assign(size(Nv)), k
            call hungarian(1.0_dp - abs(Cm),Assign)
            Cm = Cm(:,Assign)
            Nv = Nv(Assign)
            do k = 1, size(Nv)
                if (Cm(k,k) < 0.0_dp) Cm(:,k) = -Cm(:,k)
            end do
        end subroutine reorder

        logical function occupation_check()
            integer :: p, h
            occupation_check = .false.
            do p = 1, N_Particle%p
                do h = 1, N_Hole%p
                    if (pN(Hole%p(h)%a) - pN(Particle%p(p)%a) < 0.0_dp) return
                end do
            end do
            do p = 1, N_Particle%n
                do h = 1, N_Hole%n
                    if (nN(Hole%n(h)%a) - nN(Particle%n(p)%a) < 0.0_dp) return
                end do
            end do
            occupation_check = .true.
        end function occupation_check

        ! Amplitudes in the new natural-orbital basis ... X_ph = sum_qg C_qp C_gh X_qg ...
        subroutine transform_amplitudes(Cp,Cn)
            real(dp), intent(in) :: Cp(:,:), Cn(:,:)
            integer :: JJ, PP, N_ph, nu, p, h, q, g, ph, qg
            complex(dp), allocatable :: Xn(:,:), Yn(:,:)
            complex(dp) :: XS, YS
            !$omp parallel do collapse(2) schedule(dynamic) private(N_ph,Xn,Yn,nu,h,p,ph,XS,YS,g,q,qg)
            do JJ = 0, J_max
                do PP = 1, 2
                    N_ph = N_nu(JJ,PP)
                    allocate(Xn(N_ph,N_ph), Yn(N_ph,N_ph))
                    Xn = (0.0_dp, 0.0_dp)
                    Yn = (0.0_dp, 0.0_dp)
                    do nu = 1, N_ph
                        do h = 1, N_Hole%p
                            do p = 1, N_Particle%p
                                if (.not. ((mod(Particle%p(p)%l + Hole%p(h)%l,2) + 1) == PP .and. abs(Particle%p(p)%j - &
                                    Hole%p(h)%j)/2 <= JJ .and. JJ <= (Particle%p(p)%j + Hole%p(h)%j)/2)) cycle
                                ph = ParticleHole%p(JJ,PP,p,h)
                                XS = (0.0_dp, 0.0_dp)
                                YS = (0.0_dp, 0.0_dp)
                                do g = 1, N_Hole%p
                                    if (.not. (Hole%p(g)%l == Hole%p(h)%l .and. Hole%p(g)%j == Hole%p(h)%j)) cycle
                                    do q = 1, N_Particle%p
                                        if (.not. (Particle%p(q)%l == Particle%p(p)%l .and. Particle%p(q)%j == Particle%p(p)%j)) cycle
                                        qg = ParticleHole%p(JJ,PP,q,g)
                                        YS = YS + cmplx(Cp(Particle%p(q)%a,Particle%p(p)%a) * Cp(Hole%p(g)%a,Hole%p(h)%a),0.0_dp,kind=dp) * &
                                             Y_RPA(JJ,PP)%v(qg,nu)
                                        XS = XS + cmplx(Cp(Particle%p(q)%a,Particle%p(p)%a) * Cp(Hole%p(g)%a,Hole%p(h)%a),0.0_dp,kind=dp) * &
                                             X_RPA(JJ,PP)%v(qg,nu)
                                    end do
                                end do
                                Xn(ph,nu) = Xn(ph,nu) + XS
                                Yn(ph,nu) = Yn(ph,nu) + YS
                            end do
                        end do
                        do h = 1, N_Hole%n
                            do p = 1, N_Particle%n
                                if (.not. ((mod(Particle%n(p)%l + Hole%n(h)%l,2) + 1) == PP .and. abs(Particle%n(p)%j - &
                                    Hole%n(h)%j)/2 <= JJ .and. JJ <= (Particle%n(p)%j + Hole%n(h)%j)/2)) cycle
                                ph = ParticleHole%n(JJ,PP,p,h)
                                XS = (0.0_dp, 0.0_dp)
                                YS = (0.0_dp, 0.0_dp)
                                do g = 1, N_Hole%n
                                    if (.not. (Hole%n(g)%l == Hole%n(h)%l .and. Hole%n(g)%j == Hole%n(h)%j)) cycle
                                    do q = 1, N_Particle%n
                                        if (.not. (Particle%n(q)%l == Particle%n(p)%l .and. Particle%n(q)%j == Particle%n(p)%j)) cycle
                                        qg = ParticleHole%n(JJ,PP,q,g)
                                        YS = YS + cmplx(Cn(Particle%n(q)%a,Particle%n(p)%a) * Cn(Hole%n(g)%a,Hole%n(h)%a),0.0_dp,kind=dp) * &
                                             Y_RPA(JJ,PP)%v(qg,nu)
                                        XS = XS + cmplx(Cn(Particle%n(q)%a,Particle%n(p)%a) * Cn(Hole%n(g)%a,Hole%n(h)%a),0.0_dp,kind=dp) * &
                                             X_RPA(JJ,PP)%v(qg,nu)
                                    end do
                                end do
                                Xn(ph,nu) = Xn(ph,nu) + XS
                                Yn(ph,nu) = Yn(ph,nu) + YS
                            end do
                        end do
                    end do
                    X_RPA(JJ,PP)%v = Xn
                    Y_RPA(JJ,PP)%v = Yn
                    deallocate(Xn, Yn)
                end do
            end do
            !$omp end parallel do
        end subroutine transform_amplitudes

    end subroutine HF_RRPA_OBDM_iteration

    ! RRPA correlation energy ... sum over p-h components ...
    function HF_RRPA_energy(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA,Quiet) result(E_0)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: Y_RPA(0:,:)
        logical, intent(in), optional :: Quiet
        real(dp) :: E_0
        logical :: Loud

        Loud = .true.
        if (present(Quiet)) Loud = .not. Quiet
        if (Loud) write(*,'(/,a)') 'Calculating RRPA ground-state correlation energy ...'
        E_0 = bosonic_sum(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)
        if (Loud) then
            write(*,'(a)') TAB // 'ERPA ground-state correlation energy E_0^RRPA is ...     E_0^RRPA = ' // jstr(round_digits(E_0,6)) // TAB // 'MeV'
            write(*,'(a)') TAB // 'tTo get the total ground-state energy add the mean-field energy ...'
        end if
    end function HF_RRPA_energy

    subroutine HF_RRPA_energy_bosonic(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: Y_RPA(0:,:)
        real(dp) :: E_0
        write(*,'(/,a)') 'Calculating RRPA correlation energy ... Using the bosonic algebra formula ...'
        E_0 = bosonic_sum(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA)
        write(*,'(a)') TAB // 'RRPA ground-state correlation energy E_0^RRPA is ...     E_0 = ' // jstr(round_digits(E_0,6)) // TAB // 'MeV'
        write(*,'(a)') TAB // TAB // 'To get the total ground-state energy add the mean-field energy ...'
    end subroutine HF_RRPA_energy_bosonic

    real(dp) function bosonic_sum(Params,N_nu,N_Particle,Particle,N_Hole,Hole,E_RPA,Y_RPA) result(E_0)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: Y_RPA(0:,:)
        type(pnarray4_t) :: ParticleHole
        integer :: J_max, h, p, J, P_ph, nu, ph
        J_max = 2*Params%Calc%Nmax + 1
        ParticleHole = orbitals_ph_list(J_max,N_Particle,Particle,N_Hole,Hole)
        E_0 = 0.0_dp
        do h = 1, N_Hole%p
            do p = 1, N_Particle%p
                P_ph = mod(Particle%p(p)%l + Hole%p(h)%l,2) + 1
                do J = abs(Particle%p(p)%j - Hole%p(h)%j)/2, (Particle%p(p)%j + Hole%p(h)%j)/2
                    ph = ParticleHole%p(J,P_ph,p,h)
                    do nu = 1, N_nu(J,P_ph)
                        E_0 = E_0 - real(2*J + 1,dp) * real(E_RPA(J,P_ph)%v(nu)) * abs(Y_RPA(J,P_ph)%v(ph,nu))**2
                    end do
                end do
            end do
        end do
        do h = 1, N_Hole%n
            do p = 1, N_Particle%n
                P_ph = mod(Particle%n(p)%l + Hole%n(h)%l,2) + 1
                do J = abs(Particle%n(p)%j - Hole%n(h)%j)/2, (Particle%n(p)%j + Hole%n(h)%j)/2
                    ph = ParticleHole%n(J,P_ph,p,h)
                    do nu = 1, N_nu(J,P_ph)
                        E_0 = E_0 - real(2*J + 1,dp) * real(E_RPA(J,P_ph)%v(nu)) * abs(Y_RPA(J,P_ph)%v(ph,nu))**2
                    end do
                end do
            end do
        end do
    end function bosonic_sum

    ! Alternative 2-body density based correlation energy (printed only, as in the original) ...
    subroutine HF_RRPA_energy_density(Params,N_nu,N_Particle,Particle,N_Hole,Hole,A,E_RPA,Y_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(rmat_t), intent(in) :: A(0:,:)
        type(cvec_t), intent(in) :: E_RPA(0:,:)
        type(cmat_t), intent(in) :: Y_RPA(0:,:)
        type(pnarray4_t) :: ParticleHole
        complex(dp) :: E_0
        integer :: J_max

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Calculating RRPA correlation energy ... Using the alternative 2-body density based formula ...'
        ParticleHole = orbitals_ph_list(J_max,N_Particle,Particle,N_Hole,Hole)
        E_0 = (0.0_dp, 0.0_dp)
        call accumulate(Particle%p,Hole%p,ParticleHole%p)
        call accumulate(Particle%n,Hole%n,ParticleHole%n)
        write(*,'(a)') TAB // 'RRPA ground-state correlation energy E_0^RRPA is ...     E_0^RRPA = ' // &
                       jstr(round_digits(real(E_0),6)) // TAB // 'MeV'
        write(*,'(a)') TAB // TAB // 'To get the total ground-state energy add the mean-field energy ...'

    contains

        subroutine accumulate(Part,Hol,PH_tab)
            type(spstate_t), intent(in) :: Part(:), Hol(:)
            integer, intent(in) :: PH_tab(0:,:,:,:)
            integer :: h, p, g, q, J, P_ph, nu, ph, qg
            real(dp) :: Amp_J
            do h = 1, size(Hol)
                do p = 1, size(Part)
                    P_ph = mod(Part(p)%l + Hol(h)%l,2) + 1
                    do J = abs(Part(p)%j - Hol(h)%j)/2, (Part(p)%j + Hol(h)%j)/2
                        ph = PH_tab(J,P_ph,p,h)
                        Amp_J = 0.5_dp * real(2*J + 1,dp)
                        do nu = 1, N_nu(J,P_ph)
                            E_0 = E_0 - Amp_J * real(E_RPA(J,P_ph)%v(nu)) * abs(Y_RPA(J,P_ph)%v(ph,nu))**2
                            do g = 1, size(Hol)
                                do q = 1, size(Part)
                                    if ((mod(Part(p)%l + Hol(h)%l,2) == mod(Part(q)%l + Hol(g)%l,2)) .and. &
                                        (abs(Part(q)%j - Hol(g)%j) <= 2*J) .and. (2*J <= (Part(q)%j + Hol(g)%j))) then
                                        qg = PH_tab(J,P_ph,q,g)
                                        E_0 = E_0 - Amp_J * A(J,P_ph)%v(ph,qg) * Y_RPA(J,P_ph)%v(ph,nu) * Y_RPA(J,P_ph)%v(qg,nu)
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end subroutine accumulate

    end subroutine HF_RRPA_energy_density

    ! Reduced multipole matrix elements of RRPA phonons ... occupation-weighted p-h components ...
    subroutine HF_RRPA_rM(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Rho,TrOp,X_RPA,Y_RPA,rM)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        type(ivec_t), intent(in) :: Orb_Phonon(0:,:)
        type(phstate_t), intent(in) :: Phonon(:)
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(o1b_t), intent(in) :: Rho
        type(tr1b_t), intent(in) :: TrOp
        type(cmat_t), intent(in) :: X_RPA(0:,:), Y_RPA(0:,:)
        type(reduced_multipole_t), intent(out) :: rM

        write(*,'(/,a)') 'Calculating reduced transition matrix elements rM ...'
        call multipole(0,1,1,TrOp%E0,rM%E0)
        call multipole(1,2,2,TrOp%E1,rM%E1)
        call multipole(1,2,3,TrOp%E1_VC,rM%E1_VC)
        call multipole(1,2,3,TrOp%E1_VS,rM%E1_VS)
        call multipole(1,2,3,TrOp%E1_TC,rM%E1_TC)
        call multipole(1,2,3,TrOp%E1_TS,rM%E1_TS)
        call multipole(1,2,3,TrOp%E1_sTC,rM%E1_sTC)
        call multipole(1,2,3,TrOp%E1_sTS,rM%E1_sTS)
        call multipole(1,2,3,TrOp%E1_C,rM%E1_C)
        call multipole(1,2,3,TrOp%E1_sC,rM%E1_sC)
        call multipole(2,1,1,TrOp%E2,rM%E2)
        call multipole(3,2,1,TrOp%E3,rM%E3)
        write(*,'(a)') TAB // 'Reduced 1-body transition matrix elements rM calculated ...'

    contains

        ! Mode 1: O w (X + Y) ... mode 2 (E1): O w (-X + Y) ... mode 3 (other E1 operators): -O w (-X + Y) ...
        subroutine multipole(J,P,Mode,O,M)
            integer, intent(in) :: J, P, Mode
            type(o1b_t), intent(in) :: O
            type(pncvector_t), intent(out) :: M
            integer :: N_ph, nu, ph, i_ph, a_p, a_h
            complex(dp) :: pS, nS, Amp
            real(dp) :: Op
            N_ph = N_nu(J,P)
            allocate(M%p(N_ph), M%n(N_ph))
            !$omp parallel do private(pS,nS,ph,i_ph,a_p,a_h,Op,Amp)
            do nu = 1, N_ph
                pS = (0.0_dp, 0.0_dp)
                nS = (0.0_dp, 0.0_dp)
                do ph = 1, N_ph
                    i_ph = Orb_Phonon(J,P)%v(ph)
                    call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                    if (Phonon(i_ph)%tz == -1) then
                        Op = O%p(a_p,a_h) * sqrt(Rho%p(a_h,a_h) - Rho%p(a_p,a_p))
                    else
                        Op = O%n(a_p,a_h) * sqrt(Rho%n(a_h,a_h) - Rho%n(a_p,a_p))
                    end if
                    select case (Mode)
                    case (1)
                        Amp = cmplx(Op,0.0_dp,kind=dp) * (X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                    case (2)
                        Amp = cmplx(Op,0.0_dp,kind=dp) * (-1.0_dp * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                    case default
                        Amp = cmplx(-1.0_dp * Op,0.0_dp,kind=dp) * (-1.0_dp * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu))
                    end select
                    if (Phonon(i_ph)%tz == -1) then
                        pS = pS + Amp
                    else
                        nS = nS + Amp
                    end if
                end do
                M%p(nu) = pS
                M%n(nu) = nS
            end do
            !$omp end parallel do
        end subroutine multipole

    end subroutine HF_RRPA_rM

end module nj_hf_rrpa
