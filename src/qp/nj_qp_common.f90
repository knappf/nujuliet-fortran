! Common pieces of the BCS & HFB solvers ... mean-field + pairing energy of BCS_energy & HFB_energy ...
module nj_qp_common
    use nj_kinds
    use nj_types
    use nj_params
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_v3b, only : V3b_no2b
    implicit none
    private

    public :: pairing_energy

contains

    ! E_MF (mean-field, kinetic + NN + NNN) & pairing energies E_2N, E_3N ...
        ! Kappa_ba = .true. uses Kappa(b,a) as HFB_energy, .false. uses Kappa(a,b) as BCS_energy ...
    subroutine pairing_energy(Params,Rho,Kappa,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,E_MF,E_2N,E_3N,Kappa_ba)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, Kappa, T
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        real(dp), intent(out) :: E_MF, E_2N, E_3N
        logical, intent(in) :: Kappa_ba
        integer :: a_max, N_2max, N_3max, a, b, d, e, c, f, J, P_ab, P_abc
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_d, l_d, j_d, n_e, l_e, j_e, n_c, l_c, j_c, n_f, l_f, j_f
        real(dp) :: A_mass, p2, p3, Sum, Hat, pRho_ad, nRho_ad, pRho_be, nRho_be, pRho_cf, nRho_cf
        real(dp) :: ME001, ME101, ME011, ME111, ME113, j_a_hat, j_d_hat, pK_ab, nK_ab, pK_de, nK_de
        character(len=32) :: CMS

        A_mass = real(Params%Calc%A,dp)
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        p2 = Params%Int%cP2N
        p3 = Params%Int%cP3N
        CMS = Params%Calc%CMS

        ! Mean-field energy ...
        E_MF = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:E_MF) private(Sum,n_a,l_a,j_a,b,n_b,l_b,j_b,P_ab,d,l_d,j_d,n_d, &
        !$omp pRho_ad,nRho_ad,e,n_e,l_e,j_e,pRho_be,nRho_be,J,Hat,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f,pRho_cf,nRho_cf, &
        !$omp ME001,ME101,ME011,ME111,ME113)
        do a = 1, a_max
            Sum = 0.0_dp
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                    j_b = Orb(b)%j
                    P_ab = mod(l_a + l_b,2) + 1
                    do d = 1, a_max
                        l_d = Orb(d)%l
                        j_d = Orb(d)%j
                        if (.not. ((l_a == l_d) .and. (j_a == j_d))) cycle
                        n_d = Orb(d)%n
                        pRho_ad = Rho%p(a,d)
                        nRho_ad = Rho%n(a,d)
                        do e = 1, a_max
                            n_e = Orb(e)%n
                            l_e = Orb(e)%l
                            j_e = Orb(e)%j
                            if (.not. ((l_b == l_e .and. j_b == j_e) .and. ((2*(n_d + n_e) + l_d + l_e) <= N_2max))) cycle
                            pRho_be = Rho%p(b,e)
                            nRho_be = Rho%n(b,e)
                            do J = abs(j_a - j_b)/2, (j_a + j_b)/2
                                if (mod(l_a + l_b, 2) == mod(l_d + l_e, 2)) then
                                    Hat = real(2*J + 1,dp)
                                    Sum = Sum + 0.5_dp * Hat * pRho_ad * pRho_be * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                                    Sum = Sum + 0.5_dp * Hat * nRho_ad * nRho_be * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                                    Sum = Sum + Hat * pRho_ad * nRho_be * O2b_pn(a,b,d,e,J,P_ab,V_NN,Orb_NN)
                                end if
                                do c = 1, a_max
                                    n_c = Orb(c)%n
                                    l_c = Orb(c)%l
                                    if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                                    j_c = Orb(c)%j
                                    P_abc = mod(l_a + l_b + l_c, 2) + 1
                                    do f = 1, a_max
                                        n_f = Orb(f)%n
                                        l_f = Orb(f)%l
                                        j_f = Orb(f)%j
                                        if (.not. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                                   j_c == j_f .and. P_abc == (mod(l_d + l_e + l_f, 2) + 1))) cycle
                                        pRho_cf = Rho%p(c,f)
                                        nRho_cf = Rho%n(c,f)
                                        ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                        ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P_abc,V_NNN,Orb,Orb_NNN)
                                        Sum = Sum + 1.0_dp/6.0_dp * (ME113 * pRho_ad * pRho_be * pRho_cf + (2.0_dp * ME111 + ME113) * &
                                              pRho_ad * pRho_be * nRho_cf + (1.5_dp * ME001 + sqrt(3.0_dp/4.0_dp) * ME101 + &
                                              sqrt(3.0_dp/4.0_dp) * ME011 + 0.5_dp * ME111 + ME113) * pRho_ad * nRho_be * nRho_cf + &
                                              ME113 * nRho_ad * nRho_be * nRho_cf)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
                if (CMS == 'CMS1+2B') then
                    Sum = Sum + (Rho%p(a,b) + Rho%n(a,b)) * T%p(a,b) * real(j_a + 1,dp) * (1.0_dp - 1.0_dp / A_mass)
                else if (CMS /= 'CMS2B') then
                    Sum = Sum + (Rho%p(a,b) + Rho%n(a,b)) * T%p(a,b) * real(j_a + 1,dp)
                end if
            end do
            E_MF = E_MF + Sum
        end do
        !$omp end parallel do

        ! Pairing energy ...
        E_2N = 0.0_dp
        E_3N = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:E_2N,E_3N) private(n_a,l_a,j_a,j_a_hat,b,n_b,l_b,j_b,pK_ab,nK_ab, &
        !$omp d,n_d,l_d,j_d,j_d_hat,e,n_e,l_e,j_e,pK_de,nK_de,c,n_c,l_c,j_c,P_abc,f,n_f,l_f,j_f,pRho_cf,nRho_cf,ME111,ME113)
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            j_a_hat = sqrt(real(j_a + 1,dp))
            do b = 1, a_max
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                if (.not. ((2*(n_a + n_b) + l_a + l_b) <= N_2max .and. l_a == l_b)) cycle
                j_b = Orb(b)%j
                if (j_a /= j_b) cycle
                if (Kappa_ba) then
                    pK_ab = Kappa%p(b,a)
                    nK_ab = Kappa%n(b,a)
                else
                    pK_ab = Kappa%p(a,b)
                    nK_ab = Kappa%n(a,b)
                end if
                do d = 1, a_max
                    n_d = Orb(d)%n
                    l_d = Orb(d)%l
                    j_d = Orb(d)%j
                    j_d_hat = sqrt(real(j_d + 1,dp))
                    do e = 1, a_max
                        n_e = Orb(e)%n
                        l_e = Orb(e)%l
                        if (.not. (((2*(n_d + n_e) + l_d + l_e) <= N_2max) .and. l_d == l_e)) cycle
                        j_e = Orb(e)%j
                        if (j_d /= j_e) cycle
                        pK_de = Kappa%p(d,e)
                        nK_de = Kappa%n(d,e)
                        E_2N = E_2N + p2 * 0.25_dp * j_a_hat * j_d_hat * pK_ab * pK_de * O2b_pp(a,b,d,e,0,1,V_NN,Orb,Orb_NN)
                        E_2N = E_2N + p2 * 0.25_dp * j_a_hat * j_d_hat * nK_ab * nK_de * O2b_nn(a,b,d,e,0,1,V_NN,Orb,Orb_NN)
                        do c = 1, a_max
                            n_c = Orb(c)%n
                            l_c = Orb(c)%l
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) > N_3max) cycle
                            j_c = Orb(c)%j
                            P_abc = mod(l_a + l_b + l_c, 2) + 1
                            do f = 1, a_max
                                n_f = Orb(f)%n
                                l_f = Orb(f)%l
                                j_f = Orb(f)%j
                                if (.not. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. l_c == l_f .and. &
                                           j_c == j_f .and. P_abc == (mod(l_d + l_e + l_f, 2) + 1))) cycle
                                pRho_cf = Rho%p(c,f)
                                nRho_cf = Rho%n(c,f)
                                ME111 = V3b_no2b(a,b,c,1,d,e,f,1,0,1,P_abc,V_NNN,Orb,Orb_NNN)
                                ME113 = V3b_no2b(a,b,c,1,d,e,f,1,0,3,P_abc,V_NNN,Orb,Orb_NNN)
                                E_3N = E_3N + p3 * 0.25_dp * j_a_hat * j_d_hat * (pK_ab * pK_de * pRho_cf * ME113 + &
                                       pK_ab * pK_de * nRho_cf * 1.0_dp / 3.0_dp * (2.0_dp * ME111 + ME113) + &
                                       nK_ab * nK_de * nRho_cf * ME113 + nK_ab * nK_de * pRho_cf * 1.0_dp / 3.0_dp * (2.0_dp * ME111 + ME113))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine pairing_energy

end module nj_qp_common
