! Density-dependent mean-field & residual Hamiltonians from the bare V_NN & the NO2B V_NNN ... counterpart of
    ! src/common/Hamiltonian.jl (NuJuliet dev, new there) ...
    !   H1b      ... 1-body mean field h_N[Rho] (as HF_allocate of the HF solver),
    !   H2b      ... density-dependent residual NN interaction H_NN[Rho] (NO2B approximation of V_NNN; the fast equivalent of
    !                V2b_residual_no2b, which it equals up to the summation order),
    !   H1b_H2b  ... H2b & the Sigma-dressed mean field h_N[Rho,Sigma] (the V_NNN-Sigma contraction; the V_NN-Sigma term is
    !                commented out in dev & is not ported) ...
    ! Orb are the calculation-space orbitals (only the first a_max(Calc.Nmax) entries are used, so a longer Int.Nmax list of the
    ! legacy solvers works too), Rho is the LHO-basis OBDM (only its lj-diagonal blocks enter), V_NN & V_NNN are bare ...
    ! Not ported (unreachable in dev): the 9-argument H1b(...,Rho,Sigma) with the V_NN-Sigma term (H1b_dense), E1b & E2b ...
    ! The V_NNN blocks are indexed directly with the row-wise packing of V3b_no2b_index (Ind = Ket + Bra(Bra-1)/2) ...
module nj_hamiltonian
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_o1b, only : o1b_alloc
    use nj_o2b, only : O2b_index, O2b_pp, O2b_pn, O2b_nn, packed_index
    use nj_v3b, only : V3b_no2b, V3b_jindex
    implicit none
    private

    public :: H1b, H2b, H1b_H2b

contains

    ! Index of (Bra|V_NNN|Ket) in the stored lower triangle of a V_NNN block ... H1b_V3b_index / H2b_V3b_index ...
    pure integer(i8) function v3b_idx(Bra,Ket)
        integer, intent(in) :: Bra, Ket
        if (Bra >= Ket) then
            v3b_idx = int(Ket,i8) + (int(Bra,i8) * int(Bra - 1,i8))/2_i8
        else
            v3b_idx = int(Bra,i8) + (int(Ket,i8) * int(Ket - 1,i8))/2_i8
        end if
    end function v3b_idx

    ! All ordered pairs ab (ab_all) & pairs with b <= a (ab_low) with l_a = l_b & j_a = j_b, a outer & b inner ... H1b_lj_pairs ...
    subroutine lj_pairs(a_max,Orb,ab_all,N_all,ab_low,N_low)
        integer, intent(in) :: a_max
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable, intent(out) :: ab_all(:,:), ab_low(:,:)
        integer, intent(out) :: N_all, N_low
        integer :: a, b, Pass
        allocate(ab_all(2,0), ab_low(2,0))
        do Pass = 1, 2
            N_all = 0
            N_low = 0
            do a = 1, a_max
                do b = 1, a_max
                    if (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j) then
                        N_all = N_all + 1
                        if (Pass == 2) ab_all(:,N_all) = [a, b]
                        if (b <= a) then
                            N_low = N_low + 1
                            if (Pass == 2) ab_low(:,N_low) = [a, b]
                        end if
                    end if
                end do
            end do
            if (Pass == 1) then
                deallocate(ab_all, ab_low)
                allocate(ab_all(2,N_all), ab_low(2,N_low))
            end if
        end do
    end subroutine lj_pairs

    ! 1-body mean field h_N[Rho] ... counterpart of H1b (8 arguments) ...
    function H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho) result(h_N)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: T, Rho
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t) :: h_N
        write(*,'(/,a)') 'Evaluating 1-body mean-field Hamiltonian h_N ...'
        call h1b_rho(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho,h_N)
        write(*,'(a)') TAB // 'Completed the calculation of 1-body mean-field Hamiltonian h_N ...'
    end function H1b

    ! Rho part of h_N (2-body V_NN & 3-body V_NNN contractions with Rho) + kinetic & CM terms ... shared by H1b & H1b_H2b ...
        ! Threads over the pairs ad (a >= d, lj-equal), each summing its be pairs in the order of dev with one task ...
    subroutine h1b_rho(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho,h_N)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: T, Rho
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t), intent(out) :: h_N
        integer, allocatable :: ad(:,:), cf(:,:)
        integer :: N_2max, N_3max, a_max, ad_count, cf_count, ad_i, be_i, cf_i, a, b, c, d, e, f, J, P_ab, P_abc
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_d, l_d, j_d, n_e, l_e, j_e, n_c, l_c, n_f, l_f
        real(dp) :: A_mass, is_j_a_hat, J_j_a_hat, pRho_be, nRho_be, pRho_cf, nRho_cf, pSum, nSum, pSum_local, nSum_local, sr13
        real(dp) :: ME001, ME101, ME011, ME111, ME113
        character(len=32) :: CMS

        A_mass = real(Params%Calc%A,dp)
        N_2max = Params%Calc%N2max
        N_3max = Params%Calc%N3max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        CMS = Params%Calc%CMS
        sr13 = sqrt(1.0_dp/3.0_dp)

        ! ad ... pairs a >= d, be & cf ... all ordered pairs (with equal l & j) ...
        call lj_pairs(a_max,Orb,cf,cf_count,ad,ad_count)
        call o1b_alloc(h_N,a_max)

        !$omp parallel do schedule(dynamic) private(a,d,b,e,c,f,be_i,cf_i,J,P_ab,P_abc,n_a,l_a,j_a,n_b,l_b,j_b,n_d,l_d,j_d, &
        !$omp n_e,l_e,j_e,n_c,l_c,n_f,l_f,is_j_a_hat,J_j_a_hat,pRho_be,nRho_be,pRho_cf,nRho_cf,pSum,nSum,pSum_local,nSum_local, &
        !$omp ME001,ME101,ME011,ME111,ME113)
        do ad_i = 1, ad_count
            a = ad(1,ad_i)
            d = ad(2,ad_i)
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            n_d = Orb(d)%n
            l_d = Orb(d)%l
            j_d = Orb(d)%j
            is_j_a_hat = 1.0_dp / real(j_a + 1,dp)
            pSum = 0.0_dp
            nSum = 0.0_dp
            do be_i = 1, cf_count
                b = cf(1,be_i)
                e = cf(2,be_i)
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                n_e = Orb(e)%n
                l_e = Orb(e)%l
                j_e = Orb(e)%j
                P_ab = parity_of(l_a + l_b)
                pSum_local = 0.0_dp
                nSum_local = 0.0_dp
                if ((2*(n_a + n_b) + l_a + l_b) <= N_2max .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max) then
                    pRho_be = Rho%p(b,e)
                    nRho_be = Rho%n(b,e)
                    do J = abs(j_d - j_e)/2, (j_d + j_e)/2
                        ! 2-body NN interaction part ...
                        if (P_ab == parity_of(l_d + l_e)) then
                            J_j_a_hat = real(2*J + 1,dp) * is_j_a_hat
                            pSum_local = pSum_local + J_j_a_hat * pRho_be * O2b_pp(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                            pSum_local = pSum_local + J_j_a_hat * nRho_be * O2b_pn(a,b,d,e,J,P_ab,V_NN,Orb_NN)
                            nSum_local = nSum_local + J_j_a_hat * nRho_be * O2b_nn(a,b,d,e,J,P_ab,V_NN,Orb,Orb_NN)
                            nSum_local = nSum_local + J_j_a_hat * pRho_be * O2b_pn(b,a,e,d,J,P_ab,V_NN,Orb_NN)
                        end if
                        ! 3-body NNN interaction part ... only cf with l_c = l_f & j_c = j_f ...
                        do cf_i = 1, cf_count
                            c = cf(1,cf_i)
                            f = cf(2,cf_i)
                            n_c = Orb(c)%n
                            l_c = Orb(c)%l
                            n_f = Orb(f)%n
                            l_f = Orb(f)%l
                            P_abc = parity_of(l_a + l_b + l_c)
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max .and. &
                                (2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max .and. P_abc == parity_of(l_d + l_e + l_f)) then
                                pRho_cf = Rho%p(c,f)
                                nRho_cf = Rho%n(c,f)
                                ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P_abc,V_NNN,Orb,Orb_NNN)
                                ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P_abc,V_NNN,Orb,Orb_NNN)
                                pSum_local = pSum_local + is_j_a_hat * (0.5_dp*ME113*pRho_be*pRho_cf + 0.25_dp * (ME001 + &
                                             sr13*ME101 + sr13*ME011 + ME111/3.0_dp + 2.0_dp/3.0_dp*ME113)*nRho_be*nRho_cf + &
                                             1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*pRho_be*nRho_cf)
                                nSum_local = nSum_local + is_j_a_hat * (0.5_dp*ME113*nRho_be*nRho_cf + 0.25_dp * (ME001 + &
                                             sr13*ME101 + sr13*ME011 + ME111/3.0_dp + 2.0_dp/3.0_dp*ME113)*pRho_be*pRho_cf + &
                                             1.0_dp/3.0_dp * (2.0_dp*ME111 + ME113)*nRho_be*pRho_cf)
                            end if
                        end do
                    end do
                end if
                pSum = pSum + pSum_local
                nSum = nSum + nSum_local
            end do

            ! Kinetic energy & CM correction ... CMS1+2B (1 - 1/A) T, no CM correction T, CMS2B nothing ...
            if (CMS == 'CMS1+2B') then
                pSum = pSum + T%p(a,d) * (1.0_dp - 1.0_dp / A_mass)
                nSum = nSum + T%n(a,d) * (1.0_dp - 1.0_dp / A_mass)
            else if (CMS /= 'CMS2B') then
                pSum = pSum + T%p(a,d)
                nSum = nSum + T%n(a,d)
            end if
            h_N%p(a,d) = pSum
            h_N%p(d,a) = pSum
            h_N%n(a,d) = nSum
            h_N%n(d,a) = nSum
        end do
        !$omp end parallel do
    end subroutine h1b_rho

    ! Density-dependent residual 2-body Hamiltonian H_NN[Rho] (NO2B approximation of V_NNN) ... counterpart of H2b ...
        ! Threads over the pp pairs Bra = (a >= b) of each (J,P) ... every Bra writes only its own elements (Ket <= Bra), the pn
        ! element of each ordering is reached from exactly one (Bra,Ket) ...
    subroutine H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho,H_NN,Report)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(out) :: H_NN
        logical, intent(in), optional :: Report
        integer, allocatable :: cf(:,:), ad(:,:), j_index(:)
        integer :: a_max, J_max, cf_count, ad_count, J, P, N_pp, Bra, Ket, a, b, c, d, e, f, cf_i, ab, de, J_ab, J_de, l_c, P3B, Ind
        integer :: abc_0, def_0, abc_1, def_1, abc_3, def_3, l_a, l_b
        real(dp) :: cRes, is_J_hat, Phase_ab, Phase_de, isr12, ppSum, nnSum, pnSum_abde, pnSum_bade, pnSum_abed, pnSum_baed
        real(dp) :: ME001, ME101, ME011, ME111, ME113, pRho_fc, nRho_fc
        logical :: Rep

        Rep = .true.
        if (present(Report)) Rep = Report
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        J_max = Params%Calc%N2max + 1
        cRes = Params%Int%cRes
        isr12 = 1.0_dp / sqrt(12.0_dp)

        if (Rep) write(*,'(/,a)') 'Calculating density dependent residual 2-body interaction Hamiltonian H_NN ... (NO2B approximation of V_NNN)'

        ! All ordered lj-equal pairs cf ... only these are connected by the lj-diagonal Rho ...
        call lj_pairs(a_max,Orb,cf,cf_count,ad,ad_count)
        allocate(j_index(a_max))
        do c = 1, a_max
            j_index(c) = V3b_jindex(Orb(c)%l,Orb(c)%j)
        end do

        ! Start from the 2-body interaction V_NN ...
        H_NN = V_NN

        do J = 0, J_max
            do P = 1, 2
                if (Rep) then
                    if (P == 1) then
                        write(*,'(a)') TAB // 'Calculating   ...   J = ' // istr(J) // '/' // istr(J_max) // TAB // 'P = +'
                    else
                        write(*,'(a)') TAB // 'Calculating   ...   J = ' // istr(J) // '/' // istr(J_max) // TAB // 'P = -'
                    end if
                end if
                N_pp = Orb_NN%N(2,P,J)
                is_J_hat = 1.0_dp / real(2*J + 1,dp)

                !$omp parallel do schedule(dynamic) private(a,b,l_a,l_b,Phase_ab,ab,J_ab,Ket,d,e,Phase_de,de,J_de, &
                !$omp ppSum,nnSum,pnSum_abde,pnSum_bade,pnSum_abed,pnSum_baed,cf_i,c,f,abc_0,def_0,abc_1,def_1,abc_3,def_3, &
                !$omp l_c,P3B,ME001,ME101,ME011,ME111,ME113,pRho_fc,nRho_fc,Ind)
                do Bra = 1, N_pp
                    a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                    b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                    l_a = Orb(a)%l
                    l_b = Orb(b)%l
                    ! Phase of V_NNN for (b a T_ab = 0| = (-1)^((j_a + j_b)/2 - J) (a b T_ab = 0| ... opposite for T_ab = 1 ...
                    Phase_ab = phase((Orb(a)%j + Orb(b)%j)/2 - J)
                    ab = a*(a - 1)/2 + b
                    J_ab = J - abs(Orb(a)%j - Orb(b)%j)/2 + 1
                    do Ket = 1, Bra
                        d = Orb_NN%Ind(2,P,J)%v(1,Ket)
                        e = Orb_NN%Ind(2,P,J)%v(2,Ket)
                        Phase_de = phase((Orb(d)%j + Orb(e)%j)/2 - J)
                        de = d*(d - 1)/2 + e
                        J_de = J - abs(Orb(d)%j - Orb(e)%j)/2 + 1

                        ppSum = 0.0_dp
                        nnSum = 0.0_dp
                        pnSum_abde = 0.0_dp
                        pnSum_bade = 0.0_dp
                        pnSum_abed = 0.0_dp
                        pnSum_baed = 0.0_dp
                        do cf_i = 1, cf_count
                            c = cf(1,cf_i)
                            f = cf(2,cf_i)
                            ! Index of (ab T_ab c| & |de T_de f) in V_NNN ... 0 if not allowed by N_3max ...
                            abc_0 = Orb_NNN%Map(ab)%m(1,J_ab,c)
                            def_0 = Orb_NNN%Map(de)%m(1,J_de,f)
                            if (abc_0 == 0 .or. def_0 == 0) cycle
                            abc_1 = Orb_NNN%Map(ab)%m(2,J_ab,c)
                            def_1 = Orb_NNN%Map(de)%m(2,J_de,f)
                            abc_3 = Orb_NNN%Map(ab)%m(3,J_ab,c)
                            def_3 = Orb_NNN%Map(de)%m(3,J_de,f)

                            l_c = Orb(c)%l
                            P3B = parity_of(l_a + l_b + l_c)
                            ME001 = real(V_NNN%V(P3B,l_c,J,1,j_index(c))%v(v3b_idx(abc_0,def_0)),dp)
                            ME101 = real(V_NNN%V(P3B,l_c,J,1,j_index(c))%v(v3b_idx(abc_1,def_0)),dp)
                            ME011 = real(V_NNN%V(P3B,l_c,J,1,j_index(c))%v(v3b_idx(abc_0,def_1)),dp)
                            ME111 = real(V_NNN%V(P3B,l_c,J,1,j_index(c))%v(v3b_idx(abc_1,def_1)),dp)
                            ME113 = real(V_NNN%V(P3B,l_c,J,2,j_index(c))%v(v3b_idx(abc_3,def_3)),dp)

                            pRho_fc = Rho%p(f,c)
                            nRho_fc = Rho%n(f,c)

                            ppSum = ppSum + ME113 * pRho_fc + (2.0_dp * ME111 + ME113) / 3.0_dp * nRho_fc
                            nnSum = nnSum + ME113 * nRho_fc + (2.0_dp * ME111 + ME113) / 3.0_dp * pRho_fc

                            pnSum_abde = pnSum_abde + (0.5_dp * ME001 + ME111 / 6.0_dp + ME113 / 3.0_dp) * (pRho_fc + nRho_fc) + &
                                         ( ME101 + ME011) / sqrt(12.0_dp) * (nRho_fc - pRho_fc)
                            pnSum_bade = pnSum_bade + (0.5_dp * ME001 - ME111 / 6.0_dp - ME113 / 3.0_dp) * (pRho_fc + nRho_fc) + &
                                         (-ME101 + ME011) / sqrt(12.0_dp) * (nRho_fc - pRho_fc)
                            pnSum_abed = pnSum_abed + (0.5_dp * ME001 - ME111 / 6.0_dp - ME113 / 3.0_dp) * (pRho_fc + nRho_fc) + &
                                         ( ME101 - ME011) / sqrt(12.0_dp) * (nRho_fc - pRho_fc)
                            pnSum_baed = pnSum_baed + (0.5_dp * ME001 + ME111 / 6.0_dp + ME113 / 3.0_dp) * (pRho_fc + nRho_fc) - &
                                         ( ME101 + ME011) / sqrt(12.0_dp) * (nRho_fc - pRho_fc)
                        end do

                        ! pp & nn elements ...
                        Ind = O2b_index(a,b,d,e,J,P,1,Orb_NN)
                        H_NN%pp(P,J)%v(Ind) = cRes * (H_NN%pp(P,J)%v(Ind) + is_J_hat * ppSum)
                        H_NN%nn(P,J)%v(Ind) = cRes * (H_NN%nn(P,J)%v(Ind) + is_J_hat * nnSum)

                        ! pn elements of all orderings ... (ba|V|de) & (ab|V|ed) are the same element for ab = de ...
                        Ind = O2b_index(a,b,d,e,J,P,0,Orb_NN)
                        H_NN%pn(P,J)%v(Ind) = cRes * (H_NN%pn(P,J)%v(Ind) + is_J_hat * pnSum_abde)
                        if (a /= b) then
                            Ind = O2b_index(b,a,d,e,J,P,0,Orb_NN)
                            H_NN%pn(P,J)%v(Ind) = cRes * (H_NN%pn(P,J)%v(Ind) + is_J_hat * Phase_ab * pnSum_bade)
                        end if
                        if (d /= e .and. Ket /= Bra) then
                            Ind = O2b_index(a,b,e,d,J,P,0,Orb_NN)
                            H_NN%pn(P,J)%v(Ind) = cRes * (H_NN%pn(P,J)%v(Ind) + is_J_hat * Phase_de * pnSum_abed)
                        end if
                        if (a /= b .and. d /= e) then
                            Ind = O2b_index(b,a,e,d,J,P,0,Orb_NN)
                            H_NN%pn(P,J)%v(Ind) = cRes * (H_NN%pn(P,J)%v(Ind) + is_J_hat * Phase_ab * Phase_de * pnSum_baed)
                        end if
                    end do
                end do
                !$omp end parallel do
            end do
        end do

        if (Rep) write(*,'(a)') TAB // 'Completed the calculation of density dependent residual 2-body interaction Hamiltonian H_NN ...'
    end subroutine H2b

    ! H_NN[Rho] & the Sigma-dressed mean field h_N[Rho,Sigma] ... counterpart of H1b_H2b ...
        ! h_cf += 1/4 1/(j_c+1) sum_J sum_ab sum_de (ab c|V_NNN|de f) Sigma_ab,de over all ordered pairs ab & de, written per V_NNN
        ! block over the pairs of triples (a b c) >= (d e f), the orderings ba & ed through their phases ... Sigma is read from its
        ! packed lower triangles (dev unpacks each (J,P) block into a dense matrix, same values) ...
    subroutine H1b_H2b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho,Sigma,h_N,H_NN,Report)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: T, Rho
        type(o2b_t), intent(in) :: V_NN, Sigma
        type(v3b_t), intent(in) :: V_NNN
        type(o1b_t), intent(out) :: h_N
        type(o2b_t), intent(out) :: H_NN
        logical, intent(in), optional :: Report
        type(imat_t), allocatable :: Triples(:,:,:,:)
        integer, allocatable :: N_trip(:,:,:,:), cf(:,:), ad(:,:)
        real(dp), allocatable :: ph(:,:), nh(:,:)
        integer :: N_max, a_max, l_max, J_max, cf_count, ad_count, Pass, a, b, c, d, e, f, ab, de, J, J_i, J_min, P, P2B, P3B
        integer :: l_c, j_ind, N_ab, N_ab_pn, abc, def, abc_0, abc_1, abc_3, def_0, def_1, def_3, ab_pn, ba_pn, de_pn, ed_pn, cf_i
        real(dp) :: Phase_ab, Phase_de, ME001, ME101, ME011, ME111, ME113, ppS, nnS, s_abde, s_bade, s_abed, s_baed
        real(dp) :: S00, S10, S01, S11, pSum, nSum, sr3
        logical :: Rep

        Rep = .true.
        if (present(Report)) Rep = Report
        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        l_max = N_max
        J_max = Params%Calc%N2max + 1
        sr3 = sqrt(3.0_dp)

        if (Rep) write(*,'(/,a)') 'Evaluating 2-body Hamiltonian H_NN ...'
        call H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho,H_NN,Rep)

        if (Rep) write(*,'(/,a)') 'Evaluating 1-body mean-field Hamiltonian h_N ...'
        call h1b_rho(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho,h_N)
        call lj_pairs(a_max,Orb,cf,cf_count,ad,ad_count)

        ! 3-body triples (a b c), a >= b, of every V_NNN block (P,l_c,J,j_index) ... columns (a, b, c, ab, i_0, i_1, i_3) with
            ! ab the pp index of (a,b) in Orb_NN & i_0, i_1, i_3 the indices of (ab T_ab=0 c), (ab T_ab=1 c) with T = 1/2 &
            ! (ab T_ab=1 c) with T = 3/2 ... H1b_NNN_triples ...
        allocate(N_trip(2,0:l_max,0:J_max,2), Triples(2,0:l_max,0:J_max,2))
        do Pass = 1, 2
            N_trip = 0
            do a = 1, a_max
                do b = 1, a
                    ab = a*(a - 1)/2 + b
                    J_min = abs(Orb(a)%j - Orb(b)%j)/2
                    do J_i = 1, size(Orb_NNN%Map(ab)%m,2)
                        J = J_min + J_i - 1
                        do c = 1, a_max
                            if (Orb_NNN%Map(ab)%m(1,J_i,c) == 0) cycle
                            l_c = Orb(c)%l
                            P = parity_of(Orb(a)%l + Orb(b)%l + l_c)
                            j_ind = V3b_jindex(l_c,Orb(c)%j)
                            N_trip(P,l_c,J,j_ind) = N_trip(P,l_c,J,j_ind) + 1
                            if (Pass == 2) Triples(P,l_c,J,j_ind)%v(:,N_trip(P,l_c,J,j_ind)) = &
                                [a, b, c, Orb_NN%Dic(a,b,J,1), int(Orb_NNN%Map(ab)%m(1,J_i,c)), &
                                 int(Orb_NNN%Map(ab)%m(2,J_i,c)), int(Orb_NNN%Map(ab)%m(3,J_i,c))]
                        end do
                    end do
                end do
            end do
            if (Pass == 1) then
                do j_ind = 1, 2
                    do J = 0, J_max
                        do l_c = 0, l_max
                            do P = 1, 2
                                allocate(Triples(P,l_c,J,j_ind)%v(7,N_trip(P,l_c,J,j_ind)))
                            end do
                        end do
                    end do
                end do
            end if
        end do

        allocate(ph(a_max,a_max), nh(a_max,a_max))
        ph = 0.0_dp
        nh = 0.0_dp
        do J = 0, J_max
            do P2B = 1, 2
                N_ab = Orb_NN%N(2,P2B,J)
                N_ab_pn = Orb_NN%N(1,P2B,J)
                if (N_ab == 0) cycle
                do l_c = 0, N_max
                    P3B = parity_of(P2B - 1 + l_c)
                    do j_ind = 1, 2
                        if (N_trip(P3B,l_c,J,j_ind) == 0) cycle
                        associate (Tr => Triples(P3B,l_c,J,j_ind)%v, V_T1 => V_NNN%V(P3B,l_c,J,1,j_ind)%v, &
                                   V_T3 => V_NNN%V(P3B,l_c,J,2,j_ind)%v)
                        !$omp parallel do schedule(dynamic,4) reduction(+:ph,nh) private(a,b,c,ab,abc_0,abc_1,abc_3,ab_pn,ba_pn, &
                        !$omp Phase_ab,def,d,e,f,de,def_0,def_1,def_3,de_pn,ed_pn,Phase_de,ME001,ME101,ME011,ME111,ME113, &
                        !$omp ppS,nnS,s_abde,s_bade,s_abed,s_baed,S00,S10,S01,S11,pSum,nSum)
                        do abc = 1, N_trip(P3B,l_c,J,j_ind)
                            a = Tr(1,abc)
                            b = Tr(2,abc)
                            c = Tr(3,abc)
                            ab = Tr(4,abc)
                            abc_0 = Tr(5,abc)
                            abc_1 = Tr(6,abc)
                            abc_3 = Tr(7,abc)
                            ab_pn = Orb_NN%Dic(a,b,J,0)
                            ba_pn = Orb_NN%Dic(b,a,J,0)
                            Phase_ab = phase((Orb(a)%j + Orb(b)%j)/2 - J)
                            do def = 1, abc
                                d = Tr(1,def)
                                e = Tr(2,def)
                                f = Tr(3,def)
                                de = Tr(4,def)
                                def_0 = Tr(5,def)
                                def_1 = Tr(6,def)
                                def_3 = Tr(7,def)
                                de_pn = Orb_NN%Dic(d,e,J,0)
                                ed_pn = Orb_NN%Dic(e,d,J,0)
                                Phase_de = phase((Orb(d)%j + Orb(e)%j)/2 - J)

                                ME001 = real(V_T1(v3b_idx(abc_0,def_0)),dp)
                                ME101 = real(V_T1(v3b_idx(abc_1,def_0)),dp)
                                ME011 = real(V_T1(v3b_idx(abc_0,def_1)),dp)
                                ME111 = real(V_T1(v3b_idx(abc_1,def_1)),dp)
                                ME113 = real(V_T3(v3b_idx(abc_3,def_3)),dp)

                                ! pp & nn Sigma ... the phases of ba (ed) in V_NNN & Sigma cancel, every ordering adds the same ...
                                ppS = Sigma%pp(P2B,J)%v(packed_index(max(de,ab),min(de,ab),N_ab))
                                nnS = Sigma%nn(P2B,J)%v(packed_index(max(de,ab),min(de,ab),N_ab))
                                if (a /= b) then
                                    ppS = 2.0_dp * ppS
                                    nnS = 2.0_dp * nnS
                                end if
                                if (d /= e) then
                                    ppS = 2.0_dp * ppS
                                    nnS = 2.0_dp * nnS
                                end if

                                ! pn Sigma of all orderings ab, ba & de, ed ...
                                s_abde = Sigma%pn(P2B,J)%v(packed_index(max(de_pn,ab_pn),min(de_pn,ab_pn),N_ab_pn))
                                s_bade = 0.0_dp
                                s_abed = 0.0_dp
                                s_baed = 0.0_dp
                                if (a /= b) s_bade = Phase_ab * Sigma%pn(P2B,J)%v(packed_index(max(de_pn,ba_pn),min(de_pn,ba_pn),N_ab_pn))
                                if (d /= e) s_abed = Phase_de * Sigma%pn(P2B,J)%v(packed_index(max(ed_pn,ab_pn),min(ed_pn,ab_pn),N_ab_pn))
                                if (a /= b .and. d /= e) s_baed = Phase_ab * Phase_de * &
                                    Sigma%pn(P2B,J)%v(packed_index(max(ed_pn,ba_pn),min(ed_pn,ba_pn),N_ab_pn))

                                ! ... projected on the isospin T_ab & T_de of the pairs ...
                                S00 = s_abde + s_bade + s_abed + s_baed
                                S10 = s_abde - s_bade + s_abed - s_baed
                                S01 = s_abde + s_bade - s_abed - s_baed
                                S11 = s_abde - s_bade - s_abed + s_baed

                                ! Proton field (c, f protons) ... <1 0 1/2 +1/2|1/2 +1/2> = -sqrt(1/3) ...
                                pSum = ME113 * ppS + (ME113 + 2.0_dp * ME111) / 3.0_dp * nnS + &
                                       2.0_dp * (ME001 * S00 - (ME101 * S10 + ME011 * S01) / sr3 + (ME111 + 2.0_dp * ME113) / 3.0_dp * S11)
                                ! Neutron field (c, f neutrons) ... <1 0 1/2 -1/2|1/2 -1/2> = +sqrt(1/3) ...
                                nSum = ME113 * nnS + (ME113 + 2.0_dp * ME111) / 3.0_dp * ppS + &
                                       2.0_dp * (ME001 * S00 + (ME101 * S10 + ME011 * S01) / sr3 + (ME111 + 2.0_dp * ME113) / 3.0_dp * S11)

                                ! (abc|V|def) adds to h_cf ... (def|V|abc), not stored, adds the same to h_fc ...
                                ph(c,f) = ph(c,f) + pSum
                                nh(c,f) = nh(c,f) + nSum
                                if (def /= abc) then
                                    ph(f,c) = ph(f,c) + pSum
                                    nh(f,c) = nh(f,c) + nSum
                                end if
                            end do
                        end do
                        !$omp end parallel do
                        end associate
                    end do
                end do
            end do
        end do

        do cf_i = 1, cf_count
            c = cf(1,cf_i)
            f = cf(2,cf_i)
            h_N%p(c,f) = h_N%p(c,f) + 0.25_dp / real(Orb(c)%j + 1,dp) * ph(c,f)
            h_N%n(c,f) = h_N%n(c,f) + 0.25_dp / real(Orb(c)%j + 1,dp) * nh(c,f)
        end do

        if (Rep) write(*,'(a)') TAB // 'Completed the calculation of 1-body mean-field Hamiltonian h_N ...'
    end subroutine H1b_H2b

end module nj_hamiltonian
