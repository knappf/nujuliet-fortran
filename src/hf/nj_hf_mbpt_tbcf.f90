! HF-MBPT 2-body correlation function (cumulant) Sigma ... counterpart of HF_MBPT_TBCF of src/HF/HF_MBPT.jl
    ! (NuJuliet dev, lines 1540-2181) ...
    ! Sigma_pphh = t_pphh (1st + 2nd order amplitudes: 2nd-order, 3rd-order pp, hh & ph diagrams) built in the HF basis in the
    ! O2B packed storage & transformed to the LHO basis with C^T ... normalisation as RPA_sigma & V2b_correlation_energy
    ! (1/4 V.Sigma = 2 (E2 + E3)) ...
    ! Fix F1 (documented deviation, mirrored in the patched reference): dev classifies the particles/holes of the nn channel
    ! (T = +1) of the 2nd-order, 3rd pp & 3rd hh passes with the PROTON occupations pO (HF_MBPT.jl:1905-1906, 1995-1996,
    ! 2085-2086); the port uses the neutron occupations nO, as dev's own ph pass does ...
    ! Port: dev runs four passes over all (J,P,Bra,Ket,T) adding each diagram to Sigma (initially 0); the port evaluates the
    ! four diagrams of one element in one sweep & sums them in the same order (s2 + s3pp) + s3hh) + s3ph -> identical values ...
module nj_hf_mbpt_tbcf
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    use nj_o1b, only : O1b_transpose
    use nj_o2b, only : O2b_initialize, O2b_transformation, packed_index
    use nj_hf_mbpt_common
    implicit none
    private

    public :: HF_MBPT_TBCF

contains

    ! dev HF_MBPT_TBCF(Params,C,h_N,V_NN,Orb,Orb_HF_MBPT,Orb_NN;Report=true) ... V_NN is the residual NN interaction in the
        ! HF basis, Sigma is returned in the LHO basis ...
    subroutine HF_MBPT_TBCF(Params,C,h_N,V_NN,Orb,Orb_HF_MBPT,Orb_NN,Sigma,Report)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: C, h_N
        type(o2b_t), intent(in) :: V_NN
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb1b_hf_mbpt_t), intent(in) :: Orb_HF_MBPT
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(out) :: Sigma
        logical, intent(in), optional :: Report
        type(pnsvector_t) :: Particle, Hole
        type(ispstate_t), allocatable :: Orb_particle(:), Orb_hole(:)
        type(orb2b_t) :: Orb_NN_loc
        logical :: Rep
        integer :: J_max, J, P, N_T0, N_T1, Bra, Ket, i_a, i_b, i_c, i_d, Ind
        integer(i8) :: t0

        Rep = .true.
        if (present(Report)) Rep = Report
        J_max = Params%Calc%N2max + 1

        call HF_MBPT_states(Orb_HF_MBPT,h_N,Particle,Hole)
        call HF_MBPT_isospin_states(Orb_HF_MBPT,h_N,Orb_particle,Orb_hole)

        ! Zero Sigma with the calculation-space channels (dev O2b_initialize(Params,Orb;Make_Orb_NN=false)) ...
            ! Port guard: its channels must coincide with those of Orb_NN, which index Sigma below ...
        call O2b_initialize(Params,Orb,Sigma,Orb_NN_loc)
        if (Orb_NN%J_max < J_max) call fatal('HF_MBPT_TBCF: Orb_NN has fewer J channels than N2max + 1')
        do J = 0, J_max
            do P = 1, 2
                if (Orb_NN%N(1,P,J) /= Orb_NN_loc%N(1,P,J) .or. Orb_NN%N(2,P,J) /= Orb_NN_loc%N(2,P,J)) &
                    call fatal('HF_MBPT_TBCF: 2-body orbitals of V_NN & of the calculation space differ')
            end do
        end do

        if (Rep) then
            write(*,'(/,a)') 'Calculating 2nd & 3rd order HF-MBPT 2-body correlation function matrix Sigma ...'
            write(*,'(a)') TAB // TAB // 'Calculating 2nd order particle-hole diagram ...'
            write(*,'(a)') TAB // TAB // 'Calculating 3rd order particle-particle diagram ...'
            write(*,'(a)') TAB // TAB // 'Calculating 3rd order hole-hole diagram ...'
            write(*,'(a)') TAB // TAB // 'Calculating 3rd order particle-hole diagram ...'
        end if

        call timer_start(t0)
        do J = 0, J_max
            do P = 1, 2
                N_T0 = Orb_NN%N(1,P,J)
                N_T1 = Orb_NN%N(2,P,J)

                ! proton-proton (T = -1) ... every (Bra,Ket) element is written by its Bra task only ...
                !$omp parallel do schedule(dynamic) private(Ket,i_a,i_b,i_c,i_d,Ind)
                do Bra = 1, N_T1
                    do Ket = 1, Bra
                        Ind = packed_index(Bra,Ket,N_T1)
                        i_a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                        i_b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                        i_c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                        i_d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                        Sigma%pp(P,J)%v(Ind) = sigma_element(i_a,-1,i_b,-1,i_c,-1,i_d,-1,J,P)
                    end do
                end do
                !$omp end parallel do

                ! proton-neutron (T = 0) ...
                !$omp parallel do schedule(dynamic) private(Ket,i_a,i_b,i_c,i_d,Ind)
                do Bra = 1, N_T0
                    do Ket = 1, Bra
                        Ind = packed_index(Bra,Ket,N_T0)
                        i_a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                        i_b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                        i_c = Orb_NN%Ind(1,P,J)%v(1,Ket)
                        i_d = Orb_NN%Ind(1,P,J)%v(2,Ket)
                        Sigma%pn(P,J)%v(Ind) = sigma_element(i_a,-1,i_b,1,i_c,-1,i_d,1,J,P)
                    end do
                end do
                !$omp end parallel do

                ! neutron-neutron (T = +1) ... occupations nO in all four diagrams (fix F1) ...
                !$omp parallel do schedule(dynamic) private(Ket,i_a,i_b,i_c,i_d,Ind)
                do Bra = 1, N_T1
                    do Ket = 1, Bra
                        Ind = packed_index(Bra,Ket,N_T1)
                        i_a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                        i_b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                        i_c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                        i_d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                        Sigma%nn(P,J)%v(Ind) = sigma_element(i_a,1,i_b,1,i_c,1,i_d,1,J,P)
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        if (Rep) call timer_report(t0)

        ! HF -> LHO basis ...
        call timer_start(t0)
        call O2b_transformation(Params,Orb,Orb_NN,Sigma,O1b_transpose(C),Report=Rep)
        if (Rep) call timer_report(t0)

        if (Rep) write(*,'(/,a)') TAB // 'Finished allocation & basis transformation of the HF-MBPT 2-body correlation function ...'

    contains

        ! dev particle(O) & hole(O) (local functions of HF_MBPT_TBCF) ...
        pure logical function is_particle(O)
            real(dp), intent(in) :: O
            is_particle = (O <= 1.0e-3_dp)
        end function is_particle

        pure logical function is_hole(O)
            real(dp), intent(in) :: O
            is_hole = (O >= 0.999_dp)
        end function is_hole

        pure real(dp) function occ(a,t)
            integer, intent(in) :: a, t
            if (t == -1) then
                occ = Orb(a)%pO
            else
                occ = Orb(a)%nO
            end if
        end function occ

        ! Sum of the four diagrams of one stored element <ab|Sigma|cd>_J (pphh, or hhpp evaluated as the swapped pphh) ...
        real(dp) function sigma_element(a,t_a,b,t_b,c,t_c,d,t_d,J,P) result(S)
            integer, intent(in) :: a, t_a, b, t_b, c, t_c, d, t_d, J, P
            real(dp) :: O_a, O_b, O_c, O_d
            O_a = occ(a,t_a)
            O_b = occ(b,t_b)
            O_c = occ(c,t_c)
            O_d = occ(d,t_d)
            S = 0.0_dp
            if (is_particle(O_a) .and. is_particle(O_b) .and. is_hole(O_c) .and. is_hole(O_d)) then
                S = sigma2(a,t_a,b,t_b,c,t_c,d,t_d,J,P)
                S = S + sigma3_pp(a,t_a,b,t_b,c,t_c,d,t_d,J,P)
                S = S + sigma3_hh(a,t_a,b,t_b,c,t_c,d,t_d,J,P)
                S = S + sigma3_ph(a,t_a,b,t_b,c,t_c,d,t_d,J)
            else if (is_hole(O_a) .and. is_hole(O_b) .and. is_particle(O_c) .and. is_particle(O_d)) then
                S = sigma2(c,t_c,d,t_d,a,t_a,b,t_b,J,P)
                S = S + sigma3_pp(c,t_c,d,t_d,a,t_a,b,t_b,J,P)
                S = S + sigma3_hh(c,t_c,d,t_d,a,t_a,b,t_b,J,P)
                S = S + sigma3_ph(c,t_c,d,t_d,a,t_a,b,t_b,J)
            end if
        end function sigma_element

        ! dev HF_MBPT2_sigma ... - V_abcd / (E_a + E_b - E_c - E_d) ...
        real(dp) function sigma2(a,t_a,b,t_b,c,t_c,d,t_d,J,P) result(sigma)
            integer, intent(in) :: a, t_a, b, t_b, c, t_c, d, t_d, J, P
            real(dp) :: E_a, E_b, E_c, E_d, ME
            E_a = HF_MBPT_h1b(a,t_a,h_N)
            E_b = HF_MBPT_h1b(b,t_b,h_N)
            E_c = HF_MBPT_h1b(c,t_c,h_N)
            E_d = HF_MBPT_h1b(d,t_d,h_N)
            sigma = 0.0_dp
            if (abs(t_a + t_b) == 2 .or. (t_a + t_b) == 0) then
                ME = HF_MBPT_V2b(a,b,c,d,J,P,t_a + t_b,V_NN,Orb,Orb_NN) / (E_a + E_b - E_c - E_d)
                sigma = sigma - ME
            end if
        end function sigma2

        ! dev HF_MBPT3_sigma_pp ... f V_ab,rs V_rs,cd / (E_a+E_b-E_c-E_d) / (E_r+E_s-E_c-E_d), r,s particles, f = 0.5 (pp, nn) / 1 (pn) ...
        real(dp) function sigma3_pp(a,t_a,b,t_b,c,t_c,d,t_d,J,P) result(sigma)
            integer, intent(in) :: a, t_a, b, t_b, c, t_c, d, t_d, J, P
            sigma = sigma3_ladder(a,t_a,b,t_b,c,t_c,d,t_d,J,P,.true.)
        end function sigma3_pp

        ! dev HF_MBPT3_sigma_hh ... f V_ab,rs V_rs,cd / (E_a+E_b-E_c-E_d) / (E_a+E_b-E_r-E_s), r,s holes ...
        real(dp) function sigma3_hh(a,t_a,b,t_b,c,t_c,d,t_d,J,P) result(sigma)
            integer, intent(in) :: a, t_a, b, t_b, c, t_c, d, t_d, J, P
            sigma = sigma3_ladder(a,t_a,b,t_b,c,t_c,d,t_d,J,P,.false.)
        end function sigma3_hh

        real(dp) function sigma3_ladder(a,t_a,b,t_b,c,t_c,d,t_d,J,P,pp) result(sigma)
            integer, intent(in) :: a, t_a, b, t_b, c, t_c, d, t_d, J, P
            logical, intent(in) :: pp
            real(dp) :: E_a, E_b, E_c, E_d, E_r, E_s, ME, f
            integer :: T_ab, ir, is
            type(spstate_t) :: sr, ss

            E_a = HF_MBPT_h1b(a,t_a,h_N)
            E_b = HF_MBPT_h1b(b,t_b,h_N)
            E_c = HF_MBPT_h1b(c,t_c,h_N)
            E_d = HF_MBPT_h1b(d,t_d,h_N)
            sigma = 0.0_dp
            T_ab = t_a + t_b
            if (T_ab == 0) then
                f = 1.0_dp
            else
                f = 0.5_dp
            end if
            if (T_ab /= -2 .and. T_ab /= 0 .and. T_ab /= 2) return

            do ir = 1, n_list(pp,T_ab,1)
                sr = state(pp,T_ab,1,ir)
                E_r = sr%E
                do is = 1, n_list(pp,T_ab,2)
                    ss = state(pp,T_ab,2,is)
                    E_s = ss%E
                    if ((.not. triangle(sr%j,ss%j,2*J)) .or. (mod(sr%l + ss%l,2) + 1 /= P)) cycle
                    if (pp) then
                        ME = f * HF_MBPT_V2b(a,b,sr%a,ss%a,J,P,T_ab,V_NN,Orb,Orb_NN) * &
                             HF_MBPT_V2b(sr%a,ss%a,c,d,J,P,T_ab,V_NN,Orb,Orb_NN) / (E_a + E_b - E_c - E_d) / (E_r + E_s - E_c - E_d)
                    else
                        ME = f * HF_MBPT_V2b(a,b,sr%a,ss%a,J,P,T_ab,V_NN,Orb,Orb_NN) * &
                             HF_MBPT_V2b(sr%a,ss%a,c,d,J,P,T_ab,V_NN,Orb,Orb_NN) / (E_a + E_b - E_c - E_d) / (E_a + E_b - E_r - E_s)
                    end if
                    sigma = sigma + ME
                end do
            end do
        end function sigma3_ladder

        ! Length of the summation list of orbital k (1 = r, 2 = s) for the channel T_ab ... particles (pp) or holes ...
            ! pn: r proton, s neutron ...
        pure integer function n_list(pp,T_ab,k)
            logical, intent(in) :: pp
            integer, intent(in) :: T_ab, k
            logical :: proton
            proton = (T_ab == -2) .or. (T_ab == 0 .and. k == 1)
            if (pp) then
                if (proton) then
                    n_list = size(Particle%p)
                else
                    n_list = size(Particle%n)
                end if
            else
                if (proton) then
                    n_list = size(Hole%p)
                else
                    n_list = size(Hole%n)
                end if
            end if
        end function n_list

        pure type(spstate_t) function state(pp,T_ab,k,i)
            logical, intent(in) :: pp
            integer, intent(in) :: T_ab, k, i
            logical :: proton
            proton = (T_ab == -2) .or. (T_ab == 0 .and. k == 1)
            if (pp) then
                if (proton) then
                    state = Particle%p(i)
                else
                    state = Particle%n(i)
                end if
            else
                if (proton) then
                    state = Hole%p(i)
                else
                    state = Hole%n(i)
                end if
            end if
        end function state

        ! dev HF_MBPT3_sigma_ph_direct ...
        real(dp) function sigma3_ph_direct(sa,sb,sc,sd,J) result(sigma)
            type(ispstate_t), intent(in) :: sa, sb, sc, sd
            integer, intent(in) :: J
            integer :: I, I_min, I_max, ir, jf
            type(ispstate_t) :: sr, sf
            real(dp) :: Amp_abcd_I, V_fbcr, V_arfd, ME

            sigma = 0.0_dp
            I_min = max(abs(sa%j - sd%j), abs(sb%j - sc%j))/2
            I_max = min(sa%j + sd%j, sb%j + sc%j)/2
            do I = I_min, I_max
                Amp_abcd_I = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*J,sc%j,sd%j,2*I) / (sa%E + sb%E - sc%E - sd%E)
                do ir = 1, size(Orb_particle)
                    sr = Orb_particle(ir)
                    do jf = 1, size(Orb_hole)
                        sf = Orb_hole(jf)
                        if ((mod(sf%l + sb%l,2) /= mod(sc%l + sr%l,2)) .or. (mod(sa%l + sr%l,2) /= mod(sf%l + sd%l,2)) .or. &
                            (sb%t + sf%t /= sc%t + sr%t) .or. (sa%t + sr%t /= sd%t + sf%t) .or. (.not. triangle(sr%j,sf%j,2*I))) cycle
                        V_fbcr = HF_MBPT3_V2b_Pandya(sf,sb,sc,sr,I,Orb,Orb_NN,V_NN)
                        V_arfd = HF_MBPT3_V2b_Pandya(sa,sr,sf,sd,I,Orb,Orb_NN,V_NN)
                        ME = Amp_abcd_I * V_fbcr * V_arfd / (sa%E + sr%E - sd%E - sf%E)
                        sigma = sigma - ME
                    end do
                end do
            end do
        end function sigma3_ph_direct

        ! dev HF_MBPT3_sigma_ph ... direct term & its three exchanges with Phase_ab, Phase_cd ...
        real(dp) function sigma3_ph(a,t_a,b,t_b,c,t_c,d,t_d,J) result(sigma)
            integer, intent(in) :: a, t_a, b, t_b, c, t_c, d, t_d, J
            type(ispstate_t) :: sa, sb, sc, sd
            real(dp) :: Phase_ab, Phase_cd

            sa = ispstate_t(a,Orb(a)%l,Orb(a)%j,t_a,HF_MBPT_h1b(a,t_a,h_N))
            sb = ispstate_t(b,Orb(b)%l,Orb(b)%j,t_b,HF_MBPT_h1b(b,t_b,h_N))
            sc = ispstate_t(c,Orb(c)%l,Orb(c)%j,t_c,HF_MBPT_h1b(c,t_c,h_N))
            sd = ispstate_t(d,Orb(d)%l,Orb(d)%j,t_d,HF_MBPT_h1b(d,t_d,h_N))

            Phase_ab = phase((sa%j + sb%j)/2 - J)
            Phase_cd = phase((sc%j + sd%j)/2 - J)
            if (t_a == t_b) Phase_ab = -Phase_ab
            if (t_c == t_d) Phase_cd = -Phase_cd

            sigma = sigma3_ph_direct(sa,sb,sc,sd,J) + &
                    sigma3_ph_direct(sb,sa,sc,sd,J) * Phase_ab + &
                    sigma3_ph_direct(sa,sb,sd,sc,J) * Phase_cd + &
                    sigma3_ph_direct(sb,sa,sd,sc,J) * Phase_ab * Phase_cd
        end function sigma3_ph

    end subroutine HF_MBPT_TBCF

end module nj_hf_mbpt_tbcf
