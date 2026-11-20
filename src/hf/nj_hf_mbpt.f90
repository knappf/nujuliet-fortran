! HF many-body perturbation theory (NuJuliet dev) ... counterpart of src/HF/HF_MBPT.jl ...
    !   HF_MBPT                   ... particle/hole orbitals -> MBPT(2)+(3) ground-state energy -> MBPT OBDM & natural orbitals,
    !   HF_MBPT_energy            ... E^(2), E^(3)_pp, E^(3)_ph, E^(3)_hh (summary block in io/<Path>/HF/HF_summary.dat),
    !   HF_MBPT_OBDM              ... 1p1h & 2p2h corrections of the HF-basis OBDM, natural orbitals (OBDM_natural_orbitals),
    !                                 HF/densities/HF_MBPT_{radial_densities,occupation}.dat & (Export) HF/bin/{orb,C,rho}_temp.bin,
    !   HF_MBPT_occupation_check  ... MBPT needs an idempotent HF density,
    !   HF_MBPT_export            ... correlated reference io/<Path>/HF-MBPT/bin/{orb,C,H1b,rho,H2b[,sigma]}.bin ...
    ! Helpers (orbital lists, Pandya elements) are in nj_hf_mbpt_common, the 2-body correlation function in nj_hf_mbpt_tbcf ...
    ! Isospin convention of HF_MBPT.jl: t = -1 proton, t = +1 neutron; single-particle energies E_a = h_N(a,a) ...
    ! Threads: every MBPT energy sector is one parallel loop over the (flattened) outer pair list with reduction(+); the OBDM
    ! loops reduce into private copies of dpRho/dnRho (dev writes elements owned by its loop index, same sums) ...
module nj_hf_mbpt
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_orbitals, only : orbitals_export, orbitals_import, OBDM_natural_orbitals
    use nj_o1b, only : O1b_import, O1b_export, O1b_transformation, O1b_transpose
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn, O2b_transformation, O2b_export
    use nj_obdm, only : rho_export
    use nj_hamiltonian, only : H1b, H2b, H1b_H2b
    use nj_hf_mbpt_common
    use nj_hf_mbpt_tbcf, only : HF_MBPT_TBCF
    implicit none
    private

    public :: HF_MBPT, HF_MBPT_orbitals_make, HF_MBPT_energy, HF_MBPT_OBDM, HF_MBPT_occupation, HF_MBPT_occupation_check
    public :: HF_MBPT_TBCF, HF_MBPT_export

contains

    ! dev HF_MBPT(Params,Orb,Orb_NN,h_N,V_NN) ... V_NN is the residual NN interaction in the HF basis ...
    subroutine HF_MBPT(Params,Orb,Orb_NN,h_N,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: h_N
        type(o2b_t), intent(in) :: V_NN
        type(orb1b_hf_mbpt_t) :: Orb_HF_MBPT
        integer(i8) :: t0

        Orb_HF_MBPT = HF_MBPT_orbitals_make(Params,Orb)

        call timer_start(t0)
        call HF_MBPT_energy(Params,Orb,Orb_HF_MBPT,Orb_NN,h_N,V_NN)
        call timer_report(t0)

        call timer_start(t0)
        call HF_MBPT_OBDM(Params,Orb,Orb_HF_MBPT,Orb_NN,h_N,V_NN)
        call timer_report(t0)
    end subroutine HF_MBPT

    ! MBPT(2) & MBPT(3) corrections to the ground-state energy ... dev HF_MBPT_energy (154-753) ...
    subroutine HF_MBPT_energy(Params,Orb,Orb_HF_MBPT,Orb_NN,h_N,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb1b_hf_mbpt_t), intent(in) :: Orb_HF_MBPT
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: h_N
        type(o2b_t), intent(in) :: V_NN
        type(pnsvector_t) :: Particle, Hole
        type(ispstate_t), allocatable :: Orb_particle(:), Orb_hole(:)
        real(dp) :: dE, dE2, dE3_pp, dE3_ph, dE3_hh, A
        integer :: u
        integer(i8) :: t0

        call HF_MBPT_states(Orb_HF_MBPT,h_N,Particle,Hole)

        ! HF-MBPT(2) ... pp, nn & pn ...
        write(*,'(/,a)') 'Calculating HF-MBPT(2) ground-state energy correction ... '
        dE2 = 0.0_dp
        dE2 = dE2 + mbpt2(Particle%p,Particle%p,Hole%p,Hole%p,-2,0.25_dp)
        dE2 = dE2 + mbpt2(Particle%n,Particle%n,Hole%n,Hole%n,2,0.25_dp)
        dE2 = dE2 + mbpt2(Particle%p,Particle%n,Hole%p,Hole%n,0,1.0_dp)

        ! HF-MBPT(3) ...
        write(*,'(/,a)') 'Calculating HF-MBPT(3) ground-state energy correction ... '

        ! particle-particle sector ... pp, pn, nn ...
        dE3_pp = 0.0_dp
        call timer_start(t0)
        dE3_pp = dE3_pp + mbpt3_ladder(Particle%p,Particle%p,Hole%p,Hole%p,-2,0.125_dp,.true.)
        call timer_report(t0)
        call timer_start(t0)
        dE3_pp = dE3_pp + mbpt3_ladder(Particle%p,Particle%n,Hole%p,Hole%n,0,1.0_dp,.true.)
        call timer_report(t0)
        call timer_start(t0)
        dE3_pp = dE3_pp + mbpt3_ladder(Particle%n,Particle%n,Hole%n,Hole%n,2,0.125_dp,.true.)
        call timer_report(t0)

        ! particle-hole sector ...
        call HF_MBPT_isospin_states(Orb_HF_MBPT,h_N,Orb_particle,Orb_hole)
        call timer_start(t0)
        dE3_ph = mbpt3_ph(Orb_particle,Orb_hole)
        call timer_report(t0)

        ! hole-hole sector ... pp, pn, nn ...
        dE3_hh = 0.0_dp
        call timer_start(t0)
        dE3_hh = dE3_hh + mbpt3_ladder(Particle%p,Particle%p,Hole%p,Hole%p,-2,0.125_dp,.false.)
        call timer_report(t0)
        call timer_start(t0)
        dE3_hh = dE3_hh + mbpt3_ladder(Particle%p,Particle%n,Hole%p,Hole%n,0,1.0_dp,.false.)
        call timer_report(t0)
        call timer_start(t0)
        dE3_hh = dE3_hh + mbpt3_ladder(Particle%n,Particle%n,Hole%n,Hole%n,2,0.125_dp,.false.)
        call timer_report(t0)

        dE = dE2 + dE3_pp + dE3_ph + dE3_hh

        write(*,'(a)') TAB // 'HF-MBPT(2) energy    ...   E^(2)    = ' // jstr(round_sig(dE2,9)) // ' MeV'
        write(*,'(a)') TAB // 'HF-MBPT(3) energy    ...   E^(3)    = ' // jstr(round_sig(dE3_pp + dE3_ph + dE3_hh,9)) // ' MeV'
        write(*,'(a)') TAB // TAB // '                           E^(3)_pp = ' // jstr(round_sig(dE3_pp,9)) // ' MeV'
        write(*,'(a)') TAB // TAB // '                           E^(3)_ph = ' // jstr(round_sig(dE3_ph,9)) // ' MeV'
        write(*,'(a)') TAB // TAB // '                           E^(3)_hh = ' // jstr(round_sig(dE3_hh,9)) // ' MeV'
        write(*,'(a,/)') TAB // 'HF-MBPT(3) energy    ...   E^(tot)  = ' // jstr(round_sig(dE,9)) // ' MeV'

        ! Summary block (dev 734-750; the Julia text ends without a newline, rho_export continues with an empty line) ...
        A = real(Params%Calc%A,dp)
        open(newunit=u, file=io_path(Params,'/HF/HF_summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Spherical Hartree-Fock Many-Body Perturbation Theory solution review:'
        write(u,'(/,a)') line('E_0^(2)        = ',dE2,'LO HF-MBPT(2) correction to ground-state energy')
        write(u,'(a)') line('E_0^(3)        = ',dE3_pp + dE3_ph + dE3_hh,'NLO HF-MBPT(3) correction to ground-state energy')
        write(u,'(a)') line('E_0^(3)_pp     = ',dE3_pp,'NLO HF-MBPT(3) pp correction to ground-state energy')
        write(u,'(a)') line('E_0^(3)_ph     = ',dE3_ph,'NLO HF-MBPT(3) ph correction to ground-state energy')
        write(u,'(a)') line('E_0^(3)_hh     = ',dE3_hh,'NLO HF-MBPT(3) hh correction to ground-state energy')
        write(u,'(a)') repeat('_',124)
        write(u,'(/,a)') line('E_0^(tot)      = ',dE,'NLO HF-MBPT(2+3) correction to ground-state energy')
        write(u,'(//,a)') line('E_0^(2) / A    = ',dE2 / A,'NLO HF-MBPT(2) correction to ground-state energy per nucleon')
        ! (sic) no "per nucleon" in dev ...
        write(u,'(a)') line('E_0^(3) / A    = ',(dE3_hh + dE3_ph + dE3_pp) / A,'NLO HF-MBPT(3) correction to ground-state energy')
        write(u,'(a)') line('E_0^(3)_pp / A = ',dE3_pp / A,'NLO HF-MBPT(3) pp correction to ground-state energy per nucleon')
        write(u,'(a)') line('E_0^(3)_ph / A = ',dE3_ph / A,'NLO HF-MBPT(3) ph correction to ground-state energy per nucleon')
        write(u,'(a)') line('E_0^(3)_hh / A = ',dE3_hh / A,'NLO HF-MBPT(3) hh correction to ground-state energy per nucleon')
        write(u,'(a)') repeat('_',124)
        write(u,'(/,a)') line('E_0^(tot) / A  = ',dE / A,'NLO HF-MBPT(2+3) correction to ground-state energy per nucleon')
        close(u)

    contains

        ! "%15.8f \t MeV \t\t ... \t <Text>" ...
        function line(Label,x,Text) result(s)
            character(len=*), intent(in) :: Label, Text
            real(dp), intent(in) :: x
            character(len=:), allocatable :: s
            s = Label // cfmt_f(x,15,8) // ' ' // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // ' ' // Text
        end function line

        ! MBPT(2) sector (dev 169-265): f (2J+1) |<pq|V|hg>_J|^2 / (E_h+E_g-E_p-E_q), p in P1, q in P2, h in H1, g in H2 ...
        real(dp) function mbpt2(P1,P2,H1,H2,T_ab,f) result(dE_s)
            type(spstate_t), intent(in) :: P1(:), P2(:), H1(:), H2(:)
            integer, intent(in) :: T_ab
            real(dp), intent(in) :: f
            type(spstate_t) :: sp, sq, sh, sg
            integer :: p, q, h, g, J, Par
            real(dp) :: Sum, E_pqhg

            dE_s = 0.0_dp
            !$omp parallel do schedule(dynamic) reduction(+:dE_s) private(sp,q,sq,Par,J,h,sh,g,sg,E_pqhg,Sum)
            do p = 1, size(P1)
                sp = P1(p)
                Sum = 0.0_dp
                do q = 1, size(P2)
                    sq = P2(q)
                    Par = mod(sp%l + sq%l,2) + 1
                    do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                        do h = 1, size(H1)
                            sh = H1(h)
                            do g = 1, size(H2)
                                sg = H2(g)
                                if ((Par /= (mod(sh%l + sg%l,2) + 1)) .or. (abs(sh%j - sg%j) > 2*J) .or. (2*J > (sh%j + sg%j))) cycle
                                E_pqhg = sh%E + sg%E - sp%E - sq%E
                                Sum = Sum + f * real(2*J + 1,dp) * &
                                      abs(HF_MBPT_V2b(sp%a,sq%a,sh%a,sg%a,J,Par,T_ab,V_NN,Orb,Orb_NN))**2 / E_pqhg
                            end do
                        end do
                    end do
                end do
                dE_s = dE_s + Sum
            end do
            !$omp end parallel do
        end function mbpt2

        ! MBPT(3) pp (Ladder = true, dev 270-452) & hh (Ladder = false, dev 540-716) sectors, parallel over the ordered
            ! particle pairs (p,q) in P1 x P2 ...
            ! pp: f (2J+1) V_hg,pq V_pq,rs V_rs,hg / E_pqhg / E_rshg, r,s particles ...
            ! hh: f (2J+1) V_pq,hg V_hg,fe V_fe,pq / E_pqhg / E_pqfe, f,e holes ...
        real(dp) function mbpt3_ladder(P1,P2,H1,H2,T_ab,f,Ladder) result(dE_s)
            type(spstate_t), intent(in) :: P1(:), P2(:), H1(:), H2(:)
            integer, intent(in) :: T_ab
            real(dp), intent(in) :: f
            logical, intent(in) :: Ladder
            type(spstate_t) :: sp, sq, sh, sg, sr, ss
            integer :: k, N1, ip, iq, h, g, r, s, J, P_pq, N_r, N_s
            real(dp) :: Sum, E_pqhg, E_2, J_hat, ME

            N1 = size(P1)
            if (Ladder) then
                N_r = size(P1)
                N_s = size(P2)
            else
                N_r = size(H1)
                N_s = size(H2)
            end if
            dE_s = 0.0_dp
            !$omp parallel do schedule(dynamic) reduction(+:dE_s) &
            !$omp private(ip,iq,sp,sq,P_pq,h,sh,g,sg,E_pqhg,r,sr,s,ss,E_2,J,J_hat,ME,Sum)
            do k = 1, size(P1)*size(P2)
                ip = mod(k - 1,N1) + 1
                iq = (k - 1)/N1 + 1
                sp = P1(ip)
                sq = P2(iq)
                P_pq = mod(sp%l + sq%l,2) + 1
                Sum = 0.0_dp
                do h = 1, size(H1)
                    sh = H1(h)
                    do g = 1, size(H2)
                        sg = H2(g)
                        if ((mod(sh%l + sg%l,2) + 1) /= P_pq) cycle
                        E_pqhg = sh%E + sg%E - sp%E - sq%E
                        do r = 1, N_r
                            if (Ladder) then
                                sr = P1(r)
                            else
                                sr = H1(r)
                            end if
                            do s = 1, N_s
                                if (Ladder) then
                                    ss = P2(s)
                                else
                                    ss = H2(s)
                                end if
                                if ((mod(sr%l + ss%l,2) + 1) /= P_pq) cycle
                                if (Ladder) then
                                    E_2 = sh%E + sg%E - sr%E - ss%E
                                else
                                    E_2 = sr%E + ss%E - sp%E - sq%E
                                end if
                                do J = max(abs(sh%j - sg%j),abs(sp%j - sq%j),abs(sr%j - ss%j))/2, &
                                       min(sh%j + sg%j,sp%j + sq%j,sr%j + ss%j)/2
                                    J_hat = real(2*J + 1,dp)
                                    if (Ladder) then
                                        ME = f * J_hat * HF_MBPT_V2b(sh%a,sg%a,sp%a,sq%a,J,P_pq,T_ab,V_NN,Orb,Orb_NN) * &
                                             HF_MBPT_V2b(sp%a,sq%a,sr%a,ss%a,J,P_pq,T_ab,V_NN,Orb,Orb_NN) * &
                                             HF_MBPT_V2b(sr%a,ss%a,sh%a,sg%a,J,P_pq,T_ab,V_NN,Orb,Orb_NN) / E_pqhg / E_2
                                    else
                                        ME = f * J_hat * HF_MBPT_V2b(sp%a,sq%a,sh%a,sg%a,J,P_pq,T_ab,V_NN,Orb,Orb_NN) * &
                                             HF_MBPT_V2b(sh%a,sg%a,sr%a,ss%a,J,P_pq,T_ab,V_NN,Orb,Orb_NN) * &
                                             HF_MBPT_V2b(sr%a,ss%a,sp%a,sq%a,J,P_pq,T_ab,V_NN,Orb,Orb_NN) / E_pqhg / E_2
                                    end if
                                    Sum = Sum + ME
                                end do
                            end do
                        end do
                    end do
                end do
                dE_s = dE_s + Sum
            end do
            !$omp end parallel do
        end function mbpt3_ladder

        ! MBPT(3) ph sector (dev 454-538), parallel over the ordered pairs (p,q) of the combined particle list ...
        real(dp) function mbpt3_ph(Orb_particle,Orb_hole) result(dE_s)
            type(ispstate_t), intent(in) :: Orb_particle(:), Orb_hole(:)
            type(ispstate_t) :: sp, sq, sh, sg, sr, sf
            integer :: k, N_p, ip, iq, h, g, r, f, J
            real(dp) :: Sum, J_hat, V_pqhg, V_hfrq, V_fgpr, E_pqhg, E_prgf, ME

            N_p = size(Orb_particle)
            dE_s = 0.0_dp
            !$omp parallel do schedule(dynamic) reduction(+:dE_s) &
            !$omp private(ip,iq,sp,sq,h,sh,g,sg,J,J_hat,V_pqhg,E_pqhg,r,sr,f,sf,V_hfrq,V_fgpr,E_prgf,ME,Sum)
            do k = 1, N_p*N_p
                ip = mod(k - 1,N_p) + 1
                iq = (k - 1)/N_p + 1
                sp = Orb_particle(ip)
                sq = Orb_particle(iq)
                Sum = 0.0_dp
                do h = 1, size(Orb_hole)
                    sh = Orb_hole(h)
                    do g = 1, size(Orb_hole)
                        sg = Orb_hole(g)
                        if ((sp%t + sq%t) /= (sh%t + sg%t)) cycle
                        do J = max(abs(sp%j - sg%j),abs(sh%j - sq%j))/2, min(sp%j + sg%j,sh%j + sq%j)/2
                            J_hat = real(2*J + 1,dp)
                            V_pqhg = HF_MBPT3_V2b_Pandya(sp,sq,sh,sg,J,Orb,Orb_NN,V_NN)
                            if (abs(V_pqhg) < 1.0e-12_dp) cycle
                            E_pqhg = sh%E + sg%E - sp%E - sq%E
                            do r = 1, N_p
                                sr = Orb_particle(r)
                                do f = 1, size(Orb_hole)
                                    sf = Orb_hole(f)
                                    if (((sh%t + sr%t) /= (sq%t + sf%t)) .or. ((sg%t + sf%t) /= (sp%t + sr%t)) .or. &
                                        (.not. triangle(sf%j,sr%j,2*J))) cycle
                                    V_hfrq = HF_MBPT3_V2b_Pandya(sh,sr,sf,sq,J,Orb,Orb_NN,V_NN)
                                    if (abs(V_hfrq) < 1.0e-12_dp) cycle
                                    V_fgpr = HF_MBPT3_V2b_Pandya(sf,sg,sp,sr,J,Orb,Orb_NN,V_NN)
                                    if (abs(V_fgpr) < 1.0e-12_dp) cycle
                                    E_prgf = sg%E + sf%E - sp%E - sr%E
                                    ME = J_hat * V_pqhg * V_hfrq * V_fgpr / E_pqhg / E_prgf
                                    Sum = Sum + ME
                                end do
                            end do
                        end do
                    end do
                end do
                dE_s = dE_s + Sum
            end do
            !$omp end parallel do
        end function mbpt3_ph

    end subroutine HF_MBPT_energy

    ! MBPT corrections of the HF-basis OBDM, natural orbitals & radial densities ... dev HF_MBPT_OBDM (755-1538) ...
        ! The dRho loops (1p1h pp/pn & 2p2h pp/pn for protons, same for neutrons) are those of the base version (formulas,
        ! filters & index orders identical in dev) ...
    subroutine HF_MBPT_OBDM(Params,Orb,Orb_HF_MBPT,Orb_NN,h_N,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb1b_hf_mbpt_t), intent(in) :: Orb_HF_MBPT
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: h_N
        type(o2b_t), intent(in) :: V_NN
        type(pnsvector_t) :: Particle, Hole
        type(pninteger_t) :: N_Particle, N_Hole
        real(dp), allocatable :: dpRho(:,:), dnRho(:,:)
        type(o1b_t) :: C_HF, Rho_HF, Rho, D, Rho_NAT, C_NAT, Rho_LHO
        type(orb1b_t), allocatable :: Orb_NAT(:)
        integer :: a_max, p, h, f, J, r, q, g, a, c, P_pf, P_hf, P_pr, P_hr, P_pq, P_hg
        type(spstate_t) :: sp, sh, sf, sr, sq, sg, sa, sc
        real(dp) :: ME

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(/,a)') 'Calculating HF-MBPT OBDM corrections ...'

        call HF_MBPT_states(Orb_HF_MBPT,h_N,Particle,Hole)
        N_Particle = Orb_HF_MBPT%N_p
        N_Hole = Orb_HF_MBPT%N_h

        allocate(dpRho(a_max,a_max), dnRho(a_max,a_max))
        dpRho = 0.0_dp
        dnRho = 0.0_dp

        ! pRho 1p-1h pp ...
        !$omp parallel do schedule(dynamic) reduction(+:dpRho) private(sp,h,sh,f,sf,P_pf,P_hf,J,r,sr,q,sq,ME,g,sg,P_pr,P_hr)
        do p = 1, N_Particle%p
            sp = Particle%p(p)
            do h = 1, N_Hole%p
                sh = Hole%p(h)
                if (sp%j == sh%j .and. sp%l == sh%l) then
                    do f = 1, N_Hole%p
                        sf = Hole%p(f)
                        P_pf = mod(sp%l + sf%l,2) + 1
                        P_hf = mod(sh%l + sf%l,2) + 1
                        do J = abs(sp%j - sf%j)/2, (sp%j + sf%j)/2
                            do r = 1, N_Particle%p
                                sr = Particle%p(r)
                                do q = 1, N_Particle%p
                                    sq = Particle%p(q)
                                    if ((abs(sq%j - sr%j) <= 2*J) .and. (2*J <= (sq%j + sr%j)) .and. (mod(sq%l + sr%l, 2) == mod(sh%l + sf%l, 2))) then
                                        ME = 0.5_dp * real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pp(sp%a,sf%a,sq%a,sr%a,J,P_pf,V_NN,Orb,Orb_NN) * &
                                             O2b_pp(sq%a,sr%a,sh%a,sf%a,J,P_hf,V_NN,Orb,Orb_NN) / (sh%E - sp%E) / (sh%E + sf%E - sq%E - sr%E)
                                        dpRho(sp%a,sh%a) = dpRho(sp%a,sh%a) + ME
                                        dpRho(sh%a,sp%a) = dpRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do

                    do r = 1, N_Particle%p
                        sr = Particle%p(r)
                        P_pr = mod(sp%l + sr%l,2) + 1
                        P_hr = mod(sp%l + sr%l,2) + 1
                        do J = abs(sh%j - sr%j)/2, (sh%j + sr%j)/2
                            do f = 1, N_Hole%p
                                sf = Hole%p(f)
                                do g = 1, N_Hole%p
                                    sg = Hole%p(g)
                                    if ((abs(sg%j - sf%j) <= 2*J) .and. (2*J <= (sg%j + sf%j)) .and. (mod(sp%l + sr%l, 2) == mod(sg%l + sf%l, 2))) then
                                        ME = -0.5_dp * real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pp(sp%a,sr%a,sg%a,sf%a,J,P_pr,V_NN,Orb,Orb_NN) * &
                                             O2b_pp(sg%a,sf%a,sh%a,sr%a,J,P_hr,V_NN,Orb,Orb_NN) / (sh%E - sp%E) / (sg%E + sf%E - sp%E - sr%E)
                                        dpRho(sp%a,sh%a) = dpRho(sp%a,sh%a) + ME
                                        dpRho(sh%a,sp%a) = dpRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do
                end if
            end do
        end do
        !$omp end parallel do

        ! pRho 1p-1h pn ...
        !$omp parallel do schedule(dynamic) reduction(+:dpRho) private(sp,h,sh,f,sf,P_pf,P_hf,J,r,sr,q,sq,ME,g,sg,P_pr,P_hr)
        do p = 1, N_Particle%p
            sp = Particle%p(p)
            do h = 1, N_Hole%p
                sh = Hole%p(h)
                if (sp%j == sh%j .and. sp%l == sh%l) then
                    do f = 1, N_Hole%n
                        sf = Hole%n(f)
                        P_pf = mod(sp%l + sf%l,2) + 1
                        P_hf = mod(sh%l + sf%l,2) + 1
                        do J = abs(sp%j - sf%j)/2, (sp%j + sf%j)/2
                            do r = 1, N_Particle%n
                                sr = Particle%n(r)
                                do q = 1, N_Particle%p
                                    sq = Particle%p(q)
                                    if ((abs(sq%j - sr%j) <= 2*J) .and. (2*J <= (sq%j + sr%j)) .and. (mod(sq%l + sr%l, 2) == mod(sh%l + sf%l, 2))) then
                                        ME = real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pn(sp%a,sf%a,sq%a,sr%a,J,P_pf,V_NN,Orb_NN) * &
                                             O2b_pn(sq%a,sr%a,sh%a,sf%a,J,P_hf,V_NN,Orb_NN) / (sh%E - sp%E) / (sh%E + sf%E - sq%E - sr%E)
                                        dpRho(sp%a,sh%a) = dpRho(sp%a,sh%a) + ME
                                        dpRho(sh%a,sp%a) = dpRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do

                    do r = 1, N_Particle%n
                        sr = Particle%n(r)
                        P_pr = mod(sp%l + sr%l,2) + 1
                        P_hr = mod(sh%l + sr%l,2) + 1
                        do J = abs(sh%j - sr%j)/2, (sh%j + sr%j)/2
                            do f = 1, N_Hole%n
                                sf = Hole%n(f)
                                do g = 1, N_Hole%p
                                    sg = Hole%p(g)
                                    if ((abs(sg%j - sf%j) <= 2*J) .and. (2*J <= (sg%j + sf%j)) .and. (mod(sp%l + sr%l, 2) == mod(sg%l + sf%l, 2))) then
                                        ME = -real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pn(sp%a,sr%a,sg%a,sf%a,J,P_pr,V_NN,Orb_NN) * &
                                             O2b_pn(sg%a,sf%a,sh%a,sr%a,J,P_hr,V_NN,Orb_NN) / (sh%E - sp%E) / (sg%E + sf%E - sp%E - sr%E)
                                        dpRho(sp%a,sh%a) = dpRho(sp%a,sh%a) + ME
                                        dpRho(sh%a,sp%a) = dpRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do
                end if
            end do
        end do
        !$omp end parallel do

        ! pRho 2p-2h pp ... particle block ...
        !$omp parallel do schedule(dynamic) reduction(+:dpRho) private(sp,q,sq,P_pq,J,h,sh,g,sg,P_hg,a,sa,ME)
        do p = 1, N_Particle%p
            sp = Particle%p(p)
            do q = 1, N_Particle%p
                sq = Particle%p(q)
                P_pq = mod(sp%l + sq%l,2) + 1
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%p
                        sh = Hole%p(h)
                        do g = 1, N_Hole%p
                            sg = Hole%p(g)
                            P_hg = mod(sh%l + sg%l,2) + 1
                            if ((abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do a = 1, N_Particle%p
                                    sa = Particle%p(a)
                                    if ((mod(sh%l + sg%l, 2) == mod(sq%l + sa%l, 2)) .and. (sa%j == sp%j) .and. (sa%l == sp%l)) then
                                        ME = 0.5_dp * real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pp(sp%a,sq%a,sh%a,sg%a,J,P_pq,V_NN,Orb,Orb_NN) * &
                                             O2b_pp(sh%a,sg%a,sa%a,sq%a,J,P_hg,V_NN,Orb,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sh%E + sg%E - sa%E - sq%E)
                                        dpRho(sp%a,sa%a) = dpRho(sp%a,sa%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! pRho 2p-2h pp ... hole block ...
        !$omp parallel do schedule(dynamic) reduction(+:dpRho) private(sh,g,sg,P_hg,J,p,sp,q,sq,P_pq,c,sc,ME)
        do h = 1, N_Hole%p
            sh = Hole%p(h)
            do g = 1, N_Hole%p
                sg = Hole%p(g)
                P_hg = mod(sh%l + sg%l,2) + 1
                do J = abs(sh%j - sg%j)/2, (sh%j + sg%j)/2
                    do p = 1, N_Particle%p
                        sp = Particle%p(p)
                        do q = 1, N_Particle%p
                            sq = Particle%p(q)
                            P_pq = mod(sp%l + sq%l,2) + 1
                            if ((abs(sp%j - sq%j) <= 2*J) .and. (2*J <= (sp%j + sq%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do c = 1, N_Hole%p
                                    sc = Hole%p(c)
                                    if ((mod(sg%l + sc%l, 2) == mod(sq%l + sp%l, 2)) .and. (sc%j == sh%j) .and. (sc%l == sh%l)) then
                                        ME = -0.5_dp * real(2*J + 1,dp) / real(sh%j + 1,dp) * O2b_pp(sp%a,sq%a,sh%a,sg%a,J,P_hg,V_NN,Orb,Orb_NN) * &
                                             O2b_pp(sc%a,sg%a,sp%a,sq%a,J,P_pq,V_NN,Orb,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sc%E + sg%E - sp%E - sq%E)
                                        dpRho(sc%a,sh%a) = dpRho(sc%a,sh%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! pRho 2p-2h pn ... particle block ...
        !$omp parallel do schedule(dynamic) reduction(+:dpRho) private(sp,q,sq,P_pq,J,h,sh,g,sg,P_hg,a,sa,ME)
        do p = 1, N_Particle%p
            sp = Particle%p(p)
            do q = 1, N_Particle%n
                sq = Particle%n(q)
                P_pq = mod(sp%l + sq%l,2) + 1
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%p
                        sh = Hole%p(h)
                        do g = 1, N_Hole%n
                            sg = Hole%n(g)
                            P_hg = mod(sh%l + sg%l,2) + 1
                            if ((abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do a = 1, N_Particle%p
                                    sa = Particle%p(a)
                                    if ((mod(sh%l + sg%l, 2) == mod(sq%l + sa%l, 2)) .and. (sa%j == sp%j) .and. (sa%l == sp%l)) then
                                        ME = real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pn(sp%a,sq%a,sh%a,sg%a,J,P_pq,V_NN,Orb_NN) * &
                                             O2b_pn(sh%a,sg%a,sa%a,sq%a,J,P_hg,V_NN,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sh%E + sg%E - sa%E - sq%E)
                                        dpRho(sp%a,sa%a) = dpRho(sp%a,sa%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! pRho 2p-2h pn ... hole block ...
        !$omp parallel do schedule(dynamic) reduction(+:dpRho) private(sh,g,sg,J,p,sp,q,sq,P_pq,c,sc,ME)
        do h = 1, N_Hole%p
            sh = Hole%p(h)
            do g = 1, N_Hole%n
                sg = Hole%n(g)
                do J = abs(sh%j - sg%j)/2, (sh%j + sg%j)/2
                    do p = 1, N_Particle%p
                        sp = Particle%p(p)
                        do q = 1, N_Particle%n
                            sq = Particle%n(q)
                            P_pq = mod(sp%l + sq%l,2) + 1
                            if ((abs(sp%j - sq%j) <= 2*J) .and. (2*J <= (sp%j + sq%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do c = 1, N_Hole%p
                                    sc = Hole%p(c)
                                    if ((mod(sg%l + sc%l, 2) == mod(sq%l + sp%l, 2)) .and. (sc%j == sh%j) .and. (sc%l == sh%l)) then
                                        ME = -real(2*J + 1,dp) / real(sh%j + 1,dp) * O2b_pn(sp%a,sq%a,sh%a,sg%a,J,P_pq,V_NN,Orb_NN) * &
                                             O2b_pn(sc%a,sg%a,sp%a,sq%a,J,P_pq,V_NN,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sc%E + sg%E - sp%E - sq%E)
                                        dpRho(sc%a,sh%a) = dpRho(sc%a,sh%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! nRho 1p-1h nn ...
        !$omp parallel do schedule(dynamic) reduction(+:dnRho) private(sp,h,sh,f,sf,P_pf,P_hf,J,r,sr,q,sq,ME,g,sg,P_pr,P_hr)
        do p = 1, N_Particle%n
            sp = Particle%n(p)
            do h = 1, N_Hole%n
                sh = Hole%n(h)
                if (sp%j == sh%j .and. sp%l == sh%l) then
                    do f = 1, N_Hole%n
                        sf = Hole%n(f)
                        P_pf = mod(sp%l + sf%l,2) + 1
                        P_hf = mod(sh%l + sf%l,2) + 1
                        do J = abs(sp%j - sf%j)/2, (sp%j + sf%j)/2
                            do r = 1, N_Particle%n
                                sr = Particle%n(r)
                                do q = 1, N_Particle%n
                                    sq = Particle%n(q)
                                    if ((abs(sq%j - sr%j) <= 2*J) .and. (2*J <= (sq%j + sr%j)) .and. (mod(sq%l + sr%l, 2) == mod(sh%l + sf%l, 2))) then
                                        ME = 0.5_dp * real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_nn(sp%a,sf%a,sq%a,sr%a,J,P_pf,V_NN,Orb,Orb_NN) * &
                                             O2b_nn(sq%a,sr%a,sh%a,sf%a,J,P_hf,V_NN,Orb,Orb_NN) / (sh%E - sp%E) / (sh%E + sf%E - sq%E - sr%E)
                                        dnRho(sp%a,sh%a) = dnRho(sp%a,sh%a) + ME
                                        dnRho(sh%a,sp%a) = dnRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do

                    do r = 1, N_Particle%n
                        sr = Particle%n(r)
                        P_pr = mod(sp%l + sr%l,2) + 1
                        P_hr = mod(sh%l + sr%l,2) + 1
                        do J = abs(sh%j - sr%j)/2, (sh%j + sr%j)/2
                            do f = 1, N_Hole%n
                                sf = Hole%n(f)
                                do g = 1, N_Hole%n
                                    sg = Hole%n(g)
                                    if ((abs(sg%j - sf%j) <= 2*J) .and. (2*J <= (sg%j + sf%j)) .and. (mod(sp%l + sr%l, 2) == mod(sg%l + sf%l, 2))) then
                                        ME = -0.5_dp * real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_nn(sp%a,sr%a,sg%a,sf%a,J,P_pr,V_NN,Orb,Orb_NN) * &
                                             O2b_nn(sg%a,sf%a,sh%a,sr%a,J,P_hr,V_NN,Orb,Orb_NN) / (sh%E - sp%E) / (sg%E + sf%E - sp%E - sr%E)
                                        dnRho(sp%a,sh%a) = dnRho(sp%a,sh%a) + ME
                                        dnRho(sh%a,sp%a) = dnRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do
                end if
            end do
        end do
        !$omp end parallel do

        ! nRho 1p-1h pn ...
        !$omp parallel do schedule(dynamic) reduction(+:dnRho) private(sp,h,sh,f,sf,P_pf,P_hf,J,r,sr,q,sq,ME,g,sg,P_pr,P_hr)
        do p = 1, N_Particle%n
            sp = Particle%n(p)
            do h = 1, N_Hole%n
                sh = Hole%n(h)
                if (sp%j == sh%j .and. sp%l == sh%l) then
                    do f = 1, N_Hole%p
                        sf = Hole%p(f)
                        P_pf = mod(sp%l + sf%l,2) + 1
                        P_hf = mod(sh%l + sf%l,2) + 1
                        do J = abs(sp%j - sf%j)/2, (sp%j + sf%j)/2
                            do r = 1, N_Particle%p
                                sr = Particle%p(r)
                                do q = 1, N_Particle%n
                                    sq = Particle%n(q)
                                    if ((abs(sq%j - sr%j) <= 2*J) .and. (2*J <= (sq%j + sr%j)) .and. (mod(sq%l + sr%l, 2) == mod(sh%l + sf%l, 2))) then
                                        ME = real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pn(sf%a,sp%a,sr%a,sq%a,J,P_pf,V_NN,Orb_NN) * &
                                             O2b_pn(sr%a,sq%a,sf%a,sh%a,J,P_hf,V_NN,Orb_NN) / (sh%E - sp%E) / (sh%E + sf%E - sq%E - sr%E)
                                        dnRho(sp%a,sh%a) = dnRho(sp%a,sh%a) + ME
                                        dnRho(sh%a,sp%a) = dnRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do

                    do r = 1, N_Particle%p
                        sr = Particle%p(r)
                        P_pr = mod(sp%l + sr%l,2) + 1
                        P_hr = mod(sh%l + sr%l,2) + 1
                        do J = abs(sh%j - sr%j)/2, (sh%j + sr%j)/2
                            do f = 1, N_Hole%p
                                sf = Hole%p(f)
                                do g = 1, N_Hole%n
                                    sg = Hole%n(g)
                                    if ((abs(sg%j - sf%j) <= 2*J) .and. (2*J <= (sg%j + sf%j)) .and. (mod(sp%l + sr%l, 2) == mod(sg%l + sf%l, 2))) then
                                        ME = -real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pn(sr%a,sp%a,sf%a,sg%a,J,P_pr,V_NN,Orb_NN) * &
                                             O2b_pn(sf%a,sg%a,sr%a,sh%a,J,P_hr,V_NN,Orb_NN) / (sh%E - sp%E) / (sg%E + sf%E - sp%E - sr%E)
                                        dnRho(sp%a,sh%a) = dnRho(sp%a,sh%a) + ME
                                        dnRho(sh%a,sp%a) = dnRho(sh%a,sp%a) + ME
                                    end if
                                end do
                            end do
                        end do
                    end do
                end if
            end do
        end do
        !$omp end parallel do

        ! nRho 2p-2h nn ... particle block ...
        !$omp parallel do schedule(dynamic) reduction(+:dnRho) private(sp,q,sq,P_pq,J,h,sh,g,sg,P_hg,a,sa,ME)
        do p = 1, N_Particle%n
            sp = Particle%n(p)
            do q = 1, N_Particle%n
                sq = Particle%n(q)
                P_pq = mod(sp%l + sq%l,2) + 1
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%n
                        sh = Hole%n(h)
                        do g = 1, N_Hole%n
                            sg = Hole%n(g)
                            P_hg = mod(sh%l + sg%l,2) + 1
                            if ((abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do a = 1, N_Particle%n
                                    sa = Particle%n(a)
                                    if ((mod(sh%l + sg%l, 2) == mod(sq%l + sa%l, 2)) .and. (sa%j == sp%j) .and. (sa%l == sp%l)) then
                                        ME = 0.5_dp * real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_nn(sp%a,sq%a,sh%a,sg%a,J,P_pq,V_NN,Orb,Orb_NN) * &
                                             O2b_nn(sh%a,sg%a,sa%a,sq%a,J,P_hg,V_NN,Orb,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sh%E + sg%E - sa%E - sq%E)
                                        dnRho(sp%a,sa%a) = dnRho(sp%a,sa%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! nRho 2p-2h nn ... hole block ...
        !$omp parallel do schedule(dynamic) reduction(+:dnRho) private(sh,g,sg,J,p,sp,q,sq,P_pq,c,sc,ME)
        do h = 1, N_Hole%n
            sh = Hole%n(h)
            do g = 1, N_Hole%n
                sg = Hole%n(g)
                do J = abs(sh%j - sg%j)/2, (sh%j + sg%j)/2
                    do p = 1, N_Particle%n
                        sp = Particle%n(p)
                        do q = 1, N_Particle%n
                            sq = Particle%n(q)
                            P_pq = mod(sp%l + sq%l,2) + 1
                            if ((abs(sp%j - sq%j) <= 2*J) .and. (2*J <= (sp%j + sq%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do c = 1, N_Hole%n
                                    sc = Hole%n(c)
                                    if ((mod(sg%l + sc%l, 2) == mod(sq%l + sp%l, 2)) .and. (sc%j == sh%j) .and. (sc%l == sh%l)) then
                                        ME = -0.5_dp * real(2*J + 1,dp) / real(sh%j + 1,dp) * O2b_nn(sp%a,sq%a,sh%a,sg%a,J,P_pq,V_NN,Orb,Orb_NN) * &
                                             O2b_nn(sc%a,sg%a,sp%a,sq%a,J,P_pq,V_NN,Orb,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sc%E + sg%E - sp%E - sq%E)
                                        dnRho(sc%a,sh%a) = dnRho(sc%a,sh%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! nRho 2p-2h pn ... particle block ...
        !$omp parallel do schedule(dynamic) reduction(+:dnRho) private(sp,q,sq,J,h,sh,g,sg,P_hg,a,sa,ME)
        do p = 1, N_Particle%n
            sp = Particle%n(p)
            do q = 1, N_Particle%p
                sq = Particle%p(q)
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%n
                        sh = Hole%n(h)
                        do g = 1, N_Hole%p
                            sg = Hole%p(g)
                            P_hg = mod(sh%l + sg%l,2) + 1
                            if ((abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do a = 1, N_Particle%n
                                    sa = Particle%n(a)
                                    if ((mod(sh%l + sg%l, 2) == mod(sq%l + sa%l, 2)) .and. (sa%j == sp%j) .and. (sa%l == sp%l)) then
                                        ME = real(2*J + 1,dp) / real(sp%j + 1,dp) * O2b_pn(sq%a,sp%a,sg%a,sh%a,J,P_hg,V_NN,Orb_NN) * &
                                             O2b_pn(sg%a,sh%a,sq%a,sa%a,J,P_hg,V_NN,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sh%E + sg%E - sa%E - sq%E)
                                        dnRho(sp%a,sa%a) = dnRho(sp%a,sa%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! nRho 2p-2h pn ... hole block ...
        !$omp parallel do schedule(dynamic) reduction(+:dnRho) private(sh,g,sg,J,p,sp,q,sq,P_pq,c,sc,ME)
        do h = 1, N_Hole%n
            sh = Hole%n(h)
            do g = 1, N_Hole%p
                sg = Hole%p(g)
                do J = abs(sh%j - sg%j)/2, (sh%j + sg%j)/2
                    do p = 1, N_Particle%n
                        sp = Particle%n(p)
                        do q = 1, N_Particle%p
                            sq = Particle%p(q)
                            P_pq = mod(sp%l + sq%l,2) + 1
                            if ((abs(sp%j - sq%j) <= 2*J) .and. (2*J <= (sp%j + sq%j)) .and. (mod(sh%l + sg%l, 2) == mod(sp%l + sq%l, 2))) then
                                do c = 1, N_Hole%n
                                    sc = Hole%n(c)
                                    if ((mod(sg%l + sc%l, 2) == mod(sq%l + sp%l, 2)) .and. (sc%j == sh%j) .and. (sc%l == sh%l)) then
                                        ME = -real(2*J + 1,dp) / real(sh%j + 1,dp) * O2b_pn(sq%a,sp%a,sg%a,sh%a,J,P_pq,V_NN,Orb_NN) * &
                                             O2b_pn(sg%a,sc%a,sq%a,sp%a,J,P_pq,V_NN,Orb_NN) / (sh%E + sg%E - sp%E - sq%E) / &
                                             (sc%E + sg%E - sp%E - sq%E)
                                        dnRho(sc%a,sh%a) = dnRho(sc%a,sh%a) + ME
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! HF transformation matrix C (LHO -> HF) & HF-basis OBDM ...
        C_HF = O1b_import(Params,Orb,io_path(Params,'/HF/bin/C.bin'))
        Rho_HF = O1b_import(Params,Orb,io_path(Params,'/HF/bin/rho.bin'))

        ! HF-MBPT corrections ...
        Rho%p = Rho_HF%p + dpRho
        Rho%n = Rho_HF%n + dnRho

        ! HF-MBPT natural orbitals (Hungarian ordering, positive diagonal phases) ...
        call OBDM_natural_orbitals(Params,Orb,Rho,Orb_NAT,D,Rho_NAT)

        ! Radial HF-MBPT densities & radii ... Rho in the LHO basis (C D Rho_NAT D^T C^T) & the LHO -> NAT matrix C D ...
        C_NAT%p = mm(C_HF%p,D%p)
        C_NAT%n = mm(C_HF%n,D%n)
        Rho_LHO%p = mmt(mm(mm(C_HF%p,D%p),Rho_NAT%p),mm(C_HF%p,D%p))
        Rho_LHO%n = mmt(mm(mm(C_HF%n,D%n),Rho_NAT%n),mm(C_HF%n,D%n))
        call rho_export(Params,Orb,io_path(Params,'/HF/HF_summary.dat'), &
                        io_path(Params,'/HF/densities/HF_MBPT_radial_densities.dat'),Rho_LHO,C_NAT)

        call HF_MBPT_occupation(Params,Rho_HF,Rho_NAT)

        ! Natural orbitals, C_NAT (LHO -> NAT) & the NAT-basis OBDM for HF_MBPT_export (temporary files as in dev) ...
        if (Params%Calc%HF%MBPT%Export) then
            call orbitals_export(Params,Orb_NAT,io_path(Params,'/HF/bin/orb_temp.bin'))
            call O1b_export(Params,Orb_NAT,C_NAT,io_path(Params,'/HF/bin/C_temp.bin'))
            call O1b_export(Params,Orb_NAT,Rho_NAT,io_path(Params,'/HF/bin/rho_temp.bin'))
        end if
    end subroutine HF_MBPT_OBDM

    ! HF & HF-MBPT (natural-orbital) occupations ... io/<Path>/HF/densities/HF_MBPT_occupation.dat ...
    subroutine HF_MBPT_occupation(Params,Rho_HF,Rho_NAT)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho_HF, Rho_NAT
        integer :: a_max, a, u

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        open(newunit=u, file=io_path(Params,'/HF/densities/HF_MBPT_occupation.dat'), status='replace', action='write')
        write(u,'(a)') cfmt_s('a',6,.true.) // ' ' // cfmt_s('pN_HF',12,.true.) // ' ' // cfmt_s('pN_PT',12,.true.) // ' ' // &
                       cfmt_s('d_pN',12,.true.) // ' ' // cfmt_s('nN_HF',12,.true.) // ' ' // cfmt_s('nN_PT',12,.true.) // ' ' // &
                       cfmt_s('d_nN',12,.true.)
        do a = 1, a_max
            write(u,'(a)') cfmt_d(a,6,.true.) // ' ' // cfmt_f(Rho_HF%p(a,a),12,6,.true.) // ' ' // &
                           cfmt_f(Rho_NAT%p(a,a),12,6,.true.) // ' ' // cfmt_f(Rho_HF%p(a,a) - Rho_NAT%p(a,a),12,6,.true.) // ' ' // &
                           cfmt_f(Rho_HF%n(a,a),12,6,.true.) // ' ' // cfmt_f(Rho_NAT%n(a,a),12,6,.true.) // ' ' // &
                           cfmt_f(Rho_HF%n(a,a) - Rho_NAT%n(a,a),12,6,.true.)
        end do
        close(u)
    end subroutine HF_MBPT_occupation

    ! dev HF_MBPT_occupation_check ... HF-MBPT needs an idempotent HF density (diagonal of the HF-basis Rho) ...
    logical function HF_MBPT_occupation_check(Params,Rho_HF) result(OK)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho_HF
        real(dp), parameter :: eps = 1.0e-5_dp
        real(dp) :: pRho_a, nRho_a
        integer :: a, a_max

        a_max = (Params%Calc%Nmax + 1) * (Params%Calc%Nmax + 2) / 2
        OK = .true.
        do a = 1, a_max
            pRho_a = Rho_HF%p(a,a)
            nRho_a = Rho_HF%n(a,a)
            if ((abs(pRho_a**2 - pRho_a) > eps) .or. (abs(nRho_a**2 - nRho_a) > eps)) then
                write(*,'(/,a)') 'HF-MBPT calculation will not be performed due to non-idempotent HF solution ...'
                OK = .false.
                return
            end if
        end do
    end function HF_MBPT_occupation_check

    ! dev HF_MBPT_export(Params,Orb,Orb_NN,Orb_NNN,C,T,h_N,V_NN,V_NN_res,V_NNN) (2223-2327) ... V_NN bare, V_NN_res residual
        ! (HF basis) ... correlated reference in the HF-MBPT natural-orbital basis: io/<Path>/HF-MBPT/bin/{orb,C,H1b,rho,H2b}.bin
        ! & sigma.bin (TBCF) ...
    subroutine HF_MBPT_export(Params,Orb,Orb_NN,Orb_NNN,C,T,h_N,V_NN,V_NN_res,V_NNN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: C, T, h_N
        type(o2b_t), intent(in) :: V_NN, V_NN_res
        type(v3b_t), intent(in) :: V_NNN
        type(orb1b_hf_mbpt_t) :: Orb_HF_MBPT
        type(orb1b_t), allocatable :: Orb_NAT(:)
        type(o1b_t) :: C_NAT, Rho_NAT, h_N_NAT, C_HF_NAT
        type(o2b_t) :: Sigma_NAT, V_NN_NAT
        logical :: TBCF

        TBCF = Params%Calc%HF%MBPT%TBCF

        Orb_HF_MBPT = HF_MBPT_orbitals_make(Params,Orb)

        write(*,'(/,a)') 'Preparing export of perturbatively correlated HF-MBPT ground state ...'

        ! Fresh io/<Path>/HF-MBPT ...
        if (is_dir(io_path(Params,'/HF-MBPT'))) call remove_dir(io_path(Params,'/HF-MBPT'))
        call make_dir(io_path(Params,'/HF-MBPT'))
        call make_dir(io_path(Params,'/HF-MBPT/bin'))

        write(*,'(/,a)') TAB // 'Preparing HF-MBPT natural orbitals & 1-body density matrix ...'

        ! Natural orbitals from the temporary files of HF_MBPT_OBDM (read with the (l,j) blocks of Orb, as dev) ...
        Orb_NAT = orbitals_import(Params,io_path(Params,'/HF/bin/orb_temp.bin'))
        C_NAT = O1b_import(Params,Orb,io_path(Params,'/HF/bin/C_temp.bin'))
        Rho_NAT = O1b_import(Params,Orb,io_path(Params,'/HF/bin/rho_temp.bin'))

        ! Rho -> LHO basis (C_NAT Rho C_NAT^T) ...
        Rho_NAT = O1b_transformation(Rho_NAT,O1b_transpose(C_NAT))

        write(*,'(/,a)') TAB // TAB // 'HF-MBPT natural orbitals & 1-body density matrix prepared ...'

        ! 2-body correlation function ... evaluated in the HF basis & returned in the LHO basis ...
        if (TBCF) then
            write(*,'(/,a)') TAB // 'Constructing HF-MBPT 2-body correlation function ...'
            call HF_MBPT_TBCF(Params,C,h_N,V_NN_res,Orb,Orb_HF_MBPT,Orb_NN,Sigma_NAT,Report=.false.)
            write(*,'(/,a)') TAB // TAB // 'HF-MBPT 2-body correlation function constructed ...'
        end if

        write(*,'(/,a)') TAB // 'Constructing HF-MBPT 1- & 2-body Hamiltonian ...'

        if (Params%Calc%HF%MBPT%ScH) then
            ! Polarized Hamiltonian of the correlated ground state ...
            write(*,'(a)') TAB // TAB // TAB // 'Constructing polarized HF-MBPT Hamiltonian in the NAT basis...'
            if (TBCF) then
                call H1b_H2b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho_NAT,Sigma_NAT,h_N_NAT,V_NN_NAT,Report=.false.)
            else
                h_N_NAT = H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho_NAT)
                call H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho_NAT,V_NN_NAT,Report=.false.)
            end if
            h_N_NAT = O1b_transformation(h_N_NAT,C_NAT)
            call O2b_transformation(Params,Orb,Orb_NN,V_NN_NAT,C_NAT,Report=.false.)
        else
            ! "HF Hamiltonian" in the NAT basis ... h from the correlated Rho_NAT as in dev (UPSTREAM_BUGS), V = V_NN_res
                ! transformed HF -> NAT with C^T C_NAT ...
            write(*,'(a)') TAB // TAB // TAB // 'Constructing HF Hamiltonian in the NAT basis ...'
            h_N_NAT = H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho_NAT)
            V_NN_NAT = V_NN_res
            h_N_NAT = O1b_transformation(h_N_NAT,C_NAT)
            C_HF_NAT%p = mtm(C%p,C_NAT%p)
            C_HF_NAT%n = mtm(C%n,C_NAT%n)
            call O2b_transformation(Params,Orb,Orb_NN,V_NN_NAT,C_HF_NAT,Report=.false.)
        end if

        write(*,'(/,a)') TAB // TAB // 'Finished construction of HF-MBPT 1- & 2-body Hamiltonian ...'

        ! Rho & Sigma -> NAT basis ...
        Rho_NAT = O1b_transformation(Rho_NAT,C_NAT)
        if (TBCF) call O2b_transformation(Params,Orb,Orb_NN,Sigma_NAT,C_NAT,Report=.false.)

        ! Remove the temporary files ...
        call remove_file(io_path(Params,'/HF/bin/orb_temp.bin'))
        call remove_file(io_path(Params,'/HF/bin/C_temp.bin'))
        call remove_file(io_path(Params,'/HF/bin/rho_temp.bin'))

        ! Export ...
        call orbitals_export(Params,Orb_NAT,io_path(Params,'/HF-MBPT/bin/orb.bin'))
        call O1b_export(Params,Orb_NAT,C_NAT,io_path(Params,'/HF-MBPT/bin/C.bin'))
        call O1b_export(Params,Orb_NAT,h_N_NAT,io_path(Params,'/HF-MBPT/bin/H1b.bin'))
        call O1b_export(Params,Orb_NAT,Rho_NAT,io_path(Params,'/HF-MBPT/bin/rho.bin'))
        call O2b_export(Params,Orb_NN,V_NN_NAT,io_path(Params,'/HF-MBPT/bin/H2b.bin'))
        if (TBCF) call O2b_export(Params,Orb_NN,Sigma_NAT,io_path(Params,'/HF-MBPT/bin/sigma.bin'))

        write(*,'(/,a)') TAB // TAB // 'Finished performing export & garbage collection for HF-MBPT perturbatively correlated ' // &
                         'ground state ...'
    end subroutine HF_MBPT_export

end module nj_hf_mbpt
