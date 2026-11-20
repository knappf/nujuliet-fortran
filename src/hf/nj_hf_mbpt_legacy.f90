! LEGACY (removed in cleanup stage) ... base-version HF-MBPT(2) of the legacy HF (HF_LEGACY), reads io/<Path>/Bin/* ...
    ! counterpart of base Source/HF/HF_MBPT.jl ... the dev HF-MBPT is nj_hf_mbpt ...
module nj_hf_mbpt_legacy
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_orbitals, only : orbitals_ph_make
    use nj_o1b, only : O1b_import
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_obdm, only : OBDM_export
    use nj_hf_legacy, only : HF_density_operator
    implicit none
    private

    public :: HF_MBPT, HF_MBPT_energy_correction, HF_MBPT_OBDM, HF_MBPT_OBDM_correction, HF_MBPT_reorder

contains

    subroutine HF_MBPT(Params,Orb,Orb_NN,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: V_NN
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        integer(i8) :: t0

        call orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)

        call timer_start(t0)
        call HF_MBPT_energy_correction(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN)
        call timer_report(t0)

        call timer_start(t0)
        call HF_MBPT_OBDM_correction(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN)
        call timer_report(t0)
    end subroutine HF_MBPT

    ! Leading-order MBPT(2) correction to the ground-state energy ...
    subroutine HF_MBPT_energy_correction(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(o2b_t), intent(in) :: V_NN
        real(dp) :: dE, Sum
        integer :: p, q, h, g, J, Par, u
        type(spstate_t) :: sp, sq, sh, sg

        write(*,'(/,a)') 'Calculating HF-MBPT(2) ground-state energy correction ... '

        ! proton-proton contribution ...
        Sum = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:Sum) private(sp,q,sq,Par,J,h,sh,g,sg)
        do p = 1, N_Particle%p
            sp = Particle%p(p)
            do q = 1, N_Particle%p
                sq = Particle%p(q)
                Par = mod(sp%l + sq%l, 2) + 1
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%p
                        sh = Hole%p(h)
                        do g = 1, N_Hole%p
                            sg = Hole%p(g)
                            if ((Par == (mod(sh%l + sg%l, 2) + 1)) .and. (abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j))) then
                                Sum = Sum + real(2*J + 1,dp) / 4.0_dp * abs(O2b_pp(sp%a,sq%a,sh%a,sg%a,J,Par,V_NN,Orb,Orb_NN))**2 / &
                                      (sh%E + sg%E - sp%E - sq%E)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        dE = Sum

        ! neutron-neutron contribution ...
        Sum = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:Sum) private(sp,q,sq,Par,J,h,sh,g,sg)
        do p = 1, N_Particle%n
            sp = Particle%n(p)
            do q = 1, N_Particle%n
                sq = Particle%n(q)
                Par = mod(sp%l + sq%l, 2) + 1
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%n
                        sh = Hole%n(h)
                        do g = 1, N_Hole%n
                            sg = Hole%n(g)
                            if ((Par == (mod(sh%l + sg%l, 2) + 1)) .and. (abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j))) then
                                Sum = Sum + real(2*J + 1,dp) / 4.0_dp * abs(O2b_nn(sp%a,sq%a,sh%a,sg%a,J,Par,V_NN,Orb,Orb_NN))**2 / &
                                      (sh%E + sg%E - sp%E - sq%E)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        dE = dE + Sum

        ! proton-neutron contribution ...
        Sum = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:Sum) private(sp,q,sq,Par,J,h,sh,g,sg)
        do p = 1, N_Particle%p
            sp = Particle%p(p)
            do q = 1, N_Particle%n
                sq = Particle%n(q)
                Par = mod(sp%l + sq%l, 2) + 1
                do J = abs(sp%j - sq%j)/2, (sp%j + sq%j)/2
                    do h = 1, N_Hole%p
                        sh = Hole%p(h)
                        do g = 1, N_Hole%n
                            sg = Hole%n(g)
                            if ((Par == (mod(sh%l + sg%l, 2) + 1)) .and. (abs(sh%j - sg%j) <= 2*J) .and. (2*J <= (sh%j + sg%j))) then
                                Sum = Sum + real(2*J + 1,dp) * abs(O2b_pn(sp%a,sq%a,sh%a,sg%a,J,Par,V_NN,Orb_NN))**2 / &
                                      (sh%E + sg%E - sp%E - sq%E)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        dE = dE + Sum

        write(*,'(a)') TAB // 'HF-MBPT(2) energy    ...   E^(2) = ' // jstr(round_sig(dE,9)) // ' MeV'

        open(newunit=u, file=io_path(Params,'/HF/HF_Summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Spherical Hartree-Fock Leading Order Many-Body Perturbation Theory solution review:'
        write(u,'(/,a)',advance='no') 'E_0^(2)     = ' // cfmt_f(dE,15,8) // ' ' // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // &
                                      ' LO HF-MBPT(2) correction to ground-state energy'
        write(u,'(/,a)') 'E_0^(2) / A = ' // cfmt_f(dE / real(Params%Calc%A,dp),15,8) // ' ' // TAB // ' MeV ' // TAB // &
                                      TAB // ' ... ' // TAB // ' LO HF-MBPT(2) correction to ground-state energy per nucleon'
        close(u)
    end subroutine HF_MBPT_energy_correction

    ! MBPT(3) leading-order corrections to the 1-body density matrix in the HF basis ...
    subroutine HF_MBPT_OBDM(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN,dpRho,dnRho)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(o2b_t), intent(in) :: V_NN
        real(dp), allocatable, intent(out) :: dpRho(:,:), dnRho(:,:)
        integer :: a_max, p, h, f, J, r, q, g, a, c, P_pf, P_hf, P_pr, P_hr, P_pq, P_hg
        type(spstate_t) :: sp, sh, sf, sr, sq, sg, sa, sc
        real(dp) :: ME

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
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
    end subroutine HF_MBPT_OBDM

    subroutine HF_MBPT_OBDM_correction(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(pninteger_t), intent(in) :: N_Particle, N_Hole
        type(pnsvector_t), intent(in) :: Particle, Hole
        type(o2b_t), intent(in) :: V_NN
        real(dp), allocatable :: dpRho(:,:), dnRho(:,:), pRho(:,:), nRho(:,:)
        type(o1b_t) :: C, Rho_HF, Rho, D, Rho_LHO, C_NAT
        type(pnvector_t) :: Occ
        integer :: a_max, a, b

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(/,a)') 'Evaluating density LO corrections ...'

        call HF_MBPT_OBDM(Params,Orb,Orb_NN,N_Particle,Particle,N_Hole,Hole,V_NN,dpRho,dnRho)

        ! Transformation matrix C ... LHO -> HF ...
        C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))

        ! HF OBDM in the HF basis ...
        Rho_HF = HF_density_operator(a_max,C,Orb)
        Rho_HF%p = mtmm(C%p,Rho_HF%p)
        Rho_HF%n = mtmm(C%n,Rho_HF%n)

        ! Include HF-MBPT(3) corrections ...
        pRho = Rho_HF%p + dpRho
        nRho = Rho_HF%n + dnRho

        ! Symmetrize, clean the (l,j) block structure & lift degeneracies deterministically ...
        pRho = 0.5_dp * (pRho + transpose(pRho))
        nRho = 0.5_dp * (nRho + transpose(nRho))
        do a = 1, a_max
            do b = 1, a_max
                if ((Orb(a)%j /= Orb(b)%j) .or. (Orb(a)%l /= Orb(b)%l)) then
                    pRho(a,b) = 0.0_dp
                    nRho(a,b) = 0.0_dp
                end if
            end do
        end do
        do a = 1, a_max
            pRho(a,a) = pRho(a,a) + 1.0e-10_dp * real(a,dp)
            nRho(a,a) = nRho(a,a) + 1.0e-10_dp * real(a,dp)
        end do

        allocate(Occ%p(a_max), Occ%n(a_max), D%p(a_max,a_max), D%n(a_max,a_max))
        call eig_sym(pRho,Occ%p,D%p)
        call eig_sym(nRho,Occ%n,D%n)

        ! Reorder D & make the diagonal OBDM in the HF-NAT basis ...
        call HF_MBPT_reorder(Params,Occ,D,Rho)

        ! Radial HF-MBPT(3) densities & radii ...
        C_NAT%p = mm(C%p,D%p)
        C_NAT%n = mm(C%n,D%n)
        Rho_LHO%p = mmt(mm(C_NAT%p,Rho%p),C_NAT%p)
        Rho_LHO%n = mmt(mm(C_NAT%n,Rho%n),C_NAT%n)
        call OBDM_export(Params,Orb,io_path(Params,'/HF/HF_Summary.dat'), &
                         io_path(Params,'/HF/Densities/HFMBPT_Radial_Densities.dat'),Rho_LHO,C_NAT)

        call HF_MBPT_occupation(Params,Rho_HF,Rho)
    end subroutine HF_MBPT_OBDM_correction

    ! Orders natural orbitals by their largest overlap with the HF orbitals ... returns diag(n) ...
    subroutine HF_MBPT_reorder(Params,Occ,C,Rho)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(inout) :: Occ
        type(o1b_t), intent(inout) :: C
        type(o1b_t), intent(out) :: Rho
        integer :: a_max, a, b, pIndex, nIndex
        integer, allocatable :: pOrdering(:), nOrdering(:)
        real(dp) :: pMaximum, nMaximum, pOverlap, nOverlap

        a_max = (Params%Calc%Nmax + 1) * (Params%Calc%Nmax + 2) / 2
        allocate(pOrdering(a_max), nOrdering(a_max))
        do a = 1, a_max
            pIndex = 0
            nIndex = 0
            pMaximum = 0.0_dp
            nMaximum = 0.0_dp
            do b = 1, a_max
                ! Overlap of the unit vector e_a with the column b ...
                pOverlap = abs(C%p(a,b))
                nOverlap = abs(C%n(a,b))
                if (pOverlap > pMaximum) then
                    pIndex = b
                    pMaximum = pOverlap
                end if
                if (nOverlap > nMaximum) then
                    nIndex = b
                    nMaximum = nOverlap
                end if
            end do
            pOrdering(a) = pIndex
            nOrdering(a) = nIndex
        end do

        Occ%p = Occ%p(pOrdering)
        C%p = C%p(:,pOrdering)
        Occ%n = Occ%n(nOrdering)
        C%n = C%n(:,nOrdering)

        Rho%p = diagm(Occ%p)
        Rho%n = diagm(Occ%n)
    end subroutine HF_MBPT_reorder

    subroutine HF_MBPT_occupation(Params,Rho_HF,Rho_NAT)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho_HF, Rho_NAT
        integer :: a_max, a, u

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        open(newunit=u, file=io_path(Params,'/HF/Densities/HFMBPT_Occupation.dat'), status='replace', action='write')
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

end module nj_hf_mbpt_legacy
