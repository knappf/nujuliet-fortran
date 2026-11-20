! Spherical Equal-Filling Approximation HF ... counterpart of HF/HF_EFA.jl ...
module nj_hf_efa
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_linalg
    use nj_o1b, only : T1b, T1b_energy
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn, O2b_transformation, O2b_move
    use nj_v2b, only : V2b_read, V2b_residual_no2b
    use nj_obdm, only : OBDM_export
    use nj_hf, only : HF_energy
    implicit none
    private

    public :: HF_EFA

contains

    subroutine HF_EFA(Params,Orb,Orb_NNN,SPE,C,V_NNN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb3b_t), intent(in) :: Orb_NNN
        type(pnvector_t), intent(in) :: SPE
        type(o1b_t), intent(in) :: C
        type(v3b_t), intent(in) :: V_NNN
        type(params_t) :: Params_EFA
        type(o1b_t) :: T, Rho, Rho_HF
        type(o2b_t) :: V_NN, V_NN_res
        type(orb2b_t) :: Orb_NN
        real(dp) :: E_HF, T_HF
        integer(i8) :: t0

        write(*,'(/,a)') 'Starting Spherical Equal-Filling Approximation Hartree-Fock calculation ... '
        write(*,'(/,a)') 'The HF-EFA target nucleus is set to ...'
        write(*,'(a,/)') 'A = ' // istr(Params%Calc%HF%EFA%A) // ', Z = ' // istr(Params%Calc%HF%EFA%Z)
        write(*,'(/,a)') 'The HF-EFA reference nucleus is ...'
        write(*,'(a)') 'A = ' // istr(Params%Calc%A) // ', Z = ' // istr(Params%Calc%Z)

        ! Parameters of the HF-EFA calculation ... other solver settings back at their defaults as in the original ...
        Params_EFA%Int = Params%Int
        Params_EFA%Calc%A = Params%Calc%HF%EFA%A
        Params_EFA%Calc%Z = Params%Calc%HF%EFA%Z
        Params_EFA%Calc%hw = Params%Calc%hw
        Params_EFA%Calc%Nmax = Params%Calc%Nmax
        Params_EFA%Calc%N2max = Params%Calc%N2max
        Params_EFA%Calc%N3max = Params%Calc%N3max
        Params_EFA%Calc%CMS = Params%Calc%CMS
        Params_EFA%Calc%Path = Params%Calc%Path
        Params_EFA%Calc%HF = Params%Calc%HF

        T = T1b(Params_EFA,Orb)
        call V2b_read(Params_EFA,V_NN,Orb_NN)
        Rho = HF_EFA_density_operator(Params_EFA,SPE,C,Orb)
        E_HF = HF_energy(Params_EFA,Rho,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN)
        T_HF = T1b_energy(Params_EFA,Rho,Orb,T)
        call HF_EFA_export_summary(Params_EFA,E_HF,T_HF)

        call OBDM_export(Params_EFA,Orb,io_path(Params,'/HF-EFA/HF-EFA_Summary.dat'), &
                         io_path(Params,'/HF-EFA/Densities/HF_Radial_Densities.dat'),Rho,C)

        Rho_HF%p = mtmm(C%p,Rho%p)
        Rho_HF%n = mtmm(C%n,Rho%n)
        call HF_EFA_export_SPS_summary(Params_EFA,SPE,Rho_HF,Orb)

        if (Params%Calc%HF%EFA%BMF) then
            call timer_start(t0)
            call V2b_residual_no2b(Params_EFA,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)
            call O2b_move(V_NN_res,V_NN)
            call timer_report(t0)
            ! Transformation to the HF basis ... in place, as the Julia O2b_transformation mutates V_NN too ...
            call timer_start(t0)
            call O2b_transformation(Params_EFA,Orb,Orb_NN,V_NN,C)
            call timer_report(t0)
            call HF_EFA_MBPT(Params_EFA,Orb,Orb_NN,SPE,Rho_HF,V_NN)
        end if

        write(*,'(/,a)') 'The Spherical Equal-Filling Approximation Hartree-Fock calculation finished ... '
    end subroutine HF_EFA

    function HF_EFA_density_operator(Params,SPE,C,Orb) result(Rho)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SPE
        type(o1b_t), intent(in) :: C
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t) :: Rho
        real(dp), allocatable :: pOcc(:), nOcc(:)
        integer, allocatable :: Perm(:)
        integer :: Z, N, a_max, Z_count, N_count, i, a, m_a, k, l
        real(dp) :: pSum, nSum

        Z = Params%Calc%HF%EFA%Z
        N = Params%Calc%HF%EFA%A - Z
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        allocate(pOcc(a_max), nOcc(a_max), Perm(a_max))
        pOcc = 0.0_dp
        nOcc = 0.0_dp
        allocate(Rho%p(a_max,a_max), Rho%n(a_max,a_max))
        Rho%p = 0.0_dp
        Rho%n = 0.0_dp

        Z_count = 0
        call sortperm(SPE%p,Perm)
        do i = 1, a_max
            a = Perm(i)
            do m_a = -Orb(a)%j, Orb(a)%j, 2
                if ((Z_count + 1) <= Z) then
                    pOcc(a) = pOcc(a) + 1.0_dp / real(Orb(a)%j + 1,dp)
                    Z_count = Z_count + 1
                end if
            end do
            if (Z_count == Z) exit
        end do

        N_count = 0
        call sortperm(SPE%n,Perm)
        do i = 1, a_max
            a = Perm(i)
            do m_a = -Orb(a)%j, Orb(a)%j, 2
                if ((N_count + 1) <= N) then
                    nOcc(a) = nOcc(a) + 1.0_dp / real(Orb(a)%j + 1,dp)
                    N_count = N_count + 1
                end if
            end do
            if (N_count == N) exit
        end do

        do k = 1, a_max
            do l = 1, a_max
                if (Orb(k)%l == Orb(l)%l .and. Orb(k)%j == Orb(l)%j) then
                    pSum = 0.0_dp
                    nSum = 0.0_dp
                    do a = 1, a_max
                        if (Orb(a)%l == Orb(k)%l .and. Orb(a)%j == Orb(k)%j) then
                            pSum = pSum + C%p(k,a) * C%p(l,a) * pOcc(a)
                            nSum = nSum + C%n(k,a) * C%n(l,a) * nOcc(a)
                        end if
                    end do
                    Rho%p(k,l) = pSum
                    Rho%n(k,l) = nSum
                end if
            end do
        end do
    end function HF_EFA_density_operator

    subroutine HF_EFA_MBPT(Params,Orb,Orb_NN,SPE,Rho,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(pnvector_t), intent(in) :: SPE
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(in) :: V_NN
        integer, allocatable :: pPart(:), pHole(:), nPart(:), nHole(:)
        integer :: a_max, i, ip, iq, ih, ig, p, q, h, g, J, Par, u
        real(dp) :: dE, Sum, Denom

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call indices(Rho%p,pPart,pHole)
        call indices(Rho%n,nPart,nHole)

        write(*,'(/,a)') 'Calculating HF-EFA-MBPT(2) ground-state correlation energy ... '

        ! proton-proton contribution ...
        Sum = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:Sum) private(p,iq,q,Par,J,ih,h,ig,g,Denom)
        do ip = 1, size(pPart)
            p = pPart(ip)
            do iq = 1, size(pPart)
                q = pPart(iq)
                Par = mod(Orb(p)%l + Orb(q)%l, 2) + 1
                do J = abs(Orb(p)%j - Orb(q)%j)/2, (Orb(p)%j + Orb(q)%j)/2
                    do ih = 1, size(pHole)
                        h = pHole(ih)
                        do ig = 1, size(pHole)
                            g = pHole(ig)
                            Denom = SPE%p(h) + SPE%p(g) - SPE%p(p) - SPE%p(q)
                            if ((Par == (mod(Orb(h)%l + Orb(g)%l, 2) + 1)) .and. (abs(Orb(h)%j - Orb(g)%j) <= 2*J) .and. &
                                (2*J <= (Orb(h)%j + Orb(g)%j)) .and. abs(Denom) > 1.0e-3_dp) then
                                Sum = Sum + real(2*J + 1,dp) / 4.0_dp * abs(O2b_pp(p,q,h,g,J,Par,V_NN,Orb,Orb_NN))**2 / Denom * &
                                      (1.0_dp - Rho%p(p,p)) * (1.0_dp - Rho%p(q,q)) * Rho%p(h,h) * Rho%p(g,g)
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
        !$omp parallel do schedule(dynamic) reduction(+:Sum) private(p,iq,q,Par,J,ih,h,ig,g,Denom)
        do ip = 1, size(nPart)
            p = nPart(ip)
            do iq = 1, size(nPart)
                q = nPart(iq)
                Par = mod(Orb(p)%l + Orb(q)%l, 2) + 1
                do J = abs(Orb(p)%j - Orb(q)%j)/2, (Orb(p)%j + Orb(q)%j)/2
                    do ih = 1, size(nHole)
                        h = nHole(ih)
                        do ig = 1, size(nHole)
                            g = nHole(ig)
                            Denom = SPE%n(h) + SPE%n(g) - SPE%n(p) - SPE%n(q)
                            if ((Par == (mod(Orb(h)%l + Orb(g)%l, 2) + 1)) .and. (abs(Orb(h)%j - Orb(g)%j) <= 2*J) .and. &
                                (2*J <= (Orb(h)%j + Orb(g)%j)) .and. abs(Denom) > 1.0e-3_dp) then
                                Sum = Sum + real(2*J + 1,dp) / 4.0_dp * abs(O2b_nn(p,q,h,g,J,Par,V_NN,Orb,Orb_NN))**2 / Denom * &
                                      (1.0_dp - Rho%n(p,p)) * (1.0_dp - Rho%n(q,q)) * Rho%n(h,h) * Rho%n(g,g)
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
        !$omp parallel do schedule(dynamic) reduction(+:Sum) private(p,iq,q,Par,J,ih,h,ig,g,Denom)
        do ip = 1, size(pPart)
            p = pPart(ip)
            do iq = 1, size(nPart)
                q = nPart(iq)
                Par = mod(Orb(p)%l + Orb(q)%l, 2) + 1
                do J = abs(Orb(p)%j - Orb(q)%j)/2, (Orb(p)%j + Orb(q)%j)/2
                    do ih = 1, size(pHole)
                        h = pHole(ih)
                        do ig = 1, size(nHole)
                            g = nHole(ig)
                            Denom = SPE%p(h) + SPE%n(g) - SPE%p(p) - SPE%n(q)
                            if ((Par == (mod(Orb(h)%l + Orb(g)%l, 2) + 1)) .and. (abs(Orb(h)%j - Orb(g)%j) <= 2*J) .and. &
                                (2*J <= (Orb(h)%j + Orb(g)%j)) .and. abs(Denom) > 1.0e-3_dp) then
                                Sum = Sum + real(2*J + 1,dp) * abs(O2b_pn(p,q,h,g,J,Par,V_NN,Orb_NN))**2 / Denom * &
                                      (1.0_dp - Rho%p(p,p)) * (1.0_dp - Rho%n(q,q)) * Rho%p(h,h) * Rho%n(g,g)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        dE = dE + Sum

        write(*,'(a)') TAB // 'HF-EFA-MBPT(2) energy    ...   E^(2) = ' // jstr(round_sig(dE,9)) // ' MeV'

        open(newunit=u, file=io_path(Params,'/HF-EFA/HF-EFA_Summary.dat'), position='append', action='write')
        write(u,'(/,a)') 'Spherical Equal-Filling Approximation Hartree-Fock Leading Order Many-Body Perturbation Theory solution review:'
        write(u,'(/,a)') 'E_0^(2) = ' // jstr(round_sig(dE,9)) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // &
                         ' LO HF-EFA-MBPT(2) correction to ground-state energy'
        write(u,'(/,a)') 'E_0^(2) / A = ' // jstr(round_sig(dE / real(Params%Calc%HF%EFA%A,dp),9)) // TAB // ' MeV ' // TAB // TAB // &
                         ' ... ' // TAB // ' LO HF-EFA-MBPT(2) correction to ground-state energy per nucleon'
        close(u)

    contains

        ! Particle (occupation /= 1) & hole (occupation /= 0) lists ...
        subroutine indices(R,Part,Hol)
            real(dp), intent(in) :: R(:,:)
            integer, allocatable, intent(out) :: Part(:), Hol(:)
            integer :: a, np, nh
            allocate(Part(count([(abs(R(a,a) - 1.0_dp) > 1.0e-3_dp, a = 1, a_max)])))
            allocate(Hol(count([(abs(R(a,a)) > 1.0e-3_dp, a = 1, a_max)])))
            np = 0
            nh = 0
            do a = 1, a_max
                if (abs(R(a,a) - 1.0_dp) > 1.0e-3_dp) then
                    np = np + 1
                    Part(np) = a
                end if
                if (abs(R(a,a)) > 1.0e-3_dp) then
                    nh = nh + 1
                    Hol(nh) = a
                end if
            end do
        end subroutine indices

    end subroutine HF_EFA_MBPT

    subroutine HF_EFA_export_summary(Params,E_HF,T_HF)
        type(params_t), intent(in) :: Params
        real(dp), intent(in) :: E_HF, T_HF
        integer :: u, N_max, A

        A = Params%Calc%HF%EFA%A
        N_max = Params%Calc%Nmax
        write(*,'(/,a)') 'Exporting HF-EFA calculation summary ...'
        open(newunit=u, file=io_path(Params,'/HF-EFA/HF-EFA_Summary.dat'), position='append', action='write')
        write(u,'(a)') 'Nuclid data:    A = ' // istr(A) // ', Z = ' // istr(Params%Calc%HF%EFA%Z)
        write(u,'(/,a)') 'Calculation data:    hw = ' // jstr(Params%Calc%hw) // ' , N_max = ' // istr(N_max) // &
                         ' , N_2max = ' // istr(Params%Calc%N2max) // ' , N_3max = ' // istr(Params%Calc%N3max) // &
                         ', J-basis size = ' // istr((N_max+1)*(N_max+2)/2) // ', M-basis size = ' // &
                         istr((N_max+1)*(N_max+2)*(N_max+3)/6)
        write(u,'(a)') TAB // TAB // TAB // "Center-of-Mass System correction is set to:   CMS = ''" // trim(Params%Calc%CMS) // "''"
        if (Params%Calc%CMS == 'CMS1+2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Combined 1-body + 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else if (Params%Calc%CMS == 'CMS2B') then
            write(u,'(/,a)') TAB // TAB // TAB // 'Only pure 2-body Center-of-Mass System (CMS) motion correction is set ...'
        else
            write(u,'(/,a)') TAB // TAB // TAB // 'No Center-of-Mass System (CMS) motion correction is set ...'
        end if
        write(u,'(/,a)') 'Spherical Equal-Filling Approximation Hartree-Fock solution review:'
        write(u,'(/,a)') 'E_HF = ' // jstr(round_sig(E_HF,9)) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // &
                         ' Mean-field ground state energy'
        write(u,'(a)') 'T_HF = ' // jstr(round_sig(T_HF,9)) // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // &
                       ' Mean-field total kinetic energy'
        write(u,'(/,a)') 'E_HF / A = ' // jstr(round_sig(E_HF / real(A,dp),9)) // TAB // ' MeV / Nucleon ' // TAB // TAB // &
                         ' ... ' // TAB // ' Mean-field ground state energy'
        write(u,'(a)') 'T_HF / A = ' // jstr(round_sig(T_HF / real(A,dp),9)) // TAB // ' MeV / Nucleon ' // TAB // TAB // &
                       ' ... ' // TAB // ' Mean-field total kinetic energy'
        close(u)
    end subroutine HF_EFA_export_summary

    subroutine HF_EFA_export_SPS_summary(Params,SPE,Rho,Orb)
        type(params_t), intent(in) :: Params
        type(pnvector_t), intent(in) :: SPE
        type(o1b_t), intent(in) :: Rho
        type(orb1b_t), intent(in) :: Orb(:)
        integer :: N_max, a_max, Ind_a, a, u, iso
        integer, allocatable :: Sort_Orbs(:), nr(:,:)
        real(dp) :: O_a, E_a

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(Sort_Orbs(a_max), nr(0:N_max,0:max(Params%Calc%N2max+1,2*N_max+1)))

        write(*,'(/,a)') 'Exporting HF-EFA calculation summary ...'
        open(newunit=u, file=io_path(Params,'/HF-EFA/HF-EFA_Summary.dat'), position='append', action='write')
        do iso = 1, 2
            nr = 0
            if (iso == 1) then
                write(u,'(/,a)') 'Proton single-particle states'
                call sortperm(SPE%p,Sort_Orbs)
            else
                write(u,'(/,a)') 'Neutron single-particle states'
                call sortperm(SPE%n,Sort_Orbs)
            end if
            write(u,'(a)') '____________________________________________________________________'
            write(u,'(a)') cfmt_s('n',4) // ' ' // cfmt_s('l',4) // ' ' // cfmt_s('2j',4) // ' ' // cfmt_s('O_a',8) // ' ' // &
                           cfmt_s('E',14)
            do Ind_a = 1, a_max
                a = Sort_Orbs(Ind_a)
                if (iso == 1) then
                    O_a = abs(Rho%p(a,a))
                    E_a = SPE%p(a)
                else
                    O_a = abs(Rho%n(a,a))
                    E_a = SPE%n(a)
                end if
                write(u,'(a)') cfmt_d(nr(Orb(a)%l,Orb(a)%j),4) // ' ' // cfmt_d(Orb(a)%l,4) // ' ' // cfmt_d(Orb(a)%j,4) // &
                               ' ' // cfmt_f(O_a,8,3) // ' ' // cfmt_f(E_a,14,5) // TAB // 'MeV'
                nr(Orb(a)%l,Orb(a)%j) = nr(Orb(a)%l,Orb(a)%j) + 1
            end do
            write(u,'(a)') '____________________________________________________________________'
        end do
        close(u)
    end subroutine HF_EFA_export_SPS_summary

end module nj_hf_efa
