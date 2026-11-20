! HF-MBPT helpers shared by the energy/OBDM & TBCF modules ... counterpart of HF_MBPT_orbitals_make, HF_MBPT_h1b &
    ! HF_MBPT3_V2b_Pandya of src/HF/HF_MBPT.jl (NuJuliet dev, lines 14-152) ...
    ! Isospin convention of HF_MBPT.jl: t = -1 proton, t = +1 neutron ... single-particle energies are the diagonal of the
    ! HF-basis single-particle Hamiltonian h_N (diag(SPE) from HF_solve) ...
module nj_hf_mbpt_common
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    implicit none
    private

    ! Single-particle state with its isospin label ... the (Orb1B, t) tuples of the dev particle & hole lists ...
    type, public :: ispstate_t
        integer :: a = 0, l = 0, j = 0, t = 0
        real(dp) :: E = 0.0_dp
    end type ispstate_t

    public :: HF_MBPT_orbitals_make, HF_MBPT_states, HF_MBPT_isospin_states, HF_MBPT_h1b, HF_MBPT3_V2b_Pandya, HF_MBPT_V2b

contains

    ! Particle & hole orbitals of the HF reference ... dev HF_MBPT_orbitals_make (14-76): proton hole if |pO| > 1e-8,
        ! proton particle if |1 - pO| > 1e-8 (same for neutrons), in index order ...
    function HF_MBPT_orbitals_make(Params,Orb) result(Orb_HF_MBPT)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb1b_hf_mbpt_t) :: Orb_HF_MBPT
        integer :: a_max, a, Pass
        integer :: pN_particle, pN_hole, nN_particle, nN_hole

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2

        ! Pass 1 counts, pass 2 fills ...
        do Pass = 1, 2
            pN_particle = 0
            pN_hole = 0
            nN_particle = 0
            nN_hole = 0
            do a = 1, a_max
                if (abs(1.0_dp - Orb(a)%pO) > 1.0e-8_dp) then
                    pN_particle = pN_particle + 1
                    if (Pass == 2) Orb_HF_MBPT%p%p(pN_particle) = Orb(a)
                end if
                if (abs(Orb(a)%pO) > 1.0e-8_dp) then
                    pN_hole = pN_hole + 1
                    if (Pass == 2) Orb_HF_MBPT%h%p(pN_hole) = Orb(a)
                end if
                if (abs(1.0_dp - Orb(a)%nO) > 1.0e-8_dp) then
                    nN_particle = nN_particle + 1
                    if (Pass == 2) Orb_HF_MBPT%p%n(nN_particle) = Orb(a)
                end if
                if (abs(Orb(a)%nO) > 1.0e-8_dp) then
                    nN_hole = nN_hole + 1
                    if (Pass == 2) Orb_HF_MBPT%h%n(nN_hole) = Orb(a)
                end if
            end do
            if (Pass == 1) then
                allocate(Orb_HF_MBPT%p%p(pN_particle), Orb_HF_MBPT%h%p(pN_hole))
                allocate(Orb_HF_MBPT%p%n(nN_particle), Orb_HF_MBPT%h%n(nN_hole))
            end if
        end do
        Orb_HF_MBPT%N_p = pninteger_t(pN_particle,nN_particle)
        Orb_HF_MBPT%N_h = pninteger_t(pN_hole,nN_hole)
    end function HF_MBPT_orbitals_make

    ! dev HF_MBPT_h1b(a,t_a,h_N) ... single-particle energy of orbital a with isospin t_a ...
    pure real(dp) function HF_MBPT_h1b(a,t_a,h_N)
        integer, intent(in) :: a, t_a
        type(o1b_t), intent(in) :: h_N
        if (t_a == -1) then
            HF_MBPT_h1b = h_N%p(a,a)
        else
            HF_MBPT_h1b = h_N%n(a,a)
        end if
    end function HF_MBPT_h1b

    ! Port: flat (a,l,j,E) particle & hole lists per isospin (the Orb1B lists with the h_N energies attached) ...
    subroutine HF_MBPT_states(Orb_HF_MBPT,h_N,Particle,Hole)
        type(orb1b_hf_mbpt_t), intent(in) :: Orb_HF_MBPT
        type(o1b_t), intent(in) :: h_N
        type(pnsvector_t), intent(out) :: Particle, Hole
        call fill(Orb_HF_MBPT%p%p,-1,Particle%p)
        call fill(Orb_HF_MBPT%p%n,+1,Particle%n)
        call fill(Orb_HF_MBPT%h%p,-1,Hole%p)
        call fill(Orb_HF_MBPT%h%n,+1,Hole%n)
    contains
        subroutine fill(O,t,S)
            type(orb1b_t), intent(in) :: O(:)
            integer, intent(in) :: t
            type(spstate_t), allocatable, intent(out) :: S(:)
            integer :: k
            allocate(S(size(O)))
            do k = 1, size(O)
                S(k) = spstate_t(O(k)%a,O(k)%l,O(k)%j,HF_MBPT_h1b(O(k)%a,t,h_N))
            end do
        end subroutine fill
    end subroutine HF_MBPT_states

    ! Combined isospin-labelled lists [(Orb,-1) for protons; (Orb,+1) for neutrons] of the ph sector & TBCF ...
    subroutine HF_MBPT_isospin_states(Orb_HF_MBPT,h_N,Orb_particle,Orb_hole)
        type(orb1b_hf_mbpt_t), intent(in) :: Orb_HF_MBPT
        type(o1b_t), intent(in) :: h_N
        type(ispstate_t), allocatable, intent(out) :: Orb_particle(:), Orb_hole(:)
        allocate(Orb_particle(size(Orb_HF_MBPT%p%p) + size(Orb_HF_MBPT%p%n)))
        allocate(Orb_hole(size(Orb_HF_MBPT%h%p) + size(Orb_HF_MBPT%h%n)))
        call fill(Orb_HF_MBPT%p%p,Orb_HF_MBPT%p%n,Orb_particle)
        call fill(Orb_HF_MBPT%h%p,Orb_HF_MBPT%h%n,Orb_hole)
    contains
        subroutine fill(Op,On,S)
            type(orb1b_t), intent(in) :: Op(:), On(:)
            type(ispstate_t), intent(inout) :: S(:)
            integer :: k
            do k = 1, size(Op)
                S(k) = ispstate_t(Op(k)%a,Op(k)%l,Op(k)%j,-1,HF_MBPT_h1b(Op(k)%a,-1,h_N))
            end do
            do k = 1, size(On)
                S(size(Op) + k) = ispstate_t(On(k)%a,On(k)%l,On(k)%j,+1,HF_MBPT_h1b(On(k)%a,+1,h_N))
            end do
        end subroutine fill
    end subroutine HF_MBPT_isospin_states

    ! Antisymmetrized NN matrix element <ab|V|cd>_J of the channel given by t_a + t_b (-2 pp, 0 pn, +2 nn) ...
    pure real(dp) function HF_MBPT_V2b(a,b,c,d,J,P,T_ab,V_NN,Orb,Orb_NN) result(V)
        integer, intent(in) :: a, b, c, d, J, P, T_ab
        type(o2b_t), intent(in) :: V_NN
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        select case (T_ab)
        case (-2)
            V = O2b_pp(a,b,c,d,J,P,V_NN,Orb,Orb_NN)
        case (0)
            V = O2b_pn(a,b,c,d,J,P,V_NN,Orb_NN)
        case (2)
            V = O2b_nn(a,b,c,d,J,P,V_NN,Orb,Orb_NN)
        case default
            V = 0.0_dp
        end select
    end function HF_MBPT_V2b

    ! Cross-coupled (Pandya) matrix element ... dev HF_MBPT3_V2b_Pandya (87-152):
        ! V = - sum_I (2I+1) {j_a j_b I; j_c j_d J} X_I with X_I the pp/nn/pn element <ab|V|cd>_I (pn with the exchange
        ! phases of the (+1,-1) orderings) ... for |t_a + t_b| = 2 t_c, t_d are not checked (callers conserve charge) ...
    pure real(dp) function HF_MBPT3_V2b_Pandya(sa,sb,sc,sd,J,Orb,Orb_NN,V_NN) result(V)
        type(ispstate_t), intent(in) :: sa, sb, sc, sd
        integer, intent(in) :: J
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: V_NN
        integer :: P_ab, I_min, I_max, I
        real(dp) :: ME

        V = 0.0_dp
        P_ab = mod(sa%l + sb%l,2) + 1
        if ((mod(sc%l + sd%l,2) + 1) /= P_ab) return

        I_min = max(abs(sa%j - sb%j), abs(sc%j - sd%j))/2
        I_max = min(sa%j + sb%j, sc%j + sd%j)/2

        if ((sa%t + sb%t) == -2) then
            do I = I_min, I_max
                ME = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*I,sc%j,sd%j,2*J) * O2b_pp(sa%a,sb%a,sc%a,sd%a,I,P_ab,V_NN,Orb,Orb_NN)
                V = V - ME
            end do
            return
        end if

        if ((sa%t + sb%t) == 2) then
            do I = I_min, I_max
                ME = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*I,sc%j,sd%j,2*J) * O2b_nn(sa%a,sb%a,sc%a,sd%a,I,P_ab,V_NN,Orb,Orb_NN)
                V = V - ME
            end do
            return
        end if

        if ((sa%t + sb%t) == 0) then
            if (sa%t == -1 .and. sb%t == 1 .and. sc%t == -1 .and. sd%t == 1) then
                do I = I_min, I_max
                    ME = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*I,sc%j,sd%j,2*J) * O2b_pn(sa%a,sb%a,sc%a,sd%a,I,P_ab,V_NN,Orb_NN)
                    V = V - ME
                end do
                return
            end if
            if (sa%t == -1 .and. sb%t == 1 .and. sc%t == 1 .and. sd%t == -1) then
                do I = I_min, I_max
                    ME = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*I,sc%j,sd%j,2*J) * O2b_pn(sa%a,sb%a,sd%a,sc%a,I,P_ab,V_NN,Orb_NN) * &
                         phase((sc%j + sd%j)/2 - I)
                    V = V - ME
                end do
                return
            end if
            if (sa%t == 1 .and. sb%t == -1 .and. sc%t == -1 .and. sd%t == 1) then
                do I = I_min, I_max
                    ME = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*I,sc%j,sd%j,2*J) * O2b_pn(sb%a,sa%a,sc%a,sd%a,I,P_ab,V_NN,Orb_NN) * &
                         phase((sa%j + sb%j)/2 - I)
                    V = V - ME
                end do
                return
            end if
            if (sa%t == 1 .and. sb%t == -1 .and. sc%t == 1 .and. sd%t == -1) then
                do I = I_min, I_max
                    ME = real(2*I + 1,dp) * f6j(sa%j,sb%j,2*I,sc%j,sd%j,2*J) * O2b_pn(sb%a,sa%a,sd%a,sc%a,I,P_ab,V_NN,Orb_NN) * &
                         phase((sa%j + sb%j + sc%j + sd%j)/2)
                    V = V - ME
                end do
                return
            end if
        end if
    end function HF_MBPT3_V2b_Pandya

end module nj_hf_mbpt_common
