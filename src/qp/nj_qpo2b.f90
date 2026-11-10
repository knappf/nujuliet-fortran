! 2-body operators in the quasiparticle representation ... counterpart of MatrixElements/qpO2B.jl, qpO2B_40.jl & qpO2B_22.jl ...
    ! Only the components used by the solvers are computed (O^40 & O^22 without the 1111 part); O^31 & pn1111 stay zero
    ! exactly as in the original (their allocation calls are commented out there) ...
module nj_qpo2b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_orbitals, only : orbitals_lj_make
    use nj_o2b
    implicit none
    private

    public :: qpO2b_initialize, qpO2b, qpO2b_export, qpO2b_import
    public :: qpO2b_40_pp, qpO2b_40_pn, qpO2b_40_nn, qpO2b_22_pp, qpO2b_2002_pn, qpO2b_22_nn
    public :: qp_BMBPT_energy

    ! Specification of the 4-index quasiparticle transformations ... source term, U/V selectors & signs ...
    integer, parameter :: Src_V = 1, Src_Fcd = 2, Src_Fdc = 3
    integer, parameter :: Sel_U = 1, Sel_V = 2

contains

    ! Zero quasiparticle operator ... layouts as qpO2b_initialize in qpO2B.jl ...
    subroutine qpO2b_initialize(Params,Orb,Q_NN,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(qpo2b_t), intent(out) :: Q_NN
        type(orb2b_t), intent(out) :: Orb_NN
        type(o2b_t) :: Dummy
        integer :: J_max, P, J, N0, N1

        call O2b_initialize(Params,Orb,Dummy,Orb_NN)
        J_max = Params%Calc%N2max + 1
        allocate(Q_NN%qp40%pp(2,0:J_max), Q_NN%qp40%pn(2,0:J_max), Q_NN%qp40%nn(2,0:J_max))
        allocate(Q_NN%qp31%pp(2,0:J_max), Q_NN%qp31%pn2011(2,0:J_max), Q_NN%qp31%pn1120(2,0:J_max), Q_NN%qp31%nn(2,0:J_max))
        allocate(Q_NN%qp22%pp(2,0:J_max), Q_NN%qp22%pn2002(2,0:J_max), Q_NN%qp22%pn1111(2,0:J_max), Q_NN%qp22%nn(2,0:J_max))
        do P = 1, 2
            do J = 0, J_max
                N0 = Orb_NN%N(1,P,J)
                N1 = Orb_NN%N(2,P,J)
                allocate(Q_NN%qp40%pp(P,J)%v(N1*(N1+1)/2), Q_NN%qp40%pn(P,J)%v(N1,N1), Q_NN%qp40%nn(P,J)%v(N1*(N1+1)/2))
                allocate(Q_NN%qp31%pp(P,J)%v(N0,N0), Q_NN%qp31%pn2011(P,J)%v(N0,N0), Q_NN%qp31%pn1120(P,J)%v(N0,N0), &
                         Q_NN%qp31%nn(P,J)%v(N0,N0))
                allocate(Q_NN%qp22%pp(P,J)%v(N1*(N1+1)/2), Q_NN%qp22%pn2002(P,J)%v(N1,N1), &
                         Q_NN%qp22%pn1111(P,J)%v(N0*(N0+1)/2), Q_NN%qp22%nn(P,J)%v(N1*(N1+1)/2))
                Q_NN%qp40%pp(P,J)%v = 0.0_dp
                Q_NN%qp40%pn(P,J)%v = 0.0_dp
                Q_NN%qp40%nn(P,J)%v = 0.0_dp
                Q_NN%qp31%pp(P,J)%v = 0.0_dp
                Q_NN%qp31%pn2011(P,J)%v = 0.0_dp
                Q_NN%qp31%pn1120(P,J)%v = 0.0_dp
                Q_NN%qp31%nn(P,J)%v = 0.0_dp
                Q_NN%qp22%pp(P,J)%v = 0.0_dp
                Q_NN%qp22%pn2002(P,J)%v = 0.0_dp
                Q_NN%qp22%pn1111(P,J)%v = 0.0_dp
                Q_NN%qp22%nn(P,J)%v = 0.0_dp
            end do
        end do
    end subroutine qpO2b_initialize

    ! Quasiparticle representation of the 2-body operator O_NN ... O_NN is transformed to the basis C in place ...
    subroutine qpO2b(Params,Orb,Orb_NN,O_NN,C,U,V,Q_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(inout) :: O_NN
        type(o1b_t), intent(in) :: C, U, V
        type(qpo2b_t), intent(out) :: Q_NN
        type(orb2b_t) :: Orb_NN_q
        type(o2b_temp_t) :: F_NN
        type(orb2b_temp_t) :: Orb_NN_t
        type(ivec_t), allocatable :: Orb_lj(:,:)
        real(dp), allocatable :: UV(:,:,:,:)
        integer :: a_max, J_t
        integer(i8) :: t0

        write(*,'(/,a)') 'Preparing the given 2-body NN operator O expressed in the J-scheme quasiparticle representation ...'
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2

        call O2b_transformation(Params,Orb,Orb_NN,O_NN,C)
        call qpO2b_initialize(Params,Orb,Q_NN,Orb_NN_q)

        call timer_start(t0)
        call O2b_Pandya_transformation(Params,Orb,Orb_NN,O_NN,F_NN)
        call timer_report(t0)

        ! UV(:,:,Sel,iso) ... U & V amplitudes for protons (1) & neutrons (2) ...
        allocate(UV(a_max,a_max,2,2))
        UV(:,:,Sel_U,1) = U%p
        UV(:,:,Sel_V,1) = V%p
        UV(:,:,Sel_U,2) = U%n
        UV(:,:,Sel_V,2) = V%n

        call qp_temp_orbitals(Params,Orb,Orb_NN_t,J_t)
        call orbitals_lj_make(Params,Orb,Orb_lj)

        call timer_start(t0)
        write(*,'(/,a)') TAB // 'Allocating components O^(40) of the given NN operator 2-body in the quasiparticle representation ...'
        write(*,'(/,a)') TAB // TAB // 'Allocating the components of the 2-body proton-proton operator O^(40) ...'
        call engine(1, 6, [Src_V,Src_V,Src_Fdc,Src_Fdc,Src_Fcd,Src_Fcd], &
                    reshape([Sel_U,Sel_U,Sel_V,Sel_V, Sel_V,Sel_V,Sel_U,Sel_U, Sel_U,Sel_V,Sel_U,Sel_V, &
                             Sel_V,Sel_U,Sel_V,Sel_U, Sel_U,Sel_V,Sel_V,Sel_U, Sel_V,Sel_U,Sel_U,Sel_V],[4,6]), &
                    [1,1,1,1], [-1.0_dp,-1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp], 1)
        write(*,'(/,a)') TAB // TAB // 'Allocating the components of the 2-body proton-neutron operator O^(40) ...'
        call engine(2, 4, [Src_Fdc,Src_Fdc,Src_Fcd,Src_Fcd,0,0], &
                    reshape([Sel_U,Sel_V,Sel_U,Sel_V, Sel_V,Sel_U,Sel_V,Sel_U, Sel_U,Sel_V,Sel_V,Sel_U, &
                             Sel_V,Sel_U,Sel_U,Sel_V, 0,0,0,0, 0,0,0,0],[4,6]), &
                    [1,1,2,2], [1.0_dp,1.0_dp,-1.0_dp,-1.0_dp,0.0_dp,0.0_dp], 2)
        write(*,'(/,a)') TAB // TAB // 'Allocating the components of the 2-body neutron-neutron operator O^(40) ...'
        call engine(3, 6, [Src_V,Src_V,Src_Fdc,Src_Fdc,Src_Fcd,Src_Fcd], &
                    reshape([Sel_U,Sel_U,Sel_V,Sel_V, Sel_V,Sel_V,Sel_U,Sel_U, Sel_U,Sel_V,Sel_U,Sel_V, &
                             Sel_V,Sel_U,Sel_V,Sel_U, Sel_U,Sel_V,Sel_V,Sel_U, Sel_V,Sel_U,Sel_U,Sel_V],[4,6]), &
                    [2,2,2,2], [-1.0_dp,-1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp], 3)
        write(*,'(/,a)') TAB // 'Components of O^(40) of the given 2-body NN quasiparticle operator succesfully allocated ...'

        write(*,'(/,a)') TAB // 'Allocating components O^(22) of the given NN operator 2-body in the quasiparticle representation ...'
        write(*,'(/,a)') TAB // TAB // 'Allocating the components of the 2-body proton-proton operator O^(22) ...'
        call engine(1, 6, [Src_V,Src_V,Src_Fcd,Src_Fcd,Src_Fdc,Src_Fdc], &
                    reshape([Sel_U,Sel_U,Sel_U,Sel_U, Sel_V,Sel_V,Sel_V,Sel_V, Sel_U,Sel_V,Sel_U,Sel_V, &
                             Sel_V,Sel_U,Sel_V,Sel_U, Sel_U,Sel_V,Sel_V,Sel_U, Sel_V,Sel_U,Sel_U,Sel_V],[4,6]), &
                    [1,1,1,1], [1.0_dp,1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp], 4)
        write(*,'(/,a)') TAB // TAB // 'Allocating the components of the 2-body proton-neutron operator O^(22) ...'
        call engine(2, 4, [Src_Fcd,Src_Fcd,Src_Fdc,Src_Fdc,0,0], &
                    reshape([Sel_U,Sel_V,Sel_U,Sel_V, Sel_V,Sel_U,Sel_V,Sel_U, Sel_U,Sel_V,Sel_V,Sel_U, &
                             Sel_V,Sel_U,Sel_U,Sel_V, 0,0,0,0, 0,0,0,0],[4,6]), &
                    [1,1,2,2], [1.0_dp,1.0_dp,-1.0_dp,-1.0_dp,0.0_dp,0.0_dp], 5)
        write(*,'(/,a)') TAB // TAB // 'Allocating the components of the 2-body neutron-neutron operator O^(22) ...'
        call engine(3, 6, [Src_V,Src_V,Src_Fcd,Src_Fcd,Src_Fdc,Src_Fdc], &
                    reshape([Sel_U,Sel_U,Sel_U,Sel_U, Sel_V,Sel_V,Sel_V,Sel_V, Sel_U,Sel_V,Sel_U,Sel_V, &
                             Sel_V,Sel_U,Sel_V,Sel_U, Sel_U,Sel_V,Sel_V,Sel_U, Sel_V,Sel_U,Sel_U,Sel_V],[4,6]), &
                    [2,2,2,2], [1.0_dp,1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp], 6)
        write(*,'(/,a)') TAB // 'Components of O^(22) of the given 2-body NN quasiparticle operator succesfully allocated ...'
        call timer_report(t0)

        write(*,'(/,a,/)') 'The given 2-body NN operator O expressed in the quasiparticke represation with components (kl) ...'

    contains

        ! One quasiparticle transformation ... Fch = 1 pp, 2 pn, 3 nn (F_NN channel & direct V) ...
            ! Sel(s,k) is the U/V selector of index s for component k, Iso(s) the isospin of index s ...
            ! Out = 1 qp40%pp, 2 qp40%pn, 3 qp40%nn, 4 qp22%pp, 5 qp22%pn2002, 6 qp22%nn ...
        subroutine engine(Fch,K,Src,Sel,Iso,Sgn,Out)
            integer, intent(in) :: Fch, K, Src(6), Sel(4,6), Iso(4), Out
            real(dp), intent(in) :: Sgn(6)
            type(rmat_t), allocatable :: T1(:,:,:), T2(:,:,:)
            integer :: J, P, N, Bra, Ket, BraKet, a, b, c, d, i, ii, kk, Ket_dc, Bra_ib, Bra_ai, Ket_id, Bra_ab, Ket_ci
            integer :: N1, Ind, Ket_max
            real(dp) :: Sums(6), Ph, V_ibcd, F_ibcd, F_ibdc, Src_val(3), Tot

            allocate(T1(6,2,0:J_t), T2(6,2,0:J_t))
            do J = 0, J_t
                do P = 1, 2
                    N = Orb_NN_t%N(P,J)
                    do kk = 1, K
                        allocate(T1(kk,P,J)%v(N,N), T2(kk,P,J)%v(N,N))
                    end do
                end do
            end do

            ! 1st index ...
            do J = 0, J_t
                do P = 1, 2
                    N = Orb_NN_t%N(P,J)
                    !$omp parallel do schedule(dynamic) private(Bra,Ket,a,b,c,d,Ph,Ket_dc,Sums,ii,i,Bra_ib,V_ibcd,F_ibdc,F_ibcd,Src_val,kk)
                    do BraKet = 1, N*N
                        Bra = (BraKet - 1)/N + 1
                        Ket = mod(BraKet - 1,N) + 1
                        a = Orb_NN_t%Ind(P,J)%v(1,Bra)
                        b = Orb_NN_t%Ind(P,J)%v(2,Bra)
                        c = Orb_NN_t%Ind(P,J)%v(1,Ket)
                        d = Orb_NN_t%Ind(P,J)%v(2,Ket)
                        Ph = phase(J + (Orb(c)%j + Orb(d)%j)/2)
                        Ket_dc = Orb_NN_t%Dic(d,c,J)
                        Sums = 0.0_dp
                        do ii = 1, size(Orb_lj(Orb(a)%l,Orb(a)%j)%v)
                            i = Orb_lj(Orb(a)%l,Orb(a)%j)%v(ii)
                            Bra_ib = Orb_NN_t%Dic(i,b,J)
                            V_ibcd = 0.0_dp
                            if (Fch == 1) then
                                if (any(Src(1:K) == Src_V)) V_ibcd = O2b_pp(i,b,c,d,J,P,O_NN,Orb,Orb_NN)
                                F_ibdc = F_NN%pp(P,J)%v(Bra_ib,Ket_dc)
                                F_ibcd = F_NN%pp(P,J)%v(Bra_ib,Ket)
                            else if (Fch == 2) then
                                F_ibdc = F_NN%pn(P,J)%v(Bra_ib,Ket_dc)
                                F_ibcd = F_NN%pn(P,J)%v(Bra_ib,Ket)
                            else
                                if (any(Src(1:K) == Src_V)) V_ibcd = O2b_nn(i,b,c,d,J,P,O_NN,Orb,Orb_NN)
                                F_ibdc = F_NN%nn(P,J)%v(Bra_ib,Ket_dc)
                                F_ibcd = F_NN%nn(P,J)%v(Bra_ib,Ket)
                            end if
                            Src_val = [V_ibcd, F_ibcd, Ph * F_ibdc]
                            do kk = 1, K
                                Sums(kk) = Sums(kk) + UV(i,a,Sel(1,kk),Iso(1)) * Src_val(Src(kk))
                            end do
                        end do
                        do kk = 1, K
                            T1(kk,P,J)%v(Bra,Ket) = Sums(kk)
                        end do
                    end do
                    !$omp end parallel do
                end do
            end do

            ! 2nd index ...
            do J = 0, J_t
                do P = 1, 2
                    N = Orb_NN_t%N(P,J)
                    !$omp parallel do schedule(dynamic) private(a,b,Ket,Sums,ii,i,Bra_ai,kk)
                    do Bra = 1, N
                        a = Orb_NN_t%Ind(P,J)%v(1,Bra)
                        b = Orb_NN_t%Ind(P,J)%v(2,Bra)
                        do Ket = 1, N
                            Sums = 0.0_dp
                            do ii = 1, size(Orb_lj(Orb(b)%l,Orb(b)%j)%v)
                                i = Orb_lj(Orb(b)%l,Orb(b)%j)%v(ii)
                                Bra_ai = Orb_NN_t%Dic(a,i,J)
                                do kk = 1, K
                                    Sums(kk) = Sums(kk) + UV(i,b,Sel(2,kk),Iso(2)) * T1(kk,P,J)%v(Bra_ai,Ket)
                                end do
                            end do
                            do kk = 1, K
                                T2(kk,P,J)%v(Bra,Ket) = Sums(kk)
                            end do
                        end do
                    end do
                    !$omp end parallel do
                end do
            end do

            ! 3rd index ... back into T1 ...
            do J = 0, J_t
                do P = 1, 2
                    N = Orb_NN_t%N(P,J)
                    !$omp parallel do schedule(dynamic) private(Ket,c,d,Sums,ii,i,Ket_id,kk)
                    do Bra = 1, N
                        do Ket = 1, N
                            c = Orb_NN_t%Ind(P,J)%v(1,Ket)
                            d = Orb_NN_t%Ind(P,J)%v(2,Ket)
                            Sums = 0.0_dp
                            do ii = 1, size(Orb_lj(Orb(c)%l,Orb(c)%j)%v)
                                i = Orb_lj(Orb(c)%l,Orb(c)%j)%v(ii)
                                Ket_id = Orb_NN_t%Dic(i,d,J)
                                do kk = 1, K
                                    Sums(kk) = Sums(kk) + UV(i,c,Sel(3,kk),Iso(3)) * T2(kk,P,J)%v(Bra,Ket_id)
                                end do
                            end do
                            do kk = 1, K
                                T1(kk,P,J)%v(Bra,Ket) = Sums(kk)
                            end do
                        end do
                    end do
                    !$omp end parallel do
                end do
            end do

            ! 4th index & allocation of the quasiparticle components ...
            do J = 0, J_t
                do P = 1, 2
                    N1 = Orb_NN%N(2,P,J)
                    !$omp parallel do schedule(dynamic) private(a,b,Bra_ab,Ket_max,Ket,c,d,Tot,ii,i,Ket_ci,kk,Ind)
                    do Bra = 1, N1
                        a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                        b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                        Bra_ab = Orb_NN_t%Dic(a,b,J)
                        Ket_max = Bra
                        if (Out == 2 .or. Out == 5) Ket_max = N1
                        do Ket = 1, Ket_max
                            c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                            d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                            Tot = 0.0_dp
                            do ii = 1, size(Orb_lj(Orb(d)%l,Orb(d)%j)%v)
                                i = Orb_lj(Orb(d)%l,Orb(d)%j)%v(ii)
                                Ket_ci = Orb_NN_t%Dic(c,i,J)
                                do kk = 1, K
                                    Tot = Tot + Sgn(kk) * (UV(i,d,Sel(4,kk),Iso(4)) * T1(kk,P,J)%v(Bra_ab,Ket_ci))
                                end do
                            end do
                            Ind = packed_index(Bra,Ket,N1)
                            select case (Out)
                            case (1)
                                Q_NN%qp40%pp(P,J)%v(Ind) = Tot
                            case (2)
                                Q_NN%qp40%pn(P,J)%v(Bra,Ket) = Tot
                            case (3)
                                Q_NN%qp40%nn(P,J)%v(Ind) = Tot
                            case (4)
                                Q_NN%qp22%pp(P,J)%v(Ind) = Tot
                            case (5)
                                Q_NN%qp22%pn2002(P,J)%v(Bra,Ket) = Tot
                            case (6)
                                Q_NN%qp22%nn(P,J)%v(Ind) = Tot
                            end select
                        end do
                    end do
                    !$omp end parallel do
                end do
            end do
        end subroutine engine

    end subroutine qpO2b

    ! Temporary orbitals of all ordered pairs (a,b) with N_2max = 2*N_max ... qpO2b_40_initialize / qpO2b_22_initialize ...
        ! J_t is the highest channel iterated over (J loops of the original run to N2max + 1) ...
    subroutine qp_temp_orbitals(Params,Orb,Orb_NN_t,J_t)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_temp_t), intent(out) :: Orb_NN_t
        integer, intent(out) :: J_t
        integer :: N_max, a_max, J_max, a, b, P, J, Pass
        integer, allocatable :: N_Orb(:,:)

        N_max = Params%Calc%Nmax
        a_max = (N_max + 1)*(N_max + 2)/2
        J_max = 2*N_max + 1
        J_t = min(Params%Calc%N2max + 1, J_max)
        allocate(N_Orb(2,0:J_max), Orb_NN_t%Ind(2,0:J_max), Orb_NN_t%Dic(a_max,a_max,0:J_max))
        Orb_NN_t%Dic = 0
        do Pass = 1, 2
            N_Orb = 0
            do a = 1, a_max
                do b = 1, a_max
                    if ((2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l) <= 2*N_max) then
                        P = mod(Orb(a)%l + Orb(b)%l,2) + 1
                        do J = abs(Orb(a)%j - Orb(b)%j)/2, (Orb(a)%j + Orb(b)%j)/2
                            N_Orb(P,J) = N_Orb(P,J) + 1
                            if (Pass == 2) then
                                Orb_NN_t%Dic(a,b,J) = N_Orb(P,J)
                                Orb_NN_t%Ind(P,J)%v(1,N_Orb(P,J)) = a
                                Orb_NN_t%Ind(P,J)%v(2,N_Orb(P,J)) = b
                            end if
                        end do
                    end if
                end do
            end do
            if (Pass == 1) then
                do P = 1, 2
                    do J = 0, J_max
                        allocate(Orb_NN_t%Ind(P,J)%v(2,N_Orb(P,J)))
                    end do
                end do
            end if
        end do
        Orb_NN_t%N = N_Orb
    end subroutine qp_temp_orbitals

    ! Exchange ordering of T = 1 pair indices with phases ... same convention as O2b_pp ...
    pure subroutine t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        integer, intent(in) :: a, b, c, d, J
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(out) :: a_s, b_s, c_s, d_s
        real(dp), intent(out) :: Amp
        real(dp) :: Sign_ab, Sign_cd
        a_s = a
        b_s = b
        c_s = c
        d_s = d
        Sign_ab = 1.0_dp
        Sign_cd = 1.0_dp
        if (a < b) then
            a_s = b
            b_s = a
            Sign_ab = phase(J + 1 + (Orb(a)%j + Orb(b)%j)/2)
        end if
        if (c < d) then
            c_s = d
            d_s = c
            Sign_cd = phase(J + 1 + (Orb(c)%j + Orb(d)%j)/2)
        end if
        Amp = Sign_ab * Sign_cd
        if ((a < b) .and. (c < d)) Amp = phase((Orb(a)%j + Orb(b)%j + Orb(c)%j + Orb(d)%j)/2)
    end subroutine t1_order

    pure real(dp) function qpO2b_40_pp(a,b,c,d,J,P,Q,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(qpo2b_t), intent(in) :: Q
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        qpO2b_40_pp = Amp * Q%qp40%pp(P,J)%v(O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN))
    end function qpO2b_40_pp

    pure real(dp) function qpO2b_40_nn(a,b,c,d,J,P,Q,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(qpo2b_t), intent(in) :: Q
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        qpO2b_40_nn = Amp * Q%qp40%nn(P,J)%v(O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN))
    end function qpO2b_40_nn

    pure real(dp) function qpO2b_40_pn(a,b,c,d,J,P,Q,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(qpo2b_t), intent(in) :: Q
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        qpO2b_40_pn = Amp * Q%qp40%pn(P,J)%v(Orb_NN%Dic(a_s,b_s,J,1),Orb_NN%Dic(c_s,d_s,J,1))
    end function qpO2b_40_pn

    pure real(dp) function qpO2b_22_pp(a,b,c,d,J,P,Q,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(qpo2b_t), intent(in) :: Q
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        qpO2b_22_pp = Amp * Q%qp22%pp(P,J)%v(O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN))
    end function qpO2b_22_pp

    pure real(dp) function qpO2b_22_nn(a,b,c,d,J,P,Q,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(qpo2b_t), intent(in) :: Q
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        qpO2b_22_nn = Amp * Q%qp22%nn(P,J)%v(O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN))
    end function qpO2b_22_nn

    pure real(dp) function qpO2b_2002_pn(a,b,c,d,J,P,Q,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(qpo2b_t), intent(in) :: Q
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        qpO2b_2002_pn = Amp * Q%qp22%pn2002(P,J)%v(Orb_NN%Dic(a_s,b_s,J,1),Orb_NN%Dic(c_s,d_s,J,1))
    end function qpO2b_2002_pn

    ! Binary layout: 11 blocks (H = 1..11), each over J = 0..2*Nmax+1 & P = 1,2 ... packed lower triangles or full matrices ...
    subroutine qpO2b_export(Params,Orb_NN,Q,Export_Path)
        type(params_t), intent(in) :: Params
        type(orb2b_t), intent(in) :: Orb_NN
        type(qpo2b_t), intent(in) :: Q
        character(len=*), intent(in) :: Export_Path
        type(qpo2b_t) :: Q_copy
        integer :: u

        write(*,'(/,a)') 'Exporting the given quasiparticle 2-body operator in the binary format ...'
        Q_copy = Q
        open(newunit=u, file=Export_Path, access='stream', form='unformatted', status='replace', action='write')
        call stream_qp(u,Params,Orb_NN,Q_copy,.true.)
        close(u)
        write(*,'(/,a)') 'Given 2-body quasiparticle operator succesfully exported in the binary format ...'
    end subroutine qpO2b_export

    subroutine qpO2b_import(Params,Orb,Import_Path,Q,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: Import_Path
        type(qpo2b_t), intent(out) :: Q
        type(orb2b_t), intent(out) :: Orb_NN
        integer :: u

        call qpO2b_initialize(Params,Orb,Q,Orb_NN)
        write(*,'(/,a)') 'Reading given 2-body quasiparticle NN operator ...'
        open(newunit=u, file=Import_Path, access='stream', form='unformatted', status='old', action='read')
        call stream_qp(u,Params,Orb_NN,Q,.false.)
        close(u)
        write(*,'(/,a)') 'Given 2-body quasiparticle NN operator succesfully imported ...'
    end subroutine qpO2b_import

    subroutine stream_qp(u,Params,Orb_NN,Q,Write_Mode)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        type(orb2b_t), intent(in) :: Orb_NN
        type(qpo2b_t), intent(inout) :: Q
        logical, intent(in) :: Write_Mode
        integer :: H, J, P, J_e

        J_e = min(2*Params%Calc%Nmax + 1, Orb_NN%J_max)
        do H = 1, 11
            do J = 0, J_e
                do P = 1, 2
                    select case (H)
                    case (1)
                        call packed(Q%qp40%pp(P,J)%v, Orb_NN%N(2,P,J))
                    case (2)
                        call full(Q%qp40%pn(P,J)%v)
                    case (3)
                        call packed(Q%qp40%nn(P,J)%v, Orb_NN%N(2,P,J))
                    case (4)
                        call full(Q%qp31%pp(P,J)%v)
                    case (5)
                        call full(Q%qp31%pn2011(P,J)%v)
                    case (6)
                        call full(Q%qp31%pn1120(P,J)%v)
                    case (7)
                        call full(Q%qp31%nn(P,J)%v)
                    case (8)
                        call packed(Q%qp22%pp(P,J)%v, Orb_NN%N(2,P,J))
                    case (9)
                        call full(Q%qp22%pn2002(P,J)%v)
                    case (10)
                        call packed(Q%qp22%pn1111(P,J)%v, Orb_NN%N(1,P,J))
                    case (11)
                        call packed(Q%qp22%nn(P,J)%v, Orb_NN%N(2,P,J))
                    end select
                end do
            end do
        end do
    contains
        ! Bra-major order of the lower triangle ...
        subroutine packed(X,N)
            real(dp), intent(inout) :: X(:)
            integer, intent(in) :: N
            real(dp), allocatable :: Buf(:)
            integer :: Bra, Ket
            allocate(Buf(N*(N+1)/2))
            if (Write_Mode) then
                do Bra = 1, N
                    do Ket = 1, Bra
                        Buf(Bra*(Bra-1)/2 + Ket) = X(packed_index(Bra,Ket,N))
                    end do
                end do
                write(u) Buf
            else
                read(u) Buf
                do Bra = 1, N
                    do Ket = 1, Bra
                        X(packed_index(Bra,Ket,N)) = Buf(Bra*(Bra-1)/2 + Ket)
                    end do
                end do
            end if
        end subroutine packed
        ! Bra-major (row-major) order of a full matrix ...
        subroutine full(X)
            real(dp), intent(inout) :: X(:,:)
            if (Write_Mode) then
                write(u) transpose(X)
            else
                read(u) X
                X = transpose(X)
            end if
        end subroutine full
    end subroutine stream_qp

    ! Leading-order BMBPT(2) correlation energy ... HF_BCS_BMBPT_energy & HFB_BMBPT_energy ...
    subroutine qp_BMBPT_energy(Params,Orb,Orb_NN,H_N,H_NN,Label,Summary_File)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(qpo1b_t), intent(in) :: H_N
        type(qpo2b_t), intent(in) :: H_NN
        character(len=*), intent(in) :: Label, Summary_File
        integer :: a_max, ab, a, b, T_ab, c, d, J, Par, u
        real(dp) :: dE, E_a, E_b, E_c, E_d, Denom

        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        write(*,'(/,a)') 'Calculating ' // Label // '-BMBPT(2) ground-state energy correction ... '

        dE = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:dE) private(T_ab,a,b,Par,E_a,E_b,c,d,E_c,E_d,Denom,J)
        do ab = 1, 3*a_max*a_max
            T_ab = (ab - 1)/(a_max*a_max) - 1
            a = mod((ab - 1)/a_max, a_max) + 1
            b = mod(ab - 1, a_max) + 1
            Par = mod(Orb(a)%l + Orb(b)%l,2) + 1
            if (T_ab == 1) then
                E_a = H_N%qp11%n(a,a)
                E_b = H_N%qp11%n(b,b)
            else
                E_a = H_N%qp11%p(a,a)
                E_b = H_N%qp11%p(b,b)
            end if
            do c = 1, a_max
                do d = 1, a_max
                    if ((mod(Orb(c)%l + Orb(d)%l,2) + 1) /= Par) cycle
                    if (T_ab == -1) then
                        E_c = H_N%qp11%p(c,c)
                        E_d = H_N%qp11%p(d,d)
                        Denom = -1.0_dp / (E_a + E_b + E_c + E_d) / 24.0_dp
                    else if (T_ab == 0) then
                        E_c = H_N%qp11%n(c,c)
                        E_d = H_N%qp11%n(d,d)
                        Denom = -1.0_dp / (E_a + E_b + E_c + E_d) / 4.0_dp
                    else
                        E_c = H_N%qp11%n(c,c)
                        E_d = H_N%qp11%n(d,d)
                        Denom = -1.0_dp / (E_a + E_b + E_c + E_d) / 24.0_dp
                    end if
                    do J = max(abs(Orb(a)%j - Orb(b)%j)/2, abs(Orb(c)%j - Orb(d)%j)/2), &
                           min((Orb(a)%j + Orb(b)%j)/2, (Orb(c)%j + Orb(d)%j)/2)
                        if (T_ab == -1) then
                            dE = dE + real(2*J + 1,dp) * Denom * abs(qpO2b_40_pp(a,b,c,d,J,Par,H_NN,Orb,Orb_NN))**2
                        else if (T_ab == 0) then
                            dE = dE + real(2*J + 1,dp) * Denom * abs(qpO2b_40_pn(a,b,c,d,J,Par,H_NN,Orb,Orb_NN))**2
                        else
                            dE = dE + real(2*J + 1,dp) * Denom * abs(qpO2b_40_nn(a,b,c,d,J,Par,H_NN,Orb,Orb_NN))**2
                        end if
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        write(*,'(a,/)') TAB // Label // '-BMBPT(2) energy    ...   E^(2) = ' // jstr(round_sig(dE,9)) // ' MeV'

        open(newunit=u, file=Summary_File, position='append', action='write')
        if (Label == 'HFB') then
            write(u,'(/,a)') 'Spherical Hartree-Fock-Bogoliubov Leading Order Bogoliubov Many-Body Perturbation Theory (2) solution review:'
        else
            write(u,'(/,a)') 'Spherical Hartree-Fock-BCS Leading Order Bogoliubov Many-Body Perturbation Theory (2) solution review:'
        end if
        write(u,'(/,a)',advance='no') 'E_0^(2)     = ' // cfmt_f(dE,15,8) // ' ' // TAB // ' MeV ' // TAB // TAB // ' ... ' // TAB // &
                                      ' LO ' // Label // '-BMBPT(2) correction to ground-state energy'
        write(u,'(/,a)') 'E_0^(2) / A = ' // cfmt_f(dE / real(Params%Calc%A,dp),15,8) // ' ' // TAB // ' MeV / Nucleon ' // TAB // &
                         TAB // ' ... ' // TAB // ' LO ' // Label // '-BMBPT(2) correction to ground-state energy'
        close(u)
    end subroutine qp_BMBPT_energy

end module nj_qpo2b
