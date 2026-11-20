! 2-body NN interaction ... counterpart of src/common/V2b.jl (NuJuliet dev; MatrixElements/V2B.jl before) ...
module nj_v2b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_o1b, only : T1b, T2b
    use nj_o2b
    use nj_v3b, only : V3b_no2b
    use nj_orbitals, only : orbitals_make
    implicit none
    private

    public :: V2b_read, V2b_count, V2b_residual_no2b, V2b_parse, V2b_correlation_energy

    ! Buffer size (# of Float64 matrix elements read/written at once) ... a multiple of 4 ...
    integer, parameter :: Buffer_size = 1000000

contains

    ! Reads the NuHamil NN file (Float64 stream of T=0, nn, T=1 pn, pp per J) & adds the 2-body CM correction ...
        ! The interaction-space orbitals (Int.Nmax) are made here, as in dev ... the file is read sequentially in blocks of
        ! at most Buffer_size elements, & the reading stops after the last chunk a with N_a <= Calc.Nmax (nothing later is
        ! stored; dev reads the whole file) ...
    subroutine V2b_read(Params,V_NN,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), allocatable :: Orb(:)
        type(o2b_t), intent(out) :: V_NN
        type(orb2b_t), intent(out) :: Orb_NN
        type(o1b_t) :: T
        integer :: N_max_Int, N_2max_Int, N_max_Calc, N_2max_Calc, A_mass, a_max, J_max
        integer :: u, a, b, c, d, d_max, J, P, Ind, k, Bra, Ket, N_t0, N_t1, N_T1_00
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_c, l_c, j_c, n_d, l_d, j_d
        integer :: NN_a, NN_b, NN_c, NN_d, NN_ab, NN_cd, NN_min, NN_max
        integer(i8), allocatable :: N_Chunk(:)
        integer(i8) :: N_left
        integer :: Buf_N
        real(dp), allocatable :: Buf(:)
        real(dp) :: hw, Amp, ME_00, ME_nn, ME_10, ME_pp, Ph, T_sym, T_antisym, ME, Amp_2, pG, nG
        character(len=32) :: CMS

        hw = Params%Int%hw
        N_max_Int = Params%Int%Nmax
        N_2max_Int = Params%Int%N2max
        N_max_Calc = Params%Calc%Nmax
        N_2max_Calc = Params%Calc%N2max
        A_mass = Params%Calc%A
        CMS = Params%Calc%CMS
        a_max = (N_max_Int + 1)*(N_max_Int + 2)/2
        J_max = N_2max_Calc + 1

        ! Orbitals of the interaction space ...
        Orb = orbitals_make(Params,Use_Int=.true.)

        ! Initialize 1-body kinetic operator ...
        T = T1b(Params,Orb)

        ! Initialize NN interaction array ...
        write(*,'(/,a)') 'Preparing 2-body NN interaction array...'
        call O2b_initialize(Params,Orb,V_NN,Orb_NN)

        ! Count the matrix elements of each chunk a ...
        N_Chunk = V2b_count(N_max_Int,N_2max_Int,Orb)

        write(*,'(/,a)') 'Reading 2-body NN interaction file...'
        open(newunit=u, file=trim(Params%Int%NN_File), access='stream', form='unformatted', status='old', action='read')
        allocate(Buf(Buffer_size))
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            NN_a = 2*n_a + l_a
            ! Orbitals are ordered by N ... no later chunk holds calculation-space elements ...
            if (NN_a > N_max_Calc) exit
            ! The chunk of the orbital a is read by blocks of at most Buffer_size matrix elements ...
            N_left = N_Chunk(a)
            Buf_N = 0
            k = 0
            do b = 1, a
                n_b = Orb(b)%n
                l_b = Orb(b)%l
                j_b = Orb(b)%j
                NN_b = 2*n_b + l_b
                NN_ab = NN_a + NN_b
                P = mod(l_a + l_b,2) + 1
                if (NN_ab > N_2max_Int) cycle
                do c = 1, a
                    n_c = Orb(c)%n
                    l_c = Orb(c)%l
                    j_c = Orb(c)%j
                    NN_c = 2*n_c + l_c
                    d_max = c
                    if (a == c) d_max = b
                    do d = 1, d_max
                        n_d = Orb(d)%n
                        l_d = Orb(d)%l
                        j_d = Orb(d)%j
                        NN_d = 2*n_d + l_d
                        NN_cd = NN_c + NN_d
                        if (.not. ((NN_cd <= N_2max_Int) .and. (mod(l_a + l_b,2) == mod(l_c + l_d,2)))) cycle
                        do J = max(abs(j_a - j_b), abs(j_c - j_d))/2, min(j_a + j_b, j_c + j_d)/2
                            ! Read the next block of the chunk, if the buffer was used up ...
                            if (k == Buf_N) then
                                Buf_N = int(min(N_left, int(Buffer_size,i8)))
                                read(u) Buf(1:Buf_N)
                                N_left = N_left - Buf_N
                                k = 0
                            end if
                            ME_00 = Buf(k+1)
                            ME_nn = Buf(k+2)
                            ME_10 = Buf(k+3)
                            ME_pp = Buf(k+4)
                            k = k + 4

                            if ((NN_a <= N_max_Calc) .and. (NN_b <= N_max_Calc) .and. (NN_c <= N_max_Calc) .and. &
                                (NN_d <= N_max_Calc) .and. (NN_ab <= N_2max_Calc) .and. (NN_cd <= N_2max_Calc)) then
                                Amp = 1.0_dp
                                if (a == b) Amp = Amp / sqrt(2.0_dp)
                                if (c == d) Amp = Amp / sqrt(2.0_dp)

                                if ((abs(Amp - 1.0_dp) < 1.0e-4_dp) .or. (mod(J,2) == 0)) then
                                    Ind = O2b_index(a,b,c,d,J,P,1,Orb_NN)
                                    V_NN%pp(P,J)%v(Ind) = V_NN%pp(P,J)%v(Ind) + ME_pp
                                    V_NN%nn(P,J)%v(Ind) = V_NN%nn(P,J)%v(Ind) + ME_nn
                                end if

                                ! T = 1 part of the pn interaction ...
                                if ((abs(Amp - 1.0_dp) < 1.0e-4_dp) .or. (mod(J,2) == 0)) then
                                    call add_pn(a,b,c,d,0.5_dp * ME_10)
                                    if (c /= d) then
                                        Ph = phase((j_c + j_d)/2 + J + 1)
                                        call add_pn(a,b,d,c,0.5_dp * Ph * ME_10)
                                    end if
                                    if ((a /= b) .and. (c /= d)) then
                                        Ph = phase((j_a + j_b + j_c + j_d)/2)
                                        call add_pn(b,a,d,c,0.5_dp * Ph * ME_10)
                                    end if
                                    if ((a /= b) .and. ((a /= c) .or. (b /= d))) then
                                        Ph = phase((j_a + j_b)/2 + J + 1)
                                        call add_pn(b,a,c,d,0.5_dp * Ph * ME_10)
                                    end if
                                end if

                                ! T = 0 part of the pn interaction ...
                                if ((abs(Amp - 1.0_dp) < 1.0e-4_dp) .or. (mod(J,2) == 1)) then
                                    call add_pn(a,b,c,d,0.5_dp * ME_00)
                                    if (c /= d) then
                                        Ph = phase((j_c + j_d)/2 + J)
                                        call add_pn(a,b,d,c,0.5_dp * Ph * ME_00)
                                    end if
                                    if ((a /= b) .and. (c /= d)) then
                                        Ph = phase((j_a + j_b + j_c + j_d)/2)
                                        call add_pn(b,a,d,c,0.5_dp * Ph * ME_00)
                                    end if
                                    if ((a /= b) .and. ((a /= c) .or. (b /= d))) then
                                        Ph = phase((j_a + j_b)/2 + J)
                                        call add_pn(b,a,c,d,0.5_dp * Ph * ME_00)
                                    end if
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
        deallocate(Buf)
        close(u)
        write(*,'(/,a)') 'NN 2-body interaction loaded...'

        ! 2-body CM correction ...
        if (CMS == 'CMS1+2B' .or. CMS == 'CMS2B') then
            write(*,'(/,a)') 'Loading 2-body CMS correction...'
            do J = 0, J_max
                do P = 1, 2
                    N_t0 = Orb_NN%N(1,P,J)
                    N_t1 = Orb_NN%N(2,P,J)
                    !$omp parallel do schedule(dynamic) private(Ket,Ind,a,b,c,d,j_c,j_d,T_sym,T_antisym,Amp,Amp_2,ME)
                    do Bra = 1, max(N_t0,N_t1)
                        do Ket = 1, Bra
                            if (Bra <= N_t0) then
                                Ind = packed_index(Bra,Ket,N_t0)
                                a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                                b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                                c = Orb_NN%Ind(1,P,J)%v(1,Ket)
                                d = Orb_NN%Ind(1,P,J)%v(2,Ket)
                                if (CMS == 'CMS1+2B') then
                                    T_sym = T2b(Orb,a,b,c,d,J) * hw / real(A_mass,dp)
                                else
                                    ! The original indexes the O1B struct as T[a,c] here ... T.p == T.n is used ...
                                    T_sym = (hw * T2b(Orb,a,b,c,d,J) + T%p(a,c) * real(kronecker_delta(b,d),dp) + &
                                             T%p(b,d) * real(kronecker_delta(a,c),dp)) / real(A_mass,dp)
                                end if
                                V_NN%pn(P,J)%v(Ind) = V_NN%pn(P,J)%v(Ind) + T_sym
                            end if

                            if (Bra <= N_t1) then
                                Ind = packed_index(Bra,Ket,N_t1)
                                a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                                b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                                c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                                d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                                j_c = Orb(c)%j
                                j_d = Orb(d)%j
                                if (CMS == 'CMS1+2B') then
                                    T_antisym = 1.0_dp / sqrt(real(1 + kronecker_delta(a,b),dp) * real(1 + kronecker_delta(c,d),dp)) * hw * &
                                                (T2b(Orb,a,b,c,d,J) - phase((j_c + j_d)/2 - J) * T2b(Orb,a,b,d,c,J)) / real(A_mass,dp)
                                else
                                    Amp = hw / sqrt(real(1 + kronecker_delta(a,b),dp) * real(1 + kronecker_delta(c,d),dp))
                                    Amp_2 = 1.0_dp / sqrt(real(1 + kronecker_delta(a,b),dp) * real(1 + kronecker_delta(c,d),dp))
                                    T_antisym = (Amp * T2b(Orb,a,b,c,d,J) + Amp_2 * (kronecker_delta(b,d) * T%p(a,c) + &
                                                kronecker_delta(a,c) * T%p(b,d)) - phase((j_c + j_d)/2 - J) * &
                                                (Amp * T2b(Orb,a,b,d,c,J) + Amp_2 * (real(kronecker_delta(b,c),dp) * &
                                                T%p(a,d) + real(kronecker_delta(a,d),dp) * T%p(b,c)))) / real(A_mass,dp)
                                end if
                                ME = T_antisym * sqrt(real(1 + kronecker_delta(a,b),dp) * real(1 + kronecker_delta(c,d),dp))
                                V_NN%pp(P,J)%v(Ind) = V_NN%pp(P,J)%v(Ind) + ME
                                V_NN%nn(P,J)%v(Ind) = V_NN%nn(P,J)%v(Ind) + ME
                            end if
                        end do
                    end do
                    !$omp end parallel do
                end do
            end do
            write(*,'(/,a)') 'CMS 2-body correction added...'
        end if

        ! Lawson term beta N_cm (CM_Lawson) ... pair operator beta/A [b+(1).b(2) + b+(2).b(1) + (N_1 + N_2)/(A-1)], the 1-body
            ! part folded in by 1/(A-1) (exact for A particles) ... b+(1).b(2) + b+(2).b(1) = r1.r2/b^2 + p1.p2 b^2/hbar^2 = 2 p1.p2
            ! b^2/hbar^2 for e_a + e_b = e_c + e_d (one quantum up, one down), 0 otherwise ... T2b = -<p1.p2> b^2/hbar^2 (the sign of
            ! the CMS1+2B correction -p1.p2/(m A)) ...
        if (abs(Params%Calc%CM_Lawson) > 0.0_dp) then
            write(*,'(/,a)') 'Adding the Lawson term beta N_cm, beta = ' // jstr(Params%Calc%CM_Lawson) // ' MeV ...'
            do J = 0, J_max
                do P = 1, 2
                    N_t0 = Orb_NN%N(1,P,J)
                    N_t1 = Orb_NN%N(2,P,J)
                    !$omp parallel do schedule(dynamic) private(Ket,Ind,a,b,c,d,j_c,j_d,ME)
                    do Bra = 1, max(N_t0,N_t1)
                        do Ket = 1, Bra
                            if (Bra <= N_t0) then
                                Ind = packed_index(Bra,Ket,N_t0)
                                a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                                b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                                c = Orb_NN%Ind(1,P,J)%v(1,Ket)
                                d = Orb_NN%Ind(1,P,J)%v(2,Ket)
                                V_NN%pn(P,J)%v(Ind) = V_NN%pn(P,J)%v(Ind) + Params%Calc%CM_Lawson / real(A_mass,dp) * &
                                                      xcm(a,b,c,d)
                            end if
                            if (Bra <= N_t1) then
                                Ind = packed_index(Bra,Ket,N_t1)
                                a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                                b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                                c = Orb_NN%Ind(2,P,J)%v(1,Ket)
                                d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                                j_c = Orb(c)%j
                                j_d = Orb(d)%j
                                ME = Params%Calc%CM_Lawson / real(A_mass,dp) * (xcm(a,b,c,d) - phase((j_c + j_d)/2 - J) * &
                                     xcm(a,b,d,c))
                                V_NN%pp(P,J)%v(Ind) = V_NN%pp(P,J)%v(Ind) + ME
                                V_NN%nn(P,J)%v(Ind) = V_NN%nn(P,J)%v(Ind) + ME
                            end if
                        end do
                    end do
                    !$omp end parallel do
                end do
            end do
            write(*,'(/,a)') 'Lawson term added...'
        end if

        ! Enhance contact V_aabb^J=0 matrix elements ...
        pG = Params%Int%pG2N
        nG = Params%Int%nG2N
        NN_min = -1
        NN_max = 25
        if (abs(pG) > 1.0e-8_dp .or. abs(nG) > 1.0e-8_dp) then
            N_T1_00 = Orb_NN%N(2,1,0)
            do Bra = 1, N_T1_00
                a = Orb_NN%Ind(2,1,0)%v(1,Bra)
                b = Orb_NN%Ind(2,1,0)%v(2,Bra)
                if (a == b) then
                    l_a = Orb(a)%l
                    j_a = Orb(a)%j
                    NN_a = 2*Orb(a)%n + l_a
                    if (NN_a < NN_min .or. NN_a > NN_max) cycle
                    do Ket = 1, Bra
                        c = Orb_NN%Ind(2,1,0)%v(1,Ket)
                        d = Orb_NN%Ind(2,1,0)%v(2,Ket)
                        if (c == d) then
                            l_c = Orb(c)%l
                            j_c = Orb(c)%j
                            NN_c = 2*Orb(c)%n + l_c
                            if (NN_c < NN_min .or. NN_c > NN_max) exit
                            Ind = packed_index(Bra,Ket,N_T1_00)
                            Amp = -0.5_dp * sqrt(real((j_a + 1) * (j_c + 1),dp)) * phase(l_a + l_c)
                            V_NN%pp(1,0)%v(Ind) = V_NN%pp(1,0)%v(Ind) + Amp * pG
                            V_NN%nn(1,0)%v(Ind) = V_NN%nn(1,0)%v(Ind) + Amp * nG
                        end if
                    end do
                end if
            end do
        end if

    contains

        ! direct pair element of b+(1).b(2) + b+(2).b(1) + (N_1 + N_2)/(A-1), coupled to J (not antisymmetrized) ...
        real(dp) function xcm(i1,i2,i3,i4)
            integer, intent(in) :: i1, i2, i3, i4
            integer :: e1, e2, e3, e4
            e1 = 2*Orb(i1)%n + Orb(i1)%l
            e2 = 2*Orb(i2)%n + Orb(i2)%l
            e3 = 2*Orb(i3)%n + Orb(i3)%l
            e4 = 2*Orb(i4)%n + Orb(i4)%l
            xcm = 0.0_dp
            if (e1 + e2 == e3 + e4) xcm = -2.0_dp * T2b(Orb,i1,i2,i3,i4,J)
            if (i1 == i3 .and. i2 == i4) xcm = xcm + real(e1 + e2,dp) / real(A_mass - 1,dp)
        end function xcm

        subroutine add_pn(i1,i2,i3,i4,Val)
            integer, intent(in) :: i1, i2, i3, i4
            real(dp), intent(in) :: Val
            integer :: Ind_pn
            Ind_pn = O2b_index(i1,i2,i3,i4,J,P,0,Orb_NN)
            V_NN%pn(P,J)%v(Ind_pn) = V_NN%pn(P,J)%v(Ind_pn) + Val
        end subroutine add_pn

    end subroutine V2b_read

    ! Number of Float64 entries stored for each leading orbital a ...
    function V2b_count(N_max,N_2max,Orb) result(N_Chunk)
        integer, intent(in) :: N_max, N_2max
        type(orb1b_t), intent(in) :: Orb(:)
        integer(i8), allocatable :: N_Chunk(:)
        integer :: a_max, a, b, c, d, d_max, J
        integer(i8) :: N

        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(N_Chunk(a_max))
        do a = 1, a_max
            N = 0
            do b = 1, a
                if ((2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l) <= N_2max) then
                    do c = 1, a
                        d_max = c
                        if (a == c) d_max = b
                        do d = 1, d_max
                            if (((2*(Orb(c)%n + Orb(d)%n) + Orb(c)%l + Orb(d)%l) <= N_2max) .and. &
                                (mod(Orb(a)%l + Orb(b)%l,2) == mod(Orb(c)%l + Orb(d)%l,2))) then
                                do J = max(abs(Orb(a)%j - Orb(b)%j), abs(Orb(c)%j - Orb(d)%j))/2, &
                                       min(Orb(a)%j + Orb(b)%j, Orb(c)%j + Orb(d)%j)/2
                                    N = N + 4
                                end do
                            end if
                        end do
                    end do
                end if
            end do
            N_Chunk(a) = N
        end do
    end function V2b_count

    ! Density-dependent residual NN interaction ... NO2B contraction of the NNN interaction with Rho ... V_res = (V_NN + Sum)*cRes ...
        ! Non-mutating as in dev (V_NN_res = deepcopy(V_NN)) ... the bare V_NN is kept ...
    subroutine V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(in) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
        type(o2b_t), intent(out) :: V_NN_res
        integer :: N_max, N_3max, a_max, J_max, J, P, N_T0, N_T1, Bra_ind, Bra, Ket, Ind, a, b, c, d, e, f, P3B
        integer :: n_a, l_a, n_b, l_b, n_c, l_c, j_c, n_d, l_d, n_e, l_e, n_f, l_f, j_f
        real(dp) :: cRes, isr12, i6, i3, Hat, pnSum, ppSum, nnSum
        real(dp) :: ME001, ME101, ME011, ME111, ME113

        N_max = Params%Calc%Nmax
        N_3max = Params%Calc%N3max
        a_max = (N_max + 1)*(N_max + 2)/2
        J_max = Params%Calc%N2max + 1
        cRes = Params%Int%cRes
        isr12 = 1.0_dp / sqrt(12.0_dp)
        i6 = 1.0_dp / 6.0_dp
        i3 = 1.0_dp / 3.0_dp

        V_NN_res = V_NN

        if (Verbose) write(*,'(/,a,/)') 'Starting calculation of residual NN interaction ...'
        if (Verbose) write(*,'(/,a)') 'Making density dependent residual 2-body interaction...'
        do J = 0, J_max
            do P = 1, 2
                if (P == 1) then
                    if (Verbose) write(*,'(a)') TAB // 'Calculating   ...   J = ' // istr(J) // '/' // istr(Params%Calc%N2max+1) // TAB // 'P = +'
                else
                    if (Verbose) write(*,'(a)') TAB // 'Calculating   ...   J = ' // istr(J) // '/' // istr(Params%Calc%N2max+1) // TAB // 'P = -'
                end if
                N_T0 = Orb_NN%N(1,P,J)
                N_T1 = Orb_NN%N(2,P,J)
                Hat = 1.0_dp / real(2*J + 1,dp)

                !$omp parallel do schedule(dynamic) private(Bra,a,b,n_a,l_a,n_b,l_b,Ket,Ind,d,e,n_d,l_d,n_e,l_e, &
                !$omp pnSum,ppSum,nnSum,c,n_c,l_c,j_c,P3B,f,n_f,l_f,j_f,ME001,ME101,ME011,ME111,ME113)
                do Bra_ind = 1, N_T0 + N_T1
                    ! pnV loop ...
                    if (Bra_ind <= N_T0) then
                        Bra = Bra_ind
                        a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                        b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                        n_a = Orb(a)%n
                        l_a = Orb(a)%l
                        n_b = Orb(b)%n
                        l_b = Orb(b)%l
                        do Ket = 1, Bra
                            Ind = packed_index(Bra,Ket,N_T0)
                            d = Orb_NN%Ind(1,P,J)%v(1,Ket)
                            e = Orb_NN%Ind(1,P,J)%v(2,Ket)
                            n_d = Orb(d)%n
                            l_d = Orb(d)%l
                            n_e = Orb(e)%n
                            l_e = Orb(e)%l
                            pnSum = 0.0_dp
                            do c = 1, a_max
                                n_c = Orb(c)%n
                                l_c = Orb(c)%l
                                j_c = Orb(c)%j
                                if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max) then
                                    P3B = mod(l_a + l_b + l_c, 2) + 1
                                    do f = 1, a_max
                                        n_f = Orb(f)%n
                                        l_f = Orb(f)%l
                                        j_f = Orb(f)%j
                                        if ((l_c == l_f) .and. (j_c == j_f) .and. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)) then
                                            ME001 = V3b_no2b(a,b,c,0,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                            ME101 = V3b_no2b(a,b,c,1,d,e,f,0,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                            ME011 = V3b_no2b(a,b,c,0,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                            ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                            ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)
                                            pnSum = pnSum + Hat * ((0.5_dp * ME001 - ME101 * isr12 - ME011 * isr12 + &
                                                    ME111 * i6 + ME113 * i3) * Rho%p(f,c) + &
                                                    (0.5_dp * ME001 + ME101 * isr12 + ME011 * isr12 + &
                                                    ME111 * i6 + ME113 * i3) * Rho%n(f,c))
                                        end if
                                    end do
                                end if
                            end do
                            V_NN_res%pn(P,J)%v(Ind) = V_NN_res%pn(P,J)%v(Ind) + pnSum
                            V_NN_res%pn(P,J)%v(Ind) = V_NN_res%pn(P,J)%v(Ind) * cRes
                        end do
                    end if

                    ! ppV & nnV loop ...
                    if (Bra_ind > N_T0) then
                        Bra = Bra_ind - N_T0
                        a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                        b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                        n_a = Orb(a)%n
                        l_a = Orb(a)%l
                        n_b = Orb(b)%n
                        l_b = Orb(b)%l
                        do Ket = 1, Bra
                            Ind = packed_index(Bra,Ket,N_T1)
                            d = Orb_NN%Ind(2,P,J)%v(1,Ket)
                            e = Orb_NN%Ind(2,P,J)%v(2,Ket)
                            n_d = Orb(d)%n
                            l_d = Orb(d)%l
                            n_e = Orb(e)%n
                            l_e = Orb(e)%l
                            ppSum = 0.0_dp
                            nnSum = 0.0_dp
                            do c = 1, a_max
                                n_c = Orb(c)%n
                                l_c = Orb(c)%l
                                j_c = Orb(c)%j
                                if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max) then
                                    P3B = mod(l_a + l_b + l_c, 2) + 1
                                    do f = 1, a_max
                                        n_f = Orb(f)%n
                                        l_f = Orb(f)%l
                                        j_f = Orb(f)%j
                                        if ((l_c == l_f) .and. (j_c == j_f) .and. ((2*(n_d + n_e + n_f) + l_d + l_e + l_f) <= N_3max)) then
                                            ME111 = V3b_no2b(a,b,c,1,d,e,f,1,J,1,P3B,V_NNN,Orb,Orb_NNN)
                                            ME113 = V3b_no2b(a,b,c,1,d,e,f,1,J,3,P3B,V_NNN,Orb,Orb_NNN)
                                            ppSum = ppSum + Hat * (ME113 * Rho%p(f,c) + (2.0_dp * ME111 + ME113) * Rho%n(f,c) * i3)
                                            nnSum = nnSum + Hat * (ME113 * Rho%n(f,c) + (2.0_dp * ME111 + ME113) * Rho%p(f,c) * i3)
                                        end if
                                    end do
                                end if
                            end do
                            V_NN_res%pp(P,J)%v(Ind) = V_NN_res%pp(P,J)%v(Ind) + ppSum
                            V_NN_res%nn(P,J)%v(Ind) = V_NN_res%nn(P,J)%v(Ind) + nnSum
                            V_NN_res%pp(P,J)%v(Ind) = V_NN_res%pp(P,J)%v(Ind) * cRes
                            V_NN_res%nn(P,J)%v(Ind) = V_NN_res%nn(P,J)%v(Ind) * cRes
                        end do
                    end if
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine V2b_residual_no2b

    ! Cuts Int.NN_File (space Int.Nmax/N2max) down to (N_max_Out,N_2max_Out) & writes it in the same NuHamil Float64 layout to
        ! Out_File ... counterpart of V2b_parse ... chunks a are processed in parallel (largest first), each thread with its own
        ! units, writing disjoint byte ranges of the pre-sized output file ...
    subroutine V2b_parse(Params,N_max_Out,N_2max_Out,Out_File)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_max_Out, N_2max_Out
        character(len=*), intent(in) :: Out_File
        type(orb1b_t), allocatable :: Orb(:)
        real(dp), allocatable :: Buffer_In(:), Buffer_Out(:), Sizes(:)
        integer(i8), allocatable :: N_Chunk_In(:), N_Chunk_Out(:), N_Skip_In(:), N_Skip_Out(:)
        integer, allocatable :: Perm(:)
        integer(i8) :: Size_In, Size_Out, File_size, Chunk_i, Buffer_0, Out_N, N_In, N_Out, Pos_Out
        integer :: N_max_In, N_2max_In, a_max_Out, a, b, c, d, d_max, q, u_in, u_out, Buffer_N, Out_i, i_J, N_J, Err
        integer :: J_min, J_max, NN_a, NN_b, NN_c, NN_d, NN_ab, NN_cd
        logical :: Failed
        character(len=:), allocatable :: In_File

        In_File = trim(Params%Int%NN_File)
        N_max_In = Params%Int%Nmax
        N_2max_In = Params%Int%N2max
        a_max_Out = (N_max_Out + 1)*(N_max_Out + 2)/2

        ! The output space has to be contained in the interaction space ...
        if ((N_max_Out > N_max_In) .or. (N_2max_Out > N_2max_In)) &
            call fatal('Output space (' // istr(N_max_Out) // ',' // istr(N_2max_Out) // &
                       ') is not contained in the interaction space (' // istr(N_max_In) // ',' // istr(N_2max_In) // ') ...')
        if (min(N_max_Out,N_2max_Out) < 0) call fatal('V2b_parse: negative output truncation ...')
        if (abspath_of(In_File) == abspath_of(Out_File)) call fatal('Output file has to differ from the interaction file ...')

        Orb = orbitals_make(Params,Use_Int=.true.)

        write(*,'(/,a)') 'Precounting # of 2-body NN interaction matrix elements ...'
        N_Chunk_In = V2b_count(N_max_In,N_2max_In,Orb)
        N_Chunk_Out = V2b_count(N_max_Out,N_2max_Out,Orb)
        allocate(N_Skip_In(size(N_Chunk_In)), N_Skip_Out(size(N_Chunk_Out)))
        Size_In = 0
        do a = 1, size(N_Chunk_In)
            N_Skip_In(a) = Size_In
            Size_In = Size_In + 8*N_Chunk_In(a)
        end do
        Size_Out = 0
        do a = 1, size(N_Chunk_Out)
            N_Skip_Out(a) = Size_Out
            Size_Out = Size_Out + 8*N_Chunk_Out(a)
        end do

        ! Size of the interaction file ...
        open(newunit=u_in, file=In_File, access='stream', form='unformatted', status='old', action='read', iostat=Err)
        if (Err /= 0) call fatal('Cannot open the NN interaction file ' // In_File)
        inquire(unit=u_in, size=File_size)
        close(u_in)
        if (File_size /= Size_In) call fatal('Size of ' // In_File // ' (' // istr(File_size) // &
                                             ' B) does not match the interaction space (' // istr(Size_In) // ' B) ...')

        call presize_file(Out_File,Size_Out)

        allocate(Sizes(a_max_Out), Perm(a_max_Out))
        Sizes = real(N_Chunk_Out,dp)
        call sortperm_rev(Sizes,Perm)

        write(*,'(/,a)') 'Exporting 2-body NN interaction file into N_max = ' // istr(N_max_Out) // ' , N_2max = ' // &
                         istr(N_2max_Out) // ' ...'
        Failed = .false.
        !$omp parallel default(shared) private(u_in,u_out,Buffer_In,Buffer_Out,q,a,b,c,d,d_max,i_J,N_J,J_min,J_max, &
        !$omp NN_a,NN_b,NN_c,NN_d,NN_ab,NN_cd,Chunk_i,Buffer_0,Buffer_N,Out_i,Out_N,N_In,N_Out,Pos_Out)
        open(newunit=u_in, file=In_File, access='stream', form='unformatted', status='old', action='read')
        open(newunit=u_out, file=trim(Out_File), access='stream', form='unformatted', status='old', action='readwrite')
        allocate(Buffer_In(Buffer_size), Buffer_Out(Buffer_size))
        !$omp do schedule(dynamic,1)
        do q = 1, a_max_Out
            a = Perm(q)
            N_In = N_Chunk_In(a)
            N_Out = N_Chunk_Out(a)
            NN_a = 2*Orb(a)%n + Orb(a)%l
            Chunk_i = 0
            Buffer_0 = 0
            Buffer_N = 0
            Out_i = 0
            Out_N = 0
            Pos_Out = N_Skip_Out(a) + 1
            do b = 1, a
                NN_b = 2*Orb(b)%n + Orb(b)%l
                NN_ab = NN_a + NN_b
                if (NN_ab > N_2max_In) cycle
                do c = 1, a
                    NN_c = 2*Orb(c)%n + Orb(c)%l
                    d_max = c
                    if (a == c) d_max = b
                    do d = 1, d_max
                        NN_d = 2*Orb(d)%n + Orb(d)%l
                        NN_cd = NN_c + NN_d
                        if (.not. ((NN_cd <= N_2max_In) .and. (mod(Orb(a)%l + Orb(b)%l,2) == mod(Orb(c)%l + Orb(d)%l,2)))) cycle
                        J_min = max(abs(Orb(a)%j - Orb(b)%j), abs(Orb(c)%j - Orb(d)%j))/2
                        J_max = min(Orb(a)%j + Orb(b)%j, Orb(c)%j + Orb(d)%j)/2
                        N_J = 4*max(J_max - J_min + 1,0)
                        if ((NN_a <= N_max_Out) .and. (NN_b <= N_max_Out) .and. (NN_c <= N_max_Out) .and. (NN_d <= N_max_Out) .and. &
                            (NN_ab <= N_2max_Out) .and. (NN_cd <= N_2max_Out)) then
                            ! Copy the 4 matrix elements (00,nn,10,pp) per J ...
                            do i_J = 1, N_J
                                Chunk_i = Chunk_i + 1
                                if (Chunk_i > Buffer_0 + Buffer_N) then
                                    Buffer_0 = Chunk_i - 1
                                    Buffer_N = int(min(int(Buffer_size,i8), N_In - Buffer_0))
                                    read(u_in, pos=N_Skip_In(a) + 8_i8*Buffer_0 + 1_i8) Buffer_In(1:Buffer_N)
                                end if
                                if (Out_i == Buffer_size) then
                                    write(u_out, pos=Pos_Out) Buffer_Out(1:Out_i)
                                    Pos_Out = Pos_Out + 8_i8*Out_i
                                    Out_i = 0
                                end if
                                Out_i = Out_i + 1
                                Out_N = Out_N + 1
                                Buffer_Out(Out_i) = Buffer_In(Chunk_i - Buffer_0)
                            end do
                        else
                            ! Skip the elements outside the output space ...
                            Chunk_i = Chunk_i + N_J
                        end if
                    end do
                end do
            end do
            if (Out_i > 0) write(u_out, pos=Pos_Out) Buffer_Out(1:Out_i)
            if ((Chunk_i /= N_In) .or. (Out_N /= N_Out)) then
                write(*,'(a)') 'ERROR: Wrong # of matrix elements in the chunk a = ' // istr(a) // ' ...'
                Failed = .true.
            end if
        end do
        !$omp end do
        close(u_in)
        close(u_out)
        deallocate(Buffer_In, Buffer_Out)
        !$omp end parallel
        if (Failed) call fatal('V2b_parse: wrong # of matrix elements ...')

        write(*,'(/,a)') 'Finished export of 2-body NN interaction into ' // trim(Out_File) // ' ...'
    end subroutine V2b_parse

    ! Correlation energy by contraction with the 2-body correlation function Sigma ... counterpart of V2b_correlation_energy ...
        ! E = sum_{abcd,J} (2J+1) [1/4 V_pp Sigma_pp + 1/4 V_nn Sigma_nn + V_pn Sigma_pn] over all ordered pairs ab & cd ...
    real(dp) function V2b_correlation_energy(Params,Orb,Orb_NN,Sigma,V_NN) result(E_corr)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: Sigma, V_NN
        integer :: N_2max, a_max, a, b, c, d, P, J, J_min, J_max
        real(dp) :: Sum, J_hat

        N_2max = Params%Calc%N2max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2

        write(*,'(/,a)') 'Computing the correlation ground-state energy by contraction with the provided 2-body correlation function ...'
        E_corr = 0.0_dp
        !$omp parallel do schedule(dynamic) reduction(+:E_corr) private(b,c,d,P,J,J_min,J_max,Sum,J_hat)
        do a = 1, a_max
            Sum = 0.0_dp
            do b = 1, a_max
                P = mod(Orb(a)%l + Orb(b)%l,2) + 1
                if ((2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l) > N_2max) cycle
                do c = 1, a_max
                    do d = 1, a_max
                        if (((2*(Orb(c)%n + Orb(d)%n) + Orb(c)%l + Orb(d)%l) > N_2max) .or. &
                            (mod(Orb(a)%l + Orb(b)%l,2) /= mod(Orb(c)%l + Orb(d)%l,2))) cycle
                        J_min = max(abs(Orb(a)%j - Orb(b)%j)/2, abs(Orb(c)%j - Orb(d)%j)/2)
                        J_max = min((Orb(a)%j + Orb(b)%j)/2, (Orb(c)%j + Orb(d)%j)/2)
                        do J = J_min, J_max
                            J_hat = real(2*J + 1,dp)
                            Sum = Sum + 0.25_dp * J_hat * O2b_pp(a,b,c,d,J,P,V_NN,Orb,Orb_NN) * O2b_pp(a,b,c,d,J,P,Sigma,Orb,Orb_NN)
                            Sum = Sum + 0.25_dp * J_hat * O2b_nn(a,b,c,d,J,P,V_NN,Orb,Orb_NN) * O2b_nn(a,b,c,d,J,P,Sigma,Orb,Orb_NN)
                            Sum = Sum + J_hat * O2b_pn(a,b,c,d,J,P,V_NN,Orb_NN) * O2b_pn(a,b,c,d,J,P,Sigma,Orb_NN)
                        end do
                    end do
                end do
            end do
            E_corr = E_corr + Sum
        end do
        !$omp end parallel do
        write(*,'(a,/)') TAB // 'Performed calculation of correlation energy by contraction with 2-body correlation function ...'
    end function V2b_correlation_energy

end module nj_v2b
