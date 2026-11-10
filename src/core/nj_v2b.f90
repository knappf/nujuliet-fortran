! 2-body NN interaction ... counterpart of MatrixElements/V2B.jl ...
module nj_v2b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_o1b, only : T1b, T2b
    use nj_o2b
    use nj_v3b, only : V3b_no2b
    implicit none
    private

    public :: V2b_read, V2b_count, V2b_residual_no2b

contains

    ! Reads the NuHamil NN file (Float64 stream of T=0, nn, T=1 pn, pp per J) & adds the 2-body CM correction ...
    subroutine V2b_read(Params,Orb,V_NN,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o2b_t), intent(out) :: V_NN
        type(orb2b_t), intent(out) :: Orb_NN
        type(o1b_t) :: T
        integer :: N_max_Int, N_2max_Int, N_max_Calc, N_2max_Calc, A_mass, a_max, J_max
        integer :: u, a, b, c, d, d_max, J, P, Ind, k, Bra, Ket, N_t0, N_t1, N_T1_00
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_c, l_c, j_c, n_d, l_d, j_d
        integer :: NN_a, NN_b, NN_c, NN_d, NN_ab, NN_cd, NN_min, NN_max
        integer(i8), allocatable :: N_Chunk(:)
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

        ! Initialize 1-body kinetic operator ...
        T = T1b(Params,Orb)

        ! Initialize NN interaction array ...
        write(*,'(/,a)') 'Preparing 2-body NN interaction array...'
        call O2b_initialize(Params,Orb,V_NN,Orb_NN)

        ! Count the matrix elements of each chunk a ...
        N_Chunk = V2b_count(N_max_Int,N_2max_Int,Orb)

        write(*,'(/,a)') 'Reading 2-body NN interaction file...'
        open(newunit=u, file=trim(Params%Int%NN_File), access='stream', form='unformatted', status='old', action='read')
        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            NN_a = 2*n_a + l_a
            allocate(Buf(N_Chunk(a)))
            if (N_Chunk(a) > 0) read(u) Buf
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
            deallocate(Buf)
        end do
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

    ! Density-dependent residual NN interaction ... NO2B contraction of the NNN interaction with Rho ...
    subroutine V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: Rho
        type(o2b_t), intent(inout) :: V_NN
        type(v3b_t), intent(in) :: V_NNN
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
                            V_NN%pn(P,J)%v(Ind) = V_NN%pn(P,J)%v(Ind) + pnSum
                            V_NN%pn(P,J)%v(Ind) = V_NN%pn(P,J)%v(Ind) * cRes
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
                            V_NN%pp(P,J)%v(Ind) = V_NN%pp(P,J)%v(Ind) + ppSum
                            V_NN%nn(P,J)%v(Ind) = V_NN%nn(P,J)%v(Ind) + nnSum
                            V_NN%pp(P,J)%v(Ind) = V_NN%pp(P,J)%v(Ind) * cRes
                            V_NN%nn(P,J)%v(Ind) = V_NN%nn(P,J)%v(Ind) * cRes
                        end do
                    end if
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine V2b_residual_no2b

end module nj_v2b
