! 3-body NNN interaction in the normal-ordered 2-body (NO2B) storage ... counterpart of MatrixElements/V3B_NO2B.jl ...
    ! Blocks V(P,l_c,J,T_ind,j_index) hold packed lower triangles in Float32 ...
    ! Orb_NNN%Dic maps the packed key of (a,b,c,J,P,T_ab,T) to the index within its block ...
module nj_v3b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_hashmap
    implicit none
    private

    public :: V3b_no2b_initialize, V3b_no2b_key, V3b_no2b_index, V3b_no2b, V3b_no2b_read

contains

    subroutine V3b_no2b_initialize(N_max,N_2max,N_3max,Orb,V_NNN,Orb_NNN)
        integer, intent(in) :: N_max, N_2max, N_3max
        type(orb1b_t), intent(in) :: Orb(:)
        type(v3b_t), intent(out) :: V_NNN
        type(orb3b_t), intent(out) :: Orb_NNN
        integer :: a_max, l_max, J_max, a, b, c, P, J, T, T_ab, T_tot, j_ind, l, N_count
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_c, l_c, j_c
        integer(i8) :: Key

        a_max = (N_max + 1)*(N_max + 2)/2
        l_max = N_max
        J_max = N_2max + 1

        allocate(Orb_NNN%N(2,0:l_max,0:J_max,2,2))
        Orb_NNN%N = 0
        call Orb_NNN%Dic%init(4096)

        do a = 1, a_max
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            if ((2*n_a + l_a) <= N_max) then
                do b = 1, a
                    n_b = Orb(b)%n
                    l_b = Orb(b)%l
                    j_b = Orb(b)%j
                    if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                        do c = 1, a_max
                            n_c = Orb(c)%n
                            l_c = Orb(c)%l
                            j_c = Orb(c)%j
                            P = mod(l_a + l_b + l_c,2) + 1
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max .and. &
                                (2*(n_a + n_c) + l_a + l_c) <= N_2max .and. &
                                (2*(n_b + n_c) + l_b + l_c) <= N_2max) then
                                do J = abs(j_a - j_b)/2, (j_a + j_b)/2
                                    do T = 1, 2
                                        do T_ab = T - 1, 1
                                            T_tot = 2*T - 1
                                            Key = V3b_no2b_key(a,b,n_c,l_c,j_c,J,P,T_ab,T_tot)
                                            j_ind = jindex(l_c,j_c)
                                            Orb_NNN%N(P,l_c,J,T,j_ind) = Orb_NNN%N(P,l_c,J,T,j_ind) + 1
                                            call Orb_NNN%Dic%set(Key, Orb_NNN%N(P,l_c,J,T,j_ind))
                                        end do
                                    end do
                                end do
                            end if
                        end do
                    end if
                end do
            end if
        end do

        allocate(V_NNN%V(2,0:l_max,0:J_max,2,2))
        do P = 1, 2
            do l = 0, l_max
                do J = 0, J_max
                    do T = 1, 2
                        do j_ind = 1, 2
                            N_count = Orb_NNN%N(P,l,J,T,j_ind)
                            allocate(V_NNN%V(P,l,J,T,j_ind)%v(N_count*(N_count+1)/2))
                            V_NNN%V(P,l,J,T,j_ind)%v = 0.0_sp
                        end do
                    end do
                end do
            end do
        end do
    end subroutine V3b_no2b_initialize

    ! j index within the (l, j = l -+ 1/2) pair ... 1 for j = |2l-1|, 2 for j = 2l+1 ...
    pure integer function jindex(l_c,j_c)
        integer, intent(in) :: l_c, j_c
        jindex = merge(1, 2, j_c == abs(2*l_c - 1))
    end function jindex

    pure integer(i8) function V3b_no2b_key(a,b,n_c,l_c,j_c,J,P,T_ab,T)
        integer, intent(in) :: a, b, n_c, l_c, j_c, J, P, T_ab, T
        V3b_no2b_key = ior(ior(ior(ior(ior(ior(ior(ior(int(J,i8), shiftl(int(P,i8),6)), shiftl(int(T_ab,i8),8)), &
                       shiftl(int(T,i8),10)), shiftl(int(a,i8),12)), shiftl(int(b,i8),20)), shiftl(int(n_c,i8),28)), &
                       shiftl(int(l_c,i8),34)), shiftl(int(j_c,i8),40))
    end function V3b_no2b_key

    integer function V3b_no2b_index(a,b,c,T_ab,d,e,f,T_de,J,T,P,Orb,Orb_NNN)
        integer, intent(in) :: a, b, c, T_ab, d, e, f, T_de, J, T, P
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb3b_t), intent(in) :: Orb_NNN
        integer :: Bra, Ket, Tmp, N, n_c, l_c, j_c
        n_c = Orb(c)%n
        l_c = Orb(c)%l
        j_c = Orb(c)%j
        ! Non-polymorphic hashmap_lookup, not the type-bound Dic%get: in the hot threaded loops ifx built the class
        ! descriptor in static memory on every call (all threads writing the same cache lines) ... 40Ca, 9 shells,
        ! 16 threads: V2b_residual_no2b 51.8 s -> 2.3 s, HF_RRPA 2796 s -> 166 s, outputs identical to 1e-12 ...
        ! a missing key (index -1) is caught below, as Julia's KeyError ...
        Bra = hashmap_lookup(Orb_NNN%Dic, V3b_no2b_key(a,b,n_c,l_c,j_c,J,P,T_ab,T))
        Ket = hashmap_lookup(Orb_NNN%Dic, V3b_no2b_key(d,e,Orb(f)%n,l_c,j_c,J,P,T_de,T))
        if (min(Bra,Ket) < 0) call V3b_no2b_missing_key(a,b,c,T_ab,d,e,f,T_de,J,T,P)
        if (Bra < Ket) then
            Tmp = Bra
            Bra = Ket
            Ket = Tmp
        end if
        N = Orb_NNN%N(P,l_c,J,(T+1)/2,jindex(l_c,j_c))
        V3b_no2b_index = Bra + (Ket - 1)*N - (Ket * (Ket - 1))/2
    end function V3b_no2b_index

    ! Cold error path of V3b_no2b_index ... kept out of line so the lookup stays small ...
    subroutine V3b_no2b_missing_key(a,b,c,T_ab,d,e,f,T_de,J,T,P)
        integer, intent(in) :: a, b, c, T_ab, d, e, f, T_de, J, T, P
        write(*,'(a,11(1x,I0))') 'ERROR: missing key in NNN hash map lookup ... a,b,c,T_ab,d,e,f,T_de,J,T,P =', &
                                 a, b, c, T_ab, d, e, f, T_de, J, T, P
        error stop 1
    end subroutine V3b_no2b_missing_key

    ! NO2B matrix element <ab(J_de T_ab) c | V | de(J_de T_de) f> with total isospin T (doubled) ...
    real(dp) function V3b_no2b(a,b,c,T_ab,d,e,f,T_de,J_de,T,P,V_NNN,Orb,Orb_NNN)
        integer, intent(in) :: a, b, c, T_ab, d, e, f, T_de, J_de, T, P
        type(v3b_t), intent(in) :: V_NNN
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb3b_t), intent(in) :: Orb_NNN
        integer :: Ind
        real(dp) :: Ph
        if (a >= b) then
            if (d >= e) then
                Ind = V3b_no2b_index(a,b,c,T_ab,d,e,f,T_de,J_de,T,P,Orb,Orb_NNN)
                Ph = 1.0_dp
            else
                Ind = V3b_no2b_index(a,b,c,T_ab,e,d,f,T_de,J_de,T,P,Orb,Orb_NNN)
                Ph = phase((Orb(d)%j + Orb(e)%j)/2 - J_de - T_de)
            end if
        else
            if (d >= e) then
                Ind = V3b_no2b_index(b,a,c,T_ab,d,e,f,T_de,J_de,T,P,Orb,Orb_NNN)
                Ph = phase((Orb(a)%j + Orb(b)%j)/2 - J_de - T_ab)
            else
                Ind = V3b_no2b_index(b,a,c,T_ab,e,d,f,T_de,J_de,T,P,Orb,Orb_NNN)
                Ph = phase((Orb(a)%j + Orb(b)%j + Orb(d)%j + Orb(e)%j)/2 - T_de - T_ab)
            end if
        end if
        V3b_no2b = Ph * real(V_NNN%V(P,Orb(c)%l,J_de,(T+1)/2,jindex(Orb(c)%l,Orb(c)%j))%v(Ind),dp)
    end function V3b_no2b

    ! Number of Float32 matrix elements of the (a,b) chunk of the NO2B file ... counterpart of V3b_no2b_count ...
    integer(i8) function V3b_no2b_chunk_count(a,b,N_max,N_2max,N_3max,Orb) result(N)
        integer, intent(in) :: a, b, N_max, N_2max, N_3max
        type(orb1b_t), intent(in) :: Orb(:)
        integer :: a_max, c, d, e, f, NN_a, NN_b, NN_c, NN_d, NN_e, NN_f, P_abc, P_def, N_JT
        a_max = (N_max + 1)*(N_max + 2)/2
        NN_a = 2*Orb(a)%n + Orb(a)%l
        NN_b = 2*Orb(b)%n + Orb(b)%l
        N = 0
        do c = 1, a_max
            NN_c = 2*Orb(c)%n + Orb(c)%l
            if (.not. ((NN_a + NN_c <= N_2max) .and. (NN_b + NN_c <= N_2max) .and. (NN_a + NN_b + NN_c <= N_3max))) cycle
            P_abc = mod(Orb(a)%l + Orb(b)%l + Orb(c)%l,2) + 1
            do d = 1, a
                NN_d = 2*Orb(d)%n + Orb(d)%l
                do e = 1, d
                    NN_e = 2*Orb(e)%n + Orb(e)%l
                    if (NN_d + NN_e > N_2max) cycle
                    ! (T_ab,T_de,T) combinations per J: 1 + 1 + 1 + 2 = 5 ...
                    N_JT = 5*max(0, min(Orb(a)%j + Orb(b)%j, Orb(d)%j + Orb(e)%j)/2 &
                                    - max(abs(Orb(a)%j - Orb(b)%j), abs(Orb(d)%j - Orb(e)%j))/2 + 1)
                    do f = 1, a_max
                        NN_f = 2*Orb(f)%n + Orb(f)%l
                        if (.not. ((Orb(c)%j == Orb(f)%j) .and. (Orb(c)%l == Orb(f)%l) .and. (NN_d + NN_f <= N_2max) .and. &
                                   (NN_e + NN_f <= N_2max) .and. (NN_d + NN_e + NN_f <= N_3max))) cycle
                        P_def = mod(Orb(d)%l + Orb(e)%l + Orb(f)%l,2) + 1
                        if (P_abc /= P_def) cycle
                        N = N + N_JT
                    end do
                end do
            end do
        end do
    end function V3b_no2b_chunk_count

    ! Reads the NuHamil NO2B 3-body file (Float32 stream) ... as in the original, the file is split into
    ! (a,b) chunks whose offsets are precounted, and the chunks are read in parallel. Chunks beyond the last
    ! one that can hold calculation-space elements (orbitals are ordered by N) are not read at all ...
    subroutine V3b_no2b_read(Params,Orb,V_NNN,Orb_NNN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(v3b_t), intent(out) :: V_NNN
        type(orb3b_t), intent(out) :: Orb_NNN
        integer, parameter :: Buffer_size = 1000000
        real(sp), allocatable :: Buffer(:)
        integer(i8), allocatable :: N_Chunk(:), N_Chunk_skip(:)
        integer, allocatable :: ab(:,:)
        integer(i8) :: N_left, File_size, Pos
        integer :: N_max_Int, N_2max_Int, N_3max_Int, a_max, N_max_Calc, N_2max_Calc, N_3max_Calc
        integer :: u, a, b, c, d, e, f, J, T_ab, T_de, T, P, Ind, Buffer_i, Buffer_N, ab_i, ab_max, ab_last
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_c, l_c, j_c, n_d, l_d, j_d, n_e, l_e, j_e, n_f, l_f, j_f
        integer :: NN_a, NN_b, NN_c, NN_d, NN_e, NN_f, NN_ab, NN_ac, NN_bc, NN_abc, NN_de, NN_df, NN_ef, NN_def
        integer :: P_abc, P_def
        real(sp) :: Vme

        N_max_Int = Params%Int%Nmax
        N_2max_Int = Params%Int%N2max
        N_3max_Int = Params%Int%N3max
        a_max = (N_max_Int + 1)*(N_max_Int + 2)/2
        N_max_Calc = Params%Calc%Nmax
        N_2max_Calc = Params%Calc%N2max
        N_3max_Calc = Params%Calc%N3max

        write(*,'(/,a)') 'Preparing 3-body NNN interaction array ...'
        call V3b_no2b_initialize(N_max_Calc,N_2max_Calc,N_3max_Calc,Orb,V_NNN,Orb_NNN)
        ! No 3N interaction (NNN_File = '') ... the NO2B elements stay zero ...
        if (len_trim(Params%Int%NNN_File) == 0) then
            write(*,'(a)') TAB // 'NNN_File empty ... no 3-body interaction (NO2B elements zero)'
            return
        end if

        ! (a,b) chunks of the file ... counterpart of V3b_no2b_allocate ...
        ab_max = 0
        do a = 1, a_max
            if (2*Orb(a)%n + Orb(a)%l > N_max_Int) cycle
            do b = 1, a
                if (2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l <= N_2max_Int) ab_max = ab_max + 1
            end do
        end do
        allocate(ab(ab_max,2))
        ab_max = 0
        ab_last = 0
        do a = 1, a_max
            if (2*Orb(a)%n + Orb(a)%l > N_max_Int) cycle
            do b = 1, a
                if (2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l <= N_2max_Int) then
                    ab_max = ab_max + 1
                    ab(ab_max,:) = [a, b]
                    if (2*Orb(a)%n + Orb(a)%l <= N_max_Calc) ab_last = ab_max
                end if
            end do
        end do

        write(*,'(/,a)') 'Precounting # of 3-body NNN interaction matrix elements ...'
        allocate(N_Chunk(ab_last), N_Chunk_skip(ab_last))
        !$omp parallel do schedule(dynamic)
        do ab_i = 1, ab_last
            N_Chunk(ab_i) = V3b_no2b_chunk_count(ab(ab_i,1),ab(ab_i,2),N_max_Int,N_2max_Int,N_3max_Int,Orb)
        end do
        !$omp end parallel do
        N_Chunk_skip(1) = 0
        do ab_i = 2, ab_last
            N_Chunk_skip(ab_i) = N_Chunk_skip(ab_i-1) + N_Chunk(ab_i-1)
        end do

        write(*,'(/,a)') 'Reading 3-body NO2B NNN interaction file. ..'
        open(newunit=u, file=trim(Params%Int%NNN_File), access='stream', form='unformatted', status='old', action='read')
        inquire(unit=u, size=File_size)
        close(u)
        if (ab_last > 0) then
            N_left = 4*(N_Chunk_skip(ab_last) + N_Chunk(ab_last))
            if ((N_left > File_size) .or. ((ab_last == ab_max) .and. (N_left /= File_size))) &
                call fatal('Size of the NNN interaction file does not match the interaction space ' // trim(Params%Int%NNN_File))
        end if

        !$omp parallel default(shared) private(u,Buffer,N_left,Pos,a,b,c,d,e,f,J,T_ab,T_de,T,P,Ind,Buffer_i,Buffer_N, &
        !$omp n_a,l_a,j_a,n_b,l_b,j_b,n_c,l_c,j_c,n_d,l_d,j_d,n_e,l_e,j_e,n_f,l_f,j_f, &
        !$omp NN_a,NN_b,NN_c,NN_d,NN_e,NN_f,NN_ab,NN_ac,NN_bc,NN_abc,NN_de,NN_df,NN_ef,NN_def,P_abc,P_def,Vme)
        open(newunit=u, file=trim(Params%Int%NNN_File), access='stream', form='unformatted', status='old', action='read')
        allocate(Buffer(Buffer_size))
        !$omp do schedule(dynamic)
        do ab_i = 1, ab_last
            a = ab(ab_i,1)
            b = ab(ab_i,2)
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            NN_a = 2*n_a + l_a
            n_b = Orb(b)%n
            l_b = Orb(b)%l
            j_b = Orb(b)%j
            NN_b = 2*n_b + l_b
            NN_ab = NN_a + NN_b
            ! Nothing of this chunk belongs to the calculation space ...
            if ((NN_a > N_max_Calc) .or. (NN_b > N_max_Calc) .or. (NN_ab > N_2max_Calc)) cycle
            N_left = N_Chunk(ab_i)
            Pos = 4*N_Chunk_skip(ab_i) + 1
            Buffer_i = 1
            Buffer_N = 0
            do c = 1, a_max
                n_c = Orb(c)%n
                l_c = Orb(c)%l
                j_c = Orb(c)%j
                NN_c = 2*n_c + l_c
                NN_ac = NN_a + NN_c
                NN_bc = NN_b + NN_c
                NN_abc = NN_a + NN_b + NN_c
                P = mod(l_a + l_b + l_c,2) + 1
                if (.not. ((NN_ac <= N_2max_Int) .and. (NN_bc <= N_2max_Int) .and. (NN_abc <= N_3max_Int))) cycle
                P_abc = P
                do d = 1, a
                    n_d = Orb(d)%n
                    l_d = Orb(d)%l
                    j_d = Orb(d)%j
                    NN_d = 2*n_d + l_d
                    do e = 1, d
                        n_e = Orb(e)%n
                        l_e = Orb(e)%l
                        j_e = Orb(e)%j
                        NN_e = 2*n_e + l_e
                        NN_de = NN_d + NN_e
                        if (NN_de > N_2max_Int) cycle
                        do f = 1, a_max
                            n_f = Orb(f)%n
                            l_f = Orb(f)%l
                            j_f = Orb(f)%j
                            NN_f = 2*n_f + l_f
                            NN_df = NN_d + NN_f
                            NN_ef = NN_e + NN_f
                            NN_def = NN_d + NN_e + NN_f
                            if (.not. ((j_c == j_f) .and. (l_c == l_f) .and. (NN_df <= N_2max_Int) .and. &
                                       (NN_ef <= N_2max_Int) .and. (NN_def <= N_3max_Int))) cycle
                            P_def = mod(l_d + l_e + l_f,2) + 1
                            if (P_abc /= P_def) cycle
                            do J = max(abs(j_a - j_b),abs(j_d - j_e))/2, min((j_a + j_b),(j_d + j_e))/2
                                do T_ab = 0, 1
                                    do T_de = 0, 1
                                        do T = max(abs(2*T_ab-1),abs(2*T_de-1)), min((2*T_ab+1),(2*T_de+1)), 2
                                            if (Buffer_i > Buffer_N) then
                                                Buffer_N = int(min(N_left, int(Buffer_size,i8)))
                                                read(u, pos=Pos) Buffer(1:Buffer_N)
                                                Pos = Pos + 4_i8*Buffer_N
                                                N_left = N_left - Buffer_N
                                                Buffer_i = 1
                                            end if
                                            Vme = Buffer(Buffer_i)
                                            Buffer_i = Buffer_i + 1

                                            if ((NN_c <= N_max_Calc) .and. (NN_d <= N_max_Calc) .and. &
                                                (NN_e <= N_max_Calc) .and. (NN_f <= N_max_Calc) .and. &
                                                (NN_ac <= N_2max_Calc) .and. (NN_bc <= N_2max_Calc) .and. &
                                                (NN_de <= N_2max_Calc) .and. (NN_df <= N_2max_Calc) .and. (NN_ef <= N_2max_Calc) .and. &
                                                (NN_abc <= N_3max_Calc) .and. (NN_def <= N_3max_Calc)) then
                                                Ind = V3b_no2b_index(a,b,c,T_ab,d,e,f,T_de,J,T,P,Orb,Orb_NNN)
                                                V_NNN%V(P,l_c,J,(T+1)/2,jindex(l_c,j_c))%v(Ind) = Vme
                                            end if
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do
        close(u)
        deallocate(Buffer)
        !$omp end parallel

        write(*,'(/,a)') 'Finished loading of 3-body NO2B NNN interaction ...'
    end subroutine V3b_no2b_read

end module nj_v3b
