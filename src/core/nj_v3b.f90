! 3-body NNN interaction in the normal-ordered 2-body (NO2B) storage ... counterpart of src/common/V3b_no2b.jl (NuJuliet dev;
    ! MatrixElements/V3B_NO2B.jl before) ...
    ! Blocks V(P,l_c,J,T_ind,j_index) hold packed lower triangles in Float32, stored row by row (Ind = Ket + Bra(Bra-1)/2,
    ! Bra >= Ket) as in dev (base NuJuliet stored them column by column) ... the 2+1-body states (ab T_ab c) of a block are
    ! numbered by the dense map Orb_NNN%Map(ab)%m(T_ab+T_ind, J-J_min+1, c) (dev Orb3B.Map, the hash map of the base port is
    ! gone) ... packed indices are 64-bit integers (a block of N states has N(N+1)/2 elements, > 2^31 for N > 65535) ...
module nj_v3b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_orbitals, only : orbitals_make
    implicit none
    private

    public :: V3b_no2b_initialize, V3b_no2b_lookup, V3b_no2b_index, V3b_no2b, V3b_no2b_read, V3b_no2b_parse, V3b_jindex

    ! Buffer size (# of Float32 matrix elements read/written at once by each thread) ...
    integer, parameter :: Buffer_size = 1000000

contains

    ! Allocates the zero NO2B interaction blocks & the dense 2+1-body lookup ... counterpart of V3b_no2b_initialize ...
    subroutine V3b_no2b_initialize(N_max,N_2max,N_3max,Orb,V_NNN,Orb_NNN)
        integer, intent(in) :: N_max, N_2max, N_3max
        type(orb1b_t), intent(in) :: Orb(:)
        type(v3b_t), intent(out) :: V_NNN
        type(orb3b_t), intent(out) :: Orb_NNN
        integer :: a_max, l_max, J_max, a, b, c, P, J, J_min, T, T_ab, j_ind, l, N_ab, J_count, ab
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_c, l_c, j_c
        integer(i8) :: N_count

        a_max = (N_max + 1)*(N_max + 2)/2
        l_max = N_max
        J_max = N_2max + 1

        ! Dense lookup Map(ab)%m(TT, J-J_min+1, c) for every b <= a ... J_count = 0 if N_ab > N_2max ...
        allocate(Orb_NNN%Map(a_max*(a_max + 1)/2))
        do a = 1, a_max
            do b = 1, a
                N_ab = 2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l
                J_count = 0
                if (N_ab <= N_2max) J_count = (Orb(a)%j + Orb(b)%j)/2 - abs(Orb(a)%j - Orb(b)%j)/2 + 1
                allocate(Orb_NNN%Map(a*(a - 1)/2 + b)%m(3,J_count,a_max))
                Orb_NNN%Map(a*(a - 1)/2 + b)%m = 0
            end do
        end do

        ! Orbital index of (n,l,j) ...
        allocate(Orb_NNN%nlj(0:N_max/2,0:l_max,2))
        Orb_NNN%nlj = 0
        do a = 1, a_max
            Orb_NNN%nlj(Orb(a)%n,Orb(a)%l,V3b_jindex(Orb(a)%l,Orb(a)%j)) = a
        end do

        allocate(Orb_NNN%N(2,0:l_max,0:J_max,2,2))
        Orb_NNN%N = 0

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
                        ab = a*(a - 1)/2 + b
                        J_min = abs(j_a - j_b)/2
                        do c = 1, a_max
                            n_c = Orb(c)%n
                            l_c = Orb(c)%l
                            j_c = Orb(c)%j
                            P = mod(l_a + l_b + l_c,2) + 1
                            if ((2*(n_a + n_b + n_c) + l_a + l_b + l_c) <= N_3max .and. &
                                (2*(n_a + n_c) + l_a + l_c) <= N_2max .and. &
                                (2*(n_b + n_c) + l_b + l_c) <= N_2max) then
                                j_ind = V3b_jindex(l_c,j_c)
                                do J = abs(j_a - j_b)/2, (j_a + j_b)/2
                                    do T = 1, 2
                                        do T_ab = T - 1, 1
                                            Orb_NNN%N(P,l_c,J,T,j_ind) = Orb_NNN%N(P,l_c,J,T,j_ind) + 1
                                            Orb_NNN%Map(ab)%m(T_ab + T,J - J_min + 1,c) = int(Orb_NNN%N(P,l_c,J,T,j_ind),kind(Orb_NNN%Map(ab)%m))
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
                            allocate(V_NNN%V(P,l,J,T,j_ind)%v(N_count*(N_count + 1_i8)/2_i8))
                            V_NNN%V(P,l,J,T,j_ind)%v = 0.0_sp
                        end do
                    end do
                end do
            end do
        end do
    end subroutine V3b_no2b_initialize

    ! j index within the (l, j = l -+ 1/2) pair ... 1 for j = |2l-1|, 2 for j = 2l+1 ...
    pure integer function V3b_jindex(l_c,j_c)
        integer, intent(in) :: l_c, j_c
        V3b_jindex = merge(1, 2, j_c == abs(2*l_c - 1))
    end function V3b_jindex

    ! Index of the 2+1-body state (ab T_ab c) with total J_ab = J & T = T_ind (1 for 1/2, 2 for 3/2) in its V_NNN block ...
        ! 0 if not stored (b > a, (T_ab,T) = (0,3/2), or outside the space) ... counterpart of V3b_no2b_lookup ...
    pure integer function V3b_no2b_lookup(a,b,c,T_ab,J,T_ind,Orb,Orb_NNN) result(Ind)
        integer, intent(in) :: a, b, c, T_ab, J, T_ind
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb3b_t), intent(in) :: Orb_NNN
        integer :: ab, J_min
        Ind = 0
        ab = a*(a - 1)/2 + b
        if (b > a .or. b < 1 .or. ab > size(Orb_NNN%Map) .or. T_ab < T_ind - 1 .or. T_ab > 1) return
        J_min = abs(Orb(a)%j - Orb(b)%j)/2
        if (J < J_min .or. (J - J_min + 1) > size(Orb_NNN%Map(ab)%m,2) .or. c < 1 .or. c > size(Orb_NNN%Map(ab)%m,3)) return
        Ind = Orb_NNN%Map(ab)%m(T_ab + T_ind,J - J_min + 1,c)
    end function V3b_no2b_lookup

    ! Packed index of (ab T_ab c|V|de T_de f) in the block V(P,l_c,J,(T+1)/2,j_index) ... a >= b & d >= e ...
        ! The dense lookups are plain array reads (no type-bound procedures on class objects in the hot threaded loops: with
        ! ifx they serialized the threads in the old hash map version, see commit 33d40d7) ... a missing element (Julia
        ! KeyError) stops in the out-of-line V3b_no2b_missing_key ...
    integer(i8) function V3b_no2b_index(a,b,c,T_ab,d,e,f,T_de,J,T,P,Orb,Orb_NNN) result(Ind)
        integer, intent(in) :: a, b, c, T_ab, d, e, f, T_de, J, T, P
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb3b_t), intent(in) :: Orb_NNN
        integer :: T_ind, Bra, Ket, Tmp, l_c, j_c, n_f, f_c, P_abc, P_def
        T_ind = (T + 1)/2
        Bra = V3b_no2b_lookup(a,b,c,T_ab,J,T_ind,Orb,Orb_NNN)
        if (Bra == 0) call V3b_no2b_missing_key(a,b,c,T_ab,d,e,f,T_de,J,T,P)
        l_c = Orb(c)%l
        j_c = Orb(c)%j
        ! The Ket orbital is (n_f,l_c,j_c) ... f itself unless (l_f,j_f) /= (l_c,j_c) ...
        if (Orb(f)%l == l_c .and. Orb(f)%j == j_c) then
            f_c = f
        else
            n_f = Orb(f)%n
            f_c = 0
            if (n_f <= ubound(Orb_NNN%nlj,1) .and. l_c <= ubound(Orb_NNN%nlj,2)) f_c = Orb_NNN%nlj(n_f,l_c,V3b_jindex(l_c,j_c))
        end if
        Ket = 0
        if (f_c > 0) Ket = V3b_no2b_lookup(d,e,f_c,T_de,J,T_ind,Orb,Orb_NNN)
        P_abc = mod(Orb(a)%l + Orb(b)%l + l_c,2) + 1
        P_def = mod(Orb(d)%l + Orb(e)%l + l_c,2) + 1
        if (Ket == 0 .or. P /= P_abc .or. P /= P_def) call V3b_no2b_missing_key(a,b,c,T_ab,d,e,f,T_de,J,T,P)
        if (Bra < Ket) then
            Tmp = Bra
            Bra = Ket
            Ket = Tmp
        end if
        ! Packed lower triangle stored row by row (Bra >= Ket) ...
        Ind = int(Ket,i8) + (int(Bra,i8) * int(Bra - 1,i8))/2_i8
    end function V3b_no2b_index

    ! Cold error path of V3b_no2b_index ... kept out of line so the lookup stays small ...
    subroutine V3b_no2b_missing_key(a,b,c,T_ab,d,e,f,T_de,J,T,P)
        integer, intent(in) :: a, b, c, T_ab, d, e, f, T_de, J, T, P
        write(*,'(a,11(1x,I0))') 'ERROR: missing NNN matrix element (KeyError) ... a,b,c,T_ab,d,e,f,T_de,J,T,P =', &
                                 a, b, c, T_ab, d, e, f, T_de, J, T, P
        error stop 1
    end subroutine V3b_no2b_missing_key

    ! NO2B matrix element <ab(J_de T_ab) c | V | de(J_de T_de) f> with total isospin T (doubled) ...
    real(dp) function V3b_no2b(a,b,c,T_ab,d,e,f,T_de,J_de,T,P,V_NNN,Orb,Orb_NNN)
        integer, intent(in) :: a, b, c, T_ab, d, e, f, T_de, J_de, T, P
        type(v3b_t), intent(in) :: V_NNN
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb3b_t), intent(in) :: Orb_NNN
        integer(i8) :: Ind
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
        V3b_no2b = Ph * real(V_NNN%V(P,Orb(c)%l,J_de,(T+1)/2,V3b_jindex(Orb(c)%l,Orb(c)%j))%v(Ind),dp)
    end function V3b_no2b

    ! (a,b) chunks of the NO2B file of the space (N_max,N_2max) ... counterpart of V3b_no2b_allocate ...
    subroutine V3b_no2b_chunks(N_max,N_2max,Orb,ab,ab_max)
        integer, intent(in) :: N_max, N_2max
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable, intent(out) :: ab(:,:)
        integer, intent(out) :: ab_max
        integer :: a_max, a, b, Pass
        a_max = (N_max + 1)*(N_max + 2)/2
        allocate(ab(0,2))
        do Pass = 1, 2
            ab_max = 0
            do a = 1, a_max
                if (2*Orb(a)%n + Orb(a)%l > N_max) cycle
                do b = 1, a
                    if (2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l <= N_2max) then
                        ab_max = ab_max + 1
                        if (Pass == 2) ab(ab_max,:) = [a, b]
                    end if
                end do
            end do
            if (Pass == 1) then
                deallocate(ab)
                allocate(ab(ab_max,2))
            end if
        end do
    end subroutine V3b_no2b_chunks

    ! Orbitals of each (l,j) class in ascending order ... f_list(1:f_count(k),k), k = 2l + j_index - 1 ... the f loops of the
        ! file layout only reach orbitals with (l_f,j_f) = (l_c,j_c), so running over the class keeps the element order ...
    subroutine V3b_lj_classes(a_max,Orb,f_count,f_list)
        integer, intent(in) :: a_max
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable, intent(out) :: f_count(:), f_list(:,:)
        integer :: f, k
        allocate(f_count(0:2*maxval(Orb(1:a_max)%l) + 1))
        f_count = 0
        do f = 1, a_max
            k = 2*Orb(f)%l + V3b_jindex(Orb(f)%l,Orb(f)%j) - 1
            f_count(k) = f_count(k) + 1
        end do
        allocate(f_list(maxval(f_count),0:ubound(f_count,1)))
        f_count = 0
        do f = 1, a_max
            k = 2*Orb(f)%l + V3b_jindex(Orb(f)%l,Orb(f)%j) - 1
            f_count(k) = f_count(k) + 1
            f_list(f_count(k),k) = f
        end do
    end subroutine V3b_lj_classes

    ! Number of Float32 matrix elements of the (a,b) chunk of the NO2B file ... counterpart of V3b_no2b_count (per chunk) ...
    integer(i8) function V3b_no2b_chunk_count(a,b,N_max,N_2max,N_3max,Orb,f_count,f_list) result(N)
        integer, intent(in) :: a, b, N_max, N_2max, N_3max
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: f_count(0:), f_list(:,0:)
        integer :: a_max, c, d, e, f, k, f_i, NN_a, NN_b, NN_c, NN_d, NN_e, NN_f, P_abc, P_def, N_J
        a_max = (N_max + 1)*(N_max + 2)/2
        NN_a = 2*Orb(a)%n + Orb(a)%l
        NN_b = 2*Orb(b)%n + Orb(b)%l
        N = 0
        do c = 1, a_max
            NN_c = 2*Orb(c)%n + Orb(c)%l
            if (.not. ((NN_a + NN_c <= N_2max) .and. (NN_b + NN_c <= N_2max) .and. (NN_a + NN_b + NN_c <= N_3max))) cycle
            P_abc = mod(Orb(a)%l + Orb(b)%l + Orb(c)%l,2) + 1
            k = 2*Orb(c)%l + V3b_jindex(Orb(c)%l,Orb(c)%j) - 1
            do d = 1, a
                NN_d = 2*Orb(d)%n + Orb(d)%l
                do e = 1, d
                    NN_e = 2*Orb(e)%n + Orb(e)%l
                    if (NN_d + NN_e > N_2max) cycle
                    P_def = mod(Orb(d)%l + Orb(e)%l + Orb(c)%l,2) + 1
                    if (P_abc /= P_def) cycle
                    ! (T_ab,T_de,T) combinations per J: (0,0,1/2),(0,1,1/2),(1,0,1/2),(1,1,1/2),(1,1,3/2) ...
                    N_J = 5*max(0, min(Orb(a)%j + Orb(b)%j, Orb(d)%j + Orb(e)%j)/2 &
                                   - max(abs(Orb(a)%j - Orb(b)%j), abs(Orb(d)%j - Orb(e)%j))/2 + 1)
                    do f_i = 1, f_count(k)
                        f = f_list(f_i,k)
                        if (f > a_max) exit
                        NN_f = 2*Orb(f)%n + Orb(f)%l
                        if ((NN_d + NN_f <= N_2max) .and. (NN_e + NN_f <= N_2max) .and. (NN_d + NN_e + NN_f <= N_3max)) N = N + N_J
                    end do
                end do
            end do
        end do
    end function V3b_no2b_chunk_count

    ! Reads the NuHamil NO2B 3-body file (Float32 stream) ... counterpart of V3b_no2b_read & V3b_no2b_read_chunk! ...
        ! The interaction-space orbitals are made here (Int.Nmax) ... the (a,b) chunks needed in the calculation space are
        ! read in parallel, largest first, each thread streaming its chunks through its own buffer; the (c,d,e,f) blocks outside
        ! the calculation space are skipped without reading ...
        ! Port extensions: NNN_File = '' gives a zero NNN interaction, & the file size is checked. Only the chunk offsets up
        ! to the last chunk with N_a <= Calc.Nmax are counted (orbitals are ordered by N) & the c loop of a chunk stops at
        ! the last calculation-space orbital (nothing after it is stored) ... the stored values are those of dev ...
    subroutine V3b_no2b_read(Params,V_NNN,Orb_NNN)
        type(params_t), intent(in) :: Params
        type(v3b_t), intent(out) :: V_NNN
        type(orb3b_t), intent(out) :: Orb_NNN
        type(orb1b_t), allocatable :: Orb(:)
        real(sp), allocatable :: Buffer(:)
        real(dp), allocatable :: Sizes(:)
        integer(i8), allocatable :: N_Chunk(:), N_Chunk_skip(:)
        integer, allocatable :: ab(:,:), Queue(:), Perm(:), f_count(:), f_list(:,:)
        integer(i8) :: File_size, N_Total, Chunk_i, Buffer_0, N_skip, N_ab_chunk
        integer :: N_max_Int, N_2max_Int, N_3max_Int, a_max, a_max_Calc, N_max_Calc, N_2max_Calc, N_3max_Calc
        integer :: u, a, b, c, d, e, f, f_i, k, J, J_min, J_max, T_ab, T_de, T, P, Buffer_N, ab_i, ab_max, ab_last, N_Queue, q
        integer :: n_a, l_a, j_a, n_b, l_b, j_b, n_c, l_c, j_c, n_d, l_d, j_d, n_e, l_e, j_e, n_f, j_ind
        integer :: NN_a, NN_b, NN_c, NN_d, NN_e, NN_f, NN_ab, NN_ac, NN_bc, NN_abc, NN_de, NN_df, NN_ef, NN_def
        integer :: P_abc, P_def

        N_max_Int = Params%Int%Nmax
        N_2max_Int = Params%Int%N2max
        N_3max_Int = Params%Int%N3max
        a_max = (N_max_Int + 1)*(N_max_Int + 2)/2
        N_max_Calc = Params%Calc%Nmax
        N_2max_Calc = Params%Calc%N2max
        N_3max_Calc = Params%Calc%N3max
        a_max_Calc = (N_max_Calc + 1)*(N_max_Calc + 2)/2

        ! Orbitals of the interaction space ...
        Orb = orbitals_make(Params,Use_Int=.true.)

        write(*,'(/,a)') 'Preparing 3-body NNN interaction array ...'
        call V3b_no2b_initialize(N_max_Calc,N_2max_Calc,N_3max_Calc,Orb,V_NNN,Orb_NNN)
        ! No 3N interaction (NNN_File = '') ... the NO2B elements stay zero ...
        if (len_trim(Params%Int%NNN_File) == 0) then
            write(*,'(a)') TAB // 'NNN_File empty ... no 3-body interaction (NO2B elements zero)'
            return
        end if

        ! (a,b) chunks of the file ... only those up to the last one with N_a <= Calc.Nmax are needed ...
        call V3b_no2b_chunks(N_max_Int,N_2max_Int,Orb,ab,ab_max)
        ab_last = 0
        do ab_i = 1, ab_max
            if (2*Orb(ab(ab_i,1))%n + Orb(ab(ab_i,1))%l <= N_max_Calc) ab_last = ab_i
        end do
        call V3b_lj_classes(a_max,Orb,f_count,f_list)

        write(*,'(/,a)') 'Precounting # of 3-body NNN interaction matrix elements ...'
        allocate(N_Chunk(ab_last), N_Chunk_skip(ab_last))
        !$omp parallel do schedule(dynamic)
        do ab_i = 1, ab_last
            N_Chunk(ab_i) = V3b_no2b_chunk_count(ab(ab_i,1),ab(ab_i,2),N_max_Int,N_2max_Int,N_3max_Int,Orb,f_count,f_list)
        end do
        !$omp end parallel do
        N_Total = 0
        do ab_i = 1, ab_last
            N_Chunk_skip(ab_i) = N_Total
            N_Total = N_Total + N_Chunk(ab_i)
        end do

        open(newunit=u, file=trim(Params%Int%NNN_File), access='stream', form='unformatted', status='old', action='read')
        inquire(unit=u, size=File_size)
        close(u)
        if (ab_last > 0) then
            if ((4*N_Total > File_size) .or. ((ab_last == ab_max) .and. (4*N_Total /= File_size))) &
                call fatal('Size of the NNN interaction file does not match the interaction space ' // trim(Params%Int%NNN_File))
        end if

        ! Queue of the chunks needed in the calculation space ... the largest first to balance the threads ...
        allocate(Sizes(ab_last), Perm(ab_last), Queue(ab_last))
        Sizes = real(N_Chunk,dp)
        call sortperm_rev(Sizes,Perm)
        N_Queue = 0
        do q = 1, ab_last
            ab_i = Perm(q)
            NN_a = 2*Orb(ab(ab_i,1))%n + Orb(ab(ab_i,1))%l
            NN_b = 2*Orb(ab(ab_i,2))%n + Orb(ab(ab_i,2))%l
            if ((NN_a <= N_max_Calc) .and. (NN_b <= N_max_Calc) .and. (NN_a + NN_b <= N_2max_Calc)) then
                N_Queue = N_Queue + 1
                Queue(N_Queue) = ab_i
            end if
        end do

        write(*,'(/,a)') 'Reading 3-body NO2B NNN interaction file. ..'
        !$omp parallel default(shared) private(u,Buffer,Chunk_i,Buffer_0,Buffer_N,N_skip,N_ab_chunk,ab_i,a,b,c,d,e,f,f_i,k, &
        !$omp J,J_min,J_max,T_ab,T_de,T,P,j_ind,n_a,l_a,j_a,n_b,l_b,j_b,n_c,l_c,j_c,n_d,l_d,j_d,n_e,l_e,j_e,n_f, &
        !$omp NN_a,NN_b,NN_c,NN_d,NN_e,NN_f,NN_ab,NN_ac,NN_bc,NN_abc,NN_de,NN_df,NN_ef,NN_def,P_abc,P_def)
        open(newunit=u, file=trim(Params%Int%NNN_File), access='stream', form='unformatted', status='old', action='read')
        allocate(Buffer(Buffer_size))
        !$omp do schedule(dynamic,1)
        do q = 1, N_Queue
            ab_i = Queue(q)
            a = ab(ab_i,1)
            b = ab(ab_i,2)
            N_skip = N_Chunk_skip(ab_i)
            N_ab_chunk = N_Chunk(ab_i)
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            NN_a = 2*n_a + l_a
            n_b = Orb(b)%n
            l_b = Orb(b)%l
            j_b = Orb(b)%j
            NN_b = 2*n_b + l_b
            NN_ab = NN_a + NN_b
            ! Chunk_i ... index of the element within the chunk, the Buffer holds the elements (Buffer_0+1):(Buffer_0+Buffer_N) ...
            Chunk_i = 0
            Buffer_0 = 0
            Buffer_N = 0
            ! Nothing beyond the last calculation-space orbital c is stored ...
            do c = 1, a_max_Calc
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
                j_ind = V3b_jindex(l_c,j_c)
                k = 2*l_c + j_ind - 1
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
                        P_def = mod(l_d + l_e + l_c,2) + 1
                        if (P_abc /= P_def) cycle
                        J_min = max(abs(j_a - j_b),abs(j_d - j_e))/2
                        J_max = min((j_a + j_b),(j_d + j_e))/2
                        do f_i = 1, f_count(k)
                            f = f_list(f_i,k)
                            n_f = Orb(f)%n
                            NN_f = 2*n_f + l_c
                            NN_df = NN_d + NN_f
                            NN_ef = NN_e + NN_f
                            NN_def = NN_d + NN_e + NN_f
                            if (.not. ((NN_df <= N_2max_Int) .and. (NN_ef <= N_2max_Int) .and. (NN_def <= N_3max_Int))) cycle
                            if ((NN_c <= N_max_Calc) .and. (NN_d <= N_max_Calc) .and. &
                                (NN_e <= N_max_Calc) .and. (NN_f <= N_max_Calc) .and. &
                                (NN_ac <= N_2max_Calc) .and. (NN_bc <= N_2max_Calc) .and. &
                                (NN_de <= N_2max_Calc) .and. (NN_df <= N_2max_Calc) .and. (NN_ef <= N_2max_Calc) .and. &
                                (NN_abc <= N_3max_Calc) .and. (NN_def <= N_3max_Calc)) then
                                do J = J_min, J_max
                                    do T_ab = 0, 1
                                        do T_de = 0, 1
                                            do T = max(abs(2*T_ab-1),abs(2*T_de-1)), min((2*T_ab+1),(2*T_de+1)), 2
                                                Chunk_i = Chunk_i + 1
                                                ! Refill the buffer from the file once Chunk_i is out of it ...
                                                if (Chunk_i > Buffer_0 + Buffer_N) then
                                                    Buffer_0 = Chunk_i - 1
                                                    Buffer_N = int(min(int(Buffer_size,i8), N_ab_chunk - Buffer_0))
                                                    read(u, pos=4_i8*(N_skip + Buffer_0) + 1_i8) Buffer(1:Buffer_N)
                                                end if
                                                V_NNN%V(P,l_c,J,(T+1)/2,j_ind)%v(V3b_no2b_index(a,b,c,T_ab,d,e,f,T_de,J,T,P,Orb,Orb_NNN)) = &
                                                    Buffer(Chunk_i - Buffer_0)
                                            end do
                                        end do
                                    end do
                                end do
                            else
                                ! Skip the elements outside the calculation space ... 5 (T_ab,T_de,T) combinations per J ...
                                Chunk_i = Chunk_i + 5*max(J_max - J_min + 1,0)
                            end if
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

    ! Cuts the NO2B file Int.NNN_File (space Int.Nmax/N2max/N3max) down to (N_max_Out,N_2max_Out,N_3max_Out) & writes it in the same
        ! NuHamil Float32 layout to Out_File ... counterpart of V3b_no2b_parse ... chunks are processed in parallel (largest first),
        ! each thread with its own units, writing disjoint byte ranges of the pre-sized output file ...
    subroutine V3b_no2b_parse(Params,N_max_Out,N_2max_Out,N_3max_Out,Out_File)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_max_Out, N_2max_Out, N_3max_Out
        character(len=*), intent(in) :: Out_File
        type(orb1b_t), allocatable :: Orb(:)
        real(sp), allocatable :: Buffer_In(:), Buffer_Out(:)
        real(dp), allocatable :: Sizes(:)
        integer(i8), allocatable :: N_Chunk_In(:), N_Chunk_Out(:), N_Skip_In(:), N_Skip_Out(:)
        integer, allocatable :: ab_In(:,:), ab_Out(:,:), ab_In_index(:), Perm(:), f_count(:), f_list(:,:)
        integer(i8) :: Size_In, Size_Out, File_size, Chunk_i, Buffer_0, Out_N, N_In, N_Out, Sum_In, Sum_Out, Pos_Out
        integer :: N_max_In, N_2max_In, N_3max_In, a_max, ab_max_In, ab_max_Out, ab_i, ab_o, q, u_in, u_out, Buffer_N, Out_i
        integer :: a, b, c, d, e, f, f_i, k, J_min, J_max, i_J, N_J, Err
        integer :: NN_a, NN_b, NN_c, NN_d, NN_e, NN_f, NN_ab, NN_ac, NN_bc, NN_abc, NN_de, NN_df, NN_ef, NN_def, P_abc, P_def
        logical :: Failed
        character(len=:), allocatable :: In_File

        In_File = trim(Params%Int%NNN_File)
        N_max_In = Params%Int%Nmax
        N_2max_In = Params%Int%N2max
        N_3max_In = Params%Int%N3max
        a_max = (N_max_In + 1)*(N_max_In + 2)/2

        ! The output space has to be contained in the interaction space ...
        if ((N_max_Out > N_max_In) .or. (N_2max_Out > N_2max_In) .or. (N_3max_Out > N_3max_In)) &
            call fatal('Output space (' // istr(N_max_Out) // ',' // istr(N_2max_Out) // ',' // istr(N_3max_Out) // &
                       ') is not contained in the interaction space (' // istr(N_max_In) // ',' // istr(N_2max_In) // ',' // &
                       istr(N_3max_In) // ') ...')
        if (min(N_max_Out,N_2max_Out,N_3max_Out) < 0) call fatal('V3b_no2b_parse: negative output truncation ...')
        if (len_trim(In_File) == 0) call fatal('V3b_no2b_parse: NNN_File is empty ...')
        if (abspath_of(In_File) == abspath_of(Out_File)) call fatal('Output file has to differ from the interaction file ...')

        Orb = orbitals_make(Params,Use_Int=.true.)
        call V3b_lj_classes(a_max,Orb,f_count,f_list)

        write(*,'(/,a)') 'Precounting # of 3-body NNN interaction matrix elements ...'
        call V3b_no2b_chunks(N_max_In,N_2max_In,Orb,ab_In,ab_max_In)
        call V3b_no2b_chunks(N_max_Out,N_2max_Out,Orb,ab_Out,ab_max_Out)
        allocate(N_Chunk_In(ab_max_In), N_Skip_In(ab_max_In), N_Chunk_Out(ab_max_Out), N_Skip_Out(ab_max_Out))
        !$omp parallel do schedule(dynamic)
        do ab_i = 1, ab_max_In
            N_Chunk_In(ab_i) = V3b_no2b_chunk_count(ab_In(ab_i,1),ab_In(ab_i,2),N_max_In,N_2max_In,N_3max_In,Orb,f_count,f_list)
        end do
        !$omp end parallel do
        !$omp parallel do schedule(dynamic)
        do ab_i = 1, ab_max_Out
            N_Chunk_Out(ab_i) = V3b_no2b_chunk_count(ab_Out(ab_i,1),ab_Out(ab_i,2),N_max_Out,N_2max_Out,N_3max_Out,Orb,f_count,f_list)
        end do
        !$omp end parallel do
        Sum_In = 0
        do ab_i = 1, ab_max_In
            N_Skip_In(ab_i) = 4*Sum_In
            Sum_In = Sum_In + N_Chunk_In(ab_i)
        end do
        Sum_Out = 0
        do ab_i = 1, ab_max_Out
            N_Skip_Out(ab_i) = 4*Sum_Out
            Sum_Out = Sum_Out + N_Chunk_Out(ab_i)
        end do
        Size_In = 4*Sum_In
        Size_Out = 4*Sum_Out

        ! Size of the interaction file ...
        open(newunit=u_in, file=In_File, access='stream', form='unformatted', status='old', action='read', iostat=Err)
        if (Err /= 0) call fatal('Cannot open the NNN interaction file ' // In_File)
        inquire(unit=u_in, size=File_size)
        close(u_in)
        if (File_size /= Size_In) call fatal('Size of ' // In_File // ' (' // istr(File_size) // &
                                             ' B) does not match the interaction space (' // istr(Size_In) // ' B) ...')

        ! Input chunk of each output chunk ... the orbitals of the smaller space are the first ones, so the chunks follow in
            ! the same order ...
        allocate(ab_In_index(ab_max_Out))
        ab_i = 1
        do ab_o = 1, ab_max_Out
            do while ((ab_In(ab_i,1) /= ab_Out(ab_o,1)) .or. (ab_In(ab_i,2) /= ab_Out(ab_o,2)))
                ab_i = ab_i + 1
            end do
            ab_In_index(ab_o) = ab_i
        end do

        call presize_file(Out_File,Size_Out)

        allocate(Sizes(ab_max_Out), Perm(ab_max_Out))
        Sizes = real(N_Chunk_Out,dp)
        call sortperm_rev(Sizes,Perm)

        write(*,'(/,a)') 'Exporting 3-body NO2B NNN interaction file into N_max = ' // istr(N_max_Out) // ' , N_2max = ' // &
                         istr(N_2max_Out) // ' , N_3max = ' // istr(N_3max_Out) // ' ...'
        Failed = .false.
        !$omp parallel default(shared) private(u_in,u_out,Buffer_In,Buffer_Out,q,ab_o,ab_i,a,b,c,d,e,f,f_i,k,J_min,J_max,i_J, &
        !$omp N_J,NN_a,NN_b,NN_c,NN_d,NN_e,NN_f,NN_ab,NN_ac,NN_bc,NN_abc,NN_de,NN_df,NN_ef,NN_def,P_abc,P_def, &
        !$omp Chunk_i,Buffer_0,Buffer_N,Out_i,Out_N,N_In,N_Out,Pos_Out)
        ! (Params%Int%NNN_File, not the deferred-length In_File: gfortran passes a garbage name in the parallel region) ...
        open(newunit=u_in, file=trim(Params%Int%NNN_File), access='stream', form='unformatted', status='old', action='read')
        open(newunit=u_out, file=trim(Out_File), access='stream', form='unformatted', status='old', action='readwrite')
        allocate(Buffer_In(Buffer_size), Buffer_Out(Buffer_size))
        !$omp do schedule(dynamic,1)
        do q = 1, ab_max_Out
            ab_o = Perm(q)
            ab_i = ab_In_index(ab_o)
            a = ab_Out(ab_o,1)
            b = ab_Out(ab_o,2)
            N_In = N_Chunk_In(ab_i)
            N_Out = N_Chunk_Out(ab_o)
            NN_a = 2*Orb(a)%n + Orb(a)%l
            NN_b = 2*Orb(b)%n + Orb(b)%l
            NN_ab = NN_a + NN_b
            Chunk_i = 0
            Buffer_0 = 0
            Buffer_N = 0
            Out_i = 0
            Out_N = 0
            Pos_Out = N_Skip_Out(ab_o) + 1
            do c = 1, a_max
                NN_c = 2*Orb(c)%n + Orb(c)%l
                NN_ac = NN_a + NN_c
                NN_bc = NN_b + NN_c
                NN_abc = NN_a + NN_b + NN_c
                if (.not. ((NN_ac <= N_2max_In) .and. (NN_bc <= N_2max_In) .and. (NN_abc <= N_3max_In))) cycle
                P_abc = mod(Orb(a)%l + Orb(b)%l + Orb(c)%l,2) + 1
                k = 2*Orb(c)%l + V3b_jindex(Orb(c)%l,Orb(c)%j) - 1
                do d = 1, a
                    NN_d = 2*Orb(d)%n + Orb(d)%l
                    do e = 1, d
                        NN_e = 2*Orb(e)%n + Orb(e)%l
                        NN_de = NN_d + NN_e
                        if (NN_de > N_2max_In) cycle
                        P_def = mod(Orb(d)%l + Orb(e)%l + Orb(c)%l,2) + 1
                        if (P_abc /= P_def) cycle
                        J_min = max(abs(Orb(a)%j - Orb(b)%j),abs(Orb(d)%j - Orb(e)%j))/2
                        J_max = min(Orb(a)%j + Orb(b)%j,Orb(d)%j + Orb(e)%j)/2
                        N_J = 5*max(J_max - J_min + 1,0)
                        do f_i = 1, f_count(k)
                            f = f_list(f_i,k)
                            NN_f = 2*Orb(f)%n + Orb(f)%l
                            NN_df = NN_d + NN_f
                            NN_ef = NN_e + NN_f
                            NN_def = NN_d + NN_e + NN_f
                            if (.not. ((NN_df <= N_2max_In) .and. (NN_ef <= N_2max_In) .and. (NN_def <= N_3max_In))) cycle
                            if ((NN_a <= N_max_Out) .and. (NN_b <= N_max_Out) .and. (NN_c <= N_max_Out) .and. &
                                (NN_d <= N_max_Out) .and. (NN_e <= N_max_Out) .and. (NN_f <= N_max_Out) .and. &
                                (NN_ab <= N_2max_Out) .and. (NN_ac <= N_2max_Out) .and. (NN_bc <= N_2max_Out) .and. &
                                (NN_de <= N_2max_Out) .and. (NN_df <= N_2max_Out) .and. (NN_ef <= N_2max_Out) .and. &
                                (NN_abc <= N_3max_Out) .and. (NN_def <= N_3max_Out)) then
                                ! Copy the 5 (T_ab,T_de,T) combinations per J ...
                                do i_J = 1, N_J
                                    Chunk_i = Chunk_i + 1
                                    if (Chunk_i > Buffer_0 + Buffer_N) then
                                        Buffer_0 = Chunk_i - 1
                                        Buffer_N = int(min(int(Buffer_size,i8), N_In - Buffer_0))
                                        read(u_in, pos=N_Skip_In(ab_i) + 4_i8*Buffer_0 + 1_i8) Buffer_In(1:Buffer_N)
                                    end if
                                    if (Out_i == Buffer_size) then
                                        write(u_out, pos=Pos_Out) Buffer_Out(1:Out_i)
                                        Pos_Out = Pos_Out + 4_i8*Out_i
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
            end do
            if (Out_i > 0) write(u_out, pos=Pos_Out) Buffer_Out(1:Out_i)
            if ((Chunk_i /= N_In) .or. (Out_N /= N_Out)) then
                write(*,'(a)') 'ERROR: Wrong # of matrix elements in the chunk (a,b) = (' // istr(a) // ',' // istr(b) // ') ...'
                Failed = .true.
            end if
        end do
        !$omp end do
        close(u_in)
        close(u_out)
        deallocate(Buffer_In, Buffer_Out)
        !$omp end parallel
        if (Failed) call fatal('V3b_no2b_parse: wrong # of matrix elements ...')

        write(*,'(/,a)') 'Finished export of 3-body NO2B NNN interaction into ' // trim(Out_File) // ' ...'
    end subroutine V3b_no2b_parse

end module nj_v3b
