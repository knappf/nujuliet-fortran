! 2-body NN operators ... counterpart of src/common/O2b.jl (NuJuliet dev; MatrixElements/O2B.jl before) ...
    ! Each (P,J) channel stores the lower triangle Bra >= Ket packed as Ind = Bra + (Ket-1)*N - Ket*(Ket-1)/2 ...
    ! T = 0 (pn) channels use all ordered pairs (a,b), T = 1 (pp & nn) channels pairs with a >= b ...
module nj_o2b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    implicit none
    private

    public :: O2b_initialize, O2b_index, O2b_pp, O2b_pn, O2b_nn, O2b_export, O2b_export_HRF, O2b_import
    public :: O2b_temp_initialize, O2b_temp_index, O2b_transformation, O2b_Pandya_transformation
    public :: packed_index, Orb_PreComp, O2b_move
    public :: O2b_temp, O2b_Pandya_initialize, O2b_Pandya_pp, O2b_Pandya_pnpn, O2b_Pandya_pnnp, O2b_Pandya_nn

    ! O2b_Pandya_transformation(Params,Orb,Orb_NN,O_NN,F_NN,Orb_NN_t) with type(o2b_pandya_t) :: F_NN (dev, 4 blocks) ...
        ! the 5-argument form with type(o2b_temp_t) :: F_NN is LEGACY (removed in cleanup stage) ...
    interface O2b_Pandya_transformation
        module procedure O2b_Pandya_transformation_temp, O2b_Pandya_transformation_4
    end interface O2b_Pandya_transformation

contains

    pure integer function packed_index(Bra,Ket,N)
        integer, intent(in) :: Bra, Ket, N
        packed_index = Bra + (Ket - 1) * N - (Ket * (Ket - 1)) / 2
    end function packed_index

    ! Allocates a zero 2-body operator ... optionally also the 2-body orbitals Orb_NN ...
    subroutine O2b_initialize(Params,Orb,O_NN,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o2b_t), intent(out) :: O_NN
        type(orb2b_t), intent(out), optional :: Orb_NN
        integer :: N_max, N_2max, a_max, J_max, a, b, P, J, Pass, Length1, Length2
        integer, allocatable :: N_Orb_NN(:,:,:)
        type(orb2b_t) :: Orb_loc

        N_max = Params%Calc%Nmax
        N_2max = Params%Calc%N2max
        a_max = (N_max + 1)*(N_max + 2)/2
        J_max = N_2max + 1

        allocate(N_Orb_NN(2,2,0:J_max))
        Orb_loc%a_max = a_max
        Orb_loc%J_max = J_max
        allocate(Orb_loc%Ind(2,2,0:J_max), Orb_loc%Dic(a_max,a_max,0:J_max,0:1))
        Orb_loc%Dic = 0

        ! Pass 1 counts the 2-body orbitals, pass 2 fills the indices ...
        do Pass = 1, 2
            N_Orb_NN = 0
            do a = 1, a_max
                if ((2*Orb(a)%n + Orb(a)%l) <= N_max) then
                    do b = 1, a_max
                        if ((2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l) <= N_2max) then
                            P = mod(Orb(a)%l + Orb(b)%l, 2) + 1
                            do J = abs(Orb(a)%j - Orb(b)%j)/2, (Orb(a)%j + Orb(b)%j)/2
                                if (b <= a) then
                                    N_Orb_NN(2,P,J) = N_Orb_NN(2,P,J) + 1
                                    if (Pass == 2) then
                                        Orb_loc%Dic(a,b,J,1) = N_Orb_NN(2,P,J)
                                        Orb_loc%Ind(2,P,J)%v(1,N_Orb_NN(2,P,J)) = a
                                        Orb_loc%Ind(2,P,J)%v(2,N_Orb_NN(2,P,J)) = b
                                    end if
                                end if
                                N_Orb_NN(1,P,J) = N_Orb_NN(1,P,J) + 1
                                if (Pass == 2) then
                                    Orb_loc%Dic(a,b,J,0) = N_Orb_NN(1,P,J)
                                    Orb_loc%Ind(1,P,J)%v(1,N_Orb_NN(1,P,J)) = a
                                    Orb_loc%Ind(1,P,J)%v(2,N_Orb_NN(1,P,J)) = b
                                end if
                            end do
                        end if
                    end do
                end if
            end do
            if (Pass == 1) then
                do P = 1, 2
                    do J = 0, J_max
                        allocate(Orb_loc%Ind(1,P,J)%v(2,N_Orb_NN(1,P,J)))
                        allocate(Orb_loc%Ind(2,P,J)%v(2,N_Orb_NN(2,P,J)))
                    end do
                end do
            end if
        end do
        Orb_loc%N = N_Orb_NN

        ! Zero entries of O ...
        allocate(O_NN%pp(2,0:J_max), O_NN%pn(2,0:J_max), O_NN%nn(2,0:J_max))
        do P = 1, 2
            do J = 0, J_max
                Length1 = N_Orb_NN(1,P,J)*(N_Orb_NN(1,P,J) + 1)/2
                Length2 = N_Orb_NN(2,P,J)*(N_Orb_NN(2,P,J) + 1)/2
                allocate(O_NN%pp(P,J)%v(Length2), O_NN%pn(P,J)%v(Length1), O_NN%nn(P,J)%v(Length2))
                O_NN%pp(P,J)%v = 0.0_dp
                O_NN%pn(P,J)%v = 0.0_dp
                O_NN%nn(P,J)%v = 0.0_dp
            end do
        end do

        if (present(Orb_NN)) call move_orb2b(Orb_loc,Orb_NN)
    end subroutine O2b_initialize

    ! Moves the blocks of Src into Dst without copying (Src ends deallocated) ...
    subroutine O2b_move(Src,Dst)
        type(o2b_t), intent(inout) :: Src
        type(o2b_t), intent(inout) :: Dst
        call move_alloc(Src%pp, Dst%pp)
        call move_alloc(Src%pn, Dst%pn)
        call move_alloc(Src%nn, Dst%nn)
    end subroutine O2b_move

    subroutine move_orb2b(Src,Dst)
        type(orb2b_t), intent(inout) :: Src
        type(orb2b_t), intent(out) :: Dst
        Dst%a_max = Src%a_max
        Dst%J_max = Src%J_max
        call move_alloc(Src%N, Dst%N)
        call move_alloc(Src%Ind, Dst%Ind)
        call move_alloc(Src%Dic, Dst%Dic)
    end subroutine move_orb2b

    ! Packed index of <ab|O|cd> in channel (J,P,T) ...
    pure integer function O2b_index(a,b,c,d,J,P,T,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P, T
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: N, Bra, Ket
        N = Orb_NN%N(T+1,P,J)
        Bra = Orb_NN%Dic(a,b,J,T)
        Ket = Orb_NN%Dic(c,d,J,T)
        O2b_index = max(Bra,Ket) + (min(Bra,Ket) - 1) * N - (min(Bra,Ket) * (min(Bra,Ket) - 1)) / 2
    end function O2b_index

    ! Antisymmetrized pp matrix element with exchange phases for a < b or c < d ...
    pure real(dp) function O2b_pp(a,b,c,d,J,P,O,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_t), intent(in) :: O
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call o2b_t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        O2b_pp = Amp * O%pp(P,J)%v(O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN))
    end function O2b_pp

    pure real(dp) function O2b_pn(a,b,c,d,J,P,O,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_t), intent(in) :: O
        type(orb2b_t), intent(in) :: Orb_NN
        O2b_pn = O%pn(P,J)%v(O2b_index(a,b,c,d,J,P,0,Orb_NN))
    end function O2b_pn

    pure real(dp) function O2b_nn(a,b,c,d,J,P,O,Orb,Orb_NN)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_t), intent(in) :: O
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        integer :: a_s, b_s, c_s, d_s
        real(dp) :: Amp
        call o2b_t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        O2b_nn = Amp * O%nn(P,J)%v(O2b_index(a_s,b_s,c_s,d_s,J,P,1,Orb_NN))
    end function O2b_nn

    pure subroutine o2b_t1_order(a,b,c,d,J,Orb,a_s,b_s,c_s,d_s,Amp)
        integer, intent(in) :: a, b, c, d, J
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(out) :: a_s, b_s, c_s, d_s
        real(dp), intent(out) :: Amp
        real(dp) :: Sign_ab, Sign_cd
        if (a < b) then
            a_s = b
            b_s = a
            Sign_ab = phase(J + 1 + (Orb(a)%j + Orb(b)%j)/2)
        else
            a_s = a
            b_s = b
            Sign_ab = 1.0_dp
        end if
        if (c < d) then
            c_s = d
            d_s = c
            Sign_cd = phase(J + 1 + (Orb(c)%j + Orb(d)%j)/2)
        else
            c_s = c
            d_s = d
            Sign_cd = 1.0_dp
        end if
        Amp = Sign_ab * Sign_cd
        if ((a < b) .and. (c < d)) Amp = phase((Orb(a)%j + Orb(b)%j + Orb(c)%j + Orb(d)%j)/2)
    end subroutine o2b_t1_order

    ! Highest J channel written by O2b_export & read by O2b_import ... 2*N_max + 1 as in the original ...
        ! limited to the allocated channels so that N2max < 2*Nmax does not run out of bounds ...
    pure integer function export_J_max(Params,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb2b_t), intent(in) :: Orb_NN
        export_J_max = min(2*Params%Calc%Nmax + 1, Orb_NN%J_max)
    end function export_J_max

    ! Binary export ... pp (lower triangle), pn (full square loop of the packed formula, as in the original), nn ...
    subroutine O2b_export(Params,Orb_NN,O_NN,Export_Path)
        type(params_t), intent(in) :: Params
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: O_NN
        character(len=*), intent(in) :: Export_Path
        integer :: u, J, P, Bra, Ket, N
        real(dp), allocatable :: Buf(:)

        write(*,'(/,a)') 'Exporting given 2-body operator in internal binary format ...'
        open(newunit=u, file=Export_Path, access='stream', form='unformatted', status='replace', action='write')
        do J = 0, export_J_max(Params,Orb_NN)
            do P = 1, 2
                N = Orb_NN%N(2,P,J)
                allocate(Buf(N*(N+1)/2))
                do Bra = 1, N
                    do Ket = 1, Bra
                        Buf(Bra*(Bra-1)/2 + Ket) = O_NN%pp(P,J)%v(packed_index(Bra,Ket,N))
                    end do
                end do
                write(u) Buf
                deallocate(Buf)
            end do
        end do
        do J = 0, export_J_max(Params,Orb_NN)
            do P = 1, 2
                N = Orb_NN%N(1,P,J)
                allocate(Buf(N*N))
                do Bra = 1, N
                    do Ket = 1, N
                        Buf((Bra-1)*N + Ket) = O_NN%pn(P,J)%v(packed_index(max(Bra,Ket),min(Bra,Ket),N))   ! full symmetric pn matrix
                    end do
                end do
                write(u) Buf
                deallocate(Buf)
            end do
        end do
        do J = 0, export_J_max(Params,Orb_NN)
            do P = 1, 2
                N = Orb_NN%N(2,P,J)
                allocate(Buf(N*(N+1)/2))
                do Bra = 1, N
                    do Ket = 1, Bra
                        Buf(Bra*(Bra-1)/2 + Ket) = O_NN%nn(P,J)%v(packed_index(Bra,Ket,N))
                    end do
                end do
                write(u) Buf
                deallocate(Buf)
            end do
        end do
        close(u)
        write(*,'(/,a)') "Given 2-body NN operator succesfully exported in the binary format ... ''" // Export_Path // "''"
    end subroutine O2b_export

    subroutine O2b_export_HRF(Params,Orb_NN,O_NN,Export_Path)
        type(params_t), intent(in) :: Params
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: O_NN
        character(len=*), intent(in) :: Export_Path
        integer :: u, J, P, Bra, Ket, N, a, b, c, d

        write(*,'(/,a)') 'Exporting given 2-body operator in the Human-Readable-Format ...'
        open(newunit=u, file=Export_Path, status='replace', action='write')
        write(u,'(a)') 't' // TAB // 'a' // TAB // 'b' // TAB // 'c' // TAB // 'd' // TAB // 'J' // TAB // 'V'
        do J = 0, export_J_max(Params,Orb_NN)
            do P = 1, 2
                N = Orb_NN%N(2,P,J)
                do Bra = 1, N
                    a = 2*Orb_NN%Ind(2,P,J)%v(1,Bra) - 1
                    b = 2*Orb_NN%Ind(2,P,J)%v(2,Bra) - 1
                    do Ket = 1, Bra
                        c = 2*Orb_NN%Ind(2,P,J)%v(1,Ket) - 1
                        d = 2*Orb_NN%Ind(2,P,J)%v(2,Ket) - 1
                        call row('-1', O_NN%pp(P,J)%v(packed_index(Bra,Ket,N)))
                    end do
                end do
            end do
        end do
        do J = 0, export_J_max(Params,Orb_NN)
            do P = 1, 2
                N = Orb_NN%N(1,P,J)
                do Bra = 1, N
                    a = 2*Orb_NN%Ind(1,P,J)%v(1,Bra) - 1
                    b = 2*Orb_NN%Ind(1,P,J)%v(2,Bra)
                    do Ket = 1, N
                        c = 2*Orb_NN%Ind(1,P,J)%v(1,Ket) - 1
                        d = 2*Orb_NN%Ind(1,P,J)%v(2,Ket)
                        call row('0', O_NN%pn(P,J)%v(packed_index(max(Bra,Ket),min(Bra,Ket),N)))
                    end do
                end do
            end do
        end do
        do J = 0, export_J_max(Params,Orb_NN)
            do P = 1, 2
                N = Orb_NN%N(2,P,J)
                do Bra = 1, N
                    a = 2*Orb_NN%Ind(2,P,J)%v(1,Bra)
                    b = 2*Orb_NN%Ind(2,P,J)%v(2,Bra)
                    do Ket = 1, Bra
                        c = 2*Orb_NN%Ind(2,P,J)%v(1,Ket)
                        d = 2*Orb_NN%Ind(2,P,J)%v(2,Ket)
                        call row('1', O_NN%nn(P,J)%v(packed_index(Bra,Ket,N)))
                    end do
                end do
            end do
        end do
        close(u)
        write(*,'(/,a)') "Given 2-body NN operator succesfully exported in the Human-Readable-Format ... ''" // Export_Path // "''"
    contains
        subroutine row(t,ME)
            character(len=*), intent(in) :: t
            real(dp), intent(in) :: ME
            write(u,'(a)') t // TAB // istr(a) // TAB // istr(b) // TAB // istr(c) // TAB // istr(d) // TAB // istr(J) // &
                           TAB // jstr(anint(ME * 1.0e10_dp) / 1.0e10_dp)
        end subroutine row
    end subroutine O2b_export_HRF

    subroutine O2b_import(Params,Orb,Import_Path,O_NN,Orb_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: Import_Path
        type(o2b_t), intent(out) :: O_NN
        type(orb2b_t), intent(out), optional :: Orb_NN
        type(orb2b_t) :: Orb_loc
        integer :: u, J, P, Bra, Ket, N, Ios
        character(len=256) :: Msg
        real(dp), allocatable :: Buf(:)

        write(*,'(/,a)') 'Reading given 2-body NN operator ...'
        call O2b_initialize(Params,Orb,O_NN,Orb_loc)
        open(newunit=u, file=Import_Path, access='stream', form='unformatted', status='old', action='read', &
             iostat=Ios, iomsg=Msg)
        if (Ios /= 0) call fatal('O2b_import: cannot open ' // trim(Import_Path) // ' ... ' // trim(Msg))
        do J = 0, export_J_max(Params,Orb_loc)
            do P = 1, 2
                N = Orb_loc%N(2,P,J)
                allocate(Buf(N*(N+1)/2))
                read(u) Buf
                do Bra = 1, N
                    do Ket = 1, Bra
                        O_NN%pp(P,J)%v(packed_index(Bra,Ket,N)) = Buf(Bra*(Bra-1)/2 + Ket)
                    end do
                end do
                deallocate(Buf)
            end do
        end do
        do J = 0, export_J_max(Params,Orb_loc)
            do P = 1, 2
                N = Orb_loc%N(1,P,J)
                allocate(Buf(N*N))
                read(u) Buf
                ! lower triangle only ... the upper one of files written before the packed-index fix holds other elements
                do Bra = 1, N
                    do Ket = 1, Bra
                        O_NN%pn(P,J)%v(packed_index(Bra,Ket,N)) = Buf((Bra-1)*N + Ket)
                    end do
                end do
                deallocate(Buf)
            end do
        end do
        do J = 0, export_J_max(Params,Orb_loc)
            do P = 1, 2
                N = Orb_loc%N(2,P,J)
                allocate(Buf(N*(N+1)/2))
                read(u) Buf
                do Bra = 1, N
                    do Ket = 1, Bra
                        O_NN%nn(P,J)%v(packed_index(Bra,Ket,N)) = Buf(Bra*(Bra-1)/2 + Ket)
                    end do
                end do
                deallocate(Buf)
            end do
        end do
        close(u)
        write(*,'(/,a)') 'Given 2-body NN operator succesfully read ...'
        if (present(Orb_NN)) call move_orb2b(Orb_loc,Orb_NN)
    end subroutine O2b_import

    ! Temporary full-matrix 2-body operator over all ordered pairs (a,b) ...
    subroutine O2b_temp_initialize(Params,Orb,O_NN_t,Orb_NN_t)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o2b_temp_t), intent(out) :: O_NN_t
        type(orb2b_temp_t), intent(out), optional :: Orb_NN_t
        integer :: J_max, P, J, Length
        type(orb2b_temp_t) :: Orb_loc

        call orb2b_temp_make(Params,Orb,Orb_loc)
        J_max = Params%Calc%N2max + 1

        allocate(O_NN_t%pp(2,0:J_max), O_NN_t%pn(2,0:J_max), O_NN_t%nn(2,0:J_max))
        do P = 1, 2
            do J = 0, J_max
                Length = Orb_loc%N(P,J)
                allocate(O_NN_t%pp(P,J)%v(Length,Length), O_NN_t%pn(P,J)%v(Length,Length), O_NN_t%nn(P,J)%v(Length,Length))
                O_NN_t%pp(P,J)%v = 0.0_dp
                O_NN_t%pn(P,J)%v = 0.0_dp
                O_NN_t%nn(P,J)%v = 0.0_dp
            end do
        end do

        if (present(Orb_NN_t)) then
            call move_alloc(Orb_loc%N, Orb_NN_t%N)
            call move_alloc(Orb_loc%Ind, Orb_NN_t%Ind)
            call move_alloc(Orb_loc%Dic, Orb_NN_t%Dic)
        end if
    end subroutine O2b_temp_initialize

    ! Pairs of the temporary full-matrix operators (O2B_Temp & O2B_Pandya) ... a with 2n_a + l_a <= Nmax, all b with
        ! 2(n_a + n_b) + l_a + l_b <= N2max, J_max = N2max + 1 ... Dic(a,b,J) = Julia Map[a,b,J+1] ...
    subroutine orb2b_temp_make(Params,Orb,Orb_loc)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_temp_t), intent(out) :: Orb_loc
        integer :: N_max, N_2max, a_max, J_max, a, b, P, J, Pass
        integer, allocatable :: N_Orb_NN(:,:)

        N_max = Params%Calc%Nmax
        N_2max = Params%Calc%N2max
        a_max = (N_max + 1)*(N_max + 2)/2
        J_max = N_2max + 1
        allocate(N_Orb_NN(2,0:J_max), Orb_loc%Ind(2,0:J_max), Orb_loc%Dic(a_max,a_max,0:J_max))
        Orb_loc%Dic = 0

        do Pass = 1, 2
            N_Orb_NN = 0
            do a = 1, a_max
                if ((2*Orb(a)%n + Orb(a)%l) <= N_max) then
                    do b = 1, a_max
                        if ((2*(Orb(a)%n + Orb(b)%n) + Orb(a)%l + Orb(b)%l) <= N_2max) then
                            P = mod(Orb(a)%l + Orb(b)%l, 2) + 1
                            do J = abs(Orb(a)%j - Orb(b)%j)/2, (Orb(a)%j + Orb(b)%j)/2
                                N_Orb_NN(P,J) = N_Orb_NN(P,J) + 1
                                if (Pass == 2) then
                                    Orb_loc%Dic(a,b,J) = N_Orb_NN(P,J)
                                    Orb_loc%Ind(P,J)%v(1,N_Orb_NN(P,J)) = a
                                    Orb_loc%Ind(P,J)%v(2,N_Orb_NN(P,J)) = b
                                end if
                            end do
                        end if
                    end do
                end if
            end do
            if (Pass == 1) then
                do P = 1, 2
                    do J = 0, J_max
                        allocate(Orb_loc%Ind(P,J)%v(2,N_Orb_NN(P,J)))
                    end do
                end do
            end if
        end do
        call move_alloc(N_Orb_NN, Orb_loc%N)
    end subroutine orb2b_temp_make

    pure integer function O2b_temp_index(a,b,J,Orb_NN_t)
        integer, intent(in) :: a, b, J
        type(orb2b_temp_t), intent(in) :: Orb_NN_t
        O2b_temp_index = Orb_NN_t%Dic(a,b,J)
    end function O2b_temp_index

    ! Element of a temporary full-matrix operator ... counterpart of O2b_temp (O2b.jl:601) ... T = -1 pp, 0 pn, +1 nn
        ! (unlike O2b_index, where T = 0 is pn & T = 1 like particles) ... the pairs must exist (dev throws, the port
        ! does not check) ...
    pure real(dp) function O2b_temp(a,b,c,d,J,P,T,O,Orb_NN_t)
        integer, intent(in) :: a, b, c, d, J, P, T
        type(o2b_temp_t), intent(in) :: O
        type(orb2b_temp_t), intent(in) :: Orb_NN_t
        integer :: Bra, Ket
        Bra = Orb_NN_t%Dic(a,b,J)
        Ket = Orb_NN_t%Dic(c,d,J)
        select case (T)
        case (-1)
            O2b_temp = O%pp(P,J)%v(Bra,Ket)
        case (0)
            O2b_temp = O%pn(P,J)%v(Bra,Ket)
        case default
            O2b_temp = O%nn(P,J)%v(Bra,Ket)
        end select
    end function O2b_temp

    ! Empty Pandya-transformed operator over the pairs of orb2b_temp_make ... counterpart of O2b_Pandya_initialize ...
    subroutine O2b_Pandya_initialize(Params,Orb,F_NN,Orb_NN_t)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o2b_pandya_t), intent(out) :: F_NN
        type(orb2b_temp_t), intent(out), optional :: Orb_NN_t
        integer :: J_max, P, J, Length
        type(orb2b_temp_t) :: Orb_loc

        call orb2b_temp_make(Params,Orb,Orb_loc)
        J_max = Params%Calc%N2max + 1
        allocate(F_NN%pp(2,0:J_max), F_NN%pnnp(2,0:J_max), F_NN%pnpn(2,0:J_max), F_NN%nn(2,0:J_max))
        do P = 1, 2
            do J = 0, J_max
                Length = Orb_loc%N(P,J)
                allocate(F_NN%pp(P,J)%v(Length,Length), F_NN%pnnp(P,J)%v(Length,Length), &
                         F_NN%pnpn(P,J)%v(Length,Length), F_NN%nn(P,J)%v(Length,Length))
                F_NN%pp(P,J)%v = 0.0_dp
                F_NN%pnnp(P,J)%v = 0.0_dp
                F_NN%pnpn(P,J)%v = 0.0_dp
                F_NN%nn(P,J)%v = 0.0_dp
            end do
        end do
        if (present(Orb_NN_t)) then
            call move_alloc(Orb_loc%N, Orb_NN_t%N)
            call move_alloc(Orb_loc%Ind, Orb_NN_t%Ind)
            call move_alloc(Orb_loc%Dic, Orb_NN_t%Dic)
        end if
    end subroutine O2b_Pandya_initialize

    ! Accessors O2b_Pandya_pp / pnpn / pnnp / nn ...
    pure real(dp) function O2b_Pandya_pp(a,b,c,d,J,P,F,Orb_NN_t)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_pandya_t), intent(in) :: F
        type(orb2b_temp_t), intent(in) :: Orb_NN_t
        O2b_Pandya_pp = F%pp(P,J)%v(Orb_NN_t%Dic(a,b,J),Orb_NN_t%Dic(c,d,J))
    end function O2b_Pandya_pp

    pure real(dp) function O2b_Pandya_pnpn(a,b,c,d,J,P,F,Orb_NN_t)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_pandya_t), intent(in) :: F
        type(orb2b_temp_t), intent(in) :: Orb_NN_t
        O2b_Pandya_pnpn = F%pnpn(P,J)%v(Orb_NN_t%Dic(a,b,J),Orb_NN_t%Dic(c,d,J))
    end function O2b_Pandya_pnpn

    pure real(dp) function O2b_Pandya_pnnp(a,b,c,d,J,P,F,Orb_NN_t)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_pandya_t), intent(in) :: F
        type(orb2b_temp_t), intent(in) :: Orb_NN_t
        O2b_Pandya_pnnp = F%pnnp(P,J)%v(Orb_NN_t%Dic(a,b,J),Orb_NN_t%Dic(c,d,J))
    end function O2b_Pandya_pnnp

    pure real(dp) function O2b_Pandya_nn(a,b,c,d,J,P,F,Orb_NN_t)
        integer, intent(in) :: a, b, c, d, J, P
        type(o2b_pandya_t), intent(in) :: F
        type(orb2b_temp_t), intent(in) :: Orb_NN_t
        O2b_Pandya_nn = F%nn(P,J)%v(Orb_NN_t%Dic(a,b,J),Orb_NN_t%Dic(c,d,J))
    end function O2b_Pandya_nn

    ! Orbitals sharing l & 2j ...
    function Orb_PreComp(a_max,j,l,Orb) result(Orb_x)
        integer, intent(in) :: a_max, j, l
        type(orb1b_t), intent(in) :: Orb(:)
        integer, allocatable :: Orb_x(:)
        integer :: x, n
        allocate(Orb_x(count(Orb(1:a_max)%l == l .and. Orb(1:a_max)%j == j)))
        n = 0
        do x = 1, a_max
            if (Orb(x)%l == l .and. Orb(x)%j == j) then
                n = n + 1
                Orb_x(n) = Orb(x)%a
            end if
        end do
    end function Orb_PreComp

    ! Transformation of a 2-body operator O_NN into the basis given by C ... index by index ...
    subroutine O2b_transformation(Params,Orb,Orb_NN,O_NN,C,Report)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(inout) :: O_NN
        type(o1b_t), intent(in) :: C
        logical, intent(in), optional :: Report     ! dev keyword Report (default true) ... false suppresses the prints
        logical :: Rep
        type(o2b_temp_t) :: O_NN_t1, O_NN_t2
        type(orb2b_temp_t) :: Orb_NN_t
        integer :: a_max, J_max, J, P, N, Bra, Ket, a, b, cc, d, k, i, Bra_ai, Ket_id, Ind, N_T0, N_T1, Bra_ab, Ket_ci
        integer, allocatable :: Orb_x(:)
        real(dp) :: OppSum, OpnSum, OnnSum
        integer(i8) :: t0

        Rep = Verbose
        if (present(Report)) Rep = Rep .and. Report
        if (Rep) write(*,'(/,a)') 'Performing the Transformation of 2-body operator O_NN into the target basis due to C ...'
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        J_max = Params%Calc%N2max + 1

        ! Transformation of the 1st index ...
        call timer_start(t0)
        call O2b_temp_initialize(Params,Orb,O_NN_t1,Orb_NN_t)
        if (Rep) write(*,'(/,a)') 'Performing transformation of the given 2-body operator in the 1st index...'
        do J = 0, J_max
            do P = 1, 2
                N = Orb_NN_t%N(P,J)
                !$omp parallel do schedule(dynamic) private(a,b,Orb_x,Ket,cc,d,OppSum,OpnSum,OnnSum,k,i)
                do Bra = 1, N
                    a = Orb_NN_t%Ind(P,J)%v(1,Bra)
                    b = Orb_NN_t%Ind(P,J)%v(2,Bra)
                    Orb_x = Orb_PreComp(a_max,Orb(a)%j,Orb(a)%l,Orb)
                    do Ket = 1, N
                        cc = Orb_NN_t%Ind(P,J)%v(1,Ket)
                        d = Orb_NN_t%Ind(P,J)%v(2,Ket)
                        OppSum = 0.0_dp
                        OpnSum = 0.0_dp
                        OnnSum = 0.0_dp
                        do k = 1, size(Orb_x)
                            i = Orb_x(k)
                            OppSum = OppSum + C%p(i,a) * O2b_pp(i,b,cc,d,J,P,O_NN,Orb,Orb_NN)
                            OpnSum = OpnSum + C%p(i,a) * O2b_pn(i,b,cc,d,J,P,O_NN,Orb_NN)
                            OnnSum = OnnSum + C%n(i,a) * O2b_nn(i,b,cc,d,J,P,O_NN,Orb,Orb_NN)
                        end do
                        O_NN_t1%pp(P,J)%v(Bra,Ket) = OppSum
                        O_NN_t1%pn(P,J)%v(Bra,Ket) = OpnSum
                        O_NN_t1%nn(P,J)%v(Bra,Ket) = OnnSum
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        if (Rep) call timer_report(t0)

        ! Transformation of the 2nd index ...
        call timer_start(t0)
        call O2b_temp_initialize(Params,Orb,O_NN_t2)
        if (Rep) write(*,'(/,a)') 'Performing transformation of the given 2-body operator O_NN in the 2nd index...'
        do J = 0, J_max
            do P = 1, 2
                N = Orb_NN_t%N(P,J)
                !$omp parallel do schedule(dynamic) private(a,b,Orb_x,Ket,OppSum,OpnSum,OnnSum,k,i,Bra_ai)
                do Bra = 1, N
                    a = Orb_NN_t%Ind(P,J)%v(1,Bra)
                    b = Orb_NN_t%Ind(P,J)%v(2,Bra)
                    Orb_x = Orb_PreComp(a_max,Orb(b)%j,Orb(b)%l,Orb)
                    do Ket = 1, N
                        OppSum = 0.0_dp
                        OpnSum = 0.0_dp
                        OnnSum = 0.0_dp
                        do k = 1, size(Orb_x)
                            i = Orb_x(k)
                            Bra_ai = O2b_temp_index(a,i,J,Orb_NN_t)
                            OppSum = OppSum + C%p(i,b) * O_NN_t1%pp(P,J)%v(Bra_ai,Ket)
                            OpnSum = OpnSum + C%n(i,b) * O_NN_t1%pn(P,J)%v(Bra_ai,Ket)
                            OnnSum = OnnSum + C%n(i,b) * O_NN_t1%nn(P,J)%v(Bra_ai,Ket)
                        end do
                        O_NN_t2%pp(P,J)%v(Bra,Ket) = OppSum
                        O_NN_t2%pn(P,J)%v(Bra,Ket) = OpnSum
                        O_NN_t2%nn(P,J)%v(Bra,Ket) = OnnSum
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        if (Rep) call timer_report(t0)

        ! Transformation of the 3rd index ... result stored back into O_NN_t1 ...
        call timer_start(t0)
        if (Rep) write(*,'(/,a)') 'Performing transformation of the given 2-body operator O_NN in the 3rd index...'
        do J = 0, J_max
            do P = 1, 2
                N = Orb_NN_t%N(P,J)
                !$omp parallel do schedule(dynamic) private(cc,d,Orb_x,Bra,OppSum,OpnSum,OnnSum,k,i,Ket_id)
                do Ket = 1, N
                    cc = Orb_NN_t%Ind(P,J)%v(1,Ket)
                    d = Orb_NN_t%Ind(P,J)%v(2,Ket)
                    Orb_x = Orb_PreComp(a_max,Orb(cc)%j,Orb(cc)%l,Orb)
                    do Bra = 1, N
                        OppSum = 0.0_dp
                        OpnSum = 0.0_dp
                        OnnSum = 0.0_dp
                        do k = 1, size(Orb_x)
                            i = Orb_x(k)
                            Ket_id = O2b_temp_index(i,d,J,Orb_NN_t)
                            OppSum = OppSum + C%p(i,cc) * O_NN_t2%pp(P,J)%v(Bra,Ket_id)
                            OpnSum = OpnSum + C%p(i,cc) * O_NN_t2%pn(P,J)%v(Bra,Ket_id)
                            OnnSum = OnnSum + C%n(i,cc) * O_NN_t2%nn(P,J)%v(Bra,Ket_id)
                        end do
                        O_NN_t1%pp(P,J)%v(Bra,Ket) = OppSum
                        O_NN_t1%pn(P,J)%v(Bra,Ket) = OpnSum
                        O_NN_t1%nn(P,J)%v(Bra,Ket) = OnnSum
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        if (Rep) call timer_report(t0)

        ! Transformation of the 4th index ... back into the packed storage of O_NN ...
        call timer_start(t0)
        if (Rep) write(*,'(/,a)') 'Performing transformation of the given 2-body operator O_NN in the 4th index...'
        do J = 0, J_max
            do P = 1, 2
                N_T0 = Orb_NN%N(1,P,J)
                N_T1 = Orb_NN%N(2,P,J)
                !$omp parallel do schedule(dynamic) private(Ket,Ind,a,b,cc,d,Bra_ab,Orb_x,OppSum,OpnSum,OnnSum,k,i,Ket_ci)
                do Bra = 1, max(N_T0, N_T1)
                    do Ket = 1, Bra
                        ! Case of pn ... T = 0 ...
                        if (Bra <= N_T0) then
                            Ind = packed_index(Bra,Ket,N_T0)
                            a = Orb_NN%Ind(1,P,J)%v(1,Bra)
                            b = Orb_NN%Ind(1,P,J)%v(2,Bra)
                            cc = Orb_NN%Ind(1,P,J)%v(1,Ket)
                            d = Orb_NN%Ind(1,P,J)%v(2,Ket)
                            Bra_ab = O2b_temp_index(a,b,J,Orb_NN_t)
                            Orb_x = Orb_PreComp(a_max,Orb(d)%j,Orb(d)%l,Orb)
                            OpnSum = 0.0_dp
                            do k = 1, size(Orb_x)
                                i = Orb_x(k)
                                Ket_ci = O2b_temp_index(cc,i,J,Orb_NN_t)
                                OpnSum = OpnSum + C%n(i,d) * O_NN_t1%pn(P,J)%v(Bra_ab,Ket_ci)
                            end do
                            O_NN%pn(P,J)%v(Ind) = OpnSum
                        end if
                        ! Case of pp & nn ... T = 1 ...
                        if (Bra <= N_T1) then
                            Ind = packed_index(Bra,Ket,N_T1)
                            a = Orb_NN%Ind(2,P,J)%v(1,Bra)
                            b = Orb_NN%Ind(2,P,J)%v(2,Bra)
                            cc = Orb_NN%Ind(2,P,J)%v(1,Ket)
                            d = Orb_NN%Ind(2,P,J)%v(2,Ket)
                            Bra_ab = O2b_temp_index(a,b,J,Orb_NN_t)
                            Orb_x = Orb_PreComp(a_max,Orb(d)%j,Orb(d)%l,Orb)
                            OppSum = 0.0_dp
                            OnnSum = 0.0_dp
                            do k = 1, size(Orb_x)
                                i = Orb_x(k)
                                Ket_ci = O2b_temp_index(cc,i,J,Orb_NN_t)
                                OppSum = OppSum + C%p(i,d) * O_NN_t1%pp(P,J)%v(Bra_ab,Ket_ci)
                                OnnSum = OnnSum + C%n(i,d) * O_NN_t1%nn(P,J)%v(Bra_ab,Ket_ci)
                            end do
                            O_NN%pp(P,J)%v(Ind) = OppSum
                            O_NN%nn(P,J)%v(Ind) = OnnSum
                        end if
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        if (Rep) call timer_report(t0)

        if (Rep) write(*,'(/,a,/)') 'Given 2-body operator O_NN transformed ...'
    end subroutine O2b_transformation

    ! LEGACY (removed in cleanup stage) ... Pandya transformation into the 3-block o2b_temp_t (pn = dev pnnp) of the
        ! not-yet-ported qp chain ...
    subroutine O2b_Pandya_transformation_temp(Params,Orb,Orb_NN,O_NN,F_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: O_NN
        type(o2b_temp_t), intent(out) :: F_NN
        type(orb2b_temp_t) :: Orb_NN_t
        integer :: J_max, J, P, N, Bra, Ket, a, b, c, d

        J_max = Params%Calc%N2max + 1
        call O2b_temp_initialize(Params,Orb,F_NN,Orb_NN_t)
        write(*,'(/,a)') 'Performing the Pandya transformation of the given 2-body operator ..'
        do J = 0, J_max
            do P = 1, 2
                N = Orb_NN_t%N(P,J)
                !$omp parallel do schedule(dynamic) private(a,b,Ket,c,d)
                do Bra = 1, N
                    a = Orb_NN_t%Ind(P,J)%v(1,Bra)
                    b = Orb_NN_t%Ind(P,J)%v(2,Bra)
                    do Ket = 1, N
                        c = Orb_NN_t%Ind(P,J)%v(1,Ket)
                        d = Orb_NN_t%Ind(P,J)%v(2,Ket)
                        F_NN%pp(P,J)%v(Bra,Ket) = pandya_me(1,a,b,c,d,J,Orb,Orb_NN,O_NN)
                        F_NN%pn(P,J)%v(Bra,Ket) = pandya_me(0,a,b,c,d,J,Orb,Orb_NN,O_NN)
                        F_NN%nn(P,J)%v(Bra,Ket) = pandya_me(2,a,b,c,d,J,Orb,Orb_NN,O_NN)
                    end do
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine O2b_Pandya_transformation_temp

    ! Pandya (particle-hole) transformation of a 2-body operator into full matrices ... counterpart of
        ! O2b_Pandya_transformation (O2b.jl:1036) with the four blocks pp, pnnp, pnpn & nn ... J outer, P inner ...
    subroutine O2b_Pandya_transformation_4(Params,Orb,Orb_NN,O_NN,F_NN,Orb_NN_t)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: O_NN
        type(o2b_pandya_t), intent(out) :: F_NN
        type(orb2b_temp_t), intent(out) :: Orb_NN_t
        integer :: J_max, J, P, N, Bra, Ket, a, b, c, d

        J_max = Params%Calc%N2max + 1
        call O2b_Pandya_initialize(Params,Orb,F_NN,Orb_NN_t)
        write(*,'(/,a)') 'Performing the Pandya transformation of the given 2-body operator ..'
        do J = 0, J_max
            do P = 1, 2
                N = Orb_NN_t%N(P,J)
                !$omp parallel do schedule(dynamic) private(a,b,Ket,c,d)
                do Bra = 1, N
                    a = Orb_NN_t%Ind(P,J)%v(1,Bra)
                    b = Orb_NN_t%Ind(P,J)%v(2,Bra)
                    do Ket = 1, N
                        c = Orb_NN_t%Ind(P,J)%v(1,Ket)
                        d = Orb_NN_t%Ind(P,J)%v(2,Ket)
                        F_NN%pp(P,J)%v(Bra,Ket) = pandya_me(1,a,b,c,d,J,Orb,Orb_NN,O_NN)
                        F_NN%pnpn(P,J)%v(Bra,Ket) = pandya_me(3,a,b,c,d,J,Orb,Orb_NN,O_NN)
                        F_NN%pnnp(P,J)%v(Bra,Ket) = pandya_me(0,a,b,c,d,J,Orb,Orb_NN,O_NN)
                        F_NN%nn(P,J)%v(Bra,Ket) = pandya_me(2,a,b,c,d,J,Orb,Orb_NN,O_NN)
                    end do
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine O2b_Pandya_transformation_4

    ! One Pandya-transformed element ... Kind = 1 (pp), 0 (pnnp, the former pn), 2 (nn), 3 (pnpn: phase(1) = -1 & the
        ! exchanged ket (c,b) of O2b_pn) ... the factors multiply left to right as in dev (phase * hat * w6j * O) ...
    real(dp) function pandya_me(Kind,a,b,c,d,J,Orb,Orb_NN,O_NN) result(F)
        integer, intent(in) :: Kind, a, b, c, d, J
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: O_NN
        integer :: Q, I, j_a, j_b, j_c, j_d
        real(dp) :: ME
        Q = mod(Orb(a)%l + Orb(d)%l,2) + 1
        F = 0.0_dp
        if (Q == (mod(Orb(b)%l + Orb(c)%l,2) + 1)) then
            j_a = Orb(a)%j
            j_b = Orb(b)%j
            j_c = Orb(c)%j
            j_d = Orb(d)%j
            do I = max(abs(j_a - j_d)/2, abs(j_b - j_c)/2), min((j_a + j_d)/2, (j_b + j_c)/2)
                select case (Kind)
                case (1)
                    ME = O2b_pp(a,d,b,c,I,Q,O_NN,Orb,Orb_NN)
                case (0)
                    ME = O2b_pn(a,d,b,c,I,Q,O_NN,Orb_NN)
                case (3)
                    ME = O2b_pn(a,d,c,b,I,Q,O_NN,Orb_NN)
                case default
                    ME = O2b_nn(a,d,b,c,I,Q,O_NN,Orb,Orb_NN)
                end select
                if (Kind == 3) then
                    F = F + phase(1) * real(2*I + 1,dp) * f6j(j_a,j_b,2*J,j_c,j_d,2*I) * ME
                else
                    F = F + phase((j_b + j_c)/2 + I) * real(2*I + 1,dp) * f6j(j_a,j_b,2*J,j_c,j_d,2*I) * ME
                end if
            end do
        end if
    end function pandya_me

end module nj_o2b
