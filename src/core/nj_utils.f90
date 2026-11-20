! General utilities ... Julia-compatible number formatting, sorting, file system helpers & timers ...
module nj_utils
    use nj_kinds
    implicit none
    private

    public :: kronecker_delta, phase, iphase, parity_of, triangle, hat_of
    public :: jstr, istr, round_sig, round_digits, round_even, cfmt_f, cfmt_e, cfmt_d, cfmt_s
    public :: sortperm, sortperm_rev
    public :: make_dir, remove_dir, remove_file, is_dir
    public :: timer_start, timer_report
    public :: fatal
    public :: abspath_of, presize_file

    ! Integer to string for default & 64-bit integers ...
    interface istr
        module procedure istr_i4, istr_i8
    end interface istr

    ! Terminal output switch ... used to silence routines called repeatedly inside iterations,
        ! as the original does by redirecting stdout to /dev/null ...
    logical, public :: Verbose = .true.

contains

    pure integer function kronecker_delta(a,b)
        integer, intent(in) :: a, b
        kronecker_delta = merge(1, 0, a == b)
    end function kronecker_delta

    ! (-1)^Arg as Float64 ...
    pure real(dp) function phase(Arg)
        integer, intent(in) :: Arg
        phase = merge(-1.0_dp, 1.0_dp, modulo(Arg,2) == 1)
    end function phase

    ! (-1)^n as integer ...
    pure integer function iphase(n)
        integer, intent(in) :: n
        iphase = merge(-1, 1, modulo(n,2) == 1)
    end function iphase

    ! Julia parity(L) = rem(L,2) + 1 ... 1 = positive, 2 = negative parity of l_a + l_b ...
        ! Named parity_of because Parity is a local & namelist variable elsewhere (a use-associated name cannot be redeclared) ...
    pure integer function parity_of(L)
        integer, intent(in) :: L
        parity_of = mod(L,2) + 1
    end function parity_of

    ! Julia triangle(j1,j2,j3) ... |j1 - j2| <= j3 <= j1 + j2 for doubled angular momenta (no parity-of-sum test) ...
    pure logical function triangle(j_1,j_2,j_3)
        integer, intent(in) :: j_1, j_2, j_3
        triangle = (abs(j_1 - j_2) <= j_3) .and. (j_3 <= (j_1 + j_2))
    end function triangle

    ! Julia hat(j) = Float64(j + 1) for a doubled j (hat_of(2*I) = 2I + 1) ... named hat_of because Hat is a local
        ! variable of several solvers ...
    pure real(dp) function hat_of(j)
        integer, intent(in) :: j
        hat_of = real(j + 1,dp)
    end function hat_of

    subroutine fatal(Msg)
        character(len=*), intent(in) :: Msg
        write(*,'(a)') 'ERROR: ' // trim(Msg)
        error stop 1
    end subroutine fatal

    ! Integer to string ... like Julia string(::Int) ...
    function istr_i4(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: Buf
        write(Buf,'(I0)') i
        s = trim(Buf)
    end function istr_i4

    function istr_i8(i) result(s)
        integer(i8), intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: Buf
        write(Buf,'(I0)') i
        s = trim(Buf)
    end function istr_i8

    ! Absolute, normalized path ... like Julia abspath (relative paths are taken w.r.t. $PWD; '.', '..' & '//' resolved,
        ! symbolic links are not) ...
    function abspath_of(Path) result(s)
        character(len=*), intent(in) :: Path
        character(len=:), allocatable :: s, Full, Comp
        character(len=4096) :: Cwd
        integer :: i, j, n, Stat
        character(len=4096), allocatable :: Parts(:)
        if (len_trim(Path) > 0 .and. Path(1:1) == '/') then
            Full = trim(Path)
        else
            call get_environment_variable('PWD', Cwd, status=Stat)
            if (Stat /= 0) Cwd = '.'
            Full = trim(Cwd) // '/' // trim(Path)
        end if
        allocate(Parts(len(Full)/2 + 1))
        n = 0
        i = 1
        do while (i <= len(Full))
            j = index(Full(i:), '/')
            if (j == 0) then
                Comp = Full(i:)
                i = len(Full) + 1
            else
                Comp = Full(i:i+j-2)
                i = i + j
            end if
            if (len(Comp) == 0 .or. Comp == '.') cycle
            if (Comp == '..') then
                n = max(n - 1, 0)
            else
                n = n + 1
                Parts(n) = Comp
            end if
        end do
        s = ''
        do i = 1, n
            s = s // '/' // trim(Parts(i))
        end do
        if (n == 0) s = '/'
    end function abspath_of

    ! Creates (or replaces) a binary file of Size bytes (zero-filled) ... Julia open(File,"w") + truncate(io,Size) ...
    subroutine presize_file(File,Size)
        character(len=*), intent(in) :: File
        integer(i8), intent(in) :: Size
        integer :: u, Ios
        character(len=512) :: Msg
        open(newunit=u, file=trim(File), access='stream', form='unformatted', status='replace', action='write', &
             iostat=Ios, iomsg=Msg)
        if (Ios /= 0) call fatal('Cannot create ' // trim(File) // ' ... ' // trim(Msg))
        if (Size > 0) write(u, pos=Size) achar(0)
        close(u)
    end subroutine presize_file

    ! Shortest round-trip representation of a Float64 ... mimics Julia string(::Float64) ...
        ! Plain decimal notation for 1e-4 <= |x| < 1e6, scientific "d.ddde±x" otherwise ...
    function jstr(x) result(s)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: s
        character(len=64) :: Buf
        character(len=32) :: Digits
        real(dp) :: y
        integer :: d, e10, n, Ios, Pos_e

        if (x /= x) then
            s = 'NaN'
            return
        end if
        if (abs(x) > huge(x)) then
            s = merge('-Inf', ' Inf', x < 0.0_dp)
            s = trim(adjustl(s))
            return
        end if
        if (x == 0.0_dp) then
            if (sign(1.0_dp,x) < 0.0_dp) then
                s = '-0.0'
            else
                s = '0.0'
            end if
            return
        end if

        ! Find the minimal number of significant digits that round-trips ...
        do d = 1, 17
            write(Buf,'(ES40.' // istr(d-1) // 'E3)') x
            read(Buf,*,iostat=Ios) y
            if (y == x) exit
        end do
        Buf = adjustl(Buf)
        Pos_e = index(Buf,'E')
        read(Buf(Pos_e+1:),*) e10

        ! Collect the significant digits without sign & decimal point ...
        Digits = ''
        n = 0
        do d = 1, Pos_e - 1
            if (Buf(d:d) >= '0' .and. Buf(d:d) <= '9') then
                n = n + 1
                Digits(n:n) = Buf(d:d)
            end if
        end do
        ! Strip trailing zeros ... keep at least one digit ...
        do while (n > 1 .and. Digits(n:n) == '0')
            n = n - 1
        end do

        if (e10 >= -5+1 .and. e10 < 6) then
            if (e10 >= 0) then
                if (n <= e10 + 1) then
                    s = Digits(1:n) // repeat('0', e10 + 1 - n) // '.0'
                else
                    s = Digits(1:e10+1) // '.' // Digits(e10+2:n)
                end if
            else
                s = '0.' // repeat('0', -e10 - 1) // Digits(1:n)
            end if
        else
            if (n == 1) then
                s = Digits(1:1) // '.0e' // istr(e10)
            else
                s = Digits(1:1) // '.' // Digits(2:n) // 'e' // istr(e10)
            end if
        end if
        if (x < 0.0_dp) s = '-' // s
    end function jstr

    ! Julia round(x, sigdigits = n) ...
    function round_sig(x,n) result(y)
        real(dp), intent(in) :: x
        integer, intent(in) :: n
        real(dp) :: y
        character(len=64) :: Buf
        if (x == 0.0_dp .or. x /= x) then
            y = x
            return
        end if
        write(Buf,'(ES40.' // istr(n-1) // 'E3)') x
        read(Buf,*) y
    end function round_sig

    ! Julia round(x) for integers ... RoundNearest, ties to even ...
    pure integer function round_even(x)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = floor(x)
        if (x - f > 0.5_dp) then
            round_even = int(f) + 1
        else if (x - f < 0.5_dp) then
            round_even = int(f)
        else
            round_even = int(f) + modulo(int(f),2)
        end if
    end function round_even

    ! Julia round(x, digits = d) ... scaling by 10^d with ties to even ...
    function round_digits(x,d) result(y)
        real(dp), intent(in) :: x
        integer, intent(in) :: d
        real(dp) :: y, s, f, t
        if (x /= x .or. x == 0.0_dp) then
            y = x
            return
        end if
        s = 10.0_dp**d
        t = x * s
        f = floor(t)
        if (t - f > 0.5_dp) then
            f = f + 1.0_dp
        else if (t - f == 0.5_dp) then
            if (modulo(f,2.0_dp) /= 0.0_dp) f = f + 1.0_dp
        end if
        y = f / s
    end function round_digits

    ! C-like "%W.Pf" / "%-W.Pf" formatting ... never prints asterisks ...
    function cfmt_f(x,W,P,Left) result(s)
        real(dp), intent(in) :: x
        integer, intent(in) :: W, P
        logical, intent(in), optional :: Left
        character(len=:), allocatable :: s
        character(len=128) :: Buf
        logical :: Lj
        Lj = .false.
        if (present(Left)) Lj = Left
        write(Buf,'(F127.' // istr(P) // ')') x
        s = trim(adjustl(Buf))
        if (s(1:1) == '.') s = '0' // s
        if (len(s) > 1) then
            if (s(1:2) == '-.') s = '-0' // s(2:)
        end if
        s = pad(s,W,Lj)
    end function cfmt_f

    ! C-like "%W.Pe" formatting ... lowercase e, at least 2 exponent digits ...
    function cfmt_e(x,W,P,Left) result(s)
        real(dp), intent(in) :: x
        integer, intent(in) :: W, P
        logical, intent(in), optional :: Left
        character(len=:), allocatable :: s
        character(len=128) :: Buf
        integer :: Pos_e, e10
        logical :: Lj
        Lj = .false.
        if (present(Left)) Lj = Left
        write(Buf,'(ES127.' // istr(P) // 'E3)') x
        Buf = adjustl(Buf)
        Pos_e = index(Buf,'E')
        read(Buf(Pos_e+1:),*) e10
        s = Buf(1:Pos_e-1) // 'e'
        if (e10 < 0) then
            s = s // '-'
        else
            s = s // '+'
        end if
        if (abs(e10) < 10) then
            s = s // '0' // istr(abs(e10))
        else
            s = s // istr(abs(e10))
        end if
        s = pad(s,W,Lj)
    end function cfmt_e

    ! C-like "%Wd" / "%-Wd" formatting ...
    function cfmt_d(i,W,Left) result(s)
        integer, intent(in) :: i, W
        logical, intent(in), optional :: Left
        character(len=:), allocatable :: s
        logical :: Lj
        Lj = .false.
        if (present(Left)) Lj = Left
        s = pad(istr(i),W,Lj)
    end function cfmt_d

    ! C-like "%Ws" / "%-Ws" formatting ...
    function cfmt_s(Str,W,Left) result(s)
        character(len=*), intent(in) :: Str
        integer, intent(in) :: W
        logical, intent(in), optional :: Left
        character(len=:), allocatable :: s
        logical :: Lj
        Lj = .false.
        if (present(Left)) Lj = Left
        s = pad(Str,W,Lj)
    end function cfmt_s

    function pad(Str,W,Left) result(s)
        character(len=*), intent(in) :: Str
        integer, intent(in) :: W
        logical, intent(in) :: Left
        character(len=:), allocatable :: s
        if (len(Str) >= W) then
            s = Str
        else if (Left) then
            s = Str // repeat(' ', W - len(Str))
        else
            s = repeat(' ', W - len(Str)) // Str
        end if
    end function pad

    ! Stable ascending sortperm ... as Julia sortperm(x) ... (merge sort) ...
    subroutine sortperm(x,Perm)
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: Perm(:)
        integer, allocatable :: Tmp(:)
        integer :: i, n
        n = size(x)
        do i = 1, n
            Perm(i) = i
        end do
        allocate(Tmp(n))
        call merge_sort(1,n)
    contains
        recursive subroutine merge_sort(Lo,Hi)
            integer, intent(in) :: Lo, Hi
            integer :: Mid, i, j, k
            if (Hi <= Lo) return
            Mid = (Lo + Hi) / 2
            call merge_sort(Lo,Mid)
            call merge_sort(Mid+1,Hi)
            i = Lo
            j = Mid + 1
            k = Lo
            do while (i <= Mid .and. j <= Hi)
                ! Take from the left run on ties ... keeps the sort stable ...
                if (x(Perm(j)) < x(Perm(i))) then
                    Tmp(k) = Perm(j)
                    j = j + 1
                else
                    Tmp(k) = Perm(i)
                    i = i + 1
                end if
                k = k + 1
            end do
            do while (i <= Mid)
                Tmp(k) = Perm(i)
                i = i + 1
                k = k + 1
            end do
            do while (j <= Hi)
                Tmp(k) = Perm(j)
                j = j + 1
                k = k + 1
            end do
            Perm(Lo:Hi) = Tmp(Lo:Hi)
        end subroutine merge_sort
    end subroutine sortperm

    ! Stable descending sortperm ... as Julia sortperm(x, rev = true) ...
    subroutine sortperm_rev(x,Perm)
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: Perm(:)
        call sortperm(-x,Perm)
    end subroutine sortperm_rev

    logical function is_dir(Path)
        character(len=*), intent(in) :: Path
        integer :: Stat
        call execute_command_line('test -d "' // trim(Path) // '"', exitstat=Stat)
        is_dir = (Stat == 0)
    end function is_dir

    subroutine make_dir(Path)
        character(len=*), intent(in) :: Path
        integer :: Stat
        call execute_command_line('mkdir -p "' // trim(Path) // '"', exitstat=Stat)
        if (Stat /= 0) call fatal('Cannot create directory ' // trim(Path))
    end subroutine make_dir

    ! Deletes a file if it exists (Julia rm(path)) ...
    subroutine remove_file(Path)
        character(len=*), intent(in) :: Path
        integer :: u, Ios
        logical :: Exists
        inquire(file=trim(Path), exist=Exists)
        if (.not. Exists) return
        open(newunit=u, file=trim(Path), status='old', iostat=Ios)
        if (Ios == 0) close(u, status='delete')
    end subroutine remove_file

    subroutine remove_dir(Path)
        character(len=*), intent(in) :: Path
        integer :: Stat
        if (len_trim(Path) == 0) return
        call execute_command_line('rm -rf "' // trim(Path) // '"', exitstat=Stat)
    end subroutine remove_dir

    ! Wall-clock timers ... a replacement of Julia's @time ...
    subroutine timer_start(t0)
        integer(i8), intent(out) :: t0
        call system_clock(t0)
    end subroutine timer_start

    subroutine timer_report(t0,Label)
        integer(i8), intent(in) :: t0
        character(len=*), intent(in), optional :: Label
        integer(i8) :: t1, Rate
        call system_clock(t1,Rate)
        if (.not. Verbose) return
        if (present(Label)) then
            write(*,'(a)') '  ' // cfmt_f(real(t1 - t0,dp) / real(Rate,dp),0,6) // ' seconds ... ' // Label
        else
            write(*,'(a)') '  ' // cfmt_f(real(t1 - t0,dp) / real(Rate,dp),0,6) // ' seconds'
        end if
    end subroutine timer_report

end module nj_utils
