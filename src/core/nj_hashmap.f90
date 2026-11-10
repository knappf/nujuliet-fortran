! Open-addressing hash map Int64 -> Int32 ... replacement for Julia's Dict{UInt64,Int64} ...
    ! Keys must be non-negative ... the empty slot marker is -1 ...
module nj_hashmap
    use nj_kinds
    implicit none
    private

    public :: hashmap_t, hashmap_lookup

    integer(i8), parameter :: Empty = -1_i8

    type :: hashmap_t
        integer(i8), allocatable :: Keys(:)
        integer, allocatable :: Vals(:)
        integer :: Count = 0
        integer(i8) :: Mask = 0
    contains
        procedure :: init => hashmap_init
        procedure :: set => hashmap_set
        procedure :: get => hashmap_get
        procedure :: find => hashmap_find
    end type hashmap_t

contains

    pure integer(i8) function mix(Key)
        integer(i8), intent(in) :: Key
        integer(i8) :: z
        ! xorshift64 mixing rounds ... bit operations only, no signed overflow ...
        z = ieor(Key, 88172645463325252_i8)
        z = ieor(z, shiftl(z,13))
        z = ieor(z, shiftr(z,7))
        z = ieor(z, shiftl(z,17))
        z = ieor(z, shiftr(z,29))
        z = ieor(z, shiftl(z,11))
        mix = ieor(z, shiftr(z,23))
    end function mix

    subroutine hashmap_init(Map,Capacity)
        class(hashmap_t), intent(inout) :: Map
        integer, intent(in) :: Capacity
        integer(i8) :: n
        n = 16
        do while (n < 2_i8 * int(Capacity,i8) + 16)
            n = n * 2
        end do
        if (allocated(Map%Keys)) deallocate(Map%Keys, Map%Vals)
        allocate(Map%Keys(0:n-1), Map%Vals(0:n-1))
        Map%Keys = Empty
        Map%Vals = 0
        Map%Count = 0
        Map%Mask = n - 1
    end subroutine hashmap_init

    subroutine grow(Map)
        class(hashmap_t), intent(inout) :: Map
        integer(i8), allocatable :: Old_Keys(:)
        integer, allocatable :: Old_Vals(:)
        integer(i8) :: i
        call move_alloc(Map%Keys, Old_Keys)
        call move_alloc(Map%Vals, Old_Vals)
        call hashmap_init(Map, int(size(Old_Keys,kind=i8)))
        do i = 0, size(Old_Keys,kind=i8) - 1
            if (Old_Keys(i) /= Empty) call hashmap_set(Map, Old_Keys(i), Old_Vals(i))
        end do
    end subroutine grow

    subroutine hashmap_set(Map,Key,Val)
        class(hashmap_t), intent(inout) :: Map
        integer(i8), intent(in) :: Key
        integer, intent(in) :: Val
        integer(i8) :: i
        if (.not. allocated(Map%Keys)) call hashmap_init(Map,1024)
        if (2_i8 * int(Map%Count + 1,i8) > Map%Mask) call grow(Map)
        i = iand(mix(Key), Map%Mask)
        do
            if (Map%Keys(i) == Empty) then
                Map%Keys(i) = Key
                Map%Vals(i) = Val
                Map%Count = Map%Count + 1
                return
            else if (Map%Keys(i) == Key) then
                Map%Vals(i) = Val
                return
            end if
            i = iand(i + 1, Map%Mask)
        end do
    end subroutine hashmap_set

    ! Returns the value stored for Key, or -1 if the Key is absent ...
    pure integer function hashmap_find(Map,Key)
        class(hashmap_t), intent(in) :: Map
        integer(i8), intent(in) :: Key
        integer(i8) :: i
        hashmap_find = -1
        if (.not. allocated(Map%Keys)) return
        i = iand(mix(Key), Map%Mask)
        do
            if (Map%Keys(i) == Key) then
                hashmap_find = Map%Vals(i)
                return
            else if (Map%Keys(i) == Empty) then
                return
            end if
            i = iand(i + 1, Map%Mask)
        end do
    end function hashmap_find

    ! As hashmap_find with a non-polymorphic argument ... for hot threaded loops: a type-bound call on class(hashmap_t)
    ! made ifx build the class descriptor in static memory on every call (all threads writing the same cache lines) ...
    pure integer function hashmap_lookup(Map,Key)
        type(hashmap_t), intent(in) :: Map
        integer(i8), intent(in) :: Key
        integer(i8) :: i
        hashmap_lookup = -1
        if (.not. allocated(Map%Keys)) return
        i = iand(mix(Key), Map%Mask)
        do
            if (Map%Keys(i) == Key) then
                hashmap_lookup = Map%Vals(i)
                return
            else if (Map%Keys(i) == Empty) then
                return
            end if
            i = iand(i + 1, Map%Mask)
        end do
    end function hashmap_lookup

    ! As find ... but a missing key is a fatal error ... same as Julia's KeyError ...
    integer function hashmap_get(Map,Key)
        class(hashmap_t), intent(in) :: Map
        integer(i8), intent(in) :: Key
        hashmap_get = hashmap_find(Map,Key)
        if (hashmap_get < 0) then
            write(*,'(a,I0)') 'ERROR: missing key in hash map lookup ... key = ', Key
            error stop 1
        end if
    end function hashmap_get

end module nj_hashmap
