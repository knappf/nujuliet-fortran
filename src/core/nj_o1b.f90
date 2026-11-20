! 1-body operators, kinetic energy & 2-body kinetic CM terms ... counterpart of O1B.jl, T1B.jl & T2B.jl ...
module nj_o1b
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    use nj_linalg, only : mtm, mtmm
    implicit none
    private

    public :: O1b_export, O1b_import, qpO1B_export, qpO1B_import, qpN1B_make
    public :: T1b, T1b_energy, T2b, rGrad, o1b_alloc
    public :: O1b_transformation, O1b_transpose, O1b_diagonal, O1b_deepcopy, O1b_zeroes, O1b_ones, O1b_unit_matrix

contains

    subroutine o1b_alloc(O,a_max)
        type(o1b_t), intent(inout) :: O
        integer, intent(in) :: a_max
        if (allocated(O%p)) deallocate(O%p)
        if (allocated(O%n)) deallocate(O%n)
        allocate(O%p(a_max,a_max), O%n(a_max,a_max))
        O%p = 0.0_dp
        O%n = 0.0_dp
    end subroutine o1b_alloc

    ! Helpers of src/common/O1b.jl (NuJuliet dev) ... O1b_transformation(O,C) = C^T O C per isospin (as C' * O * C) ...
    function O1b_transformation(O,C) result(R)
        type(o1b_t), intent(in) :: O, C
        type(o1b_t) :: R
        R%p = mtmm(C%p,O%p)
        R%n = mtmm(C%n,O%n)
    end function O1b_transformation

    function O1b_transpose(O) result(R)
        type(o1b_t), intent(in) :: O
        type(o1b_t) :: R
        R%p = transpose(O%p)
        R%n = transpose(O%n)
    end function O1b_transpose

    ! diagm(diag(O)) per isospin ...
    function O1b_diagonal(O) result(R)
        type(o1b_t), intent(in) :: O
        type(o1b_t) :: R
        integer :: a
        call o1b_alloc(R,size(O%p,1))
        do a = 1, size(O%p,1)
            R%p(a,a) = O%p(a,a)
            R%n(a,a) = O%n(a,a)
        end do
    end function O1b_diagonal

    ! Julia deepcopy ... a plain assignment copies the allocatable components ...
    function O1b_deepcopy(O) result(R)
        type(o1b_t), intent(in) :: O
        type(o1b_t) :: R
        R = O
    end function O1b_deepcopy

    ! a_max(Calc.Nmax) square operators filled with 0, 1 & the unit matrix ...
    function O1b_zeroes(Params) result(R)
        type(params_t), intent(in) :: Params
        type(o1b_t) :: R
        call o1b_alloc(R,(Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2)
    end function O1b_zeroes

    function O1b_ones(Params) result(R)
        type(params_t), intent(in) :: Params
        type(o1b_t) :: R
        call o1b_alloc(R,(Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2)
        R%p = 1.0_dp
        R%n = 1.0_dp
    end function O1b_ones

    function O1b_unit_matrix(Params) result(R)
        type(params_t), intent(in) :: Params
        type(o1b_t) :: R
        integer :: a
        call o1b_alloc(R,(Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2)
        do a = 1, size(R%p,1)
            R%p(a,a) = 1.0_dp
            R%n(a,a) = 1.0_dp
        end do
    end function O1b_unit_matrix

    ! Binary export of the (l,j)-block-diagonal entries of a 1-body operator ... raw Float64 stream ...
    subroutine O1b_export(Params,Orb,A,A_Export_Path)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t), intent(in) :: A
        character(len=*), intent(in) :: A_Export_Path
        integer :: u
        open(newunit=u, file=A_Export_Path, access='stream', form='unformatted', status='replace', action='write')
        call write_block(u,Params,Orb,A%p)
        call write_block(u,Params,Orb,A%n)
        close(u)
    end subroutine O1b_export

    function O1b_import(Params,Orb,Import_Path) result(A)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: Import_Path
        type(o1b_t) :: A
        integer :: u, a_max, Ios
        character(len=256) :: Msg
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(A,a_max)
        open(newunit=u, file=Import_Path, access='stream', form='unformatted', status='old', action='read', &
             iostat=Ios, iomsg=Msg)
        if (Ios /= 0) call fatal('O1b_import: cannot open ' // trim(Import_Path) // ' ... ' // trim(Msg))
        call read_block(u,Params,Orb,A%p)
        call read_block(u,Params,Orb,A%n)
        close(u)
    end function O1b_import

    subroutine qpO1B_export(Params,Orb,A,A_Export_Path)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(qpo1b_t), intent(in) :: A
        character(len=*), intent(in) :: A_Export_Path
        integer :: u
        open(newunit=u, file=A_Export_Path, access='stream', form='unformatted', status='replace', action='write')
        call write_block(u,Params,Orb,A%qp11%p)
        call write_block(u,Params,Orb,A%qp11%n)
        call write_block(u,Params,Orb,A%qp20%p)
        call write_block(u,Params,Orb,A%qp20%n)
        close(u)
    end subroutine qpO1B_export

    function qpO1B_import(Params,Orb,Import_Path) result(A)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        character(len=*), intent(in) :: Import_Path
        type(qpo1b_t) :: A
        integer :: u, a_max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        call o1b_alloc(A%qp11,a_max)
        call o1b_alloc(A%qp20,a_max)
        open(newunit=u, file=Import_Path, access='stream', form='unformatted', status='old', action='read')
        call read_block(u,Params,Orb,A%qp11%p)
        call read_block(u,Params,Orb,A%qp11%n)
        call read_block(u,Params,Orb,A%qp20%p)
        call read_block(u,Params,Orb,A%qp20%n)
        close(u)
    end function qpO1B_import

    subroutine write_block(u,Params,Orb,M)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(in) :: M(:,:)
        integer :: a, b, a_max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j) write(u) M(a,b)
            end do
        end do
    end subroutine write_block

    subroutine read_block(u,Params,Orb,M)
        integer, intent(in) :: u
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp), intent(inout) :: M(:,:)
        integer :: a, b, a_max
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        do a = 1, a_max
            do b = 1, a_max
                if (Orb(a)%l == Orb(b)%l .and. Orb(a)%j == Orb(b)%j) read(u) M(a,b)
            end do
        end do
    end subroutine read_block

    function qpN1B_make(Params,U,V) result(N_qp)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: U, V
        type(qpo1b_t) :: N_qp
        N_qp%qp11%p = mtm(U%p,U%p) - mtm(V%p,V%p)
        N_qp%qp11%n = mtm(U%n,U%n) - mtm(V%n,V%n)
        N_qp%qp20%p = mtm(U%p,V%p)
        N_qp%qp20%n = mtm(U%n,V%n)
    end function qpN1B_make

    ! 1-body kinetic operator in the LHO basis ...
        ! Fix F4: the size is that of the Orb passed in (dev T1b.jl:3 uses a_max(Int.Nmax), out of bounds for a
        ! Calc.Nmax-sized Orb when Int.Nmax > Calc.Nmax) ...
    function T1b(Params,Orb) result(T)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(o1b_t) :: T
        integer :: a_max, k, l, n_k, l_k, j_k, n_l, l_l, j_l
        real(dp) :: hw, NN_k
        real(dp), allocatable :: TM(:,:)

        hw = Params%Calc%hw
        a_max = size(Orb)
        allocate(TM(a_max,a_max))
        TM = 0.0_dp
        do k = 1, a_max
            n_k = Orb(k)%n
            l_k = Orb(k)%l
            j_k = Orb(k)%j
            NN_k = real(2*n_k + l_k,dp)
            TM(k,k) = TM(k,k) + 0.5_dp * hw * (NN_k + 1.5_dp)
            do l = 1, a_max
                n_l = Orb(l)%n
                l_l = Orb(l)%l
                j_l = Orb(l)%j
                if (l_k == l_l .and. j_k == j_l) then
                    if (n_k == (n_l + 1)) TM(k,l) = TM(k,l) + 0.5_dp * hw * sqrt(0.5_dp * real(n_k * (2*n_k + 2*l_k + 1),dp))
                    if ((n_k + 1) == n_l) TM(k,l) = TM(k,l) + 0.5_dp * hw * sqrt(0.5_dp * real(n_l * (2*n_l + 2*l_l + 1),dp))
                end if
            end do
        end do
        T%p = TM
        T%n = TM
    end function T1b

    ! Total kinetic energy of the reference state including the CM corrections ...
    function T1b_energy(Params,Rho,Orb,T) result(t_tot)
        type(params_t), intent(in) :: Params
        type(o1b_t), intent(in) :: Rho, T
        type(orb1b_t), intent(in) :: Orb(:)
        real(dp) :: t_tot
        integer :: A_mass, N_max, N_2max, a_max, a, d, b, e, J
        integer :: n_a, l_a, j_a, n_d, l_d, j_d, n_b, l_b, j_b, n_e, l_e, j_e
        real(dp) :: hw, Sum, T_NN_Amp, TNN_sym, TNN_antisym, Amp, Amp_2
        character(len=32) :: CMS

        A_mass = Params%Calc%A
        hw = Params%Calc%hw
        N_max = Params%Calc%Nmax
        N_2max = Params%Calc%N2max
        CMS = Params%Calc%CMS
        a_max = (N_max + 1)*(N_max + 2)/2

        write(*,'(/,a)') 'Calculating the total kinetic energy of given ground-state reference state ...'
        t_tot = 0.0_dp
        !$omp parallel do schedule(static) reduction(+:t_tot) &
        !$omp private(Sum,n_a,l_a,j_a,d,n_d,l_d,j_d,b,n_b,l_b,j_b,e,n_e,l_e,j_e,J,T_NN_Amp,TNN_sym,TNN_antisym,Amp,Amp_2)
        do a = 1, a_max
            Sum = 0.0_dp
            n_a = Orb(a)%n
            l_a = Orb(a)%l
            j_a = Orb(a)%j
            do d = 1, a_max
                n_d = Orb(d)%n
                l_d = Orb(d)%l
                j_d = Orb(d)%j
                if (l_a == l_d .and. j_a == j_d) then
                    ! 2-body CM correction to the kinetic energy ...
                    do b = 1, a_max
                        n_b = Orb(b)%n
                        l_b = Orb(b)%l
                        j_b = Orb(b)%j
                        if ((2*(n_a + n_b) + l_a + l_b) <= N_2max) then
                            do e = 1, a_max
                                n_e = Orb(e)%n
                                l_e = Orb(e)%l
                                j_e = Orb(e)%j
                                if (l_b == l_e .and. j_b == j_e .and. (2*(n_d + n_e) + l_d + l_e) <= N_2max) then
                                    if (mod(l_a + l_b,2) == mod(l_d + l_e,2)) then
                                        do J = abs(j_a - j_b)/2, (j_a + j_b)/2
                                            T_NN_Amp = real(2*J + 1,dp) / real(A_mass,dp)
                                            if (CMS == 'CMS1+2B') then
                                                TNN_sym = T2b(Orb,a,b,d,e,J) * hw
                                                TNN_antisym = 1.0_dp / sqrt(real((1 + kronecker_delta(a,b))*(1 + kronecker_delta(d,e)),dp)) * &
                                                              (T2b(Orb,a,b,d,e,J) - phase((j_d + j_e)/2 - J) * T2b(Orb,a,b,e,d,J)) * hw
                                            else if (CMS == 'CMS2B') then
                                                Amp = hw / sqrt(real(1 + kronecker_delta(a,b),dp) * real(1 + kronecker_delta(d,e),dp))
                                                Amp_2 = 1.0_dp / sqrt(real(1 + kronecker_delta(a,b),dp) * real(1 + kronecker_delta(d,e),dp))
                                                TNN_sym = hw * T2b(Orb,a,b,d,e,J) + T%p(a,d) * kronecker_delta(b,e) + &
                                                          T%p(b,e) * kronecker_delta(a,d)
                                                TNN_antisym = (Amp * T2b(Orb,a,b,d,e,J) + Amp_2 * (kronecker_delta(b,e) * T%p(a,d) + &
                                                              kronecker_delta(a,d) * T%p(b,e)) - phase((j_d + j_e)/2 - J) * &
                                                              (Amp * T2b(Orb,a,b,e,d,J) + Amp_2 * (real(kronecker_delta(b,d),dp) * &
                                                              T%p(a,e) + real(kronecker_delta(a,e),dp) * T%p(b,d)))) / real(A_mass,dp)
                                            else
                                                TNN_sym = 0.0_dp
                                                TNN_antisym = 0.0_dp
                                            end if
                                            Sum = Sum + T_NN_Amp * TNN_antisym * Rho%p(b,e)
                                            Sum = Sum + T_NN_Amp * TNN_sym * Rho%n(b,e)
                                            Sum = Sum + T_NN_Amp * TNN_antisym * Rho%n(b,e)
                                            Sum = Sum + T_NN_Amp * TNN_sym * Rho%p(b,e)
                                        end do
                                    end if
                                end if
                            end do
                        end if
                    end do

                    ! 1-body kinetic operator & CMS correction ...
                    if (CMS == 'CMS1+2B') then
                        Sum = Sum + T%p(a,d) * (1.0_dp - 1.0_dp / real(A_mass,dp)) * Rho%p(a,d) * real(Orb(a)%j + 1,dp)
                        Sum = Sum + T%n(a,d) * (1.0_dp - 1.0_dp / real(A_mass,dp)) * Rho%n(a,d) * real(Orb(a)%j + 1,dp)
                    else if (CMS /= 'CMS2B') then
                        Sum = Sum + T%p(a,d) * Rho%p(a,d)
                        Sum = Sum + T%n(a,d) * Rho%n(a,d)
                    end if
                end if
            end do
            t_tot = t_tot + Sum
        end do
        !$omp end parallel do

        write(*,'(a)') TAB // 'Total kinetic energy reads:    <T> = ' // jstr(t_tot) // ' MeV'
    end function T1b_energy

    ! 2-body kinetic matrix element (gradient product) ...
    pure real(dp) function T2b(Orb,a,b,c,d,J)
        type(orb1b_t), intent(in) :: Orb(:)
        integer, intent(in) :: a, b, c, d, J
        real(dp) :: Grad1, Grad2
        Grad1 = rGrad(a,c,Orb)
        Grad2 = rGrad(b,d,Orb)
        T2b = phase((Orb(b)%j + Orb(c)%j)/2 + J) * Grad1 * Grad2 * f6j(Orb(a)%j, Orb(b)%j, 2*J, Orb(d)%j, Orb(c)%j, 2)
    end function T2b

    ! Reduced matrix element of the gradient operator in the LHO basis ...
    pure real(dp) function rGrad(a,b,Orb)
        integer, intent(in) :: a, b
        type(orb1b_t), intent(in) :: Orb(:)
        integer :: n_a, l_a, j_a, n_b, l_b, j_b
        real(dp) :: Amp
        n_a = Orb(a)%n
        l_a = Orb(a)%l
        j_a = Orb(a)%j
        n_b = Orb(b)%n
        l_b = Orb(b)%l
        j_b = Orb(b)%j
        Amp = f6j(j_a,j_b,2,2*l_b,2*l_a,1) * sqrt(real((j_a + 1) * (j_b + 1),dp)) * phase(l_a + (j_b + 1)/2)
        rGrad = Amp * (sqrt(real(l_b+1,dp)*(real(n_b+l_b,dp)+3.0_dp/2.0_dp)) * real(kronecker_delta(l_a,l_b+1)*kronecker_delta(n_a,n_b),dp) + &
                sqrt(real((l_b+1)*n_b,dp)) * real(kronecker_delta(l_a,l_b+1)*kronecker_delta(n_a,n_b-1),dp) + &
                sqrt(real(l_b,dp)*(real(n_b+l_b,dp)+1.0_dp/2.0_dp)) * real(kronecker_delta(l_a,l_b-1)*kronecker_delta(n_a,n_b),dp) + &
                sqrt(real(l_b*(n_b+1),dp)) * real(kronecker_delta(l_a,l_b-1)*kronecker_delta(n_a,n_b+1),dp))
    end function rGrad

end module nj_o1b
