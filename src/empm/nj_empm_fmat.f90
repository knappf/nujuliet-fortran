! Particle-hole (Pandya) transformed interaction F of the EMPM chain ... counterpart of fmat_unfi.f90 of the original EMPM codes ...
    ! F^L(i,j,k,l) = sum_J (2J+1) (-1)^(j_j+j_k-L-J) {j_i j_j L; j_l j_k J} V^J(ik,jl) ... V antisymmetrized & not normalized ...
    ! Channel (P,L) holds F as a symmetric matrix over the ordered same-isospin level pairs (i,j) of the restricted space
        ! with parity P of l_i + l_j & |j_i - j_j| <= L <= j_i + j_j ... pp, pn, np & nn blocks together ...
    ! Symmetries ... F^L(i,j,k,l) = F^L(k,l,i,j) for all isospin combinations ... F_np of the original is the transposed pn block ...
module nj_empm_fmat
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : f6j
    use nj_o2b, only : O2b_pp, O2b_pn, O2b_nn
    use nj_empm_space
    implicit none
    private

    ! Pair channels ... N(P,0:L_max) pairs, Ind(P,L)%v(1:2,Pair) levels of the pair, Dic(i,j,L) pair index (0 = none)
        ! & V(P,L)%v(Bra,Ket) the F-matrix elements ...
    type, public :: fmat_t
        integer :: L_max = 0
        integer, allocatable :: N(:,:)
        type(imat_t), allocatable :: Ind(:,:)
        integer, allocatable :: Dic(:,:,:)
        type(rmat_t), allocatable :: V(:,:)
    end type fmat_t

    public :: Fmat_initialize, Fmat_make, Fmat_ME, Fmat_export, Fmat_import

contains

    ! Zero F-matrix with the pair channels of the restricted space ...
    subroutine Fmat_initialize(Space,F_ph)
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(out) :: F_ph
        integer :: N_lev, i1, i2, Lam, P, Pass

        N_lev = Space%N_lev
        F_ph%L_max = maxval(Space%Lev(1:N_lev)%j)
        allocate(F_ph%N(2,0:F_ph%L_max), F_ph%Ind(2,0:F_ph%L_max), F_ph%V(2,0:F_ph%L_max))
        allocate(F_ph%Dic(N_lev,N_lev,0:F_ph%L_max))
        F_ph%Dic = 0

        ! Pass 1 counts the pairs, pass 2 fills the indices ...
        do Pass = 1, 2
            F_ph%N = 0
            do i1 = 1, N_lev
                do i2 = 1, N_lev
                    if (Space%Lev(i1)%t /= Space%Lev(i2)%t) cycle
                    P = mod(Space%Lev(i1)%l + Space%Lev(i2)%l,2) + 1
                    do Lam = abs(Space%Lev(i1)%j - Space%Lev(i2)%j)/2, (Space%Lev(i1)%j + Space%Lev(i2)%j)/2
                        F_ph%N(P,Lam) = F_ph%N(P,Lam) + 1
                        if (Pass == 2) then
                            F_ph%Dic(i1,i2,Lam) = F_ph%N(P,Lam)
                            F_ph%Ind(P,Lam)%v(1,F_ph%N(P,Lam)) = i1
                            F_ph%Ind(P,Lam)%v(2,F_ph%N(P,Lam)) = i2
                        end if
                    end do
                end do
            end do
            if (Pass == 1) then
                do Lam = 0, F_ph%L_max
                    do P = 1, 2
                        allocate(F_ph%Ind(P,Lam)%v(2,F_ph%N(P,Lam)), F_ph%V(P,Lam)%v(F_ph%N(P,Lam),F_ph%N(P,Lam)))
                        F_ph%V(P,Lam)%v = 0.0_dp
                    end do
                end do
            end if
        end do
    end subroutine Fmat_initialize

    ! F-matrix from the residual 2-body interaction in the HF basis ... lower triangle computed, upper one mirrored ...
    subroutine Fmat_make(Params,Orb,Orb_NN,V_NN,Space,F_ph)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o2b_t), intent(in) :: V_NN
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(out) :: F_ph
        integer :: Lam, P, N, Bra, Ket, e_max, i

        ! Every 2-body orbital of the restricted space has to be stored in V_NN ...
        e_max = 0
        do i = 1, Space%N_lev
            e_max = max(e_max, 2*Orb(Space%Lev(i)%a)%n + Orb(Space%Lev(i)%a)%l)
        end do
        if (2*e_max > Params%Calc%N2max) call fatal('EMPM restricted space needs N2max >= ' // istr(2*e_max))

        write(*,'(/,a)') 'Calculating the F-matrix (Pandya transformed interaction) in the restricted EMPM space ...'
        call Fmat_initialize(Space,F_ph)
        do Lam = 0, F_ph%L_max
            do P = 1, 2
                N = F_ph%N(P,Lam)
                !$omp parallel do schedule(dynamic) private(Ket)
                do Bra = 1, N
                    do Ket = 1, Bra
                        F_ph%V(P,Lam)%v(Bra,Ket) = F_value(F_ph%Ind(P,Lam)%v(1,Bra),F_ph%Ind(P,Lam)%v(2,Bra), &
                                                           F_ph%Ind(P,Lam)%v(1,Ket),F_ph%Ind(P,Lam)%v(2,Ket),Lam)
                        F_ph%V(P,Lam)%v(Ket,Bra) = F_ph%V(P,Lam)%v(Bra,Ket)
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        write(*,'(a)') TAB // 'F-matrix calculated ... L_max = ' // istr(F_ph%L_max) // ' , number of level pairs = ' // &
                       istr(sum(F_ph%N))

    contains

        real(dp) function F_value(i1,i2,i3,i4,Lam)
            integer, intent(in) :: i1, i2, i3, i4, Lam
            integer :: J, Q, j1, j2, j3, j4
            F_value = 0.0_dp
            if (mod(Space%Lev(i1)%l + Space%Lev(i3)%l,2) /= mod(Space%Lev(i2)%l + Space%Lev(i4)%l,2)) return
            Q = mod(Space%Lev(i1)%l + Space%Lev(i3)%l,2) + 1
            j1 = Space%Lev(i1)%j
            j2 = Space%Lev(i2)%j
            j3 = Space%Lev(i3)%j
            j4 = Space%Lev(i4)%j
            do J = max(abs(j1 - j3), abs(j2 - j4))/2, min(j1 + j3, j2 + j4)/2
                F_value = F_value + phase((j2 + j3)/2 - Lam - J) * real(2*J + 1,dp) * f6j(j1,j2,2*Lam,j4,j3,2*J) * &
                                    V_ME(i1,i3,i2,i4,J,Q)
            end do
        end function F_value

        ! V^J(ab,cd) ... a & c of one isospin, b & d of the other one ... np case by exchange of both pairs ...
        real(dp) function V_ME(i1,i2,i3,i4,J,Q)
            integer, intent(in) :: i1, i2, i3, i4, J, Q
            integer :: a, b, c, d
            a = Space%Lev(i1)%a
            b = Space%Lev(i2)%a
            c = Space%Lev(i3)%a
            d = Space%Lev(i4)%a
            if (Space%Lev(i1)%t == -1 .and. Space%Lev(i2)%t == -1) then
                V_ME = O2b_pp(a,b,c,d,J,Q,V_NN,Orb,Orb_NN)
            else if (Space%Lev(i1)%t == -1) then
                V_ME = O2b_pn(a,b,c,d,J,Q,V_NN,Orb_NN)
            else if (Space%Lev(i2)%t == -1) then
                V_ME = phase((Orb(a)%j + Orb(b)%j + Orb(c)%j + Orb(d)%j)/2) * O2b_pn(b,a,d,c,J,Q,V_NN,Orb_NN)
            else
                V_ME = O2b_nn(a,b,c,d,J,Q,V_NN,Orb,Orb_NN)
            end if
        end function V_ME

    end subroutine Fmat_make

    ! F^L(i,j,k,l) for levels i, j, k, l of the restricted space ... zero outside of the stored channels ...
    pure real(dp) function Fmat_ME(F_ph,Space,i1,i2,i3,i4,Lam)
        type(fmat_t), intent(in) :: F_ph
        type(empm_space_t), intent(in) :: Space
        integer, intent(in) :: i1, i2, i3, i4, Lam
        integer :: Bra, Ket, P
        Fmat_ME = 0.0_dp
        if (Lam > F_ph%L_max) return
        Bra = F_ph%Dic(i1,i2,Lam)
        Ket = F_ph%Dic(i3,i4,Lam)
        if (Bra == 0 .or. Ket == 0) return
        P = mod(Space%Lev(i1)%l + Space%Lev(i2)%l,2) + 1
        if (P /= mod(Space%Lev(i3)%l + Space%Lev(i4)%l,2) + 1) return
        Fmat_ME = F_ph%V(P,Lam)%v(Bra,Ket)
    end function Fmat_ME

    ! Raw Float64 stream of the channel matrices, channel by channel (L outer, P inner), empty channels skipped ...
        ! the pair structure is rebuilt from the restricted space by Fmat_initialize ...
    subroutine Fmat_export(F_ph,Export_Path)
        type(fmat_t), intent(in) :: F_ph
        character(len=*), intent(in) :: Export_Path
        integer :: u, Lam, P

        write(*,'(/,a)') 'Exporting the F-matrix in internal binary format ...'
        open(newunit=u, file=Export_Path, access='stream', form='unformatted', status='replace', action='write')
        do Lam = 0, F_ph%L_max
            do P = 1, 2
                if (F_ph%N(P,Lam) == 0) cycle
                write(u) F_ph%V(P,Lam)%v
            end do
        end do
        close(u)
        write(*,'(a)') TAB // "F-matrix succesfully exported ... ''" // Export_Path // "''"
    end subroutine Fmat_export

    subroutine Fmat_import(Space,Import_Path,F_ph)
        type(empm_space_t), intent(in) :: Space
        character(len=*), intent(in) :: Import_Path
        type(fmat_t), intent(out) :: F_ph
        integer :: u, Lam, P

        write(*,'(/,a)') 'Reading the F-matrix ...'
        call Fmat_initialize(Space,F_ph)
        open(newunit=u, file=Import_Path, access='stream', form='unformatted', status='old', action='read')
        do Lam = 0, F_ph%L_max
            do P = 1, 2
                if (F_ph%N(P,Lam) == 0) cycle
                read(u) F_ph%V(P,Lam)%v
            end do
        end do
        close(u)
        write(*,'(a)') TAB // 'F-matrix succesfully read ...'
    end subroutine Fmat_import

end module nj_empm_fmat
