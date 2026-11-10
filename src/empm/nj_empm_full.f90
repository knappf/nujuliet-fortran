! Couplings between the 0-, 1- & 2-phonon subspaces & the full EMPM Hamiltonian of one J^P channel ... counterpart of
    ! nondiag12 (vintn12), nondiag (vintn2) & fullham (nondiag1r_v2, nondiag2, fulhamt) of the original EMPM codes,
    ! with the redefined interaction computed on the fly ...
    ! For the full 2-phonon basis k = |alpha lambda; J> & the 1-phonon states a of J^P ...
        ! V(a;k)  = (-1)^(Ja+Jal+Jl)/(2Ja+1) sum_s rho_(a alpha)(Jl; s) W_l(s)
        ! W_l(s)  = sum_(ph of t) kappa c^l_(ph) F^Jl(p,h,i1,i2) , s = (i1,i2), kappa = 1/2 for t the isospin of s ...
        ! V0(k)   = (-1)^(Jal+Jl) sqrt(2Jl+1) d(Jal,Jl) [1/4 c^al_p F c^l_p + 1/4 c^al_n F c^l_n + c^al_p F c^l_n] (0+ only)
    ! The correlated 2-phonon states nu enter through X^nu = sqrt(2J+1) D^T C^nu, i.e. H_12(nu,a) = sqrt(2J+1) sum_i C^nu_i P(i,a)
        ! with P(i,a) = sum_k D_i(k) V(a;k) accumulated from the D rows of EMPM_eqm2_build ...
module nj_empm_full
    use nj_kinds
    use nj_types
    use nj_utils
    use nj_empm_space
    use nj_empm_fmat, only : fmat_t, Fmat_ME
    use nj_empm_phonons, only : empm_phonons_t
    use nj_empm_dens
    use nj_empm_eqm2, only : empm_basis2_t, empm_solution2_t, EMPM_syevd
    use nj_linalg, only : mtm
    implicit none
    private

    interface
        subroutine dgemm(Transa,Transb,m,n,k,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
            import :: dp
            character :: Transa, Transb
            integer :: m, n, k, Lda, Ldb, Ldc
            real(dp) :: Alpha, Beta, A(Lda,*), B(Ldb,*), C(Ldc,*)
        end subroutine dgemm
    end interface

    ! Full Hamiltonian of a channel ... N0 (0 or 1), N1 1-phonon states A1 (phonon indices, spurious one excluded),
        ! N2 physical 2-phonon states, eigenvalues E & eigenvectors Z(N,N) with the basis order [0 ; 1 ; 2] ...
    type, public :: empm_full_t
        integer :: N0 = 0, N1 = 0, N2 = 0, N = 0
        integer, allocatable :: A1(:)
        real(dp), allocatable :: E(:), Z(:,:)
    end type empm_full_t

    public :: EMPM_coupling_columns, EMPM_coupling_V, EMPM_full_solve, EMPM_phonon_W, EMPM_coupling_Vn

contains

    ! 1-phonon states of J^P (without the spurious phonon) & whether the ground state couples (0+) ...
    subroutine EMPM_coupling_columns(Ph,J,P,A1,N0)
        type(empm_phonons_t), intent(in) :: Ph
        integer, intent(in) :: J, P
        integer, allocatable, intent(out) :: A1(:)
        integer, intent(out) :: N0
        integer :: i
        A1 = pack([(i, i = 1, Ph%N)], Ph%J == J .and. Ph%P == P .and. [(i, i = 1, Ph%N)] /= Ph%i_cm)
        N0 = merge(1, 0, J == 0 .and. P == 1)
    end subroutine EMPM_coupling_columns

    ! Vk(k,c) over the full basis ... columns c = 1 .. size(A1) the 1-phonon states, column size(A1)+1 the ground state (N0 = 1) ...
    subroutine EMPM_coupling_V(Space,F_ph,Pairs,Ph,Bs,A1,N0,Vk)
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_basis2_t), intent(in) :: Bs
        integer, intent(in) :: A1(:), N0
        real(dp), allocatable, intent(out) :: Vk(:,:)
        real(dp), allocatable :: W(:,:), G(:,:)
        type(phonon_block_t) :: Rho_a
        real(dp), allocatable :: Wc(:,:)
        integer :: c, a, Ja, Jb, Pb, Jc, Pc, nb, nc, ib, ic, b, l, k, Pp, Ka
        real(dp) :: Fac

        allocate(Vk(Bs%N_full,size(A1)+N0))
        Vk = 0.0_dp
        call EMPM_phonon_W(Space,F_ph,Pairs,Ph,W)

        ! 1-phonon columns ... in parallel, each column c written by one thread ...
        !$omp parallel do schedule(dynamic) private(a,Ja,Rho_a,Jb,Pb,nb,Jc,Pc,nc,G,Wc,Fac,ic,l,ib,b,k,Pp,Ka)
        do c = 1, size(A1)
            a = A1(c)
            Ja = Ph%J(a)
            call EMPM_density(Pairs,Ph,a,Rho_a)
            do Jb = 0, Ph%J_max
                do Pb = 1, 2
                    nb = Ph%N_ch(Jb,Pb)
                    if (nb == 0) cycle
                    Pp = Rho_a%Pi(Jb,Pb)
                    do Jc = abs(Bs%J - Jb), min(Bs%J + Jb, Ph%J_max)
                        if (Jc < Rho_a%Lo(Jb,Pb) .or. Jc > Rho_a%Hi(Jb,Pb)) cycle
                        do Pc = 1, 2
                            nc = Ph%N_ch(Jc,Pc)
                            if (nc == 0) cycle
                            if (merge(1,-1,Pb == 1)*merge(1,-1,Pc == 1) /= merge(1,-1,Bs%P == 1)) cycle
                            Ka = Pairs%Na(Jc,Pp)
                            if (Ka == 0) cycle
                            allocate(G(nb,nc))
                            Wc = W(Pairs%Al(Jc,Pp)%v,Ph%First(Jc,Pc):Ph%First(Jc,Pc)+nc-1)
                            call dgemm('T','N',nb,nc,Ka,1.0_dp,Rho_a%B(Jb,Pb)%M(Jc)%v,Ka,Wc,Ka,0.0_dp,G,nb)
                            Fac = phase(Ja + Jb + Jc) / real(2*Ja + 1,dp)
                            do ic = 1, nc
                                l = Ph%First(Jc,Pc) + ic - 1
                                do ib = 1, nb
                                    b = Ph%First(Jb,Pb) + ib - 1
                                    k = Bs%K(b,l)
                                    if (k > 0) Vk(k,c) = Fac * G(ib,ic)
                                end do
                            end do
                            deallocate(G)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! Ground-state column (0+) ...
        if (N0 == 1) then
            !$omp parallel do schedule(dynamic) private(b,l)
            do k = 1, Bs%N_full
                b = Bs%Al(k)
                l = Bs%La(k)
                if (Ph%J(b) /= Ph%J(l)) cycle
                Vk(k,size(A1)+1) = phase(Ph%J(b) + Ph%J(l)) * sqrt(real(2*Ph%J(l) + 1,dp)) * V0(b,l)
            end do
            !$omp end parallel do
        end if

    contains

        ! Ground-state redefined interaction of vintn2 ... pp & nn with 1/4, pn once (proton ph of g, neutron ph of l) ...
        real(dp) function V0(g,l)
            integer, intent(in) :: g, l
            integer :: t, u, p, h, q, r
            real(dp) :: Kap
            V0 = 0.0_dp
            do t = 1, 2
                do u = t, 2
                    Kap = merge(0.25_dp, 1.0_dp, t == u)
                    do h = 1, Pairs%Nk(t+2)
                        do p = 1, Pairs%Nk(t)
                            if (amp(t,p,h,g) == 0.0_dp) cycle
                            do r = 1, Pairs%Nk(u+2)
                                do q = 1, Pairs%Nk(u)
                                    if (amp(u,q,r,l) == 0.0_dp) cycle
                                    V0 = V0 + Kap * amp(t,p,h,g) * amp(u,q,r,l) * &
                                              Fmat_ME(F_ph,Space,Pairs%Pos(t)%v(p),Pairs%Pos(t+2)%v(h),Pairs%Pos(u)%v(q), &
                                                      Pairs%Pos(u+2)%v(r),Ph%J(l))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end function V0

        pure real(dp) function amp(t,p,h,i)
            integer, intent(in) :: t, p, h, i
            if (t == 1) then
                amp = Ph%Cp(p,h,i)
            else
                amp = Ph%Cn(p,h,i)
            end if
        end function amp

    end subroutine EMPM_coupling_V

    ! W(s,l) = sum_(p h, t) kappa c^l_(p h) F^(J_l)(p,h; s) for all phonons l & pairs s ... the phonon side of the couplings
        ! between the n- & (n+1)-phonon subspaces (kappa = 1/2 for the isospin of the pair s) ...
    subroutine EMPM_phonon_W(Space,F_ph,Pairs,Ph,W)
        type(empm_space_t), intent(in) :: Space
        type(fmat_t), intent(in) :: F_ph
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        real(dp), allocatable, intent(out) :: W(:,:)
        integer :: l, k, i1, i2, s, t, p, h
        real(dp) :: Sum, Kap
        allocate(W(Pairs%S,Ph%N))
        !$omp parallel do schedule(dynamic) private(k,i1,i2,s,t,p,h,Sum,Kap)
        do l = 1, Ph%N
            do k = 1, 4
                do i2 = 1, Pairs%Nk(k)
                    do i1 = 1, Pairs%Nk(k)
                        s = Pairs%Off(k) + i1 + (i2 - 1)*Pairs%Nk(k)
                        Sum = 0.0_dp
                        do t = 1, 2
                            Kap = merge(0.5_dp, 1.0_dp, t == 2 - mod(k,2))
                            do h = 1, Pairs%Nk(t+2)
                                do p = 1, Pairs%Nk(t)
                                    if (amp(t,p,h,l) == 0.0_dp) cycle
                                    Sum = Sum + Kap * amp(t,p,h,l) * Fmat_ME(F_ph,Space,Pairs%Pos(t)%v(p),Pairs%Pos(t+2)%v(h), &
                                                                              Pairs%Pos(k)%v(i1),Pairs%Pos(k)%v(i2),Ph%J(l))
                                end do
                            end do
                        end do
                        W(s,l) = Sum
                    end do
                end do
            end do
        end do
        !$omp end parallel do

    contains

        pure real(dp) function amp(t,p,h,i)
            integer, intent(in) :: t, p, h, i
            if (t == 1) then
                amp = Ph%Cp(p,h,i)
            else
                amp = Ph%Cn(p,h,i)
            end if
        end function amp

    end subroutine EMPM_phonon_W

    ! Coupling of the n-phonon state a to the (n+1)-phonon basis B |alpha lambda> ... V(k) = (-1)^(Ja+J_alpha+J_lambda)/(2Ja+1)
        ! sum_s rho_(a alpha)(J_lambda; s) W(s,lambda), Rho_a the n-phonon densities of a (blocks over the channels of the
        ! n-phonon states S), as EMPM_coupling_V (n = 1) & vintn23 (n = 2) ... serial ...
    subroutine EMPM_coupling_Vn(Pairs,Ph,W,Rho_a,Ja,S_J_max,S_N_ch,S_First,B,V)
        type(empm_pairs_t), intent(in) :: Pairs
        type(empm_phonons_t), intent(in) :: Ph
        real(dp), intent(in) :: W(:,:)
        type(phonon_block_t), intent(in) :: Rho_a
        integer, intent(in) :: Ja, S_J_max, S_N_ch(0:,:), S_First(0:,:)
        type(empm_basis2_t), intent(in) :: B
        real(dp), intent(out) :: V(:)
        real(dp), allocatable :: G(:,:), Wc(:,:)
        integer :: Jb, Pb, Jc, Pc, nb, nc, ib, ic, b2, l, k, Pp, Ka
        real(dp) :: Fac

        V = 0.0_dp
        do Jb = 0, S_J_max
            do Pb = 1, 2
                nb = S_N_ch(Jb,Pb)
                if (nb == 0) cycle
                Pp = Rho_a%Pi(Jb,Pb)
                do Jc = abs(B%J - Jb), min(B%J + Jb, Ph%J_max)
                    if (Jc < Rho_a%Lo(Jb,Pb) .or. Jc > Rho_a%Hi(Jb,Pb)) cycle
                    do Pc = 1, 2
                        nc = Ph%N_ch(Jc,Pc)
                        if (nc == 0) cycle
                        if (merge(1,-1,Pb == 1)*merge(1,-1,Pc == 1) /= merge(1,-1,B%P == 1)) cycle
                        Ka = Pairs%Na(Jc,Pp)
                        if (Ka == 0) cycle
                        allocate(G(nb,nc))
                        Wc = W(Pairs%Al(Jc,Pp)%v,Ph%First(Jc,Pc):Ph%First(Jc,Pc)+nc-1)
                        call dgemm('T','N',nb,nc,Ka,1.0_dp,Rho_a%B(Jb,Pb)%M(Jc)%v,Ka,Wc,Ka,0.0_dp,G,nb)
                        Fac = phase(Ja + Jb + Jc) / real(2*Ja + 1,dp)
                        do ic = 1, nc
                            l = Ph%First(Jc,Pc) + ic - 1
                            do ib = 1, nb
                                b2 = S_First(Jb,Pb) + ib - 1
                                k = B%K(b2,l)
                                if (k > 0) V(k) = Fac * G(ib,ic)
                            end do
                        end do
                        deallocate(G)
                    end do
                end do
            end do
        end do
    end subroutine EMPM_coupling_Vn

    ! Full Hamiltonian [0 ; 1 ; 2] of the channel & its diagonalization ... physical 2-phonon solutions only (spurious excluded),
        ! Pk(i,c) = sum_k D_i(k) Vk(k,c) for the truncated rows i ...
    subroutine EMPM_full_solve(Ph,Bs,Sol,A1,N0,Pk,Full)
        type(empm_phonons_t), intent(in) :: Ph
        type(empm_basis2_t), intent(in) :: Bs
        type(empm_solution2_t), intent(in) :: Sol
        integer, intent(in) :: A1(:), N0
        real(dp), intent(in) :: Pk(:,:)
        type(empm_full_t), intent(out) :: Full
        real(dp), allocatable :: H12(:,:), Hf(:,:)
        integer :: N, n1, n2, i

        n1 = size(A1)
        n2 = Sol%N_ind - Sol%N_s
        Full%N0 = N0
        Full%N1 = n1
        Full%N2 = n2
        Full%A1 = A1
        N = N0 + n1 + n2
        Full%N = N

        ! H_12(nu,c) = sqrt(2J+1) sum_(independent rows) C(:,nu) P(Ind,c) ...
        H12 = sqrt(real(2*Bs%J + 1,dp)) * mtm(Sol%C(:,1:n2),Pk(Sol%Ind,:))

        allocate(Hf(N,N))
        Hf = 0.0_dp
        do i = 1, n1
            Hf(N0+i,N0+i) = Ph%E(A1(i))
        end do
        do i = 1, n2
            Hf(N0+n1+i,N0+n1+i) = Sol%E(i)
        end do
        Hf(N0+n1+1:N,N0+1:N0+n1) = H12(:,1:n1)
        Hf(N0+1:N0+n1,N0+n1+1:N) = transpose(H12(:,1:n1))
        if (N0 == 1) then
            Hf(N0+n1+1:N,1) = H12(:,n1+1)
            Hf(1,N0+n1+1:N) = H12(:,n1+1)
        end if

        call EMPM_syevd(Hf,Full%E)
        call move_alloc(Hf,Full%Z)
    end subroutine EMPM_full_solve

end module nj_empm_full
