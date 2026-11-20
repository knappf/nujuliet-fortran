! Stage-4 unit check of the HF-ScRPA building blocks against the Julia dumps of scrpa_blocks.jl (see run_unit.sh) ...
    ! The Julia initial RPA amplitudes (dump/XY.bin) are the common input, so no eigenvector gauge enters ...
program scrpa_blocks
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_input
    use nj_wigner, only : wigner_init_float
    use nj_orbitals, only : orbitals_import, orbitals_1p1h_make, OBDM_natural_orbitals, a_max_of
    use nj_o1b, only : O1b_import, O1b_transformation, O1b_unit_matrix, T1b
    use nj_o2b, only : O2b_import
    use nj_v2b, only : V2b_read
    use nj_v3b, only : V3b_no2b_read
    use nj_tr1b
    use nj_rpa, only : RPA_allocate, RPA_diagonalize, RPA_energy_QBA, RPA_energy, RPA_sigma, rpa_energy_t, RPA_sigma_mix
    use nj_rpa_rho
    use nj_rpa_hamiltonian
    implicit none
    type(run_input_t) :: Run
    type(params_t) :: Params, Pr
    type(orb1b_t), allocatable :: Orb_ref(:), Orb(:), Orb_tmp(:)
    type(orb2b_1p1h_t) :: Orb_1p1h
    type(orb1b_rpa_obdm_t) :: Orb_ph
    type(orb2b_t) :: Orb_NN
    type(orb3b_t) :: Orb_NNN
    type(o1b_t) :: T, C_LHO_HF, Rho_ref, Rho1, Cn, Rd, Rho, C_HF_NAT, C_LHO_NAT, hN, Jref, Rho_ref_nat, Ref(3)
    type(o2b_t) :: V_NN_bare, Sigma_in, Sigma0, VNN, Vref, Sig1, Smix
    type(v3b_t) :: V_NNN_bare
    type(cmat_t), allocatable :: X(:,:), Y(:,:), X2(:,:), Y2(:,:), Xc(:,:), Yc(:,:), X2c(:,:), Y2c(:,:), X1(:,:), Y1(:,:)
    type(cvec_t), allocatable :: E(:,:), E1(:,:)
    type(rmat_t), allocatable :: A(:,:), B(:,:), N(:,:)
    type(tr1b_t) :: TrOp_LHO, TrOp1
    type(rpa_energy_t) :: Ec
    logical :: Stab
    real(dp) :: Ej(7), d, E_QBA0, E_QBA1, worst
    integer :: u, J, P, Np, a_max, k
    character(len=16) :: Tags(3) = [character(len=16) :: 'ScH', 'ScH_TBCF', 'HFref']
    character(len=*), parameter :: IO = 'io/He4/HF'

    worst = 0.0_dp
    call read_input('../../fortran/bundled/08_He4_ScRPA.nml',Run)
    Params = Run%Params
    call wigner_init_float(75,'Jmax',9)
    a_max = a_max_of(Params%Calc%Nmax)

    Orb_ref = orbitals_import(Params,IO // '/bin/orb.bin')
    T = T1b(Params,Orb_ref)
    call V2b_read(Params,V_NN_bare,Orb_NN)
    call V3b_no2b_read(Params,V_NNN_bare,Orb_NNN)
    C_LHO_HF = O1b_import(Params,Orb_ref,IO // '/bin/C.bin')
    Rho_ref = O1b_import(Params,Orb_ref,IO // '/bin/rho.bin')
    Orb_1p1h = orbitals_1p1h_make(Params,Orb_ref,Params%Calc%RPA%Tol_par,Params%Calc%RPA%Tol_hol)
    Orb_ph = RPA_OBDM_orbitals_make(Params,Orb_ref)
    TrOp_LHO = Tr1b_initialize(Params,Orb_ref,1.0_dp)
    write(*,'(a,4i4)') 'Orb_ph N_p/N_h (p,n): ', Orb_ph%N_p%p, Orb_ph%N_h%p, Orb_ph%N_p%n, Orb_ph%N_h%n

    ! Julia amplitudes ...
    allocate(X(2,0:Orb_1p1h%J_max), Y(2,0:Orb_1p1h%J_max), E(2,0:Orb_1p1h%J_max))
    open(newunit=u, file='dump/XY.bin', access='stream', form='unformatted', status='old')
    do J = 0, Orb_1p1h%J_max
        do P = 1, 2
            block
                integer(i8) :: n8
                real(dp), allocatable :: R1(:,:), R2(:,:), e1v(:), e2v(:)
                read(u) n8
                Np = int(n8)
                if (Np /= Orb_1p1h%N(P,J)) stop 'XY.bin: block size mismatch'
                allocate(R1(Np,Np), R2(Np,Np), e1v(Np), e2v(Np))
                if (Np > 0) then
                    read(u) R1, R2
                    X(P,J)%v = cmplx(R1,R2,kind=dp)
                    read(u) R1, R2
                    Y(P,J)%v = cmplx(R1,R2,kind=dp)
                    read(u) e1v, e2v
                    E(P,J)%v = cmplx(e1v,e2v,kind=dp)
                else
                    allocate(X(P,J)%v(0,0), Y(P,J)%v(0,0), E(P,J)%v(0))
                end if
            end block
        end do
    end do
    close(u)

    ! oY4, natural orbitals, transformation ...
    Rho1 = RPA_OBDM_oY4(Params,Orb_ref,Orb_ph,Orb_1p1h,Rho_ref,X,Y)
    call cmp1('oY4',Rho1,read1('dump/oY4.bin',1))
    call OBDM_natural_orbitals(Params,Orb_ref,Rho1,Orb_tmp,Cn,Rd)
    call cmp1('OBDM_natural_orbitals C',Cn,read1('dump/nat.bin',1))
    call RPA_OBDM_transform_amplitudes(Orb_1p1h,X,Y,Cn,X2,Y2)
    call cmp_amp('transform_amplitudes',X2,Y2)

    ! complex path: column phases d_nu in {i, -i} (D^2 = -1) leave dev's M = -1/2 X^T W X^* (transpose, not adjoint:
    ! M -> D M D^*, G -> Y Y^+ + Y D^2 M D^*2 Y^+) & thus oY4 invariant; the transformation commutes with any phase ...
    allocate(Xc(2,0:Orb_1p1h%J_max), Yc(2,0:Orb_1p1h%J_max))
    do J = 0, Orb_1p1h%J_max
        do P = 1, 2
            Xc(P,J)%v = X(P,J)%v
            Yc(P,J)%v = Y(P,J)%v
            do k = 1, size(X(P,J)%v,2)
                Xc(P,J)%v(:,k) = X(P,J)%v(:,k) * ph4(k + J)
                Yc(P,J)%v(:,k) = Y(P,J)%v(:,k) * ph4(k + J)
            end do
        end do
    end do
    call cmp1('oY4 complex path vs real path',RPA_OBDM_oY4(Params,Orb_ref,Orb_ph,Orb_1p1h,Rho_ref,Xc,Yc),Rho1)
    ! generic phases: the complex path against dev's (phase-dependent) formula ...
    do J = 0, Orb_1p1h%J_max
        do P = 1, 2
            do k = 1, size(X(P,J)%v,2)
                Xc(P,J)%v(:,k) = X(P,J)%v(:,k) * exp(cmplx(0.0_dp,0.37_dp*k + J,kind=dp))
                Yc(P,J)%v(:,k) = Y(P,J)%v(:,k) * exp(cmplx(0.0_dp,0.37_dp*k + J,kind=dp))
            end do
        end do
    end do
    call cmp1('oY4 complex path vs Julia (generic phases)',RPA_OBDM_oY4(Params,Orb_ref,Orb_ph,Orb_1p1h,Rho_ref,Xc,Yc), &
              read1('dump/oY4c.bin',1))
    call RPA_OBDM_transform_amplitudes(Orb_1p1h,Xc,Yc,Cn,X2c,Y2c)
    d = 0.0_dp
    do J = 0, Orb_1p1h%J_max
        do P = 1, 2
            do k = 1, size(X(P,J)%v,2)
                d = max(d, maxval(abs(X2c(P,J)%v(:,k) - X2(P,J)%v(:,k) * exp(cmplx(0.0_dp,0.37_dp*k + J,kind=dp)))), &
                           maxval(abs(Y2c(P,J)%v(:,k) - Y2(P,J)%v(:,k) * exp(cmplx(0.0_dp,0.37_dp*k + J,kind=dp)))))
            end do
        end do
    end do
    call report('transform complex path',d)

    ! HF_ScRPA_OBDM_solve (default, OBDM_diagonal, alpha = 1) ...
    call solve_case('',Params)
    Pr = Params
    Pr%Calc%RPA%OBDM_diagonal = .true.
    call solve_case('_diag',Pr)
    Pr = Params
    Pr%Calc%RPA%OBDM_alpha = 1.0_dp
    call solve_case('_a1',Pr)
    ! default solution kept for the following tests ...
    Orb = Orb_ref
    Rho = Rho_ref
    C_HF_NAT = O1b_unit_matrix(Params)
    call HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT,C_LHO_NAT)

    ! RPA_sigma (F2 in both codes) & RPA_Hamiltonian variants ...
    call O2b_import(Params,Orb,'dump/Sigma_in.bin',Sigma_in)
    call cmp2('RPA_sigma',RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y),Sigma_in)
    Smix = Sigma_in
    call RPA_sigma_mix(Smix,Sigma0_like(Sigma_in),0.25_dp)
    call cmp2('RPA_sigma_mix (0.25*S + 0.75*0)',Smix,scaled(Sigma_in,0.25_dp))
    do k = 1, 3
        Pr = Params
        if (k == 2) Pr%Calc%RPA%TBCF = .true.
        if (k == 3) Pr%Calc%RPA%ScH = .false.
        call RPA_Hamiltonian(Pr,Orb,Orb_NN,Orb_NNN,C_LHO_NAT,Rho,Sigma_in,T,V_NN_bare,V_NNN_bare,hN,VNN)
        Jref = O1b_import(Pr,Orb,'dump/H1b_' // trim(Tags(k)) // '.bin')
        call cmp1lj('RPA_Hamiltonian h_N ' // trim(Tags(k)),hN,Jref)
        call O2b_import(Pr,Orb,'dump/H2b_' // trim(Tags(k)) // '.bin',Vref)
        call cmp2('RPA_Hamiltonian V_NN ' // trim(Tags(k)),VNN,Vref)
    end do

    ! One ScRPA step & the energies ...
    call RPA_Hamiltonian(Params,Orb,Orb_NN,Orb_NNN,C_LHO_NAT,Rho,Sigma0,T,V_NN_bare,V_NNN_bare,hN,VNN)
    TrOp1 = Tr1b_transformation(Params,Orb,C_LHO_NAT,TrOp_LHO)
    call RPA_allocate(Params,Orb,Orb_1p1h,Orb_NN,hN,VNN,Rho,Sigma0,.false.,A,B,N)
    call RPA_diagonalize(Params,Orb,Orb_1p1h,A,B,N,TrOp1,E1,X1,Y1,Stab)
    E_QBA0 = RPA_energy_QBA(Params,Orb_ref,Orb_1p1h,E,X,Y)
    E_QBA1 = RPA_energy_QBA(Params,Orb,Orb_1p1h,E1,X1,Y1)
    Sig1 = RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X1,Y1)
    Rho_ref_nat = O1b_transformation(Rho_ref,C_HF_NAT)
    Ec = RPA_energy(Params,Orb,Orb_1p1h,Orb_NN,E1,X1,Y1,Rho,Rho_ref_nat,hN,Sig1,VNN,.false.,Sigma0,Density=.true.)
    open(newunit=u, file='dump/energy.txt', status='old')
    read(u,*) Ej
    close(u)
    call report('E_QBA initial (rel)',abs(E_QBA0 - Ej(1))/abs(Ej(1)))
    call report('E_QBA step 1 (rel)',abs(E_QBA1 - Ej(2))/abs(Ej(2)))
    call report('RPA_energy E_QBA (rel)',abs(Ec%E_QBA - Ej(3))/abs(Ej(3)))
    call report('RPA_energy E_sigma (rel)',abs(Ec%E_sigma - Ej(4))/abs(Ej(4)))
    call report('RPA_energy E_rho (rel)',abs(Ec%E_rho - Ej(5))/abs(Ej(5)))
    call report('RPA_energy E_dens (rel)',abs(Ec%E_dens - Ej(7))/abs(Ej(7)))
    write(*,'(a,es10.3)') 'UNIT_S4 worst deviation = ', worst

contains

    complex(dp) function ph4(m)
        integer, intent(in) :: m
        ph4 = merge((0.0_dp,1.0_dp), (0.0_dp,-1.0_dp), mod(m,2) == 0)
    end function ph4

    subroutine report(Label,Dev)
        character(len=*), intent(in) :: Label
        real(dp), intent(in) :: Dev
        write(*,'(a,t50,es10.3)') 'CHECK ' // Label, Dev
        worst = max(worst,Dev)
    end subroutine report

    function read1(File,i_rec) result(O)
        character(len=*), intent(in) :: File
        integer, intent(in) :: i_rec
        type(o1b_t) :: O
        integer :: uu, r
        allocate(O%p(a_max,a_max), O%n(a_max,a_max))
        open(newunit=uu, file=File, access='stream', form='unformatted', status='old')
        do r = 1, i_rec
            read(uu) O%p, O%n
        end do
        close(uu)
    end function read1

    subroutine cmp1(Label,O,R)
        character(len=*), intent(in) :: Label
        type(o1b_t), intent(in) :: O, R
        call report(Label,max(maxval(abs(O%p - R%p)),maxval(abs(O%n - R%n))))
    end subroutine cmp1

    ! only the (l,j)-equal entries (the O1b file layout) ...
    subroutine cmp1lj(Label,O,R)
        character(len=*), intent(in) :: Label
        type(o1b_t), intent(in) :: O, R
        integer :: a, b
        real(dp) :: dd
        dd = 0.0_dp
        do a = 1, a_max
            do b = 1, a_max
                if (Orb_ref(a)%l /= Orb_ref(b)%l .or. Orb_ref(a)%j /= Orb_ref(b)%j) cycle
                dd = max(dd,abs(O%p(a,b) - R%p(a,b)),abs(O%n(a,b) - R%n(a,b)))
            end do
        end do
        call report(Label,dd)
    end subroutine cmp1lj

    subroutine cmp2(Label,O,R)
        character(len=*), intent(in) :: Label
        type(o2b_t), intent(in) :: O, R
        integer :: i, jj
        real(dp) :: dd
        dd = 0.0_dp
        do jj = lbound(O%pp,2), ubound(O%pp,2)
            do i = 1, 2
                if (allocated(O%pp(i,jj)%v)) dd = max(dd,maxval(abs(O%pp(i,jj)%v - R%pp(i,jj)%v),dim=1,mask=.true.))
                if (allocated(O%pn(i,jj)%v)) dd = max(dd,maxval(abs(O%pn(i,jj)%v - R%pn(i,jj)%v),dim=1,mask=.true.))
                if (allocated(O%nn(i,jj)%v)) dd = max(dd,maxval(abs(O%nn(i,jj)%v - R%nn(i,jj)%v),dim=1,mask=.true.))
            end do
        end do
        call report(Label,dd)
    end subroutine cmp2

    function Sigma0_like(S) result(Z)
        type(o2b_t), intent(in) :: S
        type(o2b_t) :: Z
        Z = scaled(S,0.0_dp)
    end function Sigma0_like

    function scaled(S,f) result(Z)
        type(o2b_t), intent(in) :: S
        real(dp), intent(in) :: f
        type(o2b_t) :: Z
        integer :: i, jj
        Z = S
        do jj = lbound(Z%pp,2), ubound(Z%pp,2)
            do i = 1, 2
                if (allocated(Z%pp(i,jj)%v)) Z%pp(i,jj)%v = f * Z%pp(i,jj)%v
                if (allocated(Z%pn(i,jj)%v)) Z%pn(i,jj)%v = f * Z%pn(i,jj)%v
                if (allocated(Z%nn(i,jj)%v)) Z%nn(i,jj)%v = f * Z%nn(i,jj)%v
            end do
        end do
    end function scaled

    subroutine cmp_amp(Label,Xa,Ya)
        character(len=*), intent(in) :: Label
        type(cmat_t), intent(in) :: Xa(:,0:), Ya(:,0:)
        integer :: uu, JJ, PP, n
        real(dp), allocatable :: R1(:,:), R2(:,:)
        real(dp) :: dd
        dd = 0.0_dp
        open(newunit=uu, file='dump/tr.bin', access='stream', form='unformatted', status='old')
        do JJ = 0, Orb_1p1h%J_max
            do PP = 1, 2
                n = Orb_1p1h%N(PP,JJ)
                if (n == 0) cycle
                allocate(R1(n,n), R2(n,n))
                read(uu) R1, R2
                dd = max(dd,maxval(abs(Xa(PP,JJ)%v - cmplx(R1,R2,kind=dp))))
                read(uu) R1, R2
                dd = max(dd,maxval(abs(Ya(PP,JJ)%v - cmplx(R1,R2,kind=dp))))
                deallocate(R1, R2)
            end do
        end do
        close(uu)
        call report(Label,dd)
    end subroutine cmp_amp

    subroutine solve_case(Tag,Ps)
        character(len=*), intent(in) :: Tag
        type(params_t), intent(in) :: Ps
        type(orb1b_t), allocatable :: Ob(:)
        type(o1b_t) :: Rh, CHN, CLN
        integer :: uu, a
        real(dp) :: dd
        Ob = Orb_ref
        Rh = Rho_ref
        CHN = O1b_unit_matrix(Ps)
        call HF_ScRPA_OBDM_solve(Ps,Ob,Orb_ph,Orb_1p1h,X,Y,Rh,CHN,CLN)
        do k = 1, 3
            allocate(Ref(k)%p(a_max,a_max), Ref(k)%n(a_max,a_max))
        end do
        open(newunit=uu, file='dump/solve' // Tag // '.bin', access='stream', form='unformatted', status='old')
        read(uu) Ref(1)%p, Ref(1)%n, Ref(2)%p, Ref(2)%n, Ref(3)%p, Ref(3)%n
        close(uu)
        call cmp1('OBDM_solve' // Tag // ' Rho',Rh,Ref(1))
        call cmp1('OBDM_solve' // Tag // ' C_LHO_NAT',CLN,Ref(2))
        call cmp1('OBDM_solve' // Tag // ' C_HF_NAT',CHN,Ref(3))
        dd = 0.0_dp
        do a = 1, a_max
            dd = max(dd,abs(Ob(a)%pO - Ref(1)%p(a,a)),abs(Ob(a)%nO - Ref(1)%n(a,a)))
        end do
        call report('OBDM_solve' // Tag // ' Orb occupations',dd)
        do k = 1, 3
            deallocate(Ref(k)%p, Ref(k)%n)
        end do
    end subroutine solve_case

end program scrpa_blocks
