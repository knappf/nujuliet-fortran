! Post-processing of stored HF-TDA & HF-RPA solutions ... radial transition densities & transition current maps ...
    ! Counterpart of HF_RPA_import_binary, HF_RPA_transition_densities & HF_RPA_transition_currents (HF_RPA_External.jl) ...
    ! The same routines serve HF-RRPA (Kind = 'RRPA') ... in the original these are near-identical copies ...
module nj_hf_rpa_external
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_wigner, only : wigner_init_float, fCG
    use nj_functions, only : Psi_rad_LHO, differentiate_Oh6, J1b_convective, J1b_spin
    use nj_orbitals, only : orbitals_make, orbitals_ph_make, orbitals_one_phonon_make
    use nj_o1b, only : O1b_import
    use nj_tr1b
    use nj_hf_rpa, only : HF_RPA_phonon_count, ph_orbitals
    implicit none
    private

    public :: HF_RPA_import_binary, parse_JP, phonon_transition_densities, phonon_transition_currents

contains

    ! Reads TDA (optional) & RPA energies & amplitudes written by the binary export ...
        ! Prefix = 'RPA' or 'RRPA' ... TDA solutions are read only when Read_TDA = .true. ...
    subroutine HF_RPA_import_binary(Params,N_nu,Prefix,Read_TDA,X_TDA,E_TDA,X_RPA,Y_RPA,E_RPA)
        type(params_t), intent(in) :: Params
        integer, intent(in) :: N_nu(0:,:)
        character(len=*), intent(in) :: Prefix
        logical, intent(in) :: Read_TDA
        type(rmat_t), allocatable, intent(out) :: X_TDA(:,:)
        type(rvec_t), allocatable, intent(out) :: E_TDA(:,:)
        type(cmat_t), allocatable, intent(out) :: X_RPA(:,:), Y_RPA(:,:)
        type(cvec_t), allocatable, intent(out) :: E_RPA(:,:)
        integer :: J_max, J, P, N_ph, u1, u2, u3, u4, u5

        J_max = 2*Params%Calc%Nmax + 1
        write(*,'(/,a)') 'Importing TDA & ' // Prefix // ' solutions from binary format ...'
        allocate(X_TDA(0:J_max,2), E_TDA(0:J_max,2), X_RPA(0:J_max,2), Y_RPA(0:J_max,2), E_RPA(0:J_max,2))
        if (Read_TDA) then
            open(newunit=u1, file=io_path(Params,'/Bin/TDA_X.bin'), access='stream', form='unformatted', status='old', action='read')
            open(newunit=u2, file=io_path(Params,'/Bin/TDA_E.bin'), access='stream', form='unformatted', status='old', action='read')
        end if
        open(newunit=u3, file=io_path(Params,'/Bin/' // Prefix // '_X.bin'), access='stream', form='unformatted', status='old', action='read')
        open(newunit=u4, file=io_path(Params,'/Bin/' // Prefix // '_Y.bin'), access='stream', form='unformatted', status='old', action='read')
        open(newunit=u5, file=io_path(Params,'/Bin/' // Prefix // '_E.bin'), access='stream', form='unformatted', status='old', action='read')
        do J = 0, J_max
            do P = 1, 2
                N_ph = N_nu(J,P)
                allocate(X_TDA(J,P)%v(N_ph,N_ph), E_TDA(J,P)%v(N_ph), X_RPA(J,P)%v(N_ph,N_ph), Y_RPA(J,P)%v(N_ph,N_ph), &
                         E_RPA(J,P)%v(N_ph))
                if (N_ph == 0) cycle
                if (Read_TDA) then
                    read(u1) X_TDA(J,P)%v
                    read(u2) E_TDA(J,P)%v
                end if
                read(u3) X_RPA(J,P)%v
                read(u4) Y_RPA(J,P)%v
                read(u5) E_RPA(J,P)%v
            end do
        end do
        if (Read_TDA) then
            close(u1)
            close(u2)
        end if
        close(u3)
        close(u4)
        close(u5)
        write(*,'(a)') TAB // 'TDA & ' // Prefix // ' solutions have been succesufully imported ...'
    end subroutine HF_RPA_import_binary

    ! "1-" -> J = 1, P = 2 ...
    subroutine parse_JP(JP,J,P)
        character(len=*), intent(in) :: JP
        integer, intent(out) :: J, P
        character(len=:), allocatable :: S
        integer :: Ios
        S = trim(adjustl(JP))
        if (len(S) < 2) call fatal('Invalid JP = "' // S // '"; use J+ or J-.')
        if (S(len(S):len(S)) == '+') then
            P = 1
        else if (S(len(S):len(S)) == '-') then
            P = 2
        else
            call fatal('Invalid JP = "' // S // '"; use J+ or J-.')
        end if
        read(S(1:len(S)-1),*,iostat=Ios) J
        if (Ios /= 0) call fatal('Invalid JP = "' // S // '"; use J+ or J-.')
    end subroutine parse_JP

    ! Radial transition densities of selected phonons ... HF_RPA_transition_densities ...
        ! Kind = 'RPA' (TDA & RPA files) or 'RRPA' (RRPA files only) ...
    subroutine phonon_transition_densities(Params,Kind,JP,nu_list,File_Name)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind, JP, File_Name
        integer, intent(in) :: nu_list(:)
        type(orb1b_t), allocatable :: Orb(:)
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        type(phstate_t), allocatable :: Phonon(:)
        integer, allocatable :: N_nu(:,:)
        type(ivec_t), allocatable :: Orb_Phonon(:,:)
        type(o1b_t) :: C, Rho
        type(tr1b_t) :: TrOp
        type(rmat_t), allocatable :: X_TDA(:,:)
        type(rvec_t), allocatable :: E_TDA(:,:)
        type(cmat_t), allocatable :: X_RPA(:,:), Y_RPA(:,:)
        type(cvec_t), allocatable :: E_RPA(:,:)
        complex(dp), allocatable :: pM_TDA(:), nM_TDA(:), isM_TDA(:), ivM_TDA(:), pM_RPA(:), nM_RPA(:), isM_RPA(:), ivM_RPA(:)
        real(dp), allocatable :: r_grid(:), pOrb_grid(:,:), nOrb_grid(:,:), pRho_TDA(:,:), nRho_TDA(:,:), pRho_RPA(:,:), nRho_RPA(:,:)
        real(dp) :: A_mass, Z, hw, nu_proton, nu_neutron, r_min, r_max, Amp, r
        complex(dp) :: Gauge_TDA, MT, MR, pTs, nTs, isTs, ivTs, pRs, nRs, isRs, ivRs, Op_ME
        real(dp) :: Gauge_RPA, pTsum, nTsum, pRsum, nRsum
        integer, parameter :: N_grid = 2**13 + 1
        integer :: J, P, N_phonon, N_ph, N_Ph_tot, a_max, i, a, k, nu_ind, nu, ph, i_ph, a_p, a_h
        logical :: With_TDA, Averaging
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        With_TDA = (Kind == 'RPA')
        A_mass = real(Params%Calc%A,dp)
        Z = real(Params%Calc%Z,dp)
        hw = Params%Calc%hw
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        nu_proton = 0.5_dp * m_p * hw / HbarC**2
        nu_neutron = 0.5_dp * m_n * hw / HbarC**2
        r_min = 0.0_dp + 1.0e-8_dp
        r_max = 2.5_dp * 1.2_dp * A_mass**(1.0_dp/3.0_dp)
        call parse_JP(JP,J,P)

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        call orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole,N_Ph_tot,Phonon)
        call HF_RPA_phonon_count(Params,N_Ph_tot,Phonon,N_nu,Orb_Phonon)
        call timer_start(t0)
        if (With_TDA) then
            C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        else
            C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_RRPA.bin'))
            Rho = O1b_import(Params,Orb,io_path(Params,'/Bin/Rho_RRPA.bin'))
        end if
        call timer_report(t0)
        call timer_start(t0)
        call HF_RPA_import_binary(Params,N_nu,Kind,With_TDA,X_TDA,E_TDA,X_RPA,Y_RPA,E_RPA)
        call timer_report(t0)

        N_phonon = size(nu_list)
        N_ph = N_nu(J,P)
        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        allocate(pM_TDA(N_phonon), nM_TDA(N_phonon), isM_TDA(N_phonon), ivM_TDA(N_phonon))
        allocate(pM_RPA(N_phonon), nM_RPA(N_phonon), isM_RPA(N_phonon), ivM_RPA(N_phonon))
        Averaging = (J == 0 .and. P == 1) .or. (J == 1 .and. P == 2) .or. (J == 2 .and. P == 1) .or. (J == 3 .and. P == 2)

        ! Transition matrix elements of the selected phonons ...
        if (Averaging) then
            do nu_ind = 1, N_phonon
                nu = nu_list(nu_ind)
                pTs = (0.0_dp,0.0_dp); nTs = (0.0_dp,0.0_dp); isTs = (0.0_dp,0.0_dp); ivTs = (0.0_dp,0.0_dp)
                pRs = (0.0_dp,0.0_dp); nRs = (0.0_dp,0.0_dp); isRs = (0.0_dp,0.0_dp); ivRs = (0.0_dp,0.0_dp)
                do ph = 1, N_ph
                    i_ph = Orb_Phonon(J,P)%v(ph)
                    call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                    if (With_TDA) then
                        MT = phase(J) * X_TDA(J,P)%v(ph,nu) / sqrt(real(2*J + 1,dp))
                    else
                        MT = (0.0_dp, 0.0_dp)
                    end if
                    MR = (phase(J) * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu)) / sqrt(real(2*J + 1,dp))
                    if (.not. With_TDA) MR = occ_weight(Phonon(i_ph)%tz,a_p,a_h) * MR
                    if (J == 1 .and. P == 2) then
                        if (Phonon(i_ph)%tz == -1) then
                            pTs = pTs + TrOp%E1%p(a_p,a_h) * MT
                            isTs = isTs + 0.5_dp * TrOp%E1_C%p(a_p,a_h) * MT
                            ivTs = ivTs + TrOp%E1%p(a_p,a_h) * MT * (A_mass - Z) / A_mass
                            pRs = pRs + TrOp%E1%p(a_p,a_h) * MR
                            isRs = isRs + 0.5_dp * TrOp%E1_C%p(a_p,a_h) * MR
                            ivRs = ivRs + TrOp%E1%p(a_p,a_h) * MR * (A_mass - Z) / A_mass
                        else
                            nTs = nTs + TrOp%E1%n(a_p,a_h) * MT
                            isTs = isTs + 0.5_dp * TrOp%E1_C%n(a_p,a_h) * MT
                            ivTs = ivTs - TrOp%E1%n(a_p,a_h) * MT * Z / A_mass
                            nRs = nRs + TrOp%E1%n(a_p,a_h) * MR
                            isRs = isRs + 0.5_dp * TrOp%E1_C%n(a_p,a_h) * MR
                            ivRs = ivRs - TrOp%E1%n(a_p,a_h) * MR * Z / A_mass
                        end if
                    else
                        Op_ME = multipole_ME(J,P,Phonon(i_ph)%tz,a_p,a_h)
                        MT = Op_ME * MT
                        MR = Op_ME * MR
                        if (Phonon(i_ph)%tz == -1) then
                            pTs = pTs + MT
                            isTs = isTs + 0.5_dp * MT
                            ivTs = ivTs + 0.5_dp * MT
                            pRs = pRs + MR
                            isRs = isRs + 0.5_dp * MR
                            ivRs = ivRs + 0.5_dp * MR
                        else
                            nTs = nTs + MT
                            isTs = isTs + 0.5_dp * MT
                            ivTs = ivTs - 0.5_dp * MT
                            nRs = nRs + MR
                            isRs = isRs + 0.5_dp * MR
                            ivRs = ivRs - 0.5_dp * MR
                        end if
                    end if
                end do
                pM_TDA(nu_ind) = pTs; nM_TDA(nu_ind) = nTs; isM_TDA(nu_ind) = isTs; ivM_TDA(nu_ind) = ivTs
                pM_RPA(nu_ind) = pRs; nM_RPA(nu_ind) = nRs; isM_RPA(nu_ind) = isRs; ivM_RPA(nu_ind) = ivRs
            end do
        end if

        ! Radial grid & reference basis orbitals on the grid ...
        r_grid = lin_range(r_min,r_max,N_grid)
        call orbitals_on_grid(r_grid,C,pOrb_grid,nOrb_grid)

        allocate(pRho_TDA(N_phonon,N_grid), nRho_TDA(N_phonon,N_grid), pRho_RPA(N_phonon,N_grid), nRho_RPA(N_phonon,N_grid))
        !$omp parallel do private(r,nu_ind,nu,pTsum,nTsum,pRsum,nRsum,Gauge_TDA,Gauge_RPA,ph,i_ph,a_p,a_h,Amp)
        do i = 1, N_grid
            r = r_grid(i)
            do nu_ind = 1, N_phonon
                nu = nu_list(nu_ind)
                pTsum = 0.0_dp; nTsum = 0.0_dp; pRsum = 0.0_dp; nRsum = 0.0_dp
                Gauge_TDA = (1.0_dp, 0.0_dp)
                Gauge_RPA = 1.0_dp
                if (Averaging) then
                    if (With_TDA) Gauge_TDA = pM_TDA(nu_ind) / abs(pM_TDA(nu_ind))
                    Gauge_RPA = real(pM_RPA(nu_ind)) / abs(pM_RPA(nu_ind))
                end if
                do ph = 1, N_ph
                    i_ph = Orb_Phonon(J,P)%v(ph)
                    call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                    if (Phonon(i_ph)%tz == -1) then
                        if (With_TDA) pTsum = pTsum + real(Gauge_TDA * phase(J) * X_TDA(J,P)%v(ph,nu)) * &
                                              pOrb_grid(a_p,i) * pOrb_grid(a_h,i)
                        pRsum = pRsum + Gauge_RPA * real(phase(J) * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu)) * &
                                pOrb_grid(a_p,i) * pOrb_grid(a_h,i)
                    else
                        if (With_TDA) nTsum = nTsum + real(Gauge_TDA * phase(J) * X_TDA(J,P)%v(ph,nu)) * &
                                              nOrb_grid(a_p,i) * nOrb_grid(a_h,i)
                        nRsum = nRsum + Gauge_RPA * real(phase(J) * X_RPA(J,P)%v(ph,nu) + Y_RPA(J,P)%v(ph,nu)) * &
                                nOrb_grid(a_p,i) * nOrb_grid(a_h,i)
                    end if
                end do
                Amp = r**2 / real(2*J + 1,dp)
                pRho_TDA(nu_ind,i) = Amp * pTsum
                nRho_TDA(nu_ind,i) = Amp * nTsum
                pRho_RPA(nu_ind,i) = Amp * pRsum
                nRho_RPA(nu_ind,i) = Amp * nRsum
            end do
        end do
        !$omp end parallel do

        ! Individual phonon densities ...
        do nu_ind = 1, N_phonon
            nu = nu_list(nu_ind)
            if (With_TDA) call write_density(io_path(Params,'/RPA/Densities/TDA_' // trim(File_Name) // '_JP' // trim(JP) // &
                                             '_nu' // istr(nu) // '_E' // jstr(round_digits(E_TDA(J,P)%v(nu),3)) // '.dat'), &
                                             pRho_TDA(nu_ind,:), nRho_TDA(nu_ind,:))
            call write_density(io_path(Params,'/' // Kind // '/Densities/' // Kind // '_' // trim(File_Name) // '_JP' // trim(JP) // &
                               '_nu' // istr(nu) // '_E' // jstr(round_digits(real(E_RPA(J,P)%v(nu)),3)) // '.dat'), &
                               pRho_RPA(nu_ind,:), nRho_RPA(nu_ind,:))
        end do
        write(*,'(a)') '"Export of individual phonon radial transition densities successfuly completed ..."'

        if (.not. Averaging) return

        ! Averaged densities weighted by the normalized transition matrix elements ...
        if (With_TDA) then
            call renormalize(pM_TDA)
            call renormalize(nM_TDA)
            call renormalize(isM_TDA)
            call renormalize(ivM_TDA)
        end if
        call renormalize(pM_RPA)
        call renormalize(nM_RPA)
        call renormalize(isM_RPA)
        call renormalize(ivM_RPA)
        if (With_TDA) call write_averaged(io_path(Params,'/RPA/Densities/TDA_' // trim(File_Name) // '_JP' // trim(JP) // &
                                          '_Averaged.dat'), pRho_TDA, nRho_TDA, pM_TDA, nM_TDA, isM_TDA, ivM_TDA)
        call write_averaged(io_path(Params,'/' // Kind // '/Densities/' // Kind // '_' // trim(File_Name) // '_JP' // trim(JP) // &
                            '_Averaged.dat'), pRho_RPA, nRho_RPA, pM_RPA, nM_RPA, isM_RPA, ivM_RPA)
        write(*,'(a)') '"Export of averaged phonon radial transitions densities successfuly completed ..."'

    contains

        ! RRPA weight sqrt(rho_hh - rho_pp) of a particle-hole component ...
        real(dp) function occ_weight(tz,a_p,a_h)
            integer, intent(in) :: tz, a_p, a_h
            if (tz == -1) then
                occ_weight = sqrt(Rho%p(a_h,a_h) - Rho%p(a_p,a_p))
            else
                occ_weight = sqrt(Rho%n(a_h,a_h) - Rho%n(a_p,a_p))
            end if
        end function occ_weight

        complex(dp) function multipole_ME(J,P,tz,a_p,a_h)
            integer, intent(in) :: J, P, tz, a_p, a_h
            real(dp) :: v
            v = 1.0_dp
            if (J == 0 .and. P == 1) then
                v = merge(TrOp%E0%p(a_p,a_h), TrOp%E0%n(a_p,a_h), tz == -1)
            else if (J == 2 .and. P == 1) then
                v = merge(TrOp%E2%p(a_p,a_h), TrOp%E2%n(a_p,a_h), tz == -1)
            else if (J == 3 .and. P == 2) then
                v = merge(TrOp%E3%p(a_p,a_h), TrOp%E3%n(a_p,a_h), tz == -1)
            end if
            multipole_ME = cmplx(v, 0.0_dp, kind=dp)
        end function multipole_ME

        ! M ./ sqrt(sum(M.^2)) with abs() ... complex sum of squares as in the original ...
        subroutine renormalize(M)
            complex(dp), intent(inout) :: M(:)
            M = abs(M) / sqrt(sum(M**2))
        end subroutine renormalize

        subroutine write_density(Path,pRho,nRho)
            character(len=*), intent(in) :: Path
            real(dp), intent(in) :: pRho(:), nRho(:)
            integer :: u, i
            open(newunit=u, file=Path, status='replace', action='write')
            write(u,'(a)') hdr5()
            do i = 1, N_grid
                write(u,'(a)') row5(r_grid(i), pRho(i), nRho(i), 0.5_dp * (pRho(i) + nRho(i)), 0.5_dp * (pRho(i) - nRho(i)))
            end do
            close(u)
        end subroutine write_density

        subroutine write_averaged(Path,pR,nR,pM,nM,isM,ivM)
            character(len=*), intent(in) :: Path
            real(dp), intent(in) :: pR(:,:), nR(:,:)
            complex(dp), intent(in) :: pM(:), nM(:), isM(:), ivM(:)
            integer :: u, i, k
            complex(dp) :: pS, nS, isS, ivS
            open(newunit=u, file=Path, status='replace', action='write')
            write(u,'(a)') hdr5()
            do i = 1, N_grid
                pS = (0.0_dp,0.0_dp); nS = (0.0_dp,0.0_dp); isS = (0.0_dp,0.0_dp); ivS = (0.0_dp,0.0_dp)
                do k = 1, N_phonon
                    pS = pS + pM(k) * pR(k,i)
                    nS = nS + nM(k) * nR(k,i)
                    isS = isS + 0.5_dp * isM(k) * (pR(k,i) + nR(k,i))
                    ivS = ivS + 0.5_dp * ivM(k) * (pR(k,i) - nR(k,i))
                end do
                write(u,'(a)') row5(r_grid(i), real(pS), real(nS), real(isS), real(ivS))
            end do
            close(u)
        end subroutine write_averaged

        function hdr5() result(s)
            character(len=:), allocatable :: s
            s = cfmt_s('r',20,.true.) // ' ' // cfmt_s('rho_p',20,.true.) // ' ' // cfmt_s('rho_n',20,.true.) // ' ' // &
                cfmt_s('rho_is',20,.true.) // ' ' // cfmt_s('rho_iv',20,.true.)
        end function hdr5

        function row5(x1,x2,x3,x4,x5) result(s)
            real(dp), intent(in) :: x1, x2, x3, x4, x5
            character(len=:), allocatable :: s
            s = cfmt_f(x1,20,8,.true.) // ' ' // cfmt_f(x2,20,8,.true.) // ' ' // cfmt_f(x3,20,8,.true.) // ' ' // &
                cfmt_f(x4,20,8,.true.) // ' ' // cfmt_f(x5,20,8,.true.)
        end function row5

        subroutine orbitals_on_grid(rg,Cm,pG,nG)
            real(dp), intent(in) :: rg(:)
            type(o1b_t), intent(in) :: Cm
            real(dp), allocatable, intent(out) :: pG(:,:), nG(:,:)
            integer :: ii, aa, kk
            real(dp) :: pS, nS
            allocate(pG(a_max,size(rg)), nG(a_max,size(rg)))
            !$omp parallel do private(aa,pS,nS,kk)
            do ii = 1, size(rg)
                do aa = 1, a_max
                    pS = 0.0_dp
                    nS = 0.0_dp
                    do kk = 1, a_max
                        if (Orb(aa)%l /= Orb(kk)%l .or. Orb(aa)%j /= Orb(kk)%j) cycle
                        pS = pS + Cm%p(kk,aa) * Psi_rad_LHO(rg(ii),Orb(kk)%n,Orb(kk)%l,nu_proton)
                        nS = nS + Cm%n(kk,aa) * Psi_rad_LHO(rg(ii),Orb(kk)%n,Orb(kk)%l,nu_neutron)
                    end do
                    pG(aa,ii) = pS
                    nG(aa,ii) = nS
                end do
            end do
            !$omp end parallel do
        end subroutine orbitals_on_grid

    end subroutine phonon_transition_densities

    ! Julia collect(range(a, stop = b, length = N)) ...
    function lin_range(a,b,N) result(g)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: N
        real(dp), allocatable :: g(:)
        integer :: i
        allocate(g(N))
        do i = 1, N
            g(i) = (a * real(N - i,dp) + b * real(i - 1,dp)) / real(N - 1,dp)
        end do
    end function lin_range

    ! Convective & spin transition current maps in the x-z plane ... HF_RPA_transition_currents ...
    subroutine phonon_transition_currents(Params,Kind,JP,nu_list,File_Name,M_in,xN_grid,zN_grid)
        type(params_t), intent(in) :: Params
        character(len=*), intent(in) :: Kind, JP, File_Name
        integer, intent(in) :: nu_list(:), M_in, xN_grid, zN_grid
        type(orb1b_t), allocatable :: Orb(:)
        type(pninteger_t) :: N_Particle, N_Hole
        type(pnsvector_t) :: Particle, Hole
        type(phstate_t), allocatable :: Phonon(:)
        integer, allocatable :: N_nu(:,:), i_map(:,:)
        type(ivec_t), allocatable :: Orb_Phonon(:,:)
        type(o1b_t) :: C, Rho
        type(tr1b_t) :: TrOp
        type(rmat_t), allocatable :: X_TDA(:,:)
        type(rvec_t), allocatable :: E_TDA(:,:)
        type(cmat_t), allocatable :: X_RPA(:,:), Y_RPA(:,:)
        type(cvec_t), allocatable :: E_RPA(:,:)
        real(dp), allocatable :: pM_TDA(:), nM_TDA(:), isM_TDA(:), ivM_TDA(:)
        complex(dp), allocatable :: pM_RPA(:), nM_RPA(:), isM_RPA(:), ivM_RPA(:)
        real(dp), allocatable :: r_grid(:), theta_grid(:), pR_grid(:,:), nR_grid(:,:), pdR_grid(:,:), ndR_grid(:,:), x_map(:,:)
        ! Single phonon component maps (component, ph, map point) ... ph & hp ordering, convective & spin ...
        real(dp), allocatable :: J_ph_conv(:,:,:), J_ph_spin(:,:,:), J_hp_conv(:,:,:), J_hp_spin(:,:,:)
        ! Phonon level maps (component, level, map point) ...
        real(dp), allocatable :: pT_conv(:,:,:), nT_conv(:,:,:), pT_spin(:,:,:), nT_spin(:,:,:)
        real(dp), allocatable :: pR_conv(:,:,:), nR_conv(:,:,:), pR_spin(:,:,:), nR_spin(:,:,:)
        real(dp), parameter :: g_n = -3.82608545_dp, g_p = 5.585694713_dp
        real(dp) :: A_mass, Z, hw, nu_proton, nu_neutron, r_min, r_max, x_min, x_max, z_min, z_max, dx, dz, dr, dtheta
        real(dp) :: x, zc, r, theta, MT, Gt, Wt
        complex(dp) :: MR, Gr
        integer, parameter :: rN_grid = 2**13 + 1, tN_grid = 256
        integer :: J, P, M, N_phonon, N_ph, N_Ph_tot, a_max, i, a, ind_x, ind_z, N_map, k_map, n_r, m_t, nu_ind, nu, ph, i_ph
        integer :: a_p, a_h
        logical :: With_TDA, Averaging
        integer(i8) :: t0

        call wigner_init_float(75,'Jmax',9)
        With_TDA = (Kind == 'RPA')
        A_mass = real(Params%Calc%A,dp)
        Z = real(Params%Calc%Z,dp)
        hw = Params%Calc%hw
        a_max = (Params%Calc%Nmax + 1)*(Params%Calc%Nmax + 2)/2
        nu_proton = 0.5_dp * m_p * hw / HbarC**2
        nu_neutron = 0.5_dp * m_n * hw / HbarC**2
        r_min = 0.01_dp
        r_max = 2.5_dp * 1.2_dp * A_mass**(1.0_dp/3.0_dp)
        call parse_JP(JP,J,P)
        M = M_in
        if (M == -1) M = 0

        Orb = orbitals_make(Params,Use_Int=.true.)        ! LEGACY Int.Nmax orbitals (old solvers, removed in cleanup stage)
        call orbitals_ph_make(Params,Orb,N_Particle,Particle,N_Hole,Hole)
        call orbitals_one_phonon_make(N_Particle,Particle,N_Hole,Hole,N_Ph_tot,Phonon)
        call HF_RPA_phonon_count(Params,N_Ph_tot,Phonon,N_nu,Orb_Phonon)
        call timer_start(t0)
        if (With_TDA) then
            C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_HF.bin'))
        else
            C = O1b_import(Params,Orb,io_path(Params,'/Bin/C_RRPA.bin'))
            Rho = O1b_import(Params,Orb,io_path(Params,'/Bin/Rho_RRPA.bin'))
        end if
        call timer_report(t0)
        call timer_start(t0)
        call HF_RPA_import_binary(Params,N_nu,Kind,With_TDA,X_TDA,E_TDA,X_RPA,Y_RPA,E_RPA)
        call timer_report(t0)

        N_phonon = size(nu_list)
        N_ph = N_nu(J,P)
        TrOp = Tr1b_initialize(Params,Orb,1.0_dp)
        TrOp = Tr1b_transformation(Params,Orb,C,TrOp)

        allocate(pM_TDA(N_phonon), nM_TDA(N_phonon), isM_TDA(N_phonon), ivM_TDA(N_phonon))
        allocate(pM_RPA(N_phonon), nM_RPA(N_phonon), isM_RPA(N_phonon), ivM_RPA(N_phonon))
        pM_TDA = 0.0_dp; nM_TDA = 0.0_dp; isM_TDA = 0.0_dp; ivM_TDA = 0.0_dp
        pM_RPA = (0.0_dp,0.0_dp); nM_RPA = (0.0_dp,0.0_dp); isM_RPA = (0.0_dp,0.0_dp); ivM_RPA = (0.0_dp,0.0_dp)
        Averaging = (J == 0 .and. P == 1) .or. (J == 1 .and. P == 2) .or. (J == 2 .and. P == 1) .or. (J == 3 .and. P == 2)

        if (Averaging) then
            do nu_ind = 1, N_phonon
                nu = nu_list(nu_ind)
                do ph = 1, N_ph
                    i_ph = Orb_Phonon(J,P)%v(ph)
                    call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                    if (With_TDA) then
                        MT = phase(J) / sqrt(real(2*J + 1,dp)) * X_TDA(J,P)%v(ph,nu)
                    else
                        MT = 0.0_dp
                    end if
                    MR = (phase(J) * X_RPA(J,P)%v(ph,nu) - Y_RPA(J,P)%v(ph,nu)) / sqrt(real(2*J + 1,dp))
                    if (.not. With_TDA) MR = occ_weight(Phonon(i_ph)%tz,a_p,a_h) * MR
                    if (J == 1 .and. P == 2) then
                        if (Phonon(i_ph)%tz == -1) then
                            pM_TDA(nu_ind) = pM_TDA(nu_ind) + TrOp%E1%p(a_p,a_h) * MT
                            isM_TDA(nu_ind) = isM_TDA(nu_ind) + 0.5_dp * TrOp%E1_C%p(a_p,a_h) * MT
                            ivM_TDA(nu_ind) = ivM_TDA(nu_ind) + TrOp%E1%p(a_p,a_h) * MT * (A_mass - Z) / A_mass
                            pM_RPA(nu_ind) = pM_RPA(nu_ind) + TrOp%E1%p(a_p,a_h) * MR
                            isM_RPA(nu_ind) = isM_RPA(nu_ind) + 0.5_dp * TrOp%E1_C%p(a_p,a_h) * MR
                            ivM_RPA(nu_ind) = ivM_RPA(nu_ind) + TrOp%E1%p(a_p,a_h) * MR * (A_mass - Z) / A_mass
                        else
                            nM_TDA(nu_ind) = nM_TDA(nu_ind) + TrOp%E1%n(a_p,a_h) * MT
                            isM_TDA(nu_ind) = isM_TDA(nu_ind) + 0.5_dp * TrOp%E1_C%n(a_p,a_h) * MT
                            ivM_TDA(nu_ind) = ivM_TDA(nu_ind) - TrOp%E1%n(a_p,a_h) * MT * Z / A_mass
                            nM_RPA(nu_ind) = nM_RPA(nu_ind) + TrOp%E1%n(a_p,a_h) * MR
                            isM_RPA(nu_ind) = isM_RPA(nu_ind) + 0.5_dp * TrOp%E1_C%n(a_p,a_h) * MR
                            ivM_RPA(nu_ind) = ivM_RPA(nu_ind) - TrOp%E1%n(a_p,a_h) * MR * Z / A_mass
                        end if
                    else
                        MT = op_ME(J,P,Phonon(i_ph)%tz,a_p,a_h) * MT
                        MR = op_ME(J,P,Phonon(i_ph)%tz,a_p,a_h) * MR
                        if (Phonon(i_ph)%tz == -1) then
                            pM_TDA(nu_ind) = pM_TDA(nu_ind) + MT
                            isM_TDA(nu_ind) = isM_TDA(nu_ind) + 0.5_dp * MT
                            ivM_TDA(nu_ind) = ivM_TDA(nu_ind) + 0.5_dp * MT
                            pM_RPA(nu_ind) = pM_RPA(nu_ind) + MR
                            isM_RPA(nu_ind) = isM_RPA(nu_ind) + 0.5_dp * MR
                            ivM_RPA(nu_ind) = ivM_RPA(nu_ind) + 0.5_dp * MR
                        else
                            nM_TDA(nu_ind) = nM_TDA(nu_ind) + MT
                            isM_TDA(nu_ind) = isM_TDA(nu_ind) + 0.5_dp * MT
                            ivM_TDA(nu_ind) = ivM_TDA(nu_ind) - 0.5_dp * MT
                            nM_RPA(nu_ind) = nM_RPA(nu_ind) + MR
                            isM_RPA(nu_ind) = isM_RPA(nu_ind) + 0.5_dp * MR
                            ivM_RPA(nu_ind) = ivM_RPA(nu_ind) - 0.5_dp * MR
                        end if
                    end if
                end do
            end do
        end if

        ! Radial & angular grids, orbitals & their radial derivatives ...
        r_grid = lin_range(r_min,r_max,rN_grid)
        theta_grid = lin_range(0.0_dp,2.0_dp*pi,tN_grid)
        allocate(pR_grid(a_max,rN_grid), nR_grid(a_max,rN_grid), pdR_grid(a_max,rN_grid), ndR_grid(a_max,rN_grid))
        !$omp parallel do private(a)
        do i = 1, rN_grid
            do a = 1, a_max
                pR_grid(a,i) = radial_orbital(a,r_grid(i),C%p,nu_proton)
                nR_grid(a,i) = radial_orbital(a,r_grid(i),C%n,nu_neutron)
            end do
        end do
        !$omp end parallel do
        do a = 1, a_max
            pdR_grid(a,:) = differentiate_Oh6(r_grid,pR_grid(a,:))
            ndR_grid(a,:) = differentiate_Oh6(r_grid,nR_grid(a,:))
        end do

        ! Cartesian map points snapped to the radial & angular grids ...
        x_min = -0.8_dp * r_max
        x_max = 0.8_dp * r_max
        z_min = -0.8_dp * r_max
        z_max = 0.8_dp * r_max
        dx = (x_max - x_min) / real(xN_grid,dp)
        dz = (z_max - z_min) / real(zN_grid,dp)
        N_map = xN_grid * zN_grid
        allocate(i_map(2,N_map), x_map(2,N_map))
        dr = (r_max - r_min) / real(rN_grid - 1,dp)
        dtheta = 2.0_dp * pi / real(tN_grid - 1,dp)
        k_map = 0
        do ind_x = 1, xN_grid
            x = x_min + real(ind_x,dp) * dx
            do ind_z = 1, zN_grid
                zc = z_min + real(ind_z,dp) * dz
                r = sqrt(x**2 + zc**2)
                theta = atan2(x,zc)
                if (theta < 0.0_dp) theta = theta + 2.0_dp*pi
                n_r = min(max(round_even((r - r_min) / dr) + 1, 1), rN_grid)
                m_t = min(max(round_even((theta - 0.0_dp) / dtheta) + 1, 1), tN_grid)
                k_map = k_map + 1
                i_map(1,k_map) = n_r
                i_map(2,k_map) = m_t
                x_map(1,k_map) = r_grid(n_r)
                x_map(2,k_map) = theta_grid(m_t)
            end do
        end do

        ! 1-body transition current matrix elements on the map ...
        allocate(J_ph_conv(3,N_ph,N_map), J_ph_spin(3,N_ph,N_map), J_hp_conv(3,N_ph,N_map), J_hp_spin(3,N_ph,N_map))
        call timer_start(t0)
        !$omp parallel do schedule(dynamic)
        do k_map = 1, N_map
            call map_point(k_map)
        end do
        !$omp end parallel do
        call timer_report(t0)

        ! Transition currents of the selected phonon levels ...
        allocate(pT_conv(3,N_phonon,N_map), nT_conv(3,N_phonon,N_map), pT_spin(3,N_phonon,N_map), nT_spin(3,N_phonon,N_map))
        allocate(pR_conv(3,N_phonon,N_map), nR_conv(3,N_phonon,N_map), pR_spin(3,N_phonon,N_map), nR_spin(3,N_phonon,N_map))
        pT_conv = 0.0_dp; nT_conv = 0.0_dp; pT_spin = 0.0_dp; nT_spin = 0.0_dp
        pR_conv = 0.0_dp; nR_conv = 0.0_dp; pR_spin = 0.0_dp; nR_spin = 0.0_dp
        !$omp parallel do private(nu_ind,nu,Gt,Gr,ph,i_ph,MT,a_p,a_h,Wt)
        do k_map = 1, N_map
            do nu_ind = 1, N_phonon
                nu = nu_list(nu_ind)
                Gt = 1.0_dp
                Gr = (1.0_dp, 0.0_dp)
                if (Averaging) then
                    if (With_TDA) Gt = pM_TDA(nu_ind) / abs(pM_TDA(nu_ind))
                    Gr = pM_RPA(nu_ind) / abs(pM_RPA(nu_ind))
                end if
                do ph = 1, N_ph
                    i_ph = Orb_Phonon(J,P)%v(ph)
                    if (With_TDA) then
                        MT = Gt * X_TDA(J,P)%v(ph,nu)
                    else
                        MT = 0.0_dp
                    end if
                    call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                    Wt = 1.0_dp
                    if (.not. With_TDA) Wt = occ_weight(Phonon(i_ph)%tz,a_p,a_h)
                    call add_level(Phonon(i_ph)%tz == -1, ph, k_map, nu_ind, MT, Wt * real(Gr * X_RPA(J,P)%v(ph,nu)), &
                                   Wt * real(Gr * Y_RPA(J,P)%v(ph,nu)))
                end do
            end do
        end do
        !$omp end parallel do

        ! Physical units [c] ...
        pT_conv = pT_conv * HbarC / m_p
        nT_conv = nT_conv * HbarC / m_p
        pT_spin = pT_spin * HbarC / m_p
        nT_spin = nT_spin * HbarC / m_p
        pR_conv = pR_conv * HbarC / m_p
        nR_conv = nR_conv * HbarC / m_p
        pR_spin = pR_spin * HbarC / m_p
        nR_spin = nR_spin * HbarC / m_p

        ! Individual phonon level maps ...
        do nu_ind = 1, N_phonon
            nu = nu_list(nu_ind)
            if (With_TDA) then
                call write_map(level_path('TDA','Convective',jstr(round_digits(E_TDA(J,P)%v(nu),3)),nu), pT_conv, nT_conv, nu_ind)
                call write_map(level_path('TDA','Spin',jstr(round_digits(E_TDA(J,P)%v(nu),3)),nu), pT_spin, nT_spin, nu_ind)
            end if
            call write_map(level_path(Kind,'Convective',jstr(round_digits(real(E_RPA(J,P)%v(nu)),3)),nu), pR_conv, nR_conv, nu_ind)
            call write_map(level_path(Kind,'Spin',jstr(round_digits(real(E_RPA(J,P)%v(nu)),3)),nu), pR_spin, nR_spin, nu_ind)
        end do
        write(*,'(a)') '"Export of individual phonon level transition current maps successfuly completed ..."'

        if (.not. Averaging) return

        if (With_TDA) then
            pM_TDA = abs(pM_TDA) / sqrt(sum(pM_TDA**2))
            nM_TDA = abs(nM_TDA) / sqrt(sum(nM_TDA**2))
            isM_TDA = abs(isM_TDA) / sqrt(sum(isM_TDA**2))
            ivM_TDA = abs(ivM_TDA) / sqrt(sum(ivM_TDA**2))
        end if
        pM_RPA = abs(pM_RPA) / sqrt(sum(pM_RPA**2))
        nM_RPA = abs(nM_RPA) / sqrt(sum(nM_RPA**2))
        isM_RPA = abs(isM_RPA) / sqrt(sum(isM_RPA**2))
        ivM_RPA = abs(ivM_RPA) / sqrt(sum(ivM_RPA**2))

        if (With_TDA) then
            call write_averaged_map(avg_path('TDA','Convective'), pT_conv, nT_conv, pM_TDA, nM_TDA, isM_TDA, ivM_TDA)
            call write_averaged_map(avg_path('TDA','Spin'), pT_spin, nT_spin, pM_TDA, nM_TDA, isM_TDA, ivM_TDA)
        end if
        call write_averaged_map(avg_path(Kind,'Convective'), pR_conv, nR_conv, real(pM_RPA), real(nM_RPA), real(isM_RPA), &
                                real(ivM_RPA))
        call write_averaged_map(avg_path(Kind,'Spin'), pR_spin, nR_spin, real(pM_RPA), real(nM_RPA), real(isM_RPA), real(ivM_RPA))
        write(*,'(a)') '"Export of averaged phonon level transition current maps successfuly completed ..."'

    contains

        ! RRPA weight sqrt(rho_hh - rho_pp) of a particle-hole component ...
        real(dp) function occ_weight(tz,a_p,a_h)
            integer, intent(in) :: tz, a_p, a_h
            if (tz == -1) then
                occ_weight = sqrt(Rho%p(a_h,a_h) - Rho%p(a_p,a_p))
            else
                occ_weight = sqrt(Rho%n(a_h,a_h) - Rho%n(a_p,a_p))
            end if
        end function occ_weight

        real(dp) function op_ME(J,P,tz,a_p,a_h)
            integer, intent(in) :: J, P, tz, a_p, a_h
            op_ME = 1.0_dp
            if (J == 0 .and. P == 1) then
                op_ME = merge(TrOp%E0%p(a_p,a_h), TrOp%E0%n(a_p,a_h), tz == -1)
            else if (J == 2 .and. P == 1) then
                op_ME = merge(TrOp%E2%p(a_p,a_h), TrOp%E2%n(a_p,a_h), tz == -1)
            else if (J == 3 .and. P == 2) then
                op_ME = merge(TrOp%E3%p(a_p,a_h), TrOp%E3%n(a_p,a_h), tz == -1)
            end if
        end function op_ME

        real(dp) function radial_orbital(a,r,Cm,nu_osc)
            integer, intent(in) :: a
            real(dp), intent(in) :: r, Cm(:,:), nu_osc
            integer :: kk
            radial_orbital = 0.0_dp
            do kk = 1, a_max
                if (Orb(a)%l /= Orb(kk)%l .or. Orb(a)%j /= Orb(kk)%j) cycle
                radial_orbital = radial_orbital + Cm(kk,a) * Psi_rad_LHO(r,Orb(kk)%n,Orb(kk)%l,nu_osc)
            end do
        end function radial_orbital

        ! Particle-hole & hole-particle 1-body currents of every phonon component at map point k ...
        subroutine map_point(k)
            integer, intent(in) :: k
            integer :: n_r, m_t, ph, i_ph, a_p, a_h, l_p, j_p, l_h, j_h, m_p, m_h
            real(dp) :: r, theta, g, R_p, R_h, dR_p, dR_h, Amp
            complex(dp) :: jx, jy, jz
            real(dp) :: Jc_ph(3), Js_ph(3), Jc_hp(3), Js_hp(3)

            n_r = i_map(1,k)
            m_t = i_map(2,k)
            r = r_grid(n_r)
            theta = theta_grid(m_t)
            do ph = 1, N_ph
                i_ph = Orb_Phonon(J,P)%v(ph)
                call ph_orbitals(i_ph,Phonon,Particle,Hole,a_p,a_h)
                l_p = Orb(a_p)%l
                j_p = Orb(a_p)%j
                l_h = Orb(a_h)%l
                j_h = Orb(a_h)%j
                Jc_ph = 0.0_dp
                Js_ph = 0.0_dp
                Jc_hp = 0.0_dp
                Js_hp = 0.0_dp
                if (Phonon(i_ph)%tz == -1) then
                    R_p = pR_grid(a_p,n_r)
                    R_h = pR_grid(a_h,n_r)
                    dR_p = pdR_grid(a_p,n_r)
                    dR_h = pdR_grid(a_h,n_r)
                    g = 0.5_dp * g_p
                else
                    R_p = nR_grid(a_p,n_r)
                    R_h = nR_grid(a_h,n_r)
                    dR_p = ndR_grid(a_p,n_r)
                    dR_h = ndR_grid(a_h,n_r)
                    g = 0.5_dp * g_n
                end if
                do m_p = -j_p, j_p, 2
                    do m_h = -j_h, j_h, 2
                        ! Particle-hole contributions ...
                        if ((m_p + m_h) == 2*M) then
                            Amp = fCG(j_p,j_h,2*J,m_p,m_h,2*M) * phase((j_h - m_h)/2)
                            call J1b_convective(r,theta,0.0_dp,l_p,j_p,m_p,R_p,dR_p,l_h,j_h,-m_h,R_h,dR_h,jx,jy,jz)
                            Jc_ph = Jc_ph + Amp * [real(jx), real(jy), real(jz)]
                            call J1b_spin(r,theta,0.0_dp,l_p,j_p,m_p,R_p,dR_p,l_h,j_h,-m_h,R_h,dR_h,jx,jy,jz)
                            Js_ph = Js_ph + Amp * [real(g*jx), real(g*jy), real(g*jz)]
                        end if
                        ! Hole-particle contributions ...
                        if ((m_p + m_h) == -2*M) then
                            Amp = fCG(j_p,j_h,2*J,m_p,m_h,-2*M) * phase(J + M + (j_h + m_h)/2)
                            call J1b_convective(r,theta,0.0_dp,l_h,j_h,-m_h,R_h,dR_h,l_p,j_p,m_p,R_p,dR_p,jx,jy,jz)
                            Jc_hp = Jc_hp + Amp * [real(jx), real(jy), real(jz)]
                            call J1b_spin(r,theta,0.0_dp,l_h,j_h,-m_h,R_h,dR_h,l_p,j_p,m_p,R_p,dR_p,jx,jy,jz)
                            Js_hp = Js_hp + Amp * [real(g*jx), real(g*jy), real(g*jz)]
                        end if
                    end do
                end do
                J_ph_conv(:,ph,k) = Jc_ph
                J_ph_spin(:,ph,k) = Js_ph
                J_hp_conv(:,ph,k) = Jc_hp
                J_hp_spin(:,ph,k) = Js_hp
            end do
        end subroutine map_point

        subroutine add_level(Proton,ph,k,lv,xT,xR,yR)
            logical, intent(in) :: Proton
            integer, intent(in) :: ph, k, lv
            real(dp), intent(in) :: xT, xR, yR
            if (Proton) then
                pT_conv(:,lv,k) = pT_conv(:,lv,k) + xT * J_ph_conv(:,ph,k)
                pT_spin(:,lv,k) = pT_spin(:,lv,k) + xT * J_ph_spin(:,ph,k)
                pR_conv(:,lv,k) = pR_conv(:,lv,k) + (xR * J_ph_conv(:,ph,k) - yR * J_hp_conv(:,ph,k))
                pR_spin(:,lv,k) = pR_spin(:,lv,k) + (xR * J_ph_spin(:,ph,k) - yR * J_hp_spin(:,ph,k))
            else
                nT_conv(:,lv,k) = nT_conv(:,lv,k) + xT * J_ph_conv(:,ph,k)
                nT_spin(:,lv,k) = nT_spin(:,lv,k) + xT * J_ph_spin(:,ph,k)
                nR_conv(:,lv,k) = nR_conv(:,lv,k) + (xR * J_ph_conv(:,ph,k) - yR * J_hp_conv(:,ph,k))
                nR_spin(:,lv,k) = nR_spin(:,lv,k) + (xR * J_ph_spin(:,ph,k) - yR * J_hp_spin(:,ph,k))
            end if
        end subroutine add_level

        function level_path(Pre,Typ,sE,nu) result(s)
            character(len=*), intent(in) :: Pre, Typ, sE
            integer, intent(in) :: nu
            character(len=:), allocatable :: s
            character(len=:), allocatable :: Dir
            Dir = merge('RPA ', Kind // repeat(' ',max(0,4-len(Kind))), Pre == 'TDA')
            s = io_path(Params,'/' // trim(Dir) // '/Densities/' // Pre // '_' // Typ // '_' // trim(File_Name) // '_JP' // &
                trim(JP) // '_nu' // istr(nu) // '_E' // sE // '.dat')
        end function level_path

        function avg_path(Pre,Typ) result(s)
            character(len=*), intent(in) :: Pre, Typ
            character(len=:), allocatable :: s
            character(len=:), allocatable :: Dir
            Dir = merge('RPA ', Kind // repeat(' ',max(0,4-len(Kind))), Pre == 'TDA')
            s = io_path(Params,'/' // trim(Dir) // '/Densities/' // Pre // '_' // Typ // '_' // trim(File_Name) // '_JP' // &
                trim(JP) // '_Averaged.dat')
        end function avg_path

        subroutine write_map(Path,pJ,nJ,lv)
            character(len=*), intent(in) :: Path
            real(dp), intent(in) :: pJ(:,:,:), nJ(:,:,:)
            integer, intent(in) :: lv
            integer :: u, k
            real(dp) :: xx, zz, pJ_x, pJ_z, nJ_x, nJ_z
            character(len=8), parameter :: H(10) = ['x    ', 'z    ', 'pJ_x ', 'pJ_z ', 'nJ_x ', 'nJ_z ', 'isJ_x', 'isJ_z', &
                                                    'ivJ_x', 'ivJ_z']
            open(newunit=u, file=Path, status='replace', action='write')
            write(u,'(a)') header(H)
            do k = 1, N_map
                xx = x_map(1,k) * sin(x_map(2,k))
                zz = x_map(1,k) * cos(x_map(2,k))
                pJ_x = pJ(1,lv,k)
                pJ_z = pJ(3,lv,k)
                nJ_x = nJ(1,lv,k)
                nJ_z = nJ(3,lv,k)
                write(u,'(a)') rowf([xx, zz, pJ_x, pJ_z, nJ_x, nJ_z, 0.5_dp*(pJ_x + nJ_x), 0.5_dp*(pJ_z + nJ_z), &
                                     0.5_dp*(pJ_x - nJ_x), 0.5_dp*(pJ_z - nJ_z)])
            end do
            close(u)
        end subroutine write_map

        subroutine write_averaged_map(Path,pJ,nJ,pM,nM,isM,ivM)
            character(len=*), intent(in) :: Path
            real(dp), intent(in) :: pJ(:,:,:), nJ(:,:,:), pM(:), nM(:), isM(:), ivM(:)
            integer :: u, k, lv
            real(dp) :: xx, zz, S(12), pJ_x, pJ_z, nJ_x, nJ_z
            character(len=8), parameter :: H(14) = ['x    ', 'z    ', 'ppJ_x', 'ppJ_z', 'pnJ_x', 'pnJ_z', 'npJ_x', 'npJ_z', &
                                                    'nnJ_x', 'nnJ_z', 'isJ_x', 'isJ_z', 'ivJ_x', 'ivJ_z']
            open(newunit=u, file=Path, status='replace', action='write')
            write(u,'(a)') header(H)
            do k = 1, N_map
                xx = x_map(1,k) * sin(x_map(2,k))
                zz = x_map(1,k) * cos(x_map(2,k))
                S = 0.0_dp
                do lv = 1, N_phonon
                    pJ_x = pJ(1,lv,k)
                    pJ_z = pJ(3,lv,k)
                    nJ_x = nJ(1,lv,k)
                    nJ_z = nJ(3,lv,k)
                    S(1) = S(1) + pM(lv) * pJ_x
                    S(2) = S(2) + pM(lv) * pJ_z
                    S(3) = S(3) + pM(lv) * nJ_x
                    S(4) = S(4) + pM(lv) * nJ_z
                    S(5) = S(5) + nM(lv) * pJ_x
                    S(6) = S(6) + nM(lv) * pJ_z
                    S(7) = S(7) + nM(lv) * nJ_x
                    S(8) = S(8) + nM(lv) * nJ_z
                    S(9) = S(9) + 0.5_dp * isM(lv) * (pJ_x + nJ_x)
                    S(10) = S(10) + 0.5_dp * isM(lv) * (pJ_z + nJ_z)
                    S(11) = S(11) + 0.5_dp * ivM(lv) * (pJ_x - nJ_x)
                    S(12) = S(12) + 0.5_dp * ivM(lv) * (pJ_z - nJ_z)
                end do
                write(u,'(a)') rowf([xx, zz, S])
            end do
            close(u)
        end subroutine write_averaged_map

        function header(H) result(s)
            character(len=*), intent(in) :: H(:)
            character(len=:), allocatable :: s
            integer :: k
            s = cfmt_s(trim(H(1)),20,.true.)
            do k = 2, size(H)
                s = s // ' ' // cfmt_s(trim(H(k)),20,.true.)
            end do
        end function header

        function rowf(v) result(s)
            real(dp), intent(in) :: v(:)
            character(len=:), allocatable :: s
            integer :: k
            s = cfmt_f(v(1),20,8,.true.)
            do k = 2, size(v)
                s = s // ' ' // cfmt_f(v(k),20,8,.true.)
            end do
        end function rowf

    end subroutine phonon_transition_currents

end module nj_hf_rpa_external
