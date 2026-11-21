! Hamiltonian of the HF-ScRPA iteration in the current natural-orbital (NAT) basis (NuJuliet dev, stage 4) ... counterpart
    ! of src/RPA/RPA_Hamiltonian.jl (RPA_Hamiltonian, RPA_Hamiltonian_HF_reference) ... replaces the base-version
    ! HF_RRPA_h1b / HF_RRPA_V2b / HF_RRPA_V2b_transform (HF-RRPA, removed) ... uses the core H1b / H2b / H1b_H2b (nj_hamiltonian) ...
    ! Deviation (B6, no change of results beyond rounding): for ScH & TBCF dev transforms the caller's Sigma in place to the LHO
    ! basis & back (NAT -> LHO -> NAT); the port transforms a copy, the caller's Sigma is unchanged ...
module nj_rpa_hamiltonian
    use nj_kinds
    use nj_types
    use nj_params
    use nj_utils
    use nj_o1b, only : O1b_import, O1b_transformation, O1b_transpose
    use nj_o2b, only : O2b_import, O2b_transformation
    use nj_hamiltonian, only : H1b, H2b, H1b_H2b
    implicit none
    private

    public :: RPA_Hamiltonian, RPA_Hamiltonian_HF_reference

contains

    ! h_N & V_NN in the NAT basis ... C = C_LHO_NAT, Rho & Sigma in the NAT basis, T, V_NN_bare & V_NNN_bare bare (LHO) ...
        !   ScH = true:  Rho_L = C Rho C^T, h_L = H1b(Rho_L) [TBCF: H1b_H2b(Rho_L, Sigma_L = C Sigma C^T) with the V3N.Sigma term],
        !                V_L = H2b(Rho_L), h_N = C^T h_L C, V_NN = V_L transformed with C ...
        !   ScH = false: the reference H1b.bin / H2b.bin rotated HF -> NAT (RPA_Hamiltonian_HF_reference); Sigma, T & the bare
        !                interactions are not used (they may be unallocated) ...
    subroutine RPA_Hamiltonian(Params,Orb,Orb_NN,Orb_NNN,C,Rho,Sigma,T,V_NN_bare,V_NNN_bare,h_N,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(orb3b_t), intent(in) :: Orb_NNN
        type(o1b_t), intent(in) :: C, Rho, T
        type(o2b_t), intent(in) :: Sigma, V_NN_bare
        type(v3b_t), intent(in) :: V_NNN_bare
        type(o1b_t), intent(out) :: h_N
        type(o2b_t), intent(out) :: V_NN
        type(o1b_t) :: Rho_L, h_L, C_T
        type(o2b_t) :: Sigma_L

        write(*,'(/,a)') 'Constructing Hamiltonian for RPA calculation ...'
        if (Params%Calc%RPA%ScH) then
            C_T = O1b_transpose(C)
            Rho_L = O1b_transformation(Rho,C_T)
            if (Params%Calc%RPA%TBCF) then
                Sigma_L = Sigma
                call O2b_transformation(Params,Orb,Orb_NN,Sigma_L,C_T,Report=.false.)
                call H1b_H2b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN_bare,V_NNN_bare,Rho_L,Sigma_L,h_L,V_NN,Report=.false.)
            else
                h_L = H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN_bare,V_NNN_bare,Rho_L)
                call H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN_bare,V_NNN_bare,Rho_L,V_NN,Report=.false.)
            end if
            h_N = O1b_transformation(h_L,C)
            call O2b_transformation(Params,Orb,Orb_NN,V_NN,C,Report=.false.)
        else
            call RPA_Hamiltonian_HF_reference(Params,Orb,Orb_NN,C,h_N,V_NN)
        end if
        write(*,'(a)') TAB // 'Finished construction of RPA Hamiltonian ...'
    end subroutine RPA_Hamiltonian

    ! Reference (HF) Hamiltonian io/<Path>/<RPA.Reference>/bin/{H1b,H2b}.bin rotated to the NAT basis with
        ! C_HF_NAT = C_LHO_HF^T C_LHO_NAT (C_LHO_HF from C.bin) ... counterpart of RPA_Hamiltonian_HF_reference ...
    subroutine RPA_Hamiltonian_HF_reference(Params,Orb,Orb_NN,C,h_N,V_NN)
        type(params_t), intent(in) :: Params
        type(orb1b_t), intent(in) :: Orb(:)
        type(orb2b_t), intent(in) :: Orb_NN
        type(o1b_t), intent(in) :: C
        type(o1b_t), intent(out) :: h_N
        type(o2b_t), intent(out) :: V_NN
        type(o1b_t) :: H_ref, C_HF, C_HF_NAT
        character(len=:), allocatable :: IO

        IO = io_dir(Params,Params%Calc%RPA%Reference)
        H_ref = O1b_import(Params,Orb,IO // '/bin/H1b.bin')
        call O2b_import(Params,Orb,IO // '/bin/H2b.bin',V_NN)
        C_HF = O1b_import(Params,Orb,IO // '/bin/C.bin')
        C_HF_NAT%p = matmul(transpose(C_HF%p),C%p)
        C_HF_NAT%n = matmul(transpose(C_HF%n),C%n)
        h_N = O1b_transformation(H_ref,C_HF_NAT)
        call O2b_transformation(Params,Orb,Orb_NN,V_NN,C_HF_NAT,Report=.false.)
    end subroutine RPA_Hamiltonian_HF_reference

end module nj_rpa_hamiltonian
