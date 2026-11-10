! Data structures of NuJuliet ... counterparts of Import_Structures.jl ...
    ! Conventions of the port:
    !   Julia X[P,J+1] (channel of parity P & angular momentum J)  ->  Fortran X(P,J) with J starting at 0 ...
    !   Julia Orb_NN.Ind[t,P,J+1][Bra][k]                           ->  Fortran Orb_NN%Ind(t,P,J)%v(k,Bra) ...
    !   Julia Dict lookups of 2-body orbitals are replaced by a dense index table Orb_NN%Dic(a,b,J,T) ...
module nj_types
    use nj_kinds
    use nj_hashmap
    implicit none
    private

    ! Generic jagged-array containers ...
    type, public :: ivec_t
        integer, allocatable :: v(:)
    end type ivec_t

    type, public :: imat_t
        integer, allocatable :: v(:,:)
    end type imat_t

    type, public :: rvec_t
        real(dp), allocatable :: v(:)
    end type rvec_t

    type, public :: rmat_t
        real(dp), allocatable :: v(:,:)
    end type rmat_t

    type, public :: svec_t
        real(sp), allocatable :: v(:)
    end type svec_t

    type, public :: cvec_t
        complex(dp), allocatable :: v(:)
    end type cvec_t

    type, public :: cmat_t
        complex(dp), allocatable :: v(:,:)
    end type cmat_t

    ! Orbitals & single-particle states ...
    type, public :: orb1b_t
        integer :: a = 0, n = 0, l = 0, j = 0, pO = 0, nO = 0
    end type orb1b_t

    type, public :: spstate_t
        integer :: a = 0, l = 0, j = 0
        real(dp) :: E = 0.0_dp
    end type spstate_t

    type, public :: spstate_vec_t
        type(spstate_t), allocatable :: v(:)
    end type spstate_vec_t

    type, public :: pnsvector_t
        type(spstate_t), allocatable :: p(:), n(:)
    end type pnsvector_t

    type, public :: pninteger_t
        integer :: p = 0, n = 0
    end type pninteger_t

    type, public :: pnfloat_t
        real(dp) :: p = 0.0_dp, n = 0.0_dp
    end type pnfloat_t

    type, public :: pnvector_t
        real(dp), allocatable :: p(:), n(:)
    end type pnvector_t

    type, public :: pncvector_t
        complex(dp), allocatable :: p(:), n(:)
    end type pncvector_t

    ! pnArray of particle-hole phonon indices ... (J+1, P, p, h) in Julia -> (0:J, P, p, h) here ...
    type, public :: pnarray4_t
        integer, allocatable :: p(:,:,:,:), n(:,:,:,:)
    end type pnarray4_t

    type, public :: phstate_t
        integer :: p = 0, h = 0, tz = 0, J = 0, Par = 0
    end type phstate_t

    type, public :: ind2qp_t
        integer :: a = 0, b = 0, T = 0
    end type ind2qp_t

    type, public :: ind2qp_vec_t
        type(ind2qp_t), allocatable :: v(:)
    end type ind2qp_vec_t

    ! 2-quasiparticle orbitals ... N(P,0:J) & i(P,0:J)%v(:) ...
    type, public :: qporb2b_t
        integer, allocatable :: N(:,:)
        type(ind2qp_vec_t), allocatable :: i(:,:)
    end type qporb2b_t

    ! 1-body operators ...
    type, public :: o1b_t
        real(dp), allocatable :: p(:,:), n(:,:)
    end type o1b_t

    type, public :: qpo1b_t
        type(o1b_t) :: qp11, qp20
    end type qpo1b_t

    ! Standard 2-body NN operator ... packed lower triangles per channel (P, 0:J_max) ...
    type, public :: o2b_t
        type(rvec_t), allocatable :: pp(:,:), pn(:,:), nn(:,:)
    end type o2b_t

    ! 2-body orbitals ... N(T+1,P,0:J), Ind(T+1,P,0:J)%v(1:2,Bra) & Dic(a,b,J,T) giving the 2-body index ...
    type, public :: orb2b_t
        integer :: a_max = 0, J_max = 0
        integer, allocatable :: N(:,:,:)
        type(imat_t), allocatable :: Ind(:,:,:)
        integer, allocatable :: Dic(:,:,:,:)
    end type orb2b_t

    ! 3-body NO2B orbitals ... hash map of 2+1-body keys & counts N(P,0:l,0:J,T,j_index) ...
    type, public :: orb3b_t
        type(hashmap_t) :: Dic
        integer, allocatable :: N(:,:,:,:,:)
    end type orb3b_t

    ! 3-body NO2B interaction in single precision ... V(P,0:l,0:J,T,j_index)%v(:) ...
    type, public :: v3b_t
        type(svec_t), allocatable :: V(:,:,:,:,:)
    end type v3b_t

    ! Temporary full-matrix 2-body operators used in basis transformations ...
    type, public :: o2b_temp_t
        type(rmat_t), allocatable :: pp(:,:), pn(:,:), nn(:,:)
    end type o2b_temp_t

    type, public :: orb2b_temp_t
        integer, allocatable :: N(:,:)
        type(imat_t), allocatable :: Ind(:,:)
        integer, allocatable :: Dic(:,:,:)
    end type orb2b_temp_t

    ! Quasiparticle 2-body operators ...
    type, public :: qpo40_t
        type(rvec_t), allocatable :: pp(:,:), nn(:,:)
        type(rmat_t), allocatable :: pn(:,:)
    end type qpo40_t

    type, public :: qpo31_t
        type(rmat_t), allocatable :: pp(:,:), pn2011(:,:), pn1120(:,:), nn(:,:)
    end type qpo31_t

    type, public :: qpo22_t
        type(rvec_t), allocatable :: pp(:,:), pn1111(:,:), nn(:,:)
        type(rmat_t), allocatable :: pn2002(:,:)
    end type qpo22_t

    type, public :: qpo2b_t
        type(qpo40_t) :: qp40
        type(qpo31_t) :: qp31
        type(qpo22_t) :: qp22
    end type qpo2b_t

    ! Reduced multipole transition matrix elements per phonon ... Julia struct ReducedMultipole ...
    type, public :: reduced_multipole_t
        type(pncvector_t) :: E0, E1, E2, E3, E1_VC, E1_VS, E1_TC, E1_TS, E1_sTC, E1_sTS, E1_C, E1_sC
    end type reduced_multipole_t

    ! Physical (ph), isoscalar (is) & isovector (iv) reduced transition intensities ... Julia struct Transition ...
    type, public :: transition_t
        real(dp), allocatable :: ph(:), is(:), iv(:)
    end type transition_t

    ! Julia struct ReducedTransition ... Julia field E1_c (compressional r^3 Y1 strength) is E1_comp here
        ! because Fortran cannot distinguish it from E1_C ...
    type, public :: reduced_transition_t
        type(transition_t) :: E0, E1, E2, E3, E1_V, E1_VC, E1_VS, E1_T, E1_TC, E1_TS, E1_C, E1_CC, E1_CS, E1_comp, E1_NLO_LWA
    end type reduced_transition_t

    ! HF-RRPA auxiliary structure ...
    type, public :: obdm_iteration_pairs_t
        integer :: N = 0
        integer, allocatable :: T(:), a(:), b(:)
    end type obdm_iteration_pairs_t

end module nj_types
