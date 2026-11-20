# Dev-port spec, group "common A": basics, structures, operators

Base: Julia NuJuliet `2595e82` (`Source/...`), which the Fortran port `dev-port` mirrors.
Target: Julia dev `f3f7a33` (`NuJuliet_dev/src/...`).
Scope: `src/common/{functions.jl, orb.jl, O1b.jl, O2b.jl, T1b.jl, T2b.jl, rho.jl, rO.jl}`, `src/import_structures.jl`,
`src/import_NuJuliet.jl`. Fortran: `src/core/{nj_functions,nj_orbitals,nj_o1b,nj_o2b,nj_obdm,nj_types,nj_params,nj_utils,nj_linalg}.f90`,
`src/app/nj_input.f90`.

All line numbers refer to the dev tree unless prefixed `base:` or `F:` (Fortran).

---------------------------------------------------------------------------------------------------------------------

## 0. Summary of what matters

| # | Change | Effect on results | Fortran work |
|---|---|---|---|
| 1 | `Orb1B.pO/nO` are now `Float64`; `orbitals_make` fills the last open shell **uniformly with fractional occupation** (equal-filling) instead of leaving it empty | **yes**, for every open-shell nucleus (closed shells unchanged) | type change + rewrite of `orbitals_make`, every `%pO == 1` test in other groups |
| 2 | `orbitals_make` builds orbitals up to `Calc.Nmax` by default (`Int=true` keyword → `Int.Nmax`) | array length of `Orb` changes (first `a_max(Calc.Nmax)` entries identical) | optional argument; V2b/V3b readers use the `Int` variant |
| 3 | New 1p-1h orbital list `orbitals_1p1h_make` (occupation thresholds `Tol_par`, `Tol_hol`) replaces `orbitals_ph_make/ph_list/one_phonon_make` (kept as dead "Old remove soon" code) | defines the TDA/RPA space | new type + routine |
| 4 | New `orbitals_update`, binary `orbitals_export(...,path)` / `orbitals_import`, `OBDM_natural_orbitals` | new files `*/bin/orb.bin` | new routines |
| 5 | `J1b_convective`: two bug fixes (sqrt(l-1)→sqrt(l+1), second terms no longer inside `if l>0`) | transition currents | 10 lines |
| 6 | `rho_radial_summary`: "normalized" radii are now `sqrt(R2)/A^(1/3)` instead of `sqrt(R2/A^(1/3))` | summary files | 7 lines |
| 7 | `O2B_Pandya` type with 4 blocks (`pp, pnnp, pnpn, nn`); new `pnpn` block | qp (QTDA/QRPA) | new type, extend routine |
| 8 | New `O2b_temp` accessor, `O2b_initialize_empty`, `O1b_*` helpers, `parity/triangle/hat`, `hungarian` in-house | none (API) | small helpers |
| 9 | `OBDM_*` renamed `rho_*`; `qpO1B_*`→`qpO1b_*`, `GeneralizedLaguerre`→`generalized_laguerre`, `KroneckerDelta` removed, `fCG/f6j/f9j` wrapped as `cg/w6j/w9j` | none | rename only (case-only renames are no-ops in Fortran) |
| 10 | Parameter structs: `HF_EFA_Parameters`, `RRPA_Parameters` removed; new `HF_MBPT_Parameters`, `TDA_Parameters`, rich `RPA_Parameters`, `pnQTDA_Parameters`; `Reference` fields; HFB `Broy_Qmax`, `Broy_Bmax` removed, `Broy_Amax` 0.9→0.8 | defaults | `nj_params`, `nj_input` rewrite |
| 11 | Output root `IO/` → `io/`, per-solver subfolders `io/<Path>/<Name>/{bin,densities,...}` | layout | `io_path` helper |

---------------------------------------------------------------------------------------------------------------------

## 1. `src/import_NuJuliet.jl` (module layout & exported API)

Counterpart: base `Source/Import_NuJuliet_Module.jl`. Fortran: `src/app/nujuliet.f90` + dispatch modules + Makefile.

Changes:
* `using` list: `Hungarian` and `Mmap` dropped (`hungarian` is now implemented in `common/functions.jl`).
* Include order (common): `import_structures, functions, orb, O1b, T1b, Tr1b, rho, T2b, O2b, V2b, V3b_no2b, qpO2b, qpO2b_40,
  qpO2b_31, qpO2b_22, Hamiltonian` (`Hamiltonian.jl` new, `OBDM.jl`→`rho.jl`). **`rO.jl` (ex `rME.jl`) is still NOT included**,
  so it stays unreachable (it also redefines `rGrad`, which `T2b.jl` defines) → keep it unported.
* Solver blocks: `HF` (+`HF_MBPT`, no `HF_EFA`), `TDA` (new), `RPA` (replaces `HF_RPA`; also `HF_ScRPA`), `BCS`, `HFB`, `QTDA`,
  `pnQTDA` (new), `QRPA`. `HF_RRPA` and `HF_EFA` are gone; `Temp/` is gone.
* All names are still exported by the `names(@__MODULE__)` loop (no explicit API).

Public entry points used by `scripts/sample_script_run.jl` (= the user-facing API to mirror in `&run Solvers`):

| Julia dev call | proposed Fortran `&run` action | base/Fortran now |
|---|---|---|
| `HF(Params)` | `HF` | `HF` (+ `HF_EFA` when `&hf_efa A,Z` set → remove) |
| `TDA(Params)` | `TDA` | part of `HF_RPA` |
| `RPA(Params)` | `RPA` | `HF_RPA` |
| `HF_ScRPA(Params)` (`src/RPA/RPA.jl:72`) | `HF_ScRPA` | (`HF_RRPA` is the removed predecessor) |
| `BCS`, `HFB`, `QTDA`, `QRPA` | same | same |
| `pnQTDA(Params)` | `pnQTDA` | none |
| `TDA_transition_densities(Params; JP, nu_list, File_Name="Phonon_Densities")` | `TDA_transition_densities` | `HF_RPA_transition_densities` |
| `RPA_transition_densities(...)` | `RPA_transition_densities` | `HF_RPA_transition_densities` |
| `TDA/RPA_transition_currents(Params; JP, nu_list, File_Name="Phonon_Currents", M=-1, xN_grid=40, zN_grid=40)` (M=-1 means 0) | `TDA_transition_currents`, `RPA_transition_currents` | `HF_RPA_transition_currents` |
| `V2b_parse(Params,N_max_Out,N_2max_Out,Out_File)`, `V3b_no2b_parse(...)` | (V2b/V3b group decides) | not ported |

Fortran-only actions (`EMPM_*`, `NCSM`) stay.

Fortran: keep `nujuliet.f90` structure; the dispatch names change (owned by the solver groups). Old names may be kept as
aliases with a deprecation message (decision for the owner, see open questions).

---------------------------------------------------------------------------------------------------------------------

## 2. `src/import_structures.jl` → `nj_types.f90`, `nj_params.f90`

### 2.1 Parameter structs (complete, dev)

```julia
Interaction_Parameters: NN_File="" NNN_File="" hw=0.0 Nmax=0 N2max=0 N3max=0 cP2N=1.0 cP3N=1.0 cRes=1.0 pG2N=0.0 nG2N=0.0   # unchanged
HF_MBPT_Parameters:     Export=false TBCF=false ScH=false                                                                 # NEW
HF_Parameters:          Tol=1e-7 Imax=100 V_HF=false BMF=true HRF=false ESH=false MBPT=HF_MBPT_Parameters()               # EFA removed; V_HF, ESH, MBPT new
TDA_Parameters:         Name="TDA" Reference="HF" Reference_TBCF=false Ortho=true Tol_hol=0.7 Tol_par=0.1                 # NEW
RPA_Parameters:         Name="RPA" Reference="HF" Reference_TBCF=false RPA_Orb=false Ortho=true Transitions="Linear"
                        Tol=1e-5 Tol_hol=0.7 Tol_par=0.1 Imax=150 OBDM_diagonal=false OBDM_alpha=0.6 ScH=true
                        TBCF=false TBCF_alpha=0.8 TBCF_beta=5.0 DEV=false                                                 # was only Ortho=true
(RRPA_Parameters removed)
BCS_Parameters:         ScBCS=false Tol=1e-7 Imax=500 q=0.05 A0=0 Z0=0 pD0=0.5 nD0=0.5 BMF=true HRF=false                 # unchanged
HFB_Parameters:         Tol=1e-7 Imax=150 Pairing="Full" dLmax=0.5 pL0=-3.0 nL0=-3.0 pa0=0.25 na0=0.25 Broy_m=8
                        Broy_Amax=0.8 Broy_Tmax=0.9 LNT=false LNRes=true BMF=true HRF=false                               # Broy_Qmax, Broy_Bmax removed; Broy_Amax 0.9 -> 0.8
QTDA_Parameters:        Ortho=true Reference="HFB"                                                                        # Reference new
QRPA_Parameters:        Ortho=true Reference="HFB"                                                                        # Reference new
pnQTDA_Parameters:      Reference="HFB"                                                                                   # NEW
Calculation_Parameters: A=0 Z=0 hw=0.0 Nmax=0 N2max=0 N3max=0 CMS="CMS1+2B"
                        Path = "A$(A)_Z$(Z)_hw$(hw)_Nmax$(Nmax)_N2max$(N2max)_N3max$(N3max)_$(CMS)"   (unchanged formula)
                        HF TDA RPA BCS HFB QTDA QRPA pnQTDA                                                               # RRPA removed; TDA, pnQTDA added
Parameters(Int, Calc)
```

Where the new fields are used in dev (for the solver groups; semantic owners are those groups):

* `HF.V_HF` → `HF_solver.jl:38` export of the non-local HF radial potential (`HF_radial_potential`, file
  `io/<Path>/HF/densities/HF_radial_potential.dat`).
* `HF.ESH` → `HF_orbitals.jl:6` ("ShellClosure", enforced shell closure in the HF occupation filling).
* `HF.MBPT.Export` (`HF_solver.jl:59`, `HF_MBPT.jl:1521`), `.TBCF` (`HF_MBPT.jl:2252`), `.ScH` (`HF_MBPT.jl:2265`).
* `HF.HRF`, `BCS.HRF`, `HFB.HRF`: **no reader left in dev** (grep `\.HRF` → 0 hits). Keep them in the structs/namelists
  (inputs stay valid) but they have no effect; the HRF exports are removed.
* `TDA.Name`/`RPA.Name`: output subfolder `io/<Path>/<Name>`; `Reference`: input subfolder `io/<Path>/<Reference>`
  (e.g. "HF", "HF-MBPT", "BCS", "HFB"); `Reference_TBCF`: read `io/<Path>/<Reference>/bin/sigma.bin`.
* `Tol_hol`, `Tol_par`: passed to `orbitals_1p1h_make` (§3.2) and `RPA_sigma.jl:6`.
* `RPA.Transitions ∈ {"Exact","Linear"}` (error otherwise, `RPA_external.jl:123`), `RPA_Orb`, `OBDM_diagonal`,
  `OBDM_alpha`, `ScH`, `TBCF`, `TBCF_alpha`, `TBCF_beta`, `Tol`, `Imax` (ScRPA iteration), `DEV` (developer switch in
  `RPA_energy.jl:40`, `RPA_allocate.jl:176,327`).
* `QTDA/QRPA/pnQTDA.Reference`: input folder of `C.bin`, `U.bin`, `V.bin` (`qpH1b.bin`, `qpH2b.bin` for pnQTDA).

### 2.2 Namelist proposal (`nj_input.f90`)

Fortran-only extensions are kept unchanged: `&calculation CM_Lawson`, `&empm`, `&ncsm`.

```fortran
&interaction  NN_File='' NNN_File='' hw=0 Nmax=0 N2max=0 N3max=0 cP2N=1 cP3N=1 cRes=1 pG2N=0 nG2N=0 /     ! unchanged
&calculation  A=0 Z=0 hw=0 Nmax=0 N2max=0 N3max=0 CMS='CMS1+2B' CM_Lawson=0 Path='' /                   ! unchanged (Path '' -> default_path)
&hf           Tol=1e-7 Imax=100 V_HF=F BMF=T HRF=F ESH=F /                                              ! V_HF, ESH new; HRF accepted, no effect
&hf_mbpt      Export=F TBCF=F ScH=F /                                                                   ! NEW group (Julia HF.MBPT)
&tda          Name='TDA' Reference='HF' Reference_TBCF=F Ortho=T Tol_hol=0.7 Tol_par=0.1 /              ! NEW group
&rpa          Name='RPA' Reference='HF' Reference_TBCF=F RPA_Orb=F Ortho=T Transitions='Linear'
              Tol=1e-5 Tol_hol=0.7 Tol_par=0.1 Imax=150 OBDM_diagonal=F OBDM_alpha=0.6 ScH=T
              TBCF=F TBCF_alpha=0.8 TBCF_beta=5.0 DEV=F /                                               ! extended
&bcs          ScBCS=F Tol=1e-7 Imax=500 q=0.05 A0=0 Z0=0 pD0=0.5 nD0=0.5 BMF=T HRF=F /                  ! unchanged (HRF no effect)
&hfb          Tol=1e-7 Imax=150 Pairing='Full' dLmax=0.5 pL0=-3 nL0=-3 pa0=0.25 na0=0.25 Broy_m=8
              Broy_Amax=0.8 Broy_Tmax=0.9 LNT=F LNRes=T BMF=T HRF=F /                                   ! Broy_Qmax, Broy_Bmax removed; Amax default 0.8
&qtda         Ortho=T Reference='HFB' /
&qrpa         Ortho=T Reference='HFB' /
&pnqtda       Reference='HFB' /                                                                         ! NEW group
&run          Solvers='HF','TDA','RPA',... /
&transition_densities JP='1-' nu_list=... File_Name='Phonon_Densities' /                                ! unchanged
&transition_currents  JP='1-' M=0 nu_list=... File_Name='Phonon_Currents' xN_grid=40 zN_grid=40 /       ! unchanged (Julia M=-1 == 0)
```

Removed groups: `&hf_efa`, `&rrpa`. Reading them must not silently succeed with stale meaning; proposal: if the file
contains `&hf_efa` or `&rrpa`, stop with "removed in NuJuliet dev" (detect by a raw text scan, because a missing namelist
group is not an error).

Removed `&hfb` items `Broy_Qmax`, `Broy_Bmax`: a namelist read of an unknown item fails with an iostat error, which is
the desired behaviour (clear message) — or keep them as accepted-and-ignored dummies (owner decision).

**Implementation constraint (name clashes).** `read_input` currently declares one local variable per namelist item and
shares them between groups. The new items create type clashes:
`Transitions` is `logical` in `&ncsm` and `character` in `&rpa`; `Name`/`Reference` are fine; `Tol`, `Imax`, `Ortho`,
`TBCF`, `ScH` are shared with the same type. Restructure `read_input` into one internal (contained) subroutine per group,
each with its own local variables, defaults taken from a default-initialized `*_params_t`, e.g.

```fortran
subroutine read_rpa(u, P)          ! contained in read_input; locals shadow nothing in the host
    type(rpa_params_t), intent(inout) :: P
    character(len=64) :: Name, Reference, Transitions
    logical :: Reference_TBCF, RPA_Orb, Ortho, OBDM_diagonal, ScH, TBCF, DEV
    real(dp) :: Tol, Tol_hol, Tol_par, OBDM_alpha, TBCF_alpha, TBCF_beta
    integer :: Imax, Ios
    namelist /rpa/ Name, Reference, Reference_TBCF, RPA_Orb, Ortho, Transitions, Tol, Tol_hol, Tol_par, Imax, &
                   OBDM_diagonal, OBDM_alpha, ScH, TBCF, TBCF_alpha, TBCF_beta, DEV
    Name = P%Name; ...                 ! copy defaults
    rewind(u); read(u, nml=rpa, iostat=Ios, iomsg=Msg); call check(Ios,'rpa')
    P%Name = Name; ...
end subroutine
```

Also validate `RPA%Transitions ∈ {'Exact','Linear'}`, `TDA%Name /= TDA%Reference`, `RPA%Name /= RPA%Reference`
(dev refuses these, `TDA.jl:44`).

### 2.3 `nj_params.f90` (proposed types)

```fortran
type :: hf_mbpt_params_t;  logical :: Export=.false., TBCF=.false., ScH=.false.; end type
type :: hf_params_t
    real(dp) :: Tol = 1.0e-7_dp;  integer :: Imax = 100
    logical  :: V_HF = .false., BMF = .true., HRF = .false., ESH = .false.
    type(hf_mbpt_params_t) :: MBPT
end type                                                       ! EFA component removed
type :: tda_params_t
    character(len=64) :: Name = 'TDA', Reference = 'HF'
    logical :: Reference_TBCF = .false., Ortho = .true.
    real(dp) :: Tol_hol = 0.7_dp, Tol_par = 0.1_dp
end type
type :: rpa_params_t
    character(len=64) :: Name = 'RPA', Reference = 'HF'
    logical :: Reference_TBCF = .false., RPA_Orb = .false., Ortho = .true.
    character(len=16) :: Transitions = 'Linear'
    real(dp) :: Tol = 1.0e-5_dp, Tol_hol = 0.7_dp, Tol_par = 0.1_dp
    integer :: Imax = 150
    logical :: OBDM_diagonal = .false.;  real(dp) :: OBDM_alpha = 0.6_dp
    logical :: ScH = .true., TBCF = .false.
    real(dp) :: TBCF_alpha = 0.8_dp, TBCF_beta = 5.0_dp
    logical :: DEV = .false.
end type
type :: hfb_params_t   ! remove Broy_Qmax, Broy_Bmax; Broy_Amax = 0.8_dp
type :: qtda_params_t;   logical :: Ortho = .true.; character(len=64) :: Reference = 'HFB'; end type
type :: qrpa_params_t;   logical :: Ortho = .true.; character(len=64) :: Reference = 'HFB'; end type
type :: pnqtda_params_t; character(len=64) :: Reference = 'HFB'; end type
type :: calculation_params_t   ! components: HF, TDA, RPA, BCS, HFB, QTDA, QRPA, pnQTDA, EMPM, NCSM (RRPA removed), CM_Lawson kept
```

Delete `hf_efa_params_t`, `rrpa_params_t`. Note: `nj_input` builds `rrpa_params_t(...)`, `bcs_params_t(...)`,
`hfb_params_t(...)` with positional structure constructors — every one of those must be updated when fields change
(positional constructors silently mis-assign if two adjacent fields have the same type: use keyword constructors or
component assignment).

Path helpers (shared by all groups):

```fortran
function io_path(Params,Suffix)      ! 'io/' // trim(Params%Calc%Path) // Suffix      (was 'IO/')
function io_dir(Params,Sub)          ! 'io/' // trim(Path) // '/' // trim(Sub)  e.g. io_dir(P,P%Calc%RPA%Name)
```

Dev creates directories with plain `mkdir` (fails if the parent does not exist); the port's `make_dir` (mkdir -p) is
fine. Fortran-only consumers of the old layout (`src/empm`, `src/ncsm`) read `Bin/C_HF.bin`, `Bin/h_HF.bin`,
`Bin/Rho_HF.bin`, `Bin/V2B_HF.bin`, `HF/HF_Summary.dat`; dev HF writes `io/<Path>/HF/bin/{orb,C,rho,H1b,H2b}.bin` and
`io/<Path>/HF/HF_summary.dat` (lower-case "summary"). The mapping and the file contents/bases belong to the HF group's
spec; EMPM/NCSM must be adapted after that (see open questions).

### 2.4 Structure types (`nj_types.f90`)

| dev struct | base | Fortran now | required change |
|---|---|---|---|
| `Orb1B(a,n,l,j,pO::Float64,nO::Float64)` | `pO,nO::Int64` | `orb1b_t` integer `pO,nO` | **`real(dp) :: pO = 0, nO = 0`**. All users of `%pO == 1` (`nj_hf.f90:252,254,516,520`, `nj_bcs.f90:532-539,893,898`, `nj_hf_rpa.f90:552`, RRPA files) change with their groups; `istr(Orb%pO)` in `orbitals_export` (HRF) disappears. |
| `SPState`, `pnSVector`, `pnInteger`, `pnFloat`, `pnVector`, `pnMatrix`, `pnArray`, `pnCVector`, `PhState` | same | exist | none (`SPState`, `pnSVector` only in the dead "old" orbital code; `PhState` still used by pnQTDA) |
| `pnOrbVector(p,n::Vector{Orb1B})` | new | – | `type :: pnorbvector_t; type(orb1b_t), allocatable :: p(:), n(:)` (used by HF-MBPT) |
| `Ind2B(a,b,T)` | `Ind2qp` | `ind2qp_t` | rename to `ind2b_t` (same fields; T = -1 proton, +1 neutron) |
| `Orb2B_1p1h(N::Matrix{Int}, Ind::Matrix{Vector{Ind2B}}, Dic::Dict{UInt64,Int})` | new | – | `type :: orb2b_1p1h_t; integer, allocatable :: N(:,:) ! (2,0:J_max); type(ind2b_vec_t), allocatable :: Ind(:,:) ! (2,0:J_max); integer, allocatable :: Dic(:,:,:,:) ! (a,b,0:J_max,1:2) T=-1->1, +1->2` |
| `Orb2B_2qp(N, Ind, Dic)` | `qpOrb2B(N, i)` | `qporb2b_t(N, i)` | rename to `orb2b_2qp_t`, field `i`→`Ind`, add dense `Dic(a,b,0:J_max,1:2)` (dev never reads it: `orb2b_2qp_index` has no caller; optional) |
| `Orb2B(Map::Array{Int32,4}, N::Array{Int,3}, Ind::Array{Vector{Vector{Int}},3})` | `Dic::Dict` | `orb2b_t%Dic(a,b,0:J,0:1)` dense | none: the port already did this. Dev's `Map[a,b,J+1,T+1]` = F `Dic(a,b,J,T)`. Optionally rename `Dic`→`Map`. |
| `Orb2B_Temp(Map::Array{Int32,3}, N, Ind)` | `Dic` | `orb2b_temp_t%Dic(a,b,0:J)` dense | none |
| `Orb3B(Map::Vector{Array{Int32,3}}, N::Array{Vector{Int},4}, nlj::Array{Int,3})` | `Dic::Dict, N` | hash map `orb3b_t` | owned by the V3b group (dense per-pair map; see its spec) |
| `O2B_Pandya(pp, pnnp, pnpn, nn :: Matrix{Matrix{Float64}})` | new (base used `O2B_Temp` with `pp, pn, nn`) | `o2b_temp_t` | **new `type :: o2b_pandya_t; type(rmat_t), allocatable :: pp(:,:), pnnp(:,:), pnpn(:,:), nn(:,:)  ! (P,0:J_max)`** |
| `Tr1B`, `qpTr1B` | 15 fields | `nj_tr1b` | two new fields `F`, `GT` (Fermi, Gamow-Teller; pnQTDA) — owned by the Tr1b group |
| `Orb1B_HF_MBPT(p::pnOrbVector, h::pnOrbVector, N_p::pnInteger, N_h::pnInteger)` | new | – | HF-MBPT group (shared type, define in `nj_types`) |
| `Orb_HF_MBPT3(a,l,j,t,E)` | new | – | defined but unused in dev → skip |
| `OBDM_Iteration_Pairs` | same | `obdm_iteration_pairs_t` | unused in dev (RRPA gone) → delete with RRPA |
| `HFB_Broyden*` | same | not ported (port note B.1) | none |
| `qpO40/31/22/qpO2B/O2B_31_Temp/O2B/O2B_Temp/O1B/qpO1B/ReducedMultipole/Transition/ReducedTransition` | same | exist | none |
| `Orb1B_RPA_OBDM` (defined in `src/RPA/RPA_orbitals.jl:1`, not in import_structures) | new | – | RPA group |

---------------------------------------------------------------------------------------------------------------------

## 3. `src/common/orb.jl` (base `MatrixElements/Orb.jl`) → `nj_orbitals.f90`

### 3.1 `orbitals_make(Params; Int::Bool=false)` (orb.jl:1-89) — F `orbitals_make`

Changes vs base:
1. `N_max = Params.Calc.Nmax`, or `Params.Int.Nmax` when `Int=true` (base: always `Params.Int.Nmax`, base Orb.jl:3).
   Callers with `Int=true`: `V2b.jl:16,326`, `V3b_no2b.jl:134,396` only. Every solver now gets an `Orb` of length
   `a_max(Calc.Nmax)`, and dev code uses `length(Orb)` as a_max in places (`HF_density_operator.jl:3`,
   `radial_moment_matrix_LHO`), so the length matters.
2. Occupations are `Float64` with equal filling of the first shell that does not fit. Exact dev code (orb.jl:31-79):

```julia
if abs(A) > 1e-8 || abs(Z) > 1e-8          # base: A != 0 && Z != 0
    Z_occ, N_occ = 0.0, 0.0;  Z_fill, N_fill = true, true;  Z_mask, N_mask = falses(a_max), falses(a_max)
    while N_fill
        for N = 0:N_max, n in 0:div(N,2)            # l = N - 2n (largest l first)
            l = N - 2*n
            for j = (2*l + 1):-2:abs(2*l-1)          # j = l+1/2 first
                j_hat = Float64(j + 1)
                for a = 1:a_max
                    if Orb_arr[a,2]==n && Orb_arr[a,3]==l && Orb_arr[a,4]==j && N_mask[a]==false
                        Occ = min(1.0, (A - Z - N_occ) / j_hat)
                        Obc_occ[a,2] = Occ;  N_occ += j_hat * Occ;  N_mask[a] = true
                        if abs(A - Z - N_occ) < 1e-8;  N_fill = false;  end
    ...  # identical block for protons with Z, Z_occ, Z_mask, column 1
```

   Notes for an exact port:
   * `A`, `Z` are `Float64`; keep the arithmetic order (`(A - Z - N_occ) / j_hat`, `N_occ += j_hat*Occ`) so that the
     tiny residual values (±1e-16) given to the orbitals after the last partially filled one are reproduced bit for bit
     (they are not exactly zero in general; `min(1.0, x)` with x ≈ ±eps).
   * After `N_fill = false` the `for` loops still run to the end of the pass; all later orbitals get
     `Occ = min(1, residual/j_hat)` ≈ 0. Every orbital is visited exactly once (mask), so a second `while` pass only
     happens if the whole basis is filled without reaching A−Z → **infinite loop** in dev if A−Z (or Z) exceeds the
     basis capacity. The port should `fatal()` instead (deviation).
   * Behavioural difference vs base for open shells (e.g. 18O, N=10): base filled whole j-shells only and — because
     it kept looping after the first shell that did not fit — could fill a later, smaller shell (base filled 1s1/2
     before 0d5/2 for N=10); dev puts 2/6 into 0d5/2. Closed-shell nuclei: identical (0.0/1.0).
   * HF itself re-derives the occupations from the HF energies in `HF_orbitals_make` (`src/HF/HF_orbitals.jl`, HF group,
     uses `HF.ESH`); this is the replacement of HF-EFA.
3. Fortran signature proposal: `function orbitals_make(Params, Use_Int) result(Orb)` with `logical, optional ::
   Use_Int` (do not call the dummy `Int`: it shadows the intrinsic `int`). V2b/V3b readers call it with
   `Use_Int=.true.` internally (dev `V2b_read(Params)` no longer receives `Orb`; V2b group).
4. Banner text changed: "Constructing Linear-Harmonic Oscillator basis single-particle orbitals ..." /
   "Linear-Harmonic Oscillator basis single-particle orbitals initialized ...".

Size: ~70 lines rewritten.

### 3.2 `orbitals_1p1h_make(Params, Orb; Tol_par=1e-4, Tol_hol=0.999)` (orb.jl:333-442) — new; F: none

Replaces `orbitals_ph_make` + `orbitals_ph_list` + `orbitals_one_phonon_make` for TDA/RPA (those three are still in
orb.jl:257-330 under "Old remove soon", `orbitals_ph_make` itself is deleted; none has a caller in dev).
Callers always pass `Tol_par = Params.Calc.{TDA,RPA}.Tol_par`, `Tol_hol = ...Tol_hol` (defaults 0.1 / 0.7).

```julia
J_max = 2*N_max + 1          # N_max = Calc.Nmax
for J in 0:J_max, P in 1:2                     # (count pass and fill pass use the same order)
    for Orb_a in Orb                           # a = particle candidate, all orbitals
        for Orb_b in Orb                       # b = hole candidate, all orbitals
            if parity(l_a + l_b) != P || !triangle(j_a, j_b, 2J); continue; end
            if pO_b > Tol_hol && pO_a < Tol_par && (pO_b - pO_a) > 1e-8      # proton pair first
                N += 1; Ind[P,J+1][N] = Ind2B(a,b,-1); Dic[orb2b_1p1h_key(a,b,J,P,-1)] = N
            if nO_b > Tol_hol && nO_a < Tol_par && (nO_b - nO_a) > 1e-8      # then neutron pair
                N += 1; Ind[P,J+1][N] = Ind2B(a,b,+1); Dic[orb2b_1p1h_key(a,b,J,P,+1)] = N
```

* Storage: `N[P,J+1]` (`N_ph` is `zeros(Int64,2,J_max+1)`), `Ind[P,J+1]::Vector{Ind2B}`; ordering inside a channel:
  particle index a outer, hole index b inner, proton entry before neutron entry for the same (a,b). Note the order
  differs from the base phonon ordering (base: all proton (p,h) pairs first, then all neutron pairs, p over the list
  of particles).
* Key `orb2b_1p1h_key(a,b,J,P,T) = P | (T+1)<<2 | J<<4 | a<<10 | b<<18` (a, b < 256, J < 64);
  `orb2b_1p1h_index(a,b,J,P,T,Orb_1p1h)` returns `get(Dic,key,0)` — **0 if absent (no error)**.
  Fortran: dense `Dic(a,b,0:J_max,iT)` with `iT = (T+3)/2` (P is implied by l_a+l_b), 0 = absent;
  `pure integer function orb2b_1p1h_index(a,b,J,P,T,Orb_1p1h)`.
* With fractional occupations an open shell (e.g. occupation 1/3) is neither particle (`< Tol_par`) nor hole
  (`> Tol_hol`) with the default tolerances → excluded from the 1p1h space. This is dev behaviour; document it.
* `parity(L) = rem(L,2) + 1`, `triangle(j1,j2,j3) = |j1-j2| <= j3 <= j1+j2` (doubled j; no parity-of-sum test).

Size: ~90 lines (type + routine + index function).

### 3.3 `orbitals_2qp_make` (orb.jl:444-549) — F `orbitals_2qp_make`

Same counting/ordering as base (pp block a=1..a_max, b=1..a, J range, `a≠b or J even`; then nn block); `P = parity(...)`
identical to base `rem(l_a+l_b,2)+1`. New: returns `Orb2B_2qp(N, Ind, Dic)` with `Ind2B` entries and key
`orb2b_2qp_key` (same packing as 1p1h). Fortran: rename type/field (`i`→`Ind`) and optionally fill `Dic`. ~15 lines.

### 3.4 `orbitals_update(Orb, Rho::O1B)` (orb.jl:91-109) — new

Sets `pO_a = Rho.p[a,a]`, `nO_a = Rho.n[a,a]` (by `Orb_a.a`), mutates and returns `Orb`.
(`a_max = size(Orb)` is a tuple; harmless in Julia.) Fortran: `subroutine orbitals_update(Orb,Rho)` with
`intent(inout)`, ~12 lines. Callers: `RPA_solver.jl:400`, `RPA_export.jl:141`.

### 3.5 Binary orbital table: `orbitals_export(Params, Orb, export_path)` / `orbitals_import(Params, import_path)` (orb.jl:111-156) — new

Layout (stream, little endian, no header), for i = 1..a_max(Calc.Nmax):
`Int64 a, Int64 n, Int64 l, Int64 j, Float64 pO, Float64 nO` (48 bytes per orbital).
Files: `io/<Path>/HF/bin/orb.bin`, `.../HF/bin/orb_temp.bin`, `.../HF-MBPT/bin/orb.bin`, `io/<Path>/<RPA.Name>/bin/orb.bin`
(writers/readers belong to the solver groups). Fortran: `integer(i8)` writes, `access='stream'`; ~35 lines.
The Fortran name `orbitals_export` is taken by the HRF text export — see §3.6.

### 3.6 `orbitals_export(Params, Orb)` (HRF, orb.jl:196-255) and `orbitals_lj_make`

* HRF `orbitals_export(Params,Orb)` is unchanged but has **no caller** in dev (HRF switches are dead). It still reads
  the never-written `Bin/pE_HF.bin` (port note 1). Proposal: drop the Fortran HRF version; give the binary version the
  generic name `orbitals_export(Params,Orb,Path)` (or `orbitals_export_bin`). Port note 1 becomes obsolete.
* `orbitals_lj_make` unchanged (now also used by `O2b_transformation_index1..4` instead of `Orb_PreComp`).

### 3.7 `OBDM_natural_orbitals(Params, Orb, Rho::O1B)` (orb.jl:551-650) — new; F: closest is the RRPA `reorder` in `nj_hf_rrpa.f90:806`

Returns `(Orb_new, C::O1B, Rho_diag::O1B)`. Exact steps:
1. `pRho = 0.5(Rho.p + Rho.p')`, same for n.
2. Zero every element with `Orb[a].j != Orb[b].j || Orb[a].l != Orb[b].l`.
3. `pRho[a,a] += 1e-10*a`, `nRho[a,a] += 1e-10*a` (deterministic splitting; the shift stays in the returned occupations).
4. `pOcc, pC = eigen(Symmetric(pRho))` (ascending; F `eig_sym`), same for n.
5. `pReordering, _ = hungarian(ones(a_max,a_max) .- abs.(pC))`; `pC = pC[:,pReordering]`, `pOcc = pOcc[pReordering]`;
   same for n (F: existing `hungarian(Cost,Assign)` — same JV algorithm and tie-breaking, see §4.3).
6. Sign fix: `if pC[a,a] < 0: pC[:,a] *= -1` (same for n).
7. `Rho_diag = diagm(pOcc), diagm(nOcc)`.
8. New orbitals: for each a: `l_a = round(Int, Σ_b |pC[b,a]|^2 l_b)`, `j_a = round(Int, Σ_b |pC[b,a]|^2 j_b)`
   (**proton C is used for both isospins**); `n_a = #{b < a : Orb[b].l == l_a && Orb[b].j == j_a}` (counted on the old
   `Orb` list); `Orb_new[a] = Orb1B(a, n_a, l_a, j_a, pOcc[a], nOcc[a])`.
   Fortran: `nint` is fine (values are near-integers).
Callers: `HF_MBPT.jl:1510`, `RPA_export.jl:69`, `RPA_rho.jl:555`. Size ~70 lines. Put it in `nj_orbitals` (or `nj_obdm`).

### 3.8 Dead code to drop in Fortran
`orbitals_ph_make` (deleted in dev), `orbitals_ph_list`, `orbitals_one_phonon_make` (dead in dev), `read_h_HF_diag`
(only used by these and HRF export). Remove once the TDA/RPA/RRPA groups no longer call them (~ -150 lines).

---------------------------------------------------------------------------------------------------------------------

## 4. `src/common/functions.jl` (base `Functions/Functions.jl`, CRLF → LF) → `nj_functions`, `nj_utils`, `nj_linalg`

### 4.1 Additions / removals / renames

| dev | base | Fortran action |
|---|---|---|
| `parity(L) = rem(L,2)+1` (:9) | new | add `pure integer function parity_of(L)` in `nj_utils` (**not** `parity`: `Parity` is a local/namelist variable in `nj_input`/`nj_ncsm`, and a use-associated name clashes with a local declaration) |
| `cg`, `w3j`, `w6j`, `w9j` (:16-30) = `fCG`, `f3j`, `f6j`, `f9j` | new wrappers | none: keep calling `fCG/f6j/f9j` (`CG` is also a local variable name in `nj_functions`); `w3j` has no caller |
| `triangle(j1,j2,j3)` (:32) | new | add `pure logical function triangle(j1,j2,j3)` (check no local `triangle` variables: none found) |
| `hat(j) = Float64(j+1)` (:40) | new | add `pure real(dp) function hat(j)`; note the argument is the doubled j (`hat(2*I) = 2I+1`) |
| `KroneckerDelta` | removed (kept `kronecker_delta`) | none |
| `GeneralizedLaguerre` → `generalized_laguerre` | rename | none (F already lower case) |
| `radial_moment_matrix_LHO`: local `HbarC` → `hc` | cosmetic | none |
| `radial_phase_test` | removed | none (not ported) |
| `hungarian(costMat)` (:711-809) | new in-house (was package Hungarian.jl) | F `nj_linalg::hungarian` already implements the same algorithm; see 4.3 |
| `_assigned_vars!` / `@quiet` | commented out | none |
| `phonon_weights!(M)` (:860) | new | add `subroutine phonon_weights(M)` (complex or real vector): `Norm = sqrt(sum |M|^2)`; `M = |M|/Norm` if Norm>0 else `|M|` (used 20x by TDA/RPA external routines) |
| `logm`, `lorentzian`, `integrate_quadrature`, `legendre_polynomial`, `pretty_summarysize` | unchanged, still no caller | not ported |
| `JP_initialize`, `differentiate_Oh6`, `Psi_rad_LHO`, `integrate_trap`, `radial_moment_LHO`, `doublefactorial`, `associated_legendre_polynomial`, `normalize_spherical_angles`, `scalar_/vector_spherical_harmonic_function`, `J1b_spin` | unchanged apart from `fCG/f9j`→`cg/w9j` | none |

### 4.2 `J1b_convective` bug fixes (functions.jl:504-623) — F `nj_functions.f90:282-335`

In each of the four brackets the "2nd term" (the `Y_{l, l+1}` vector harmonic) moved **out of** the `if l_a > 0` /
`if l_b > 0` block, `Y` is recomputed before it, and the factor `sqrt(l-1)` became `sqrt(l+1)`:

```julia
# functions.jl:529-537 (bra, first bracket, 2nd term; same pattern at :555, :586, :611)
Y = scalar_spherical_harmonic_function(l_b,div(m_b-1,2),theta,phi)
(Y_x,Y_y,Y_z) = vector_spherical_harmonic_function(l_a,l_a+1,-div(m_a-1,2),theta,phi)
J = -sqrt(Float64(l_a + 1) / Float64(2*l_a + 1)) * (dR_a - Float64(l_a) / r * R_a) * R_b * Y * Amp_1 * phase(div(m_a-1,2))
# second bracket, 2nd term (:555-559)
Y = phase(div(m_a-1,2)) * scalar_spherical_harmonic_function(l_a,-div(m_a-1,2),theta,phi)
(Y_x,Y_y,Y_z) = vector_spherical_harmonic_function(l_b,l_b+1,div(m_b-1,2),theta,phi)
J = sqrt(Float64(l_b + 1) / Float64(2*l_b + 1)) * (dR_b - Float64(l_b) / r * R_b) * R_a * Y * Amp_1
```

(base had `sqrt(l-1)` — NaN for l=0 was masked by the `if l>0`; for l≥1 the old result was simply wrong).
Fortran `add_terms`: move the two `l+1` terms after the respective `end if`, recompute `Y` before each, replace
`real(l_a - 1,dp)`/`real(l_b - 1,dp)` by `real(l_a + 1,dp)`/`real(l_b + 1,dp)`. The 1st terms (`Y_{l,l-1}`) stay
inside the `if`. ~10 lines. Affects `*_transition_currents` output only.

### 4.3 `hungarian`
Dev algorithm = shortest augmenting path with potentials (JV), identical loop structure and strict `<` tie-breaking
to F `nj_linalg::hungarian` (square case). Dev additionally supports rectangular matrices (transposes if rows > cols,
`assignment[n] = 0` for unassigned rows) and maps non-finite entries to a penalty
`c_max + n*(c_max - c_min) + 1`. All dev callers (`OBDM_natural_orbitals`, `HFB_ordering.jl:149`) pass square finite
matrices → no Fortran change needed; optionally add a `size(Cost,1) == size(Cost,2)` check.

---------------------------------------------------------------------------------------------------------------------

## 5. `src/common/O1b.jl` (base `O1B.jl`) → `nj_o1b.f90`

New helpers (O1b.jl:1-39), add to `nj_o1b` (~50 lines):

```julia
O1b_transformation(O,C) = O1B(C.p' * O.p * C.p, C.n' * O.n * C.n)      # note: C' O C
O1b_transpose(O)        = O1B(transpose(O.p), transpose(O.n))
O1b_diagonal(O)         = O1B(diagm(diag(O.p)), diagm(diag(O.n)))
O1b_deepcopy(O)
O1b_zeroes(Params), O1b_ones(Params), O1b_unit_matrix(Params)          # a_max(Calc.Nmax) square
```
Fortran: `function O1b_transformation(O,C) result(R)` via `mtmm`; others trivial (`O1b_deepcopy` = assignment).
Renames `qpO1B_export/import`→`qpO1b_*`, `qpN1B_make`→`qpN1b_make`: case-only, nothing to do in Fortran.
`O1b_export/O1b_import` file layout unchanged (per isospin, (a,b) with equal (l,j), a outer).

## 6. `src/common/T1b.jl`, `T2b.jl` → `nj_o1b.f90`
`T1b.jl`: byte-identical. `T2b.jl`: `f6j`→`w6j` only. No Fortran change.

## 7. `src/common/rO.jl` (base `rME.jl`)
`f6j/fCG`→`w6j/cg` only, and still **not included** by the module. Not ported (unchanged decision).

## 8. `src/common/rho.jl` (base `OBDM.jl`) → `nj_obdm.f90`

1. Renames: `OBDM_export`→`rho_export`, `OBDM_radial_cms`→`rho_radial_cms`, `OBDM_radial_density_grid`→
   `rho_radial_density_grid`, `OBDM_radial_chdensity_grid`→`rho_radial_chdensity_grid`, `OBDM_radial_radii`→
   `rho_radial_radii`, `OBDM_radial_density_export`→`rho_radial_density_export`, `OBDM_radial_summary`→
   `rho_radial_summary`, `OBDM_rN`→`rho_rN`, `OBDM_chR2`→`rho_chR2`. Signatures unchanged. Fortran: rename (module
   `nj_obdm` → `nj_rho` optional); keep old names as thin aliases only if other groups lag behind.
2. **Formula change** in the "normalized by the scale factor 1/A^1/3" block of `rho_radial_summary` (rho.jl:276-282):
   `sqrt(abs(X)) / A^(1/3)` instead of base `sqrt(abs(X) / A^(1/3))`, e.g.
   ```julia
   @printf(Summary, "\nr_p   = %12.6f fm ...", sqrt(abs(pR2 + R_CMS[1] + R_CMS[2] + 0.5 * (197.326980 / 938.272013)^2)) / A^(1/3))
   ```
   Fortran `OBDM_radial_summary` (F nj_obdm.f90:296-309): in Pass 2 replace `sqrt(abs(X) / A13)` by
   `sqrt(abs(X)) / A13` for all 7 lines (Pass 1 has `A13 = 1`, unchanged). 7 lines.
3. Everything else (CMS formulas, grids `2^13+1` points, folding, `rho_chR2` passing `Rho_grid[2]` as charge density)
   unchanged.

## 9. `src/common/O2b.jl` (base `O2B.jl`) → `nj_o2b.f90`

1. Dense lookup `Orb2B.Map[a,b,J+1,T+1]` (`Int32`) instead of the `Dict` — the port already uses the dense table
   `Orb_NN%Dic(a,b,J,T)`. Dev `O2b_lookup` throws `KeyError` for out-of-range or 0 entries; F `O2b_index` returns
   garbage silently. Optional: a debug-build check.
2. `O2b_initialize_empty()` (O2b.jl:104): O2B with 1x1 undefined channel matrices; used as a "no Sigma" placeholder
   (`RPA_solver.jl`, `TDA_external.jl`). Fortran: an unallocated `o2b_t` plus `allocated()` tests — no routine needed
   (document for the TDA/RPA groups).
3. `O2b_pp/O2b_nn`: `Sign_ab = 1` → `1.0` (type stability only; no numeric change).
4. `O2b_export`/`O2b_import`: block-buffered I/O, **identical file layout** (pp lower triangles, pn full N×N squares
   with packed index aliasing, nn lower triangles; J = 0..2*Nmax+1, P = 1,2). Still loops `J_max = 2*N_max+1` while
   the arrays have `N2max+2` channels → **port note 3 still needed** (bug not fixed upstream).
5. `O2b_temp_index` now throws on missing entries (F returns 0). New accessor (O2b.jl:601):
   ```julia
   O2b_temp(a,b,c,d,J,P,T,O::O2B_Temp,Orb_NN_temp) = (T==-1 ? O.pp : T==0 ? O.pn : O.nn)[P,J+1][Bra,Ket]
   ```
   with `Bra = O2b_temp_index(a,b,J,P)`, `Ket = O2b_temp_index(c,d,J,P)`. **T convention here: −1 pp, 0 pn, +1 nn**
   (unlike `O2b_index`, where T=0 pn and T=1 like-particles). Used heavily by `TDA_allocate.jl`, `RPA_allocate.jl`.
   Fortran: `pure real(dp) function O2b_temp(a,b,c,d,J,P,T,O,Orb_NN_t)` in `nj_o2b` (~15 lines).
6. `O2b_transformation(...; Report=true)`: `Report=false` suppresses all prints and the `@time` lines were removed.
   Fortran: optional `logical :: Report` (default `.true.`), or keep the global `Verbose` switch (port note 8). The
   4 index steps now use `orbitals_lj_make` instead of `Orb_PreComp` (removed): no numeric change, the summation order
   over `Orb_x` (ascending orbital index) is the same. Fortran: may replace `Orb_PreComp` likewise (performance only).
7. **Pandya**: new `O2b_Pandya_initialize` (same channel/orbital structure as `O2b_temp_initialize`: a with
   2n_a+l_a ≤ Nmax, all b with 2(n_a+n_b)+l_a+l_b ≤ N2max, J_max = N2max+1) returning `O2B_Pandya` + `Orb2B_Temp`, and
   accessors `O2b_Pandya_pp/pnpn/pnnp/nn(a,b,c,d,J,P,F,Orb_NN_temp)`. `O2b_Pandya_transformation` (O2b.jl:1036) now
   fills four blocks (with `hat(2I) = 2I+1`, `phase`, `w6j = f6j`):
   ```julia
   ppF   = Σ_I phase(div(j_b+j_c,2)+I) * hat(2I) * w6j(j_a,j_b,2J,j_c,j_d,2I) * O2b_pp(a,d,b,c,I,Q,...)
   pnpnF = Σ_I phase(1)                * hat(2I) * w6j(j_a,j_b,2J,j_c,j_d,2I) * O2b_pn(a,d,c,b,I,Q,...)   # NEW: sign -1, ket (c,b)
   pnnpF = Σ_I phase(div(j_b+j_c,2)+I) * hat(2I) * w6j(j_a,j_b,2J,j_c,j_d,2I) * O2b_pn(a,d,b,c,I,Q,...)   # = base "pn"
   nnF   = Σ_I phase(div(j_b+j_c,2)+I) * hat(2I) * w6j(j_a,j_b,2J,j_c,j_d,2I) * O2b_nn(a,d,b,c,I,Q,...)
   # Q = parity(l_a+l_d), requires Q == parity(l_b+l_c); I from max(|j_a-j_d|,|j_b-j_c|)/2 to min(j_a+j_d,j_b+j_c)/2
   ```
   (O2b.jl:1061-1074). Loop: `for JP in JP_initialize(N2max+1)` (J outer, P inner), threads over Bra.
   Fortran: new type `o2b_pandya_t`, `O2b_Pandya_transformation(Params,Orb,Orb_NN,O_NN,F_NN)` with
   `type(o2b_pandya_t), intent(out) :: F_NN`, `pandya` kinds pp / pnpn / pnnp / nn; the existing `pn` result
   becomes `pnnp`. Caller `nj_qpo2b.f90:85` and all qp users of `F_NN%pn` must switch to `%pnnp` (and use `%pnpn`
   where `qpO2b_22/40` dev uses it — qp group). ~40 lines.

---------------------------------------------------------------------------------------------------------------------

## 10. Shared-module interface changes (one line each)

1. `nj_types`: `orb1b_t%pO`, `%nO` → `real(dp)` (default 0).
2. `nj_types`: rename `ind2qp_t`/`ind2qp_vec_t` → `ind2b_t`/`ind2b_vec_t` (fields a, b, T unchanged).
3. `nj_types`: new `orb2b_1p1h_t { integer N(2,0:J_max); type(ind2b_vec_t) Ind(2,0:J_max); integer Dic(a_max,a_max,0:J_max,2) }`.
4. `nj_types`: `qporb2b_t` → `orb2b_2qp_t`, component `i` → `Ind`, new optional `Dic(a_max,a_max,0:J_max,2)`.
5. `nj_types`: new `o2b_pandya_t { type(rmat_t), allocatable :: pp(:,:), pnnp(:,:), pnpn(:,:), nn(:,:) }` indexed (P,0:J_max).
6. `nj_types`: new `pnorbvector_t { type(orb1b_t), allocatable :: p(:), n(:) }` and `orb1b_hf_mbpt_t { type(pnorbvector_t) :: p, h; type(pninteger_t) :: N_p, N_h }` (HF-MBPT group uses them).
7. `nj_types`: delete `obdm_iteration_pairs_t` together with RRPA; `pnarray4_t`, `spstate_t`, `pnsvector_t` become dead once the old ph-orbital routines are gone.
8. `nj_params`: delete `hf_efa_params_t`, `rrpa_params_t`, `calculation_params_t%RRPA`, `hf_params_t%EFA`.
9. `nj_params`: `hf_params_t` add `V_HF=.false.`, `ESH=.false.`, `MBPT` (`hf_mbpt_params_t`: Export, TBCF, ScH all `.false.`).
10. `nj_params`: new `tda_params_t` (Name='TDA', Reference='HF', Reference_TBCF=F, Ortho=T, Tol_hol=0.7, Tol_par=0.1) and `calculation_params_t%TDA`.
11. `nj_params`: `rpa_params_t` extended to the 17 dev fields (§2.3).
12. `nj_params`: `hfb_params_t` drop `Broy_Qmax`, `Broy_Bmax`; `Broy_Amax` default 0.8.
13. `nj_params`: `qtda_params_t`, `qrpa_params_t` add `Reference='HFB'`; new `pnqtda_params_t(Reference='HFB')` and `calculation_params_t%pnQTDA`.
14. `nj_params`: `io_path` prefix `'IO/'` → `'io/'`; new `io_dir(Params,Sub)`.
15. `nj_input`: groups `&hf_mbpt`, `&tda`, `&pnqtda` new; `&rpa`, `&qtda`, `&qrpa`, `&hf`, `&hfb` changed; `&hf_efa`, `&rrpa` removed; restructure into per-group contained readers (type clash `Transitions`).
16. `nj_orbitals`: `orbitals_make(Params[,Use_Int])` (Calc.Nmax by default, fractional filling); new `orbitals_update`, `orbitals_export(Params,Orb,Path)` (binary), `orbitals_import(Params,Path)`, `orbitals_1p1h_make(Params,Orb,Tol_par,Tol_hol)`, `orb2b_1p1h_index`, `OBDM_natural_orbitals(Params,Orb,Rho,Orb_new,C,Rho_diag)`; `orbitals_2qp_make` returns `orb2b_2qp_t`; remove `orbitals_ph_make`, `orbitals_ph_list`, `orbitals_one_phonon_make`, HRF `orbitals_export(Params,Orb)`.
17. `nj_o1b`: add `O1b_transformation(O,C)` (= CᵀOC), `O1b_transpose`, `O1b_diagonal`, `O1b_zeroes/ones/unit_matrix(Params)`.
18. `nj_o2b`: add `O2b_temp(a,b,c,d,J,P,T,O,Orb_NN_t)` (T=-1/0/+1); `O2b_Pandya_transformation` output type `o2b_pandya_t` with `pnpn`; optional `Report` argument of `O2b_transformation`.
19. `nj_obdm`: rename `OBDM_*` → `rho_*` (9 routines), fix the A^(1/3) normalization.
20. `nj_functions`: fix `J1b_convective`; add `phonon_weights`.
21. `nj_utils`: add `parity_of(L)`, `triangle(j1,j2,j3)`, `hat(j)`.
22. `nj_v2b`/`nj_v3b` (owned by V2b/V3b group): readers build their own `Int.Nmax` orbitals (`orbitals_make(Params,Use_Int=.true.)`); the solver `Orb` is `Calc.Nmax`-sized.

---------------------------------------------------------------------------------------------------------------------

## 11. Julia bugs / races noticed (this group) and status of port notes

* `orbitals_make` (dev orb.jl:36/58): infinite `while` loop when A−Z or Z exceeds the basis capacity. Port: `fatal()`.
* `orbitals_make`: orbitals after the last (partially) filled shell get `min(1,(A−Z−N_occ)/j_hat)`, a ±1e-16 residual
  instead of exactly 0. Harmless; reproduce for bit parity (they enter `Rho = Σ O_a φφᵀ`).
* `OBDM_natural_orbitals`: `l_a, j_a` of the **neutron** natural orbitals are taken from the proton eigenvectors
  `pC`. Harmless while both blocks are (l,j)-block diagonal and Hungarian keeps column a in the block of orbital a;
  keep as in dev (or derive from nC — owner decision).
* `OBDM_natural_orbitals`: returned occupations include the `+1e-10*a` splitting.
* `orbitals_1p1h_make`: open-shell orbitals with 0.1 ≤ occupation ≤ 0.7 are silently dropped from TDA/RPA.
* `O2b_export/O2b_import/O2b_count`: still `J_max = 2*N_max+1` (out-of-bounds when N2max < 2 Nmax) → port note 3 stays.
* `V2b_read` CMS2B `T[a,c]` indexing of an `O1B` is still present (`V2b.jl:178,201`) → port note 2 stays (V2b group).
* `orbitals_export(Params,Orb)` HRF still reads non-existent `Bin/pE_HF.bin` but is dead code → port note 1 obsolete.
* Port note 4 (summary tables) and 5 (HF_allocate parallel loop) belong to the HF group.
* `J1b_convective` sqrt(l−1) bug: **fixed upstream** (§4.2); the port must take the fix.
* `rho_radial_summary` normalization: fixed upstream (§8).
* The Fortran-side `hungarian` already matches the new in-house Julia `hungarian` (square matrices).

---------------------------------------------------------------------------------------------------------------------

## 12. Effort (Fortran lines, rough)

| module | lines |
|---|---|
| nj_types | +80 / −20 |
| nj_params | +90 / −40 |
| nj_input (restructure + new groups) | +250 / −150 |
| nj_orbitals | +330 / −170 |
| nj_o1b | +50 |
| nj_o2b | +60 |
| nj_obdm (rename + fix) | ±40 |
| nj_functions, nj_utils | +40 |
| **total** | **≈ 950 touched (≈ 600 net new)** |
