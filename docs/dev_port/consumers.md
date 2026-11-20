# Dev port spec: consumers (EMPM, NCSM, app/dispatch layer)

Group: EMPM/NCSM/app consumers of the changed outputs.
Base Julia: `NuJuliet` @ 2595e82 (`Source/`), dev Julia: `NuJuliet_dev` @ f3f7a33 (`src/`).
Fortran files in scope: `src/empm/*` (11 modules + `nj_dispatch_empm.f90` + `stub/`), `src/ncsm/nj_ncsm.f90`,
`src/app/nujuliet.f90`, `src/app/nj_input.f90`, `src/{rpa,qp,empm}/nj_dispatch_*.f90` and `src/*/stub/nj_dispatch_*.f90`.

EMPM and NCSM have no Julia counterpart (they are Fortran-only extensions). Dev changes do not touch their physics. What
changes is everything they **consume**: the HF output files (names, directory, the source of the occupations), the
`Orb1B` type, the chain-A helper routines they borrow, and the parameter and run layer. The sections below list each
dependency, what dev does now, and the change to make in the port.

---------------------------------------------------------------------------------------------------------------------

## 1. Dev facts the consumers depend on (with file:line)

### 1.1 Output root and HF layout

* The root directory is lowercase **`io/`** everywhere in dev (`src/HF/HF.jl:36-42`). The base used `IO/`.
* `HF(Params)` (`src/HF/HF.jl:36-42`) removes `io/<Path>` recursively, then creates only
  `io/<Path>`, `io/<Path>/HF`, `io/<Path>/HF/bin`, `io/<Path>/HF/densities`.
  **There is no `io/<Path>/Bin` any more**, and no `HF-EFA` directory.
* Files written by HF (`src/HF/HF_export.jl:147-168`, `src/HF/HF_solver.jl:27-52`):

| base file (`IO/<Path>/...`) | dev file (`io/<Path>/...`) | content / layout |
|---|---|---|
| — (none) | `HF/bin/orb.bin` | **new**: `a_max` records of `Int64 a, n, l, j; Float64 pO, nO` (48 bytes each), `orbitals_export` `src/common/orb.jl:111-131` |
| `Bin/C_HF.bin` | `HF/bin/C.bin` | unchanged `O1b_export` layout: protons then neutrons, `for a, for b if l_a==l_b && j_a==j_b: Float64 A[a,b]` (`src/common/O1b.jl:41-70`) |
| `Bin/Rho_HF.bin` | `HF/bin/rho.bin` | unchanged layout; still the OBDM **in the HF basis**, `C' Rho C` (`HF_solver.jl:35`) |
| `Bin/h_HF.bin` | `HF/bin/H1b.bin` | unchanged layout; still `diag(SPE)` (`HF_solver.jl:135` `return O1B(diagm(SPE.p),diagm(SPE.n)),...`) |
| `Bin/V2B_HF.bin` (only if `BMF`) | `HF/bin/H2b.bin` (only if `HF.BMF == true`) | unchanged `O2b_export` layout (`src/common/O2b.jl:195-257`, the pn block is still the `N_T0^2` square loop with the packed index) |
| `Bin/V2B_HF.dat`, `Orbitals_LHO.dat`, `Orbitals_HF.dat` (`HRF`) | **removed** | `HRF` is no longer read anywhere in dev (`grep '\.HRF'` finds nothing); `O2b_export_HRF` and `orbitals_export(Params,Orb)` remain as dead code |
| `HF/HF_Summary.dat` | `HF/HF_summary.dat` | same `E_HF = %15.8f` line (`HF_export.jl:82`), so `read_ehf` still parses it |
| `HF/Densities/HF_Radial_Densities.dat` | `HF/densities/HF_radial_densities.dat` | |
| — | `HF-MBPT/bin/{orb,C,H1b,rho,H2b[,sigma]}.bin` | when `HF.MBPT.Export` is set: a second **reference** in natural orbitals (`src/HF/HF_MBPT.jl:2233-2313`), where **H1b is not diagonal** |

* **The occupations change meaning** (`src/HF/HF_orbitals.jl:1-140`). `HF_orbitals_make` re-derives the occupations
  in every iteration by Aufbau on the HF single-particle energies, and `Orb .= Orb_new` (line 137) writes them back into
  `Orb`. With `HF.ESH = false` (the default) the last level may be **fractionally** occupied:
  `Occ = min(1.0, (Z - Z_occ) / j_hat)` (line 85). `ESH = true` forces closed sub-shells or throws.
  The base kept the naive LHO filling of `orbitals_make` (integers) and only reordered the columns of C.
  `orb.bin` is therefore the **only** correct source of the occupations of a dev HF reference. Re-running
  `orbitals_make` gives the LHO guess, which can differ from it (level inversions, open shells).
* `Orb1B.pO, nO` are `Float64` in dev (`src/import_structures.jl` `struct Orb1B ... pO::Float64 nO::Float64`). In
  the base they were `Int64`.
* `orbitals_make(Params; Int=false)` (`src/common/orb.jl:1-89`) now uses **`Params.Calc.Nmax`** by default. The base used
  `Params.Int.Nmax`. Only `V2b_read`/`V2b_parse` pass `Int = true`. It also fills fractionally: `Occ = min(1.0,(A-Z-N_occ)/j_hat)`.
  The fill order is still N ascending, n ascending (l descending), j = 2l+1 first.
* The column ordering of C (HF orbital `a` has the `(l,j)` of LHO orbital `a`) is the same in base and dev.

### 1.2 TDA/RPA reference interface (the dev pattern that EMPM/NCSM should copy)

`TDA_solver_import` (`src/TDA/TDA_solver.jl:33-67`) and `RPA_solver` (`src/RPA/RPA_solver.jl:1-30`) read a reference
directory `io/<Path>/<Reference>/bin/` (`Reference = "HF"` by default; `"HF-MBPT"` is also valid):

```julia
IO = "io/" * Params.Calc.Path * "/" * Params.Calc.TDA.Reference
Orb = orbitals_import(Params,IO * "/bin/orb.bin")                       # TDA_solver.jl:38
C   = O1b_import(Params,Orb,IO * "/bin/C.bin");  Rho = O1b_import(...,"/bin/rho.bin")
Orb_1p1h = orbitals_1p1h_make(Params,Orb;Tol_par=...,Tol_hol=...)       # p/h from occupations
h_N = O1b_import(Params,Orb,IO * "/bin/H1b.bin")
V_NN, Orb_NN = O2b_import(Params,Orb,IO * "/bin/H2b.bin",Make_Orb_NN = true)
```

`TDA(Params)` refuses to run if `io/<Path>/<Reference>` is missing (`TDA.jl:26-30`), and refuses `Name == Reference` or
a non-TDA directory that already exists (`TDA.jl:43-55`, marker file `.TDA_output`, although nothing in dev writes that
marker). Output goes to `io/<Path>/<Name>/{bin,densities,transitions/{E0,E1,E2,E3},spectra,amplitudes}` (`TDA.jl:57-67`,
`RPA.jl:31-60`).

| TDA output (dev) | RPA output (dev) |
|---|---|
| `<Name>/TDA_summary.dat` | `<Name>/RPA_summary.dat` |
| `<Name>/bin/TDA_X.bin`, `TDA_E.bin` (`for J in 0:J_max, P in 1:2`: `write(permutedims(X[P,J+1]))`, `E[P,J+1]`; `TDA_export.jl:140-171`) | `<Name>/bin/RPA_{X,Y,E}.bin`, `orb.bin`, `C.bin`, `rho.bin` (+`H1b.bin` with `RPA_Orb`), `sigma.bin` (if `TBCF`) |
| `<Name>/spectra/TDA.dat`, `TDA_plot.dat` | `<Name>/spectra/RPA.dat`, `RPA_plot.dat` |
| `<Name>/amplitudes/TDA_amplitudes{,_ordered}.dat` | `<Name>/amplitudes/RPA_amplitudes{,_ordered}.dat` |
| `<Name>/transitions/E0/TDA_E0.dat`, `E1/TDA_E1{,_V,_V_conv,_V_spin,_T,_T_conv,_T_spin,_C,_C_conv,_C_spin,_c,_NLO_LWA}.dat`, `E2/TDA_E2.dat`, `E3/TDA_E3.dat` | same with `RPA_` |
| externals: `<Name>/densities/TDA_<File_Name>_JP<JP>_nu<nu>_E<E>.dat`, `..._Averaged.dat`, `<Name>/levels/TDA_<File_Name>_...` | same with `RPA_` |

The 1p-1h basis of dev TDA (`orbitals_1p1h_make`, `src/common/orb.jl:348-445`) has `a` (particle) as the outer loop
and `b` (hole) as the inner loop over **all** orbitals, with the proton and neutron entries interleaved per `(a,b)`. It
selects `pO_b > Tol_hol && pO_a < Tol_par && pO_b - pO_a > 1e-8`. The TDA defaults are `Tol_hol = 0.7`, `Tol_par = 0.1`.

### 1.3 Spurious 1- treatment in dev TDA (`src/TDA/TDA_spurious.jl`, `src/TDA/TDA_diagonalize.jl:50-93`)

* Spurious vector: `Spur[i] = (O_h - O_p) * TrOp.E1.{p,n}[a_p,a_h]`, noise cut **`abs < 1e-14`** (base: `1e-6`), then
  normalized (`TDA_spurious.jl:11-36`).
* The orthonormal complement is a closed-form Gram–Schmidt (`TDA_spurious_ortho!`, `TDA_spurious.jl:41-107`). It
  replaces the base Householder + MGS `Spur_Ortho`. The dropped column is the first `a`, counted from `N_ph` down,
  with `v2_a <= 0 || sqrt(u2_a/v2_a) < 1e-4`.
* Physical block: `eigen!(Symmetric(A_ort),Symmetric(N_ort))`. Spurious level:
  **`E[1] = dot(Spur,A,Spur)/dot(Spur,N.*Spur)`**, `X[:,1] = Spur/sqrt(N_spur)` (`TDA_diagonalize.jl:85-93`).
  The base and the current Fortran EMPM set `E(1) = 0`.
* Other channels: `H = N^-1/2 A N^-1/2`, `eigen!(Symmetric(H))`, `X = N^-1/2 R` (lines 35-48). With integer occupations
  `N = 1`, so this is the plain symmetric problem.

### 1.4 Public solver API of dev (`scripts/sample_script_run.jl`, `src/*/<X>.jl`)

`HF`, `TDA`, `RPA`, `HF_ScRPA`, `BCS`, `HFB`, `QTDA`, `pnQTDA`, `QRPA`, and for `X ∈ {TDA,RPA}` (keywords):
`X_transition_densities(Params; JP, nu_list, File_Name="Phonon_Densities")`,
`X_transition_operator_densities(Params; JP, nu_list, File_Name="Phonon_Transition_Operator_Densities")`,
`X_transition_analysis(Params; JP, nu_list, File_Name="Phonon_Transtion_Analysis" [sic], Cut=0.1)`,
`X_transition_currents(Params; JP, nu_list, File_Name="Phonon_Currents", M=-1 (=> 0), xN_grid=40, zN_grid=40)`,
`X_transitions_recalculate(Params)`. There are also the utilities
`V2b_parse(Params,N_max_Out,N_2max_Out,Out_File)` and `V3b_no2b_parse(Params,N_max_Out,N_2max_Out,N_3max_Out,Out_File)`.
**Removed**: `HF_RPA*`, `HF_RRPA*` (all of them), HF-EFA.

### 1.5 Parameter structs (`src/import_structures.jl`) — the source of truth for `nj_input`

```
Interaction_Parameters  NN_File="" NNN_File="" hw=0 Nmax=0 N2max=0 N3max=0 cP2N=1 cP3N=1 cRes=1 pG2N=0 nG2N=0   (unchanged)
HF_MBPT_Parameters      Export=false TBCF=false ScH=false                                                          (new)
HF_Parameters           Tol=1e-7 Imax=100 V_HF=false BMF=true HRF=false ESH=false MBPT=HF_MBPT_Parameters()     (EFA removed; V_HF, ESH, MBPT new; HRF dead)
TDA_Parameters          Name="TDA" Reference="HF" Reference_TBCF=false Ortho=true Tol_hol=0.7 Tol_par=0.1          (new)
RPA_Parameters          Name="RPA" Reference="HF" Reference_TBCF=false RPA_Orb=false Ortho=true Transitions="Linear"
                        Tol=1e-5 Tol_hol=0.7 Tol_par=0.1 Imax=150 OBDM_diagonal=false OBDM_alpha=0.6 ScH=true
                        TBCF=false TBCF_alpha=0.8 TBCF_beta=5.0 DEV=false                                          (was Ortho only)
RRPA_Parameters         removed
BCS_Parameters          ScBCS=false Tol=1e-7 Imax=500 q=0.05 A0=0 Z0=0 pD0=0.5 nD0=0.5 BMF=true HRF=false         (unchanged; HRF dead)
HFB_Parameters          Tol=1e-7 Imax=150 Pairing="Full" dLmax=0.5 pL0=-3 nL0=-3 pa0=0.25 na0=0.25 Broy_m=8
                        Broy_Amax=0.8 Broy_Tmax=0.9 LNT=false LNRes=true BMF=true HRF=false
                        (Broy_Qmax, Broy_Bmax removed; Broy_Amax default 0.9 -> 0.8)
QTDA_Parameters         Ortho=true Reference="HFB"     (Reference new)
QRPA_Parameters         Ortho=true Reference="HFB"     (Reference new)
pnQTDA_Parameters       Reference="HFB"                (new)
Calculation_Parameters  A Z hw Nmax N2max N3max CMS="CMS1+2B" Path=<default> HF TDA RPA BCS HFB QTDA QRPA pnQTDA
```

Note that the sample script still writes `OBDM_digonal` [sic] in a commented line and `OBDM_diagonal` in the live
call. The field is `OBDM_diagonal`.

---------------------------------------------------------------------------------------------------------------------

## 2. Shared-module interface changes needed by this group (several groups touch these)

1. **`nj_params::io_path`**: change the root from `'IO/'` to `'io/'`. Add
   `ref_path(Params,Ref,Suffix) = 'io/'//trim(Path)//'/'//trim(Ref)//Suffix`, used by every reference reader
   (TDA, RPA, EMPM, NCSM). Update the comment `"IO/" * Params.Calc.Path * Suffix`.
2. **`nj_types::orb1b_t`**: `pO, nO` go from `integer` to `real(dp)`. Every `Orb%pO == 1` or `== 0` test in the code base
   must use a tolerance (in this group: `nj_orbitals::orbitals_ph_make`, NCSM none). The constructor calls
   `orb1b_t(a,n,l,j,int,int)` in `nj_orbitals::orbitals_make` must pass reals.
3. **`nj_orbitals`**:
   * `orbitals_make(Params[,Int])`: use `Params%Calc%Nmax` by default and `Params%Int%Nmax` only with `Int = .true.`.
     Use the dev fractional filling (port `src/common/orb.jl:1-89` 1:1).
   * New `orbitals_export(Params,Orb,Path)` / `orbitals_import(Params,Path)`: binary layout
     `integer(i8) a,n,l,j; real(dp) pO,nO` per orbital, `a_max(Calc%Nmax)` records.
     **Use `integer(i8)`**: default Fortran integers are 4 bytes and would break byte compatibility with Julia.
   * The current `orbitals_export(Params,Orb)` (HRF, `Orbitals_LHO.dat`/`Orbitals_HF.dat`, which reads `Bin/h_HF.bin`)
     is dead in dev. Rename it (e.g. `orbitals_export_hrf`) or drop it with HRF (see open questions). PORTING_NOTES
     deviation 1 becomes obsolete if HRF is dropped.
   * `read_h_HF_diag(Params,Orb,pSPE,nSPE)` → `read_h1b_diag(Params,Orb,File,pSPE,nSPE)` with the file passed in
     (`ref_path(Params,Ref,'/bin/H1b.bin')`). It must fail if the off-diagonal elements are non-zero (`> 1e-10`), because
     the EMPM and the NCSM 'EMPM' Hamiltonian assume a diagonal h.
   * `orbitals_ph_make(Params,Orb,...)`: take `Orb` **from `orb.bin`** and the `H1b` file path as an argument. Holes are
     `pO >= 1 - 1e-8`, particles `pO <= 1e-8`. Any other occupation is a fatal error
     ("EMPM/NCSM need integer occupations: closed (sub)shell nucleus or &hf ESH=.true.").
4. **Spurious orthogonalization**: put the dev `TDA_spurious_ortho!` (closed-form GS) in a **core** module
   (e.g. `nj_linalg::spurious_ortho(Spur,U) result(Spur_ind)`), so that the new TDA module (chain A) and EMPM share it
   without EMPM depending on chain A. Then delete `nj_hf_rpa::Spur_Ortho`.
5. **`nj_tr1b::tr1b_t`**: dev adds the `F` and `GT` components to the operator and to `Tr1b_transformation`.
   EMPM and NCSM only use E0..E3 and E1; the E0..E3 matrix elements are unchanged (verified by diffing
   `Tr1b_EX_initialize`, only `fCG → cg`). No consumer change apart from recompiling.
6. **`nj_v2b::V2b_read`**: dev's signature is `V2b_read(Params)`, which builds `orbitals_make(Params,Int=true)` internally.
   If the core group adopts it, change `nj_ncsm::hamiltonian_nn` (line ~1686) to `call V2b_read(Params,V,Orb_NN)` and keep
   its own `Orb = orbitals_make(Params)` (Calc%Nmax, the same size as everything else in NCSM). The dev CMS2B bug
   `T[a,c]` on an `O1B` is **still present** (`src/common/V2b.jl:177,193-195`), so PORTING_NOTES deviation 2 is still needed.
7. **`nj_o2b::orb2b_t`**: dev renamed the dense lookup to `Map[a,b,J+1,T+1]` (Int32). The Fortran type already has the
   dense `Dic(a,b,J,T)`. If the core group renames it to `Map`, update `nj_ncsm::coupled` (uses `Orb_NN%Dic`,
   `%J_max`, `%a_max`) and `nj_empm_fmat` (via `O2b_pp/pn/nn` only, so no change there).
8. **`nj_params` structs** (needed by `nj_input`, see §5): `hf_params_t` (−EFA, +V_HF, +ESH, +MBPT), new
   `hf_mbpt_params_t`, new `tda_params_t`, expanded `rpa_params_t`, removed `rrpa_params_t` and `hf_efa_params_t`,
   `hfb_params_t` (−Broy_Qmax, −Broy_Bmax, Broy_Amax = 0.8), `qtda/qrpa` +Reference = 'HFB', new `pnqtda_params_t`, and
   `calculation_params_t` with RRPA removed and TDA and pnQTDA added. Fortran-only: `empm_params_t%Reference = 'HF'` (new, §3),
   and `calculation_params_t%CM_Lawson` (keep).
9. **`nj_utils`**: nothing new. `make_dir` must create parents (mkdir -p), which it already does for `'/EMPM/Spectra'`.

---------------------------------------------------------------------------------------------------------------------

## 3. EMPM (`src/empm/*`) — dependencies and the adaptation plan

### 3.1 Dependency inventory (current code)

| Fortran routine (file:line) | consumes | dev replacement |
|---|---|---|
| `EMPM_TDA` (`nj_empm_driver.f90:61-76`) | `is_dir(io_path(''))`; `orbitals_make(Params)`; `IO/<Path>/Bin/C_HF.bin`; `Bin/V2B_HF.bin` | check `ref_path(Ref,'/bin')`; `orbitals_import(ref/bin/orb.bin)`; `ref/bin/C.bin`; `ref/bin/H2b.bin` |
| `EMPM_space_make` → `orbitals_ph_make` → `read_h_HF_diag` (`nj_empm_space.f90:117`, `nj_orbitals.f90:157-217`) | `Bin/h_HF.bin` diagonal; integer `Orb%pO/nO` from `orbitals_make` | `ref/bin/H1b.bin` (diagonal check); occupations from `orb.bin` with a tolerance |
| `EMPM_2PH/3PH/4PH/FULL/TRANS/SF` (`nj_empm_driver.f90:132,355,624,715,1536,1846`) | `orbitals_make(Params)` (only for the space/window and `Tr1b`), `Bin/EMPM_F.bin`, `Bin/EMPM_TDA_{E,X}.bin` | `orbitals_import(ref/bin/orb.bin)`; EMPM binaries move (§3.3) |
| `EMPM_cm` (`nj_empm_cm.f90:93,212`) | `Bin/C_HF.bin` | `ref/bin/C.bin` |
| `EMPM_TDA_phonons` (`nj_empm_tda.f90:12,40`) | `nj_hf_rpa::HF_RPA_phonon_count` (chain A, removed) | local copy `EMPM_phonon_count` (27 lines, verbatim from `nj_hf_rpa.f90:36-62`) |
| `EMPM_TDA_diagonalize` (`nj_empm_tda.f90:12,148`) | `nj_hf_rpa::Spur_Ortho` (chain A, removed); `E(1) = 0`; noise cut `1e-6` | core `spurious_ortho` (§2.4); `E(1) = <spur|A|spur>`; cut `1e-14` (see 3.2) |
| `EMPM_TDA_export` (`nj_empm_export.f90:8,46,54,114`) | `nj_hf_rpa_export::{rpa_summary_header, dimensions_block, spectrum_block_end}` (chain A, removed or rewritten) | local copies in `nj_empm_export` (≈ 60 lines) or the new TDA summary writer (open question 4) |
| `EMPM_TDA_strengths` (`nj_empm_tda.f90:208`) | `IO/<Path>/EMPM/Transitions/TDA_E<L>.dat` | path per §3.3; formulas already match dev `TDA_rB` (`B_ph = |M_p|^2`, `B_is = |M_p+M_n|^2/4`, `B_iv(E1) = |e_p M_p - e_n M_n|^2` with `e_p = N/A`, `e_n = Z/A`, otherwise `|M_p-M_n|^2/4`; `src/TDA/TDA_transitions.jl:219-290`) |
| `Tr1b_initialize(Params,Orb,1.0_dp)` (driver:85, cm:211) | `chR2 = 1` | unchanged (E0..E3 and E1 do not depend on chR2) |
| `states_dir` (`nj_empm_driver.f90:2305-2323`) | default `io_path('/Bin')` | default `<EMPM dir>/bin` |
| `files.mk` comment | "uses chain A (nj_hf_rpa, nj_hf_rpa_export)" | after §3.2, EMPM depends on core only, so update the comment and drop the stub remark |

### 3.2 Physics/convention changes to apply (to keep "EMPM_TDA in the full space == TDA")

1. **Reference.** Add the Fortran-only `&empm Reference = 'HF'` (string, default `'HF'`). This mirrors `TDA_Parameters.Reference`.
   All six readers use `ref_path(Params,Reference,'/bin/{orb,C,H1b,H2b}.bin')`. If `H2b.bin` is missing, fail with
   "run HF with BMF = .true.". If H1b is non-diagonal (an `HF-MBPT` reference), fail, because the F-matrix
   formalism and `A = (e_p - e_h) δ - F` assume HF.
2. **Occupations.** Holes and particles come from `orb.bin`. In the full space the hole set then equals dev TDA's set with
   `Tol_hol/Tol_par` for integer occupations. For nuclei where HF Aufbau inverts levels relative to the naive LHO
   filling, **the EMPM space changes** (correctly).
3. **Spurious state** (`EMPM_TDA_diagonalize::TDA_channel` / `spurious_state`):
   * vector `ph_CMS(k) = (O_h - O_p) * E1(a_p,a_h)`, which equals `E1(a_p,a_h)` for integer occupations, cut at
     `abs < 1e-14` (was `1e-6`). This changes the spurious vector by components in `[1e-14,1e-6)` and so shifts
     energies only at roughly the 1e-10 level. Note it when comparing with old runs.
   * Use the core `spurious_ortho` for U. The eigenvalues are unchanged; eigenvectors can differ only within
     degenerate subspaces and by sign.
   * Physical block: since A is symmetric, the symmetric solver (`eig_sym`) on `U_sub^T A U_sub` can replace
     `eig_general`. The sort order (ascending) stays the same.
   * **`E(1) = dot(Spur, A Spur)`** instead of `0`, written to `EMPM_TDA_E.bin`, `TDA.dat`, `TDA_E1.dat` and the log.
     `EMPM_phonons_import` already recomputes `E(i_cm) = R(i_cm) = <cm|A|cm>` from A (`nj_empm_phonons.f90:18-20`),
     so the 2PH/3PH/4PH/FULL/TRANS/SF results do **not** change.
4. Keep the EMPM 1p1h ordering (protons first, particle outer, hole inner over the window levels). It is internal to the
   EMPM binaries and does not need to match `orbitals_1p1h_make`.

### 3.3 EMPM output layout (proposal; dev has no EMPM, so only the conventions apply)

The old `IO/<Path>/Bin` no longer exists after a dev-style HF, so EMPM must create and own its own `bin/`.
Proposal: `io/<Path>/EMPM/` as the EMPM "Name" directory, with lowercase standard subdirectories as in dev:

| current | proposed |
|---|---|
| `IO/<Path>/Bin/EMPM_F.bin` | `io/<Path>/EMPM/bin/EMPM_F.bin` |
| `Bin/EMPM_TDA_E.bin`, `Bin/EMPM_TDA_X.bin` | `EMPM/bin/EMPM_TDA_E.bin`, `EMPM/bin/EMPM_TDA_X.bin` (layout unchanged) |
| `Bin/EMPM_2PH_<JP>_{E,C}.bin`, `Bin/EMPM_FULL_<JP>_{E,Z}.bin` | `EMPM/bin/...` (same names) |
| `Bin/EMPM_FULL<n>_<JP>_{E,Z,Ncm}_<lab>.bin`, `Bin/EMPM_CMB_*`, `Bin/EMPM_CMH_*` | `EMPM/bin/...` |
| `Bin/EMPM_S<n>_<JP>.bin` (default `States_Dir`) | `EMPM/bin/EMPM_S<n>_<JP>.bin` (`States_Dir` override unchanged) |
| `Bin/EMPM_SF*` scratch (driver:2313 `Dir = io_path('/Bin')`) | `EMPM/bin` |
| `EMPM/Spectra/TDA.dat` | `EMPM/spectra/TDA.dat` |
| `EMPM/Transitions/TDA_E<L>.dat`, `EMPM/Transitions/Full<n>_<lab>_E1*.dat` | `EMPM/transitions/...` (optionally `transitions/E<L>/` like dev, see open question 3) |
| `EMPM/TDA_Summary.dat` | `EMPM/TDA_summary.dat` (dev capitalization) |
| `EMPM/Levels.dat`, `EMPM/TDA_Structure.dat`, `EMPM/{2,3,4}Phonon/*`, `EMPM/Full/<n>Phonon/*`, `EMPM/SF/*` | keep the names (Fortran-only; renaming is an owner decision, open question 3) |

Implementation: add `empm_path(Params,Suffix) = io_path(Params,'/EMPM'//Suffix)` in `nj_empm_space` and replace every
`io_path(Params,'/Bin/EMPM_...')` and `'/EMPM/...'` (about 45 call sites, listed by
`grep -n "io_path(Params,'/Bin\|io_path(Params,'/EMPM" src/empm/*.f90`). `EMPM_TDA` does
`remove_dir(EMPM)`, then `make_dir(EMPM/bin)` and `make_dir(EMPM/spectra)`. Because the binaries now live under
`EMPM/`, re-running `EMPM_TDA` also removes the 2PH..FULL binaries and the default state files. That is the right
behavior, because they are invalid after a new TDA. External `States_Dir` files survive as before.

### 3.4 Routine-by-routine edits and sizes

| file | change | ~lines |
|---|---|---|
| `nj_empm_driver.f90` | 7 entry points: `Orb = orbitals_import(...)`, the reference check, the `C.bin`/`H2b.bin` paths, all `/Bin/` paths through `empm_path`, `states_dir` default, error texts ("run HF first" means the reference is missing) | 90 |
| `nj_empm_space.f90` | `EMPM_space_make` passes the H1b path and the reference `Orb` to `orbitals_ph_make`; new `empm_path`; `Levels.dat` path | 20 |
| `nj_empm_tda.f90` | drop `use nj_hf_rpa`; local `EMPM_phonon_count`; spurious changes of §3.2.3; transitions path | 60 |
| `nj_empm_export.f90` | drop `use nj_hf_rpa_export`; local summary/dimension/spectrum helpers (copied, with the header text of dev `TDA_summary`); paths | 70 |
| `nj_empm_phonons.f90` | `EMPM_TDA_{E,X}.bin` paths | 4 |
| `nj_empm_cm.f90` | two `C.bin` paths and the `Orb` source | 6 |
| `files.mk`, `stub/` | comment only; the EMPM stub dispatcher stays | 3 |
| **EMPM total** | | **≈ 250** |

---------------------------------------------------------------------------------------------------------------------

## 4. NCSM (`src/ncsm/nj_ncsm.f90`)

### 4.1 Dependencies

| routine (line) | consumes now | change |
|---|---|---|
| `NCSM` (~286-300) | `make_dir(io_path('/NCSM'))`; `Orb = orbitals_make(Params)` (Int%Nmax size) | `io/<Path>/NCSM/`. `orbitals_make` now has Calc%Nmax size (dev default), consistent with the `a_max(Calc)` loops. For `Hamiltonian='EMPM'` or `Basis='HF'`, use `orbitals_import(ref/bin/orb.bin)` |
| `hamiltonian_nn` (~1686) | `V2b_read(Params,Orb,V,Orb_NN)`, `T1b(Params,Orb)` | follow the `V2b_read` signature of §2.6. Dev `T1b` uses `Params.Int.Nmax` (`src/common/T1b.jl:3`), unchanged |
| `hamiltonian_empm` (1718-1723) | `is_dir('/Bin')`, `Bin/V2B_HF.bin`, `Bin/h_HF.bin`, `Bin/C_HF.bin`, `Bin/Rho_HF.bin`, `HF/HF_Summary.dat` | `ref/bin/H2b.bin`, `ref/bin/H1b.bin` (diagonal check), `ref/bin/C.bin`, `ref/bin/rho.bin`, `io/<Path>/HF/HF_summary.dat` (`read_ehf` parses `E_HF =`, unchanged). Reference = `&empm Reference` (shared with the EMPM; it must be HF-like) |
| occupations `occ(a,t) = Rho%p(a,a)` (1729-1731), `RefO = occ > 0.5` (1745) | rho.bin diagonal | add a check that `occ ∈ {0,1}` within 1e-8, otherwise fail: the plain-form construction and the SD reference need a determinant (an open-subshell dev HF with `ESH = .false.` gives fractions) |
| `window_make` (697-715) → `EMPM_space_make` | as in EMPM | inherits §3 |
| `trans_init` (741) | `Bin/C_HF.bin` | `ref/bin/C.bin` |
| outputs `NCSM/NCSM_<Ham>_...dat`, `_E<L>.dat`, `_SF_*`, `_ckpt.bin`, `_vecs.bin` (574-584, 947, 1155-1157) | `IO/<Path>/NCSM/` | `io/<Path>/NCSM/` through `io_path`; no rename needed |

### 4.2 Size: ≈ 50 lines (paths, two occupation checks, the `orbitals_import` branch, the `V2b_read` call).

---------------------------------------------------------------------------------------------------------------------

## 5. App layer: `nj_input.f90`, `nujuliet.f90`, dispatchers

### 5.1 Namelists (`nj_input::read_input`)

Mirror the dev structs (§1.5). Each namelist group stays optional; defaults come from the `nj_params` types.

| namelist | variables (dev defaults) | status |
|---|---|---|
| `&interaction` | NN_File, NNN_File, hw, Nmax, N2max, N3max, cP2N, cP3N, cRes, pG2N, nG2N | unchanged |
| `&calculation` | A, Z, hw, Nmax, N2max, N3max, CMS='CMS1+2B', Path='' (→ default name), CM_Lawson=0 (Fortran-only) | unchanged |
| `&hf` | Tol=1e-7, Imax=100, V_HF=F, BMF=T, HRF=F, ESH=F | **+V_HF, +ESH** (HRF kept as input but dead in dev, open question 1) |
| `&hf_mbpt` (new) | Export=F, TBCF=F, ScH=F | **new** (`HF.MBPT`); a namelist cannot nest a derived type cleanly, so use a separate group |
| `&hf_efa` | — | **remove** (reading an old file with `&hf_efa` then hits end-of-file, which is silently ignored by `check`; consider a warning) |
| `&tda` (new) | Name='TDA', Reference='HF', Reference_TBCF=F, Ortho=T, Tol_hol=0.7, Tol_par=0.1 | **new** |
| `&rpa` | Name='RPA', Reference='HF', Reference_TBCF=F, RPA_Orb=F, Ortho=T, Transitions='Linear', Tol=1e-5, Tol_hol=0.7, Tol_par=0.1, Imax=150, OBDM_diagonal=F, OBDM_alpha=0.6, ScH=T, TBCF=F, TBCF_alpha=0.8, TBCF_beta=5.0, DEV=F | **expanded** |
| `&rrpa` | — | **remove** |
| `&bcs` | unchanged | — |
| `&hfb` | Tol, Imax, Pairing, dLmax, pL0, nL0, pa0, na0, Broy_m=8, Broy_Amax=**0.8**, Broy_Tmax=0.9, LNT, LNRes, BMF, HRF | **−Broy_Qmax, −Broy_Bmax** (a file that still sets them fails with a namelist error, so document it) |
| `&qtda`, `&qrpa` | Ortho=T, Reference='HFB' | **+Reference** |
| `&pnqtda` (new) | Reference='HFB' | **new** |
| `&empm` | + **Reference='HF'** (Fortran-only) | +1 |
| `&ncsm` | unchanged | — |
| `&transition_densities` | JP='1-', nu_list, File_Name='Phonon_Densities' | unchanged; used by `TDA_/RPA_transition_densities` |
| `&transition_currents` | JP, **M=-1** (dev default; −1 means 0), nu_list, File_Name='Phonon_Currents', xN_grid=40, zN_grid=40 | default of M changes from 0 to −1 (same result) |
| `&transition_operator_densities` (new) | JP, nu_list, File_Name='Phonon_Transition_Operator_Densities' | **new** (if the RPA group ports the routine) |
| `&transition_analysis` (new) | JP, nu_list, File_Name='Phonon_Transtion_Analysis' [sic, dev spelling], Cut=0.1 | **new** (same condition) |
| `&parse` (new) | N_max_Out, N_2max_Out, N_3max_Out, NN_Out_File='io/NN_parsed.bin', NNN_Out_File='io/NNN_parsed.bin' | **new**, for the `V2b_parse`/`V3b_no2b_parse` actions |
| `&run` | Solvers(64) | action list, §5.2 |

Character lengths: `Name` and `Reference` are `character(len=64)`. They are path components, so validate that they
contain no `/` and are non-empty.

Several TDA/RPA runs with different `Name`/`Reference` (e.g. HF-TDA and HF-ScRPA in the dev sample) need one input file
each, because a namelist group is read once per file. `nujuliet a.nml b.nml` already processes files in sequence. Only the
first file should contain `HF`, since HF wipes `io/<Path>`. Document this in the README.

`export_request_t` (in `nj_input`) gains `Cut` and the defaults `M = -1` and `File_Name` per request kind. Requests for
the four external kinds become `Densities`, `OpDensities`, `Analysis`, `Currents`.

### 5.2 `&run` actions and dispatch

| action | handler | old action |
|---|---|---|
| `HF` | `nj_hf_driver::HF` (now also V_HF, MBPT export) | `HF` |
| `TDA` | new chain-A TDA driver | part of `HF_RPA` |
| `RPA` | new chain-A RPA driver | part of `HF_RPA` |
| `HF_ScRPA` | new chain-A driver | — (base had the Temp/HF_ScRPA prototype only) |
| `TDA_transition_densities`, `TDA_transition_operator_densities`, `TDA_transition_analysis`, `TDA_transition_currents`, `TDA_transitions_recalculate` | chain-A externals | — |
| `RPA_transition_densities`, `RPA_transition_operator_densities`, `RPA_transition_analysis`, `RPA_transition_currents`, `RPA_transitions_recalculate` | chain-A externals | `HF_RPA_transition_{densities,currents}` |
| `BCS`, `HFB`, `QTDA`, `QRPA` | chain B | same |
| `pnQTDA` | chain B (new) | — |
| `V2b_parse`, `V3b_no2b_parse` | core utility (new; `&parse`) | — |
| `EMPM_TDA`, `EMPM_2PH`, `EMPM_3PH`, `EMPM_4PH`, `EMPM_FULL`, `EMPM_TRANS`, `EMPM_SF`, `EMPM` | `nj_dispatch_empm` | unchanged |
| `NCSM` | `nj_ncsm::NCSM` | unchanged |
| removed: `HF_RPA`, `HF_RPA_transition_densities`, `HF_RPA_transition_currents`, `HF_RRPA`, `HF_RRPA_transition_densities`, `HF_RRPA_transition_currents` | — | **fatal with a hint**: `'HF_RPA' was removed upstream, use 'TDA', 'RPA'` (and the same for HF_RRPA: "removed upstream, no replacement") |

Changes per file:

* `src/rpa/nj_dispatch_rpa.f90` (rename the module/dir if the RPA group renames chain A): replace the six cases with the
  13 TDA/RPA cases above. Pass `Inp%<Request>` arguments as today (`nu_list(1:N_nu)`, `M`, grids, `Cut`).
  Use `M_proj = merge(0, M, M == -1)` in the callee, as dev does (`RPA_external.jl:775`). About 50 lines.
* `src/rpa/stub/nj_dispatch_rpa.f90`: update the comment only. About 2 lines.
* `src/qp/nj_dispatch_qp.f90`: add `case ('pnQTDA')`. The stub is unchanged. About 4 lines.
* `src/empm/nj_dispatch_empm.f90`: unchanged apart from the comment. The `'EMPM'` combined action stays.
* `src/app/nujuliet.f90`: print `io/` instead of `IO/` (line 27). Add the removed-action hints before
  `Unknown solver`, and the `V2b_parse`/`V3b_no2b_parse` cases (or put them in a small `nj_dispatch_core`). Change the usage
  text "run from the directory containing io/". About 25 lines.
* `src/app/nj_input.f90`: about 250 lines changed or added (new groups, removed groups, the struct constructors
  `rrpa_params_t(...)`/`hfb_params_t(...)` with positional arguments must be rewritten. **Prefer component-wise assignment**
  over positional structure constructors, because these break silently when a field is inserted).

---------------------------------------------------------------------------------------------------------------------

## 6. Julia bugs/races seen in dev (relevant to this group) and the port-note status

| item | dev location | status / action |
|---|---|---|
| `RPA(Params)` and `HF_ScRPA` remove `io/<Path>/<Name>` with **no** `Name != Reference` guard (unlike `TDA`). `Name = "HF"` would delete the HF reference | `src/RPA/RPA.jl:31-35, 102-106` | **new deviation**: in the port, refuse `Name == Reference` as in TDA (`TDA.jl:43-46`). Add it to PORTING_NOTES |
| `RPA(...)`: `if !isdir("io/"*Path) mkdir(...); return; end` creates the dir and silently returns. This is unreachable because the reference check comes first | `RPA.jl:26-29`, `97-100` | harmless; port as a no-op |
| `RPA_transition_analysis` writes `IO_RPA/levels/...` without creating `levels/` (the TDA version does `mkdir(IO_TDA*"/levels")`, `TDA_external.jl:807`) | `src/RPA/RPA_external.jl:~727` | **new deviation**: `make_dir(<Name>/levels)` |
| `TDA(Params)` accepts an existing directory only if it holds `.TDA_output`, but nothing writes that marker, so a re-run with the same Name fails unless the dir is empty | `TDA.jl:48-55` | owner decision (open question 6). Faithful: same check plus writing `.TDA_output` after `mkdir` |
| `V2b_read` `CMS == "CMS2B"` indexes `T[a,c]` on an `O1B` | `src/common/V2b.jl:177,193-195` | **still a bug**, so PORTING_NOTES deviation 2 is still needed |
| `orbitals_export(Params,Orb)` (HRF) reads `IO/<Path>/Bin/pE_HF.bin`; it is now dead code | `src/common/orb.jl:196-215` | deviation 1 is obsolete if HRF is dropped; if HRF is kept as a Fortran extension, read `HF/bin/H1b.bin` |
| `O2b_export/import` pn block: the square loop with the packed index (non-injective, last write wins) | `src/common/O2b.jl:214-226`, `367-384` | unchanged from base. Per CLAUDE.md the port writes the true symmetric pn matrix, so Julia (base and dev) **cannot read a pn block of `H2b.bin`/`sigma.bin` written by the port**, while the port reads Julia files. No change is needed for EMPM/NCSM, which only read. A mixed Fortran-HF to Julia-TDA chain would need the port to write the Julia aliasing (core group decision) |
| `O2b_export` writes channels `J ∈ 0:2Nmax+1` while `O2b_initialize` allocates `0:N2max+1` | `O2b.jl:198`, `O2b.jl:5` | unchanged from base, so PORTING_NOTES deviation 3 is still needed |
| HF `HRF`, BCS `HRF`, HFB `HRF` fields are read nowhere | `import_structures.jl` | dead parameters; keep them in the namelist for compatibility, or drop them (open question 1) |
| HF port notes 6, 7, 9 (HF_RPA_OBDM race, RRPA `OBDM == "Full"`, the `norm(Y)^2 < 0.5` cut) | base HF_RPA/HF_RRPA | for RRPA: **obsolete**, the code was removed upstream. For HF_RPA_OBDM: the RPA group must check `RPA_OBDM_oY2` (`src/RPA/RPA_rho.jl:1-90`) |

---------------------------------------------------------------------------------------------------------------------

## 7. README.md sections to rewrite (Fortran repo)

1. **Title table (lines 7-21)**: replace the HF_RPA/HF_RRPA rows with the action list of §5.2, and the "Julia call" column
   with the dev names (`TDA(Params)`, `RPA(Params)`, `HF_ScRPA(Params)`, `pnQTDA(Params)`, `TDA_/RPA_transition_*`,
   `V2b_parse`, `V3b_no2b_parse`). Drop "HF-EFA" from the HF row and add "HF-MBPT export, ESH".
2. **Running (40-75)**: `IO/` → `io/`; the namelist block per §5.1 (new `&hf_mbpt`, `&tda`, `&rpa`, `&pnqtda`, `&parse`,
   `&transition_operator_densities`, `&transition_analysis`; removed `&hf_efa`, `&rrpa`; HFB Broyden changes); the sample
   `&run Solvers='HF','TDA','RPA'`. Explain one input file per TDA/RPA flavour (Name/Reference), and that `Scripts/` is now
   `scripts/sample_script_run.jl` with `io/NN.bin`, `io/NNN.bin`. Mention the `inputs/` directory (referenced but not
   present in the repo).
3. **Output layout (new section)**: `io/<Path>/HF/{bin/{orb,C,rho,H1b,H2b}.bin, densities, HF_summary.dat}`,
   `HF-MBPT/bin`, `<TDA Name>/...`, `<RPA Name>/...`, `EMPM/...`, `NCSM/...`, and the `orb.bin` record layout.
4. **EMPM extension (76-98)**: inputs `HF/bin/{orb,C,H1b,H2b}.bin` of `&empm Reference`; the integer-occupation
   requirement (`ESH`); the output list per §3.3; "equals HF_RPA (TDA)" becomes "equals `TDA` with
   `Reference='HF'`", and the spurious E(1) convention.
5. **EMPM 2PH/3PH/4PH/FULL/TRANS/SF (99-484)**: every `Bin/EMPM_*` → `EMPM/bin/EMPM_*` (≈ 20 mentions, lines 101-102, 128,
   145, 206-207, 225-230, 273, 300, 384-386, 409, 446, 458), `EMPM/Spectra`/`Transitions` → lowercase.
   The test-tool references (`tests/empm_fci_check.py ... V2B_HF.dat`, `Orbitals_HF.dat`, `Bin/h_HF.bin`, lines 414-417) depend
   on the HRF output that dev removed (open question 1).
6. **NCSM (485-554)**: `IO/<Path>` → `io/<Path>`, "HF run in `io/<Path>/HF/bin`" (the reference files), the
   `V2B_HF`/`h_HF`/`C_HF` names → `H2b`/`H1b`/`C`, and the integer-occupation check. The line 550-551 comparison "vs HF_RPA
   TDA" becomes "vs TDA".
7. **Compatibility with the Julia code (555-562)**: "byte-compatible with Julia dev f3f7a33"; `io/*/HF/bin/*.bin`,
   `<Name>/bin/*.bin`; the `orb.bin` int64 layout.
8. **Verification (564-567)**: `tests/compare_outputs.py` must learn the lowercase tree (`bin/*.bin`).
9. `CLAUDE.md` (Runs → solver mapping): "HF_RPA gives both HF-TDA & HF-RPA; HF_RRPA is a separate solver" → "TDA, RPA,
   HF_ScRPA are separate actions"; the "Original Julia code" paragraph should point to `NuJuliet_dev` and
   `include("src/import_NuJuliet.jl")`, run from the directory holding `io/`.

---------------------------------------------------------------------------------------------------------------------

## 8. Open questions for the code owner

1. **HRF output** (`V2B_HF.dat`, `Orbitals_*.dat`): dev removed it, but the EMPM cross-check scripts (`tests/empm_fci_check.py`,
   `empm_cm_exact_check.py`) read it. Drop it (full parity), or keep it as a documented Fortran extension writing
   `HF/H2b.dat` and `HF/orbitals.dat` under `&hf HRF=.true.`?
2. **EMPM Reference**: allow only `'HF'` (recommended, with diagonal-H1b and integer-occupation checks), or also
   `HF-MBPT` natural orbitals (not diagonal, so the F-matrix/TDA formalism would need the full h)?
3. **EMPM directory naming**: (a) move only the binaries into `EMPM/bin` and lowercase `spectra`/`transitions`
   (recommended), (b) also lowercase `2Phonon`, `Full`, `SF`, `Levels.dat` ..., or (c) keep everything and only
   create `io/<Path>/Bin` from EMPM_TDA. Old runs and test scripts depend on the current names.
4. **EMPM TDA summary**: keep the HF_RPA-style header (moved into EMPM) or switch to the dev `TDA_summary` text
   (stability line, lowest excitation)?
5. **Spurious conventions in EMPM_TDA**: adopt dev `E(1) = <spur|A|spur>` and the `1e-14` noise cut (recommended; downstream
   unchanged, TDA outputs differ in E(1) and possibly at about 1e-10), or keep `E(1) = 0` and `1e-6` for continuity with old runs?
6. **TDA directory guard**: port the `.TDA_output` marker check faithfully (it blocks re-runs into a non-empty dir), or
   write the marker so that re-runs work?
7. **Removed actions**: fail hard on `HF_RPA`/`HF_RRPA` (recommended), or keep an alias `HF_RPA → TDA + RPA`?
8. **Open-shell HF** (`ESH = .false.`, fractional occupations): should EMPM/NCSM refuse (recommended), or round the
   occupations and warn?

---------------------------------------------------------------------------------------------------------------------

## 9. Effort summary

| area | Fortran lines (changed/added) |
|---|---|
| EMPM (`src/empm/*`) | ≈ 250 |
| NCSM | ≈ 50 |
| `nj_input.f90` | ≈ 250 |
| `nujuliet.f90` + 3 dispatchers + stubs | ≈ 85 |
| shared helpers owned by this group if nobody else does them (`ref_path`, `orbitals_import/export`, `read_h1b_diag`, the tolerant `orbitals_ph_make`) | ≈ 80 |
| **total** | **≈ 715** (plus about 150 lines of README/CLAUDE.md prose) |

Order: core (§2: `io_path`, `orb1b_t` real occupations, `orbitals_import/export`, `spurious_ortho`) → HF dev output →
`nj_input`/dispatch → EMPM → NCSM → README. EMPM/NCSM can be validated against the old binaries by symlinking
`HF/bin/{C,H1b,H2b,rho}.bin` to the old `Bin/{C_HF,h_HF,V2B_HF,Rho_HF}.bin` plus an `orb.bin` written from the naive filling.
For closed-shell nuclei with a physical HF level order, results must be bit-identical apart from §3.2.3.
