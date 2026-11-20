# Dev-port spec: quasiparticle chain (BCS, HFB, qpO2b, QTDA, QRPA) + pnQTDA

Base: NuJuliet `2595e82` (`Source/`), target: NuJuliet dev `f3f7a33` (`src/`).
Fortran: `src/qp/*` (nj_qpo2b, nj_qp_common, nj_bcs, nj_hfb, nj_qex, nj_qtda, nj_qrpa, nj_dispatch_qp).
Line numbers `dev:file:L` refer to `/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev/src/...`.

Summary of what actually changed physically/numerically (everything else is rename/paths):

| # | Change | Effect on results |
|---|---|---|
| A | `qpO2b_40` pp & nn: sign of the two direct (V) terms flipped | H^40 pp/nn change -> HFB/BCS-BMBPT(2) energy, QRPA B matrix |
| B | `qpO2b_22` pn(11,11) component now computed (needs new Pandya `pnpn`) | qpH2b.bin block H=10 no longer zero; needed by pnQTDA |
| C | HFB: real modified-Broyden mixing of (Rho,Kappa) | HFB iteration path / iteration count; converged result same within Tol |
| D | HFB: safeguarded secant/bisection for Lambda + `Lambda_H` bookkeeping + closed-shell gap Lambda | Lambda iteration; fixes base wrong H shift at L=1 of iterations > 0 |
| E | HFB (l,j)-ordering by block weights + Hungarian + Löwdin projection | never fails; mixed degenerate vectors are projected back to blocks |
| F | BCS & HFB(Pairing="BCS") use `HF_orbitals_make` (HF group) instead of `HF_orbital_ordering` | BCS reference HF now uses *target* (A,Z) fractional filling; see §1.3 |
| G | QTDA/QRPA/pnQTDA `Reference` option (`"HFB"`/`"BCS"`) selects the input folder | paths |
| H | New solver pnQTDA | new |
| I | Output layout `io/<Path>/<Solver>/{bin,densities,transitions,spectra,amplitudes}`, lowercase names | paths |
| J | qpO2b_import: whole-chunk reads + file-size check | error check only |

---------------------------------------------------------------------------------------------------

## 0. Parameters & input (shared: nj_params, nj_input)

Dev structs (`dev:import_structures.jl:61-125`):

```julia
BCS_Parameters  : unchanged (ScBCS=false, Tol=1e-7, Imax=500, q=0.05, A0=0, Z0=0, pD0=0.5, nD0=0.5, BMF=true, HRF=false)
HFB_Parameters  : Tol=1e-7, Imax=150, Pairing="Full", dLmax=0.5, pL0=-3.0, nL0=-3.0, pa0=0.25, na0=0.25,
                  Broy_m=8, Broy_Amax=0.8, Broy_Tmax=0.9, LNT=false, LNRes=true, BMF=true, HRF=false
                  # REMOVED: Broy_Qmax (=3), Broy_Bmax (=0.3);  Broy_Amax default 0.9 -> 0.8
QTDA_Parameters : Ortho=true, Reference="HFB"      # Reference NEW
QRPA_Parameters : Ortho=true, Reference="HFB"      # Reference NEW
pnQTDA_Parameters: Reference="HFB"                 # NEW struct, Calculation_Parameters.pnQTDA
```

Fortran:
* `nj_params.f90` `hfb_params_t`: delete `Broy_Qmax`, `Broy_Bmax`; `Broy_Amax = 0.8_dp`.
* `qtda_params_t`, `qrpa_params_t`: add `character(len=16) :: Reference = 'HFB'`.
* new `type pnqtda_params_t; character(len=16) :: Reference = 'HFB'; end type`, field `Calc%pnQTDA`.
* `nj_input.f90`: `namelist /hfb/` drop `Broy_Qmax, Broy_Bmax`; constructor `hfb_params_t(...)` call loses 2 args.
  `namelist /qtda/ Ortho, Reference`, `namelist /qrpa/ Ortho, Reference`, new `namelist /pnqtda/ Reference` (optional
  group, same `rewind/read/check` pattern). Validate Reference in {`HFB`,`BCS`} (Julia does not validate; it only fails
  on the missing files, see §5).
* `nj_dispatch_qp.f90`: add `case ('pnQTDA') ; call pnQTDA(Inp%Params)`.
* `io_path` (nj_params.f90:292) `'IO/'` -> `'io/'` (shared with all groups).
* HF group owns `HF_Parameters.ESH` (default false); the qp chain depends on it through `HF_orbitals_make` (§1.3).

---------------------------------------------------------------------------------------------------

## 1. BCS (`dev:BCS/*`)

| dev file | base | Fortran now |
|---|---|---|
| BCS.jl | BCS.jl | nj_bcs: `BCS`, `BCS_reference` |
| BCS_solver.jl | BCS_Solver.jl | `BCS_solver`, `HF_BCS_solve` |
| BCS_V2b.jl | BCS_Residual_Interaction.jl (`BCS_V2b_res`) | `BCS_V2b_res` |
| BCS_chemical_potential.jl | BCS_Chemical_Potential.jl (`BCS_Lambda`) | `BCS_Lambda`, `HF_BCS_initialize_chemical_potential` |
| BCS_allocate / density_operator / energy / particle_number / BMBPT | identical (renames only) | unchanged |
| BCS_export.jl | BCS_Export.jl | `BCS_summary`, `BCS_summary_SQS`, `BCS_export` |

### 1.1 Renames (no semantic change)
`BCS_V2b_res` -> `BCS_V2b`; `BCS_Lambda` -> `BCS_chemical_potential`; `OBDM_export` -> `rho_export`;
`qpO1B_export` -> `qpO1b_export`; `V2b_read(Params,Orb)` -> `V2b_read(Params)`, `V3b_no2b_read(Params,Orb)` ->
`V3b_no2b_read(Params)` (core group). Fortran may keep its routine names; renaming is optional.

### 1.2 Output layout (BCS.jl:36-42, BCS_export.jl, BCS_solver.jl:33-34)
```
rm -r io/<Path>; mkdir io/<Path>, io/<Path>/BCS, io/<Path>/BCS/bin, io/<Path>/BCS/densities
io/<Path>/BCS/BCS_summary.dat                     (was IO/<Path>/BCS/BCS_Summary.dat)
io/<Path>/BCS/densities/BCS_radial_densities.dat  (was .../BCS/Densities/BCS_Radial_Densities.dat)
io/<Path>/BCS/bin/{C.bin,U.bin,V.bin,qpH1b.bin,qpH2b.bin}   (was IO/<Path>/Bin/{C,U,V,qpH1B,qpH2B}.bin)
```
No top-level `Bin/` any more. Binary layouts unchanged (O1b_export / qpO1b_export / qpO2b_export).
`BCS_BMBPT.jl:121` appends to `io/<Path>/BCS/BCS_summary.dat`.

### 1.3 `HF_orbital_ordering` -> `HF_orbitals_make` (BCS_solver.jl:240-241) and the reference HF

```julia
#C, SPE = HF_orbital_ordering(Orb,a_max,O1B(pC,nC),pnVector(pSPE,nSPE))
Orb, SPE, C = HF_orbitals_make(Params,pnVector(pSPE,nSPE),O1B(pC,nC),Orb)      # dev:BCS/BCS_solver.jl:241
```
`HF_orbitals_make` (dev:HF/HF_orbitals.jl, HF group) = old ordering (identical 1e-8 <l>,<j> matching) **plus** it
recomputes `Orb[a].pO/nO` (now `Float64`) by filling the energy-ordered levels with the **target** A,Z of `Params`
(fractional last level, or throw if `Params.Calc.HF.ESH` and no closed sub-shell fits) and writes them **in place**
(`Orb .= Orb_new`). Consequences for BCS:

* ScBCS loop: `Orb` (BCS_solver's `Orb`) gets target occupations; no BCS routine reads `pO/nO` except
  `HF_BCS_initialize_chemical_potential` and `BCS_summary_SQS` (which builds its own `orbitals_make(Params_Ref)`),
  so only the possible ESH throw is observable.
* **Reference HF** (`HF_BCS_HF_solve`, BCS_solver.jl:96-133): unchanged call `HF_solve(Params,Orb_HF,...)` with
  `Orb_HF = orbitals_make(Params_Ref)`, but dev `HF_solve` now calls `HF_orbitals_make(Params,...)` every iteration
  (HF group). Hence the "reference" HF is self-consistent with the **target** (A,Z) filling (fractional last level
  unless ESH), not the closed-shell reference occupations of base; `Orb_HF` leaves `HF_solve` with target
  occupations and `HF_BCS_initialize_chemical_potential(Params_Ref,SPE,Orb_HF)` sees them.
* `HF_BCS_initialize_chemical_potential` (BCS_chemical_potential.jl:54-93) unchanged text but `Orb[a].pO == 1` /
  `== 0` are now Float64 tests: a fractionally occupied level is neither hole nor particle. Formula kept:
  `Lambda = 0.5*(Eps_hole - Eps_particle)` (sic, see §8).
* `BCS_summary_SQS` (BCS_export.jl:87,104): `O_a = Orb[a].pO == 1 ? 1.0 : 0.0` with Float64 occupations.

Fortran (nj_bcs):
* `HF_BCS_solve` line ~268: replace `call HF_orbital_ordering(Orb,a_max,C,SPE)` by
  `call HF_orbitals_make(Params,SPE,C,Orb)` (Orb `intent(inout)` through `HF_BCS_solve` and `BCS_solver`).
* `HF_solve(Params,Orb_HF,...)` must be the dev HF_solve (HF group), with `Orb_HF` `intent(inout)` so that the
  mutated occupations reach `HF_BCS_initialize_chemical_potential`.
* `HF_BCS_initialize_chemical_potential` & `BCS_summary_SQS`: compare `Orb(a)%pO == 1.0_dp` / `== 0.0_dp` once
  `orb1b_t%pO/nO` are `real(dp)` (shared change).
* Paths per §1.2; `BCS` creates `/BCS`, `/BCS/bin`, `/BCS/densities` (no `/Bin`).
* Port-note B.2 (`V_NN_Ref` read but unused) still applies (dev still reads `V2b_read(Params)` at line 104 and
  never uses it).
Size: ~40 lines changed.

### 1.4 `BCS_reference` (BCS.jl:73-130)
Copies the new Calc fields `TDA, QTDA, QRPA, pnQTDA` and drops `RRPA`. Fortran copies the whole `Params%Calc`, so
no change beyond the type itself.

---------------------------------------------------------------------------------------------------

## 2. HFB (`dev:HFB/*`)

| dev file | base | Fortran now (nj_hfb) |
|---|---|---|
| HFB.jl | HFB.jl | `HFB` |
| HFB_solver.jl | HFB_Solver.jl | `HFB_solver`, `HFB_solve`, `HFB_diagonalize*`, `bogoliubov_eig`, `bmz_amplitudes` |
| HFB_Broyden.jl | **new** (replaces base HFB_Broyden.jl) | `HFB_Broyden_update(Orb,Rho,Kappa)` (net-effect stub) |
| HFB_Broyden_BackUp.jl | = base HFB_Broyden.jl, **not included** by `import_NuJuliet.jl` | – (ignore) |
| HFB_chemical_potential.jl | HFB_Chemical_Potential.jl (`HFB_Lambda_secant`) | `HFB_Lambda_secant` |
| HFB_ordering.jl | HFB_Ordering.jl | `HFB_orbital_ordering`, `HFB_canonical_orbital_ordering`, `canonical_order` |
| HFB_canonical_basis.jl | identical | `HFB_canonical_basis`, `HFB_canonical_basis_BCS` (call ordering, §2.4) |
| HFB_allocate / density_operator / energy / particle_number / BMBPT / Lipkin_Nogami | renames/paths/whitespace only | unchanged |
| HFB_export.jl | paths | `HFB_solver` export block, `HFB_summary*` |

### 2.1 Output layout
```
rm -r io/<Path>; mkdir io/<Path>, io/<Path>/HFB, io/<Path>/HFB/bin, io/<Path>/HFB/densities
io/<Path>/HFB/HFB_summary.dat
io/<Path>/HFB/densities/HFB_radial_densities.dat
io/<Path>/HFB/bin/{C,U,V,Rho,Kappa}.bin, qpH1b.bin, qpH2b.bin
io/<Path>/HFB/bin/lambda_2.bin          (Lipkin-Nogami, was IO/<Path>/Bin/Lambda_2.bin; export AND import path)
```
Message text changes: `"hw = "` instead of `"HbarOmega = "` in the run header, `"Starting spherical Hartree-Fock-Bogoliubov
calculation ..."`; `HFB_energy` prints its header before allocating. `HFB_solver` sets BLAS threads to 1 around the
solve (performance only; Fortran: nothing, or `mkl_set_num_threads_local(1)` optionally).

### 2.2 Chemical potential (dev:HFB/HFB_chemical_potential.jl) — replaces `HFB_Lambda_secant`

```julia
HFB_chemical_potential_secant(Params,Lambda_1,dA_1,Lambda_2,dA_2,Lambda_lo,Lambda_hi) -> (Lambda,Lambda_lo,Lambda_hi)
   # per species: HFB_chemical_potential_secant_step(L1,dA1,L2,dA2,lo,hi,dLmax=Params.Calc.HFB.dLmax,Tol=Params.Calc.HFB.Tol)
```
`secant_step` (lines 15-62), exact:
```
if |dA_2| < Tol: return (L2, lo, hi)                     # converged, keep
if dA_2 > 0: lo = L2  else: hi = L2                      # dA = A_target - A(Lambda), A increasing in Lambda
Denom = dA_2 - dA_1; dL = 0
if Denom != 0 && L1 != L2 && isfinite(Denom): dL = -dA_2*(L2 - L1)/Denom
if !(dL*dA_2 > 0): dL = sign(dA_2)*dLmax                 # wrong direction / flat / NaN -> max step toward root
if |dL| > dLmax: dL = sign(dL)*dLmax
L = L2 + dL
if !isnan(lo) && !isnan(hi) && (L <= lo || L >= hi): L = 0.5*(lo + hi)    # bisection
return (L, lo, hi)
```
(`sign` = Julia sign, 0 for 0; NaN propagates: `!(NaN>0)` is true -> max step.) Brackets are NaN = unknown.

`HFB_chemical_potential_gap(Params,H,Lambda_H,Delta,Orb)` (lines 64-84): per species, if `maximum(abs.(Delta.x)) < Tol`
(pairing collapsed) return `HFB_chemical_potential_gap_species(H.x,Lambda_H.x,A_target,Orb)`, else NaN.
`..._gap_species` (86-128): for each (l,j) block in order of first appearance (block = all `b` with same l,j, in
orbital order) take `eigvals(Symmetric(H[Block,Block]))` (ascending), push `E + Lambda_H` with degeneracy `j+1`;
`Order = sortperm(SPE)` (stable); accumulate `A_filled += Deg[Order[i]]` for `i = 1:length-1`; if
`|A_filled - A_target| < 1e-8` return `0.5*(SPE[Order[i]] + SPE[Order[i+1]])`; if `A_filled > A_target` break; else NaN.
Note `H` here already contains `-Lambda_H*I`, hence the `+Lambda_H`.

### 2.3 `HFB_solve` changes (dev:HFB/HFB_solver.jl:84-392)

New state: `Lambda_H` (Lambda currently contained in field H), `Lambda_lo, Lambda_hi` (pnFloat NaN), `Lambda_gap`
(pnFloat NaN); Broyden state (`Broyden, BroyVec = HFB_Broyden_initialize(Params,Rho,Kappa,Orb)` after the initial guess).
Loop skeleton (only the differences vs Fortran `HFB_solve`):
```
Lambda_H = Lambda (initial)                                         # L142
outer while ...:
  for L in 1:Imax
    if L == 1: Lambda_lo = Lambda_hi = (NaN,NaN)                     # L155
    if L == 1 && Iteration == 0:   (init block as before, but)
        Lambda,lo,hi = HFB_chemical_potential_secant(...,lo,hi)      # L189
        Lambda_1 = Lambda_2; Lambda_2 = Lambda
        Lambda_H = Lambda_1                                          # L196: H = H_2 was built with old Lambda_2
        SQE,U,V,H,Delta = ..._2
    end
    if L == 1: Lambda_gap = HFB_chemical_potential_gap(Params,H,Lambda_H,Delta,Orb)     # L209-211
    if !isnan(Lambda_gap.p): Lambda.p = Lambda_gap.p                # L212-217
    if !isnan(Lambda_gap.n): Lambda.n = Lambda_gap.n
    Lambda_2 = Lambda                                               # L218 (always)
    H += (Lambda_H - Lambda) * I                                    # L222  (base: (Lambda_1 - Lambda))
    Lambda_H = Lambda                                               # L223
    diagonalize; Rho_Temp,Kappa_Temp; Z,N; Z_1=Z_2; Z_2=Z (same for N)
    Lambda,lo,hi = HFB_chemical_potential_secant(Params,Lambda_1,dA_1,Lambda_2,dA_2,lo,hi)   # L245
    Lambda_1 = Lambda_2; Lambda_2 = Lambda
    dZ,dN; if dZ<eps && dN<eps: break
  end
  Rho,Kappa = HFB_density_operator(U,V)
  Rho,Kappa,BroyVec = HFB_Broyden_update(Params,Iteration,Rho,Kappa,Broyden,BroyVec,Orb)   # L264 (Iteration before ++)
  H,Delta = HFB_allocate(Params,Lambda,Rho,Kappa,...)
  Lambda_H = Lambda                                                 # L268
  LN (Iteration > 15), diagonalize, particle numbers, dE/d2E, degeneracy logic: unchanged
```
The base used `Lambda_1 - Lambda` for the shift; at `L == 1` of iterations > 0 `Lambda_1` is not the Lambda contained
in `H` (H was rebuilt with `Lambda`), so base shifted H by a spurious `Lambda_1 - Lambda`. Dev fixes this (Fortran
reproduces the base bug in `HFB_solve` lines ~236-241: must change to the `Lambda_H` scheme).
Degeneracy branch, `E_HFB_1` with stale `Rho_1,Kappa_1` (port-note B "kept on purpose") unchanged in dev.
The final canonical basis + `HFB_orbital_ordering(...,Final_Ordering=true,C=C)` unchanged (C kwarg unused).

Fortran: rewrite inner loop of `HFB_solve` (~40 lines), replace `HFB_Lambda_secant` by
`HFB_chemical_potential_secant` + `_secant_step` + `_gap` + `_gap_species` (~110 lines; `_gap_species` uses
`eig_sym` values of the extracted block and the stable `sortperm`).

### 2.4 Orbital ordering (dev:HFB/HFB_ordering.jl) — replaces <l>,<j> 1e-3 matching

All three routines now call `HFB_orbital_matching(X,Orb) -> (Order, X_new)`:
* `HFB_orbital_ordering`: optional final sort (unchanged: if `Final_Ordering && Pairing != "BCS"`, stable
  `sortperm(norm(U[:,a])^2)` ascending applied to SQE, U, V), then `Order, UV = matching(vcat(U,V))`,
  `SQE = SQE[Order]; U = UV[1:a_max,Order]; V = UV[a_max+1:2a_max,Order]`.
* `HFB_canonical_orbital_ordering(Params,Occ,C,Orb)`: `Order, C = matching(C)`; `Occ = Occ[Order]; C = C[:,Order]`.
* `HFB_canonical_basis_orbital_ordering(Params,C,Orb)`: same without Occ (used by `HFB_canonical_basis`,
  `HFB_canonical_basis_BCS`).

`HFB_orbital_matching(X, Orb)` (lines 79-185), `a_max = length(Orb)`, `Row_max = size(X,1)` (a_max or 2a_max),
row r belongs to orbital `b = rem(r-1,a_max)+1`:
1. Blocks: `Block[a]` = index of the (l,j) block of orbital a, numbered by first appearance (1..B_max).
2. Weights `W[B,a] = sum_{r: Block[b(r)]==B} X[r,a]^2`, then `W[:,a] ./= sum(W[:,a])`.
3. `Vector_block[a] = argmax(W[:,a])` (first maximum).
4. If the per-block counts of `Vector_block` differ from the block sizes: `Cost[b,a] = 1 - W[Block[b],a]` (a_max x
   a_max), `Assignment,_ = hungarian(Cost)` (row b -> column a), `Vector_block[Assignment[b]] = Block[b]`.
5. `Order[Block_orbitals_B] = Block_vectors_B` for every B (both lists ascending: k-th vector of a block goes to the
   k-th orbital of that block, input order kept).
6. Projection: for each B, if NOT all `W[B,a] > 1 - 1e-10` for its vectors: `P = X[Block_rows,Block_vectors]`,
   `X_new[Block_rows,Block_vectors] = P * inv(sqrt(Symmetric(P'P)))` (Löwdin; Fortran: `eig_sym(P^T P)` ->
   `Z diag(w^-1/2) Z^T`), `X_new[Other_rows,Block_vectors] = 0`. Pure blocks untouched (`X_new = copy(X)`).
   Note: when one vector of a block is mixed, ALL vectors of that block are Löwdin-transformed.
7. Returns `(Order, X_new)`; callers index `X_new[:,Order]`.

`hungarian` is the drop-in implementation in `dev:common/functions.jl:711-858` (shortest augmenting path with dual
potentials, Jonker–Volgenant style, `assignment[row] = col`). Port it **verbatim** (loop order, strict `<`
comparisons, `Inf` init) so ties resolve identically; it is shared with `OBDM_natural_orbitals` (dev:common/orb.jl:593),
so it belongs in `src/core/nj_functions.f90` (shared).

Fortran (nj_hfb):
* replace `order_one` in `HFB_orbital_ordering` and `canonical_order` by `orbital_matching(X,Orb,Order,X_new)`.
  `HFB_canonical_basis` / `_BCS` currently do `C%p = C%p(:,canonical_order(C%p,Orb))` -> must also take the projected
  `X_new`. ~150 lines.
* Port-note B.5 (fatal on failed ordering) becomes obsolete: matching cannot fail.

### 2.5 `HFB_diagonalize_BCS` (HFB_solver.jl:394-460)
`HF_orbital_ordering(Orb,a_max,C,SPE)` -> `Orb, SPE, C = HF_orbitals_make(Params,SPE,C,Orb)` (line 409). The ordering
of C/SPE is the same as before; side effect: the caller's `Orb` gets target-filling occupations (in place), and
with `HF.ESH = true` and an open-shell target HF_orbitals_make throws (§8). No HFB routine reads `pO/nO` after the
initial `HFB_density_operator_initialize` (which does not use them either), so results are unaffected.
Fortran: call the HF group's `HF_orbitals_make` on a local copy `Orb_tmp = Orb` (keeps `Orb` `intent(in)` in
`HFB_solve`), or propagate `intent(inout)`; the ESH fatal must be reproduced either way.

### 2.6 Modified Broyden mixing (dev:HFB/HFB_Broyden.jl) — replaces the base stub

Port-note B.1 ("only the net effect is ported") is **obsolete**: dev mixes for real. Uses existing structs
`HFB_Broyden{Map,Key,M,m}` and `HFB_Broyden_Vector{x,y,z,X,r,s,R}_{Rho,Kappa}` (pnVector / O1B with M x m matrices).
Fortran: module-private derived type in nj_hfb (not shared), e.g.
```fortran
type :: broyden_t
  integer :: M, m
  integer, allocatable :: Key(:,:)        ! (2,M): a,b with a<=b
  real(dp), allocatable :: w(:)           ! weights
  real(dp), allocatable :: x(:,:), y(:,:), z(:,:), r(:,:), s(:,:)   ! (M,4): 1 pRho 2 nRho 3 pKappa 4 nKappa
  real(dp), allocatable :: X(:,:,:), R(:,:,:)                       ! (M,m,4)
end type
```
Exact algorithm:
* **Mapping** (26-59): `M=0; for a in 1:a_max, b in a:a_max: if l_a==l_b && j_a==j_b: M+=1; Map[a,b]=Map[b,a]=M`;
  `Key[i]=(a,b)` (a <= b). (Base stored all ordered pairs; dev only a<=b.)
* **Weights** (61-74): `w[i] = sqrt(j_a + 1)` (j doubled), times `sqrt(2)` if `a != b`.
* **Linearize** (100-116): `x[i] = 0.5*w[i]*(Rho[a,b] + Rho[b,a])` (same for Kappa, p & n). Elements outside the
  (l,j) blocks are dropped.
* **Reconstruct** (118-140): `Rho[a,b] = Rho[b,a] = x[i]/w[i]`, zeros elsewhere. **No 1e-14 noise cleaning any more**
  (base had it).
* **Norm** (142-145): `sqrt(sum over the 4 components of x.^2)` (one joint norm for pRho,nRho,pKappa,nKappa).
* **Initialize** (76-98): `x = linearize(initial Rho,Kappa)`; y,z,r,s zero; X,R zero (M x m each component).
* **history_clear**: zero all X, R.
* **history_update(Iteration)** (161-187): `k = rem(Iteration,m) + 1`; `dR = r - s` (4 comps), `dR_norm = norm(dR)`;
  if `dR_norm < 1e-14` return (no store); `X[:,k] = (x - y)/dR_norm`, `R[:,k] = dR/dR_norm` (all 4 comps).
* **gamma** (189-215): `w0 = 0.01`; `A = sum_comp R_c' R_c` (m x m), `c = sum_comp R_c' r_c`; `A[k,k] += w0^2` for all
  k; `gamma = Symmetric(A) \ c` (Julia: Bunch–Kaufman on upper triangle -> Fortran `dsysv('U',...)`); if any
  non-finite -> `gamma = 0`.
* **update(Params,Iteration,Rho,Kappa,...)** (217-306), `alpha = Broy_Amax`, `T_max = Broy_Tmax`,
  `Iteration_Broyden = LNT ? 17 : 1`:
  ```
  F = linearize(Rho,Kappa)                    # output densities of this iteration
  s = r; r = F - x
  if Iteration < Iteration_Broyden: z = F                                   # plain iteration
  else
      r_norm = norm(r); s_norm = norm(s)
      if Iteration == Iteration_Broyden && Iteration_Broyden > 1: history_clear      # LN switched on
      elseif r_norm > s_norm / T_max: history_clear; history_update(Iteration)
      else: history_update(Iteration)
      gamma = gamma(); z_c = x_c + alpha*r_c - (X_c + alpha*R_c)*gamma   for each of the 4 comps
  end
  if norm(F_pKappa) < 1e-12: z_pKappa = 0 ;  if norm(F_nKappa) < 1e-12: z_nKappa = 0     # collapsed pairing
  y = x; x = z
  Rho,Kappa = reconstruct(z)
  ```
  Note the `history_clear` at `Iteration == 17` with LNT does NOT store a new vector, so the first LN step is
  linear mixing (`gamma = 0` because R = 0).
  The degeneracy pass of `HFB_solve` calls update again with the same (non-incremented) `Iteration`.
Size ~200 Fortran lines (+ dsysv interface if `nj_linalg` lacks one: shared, see §7).

### 2.7 Lipkin–Nogami, energy, particle number, allocate, BMBPT
Only `lambda_2.bin` path (§2.1), print order, whitespace. `HFB_particle_number` initial `Z, N = 0.0, 0.0` (base
`0,0, 0,0` gave the same zeros). No Fortran change except paths. Port-note B items "HFB_diagonalize_MCA uses
unmodified Delta" and "LN leaves Delta unchanged" are still present in dev -> keep.

---------------------------------------------------------------------------------------------------

## 3. Quasiparticle 2-body operators (`dev:common/qpO2b*.jl`) — Fortran nj_qpo2b + core nj_o2b

### 3.1 Structural renames (no effect)
`Orb2B.Dic` (Dict) -> dense `Orb2B.Map[a,b,J+1,T+1]::Int32` + `O2b_lookup` (Fortran already dense `Dic(a,b,J,T)`);
`Orb2B_Temp.Map[a,b,J+1]` (Fortran `Orb_NN_t%Dic(a,b,J)` already); `qpOrb2B{N,i::Ind2qp}` ->
`Orb2B_2qp{N,Ind::Ind2B,Dic}` (Fortran `qporb2b_t%i` may stay; the Dict is unused by this group).
`Sign_ab = 1` -> `1.0` (type stability only). Removed (never ported): `qpO2b_40_allocate_canonical`,
`qpO2b_22_allocate_canonical`.

### 3.2 Pandya transformation type `O2B_Pandya{pp,pnnp,pnpn,nn}` (dev:common/O2b.jl:913-1150) — **shared core**
Old `O2B_Temp{pp,pn,nn}` from `O2b_temp_initialize`; new `O2b_Pandya_initialize` builds the same pair set (all ordered
(a,b) with 2n_a+l_a <= Nmax and 2(n_a+n_b)+l_a+l_b <= N2max, J = |j_a-j_b|/2..(j_a+j_b)/2, channel (P,J), J = 0..N2max+1).
Components, with `Q = parity(l_a+l_d)` required `== parity(l_b+l_c)`, `I = max(|j_a-j_d|,|j_b-j_c|)/2 .. min(j_a+j_d,j_b+j_c)/2`:
```julia
pp  : sum_I (-1)^((j_b+j_c)/2+I) (2I+1) {j_a j_b J; j_c j_d I} O_pp(a,d,b,c;I,Q)     # unchanged
pnnp: sum_I (-1)^((j_b+j_c)/2+I) (2I+1) {...}                O_pn(a,d,b,c;I,Q)     # = old "pn"
pnpn: sum_I (-1)            (2I+1) {...}                     O_pn(a,d,c,b;I,Q)     # NEW  (dev:O2b.jl:1059-1072)
nn  : sum_I (-1)^((j_b+j_c)/2+I) (2I+1) {...}                O_nn(a,d,b,c;I,Q)     # unchanged
```
(`phase(1)` = -1 constant, `hat(2I)` = 2I+1, `w6j` = f6j with doubled args.)
Fortran core `nj_o2b`: add `pnpn` to the Pandya result (either rename `o2b_temp_t%pn` -> `pnnp` and add `pnpn`, or a
new `o2b_pandya_t{pp,pnnp,pnpn,nn}`) and a 4th `pandya(Kind=3)` branch `-(2I+1)*f6j(...)*O2b_pn(a,d,c,b,I,Q,...)`.
Memory of F_NN grows by 1/3.

### 3.3 O^(40) pp & nn: sign of the direct terms (dev:common/qpO2b_40.jl:300, 649)
```julia
Sum += (-ME1 - ME2 + ME3 + ME4 - ME5 - ME6)   # base
Sum += ( ME1 + ME2 + ME3 + ME4 - ME5 - ME6)   # dev  (ME1,ME2 = V_ibcd terms, UUVV and VVUU)
```
Fortran `qpO2b` engine calls for Out=1 (qp40 pp) and Out=3 (qp40 nn): `Sgn = [-1,-1,1,1,-1,-1]` ->
`[1,1,1,1,-1,-1]`. pn O^40 unchanged apart from reading `F_NN.pnnp` (= old pn).

### 3.4 O^(22) pn(11,11) — NEW (dev:common/qpO2b_22.jl:477-633, called at line 14)
Indices: a,c protons; b,d neutrons; bra/ket = T0 (pn) pairs of `Orb_NN` (index 1, all ordered pairs under N2max).
4 terms, all `+`; transformation order 1st idx a (proton U/V), 2nd b (neutron), 3rd c (proton), 4th d (neutron):

| term | source at stage 1 (sum over i in (l_a,j_a) block) | stage 1 (a,p) | 2 (b,n) | 3 (c,p) | 4 (d,n) |
|---|---|---|---|---|---|
| 1 | `O_ibcd = O2b_pn(i,b,c,d,J,P)` | U | U | U | U |
| 2 | `O_ibcd` | V | V | V | V |
| 3 | `F_ibcd = F_NN.pnpn[P,J][Bra_ib,Ket]` | U | V | U | V |
| 4 | `F_ibcd` | V | U | V | U |

Final stage: `for Bra in 1:N0, Ket in 1:Bra` over `Orb_NN.Ind[1,P,J+1]`, `Ind = Bra + (Ket-1)*N0 - Ket*(Ket-1)/2`,
`Q_NN.qp22.pn1111[P,J+1][Ind] = sum_i (UVsel(i,d) * T1_k[Bra_ab,Ket_ci])` summed over k=1..4 with `+`.
Fortran engine: new `Fch = 4` (V source = `O2b_pn(i,b,c,d,J,P,O_NN,Orb_NN)`, F source `Src_Fcd` = `F_NN%pnpn`),
`Iso = [1,2,1,2]`, `Sel = (U,U,U,U),(V,V,V,V),(U,V,U,V),(V,U,V,U)`, `Sgn = [1,1,1,1]`, new `Out = 7`: 4th stage loops
`Orb_NN%N(1,P,J)` / `Orb_NN%Ind(1,P,J)` and writes `Q_NN%qp22%pn1111(P,J)%v(packed_index(Bra,Ket,N0))`.
(The 4th stage currently hard-codes T1 pairs `Orb_NN%N(2,...)`; make the pair set depend on Out.)
New accessor (dev:qpO2b_22.jl:899-902):
```julia
qpO2b_1111_pn(a,b,c,d,J,P,O_NN,Orb_NN) = O_NN.qp22.pn1111[P,J+1][packed(O2b_lookup(a,b,J,0), O2b_lookup(c,d,J,0), N0)]
```
no sign/reorder (a,c proton; b,d neutron). Fortran: `qpO2b_1111_pn` using `Orb_NN%Dic(a,b,J,0)` + `packed_index`
(max/min of Bra,Ket); fatal if a Dic entry is 0. Update the module header comment ("pn1111 stays zero").
Size: ~40 lines nj_qpo2b + ~20 lines nj_o2b.

### 3.5 Export / import (dev:common/qpO2b.jl:176-545)
Export: unchanged layout (11 blocks H=1..11, each `J=0:2Nmax+1`, `P=1:2`; symmetric blocks H=1,3,8,10,11 packed
Bra-major `Ket<=Bra`; full blocks row-major; H=10 = pn1111 now non-zero). Port-note B.3 (`min(2Nmax+1,N2max+1)`
channels) **still needed** (dev export still indexes `Orb_NN.N[...,J+1]` up to 2Nmax+1).
Import rewritten: per (H,J,P) chunk one `seek(8*skip)` + `read!` into a buffer (threads), symmetric chunks unpacked to
the packed lower triangle, full chunks transposed; **new check**: `filesize == 8*N_ME` else `error("Size of the binary
file ''<path>'' (<n> B) does not match the expected size ... (<8*N_ME> B) ...")`.
Fortran `qpO2b_import`: add `inquire(file=...,size=)` check against the expected size computed with the *same* channel
range as `stream_qp` (`J_e = min(2Nmax+1, Orb_NN%J_max)`); fatal with the Julia message. ~15 lines.
`qpO2b_31` unchanged (still unused).

---------------------------------------------------------------------------------------------------

## 4. QTDA (`dev:QTDA/*`) — Fortran nj_qtda (+ nj_qex shared helpers)

| dev file | base | Fortran |
|---|---|---|
| QTDA.jl | QTDA.jl | `QTDA` |
| QTDA_solver.jl | QTDA_Solver.jl | `QTDA_solver` |
| QTDA_allocate.jl | identical up to `Orb_2qp.i` -> `.Ind` | `qex_allocate` |
| QTDA_diagonalize.jl | identical (commented code removed) | `QTDA_diagonalize` |
| QTDA_spurious.jl / QTDA_transitions.jl | identical up to `.Ind` | `spurious_*`, `qex_rM`, `qex_rB` |
| QTDA_export.jl | paths/names | `QTDA_export`, `qex_*_export` |
| QTDA_external.jl | **new** `QTDA_import_binary(Params)` | none |

Changes:
* **Reference** (`R = Params.Calc.QTDA.Reference`): precondition (QTDA.jl:21-26) is now
  `isfile(io/<Path>/<R>/bin/C.bin) && .../U.bin && .../V.bin`, else print the same error and return.
  Inputs (QTDA_solver.jl:9-22): `io/<Path>/<R>/bin/{C,U,V,qpH1b,qpH2b}.bin`; `QTDA_transitions_export` imports
  `io/<Path>/<R>/bin/C.bin` (QTDA_export.jl:292).
* Layout: `rm -r io/<Path>/QTDA`; mkdir `QTDA`, `QTDA/bin`, `QTDA/densities`, `QTDA/transitions/{E0,E1,E2,E3}`,
  `QTDA/spectra`. Files: `QTDA/QTDA_summary.dat`, `QTDA/bin/QTDA_X.bin`, `QTDA/bin/QTDA_E.bin` (layout unchanged:
  J=0..2Nmax+1 outer, P=1,2, column-major X[qp,nu], E[nu]), `QTDA/spectra/QTDA.dat`, `QTDA/spectra/QTDA_plot.dat`
  (lowercase p), `QTDA/transitions/E*/QTDA_*.dat` (names unchanged).
* Renames: `OBDM_rN` -> `rho_rN`, `OBDM_chR2` -> `rho_chR2`, `qpN1B_make` -> `qpN1b_make`, `qpO1B_import` -> `qpO1b_import`.
* `QTDA_import_binary` (QTDA_external.jl): reads back `QTDA/bin/QTDA_{X,E}.bin` (per JP block `read!` column-major);
  only consumer would be external analysis routines (none in this group). Port only if the transition-density/current
  tools of the TDA/RPA group are extended to QTDA; otherwise optional (~30 lines).
Fortran (nj_qtda + nj_qex): replace `io_path(Params,'/Bin/...')` by `io_path(Params,'/'//trim(Ref)//'/bin/...')`
in `QTDA_solver` and `qex_transitions_export` (needs the Reference passed in: add argument `Ref`); dirs; summary file
`Label//'_summary.dat'`; transitions dir `'/transitions/'`; spectra dir `'/spectra/'`; plot file `QTDA_plot.dat`
vs `QRPA_Plot.dat` (pass the file name, they differ in case); binaries to `'/QTDA/bin/'`. ~40 lines.
Physics changes reaching QTDA through shared code: H^22 pp/nn unchanged; H^40 not used; `qpTr1b_initialize` changed
in Tr1b group (T_odd phase for time-odd E1 current operators in the qp20 part, E0/EX 11-part formula, new F/GT
operators) -> owned by the transitions group, listed in §7.

## 5. QRPA (`dev:QRPA/*`) — Fortran nj_qrpa (+ nj_qex)

| dev file | base | Fortran |
|---|---|---|
| QRPA.jl / QRPA_solver.jl | paths, Reference | `QRPA`, `QRPA_solver` |
| QRPA_allocate / diagonalize / A_Phonon / energy / OBDM / spurious / transitions | `.i` -> `.Ind`, removed comments | `qex_allocate`, `QRPA_diagonalize`, `solve_block`, `QRPA_A_phonon_particle_number`, `QRPA_energy`, `QRPA_OBDM`, `qex_rM`, `qex_rB` |
| QRPA_export.jl | paths | `QRPA_export` |
| QRPA_external.jl | **new** `QRPA_import_binary(Params)` (ComplexF64 X,Y,E) | none |

* Reference: identical scheme to QTDA with `Params.Calc.QRPA.Reference` (QRPA.jl:21-26; QRPA_solver.jl:9-19;
  QRPA_export.jl:473 C import for transitions).
* Layout: mkdir `QRPA`, `QRPA/bin`, `QRPA/densities`, `QRPA/transitions/{E0..E3}`, `QRPA/spectra`, `QRPA/amplitudes`.
  Files: `QRPA/QRPA_summary.dat`; `QRPA/bin/QRPA_{X,Y,E}.bin` (ComplexF64, layout unchanged);
  `QRPA/spectra/QRPA.dat`, `QRPA/spectra/QRPA_Plot.dat` (capital P kept!);
  `QRPA/amplitudes/QRPA_amplitudes_ordered.dat`, `QRPA/amplitudes/QRPA_amplitudes.dat`;
  `QRPA/densities/QRPA_radial_densities.dat`; `QRPA/transitions/E3/QRPAE3.dat` (typo kept), `Mev` kept.
* QRPA results change through §3.3 (B matrix uses H^40 pp/nn via `qpO2b_40_pp/nn`) — i.e. only when the reference
  qpH2b.bin was produced by the new code.
* `QRPA_import_binary`: optional, as for QTDA (~35 lines).
Fortran: paths as §4 (~35 lines). Port-note B items (complex-eigvec stop, E3 name, Mev) still valid.

---------------------------------------------------------------------------------------------------

## 6. pnQTDA — entirely new solver (`dev:pnQTDA/*`), new Fortran module `src/qp/nj_pnqtda.f90`

Charge-exchange (proton-particle/neutron-quasiparticle) QTDA on top of an HFB/BCS reference: 2qp states
`|a_p b_n; J P>`, all (a,b) pairs, A matrix = diagonal 2qp energies + H^22 pn(11,11). No B, no spurious removal,
no transitions, no densities.

### 6.1 Driver `pnQTDA(Params)` (pnQTDA.jl)
```
wigner_init_float(75,"Jmax",9)
print header (same text pattern as QTDA: "Starting spherical pnQTDA calculation with residual NO2B NN+NNN interaction
  based on available reference quasiparticle mean-field calculation", "Calculation data:", thread report, Path)
R = Params.Calc.pnQTDA.Reference
if !isfile(io/<Path>/<R>/bin/qpH1b.bin) || !isfile(io/<Path>/<R>/bin/qpH2b.bin):
    println("\nError! ... No precomputed reference quasiparticle mean-field (HFB, BCS) solution is available in given
             Input_File path ... run HFB/BCS solver first ..."); return
rm -r io/<Path>/pnQTDA (if exists); mkdir pnQTDA, pnQTDA/bin, pnQTDA/spectra
"Starting pnQTDA calculation ..." ; pnQTDA_solver ; "Spherical pnQTDA solver run completed ..."
```
### 6.2 `pnQTDA_solver` (pnQTDA_solver.jl)
`Orb = orbitals_make(Params)`; `Orb_pn2qp = orbitals_pn2qp_make(Params,Orb)`;
`H_N = qpO1b_import(io/<Path>/<R>/bin/qpH1b.bin)`; `H_NN, Orb_NN = qpO2b_import(io/<Path>/<R>/bin/qpH2b.bin)`;
`A = pnQTDA_allocate`; `E, X = pnQTDA_diagonalize`; `pnQTDA_export`. (C, U, V not read.)

### 6.3 `orbitals_pn2qp_make` (pnQTDA_orbitals.jl)
`a_max` from Nmax, `J_max = 2Nmax+1`. Order: `for a in 1:a_max (proton), for b in 1:a_max (neutron),
P = parity(l_a+l_b), for J in |j_a-j_b|/2:(j_a+j_b)/2` -> append `Ind2B(a,b,0)` to channel (P,J). No a>=b restriction,
no N2max truncation, no a==b/J-even restriction. Returns `Orb2B_2qp(N,Ind,Dic)` (Dic unused).
Fortran: reuse `qporb2b_t` with `T = 0` (~35 lines).

### 6.4 `pnQTDA_allocate` (pnQTDA_allocate.jl)
For each (J,P) (threaded over JP), `N = Orb_pn2qp.N[P,J+1]`:
```
A[i,j] (j <= i, mirrored) = delta(a,c)delta(b,d) * (H_N.qp11.p[a,a] + H_N.qp11.n[b,b]) + qpO2b_1111_pn(a,b,c,d,J,P,H_NN,Orb_NN)
```
with (a,b) = state i, (c,d) = state j. No normalization factors.
**Domain issue**: `qpO2b_1111_pn` looks up `Orb_NN` T=0 pairs, which exist only for `2(n_a+n_b)+l_a+l_b <= N2max` and
`J <= N2max+1`; pn2qp contains all pairs up to `2Nmax`. For `N2max < 2Nmax` Julia throws `KeyError` (O2b_lookup).
Fortran: fatal with a clear message in that case (or decide with the owner, §open questions).

### 6.5 `pnQTDA_diagonalize`
For each JP (serial): `E, X = eigen!(A[P,J+1], sortby=+)` — A is exactly symmetric, Julia dispatches to the
symmetric solver: Fortran `eig_sym(A,E,X)` ascending. "no spurious Goldstone levels in the pn channel".

### 6.6 `pnQTDA_export` (pnQTDA_export.jl)
1. `pnQTDA/pnQTDA_summary.dat` (lines 17-84):
```
Spherical proton-neutron Quasiparticle Tamm-Dancoff Approximation review:
\nNuclid data:    A = <A>, Z = <Z>
\nCalculation data:    hw = <jstr(hw)> , N_max = <Nmax>, s.p. J-basis size = <a_max>, s.p. M-basis size = <..>
\t\t\tCenter-of-Mass System correction is set to:   CMS = ''<CMS>''
\n\t\t\t<CMS sentence: CMS1+2B | CMS2B | none>  (same three sentences as QTDA_summary)
\nReference quasiparticle mean-field solution:   Reference = ''<R>''
\nInformation on dimensions of pnQTDA 1-phonon subspaces ...
\nTotal number of 1-phonon pn quasiparticle states in J-scheme = <sum N>
Total number of 1-phonon pn quasiparticle states in M-scheme = <sum N*(2J+1)>
\nJ\tN_phonon            then rows "J\t(N[+]+N[-])" for J=0..2Nmax+1
\nSubspace\tP = +        rows "J\tN[1,J]"
\nSubspace\tP = -        rows "J\tN[2,J]"
```
2. Binary (86-119): `pnQTDA/bin/pnQTDA_X.bin`: for J=0..2Nmax+1, P=1,2: `write(X[P,J+1])` (N x N Float64,
   column-major, i.e. qp fastest); `pnQTDA/bin/pnQTDA_E.bin`: same loop, `E[P,J+1]` (N Float64). Stream, no headers.
3. Spectrum `pnQTDA/spectra/pnQTDA.dat` (121-153): header `@printf("%-5s %-5s %-10s\n","J","P","E")`; for J, P:
   rows `@printf("%-5d %-5s %-10.6f\n", J, "+"/"-", E)`, then `println(f,"\n")` (=> two newlines) after every (J,P).
4. Plot `pnQTDA/spectra/pnQTDA_plot.dat` (155-228): N = sum N (no ground-state row), columns (J,P,E,E); stable
   `sortperm(E)`; label shift `for nu=2:N: if E4[nu] < E4[nu-1] or |E4[nu]-E4[nu-1]| < 1.0: E4[nu] = E4[nu-1] + 1.0`;
   header `"%-5s %-5s %-12s %-12s\n"` (J,P,E,label); rows `"%-5s %-5s %-12.6f %-12.6f\n"` with J integer, P "+"/"-";
   **all N rows written** (unlike QTDA/QRPA which drop the highest level).
5. `pnQTDA_import_binary` (pnQTDA_external.jl): reader of the two binaries (optional, as QTDA).
Fortran module size: ~300 lines (driver 40, solver 25, orbitals 35, allocate 40, diag 20, export 120, import 30 optional).
Can reuse `qex_header`-style printing and the QTDA spectrum writer with a flag for "no ground state / write all".

---------------------------------------------------------------------------------------------------

## 7. Shared-module interface changes needed by this group

1. `nj_types.f90` `orb1b_t`: `pO, nO` integer -> `real(dp)` (dev `Orb1B.pO/nO::Float64`). nj_bcs comparisons
   `== 1` / `== 0` -> `== 1.0_dp` / `== 0.0_dp` (HF/core group owns the type change).
2. `nj_types.f90`: Pandya type with four components `pp, pnnp, pnpn, nn` (rename `o2b_temp_t%pn` -> `pnnp` or new
   `o2b_pandya_t`); `nj_o2b.f90` `O2b_Pandya_transformation`: compute `pnpn = sum_I -(2I+1) f6j(j_a,j_b,2J,j_c,j_d,2I)
   O2b_pn(a,d,c,b,I,Q)` in addition (§3.2).
3. `nj_functions.f90` (or nj_linalg): `hungarian(Cost, Assignment)` verbatim port of dev:common/functions.jl:711-858
   (also needed by `OBDM_natural_orbitals`).
4. `nj_linalg.f90`: symmetric indefinite solve (`dsysv`, uplo='U') for the Broyden gamma if not present; inverse
   square root of an SPD matrix (via `eig_sym`) for the Löwdin step (can be local to nj_hfb).
5. HF group: `HF_orbitals_make(Params,SPE,C,Orb)` (`Orb` intent(inout), returns reordered SPE, C; fatal on ESH failure)
   replacing `HF_orbital_ordering`, and dev `HF_solve` (calls it each iteration, mutates `Orb` occupations).
6. `nj_params.f90`: `io_path` prefix `'io/'`; `hfb_params_t` (drop Broy_Qmax/Broy_Bmax, Broy_Amax = 0.8);
   `qtda_params_t%Reference`, `qrpa_params_t%Reference` (default 'HFB'); new `pnqtda_params_t` + `Calc%pnQTDA`.
7. `nj_input.f90`: `&hfb` namelist without Broy_Qmax/Broy_Bmax; `&qtda`/`&qrpa` + `Reference`; new `&pnqtda Reference`.
8. `nj_obdm.f90` renames `OBDM_export/OBDM_chR2/OBDM_rN` -> `rho_export/rho_chR2/rho_rN` (optional; owned by core group).
9. `nj_o1b.f90` renames `qpO1B_export/import`, `qpN1B_make` -> `qpO1b_*`, `qpN1b_make` (optional).
10. `nj_tr1b.f90` `qpTr1b_initialize` (owned by the transitions group, consumed by QTDA/QRPA): dev adds `T_odd`
    phase (`-1` on the `V_ka U_lb` term of qp20) for E1_VC/VS/TC/TS/sTC/sTS, changes the qp11 part to
    `O_kl U U + (-1)^((j_a+j_b)/2+X) O_lk V V` (drops the E0 c-number term), and adds `F`, `GT` components to `qpTr1B`.
11. `nj_orbitals.f90`: dev `orbitals_make` fills occupations fractionally (j = 2l+1 first within (n,l)); owned by
    HF/core group; affects BCS reference only through §1.3.

---------------------------------------------------------------------------------------------------

## 8. Julia issues seen in dev & status of port notes

New/remaining issues in dev:
1. `HFB_diagonalize_BCS` -> `HF_orbitals_make` throws for `HF.ESH = true` with an open-shell target (HFB/BCS-pairing
   is meant for open shells); and mutates HFB's `Orb` occupations as a side effect (harmless).
2. BCS reference HF now runs with the target's (fractional) occupations because dev `HF_solve` re-fills with
   `Params` (target) — `Params_Ref` only seeds the first density. Probably unintended (open question).
3. `HF_BCS_initialize_chemical_potential`: `Lambda = 0.5*(Eps_hole - Eps_particle)` (pre-existing; a mid-gap value
   would be `+`), and Float64 `== 1`/`== 0` tests now skip fractionally occupied levels.
4. pnQTDA: pair set not truncated by N2max -> `KeyError` in `qpO2b_1111_pn` when `N2max < 2Nmax`.
5. `qpO2b_export` still loops `J = 0:2Nmax+1` over arrays sized `N2max+1` (pre-existing) -> port-note B.3 still needed;
   the new import size check must use the port's channel range.
6. `BCS()` and `HFB()` both `rm -r io/<Path>`, so `Reference="BCS"` and `"HFB"` cannot coexist under one Path
   (design limitation; mirror).
7. QTDA/QRPA plot export still writes only N of the N+1 sorted entries (drops highest level); pnQTDA writes all.
8. HFB_Broyden.jl header says `Broy_Amax` default 0.9, struct default is 0.8 (doc only; use 0.8).
9. QTDA plot file `QTDA_plot.dat` vs QRPA `QRPA_Plot.dat` inconsistent case (mirror).

Port notes (PORTING_NOTES.md, chain B):
* B.1 Broyden net effect -> **obsolete** (dev implements mixing; port §2.6).
* B.2 `V_NN_Ref` unused -> still valid.
* B.3 qpO2b export/import channel range -> still needed (+ size check consistency).
* B.4 silenced output -> still valid (dev still uses `redirect_stdout` in ScBCS and degeneracy branches).
* B.5 failed orderings fatal -> **obsolete** (matching + Hungarian never fails).
* Kept-on-purpose list: MCA Delta, degeneracy `Rho_1`, LN Delta, QRPA complex eigvecs, `QRPAE3.dat`/`Mev`: all still
  present in dev -> keep.
* "Not ported" list: `HFB_diagonalize_old` and `qpO2b_*_allocate_canonical` deleted in dev (nothing to do);
  `qpO2b_1111_index`/`qpO2b_1111_pn` **must now be ported** (pnQTDA); `qpO2B_31.jl` still unused.
* New fix vs base: wrong H shift at L=1 for iterations > 0 (§2.3) — Fortran currently reproduces the base behaviour.

---------------------------------------------------------------------------------------------------

## 9. Work list (Fortran, by file) and size

| file | work | lines |
|---|---|---|
| nj_hfb.f90 | paths/dirs; Broyden type + 8 routines; chemical potential secant/gap; HFB_solve inner loop; orbital matching + Löwdin replacing `order_one`/`canonical_order`; HF_orbitals_make in HFB_diagonalize_BCS | ~520 |
| nj_bcs.f90 | paths/dirs; HF_orbitals_make; Orb intent(inout); real occupation tests | ~40 |
| nj_qpo2b.f90 | Sgn fix qp40 pp/nn; engine Fch=4/Out=7 pn1111; `qpO2b_1111_pn`; import size check | ~70 |
| nj_o2b.f90 (core) | Pandya pnpn | ~25 |
| nj_functions.f90 (core) | hungarian | ~110 |
| nj_qtda.f90, nj_qrpa.f90, nj_qex.f90 | Reference paths, io layout, file names | ~80 |
| nj_pnqtda.f90 (new) + files.mk + dispatcher | pnQTDA | ~320 |
| nj_params / nj_input | params & namelists | ~35 |
| optional `*_import_binary` (QTDA/QRPA/pnQTDA) | readers | ~90 |
| **total** | | **~1200 (1100 without optional readers)** |
