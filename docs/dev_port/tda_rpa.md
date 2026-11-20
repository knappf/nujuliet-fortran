# Dev-port spec: TDA + RPA (non self-consistent)

Group: **TDA + RPA (plain, non self-consistent)**.
Julia base: `NuJuliet` @ 2595e82 (`Source/HF_RPA/*`). Julia dev: `NuJuliet_dev` @ f3f7a33 (`src/TDA/*`, `src/RPA/*`).
Fortran now: `src/rpa/nj_hf_rpa*.f90`, `nj_dispatch_rpa.f90` (+ RRPA files to be removed).
Out of scope (other group): the HF-ScRPA iteration (`HF_ScRPA_solver`, `RPA_rho.jl` except `RPA_OBDM_oY2`,
`RPA_sigma.jl`, `RPA_Hamiltonian.jl`, `RPA_orbitals.jl`). Where the plain RPA path calls into them it is described below.

Line references `file:NN` are dev lines (`NuJuliet_dev/src/...`).

---------------------------------------------------------------------------------------------------------------------

## 0. Executive summary

* Base `HF_RPA` ran TDA **and** RPA together on the HF solution read from `IO/<Path>/Bin/{C_HF,V2B_HF,h_HF}.bin`
  and wrote everything into `IO/<Path>/RPA/` + `IO/<Path>/Bin/`.
* Dev splits this into two independent solvers, `TDA(Params)` and `RPA(Params)`, each with its own parameter struct
  (`TDA_Parameters`, `RPA_Parameters`), reading a **reference directory** `io/<Path>/<Reference>/bin/` (`orb.bin`,
  `C.bin`, `rho.bin`, `H1b.bin`, `H2b.bin`, optional `sigma.bin`) and writing to `io/<Path>/<Name>/`.
* The formalism is generalized to **fractional reference occupations** (correlated references such as `HF-MBPT`
  natural orbitals): 1p-1h space from occupation tolerances, metric `N = diag(n_h - n_p)`, a general 1-body
  Hamiltonian `h_N` (not just diagonal HF energies), optional reference 2-body correlation function `Sigma`
  (box diagrams added to A and B).
* New diagonalization (non-symmetric `H = H_p H_m` product form), new normalization with the metric `N`, new
  spurious-state treatment (closed-form Gram-Schmidt + explicit 1x1 / 2x2 spurious-subspace solution, zero-norm
  Goldstone mode set to `X = Y = 0`).
* New OBDM (`RPA_OBDM_oY2`, all states, weights `(n_h - n_p)^2`, race-free), new correlation energy (quasi-boson with
  metric, physical states only, optional density energy), new transition matrix elements (`Transitions = "Linear"` or
  `"Exact"`, explicit `M[h,p]` in the Y term, 2-body Sigma part), fixed isovector E1 V/T/C strengths (`f3f7a33`).
* Output layout, file names, binary layout (now **row-major** amplitudes), summary contents all changed.
* External (post-processing) routines rewritten around a generic transition density matrix `Gamma = D_X X + D_Y Y`;
  new ones (`*_transition_operator_densities`, `*_transition_analysis`, `*_transitions_recalculate`) exist for TDA and RPA.
* Port notes 6 and 9 are obsolete in their present form (see section 7).

---------------------------------------------------------------------------------------------------------------------

## 1. File map (dev file -> base counterpart -> Fortran now)

| dev file | base counterpart | Fortran now | action |
|---|---|---|---|
| `TDA/TDA.jl` | part of `HF_RPA/HF_RPA.jl` | `nj_hf_rpa_driver.f90::HF_RPA` | new driver `TDA` |
| `TDA/TDA_solver.jl` | part of `HF_RPA_Solver.jl` | `HF_RPA_solver` | new `TDA_solver`, `TDA_solver_import` |
| `TDA/TDA_allocate.jl` | A part of `HF_RPA_Allocate.jl` | `HF_RPA_allocate` | rewrite (+ new `TDA_allocate_TBCF`) |
| `TDA/TDA_diagonalize.jl` | TDA part of `HF_RPA_Diagonalize.jl` | `HF_RPA_diagonalize::TDA_channel` | rewrite |
| `TDA/TDA_spurious.jl` | `HF_RPA_Spurious.jl` | `Spur_Ini`, `Spur_Ortho` | rewrite (identical to RPA_spurious, share) |
| `TDA/TDA_transitions.jl` | TDA part of `HF_RPA_Transitions.jl` | `HF_RPA_rM`, `HF_RPA_rB` | rewrite (exact kernel) |
| `TDA/TDA_export.jl` | TDA parts of `HF_RPA_Export.jl` | `nj_hf_rpa_export.f90` | rewrite |
| `TDA/TDA_external.jl` | TDA parts of `HF_RPA_External.jl` | `nj_hf_rpa_external.f90` (With_TDA) | rewrite + 3 new routines |
| `RPA/RPA.jl` (`RPA`) | `HF_RPA.jl` | `HF_RPA` | new driver `RPA` (`HF_ScRPA` = other group) |
| `RPA/RPA_solver.jl::RPA_solver` | `HF_RPA_Solver.jl` | `HF_RPA_solver` | rewrite |
| `RPA/RPA_allocate.jl::RPA_allocate`, `RPA_allocate_TBCF` | `HF_RPA_Allocate.jl` | `HF_RPA_allocate` | rewrite |
| `RPA/RPA_allocate.jl::RPA_allocate_SE`, `RPA_allocate_SE_Old` | new | none | DEV-only, see 4.4 |
| `RPA/RPA_diagonalize.jl` | `HF_RPA_Diagonalize.jl` | `RPA_eigen_solve`, `RPA_normalize`, `HF_RPA_diagonalize::RPA_channel` | rewrite |
| `RPA/RPA_spurious.jl` | `HF_RPA_Spurious.jl` | `Spur_Ini`, `Spur_Ortho` | rewrite (shared with TDA) |
| `RPA/RPA_energy.jl` | `HF_RPA_Energy.jl` | `HF_RPA_energy` | rewrite |
| `RPA/RPA_rho.jl::RPA_OBDM_oY2` | `HF_RPA_OBDM.jl` | `HF_RPA_OBDM` | rewrite (rest of file: ScRPA group) |
| `RPA/RPA_transitions.jl` | `HF_RPA_Transitions.jl` | `HF_RPA_rM`, `HF_RPA_rB` | rewrite rM, patch rB |
| `RPA/RPA_export.jl` | `HF_RPA_Export.jl` | `nj_hf_rpa_export.f90` | rewrite |
| `RPA/RPA_external.jl` | `HF_RPA_External.jl` | `nj_hf_rpa_external.f90` | rewrite + 3 new routines |
| (removed) | `HF_RPA_Collectivity.jl`, `HF_RPA_Phonon_Count.jl` | `HF_RPA_phonon_count` | collectivity: drop; phonon_count: keep only for EMPM |
| (removed) | `HF_RRPA/*` | `nj_hf_rrpa*.f90` | delete (dev removed HF_RRPA) |

`RPA_sigma.jl`, `RPA_rho.jl` (oY4, sigma, HF_ScRPA_OBDM_solve), `RPA_Hamiltonian.jl`, `RPA_orbitals.jl`:
ScRPA group. Plain RPA uses only `RPA_OBDM_oY2` and (for `TBCF = true`) `RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y)`
with defaults `beta = 1.0`, `Report = true`.

---------------------------------------------------------------------------------------------------------------------

## 2. Parameters

### 2.1 `TDA_Parameters` (new; `import_structures.jl:32`)

| field | default | meaning / use |
|---|---|---|
| `Name::String` | `"TDA"` | output dir `io/<Path>/<Name>`; must differ from `Reference` |
| `Reference::String` | `"HF"` | reference dir `io/<Path>/<Reference>/bin/` |
| `Reference_TBCF::Bool` | `false` | read `<Reference>/bin/sigma.bin`, add box diagrams to A (`TDA_allocate_TBCF`) and the 2-body part of the exact transition kernel; error if file missing |
| `Ortho::Bool` | `true` | spurious 1- orthogonalization |
| `Tol_hol::Float64` | `0.7` | orbital is a hole if occupation `> Tol_hol` |
| `Tol_par::Float64` | `0.1` | orbital is a particle if occupation `< Tol_par` |

TDA transitions are always "Exact" (w.r.t. the reference ground state).

### 2.2 `RPA_Parameters` (base had only `Ortho`; `import_structures.jl:41`)

| field | default | used by plain RPA? | meaning |
|---|---|---|---|
| `Name` | `"RPA"` | yes | output dir |
| `Reference` | `"HF"` | yes | reference dir |
| `Reference_TBCF` | `false` | yes | if true **and** `<Reference>/bin/sigma.bin` exists: Sigma box diagrams in A, B (`RPA_allocate_TBCF`), Sigma_AB energy term; silently ignored if the file is missing |
| `RPA_Orb` | `false` | yes (export) | export RPA natural orbitals as a new reference basis (`orb/C/rho/H1b.bin` in NAT basis) |
| `Ortho` | `true` | yes | spurious 1- treatment |
| `Transitions` | `"Linear"` | yes | `"Linear"` = `<0|[Q,M]|0>` with Orb occupations; `"Exact"` = `<0|Q M|0>` with RPA OBDM (+ RPA Sigma if `TBCF`); anything else -> error |
| `Tol` | `1e-5` | no (ScRPA) | |
| `Tol_hol` | `0.7` | yes | 1p-1h space |
| `Tol_par` | `0.1` | yes | 1p-1h space |
| `Imax` | `150` | no (ScRPA) | |
| `OBDM_diagonal` | `false` | no (ScRPA) | |
| `OBDM_alpha` | `0.6` | no (ScRPA) | |
| `ScH` | `true` | no (ScRPA) | |
| `TBCF` | `false` | yes | compute RPA 2-body correlation function `RPA_sigma` after the diagonalization; enables density energy; 2-body part of "Exact" rM; exports `bin/sigma.bin` |
| `TBCF_alpha` | `0.8` | no (ScRPA) | |
| `TBCF_beta` | `5.0` | no (ScRPA) | |
| `DEV` | `false` | yes (only with TBCF-allocate) | adds `pf/nf` polarization contact terms in `RPA_allocate_TBCF` and calls `RPA_allocate_SE(...;Mode="OnShell")` in `RPA_allocate` and `RPA_energy` |

`RRPA_Parameters` is gone in dev (remove `rrpa_params_t`, `&rrpa`).

### 2.3 Meaning of `Reference`

`Reference` is just a sub-directory name of `io/<Path>/`. The solver reads from `io/<Path>/<Reference>/bin/`:

| file | content | basis |
|---|---|---|
| `orb.bin` | `a_max` records `Int64 a,n,l,j; Float64 pO,nO` (`orbitals_import`) | occupations = diagonal of reference OBDM |
| `C.bin` | LHO -> reference basis (`O1b_import`) | |
| `rho.bin` | reference OBDM | reference basis |
| `H1b.bin` | reference 1-body Hamiltonian `h_N` | reference basis |
| `H2b.bin` | reference residual 2-body interaction `V_NN` (`O2b_import(...,Make_Orb_NN=true)`) | reference basis |
| `sigma.bin` | reference 2-body correlation function (only if `Reference_TBCF`) | reference basis |

Producers in dev: `HF` (writes `io/<Path>/HF/bin/{orb,C,rho,H1b,H2b}.bin`, `HF_export.jl:96`, `HF_solver.jl:51`)
and `HF-MBPT` (natural orbitals, `{orb,C,H1b,rho,H2b}.bin` and `sigma.bin` if `HF.MBPT.TBCF`, `HF_MBPT.jl:2307-2313`).
An RPA output dir is **not** a complete reference (no `H2b.bin`; `H1b.bin` only with `RPA_Orb=true`).
With `HF` the occupations are 0/1, `h_N` is diagonal and the dev equations reduce to the base ones (up to ordering,
spurious treatment, normalization conventions listed below). With `HF-MBPT` occupations are fractional.

---------------------------------------------------------------------------------------------------------------------

## 3. Common building blocks (TDA and RPA)

### 3.1 1p-1h space `orbitals_1p1h_make` (`common/orb.jl:348`, kwargs `Tol_par`, `Tol_hol`)

Replaces `orbitals_ph_make` + `orbitals_one_phonon_make` + `HF_RPA_phonon_count`. For `J = 0..J_max` (`J_max = 2 Nmax+1`),
`P = 1,2`, loop **`Orb_a` (particle) outer, `Orb_b` (hole) inner over all a_max orbitals**; accept if
`parity(l_a+l_b) == P` and `triangle(j_a,j_b,2J)`; then
```
if pO_b > Tol_hol && pO_a < Tol_par && (pO_b - pO_a) > 1e-8  -> push Ind2B(a,b,-1)   (proton)
if nO_b > Tol_hol && nO_a < Tol_par && (nO_b - nO_a) > 1e-8  -> push Ind2B(a,b, 1)   (neutron)
```
(proton and neutron entries are interleaved per (a,b) pair). Structure `Orb2B_1p1h(N[P,J+1], Ind[P,J+1]::Vector{Ind2B}, Dic)`,
`Dic` key `P | (T+1)<<2 | J<<4 | a<<10 | b<<18` (used only by ScRPA `RPA_sigma`/`RPA_rho`).
**Ordering change vs base** (base: all proton pairs `p` outer / `h` inner over HF particle/hole lists, then neutrons, J innermost):
row order of X, Y, binaries, amplitude files, the spurious index all change. Base particle/hole classification was
`pO == 0` / `pO == 1`; dev uses tolerances (an orbital can be both particle and hole when `Tol_par > Tol_hol`).

### 3.2 Spurious 1- vector `TDA_spurious_ini` == `RPA_spurious_ini` (3-arg) (`RPA_spurious.jl:1`)
Channel `(J,P) = (1,2)`, over `Orb_1p1h.Ind[2,2]`:
```
Spur[i_ph] = (n_h - n_p) * TrOp.E1.t[a_p,a_h]      (t = p for T=-1, n for T=1)
|Spur[i]| < 1e-14 -> 0     (base: 1e-6, and no (n_h-n_p) factor)
Spur ./= norm(Spur)
```
The 4-arg method `RPA_spurious_ini(...,Sigma)` (`RPA_spurious.jl:40`) is dead code (undefined `J`, `P_pa`, `Orb_NN`,
malformed `! =` token) - **do not port**.

### 3.3 Orthonormal complement `TDA_spurious_ortho!` == `RPA_spurious_ortho!` (`RPA_spurious.jl:132`)
Closed-form Gram-Schmidt replacing base `Spur_Ortho` (same subspace in exact arithmetic; numerically different):
```
U .= 0; N_a = 0; Spur_ind = 0
for a = N_ph:-1:1
    N_next = N_a; N_a += Spur[a]^2
    c_a = 0; u2_a = 1
    if N_a > 0: c_a = Spur[a]/N_a; u2_a = N_next/N_a
    v2_a = 1 - Spur[a]^2
    if Spur_ind == 0 && (v2_a <= 0 || sqrt(u2_a/v2_a) < 1e-4): Spur_ind = a; continue
    norm_a = 1/sqrt(u2_a)
    U[b,a] = -norm_a*c_a*Spur[b]   for b = a+1..N_ph
    U[a,a] = sqrt(u2_a)
U[:,Spur_ind] = Spur
```
Callers then shift columns `Spur_ind+1..N_ph` one left and use `U_sub = U[:,1:N_ph-1]`.
Keep the old Fortran `Spur_Ortho` (EMPM uses it, see 8.3).

### 3.4 Loop / storage conventions
* Channel arrays are `[P,J+1]` (P outer index) in dev; recommend Fortran `(P,J)` with `J = 0..J_max` like `qporb2b_t`.
* `JP_initialize(J_max)` order = `J` outer, `P` inner; all exports and binaries use this order.
* `phase(n) = (-1)^n`, `hat(j2) = j2+1` (argument is doubled), `parity(L) = rem(L,2)+1`, `w6j = f6j` (doubled args).

---------------------------------------------------------------------------------------------------------------------

## 4. RPA (plain)

### 4.1 Driver `RPA(Params)` (`RPA/RPA.jl:1-70`)
1. `wigner_init_float(75,"Jmax",9)`, print header + thread report.
2. If `!isdir("io/<Path>/<Reference>")` -> print error, return.
3. (`if !isdir("io/<Path>") mkdir; return` - dead, unreachable after 2.)
4. `IO_RPA = io/<Path>/<Name>`; **`rm(IO_RPA, recursive=true)` if it exists (no guard)**; mkdir
   `IO_RPA`, `bin`, `densities`, `transitions/{E0,E1,E2,E3}`, `spectra`, `amplitudes`.
5. `RPA_solver(Params)`.

### 4.2 `RPA_solver` (`RPA_solver.jl:1-92`) call sequence
```
IO = io/<Path>/<Reference>
Orb  = orbitals_import(IO/bin/orb.bin)
C    = O1b_import(IO/bin/C.bin);  Rho = O1b_import(IO/bin/rho.bin);  h_N = O1b_import(IO/bin/H1b.bin)
V_NN, Orb_NN = O2b_import(IO/bin/H2b.bin, Make_Orb_NN=true)
Sigma = empty; if isfile(IO/bin/sigma.bin) && Reference_TBCF: Sigma = O2b_import(sigma.bin)
Orb_1p1h = orbitals_1p1h_make(Orb; Tol_par, Tol_hol)
TrOp = Tr1b_transformation(C, Tr1b_initialize(Orb, 1.0))
TBCF_allocate = Reference_TBCF && isfile(sigma.bin)
A, B, N = RPA_allocate(..., h_N, V_NN, Rho, Sigma, TBCF_allocate)
Sigma_AB = TBCF_allocate ? Sigma : nothing
E, X, Y, Stability = RPA_diagonalize(Orb, Orb_1p1h, A, B, N, TrOp)
Rho_ref = copy(Rho)
Rho = RPA_OBDM_oY2(Orb, Orb_1p1h, Y)                         # reference basis, diagonal
Sigma = TBCF ? RPA_sigma(Orb,Orb_1p1h,Orb_NN,X,Y) : empty      # ScRPA group routine
chR2 = rho_chR2(Orb, C, Rho)
TrOp = Tr1b_transformation(C, Tr1b_initialize(Orb, chR2))
E_corr = RPA_energy(..., E,X,Y, Rho, Rho_ref, h_N, Sigma, V_NN; Sigma_AB)   # Density defaults to TBCF
rM = RPA_rM(..., TrOp, Rho, Sigma, X, Y)
rB = RPA_rB(Orb_1p1h, rM, real.(E))
Rho = C*Rho*C'   (LHO basis)
RPA_export(Orb, Orb_1p1h, Orb_NN, Stability, E_corr, E, X, Y, rB, C, Rho, Sigma, h_N)
```
Differences vs base: no TDA, reference files instead of `C_HF/V2B_HF/h_HF`, OBDM before Sigma/energy,
energy after the TrOp re-init, rM uses Rho & Sigma.

### 4.3 `RPA_allocate` (`RPA_allocate.jl:1-185`)
Threads over JP (each writes its own block; no race). Per JP, **lower triangle** `i_qg <= i_ph`, mirrored.
For pair `ph = (p,h,T_ph)`, `qg = (q,g,T_qg)`: `dn_ph = n_h - n_p`, `dn_qg = n_g - n_q` (occupations of the isospin `T`).
1-body terms only if `T_ph == T_qg` (`h = h_N.p` or `.n`), `:79-101`:
```
A += 0.5*(dn_ph+dn_qg)*(h[p,q]*δ(h,g) - h[h,g]*δ(p,q))
B += 0.5*phase(J+(j_p+j_h)/2)*(dn_ph-dn_qg)*(h[p,g]*δ(h,q) - h[h,q]*δ(p,g))
N += δ(p,q)*δ(h,g)*dn_ph                        # N is diagonal
```
2-body A (`:107-127`), if `rem(l_p+l_g,2) == rem(l_q+l_h,2)`, `I = max(|j_p-j_g|,|j_q-j_h|)/2 .. min(j_p+j_g,j_q+j_h)/2`,
`P_pg = parity(l_p+l_g)`, `Amp0 = dn_ph*dn_qg*(2I+1)`:

| T_ph,T_qg | term |
|---|---|
| -1,-1 | `Amp0*phase((j_h+j_q)/2+I)*w6j(j_p,j_h,2J,j_q,j_g,2I)*O2b_pp(p,g,h,q,I,P_pg)` |
| -1, 1 | `Amp0*phase((j_h+j_q)/2+I)*w6j(j_p,j_h,2J,j_q,j_g,2I)*O2b_pn(p,g,h,q,I,P_pg)` |
|  1,-1 | `Amp0*phase((j_g+j_p)/2+I)*w6j(j_q,j_g,2J,j_p,j_h,2I)*O2b_pn(q,h,g,p,I,P_pg)` |
|  1, 1 | `Amp0*phase((j_h+j_q)/2+I)*w6j(j_p,j_h,2J,j_q,j_g,2I)*O2b_nn(p,g,h,q,I,P_pg)` |

2-body B (`:131-151`), if `rem(l_p+l_q,2) == rem(l_h+l_g,2)`, `I = max(|j_p-j_q|,|j_h-j_g|)/2 .. min(j_p+j_q,j_h+j_g)/2`, `P_pq`:

| T_ph,T_qg | term |
|---|---|
| -1,-1 | `Amp0*phase((j_h+j_q)/2+J+I)*w6j(j_p,j_h,2J,j_g,j_q,2I)*O2b_pp(p,q,h,g,I,P_pq)` |
| -1, 1 | same phase/6j, `O2b_pn(p,q,h,g,I,P_pq)` |
|  1,-1 | `Amp0*phase((j_g+j_p)/2+J+I)*w6j(j_q,j_g,2J,j_h,j_p,2I)*O2b_pn(q,p,g,h,I,P_pq)` |
|  1, 1 | as (-1,-1) with `O2b_nn` |

These are the base expressions times `dn_ph*dn_qg` (base A diagonal was `E_p - E_h`; base filled the full matrix
without mirroring - Fortran must mirror the lower triangle). Base `A`/`B` were `[J+1,P]`.

If `TBCF` (argument): `A, B = RPA_allocate_TBCF(...)`; if additionally `RPA.DEV`: `A = RPA_allocate_SE(...;Mode="OnShell")`.

`RPA_allocate_TBCF` (`:187-566`) - port literally:
1. `W1_NN, Orb_NN_temp = O2b_temp_initialize(Params,Orb,Make_Orb_NN=true)`; per JP of the temp space (`J = 0..N2max+1`),
   for every `Bra=(k,l)`, `Col=(b,d)` cross-coupled sums over
   `J_1 = max(|j_k-j_b|,|j_l-j_d|)/2 .. min(j_k+j_b,j_l+j_d)/2`, `Amp_J1 = (2J_1+1)*w6j(j_k,j_l,2J,j_d,j_b,2J_1)`, `Amp = phase(J_1)*Amp_J1`,
   `P_kb = parity(l_k+l_b)`, for both V (`V_NN`) and S (`Sigma`):
   `pp += Amp*O2b_pp(k,b,l,d)`, `pn += Amp*O2b_pn(k,b,l,d)`, `nn += Amp*O2b_nn(k,b,l,d)`,
   `np += phase((j_k+j_l)/2)*Amp*O2b_pn(b,k,d,l)`, `mpn += Amp_J1*O2b_pn(k,b,d,l)`.
   Then `W1.pp = ppXV ppXS' + ppXS ppXV' + pnXV pnXS' + pnXS pnXV'`, `W1.pn = mpnXV mpnXS' + mpnXS mpnXV'`,
   `W1.nn = nnXV nnXS' + nnXS nnXV' + npXV npXS' + npXS npXV'` (dgemm, `:223-236`).
2. `W2_NN`: direct `O[Bra,Ket] = <kl;J|O|mn;J>` for V and S, `W2.t = V S' + S V'` for t = pp, pn, nn (`:239-277`).
3. DEV only: `pf, nf` 1-body polarization potential (`:279-337`).
4. Per JP of the 1p-1h space, lower triangle, threads over `i_ph` (`:339-560`): DEV-only pf/nf contact terms
   (`Amp_A = -1`, `Amp_B = phase(1-J+(j_p+j_h)/2)`); then
   * A box 1: `I = max(|j_p-j_q|,|j_h-j_g|)/2..min(...)`, `Amp = phase(2+J+I+(j_h+j_g)/2)*(2I+1)*w6j(j_p,j_q,2I,j_g,j_h,2J)`
     (for T=(-1,1) phase uses `(j_h+j_q)/2`; for (1,-1) `phase(2+J+I+(j_g+j_p)/2)*w6j(j_q,j_p,2I,j_h,j_g,2J)` and
     `O2b_temp(q,p,g,h,I,P_pq,0,W1)`), `O2b_temp(p,q,h,g,I,P_pq,T,W1)` with T=-1 (pp), 0 (pn), 1 (nn).
   * A box 2: `I` over (p,g),(q,h); `Amp = phase(1+I+(j_h+j_q)/2)*(2I+1)*w6j(j_p,j_h,2J,j_q,j_g,2I)`;
     pp/nn: `0.5*Amp*O2b_temp(p,g,h,q,I,P_pg,∓1,W2)`, (-1,1): `Amp*O2b_temp(p,g,h,q,..,0,W2)`,
     (1,-1): `phase((j_p+j_h+j_q+j_g)/2)*Amp*O2b_temp(g,p,q,h,..,0,W2)`.
   * B box 1 and box 2 analogous (`:497-545`; quote exactly when porting).
   Add `A[P,J+1] .+= A_JP`, `B .+= B_JP`.

`O2b_temp(a,b,c,d,J,P,T,O,Orb_NN_temp)` (`common/O2b.jl:601`) = `O.{pp,pn,nn}[P,J+1][Map[a,b,J+1],Map[c,d,J+1]]`
for `T = -1,0,1` (throws if the pair is not in the temp space).

### 4.4 `RPA_allocate_SE` (DEV only, ~650 Julia lines, `RPA_allocate.jl:1232-1868`)
Energy-dependent 1p-1h <-> 2p-2h self-energy (modes "QP"/"OnShell"/"Fixed", regularization "DSRG"/"Imaginary"/"None",
kwargs `Gamma=2.0`, `Omega=0.0`, `dE_bin=Gamma/20`, `Chunk=32768`). Only reachable with `RPA.DEV = true` **and**
TBCF-allocation. `RPA_allocate_SE_Old` is unreachable ("for testing"). Recommendation: accept `DEV` in the namelist,
stop with "not ported" if `DEV=.true.` (owner decision, see open questions).

### 4.5 `RPA_diagonalize` (`RPA_diagonalize.jl`)
Serial over JP (shared work arrays `M1_temp`, `M2_temp` of size `maximum(Orb_1p1h.N)`). Returns
`E[P,J+1]::ComplexF64`, `X`, `Y::ComplexF64` (N_ph x N_ph), `Stability`.

**Non-orthogonalized channel** (`Ortho=false` or `(J,P) != (1,2)`), `:56-139`:
```
Nd = diag(N); s = 1 ./ sqrt.(Nd)
K  = s.^2 .* (A - B) .* s'                 # = N^-1 (A-B) N^-1/2
H  = (A*K + B*K) .* s                      # rows scaled: H = N^-1/2 (A+B) N^-1 (A-B) N^-1/2 = H_p H_m
E2, R = eigen!(H)                          # general real (geevx 'B'); symmetric H -> syevr; sorted by (re,im) of E2
E[nu] = sqrt(ComplexF64(E2[nu]))
XmY = s .* R ;  EXpY = K*R   (real and imaginary parts of R multiplied separately)
X[:,nu] = 0.5*EXpY/E + 0.5*XmY ;  Y[:,nu] = 0.5*EXpY/E - 0.5*XmY
norm_nu = sum_ph Nd[ph]*(|X|^2 - |Y|^2)
|norm_nu| > 1e-12 : X,Y *= 1/sqrt(|norm_nu|)          (negative-norm states keep negative norm, NO -i factor)
else              : X = Y = 0, print "Zero-norm RPA state detected, setting X = Y = 0 ..."
```
Base: three branches (A-B posdef / A+B posdef / complex), unit normalization of `|X|^2-|Y|^2` without metric,
`-i` factor for negative norm, `|X|^2+|Y|^2 = 1` for near-zero norm. All of that is gone.
Fortran `eig_general` already reproduces Julia's `eigen!` dispatch (symmetric -> dsyevr; else dgeevx 'B', sorted by (re,im)).

**Orthogonalized 1- channel**, `:141-361`:
```
Spur = RPA_spurious_ini(Orb,Orb_1p1h,TrOp);  U <- RPA_spurious_ortho!;  shift columns;  U_sub = U[:,1:N_ort], N_ort = N_ph-1
N_ort_JP = U_sub' * (Nd .* U_sub)
λ, V = eigen!(Symmetric(N_ort_JP))   (syevr, uplo U)
V .*= λ'.^(-1/4);  inv_sqrtN_ort = V*V'
W = U_sub*inv_sqrtN_ort                       (N_ph x N_ort)
H_A = W'*A*W ; H_B = W'*B*W ; H_p = H_A + H_B ; H_m = H_A - H_B      (computed as H_p .+= H_m; H_m .= H_p .- 2 H_m)
E2, R = eigen!(H_p*H_m);  E[2:N_ph] = sqrt(E2)
X-Y = W*R ;  E(X+Y) = W*(H_m*R)  -> X[:,2:], Y[:,2:] as above; normalization with the full-space Nd (same rule)
Spurious 2x2 subspace:
  A_s = Spur'*A*Spur ; B_s = Spur'*B*Spur ; N_s = Spur'*(Nd.*Spur) ; c = 1/sqrt(N_s)
  Hp_s = c*(A_s+B_s)*c ; Hm_s = c*(A_s-B_s)*c ; E_s = sqrt(ComplexF64(Hp_s*Hm_s))
  r = c ; p = c*Hm_s ; x = 0.5*p/E_s + 0.5*r ; y = 0.5*p/E_s - 0.5*r ; n_s = N_s*(|x|^2-|y|^2)
  print "Spurious 1- RPA Goldstone boson ... E = (E_s) MeV"
  if |E_s| > 1e-8 && |n_s| > 1e-8 : E[1] = E_s ; X[:,1] = (x/sqrt|n_s|)*Spur ; Y[:,1] = (y/sqrt|n_s|)*Spur
  if |E_s| < 1e-8 || |n_s| < 1e-8 : print zero-norm message ; E[1] = 0 ; X[:,1] = Y[:,1] = 0
```
Base: spurious stored as `E = 0`, `X = Spur/sqrt(2)`, `Y = -Spur/sqrt(2)`, physical solutions from the U-projected
A,B with the three-branch solver and no metric. Note: an imaginary `E_s` gives `n_s = 0` -> zeroed;
a real `E_s` gives a **regular normalized spurious state that enters OBDM, energy and Sigma** (see 7.2).

Stability: any `|imag(E)| > 1e-8` -> `"RPA INSTABILITY DETECTED in channel: J = $J, P = ±"` (base always printed `+`).

### 4.6 `RPA_OBDM_oY2` (`RPA_rho.jl:1-73`)
```
Rho = diag(Orb.pO), diag(Orb.nO)
per JP (threaded, own buffers pRho_JP[i_JP], nRho_JP[i_JP] of length a_max):
  for each 1p-1h pair i_ph: S = sum_{nu=1..N_ph} |Y[i_ph,nu]|^2 * (2J+1)        # ALL states, no norm cut
     w = (n_h - n_p)^2 ;  buf[h] -= w/(j_h+1)*S ;  buf[p] += w/(j_p+1)*S
Rho[a,a] += buf_JP[a]   summed over JP in JP_list order (J outer, P inner)
```
Base: `|Y|^2 < 0.5` cut, unsynchronized threaded update (port note 6). Dev is race-free; the Fortran should
accumulate per JP block then add blocks in JP order (bit-closest), serial is fine.

### 4.7 `RPA_energy` (`RPA_energy.jl`)
Signature `RPA_energy(Params,Orb,Orb_1p1h,Orb_NN,E,X,Y,Rho,Rho_ref,h_N,Sigma,V_NN;Sigma_AB=nothing,Density=Params.Calc.RPA.TBCF)`,
returns NamedTuple `(E_QBA, E_sigma, E_rho, E_dress, E_dens)` (NaN where not evaluated).
* `RPA_energy_QBA` (`:71`): per JP, `D_ph = n_h - n_p`; skip states with `|Σ D(|X|²-|Y|²)| < 1e-8` (zero norm,
  incl. a zeroed spurious state) and states with `|Im E| > 1e-8` (counted, warning printed);
  `E_QBA = -Σ_JP (2J+1) Σ_nu Re(E_nu) Σ_ph D_ph |Y[ph,nu]|²`. No explicit spurious skip (base skipped `nu=1` for 1- with Ortho).
* `Density = true`: `E_sigma, E_rho = RPA_energy_rho` with `E_sigma = V2b_correlation_energy(Orb,Orb_NN,Sigma,V_NN)`
  (`common/V2b.jl:572`, other group) and `E_rho = Σ_{a,b: l,j equal} (j_a+1)*(h.p[a,b]*(Rho.p-Rho_ref.p)[b,a] + same n)`;
  if `Sigma_AB !== nothing`: `dA,dB = RPA_allocate_TBCF(...,Sigma_AB, zero A, zero B)` (+ DEV SE),
  `E_dress = Σ (2J+1) Re(Y'dA Y + X'dB Y)` over the same physical states. `E_dens = E_sigma+E_rho+E_dress`.
* Prints (6-digit rounding) as in `:49-58`.

### 4.8 `RPA_rM` / `RPA_rM_calculate` / `RPA_rM_contract` (`RPA_transitions.jl:1-238`)
Operator lists per multipole: E0 `[E0]` (J,P)=(0,1); E1 `[E1,E1_VC,E1_VS,E1_TC,E1_TS,E1_sTC,E1_sTS,E1_C,E1_sC]` (1,2);
E2 `[E2]` (2,1); E3 `[E3]` (3,2). `rM = M_X*X + M_Y*Y` (dgemm/zgemm), coefficient matrices `M_X[i_op,i_ph]`.
Per ph pair (`T_ph`), `Phase_hp = phase((j_p-j_h)/2)`, operator of the same isospin `M_1`, opposite `M_2`:
* **Linear** (`:149-160`): `D = n_h - n_p`; `M_X1 = M_1[p,h]*D*phase(J)`; `M_Y1 = M_1[h,p]*D*Phase_hp`.
  Base: `M[p,h]*(phase(J)*X + Y)`. Identical for time-even (E-type) operators; **sign of the Y term flips for
  time-odd operators** (`E1_VC, E1_VS, E1_TC, E1_TS, E1_sTC, E1_sTS`, `M[h,p] = -(-1)^((j_p-j_h)/2) M[p,h]`) ->
  E1 vortical/toroidal strengths change w.r.t. base. Fortran must read `M[a_h,a_p]` from the full transformed matrix.
* **Exact** (`:163-221`): `G = (I - Rho)*M*Rho` per operator (Rho = RPA OBDM, reference basis);
  `M_X1 = G_1[p,h]*phase(J)`, `M_Y1 = -G_1[h,p]*Phase_hp`; if `RPA.TBCF`: for all `a,b` with `triangle(j_a,j_b,2J)`
  and `parity(l_p+l_b) == parity(l_h+l_a)`, `I = max(|j_h-j_a|,|j_p-j_b|)/2..min(j_h+j_a,j_p+j_b)/2`,
  `F = phase(J+I+(j_p+j_a)/2)*(2I+1)*w6j(j_p,j_h,2J,j_a,j_b,2I)`;
  proton ph: `K1 += O2b_pp(h,a,p,b,I,P_ha,Sigma)*F`, `K2 += O2b_pn(h,a,p,b,I,P_ha,Sigma)*F`;
  neutron ph: `K1 += O2b_nn(h,a,p,b)*F`, `K2 += O2b_pn(a,h,b,p)*phase((j_p+j_a+j_h+j_b)/2)*F`;
  `M_X{1,2} += phase(1+(j_p+j_h)/2)*M_{1,2}[a,b]*K{1,2}`, `M_Y{1,2} += phase(J+(j_p+j_h)/2)*phase((j_a-j_b)/2)*M_{1,2}[b,a]*K{1,2}`.
  `M_X2/M_Y2` (opposite isospin) are stored into the **other** isospin's coefficient arrays.
* Invalid `Transitions` -> error.

### 4.9 `RPA_rB` == `TDA_rB` (`RPA_transitions.jl:240-430`)
Same as base `HF_RPA_rB` (`hc=197.326980`, `g_p=5.586`, `g_n=-3.826`, `e_p=(A-Z)/A`, `e_n=Z/A`) **except the isovector
E1 fixes of f3f7a33** (`:351-366`):
```
rB_E1V [iv] = |0.5(VC.p-VC.n) + 0.125(g_p-g_n)(VS.p - VS.n)|^2         (base: VS.p + VS.n)
rB_E1VS[iv] = |0.125(g_p-g_n)(VS.p - VS.n)|^2
rB_E1T [iv] = |0.5(TC.p-TC.n) + 0.125(g_p-g_n)(TS.p - TS.n)|^2
rB_E1TS[iv] = |0.125(g_p-g_n)(TS.p - TS.n)|^2
rB_E1C [iv] = |0.5(TC.p-TC.n)+0.125(g_p-g_n)(TS.p-TS.n) - (0.5(VC.p-VC.n)+0.125(g_p-g_n)(VS.p-VS.n))|^2
rB_E1CS[iv] = |0.125(g_p-g_n)(TS.p-TS.n) - 0.125(g_p-g_n)(VS.p-VS.n)|^2
rB_E1_NLO_LWA[iv] = |e_p E1.p - e_n E1.n + 0.5 E/hc ((TC.p-TC.n) + 0.25(g_p-g_n)(TS.p - TS.n))|^2   (base: 0.125, +)
```
One shared Fortran routine for TDA and RPA.

### 4.10 `RPA_export` (`RPA_export.jl`) - order of calls
1. `RPA_summary` -> `RPA_summary.dat` (new file): title `Spherical Random-Phase Approximation review:`, nuclide, calc data,
   CMS lines, Ortho line (`Goldstone levels`), stability line, dimension tables (J-scheme total, M-scheme total,
   `J\tN_phonon`, `Subspace\tP = +`, `Subspace\tP = -`), then `"\n\nRPA ground-state correlation energy ..."`,
   `@printf "\nE_corr     = %12.6f \t MeV ... RPA ground-state correlation energy (quasi-boson)"` and, if `E_dens` not NaN,
   the `(density)` line (`:155-243`). Base had `E_corr` from the old formula.
2. `RPA_binary_export` -> `bin/RPA_X.bin`, `bin/RPA_Y.bin`, `bin/RPA_E.bin`: for J=0..J_max, P=1,2:
   `write(permutedims(X[P,J+1]))` = **row-major** (for ph: for nu: X[ph,nu]) ComplexF64 (re,im) 16 B each; `E` ComplexF64 vector.
   Base: `IO/<Path>/Bin/`, column-major (nu outer, ph inner). Fortran: stream-write `transpose(X)`.
3. `RPA_spectrum_export` -> `spectra/RPA.dat` (same format as base `RPA.dat`).
4. `RPA_plot_spectrum_export` -> `spectra/RPA_plot.dat` (was `RPA_Plot.dat`); writes **all N+1** levels (base wrote N).
5. `RPA_amplitudes_export` -> `amplitudes/RPA_amplitudes_ordered.dat` and `amplitudes/RPA_amplitudes.dat`
   (lowercase); `|X|^2`, `|Y|^2` are now metric-weighted `Σ (n_h-n_p)|X|^2` (base: plain norms; base energy-ordered
   file had an `N_ph` shadowing bug - fixed).
6. `rho_export(Params,Orb,Summary_File,Densities_File,Rho_LHO,C)` (common, other group) with
   `Summary_File = <Name>/RPA_summary.dat`, `Densities_File = <Name>/densities/<Name>_radial_densities.dat`.
7. `RPA_orbitals_export` (`:54-132`): `Rho_RPA = C'*Rho*C`; `_, C_NAT, _ = OBDM_natural_orbitals(Params,Orb,Rho_RPA)`
   (common/orb.jl:552, Hungarian ordering); `pO,nO = diag(C_NAT' Rho_RPA C_NAT)`, `pE,nE = diag(C_NAT' h C_NAT)`;
   appends to the summary: `"\nRPA single-particle orbitals ... natural orbitals of the RPA OBDM, E = <a|h|a> centroids ..."`,
   per Proton/Neutron table header `@printf "%4s %4s %4s %10s %14s\n" n l 2j O_a E`, rows sorted by `sortperm(E_t)`
   with `n` counted per (l,j) in energy order, footer `Particle number ...` (`%-14s %15.6f`);
   writes `densities/RPA_occupation.dat` (columns `a n l 2j pO_ref pO_RPA pO_dep nO_ref nO_RPA nO_dep`, `%-6d`/`%-12.6f`,
   `Orb_ref` re-read from `<Reference>/bin/orb.bin`).
   Binary: `RPA_Orb=false`: `bin/orb.bin` (Orb, reference occupations), `bin/C.bin` (C), `bin/rho.bin` (Rho_RPA);
   `RPA_Orb=true` (`RPA_natural_orbitals_export`): `orb.bin` (occupations of the NAT diagonal), `C.bin = C*C_NAT`,
   `rho.bin = diag`, `H1b.bin = C_NAT' h C_NAT`.
8. `TBCF`: `O2b_export(Orb_NN, Sigma, <Name>/bin/sigma.bin)` (RPA basis).
9. `RPA_transitions_export` (`:525-993`): `Rho_RPA = C'*Rho_LHO*C`; `pR2 = rho_rN(2,Rho_RPA,C,Orb)`, `pR4 = rho_rN(4,...)`
   (fixes base note "OBDM_rN with RPA density in the LHO basis"); appends moments `m_-1..m_3` (ph/is/iv for E0..E3) and TRK
   values to `RPA_summary.dat` (same text as base RPA block). **E1 moments skip `i_nu = 1` if `RPA.Ortho`** (base did not).
   Transition tables under `transitions/E0/RPA_E0.dat`, `transitions/E1/RPA_E1{,_V,_V_conv,_V_spin,_T,_T_conv,_T_spin,_C,_C_conv,_C_spin,_c,_NLO_LWA}.dat`,
   `transitions/E2/RPA_E2.dat`, `transitions/E3/RPA_E3.dat` (format `%-20s`x4 header, `%-20.8f`x4 rows, `real(E)`).

---------------------------------------------------------------------------------------------------------------------

## 5. TDA (new solver)

### 5.1 Driver `TDA(Params)` (`TDA/TDA.jl`)
Prints header incl. `TDA parameters: Name = ''..'' , Reference = ''..'' , Reference_TBCF = .. , Ortho = .. , Tol_hol = .. , Tol_par = ..`.
Checks: reference dir exists; `Reference_TBCF && !isfile(<ref>/bin/sigma.bin)` -> error/return;
`Name == Reference` -> error/return; if `IO_TDA` exists: remove only if it contains `.TDA_output` **or is empty**,
else error (the marker is never written - bug, 7.4). mkdir `IO_TDA`, `bin`, `densities`, `transitions/E0..E3`, `spectra`, `amplitudes`.

### 5.2 `TDA_solver` (`TDA_solver.jl`)
```
Orb, Orb_1p1h, C, Rho, h_N, V_NN, Orb_NN, Sigma = TDA_solver_import    # same reference files as RPA, Tol_par/Tol_hol of TDA
chR2 = rho_chR2(Orb,C,Rho)                                              # reference charge radius
TrOp = Tr1b_transformation(C, Tr1b_initialize(Orb,chR2))
A, N = TDA_allocate(...)            # A,N formulas of RPA_allocate (A and N parts only), + TDA_allocate_TBCF if Reference_TBCF
E, X = TDA_diagonalize(...)
rM = TDA_rM(...,TrOp,Rho,Sigma,X) ; rB = TDA_rB(Orb_1p1h,rM,E)
TDA_export(Orb,Orb_1p1h,E,X,rB,C,Rho)
```
`TDA_allocate` (`TDA_allocate.jl:1-137`) is the A/N part of `RPA_allocate` verbatim. `TDA_allocate_TBCF` (`:139-446`)
is the A part of `RPA_allocate_TBCF` verbatim (W1, W2, box 1 and box 2 of A; pf/nf block commented out, no DEV).
Share the Fortran code: one allocate routine with an `Want_B` switch, one TBCF routine with `Want_B` and `DEV` switches.

### 5.3 `TDA_diagonalize` (`TDA_diagonalize.jl`)
Non-ortho: `s = 1/sqrt(Nd)`, `H = s.*A.*s'`, `eigen!(Symmetric(H))` (dsyevr, ascending), `E = values`, `X = s.*vectors`.
Ortho 1-: `Spur`, `U_sub` as 3.2/3.3; `A_ort = U_sub'*A*U_sub`, `N_ort = U_sub'*(Nd.*U_sub)`;
`eigen!(Symmetric(A_ort),Symmetric(N_ort))` = LAPACK `sygvd(1,'V','U')`, sorted ascending (X'N X = 1);
`E[2:] = values`, `X[:,2:] = U_sub*vectors`; spurious: `E[1] = (Spur'A Spur)/(Spur'(Nd.*Spur))`,
`X[:,1] = Spur/sqrt(N_s)`, print `Spurious 1- TDA level ... E = (..) MeV`. Base TDA: spurious `E=0`, `X=Spur`.
Print `Negative TDA eigenvalue detected in channel: J = .., P = ±` if any `E < -1e-8`. `N_ph == 0` -> empty.

### 5.4 `TDA_rM` (`TDA_transitions.jl:1-217`) - always "Exact" w.r.t. the reference
`Pairs, D = TDA_transition_kernel(J,P,...,Rho,Sigma)`; `Gamma = D*X`; `prM = pM_pair*Gamma`, `nrM = nM_pair*Gamma`
(`pM_pair[i_op,i_pair] = M.p[a,b]` for proton pairs only, analogous neutrons). Kernel:
* `Reference_TBCF=false`: `Pairs = Orb_1p1h.Ind[P,J+1]` + extra pairs `(c,d,T_ph)` with
  `|δ(c,p) - Rho[p,c]|*|Rho[d,h]| > 1e-12` not yet listed (none for a diagonal Rho).
* `Reference_TBCF=true`: `Pairs` = all `(a,b,t)` with `t in (-1,1)` outer, then `a`, `b` over Orb, parity P, `triangle(j_a,j_b,2J)`.
* `D[(c,d,t_ph),ph] += (δ(c,p) - Rho[p,c])*Rho[d,h]*phase(J)` over `c` with (l_p,j_p), `d` with (l_h,j_h) (pairs not in the
  map skipped). For HF: `D[ph,ph] = (-1)^J` (= base TDA sign convention `M[p,h]*(-1)^J X`).
* TBCF: for each pair `(a,b)` with `T == T_ph`: `K1,K2` exactly as in 4.8 (same F, same O2b calls), `Phase_X = phase(1+(j_p+j_h)/2)`
  (fixed by f3f7a33, was `phase((j_p+j_h)/2)`), `D[(a,b,t_1),ph] += Phase_X*K1`, `D[(a,b,t_2),ph] += Phase_X*K2`.

### 5.5 `TDA_export` (`TDA_export.jl`)
1. `TDA_summary.dat`: title `Spherical Tamm-Dancoff Approximation review:`, nuclide/calc/CMS lines, `TDA parameters:` block
   (Name/Reference/Reference_TBCF/Ortho, `1p-1h space tolerances`, TBCF included/not, `Transition matrix elements ... ''Exact'' ...`),
   Ortho line (`spurious center-of-mass levels`), stability (`N_neg` = # of E < -1e-8 excluding spurious), lowest level
   `@printf "\tLowest TDA excitation energy:   E = %12.6f MeV   (J = %d, P = %s)\n"`, dimension tables. No energy, no radii.
2. `bin/TDA_X.bin` (Float64, row-major `permutedims`), `bin/TDA_E.bin` (Float64).
3. `spectra/TDA.dat` (`%-5s %-5s %-20s` / `%-5d %-5s %-20.6f`, blank lines between blocks), `spectra/TDA_plot.dat`.
4. `amplitudes/TDA_amplitudes_ordered.dat` and `amplitudes/TDA_amplitudes.dat` (energy-sorted): new columns
   `Ind E J P |X|^2 t p(n,l,2j) h(n,l,2j) max weight` with the dominant component (largest `(n_h-n_p)|X|^2`), formats `:338-350`.
5. `TDA_transitions_export`: identical to the RPA one except `Rho` is the reference OBDM already in the reference basis
   (no transform), output in `TDA_summary.dat`, E1 skip uses `TDA.Ortho`, TRK header text
   `"\nPhysical Thomas-Reiche-Kuhn (TRK) energy-weighted sum-rule values (reference values, not fulfilled exactly by TDA) ...\n"`;
   files `transitions/E*/TDA_*.dat` (same 15 names as RPA; base TDA lacked `_C_conv`, `_C_spin`, `_c`).
   `densities/` is created but stays empty.

---------------------------------------------------------------------------------------------------------------------

## 6. External / post-processing routines (`TDA_external.jl`, `RPA_external.jl`)

Public functions (exported from the module; Fortran: `&run` actions, see 8.2):

| Julia | kwargs (defaults) | reads | writes |
|---|---|---|---|
| `RPA_transition_densities` | `JP, nu_list, File_Name="Phonon_Densities"` | `<Name>/bin/{orb,C,rho,RPA_X,RPA_Y,RPA_E}.bin` (+`sigma.bin` if Exact&TBCF) | `densities/RPA_<File>_JP<JP>_nu<nu>_E<sE>.dat`, `..._JP<JP>_Averaged.dat` |
| `RPA_transition_operator_densities` | `File_Name="Phonon_Transition_Operator_Densities"` | same | same pattern (header `pM nM isM ivM`) |
| `RPA_transition_analysis` | `File_Name="Phonon_Transtion_Analysis"`, `Cut=0.1` | same | `levels/RPA_<File>_JP<JP>_nu<nu>_E<sE>.dat` (mkdir `levels`) |
| `RPA_transition_currents` | `File_Name="Phonon_Currents"`, `M=-1` (->0), `xN_grid=40`, `zN_grid=40` | same | `densities/RPA_{Convective,Spin}_<File>_JP<JP>_nu<nu>_E<sE>.dat`, `..._Averaged.dat` |
| `RPA_transitions_recalculate` | - | same | truncates `RPA_summary.dat` at `"\nReview of RPA electromagnetic transition moments"`, re-runs rM/rB/`RPA_transitions_export` |
| `TDA_*` (same five) | same | `<Reference>/bin/{orb,C,rho}.bin` (+`sigma.bin` if Reference_TBCF), `<Name>/bin/TDA_{X,E}.bin` | `densities/TDA_...`, `levels/TDA_...`, `TDA_summary.dat` |

`sE = string(round(E,digits=3))` (Fortran `jstr(round_digits(E,3))`). RPA externals error if `RPA_Orb = true`.
`RPA_external_JP("J±")` parses JP. All read binaries with `transpose(reshape(buf,N,N))` (row-major).

Core of the rewrite (`RPA_external.jl:173-330`): `RPA_transition_kernel` returns `Pairs, D_X, D_Y, f_N` on all pairs
`(a,b,t)` (t outer, a, b) of the multipole, keeping only pairs with a non-zero kernel:
* Linear: `D_X[(p,h,t),ph] += phase(J)*(n_h-n_p)`, `D_Y[(h,p,t),ph] += Phase_hp*(n_h-n_p)`.
* Exact: `D_X[(c,d,t)] += phase(J)*(δ(c,p)-Rho[p,c])*Rho[d,h]`, `D_Y[(c,d,t)] -= Phase_hp*(δ(c,h)-Rho[h,c])*Rho[d,p]`
  (+ Sigma part with `Phase_X`, `Phase_Y*phase((j_a-j_b)/2)` on `(a,b)` resp. `(b,a)` as in 4.8).
* `f_N[ph] = n_h - n_p`.
`Gamma = D_X*X[:,nu_list] + D_Y*Y[:,nu_list]` (complex). Multipoles (`RPA_external_multipoles`):
`Amp = Gamma/sqrt(2J+1)`, `pM += Op.p[a,b]*Amp`, ...; for E1 `isM += 0.5*(E1_C + E1_sC)[a,b]*Amp`
(base: `0.5*E1_C` only), `ivM += e_p*Op.p*Amp` / `-= e_n*Op.n*Amp`; other J: `±0.5*MME`.
Gauge `g = conj(pM)/|pM|` (1 if `pM == 0`) (base RPA: `sign(real(pM))`).
Radial densities (`:429-485`): `r^(2+Power)/(2J+1) * Σ_ab Re(g*Gamma_ab) * <a||Y_J||b> * R_a(r) R_b(r)` with
`<a||Y_J||b> = phase(J+(j_a-1)/2)*sqrt((j_a+1)(j_b+1))*cg(j_a,j_b,2J,1,-1,0)/sqrt(4π)` (0 if `l_a+l_b+J` odd)
(base had no angular factor). `Power = 0` (densities) or `(J==0 ? 2 : J)` (operator densities, J,P in E0..E3 only).
Grid `r = 1e-8 .. 2.5*1.2*A^(1/3)`, `2^13+1` points; `nu_p,n = 0.5*m*hw/hc^2`, `m_p=938.272013`, `m_n=939.565346`.
Averaged files only for (J,P) in {(0,+),(1,-),(2,+),(3,-)}, weights `phonon_weights!(M) = |M|/||M||` (common/functions.jl:860).
Currents (`:743-1052`): grid `r = 0.01..3.0*A^(1/3)` (`2^13+1`), `theta` 256 points on [0,2π]; map `xN x zN` on
`[-0.8 r_max, 0.8 r_max]^2` mapped to nearest (r,θ) node; pair maps `J_ab = Σ_{m_a} cg(j_a,j_b,2J,m_a,m_b,2M)*phase((j_b+m_b)/2)*J1b_*(...,-m_b,...)`
(spin part times `g/2`, `g_p = 5.585694713`, `g_n = -3.82608545`); level maps `Amp = -Re(g*Gamma_ab)/sqrt(2J+1)`;
scale conv `0.5*hc/m`, spin `hc/m`. TDA currents use `Amp = g*phase(J)*Gamma_ab` (no `-1/hat(J)`; inconsistency 7.6).
Analysis (`:647-741`): rows `a_p a_h T_ph Re X_ph Re Y_ph Re M_ph Re M_ph(-t)` for ph with `|(-1)^J x + y|^2 > Cut`,
`x = sqrt(f_N)X`; TDA version (`TDA_external.jl:759-907`) `a_p a_h T_ph X_ph M_ph M_ph(-t)`, criterion `|sqrt(f_N) X|^2 > Cut`.

These replace Fortran `phonon_transition_densities` / `phonon_transition_currents` (base algorithms); write them
anew around one kernel routine (Linear/Exact, with/without Y) shared by TDA and RPA.

---------------------------------------------------------------------------------------------------------------------

## 7. Port notes, Julia bugs and races

### 7.1 Port note 6 (OBDM race) - **fixed upstream / obsolete**
`RPA_OBDM_oY2` accumulates into per-JP buffers and sums serially. No deviation needed; delete note 6 (and its RRPA part
with HF_RRPA). Accumulate per JP then add in JP order to stay bit-close.

### 7.2 Port note 9 (spurious state in the OBDM) - **obsolete in its form, open question**
The 0.5 norm cut no longer exists; the spurious column is either a properly normalized solution of the 2x2 spurious
problem (`|E_s|, |n_s| > 1e-8`) or `X = Y = 0`. So the last-bit dependence is gone. But a spurious state with a **real
non-zero** `E_s` is now included in `RPA_OBDM_oY2`, `RPA_energy_QBA`, `RPA_sigma` (only E1 moments skip it). With a small
`E_s` its `|Y|^2 ~ p/(4 r E_s)` can be large (the "larger RPA radii" effect). Faithful port = include it.
A Fortran-only switch to exclude `nu = 1` of the orthogonalized 1- channel from OBDM/energy/Sigma would be a new
documented deviation (owner decision).

### 7.3 Other base port notes now fixed upstream
* "plot spectrum writes only N of N+1 levels" - fixed (`RPA_plot.dat`, `TDA_plot.dat` write N+1).
* "instability message always prints P = +" - fixed.
* "transitions_export uses OBDM_rN with the RPA density in the LHO basis" - fixed (`Rho_RPA = C'*Rho*C`).
* Notes 7, 8 and the two HF_RRPA "kept on purpose" items - obsolete with HF_RRPA removal.

### 7.4 New Julia bugs / oddities in dev (this group)
1. `TDA.jl:49` checks for a marker `.TDA_output` that is never written -> a second TDA run with the same `Name` refuses to
   overwrite (unless the dir is empty). Suggest: Fortran writes the marker (documented deviation) or reproduces.
2. `RPA.jl:34-38` (and `HF_ScRPA`) deletes `io/<Path>/<Name>` unconditionally, even if `Name == Reference` (e.g.
   `Name = "HF"` wipes the HF reference). No TDA-like guard. Also `if !isdir(io/<Path>) mkdir; return` is dead code.
3. `RPA_solver.jl:24,40`: `Reference_TBCF = true` with a missing `sigma.bin` is silently ignored (TDA errors).
4. `RPA_diagonalize.jl:345-360`: if `|E_s|` or `|n_s|` equals exactly `1e-8` neither branch fires and `E[1]` stays
   `undef` (garbage). Fortran: initialize `E(1)=0, X(:,1)=Y(:,1)=0` before (harmless deviation).
5. `RPA_spurious.jl:40-130` (Sigma variant of `RPA_spurious_ini`): undefined `J`, `P_pa`, `Orb_NN`, token `! =` -
   dead code, not ported.
6. TDA vs RPA current-map normalization differ: RPA `-Re(g Gamma)/hat(J)` (`RPA_external.jl:936`), TDA
   `g*(-1)^J*Gamma` (`TDA_external.jl:1283`) - port as is, flag upstream.
7. `RPA_OBDM_oY2` includes states with complex E (unstable) and zero-norm states are 0 anyway, while `RPA_energy_QBA`
   excludes complex states - comment claims "same set". Port as is.
8. Base `HF_RPA_amplitudes_export` `N_ph` shadowing - gone in dev.

### 7.5 Races
None found in the TDA/RPA files (all threaded loops write disjoint blocks/columns; `RPA_allocate_SE` changes BLAS
threads but is DEV-only).

---------------------------------------------------------------------------------------------------------------------

## 8. Fortran changes, routine by routine (with rough sizes)

### 8.1 New/rewritten modules in `src/rpa/` (suggested layout)
| module / routine | content | ~lines |
|---|---|---|
| `nj_ph_common.f90` | `ph_spurious_ini` (3.2), `ph_spurious_ortho` (3.3), `ph_metric` (`D = n_h-n_p`), shared `ph_allocate(Want_B)` (4.3 A/B/N), `ph_allocate_tbcf(Want_B,DEV)` (W1, W2, box diagrams; `pf/nf` under DEV) | 650 |
| `nj_tda.f90` | `TDA_diagonalize` (5.3), `TDA_transition_kernel`, `TDA_rM` (5.4) | 300 |
| `nj_rpa.f90` | `RPA_diagonalize` (4.5), `RPA_OBDM_oY2` (4.6), `RPA_energy` + `_QBA` + `_dress` + `_rho` (4.7), `RPA_rM_contract`/`RPA_rM` (4.8) | 650 |
| `nj_ph_transitions.f90` | shared `ph_rB` (4.9, iv fix), shared `ph_transitions_export(Kind)` (moments, TRK, 15 tables) | 450 (mostly from existing `HF_RPA_rB`, `HF_RPA_transitions_export`) |
| `nj_tda_export.f90` | `TDA_summary`, binary (row-major), spectra, plot, amplitudes with dominant component | 300 |
| `nj_rpa_export.f90` | `RPA_summary`, binary, spectra, plot (N+1), amplitudes (metric), `RPA_orbitals_export`, `RPA_natural_orbitals_export`, sigma export call | 400 |
| `nj_ph_external.f90` | `import_binary` (TDA/RPA, row-major), `external_import`, `transition_kernel` (Linear/Exact, TDA = no Y), multipoles, gauge, orbital grids, radial densities, operator densities, analysis, currents (TDA and RPA variants), `transitions_recalculate` | 1100 (replaces 875) |
| `nj_tda_driver.f90`, `nj_rpa_driver.f90` | `TDA`, `TDA_solver(_import)`, `RPA`, `RPA_solver` (4.1, 4.2, 5.1, 5.2) | 300 |
| `nj_dispatch_rpa.f90` | new actions (8.2) | 60 |
| DEV `RPA_allocate_SE` (optional) | | ~900 |

Total new/changed for this group ≈ 4300 Fortran lines (without `RPA_allocate_SE`).

### 8.2 `nj_dispatch_rpa.f90` actions
Remove `HF_RPA`, `HF_RPA_transition_densities`, `HF_RPA_transition_currents`, `HF_RRPA*`. Add:
`TDA`, `RPA`, `TDA_transition_densities`, `TDA_transition_operator_densities`, `TDA_transition_analysis`,
`TDA_transition_currents`, `TDA_transitions_recalculate`, and the same five with `RPA_` prefix
(`HF_ScRPA` belongs to the ScRPA group but lives in the same dispatcher).

### 8.3 Removals / things to keep
* Delete `nj_hf_rrpa.f90`, `nj_hf_rrpa_driver.f90`, `nj_hf_rrpa_export.f90` (~1700 lines), RRPA branch of the external code,
  `rrpa_params_t`, `&rrpa`, `obdm_iteration_pairs_t` (only used by RRPA - check with the core group).
* Delete `RPA_eigen_solve`, `RPA_eigen_solve_c`, `RPA_normalize`, `HF_RPA_energy`, `HF_RPA_OBDM`, base `HF_RPA_rM`,
  `HF_RPA_collectivity` equivalents.
* **Keep for EMPM** (Fortran-only): `HF_RPA_phonon_count` and `Spur_Ortho` (used by `src/empm/nj_empm_tda.f90:12`) and
  `rpa_summary_header`, `dimensions_block`, `spectrum_block_end` (used by `src/empm/nj_empm_export.f90:8`), together with
  the base `orbitals_ph_make`/`orbitals_one_phonon_make` they rely on. Move them into a small `nj_ph_legacy.f90`
  (or the EMPM tree) so that `nj_hf_rpa*.f90` can be deleted. EMPM also reads the HF output (`IO/<Path>/Bin/...`);
  its input paths must follow whatever the HF group decides.

### 8.4 Port details that matter for agreement with Julia
* Mirror lower triangles of A, B, N exactly as Julia (do not compute the upper triangle).
* Non-symmetric eigenproblem: call `eig_general` on `H` built as `(A*K + B*K)` row-scaled (same operation order).
* Generalized symmetric problem (TDA ortho): new `eig_sym_gen` = `dsygvd(itype=1,'V','U')` + stable ascending sort.
* Symmetric problems with `Symmetric(...)` use the upper triangle (`'U'`).
* Binary amplitude files are row-major; complex as (re,im) pairs.
* Initialize the spurious column (7.4-4).

---------------------------------------------------------------------------------------------------------------------

## 9. Interface changes needed in SHARED Fortran modules

1. `nj_params`: replace `rpa_params_t` by the full `RPA_Parameters` (2.2) - `Name`, `Reference` (`character(len=Str_Len)`),
   `Reference_TBCF`, `RPA_Orb`, `Ortho`, `Transitions` (`character(len=16)`, default `'Linear'`), `Tol`, `Tol_hol`, `Tol_par`,
   `Imax`, `OBDM_diagonal`, `OBDM_alpha`, `ScH`, `TBCF`, `TBCF_alpha`, `TBCF_beta`, `DEV`; add `tda_params_t`
   (`Name='TDA'`, `Reference='HF'`, `Reference_TBCF=.false.`, `Ortho=.true.`, `Tol_hol=0.7`, `Tol_par=0.1`) and
   `calculation_params_t%TDA`; remove `rrpa_params_t` / `%RRPA`.
2. `nj_params::io_path`: dev root is lowercase `io/`; add `io_dir(Params,Sub)` -> `'io/'//Path//'/'//Sub` for
   `<Reference>` and `<Name>` dirs (coordinate with HF group, which owns the reference layout).
3. `nj_input`: new `&tda` group (`Name, Reference, Reference_TBCF, Ortho, Tol_hol, Tol_par`) and extended `&rpa` group
   (all 17 fields). Namelist items share local variables inside `read_input`: `Transitions` already exists as a
   **logical** for `&ncsm` and `Name`/`Reference`/`Tol`/`Imax`/`Ortho`/`TBCF` collide with other groups -> read `&tda` and
   `&rpa` in their own internal subroutines (own locals) or rename the NCSM local. Remove `&rrpa`.
4. `nj_input::export_request_t` / new namelists: add `Cut` (default 0.1) and namelists `&transition_operator_densities`
   (`JP, nu_list, File_Name`, default `'Phonon_Transition_Operator_Densities'`) and `&transition_analysis`
   (`JP, nu_list, File_Name='Phonon_Transtion_Analysis', Cut`); currents `M` default `-1` meaning 0 (Julia).
5. `nj_types`: new `orb2b_1p1h_t` with `N(2,0:J_max)`, `Ind(2,0:J_max)%v(:)` of `ind2qp_t` (a,b,T) and dense
   `Dic(a,b,0:J_max,2)` (replaces the Julia UInt64 key; used by ScRPA `RPA_sigma`/`RPA_rho`). (`qporb2b_t` has the same shape
   without `Dic` and could be reused.)
6. `nj_orbitals`: new `orbitals_1p1h_make(Params,Orb,Tol_par,Tol_hol)` (3.1, defaults `1e-4`/`0.999` in Julia, always called
   with the solver's tolerances); `orbitals_import(Params,Path)` reading `orb.bin` (HF group may already add it);
   `orbitals_update(Orb,Rho)`.
7. `nj_o2b`: new pure accessor `O2b_temp(a,b,c,d,J,P,T,O_t,Orb_NN_t)` (T = -1 pp, 0 pn, 1 nn) on `o2b_temp_t`
   (`O2b_temp_initialize` exists); `O2b_initialize_empty` equivalent (unallocated `o2b_t`).
8. `nj_linalg`: new `eig_sym_gen(A,B,w,Z)` (dsygvd, itype 1, 'V', 'U', ascending stable sort).
9. `nj_tr1b` (f3f7a33 "toroidal fix", `common/Tr1b.jl:160-330`): in `Tr1b_E1_NLO_initialize` the convective gradient
   combinations (Bra `l_a+1` branch: `sqrt((l_a+1)/(2l_a+1))`, no `l_a>=1` guard, sign swaps of the `l_a-1`/`l_a+1` terms,
   same for Ket), spin normalization `16π` instead of `8π` under the sqrt, and factor `0.5` on `VC`, `TC`, `sTC` (p and n).
   Affects every E1 V/T/C strength (owned by the core group; required for parity of this group's outputs).
10. `nj_functions::J1b_convective` (f3f7a33, `common/functions.jl:504-622`): the `l+1` gradient term uses
    `sqrt((l+1)/(2l+1))` and is applied for all `l` (outside the `l > 0` block), for both brackets and both Bra/Ket terms.
11. `nj_obdm` (common `rho.jl`): `rho_chR2(Params,Orb,C,Rho)`, `rho_rN`, `rho_export(Params,Orb,Summary_File,Densities_File,Rho_LHO,C)`
    with explicit file names (dev `rho.jl`), `OBDM_natural_orbitals(Params,Orb,Rho) -> (Orb_NAT, C_NAT, Rho_NAT)` with
    Hungarian reordering (`common/orb.jl:552`); `nj_o1b`/`nj_o2b` import/export with `Orb`-based layout; `V2b_correlation_energy`
    (`common/V2b.jl:572`) for the density energy.
12. `nj_utils`/`nj_functions`: `phonon_weights` (`|M|/||M||`, unchanged if norm 0).
13. EMPM keeps using legacy `HF_RPA_phonon_count`, `Spur_Ortho`, `rpa_summary_header`, `dimensions_block`,
    `spectrum_block_end` (8.3) - must be moved before deleting `nj_hf_rpa*.f90`.
