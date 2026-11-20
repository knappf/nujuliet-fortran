# Dev port spec: HF-ScRPA (replaces HF-RRPA)

This spec covers Julia dev `f3f7a33` compared with base `2595e82`, which the Fortran port `dev-port` follows.

- **Dev scope:**
  - `src/RPA/RPA.jl` (`HF_ScRPA`)
  - `src/RPA/RPA_solver.jl` (`HF_ScRPA_solver`)
  - `src/RPA/RPA_rho.jl`
  - `src/RPA/RPA_sigma.jl`
  - `src/RPA/RPA_Hamiltonian.jl`
  - `src/RPA/RPA_energy.jl`
  - `src/RPA/RPA_orbitals.jl`
  - the helpers they call in `src/common`: `Hamiltonian.jl`, `orb.jl` (`OBDM_natural_orbitals`, `orbitals_update`, `orbitals_import`), `O1b.jl`, and `V2b.jl` (`V2b_correlation_energy`)
- **Base:** `Source/HF_RRPA/*` (13 files, 6.5k lines).
- **Fortran today:** `src/rpa/nj_hf_rrpa.f90` (1101 lines), `nj_hf_rrpa_driver.f90` (297 lines) and `nj_hf_rrpa_export.f90` (297 lines).

The following pieces are owned by the plain-RPA group: `RPA_allocate`, `RPA_allocate_TBCF`, `RPA_allocate_SE`, `RPA_diagonalize`, `RPA_rM`/`RPA_rB`, `RPA_export`, `orbitals_1p1h_make` and `Tr1b`. HF-ScRPA calls them unchanged, so this spec lists only their interfaces.

---------------------------------------------------------------------------------------------------

## 0. Executive summary: HF-ScRPA compared with the old HF-RRPA

| aspect | base HF-RRPA (`Source/HF_RRPA`, Fortran `nj_hf_rrpa*`) | dev HF-ScRPA (`RPA_solver.jl:94-418`) |
|---|---|---|
| entry / run action | `HF_RRPA(Params)`, output `IO/<Path>/RRPA/` | `HF_ScRPA(Params)`, output `io/<Path>/<RPA.Name>/` (dir is `rm -r`'d first) |
| params struct | `RRPA_Parameters(Ortho,Tol,Imax,ScV3N,ScOBH,dOBDM)` | shared `RPA_Parameters` (see section 6) |
| reference input | `orbitals_make` + `Bin/V2B_HF.bin`, `Bin/C_HF.bin`, HF occupations | `io/<Path>/<Reference>/bin/{orb.bin,C.bin,rho.bin,H1b.bin,H2b.bin}` |
| 1p-1h space | `orbitals_ph_make`/`one_phonon_make`/`HF_RRPA_phonon_count` (phonon lists) | `orbitals_1p1h_make(Orb_ref;Tol_par,Tol_hol)` (keyed `Orb2B_1p1h`), fixed during the run |
| RPA amplitude normalization | X^+X - Y^+Y = 1 (renormalized amplitudes) | X^+NX - Y^+NY = 1 with N = diag(D), D_ph = n_h - n_p (RPA group) |
| basis | HF -> "HF-RRPA" NO basis (`C_HF_RRPA`) | HF -> NAT (`C_HF_NAT`), the natural orbitals of the RPA OBDM, plus `C_LHO_NAT = C_LHO_HF*C_HF_NAT` |
| OBDM | `HF_RRPA_OBDM` (O(Y^4), scalar loops, quirks) | `RPA_OBDM_oY4` (same physics, BLAS `M = -1/2 X^T W X*`, `G = Y(1+M)Y^+`, each ordered element separately) |
| inner OBDM loop | no mixing, tol = `RRPA.Tol` (absolute mean abs diff), Imax = `RRPA.Imax`, D<0 aborts and returns `Rho_Ini` & identity, "cycled degenerate" exit | **linear mixing `OBDM_alpha`**, relative weighted-Frobenius residual < 1e-8 (hard-coded), max 150 (hard-coded), divergence guard, relative 2-cycle detection, **fallback = one explicit evaluation**, natural orbitals with block cleaning and degeneracy splitting |
| outer convergence | `eta` = mean abs OBDM change, weird exits with `E_old = 2` never updated | **relative change of E_QBA < `Tol` in 2 consecutive iterations** (and beta = 1 when TBCF) |
| mean field | `ScOBH`: h from Rho_LHO or Rho_HF; `ScV3N`: NO2B residual from Rho_LHO or rotation only | `ScH = true`: h = H1b(Rho_LHO) [+ V3N.Sigma], V = H2b(Rho_LHO); `ScH = false`: reference H1b.bin/H2b.bin rotated HF->NAT |
| 2-body correlations | none | optional **TBCF**: Sigma from X,Y (`RPA_sigma`), ramp `beta_i = tanh((i/TBCF_beta)^2)` (snaps to 1 above 0.95), linear mixing `TBCF_alpha`, enters A,B (`RPA_allocate_TBCF`), h (V3N.Sigma) and energy |
| stability | not tracked | `Stability` per iteration (complex eigenvalues), warnings, count written in the non-convergence summary |
| energy | `HF_RRPA_energy` (= -sum (2J+1) Re E abs2(Y)), extra prints | `RPA_energy` NamedTuple (E_QBA, E_sigma, E_rho, E_dress, E_dens), Density = true |
| non-convergence | always exports | **writes only `RPA_summary.dat` (termination summary) and returns without exporting anything else** |
| export | `HF_RRPA_export` (own layout) | generic `RPA_export` (same layout as plain RPA) |

Features that address the old **inner-loop stall** ("cycled degenerate solution", 208Pb RRPA):

1. **`OBDM_alpha` linear mixing (default 0.6).** A period-2 oscillation has a fixed-point-map eigenvalue near -1. Mixing turns it into 1-2*alpha = -0.2, so the 2-cycle is damped.
2. **`OBDM_natural_orbitals`.**
   - It zeroes the elements outside the (l,j) blocks.
   - It adds `1e-10*a` to the diagonal, a deterministic splitting that makes the eigenvectors of exactly degenerate occupations unique (for example the many zero particle occupations). Without it the Hungarian matching jumps between equivalent bases from one iteration to the next.
   - The base code diagonalized the raw matrix with `eigen`, with no symmetrization guard and no splitting.
3. **No bail-out on negative D.** Base `HF_RRPA_OBDM_occupation_check` returned `Rho_Ini` and the identity whenever one D_ph < 0, so the outer loop repeated the same state. In dev, pairs with D < 1e-8 are simply skipped.
4. **Bounded inner loop.** It stops at 150 iterations, on divergence or on a relative 2-cycle (`abs(delta - delta_prev) < 1e-3*delta`). It then falls back to one explicit evaluation, so it never hangs.
5. **Outer convergence on E_QBA**, not on the OBDM.
6. **Cost.** M and G are computed with BLAS once per (J,P) block. Each OBDM element then costs O(n_h * n_J) look-ups instead of O(n_h * n_J * N_ph^2) (base inner/mirror loops). This is the dominant runtime of the base on 208Pb, and the Fortran port has the same loop structure.

---------------------------------------------------------------------------------------------------

## 1. File-by-file

### 1.1 `src/RPA/RPA.jl` -> `HF_ScRPA(Params)` (lines 72-141)

- **Counterpart:** base `HF_RRPA/HF_RRPA.jl`. Fortran: `nj_hf_rrpa_driver.f90::HF_RRPA`.
- **Changes:**
  - It prints a header and checks that `io/<Path>/<RPA.Reference>` exists. If it does not, it prints an error and returns.
  - `IO_RPA = "io/<Path>/<RPA.Name>"`: `rm(IO_RPA, recursive=true)` if it exists, then mkdir `IO_RPA` with the subdirectories `bin`, `densities`, `transitions/{E0,E1,E2,E3}`, `spectra` and `amplitudes`. This is the same set as `RPA()`, lines 1-70.
  - It calls `HF_ScRPA_solver(Params)`.
  - The `isdir("io/<Path>")` branch (lines 97-100) is unreachable, because the reference directory already implies that it exists.
- **Julia hazard:** if `Name == Reference`, the reference solution is deleted before it is read (line 105). The port should refuse to run with `fatal` (see section 7).
- **Fortran:** new `subroutine HF_ScRPA(Params)` in a new module `src/rpa/nj_scrpa.f90`. It reuses the directory helper of the RPA group, because `RPA` and `HF_ScRPA` make identical trees (`rm -r` needs `execute_command_line('rm -rf ...')` or the port's existing utility). About 40 lines.

### 1.2 `src/RPA/RPA_solver.jl` -> `HF_ScRPA_solver(Params)` (lines 94-418)

- **Counterpart:** base `HF_RRPA_Solver.jl::HF_RRPA_solver` + `HF_RRPA_solve`. Fortran: `nj_hf_rrpa_driver.f90::HF_RRPA_solver/HF_RRPA_solve` (to be deleted).

**Setup** (lines 96-234)

- `IO = io/<Path>/<Reference>`, `epsilon = RPA.Tol`, `Imax = RPA.Imax`, `alpha_TBCF = RPA.TBCF_alpha`, `beta_TBCF = RPA.TBCF_beta`, `TBCF = RPA.TBCF`.
- Local function `HF_ScRPA_solver_beta(i) = tanh((i/beta_TBCF)^2)`, set to 1.0 when it is > 0.95. With the default beta_TBCF = 5, i = 1..6 give 0.0400, 0.159, 0.345, 0.556, 0.762, 0.894, and i >= 7 gives 1.
- Data read and built, in this order:
  1. `Orb_ref = orbitals_import(IO/bin/orb.bin)`. Binary: per orbital `Int64 a,n,l,j` and `Float64 pO,nO`, a_max records.
  2. `T = T1b(Params,Orb_ref)`.
  3. `V_NN_bare, Orb_NN = V2b_read(Params)`. Dev's `V2b_read` takes only Params.
  4. `V_NNN_bare, Orb_NNN = V3b_no2b_read(Params)`. Both are read even when `ScH = false`, where they are unused.
  5. `C_LHO_HF = O1b_import(IO/bin/C.bin)`, `Rho_ref = O1b_import(IO/bin/rho.bin)`, `h_N = O1b_import(IO/bin/H1b.bin)`, `V_NN = O2b_import(IO/bin/H2b.bin)`. No Make_Orb_NN; `Orb_NN` from `V2b_read` is used everywhere.
  6. `Sigma = empty`.
  7. `Orb_1p1h = orbitals_1p1h_make(Orb_ref;Tol_par,Tol_hol)`.
  8. `Orb_ph = RPA_OBDM_orbitals_make(Orb_ref)` (section 1.7).
  9. `TrOp_LHO = Tr1b_initialize(Orb_ref,1.0)` and `TrOp = Tr1b_transformation(Orb_ref,C_LHO_HF,TrOp_LHO)`.
- Initial RPA:
  - `RPA_allocate(...,h_N,V_NN,Rho_ref,Sigma,false)`, then `RPA_diagonalize(...,A,B,N,TrOp)`.
  - If unstable, print `"\n\tWarning! The initial RPA solution is unstable ..."`.
  - `E_QBA_old = RPA_energy_QBA(Orb_ref,...)`, `E_QBA_list = [E_QBA_old]`, `Stability_list = Bool[]`.
- State before the loop:
  - `Orb = Orb_ref`, `Rho = copy(Rho_ref)`, `C_HF_NAT = I`, `C_LHO_NAT = copy(C_LHO_HF)`.
  - `Iteration = 1`, `Converged = false`, `N_converged = 0`, `beta = beta(1)`.
  - `Sigma = Sigma_old = empty`.

**Loop** `while Iteration <= Imax` (lines 266-340)

```
C_HF_NAT_old = copy(C_HF_NAT)
1) if TBCF:  Sigma = RPA_sigma(Orb,Orb_1p1h,Orb_NN,X,Y; beta=beta, Report=false)        # new object, current NAT basis
             if Iteration > 1:  RPA_sigma_mix!(Sigma,Sigma_old,alpha_TBCF)               # Sigma = a*Sigma + (1-a)*Sigma_old
2) Rho, C_LHO_NAT, Orb = HF_ScRPA_OBDM_solve(Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT)     # C_HF_NAT updated IN PLACE
   if TBCF:  C_step = (C_HF_NAT_old' * C_HF_NAT)  (p & n)
             Sigma = O2b_transformation(Orb,Orb_NN,Sigma,C_step)                         # in place
             Sigma_old = deepcopy(Sigma)
4) h_N, V_NN = RPA_Hamiltonian(Orb,Orb_NN,Orb_NNN,C_LHO_NAT,Rho,Sigma,T,V_NN_bare,V_NNN_bare)
   TrOp = Tr1b_transformation(Orb,C_LHO_NAT,TrOp_LHO)
5) A,B,N = RPA_allocate(Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Rho,Sigma,TBCF)
   E,X,Y,Stability = RPA_diagonalize(Orb,Orb_1p1h,A,B,N,TrOp)
6) E_QBA = RPA_energy_QBA(Orb,Orb_1p1h,E,X,Y)
   dE_QBA = |E_QBA - E_QBA_old| / |E_QBA|   (absolute if |E_QBA| <= 1e-12);  E_QBA_old = E_QBA
   push E_QBA_list, Stability_list
   @printf("\tHF-ScRPA iteration:   Iteration number = %4d   E0 = %12.9f MeV   dE0_rel = %12.9f\n", ...)
   N_converged = (dE_QBA < Tol && (!TBCF || beta == 1.0)) ? N_converged+1 : 0
   if N_converged >= 2 && (!TBCF || beta == 1.0): Converged = true; break
   beta = beta(Iteration+1);  Iteration += 1
```

- Unlike `RPA_solver`, the dev OBDM step does **not** use `RPA_OBDM_oY2`. `Delta_E` and `Delta_E_old` are computed but never used.
- `N_iterations = Converged ? Iteration : Iteration-1`, then print `"\nHF-ScRPA iteration finished after %4d iterations ...   E0 = %12.9f MeV   dE0_rel = %12.9f\n"`.
- If not converged:
  - print `"\nHF-ScRPA iteration has NOT converged ... Run terminated ..."`;
  - call `ScRPA_terminate_summary(N_iterations,E_QBA_list,Stability_list,beta(N_iterations))`;
  - **return**, with no other export.

**Final step** (lines 359-415)

- If TBCF, recompute `h_N, V_NN = RPA_Hamiltonian(...)`. This gives the same Rho and Sigma as step 4 of the last iteration, so it is redundant but harmless.
- `A,B,N = RPA_allocate(...,TBCF)`, `Sigma_AB = TBCF ? Sigma : nothing`, `E,X,Y,Stability = RPA_diagonalize`.
- If unstable, print `"\nWarning! The final HF-ScRPA solution is unstable ..."`.
- `Sigma = RPA_sigma(Orb,Orb_1p1h,Orb_NN,X,Y)`. This runs **always** (beta = 1, Report = true), even for TBCF = false, because `E_dens` needs it.
- `chR2 = rho_chR2(Orb,C_LHO_NAT,Rho)`, `TrOp = Tr1b_initialize(Orb,chR2)`, then `Tr1b_transformation(Orb,C_LHO_NAT,TrOp)`.
- `Orb = orbitals_update(Orb,Rho)`. This is a no-op, because Orb already carries `diag(Rho)` (it mutates in place).
- `E_corr = RPA_energy(Orb,Orb_1p1h,Orb_NN,E,X,Y,Rho, O1b_transformation(Rho_ref,C_HF_NAT), h_N,Sigma,V_NN; Sigma_AB, Density=true)`.
- `rM = RPA_rM(Orb,Orb_1p1h,Orb_NN,TrOp,Rho,Sigma,X,Y)`, `rB = RPA_rB(Orb_1p1h,rM,real.(E))`.
- `Rho_LHO = C_LHO_NAT * Rho * C_LHO_NAT'` (written as `O1b_transformation(Rho, transpose(C_LHO_NAT))`).
- `RPA_export(Orb,Orb_1p1h,Orb_NN,Stability,E_corr,E,X,Y,rB,C_LHO_NAT,Rho_LHO,Sigma,h_N)`. `h_N` is in the NAT basis. `sigma.bin` is written only if `TBCF = true`.

**`ScRPA_terminate_summary`** (lines 117-174) writes `io/<Path>/<Name>/RPA_summary.dat` (mode `"w"`):

```
Spherical self-consistent Random-Phase Approximation (HF-ScRPA) review:
<blank>!!! The HF-ScRPA iteration has NOT converged ... no RPA results were exported !!!
<blank>Nuclide data:    A = <A>, Z = <Z>
<blank>Calculation data:    hw = <jstr hw> , N_max = <N>, s.p. J-basis size = <a_max>, s.p. M-basis size = <(N+1)(N+2)(N+3)/6>
\t\t\tCenter-of-Mass System correction is set to:   CMS = ''<CMS>''
<blank>Orthogonalization of spurious Goldstone levels (1-) is [NOT ]included ...
<blank>HF-ScRPA parameters:
\tReference = ''<Ref>'',   Reference_TBCF = <bool>,   ScH = <bool>
\tTol = <jstr>,   Imax = <int>,   Tol_hol = <jstr>,   Tol_par = <jstr>
\tOBDM_alpha = <jstr>,   OBDM_diagonal = <bool>
\tTBCF = <bool>,   TBCF_alpha = <jstr>,   TBCF_beta = <jstr>
<blank>Convergence information:
\tIterations performed = <N_it> (Imax = <Imax>)
 TBCF: \tCriterion ... relative change of E_QBA below Tol in two consecutive iterations, with beta = 1 ...
       \tLast beta (ramp of Sigma) = <jstr beta>
 else: \tCriterion ... relative change of E_QBA below Tol in two consecutive iterations ...
\tUnstable iterations (complex RPA eigenvalues) = <count false> of <length>
 if N_E>=2: \tLast E0 = <jstr> MeV,   last change = <jstr |dE|> MeV (relative <jstr |dE|/max(|E|,1e-12)>)
<blank>Last iterations ... (Iteration, E_QBA [MeV], relative change, stable)
 for i = max(2,N_E-9):N_E : "<i-1>\t<jstr E_i>\t<jstr |E_i-E_{i-1}|/max(|E_i|,1e-12)>\t<true|false>"
```

- Julia `string(Bool)` gives `true`/`false`. Floats use `jstr`.
- `N_E = length(E_QBA_list) = N_iterations + 1`.

**Fortran:** new module `nj_scrpa`:
- `HF_ScRPA_solver` (~230 lines)
- internal `scrpa_beta` (5 lines)
- `ScRPA_terminate_summary` (~60 lines)

### 1.3 `src/RPA/RPA_rho.jl`

Counterpart: base `HF_RRPA_OBDM.jl` (and `HF_RPA_OBDM` of HF_RPA for oY2). Fortran: `nj_hf_rrpa.f90::HF_RRPA_OBDM_iteration` with internals `obdm`, `pair_sum`, `reorder`, `occupation_check` and `transform_amplitudes`.

#### (a) `RPA_OBDM_oY2(Params,Orb,Orb_1p1h,Y)` (lines 1-73)

- Used by plain `RPA_solver`, not by ScRPA. Listed for completeness; the plain-RPA group owns it.
- `Rho = diag(Orb.pO/nO)`. For each (J,P) block, a per-block accumulator (no race):

  ```
  Sum = sum_nu |Y_ph,nu|^2 (2J+1)
  ```

  where `ph = (a_p,a_h,T)`.
  - `T = -1`: `pRho[h,h] -= (n_h-n_p)^2/(j_h+1) * Sum` and `pRho[p,p] += (n_h-n_p)^2/(j_p+1) * Sum`.
  - `T = +1`: the same with n.
- All solutions are included; there is **no norm cut**. This makes port note 9 (`0.5-1e-8` cut) obsolete for dev code paths, and the spurious state is now handled in `RPA_diagonalize`.

#### (b) `RPA_OBDM_sigma` (lines 91-266)

**Dead code in dev** (no caller). Do not port, or port later. The formula is documented in the comment block at lines 75-90.

#### (c) `RPA_OBDM_oY4(Params,Orb,Orb_ph,Orb_1p1h,Rho,X,Y,pM,nM)` (lines 283-478)

This replaces the base `HF_RRPA_OBDM`.

- **Flags.** `p_par[a]` is true if `a in Orb_ph.p.p`; `p_hol`, `n_par` and `n_hol` are defined the same way. Old occupations: `Rho_old = Rho` (its diagonal is used).
- **Per (J,P) block, largest blocks first**, with single-threaded BLAS inside (lines 343-368):

  ```
  pW[i_qg] = D_qg^2 if T_qg==-1 && p_par[q] && p_hol[g] && D_qg=pRho_old[g,g]-pRho_old[q,q] >= 1e-8, else 0   (nW analog)
  pM[P,J+1] = -0.5 * transpose(X_JP) * (pW .* conj(X_JP))      # M_nu,mu = -1/2 sum_qg X_qg,nu W_qg X*_qg,mu
  pG[P,J+1] = Y_JP * ((I + pM) * adjoint(Y_JP))                # G = Y (1+M) Y^+   (N_ph x N_ph, complex)
  ```

  `pM`/`nM` are caller-owned storage that is overwritten. `pG`/`nG` are local.
- **Task list** (lines 374-392), in the order `T=-1` then `T=+1`. For each T: all ordered particle pairs (a,b) in `Orb_ph.p.<t>` with equal l,j (outer a, inner b), then all ordered hole pairs.
- **Each task writes only `dRho[i]`** (lines 398-463):

  ```
  Type :p  (a,b particles): Sum = sum_{h in Hol} [D_ah>=1e-8 && D_bh>=1e-8]  sum_{J=|j_a-j_h|/2}^{(j_a+j_h)/2}
               i_ah = idx(a,h,J,P=parity(l_a+l_h),T), i_bh = idx(b,h,...); skip if 0
               Sum += (2J+1) * D_ah * D_bh * Re G_T[P,J+1][i_ah,i_bh]          (D_xh = Rho_T[h,h]-Rho_T[x,x])
  Type :h  (a,b holes):     Sum -= same with sum over particles p, D_pa = Rho_T[a,a]-Rho_T[p,p], i_pa=idx(p,a), i_pb=idx(p,b)
  dRho[i] = Sum / (j_a + 1)
  ```

- Result: `Rho_new = Rho_old + dRho` (copy, then `+=` at `[a,b]`).
- **Differences from base `HF_RRPA_OBDM`** (and from the Fortran `obdm`/`pair_sum`):
  1. With the dev normalization (X_dev = X_base/sqrt(D)), `D^2 |X_dev|^2 = D |X_base|^2` and `D_a D_b Y_dev Y_dev* = sqrt(D_a D_b) Y_base Y_base*`. **The formulas are identical** once the normalization of the RPA group is adopted. Port the dev form literally, because the amplitudes come from the new `RPA_diagonalize`.
  2. Base `M[nu,mu] = -1/2 sum D X_mu X*_nu` is the transpose (complex conjugate) of dev; the two are equal for real amplitudes.
  3. **Base quirk 1 is fixed.** Base computed only a >= b and, for a > b, added the "mirror" inner term **and** wrote the sum to both [a,b] and [b,a], so the off-diagonal inner (M) part was counted twice. Dev computes each ordered pair once.
  4. **Base quirk 2 is fixed.** The neutron-hole mirror term used `Y[pg,mu]` (Fortran `Quirk` flag). It is gone.
  5. **Base quirk 3 is gone.** Base applied `0.5(R+R')` symmetrization and a `1e-8` cleaning of Rho and C. Dev has neither (Hermitian G gives symmetric dRho).
  6. Pairs with D < 1e-8 are skipped (base used `sqrt(D D)`, which gives NaN for a negative product).
- Port note "HF_RRPA_OBDM, neutron-hole branch ... kept on purpose" is therefore **obsolete**. Port note 6 (serial accumulation for the race) is obsolete for this routine, because dev has no race here.

#### (d) `HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho_ini,C_HF_NAT)` (lines 528-715)

This replaces `HF_RRPA_OBDM_iteration`.

- **Constants and inputs.**
  - `epsilon = 1e-8` and `Iteration_max = 150` are **hard-coded**; they do not come from Params.
  - `Diagonal = RPA.OBDM_diagonal`, `alpha = RPA.OBDM_alpha`.
  - It re-imports `Rho_HF = rho.bin` and `C_LHO_HF = C.bin` from the reference directory on every call.
  - `Rho_ref = C_HF_NAT' Rho_HF C_HF_NAT` (reference OBDM in the current NAT basis).
- **Local functions.**
  - `OBDM_update(Rho_ref,Rho_old,X,Y) = Rho_ref + (oY4(Rho_old,X,Y) - Rho_old)`, which is `Rho_ref + dRho[Rho_old]`.
  - `OBDM_natural(Rho)`: `C = OBDM_natural_orbitals(Params,Orb,Rho)[2]`, `Rho_NAT = diag(diag(C' Rho C))`. It recomputes from the **unsplit** Rho and drops the off-diagonal part.
  - `OBDM_norm(pO,nO) = sqrt( sum_a sum_b (j_a+1) (pO[a,b]^2 + nO[a,b]^2) )`.
- **Iteration** (lines 600-680):

  ```
  X_it,Y_it = X,Y ; Rho_old = copy(Rho_ini) ; Rho_ref_it = copy(Rho_ref) ; C_tot = I ; Norm_ref = OBDM_norm(Rho_ref)
  delta = delta_odd = delta_even = 1.0 ; print "\n\tStarting RPA O(Y^4) OBDM iteration ..."
  while Iteration < 150
    Rho_new = OBDM_update(Rho_ref_it,Rho_old,X_it,Y_it)
    if Diagonal: Rho_new = diag(Rho_new)
    Rho_new = alpha*Rho_new + (1-alpha)*Rho_old
    if any non-finite: print "\n\tOBDM iteration stopped - diverging solution (Inf or NaN) ..."; break
    C_step, Rho_new = OBDM_natural(Rho_new)
    C_tot = C_tot*C_step ; Rho_old = C_step' Rho_old C_step ; Rho_ref_it = C_step' Rho_ref_it C_step
    X_it,Y_it = RPA_OBDM_transform_amplitudes(Orb_1p1h,X_it,Y_it,C_step)
    R = (Rho_new - Rho_old)/alpha ; delta_abs = OBDM_norm(R) ; Norm_cor = OBDM_norm(Rho_new - Rho_ref_it)
    delta = Norm_cor > 1e-12 ? delta_abs/Norm_cor : delta_abs
    Iteration += 1 ; @printf("\t\tOBDM iteration:   Iteration number =  %4d   delta = %12.9f\n")
    Rho_old = Rho_new
    if !isfinite(delta) || delta_abs > Norm_ref: print "\n\tOBDM iteration stopped - diverging solution ..."; break
    if delta < 1e-8: Converged = true; break
    if Iteration > 100: prev = iseven ? delta_even : delta_odd
        if |delta - prev| < 1e-3*delta: print "\n\tOBDM iteration stopped - cycled degenerate solution found ..."; break
        store delta into delta_even/odd
  ```

- **After the loop.**
  - Converged: `Rho = Rho_old` (diagonal, newest basis).
  - Not converged:
    1. print `"\n\tRPA O(Y^4) OBDM iteration has not converged - using a single explicit evaluation of the OBDM ..."`;
    2. `Rho_expl = OBDM_update(Rho_ref,Rho_ini,X,Y)` with the **input** amplitudes and the input basis, and **without mixing**;
    3. optionally take the diagonal;
    4. `C_tot, Rho = OBDM_natural(Rho_expl)`.
  - Then:
    - `Orb_NAT[a] = Orb1B(a,n,l,j,Rho.p[a,a],Rho.n[a,a])` (new vector; n,l,j from input Orb);
    - `C_HF_NAT .= C_HF_NAT*C_tot` (in place);
    - `C_LHO_NAT = C_LHO_HF*C_HF_NAT`;
    - return `(Rho, C_LHO_NAT, Orb_NAT)`.
- Orb, Orb_ph and Orb_1p1h stay fixed inside the loop; only Rho, X and Y rotate.
- **Removed base features:**
  - `occupation_check` and its "Singular solution" return of (Rho_Ini, I);
  - the "lower precision" exits at 150/200/250 iterations;
  - the `RRPA.Imax`/`RRPA.Tol` use;
  - the `dSum/(2 a_max^2)` measure.

#### (e) `RPA_OBDM_transform_amplitudes(Orb_1p1h,X,Y,C)` (lines 717-776)

Per (P,J) block, largest first, threads with BLAS = 1:

```
U[i_ph,i_qg] = C_T[a_q,a_p] * C_T[a_g,a_h]   if T_qg == T_ph (else 0)     ;  X_new = U*X ; Y_new = U*Y
```

- This is the same mapping as base `HF_RRPA_OBDM_transform_amplitudes` (`X_ph = sum_qg C_qp C_gh X_qg`) but done as one dense matmul per block over the fixed `Orb_1p1h` list.
- It returns new arrays; the inputs are not mutated. The Fortran version currently transforms in place; it must return copies, or the solver must keep X,Y.
- **Fortran:** build U, then `zgemm` (or a real dgemm on U with complex X: split real and imaginary parts, or promote U to complex).

**Fortran plan for `RPA_rho.jl`** (new module `src/rpa/nj_rpa_rho.f90`):
- `RPA_OBDM_oY4` ~150 lines, using `zgemm`. M and G use OpenMP over blocks with single-threaded BLAS, or a serial block loop with threaded BLAS; results must be thread-count independent.
- `HF_ScRPA_OBDM_solve` ~150 lines.
- `RPA_OBDM_transform_amplitudes` ~50 lines.
- `RPA_OBDM_oY2` ~50 lines (if the RPA group does not own it).
- `Orb_ph` builder ~40 lines.
- Delete `HF_RRPA_OBDM_iteration` and its internals (~450 lines).

### 1.4 `src/RPA/RPA_sigma.jl` (new; no base counterpart in `Source/` except `Source/Temp`)

- **Fortran:** none. New module `src/rpa/nj_rpa_sigma.f90`, shared with plain RPA (`TBCF = true`): `RPA_sigma` ~250 lines, three kernels ~120 lines, `RPA_sigma_mix` ~20 lines.

**`RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y; beta=1.0, Report=true) -> O2B`** (lines 1-295)

- **Storage and classification.**
  - `Sigma = O2b_initialize(Params,Orb)`: same packed layout and N2max channels as V_NN, in the convention `<ab|Sigma|cd>_A` (unnormalized antisymmetrized, no 1/sqrt2).
  - `particle(O) = O <= Tol_par` and `hole(O) = O >= Tol_hol`. These are **inclusive**, whereas `orbitals_1p1h_make` and `RPA_OBDM_orbitals_make` use strict comparisons. Occupations come from the **current** Orb (NAT occupations), while `Orb_1p1h` is fixed from the reference.
- **Loop order.** For each (J,P) in `JP_initialize(2Nmax+1)`, the work list BraKet comes in the order pp (`Bra=1..N_T1`, `Ket<=Bra`, `T=-1`), pn (`N_T0`, `T=0`), nn (`T=1`). Index: `Ind = Bra + (Ket-1)*N - Ket(Ket-1)/2`.
- **Parity gating.** `Q = (l_a+l_c)%2+1`, `Qx = (l_a+l_d)%2+1`, and the block is processed only if `Q == (l_b+l_d)%2+1`.
- **pp and nn elements:** the sum of these terms, times `beta`:
  - 2p-2h (a,b particles; c,d holes):

    ```
    0.25*[ S(a,b,c,d;Q) - Ph_cd S(a,b,d,c;Qx) - Ph_ab S(b,a,c,d;Qx) + Ph_ab Ph_cd S(b,a,d,c;Q) ]
    ```

    with `S = pphh` and `Ph_xy = (-1)^(J+(j_x+j_y)/2)`.
  - 2h-2p (a,b holes; c,d particles): the same with `pphh(c,d,a,b)`, `pphh(d,c,a,b)`, `pphh(c,d,b,a)`, `pphh(d,c,b,a)`.
  - 1p-1h, each as `K(a,b,c,d;Q) - (-1)^(J+(j_c+j_d)/2) K(a,b,d,c;Qx)`:
    - (p,h,p,h) and (p,h,h,p) with `K = phhp` (Y Y*);
    - (h,p,h,p) and (h,p,p,h) with `K = hpph` (X X*).
- **pn elements:** `pO` for a,c and `nO` for b,d; no exchange; times `beta`:
  - `pphh(a,b,c,d;Q)` [pp/hh];
  - `pphh(c,d,a,b;Q)` [hh/pp];
  - `phhp(a,b,c,d;Q, Y)` [a p, b h, c h, d p];
  - `hpph(a,b,c,d;Q, **Y**)` [a h, b p, c p, d h]. **See bug B1:** the pp/nn branches pass X here.

**Kernels** (`I` runs over `max(|j_a-j_c|,|j_b-j_d|)/2 .. min(j_a+j_c,j_b+j_d)/2` in all three):

```
pphh: T_ac=(t_a+t_c)/2, T_bd=(t_b+t_d)/2 ; i_ac=idx(a,c,I,P,T_ac), i_bd=idx(b,d,I,P,T_bd)
      Amp = (2I+1) (-1)^(J+I+(j_b+j_c)/2) {j_a j_b J; j_d j_c I} |(O_c-O_a)(O_d-O_b)| ; sigma += Amp sum_nu Re(Y[i_ac,nu] conj(X[i_bd,nu]))
phhp: i_ac=idx(a,c), i_db=idx(d,b) ; Amp = (2I+1)(-1)^(J+(j_c+j_d)/2) 6j |(O_c-O_a)(O_b-O_d)| ; sum Re(Y_ac conj(Y_db))
hpph: i_ca=idx(c,a), i_bd=idx(b,d) ; Amp = (2I+1)(-1)^(J+(j_a+j_b)/2) 6j |(O_a-O_c)(O_b-O_d)| ; sum Re(X_ca conj(X_bd))   [arg named X]
all three, after the I loop:  if a==c && b==d: sigma -= O_a O_b ;  if (a==d && t_a==t_d) && (b==c && t_b==t_c): sigma -= (-1)^(J+(j_a+j_b)/2) O_a
```

- `w6j(j1,j2,j3,j4,j5,j6)` takes doubled arguments: `{j_a/2 j_b/2 J; j_d/2 j_c/2 I}`, called as `w6j(j_a,j_b,2J,j_d,j_c,2I)`. Use `f6j` from `nj_wigner`.
- `idx = orb2b_1p1h_index` (key `P | (T+1)<<2 | J<<4 | a<<10 | b<<18`); it returns 0 if absent.
- `N_nu = Orb_1p1h.N[P,I+1]`.
- Threads over BraKet tuples; each writes its own element, so there is no race.
- **Caveat:** J runs to `2Nmax+1`, so the code assumes `N2max = 2*Nmax`. With `N2max < 2Nmax` the `Orb_NN.N` indexing goes out of bounds, as in port note 3 for `O2b_export`. The port should loop only over `0..N2max+1`.

**`RPA_sigma_mix!(Sigma,Sigma_old,gamma)`** (lines 420-434): if `|1-gamma| < 1e-12` it is a no-op. Otherwise `S = gamma*S + (1-gamma)*S_old` for each assigned block of pp, pn and nn.

### 1.5 `src/RPA/RPA_Hamiltonian.jl` (new; replaces `HF_RRPA_Hamiltonian.jl`)

- **Counterpart:** base `HF_RRPA_h1b`, `HF_RRPA_V2b`, `HF_RRPA_V2b_transform`. Fortran: `nj_hf_rrpa.f90::HF_RRPA_h1b/HF_RRPA_V2b/HF_RRPA_V2b_transform`.
- **`RPA_Hamiltonian(Params,Orb,Orb_NN,Orb_NNN,C,Rho,Sigma,T,V_NN_bare,V_NNN_bare) -> (h_N, V_NN)`.** `C` is `C_LHO_NAT`; Rho and Sigma are in NAT.
  - `ScH == true`, with `Rho_L = C Rho C^T` (`O1b_transformation(Rho,transpose(C))`):
    - `TBCF == true`:
      1. `Sigma_L = O2b_transformation(Sigma, C^T)` (**in place**: the caller's Sigma is changed);
      2. `h_L, V_L = H1b_H2b(Orb,Orb_NN,Orb_NNN,T,V_NN_bare,V_NNN_bare,Rho_L,Sigma_L)`;
      3. `h_N = C^T h_L C`, `V_NN = O2b_transformation(V_L, C)`;
      4. Sigma is transformed back with C (in place). The round trip leaves the caller's Sigma unchanged up to rounding. **The port must reproduce this or transform a copy**; a copy is cleaner, and its effect on output is about 1e-15.
    - `TBCF == false`: `h_L = H1b(..., Rho_L)`, `V_L = H2b(..., Rho_L)`, then `h_N = C^T h_L C` and `V_NN = O2b_transformation(V_L,C)`.
  - `ScH == false`: `RPA_Hamiltonian_HF_reference`.
    1. Fresh imports of `H1b.bin` and `H2b.bin` (HF basis) and `C.bin`.
    2. `C_HF_NAT' = C_LHO_HF^T C` (p,n).
    3. `h_N = C_HF_NAT'^T H C_HF_NAT'`, `V_NN = O2b_transformation(H_NN, C_HF_NAT')`.
    4. Sigma is ignored for h,V; TBCF then enters only through A, B and the energy.
- **Changes from base:**
  - `ScOBH` and `ScV3N` are merged into one `ScH`.
  - `ScH = false` now also keeps the HF residual interaction (base `ScV3N = false` rotated `V2B_HF.bin`; port note 7).
  - Base `ScOBH = false` evaluated h from Rho_HF with the bare interaction. Dev uses the stored HF H1b, which is the same for a converged HF.
  - The base symmetrizations `0.5(H+H')` and `0.5(Rho_LHO+Rho_LHO')` are gone.
  - stdout redirection (port note 8) is replaced by `Report = false`, but H1b still prints `"\nEvaluating 1-body mean-field Hamiltonian h_N ..."` and `"\tCompleted ..."`, and `RPA_Hamiltonian` prints `"\nConstructing Hamiltonian for RPA calculation ..."` / `"\tFinished construction of RPA Hamiltonian ..."`.
- **Fortran:** `RPA_Hamiltonian` + `RPA_Hamiltonian_HF_reference` ~80 lines in `nj_scrpa.f90`, calling the core H1b/H2b/H1b_H2b (section 2).

### 1.6 `src/RPA/RPA_energy.jl` (new; replaces `HF_RRPA_Energy.jl`)

Counterpart: base `HF_RRPA_energy`, `HF_RRPA_energy_density`, `HF_RRPA_energy_bosonic`. Fortran: `nj_hf_rrpa.f90::HF_RRPA_energy/_bosonic/bosonic_sum/_density`. Shared with plain RPA, so it may be owned jointly with the RPA group.

**`RPA_energy(...; Sigma_AB=nothing, Density=RPA.TBCF)`** returns the NamedTuple `(E_QBA, E_sigma, E_rho, E_dress, E_dens)`, with NaN for terms that were not evaluated (`E_dens = E_sigma + E_rho + E_dress`).

- ScRPA calls it with `Density = true`; the dressing is evaluated only if `Sigma_AB !== nothing`.
- Dressing:
  - `Zero` = zero N_ph x N_ph blocks;
  - `dA,dB = RPA_allocate_TBCF(Orb,Orb_1p1h,Orb_NN,V_NN,Sigma_AB,Zero,copy(Zero))`;
  - if `RPA.DEV`: `dA = RPA_allocate_SE(...,dA,Mode="OnShell")`;
  - `E_dress = RPA_energy_dress(...)`.
- Prints:
  ```
  "\tRPA correlation energy (Quasi-boson ... -sum w Y^+ N Y) ...                  E_QBA = <round 6>"
  ```
  and, if Density, the E_sig, E_rho, (D) and E_RPA lines verbatim from lines 49-57.

**`RPA_energy_QBA(Params,Orb,Orb_1p1h,E,X,Y)`** (lines 71-142), per (J,P):

```
D_ph = O_h - O_p (isospin of the pair, from Orb) ; for each nu: norm = sum D (|X|^2-|Y|^2)
skip if |norm| < 1e-8 ; if |Im E_nu| > 1e-8: N_drop++, skip
E_QBA -= (2J+1) Re(E_nu) sum_ph D_ph |Y_ph,nu|^2
```

If `N_drop > 0`, print the two warning lines (lines 137-138). Compared with base: the dev form is the same formula under the new normalization, but it also excludes zero-norm and complex solutions.

**`RPA_energy_dress`** (lines 150-214), same filter:

```
E_dress = sum (2J+1) Re( Y_nu' dA Y_nu + X_nu' dB Y_nu )
```

`dot(a,b) = a^H b`.

**`RPA_energy_rho`** (lines 225-245):

```
E_sigma = V2b_correlation_energy(Orb,Orb_NN,Sigma,V_NN)
E_rho   = sum_{a,b: l,j equal} (j_a+1) (h_N.p[a,b] (Rho.p-Rho_ref.p)[b,a] + n analog)
```

`V2b_correlation_energy` (`common/V2b.jl:572-634`):

```
E = sum_{a,b,c,d,J} [ 1/4 (2J+1) (V_pp Sigma_pp + V_nn Sigma_nn) + (2J+1) V_pn Sigma_pn ]
```

- a,b,c,d run over all orbitals with `2n+l` pair sums <= N2max and matching parity.
- J runs from `max(|j_a-j_b|,|j_c-j_d|)/2` to `min(...)/2`.
- It uses `O2b_pp/pn/nn` accessors. Threading uses `:static` + threadid accumulators, which is OK.

**Fortran:** new module `src/rpa/nj_rpa_energy.f90`: `RPA_energy` ~70 lines, `RPA_energy_QBA` ~50, `RPA_energy_dress` ~50, `RPA_energy_rho` ~30, `V2b_correlation_energy` ~50 (core `nj_v2b`). Return a derived type `rpa_ecorr_t (E_QBA,E_sigma,E_rho,E_dress,E_dens)` with `ieee_value(quiet_nan)` for unset terms. Delete `HF_RRPA_energy*` and `bosonic_sum` (~130 lines).

### 1.7 `src/RPA/RPA_orbitals.jl`: `Orb1B_RPA_OBDM` and `RPA_OBDM_orbitals_make` (lines 50-118)

- **Counterpart:** base `orbitals_ph_make` (particle and hole `SPState` lists). Fortran: `nj_orbitals::orbitals_ph_make`.
- **Definition:** particle if `|pO| < Tol_par`, hole if `|pO| > Tol_hol` (strict; the same for n), in Orb order.
- **Type:** `Orb1B_RPA_OBDM(p::pnOrbVector, h::pnOrbVector, N_p::pnInteger, N_h::pnInteger)`. In Fortran use `pnsvector_t`-like arrays of `orb1b_t` (needs real occupations, section 2), or index arrays. About 40 lines.

---------------------------------------------------------------------------------------------------

## 2. Shared helpers in `src/common` used by this group, and SHARED Fortran interface changes

| dev helper | dev location | Fortran now | action |
|---|---|---|---|
| `Orb1B.pO,nO::Float64` | `import_structures.jl:133` (base: Int64) | `nj_types::orb1b_t` `integer :: pO,nO` | **change to `real(dp)`**; NAT orbitals carry fractional occupations. Audit every `real(Orb%pO,dp)` and integer comparison (for example `orbitals_ph_make`, HF, EMPM). |
| `orbitals_import(Params,path)` / `orbitals_export(Params,Orb,path)` | `orb.jl:111-156` | none (`orbitals_export(Params,Orb)` is the HRF table) | new: stream `Int64 a,n,l,j; Float64 pO,nO` x a_max. Write `integer(int64)`. |
| `orbitals_update(Orb,Rho)` | `orb.jl:91-109` | none | new: `Orb(a)%pO = Rho%p(a,a)` etc., `intent(inout)`. |
| `OBDM_natural_orbitals(Params,Orb,Rho) -> (Orb_new,C,Rho_diag)` | `orb.jl:552-649` | none (`nj_hf_mbpt` has a private variant at line 670) | **new core routine** (~100 lines), see section 2.1. Shared with HF-MBPT and `RPA_export`. |
| `O1b_transformation(O,C) = C'OC`, `O1b_transpose`, `O1b_diagonal`, `O1b_deepcopy`, `O1b_unit_matrix` | `O1b.jl:1-39` | `mtmm` in `nj_linalg` | add thin wrappers to `nj_o1b` (~30 lines) |
| `O1b_import/export` | `O1b.jl:41-102` | `nj_o1b` | layout unchanged (p then n; (l,j)-equal elements, row a outer, b inner) |
| `O2b_initialize(Params,Orb)` without Orb_NN output | `O2b.jl:1` | `O2b_initialize(Params,Orb,O_NN,Orb_NN)` | optional `Orb_NN` |
| `O2b_initialize_empty()` | `O2b.jl:104` | none | represent as unallocated `o2b_t`; add `logical function o2b_is_empty` |
| `O2b_transformation(...;Report)` (mutates `O_NN`) | `O2b.jl:615` | `subroutine O2b_transformation(Params,Orb,Orb_NN,O_NN,C)` `intent(inout)` | unchanged semantics; add an optional `Report` |
| `V2b_read(Params)`, `V3b_no2b_read(Params)` | `V2b.jl:1`, `V3b_no2b.jl:381` | take `Orb` too | the HF/core group decides; ScRPA needs only `(V_NN_bare,Orb_NN)` and `(V_NNN_bare,Orb_NNN)` |
| `H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho)` | `Hamiltonian.jl:443-597` | `nj_hf_rrpa::HF_RRPA_h1b` (and HF's own h) | **move to new core module `src/core/nj_hamiltonian.f90`** as `H1b`. Identical math to base `HF_RRPA_h1b` (b,e with equal l,j; c,f with equal l,j; CMS handling `"CMS1+2B"`: `T*(1-1/A)`, `"CMS2B"`: none, else `T`). Result order: lower pairs (a >= d) written symmetrically. Also used by HF-MBPT (ScH). |
| `H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho;Report)` | `Hamiltonian.jl:599-760` | `nj_v2b::V2b_residual_no2b` (in place) | wrapper: copy V_NN, then `V2b_residual_no2b`. Same NO2B contraction, `H = cRes*(V + V3N.Rho)`. **Verify element-wise once** (dev reorganized the loops: pairs of `Orb_NN.Ind[2]`, pn orderings with `Phase_ab=(-1)^((j_a+j_b)/2-J)`, `pnSum_abde=(ME001/2+ME111/6+ME113/3)(pRho+nRho)+(ME101+ME011)/sqrt12 (nRho-pRho)`, etc., lines 726-729). |
| `H1b_H2b(...,Rho,Sigma;Report) -> (h_N,H_NN)` | `Hamiltonian.jl:995-1448` | none | new (~200 lines): `H_NN = H2b(Rho)`; `h = H1b(Rho)` part (lines 1123-1224, identical to H1b) **plus the V3N.Sigma term** (lines 1226-1336, section 2.2). The V_NN.Sigma term at lines 1338-1441 is commented out; do not port. |
| `H1b(...,Rho,Sigma)` (9 args, with V_NN.Sigma) | `Hamiltonian.jl:2-441` | none | dead code in dev; do not port |
| `E1b`, `E2b` | `Hamiltonian.jl:762, 905` | none | unused by this group |
| `V2b_correlation_energy` | `V2b.jl:572` | none | new in `nj_v2b` (section 1.6) |
| dense NNN lookup `Orb3B.Map[ab](TT,J-Jmin+1,c)` and `nlj` | `import_structures.jl:246` | `orb3b_t` hash map | the H1b_H2b Sigma term enumerates triples through `Map`; in Fortran build the same triple table once via hash look-ups (no type change required), or add `Map` to `orb3b_t` (owned by the V3b group) |
| `orbitals_1p1h_make`, `orb2b_1p1h_index`, `Orb2B_1p1h(N,Ind,Dic)` | `orb.jl:333-443` | none | RPA group; ScRPA needs `N(P,0:J)`, `Ind(P,J)%v` of `(a_p,a_h,T)` and the index lookup. Pair order for each (J,P): a outer, b inner over Orb; per (a,b) first T=-1 then T=+1; condition `O_b > Tol_hol && O_a < Tol_par && O_b-O_a > 1e-8`. |
| `JP_initialize(J_max)` | `functions.jl:191` | loops | order J outer, P inner |
| `rho_chR2(Params,Orb,C,Rho)` | `rho.jl:374` | `nj_obdm::OBDM_chR2` | rename only |
| `hungarian(cost)` | Julia Hungarian.jl | `nj_linalg::hungarian` | reuse |

### 2.1 `OBDM_natural_orbitals` (exact recipe, `orb.jl:552-649`)

1. `pRho = (Rho.p + Rho.p')/2` (n alike).
2. Zero every [a,b] with `j_a != j_b || l_a != l_b`.
3. `pRho[a,a] += 1e-10*a` (1-based a; n alike).
4. `pOcc,pC = eigen(Symmetric(pRho))`: ascending, LAPACK `syevr`. Fortran: `eig_sym`.
5. `cost = 1 - |C|`; `perm = hungarian(cost)`; `C = C[:,perm]`; `Occ = Occ[perm]`.
6. For each a: if `C[a,a] < 0`, then `C[:,a] *= -1`.
7. `Rho_out = diagm(Occ)`.
8. `Orb_new[a]`:
   - `l_a = round(sum_b |pC[b,a]|^2 l_b)` and `j_a = round(sum_b |pC[b,a]|^2 j_b)`, using the **proton** C;
   - `n_a` = the number of b < a with the same (l,j) in the **input** Orb;
   - `pO = pOcc[a]`, `nO = nOcc[a]`.
9. Return `(Orb_new, O1B(pC,nC), O1B(pRho_diag,nRho_diag))`.

`HF_ScRPA_OBDM_solve` uses only C; `RPA_export` uses only C.

### 2.2 V3N.Sigma term of `H1b_H2b` (`Hamiltonian.jl:1226-1336`, identical to lines 224-334)

**Triples.** For every (P3B, l_c, J) block of `V_NNN[P3B,l_c+1,J+1,T][j_index]`, the triples are `(a,b,c, ab=Orb_NN.Map[a,b,J+1,2], i_0, i_1, i_3)` for a >= b and c with `Map_ab[1,J_i,c] != 0`, where `j_index = 1` if `j_c == |2l_c-1|`, else 2. They are enumerated a outer, b <= a, J_i, c.

**Loop.** For `J = 0..N2max+1`, `P2B = 1,2`, `l_c = 0..Nmax`, `P3B = parity(P2B-1+l_c)`, and each triple abc with each def <= abc in the same block:

```
ME001=V_T1[idx(abc_0,def_0)], ME101=V_T1[idx(abc_1,def_0)], ME011=V_T1[idx(abc_0,def_1)], ME111=V_T1[idx(abc_1,def_1)], ME113=V_T3[idx(abc_3,def_3)]
   idx(B,K) = B>=K ? K + B(B-1)/2 : B + K(K-1)/2
ppS, nnS = Sigma_dense[de,ab] (pp/nn, a>=b packing) ; x2 if a!=b ; x2 if d!=e
pnS_abde = pnSigma[de_pn,ab_pn] ; pnS_bade = (a!=b) Ph_ab pnSigma[de_pn,ba_pn] ; pnS_abed = (d!=e) Ph_de pnSigma[ed_pn,ab_pn] ; pnS_baed = (a!=b&&d!=e) Ph_ab Ph_de pnSigma[ed_pn,ba_pn]
   Ph_ab = (-1)^((j_a+j_b)/2 - J)
S00 = abde+bade+abed+baed ; S10 = abde-bade+abed-baed ; S01 = abde+bade-abed-baed ; S11 = abde-bade-abed+baed
pSum = ME113 ppS + (ME113+2ME111)/3 nnS + 2[ ME001 S00 - (ME101 S10 + ME011 S01)/sqrt3 + (ME111+2ME113)/3 S11 ]
nSum = ME113 nnS + (ME113+2ME111)/3 ppS + 2[ ME001 S00 + (ME101 S10 + ME011 S01)/sqrt3 + (ME111+2ME113)/3 S11 ]
h_cf += pSum ; if def != abc: h_fc += pSum      (n alike)
```

**Finally:** `h_N[c,f] += 0.25/(j_c+1) * h_cf` for all (l,j)-equal ordered (c,f).

Threads use task-local `a_max x a_max` accumulators, so there is no race; in Fortran use an OpenMP reduction or per-thread arrays.

### 2.3 `nj_params` / `nj_input` (shared with the RPA group)

- Replace `rpa_params_t` (currently only `Ortho`) with the full dev `RPA_Parameters` (section 6).
- Delete `rrpa_params_t`, `Calc%RRPA` and namelist `&rrpa`.
- Namelist `&rpa` gets all fields (no case clashes).
- Add run action `HF_ScRPA` to `nj_dispatch_rpa`.
- Delete `HF_RRPA`, `HF_RRPA_transition_densities` and `HF_RRPA_transition_currents`. Their dev replacements are the generic `RPA_external` routines on `io/<Path>/<RPA.Name>`; the RPA group owns them.
- `nj_types::obdm_iteration_pairs_t` becomes unused and should be deleted. Dev still declares `OBDM_Iteration_Pairs` (dead).
- `nj_hf_rpa` helpers used only by RRPA become removable after the RPA group's rewrite: `RPA_eigen_solve_c` and the `Zero_Negative` branch of `RPA_normalize`. Check with grep.

---------------------------------------------------------------------------------------------------

## 3. What to change in the Fortran port, routine by routine

| Fortran routine (now) | action | dev source | est. lines |
|---|---|---|---|
| `nj_hf_rrpa_driver::HF_RRPA` | replace by `nj_scrpa::HF_ScRPA` | RPA.jl:72-141 | 40 |
| `nj_hf_rrpa_driver::HF_RRPA_solver`, `HF_RRPA_solve` | replace by `nj_scrpa::HF_ScRPA_solver` (+ `scrpa_beta`, `ScRPA_terminate_summary`) | RPA_solver.jl:94-418 | 300 |
| `nj_hf_rrpa::HF_RRPA_OBDM_iteration` (+ `obdm`, `pair_sum`, `reorder`, `occupation_check`, `transform_amplitudes`, `iteration_list`, `add_pairs`) | replace by `nj_rpa_rho::HF_ScRPA_OBDM_solve`, `RPA_OBDM_oY4`, `RPA_OBDM_transform_amplitudes`, `RPA_OBDM_orbitals_make` | RPA_rho.jl:283-776, RPA_orbitals.jl | 390 |
| (none) | new `nj_rpa_sigma::RPA_sigma`, kernels, `RPA_sigma_mix` | RPA_sigma.jl | 390 |
| `nj_hf_rrpa::HF_RRPA_h1b` | move to core `nj_hamiltonian::H1b` | Hamiltonian.jl:443-597 | 0 (move) + 20 |
| `nj_hf_rrpa::HF_RRPA_V2b`, `HF_RRPA_V2b_transform` | replace by `nj_hamiltonian::H2b` (copy + `V2b_residual_no2b`) and `nj_scrpa::RPA_Hamiltonian(_HF_reference)` | Hamiltonian.jl:599-760, RPA_Hamiltonian.jl | 100 |
| (none) | new `nj_hamiltonian::H1b_H2b` (V3N.Sigma) | Hamiltonian.jl:995-1336 | 200 |
| `nj_hf_rrpa::HF_RRPA_energy`, `_bosonic`, `bosonic_sum`, `_density` | replace by `nj_rpa_energy::RPA_energy`, `_QBA`, `_dress`, `_rho` + `nj_v2b::V2b_correlation_energy` | RPA_energy.jl, V2b.jl:572 | 250 |
| `nj_hf_rrpa::HF_RRPA_allocate`, `HF_RRPA0_diagonalize`, `HF_RRPA_diagonalize`, `HF_RRPA_rM` | delete; use the RPA group's `RPA_allocate(…,Rho,Sigma,TBCF)`, `RPA_diagonalize`, `RPA_rM`, `RPA_rB` | — | -600 |
| `nj_hf_rrpa_export` (entire module) | delete; use the RPA group's `RPA_export` | RPA_export.jl | -297 |
| (none) | core: `OBDM_natural_orbitals`, `orbitals_import/export(path)`, `orbitals_update`, O1b wrappers | orb.jl, O1b.jl | 180 |
| `nj_params`, `nj_input`, `nj_dispatch_rpa`, `files.mk`, `nj_types` | params/namelist/dispatch/types (section 2.3) | import_structures.jl:41-58 | 60 |

- New code is about 1900 lines, of which about 400 are in shared core.
- About 1700 lines are deleted (`nj_hf_rrpa*.f90`).
- `files.mk` order: `nj_rpa_sigma`, `nj_rpa_rho`, `nj_rpa_energy` after the RPA group's allocate/diagonalize modules, then `nj_scrpa`, then the dispatcher. `nj_hamiltonian` goes into the core list after `nj_v2b` and `nj_v3b`.

**Parallelization notes for parity:**
- `RPA_OBDM_oY4` and `transform_amplitudes`: OpenMP over (J,P) blocks (largest first, `schedule(dynamic)`) calling single-threaded BLAS (MKL: `mkl_set_num_threads_local(1)` inside the region), or a serial block loop with threaded BLAS.
- The element tasks are independent, so `parallel do` is fine.
- Sums inside a task are sequential in Julia; keep the same order (h in `Orb_ph` order, then J ascending) for bit-level parity.

---------------------------------------------------------------------------------------------------

## 4. Output layout changes

- **Output root:** `IO/<Path>/RRPA/{Densities,Transitions/E0..E3,Spectra,Amplitudes,...}` is replaced by `io/<Path>/<RPA.Name>/{bin,densities,transitions/E0..E3,spectra,amplitudes}`. All lowercase; the directory is deleted and recreated at the start.
- **Converged run:** files are written by the generic `RPA_export` (layout defined by the RPA group spec):
  - `RPA_summary.dat`
  - `bin/RPA_X.bin`, `bin/RPA_Y.bin`, `bin/RPA_E.bin`
  - `spectra/RPA.dat`, `spectra/RPA_plot.dat`
  - `amplitudes/RPA_amplitudes_ordered.dat`, `amplitudes/RPA_amplitudes.dat`
  - `densities/<Name>_radial_densities.dat` (+ radii appended to `RPA_summary.dat`)
  - `densities/RPA_occupation.dat`
  - `bin/orb.bin`, `bin/C.bin` (= C_LHO_NAT), `bin/rho.bin` (NAT basis). With `RPA_Orb = true` instead: `orb.bin`, `C.bin`, `rho.bin` and `H1b.bin` of the natural orbitals of the exported OBDM.
  - `bin/sigma.bin` only if `TBCF = true`
  - `transitions/E0..E3/RPA_*.dat`
  - For ScRPA, `RPA_export` receives `C = C_LHO_NAT`, `Rho = Rho_LHO` and `h = h_N` (NAT). In `RPA_orbitals_export`, `Rho_RPA = C' Rho_LHO C` is then already diagonal, so `C_NAT` is approximately I up to Hungarian signs.
- **Non-converged run:** only `RPA_summary.dat` with the termination summary (section 1.2); every other subdirectory stays empty. This is new behaviour; base always exported.
- **Inputs read:** `io/<Path>/<Reference>/bin/{orb.bin, C.bin, rho.bin, H1b.bin, H2b.bin}` replace `IO/<Path>/Bin/{V2B_HF.bin, C_HF.bin}`. The HF group must write these files in the dev layout.
- **Removed outputs:** base RRPA-specific files such as `Bin/C_RRPA.bin` and `Bin/Rho_RRPA.bin` (read by `nj_hf_rpa_external` with Kind = 'RRPA'), and the RRPA occupation and single-particle spectrum tables of `HF_RRPA_export_*`.

---------------------------------------------------------------------------------------------------

## 5. Removed features

- The run action `HF_RRPA`, the namelist `&rrpa` and the `RRPA_Parameters` fields `ScV3N`, `ScOBH` and `dOBDM`.
- The `HF_RRPA_transition_densities` and `HF_RRPA_transition_currents` actions. They are replaced by the generic `RPA_external` routines on the ScRPA output directory (RPA group).
- RRPA(0) start with the phonon-list machinery (`orbitals_ph_make`, `orbitals_one_phonon_make`, `HF_RRPA_phonon_count`), `HF_RRPA0_allocate/diagonalize`, the iteration-variant diagonalize with `Zero_Negative`, and the RRPA-specific spurious orthogonalization (`HF_RRPA_spurious_*`). These are now generic in `RPA_diagonalize`.
- `HF_RRPA_OBDM_occupation_check` bail-out, the "lower precision" exits, and the outer-loop `E_old` tests (port "kept on purpose" bullet: `E_old` never updated).
- `HF_RRPA_energy_density` and `HF_RRPA_energy_bosonic` prints, and `HF_RRPA_collectivity`.
- `HF_RRPA_export` with all of its files.
- `obdm_iteration_pairs_t`.

---------------------------------------------------------------------------------------------------

## 6. New or changed parameters (`RPA_Parameters`, `import_structures.jl:41-58`, shared with RPA)

| field | default | meaning for HF-ScRPA |
|---|---|---|
| `Name` | `"RPA"` | output dir `io/<Path>/<Name>`; must differ from `Reference` |
| `Reference` | `"HF"` | reference solution dir |
| `Reference_TBCF` | false | unused by ScRPA (written to the summary only) |
| `RPA_Orb` | false | `RPA_export`: export natural orbitals (+ H1b.bin) as a new reference |
| `Ortho` | true | spurious 1- orthogonalization in `RPA_diagonalize` |
| `Transitions` | `"Linear"` | RPA group (`RPA_rM`) |
| `Tol` | 1e-5 | outer criterion: relative change of E_QBA (was `RRPA.Tol = 1e-7`, OBDM measure) |
| `Tol_hol` | 0.7 | hole if O > Tol_hol (1p1h, Orb_ph) or O >= Tol_hol (sigma) |
| `Tol_par` | 0.1 | particle if O < Tol_par (1p1h, Orb_ph) or O <= Tol_par (sigma) |
| `Imax` | 150 | max outer iterations (inner loop fixed at 150) |
| `OBDM_diagonal` | false | keep only the diagonal of the new OBDM (replaces `RRPA.dOBDM`) |
| `OBDM_alpha` | 0.6 | linear mixing of the inner OBDM iteration |
| `ScH` | true | self-consistent h and V from the bare H with the current OBDM (replaces `ScOBH` and `ScV3N`) |
| `TBCF` | false | use the 2-body correlation function Sigma in the outer loop (A,B, h, energy, export sigma.bin) |
| `TBCF_alpha` | 0.8 | linear mixing of Sigma between outer iterations |
| `TBCF_beta` | 5.0 | ramp length: `beta_i = tanh((i/TBCF_beta)^2)`, set to 1 above 0.95 |
| `DEV` | false | experimental: adds `RPA_allocate_SE` (1p1h-2p2h self-energy) to A and to `E_dress` |

Hard-coded (not parameters) in `HF_ScRPA_OBDM_solve`: inner tolerance 1e-8, 150 iterations, 2-cycle test after iteration 100 with relative threshold 1e-3, and the D threshold 1e-8 in oY4.

---------------------------------------------------------------------------------------------------

## 7. Julia bugs and races in dev (this scope), and the status of the port notes

**B1, likely bug.** `RPA_sigma.jl:208-210`: the pn branch calls `RPA_sigma_hpph(..., J, Q, **Y**, Orb_1p1h)`, while the pp and nn branches (lines 152-161, 272-281) pass `X` for the same h-p-p-h configuration; the kernel argument is named X. Physically `<a_h^+ a_p a_p'^+ a_h'>` is proportional to X X*. It affects every TBCF quantity (A,B dressing, h V3N.Sigma, E_sigma/E_dens, sigma.bin), and also `E_dens` for TBCF = false, because the final Sigma is always computed. Owner decision: port faithfully (bit parity) or fix (documented deviation).

**B2, hazard.** `RPA.jl:104-105` / `HF_ScRPA`: `rm -r io/<Path>/<Name>` happens before the reference is read, so `Name == Reference` destroys the reference. The port should `fatal` in that case. This is a harmless deviation.

**B3, inconsistency.** The particle and hole tests differ: `orbitals_1p1h_make` and `RPA_OBDM_orbitals_make` are strict (and use `abs`), while `RPA_sigma` is inclusive and uses the **current** NAT occupations. `Orb_1p1h` and `Orb_ph` are frozen at the reference. If an occupation crosses a tolerance during the iteration, `RPA_sigma` classifies a pair differently from the 1p1h space. Look-ups then return 0, but the disconnected `-O_a O_b` and `-phase O_a` terms are still added. Keep it faithfully, and list it under "kept on purpose".

**B4, assumption.** `RPA_sigma`, `RPA_OBDM_oY2/oY4` and `RPA_energy*` loop J up to `2*Nmax+1` and index `Orb_NN.N`, which is sized `N2max+1` (`O2b_initialize`). This is out of bounds if `N2max < 2*Nmax` (compare port note 3). The port should loop to `min(2Nmax+1, N2max+1)` for the `Orb_NN`-indexed parts.

**B5, inefficiency.**
- The final `RPA_Hamiltonian` call for TBCF (line 365) repeats the last in-loop one with identical inputs.
- The bare V_NN and V_NNN are read even for `ScH = false`, where they are unused.
- `HF_ScRPA_OBDM_solve` re-reads rho.bin and C.bin on every outer iteration.
- None of these changes results; optional to skip in the port (owner's choice for V_NNN memory on 208Pb).

**B6, side effect.** `RPA_Hamiltonian` with `ScH && TBCF` round-trips the caller's Sigma in place (NAT -> LHO -> NAT). This is harmless up to rounding; the port may transform a copy.

**B7, global state.** `RPA_OBDM_oY4` and `RPA_OBDM_transform_amplitudes` call `BLAS.set_num_threads(1)` and restore it afterwards. This is a global setting, but there is no data race. The Fortran port needs a thread-safe equivalent.

**B8, semantics.** The non-convergence fallback of the inner loop is unmixed and starts from `Rho_ini`. If the inner loop never converges (strict 1e-8), the outer loop effectively becomes an unmixed fixed-point iteration on the OBDM. This is not a bug, but it explains a possible outer non-convergence.

**B9, ScH = false with TBCF.** The reference H1b and H2b are used and Sigma is ignored in h,V, so TBCF acts only via A,B and the energy. This is probably intended.

**Races.** None found in the scope:
- `RPA_sigma`, `oY2` and `oY4` write per-element or per-block.
- `H1b` and `H1b_H2b` use task-local accumulators.
- `V2b_correlation_energy` uses `:static` with threadid.
- `RPA_energy_*` work per block.
- The base race of port note 6 (`HF_RPA_OBDM`) is gone in dev's `RPA_OBDM_oY2`, which uses per-(J,P) accumulators. The "RRPA counterpart" part of port note 6 is obsolete.

**Port-note status:**

| port note | status |
|---|---|
| 6 (serial OBDM accumulation) | obsolete for ScRPA and oY2 (dev fixed upstream) |
| 7 (`RRPA.OBDM` field) | obsolete (code removed) |
| 8 (Verbose silencing of the RRPA loop) | replace by passing `Report = .false.` where dev does. Dev still prints the "Evaluating ... h_N", "Allocating RPA matrices", "Constructing Hamiltonian" and `@time` lines; reproduce them or route them through Verbose consistently. |
| 9 (0.5-1e-8 norm cut) | obsolete: dev has no norm cut; spurious handling is in `RPA_diagonalize` (RPA group) |
| "kept on purpose": `HF_RRPA_OBDM` neutron-hole `Y[pg,mu]` | obsolete (fixed in dev) |
| "kept on purpose": `HF_RRPA_solve` `E_old` | obsolete (code removed) |

Add B1 (if ported faithfully) and B3 as new "kept on purpose" entries.

---------------------------------------------------------------------------------------------------

## 8. Decisions for the code owner

1. B1: port the pn `hpph(…,Y)` faithfully, or fix it to X?
2. Remove the Fortran HF-RRPA completely, as dev did, or keep it temporarily as a legacy action? The 208Pb RRPA comparison runs (Julia base vs Fortran) are still in progress, and the uncommitted `nj_v3b.f90` speed-up is used by them.
3. Expose the hard-coded inner-loop constants (1e-8, 150, 1e-3 cycle threshold) as extra namelist fields (deviation from "namelists mirror Julia structs"), or keep them hard-coded?
4. Port `RPA.DEV` (`RPA_allocate_SE`, about 650 lines, experimental, RPA group) now, or reject `DEV = .true.` with `fatal`?
5. Guard `Name == Reference` with `fatal` (B2)?
6. Skip reading the bare V_NNN and V_NN when `ScH = false` (memory and time; get Orb_NN from `O2b_initialize`)?
7. Bit-parity versus speed for oY4 and transform: Julia uses complex BLAS per block. Is real-only arithmetic (amplitudes are real for stable solutions) acceptable when results agree to about 1e-13?
8. Port the dead dev code `RPA_OBDM_sigma` and `H1b(...,Sigma)` (the V_NN.Sigma variant)? The recommendation is no.
