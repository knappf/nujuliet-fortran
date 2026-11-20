# Dev port spec: HF + HF-MBPT

Julia base `2595e82` (`Source/HF/*`) → Julia dev `f3f7a33` (`src/HF/*`). Fortran port: `src/hf/*` on branch `dev-port`.
This is read-only analysis. Line numbers refer to `NuJuliet_dev/src/HF/<file>` unless stated otherwise.

Relevant dev commits that touch `src/HF` (in order): `4893439` (first dev upload), `9f80777` ("HF-MBPT(3) hole-hole fix, RPA big update,
HF options"), `fcb42d9`/`27b274c` (mass update). The "hole-hole fix" plus the later commits correct **isospin labels of the
single-particle energies** in the MBPT(3) energy (see §8.3). Port the final `f3f7a33` state only.

---------------------------------------------------------------------------------------------------------------------------

## 0. Summary of what changes

| Area | Base (ported now) | Dev (target) |
|---|---|---|
| Occupations | `Orb1B.pO/nO::Int64` from naive LHO filling, fixed for the whole run | `Float64`. Re-determined **every HF iteration** from the HF SPEs (`HF_orbitals_make`). Fractional occupation of the last (sub)shell is allowed unless `ESH = true` |
| HF density | `HF_density_operator(a_max,C,Orb)`, sum over orbitals with `pO == 1` in the same (l,j) block | `HF_density_operator(C,Orb)` = Σ_a O_a C[:,a] C[:,a]ᵀ (weights = real occupations) |
| HF ordering | `HF_orbital_ordering` (`HF_Ordering.jl`) | merged into `HF_orbitals_make` (`HF_orbitals.jl`) together with the new occupation rule |
| HF field / energy | full loops | same formulas, loops pruned by `l_max/j_max` of occupied orbitals (exact, performance only) |
| Radial potential | `HF_Radial_Potential` (dead code) | new `HF_radial_potential` (non-local V(r,r') per partial wave), enabled by `HF.V_HF` |
| HF-EFA | `HF_EFA.jl`, `HF.EFA` params, `HF-EFA/` dirs | **removed** |
| HRF export | `Bin/V2B_HF.dat` + HRF orbitals table if `HF.HRF` | **removed** (field `HRF` still exists in `HF_Parameters`, unused) |
| HF-MBPT | MBPT(2) energy + MBPT OBDM, NAT orbitals by max-overlap reorder | MBPT(2) **+ MBPT(3)** energy (pp, ph, hh), OBDM with NAT orbitals by Hungarian assignment (`OBDM_natural_orbitals`), occupation (idempotency) gate, optional export of a correlated reference (`HF-MBPT/bin/*`), optional 2-body correlation function (TBCF) up to 3rd order, optional polarised Hamiltonian (ScH) |
| Output layout | `IO/<Path>/Bin/{C_HF,Rho_HF,h_HF,V2B_HF}.bin`, `IO/<Path>/HF/HF_Summary.dat`, `HF/Densities/...` | `io/<Path>/HF/bin/{orb,C,rho,H1b,H2b}.bin`, `io/<Path>/HF/HF_summary.dat`, `HF/densities/...`, `io/<Path>/HF-MBPT/bin/...` |

---------------------------------------------------------------------------------------------------------------------------

## 1. Parameters

Dev `src/import_structures.jl:16-30`:

```julia
Base.@kwdef struct HF_MBPT_Parameters
    Export::Bool = false
    TBCF::Bool = false
    ScH::Bool = false
end

Base.@kwdef struct HF_Parameters
    Tol::Float64 = 1e-7
    Imax::Int64 = 100
    V_HF::Bool = false
    BMF::Bool = true
    HRF::Bool = false
    ESH::Bool = false
    MBPT::HF_MBPT_Parameters = HF_MBPT_Parameters()
end
```

Removed: `HF_EFA_Parameters` (`A`, `Z`, `BMF`) and `HF_Parameters.EFA`.

Meaning:
* `V_HF` – export the non-local HF radial potential (`HF_radial_potential`, §5).
* `BMF` – beyond-mean-field: build/export residual NN interaction and run HF-MBPT (unchanged meaning).
* `HRF` – still in the struct, **never read** in dev (the HRF block of `HF_solver` was deleted).
* `ESH` – "enforce (sub)shell closure": when filling by HF SPE order, skip orbitals that cannot be filled completely (§3.4).
* `MBPT.Export` – after the HF-MBPT OBDM, export a correlated HF-MBPT reference (natural orbitals, NAT-basis h, V, ρ [, σ])
  to `io/<Path>/HF-MBPT/bin/` (§8.7). Consumed by TDA/RPA with `Reference = "HF-MBPT"`.
* `MBPT.TBCF` – (only with `Export`) compute the MBPT(2)+(3) two-body correlation function σ and export `sigma.bin`;
  also selects `H1b_H2b(…,Sigma)` instead of `H1b`/`H2b` when `ScH`.
* `MBPT.ScH` – (only with `Export`) build the NAT-basis Hamiltonian self-consistently from the correlated density
  (`H1b`/`H2b`, or `H1b_H2b` with σ); otherwise use HF residual V (§8.7).

### Fortran changes (`src/core/nj_params.f90`, `src/app/nj_input.f90`)

```fortran
type, public :: hf_mbpt_params_t
    logical :: Export = .false.
    logical :: TBCF   = .false.
    logical :: ScH    = .false.
end type
type, public :: hf_params_t
    real(dp) :: Tol  = 1.0e-7_dp
    integer  :: Imax = 100
    logical  :: V_HF = .false.
    logical  :: BMF  = .true.
    logical  :: HRF  = .false.     ! accepted, no effect (as in dev)
    logical  :: ESH  = .false.
    type(hf_mbpt_params_t) :: MBPT
end type
```
* delete `hf_efa_params_t` and the `EFA` component.
* namelist `&hf` → `Tol, Imax, V_HF, BMF, HRF, ESH`; new namelist `&hf_mbpt` → `Export, TBCF, ScH` (mirrors the nested
  struct the same way `&hf_efa` mirrored `EFA`). Delete `&hf_efa` (decide whether an old input containing `&hf_efa` is an
  error or silently ignored – see open questions). Note `TBCF`/`ScH` local namelist variables will also be needed by the RPA
  group's `&rpa`; reset them from `C0%HF%MBPT` before reading `&hf_mbpt`.

---------------------------------------------------------------------------------------------------------------------------

## 2. `HF.jl` → `nj_hf_driver.f90 :: HF`

Base: `Source/HF/HF.jl`. Fortran: `src/hf/nj_hf_driver.f90 :: HF` (lines 61-87).

Changes:
* Output root `"IO/"` → `"io/"` (shared `io_path`, see §10).
* Directories (dev `HF.jl:36-43`): `rm -r io/<Path>` if it exists, then
  `mkdir io/<Path>`, `io/<Path>/HF`, `io/<Path>/HF/bin`, `io/<Path>/HF/densities`.
  No `Bin/`, no `HF/Densities`, no `HF-EFA/…`.
  (Behaviour kept from base: the whole `io/<Path>` is wiped, including other solvers' output.)
* Header/printing unchanged (`print_run_header`, still used by BCS/HFB).

Size: ~15 lines changed.

---------------------------------------------------------------------------------------------------------------------------

## 3. `HF_solver.jl` (+ `HF_orbitals.jl`) → `nj_hf_driver.f90 :: HF_solver`, `nj_hf.f90 :: HF_solve`

Base: `Source/HF/HF_Solver.jl`, `Source/HF/HF_Ordering.jl`.

### 3.1 `HF_solver` sequence (dev `HF_solver.jl:1-75`)

```
Orb   = orbitals_make(Params)                      # LHO orbitals, Float64 occupations (naive filling, may be fractional)
T     = T1b(Params,Orb)
V_NN, Orb_NN   = V2b_read(Params)                  # signature lost Orb (common group)
V_NNN, Orb_NNN = V3b_no2b_read(Params)             # idem
h_N, C, Rho, Iteration = HF_solve(...)             # NOTE: mutates Orb in place (occupations), §3.3
E_HF = HF_energy(...)
T_HF = T1b_energy(Params,Rho,Orb,T)
HF_summary(Params,E_HF,T_HF,Iteration)
rho_export(Params,Orb,"io/<Path>/HF/HF_summary.dat","io/<Path>/HF/densities/HF_radial_densities.dat",Rho,C)
HF_SPS_summary(Params,h_N,Orb)
HF_export(Params,C,O1B(C.p'*Rho.p*C.p, C.n'*Rho.n*C.n),h_N,Orb)
if HF.V_HF:  HF_radial_potential(Params,Orb,C,T,h_N)
if HF.BMF:
    V_NN_den = V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN)   # dev returns a deepcopy, V_NN untouched
    V_NN_res = O2b_transformation(Params,Orb,Orb_NN,V_NN_den,C)              # in place on V_NN_den
    O2b_export(Params,Orb_NN,V_NN_res,"io/<Path>/HF/bin/H2b.bin")
    if HF_MBPT_occupation_check(Params,O1B(C.p'*Rho.p*C.p, C.n'*Rho.n*C.n))
        HF_MBPT(Params,Orb,Orb_NN,h_N,V_NN_res)
        if HF.MBPT.Export:  HF_MBPT_export(Params,Orb,Orb_NN,Orb_NNN,C,T,h_N,V_NN,V_NN_res,V_NNN)
```

Differences vs base:
* `OBDM_export` → `rho_export` (rename only, common group; Fortran `nj_obdm::OBDM_export` keeps the same arguments).
* HF-EFA call removed. HRF export removed.
* `Bin/V2B_HF.bin` → `HF/bin/H2b.bin`.
* MBPT is gated by `HF_MBPT_occupation_check` and receives `h_N`.
* **`V_NN` must stay the bare interaction** after the residual is built, because `HF_MBPT_export` uses it (H1b/H2b).
  Fortran `V2b_residual_no2b` is `intent(inout)` on `V_NN` (base semantics). The driver must work on a copy
  (`V_NN_res = V_NN; call V2b_residual_no2b(…,V_NN_res,…); call O2b_transformation(…,V_NN_res,C)`), or nj_v2b must get a
  non-mutating variant (shared interface, §10). If memory matters, the copy is only needed when `HF.MBPT.Export` is true.

### 3.2 `HF_solve` (dev `HF_solver.jl:77-135`)

Loop unchanged except:
* initial `Rho = HF_density_operator(C,Orb)` with `C = I` and the **orbitals_make occupations** (may be fractional);
* after diagonalisation: `Orb, SPE, C = HF_orbitals_make(Params,SPE,C,Orb)` (replaces `HF_orbital_ordering`), then
  `Rho = HF_density_operator(C,Orb)`;
* printout text `"\tIteration number: %4d   dE = %12.9f MeV"` (was `"\tHF iteration number: …"`), start message
  `"Initializing the iteration of HF mean-field equations ..."`.
* Convergence `dE = (Σ|ΔSPE_p| + Σ|ΔSPE_n|)/(2 a_max)` and return `h_N = diagm(SPE)` unchanged.

Fortran `HF_solve`: `Orb` becomes `intent(inout)` (`type(orb1b_t), intent(inout) :: Orb(:)`), replace the ordering call by
`HF_orbitals_make`. `nj_bcs.f90` calls `HF_solve` (dev BCS no longer does – BCS group decides). ~15 lines.

### 3.3 In-place occupation update (important)

`HF_orbitals_make` ends with `Orb .= Orb_new` (`HF_orbitals.jl:137`), i.e. it overwrites the caller's vector. Hence after
`HF_solve` the `Orb` of `HF_solver` carries the **final HF occupations**, and everything downstream (HF_energy, T1b_energy,
rho_export, SPS table, `orb.bin`, MBPT particle/hole split, `O1b_export` block structure) uses them. The Fortran port must
reproduce this (intent(inout)).

### 3.4 `HF_orbitals_make(Params,SPE,C,Orb)` – new, dev `HF_orbitals.jl:1-140`

Step 1 – ordering: identical to base `HF_orbital_ordering` (expectation values of l and j of each eigenvector, assign to the
first free LHO slot with |l−⟨l⟩|<1e-8 and |j−⟨j⟩|<1e-8, reshuffle SPE and C columns). Keep the Fortran `fatal` on failure
(Julia leaves `Vector{Int}(undef)` garbage – port deviation still needed).

Step 2 – occupations, protons (neutrons identical with `N = A−Z`, `nSPE`, `nO`):

```julia
Z_occ = 0.0
while abs(Z - Z_occ) > 1e-8
    a_last, E_last = 1, 1e5
    for a in 1:a_max
        E = pSPE[a]
        if (E_last > E) && (pO_mask[a] == false)
            if (ShellClosure == true) && ((Z - Z_occ) < (Orb[a].j + 1))
                if a == a_max
                    throw("Expection: HF solver is not able to find closed-subshell solution ...")
                end
                continue
            end
            E_last = E; a_last = a
        end
    end
    j_hat = Float64(Orb[a_last].j + 1)
    pO_mask[a_last] = true
    Occ = min(1.0, (Z - Z_occ) / j_hat)
    pO[a_last] = Occ
    Z_occ += j_hat * Occ
end
```
(`ShellClosure = Params.Calc.HF.ESH`). I.e. repeatedly pick the lowest-SPE unfilled slot (strict `<`, ties → lowest index),
fill it with min(1, remaining/(2j+1)). With `ESH` skip slots with 2j+1 > remaining; the error is raised when the skipped slot
is `a_max` (Fortran: `call fatal(...)` with the same message text). Neutrons: same loop on `A - Z - N_occ`.
Finally `Orb_new[a] = Orb1B(a,n,l,j,pO[a],nO[a])` and copy into `Orb`.

Fortran: new `subroutine HF_orbitals_make(Params,SPE,C,Orb)` in `nj_hf.f90`, `intent(inout)` for all three. Implement as
`call HF_orbital_ordering(Orb,a_max,C,SPE)` (keep it public, BCS/HFB use it) followed by a new private
`HF_occupations_make(Params,SPE,Orb)`. ~70 lines.

Julia edge-case bugs in this loop (port verbatim or guard, see §11): fallback `a_last = 1` if no slot is eligible
(only possible when `a_max` is already masked) overwrites the occupation of slot 1; with `ESH` a valid candidate whose SPE is
above `SPE[a_max]` still triggers the throw when `a_max` is too large.

---------------------------------------------------------------------------------------------------------------------------

## 4. `HF_allocate.jl`, `HF_density_operator.jl`, `HF_energy.jl` → `nj_hf.f90`

### 4.1 `HF_density_operator(C,Orb)` (dev `HF_density_operator.jl:1-21`)

```julia
for Orb_a in Orb
    a, pO, nO = Orb_a.a, Orb_a.pO, Orb_a.nO
    pRho .+= pO * (C.p[:,a] * C.p[:,a]')
    nRho .+= nO * (C.n[:,a] * C.n[:,a]')
end
```
`a_max = length(Orb)`. No (l,j) block mask any more (equivalent when C is block diagonal, which the ordering guarantees).
Fortran: rewrite `HF_density_operator(C,Orb)` (drop the `a_max` argument, real weights; can skip orbitals with
`pO == 0 .and. nO == 0`… but skipping changes nothing numerically). Callers to update: `nj_hf_rpa_export.f90:416`
(RPA group), `nj_hf_mbpt.f90` (rewritten anyway). ~20 lines.

### 4.2 `HF_allocate` (dev `HF_allocate.jl:1-133`)

Formulas, loop bounds and CMS handling identical to base. New pruning:
* `HF_allocate_indices` additionally returns `l_max, j_max` = max l / max 2j over orbitals with `pO > 1e-8 || nO > 1e-8`
  (dev `:189-202`), computed from the **current** `Orb` (occupations of the previous iteration, consistent with `Rho`).
* skip (b,e) if `l_b > l_max || j_b > j_max` (`:36-37`); in the NNN part skip c if `l_c > l_max` or `j_c > j_max` (`:57-65`).
* NNN parity test `P_abc == parity(l_d+l_e+l_f)` was replaced by an assignment `P_abc = parity(l_d + l_e + l_f)` (`:81`).
  Since (a,d) and (b,e) always share (l,j) and l_c == l_f is tested, both parities are always equal → no numerical change.
All skipped terms are exactly zero (ρ vanishes there), so this is performance only. Fortran: add the pruning to
`HF_allocate`/`HF_allocate_indices` (keep port deviation 5: parallel over (a,d)). ~20 lines.

### 4.3 `HF_energy` (dev `HF_energy.jl:1-140`)

Same formulas. Changes: same `l_max/j_max` pruning on a, b, c; J range now
`max(|j_a−j_b|,|j_d−j_e|)/2 … min(j_a+j_b, j_d+j_e)/2` and `J_hat = hat(2J)` (= 2J+1) – terms outside are zero.
The kinetic term `(Rho.p[a,b]+Rho.n[a,b]) * T.p[a,b] * (j_a+1) [* (1-1/A)]` stays inside the b loop (after the N2max test),
unchanged (note: still uses `T.p` for neutrons, harmless since `T.p == T.n`). Old version kept in a comment block.
Fortran: optional pruning, ~15 lines. Results identical.

---------------------------------------------------------------------------------------------------------------------------

## 5. `HF_potential.jl` → new `nj_hf.f90 :: HF_radial_potential` (or new module `nj_hf_potential.f90`)

Base `HF_Radial_Potential.jl` was dead code (not ported). Dev `HF_potential.jl:1-86`, called when `HF.V_HF`.

Algorithm (`HF_radial_potential(Params,Orb,C,T,h_N)`):
1. `pV = C.p * h_N.p * C.p'` (h in LHO basis), same for n; subtract kinetic:
   `CMS1+2B`: `pV -= T.p[1:a_max,1:a_max]*(1-1/A)`; any CMS other than `CMS2B`: `pV -= T.p[1:a_max,1:a_max]`; `CMS2B`: nothing.
   (T may be larger than a_max – use the leading block.)
2. Grid: `r = collect(range(1e-8, stop = 2.5*1.2*A^(1/3), length = 2^13+1))[1:16:end]` → N_r = 513 points (same nodes as
   `rho_export`).
3. `nu_p = 0.5*m_p*hw/hc^2`, `nu_n = 0.5*m_n*hw/hc^2`, `hc = 197.326980, m_p = 938.272013, m_n = 939.565346`
   (`hw = Params.Calc.hw`). `pR[a,i] = Psi_rad_LHO(r_i, n_a, l_a, nu_p)` (Fortran `nj_functions::Psi_rad_LHO` exists).
4. Partial waves `lj_list = sort(unique([(l_a, j_a)]))` (sorted by l then 2j). For each (l,j), `Block` = slots with that (l,j)
   in index order: `pV_lj = pR[Block,:]' * pV[Block,Block] * pR[Block,:]` (N_r×N_r);
   `pV_tot += (j+1)/(4π) * pV_lj`.
5. Labels `"spdfghiklmnoqrtuvwxyz"[l+1] * string(j) * "/2"` (j = 2j value, e.g. `p3/2`).
6. File `io/<Path>/HF/densities/HF_radial_potential.dat`:
   ```
   # Non-local HF potential V(r,r') = h_HF - T_1b (CMS = <CMS>) ... r, r' in fm, V in MeV fm^-3
   # V_lj(r,r') = sum_ab R_a(r) V_ab R_b(r') per partial wave, V_tot = sum_lj (2j+1)/(4 pi) V_lj
   %-14s %-14s %-16s %-16s  ("r","r'","pV_tot","nV_tot") then " %-16s %-16s" per wave ("pV_<lab>","nV_<lab>") then "\n"
   rows: for i, for k: "%-14.8f %-14.8f %-16.8e %-16.8e" (r_i, r_k, pV_tot[i,k], nV_tot[i,k]) + " %-16.8e %-16.8e" per wave + "\n"
         blank line after each i block
   ```
   (513² = 263 169 data rows.)
Use `cfmt_f`/`cfmt_e`/`cfmt_s` with left justification. ~110 Fortran lines.

---------------------------------------------------------------------------------------------------------------------------

## 6. `HF_export.jl` → `nj_hf.f90 :: HF_summary, HF_SPS_summary, HF_export, sps_table`

* `HF_summary`: file `io/<Path>/HF/HF_summary.dat` (lower-case `summary`). One format change (dev `HF_export.jl:30`):
  `"\nNumber of Iterations = %5d, Precision = %.3e \t MeV\n"` – now ends with `\n`, so an extra empty line appears before
  `E_HF = …` (base had no trailing newline). Fortran currently writes `advance='no'` – change to a normal record.
* `HF_SPS_summary`: same table; `O_a` is now the real occupation `Orb[a].pO` / `Orb[a].nO` printed with `%8.3f` (base printed
  1.0/0.0 from the integer). Counter array still `zeros(Int64,N_max+1,J_max+1)` with `J_max = N2max+1` → port deviation 4
  (size by 2*Nmax+1) still needed. Fortran `sps_table`: replace `merge(1,0,Orb(a)%pO == 1)` by `Orb(a)%pO`; drop the
  `Kind`/EFA generality if desired.
* `HF_export(Params,C,Rho,h,Orb)` (dev `:58-79`) writes to `io/<Path>/HF/bin/`:

| file | content | writer / layout |
|---|---|---|
| `orb.bin` (new) | HF orbitals with final occupations | `orbitals_export(Params,Orb,path)`: for i=1..a_max(Calc.Nmax): `Int64 a, Int64 n, Int64 l, Int64 j, Float64 pO, Float64 nO` (40 bytes/record, little endian stream) |
| `C.bin` (was `Bin/C_HF.bin`) | LHO→HF coefficients C | `O1b_export`: protons then neutrons, for a=1..a_max, b=1..a_max with (l_a,j_a)==(l_b,j_b): Float64 `A[a,b]` – unchanged layout |
| `rho.bin` (was `Bin/Rho_HF.bin`) | ρ in the HF basis `Cᵀ ρ C` (diag = occupations) | `O1b_export` |
| `H1b.bin` (was `Bin/h_HF.bin`) | h_N = diag(SPE) in HF basis | `O1b_export` |

  And from `HF_solver`: `H2b.bin` (was `Bin/V2B_HF.bin`) = residual NN in the HF basis, `O2b_export` layout unchanged
  (pp lower triangles, pn full N_T0² block, nn lower triangles; J = 0..2Nmax+1, P = 1,2). Port deviation 3 (J range
  clipped to N2max+1) still needed. Note on the pn block see §11 (common-group item).

Fortran `HF_export`: new paths + `orbitals_export(Params,Orb,path)` binary variant (shared `nj_orbitals`, §10). ~15 lines.

---------------------------------------------------------------------------------------------------------------------------

## 7. Removed: `HF_EFA.jl` → delete `src/hf/nj_hf_efa.f90`

Delete the module (357 lines), its `Makefile` entry (`SRC_HF`), the `use nj_hf_efa` and the EFA blocks in
`nj_hf_driver.f90` (lines 17, 76-79, 144-151), `hf_efa_params_t`, the `&hf_efa` namelist. `summary_header` and `sps_table`
were shared with EFA; keep them for HF. Remove PORTING_NOTES references to `HF_EFA`.

Also delete the HRF branch of the driver (`orbitals_export(Params,Orb)` HRF table + `O2b_export_HRF(…,'Bin/V2B_HF.dat')`):
dev never writes it. Port deviation 1 (HRF orbitals reading `pE_HF.bin`) becomes obsolete for HF (still relevant only if
some other solver keeps an HRF export – BCS/HFB groups).

---------------------------------------------------------------------------------------------------------------------------

## 8. `HF_MBPT.jl` (2327 lines) → rewrite `src/hf/nj_hf_mbpt.f90`

Base: `Source/HF/HF_MBPT.jl` (1080 lines; ported: `HF_MBPT`, `HF_MBPT_energy_correction`, `HF_MBPT_OBDM`,
`HF_MBPT_OBDM_correction`, `HF_MBPT_reorder`, `HF_MBPT_occupation`; not ported: `HF_MBPT_OBDM_NOB`, `HF_MBPT_TBDM`).
Dev functions:

| dev function (lines) | base counterpart | Fortran now | action |
|---|---|---|---|
| `HF_MBPT(Params,Orb,Orb_NN,h_N,V_NN)` 1-12 | `HF_MBPT(Params,Orb,Orb_NN,V_NN)` | `HF_MBPT` | new signature (h_N), calls orbitals_make → energy → OBDM |
| `HF_MBPT_orbitals_make(Params,Orb)` 14-76 | `orbitals_ph_make` (common) | `nj_orbitals::orbitals_ph_make` | new local routine, §8.1 |
| `HF_MBPT_h1b(a,t,h_N)` 78-85 | – (SPState.E) | – | trivial helper |
| `HF_MBPT3_V2b_Pandya(...)` 87-152 | – | – | new, §8.2 |
| `HF_MBPT_energy` 154-753 | `HF_MBPT_energy_correction` (MBPT(2) only) | `HF_MBPT_energy_correction` | extend with MBPT(3), new summary, §8.3 |
| `HF_MBPT_OBDM` 755-1538 | `HF_MBPT_OBDM` + `HF_MBPT_OBDM_correction` + `HF_MBPT_reorder` | same names | merge; formulas unchanged, NAT step via `OBDM_natural_orbitals`, §8.4 |
| `HF_MBPT_TBCF` 1540-2181 | `HF_MBPT_TBDM` (not ported, different normalisation) | none | new, §8.6 |
| `HF_MBPT_occupation` 2183-2203 | same | same | path only, §8.4 |
| `HF_MBPT_occupation_check` 2205-2221 | – | – | new, §8.5 |
| `HF_MBPT_export` 2223-2327 | – | – | new, §8.7 |

Removed from dev: `HF_MBPT_OBDM_NOB`, `HF_MBPT_TBDM`, `HF_MBPT_reorder` (replaced by Hungarian NAT ordering).
Struct `Orb1B_HF_MBPT{p::pnOrbVector, h::pnOrbVector, N_p::pnInteger, N_h::pnInteger}` (particle and hole `Orb1B` lists);
`Orb_HF_MBPT3` is defined but unused.

Isospin convention inside the file: `t = -1` proton, `t = +1` neutron. Parities `P = rem(l1+l2,2)+1`.
All single-particle energies are `E_a = h_N.p[a,a]` (t=−1) / `h_N.n[a,a]` (t=+1), taken from the in-memory `h_N`
(base read the diagonal of `h_HF.bin` via `orbitals_ph_make` – same numbers).

### 8.1 `HF_MBPT_orbitals_make` (14-76)

For a = 1..a_max in index order: proton **hole** if `abs(pO) > 1e-8`, proton **particle** if `abs(1 - pO) > 1e-8` (same for
neutrons). With fractional occupation an orbital would be both; the occupation check (§8.5) prevents that case from reaching
MBPT. Fortran: local routine producing particle/hole `spstate_t` lists (a,l,j,E from h_N) – reuse `pnsvector_t`; do not use
`orbitals_ph_make` (it reads the old `Bin/h_HF.bin` and compares integer occupations). ~40 lines.

### 8.2 `HF_MBPT3_V2b_Pandya(a,b,c,d,J)` (87-152) – cross-coupled ME used by MBPT(3) ph and TBCF ph

```julia
P_ab = rem(l_a + l_b,2) + 1;  if (rem(l_c + l_d,2) + 1) != P_ab  return 0.0
I_min, I_max = div(max(abs(j_a - j_b), abs(j_c - j_d)),2), div(min(j_a + j_b, j_c + j_d),2)
V = - Σ_I (2I+1) * w6j(j_a,j_b,2I,j_c,j_d,2J) * X_I
  t_a+t_b == -2 : X_I = O2b_pp(a,b,c,d,I,P_ab)
  t_a+t_b == +2 : X_I = O2b_nn(a,b,c,d,I,P_ab)
  (-1,+1,-1,+1) : X_I = O2b_pn(a,b,c,d,I,P_ab)
  (-1,+1,+1,-1) : X_I = O2b_pn(a,b,d,c,I,P_ab) * phase(div(j_c+j_d,2) - I)
  (+1,-1,-1,+1) : X_I = O2b_pn(b,a,c,d,I,P_ab) * phase(div(j_a+j_b,2) - I)
  (+1,-1,+1,-1) : X_I = O2b_pn(b,a,d,c,I,P_ab) * phase(div(j_a+j_b+j_c+j_d,2))
  otherwise 0.0
```
(`w6j` = `f6j` with doubled args; `phase(n) = (-1)^n`.) For |t_a+t_b| = 2 the routine does not check t_c,t_d – all callers
guarantee charge conservation. ~60 lines.

### 8.3 `HF_MBPT_energy` (154-753)

**MBPT(2)** (170-265): identical to base (pp/nn `0.25*(2J+1)*|V|²/E_pqhg`, pn `(2J+1)*|V|²/E_pqhg`,
`E_pqhg = E_h+E_g−E_p−E_q`, J loop over (p,q), hole pair filtered by parity and triangle).

**MBPT(3)**, three sectors, each a sum of separately accumulated (pp, pn, nn) parts:

*pp sector* (270-452), parallel over all ordered particle pairs (p,q) of the given isospins
(`ppOrb_particle_list = [(x,y) for x in pOrb_particle, y in pOrb_particle]`, column-major comprehension order irrelevant
for the sum), holes h,g, particles r,s; parity filters `P(h,g) == P(p,q)`, `P(r,s) == P(p,q)`;
`J ∈ [max(|j_h−j_g|,|j_p−j_q|,|j_r−j_s|)/2, min(j_h+j_g, j_p+j_q, j_r+j_s)/2]`:
```julia
# pp, nn (317-319):
ME = 0.125 * (2J+1) * O2b_pp(h,g,p,q,J,P) * O2b_pp(p,q,r,s,J,P) * O2b_pp(r,s,h,g,J,P) / E_pqhg / E_rshg
# pn (378-380): factor 1.0 instead of 0.125, O2b_pn, p,r protons, q,s neutrons, h proton, g neutron
E_pqhg = E_h+E_g-E_p-E_q ;  E_rshg = E_h+E_g-E_r-E_s
```

*hh sector* (540-716), parallel over particle pairs (p,q), holes h,g,f,e; filters `P(h,g)==P(p,q)`, `P(f,e)==P(p,q)`;
J range as above with (f,e) in place of (r,s):
```julia
ME = 0.125 * (2J+1) * O2b_pp(p,q,h,g,J,P) * O2b_pp(h,g,f,e,J,P) * O2b_pp(f,e,p,q,J,P) / E_pqhg / E_pqfe   # pp, nn
#    factor 1.0 and O2b_pn for pn (p,h,f protons; q,g,e neutrons)
E_pqfe = E_f+E_e-E_p-E_q
```

*ph sector* (454-538): combined lists `Orb_particle = [(p,-1)…; (n,+1)…]`, `Orb_hole` likewise; parallel over all ordered
pairs of `Orb_particle` (p,q); loops h,g ∈ Orb_hole, skip unless `t_p+t_q == t_h+t_g`;
`J ∈ [max(|j_p−j_g|,|j_h−j_q|)/2, min(j_p+j_g, j_h+j_q)/2]`:
```julia
V_pqhg = Pandya(p,q,h,g;J);  if abs(V_pqhg) < 1e-12 continue
E_pqhg = E_h + E_g - E_p - E_q
for (r,t_r) in Orb_particle, (f,t_f) in Orb_hole
    if (t_h+t_r) != (t_q+t_f) || (t_g+t_f) != (t_p+t_r) || !triangle(j_f,j_r,2J)  continue
    V_hfrq = Pandya(h,r,f,q;J);  if abs(V_hfrq) < 1e-12 continue
    V_fgpr = Pandya(f,g,p,r;J);  if abs(V_fgpr) < 1e-12 continue
    E_prgf = E_g + E_f - E_p - E_r
    Sum += (2J+1) * V_pqhg * V_hfrq * V_fgpr / E_pqhg / E_prgf
```
(argument order of Pandya = (a,b,c,d) as written; isospin labels are passed along with each orbital.)

**"Hole-hole fix"** (commit `9f80777`) and its follow-up (`fcb42d9`/`27b274c`): in the earlier dev upload the energies inside the
hh pn and hh nn loops (and later found: pp pn `E_q, E_g, E_s` and the whole pp nn loop) were taken with `t = -1` (proton h_N)
for neutron orbitals. In `f3f7a33` every `HF_MBPT_h1b` call uses the isospin of the orbital it belongs to
(pn: p,h,r/f proton; q,g,s/e neutron). The Fortran port has no MBPT(3) yet, so simply port the final version; make sure the
isospin of each energy follows the orbital list it comes from.

Totals: `dE = dE2 + dE3_pp + dE3_ph + dE3_hh`. Stdout (722-731): `"\tHF-MBPT(2) energy    ...   E^(2)    = " * string(round(dE2, sigdigits=9)) * " MeV"`,
`E^(3)`, `E^(3)_pp`, `E^(3)_ph`, `E^(3)_hh`, `E^(tot)` lines (use `jstr(round_sig(x,9))`).

Summary block appended to `io/<Path>/HF/HF_summary.dat` (734-750), exactly:
```julia
println(Summary, "\nSpherical Hartree-Fock Many-Body Perturbation Theory solution review:")
@printf(Summary, "\nE_0^(2)        = %15.8f \t MeV \t\t ... \t LO HF-MBPT(2) correction to ground-state energy", dE2)
@printf(Summary, "\nE_0^(3)        = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(3) correction to ground-state energy", dE3_pp + dE3_ph + dE3_hh)
@printf(Summary, "\nE_0^(3)_pp     = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(3) pp correction to ground-state energy", dE3_pp)
@printf(Summary, "\nE_0^(3)_ph     = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(3) ph correction to ground-state energy", dE3_ph)
@printf(Summary, "\nE_0^(3)_hh     = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(3) hh correction to ground-state energy", dE3_hh)
println(Summary, "\n" * "_"^124)
@printf(Summary, "\nE_0^(tot)      = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(2+3) correction to ground-state energy\n\n", dE)
@printf(Summary, "\nE_0^(2) / A    = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(2) correction to ground-state energy per nucleon", dE2/A)
@printf(Summary, "\nE_0^(3) / A    = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(3) correction to ground-state energy", dE3/A)   # no "per nucleon" (sic)
@printf(Summary, "\nE_0^(3)_pp / A = ... per nucleon", dE3_pp/A) ; same for _ph, _hh
println(Summary, "\n" * "_"^124)
@printf(Summary, "\nE_0^(tot) / A  = %15.8f \t MeV \t\t ... \t NLO HF-MBPT(2+3) correction to ground-state energy per nucleon", dE/A)
```
(the underscore line is a literal of 124 `_` in the source; copy it from line 741.) Replaces the base header
"Spherical Hartree-Fock Leading Order Many-Body Perturbation Theory solution review:" with its two lines.

Fortran: extend `HF_MBPT_energy_correction` → `HF_MBPT_energy`; MBPT(3) loops: pp ×3, hh ×3, ph ×1 with OpenMP
`reduction(+:Sum)` over the outer pair index (collapse the (p,q) pair list into one loop index). ~450 lines.
Cost: pp/hh are O(N_p² N_h² N_{p|h}² J), ph O(N_p² N_h³ N_p J) with Pandya 6j sums – expect it to dominate HF-MBPT runtime.

### 8.4 `HF_MBPT_OBDM` (755-1538)

dρ loops (1p1h pp/pn, 2p2h pp/pn for protons; same for neutrons): **formulas, filters, index orders identical to base**
(verified term by term; only refactoring into `J_hat`, `E_hp`, … temporaries; energies from `h_N`). The existing Fortran loops
in `nj_hf_mbpt.f90 :: HF_MBPT_OBDM` can be kept; feed them the new particle/hole lists.
Thread safety in dev: every `@threads` loop writes rows/columns owned by its loop index (particle/hole sets disjoint) – no race;
the Fortran `reduction(+:dpRho)` is fine.

Tail (1497-1537), replaces base `HF_MBPT_OBDM_correction` + `HF_MBPT_reorder`:
```julia
C      = O1b_import(Params,Orb,"io/<Path>/HF/bin/C.bin")
Rho_HF = O1b_import(Params,Orb,"io/<Path>/HF/bin/rho.bin")      # HF-basis density (base recomputed it from C)
pRho, nRho = Rho_HF.p .+ dpRho, Rho_HF.n .+ dnRho
Orb_NAT, D, Rho_NAT = OBDM_natural_orbitals(Params,Orb,O1B(pRho,nRho))          # common/orb.jl:552
rho_export(Params,Orb,"io/<Path>/HF/HF_summary.dat","io/<Path>/HF/densities/HF_MBPT_radial_densities.dat",
           O1B(C.p*D.p*Rho_NAT.p*D.p'*C.p', …), O1B(C.p*D.p, C.n*D.n))
HF_MBPT_occupation(Params,Rho_HF,Rho_NAT)                       # io/<Path>/HF/densities/HF_MBPT_occupation.dat
if HF.MBPT.Export
    C_NAT = O1B(C.p*D.p, C.n*D.n)                               # LHO -> NAT
    orbitals_export(Params,Orb_NAT,"io/<Path>/HF/bin/orb_temp.bin")
    O1b_export(Params,Orb_NAT,C_NAT,"io/<Path>/HF/bin/C_temp.bin")
    O1b_export(Params,Orb_NAT,Rho_NAT,"io/<Path>/HF/bin/rho_temp.bin")
```
`OBDM_natural_orbitals` (common group; Fortran will need it in a shared module): symmetrise, zero off-(l,j)-block entries,
add `1e-10*a` to the diagonal, `eigen(Symmetric)`, **Hungarian** assignment on `Cost = 1 .- abs.(D)` (row = HF slot), reorder
columns and eigenvalues, flip column a if `D[a,a] < 0`, `Rho_NAT = diagm(occ)`; `Orb_NAT[a]` gets
`l,j = round(Σ_b D[b,a]² l_b)`, `round(Σ_b D[b,a]² j_b)` (from the proton D), `n` = number of b<a with that (l,j), occupations
`pOcc[a], nOcc[a]`. Base used max-|overlap| reordering without phase fixing (`HF_MBPT_reorder`) → occupations file and NAT
densities can differ in ordering/sign only; delete `HF_MBPT_reorder`. Fortran `nj_linalg::hungarian` already exists.

`HF_MBPT_occupation`: same format as base (`%-6s %-12s ×6` header `a pN_HF pN_PT d_pN nN_HF nN_PT d_nN`,
rows `%-6d %-12.6f ×6`), new path `io/<Path>/HF/densities/HF_MBPT_occupation.dat` (was `HF/Densities/HFMBPT_Occupation.dat`).
Radial densities file renamed `HFMBPT_Radial_Densities.dat` → `HF_MBPT_radial_densities.dat`.

The temp files exist only to pass data to `HF_MBPT_export` and are deleted there; the Fortran port may keep `Orb_NAT`,
`C_NAT`, `Rho_NAT` in memory instead (no observable difference unless the run aborts in between) – see open questions.
Fortran: ~80 lines changed (tail rewrite, delete reorder).

### 8.5 `HF_MBPT_occupation_check(Params,Rho_HF)` (2205-2221)

Input: HF-basis density `Cᵀ ρ C`. Returns false (and prints
`"\nHF-MBPT calculation will not be performed due to non-idempotent HF solution ..."`) if for any a
`abs(pRho_a^2 - pRho_a) > 1e-5 || abs(nRho_a^2 - nRho_a) > 1e-5` (diagonal only). Otherwise true. ~20 lines.

### 8.6 `HF_MBPT_TBCF(Params,C,h_N,V_NN,Orb,Orb_HF_MBPT,Orb_NN;Report=true)` (1540-2181) – new

Builds the two-body cumulant σ (MBPT 1st+2nd-order amplitudes) in the **HF basis** in O2B storage, then transforms it to the
**LHO basis**: `O2b_transformation(Params,Orb,Orb_NN,Sigma,O1b_transpose(C),Report=Report)` (2174). `V_NN` is the HF-basis
residual interaction (`V_NN_res`). Normalisation comment (1613-1614): σ_pphh = t_pphh, "1/4 V·Σ = 2(E2+E3)".

Structure:
* `Sigma = O2b_initialize(Params,Orb;Make_Orb_NN=false)` – zero O2B with channels J = 0..N2max+1, P=1,2, packed lower
  triangles (`Ind = Bra + (Ket-1)*N - Ket*(Ket-1)/2`, Ket ≤ Bra) for pp, pn and nn.
* `particle(O) = O <= 1e-3`, `hole(O) = O >= 0.999` on `Orb[].pO/nO`.
* For each JP of `JP_initialize(N2max+1)` (J outer 0..J_max, P inner 1..2), four passes (2nd order, 3rd pp, 3rd hh, 3rd ph),
  each parallel over the list of (Bra,Ket,T) with Ket ≤ Bra: T=−1 over `N_T1 = Orb_NN.N[2,P,J+1]` with pairs
  `Orb_NN.Ind[2,P,J+1]` (a ≥ b), T=0 over `N_T0 = Orb_NN.N[1,P,J+1]` with `Ind[1,…]` (a proton, b neutron), T=+1 over N_T1.
  Each (Bra,Ket,T) element is written by one task only (no race).
* In each pass, with (a,b)=Bra pair, (c,d)=Ket pair and isospins (−1,−1,−1,−1) / (−1,+1,−1,+1) / (+1,+1,+1,+1):
  `if particle(O_a)&&particle(O_b)&&hole(O_c)&&hole(O_d)  X(a,b,c,d)  elseif hole(a,b)&&particle(c,d)  X(c,d,a,b)`; add to
  `Sigma.{pp,pn,nn}[P,J+1][Ind]`.
  - 2nd order (1615-1634): `σ = - V_{abcd}(J) / (E_a+E_b-E_c-E_d)` (pp/pn/nn ME by t_a+t_b; no 1/2 factors in f3f7a33).
  - 3rd pp (1636-1697): `σ += f * V_{ab,rs} V_{rs,cd} / (E_a+E_b-E_c-E_d) / (E_r+E_s-E_c-E_d)`, r,s particles of the channel
    isospins (pn: r proton, s neutron), filter `triangle(j_r,j_s,2J)` and `P(r,s) == P`; `f = 0.5` (pp, nn), `1.0` (pn).
  - 3rd hh (1699-1763): same with r,s holes and second denominator `(E_a+E_b-E_r-E_s)`.
  - 3rd ph (1765-1825):
    ```julia
    direct(a,b,c,d;J) = - Σ_{I=max(|j_a-j_d|,|j_b-j_c|)/2}^{min(j_a+j_d,j_b+j_c)/2} (2I+1) w6j(j_a,j_b,2J,j_c,j_d,2I)/(E_a+E_b-E_c-E_d)
                         × Σ_{(r,t_r)∈particles, (f,t_f)∈holes} Pandya(f,b,c,r;I) * Pandya(a,r,f,d;I) / (E_a+E_r-E_d-E_f)
        # filter: rem(l_f+l_b,2)==rem(l_c+l_r,2), rem(l_a+l_r,2)==rem(l_f+l_d,2), t_b+t_f==t_c+t_r, t_a+t_r==t_d+t_f, triangle(j_r,j_f,2I)
    Phase_ab = (-1)^((j_a+j_b)/2 - J) * (t_a == t_b ? -1 : 1);  Phase_cd likewise
    σ_ph = direct(a,b,c,d) + direct(b,a,c,d)*Phase_ab + direct(a,b,d,c)*Phase_cd + direct(b,a,d,c)*Phase_ab*Phase_cd
    ```
    (in this pass the nn channel correctly uses `nO`, line 2149.)
* Only the stored lower triangle is filled; σ is symmetric by construction (swapped call for hh→pp case).

Fortran: new subroutine `HF_MBPT_TBCF(Params,C,h_N,V_NN,Orb,ph lists,Orb_NN,Sigma,Report)` returning an `o2b_t` in the LHO
basis. Needs an `O2b_initialize`-equivalent allocator in nj_o2b (shared) and `O2b_transformation` (exists, in-place) with
`transpose(C)`. ~400 lines. See §11 for the `pO`-for-neutrons bug.

### 8.7 `HF_MBPT_export(Params,Orb,Orb_NN,Orb_NNN,C,T,h_N,V_NN,V_NN_res,V_NNN)` (2223-2327) – new

Runs only if the occupation check passed and `HF.MBPT.Export`:
1. `Orb_HF_MBPT = HF_MBPT_orbitals_make(Params,Orb)`.
2. `rm -r io/<Path>/HF-MBPT` if present; `mkdir io/<Path>/HF-MBPT`, `io/<Path>/HF-MBPT/bin`.
3. Read `Orb_NAT = orbitals_import(".../HF/bin/orb_temp.bin")`, `C_NAT = O1b_import(Params,Orb,".../C_temp.bin")`,
   `Rho_NAT = O1b_import(Params,Orb,".../rho_temp.bin")`; `Rho_NAT = O1b_transformation(Rho_NAT,O1b_transpose(C_NAT))`
   → LHO basis (`O1b_transformation(O,C) = Cᵀ O C`, so this is `C_NAT ρ C_NATᵀ`).
4. If `TBCF`: `Sigma_NAT = HF_MBPT_TBCF(Params,C,h_N,V_NN_res,Orb,Orb_HF_MBPT,Orb_NN,Report=false)` (LHO basis).
5. Hamiltonian in NAT basis:
   * `ScH && TBCF`: `h, V = H1b_H2b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho_NAT(LHO),Sigma_NAT,Report=false)`,
     `h = C_NATᵀ h C_NAT`, `V = O2b_transformation(…,V,C_NAT)`.
   * `ScH && !TBCF`: `h = H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho_NAT)`, `V = H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho_NAT,Report=false)`, transform both with C_NAT.
   * `!ScH`: `h = H1b(…,Rho_NAT)` (note: still built from the **correlated** density, see §11) transformed with C_NAT;
     `V = O2b_transformation(…, deepcopy(V_NN_res), O1B(C.p'*C_NAT.p, C.n'*C_NAT.n))` (HF → NAT).
   `V_NN` here must be the **bare** NN interaction (§3.1).
6. `Rho_NAT = O1b_transformation(Rho_NAT,C_NAT)` (back to NAT basis, diagonal); if TBCF `Sigma_NAT = O2b_transformation(…,Sigma_NAT,C_NAT)`.
7. Delete the three temp files.
8. Export to `io/<Path>/HF-MBPT/bin/`:
   `orb.bin` (`orbitals_export`, Orb_NAT with fractional NAT occupations), `C.bin` (C_NAT, LHO→NAT), `H1b.bin` (h, full
   (l,j)-block matrix, not diagonal), `rho.bin` (diag NAT occupations), `H2b.bin` (V), and if TBCF `sigma.bin` (σ in NAT basis);
   O1B files written with `O1b_export(Params,Orb_NAT,…)` (same (l,j) block structure as Orb), O2B with `O2b_export`.
Consumers: dev TDA/RPA read `io/<Path>/<Reference>/bin/{orb,C,H1b,H2b,rho[,sigma]}.bin` with `Reference ∈ {"HF","HF-MBPT"}`
(`RPA_solver.jl:7-40`, `RPA_Hamiltonian.jl:37-41`, `RPA_rho.jl:575-577`).

Fortran: ~130 lines + dependencies `H1b`, `H2b`, `H1b_H2b` (dev `common/Hamiltonian.jl`, shared with RPA – must exist in a
shared Fortran module before this can be ported).

---------------------------------------------------------------------------------------------------------------------------

## 9. Files written by HF in dev (complete list)

| path (under `io/<Path>/`) | written by | when |
|---|---|---|
| `HF/HF_summary.dat` (append) | HF_summary, rho_export (radii), HF_SPS_summary, HF_MBPT_energy, rho_export (MBPT radii) – in this order | always / BMF |
| `HF/densities/HF_radial_densities.dat` | rho_export | always |
| `HF/densities/HF_radial_potential.dat` | HF_radial_potential | `V_HF` |
| `HF/bin/orb.bin`, `C.bin`, `rho.bin`, `H1b.bin` | HF_export | always |
| `HF/bin/H2b.bin` | O2b_export (residual NN, HF basis) | `BMF` |
| `HF/densities/HF_MBPT_radial_densities.dat` | HF_MBPT_OBDM → rho_export | `BMF` && check |
| `HF/densities/HF_MBPT_occupation.dat` | HF_MBPT_occupation | `BMF` && check |
| `HF/bin/{orb,C,rho}_temp.bin` | HF_MBPT_OBDM (deleted by HF_MBPT_export) | `BMF` && check && `Export` |
| `HF-MBPT/bin/{orb,C,H1b,rho,H2b}.bin` | HF_MBPT_export | `BMF` && check && `Export` |
| `HF-MBPT/bin/sigma.bin` | HF_MBPT_export | … && `TBCF` |

Removed: `Bin/C_HF.bin`, `Bin/Rho_HF.bin`, `Bin/h_HF.bin`, `Bin/V2B_HF.bin`, `Bin/V2B_HF.dat`, HRF orbital table,
`HF/HF_Summary.dat`, `HF/Densities/*`, `HF/Densities/HFMBPT_*`, `HF-EFA/*`.

---------------------------------------------------------------------------------------------------------------------------

## 10. Interface changes needed in SHARED Fortran modules

1. `nj_types :: orb1b_t` – `pO, nO` become `real(dp)` (dev `Orb1B.pO/nO::Float64`). Every `Orb(a)%pO == 1` / `== 0` test in the
   code base must be rewritten (HF: `HF_density_operator`, `sps_table`; others: `nj_orbitals::orbitals_ph_make`, RPA/EMPM/NCSM/qp).
2. `nj_params :: io_path` – prefix `'IO/'` → `'io/'` (all groups).
3. `nj_params` – `hf_params_t` (+`V_HF`, `ESH`, `MBPT`), new `hf_mbpt_params_t`, delete `hf_efa_params_t` (§1).
4. `nj_input` – `&hf` gets `V_HF, ESH`; new `&hf_mbpt (Export, TBCF, ScH)`; delete `&hf_efa` (§1).
5. `nj_orbitals` – binary `orbitals_export(Params,Orb,path)` and `orbitals_import(Params,path)` (record: 4×Int64 a,n,l,j +
   2×Float64 pO,nO; a_max(Calc.Nmax) records). The existing HRF `orbitals_export(Params,Orb)` is no longer called by HF.
   `orbitals_make(Params)` must produce real (possibly fractional) occupations with dev's filling order (common group).
   `read_h_HF_diag`/`orbitals_ph_make` read `Bin/h_HF.bin` → must read `HF/bin/H1b.bin` (or the `<Reference>/bin/H1b.bin`) if
   other modules keep using them; HF-MBPT will not use them any more.
6. `nj_orbitals` (or `nj_obdm`) – new `OBDM_natural_orbitals(Params,Orb,Rho,Orb_NAT,D,Rho_NAT)` (§8.4), uses
   `nj_linalg::hungarian` (exists). Shared with RPA (dev `RPA_rho.jl`).
7. `nj_obdm :: OBDM_export` – dev name `rho_export`; HF only needs the path arguments (no change if the common group keeps
   the Fortran name).
8. `nj_v2b :: V2b_residual_no2b` – dev returns a new object and leaves `V_NN` untouched (deepcopy). Either change to
   `V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)` with `V_NN` `intent(in)`, or have the HF driver copy.
   PORTING_NOTES "In-place mutation" paragraph must be updated (HF_EFA no longer relies on mutation).
9. `nj_v2b :: V2b_read`, `nj_v3b :: V3b_no2b_read` – dev signatures without `Orb` (common group); HF driver call sites change.
10. `nj_o2b` – a zero-initialised O2B allocator (`O2b_initialize(Params,Orb)` equivalent: channels J=0..N2max+1, P=1..2,
    packed lengths N(N+1)/2 for pp, pn and nn) for TBCF; `O2b_transformation` is used with `transpose(C)` and with
    `C_HFᵀ C_NAT` – existing `intent(inout)` semantics fit.
11. `nj_o1b` – helpers `O1b_transformation(O,C) = (Cᵀ O C)` and transpose (trivial; can be inline `mtmm`/`mmt`).
12. New shared Hamiltonian module (dev `common/Hamiltonian.jl`): `H1b(…,Rho)`, `H1b(…,Rho,Sigma)`, `H2b(…,Rho;Report)`,
    `H1b_H2b(…,Rho,Sigma;Report)` – needed by `HF_MBPT_export`; owned by whichever group ports Hamiltonian.jl (RPA uses it too).
13. `nj_hf` public API used by other groups: `HF_solve` (Orb becomes `intent(inout)` – BCS caller), `HF_density_operator`
    (signature `(C,Orb)`, real weights – `nj_hf_rpa_export.f90:416` caller), `HF_orbital_ordering` (keep public for BCS/HFB).
14. `nj_utils` – `remove_file` (delete temp files) if the temp-file route is kept.
15. Fortran-only modules that read HF output (`src/empm/*`, `src/ncsm/nj_ncsm.f90:741`, `rpa/nj_hf_rpa_external.f90`,
    `rpa/nj_hf_rrpa_driver.f90`) open `Bin/C_HF.bin`, `Bin/V2B_HF.bin`, `Bin/h_HF.bin` → must move to `HF/bin/C.bin`,
    `HF/bin/H2b.bin`, `HF/bin/H1b.bin` (and should take occupations from `HF/bin/orb.bin`, since HF occupations can now differ
    from the naive LHO filling).

---------------------------------------------------------------------------------------------------------------------------

## 11. Julia bugs / races in dev and status of the port notes

Bugs noticed in dev (`f3f7a33`):
1. **TBCF neutron channel uses proton occupations**: in the 2nd-order, 3rd-pp and 3rd-hh passes the `T == 1` (nn) branch sets
   `O_a, O_b = Orb[a].pO, Orb[b].pO; O_c, O_d = Orb[c].pO, Orb[d].pO` (`HF_MBPT.jl:1905-1906, 1995-1996, 2085-2086`);
   only the ph pass uses `nO` (2149). Wrong particle/hole classification of neutron pairs whenever proton and neutron
   occupation patterns differ (N ≠ Z, e.g. 48Ca, 208Pb). Recommend fixing in the port as a documented deviation (owner decision).
2. `HF_orbitals_make` edge cases (§3.4): fallback slot 1 if nothing eligible; spurious ESH throw when a cheaper-than-a_max
   candidate exists but has higher SPE than `a_max`. Unreachable in normal runs; port verbatim or with a `fatal`.
3. `HF_orbitals_make` ordering failure leaves `Vector{Int}(undef)` indices (as base) – keep the Fortran `fatal` (deviation).
4. `HF_MBPT_export` with `ScH = false`: comment says "take HF Hamiltonian", but `h_N_NAT = H1b(…,Rho_NAT)` is built with the
   correlated density, while V is the HF residual interaction. Possibly intended `h = C_NATᵀ (C h_N Cᵀ) C_NAT`. Owner decision;
   port verbatim by default.
5. `O2b_export` (common) writes the pn block as an N_T0² loop with the **packed lower-triangle index also for Ket > Bra**, which
   aliases other elements instead of the transpose (dev `common/O2b.jl:216-225`, same in base). The Fortran `O2b_export`
   writes the symmetric element (`packed_index(max,min)`) – a deviation noted in the repo CLAUDE.md but not in
   PORTING_NOTES (Julia cannot read a port-written pn block); affects `HF/bin/H2b.bin`,
   `HF-MBPT/bin/{H2b,sigma}.bin` only if a consumer reads the upper triangle. Common group should decide and document.
6. `T1b` (common) loops over `a_max(Int.Nmax)` but indexes `Orb` of size `a_max(Calc.Nmax)` → out of bounds when
   `Int.Nmax > Calc.Nmax` (HF passes the Calc-sized Orb). Common group.
7. Divisions by zero in MBPT energies would occur for orbitals that are both particle and hole (fractional occupation); prevented
   by `HF_MBPT_occupation_check`. No action.
8. `HF.HRF` is accepted but ignored.
9. Summary text `E_0^(3) / A … correction to ground-state energy` lacks "per nucleon" (cosmetic, keep for parity).

Races: none found in HF/HF-MBPT dev code. `HF_allocate`/`HF_energy` use per-thread accumulators (`:static`), MBPT energy
likewise; OBDM and TBCF loops write disjoint elements.

PORTING_NOTES status for this group:
* Deviation 1 (HRF orbitals reads `pE_HF.bin`): **obsolete** for HF (HRF export removed).
* Deviation 3 (`O2b_export` J range clipped): **still needed** (dev still loops J = 0..2Nmax+1).
* Deviation 4 (SPS counters sized by 2Nmax+1): **still needed** (dev still `J_max = N2max+1`).
* Deviation 5 (`HF_allocate` parallel over (a,d)): keep.
* "In-place mutation" note: `V2b_residual_no2b` no longer mutates in dev → update; `HF_orbitals_make` now relies on in-place
  mutation of `Orb` → add.
* "Not ported" list: drop `HF_MBPT_OBDM_NOB`, `HF_MBPT_TBDM`, `HF_Radial_Potential` (gone/replaced); `HF_radial_potential` and
  `HF_MBPT_TBCF` are now reachable and must be ported.
* Remove the whole HF-EFA mention.

---------------------------------------------------------------------------------------------------------------------------

## 12. Work estimate (Fortran lines, new or rewritten)

| item | lines |
|---|---|
| params/input (`hf_params_t`, `&hf`, `&hf_mbpt`, delete EFA) | 40 |
| `nj_hf_driver` (dirs, sequence, copy of V_NN, MBPT gating, V_HF) | 80 |
| `nj_hf`: HF_solve, HF_orbitals_make/occupations, density operator, pruning, summary/SPS/export paths, orb.bin | 180 |
| `HF_radial_potential` | 110 |
| `nj_hf_mbpt`: orbitals_make, h1b, Pandya | 110 |
| MBPT(3) energy + summary | 480 |
| OBDM tail (NAT via OBDM_natural_orbitals, paths, temp/in-memory hand-over), occupation check | 110 |
| TBCF | 400 |
| HF_MBPT_export | 130 |
| delete `nj_hf_efa.f90` | −357 |
| **total new/changed (excluding shared Hamiltonian/NAT modules)** | **≈1650** |
