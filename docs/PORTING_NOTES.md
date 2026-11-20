# NuJuliet Julia → Fortran porting notes

These notes record the conventions of the port and every intentional deviation
from the Julia original (https://github.com/folpradek/NuJuliet, commit 2595e82).
The port is being brought to parity with NuJuliet dev (commit f3f7a33, branch `dev-port`); the
section "NuJuliet dev port" at the end lists the changes and deviations of that work. Items of the
older sections that dev makes obsolete are marked **(obsolete in dev)**.

## Translation conventions

| Julia | Fortran |
|---|---|
| `X[P,J+1]` channel arrays (`Matrix{Vector{Float64}}`) | `X(P,J)%v` with `J` starting at 0 (`rvec_t`, `rmat_t`, …) |
| `Orb_NN.Ind[t,P,J+1][Bra][k]` | `Orb_NN%Ind(t,P,J)%v(k,Bra)` |
| `Dict{UInt64,Int64}` lookup of 2-body orbitals (`O2b_index`) (dev: dense `Orb_NN.Map[a,b,J+1,T+1]`) | dense table `Orb_NN%Dic(a,b,J,T)` |
| dev `Orb3B.Map[ab][TT,J-J_min+1,c]`, `Orb3B.nlj` (dense 3-body NO2B lookup; base: `Dict`) | `Orb_NNN%Map(ab)%m(TT,J-J_min+1,c)` (Int32), `Orb_NNN%nlj(n,l,j_index)` |
| `eigen(Symmetric(A), sortby=±)` | `eig_sym(A,w,Z[,Descending])` (LAPACK `dsyevr`, stable re-sort) |
| `eigen(A::Matrix)` | `eig_general` (symmetric → `dsyevr`, else `dgeevx('B',…)`, sorted by (real,imag)) |
| `sortperm(x)` / `sortperm(x, rev=true)` | stable `sortperm` / `sortperm_rev` |
| `fCG`, `f6j`, `f9j` (CGcoefficient.jl, doubled arguments) | `nj_wigner` — bit-identical port of the float algorithms |
| `string(x::Float64)` | `jstr(x)` (shortest round-trip, Julia layout) |
| `@printf "%W.Pf"` etc. | `cfmt_f`, `cfmt_e`, `cfmt_d`, `cfmt_s` |
| `Threads.@threads` | OpenMP `parallel do` (reductions instead of thread-local accumulators) |
| binary `write(io, Float64)` | `access='stream'` — files are byte-compatible with Julia's |

### Case-insensitivity

Fortran does not distinguish `N_a`/`n_a`, `A`/`a`, `P`/`p`, `C`/`c`, `J`/`j`.
When the Julia code uses both spellings in one scope the renames are:

* mass number `A` → `A_mass`
* oscillator shell number `N_x = 2n_x + l_x` → `NN_x`
* parity `P` vs. particle index `p` → parity `Par`
* orbital index `c` vs. matrix `C` → orbital `cc`
* the `PhState` parity field is `Par`

### In-place mutation

Julia functions that write into their arguments (e.g. `O2b_transformation`,
`HF_orbital_ordering`, dev `HF_orbitals_make`) are ported with `intent(inout)` so that callers relying on the mutation
keep working. Dev `HF_orbitals_make` ends with `Orb .= Orb_new`: `HF_solve` updates the occupations of the caller's
`Orb` in place, and everything after it in `HF_solver` (HF energy, kinetic energy, radii, SPS table, `orb.bin`,
`O1b_export` blocks, HF-MBPT particle/hole split) uses the final HF occupations; the port's `HF_solve` has `Orb`
`intent(inout)` for this. (`V2b_residual_no2b` is non-mutating in dev and in the port since group common B:
it returns `V_NN_res`; the legacy solvers move it back into their `V_NN` with `O2b_move`.)

## Intentional deviations from the original

1. **(obsolete in dev: the HRF `orbitals_export(Params,Orb)` has no caller in dev; the Fortran
   version is LEGACY and goes away in the cleanup stage)** `orbitals_export` (HRF output) read
   `Bin/pE_HF.bin`/`Bin/nE_HF.bin`, which no solver writes; the port reads the diagonal of
   `Bin/h_HF.bin` instead.
2. `V2b_read`, `CMS = "CMS2B"`: the original indexes the `O1B` struct as `T[a,c]`
   (a runtime error in Julia); the port uses `T%p(a,c)` (kinetic `T.p == T.n`).
3. `O2b_export`/`O2b_import` write/read channels up to `min(2*Nmax+1, N2max+1)` instead of
   `2*Nmax+1`, avoiding an out-of-bounds access when `N2max < 2*Nmax` (still needed in dev).
4. Single-particle summary tables size their `(l, 2j)` counters by `2*Nmax+1`.
5. `HF_allocate` parallelizes over the outer (a,d) loop instead of the inner (b,e) loop.

6. `HF_RPA_OBDM` (and the RRPA counterpart): the original updates the diagonal density
   matrix from several threads without synchronization (`Threads.@threads` over `ph`, writing
   `pRho[a_h,a_h]`), which loses updates nondeterministically. The port accumulates serially.
   Verified: Julia run with `-t 1` agrees with the port to 1e-13 in `RPA_Radial_Densities.dat`,
   the RPA radii and the E1 compression/toroidal strengths (which depend on the RPA charge
   radius); multithreaded Julia output differs by up to 2e-3 there.

7. `HF_RRPA_solve` with `ScV3N = false` tests `Params.Calc.RRPA.OBDM == "Full"`, a field that
   does not exist in `RRPA_Parameters` (a runtime error in Julia); the port performs the
   intended basis transformation of the residual interaction.
8. Routines that the original silences by redirecting stdout to `/dev/null` inside the RRPA
   iteration are silenced through the `Verbose` switch of `nj_utils`. (RRPA is removed in dev; the
   dev `H2b`/`H1b_H2b` have a `Report` keyword, ported as an optional argument.)
9. `HF_RPA_OBDM` keeps phonons with `norm(Y)^2 < 0.5`, which is meant to drop the spurious 1- state.
   That state is built as `Y = -Spur_State/sqrt(2)`, so its norm is 1/2 up to rounding and the
   original includes or excludes it depending on the last bit (O16, Nmax = 12, DN2LOGO394:
   Julia 0.4999999999999999 → included, r_p = 2.987 fm instead of 2.805 fm). The port cuts at
   `0.5 - 1e-8`, so the spurious state is always excluded (physical states have norms ≲ 1e-2).

## Original behaviour kept on purpose (possible bugs worth a look upstream)

* `HF_RRPA_OBDM`, neutron-hole branch (`Type == 4`): the mirrored term uses `Y[pg,mu]` where
  the proton-hole branch uses `Y[ph,mu]`.
* `HF_RRPA_solve`: `E_old` is never updated, so the `E_new < E_old` exit tests compare with 2.0.

* `HF_RPA_plot_spectrum_export` sorts N+1 entries (ground state + N phonons) but writes only
  the first N, dropping the highest level.
* HF-RPA instability messages always print `P = +`.
* `HF_RPA_transitions_export` evaluates `OBDM_rN` with the RPA density already transformed
  back to the LHO basis, while the HF density is passed in the HF basis.

## Not ported (unreachable from the public entry points)

`HF_MBPT_OBDM_NOB`, `HF_MBPT_TBDM` (base only, removed in dev), `HFB_diagonalize_old`, base `HF_Radial_Potential`
(commented out in base; replaced by the dev `HF_radial_potential`, which is ported), `radial_phase_test` (removed in dev), `integrate_quadrature`,
`legendre_polynomial`, `logm`, `lorentzian`, `pretty_summarysize` (dev), `w3j` (dev, no caller),
`MatrixElements/rME.jl` (dev `src/common/rO.jl`, still not included by the module),
`Source/Temp/*` (not included by the module, gone in dev).

## Chain B (BCS, HFB, QTDA, QRPA)

Sources: `src/qp/` (`nj_qpo2b`, `nj_qp_common`, `nj_bcs`, `nj_hfb`, `nj_qex`, `nj_qtda`, `nj_qrpa`).

Structure:

* The O^(40) & O^(22) quasiparticle transformations of `qpO2B_40.jl` / `qpO2B_22.jl` (six nearly
  identical 4-stage routines) are one table-driven engine in `qpO2b`; each table row
  (source term, U/V selectors of the four indices, sign) is transcribed from the Julia code.
* `BCS_energy`/`HFB_energy`, `HF_BCS_BMBPT_energy`/`HFB_BMBPT_energy` and the QTDA/QRPA parts that
  are identical in the original (A/B allocation, spurious-state routines, rM, rB, moments &
  transition tables, plot spectra) are shared routines.

Intentional deviations / simplifications (no effect on results):

1. `HFB_Broyden_update`: the trial densities of the original equal the new densities (all mixing
   terms are commented out), so the Broyden history and its `beta` coefficients never influence
   the result. The port applies only the net effect (block symmetrization + 1e-14 cleaning of
   `HFB_Broyden_reconstruct`).
2. `HF_BCS_solve` re-reads the NN interaction for the reference nucleus (`V_NN_Ref`) without ever
   using it; the port skips this read.
3. `qpO2b_export`/`qpO2b_import` use `min(2*Nmax+1, N2max+1)` channels like `O2b_export`.
4. Output silenced by `redirect_stdout` in the original (degenerate-solution energies, ScBCS
   residual interaction) is printed or suppressed with a quiet flag.
5. Failed (l,j) orderings (`HFB_orbital_ordering`, `HFB_canonical_orbital_ordering`) stop with an
   error instead of using uninitialized/zero permutation indices.

Original behaviour kept on purpose (possible bugs worth a look upstream):

* `HFB_diagonalize_MCA` builds the canonical-basis-diagonal pairing field but diagonalizes the
  second time with the unmodified `Delta`.
* `HFB_solve` degeneracy branch evaluates `E_HFB_1` with the densities `Rho_1, Kappa_1` of the
  initial secant step, not of the stalled solution.
* `HFB_Lipkin_Nogami_allocate` leaves `Delta` unchanged (its correction is commented out).
* `QRPA_diagonalize` stores the eigenvectors of (A−B)(A+B) in a real buffer: complex eigenvectors
  would raise an `InexactError` in Julia; the port stops with an explicit error in that case.
* `QRPA_transitions_export` names the E3 table `QRPAE3.dat` and spells the TRK unit `Mev`.

Not ported (unreachable): `HFB_diagonalize_old`, `HFB_canonical_basis_diagonalize`,
`qpO2b_40_allocate_canonical`, `qpO2b_22_allocate_canonical`, the whole of `qpO2B_31.jl`
(its only caller is commented out), `qpO2b_1111_index`, `qpO2b_1111_pn`.

## Verification

`tests/compare_outputs.py REF NEW` compares output folders: `Bin/*.bin` element-wise
(also modulo sign, since eigenvector signs are an arbitrary LAPACK gauge) and `*.dat`
token-wise with a numeric tolerance.

Dev port (2026-11-20): the reference is the patched Julia worktree
`/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref` (branch `port-ref` = dev f3f7a33 + one local commit
per fix F1..F5, never pushed). `tests/dev_port/` holds the Julia drivers (`julia/<case>.jl`, `run_julia.sh`), the
comparison script `compare.py` for the dev `io/<Path>/<Name>/` layout (it replaces `compare_outputs.py` for the dev
layout) and `README.md`. Fix F5 in the reference: `OBDM_natural_orbitals` derives (l,j) of the neutron natural orbital
from `nC`; since `Orb1B` carries a single (l,j), a proton/neutron mismatch is an `error(...)` there, and the Fortran
port must stop with `fatal()` in the same case (never triggered in the test cases).

## NuJuliet dev port (f3f7a33) ... group "common A" (core types, parameters, orbitals, operators)

Spec: `docs/dev_port/common_a.md`. Stage policy: additive; the old APIs that the not-yet-ported
solvers use stay until the cleanup stage and are marked `! LEGACY (removed in cleanup stage)`.

Changes (follow dev):

* `orb1b_t%pO`, `%nO` are `real(dp)`. `orbitals_make(Params[,Use_Int])` builds `a_max(Calc.Nmax)`
  orbitals (`Use_Int=.true.` = Julia keyword `Int=true`: `Int.Nmax`) with the dev filling
  `Occ = min(1, (A-Z-N_occ)/(2j+1))` (equal filling of the first shell that does not fit). The
  residual ~1e-16 occupations of the orbitals after the last filled one are reproduced bit for bit
  (the product `j_hat*Occ` goes through a `volatile` temporary so that it is not contracted into an
  FMA). The not-yet-ported solvers call `orbitals_make(Params,Use_Int=.true.)` (their old
  `Int.Nmax` length) until they are ported.
* New: `orbitals_update`, binary `orbitals_export(Params,Orb,Path)` / `orbitals_import(Params,Path)`
  (`orb.bin`, 48-byte records Int64 a,n,l,j + Float64 pO,nO), `orbitals_1p1h_make` + `orb2b_1p1h_index`
  (dense `Dic(a,b,J,iT)`, iT = (T+3)/2, P implied by l_a+l_b, 0 = absent as `get(Dic,key,0)`),
  `orb2b_2qp_index` (dev has no caller), `OBDM_natural_orbitals(Params,Orb,Rho,Orb_new,C,Rho_diag)`.
* Types: `ind2qp_t` → `ind2b_t`, `qporb2b_t` → `orb2b_2qp_t` (component `i` → `Ind`, new dense `Dic`),
  new `orb2b_1p1h_t`, `pnorbvector_t`, `orb1b_hf_mbpt_t`, `o2b_pandya_t` (pp, pnnp, pnpn, nn).
* `nj_o1b`: `O1b_transformation` (Cᵀ O C), `O1b_transpose`, `O1b_diagonal`, `O1b_deepcopy`,
  `O1b_zeroes`, `O1b_ones`, `O1b_unit_matrix`.
* `nj_o2b`: `O2b_temp(a,b,c,d,J,P,T,O,Orb_NN_t)` (T = -1 pp, 0 pn, +1 nn), `O2b_Pandya_initialize`,
  accessors `O2b_Pandya_pp/pnpn/pnnp/nn`, `O2b_Pandya_transformation(...,F_NN,Orb_NN_t)` with
  `type(o2b_pandya_t)` (the 5-argument `o2b_temp_t` form, whose `pn` = dev `pnnp`, is LEGACY for the
  qp chain), optional `Report` of `O2b_transformation` (combined with the global `Verbose`).
* `nj_obdm`: `OBDM_*` renamed `rho_*` (the old names stay as generic aliases until cleanup);
  radii "normalized by 1/A^1/3" are `sqrt(R2)/A^(1/3)` (dev fix; was `sqrt(R2/A^(1/3))`).
* `nj_functions`: `J1b_convective` dev fixes (the `Y_{l,l+1}` terms are outside the `l > 0` blocks
  with factor `sqrt(l+1)`; was `sqrt(l-1)` inside), new `phonon_weights` (real & complex).
* `nj_utils`: `parity_of(L)` (Julia `parity`), `triangle(j1,j2,j3)`, `hat_of(j)` (Julia `hat`). The
  Fortran names differ from Julia because `Parity`/`Hat` are local or namelist variables in several
  solvers (a use-associated name cannot be redeclared).
* Parameters/namelists mirror the dev structs: `&hf` (+V_HF, +ESH; HRF accepted, no reader in dev),
  new `&hf_mbpt`, `&tda`, `&pnqtda`, extended `&rpa` (17 fields), `&qtda`/`&qrpa` + Reference,
  `&hfb` without Broy_Qmax/Broy_Bmax and `Broy_Amax = 0.8`. `read_input` has one contained reader per
  group (own locals: `Transitions` is character in `&rpa`, logical in `&ncsm`).
* Output root `io/` (was `IO/`); `io_path(Params,Suffix)`, new `io_dir(Params,Sub)`.

Deviations from dev:

1. **F4** `T1b` sizes the operator by `size(Orb)` (dev `T1b.jl:3` uses `a_max(Int.Nmax)`, out of
   bounds for a `Calc.Nmax` orbital list when `Int.Nmax > Calc.Nmax`).
2. **F5** `OBDM_natural_orbitals`: l and j of the neutron natural orbitals come from the neutron
   eigenvectors `nC` (dev `orb.jl:626-631` uses `pC` for both); as one `orb1b_t` carries one l/j,
   the port (and the patched Julia reference `port-ref`) stops with an error when the proton and
   neutron values differ.
3. `orbitals_make` stops with an error when A-Z or Z exceed the capacity of the basis (dev loops
   forever in its `while` loop).
4. `read_input` refuses `Name == Reference` in `&tda`/`&rpa` (dev refuses it in the solver and would
   otherwise delete the reference folder), `Transitions` other than 'Exact'/'Linear' (dev errors in
   `RPA_external.jl:123`), and empty Name/Reference or ones containing `/`.
5. `hungarian` handles square cost matrices only (dev's in-house version also takes rectangular /
   non-finite ones; every dev caller passes a square finite matrix) and stops otherwise.
6. The legacy groups `&hf_efa`, `&rrpa` and the `&hfb` items `Broy_Qmax`, `Broy_Bmax` are still
   accepted (with a warning) until the cleanup stage, where HF-EFA and HF-RRPA are removed.
7. `O2b_temp` and the Pandya accessors do not check that the pairs exist (dev throws a KeyError);
   `orb2b_1p1h_index`/`orb2b_2qp_index` return 0 for any absent or out-of-range pair.

Dev behaviour kept on purpose (see `docs/dev_port/UPSTREAM_BUGS.md`):

* `orbitals_make` gives the orbitals after the last (partially) filled one `min(1, residual/(2j+1))`
  with a ±1e-16 residual instead of 0 (reproduced for bit parity).
* `orbitals_1p1h_make`: open-shell orbitals with `Tol_par <= O <= Tol_hol` (e.g. 1/3 with the
  defaults 0.1/0.7) are neither particles nor holes and are dropped from the TDA/RPA space.
* `OBDM_natural_orbitals` returns occupations that include the `+1e-10*a` splitting.

## NuJuliet dev port (f3f7a33) ... group "common B" (interactions, transitions, Hamiltonian)

Spec: `docs/dev_port/common_b.md`. Modules `nj_v2b`, `nj_v3b`, `nj_tr1b`, new `nj_hamiltonian`
(counterpart of dev `src/common/Hamiltonian.jl`), `&parse` + actions in the app.

Changes (follow dev):

* `V2b_read(Params,V_NN,Orb_NN)` and `V3b_no2b_read(Params,V_NNN,Orb_NNN)` build their own
  interaction-space orbitals (`orbitals_make(Params,Use_Int=.true.)`); the NN file is read in blocks of
  at most 1e6 Float64. `V2b_residual_no2b(...,V_NN [in],V_NNN,V_NN_res [out])` is non-mutating.
  New `V2b_parse`, `V3b_no2b_parse` (byte-identical to the Julia reference on the bundled case, and
  a parse to the full space reproduces the input files), `V2b_correlation_energy`.
* `orb3b_t`: the hash map is replaced by the dev dense map `Map(ab)%m(T_ab+T_ind, J-J_min+1, c)` and
  `nlj(n,l,j_index)`; `V3b_no2b_lookup` as in dev. The packed V_NNN blocks are stored **row by row**,
  `Ind = Ket + Bra(Bra-1)/2` (dev; base was column by column) because `nj_hamiltonian` indexes the
  blocks directly. `V3b_no2b_index` and the block lengths are 64-bit (`integer(i8)`), so blocks with
  more than 65535 states no longer overflow. All lookups are plain array reads (no type-bound
  procedures in the hot loops, cf. commit 33d40d7).
* `V3b_no2b_read`: only the chunks needed in the calculation space are read (largest first), and
  (c,d,e,f) blocks outside the calculation space are skipped without reading (dev
  `V3b_no2b_read_chunk!`). Port extensions kept: `NNN_File = ''` = no 3N force, file size check.
* `nj_tr1b`: `tr1b_t`/`qptr1b_t` gain `F` & `GT` (`Tr1b_BX_initialize`; F/GT mix the isospin bases in
  `Tr1b_transformation`; `qp_transform_bx` = `qpTr1b_BX_initialize`). Upstream E1 fix: convective
  combination in `rME_nabla_dot_rN_YJL` (signs & `sqrt((l+1)/(2l+1))` of the Bra terms, no `l >= 1`
  conditions), `s = sigma/2` normalization (16 pi) in `rME_nabla_cross_S_dot_rN_YJL`, `0.5*` on the
  convective VC/TC/sTC. Quasiparticle forms: `O^11 = O_kl U U + (-1)^((j_a+j_b)/2+X) O_lk V V` (no c-number
  term for E0), time-odd phase `Phase_T = -1` in the EX 20 term of the current operators. The
  downstream E1 current/toroidal/compression strengths of the legacy solvers change accordingly
  (the dev strength formulas of RPA_rB etc. are ported by the solver groups).
* `nj_hamiltonian`: `H1b` (h_N[Rho]), `H2b` (H_NN[Rho], equal to `V2b_residual_no2b` up to summation
  order, checked to 3e-15), `H1b_H2b` (Sigma-dressed h_N). Checked against the Julia reference
  (bundled interaction, 16O, Calc (3,6,9) and (2,4,5), non-diagonal test Rho, Sigma = V_NN): h_N,
  H_NN and h_N[Rho,Sigma] agree to <= 2e-13 absolute, `V2b_correlation_energy` to 1e-15 relative.
* `&parse` (N_max_Out, N_2max_Out, N_3max_Out, NN_Out_File='io/NN_parsed.bin',
  NNN_Out_File='io/NNN_parsed.bin') with the `&run` actions `V2b_parse`, `V3b_no2b_parse` and
  `PARSE` (both; the NNN part only when `NNN_File` is set).

Deviations from dev (results unchanged):

1. `V2b_read` stops reading after the last chunk `a` with `N_a <= Calc.Nmax` (orbitals are ordered by
   N, nothing later is stored); dev reads the whole file. The reading is serial (dev threads over the
   chunks a; each chunk writes only elements whose largest orbital is a, so the result is the same).
2. `V3b_no2b_read` counts the chunk offsets only up to the last chunk with `N_a <= Calc.Nmax` and stops the
   c loop of a chunk at the last calculation-space orbital (port optimizations; dev counts the whole
   interaction space). The f loops of the count/read/parse run over the orbitals of the (l_c,j_c) class
   only (same element order).
3. `H1b`/`H1b_H2b` thread over the (a,d) pairs / the triples with OpenMP reductions instead of the dev
   `4*nthreads` strided tasks (summation order only). `H1b_H2b` reads Sigma from its packed blocks
   instead of unpacking each (J,P) block into a dense matrix (same values, no N_pn^2 temporaries).
4. Parse: the output file is pre-sized by writing its last byte (Fortran cannot truncate); Julia
   `abspath` is emulated by `abspath_of` (relative to `$PWD`, `.`/`..` resolved, no symlinks); negative
   or missing output truncations are fatal.
5. `V2b_parse`/`V3b_no2b_parse` count errors stop after all chunks (one `fatal`) instead of the first
   failing task.

Still needed: deviation 2 of the original list (`T%p(a,c)` in the CMS2B correction; dev unchanged),
the Lawson term `CM_Lawson` (port extension). Race note (unchanged from base): the NNN chunks (a,b) and
(a,e) both write the elements with d = a, which is benign only if the file is symmetric (as NuHamil
files are); Julia and the port keep the value of the last writer.

Not ported (unreachable in dev): the 9-argument `H1b(...,Rho,Sigma)` with the V_NN-Sigma term
(`H1b_dense`, also commented out in `H1b_H2b`), `E1b`, `E2b`.


## NuJuliet dev port (f3f7a33) ... stage 2, HF (spec `docs/dev_port/hf.md` sections 1-7, 9-12)

Layout: `src/hf/nj_hf.f90` (dev `HF_solve`, `HF_allocate`, `HF_density_operator`, `HF_orbitals_make`, `HF_energy`,
`HF_summary`, `HF_SPS_summary`, `HF_export`, `HF_radial_potential`), `src/hf/nj_hf_driver.f90` (dev `HF`, `HF_solver`;
`&run` action `HF`), `src/hf/nj_hf_mbpt*.f90` (dev HF-MBPT, next section). Outputs as in dev:
`io/<Path>/HF/HF_summary.dat`, `HF/densities/{HF_radial_densities,HF_radial_potential}.dat`,
`HF/bin/{orb,C,rho,H1b,H2b}.bin`.

LEGACY (removed in the cleanup stage): the base-version HF is kept unchanged as `nj_hf_legacy`,
`nj_hf_mbpt_legacy`, `nj_hf_efa` and `nj_hf_legacy_driver` (`&run` action `HF_LEGACY`, output `io/<Path>/Bin/*`,
`HF/HF_Summary.dat`, HF-EFA, HRF export, base MBPT(2)); BCS, HFB and the legacy HF_RPA export use `nj_hf_legacy`.
Action `HF` runs `HF_LEGACY` (with a NOTICE) when the same `&run` list contains a legacy consumer of
`io/<Path>/Bin/*` (`HF_RPA*`, `HF_RRPA*`, `EMPM*`, `NCSM`), so existing inputs keep working until those solvers are
replaced. The legacy chain still sees the stage-1 interim limitation (closed-shell nuclei only).

Deviations:
1. Kept from base: deviations 3 (O2b_export J range), 4 (SPS counters sized by `max(N2max+1, 2Nmax+1)`), 5
   (`HF_allocate` parallel over (a,d); dev threads over (b,e) with `:static` accumulators: summation order only).
2. `HF_orbitals_make`: an HF orbital that cannot be assigned to an LHO (l,j) slot is a `fatal` (Julia leaves
   `Vector{Int}(undef)` indices). The occupation loop stops with `fatal` when no unfilled slot is eligible
   (robustness guard; dev falls back to slot 1, which is already filled, and overwrites its occupation, possible only
   with an invalid A/Z or with ESH when every remaining slot is skipped). The ESH "closed-subshell" `throw` is a
   `fatal` with the same text (and the same trigger, including the spurious case listed in `UPSTREAM_BUGS.md`).
3. Driver: the bare `V_NN` is kept for `HF_MBPT_export` only when `HF.MBPT.Export`; otherwise it is freed right
   after the residual interaction is built (memory; no output change). The dev in-memory hand-over is unchanged
   otherwise.
4. `HF_radial_potential`: the data rows are formatted with one internal write per row and C-like left
   justification (fallback to `cfmt_e` for 3-digit exponents); the text equals Julia's `@printf` except that a
   value within an ulp of a rounding tie may print a different last digit (8 of 263169 rows for 16O; compare.py
   accepts one unit of the last digit).
5. Removed from the dev HF only: HF-EFA (`&hf_efa`, `hf_efa_params_t`, `nj_hf_efa`) and the HRF export stay
   with the LEGACY HF until the cleanup stage (`HF.HRF` is accepted and ignored by the dev HF, as in dev).

Not ported (dead in dev): the commented-out old `HF_energy` loop of `HF_energy.jl`.


## NuJuliet dev port (f3f7a33) ... stage 2, HF-MBPT (spec `docs/dev_port/hf.md` section 8, dev `src/HF/HF_MBPT.jl`)

Layout: `src/hf/nj_hf_mbpt_common.f90` (`HF_MBPT_orbitals_make`, `HF_MBPT_h1b`, `HF_MBPT3_V2b_Pandya`, flat particle/hole
lists with the `h_N` energies, isospin-labelled lists `ispstate_t` = dev `(Orb1B, t)` tuples), `src/hf/nj_hf_mbpt_tbcf.f90`
(`HF_MBPT_TBCF`), `src/hf/nj_hf_mbpt.f90` (`HF_MBPT`, `HF_MBPT_energy` MBPT(2)+(3), `HF_MBPT_OBDM`,
`HF_MBPT_occupation`, `HF_MBPT_occupation_check`, `HF_MBPT_export`). Namelist `&hf_mbpt (Export, TBCF, ScH)`. Files as
dev (spec section 9): MBPT block of `HF/HF_summary.dat` + MBPT radii, `HF/densities/HF_MBPT_{radial_densities,occupation}.dat`,
`HF/bin/{orb,C,rho}_temp.bin` (Export; deleted again by the export), `HF-MBPT/bin/{orb,C,H1b,rho,H2b[,sigma]}.bin`.
The temporary-file hand-over between `HF_MBPT_OBDM` and `HF_MBPT_export` is kept as in dev (the export reads C/Rho back
with the (l,j) blocks of `Orb`, so the in-memory values equal dev's exactly).

Deviations:
1. **F1** (clear bug, fixed, mirrored in the patched reference): dev `HF_MBPT_TBCF` classifies the particles/holes of the
   nn channel (T = +1) in the 2nd-order, 3rd-order pp and 3rd-order hh passes with the **proton** occupations `pO`
   (`HF_MBPT.jl:1905-1906, 1995-1996, 2085-2086`); the port uses `nO` (as dev's own ph pass, line 2149). Matters
   whenever the proton and neutron occupation patterns differ (N != Z); verified with 28O (`mbpt_options` case).
2. `HF_MBPT_TBCF`: dev runs four passes over all (J,P,Bra,Ket,T) adding one diagram each to the zero-initialised
   Sigma; the port evaluates the four diagrams of an element in one sweep and adds them in the same order
   (((s2 + s3pp) + s3hh) + s3ph), which gives identical values. Each element is written by its own (Bra) task.
   Robustness guard: `fatal` if the 2-body channels of `Orb_NN` (V_NN) differ from those of the calculation space
   (dev indexes the `O2b_initialize` Sigma with `Orb_NN` without a check).
3. Threads: the MBPT energy sectors are one OpenMP loop over the flattened ordered (p,q) pair list (dev: `:static`
   threads over the same list) with `reduction(+)`; the OBDM loops reduce into private copies of dRho (dev threads over
   the outer orbital and writes only elements owned by it). Summation order differs from Julia at the 1e-15 level
   only (1 vs 8 threads agree to 1e-13 in every output file).
4. The base OBDM loops of `nj_hf_mbpt_legacy` are reused unchanged (dev formulas, filters and index orders are
   identical to base, spec 8.4); the dev tail (natural orbitals via `OBDM_natural_orbitals` with Hungarian ordering,
   F5) replaces the base `HF_MBPT_OBDM_correction`/`HF_MBPT_reorder`.

Kept faithfully (listed in `UPSTREAM_BUGS.md`): `HF_MBPT_export` with `ScH = false` builds H1b from the correlated
density Rho_NAT; the summary text `E_0^(3) / A ... correction to ground-state energy` lacks "per nucleon".

Not ported (dead in dev): `Orb1B_HF_MBPT3` (defined, never used), the unreachable `return 0.0` of `HF_MBPT_h1b`.

## NuJuliet dev port (f3f7a33) ... stage 3, quasiparticle chain (spec `docs/dev_port/qp.md`; track QP)

Sources: `src/qp/` (`nj_qpo2b`, `nj_bcs`, `nj_hfb`, `nj_qex`, `nj_qtda`, `nj_qrpa`, new `nj_pnqtda`, dispatcher `pnQTDA`).
The section "Chain B" above describes the base port; this section lists what changed with dev.

Legacy handling: no legacy consumer reads the BCS/HFB/QTDA/QRPA outputs (`consumers.md`; EMPM/NCSM use the HF chain
only), so the chain was replaced **in place**; there are no `*_LEGACY` actions for it. The chain no longer uses
`nj_hf_legacy` (BCS uses the dev `HF_solve`/`HF_energy`/`HF_orbitals_make` of `nj_hf`, HFB `HF_orbitals_make`).
Orbitals are the dev `orbitals_make(Params)` (Calc.Nmax) everywhere.

Dev changes ported (spec table A-J): output layout `io/<Path>/{BCS,HFB}/{bin,densities}` (no `Bin/`), lowercase file
names (`qpH1b.bin`, `qpH2b.bin`, `*_summary.dat`, `lambda_2.bin`, ...); `Reference` of `&qtda`/`&qrpa`/`&pnqtda`
selects `io/<Path>/<Reference>/bin` (precondition: C, U, V files resp. qpH1b, qpH2b for pnQTDA, dev error text);
QTDA/QRPA folders `bin, densities, transitions/E0..E3, spectra[, amplitudes]` (`QTDA_plot.dat` vs `QRPA_Plot.dat`
as dev); O^(40) pp/nn direct-term sign flip; new H^22 pn(11,11) (engine `Fch = 4`/`Out = 7` with the 4-block Pandya
`pnpn`, accessor `qpO2b_1111_pn`); `qpO2b_import` file-size check (dev message; expected size from the port's channel
range, port note B.3); BCS & HFB(`Pairing='BCS'`) via `HF_orbitals_make` (target fractional filling); HFB modified
Broyden mixing (`broyden_t`, verbatim algorithm, gamma via `dsysv('U')` like Julia's Bunch-Kaufman `Symmetric \`),
safeguarded secant/bisection for Lambda + `Lambda_H` bookkeeping + closed-shell gap Lambda, (l,j)-block matching with
Hungarian assignment and Loewdin projection (`HFB_orbital_matching`) in all three ordering routines; new solver pnQTDA.

Deviations:
1. **F3** (clear bug, fixed, mirrored in `port-ref`): `HF_BCS_initialize_chemical_potential` uses the mid-gap value
   `Lambda = 0.5*(Eps_hole + Eps_particle)` (dev: `0.5*(Eps_hole - Eps_particle)`).
2. Port notes B.1 (Broyden net effect) and B.5 (fatal on failed ordering) of chain B are **obsolete** (dev mixes for
   real; the matching cannot fail). B.2 (`V_NN_Ref` not read), B.3 (export channel range), B.4 (silenced output) still
   apply.
3. HFB chemical-potential bracket / gap: dev's NaN markers are logical flags (`Has_lo`, `Has_hi`, `Has_gap`), so the
   logic does not depend on NaN semantics under ifx fast floating point. Loewdin `P (P^T P)^(-1/2)` uses `eig_sym` and
   `Z diag(w^-1/2) Z^T` (Julia: `inv(sqrt(Symmetric(S)))`), equal to rounding.
4. `HFB_diagonalize_BCS` calls `HF_orbitals_make` on a local copy of `Orb` (dev mutates HFB's `Orb` occupations, which
   no HFB routine reads afterwards); the ESH fatal is reproduced by `HF_orbitals_make` itself.
5. Robustness guards (fatal instead of out-of-bounds / dev `KeyError`): `qpO2b` and `pnQTDA` with `N2max < 2*Nmax`
   (the direct O2b terms resp. the pn 2qp pairs need 2-body pairs up to 2*Nmax), `qpO2b_1111_pn` with a missing pair.
6. Determinism: the Lipkin-Nogami traces (`HFB_Lipkin_Nogami_Lambda_2`) are summed per orbital and added in orbital
   order instead of an OpenMP `reduction` (the LN field feeds back into the iteration; with the reduction 1 vs 8
   threads differed by ~1e-11 with eigenvector sign flips). Other reductions only affect final energies (1e-15 rel.).
7. pnQTDA: `Orb_pn2qp` reuses `orb2b_2qp_t` with `T = 0` without the (unused) lookup; `pnQTDA_allocate` is one OpenMP
   task per (J,P) block (each block written by one thread). `QTDA/QRPA/pnQTDA_import_binary` (`*_external.jl`) have no
   caller in dev and are not ported.
8. `&qtda/&qrpa/&pnqtda Reference` is validated as a folder name (common A reader) rather than restricted to HFB/BCS;
   a wrong name gives dev's "No precomputed reference" message.

Kept faithfully (`UPSTREAM_BUGS.md`): BCS reference HF with the target fractional occupations; `HFB_diagonalize_BCS`
ESH throw for open shells; QTDA/QRPA plot files drop the highest level (pnQTDA writes all); `BCS`/`HFB` both remove
`io/<Path>`; QRPA complex (A-B)(A+B) eigenvectors stop the run (dev `InexactError`, e.g. 22Ne with the bundled
interaction); MCA `Delta`, degeneracy `Rho_1`, LN `Delta`, `QRPAE3.dat`/`Mev` (chain B list).

Not ported (dead in dev): `HFB_Broyden_BackUp.jl` (not included), `HFB_canonical_basis_diagonalize`,
`qpO2b_31` (caller commented out), the `*_import_binary` readers (no caller).

## NuJuliet dev port (f3f7a33) ... stage 3, TDA & RPA part 1 (spec `docs/dev_port/tda_rpa.md` sections 2-5, 9; track RPA)

New files (`src/rpa/`): `nj_ph_common.f90` (1p1h metric, spurious vector + closed-form Gram-Schmidt, shared
`ph_allocate(Want_B)` = RPA_allocate / TDA_allocate, `ph_allocate_tbcf(Want_B,DEV)` = RPA_allocate_TBCF /
TDA_allocate_TBCF incl. the DEV pf/nf contact terms, shared `ph_rB` with the f3f7a33 isovector E1 fixes, shared
`ph_transitions_export`, summary/dimension/plot helpers), `nj_tda.f90` (TDA_diagonalize, TDA_transition_kernel, TDA_rM,
TDA_export), `nj_rpa.f90` (RPA_diagonalize, RPA_OBDM_oY2, RPA_energy + QBA/dress/rho, RPA_rM incl. `Transitions =
'Exact'` with the RPA Sigma, RPA_sigma), `nj_rpa_export.f90`, `nj_tda_rpa_driver.f90` (actions `TDA`, `RPA`). Shared
change: `nj_linalg::eig_sym_gen` (dsygvd itype 1, 'V', 'U'; the TDA orthogonalized 1- channel). Spec section 9 items
1-7, 9-12 were already in place (stages 1-2).

Deviations:
1. **F2** (clear bug, fixed, mirrored in `port-ref`): `RPA_sigma` pn branch, hole-particle-particle-hole term uses `X`
   in `RPA_sigma_hpph` (dev passes `Y`, `RPA_sigma.jl:208-210`).
2. TDA writes the marker `<Name>/.TDA_output` (dev checks for it before deleting an existing output folder but never
   writes it, so a TDA re-run with the same `Name` was refused). `compare.py` ignores the marker.
3. RPA with `Reference_TBCF = .true.` and no `<Reference>/bin/sigma.bin` stops (`fatal`); dev silently runs without the
   reference Sigma. `Name == Reference` is refused by the `&tda`/`&rpa` readers (dev RPA deletes the reference).
4. RPA spurious 2x2 problem: when `|E_s|` or `|n_s|` is exactly `1e-8` the zero-norm (else) branch is taken (dev leaves
   `E[1]`, `X[:,1]`, `Y[:,1]` undefined).
5. Determinism: every threaded loop writes disjoint elements (A/B/N rows i_ph with their mirrored columns, the W1/W2
   rows, the kernel columns, the Sigma (Bra,Ket,T) list); `RPA_OBDM_oY2` and the QBA / dressing energies accumulate
   per (J,P) block and sum the blocks in JP order (dev does the same with threads). Only `V2b_correlation_energy`
   (core, `reduction(+)`) is order dependent, at the 1e-15 relative level.
6. Stage 3 part 2 items: `RPA DEV = .true.` together with `Reference_TBCF` (the only path reaching `RPA_allocate_SE`)
   stops with "not ported yet"; `DEV = .true.` without a reference Sigma is identical to dev (the DEV switch is then
   inert, as in dev). The pf/nf contact terms of `RPA_allocate_TBCF` are already ported (`ph_allocate_tbcf(DEV)`).
7. Port notes 6 (OBDM race) and 9 (spurious state in the OBDM) apply to the LEGACY `HF_RPA` only; the dev RPA follows
   dev exactly (spurious 1- state with a real `E_s` enters OBDM, energy and Sigma; only the E1 moments skip it).
8. The LEGACY actions `HF_RPA`, `HF_RPA_transition_densities/currents`, `HF_RRPA*` and their files stay unchanged
   (`nj_hf_rpa*.f90` is needed by EMPM until stage 5, `HF_RRPA` until HF_ScRPA in stage 4).

Not ported (dead in dev): the Sigma variant of `RPA_spurious_ini` (`RPA_spurious.jl:40-130`), `RPA_OBDM_sigma`, the
commented-out pf/nf block of `TDA_allocate_TBCF`, the commented-out Sigma_NAT transformation of `RPA_export`.

## NuJuliet dev port (f3f7a33) ... stage 3, TDA & RPA part 2 (spec `docs/dev_port/tda_rpa.md` sections 4.4, 6, 7; track RPA)

New files (`src/rpa/`): `nj_rpa_se.f90` (`RPA_allocate_SE`, the experimental DEV path), `nj_ph_external.f90` (the ten
external actions `TDA_/RPA_transition_densities`, `_transition_operator_densities`, `_transition_analysis`,
`_transition_currents`, `_transitions_recalculate` of `TDA_external.jl` / `RPA_external.jl`, incl.
`RPA_transition_kernel` for `Transitions = 'Linear'` and `'Exact'` with the RPA Sigma). Shared (additive): `nj_input`
namelists `&transition_operator_densities` (`JP, nu_list, File_Name = 'Phonon_Transition_Operator_Densities'`) and
`&transition_analysis` (`JP, nu_list, File_Name = 'Phonon_Transtion_Analysis'` (dev spelling), `Cut = 0.1`),
`run_input_t%Op_Densities`, `%Analysis`, `export_request_t%Cut`. `&transition_currents` `M = 0` is the dev default
`M = -1` (= 0).

Deviations:
1. `RPA_allocate_SE` (reached from `RPA_allocate` and from the dressing energy of `RPA_energy`, both with
   `Mode = "OnShell"` and the default `Regularization = "DSRG"`, `Gamma = 2`, `dE_bin = 0.1`, `Chunk = 32768`):
   algorithmic restructuring with the same result up to round-off. Dev accumulates the binned spectral weights
   `W[i,j,b]` (N_ph^2 x N_bin per channel, 1234-1897 bins on the bundled case) and evaluates a full `Sigma(w_i)` for
   every row (N_ph^3 N_bin); the port contracts the two bin weights of every 2p-2h configuration with
   `f(Delta_b - w_i)` first, `dA = (S + S^T)/2`, `S = (A_12 o G) A_12^T` (`G[i,c] = -g_c(w_i)`), one dgemm per block of
   configurations (plus the backward `B_12` part with `f(Delta_b + w_i)` when particle and hole sets overlap). The
   pole grid, bins, linear splitting, `1/N_c`, the configuration list and the vertices are dev's. The order of
   degenerate configurations (Julia `Dict` iteration + stable sort) only changes the round-off. `expm1` is libm's.
2. Not ported (unreachable: both callers pass `Mode = "OnShell"` and the default regularization): the modes `"QP"`
   (`SE_QP`, root search, des Cloizeaux orthogonalization) and `"Fixed"`, the regularizations `"Imaginary"` and
   `"None"`, `RPA_allocate_SE_Old`. The removed "not ported yet" stop of part 1 makes `DEV = .true.` with
   `Reference_TBCF` run as in dev.
3. External actions, robustness guards: a `JP` outside `0..J_max` or a `nu` outside the stored levels of the channel
   stops with a message (dev: `BoundsError`), a missing input binary stops with its path (dev: `SystemError`);
   `levels/`, `densities/` and `transitions/E0..E3` are created with `mkdir -p` (dev: non-recursive `mkdir`, only when
   `levels/` resp. `transitions/` is missing). The RPA externals stop for `RPA_Orb = .true.` as in dev.
4. `TDA_transition_operator_densities`: dev fills an `a_max x a_max x N_grid` buffer with `Amp(a,b) r^X R_a R_b` for
   `a >= b` and mirrors it with `(-1)^((j_b - j_a)/2)`; the port evaluates `<a||Y_J||b>` directly for every pair as in
   `RPA_transition_operator_densities` (identical, `Amp(b,a) = (-1)^((j_b-j_a)/2) Amp(a,b)`; no buffer). The TDA and RPA
   densities/currents share one code path; they differ only where dev differs (real vs complex gauge, current pair-map
   phase `(-1)^((j_b -+ m_b)/2)` and the `|cg| < 1e-14` skip of RPA, level amplitude `g (-1)^J Gamma` vs
   `-Re(g Gamma)/hat(J)`, analysis columns and criterion).
5. Current maps: the map points are snapped to the `(r, theta)` grid with `round` (ties to even). At `x = 0`,
   `theta = pi` and `pi/dtheta = 127.5` is an exact tie, so the snapped node depends on the last bit of `x` and of
   `atan2`. The port computes the map coordinates with separately rounded products (`volatile`, no FMA contraction)
   and a scalar `atan2` (`!DIR$ NOVECTOR`), which reproduces Julia (with FMA contraction the `x = 0` column of the maps
   was mirrored).
6. `RPA_transitions_recalculate` passes the RPA OBDM in the RPA basis to the shared export (dev passes `C Rho C^T`,
   which `RPA_transitions_export` transforms back). The summary is cut before the first
   `"\nReview of {TDA,RPA} electromagnetic transition moments"` as in dev.
7. The external routines are race-free: kernel columns, grid points, map points and pair maps are each written by one
   thread; results with 1 and 8 threads agree to 1e-12 (`--strict-sign`).

## NuJuliet dev port (f3f7a33) ... stage 4 part 1, HF-ScRPA building blocks (spec `docs/dev_port/scrpa.md` sections 1.3-1.7, 2)

New `src/rpa/nj_rpa_rho.f90` (`RPA_OBDM_orbitals_make`, `RPA_OBDM_oY4`, `RPA_OBDM_transform_amplitudes`,
`HF_ScRPA_OBDM_solve`), `src/rpa/nj_rpa_hamiltonian.f90` (`RPA_Hamiltonian`, `RPA_Hamiltonian_HF_reference`),
`nj_rpa::RPA_allocate` / `RPA_sigma_mix`, `nj_types::orb1b_rpa_obdm_t`, `src/core/nj_blas_threads_{mkl,ref}.f90`.
The core helpers of spec section 2 (`OBDM_natural_orbitals`, `orbitals_import/export/update`, O1b wrappers, `H1b`,
`H2b`, `H1b_H2b`, `V2b_correlation_energy`) and `RPA_energy*` / `RPA_sigma` already came with stages 1 and 3 and are
reused unchanged (no `nj_rpa_energy` / `nj_rpa_sigma` modules: those routines live in `nj_rpa`).

1. `RPA_OBDM_oY4`: dev's caller-owned `pM`/`nM` work arrays are block-local in the port (no caller reads them). Per (J,P)
   block and isospin T, `M = -1/2 X_W^T diag(W) X_W^*` uses only the rows with `W > 0` (the other rows of dev's full
   product are exact zeros), and `G = Y_T (1 + M) Y_T^+` only the rows of isospin T (the element sums never read the
   pn cross block); only `Re G` is kept. Blocks whose X and Y have exactly zero imaginary parts (no complex RPA
   solution in the block) use real `dgemm`, every other block complex `zgemm` as dev (owner decision 7; unit check:
   the complex path reproduces dev's result for complex amplitudes to 5e-18). Effect: round-off only.
2. `RPA_OBDM_transform_amplitudes`: U is block-diagonal in the isospin; `U X` is done per isospin block with `dgemm` on
   Re X and (if non-zero) Im X. Round-off only. Returns new arrays (subroutine with `X_new`, `Y_new` out).
3. B7: inside the OpenMP block loops of (1) and (2) each thread sets its own MKL thread count to 1 with
   `MKL_Set_Num_Threads_Local` (module `nj_blas_threads`, selected by the Makefile; the reference-BLAS variant of
   `COMPILER=gnu` is a no-op) and restores it afterwards. Results do not depend on the thread count (no reductions:
   each OBDM element and each block is written by one task, sums in dev's order).
4. `HF_ScRPA_OBDM_solve` interface: dev returns `(Rho, C_LHO_NAT, Orb_NAT)` and updates `C_HF_NAT` in place; the port
   has `Orb` and `Rho` `intent(inout)` (in: current orbitals / `Rho_ini`; out: same n,l,j with the new occupations /
   the new diagonal OBDM), `C_HF_NAT` `inout`, `C_LHO_NAT` `out`. Constants 1e-8, 150, 100, 1e-3 hard-coded as dev
   (decision 3). `rho.bin` and `C.bin` of the reference are re-read on every call as dev.
5. `RPA_Hamiltonian` (B6): for `ScH && TBCF` dev transforms the caller's Sigma in place NAT -> LHO -> NAT; the port
   transforms a copy (the caller's Sigma is unchanged instead of round-tripped, ~1e-15).
6. `RPA_allocate` (dev `RPA_allocate.jl:1`) is now one routine shared by `RPA_solver` and HF-ScRPA (`ph_allocate`,
   `ph_allocate_tbcf`, `RPA_allocate_SE` for DEV, exactly where dev calls them); `RPA_solver` output is byte-identical
   (rpa_options, `diff -r`).
7. Not ported (dead in dev, decision 8): `RPA_OBDM_sigma`, the 9-argument `H1b(...,Sigma)` V_NN-Sigma variant.
8. Kept on purpose: B3 (inclusive particle/hole tests of `RPA_sigma` on the current NAT occupations vs the strict
   ones of `orbitals_1p1h_make` / `RPA_OBDM_orbitals_make` on the reference), and dev's transpose (not adjoint) in
   `M = -1/2 X^T W X^*` (see UPSTREAM_BUGS, stage 4).

## NuJuliet dev port (f3f7a33) ... stage 4 part 2, HF-ScRPA solver (spec `docs/dev_port/scrpa.md` sections 1.1, 1.2, 2.3)

New `src/rpa/nj_scrpa.f90` (`HF_ScRPA`, `HF_ScRPA_solver`, `scrpa_beta`, `ScRPA_terminate_summary`), `&run` action
`HF_ScRPA` (`nj_dispatch_rpa`), `files.mk` entry after `nj_tda_rpa_driver` (whose `print_header`, `print_threads`,
`make_output_dirs` are now public and shared). The ScRPA fields of `&rpa` (`Tol`, `Imax`, `OBDM_diagonal`, `OBDM_alpha`,
`ScH`, `TBCF`, `TBCF_alpha`, `TBCF_beta`, `DEV`, `RPA_Orb`, ...) existed since stage 3 (dev `RPA_Parameters`). The
legacy `HF_RRPA` action, `&rrpa` and `nj_hf_rrpa*` are untouched (decision 2; removed in the cleanup stage).

1. B2 (decision 5): `Name == Reference` is refused by the `&rpa` reader (since stage 3) and again in `HF_ScRPA` with
   `fatal` (dev deletes the reference folder before reading it).
2. Decision 6 / B5: for `ScH = false` the bare `V2b_read` / `V3b_no2b_read` are skipped (dev reads both and never uses
   them); `Orb_NN` then comes from the import of the reference `H2b.bin` (`O2b_import` builds it with the same
   `O2b_initialize` over the first a_max(Calc.Nmax) orbitals as `V2b_read`, so the table is identical). Results are
   unchanged (scrpa_options `ScRPA-noScH`, `ScRPA-noScH-TBCF` agree with Julia, which reads the files).
3. Dev's unused `Delta_E` / `Delta_E_old` are not ported. The redundant final `RPA_Hamiltonian` for TBCF (B5) is kept
   (one call, same printed log). `Sigma_old` is freed before the final step as dev sets it to `nothing`.
4. Dev `Sigma_AB = Sigma` (the mixed Sigma of A and B) and the re-binding `Sigma = RPA_sigma(...)` of the final step are
   two variables in the port (`Sigma` passed as `Sigma_AB`, `Sigma_fin` for `E_sigma`, `RPA_rM` and `sigma.bin`).
5. The external `RPA_*` actions (`RPA_transition_densities`, `_operator_densities`, `_analysis`, `_currents`,
   `RPA_transitions_recalculate`) work on a ScRPA output folder unchanged (`&rpa Name` = the ScRPA folder), as dev.
6. Observation, not a deviation: when the inner OBDM loop diverges and its unmixed explicit fallback produces non-finite
   values, the following RPA diagonalization stops the Fortran run (`ERROR: dgeevx failed`; seen for 4He with
   `Ortho = false`, `OBDM_alpha = 1`, `TBCF`, `TBCF_beta = 2`). Dev has no guard either (its `eigen` would throw on
   NaN input); this case was not run in Julia.
7. Stage-4 verification (2026-11-20): `O1b_import` and `O2b_import` (core) check the `open` with `iostat` and stop with
   `ERROR: O1b_import / O2b_import: cannot open <path> ...` (exit status 1) like `orbitals_import`, instead of the
   ifx runtime abort (`forrtl: severe (29)` + traceback). Seen for an HF_ScRPA `Reference` that exists but is not an
   HF / HF-MBPT folder (e.g. an RPA folder without `H1b.bin`); dev throws a Julia `SystemError` there. Error path only,
   all outputs unchanged. A missing `Reference` folder keeps dev's message + return (exit status 0, no output folder).
8. Stage-4 review fix (2026-11-21): the V3N.Sigma term of the core `H1b_H2b` (`nj_hamiltonian`, used by HF-MBPT with
   `TBCF`/`ScH` and by the TBCF HF-ScRPA `RPA_Hamiltonian`) no longer sums into `ph, nh` with an OpenMP
   `reduction(+)` + `schedule(dynamic,4)` (run-to-run differences up to 7e-14). The triples of each V_NNN block are
   dealt to 64 fixed interleaved chunks (`abc = k, k+64, ...`), each chunk sums serially into its own accumulator
   (`phk(:,:,k)`, `nhk(:,:,k)`, kept over all blocks; 2 x 64 x a_max^2 doubles, 8.5 MB at 13 shells) and the chunks
   are added in the order k = 1..64 at the end. The result is independent of the thread count and schedule (bitwise
   equal for 1 and 8 threads and between runs); it differs from dev's serial order only by rounding (<= 1e-13 in the
   outputs). Dev sums serially; this is a summation-order change only.
