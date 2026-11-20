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
`HF_orbital_ordering`) are ported with `intent(inout)` so that callers relying on the mutation
keep working. (`V2b_residual_no2b` is non-mutating in dev and in the port since group common B:
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

`HF_MBPT_OBDM_NOB`, `HF_MBPT_TBDM`, `HFB_diagonalize_old`, `HF_Radial_Potential`
(commented out in the original), `radial_phase_test` (removed in dev), `integrate_quadrature`,
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

