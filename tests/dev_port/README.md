# dev-port verification harness

Reference: the **patched** Julia NuJuliet dev, worktree
`/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref`, local branch `port-ref`
(= dev `f3f7a33` + one commit per Julia bug fix, never pushed):

| commit | fix |
|---|---|
| F1 | `HF_MBPT_TBCF` nn channel (T==1) classifies particles/holes with the neutron occupations `nO` |
| F2 | `RPA_sigma` pn branch passes `X` to `RPA_sigma_hpph` (as pp/nn do) |
| F3 | `HF_BCS_initialize_chemical_potential`: `Lambda = 0.5*(Eps_hole + Eps_particle)` |
| F4 | `T1b` loops over `length(Orb)` (not `Params.Int.Nmax`) |
| F5 | `OBDM_natural_orbitals`: neutron (l,j) from `nC`; `Orb1B` has one (l,j), so a p/n mismatch is an `error` |

`git -C <ref> log --oneline f3f7a33..port-ref` lists them.

## Layout

```
tests/dev_port/
  julia/common.jl        shared driver helpers (reference include, stdout flusher, step timing + output checks)
  julia/<case>.jl        Julia driver of one case
  run_julia.sh           runs one case: work/julia_t<N>/<case>/{io/,run.log,time.log,timings_julia.txt}
  compare.py             compares two io/ trees
  work/                  run outputs (gitignored)
```

## Cases

| case | interaction | content |
|---|---|---|
| `bundled` | dev `io/NN.bin`, `io/NNN.bin` (hw 16, Nmax 3, N2max 6, N3max 9) | mirrors dev `scripts/sample_script_run.jl` + option variants: `io/O16` HF → TDA → RPA (+ transition densities, currents, operator densities, analysis for 1- ν = 5,6,7), RPA `Transitions='Exact'` (`RPA-Exact`), RPA `DEV=true` (`RPA-DEV`); `io/O16_MBPT` HF with MBPT(Export, TBCF, ScH) → TDA/RPA with `Reference='HF-MBPT'`, `Reference_TBCF=true` (`TDA-MBPT`, `RPA-MBPT`); `io/He4` HF → HF_ScRPA (`HF-ScRPA`), + `TBCF=true` (`HF-ScRPA-TBCF`), + `OBDM_diagonal=true` (`HF-ScRPA-OBDMdiag`); `io/O18_BCS` BCS → QTDA → QRPA; `io/O18_HFB` HFB → QTDA → pnQTDA → QRPA; `io/parse` V2b_parse / V3b_no2b_parse to 2/4/6 |
| `ca40_5sh` | DN2LOGO394 hw 12 (files 16/32/28) | `io/Ca40` 40Ca, Nmax 4, N2max 8, N3max 8: HF + MBPT(Export) → TDA → RPA (+ 1- transition densities ν = 1..3) → HF_ScRPA (`HF-ScRPA`) |
| `o16_open` | bundled | `io/O20` 20O (0d5/2 neutrons 4/6 filled), HF `ESH=false` → TDA → RPA with default tolerances (0d5/2 dropped from the 1p1h space), `TDA-frac`/`RPA-frac` with `Tol_par=0.9`, `Tol_hol=0.5` (0d5/2 both particle and hole) + 2+ transition densities |
| `hf_options` | bundled | HF options of no other case, HF only (`BMF = false`): `io/O16_VHF` 16O `V_HF = true` (radial HF potential, CMS1+2B); `io/O18_ESH` 18O `ESH = true` (valence neutrons in 1s1/2), `V_HF = true`, `CMS = "None"` |
| `mbpt_options` | bundled | HF-MBPT export branches of no other case, HF + MBPT only: `io/O28_TBCF` 28O (N = 20, N != Z, fix F1) `Export`, `TBCF`, `ScH = false`; `io/O16_ScH` 16O `Export`, `ScH`, `TBCF = false` |
| `rpa_options` | bundled | TDA/RPA options (stage 3): `io/O16` 16O HF → TDA, RPA `Ortho = false` (`TDA-noOrtho`, `RPA-noOrtho`), RPA `TBCF`, `RPA_Orb`, `Transitions = 'Exact'` (`RPA-TBCF`: RPA_sigma with fix F2, density energy, natural-orbital export); `io/O16_MBPT` HF + MBPT(Export, TBCF, ScH) → RPA `Reference = 'HF-MBPT'`, `Reference_TBCF`, `TBCF`, `'Exact'` (`RPA-dress`: Sigma_AB dressing energy), RPA `Reference_TBCF`, `DEV = true` (`RPA-DEV-TBCF`: RPA_allocate_SE, stage 3 part 2); part 2 additions: externals on `TDA-noOrtho`/`RPA-noOrtho` (2- analysis `Cut = 0.05`, 2+ densities `File_Name = 'D2'`, 0+/3- operator densities, 2+ currents `M = 1` on 20 x 30, both `*_transitions_recalculate`), `RPA-Exact` + 1- densities/analysis/currents/recalculation, `io/O16_MBPT` `RPA-DEV-dress` (`Reference_TBCF`, `TBCF`, `DEV`, `'Exact'`: SE in A and in the dressing energy) + all RPA externals, `TDA-MBPT` + operator densities/analysis/recalculation, `RPA-DEV-frac` (`DEV`, `Tol_par = 0.9`, `Tol_hol = 0.02`: overlapping particle/hole sets) |
| `qp_options` | bundled | quasiparticle-chain options (stage 3): `io/Ne22_HFB` 22Ne HFB (open p & n shells) → QTDA `Ortho = false` → pnQTDA (QRPA stops in both codes: complex eigenvectors); `io/O20_LN` 20O HFB `LNT`, `LNRes` → QRPA; `io/O20_MCA` HFB `Pairing = 'MCA'` → QTDA; `io/O20_HBCS` HFB `Pairing = 'BCS'` → pnQTDA; `io/O20_ScBCS` BCS `ScBCS`, `A0 = 16`, `Z0 = 8` → QTDA, pnQTDA, QRPA `Ortho = false` (`Reference = 'BCS'`) |
| `scrpa_options` | bundled | HF-ScRPA options (stage 4): `io/He4` (sample tolerances of `bundled`) `ScH = false` (`ScRPA-noScH`: oscillates, non-convergence summary after Imax = 100), `ScH = false` + `TBCF`, `TBCF_alpha = 0.5`, `TBCF_beta = 3` (`ScRPA-noScH-TBCF`), `TBCF`, `TBCF_alpha = 1`, `TBCF_beta = 2`, `OBDM_alpha = 0.8`, `RPA_Orb`, `'Exact'` (`ScRPA-ramp`), `TBCF`, `DEV` (`ScRPA-DEV`), `TBCF`, `Imax = 3` (`ScRPA-noconv`: termination summary), `HF-ScRPA` + RPA externals on the ScRPA folder (1- densities, 2+ operator densities, 1- analysis `Cut = 0.05`, 1- currents `M = 0`, recalculation); `io/O16_MBPT` HF + MBPT(Export, TBCF, ScH) → HF_ScRPA `Reference = 'HF-MBPT'` (default tolerances, `ScRPA-MBPT`); `io/O20` HF `ESH = false` → HF_ScRPA `Tol_par = 0.9`, `Tol_hol = 0.5`, `Imax = 30` (`ScRPA-frac`) |
| `empm_ncsm` | bundled (4He, 16O), DN2LOGO394 hw 12 (40Ca) | stage 5, **no Julia counterpart**: EMPM / NCSM on the dev HF reference vs a legacy baseline of the frozen stage-4 binary (section "Stage-5 EMPM / NCSM checks"): `io/He4c` 4He 3 shells, complete chain EMPM_TDA → 2PH → 3PH → 4PH → FULL (`Full_N = 4`, `Full_Diag = 'both'`, `CM_Exact`, `CM_Variants`, `CM_Beta`), `Ortho = false`, + dev TDA (`Ortho = false`); `io/He4t` 3 shells: TDA, 2PH, FULL (`Full_N = 2`), EMPM_TRANS, EMPM_SF (`'full'`/`'iv'`, `'ind'`/`'is'`), NCSM (`EMPM`/HF npnh with transitions, `EMPM`/HO with `Beta_CM`, `NN`); `io/He4o` 3 shells `Ortho = true`: TDA (+ dev TDA), 2PH, 3PH (`E3max`), FULL (`CM_Keep`, `Full_Svd`), TRANS; `io/O16` 4 shells: full-space EMPM_TDA (+ dev TDA), 2PH (`E2max`, 0+ 1- 2+), NCSM ≤ 1p1h; `io/O16w` window, action `EMPM` (`Full_H`) in one `&run` with `HF`, NCSM `Window`; `io/Ca40` 4 shells: EMPM_TDA (+ dev TDA), 2PH; error paths (fractional occupations, no `H2b.bin`, HF-MBPT reference, missing reference, no EMPM_TDA, bad `Reference`; NCSM: fractional reference with Basis HO and with `Window`, HF-MBPT reference with Hamiltonian `EMPM`; stage 6: removed action `HF_RPA`, removed `&hfb Broy_Qmax`) |

The parameters of every step are in the driver scripts; the Fortran namelist inputs must reproduce them
(`Path` and `Name` give the output directories).

## Julia (ready)

```sh
cd tests/dev_port
./run_julia.sh bundled 8      # -> work/julia_t8/bundled/io
./run_julia.sh bundled 1      # thread-race check (-t 1)
./run_julia.sh ca40_5sh 8
./run_julia.sh o16_open 8
./run_julia.sh hf_options 8   # HF V_HF / ESH (stage 2)
./run_julia.sh mbpt_options 8 # HF-MBPT export branches (stage 2)
./run_julia.sh scrpa_options 8 # HF-ScRPA options (stage 4)
```

`run_julia.sh <case> [NTHREADS=8]` (max 32) wipes `work/julia_t<N>/<case>/`, creates `io/` there (dev writes
`io/<Path>/...` relative to the cwd and needs `io/` to exist) and runs
`julia --project=@v1.10 -t N,1 julia/<case>.jl` with `OPENBLAS_NUM_THREADS=N` under `/usr/bin/time -v`. The extra
interactive thread flushes the log every 10 s. Environment overrides: `JULIA` (default juliaup 1.10.11),
`NUJULIET_REF` (default the `port-ref` worktree). The user's default `@v1.10` environment already has every
package dev needs (CGcoefficient, BenchmarkTools, DelimitedFiles, ...), so no dedicated project is needed.
Each step checks that its output directory exists (dev reports many input errors with `println` + `return`), so
a silently skipped step stops the driver.

Reference wall times (elwe10, AMD EPYC 7532), 2026-11-20:

| case | threads | wall (incl. compilation) | sum of steps | peak RSS |
|---|---|---|---|---|
| bundled | 8 | 3:26 | 90 s | 2.8 GB |
| bundled | 1 | 4:50 | 168 s | 1.1 GB |
| ca40_5sh | 8 | 2:55 | 101 s (HF+MBPT 58 s, HF_ScRPA 34 s) | 2.1 GB |
| o16_open | 8 | 1:25 | 21 s | 1.8 GB |
| scrpa_options | 8 | 2:48 | 96 s (HF_ScRPA steps 1.5-17 s) | 1.4 GB |

## Fortran

Namelists `fortran/<case>/NN_<step>.nml`, one file per Julia `step` group in run order (`*_HF.nml` = the HF steps;
the TDA/RPA/... files are run by the later stages). `run_fortran.sh <case> [NTHREADS=8] [GLOB='*.nml']` runs the files
matching GLOB in order, e.g. `./run_fortran.sh bundled 8 '*_HF.nml'` for the HF steps only (HF + HF-MBPT).
`H2b.bin` in the HF basis and `H1b.bin`, `H2b.bin`, `sigma.bin` in the HF-MBPT natural-orbital basis are SIGN when an
HF / natural orbital has the opposite sign (stage-2 check with the column-sign gauge of `C.bin` applied element-wise:
agreement to <= 2.3e-13; `python3 -I gauge_check.py REF_BIN NEW_BIN NMAX N2MAX` does this check).
Quasiparticle chain (stage 3): `qpH2b.bin`, `QTDA_X.bin`, `QRPA_X/Y.bin`, `pnQTDA_X.bin` are SIGN when a quasiparticle
state (C / U,V column) has the opposite sign; `python3 -I qp_gauge_check.py REF_IO NEW_IO NMAX N2MAX BCS|HFB` (REF_IO,
NEW_IO = `io/<Path>`) derives the quasiparticle signs from C, U, V of the reference solver and compares element-wise.

### Original plan text (stage 1)

Expected pattern, namelist files with the same `Path`/`Name` values as the Julia driver, run from a work
directory that contains `io/`: `cd <workdir> && <nujuliet> <case>.nml`. `run_fortran.sh <case> [NTHREADS=8]`
(added alongside by the solver stage) does this for the namelists in `fortran/<case>/` (run in file order) and
writes `work/fortran_t<N>/<case>/io`. Set `OMP_NUM_THREADS` <= 32, `OMP_STACKSIZE=512M`, `ulimit -s unlimited`.

## Compare

```sh
python3 -I compare.py work/julia_t8/<case>/io work/fortran_t8/<case>/io [--tol 1e-8] [--quiet] [--strict-sign]
```

Status per file: `ok`, `SIGN` (equal only in absolute value), `DIFF`, `SIZE`, `MISSING` (only in REF), `EXTRA`
(only in NEW). Exit status 1 on any DIFF/SIZE/MISSING/EXTRA (and SIGN with `--strict-sign`).

* `*.bin`: little-endian Float64 (ComplexF64 = re/im pairs); `orb.bin` = 4 Int64 (exact) + 2 Float64;
  `*_parsed.bin` = Float32. Tolerance: `|a-b| <= tol*max(1,|a|,|b|)`.
* `*.dat`: token-wise; numeric tokens with tol plus one unit of the last printed digit, text exactly (a number in
  parentheses, e.g. `(relative 0.029)` of the ScRPA termination summary, is split off and compared as a number).
* Gauge freedom handled: per-column eigenvector signs of `C.bin`, `U.bin`+`V.bin` (joint), `TDA_X`,
  `RPA_X`+`RPA_Y` (joint), `QTDA_X`, `QRPA_X`+`QRPA_Y` (joint), `pnQTDA_X` (block sizes from
  `<Name>/spectra/<Solver>.dat`, validated against `*_E.bin`); lower triangle only of the pn blocks of `H2b.bin`
  and `sigma.bin` (Julia's aliased upper triangle); per-column sign of numeric columns of `*Phonon*.dat`.
* Not resolved: the sign of the HF/natural **orbitals** themselves. When LAPACK returns an HF orbital with the
  opposite sign, `H1b/H2b/sigma/qpH2b` in that basis, the X/Y rows (p,h) and the transition-analysis tables
  change sign element-wise; such files are `SIGN`. All gauge-invariant outputs (energies, spectra,
  transition strengths, densities, summaries) are compared strictly.

Self-check: Julia `-t 1` vs `-t 8` on `bundled` (different OpenBLAS thread counts) gives 409 ok, 38 SIGN,
0 DIFF: no thread race is visible at 1e-8; all differences are eigenvector/orbital sign gauges (|.| agree to
≤ 6e-12).

## Stage-4 unit check (HF-ScRPA building blocks)

`unit_s4/run_unit.sh [BUILD=build_s4] [NT=8]`: `unit_s4/scrpa_blocks.jl` (patched reference, called directly, no Julia
source changes) dumps intermediates of 4He (HF of bundled step 07, ScRPA parameters of step 08) to
`work/unit_s4/dump/` (Julia rerun only when missing or `REDO_JULIA=1`); `unit_s4/scrpa_blocks.f90` is linked against the
objects of BUILD, reads the Julia initial RPA amplitudes and checks `RPA_OBDM_oY4` (real & complex path),
`OBDM_natural_orbitals`, `RPA_OBDM_transform_amplitudes`, `HF_ScRPA_OBDM_solve` (default, OBDM_diagonal, alpha = 1),
`RPA_sigma`, `RPA_sigma_mix`, `RPA_Hamiltonian` (ScH, ScH+TBCF, HF reference) and the energies after one ScRPA step.
Prints `CHECK <item> <max abs deviation>` and `UNIT_S4 worst deviation`.

## Stage-4 HF-ScRPA checks

Iteration histories (outer `HF-ScRPA iteration` lines, inner `OBDM iteration` lines, exits) are compared by
`diff` of the `grep -E 'OBDM iteration:|HF-ScRPA iteration|iteration stopped|has not converged - using|NOT converged|finished after|unstable'`
lines of `run.log` of both codes. SIGN files of the ScRPA folders: `RPA_X/Y.bin` (eigenvector columns) and
`sigma.bin` / `H1b.bin` (NAT basis = HF basis x C_HF_NAT, so the HF proton-column signs propagate):
`python3 -I gauge_check.py <julia>/<Name>/bin <fortran>/<Name>/bin NMAX N2MAX sigma.bin rho.bin`.

## Stage-5 EMPM / NCSM checks (no Julia counterpart)

EMPM and NCSM are Fortran-only. Their reference is a **legacy baseline** of the frozen stage-4 binary
(`/scratch/local.ssd/knapp/job_25503048.pbs-m1/nujuliet-snap-d753f29/nujuliet`, commit d753f29: legacy HF writing
`io/<Path>/Bin/{C_HF,h_HF,V2B_HF,Rho_HF}.bin`, EMPM/NCSM reading them), against which the dev-form chain (dev HF writing
`io/<Path>/HF/bin/{orb,C,rho,H1b,H2b}.bin`, EMPM reading `&empm Reference`, default `'HF'`) is compared.

```sh
cd tests/dev_port
./run_empm_ncsm.sh baseline 16      # fortran/empm_ncsm_legacy/*.nml with the frozen binary -> work/empm_baseline (~45 s)
./run_empm_ncsm.sh dev 16           # fortran/empm_ncsm/*.nml with ./nujuliet (NUJULIET=...), then fortran/empm_ncsm/errors/*.nml
./run_empm_ncsm.sh compare 16       # python3 -I empm_compare.py work/empm_baseline/io work/empm_dev_t16/io
python3 -I empm_compare.py --tda work/empm_dev_t16/io     # full-space EMPM_TDA vs the dev TDA action (Reference 'HF')
```

* Each namelist file runs in its own invocation (wall time per file in `work/.../times.txt`). The legacy inputs
  `fortran/empm_ncsm_legacy/` use `HF_LEGACY` explicitly (in the frozen binary a lone `HF` is the dev HF). `HF_LEGACY`
  was removed in stage 6 (the current `./nujuliet` stops on it), so these inputs stay as they are and **need the frozen
  binary**: regenerate the baseline only with it (`NUJULIET_LEGACY` overrides the default, any binary of commit ec14d69
  or earlier works); `TAG=<suffix>` writes a dev run to `work/empm_dev<suffix>_t<N>`.
* Error cases `fortran/empm_ncsm/errors/*.nml` run alone in a fresh `io/`; the first line `! EXPECT: <text>` must appear
  in the log (fatal errors exit 1, "Error! ..." + skip exits 0 as the dev TDA); `! MKDIR: io/...` lines create stale
  directories first (`e10`, `e11`: a left-over `io/O16/EMPM/bin` without the reference). `e12`, `e13` (stage 6): a removed
  action (`HF_RPA`) and the removed `&hfb Broy_Qmax` / `Broy_Bmax` stop the file before anything runs.
* `empm_compare.py` maps the legacy names (`Bin/EMPM_*` → `EMPM/bin/EMPM_*`, `EMPM/Spectra|Transitions` →
  `EMPM/spectra|transitions`, `TDA_Summary.dat` → `TDA_summary.dat`, `Bin/{C_HF,h_HF,V2B_HF,Rho_HF}.bin` →
  `HF/bin/{C,H1b,H2b,rho}.bin`, `HF_Summary.dat` → `HF_summary.dat` (E_HF only)) and prints the max deviation per file:
  `d` (signed) and `d_abs` (|.|), relative to max(1,|a|,|b|). Eigenvector files (`EMPM_2PH_*_C`, `EMPM_FULL*_Z*`, the
  state files `EMPM_S<n>_*`) are compared per column up to a sign and by projectors within degenerate clusters (state
  files: int32 header/index arrays exactly, R by singular values). Status `GAUGE`: the independent rows chosen by the
  pivoted Cholesky among (nearly) dependent rows differ (`Ind` of `Basis_*.dat`, `Kc` of the state files), so C on those
  rows is not comparable; energies (the matching `_E.bin`, signed) and X are compared and must agree, else DIFF.
  Text files: every number is compared beyond its own printed unit (one unit of the last printed digit; integers
  exactly), and the printed `d`/`d_abs` are that measure. `SIGN` (equal only in |.|) is allowed **only** for files
  that carry a sign gauge: `HF/bin/{C,H1b,H2b}.bin`, `EMPM_TDA_X.bin`, `EMPM_F.bin`, `TDA_Structure.dat` (components
  with equal printed |X| are sorted before comparing) and the `M_p`, `M_n` columns of the NCSM `_E<L>.dat` tables (one
  sign per row); in every other file (energies, `Ncm`, spectra, B values, NCSM energy tables, rho) a sign-only change
  is DIFF. `ok*`/`SIGN*`: equal once the spurious E(1) of EMPM_TDA is excluded (legacy 0, dev <spur|A|spur>,
  consumers.md 3.2.3; `SIGN*` only for `TDA_Structure.dat`). Exit status 1 on DIFF/SIZE/MISSING, so with the known
  He4o / O16w amplitude-basis DIFFs (45 + 2 files, see `docs/dev_port/STATUS.md`) the baseline comparison exits 1.
  `--tda`: a missing spectra / bin / `TDA_E<L>.dat` (L = 0..3) file of either side is MISSING (exit 1); the deviation
  printed per file is the one used for pass/fail.
* NCSM steps (stage 5 part 2): `06`-`08` (He4t), `13` (O16), `16` (O16w `Window`, after `14_O16w_HF_EMPM`) read the
  dev reference of `&empm Reference`. Expected: He4t NCSM files byte-identical to the baseline (same HF files);
  O16 / O16w energies equal, the `M_p`, `M_n` columns of `_E<L>.dat` may flip sign (`SIGN`, eigenvector phase; the
  16O dev HF differs from the legacy HF at 1e-15) and the residual column carries roundoff.

## Stage-6 regression (cleanup must not change any dev output)

The cleanup (removal of the LEGACY code) was checked against the frozen binary of the last stage-5 commit
(`/scratch/local.ssd/knapp/job_25503048.pbs-m1/nujuliet-snap-ec14d69/nujuliet`, commit ec14d69): every case with both
binaries at the same thread count, then `diff -rq` of the work trees (`run_fortran.sh` moves nothing, so rename
`work/fortran_t8` after each binary; EMPM/NCSM via `TAG=`):

```sh
for c in bundled ca40_5sh o16_open hf_options mbpt_options rpa_options qp_options scrpa_options; do
  NUJULIET=/path/to/frozen/nujuliet ./run_fortran.sh $c 8; done; mv work/fortran_t8 work/s6_frozen_t8
for c in ...; do ./run_fortran.sh $c 8; done; mv work/fortran_t8 work/s6_final_t8
diff -rq work/s6_frozen_t8 work/s6_final_t8                       # only run.log / time.log may differ
NUJULIET=/path/to/frozen/nujuliet TAG=_s6frozen ./run_empm_ncsm.sh dev 16
TAG=_s6final ./run_empm_ncsm.sh dev 16
diff -rq work/empm_dev_s6frozen_t16 work/empm_dev_s6final_t16      # logs, times.txt & the seconds column of Summary.dat
```

Two runs of the same binary differ in the same places (timings, free-memory limits, RSS, and the last printed digit of
a few full-precision energies in the logs, from OpenMP reduction order), so this is the noise floor.
