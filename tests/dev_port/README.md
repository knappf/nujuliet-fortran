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
* `*.dat`: token-wise; numeric tokens with tol plus one unit of the last printed digit, text exactly.
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
