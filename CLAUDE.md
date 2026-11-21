# CLAUDE.md — notes for working on nujuliet-fortran

Full documentation (namelists, solvers, outputs) is in `README.md`. These notes cover the cluster environment and
run setup.

## Environment (MetaCentrum PBS nodes, e.g. elwe10 = AMD EPYC 7532, 128 cores, ~2 TB RAM)

- The repo lives in the job's local scratch (`/scratch/local.ssd/knapp/job_<id>/nujuliet-fortran`), so it is lost
  when the job ends. A previous copy of the port (with older binaries & `benchmark14/`) is in
  `/storage/brno12-cerit/home/knapp/NuJuliet_Fortran_port/`.
- In a non-interactive shell `module` is undefined: first `source /etc/profile.d/30_meta_modules.sh`.

## Building

```sh
source /etc/profile.d/30_meta_modules.sh
module load intel-oneapi-hpc-compilers     # ifx 2025.1.1 + MKL
make                                       # -> ./nujuliet (~75 s, serial make)
```

- Default `COMPILER=intel`, `-O2 -march=core-avx2 -qopenmp`, which runs on both AMD & Intel nodes (do not use
  `-xHost`: those binaries refuse to run on AMD).
- The binary carries rpaths to MKL / libiomp5 / libimf, so it runs without the module loaded.
- Fallback: `make COMPILER=gnu` (system gfortran 14.2 + reference LAPACK/BLAS).
- `./nujuliet` with no arguments prints usage, which is a quick check that the binary starts.

## Interaction files

`/auto/brno12-cerit/nfs4/home/knapp/NuHamil_files/NN+NNN files/<interaction>/` (same as
`/storage/brno12-cerit/home/knapp/NuHamil_files/...`), interactions DN2LOGO394, DN2LOGO450, N3LOtexas, NNLOsat.
For DN2LOGO394, hw = 8..20 (step 2):
- NN: `NN-DN2LOGO394_vbare2.0_hw<hw>_emax16_e2max32.me2j.bin` (3.9 GB)
- NNN: `NNN_DN2LOGO394_hw<hw>_hwt<hw>_emax16_e2max32_e3max28_jmax3_59.stream.bin` (188 GB, read over NFS, slow)

In `&interaction`, Nmax/N2max/N3max describe the **file** (16/32/28). In `&calculation`, they set the truncation
actually used (N shells -> Nmax = N-1; convention used so far: N2max = 2 Nmax, N3max = min(2 Nmax, 28)).

## Runs

- Put runs in `runs/<name>/` (gitignored) and run from that directory: outputs go to `io/<Path>/` relative to
  the cwd (default Path e.g. `A208_Z82_hw12.0_Nmax12_N2max24_N3max24_CMS1+2B`; runs made before stage 5 used `IO/`).
- Solver mapping (dev): MBPT is part of `HF` (enabled by `&hf BMF=.true.`, the default; `&hf_mbpt Export` writes the
  HF-MBPT reference); `TDA`, `RPA` and `HF_ScRPA` are separate actions (one input file per `Name`/`Reference`); EMPM and
  NCSM read the HF reference of `&empm Reference` (default `io/<Path>/HF/bin`). The base-version `HF_RPA*`, `HF_RRPA*`
  and the interim `HF_LEGACY` were removed in stage 6 (fatal with a hint; README "Removed actions and input"); old
  inputs/runs with them (e.g. `runs/Pb208_13sh_hw12*`, `tests/dev_port/fortran/empm_ncsm_legacy/`) need a frozen binary
  of commit ec14d69 or earlier (`../nujuliet-snap-ec14d69/nujuliet`, `../nujuliet-snap-d753f29/nujuliet`).
- Run script pattern (see `runs/Pb208_13sh_hw12/run.sh`): set `OMP_NUM_THREADS`/`MKL_NUM_THREADS`,
  `OMP_STACKSIZE=512M`, `ulimit -s unlimited`, wrap in `/usr/bin/time -v`, and log to `logs/`. The user starts
  production runs manually, so prepare the scripts and do not launch them.
- Reference cost before the V3b lookup fix of 2026-11-20 (benchmark14, 15 shells, hw=16, HF+HF_RPA+HF_RRPA, 32 threads): O16 ~60 h, Ca48 ~32 h,
  peak RSS ~147 GB each.

## Original Julia code (reference / benchmark)

- The port follows NuJuliet dev: `../NuJuliet_dev` (commit f3f7a33). Entry point: `include("src/import_NuJuliet.jl");
  using .NuJuliet` (the dev sample is `scripts/sample_script_run.jl`, interaction files `io/NN.bin`, `io/NNN.bin`).
  Run from the directory that holds `io/` (dev writes `io/<Path>/...` and needs `io/` to exist).
- Reference for the comparisons: the patched worktree `../NuJuliet_dev_ref` (branch `port-ref` = f3f7a33 + the Julia
  bug fixes F1-F5, never pushed), driven by `tests/dev_port/run_julia.sh` (see `tests/dev_port/README.md`). The base
  version `../NuJuliet` (c450842, `include("Source/Import_NuJuliet_Module.jl")`, `IO/`) only matters for old runs of the
  removed base-version chain (its `runs/Pb208_13sh_hw12/` has the base-version 208Pb benchmark setup).
- Julia via juliaup: `~/.julia/juliaup/julia-1.10.11+0.x64.linux.gnu/bin/julia --project=@v1.10` (the env already
  has CGcoefficient, Hungarian, etc.). Uses OpenBLAS (not MKL), so set `OPENBLAS_NUM_THREADS` along with `-t N`.
- Parameter defaults match the Fortran namelist defaults (dev structs, `src/import_structures.jl`).
- Smoke test: `tests/dev_port/run_julia.sh bundled 8` (16O/4He/18O, bundled `io/NN.bin`, 8 threads): about 3.5 min.

## Comparing with Julia

- `tests/dev_port/` is the comparison setup (see its `README.md`): `run_julia.sh <case>` (reference in
  `tests/dev_port/work/julia_t<N>/`), `run_fortran.sh <case> [NT]`, `compare.py JULIA_IO FORTRAN_IO --tol 1e-8`.
  `runs/test_Ca40_9sh/` holds the older base-version comparison (`IO/`, `HF_RPA`, `HF_RRPA`, frozen binaries only).
- Expected differences that are not bugs:
  - `HF/bin/C.bin` (and the gauge-dependent `H2b.bin`, amplitudes): the columns can differ in sign.
  - pn block of two-body files (`H2b.bin`, `sigma.bin`): Julia writes the full N x N square with the packed
    lower-triangle index, so its upper triangle holds aliased entries. The port writes the true symmetric matrix. Each
    code reads its own file back correctly, and the port can read Julia files, but **Julia cannot read a pn block
    written by the port**.
  - The three known O20_HBCS HFB files of `qp_options` (docs/dev_port/STATUS.md).
- The Julia NNN reader keeps the whole file resident: 187 GB RSS for 40Ca at 9 shells, against 0.9 GB for the port.

## Session log

- 2026-11-19: built `./nujuliet` with ifx+MKL on elwe10. Prepared `runs/Pb208_13sh_hw12/` (`Pb208.nml` + `run.sh`):
  208Pb, DN2LOGO394, hw=12, 13 shells (12/24/24), solvers HF (+MBPT), HF_RPA (TDA+RPA), HF_RRPA. The user runs it
  manually: `./run.sh [NTHREADS]`.
- 2026-11-20: the Fortran Pb208 run was stopped by the user's request after 22.9 h (peak 34 GB). HF, TDA and RPA
  had finished. RRPA did not converge: the outer `eta` stalled at 2e-6..1e-5 after about 8 iterations, and each
  inner OBDM loop ended in "cycled degenerate solution". Prepared the matching Julia benchmark in
  `../NuJuliet/runs/Pb208_13sh_hw12/`.
- 2026-11-20: **performance fix** in `V3b_no2b_index` (`src/core/nj_v3b.f90`): the type-bound `Orb_NNN%Dic%get` call
  was replaced by `hashmap_lookup`. With ifx, the polymorphic call wrote a class descriptor into shared static
  memory on every lookup, so all threads contended for the same cache lines. This lookup is the hot path of the
  HF iterations, `V2b_residual_no2b` and every RRPA iteration. 40Ca, 9 shells, 16 threads, total wall time:
  old 48:35, new 3:02, Julia 13:49. HF_RRPA alone: 2796 s, then 166 s. The outputs are identical to 1e-12. The
  old binary is kept as `nujuliet.old_dicget`.
- 2026-11-20: prepared `runs/Pb208_13sh_hw12_v2/` (same physics as `Pb208_13sh_hw12`) with a frozen binary
  `nujuliet_78c3588` (dev-port stage 2, built in the worktree `../nujuliet-snap-78c3588`). Outputs go to `io/`; 'HF'
  runs HF_LEGACY automatically because HF_RPA/HF_RRPA are in the same &run. 16O smoke test of the chain passed. The user
  launches it: `./run.sh [NTHREADS]`.
- 2026-11-20: the v2 Pb208 run was stopped by the user's request after 3:05 (peak 34 GB). HF, MBPT, TDA and RPA had
  finished in under 2 h (outputs in `runs/Pb208_13sh_hw12_v2/io/`). RRPA stalled as before: eta 2e-6..7e-6 over
  iterations 8-11, inner loops "cycled degenerate".
- 2026-11-21: dev-port stage 4 (HF-ScRPA, commit d753f29). 40Ca 9 shells: HF_ScRPA converges in 26 iterations, 14.6 s at
  32 threads (legacy HF_RRPA 154 s, stalled). Prepared `runs/Pb208_13sh_hw12_scrpa/` with the frozen binary
  `nujuliet_d753f29` (worktree `../nujuliet-snap-d753f29`): `Pb208_scrpa.nml` (HF -> HF_ScRPA) then `Pb208_rpa.nml`
  (TDA, RPA on the same HF). Estimate: HF ~2 h, ScRPA 15-25 min. The user launches it: `./run.sh [NTHREADS]` (default 64).
- 2026-11-21: dev-port stage 5: EMPM and NCSM read the dev HF reference (`&empm Reference`, default 'HF'), EMPM outputs in
  `io/<Path>/EMPM/{bin,spectra,transitions}`. Checked against a legacy baseline (`tests/dev_port/run_empm_ncsm.sh`).
- 2026-11-21: the user ran `runs/Pb208_13sh_hw12_scrpa/` (64 threads): whole chain 6:14, peak 34 GB. HF-ScRPA converged in 10
  outer iterations, E_QBA = -845.029752 MeV (E_HF = -752.626164, plain RPA E_QBA = -1133.811798, both equal to the v2 run).
- 2026-11-21: dev-port stage 6 (cleanup): legacy HF/HF-EFA/HF-RPA/HF-RRPA chain removed (~6300 lines); the removed actions
  stop with a hint. Old inputs (`runs/Pb208_13sh_hw12*`, `runs/test_Ca40_9sh/Ca40.nml`) need the frozen binaries
  (`nujuliet_78c3588`, `../nujuliet-snap-d753f29`, `../nujuliet-snap-ec14d69`). `make COMPILER=gnu` works again.

