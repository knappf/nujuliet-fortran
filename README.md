# NuJuliet — Fortran port

Fortran (2008) port of [NuJuliet](https://github.com/folpradek/NuJuliet), a spherical nuclear
many-body solver originally written in Julia. The port follows NuJuliet **dev f3f7a33** (solvers, file formats and
output files); EMPM and NCSM are Fortran-only extensions that read the dev HF output.

| Action (`&run`)                    | Julia call (dev)                    |
|------------------------------------|-------------------------------------|
| `HF`                               | `HF(Params)` (HF, `V_HF`, `ESH`, MBPT(2,3) energies, HF-MBPT export `&hf_mbpt`) |
| `TDA`, `RPA`                       | `TDA(Params)`, `RPA(Params)` (reference `io/<Path>/<Reference>`, output `io/<Path>/<Name>`) |
| `HF_ScRPA`                         | `HF_ScRPA(Params)` (HF-based self-consistent RPA, `&rpa`) |
| `TDA_transition_densities`, `RPA_transition_densities` | `TDA_/RPA_transition_densities(Params; JP, nu_list, File_Name)` |
| `TDA_transition_operator_densities`, `RPA_transition_operator_densities` | `TDA_/RPA_transition_operator_densities(...)` |
| `TDA_transition_analysis`, `RPA_transition_analysis` | `TDA_/RPA_transition_analysis(...; Cut)` |
| `TDA_transition_currents`, `RPA_transition_currents` | `TDA_/RPA_transition_currents(...; M, xN_grid, zN_grid)` |
| `TDA_transitions_recalculate`, `RPA_transitions_recalculate` | `TDA_/RPA_transitions_recalculate(Params)` |
| `BCS`, `HFB`, `QTDA`, `QRPA`, `pnQTDA` | `BCS`, `HFB`, `QTDA`, `QRPA`, `pnQTDA` |
| `V2b_parse`, `V3b_no2b_parse` (`PARSE` = both) | `V2b_parse`, `V3b_no2b_parse` (`&parse`) |
| `EMPM_TDA`                         | — EMPM: F-matrix & TDA in a restricted ph space (`Fmat`, `Tda`) |
| `EMPM_2PH`                         | — EMPM: 2-phonon equations of motion (`phon_dens1`, `phon_int`, `eqm2`) |
| `EMPM`                             | — complete EMPM: `EMPM_TDA` + `EMPM_2PH` + couplings & full Hamiltonian (`nondiag`, `fullham`) |
| `EMPM_3PH`, `EMPM_4PH`, `EMPM_FULL` | — EMPM: 3-/4-phonon equations of motion, full Hamiltonian with up to 1 .. 4 phonons |
| `EMPM_TRANS`                       | — EMPM: B(E1; 0+_gs -> 1-) of the `EMPM_FULL` states (Full_N <= 2) |
| `EMPM_SF`                          | — EMPM: E1 strength function with up to 3 phonons by Lanczos, no 3-phonon eigenproblem |
| `NCSM`                             | — no-core shell model (HO Nmax or npnh in the HF basis), optional transitions |
| `HF_LEGACY`, `HF_RPA*`, `HF_RRPA*` | LEGACY base-version solvers (`io/<Path>/Bin/*`), kept until the cleanup stage; `HF` in a `&run` with `HF_RPA*` / `HF_RRPA*` runs `HF_LEGACY` |

> **Earlier test runs (until 2026-11-06)** are in `backup_2026-11-06/`. It holds `runs/` (inputs, logs, results and
> plots; the `runs/...` paths quoted below refer to it), `tests/` (the Python cross-checks), `doc/` (the notes), `bin/` (the old binaries
> and builds), `snapshot/` (src, Makefile, README and doc as of that date) and `runs_cmtop/` (the validation runs of
> the top-level <N_cm> completion, `4He_NN4sh_cmtop`, `4He_NN5sh_cmtop`, moved there on the same day). `runs/` starts empty.
> Runs and inputs before stage 5 (2026-11-21) use the base-version layout (`IO/<Path>/Bin/...`, `HF_RPA`).

## Building

Requirements: a Fortran 2008 compiler with OpenMP (tested with gfortran 9.4) and LAPACK/BLAS.

```sh
make                      # optimized build -> ./nujuliet
make LAPACK="-lopenblas"  # use another LAPACK/BLAS
make debug                # runtime-checked build -> ./nujuliet_debug
```

## Running

Run from the directory that contains `io/` (outputs go to `io/<Path>/`, as in dev):

```sh
OMP_NUM_THREADS=16 ./nujuliet O16_HF_TDA.nml O16_RPA.nml [more inputs ...]
```

The input files run in sequence. Each describes one calculation; all namelist groups are optional and omitted
values take the defaults of the dev parameter structs (`src/import_structures.jl`).

```fortran
&interaction  NN_File='io/NN.bin', NNN_File='io/NNN.bin', hw=16.0, Nmax=3, N2max=6, N3max=9
              ! cP2N, cP3N, cRes, pG2N, nG2N /
&calculation  A=16, Z=8, hw=16.0, Nmax=3, N2max=6, N3max=9, CMS='CMS1+2B', Path='' /  ! Path='' -> default name; CM_Lawson (Fortran)
&hf           Tol=1e-7, Imax=100, V_HF=.false., BMF=.true., ESH=.false. /     ! HRF accepted, dead in dev
&hf_mbpt      Export=.false., TBCF=.false., ScH=.false. /                     ! HF.MBPT: export io/<Path>/HF-MBPT/bin
&tda          Name='TDA', Reference='HF', Reference_TBCF=.false., Ortho=.true., Tol_hol=0.7, Tol_par=0.1 /
&rpa          Name='RPA', Reference='HF', Reference_TBCF=.false., RPA_Orb=.false., Ortho=.true., Transitions='Linear',
              Tol=1e-5, Tol_hol=0.7, Tol_par=0.1, Imax=150, OBDM_diagonal=.false., OBDM_alpha=0.6, ScH=.true.,
              TBCF=.false., TBCF_alpha=0.8, TBCF_beta=5.0, DEV=.false. /     ! also HF_ScRPA
&bcs          ScBCS=.false., Tol=1e-7, Imax=500, q=0.05, A0=0, Z0=0, pD0=0.5, nD0=0.5, BMF=.true. /
&hfb          Tol=1e-7, Imax=150, Pairing='Full', dLmax=0.5, pL0=-3.0, nL0=-3.0, pa0=0.25, na0=0.25,
              Broy_m=8, Broy_Amax=0.8, Broy_Tmax=0.9, LNT=.false., LNRes=.true., BMF=.true. /
&qtda         Ortho=.true., Reference='HFB' /
&qrpa         Ortho=.true., Reference='HFB' /
&pnqtda       Reference='HFB' /
&empm         Reference='HF', pHole=0,0, pPart=0,0, nHole=0,0, nPart=0,0, Ortho=.true., E2max=1e8, JP2='', ... /  ! see EMPM
&ncsm         Hamiltonian='NN', Basis='HO', Trunc='Nmax', Nmax_List=..., ... /                                  ! see NCSM
&transition_densities           JP='1-', nu_list=5,6,7, File_Name='Phonon_Densities' /
&transition_operator_densities  JP='1-', nu_list=5,6,7, File_Name='Phonon_Transition_Operator_Densities' /
&transition_analysis            JP='1-', nu_list=5,6,7, File_Name='Phonon_Transtion_Analysis', Cut=0.1 /   ! [sic] dev name
&transition_currents            JP='1-', M=0, nu_list=5,6,7, File_Name='Phonon_Currents', xN_grid=40, zN_grid=40 /
&parse        N_max_Out=..., N_2max_Out=..., N_3max_Out=..., NN_Out_File='io/NN_parsed.bin', NNN_Out_File='io/NNN_parsed.bin' /
&run          Solvers='HF', 'TDA', 'RPA' /
```

* One input file per TDA/RPA flavour: a namelist group is read once per file, so e.g. an HF-TDA and an HF-ScRPA
  (different `Name`/`Reference`) need two files, run in one `nujuliet a.nml b.nml` call. Only the first file should
  contain `HF`: HF wipes `io/<Path>`. `Name` and `Reference` are folder names (non-empty, no `/`); `Name = Reference`
  is refused.
* `M` of `&transition_currents`: dev's default -1 means 0; the port keeps the default 0 (same result).
* LEGACY groups `&hf_efa`, `&rrpa` and the HFB `Broy_Qmax`, `Broy_Bmax` are still accepted (base-version solvers, ignored
  by the dev solvers) until the cleanup stage.
* The dev sample `scripts/sample_script_run.jl` (with `io/NN.bin`, `io/NNN.bin`, the bundled ΔN²LO_GO(394) files) is
  mirrored by the `bundled` case of `tests/dev_port` (namelists in `tests/dev_port/fortran/bundled/`). The `inputs/`
  directory quoted in older notes is not part of the repository.

## Output layout

```
io/<Path>/HF/            HF_summary.dat, densities/, bin/{orb,C,rho,H1b,H2b}.bin   (H2b.bin only with BMF = .true.)
io/<Path>/HF-MBPT/bin/   {orb,C,rho,H1b,H2b}.bin of the HF-MBPT natural orbitals (&hf_mbpt Export)
io/<Path>/<TDA Name>/    TDA_summary.dat, bin/, spectra/, transitions/, amplitudes/, densities/, levels/
io/<Path>/<RPA Name>/    RPA_summary.dat, bin/, spectra/, transitions/, amplitudes/, densities/ (also HF_ScRPA)
io/<Path>/{BCS,HFB,QTDA,QRPA,pnQTDA}/   <solver>_summary.dat, bin/, spectra/, ...
io/<Path>/EMPM/          TDA_summary.dat, Levels.dat, TDA_Structure.dat, bin/, spectra/, transitions/,
                         2Phonon/, 3Phonon/, 4Phonon/, Full/<n>Phonon/, SF/
io/<Path>/NCSM/          NCSM_<Ham>_..._dat (+ _E<L>.dat, _SF_*, checkpoints)
```

`orb.bin` holds one record per orbital a = 1 .. a_max(Calc.Nmax): `a, n, l, 2j` (int64) and the proton & neutron
occupations `pO, nO` (float64), little-endian stream. `C.bin` (HF orbitals in the HO basis), `rho.bin` (density in the
HF basis), `H1b.bin` (one-body Hamiltonian, diagonal for HF) and `H2b.bin` (normal-ordered two-body Hamiltonian, O2b
layout) are the reference read by TDA/RPA (`Reference`), EMPM and NCSM (`&empm Reference`).

## EMPM extension (`src/empm/`)

First step of the Equation of Motion Phonon Method, ported from the `Fmat` & `Tda` programs of the original
EMPM codes. `EMPM_TDA` runs on top of a finished `HF` run: it reads the reference
`io/<Path>/<Reference>/bin/{orb,C,H1b,H2b}.bin` of `&empm Reference` (default `'HF'`). The reference must be
HF-like: `H1b.bin` diagonal (off-diagonal elements <= 1e-10, so `'HF-MBPT'` is refused), `H2b.bin` present (`&hf
BMF = .true.`) and **integer occupations** in `orb.bin` (0 or 1 within 1e-8): a closed (sub)shell nucleus, or `&hf
ESH = .true.` for an open shell (the dev HF with `ESH = .false.` gives fractional occupations, which are fatal). A
missing reference prints "Error! ... No precomputed reference ..." and skips the step (as the dev TDA).

1. restricted s.p. space ... windows `(first, last)` of the HF levels of each isospin ordered by energy
   (1 = deepest level), the numbering of `input_tda_coup.dat`; a zero bound means the full range
   (holes 1 .. N_occupied, particles N_occupied+1 .. a_max); windows must not mix occupied & empty levels,
2. F-matrix `F^L(i,j,k,l) = sum_J (2J+1) (-1)^(j_j+j_k-L-J) {j_i j_j L; j_l j_k J} V^J(ik,jl)` over all level
   pairs of the restricted space (convention of `fmat_unfi.f90`),
3. TDA `A(ph,p'h') = (e_p - e_h) delta - (-1)^(j_p'+j_h'+J) F^J(p,h,h',p')`, spurious 1- state removed with the dev
   convention (as the dev `TDA`): the CM vector |spur> is kept as the first 1- phonon with E(1) = <spur|A|spur> (not 0),
   the other 1- phonons from the symmetric eigenproblem in its orthonormal complement (`Ortho = .true.`).

```fortran
&empm  pHole = 7, 16, pPart = 17, 45, nHole = 11, 22, nPart = 23, 55, Ortho = .true. /   ! 208Pb, hw_12_55_jshells
&run   Solvers = 'EMPM_TDA' /
```

Outputs in `io/<Path>/EMPM/` (`EMPM_TDA` checks the reference, then wipes this folder, so a new TDA also removes the
later EMPM binaries): `Levels.dat` (the restricted space), `TDA_summary.dat` (the dev `TDA_summary` text with the EMPM
items), `spectra/TDA.dat`, `transitions/TDA_E<L>.dat`, `TDA_Structure.dat` (phonon composition), `bin/EMPM_F.bin`
(F-matrix, channel by channel `L`, `P`), `bin/EMPM_TDA_E.bin` & `bin/EMPM_TDA_X.bin` (channel by channel `J`, `P`, as
the TDA `bin/TDA_*.bin`). Every EMPM binary (2PH, FULL, state files, CM dumps, SF scratch) is in `EMPM/bin/`. The pair &
phonon structure behind the binary files is rebuilt from the `&empm` windows (`Fmat_initialize`, `EMPM_TDA_phonons`).
With the full space (no `&empm` group) the energies, amplitudes and transitions equal those of the `TDA` action with
`Reference = 'HF'` (incl. the spurious E(1)). The later steps (`EMPM_2PH` .. `EMPM_SF`) print "run EMPM_TDA first" and
skip when `EMPM/bin` is missing; they recompute the CM energy from A, so they do not depend on E(1).

### 2-phonon step (`EMPM_2PH`)

Port of `phon_dens1`, `phon_int`, `eqm2_svd` & `eqm2_admat`. It runs on top of `EMPM_TDA` (reads `EMPM/bin/EMPM_F.bin` &
`EMPM/bin/EMPM_TDA_*.bin`, rebuilds the space from the same `&empm` windows) ... **nothing is written to disk in between**:
the 1-phonon densities rho_gb(s;i1,i2) & the redefined interaction V_gb(s;k,l) of a phonon g (for all phonons b) are
computed when a row of A & D needs them (`nj_empm_dens.f90`), as dense blocks per (J,P) channel of b, so that the A & D
rows are dgemm contractions.

* basis |alpha lambda; J> of all phonon pairs, spurious states (alpha or lambda = spurious 1- phonon) first; truncated
  rows = all spurious states & alpha >= lambda with E_alpha + E_lambda <= `E2max`; channels `JP2 = '1-', '2+'` (none = all),
* densities & interactions are stored & contracted only over the pairs (i1,i2) allowed by the multipolarity sigma
  (triangle & parity, ~16% of the pair space), densities & exchange terms are dgemm's,
* the lambdas of the rows are processed in blocks whose cached V_l & rho~_l fit into `Cache_GB` (437 GB for all 1824
  lambdas of 208Pb 1-), the rows of a block in parallel over alpha (rho_a once per task of <= 4 rows); H = A D^T over the
  full basis is built in row blocks so that at most `Mem_GB` of A & D rows are in memory (A_I, D_I and A_J, D_J of a
  later block J, rows of later blocks recomputed); both H_IJ = A_I D_J^T and H_JI = A_J D_I^T are computed, so the
  symmetry check |AD - (AD)^T| always covers the whole H (208Pb 1-, E2max = 35: Mem_GB = 1000, Cache_GB = 150 give
  ~1.2 TB peak instead of ~2.3 TB with one block, at the price of recomputing the rows of later blocks),
* then pivoted Cholesky (tol 1e-6) of D, SVD of the spurious D columns & the spurious subspace projected out exactly
  (see "CM treatment" below; the 1e8 MeV shift of `ham_geev` only in the limit), C normalized by C^T D C = 1. The spurious
  1- phonon enters with <cm|A_TDA|cm> & the couplings <nu|A_TDA|cm>, as `tda_r_overl.dat`.

```fortran
&empm  pHole = 7, 16, pPart = 17, 45, nHole = 11, 22, nPart = 23, 55, E2max = 35.0, JP2 = '1-', Mem_GB = 400.0 /
&run   Solvers = 'EMPM_2PH' /
```

Outputs: `EMPM/2Phonon/Spectrum.dat` (J, P, E, spurious flag), `EMPM/2Phonon/Summary.dat` (dimensions, timings,
|AD - (AD)^T|), `EMPM/2Phonon/Basis_<J><P>.dat` (truncated basis with spurious & independence flags),
`EMPM/bin/EMPM_2PH_<J><P>_E.bin` & `_C.bin` (energies, amplitudes C(independent row, solution)).
The log also prints the eigenvalue spectrum of the spurious & truncated D blocks per decade.

Validated against the original binaries (208Pb, small window, 1-, 2+, 0+): identical dimensions, A & D rows equal to
1e-7 relative, energies equal to ~5e-7 MeV (median). States at the truncation edge are ill-conditioned (smallest
non-zero eigenvalues of the truncated D ~1e-5) and differ up to 1e-2 MeV, the single-precision files of the original
codes amplified by ~1e5; with the original A & D the NuJuliet solver reproduces the original energies to 5e-7 MeV.

### Couplings & full Hamiltonian (`Full_H`, action `EMPM`)

`Full_H = .true.` (always on for the action `EMPM`, which also runs `EMPM_TDA` first) couples, for each channel, the
1-phonon states a of J^P (& for 0+ the HF ground state) to the 2-phonon states. The redefined interaction
V(a; alpha lambda) = (-1)^(Ja+Jal+Jl)/(2Ja+1) sum_s rho_(a alpha)(Jl; s) W_l(s), W_l(s) = sum_ph kappa c^l_ph F^Jl(p,h; s)
(`nondiag12`) & the ground-state V0 (`vintn2`) are computed on the fly over the full basis & projected on the correlated
states while the D rows are built, P = D V, H_12 = sqrt(2J+1) C^T P (`fullham`), nothing is stored on disk. The full
Hamiltonian [ground ; 1-phonon (without the spurious phonon) ; physical 2-phonon] is diagonalized by `dsyevd`.
Outputs `EMPM/2Phonon/Full_<J><P>.dat` (E with 0-, 1- & 2-phonon weights), `Full_<J><P>_1phon.dat`,
`EMPM/bin/EMPM_FULL_<J><P>_E.bin` & `_Z.bin`. Validated against `Nondiag12` + `Fullham` (small 208Pb window, 1- & 0+ with the
ground state): same spectra & weights to the single precision of the original files.

Eigensolver of the 2-phonon step: `dsygvd` (divide & conquer) replaces `dsygv`; it returns C normalized by C^T D C = 1. `D_Spectrum = .true.` prints the
eigenvalues of the whole truncated D (a dense n^3 diagnostic); the small spurious block is always reported.

Legacy options for one-to-one comparisons with the original runs: `E1max` (cut E_lambda <= E1max, 2nd line of `en_trun.dat`),
`Legacy_Order` (phonon numbering of the original `Tda`), `Spur_Rows = 'lambda'` (spurious rows only with lambda = spurious
phonon) & `Legacy_Spur` (spurious subspace selected as the original `eqm.f90`, Cholesky of D(rows 1..N_spur, columns
1..N_spur), which misses spurious directions when only the lambda rows are kept). With all four the 208Pb 1- run of May 2025
(51041 states) is reproduced (lowest 558 states within 1.2e-4 MeV); its 5.82 MeV 1- state is a spurious leftover.

With the Intel compiler set `FORT_BUFFERED=true` ... formatted output is unbuffered by default, which makes
the text exports very slow on NFS.

### 3-phonon step (`EMPM_3PH`)

Counterpart of `eqm3_admat` & `phon_dens2` (`ro7_MPI.f90`) of the original codes, on top of `EMPM_TDA`
(`nj_empm_eqm3.f90`), without any files in between:

* the physical 2-phonon states of all channels J = 0 .. 2 J_max (truncation `E2max`, `E1max`) are solved in memory and
  kept with their amplitudes C (on the full 2-phonon basis) & X = sqrt(2J+1) D(independent rows)^T C,
* 3-phonon basis |a lambda; J^P> of the 2-phonon states a & the TDA phonons lambda (all couplings), rows
  E_a + E_lambda <= `E3max`, channels `JP3` (none = all),
* 2-phonon densities on the fly, one row state a at a time for all 2-phonon states b (ro7):
  rho_ab(s) = sum_(mu nu) X^b(mu nu) T^a_(mu nu)(s) with
  T^a_(mu nu)(s) = sum_al (-1)^(j_mu+Ja+j_nu+s) sqrt(2Ja+1) {s j_mu j_al; j_nu Ja Jb} C^a(al nu) rho_(al mu)(s)
  + sum_l1 (-1)^(j_mu+Jb+j_l1+s) sqrt(2Ja+1) {s j_nu j_l1; j_mu Ja Jb} C^a(mu l1) rho_(l1 nu)(s)
  (1-phonon densities rho of all phonons cached), one dgemm with X^b per (s, channel of b),
* A & D rows (a lambda) as in the 2-phonon step with rho_(alpha beta) -> rho_ab & E_alpha -> E_a, and the metric
  exchange term D += sum_g (-1)^(j_l'+j_l+Ja+Jb) {j_l j_g Jb; j_l' J Ja} X^a(g l') X^b(g l),
* independent rows first (`Ind_First = .true.`, default): D of all truncated rows only (no A, row blocks within
  `Mem_GB`), a blocked pivoted Cholesky selects the independent rows, then A & D rows of those only and
  H = A D^T on them (r x r instead of N x N), the |AD - (AD)^T| check on that block & the eigenproblem ... same states as
  the full H (identical spectra & pivots on complete 4He 3/4-phonon), a large gain when the rows are strongly dependent
  (complete 12C 0+: 64844 rows, 10516 independent, A D^T 723 s -> 20 s, half the memory), none when the truncated rows
  are all independent (only the extra D pass),
* `Ind_First = .false.`: H = A D^T on all truncated rows (row blocks within `Mem_GB`, H_IJ & H_JI both computed), the
  |AD - (AD)^T| check & the eigenproblem as in the 2-phonon step.

The pivoted Cholesky of the independent-rows path (`chol_piv`, `nj_empm_eqm2.f90`) pivots greedily on the largest
remaining diagonal like LAPACK dpstrf (same pivot set: identical on all complete 4He 3/4-phonon channels & on complete 12C
0+), but left-looking: the Schur-complement columns of a pool of likely pivots (the rows of largest remaining diagonal)
are kept up to date with one dgemm per panel of 64 pivots, so the cost is ~N r^2 instead of N^2 r and D is not
overwritten (complete 12C 0+, N = 64844, r = 10516: 42 s vs 99 s of dpstrf on 28 threads).

**Greedy independent rows (`Ind_Greedy = .true.` in `&empm`, 3/4-phonon steps with `Ind_First`).** For large
overcomplete bases even D of all truncated rows does not fit (4He in 4 shells, 4-phonon 3+: 188k rows, D ~280 GB).
`EMPM_ind_greedy` (`nj_empm_eqm3.f90`) never forms the full D: the rows are grouped by the (n-1)-phonon state a, the
groups are visited in a fixed pseudo-random order (LCG, reproducible), and for each block of rows (>= 256, ~Target/4) only
its D rows are computed; the Schur complement on the rows selected so far, S = D_bb - D_bI D_II^-1 D_Ib (dpotrf/dpotrs/
dgemm), is reduced by the pivoted Cholesky (`chol_piv`, absolute tolerance 1e-6) and the new pivots are added to the
selection. The loop stops when the selection reaches the npnh dimension of the channel (`EMPM_npnh_dims`,
`nj_empm_space.f90`: m-scheme counts of the n-particle n-hole configurations by a subset DP, d_J = N(2J) - N(2J+2)), which
is the exact rank of D for a complete basis; the log reports the rows selected, the npnh dimension and the rows visited
(warning if they differ). Only the visited rows cost D, so the cost is ~Target instead of N_tr rows. 4He in 3 shells:
3PH (1948 states) & 4PH (3739 states) spectra identical to the full-D path; 4 shells: every channel reaches the npnh
dimension with < 1% extra rows visited (0+ 1583 of 46028 rows, 3+ 5967 of 187560), EMPM_4PH 2:53 h, 26 GB on 32 threads.

Spurious states with `Ortho` as in "CM treatment" below. Outputs in `EMPM/3Phonon/`: `Spectrum.dat` (J, P, E, spurious
flag), `Summary.dat` (2- & 3-phonon channels with dimensions, spurious states & solutions, |AD - (AD)^T|), `States2.dat` &
`States3.dat` (index, J, P, E, spurious flag) and, as input of the 4-phonon step, the state sets `EMPM/bin/EMPM_S2_<JP>.bin`
(2-phonon) & `EMPM/bin/EMPM_S3_<JP>.bin` (3-phonon), raw streams of -3 (format marker), N_full, n, N_ind, N_phys (int32),
Al(N_full), La(N_full) (int32: the (n-1)-phonon state of each basis state, numbered as in `States2.dat` for S3 or the phonon
index for S2, and the phonon lambda), E(n), Kc(N_ind) (int32, full-basis indices of the independent rows, ascending),
C(N_ind,n) on those rows (zero on all others), X(N_full,n) = sqrt(2J+1) D(independent rows)^T C & R(N_phys,n-N_phys)
(float64, column-major); C^T X = sqrt(2J+1) 1; the N_phys physical states first, then the spurious ones with their couplings
R = <physical|H|spurious>. Files of format -2 (no N_phys & R, all states physical) & of the first format (no marker,
C(N_full,n) on the full basis) are still read. `States_Dir` puts these files into another folder than `Bin` (a leading
environment variable is expanded, e.g. `States_Dir = '$SCRATCHDIR/C12'` for the local disk of a PBS job ... all steps
of a chain then need the same `States_Dir`, and the files are lost with the job's scratch unless copied back).
The rank of the pivoted Cholesky keeps only pivots above its tolerance (MKL's dpstrf returns rank 1 for a numerically zero
D). Complete 4He (3 shells, no truncation): the independent 2- & 3-phonon states equal the 2p2h (332) & 3p3h (1948)
J-states in every J^P channel, the Pauli-forbidden 3-phonon channels 8+-, 9+ come out empty. 16O (DN2LOGO394, hw = 16 MeV,
4 shells, 13824 2-phonon states, full 3-phonon bases up
to 533528): all channels with E3max = 40 MeV (10k rows) in 104 s on 28 threads, |AD - (AD)^T| <= 4e-13 in every channel
(a sign error planted in the exchange term gives 4.3), row blocks identical to one block.

### 4-phonon step (`EMPM_4PH`)

On top of a complete `EMPM_3PH` run (same output folder): the 2- & 3-phonon states are read from `EMPM/bin/EMPM_S2_*.bin` &
`EMPM/bin/EMPM_S3_*.bin`, the 2-phonon densities of all 2-phonon states are cached and the 3-phonon densities of the row states
computed on the fly with the same recoupling formula one level up (the (n-2)-phonon densities rho_(al mu) are then the
2-phonon densities). Basis |c lambda; J^P> of the 3-phonon states c & the phonons, rows E_c + E_lambda <= `E4max`, channels
`JP4`; the A & D rows, H = A D^T, the |AD - (AD)^T| check, the eigenproblem & the exports (`EMPM/4Phonon/Spectrum.dat`,
`Summary.dat`, `States4.dat`, `EMPM/bin/EMPM_S4_<JP>.bin`) as in the 3-phonon step (shared code, `nph_channels` &
`nj_empm_eqm3.f90` with the lower level as an argument). Complete 4He (3 shells): 3739 independent 4-phonon states = the
4p4h J-states in every J^P channel, |AD - (AD)^T| <= 5.4e-13, 94 s on 28 threads.

### CM treatment (`Ortho = .true.`) at all levels

The EMPM_TDA spurious 1- phonon cm (the normalized E1 vector, the other 1- phonons in its orthogonal complement) is kept
as a basis phonon with <cm|A_TDA|cm> & the couplings R(nu) = <nu|A_TDA|cm>. The A·D^T sums need complete intermediate sets,
so the spurious states are never dropped inside a level, only projected out of its solutions:
* n-phonon spurious subspace S_n = span of the states |a cm> (a over the complete (n-1)-phonon set) and |a lambda> with a
  spurious (contained in it); they are tagged in the basis (`Spur`, `Ks`), `Spur_Rows = 'all'` keeps them as rows, `'none'`
  only as columns (the physical pairs span the physical subspace; same spectrum, fewer rows);
* the physical subspace is the null space of the D elements <s|row> of the spurious states (SVD, tol 1e-4) ... with
  Vs^T = Q R the physical states are the eigenstates of P H P, the spurious ones those of S H S (D-orthogonal to all
  physical ones), and R = <physical|H|spurious> is kept (`spur_gev` in `nj_empm_eqm2.f90`, instead of the 1e8 shift);
* each step exports physical & spurious states (a complete set of its level) with R; the next level's A rows get
  d_ll' H^(n-1)_ab (diagonal E & R) + d_ab R_ll' (phonon couplings to cm), as the 2-phonon A;
* `EMPM_FULL` leaves all spurious states (cm phonon & the spurious n-phonon states) out of its basis, the couplings
  (H20, G1, X^T V, n+2 <- n) sum over the complete sets.
Complete 4He (3 shells, `runs/4He_CM`, `run_He4_CM_chain.sh`, ~5 min): |AD - (AD)^T| <= 6e-13 at every level; the physical
2-, 3- & 4-phonon spectra agree with the exact projection (npnh space minus S_n, S_n = npnh part of O+_cm |(n-1)p(n-1)h>,
O+_cm the ph part of r_mu) to the 1e-6 print precision and EMPM_FULL with up to 1 .. 4 phonons with the projected FCI to
6e-10 MeV in every J^P channel (`tests/empm_cm_chain_check.py HRF_DIR EMPM_DIR`, ~3 min; `tests/empm_cm2_check.py` for the
2-phonon level alone, `tests/empm_cm_check.py` N_cm of the states without CM treatment). Spurious dimensions (M = 0, +/-):
1p1h 0/1, 2p2h 19/32, 3p3h 388/451, 4p4h 1239/1860 (all negative-parity 4p4h states of the 3-shell space are spurious).

The definition S_n itself is an approximation. `tests/empm_ho_nmax_check.py HRF_DIR` (~2 min) rotates the FCI space to HO
determinants, keeps the Nmax = 2 HO space (complete in 3 shells, CM separation exact: integer N_cm, Lawson energies
independent of beta) and compares the exact spurious space b+_cm (Nmax-1) with the EMPM-style one (HO reference, S_n from
the ph part of r). The dimensions agree (19/1), negative parity (1 hw) is exact, but in positive parity 12 of the 19
directions overlap only with cos^2 = 0.5: b+_cm on 1p1h(0s -> 0p) = ph part (2p2h, kept in S_2) + pp part (0p -> 1s0d,
1p1h, missing), so CM-excited strength (<N_cm> up to 0.5) stays in the physical space and e.g. the 2+_1 is 26.6 MeV instead
of 34.0 MeV. `tests/empm_fci_lawson.py HRF_DIR` does N_cm / Lawson in the full 3-shell space, which is not Nmax-complete
(non-integer N_cm, beta-dependent Lawson energies), so it is not an exact reference there.

**Exact CM condition (`CM_Exact = .true.` in `&empm`, `EMPM_FULL` only, with the complete sets of `Ortho = .false.`).**
The physical states are those with b_cm psi = 0: `nj_empm_cm.f90` builds the reduced matrix elements of the rank-1
operator L = sum_i b(i) (b = sqrt(2) (r/b)_(N -> N-1), the HO lowering part of r, transformed to HF, so with ph, hp, pp & hh
parts) between all EMPM states of the levels 0 .. Full_N ... ph & hp parts from the 1- phonon strengths and X of the next
level, pp & hh parts from the rank-1 n-phonon densities (`EMPM_densityN(..., Sig_Only = 1)`) ... N_cm = G2/A sum_J'
|<J'||L||psi>|^2/(2J+1) per channel, H diagonalized on its null space (N_cm < `CM_Tol`, default 1e-4) for every truncation
n = 1 .. Full_N (b changes the phonon number by at most 1, so the rows of level n+1 make the condition exact for n < Full_N).
`Full_unproj_<JP>.dat` lists E & N_cm of the unprojected states, the log the N_cm range & gap per channel;
`CM_Dump = .true.` writes the full H & the L blocks (EMPM/bin/EMPM_CMH_*, EMPM_CMB_*) for `tests/empm_cm_exact_check.py HRF_DIR
EMPM_DIR [--fit]` (reference: exact N_cm in the m-scheme FCI). 4He, 3 shells (`runs/4He_CMX`, copy of the `runs/4He` chain,
`run_He4_CMX_FULL.sh`, ~3 min): the hp phase (-1)^(J'-J), the normalizations & G2 reproduce N_cm of all FCI eigenstates to
1.3e-8 (`--fit`); the N_cm = 0 dimensions (with 3 phonons 1+/7-, with 4 phonons 196+/176- at M = 0, null eigenvalues
~1e-15, gap to 0.13 - 0.25) and the spectra agree with the projected FCI to 4e-10 MeV in every J^P. With 1 & 2 phonons no
state satisfies b_cm psi = 0 (the 3-shell HF space is not Nmax-complete), with 3 phonons only a few.

`tests/fci_sparse_cm.py HRF_DIR [EMPM_DIR]` is the sparse version of this reference (plain-operator H transformed to HO
determinants, N_cm block diagonal in the HO quanta, so its null space comes from small blocks; also the HF-basis
truncations <= n ph and the exact intrinsic spectra of the Nmax-complete HO spaces) ... it reproduces the dense 3-shell
results and the CM_Exact EMPM to 4e-10 MeV. 4He in 4 shells (`runs/4He_4sh`, `run_He4_4sh_chain.sh`, FCI M = 0 dimensions
39806/39492, ~30 min, results in `logs/He4_4sh_FCI.log` & `4He_4sh/fci_sparse_cm_4sh.npz`): 1159/1095 CM-free states, but
again only 1/7 with <= 3 ph excitations. The complete EMPM chain (`run_He4_4sh_chain.sh`, EMPM_4PH with `Ind_Greedy`,
EMPM_FULL with `CM_Exact`) reproduces this in every J^P (0+ .. 8+, 0- .. 7-): with 4 phonons the N_cm = 0 dimensions are
equal (0+ 110, 1+ 225, 2+ 301, 3+ 239, ..., 0- 99, 1- 224, 2- 288, 3- 234, ...; null eigenvalues ~1e-15, gap to
0.10/0.064) and the spectra agree to 1.3e-9 MeV (lowest 0+ -5.1079, 1- 20.3829 MeV relative to E_HF ... without CM
removal the lowest 1- is the spurious state at -5.45 MeV); with 3 phonons 0+ 1, 0- 2, 1- 3, 2- 2 states (= 1/7 at M = 0,
but N_cm only ~1e-5 with the next value ~3e-4), with <= 2 phonons none. EMPM_FULL 3:59 h, 85 GB on 32 threads (Hf of a
J^P up to ~1000 s, N_cm up to ~500 s per channel). For this size the within-level pp/hh blocks of L are obtained by a
contracted recursion (`EMPM_cm_dlevels`: level n from the level n-1 matrix & the phonon matrix, as in `EMPM_densityN`
at sigma = 1, instead of the full rank-1 densities of the 4-phonon states per J'), the 4-phonon part by batched dgemm
with X, and the Hf couplings 2<-1, 3<-2, 4<-3 by the threaded `mtm`.

**Weaker CM conditions (`CM_Variants = .true.` with `CM_Exact`).** N_cm is accumulated separately per phonon level q of
the intermediate rows (`EMPM_cm_dims`), so for every truncation n EMPM_FULL also imposes, after the exact condition
(rows q <= n+1): `space` P_n b_cm psi = 0 (rows q <= n) and `proj` P_(n-1) b_cm psi = 0 (rows q <= n-1; physical =
complement of P_n b+_cm acting on the <= n-1 phonon space, i.e. the SVD idea of `Ortho` with the full b_cm instead of its
ph part). Spectra & exact <N_cm> of the states go to `EMPM/bin/EMPM_FULL<n>_<JP>_E_<space|proj>.bin` & `_Ncm_<...>.bin`
(the exact output is unchanged). `tests/empm_cm_exact_check.py --variants` reproduces both to 6e-10 MeV (3 shells);
`tests/empm_cm_chain_check.py --out` gives the exact <N_cm> of the SVD (`Ortho`) states; `runs/plot_cm_svd_exact.py
IO_SVD IO_EXACT` plots no CM / SVD / proj / exact per J^P (`runs/plots/He4_{3,4}sh_cm_svd_exact*`). Findings for 4He:
SVD keeps ~10-17x the exact number of states at the top level, its lowest positive-parity states carry <N_cm> ~ 0.5
(4 shells: 1+, 2+, 3+ 4-6 MeV too low) and the lowest 3- (4 shells also 4-) is a CM excitation; it barely changes with n.
`proj` keeps the SVD dimensions at n <= 3 but its low states are nearly clean (<N_cm> 0.01-0.1) and close to the exact
4-phonon result already at n = 3 (4 shells, n = 3: 0+ -5.36, 1- 20.13, 2- 18.41 vs exact -5.11, 20.38, 18.72 MeV). At
the top level n = Full_N it does not remove the states whose b_cm psi lies in level n: most of its states then have
<N_cm> > 0.5 and a few with <N_cm> ~ 5-6 fall low (4 shells, 4 phonons: lowest 4+ at 30.5, lowest 3- at 33.8 MeV).
`space` coincided with the exact condition at n = Full_N as long as the rows of level n+1 were missing; with the
top-level completion of N_cm (below) `exact` uses the complete N_cm also there, and both remain nearly as restrictive.
(The numbers in this paragraph predate the completion.) Recommendation: `proj` for truncated chains, judged by the low states with small <N_cm>; the
exact condition where the space is complete enough to contain N_cm = 0 states.

**Projection `proj` as the CM treatment (`Full_Diag = 'proj'` or `'both'`, `CM_Ncm`; switch `CM_Proj` until 2026-11-06, called 'range' before 2026-11-06, the spurious space being the range of P_n b+_cm P_(n-1)).** `CM_Proj` imposes only the condition P_(n-1) b_cm psi = 0 (no exact /
space spectra, which leave almost no states in shell-truncated spaces: 4He NN only in 5 shells, 3 phonons, 0 - 4 states
per channel); the proj spectrum goes to `Full_proj_<JP>.dat` (format above). `CM_Ncm = .true.` (default) also builds the complete N_cm matrix: exact <N_cm> of the unprojected
(`Full_unproj_<JP>.dat`) & proj states (`Full_proj_<JP>.dat`), with `CM_Keep` and `Full_Svd` of the SVD states (`Full_svd_<JP>.dat`, H on the
physical states of each level). `CM_Ncm = .false.` is proj only: P_(n-1) b_cm psi needs no b_cm within the top level
Full_N, so `EMPM_cm_dlevels` stops at level Full_N - 1 and `EMPM_cm_blocks` skips the top pp/hh block (the dominant cost:
5 shells, 3 phonons, a dense 74352^2 matrix, 19 min). 4He 4 shells, Full_N = 3: `CM_Proj` with & without `CM_Ncm` and the
new `CM_Exact` + `CM_Variants` reproduce the previous binary bit for bit (proj, space, exact, <N_cm>; `runs/4He_4sh_rtest`);
70 s (exact + variants) -> 41 s (proj + <N_cm>) -> 11 s (proj only); 4He NN only, 5 shells, Full_N = 3
(`runs/4He_NN5sh_rtest`): proj only 11 min & 38 GB (b_cm levels 0.3 s) against 47 min & 81 GB with exact + variants,
proj spectra identical (0 difference); SVD 2.7 min & 33 GB. Benchmark in 5 shells (`runs/4He_NN5sh`, NCSM <= 3p3h
in the HF basis = the same space): proj within 0.1 - 1.3 MeV of the NCSM Lawson (beta 50) levels with <N_cm> 0.005 - 0.07
for the low states; SVD (`runs/4He_NN5sh_SVD`, `CM_Keep` run) keeps as many states but its low states carry <N_cm> 0.13 -
1.5 and lie up to 8 MeV too low; the top-level caveat above applies (lowest proj 3-, 4- with <N_cm> ~ 3.8). Notes:
`backup_2026-11-06/doc/empm_cm_conditions.pdf`.

**Exact <N_cm> on the top level (top-level completion, with `CM_Ncm`, Full_N <= 3).** The rows of level Full_N + 1 are not
in the basis, but only the ph part L_ph of L reaches them from level n = Full_N, and
sum_mu L_ph,mu+ L_ph,mu = sum_mu L_ph,mu L_ph,mu+ + C, with C = sum_mu [L_ph,mu+, L_ph,mu] a scalar 1-body operator (pp & hh
parts only): c(h,h') = sum_p <p||L||h><p||L||h'>/(2j_h+1), c(p',p) = - sum_h <p'||L||h><p||L||h>/(2j_p+1). L_ph+ maps level
n to level n-1, so the missing part on level n is G2/(A(2J+1)) [sum_J' G G^T + (2J+1) C0 + sqrt(2J+1) <||C_d||>], with
G = <level n, J || L_ph || level n-1, J'> (`EMPM_cm_up`, the ph block with rows & columns exchanged), C0 = <HF|C|HF>, and the
rank-0 within-level matrices <||C_d||> from the density recursion of `EMPM_cm_dlevels` at rank 0 (`Rk = 0`, pair weights
wC = <i1||c||i2>) (`EMPM_cm_top`). This is exact when level n-1 is complete, the same assumption the EMPM makes for H.
It costs one rank-0 pass over the levels (one J' instead of three, cheaper than the rank-1 top-level block). With
`NJ_CMTOP_CHECK` set (and the proj diagonalization), the completion of level Full_N - 1 is compared with its exact part from the rows of
level Full_N. Checks on 4He NN, 4 shells, complete sets (`backup_2026-11-06/runs_cmtop/4He_NN4sh_cmtop`, `tests/check_ncm_ncsm.py`): level-2
completion vs the exact rows at 2e-15 in all channels; <N_cm> of all unprojected states up to 1 / 2 / 3 phonons vs NCSM
<= 1p1h / 2p2h / 3p3h (HF basis, Beta_CM = 0): max |dN_cm| 5e-6, the print precision (33 / 103 / 96 states). This replaces
the earlier lower bound; with it, the 'lawson' variant of `CM_Beta` is also exact at n = Full_N.
For Full_N = 4 the completion is not implemented; the warning about it is printed only if a level 5 can exist (more than
4 hole m-states). In 4He (4 hole m-states) no 5p5h states exist, so <N_cm> is exact at Full_N = 4 without it: 4 shells, up to
4 phonons, all channels 0..4 +- vs NCSM <= 4p4h, max |dN_cm| 5e-6 (`runs/4He_NN4sh`).

**Lawson term in EMPM_FULL (`CM_Beta` > 0, MeV).** With the complete N_cm matrix (built as for `CM_Ncm`) H + beta N_cm is
diagonalized on the <= n phonon space ('lawson') and, with `Full_Diag = 'proj'/'both'`, on the proj space ('projL'); E = <H> =
eigenvalue - beta <N_cm> & <N_cm> in the order of the eigenvalues go to `_E_<lab>.bin` & `_Ncm_<lab>.bin`. For n < Full_N
the N_cm matrix is exact and the result equals the NCSM with the same beta in the same space: 4He NN only, 5 shells,
Full_N = 3, beta 50 (`runs/4He_NN5sh_rtest/He4_FULL3L.nml`, `runs/check_empm_lawson.py`): n = 2 vs NCSM <= 2p2h (HF basis)
max |dE| 4.9e-9 MeV, |d<N_cm>| 4.8e-7 for the lowest 3 states of all J^P <= 4; at n = Full_N = 3 these runs predate the
top-level completion (rows of level 4 missing, N_cm a lower bound): the EMPM levels lay 2 - 110 keV below the NCSM <= 3p3h
ones (<N_cm> 0.003 - 0.009 lower). With the completion N_cm is exact at n = Full_N <= 3 as well.
'projL' is cleaner (<N_cm> 0.002 - 0.02) but 0.5 - 4 MeV higher; unlike 'proj' its lowest 3- / 4- are physical
(<N_cm> 0.01 - 0.02 instead of 3.8). 29 min, 81 GB (proj + <N_cm> + Lawson, 24 threads).

**Lawson term in the interaction (`&calculation CM_Lawson` > 0, MeV).** beta N_cm is added to V_NN when the interaction
is read (`V2b_read`): the pair operator beta/A [b+(1).b(2) + b+(2).b(1) + (N_1 + N_2)/(A-1)], the 1-body part folded in
by 1/(A-1) (exact for A particles). The whole chain (HF -> EMPM -> NCSM 'EMPM') then works with H + beta N_cm, so the
<= n phonon EMPM diagonalizes P_n (H + beta N_cm) P_n exactly, without the b_cm rows of level n + 1 and with no CM step
in EMPM_FULL; all energies include beta <N_cm>, the NCSM run in the same space (Beta_CM = 0) gives <N_cm> and <H> =
E - beta <N_cm>. The HF reference is that of H + beta N_cm (a different space than the HF of H). Checks: NCSM 'NN' in
the HO basis, Nmax 4/3, CM_Lawson = 50 equals Beta_CM = 50 (40 states, all digits); 4He NN 4 shells, EMPM up to 2 / 3
phonons vs NCSM <= 2p2h / 3p3h in the HF basis: max |dE| 5e-7 MeV (print precision) for 101 / 96 states
(`runs/4He_lawV_4sh`, `runs/check_lawV.py`).

**CM-clean phonons with the CM condition in the full Hamiltonian (`Ortho = .true.`, `CM_Keep = .true.` in `EMPM_FULL`).**
The levels are built with the CM-clean TDA phonons of `Ortho = .true.` (cm phonon = normalized E1 vector, the other 1-
phonons in its complement) and their complete physical + spurious sets; `CM_Keep` keeps the cm phonon & the spurious
n-phonon states in the `EMPM_FULL` basis with the within-level couplings R(nu) = <nu|A|cm> and Rc = <physical|H|spurious>,
so that `CM_Exact` / `CM_Variants` (proj P_(n-1) b_cm psi = 0 with the full b_cm, ph + hp + pp + hh) act on the complete
sets in this basis, instead of the level-by-level SVD projection with the ph part of b_cm only. With complete sets this is
a rotation of the `Ortho = .false.` basis: 4He in 3 shells (`runs/4He_CMK`, the `4He_CM` sets, Full_N = 4) reproduces
`runs/4He_CMX` (none, exact, space, proj) in all 54 channel x truncation spectra with equal dimensions to 7e-12 MeV.
The point of the clean phonons is the truncated case: the spurious strength sits in the cm phonon, so the |a cm>
products (the bulk of b+_cm |a>) can be kept completely in an energy-truncated basis (spurious rows always kept).

### Chain TDA -> 2PH -> 3PH -> 4PH -> FULL

Each step reads the states saved by the previous one: `EMPM_2PH` with `Save_X = .true.` writes the 2-phonon states (basis,
E, C on the full basis, X = sqrt(2J+1) D(independent rows)^T C) to `EMPM/bin/EMPM_S2_<JP>.bin` & `EMPM/2Phonon/States2.dat`
(off by default ... large for big spaces), `EMPM_3PH` reads them (it no longer solves the 2-phonon problem itself) and
writes `EMPM/bin/EMPM_S3_<JP>.bin`, `EMPM_4PH` reads S2 & S3 and writes `EMPM/bin/EMPM_S4_<JP>.bin`. Example: `runs/run_He4_HF.sh`,
`run_He4_EMPM2.sh`, `run_He4_EMPM3.sh`, `run_He4_EMPM4.sh`, `run_He4_FULL.sh` (or `run_He4_chain.sh`), 4He in 3 shells,
~3 min on 28 threads.

### Full Hamiltonian (`EMPM_FULL`), sequentially with up to 1 .. 4 phonons

From the saved state sets, per J^P the basis [ground (0+) ; 1- ; 2- ; 3- ; 4-phonon states], H built once up to `Full_N`
(1 .. 4) phonons (channels `JPF`) and its leading [0 ; .. ; n] block diagonalized for n = 1 .. `Full_N`. Diagonal the
energies of each subspace, couplings
* n+1 <- n (1-2, 2-3, 3-4): H = X^T V, V(alpha lambda; a) = (-1)^(Ja+J_alpha+J_lambda)/(2Ja+1) sum_s rho_(a alpha)(J_lambda; s)
  W(s,lambda) with the n-phonon densities of a (as `vintn12` & `vintn23`), W(s,lambda) = sum kappa c^lambda F^(J_lambda),
* n+2 <- n (0-2, 1-3, 2-4): only the pair-creation part V20 of H, which commutes with the phonon operators ... G_0 = h20 =
  <gamma|H|0> of the 0+ 2-phonon states, G_n(beta; alpha) = 1/sqrt(2J+1) sum_(k = mu nu) C^alpha_k sum_d G_(n-1)(d; mu)
  X^beta(d nu) (n = 1 is the h31 of the original `nondiag13`).

**Diagonalizations & outputs (`Full_Diag`).** `Full_Diag = 'unproj'` (default), `'proj'` or `'both'` selects the full
diagonalization per truncation n = 1 .. Full_N: H on the complete <= n phonon basis without CM projection, and/or the CM
projection P_(n-1) b_cm psi = 0 (below; needs complete sets). All variants are written with the same structure to
`EMPM/Full/<n>Phonon/Full_<lab>_<JP>.dat` (lab = `unproj`, `proj`; `svd` with `CM_Keep` and `Full_Svd = .true.` (default
.false.), or as the unprojected spectrum of an `Ortho` basis without `CM_Keep`; the diagnostics `exact`, `space`, `lawson`, `projL`): a header line (J^P,
method, dimensions n0 .. n4, number of states & of physical states) and the columns i, E (relative to E_HF), exact <N_cm>
(-1 if not computed), phys (<N_cm> < 0.5), w0 .. w4 (0- .. 4-phonon weights). <N_cm> is computed (`CM_Ncm = .true.`,
default) whenever the sets are complete. `Summary_<lab>.dat` (unproj or svd, proj) lists per channel the dimensions, the
numbers of states & physical states and the lowest physical state; binaries `EMPM/bin/EMPM_FULL<n>_<JP>_E_<lab>.bin`,
`_Ncm_<lab>.bin`, and `_Z_<lab>.bin` (n0 .. n4, eigenvectors) for unproj & proj; the 5 lowest states in the log.
`tests/proj_physical.py IO_DIR` collects the physical proj states of all J^P (`Full_proj_physical.dat`).
`runs/plot_empm_spectra.py` (with `plot_He4_spectra.sh`) plots the spectra of any subspace sets side by side (gnuplot).
Complete 4He (3 shells): with 4 phonons the whole 4-nucleon space, all 6040 eigenvalues (3060 +, 2980 -) agree with an
independent M-scheme FCI of the same normal-ordered H to 6e-10 MeV in every J^P channel (the 1e-10 rounding of `V2B_HF.dat`):
`tests/empm_fci_check.py HRF_DIR EMPM_DIR` (HF energies from `Bin/h_HF.bin`, `V2B_HF.dat` & `Orbitals_HF.dat` of an HF run with
`&hf HRF = .true.`, e.g. `runs/4He/He4_HF_HRF.nml`, J from <J^2>, ~2 min), FCI spectrum in `HRF_DIR/FCI_Spectrum.dat`.
These text files come from the HRF output of the base-version HF (`HF_LEGACY` with `HRF = .true.`); the dev HF writes no
HRF tables (`HRF` is a dead input), so the check needs `HF_LEGACY` until it is ported to `HF/bin/*.bin`.
The pn blocks of `V2B_HF.dat` & of the full N x N pn matrices of `Bin/V2B_*.bin` were written with the packed index of the
lower triangle also for Ket > Bra (wrong upper triangle) until Nov 2026; the solvers were not affected, the import now reads
the lower triangle only, so older files remain valid.

**Memory of the state sets (2026-11-10).** `EMPM_states_read` has two optional switches: `No_Xt` (no transposed copy Xt,
which only the row construction of EMPM_nPH uses) & `Lazy` (X not read, only its stream position; `EMPM_states_load_X` /
`EMPM_states_free_X` load & release X of one channel). EMPM_FULL reads the 2- & 3-phonon sets with `No_Xt` and the
4-phonon set with `Lazy`: X(J,P) is loaded before the 4 <- 3 coupling & released after the N_cm loop, X(Jt,Pt) of the
intermediate channel only around `EMPM_cm_blocks` (one at a time), each load logged with its size & time ... at most two
4-phonon channels in memory instead of all of them twice. Results unchanged: 4He 3 shells, FULL4 all J^P with
`Full_Diag = 'both'` & `CM_Ncm`, all outputs bit-identical, memory 0.7 vs 1.3 GB. 4He NN only, 5 shells (S4 3.5 TB on
disk): FULL4 had run out of memory (1.5 TB) while importing the states; with `Lazy` the peak is ~1.3 TB (3+/3-: two X of
466 GB) and a J^P channel takes up to ~5.5 h on 96 threads (N_cm & b_cm parts up to 4.7 h, H 1 h, diagonalizations 0.5 h).

**EMPM_FULL with Full_N <= 2 without the phonon cache (2026-11-12).** The cache of `EMPM_cache1_make` (R, Rt & V of every
phonon for every rank) is built only for Full_N >= 3. For Full_N <= 2 the within-level CM matrices (ranks 1 & 0 = the b_cm
and the commutator parts) come from one phonon density at a time (`EMPM_cm_dlevels_nc`, both ranks in one pass, optional
`Sel` = [J_a, P_a, J_t, P_t] for one block only), and the 2-phonon set is read `Lazy`: X(J,P) of the channel is loaded for
the 2 <- 1 coupling and X(Jt,Pt) of the intermediate channels around `EMPM_cm_blocks`. The log has RSS / peak lines
(`memory ... RSS, peak`). 4He 4 shells (Full_N 1, 2, 3; complete and E2max 60): all outputs byte-identical. 208Pb window
(13 shells, E2max 20, 2-phonon X of 311 GB on disk): 0+ & 1- with up to 2 phonons, peak 246 GB in 2 h 14 min on 16 threads
(with the cache it ran out of memory above 600 GB).

### E1 transitions of the full EMPM states (`EMPM_TRANS`)

On top of an `EMPM_FULL` run with `Full_N` = 1 or 2 (same input, `Solvers = 'EMPM_TRANS'`; complete sets, i.e. `Ortho =
.false.` or `Ortho` with `CM_Keep`). The E1 operator of protons and neutrons (1-body, HF basis) goes through the CM-operator
machinery (`EMPM_op_make`: ph strengths bL = bR, pair weights wD; `EMPM_cm_dlevels_nc` with `Sel` for 0+ -> 1-;
`EMPM_cm_blocks`), giving <1- basis||Q_t||0+ basis> on [0 ; 1-phonon ; 2-phonon]. With the eigenvectors of
`EMPM/bin/EMPM_FULL<n>_<JP>_Z_<lab>.bin` (lab = unproj or svd, proj; n = 1 .. Full_N):
M_t(k) = Z_1-(:,k)^T B_t Z_0+(:,gs), gs = the lowest 0+ with <N_cm> < 0.5; B_ph = M_p^2, B_is = (M_p + M_n)^2/4, B_iv =
(e_p M_p - e_n M_n)^2 (e_p = N/A, e_n = Z/A), as `EMPM_TDA` & `NCSM`. Output `EMPM/transitions/Full<n>_<lab>_E1.dat`
(k, E - E_HF, Ex, N_cm, physical flag, B_ph, B_is, B_iv; header with the g.s., the sums & the B_iv centroid); for n = 2 also
from |HF> (the 0+ basis state 1, Ex = E - E_HF): `Full2_<lab>_E1_fromHF.dat`. Only X of 0+ & 1- are loaded. Checks (4He
4 shells, complete): unprojected <= 1 / <= 2 phonons equal NCSM <= 1p1h / <= 2p2h (Transitions) in all 8 / 126 1- states
(energies exactly, B to 5e-9); from |HF> the sum over all <= 2-phonon 1- states equals the <= 1 sum (HF sum rule).

### E1 strength function with up to 3 phonons by Lanczos (`EMPM_SF`)

B(E1; 0+_gs -> 1-) with the g.s. up to 2 phonons and the 1- space up to 3 phonons, without the 3-phonon eigenproblem and
without the full diagonalization (`EMPM_SF` in `nj_empm_driver.f90`). It runs on top of `EMPM_TDA`, `EMPM_2PH` (all channels,
`Save_X = .true.`) and `EMPM_FULL` with `Full_N = 2`, `JPF = '0+'` (the g.s. = lowest 0+ of `EMPM/bin/EMPM_FULL2_0+_Z_unproj.bin`);
`Ortho = .false.`. Input (`&empm`): `SF_Steps` (100), `SF_Op` = `'iv'` (default), `'is'`, `'ph'`, `SF_Reorth` (.true.),
`SF_Mode`. The 1- space is [1-phonon ; 2-phonon states ; 3-phonon sector], the couplings are those of `EMPM_FULL` taken
before the projection on the 3-phonon states (H|a_2> -> sqrt(2J+1) V(:,a), H|lambda> -> h20(gamma) on |gamma lambda>,
Q_ph|a_2; 0+> -> bL(lambda) on |a lambda>), the [1 ; 2] part of Q|g> is that of `EMPM_TRANS`.

* `SF_Mode = 'full'` (default): the 3-phonon sector is the whole redundant basis |a lambda; 1-> (all 2-phonon states x
  phonons). A vector is kept by its coefficients c; the equations of motion give H|i> = sum_k A(i,k)|k>, so H acts as
  c -> A^T c and the scalar product is c^T D d. Neither D^-1 nor a Cholesky factor is needed. With complete sets D has only
  the eigenvalues 0, 5 and 9 (sum_i |i><i| = N_p(pi) N_h(pi) + N_p(nu) N_h(nu) on 3p3h, the phonons being
  charge-conserving), so P = D (14 - D)/45 projects on its range; the null-space part of c, which the metric does not see
  but which grows by a factor ~2.4 per step (4He: breakdown after 29 steps without it), is removed by c -> P c in every
  step. One product with A^T and two with D per step. Complete sets only (no `E2max`, `E3max`); A & D are dense
  (16 N3^2 bytes), so this mode is a reference implementation for small spaces.
* `SF_Mode = 'ind'`: the independent rows I of the 3-phonon basis as in `EMPM_3PH` (pivoted Cholesky of D, or
  `Ind_Greedy`; rows `E3max`), H_II = A_I D_I^T and D_II = U^T U; Lanczos on the orthonormal states with
  H~ = U^-T H_II U^-1 and the couplings U^-T D_I g. Also for truncated bases, with the usual EMPM approximations.

Outputs `EMPM/SF/SF_E1_<op>.dat` (poles: E - E_HF, Ex, B, residual), `_steps.dat` (poles after every step) and `_coef.dat`
(alpha_j, beta_j), in the format of the NCSM strength function. The log has a Hermiticity check (<x|H y> = <H x|y> on
random vectors) and the check of the D spectrum. Checks (`runs/empm_sf/`, `compare_coef.py`): 4He in 4 shells (1- space
8 + 118 + 12030 redundant, 1266 independent) and 16O in 3 shells (10 + 190 + 28832, 4552 independent; A & D 13.3 GB)
against the NCSM (g.s. <= 2p2h, 1- <= 3p3h, `Single = .false.`): S0 equal to 4e-15, alpha_1 and beta_1 to 1e-14, all 40
coefficients to 2e-6 or better (both modes; the difference grows with the step as usual for Lanczos). Truncated
(4He, `E3max = 45`, 46 independent rows): the Lanczos poles equal the 1- eigenvalues of `EMPM_3PH` + `EMPM_FULL`
(`Full_N = 3`) to the 6 printed digits.

## NCSM (`src/ncsm/`, action `NCSM`)

No-core shell model in the m-scheme HO basis with the same matrix elements, as a reference for the EMPM (light nuclei;
4He to Nmax ~ 14, dimension ~ 1e7 with the stored matrix, larger on the fly). `&ncsm`: `Nmax_List` (many-body
truncations, HO quanta above the lowest configuration, parity (-1)^(N0+Nmax), -1 terminates), `Hamiltonian`, `M2` (2M,
default A mod 2), `N_Eig`, `Lanczos_Max` (largest number of matrix-vector products), `Lanczos_Tol` (residual),
`Beta_CM` (Lawson term H + Beta_CM N_cm, MeV per CM quantum, 0 = none), `Matrix` (`'auto'` default: stored if within
`Matrix_GB`, else on the fly; `'stored'`; `'otf'`), `Single` (stored values in single precision, energies to ~1e-6
MeV), `Lanczos_Basis` (thick-restart basis, 0 = max(2 N_Eig + 50, 100)), `Matrix_GB` (0 = 70% of MemAvailable),
`Trunc` (`'Nmax'` default; `'npnh'`: all m-states e <= Calc%Nmax, at most `Nmax_List(i)` particles outside the lowest
A m-states, both parities, output `NCSM_<Ham>_<Basis>_<n>p<n>h<p|m>[_B<beta>].dat` ... the EMPM n-phonon space),
`Basis` (`'HO'` default; `'HF'` with `'EMPM'`: H kept in the HF orbitals, N_cm built from b & N transformed with C,
so that `npnh` counts excitations w.r.t. the HF reference like the EMPM). Checks (4He NN only, 5 shells,
`runs/4He_NN5sh`): the complete space (4p4h) in HO & HF gives identical E, <J^2> & <N_cm>; the 3p3h M = 0 dimensions
39209 / 38528 equal the EMPM 0..3-phonon J-scheme counts summed over J.
Hamiltonians: `'NN'` T_int + V_NN in the HO basis ((1 - 1/A) T + V_NN with the 2-body CMS term of `V2b_read`, no HF run
needed, translationally invariant) and `'EMPM'` the normal-ordered EMPM Hamiltonian E_HF + e :a+a: + 1/4 H2b :a+a+aa:
of the HF reference `io/<Path>/<&empm Reference>/bin/{orb,H1b,H2b,C,rho}.bin` (default `io/<Path>/HF/bin`; E_HF from
its `HF_summary.dat`; H1b diagonal, residual 3-body dropped; its plain form is built in the HF basis & transformed back
to HO with C), energies also relative to E_HF. As for the EMPM the reference must be HF-like (diagonal `H1b.bin`, so
`'HF-MBPT'` is refused) with integer occupations (diagonal of `rho.bin` 0 or 1 within 1e-8, else fatal: closed
(sub)shells or `&hf ESH = .true.`). With `'EMPM'` or `Basis = 'HF'` the orbitals come from the reference `orb.bin`,
otherwise from the LHO construction (Calc.Nmax). Basis: products of proton & neutron SDs
(single-particle space e <= min(Calc%Nmax, e_F + Nmax), a warning if Calc%Nmax cuts the Nmax space), m-scheme 2-body
tables per pair block (2M, parity) within the pair quanta cut (needs N2max >= the largest pair quanta). Engine (v2):
one operator H + Beta_CM N_cm (tables combined), each column generated from the like-particle species matrices & the
pn one-body jumps; `stored` keeps the lower half by columns (64-bit pointers, 8 or 12 bytes per entry, product with
thread-private partial sums), `otf` contracts the generated columns with the vector(s) directly (no matrix). Thick-restart
Lanczos (Wu & Simon: full reorthogonalization in a basis of `Lanczos_Basis` vectors, the lowest Ritz vectors kept at a
restart), so many states (e.g. hundreds at beta 0 among the CM copies) need bounded memory. <J^2> & <N_cm> (b_mu =
sqrt(2) (r_mu/b)_(N->N-1) as in `tests/empm_cm_check.py`) of the final states on the fly in blocks of 16 vectors; E =
<H> = Ritz value - Beta <N_cm>. Output per state E, J (from <J^2>), <N_cm> & the residual in
`io/<Path>/NCSM/NCSM_<Ham>_Nmax<N>[_B<beta>].dat`. Checks: 4He in 4 shells with the EMPM Hamiltonian
(`runs/4He_4sh/He4_NCSM.nml`, Nmax 12 / 11 = the complete 4-shell space, M = 0 dimensions 39806 / 39492) reproduces
`tests/fci_sparse_cm.py` (all J^P, incl. <N_cm>) to 4.7e-9 MeV; with `'NN'` (`runs/4He_NCSM`, Nmax 2 .. 10) every <N_cm>
is an integer to 1e-16 (exact CM factorization) and the intrinsic levels of Nmax reappear unchanged as CM excitations
at Nmax + 2, ...; ground state -21.597, -26.086, -27.630, -28.190, -28.328 MeV for Nmax 2, 4, 6, 8, 10. The v2 engine
reproduces the v1 results (full matrices, plain Lanczos) of `runs/4He_NCSMsh` (Nmax 6/5 & 10/9, beta 50: all 30 states,
E, <J^2> & <N_cm> identical to print precision in `stored` & `otf`, `Single` within 1e-6 MeV); Nmax 10/9 (dimensions
196438 / 96340, 1.35e8 / 4.8e7 entries in the lower half): 2.0 min & 3.5 GB on 16 threads (v1: 4.3 min & 7.6 GB).
SD hash keys (2026-11-09): combinatorial rank sum_i C(occ_i - 1, i) of the sorted occupations (< C(ns,Nt)) instead of
the base-ns number (< ns^Nt, which overflowed int64 for 16O at Nmax >= 6, ns = 240); results unchanged (16O Nmax 2 & 4
energies bit-identical). Remaining limit: one-body jump offsets (Jo, Tot) are 32-bit, 16O Nmax 8 needs 3.2e9 jumps per
species. 16O runs: `runs/16O_NCSM/README_16O.md` (NCSM), `runs/16O_EMPM5sh/README_16O_EMPM.md` (EMPM, 5 shells, <= 2 phonons).

**npnh options (2026-11-11/12).** `Ecut_npnh` (Basis 'HF'): cuts (MeV) on the unperturbed HF excitation energy of the
n-particle n-hole part, one value per n (< 0: no cut, the last value repeats); output tag `_Ec<E>` (one value) or
`_Ecx-<E2>[-<E3>]` (per rank, x = rank 1 uncut). `Window = .true.`: only the `&empm` window (pHole, pPart, nHole, nPart)
is active, the SDs are stored in hole/particle form w.r.t. the reference, the frozen core stays in E_HF & e (so the
reference is E_HF exactly); <N_cm> covers the window only and `Beta_CM` must be 0. `Parity` (-1 both, 0 +, 1 -).
`Dense_Max = n`: full LAPACK diagonalization (all states) when the dimension <= n.

**Transitions (`Transitions = .true.`, `Trunc = 'npnh'`, `M2 = 0`).** B(EL; g.s. -> k) for E0 (r^2), E1, E2, E3 from the
lowest parity + state (`Tr_From = 'gs'`, kept for the parity - run) or from the reference (`Tr_From = 'HF'`): output
`<ncsm file>_E<L>.dat` (k, E, Ex, J, B_ph, B_is, B_iv, M_p, M_n) & the sum rule (2L+1)|Q_p g|^2 in the header (B_ph
protons with bare charge, B_is, B_iv with the recoil charges N/A, Z/A). `Tr_Lanczos = n`: strength-function Lanczos of
n steps from Q|g> for B_ph, B_is, B_iv (`<file>_E<L>_SF_{ph,is,iv}.dat`: k, E, Ex, B, residual; header E_g, dim, S0).
`Tr_Ops = 'E1iv'` (list of `E<L><ph|is|iv>`, default '' = all): only the listed strength functions (6x fewer products of H
for one of them). `Tr_Lanczos = 1` already gives m0 & the centroid m1/m0 exactly, with the width sigma as the residual.
`Tr_Gs_Nmax = m` (with e.g. `Nmax_List = 2, 3`): the g.s. of the <= m p m h run is kept (by SD occupations, remapped into
the larger basis) and used as the initial state of the larger run (g.s. <= n p n h -> excited states <= (n+1)p(n+1)h).
`Tr_Ckpt = k` (default 0 = off): checkpoint & restart of the strength-function Lanczos. Every step the Lanczos vector is
written to `<file>_E<L>_SF_<op>_vecs.bin` (n x 4 bytes with `Tr_Vec = 'single'`, n x 8 with 'double'; only with
`Tr_Reorth`), every k steps alpha, beta & the current vectors to `<file>_E<L>_SF_<op>_ckpt.bin` (via a .tmp file &
rename, so a run killed while writing keeps the previous one; about 2 x 8n bytes). A run with the same input that finds a
valid checkpoint (same dimension, operator & S0) continues after its last step, with `_steps.dat` rebuilt; a larger
`Tr_Lanczos` extends a finished run (from the checkpoint at step <= Tr_Lanczos - 1). Delete both files to start from step 1.
The setup (basis, species matrices, g.s.) is redone on restart. Check: 16O 3 shells, 40 steps in one run vs 15 + restart,
coefficients equal to the run-to-run noise (alpha 7e-8 MeV, beta 6e-12 MeV), single & double vectors.
`EMPM_TDA` writes `EMPM/transitions/TDA_E<L>.dat` (E, B_ph, B_is, B_iv as the TDA). Checks: 16O 4 shells, <= 1p1h vs
TDA (`Ortho = .false.`) all E0 .. E3 states, energies exactly & B to 4e-9 (`tests/cmp_ncsm_tda_tr.py`); 208Pb
window (NN + NNN, 13 shells), NCSM <= 1p1h vs EMPM_TDA equal (B to 3e-9); 16O g.s. <= 2p2h with 1- <= 3p3h agrees with
the full <= 3p3h to 0.05 MeV & 0.03 % in B; SF Lanczos = dense.

## Compatibility with the Julia code

* Byte-compatible with NuJuliet dev f3f7a33: binary files (`io/*/HF/bin/*.bin`, `io/*/HF-MBPT/bin/*.bin`,
  `io/*/<Name>/bin/*.bin`, NuHamil interaction files) are raw little-endian streams, `orb.bin` with the int64 labels
  `a, n, l, 2j` and the float64 occupations per orbital. The port reads Julia files, so Julia and Fortran stages can be
  mixed, with one exception: the pn block of a two-body file (`H2b.bin`, `sigma.bin`) is written by Julia as the full
  N x N square with the packed lower-triangle index (aliased upper triangle), by the port as the true symmetric matrix,
  so Julia cannot read a pn block written by the port (Fortran HF -> Julia TDA/RPA).
* Text outputs have the same names, layout and number formats.
* Results agree with the Julia code to ~1e-12 or better (1e-8 in the test harness); eigenvector-based quantities can
  differ by arbitrary signs/phases (LAPACK gauge) without changing any physical output.
* EMPM and NCSM have no Julia counterpart; their file conventions are described above (`&empm Reference`, `EMPM/bin`).
* Intentional deviations, bug fixes and Julia quirks are listed in `docs/PORTING_NOTES.md`; the state of the dev port is
  in `docs/dev_port/STATUS.md`.

## Verification

`tests/dev_port/` (see its `README.md`): Julia reference runs of the dev code (`run_julia.sh`), the same cases with the
port (`run_fortran.sh`, namelists in `tests/dev_port/fortran/<case>/`) and `compare.py JULIA_IO FORTRAN_IO [--tol 1e-8]`
(walks the lowercase `io/<Path>/...` trees: binary files element-wise and modulo column signs, text files token-wise with
a numeric tolerance; `gauge_check.py` for eigenvector gauges). EMPM / NCSM: `run_empm_ncsm.sh baseline|dev|compare`
against a legacy baseline of the frozen stage-4 binary (`empm_compare.py`).
