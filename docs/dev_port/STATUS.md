
## 2026-11-20 — stage 1, verification harness

Done:
- Patched Julia reference: worktree `/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref`, local branch
  `port-ref` = f3f7a33 + commits F1 (204b5be), F2 (95bc34d), F3 (60d650a), F4 (990ef85), F5 (55c5880); reviewed, each
  minimal. F5 note: `Orb1B` has a single (l,j), so the neutron (l,j) from `nC` is checked against the proton one and
  a mismatch is an `error` (Fortran must `fatal()` likewise).
- Runs with Julia 1.10.11 and the user's `@v1.10` env (all packages present, no dedicated project needed).
- `tests/dev_port/`: drivers `julia/{common,bundled,ca40_5sh,o16_open}.jl`, `run_julia.sh`, `compare.py`, `README.md`;
  `/tests/dev_port/work/` added to `.gitignore`.
- Reference outputs: `tests/dev_port/work/julia_t8/{bundled,ca40_5sh,o16_open}/io`, `work/julia_t1/bundled/io`.
  Wall times: bundled 3:26 (t8) / 4:50 (t1), ca40_5sh 2:55 (t8), o16_open 1:25 (t8).
- Thread-race check bundled t1 vs t8: 409 ok, 38 SIGN, 0 DIFF at 1e-8 (only eigenvector/orbital sign gauges).

Open items:
- Fortran namelist inputs for the three cases (solver stages); `compare.py` cannot remove the HF-orbital sign gauge
  (files in the HF basis come out as SIGN when an orbital flips; see README).
- The interactive-thread stdout flusher occasionally duplicates a few log lines in run.log (concurrent flush);
  cosmetic only.

Not ported (dead Julia code, to be listed by the solver stages): 9-arg H1b with Sigma, E1b, E2b, RPA_OBDM_sigma,
dead RPA_spurious Sigma variant.

## 2026-11-20 — stage 1, group common A (core types, parameters, orbitals, operators)

Done (spec `common_a.md`; resumed from an interrupted attempt, whose edits were reviewed and kept):
- `nj_types`: real `orb1b_t%pO/nO`; `ind2b_t`/`ind2b_vec_t` (ex `ind2qp_*`), `orb2b_2qp_t` (ex `qporb2b_t`, `i`→`Ind`,
  dense `Dic`), new `orb2b_1p1h_t`, `pnorbvector_t`, `orb1b_hf_mbpt_t`, `o2b_pandya_t`.
- `nj_params`/`nj_input`: dev structs as namelists (`&hf` +V_HF/ESH, new `&hf_mbpt`, `&tda`, `&pnqtda`, 17-field `&rpa`,
  `&qtda`/`&qrpa` +Reference, `&hfb` dev Broyden fields); one contained reader per group; validation of
  Transitions, Name/Reference (non-empty, no '/', Name /= Reference). `io/` root, `io_path`, new `io_dir`.
- `nj_orbitals`: dev `orbitals_make(Params[,Use_Int])` (Calc.Nmax default, fractional filling, fatal instead of the
  infinite loop), `orbitals_update`, binary `orbitals_export(Params,Orb,Path)`/`orbitals_import`, `orbitals_1p1h_make` +
  `orb2b_1p1h_index`, `orb2b_2qp_index`, `OBDM_natural_orbitals` with **F5**.
- `nj_o1b`: O1b helpers; **F4** (`T1b` sized by `size(Orb)`).
- `nj_o2b`: `O2b_temp` accessor, `O2b_Pandya_initialize`, Pandya accessors, 4-block `O2b_Pandya_transformation`
  (pnpn new) next to the LEGACY 3-block form, optional `Report` of `O2b_transformation`.
- `nj_obdm`: `rho_*` names (OBDM_* kept as generic aliases), dev radius normalization `sqrt(R2)/A^(1/3)`.
- `nj_functions`: `J1b_convective` dev fixes, `phonon_weights`; `nj_utils`: `parity_of`, `triangle`, `hat_of`;
  `nj_linalg::hungarian` square check.
- Mechanical: legacy solvers (hf, rpa, qp, empm, ncsm) call `orbitals_make(Params,Use_Int=.true.)` (their old Int.Nmax
  length); qp chain uses the renamed 2qp type; messages `IO/` → `io/`.
- Fortran namelists for the harness cases: `tests/dev_port/fortran/{bundled,ca40_5sh,o16_open}/*.nml` (run in file
  order) + `tests/dev_port/run_fortran.sh`. All 14 files parse with the new reader (actions not yet implemented
  are the solver stages' job; `&transition_operator_densities`, `&transition_analysis`, `&parse` are read by later stages).

Verification:
- Build (ifx) clean.
- Regression, old solvers, 40Ca 9 shells (`runs/test_Ca40_9sh/common_a`, 16 threads, 3:16 wall) vs `fortran_new/IO`:
  77 files, same file set (under `io/` instead of `IO/`), 74 identical (max err 0), 3 differ: `HF/HF_Summary.dat`,
  `RPA/RPA_Summary.dat`, `RRPA/RRPA_Summary.dat`, only in the "normalized by 1/A^1/3" radius block (dev formula
  `sqrt(R2)/A^(1/3)`). The E1 toroidal/vortical/compression files are unchanged: `J1b_convective` only enters the
  transition-current exports, which this test does not run.
- Unit comparison against the patched Julia reference: `orbitals_make` (20O, 18O, 23Na; Nmax 3) occupations and
  `orbitals_1p1h_make` channel sizes / ordering (both tolerance sets) bit-identical; `OBDM_natural_orbitals` on a
  synthetic 16O density: identical n,l,j, occupations and C diagonals within 6e-17.

Open items:
- Spec items left to other groups by design: V2b/V3b readers building their own `Use_Int` orbitals (common B), the
  solvers switching to `orbitals_1p1h_make`/`io_dir`/`rho_*`, `&empm Reference`, `export_request_t` M=-1/Cut and the
  new external namelist groups (consumers), qp use of `o2b_pandya_t`/pnpn (qp group).
- Cleanup stage: remove `hf_efa_params_t`, `rrpa_params_t`, `&hf_efa`/`&rrpa` readers, `Broy_Qmax/Bmax` dummies,
  `orbitals_ph_make`, `orbitals_ph_list`, `orbitals_one_phonon_make`, `read_h_HF_diag`, HRF `orbitals_export(Params,Orb)`,
  the 3-block Pandya form, `OBDM_*` aliases, `obdm_iteration_pairs_t`, `Use_Int=.true.` in the solvers.
- Not ported (dead/unreachable in dev, group common A): `w3j`, `radial_phase_test` (removed), `logm`, `lorentzian`,
  `integrate_quadrature`, `legendre_polynomial`, `pretty_summarysize`, `src/common/rO.jl` (not included),
  `O2b_initialize_empty` (an unallocated `o2b_t` plays its role), `orb2b_*_key` packing (dense lookups instead).

## 2026-11-20 — stage 1, group common B (interactions, transitions, Hamiltonian)

Done (spec `common_b.md`; details & deviations in `docs/PORTING_NOTES.md`, section "common B"):
- `nj_v2b`: `V2b_read(Params,V_NN,Orb_NN)` (own Int orbitals, 1e6-element block reads, stops after the last
  Calc chunk), non-mutating `V2b_residual_no2b(...,V_NN,V_NNN,V_NN_res)`, new `V2b_parse`, `V2b_correlation_energy`.
- `nj_v3b`: dense `Orb3B.Map`/`nlj` (hash map gone, `V3b_no2b_key` deleted), `V3b_no2b_lookup`, row-wise packed blocks
  `Ind = Ket + Bra(Bra-1)/2` with 64-bit indices/lengths, `V3b_no2b_read(Params,V_NNN,Orb_NNN)` (own Int orbitals, only
  needed chunks, non-Calc (c,d,e,f) blocks skipped without reading), new `V3b_no2b_parse`.
- `nj_tr1b`: F & GT (`Tr1b_BX_initialize`, pn-mixing transformation, `qp_transform_bx`), upstream E1 convective /
  spin-normalization fix, 0.5 on VC/TC/sTC, dev qp 11/20 forms (O_lk, no E0 c-number, time-odd Phase_T).
- New `src/core/nj_hamiltonian.f90` (Makefile SRC_CORE after nj_v2b): `H1b`, `H2b`, `H1b_H2b`.
- `nj_types`: `i3_t`, `orb3b_t` = (Map, N, nlj); `nj_utils`: generic `istr` (int64), `abspath_of`, `presize_file`;
  `nj_o2b`: `O2b_move`. App: `&parse` group, `&run` actions `V2b_parse`, `V3b_no2b_parse`, `PARSE`.
- Callers adapted mechanically (hf driver, hf_efa, bcs, hfb, rrpa driver/rrpa, ncsm): no Orb argument to the readers;
  the legacy solvers move `V_NN_res` back into `V_NN` (`O2b_move`), RRPA passes its bare V_NN directly.

Verification:
- Build (ifx) clean.
- Parse (bundled, Nmax 2/4/6): `NN_parsed.bin`, `NNN_parsed.bin` byte-identical to the Julia reference
  (`tests/dev_port/work/julia_t8/bundled/io/parse`); `PARSE` to the full space (3/6/9) reproduces `NN.bin`/`NNN.bin`.
- Hamiltonian unit test vs patched Julia (bundled, 16O, Calc (3,6,9) & (2,4,5), non-diagonal Rho, Sigma = V_NN):
  H1b <= 1.4e-14, H1b_H2b <= 2.3e-13 (scale 370), H2b (all blocks) <= 3.6e-15, V2b_correlation_energy 1.3e-15 rel.;
  H2b vs V2b_residual_no2b <= 3e-15. Also run with `-check all` (no bounds errors).
- Regression, old solvers, 40Ca 9 shells (`runs/test_Ca40_9sh/common_b`, 16 threads, 2:32 wall vs 3:16 common A):
  vs `fortran_new/IO` 29 of 77 files differ: the 3 radius summaries of common A, and 26 E1 current files, all intended
  (upstream E1 fix + 0.5/1/16pi normalization): `{TDA,RPA,RRPA}_E1_{V,T}{,_conv,_spin}`, `{TDA,RPA,RRPA}_E1_C`,
  `{RPA,RRPA}_E1_C_conv`, `{TDA,RPA,RRPA}_E1_NLO_LWA` (the `_spin` strengths are exactly 0.5x). Against `common_a/io`
  exactly these 26 files differ; all others (HF/MBPT energies, spectra, E0/E1/E2/E3/E1_c/M strengths, densities,
  Bin/*) are identical (max err 0).
- Performance (16 threads, HF only, warm page cache, base = 33d40d7 rebuilt, 3 alternating runs each):
  V2b_residual_no2b 2.25 s -> 1.84 s; NNN read 4.9-5.0 s -> 0.61-0.64 s; NN read 1.6 s -> 0.36 s; HF iterations
  2.8 s -> 1.8 s. Cold-cache first run of the new binary: NNN read 8.1 s (base cold: 12.9 s).
- QP chain smoke test (bundled O18 BCS -> QTDA -> QRPA with the legacy solvers, incl. F/GT qp transforms): runs;
  not comparable with Julia yet (dev BCS changes belong to the qp group).

Open items:
- The legacy RPA/RRPA/QRPA strength formulas of the E1 current operators still use the base combination; the dev
  formulas (VC + 0.5 g_s VS etc., RPA_rB) come with the tda_rpa/qp groups.
- `H1b`, `H1b_H2b`, `V2b_correlation_energy`, F/GT have no solver caller yet (RPA/HF-MBPT/pnQTDA groups).
- `V2b_parse`/`V3b_no2b_parse` count the whole interaction space (as dev); for the full DN2LOGO394 NNN file
  (emax 16) the count is O(minutes) with 32 threads and the parse reads ~188 GB over NFS - untested at that size.
- Not ported (dead in dev, group common B): 9-arg `H1b(...,Rho,Sigma)` (`H1b_dense`), `E1b`, `E2b`.

### Known interim limitation (review of stage 1, accepted until the cleanup stage)
- `orbitals_make` now uses dev's fractional filling of the last open shell. The LEGACY solvers still use whole-shell
  logic: old HF/BCS take `pO/nO == 1` as occupied, and old `orbitals_ph_make` takes `== 0` as a particle. For
  **open-shell** nuclei (e.g. 18O) the legacy HF density misses the fractional shell (wrong particle number), and
  legacy TDA/RPA/RRPA drop that orbital from the particle space. Closed-shell runs are unchanged (40Ca regression).
  Use closed-shell nuclei with the legacy solvers; the dev solvers (stages 2-4) replace them.

## 2026-11-20 — stage 2, HF (spec `hf.md` sections 1-7, 9-12; section 8 HF-MBPT is the next step)

Done:
- `src/hf/nj_hf.f90` rewritten to dev: `HF_solve` (`Orb` intent(inout), dev messages), `HF_orbitals_make` (= public
  `HF_orbital_ordering` + dev occupation rule with `ESH`, in-place update of `Orb`), `HF_density_operator(C,Orb)` with
  real weights, `HF_allocate`/`HF_energy` with the l_max/j_max pruning and the dev J range, `HF_summary` (new name,
  trailing newline), `HF_SPS_summary` (real O_a), `HF_export` (`HF/bin/{orb,C,rho,H1b}.bin`), new
  `HF_radial_potential` (`HF.V_HF`, `HF/densities/HF_radial_potential.dat`).
- `src/hf/nj_hf_driver.f90`: dev `HF`/`HF_solver` (dirs `io/<Path>/HF/{bin,densities}`, Calc-sized `orbitals_make`,
  `rho_export`, V_HF, BMF: non-mutating residual, `HF/bin/H2b.bin`, MBPT gated by `HF_MBPT_occupation_check`,
  `HF_MBPT(…,h_N,V_NN_res)`, `HF_MBPT_export(…,V_NN,V_NN_res,…)` with the bare V_NN kept only if `MBPT.Export`).
- `src/hf/nj_hf_mbpt.f90`: new dev interface; `HF_MBPT_occupation_check` ported; `HF_MBPT` and `HF_MBPT_export` are
  STUBS (print a NOTICE) for the HF-MBPT engineer (spec section 8). The base MBPT(2) moved unchanged to
  `nj_hf_mbpt_legacy.f90`.
- LEGACY chain kept (marked "LEGACY (removed in cleanup stage)"): `nj_hf_legacy.f90` (copy of the base nj_hf),
  `nj_hf_mbpt_legacy.f90`, `nj_hf_efa.f90`, `nj_hf_legacy_driver.f90` (action `HF_LEGACY`); BCS/HFB/HF_RPA export use
  `nj_hf_legacy`. In `nujuliet.f90`, `HF` runs `HF_LEGACY` (NOTICE) when the same `&run` has `HF_RPA*`, `HF_RRPA*`,
  `EMPM*` or `NCSM`. `nj_obdm::radial_grid` made public (HF_radial_potential grid). Makefile SRC_HF updated.
- Spec section 10 items: 1-5, 7-9, 13 were already done in stage 1 (types, io/, params/namelists, orbitals_export/import,
  rho_export, readers, non-mutating residual); 13 done here (HF_solve/HF_density_operator signatures; the legacy
  callers keep the legacy versions). 6, 10, 11, 12, 14 belong to HF-MBPT / already exist (OBDM_natural_orbitals,
  O1b helpers, nj_hamiltonian); 15 (EMPM/NCSM/RRPA reading Bin/*) stays with the legacy HF until stage 5/cleanup.
- Harness: namelists split into one file per Julia step group (`fortran/<case>/NN_<step>.nml`, HF steps
  `*_HF.nml`), `run_fortran.sh <case> [NT] [GLOB]`; new case `hf_options` (Julia driver `julia/hf_options.jl` +
  `fortran/hf_options/`, reference in `work/julia_t8/hf_options`): 16O `V_HF`, 18O `ESH` + `V_HF` + `CMS='None'`.

Verification (tol 1e-8, compare.py on `io/<Path>/HF`):
- bundled O16, O16_MBPT, He4; o16_open O20 (fractional 4/6 neutron 0d5/2); ca40_5sh Ca40 (16 threads):
  `orb.bin` exact, `C.bin` <= 8.3e-16 (column signs), `rho.bin` <= 2.2e-15, `H1b.bin` <= 1.5e-14,
  `HF_radial_densities.dat` exact, `HF_summary.dat` identical (O20 whole file; the others identical up to the
  missing HF-MBPT block), `H2b.bin` SIGN only: with the sign gauge of `C.bin` applied element-wise
  (scratch script) <= 2e-14 for all four nuclei. Missing: `HF_MBPT_occupation.dat`, `HF_MBPT_radial_densities.dat`
  (MBPT stub). Iteration counts identical.
- hf_options (new Julia reference): all 14 files ok, `HF_radial_potential.dat` max err 0 (8 of 263169 rows differ
  in the last printed digit at rounding ties). 20O with `ESH = true` stops with the same "closed-subshell" error in
  both codes (UPSTREAM_BUGS).
- Debug build (`-check all`): HF steps of all four cases run without runtime errors.
- Legacy regression `runs/test_Ca40_9sh/stage2_hf` (`HF_LEGACY`, `HF_RPA`, `HF_RRPA`, 16 threads, 2:19 wall) and
  `stage2_hf_auto` (original input with `HF`, auto-routed): both trees byte-identical to `common_b/io` (`diff -rq`).
- Legacy BCS -> QTDA -> QRPA and HFB -> QTDA (bundled 18O) still run (pnQTDA is stage 3).

Open items:
- HF-MBPT (spec section 8) replaces the stubs in `nj_hf_mbpt.f90`: `HF_MBPT`, `HF_MBPT_export` (signatures fixed by the
  driver), plus TBCF with F1. With `BMF = true` the dev HF currently writes `H2b.bin` and skips MBPT with a NOTICE.
- Cleanup stage: delete `nj_hf_legacy*`, `nj_hf_mbpt_legacy`, `nj_hf_efa`, the `HF_LEGACY` action/auto-routing.
- README.md / CLAUDE.md still describe the old HF layout (`IO/<Path>/Bin`), to be updated in the cleanup stage.

Not ported (dead in dev): the commented-out old loop in `HF_energy.jl`; `HF.HRF` is read but unused (as in dev).

## 2026-11-20 — stage 2, HF-MBPT (spec `hf.md` section 8; dev `src/HF/HF_MBPT.jl`)

Done:
- `src/hf/nj_hf_mbpt.f90` rewritten (stubs replaced): `HF_MBPT`, `HF_MBPT_energy` (MBPT(2) + MBPT(3) pp/ph/hh sectors,
  stdout lines and the dev summary block in `HF/HF_summary.dat`), `HF_MBPT_OBDM` (base dRho loops, dev tail:
  `HF/bin/{C,rho}.bin` read back, `OBDM_natural_orbitals`, `rho_export` to `HF/densities/HF_MBPT_radial_densities.dat`,
  `HF_MBPT_occupation` to `HF/densities/HF_MBPT_occupation.dat`, Export: `HF/bin/{orb,C,rho}_temp.bin`),
  `HF_MBPT_occupation_check`, `HF_MBPT_export` (all four `TBCF`/`ScH` branches, `io/<Path>/HF-MBPT/bin/{orb,C,H1b,rho,
  H2b[,sigma]}.bin`, temp files removed).
- New `src/hf/nj_hf_mbpt_common.f90` (`HF_MBPT_orbitals_make`, `HF_MBPT_h1b`, `HF_MBPT3_V2b_Pandya`, particle/hole lists)
  and `src/hf/nj_hf_mbpt_tbcf.f90` (`HF_MBPT_TBCF` with **F1**: nn channel classified with `nO`). Makefile SRC_HF.
- `nj_utils::remove_file` (spec 10.14). The HF driver needed no change (signatures fixed by the HF step).
- OpenMP: MBPT energy sectors parallel over the flattened (p,q) pair list with `reduction(+)`; OBDM loops reduce into
  private dRho copies; TBCF parallel over Bra per (J,P,T), each element written by one task (no race).
- Harness: new case `mbpt_options` (`julia/mbpt_options.jl`, `fortran/mbpt_options/0{1,2}_*_HF.nml`, reference
  `work/julia_t8/mbpt_options`): 28O (N != Z, fix F1) `Export, TBCF, ScH = false`; 16O `Export, ScH, TBCF = false`.
  New tool `tests/dev_port/gauge_check.py` (element-wise comparison of rotated-basis O1b/O2b files with the column-sign
  gauge of `C.bin`). The existing `*_HF.nml` of bundled/ca40_5sh/o16_open already carried the `&hf_mbpt` settings.

Verification (compare.py tol 1e-8 on `io/<Path>/HF` and `io/<Path>/HF-MBPT`; SIGN files re-checked with gauge_check.py):
- bundled O16, He4 (MBPT, no export): all 9 HF files ok (`HF_summary.dat` incl. the MBPT block & MBPT radii exact,
  `HF_MBPT_occupation.dat`, `HF_MBPT_radial_densities.dat` exact); `H2b.bin` SIGN (HF-basis, as before).
- bundled O16_MBPT (Export, TBCF, ScH): HF 9/9; HF-MBPT `orb.bin` 4.4e-16, `C.bin` 3.3e-16 (3 columns flipped),
  `rho.bin` 8.9e-16, gauge: `H1b.bin` 1.1e-14, `H2b.bin` 2.8e-15, `sigma.bin` 6.7e-16.
- ca40_5sh Ca40 (Export, ScH = false): HF 9/9 (MBPT E^(2) = -55.4908088, E^(tot) = -60.7876419 as Julia); HF-MBPT
  `orb.bin` 2.4e-15, `C.bin` 6.2e-15, `rho.bin` 5.2e-15, gauge: `H1b.bin` 2.3e-13, `H2b.bin` 3.4e-14.
- mbpt_options O28_TBCF: HF 9/9; HF-MBPT gauge `H1b.bin` 1.3e-13, `H2b.bin` 4.3e-14, `sigma.bin` 8.8e-16 (F1 active:
  protons 0s0p, neutrons up to 0d3/2). O16_ScH: HF 9/9; gauge `H1b.bin` 1.1e-14, `H2b.bin` 2.8e-15.
- o16_open O20 (fractional 0d5/2): MBPT skipped by the occupation check in both codes; HF 7/7 as before.
- Races: bundled HF steps with 1 vs 8 threads: every file identical to 1e-13 with `--strict-sign`.
- Debug build (`-check all -warn all`): bundled O16_MBPT and mbpt_options run without runtime errors or new warnings.
- Scaling: 40Ca 9 shells (DN2LOGO394, `runs/test_Ca40_9sh/stage2_mbpt`, Export+TBCF+ScH, 16 threads): whole HF run 7.0 s,
  MBPT energy 0.49 s (ph sector 0.30 s), OBDM 0.70 s, export 2.3 s; E^(2) = -146.50856125 = legacy MBPT(2) value.
- Timing ca40_5sh (8 threads, HF+MBPT export step): Fortran 0.93 s wall (MBPT energy 0.017 s, OBDM + radii 0.31 s, export
  0.019 s) vs Julia 58.5 s (incl. JIT and the NN/NNN reads; MBPT energy 2.2 s, OBDM 5.3 s, export 2.2 s).
- Legacy chain: `runs/test_Ca40_9sh/stage2_mbpt_legacy` (original input `HF`+`HF_RPA`+`HF_RRPA`, auto-routed to
  `HF_LEGACY`, 16 threads, 2:20 wall): `diff -rq` against `stage2_hf/io` byte-identical.

Open items:
- Later stages: TDA/RPA with `Reference = 'HF-MBPT'` (+ `Reference_TBCF`) consume `HF-MBPT/bin/*`.
- The temp-file hand-over (`HF/bin/*_temp.bin`) is kept as in dev; an in-memory hand-over is possible (no output change).

Not ported (dead in dev): `Orb1B_HF_MBPT3` (unused struct), the unreachable `return 0.0` of `HF_MBPT_h1b`.
