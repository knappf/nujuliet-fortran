
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

## 2026-11-20 — stage 3, track QP (spec `qp.md`: BCS, HFB, qpO2b, QTDA, QRPA, pnQTDA)

Done (details & deviations in `docs/PORTING_NOTES.md`, section "stage 3, quasiparticle chain"):
- Legacy: no legacy consumer of the BCS/HFB/QTDA/QRPA outputs (`consumers.md`), so the chain was replaced in place
  (actions `BCS`, `HFB`, `QTDA`, `QRPA`, new `pnQTDA`); it no longer uses `nj_hf_legacy`.
- `nj_qpo2b`: O^40 pp/nn direct-term sign (dev), H^22 pn(11,11) (engine `Fch = 4`, `Out = 7`, 4-block Pandya `pnpn`
  of common A), `qpO2b_1111_pn`, import file-size check, guard `N2max >= 2*Nmax`.
- `nj_bcs`: dev layout `BCS/{bin,densities}`, dev `HF_solve`/`HF_energy`/`HF_orbitals_make`, real occupation tests,
  **F3** (mid-gap Lambda).
- `nj_hfb`: dev layout `HFB/{bin,densities}` (+ `bin/lambda_2.bin`), modified Broyden mixing (`broyden_t`, `dsysv`),
  `HFB_chemical_potential_secant/_step/_gap/_gap_species` with the `Lambda_H` scheme, `HFB_orbital_matching`
  (Hungarian of core `nj_linalg` + Loewdin) in all orderings and both canonical bases, `HF_orbitals_make` in
  `HFB_diagonalize_BCS`, deterministic Lipkin-Nogami traces.
- `nj_qtda`/`nj_qrpa`/`nj_qex`: `Reference` folders, dev layout and file names (`transitions/`, `spectra/`,
  `amplitudes/`, `*_summary.dat`, `QTDA_plot.dat`/`QRPA_Plot.dat`); F/GT qp operators and the dev E1 current formulas come
  from stage-1 `nj_tr1b` (`qpTr1b_initialize`), the strength combination of `qex_rB` is unchanged in dev (verified).
- New `src/qp/nj_pnqtda.f90` (driver, `orbitals_pn2qp_make`, allocate, diagonalize, export) + `files.mk` + dispatcher.
- Harness: new case `qp_options` (`julia/qp_options.jl`, `fortran/qp_options/0[1-5]_*.nml`, reference
  `work/julia_t8/qp_options`); new tool `tests/dev_port/qp_gauge_check.py` (quasiparticle sign gauge for qpH2b and the
  QTDA/QRPA/pnQTDA amplitudes). The existing `bundled` namelists 11/12 were already correct.

Verification (binary `nujuliet_qp`, build dir `build_qp`, full tree; compare.py tol 1e-8; SIGN files with
qp_gauge_check.py):
- bundled O18_BCS (BCS -> QTDA -> QRPA, Reference BCS): 51 files, 47 ok, 4 SIGN (qpH2b, QTDA_X, QRPA_X/Y); gauge:
  qpH2b 3.0e-15, QTDA_X 1.9e-13, QRPA_X+Y 9.4e-12. All other files (BCS summary, C/U/V, qpH1b, densities, spectra,
  transitions incl. E1 currents) <= 3.7e-12 (QRPA_E.bin), text files to the printed digit.
- bundled O18_HFB (HFB -> QTDA -> pnQTDA -> QRPA): 58 files, 53 ok, 5 SIGN; gauge: qpH2b 3.6e-11, QTDA_X 3.5e-11,
  QRPA_X+Y 5.2e-11, pnQTDA_X 3.1e-11; same HFB iteration count (13) and history as Julia.
- qp_options (167 files): Ne22_HFB (Full pairing, open p & n; QTDA Ortho=false; pnQTDA), O20_LN (LNT/LNRes, Broyden
  restart, QRPA), O20_MCA, O20_ScBCS (ScBCS, A0/Z0, QTDA, pnQTDA, QRPA Ortho=false): all ok or SIGN, gauge <= 4.0e-10.
  O20_HBCS (HFB Pairing=BCS -> pnQTDA): same 14-iteration history, energies agree to 4e-10 rel., qpH2b/pnQTDA_X gauge
  1.3e-8 / 5.0e-8 and the printed dZ dispersion 4.5e-5 vs 5.5e-5: ill-conditioned converged solution (debug vs optimized
  Fortran builds differ by the same amount; see UPSTREAM_BUGS). 22Ne QRPA stops in both codes (complex eigenvectors).
- Threads: 1 vs 8 threads, bundled O18 steps and qp_options: every file identical to 1e-12 with `--strict-sign`
  (after making the LN traces deterministic).
- Debug build (`-check all -warn all`, `build_qp_dbg`): bundled O18 steps and qp_options run without runtime errors
  (only array-temporary warnings of the qpO2b engine call); no compiler warnings.
- Timing (8 threads, bundled interaction): bundled O18 BCS/HFB chains Fortran 1.6 s wall (whole run incl. reads) vs Julia
  30.1 s (sum of steps, excl. JIT); qp_options Fortran 2.7 s vs Julia 50.4 s.

Open items:
- `README.md` (top level) still documents the base chain-B layout; to be updated in the cleanup stage.
- The BCS reference-HF occupation question (UPSTREAM_BUGS, stage 3) is open for the owner.

Not ported (dead in dev): `HFB_Broyden_BackUp.jl` (not included), `HFB_canonical_basis_diagonalize`, `qpO2b_31`
(caller commented out), `QTDA/QRPA/pnQTDA_import_binary` (`*_external.jl`, no caller).

## 2026-11-20 — stage 3, TDA & RPA part 1 (spec `tda_rpa.md` sections 2-5 + section 9; track RPA)

Done:
- New `src/rpa/nj_ph_common.f90`, `nj_tda.f90`, `nj_rpa.f90`, `nj_rpa_export.f90`, `nj_tda_rpa_driver.f90`
  (`src/rpa/files.mk`), `&run` actions `TDA` and `RPA` in `nj_dispatch_rpa.f90`. Reference import from
  `io/<Path>/<Reference>/bin/{orb,C,rho,H1b,H2b[,sigma]}.bin` (HF or HF-MBPT, `Reference_TBCF`), 1p1h space from
  `orbitals_1p1h_make` (fractional occupations, Tol_par/Tol_hol), A/B/N with the metric and `h_N`, TBCF box diagrams
  (A & B; DEV pf/nf contact terms), TDA (dsyevr / dsygvd in the spurious complement) and RPA (H_p H_m product form,
  metric normalization, zero-norm states, 2x2 spurious problem) diagonalization, `RPA_OBDM_oY2`, QBA / density /
  dressing energies, `RPA_sigma` with **F2**, TDA exact kernel (incl. reference Sigma), RPA `Linear` / `Exact` rM (incl.
  RPA Sigma), shared rB (f3f7a33 E1 fixes), dev output layout `io/<Path>/<Name>/{bin,densities,transitions,spectra,
  amplitudes}` incl. RPA radii, natural orbitals (`RPA_Orb`), `sigma.bin`, summaries, spectra, amplitudes, 15 transition
  tables per solver. Shared: `nj_linalg::eig_sym_gen`; `compare.py` ignores the port-only `.TDA_output` marker.
- Legacy kept unchanged: `HF_RPA`, `HF_RPA_transition_densities/currents`, `HF_RRPA*` (+ `HF_LEGACY` routing).
- New test case `rpa_options` (`julia/rpa_options.jl`, reference `work/julia_t8/rpa_options`, `fortran/rpa_options/`):
  Ortho = false, RPA TBCF + RPA_Orb + Exact, RPA on HF-MBPT with Reference_TBCF + TBCF (dressing energy), and
  `RPA-DEV-TBCF` (step 06, DEV with Reference_TBCF, for part 2).

Left to PART 2 (not done here): `TDA_/RPA_transition_densities`, `_transition_operator_densities`,
`_transition_analysis`, `_transition_currents`, `_transitions_recalculate` (the external actions of the `02_*` / `06_*`
namelists of bundled, ca40_5sh, o16_open are still unknown actions), `RPA_allocate_SE` (RPA `DEV = .true.` with
`Reference_TBCF` stops with "not ported yet"; step `rpa_options/06`), the `levels/` dir. `Transitions = 'Exact'` is not
separate and is done. Part 2 can reuse `ph_operators`, `ph_metric`, `ph_occ`, `TDA_transition_kernel`, the X/Y layout
(channel arrays `(P,0:J_max)`, row-major binaries) and `ph_allocate_tbcf(...,DEV=.true.,...)` (pf/nf done).

Verification (compare.py tol 1e-8, TDA/RPA steps only, external actions stripped from `&run`; X/Y binaries checked
element-wise with the HF/natural-orbital row-sign gauge of `C.bin` plus the column gauge, scratch script):
- bundled O16 `TDA`, `RPA`, `RPA-Exact`, `RPA-DEV` (HF reference) and O16_MBPT `TDA-MBPT`, `RPA-MBPT` (HF-MBPT with
  Reference_TBCF): all summaries, spectra, amplitudes, radii, occupations, `orb/C/rho.bin`, all 15 transition tables per
  solver ok (0 DIFF); `TDA_X` <= 6.4e-14, `RPA_X/Y` <= 1.7e-13 with the gauge.
- rpa_options (`TDA-noOrtho`, `RPA-noOrtho`, `RPA-TBCF`, `RPA-dress`): 0 DIFF; `RPA_X/Y` <= 2.2e-12; `sigma.bin` gauge
  2.7e-14 (RPA-TBCF, HF basis, F2 active) and 2.1e-15 (RPA-dress); `RPA_Orb` `H1b.bin` 1.3e-14, `rho.bin` 5.5e-14.
- o16_open O20 `TDA`, `RPA` (0d5/2 dropped) and `TDA-frac`, `RPA-frac` (fractional 1p1h): 0 DIFF, X/Y <= 1.9e-13.
- ca40_5sh Ca40 `TDA`, `RPA` (16 threads): 0 DIFF, X/Y <= 5.8e-13. Fortran TDA 0.02 s, RPA 0.21 s.
- Threads: bundled 1 vs 8 threads, all 180 files equal at 1e-12 with `--strict-sign`.
- Debug build (`-check all -warn all`): bundled, rpa_options, o16_open TDA/RPA steps without runtime errors.
- Legacy: `runs/test_Ca40_9sh/stage3_rpa_legacy` (Ca40.nml, HF → HF_LEGACY, HF_RPA, HF_RRPA, 16 threads, 2:32 wall)
  byte-identical to `stage2_mbpt_legacy/io` (`diff -rq`).

Open items:
- Part 2 items above. Cleanup stage: the LEGACY `nj_hf_rpa*` / `nj_hf_rrpa*` files once EMPM (stage 5) and HF_ScRPA
  (stage 4) no longer need them.

Not ported (dead in dev): Sigma variant of `RPA_spurious_ini`, `RPA_OBDM_sigma`, commented-out blocks of
`TDA_allocate_TBCF` (pf/nf) and `RPA_export` (Sigma_NAT).

## 2026-11-20 — stage 3, TDA & RPA part 2 (spec `tda_rpa.md` sections 4.4, 6, 7; track RPA)

Done (details & deviations in `docs/PORTING_NOTES.md`, section "stage 3, TDA & RPA part 2"):
- New `src/rpa/nj_ph_external.f90`: `&run` actions `TDA_/RPA_transition_densities`, `_transition_operator_densities`,
  `_transition_analysis` (`levels/`), `_transition_currents`, `_transitions_recalculate` (dev names, files and folders),
  `RPA_transition_kernel` (`Linear`, `Exact` incl. the RPA Sigma), TDA with the exact kernel incl. the reference Sigma.
  New namelists `&transition_operator_densities`, `&transition_analysis` (`Cut`) in `nj_input` (additive).
- New `src/rpa/nj_rpa_se.f90`: `RPA_allocate_SE` (DEV, Mode OnShell, DSRG), called from the RPA allocation and from the
  dressing energy as in dev; the part-1 "not ported yet" stop is removed.
- Section 7: no races in the new code (disjoint writes per thread), guards for JP / nu range and missing inputs,
  `mkdir -p` of `levels/`, `densities/`, `transitions/E0..E3` (the part-1 guards: `.TDA_output`, Name == Reference,
  spurious 1e-8 else-branch, Reference_TBCF without sigma.bin).
- Test case `rpa_options` extended (Julia driver + reference regenerated, `fortran/rpa_options/07..12`): externals on
  TDA/RPA-noOrtho (analysis of the non-EJ 2-, 2+ densities with File_Name, 0+/3- operator densities, 2+ currents with
  M = 1 on 20 x 30, both recalculations), `RPA-Exact` + externals, `RPA-DEV-dress` (Reference_TBCF, TBCF, DEV, Exact:
  SE in A and in the dressing energy) + all RPA externals with the RPA sigma.bin, `TDA-MBPT` + externals with the
  reference Sigma, `RPA-DEV-frac` (Tol_par 0.9, Tol_hol 0.02: overlapping particle/hole sets, backward B_12 coupling).
  The old `O20_MBPT` idea was dropped: dev skips HF-MBPT for the non-idempotent open-shell HF.

Verification (binary `build_rpa/nujuliet_rpa`; compare.py tol 1e-8; work dirs `tests/dev_port/work/fortran_rpa*_t*`):
- bundled `io/O16` (steps 01-04: TDA, RPA, RPA-Exact, RPA-DEV + 1- densities, currents, operator densities, analysis
  of TDA and RPA): 153 files, 141 ok, 12 SIGN (eigenvector / HF-orbital signs of the binaries and of the analysis
  files), 0 DIFF; `io/O16_MBPT` (TDA-MBPT/RPA-MBPT + densities & currents): 89 files, 82 ok, 7 SIGN.
- rpa_options (new reference): 337 files, 308 ok, 29 SIGN, 0 DIFF/MISSING/EXTRA; `RPA-DEV-TBCF`, `RPA-DEV-dress`,
  `RPA-DEV-frac` (SE path) all ok. The analysis SIGN files (13 in rpa_options, 6 in bundled) were checked with a
  scratch script: identical a_p/a_h/T rows, |values| equal to 1e-8, the M columns flip with one sign per level.
- o16_open (TDA/RPA, TDA-frac/RPA-frac + 2+ densities): 113 files, 106 ok, 7 SIGN. ca40_5sh steps 01-02 (incl. the 1-
  densities): 66 ok, 6 SIGN (only HF-ScRPA, stage 4, missing).
- Threads: bundled, rpa_options, o16_open, ca40_5sh: 1 vs 8 (ca40 16) threads all files equal at 1e-12 with
  `--strict-sign`. Debug build (`build_rpa_dbg`, `-check all -warn all`): bundled 01-06 and rpa_options without
  runtime errors, results equal to the optimized build (tol 1e-10, eigenvector signs).
- Timing ca40_5sh (8 threads, wall per action incl. start-up and imports): TDA 0.41 s, RPA 0.63 s, TDA+RPA 1- densities
  0.97 s; Julia (`timings_julia.txt`, 8 threads) TDA 3.67 s, RPA 4.45 s, densities 0.94 s.
- Legacy: `runs/test_Ca40_9sh/stage3_rpa2_legacy` (HF → HF_LEGACY, HF_RPA, HF_RRPA, 16 threads, 3:07 wall)
  byte-identical to `stage2_mbpt_legacy/io` (`diff -rq`).

Open items:
- Cleanup stage: LEGACY `nj_hf_rpa*` / `nj_hf_rrpa*` (EMPM stage 5, HF_ScRPA stage 4). The full-namelist runs of
  bundled/ca40_5sh still stop at the stage-4 HF_ScRPA steps (not this track).
- `RPA_allocate_SE` memory: the configuration list of a channel is kept in memory (8 ints + 2 reals per 2p-2h
  configuration); fine for the test cases, large for production bases (DEV path only).

Not ported (unreachable in dev): `RPA_allocate_SE` modes `"QP"` (`SE_QP`) and `"Fixed"`, regularizations
`"Imaginary"`/`"None"`, `RPA_allocate_SE_Old`.

## 2026-11-20 — stage 3, independent verification (merged tree, no code changes)

- Clean build of the whole tree (`build/`, `./nujuliet`, plain `make`): no errors, no warnings (`build_rpa`, `build_qp`,
  `nujuliet_rpa`, `nujuliet_qp` removed).
- All cases, 8 threads, vs `work/julia_t8` (compare.py 1e-8): hf_options 14 ok; mbpt_options 23 ok + 6 SIGN; o16_open
  106 ok + 7 SIGN; rpa_options 308 ok + 29 SIGN; ca40_5sh 66 ok + 6 SIGN + 28 MISSING (HF-ScRPA, stage 4); bundled 333 ok
  + 29 SIGN + 85 MISSING (He4 HF-ScRPA*, stage 4; steps 08-10 / ca40 03 not run); qp_options 150 ok + 14 SIGN + 3 DIFF
  (O20_HBCS: HFB summary dZ 4.5e-5 vs 5.5e-5, qpH2b gauge 1.3e-8, pnQTDA_X gauge 5.0e-8; Julia -t 1 vs -t 8 differs by
  the same amount there (6.6e-9 / 2.0e-8, dZ 5e-5): ill-conditioned converged HFB, not a port defect). No EXTRA.
- SIGN files re-checked with gauges: H1b/H2b/sigma (gauge_check.py) <= 2.3e-13; TDA_X/RPA_X+Y with the HF / natural
  orbital row gauge + column gauge <= 5.7e-12 (RPA-DEV-frac), all others <= 2.2e-12; transition-analysis `levels/`
  files exact to the printed digits; qp files (qp_gauge_check.py) <= 4.0e-10 (except O20_HBCS above); HFB Rho/Kappa
  (canonical-basis gauge) <= 2.1e-12, U/V <= 1.3e-11.
- bundled Fortran -t 1 vs Julia -t 1: 350 ok + 12 SIGN, 0 DIFF. Fortran 1 vs 8 threads (--strict-sign, 1e-12): all cases
  identical except rpa_options RPA-DEV-frac RPA_E/rho 1.5e-12: the stage-1/2 `reduction(+)` loops with dynamic schedule
  (e.g. `nj_hamiltonian` H1b, used for HF-MBPT `H1b.bin`) vary at 1 ulp from run to run (also with 8 threads twice);
  RPA-DEV-frac amplifies this to ~1e-12. With a fixed input the RPA step is bitwise reproducible with 8 threads.
- Legacy: Ca40.nml (HF, HF_RPA, HF_RRPA, 16 threads, 2:54) byte-identical to `stage2_mbpt_legacy/io`.
- 40Ca 9 shells dev actions (`runs/test_Ca40_9sh/stage3_dev`, HF + TDA + RPA, 16 threads, 8.6 s wall): TDA 0.18 s, RPA
  0.50 s vs legacy HF_RPA 0.71 s (TDA+RPA). Spectra equal to the legacy ones except the spurious 1- energy (dev keeps
  0.586 / 0.433 MeV, legacy 0) and the sign of Im E of the unstable 3- pair (legacy prints the channel as "3+").
- Timing (8 threads, whole runs incl. reads): Fortran bundled 9.0 s (Julia steps 90 s), ca40_5sh 2.5 s (101 s),
  o16_open 1.8 s (21 s), rpa_options 8.0 s (55 s), hf_options 1.4 s (12.6 s), mbpt_options 1.6 s (22 s), qp_options
  2.9 s (50 s). No step slower than Julia.
- Edge cases: Name == Reference fatal (TDA & RPA); TDA re-run with the same Name in the same io dir works (old folder
  replaced, identical results); a TDA into a non-TDA folder is refused; RPA Reference_TBCF without sigma.bin fatal.
  TDA Reference_TBCF without sigma.bin prints dev's error and skips the step (exit status 0, as dev); the RPA fatal
  comes after the output folder tree has been created (empty `<Name>/` tree remains).

## 2026-11-20 — stage 4 part 1, HF-ScRPA building blocks (spec `scrpa.md` sections 1.3-1.7, 2; no solver yet)

Done (details & deviations in `docs/PORTING_NOTES.md`, section "stage 4 part 1"):
- Core (spec section 2): nothing new needed: `OBDM_natural_orbitals` (exact recipe + F5), `orbitals_import/export(path)`,
  `orbitals_update`, O1b wrappers, `H1b`, `H2b`, `H1b_H2b` (V3N.Sigma), `V2b_correlation_energy` exist since stage 1;
  `RPA_energy`, `_QBA`, `_dress`, `_rho` (Density) and `RPA_sigma` (F2) since stage 3 (`nj_rpa`). Reused, not duplicated.
- New `src/rpa/nj_rpa_rho.f90`, `src/rpa/nj_rpa_hamiltonian.f90` (files.mk after `nj_rpa_export`), new
  `nj_rpa::RPA_allocate` (now also used by `RPA_solver`) and `nj_rpa::RPA_sigma_mix`, type `nj_types::orb1b_rpa_obdm_t`,
  thread-local MKL control `src/core/nj_blas_threads_mkl.f90` (`COMPILER=gnu`: `_ref.f90`, no-op; Makefile `BLASTHR_SRC`).
- Unit check `tests/dev_port/unit_s4/` (Julia dump script + Fortran test program + runner, README section).

API for part 2 (all `public`; X/Y/E in the `RPA_diagonalize` layout `(2,0:J_max)`):
- `type(orb1b_rpa_obdm_t)` (`nj_types`): `p, h` (`pnorbvector_t`, particle / hole orbitals per isospin), `N_p, N_h`.
- `nj_rpa_rho`:
  - `function RPA_OBDM_orbitals_make(Params,Orb) result(Orb_ph)` (strict `|O| < Tol_par` / `> Tol_hol`, RPA tolerances).
  - `function RPA_OBDM_oY4(Params,Orb,Orb_ph,Orb_1p1h,Rho,X,Y) result(Rho_new)` (= Rho + dRho; dev's pM/nM dropped).
  - `subroutine RPA_OBDM_transform_amplitudes(Orb_1p1h,X,Y,C,X_new,Y_new)` (`cmat_t, allocatable, intent(out)`).
  - `subroutine HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT,C_LHO_NAT)`: `Orb` inout (occupations
    replaced: dev's returned `Orb_NAT`), `Rho` inout (in `Rho_ini`, out new diagonal OBDM, new NAT basis), `C_HF_NAT`
    inout (dev in place), `C_LHO_NAT` out; X, Y intent(in). Reads `io/<Path>/<RPA.Reference>/bin/{rho,C}.bin`.
    Dev call `Rho, C_LHO_NAT, Orb = HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT)` becomes
    `call HF_ScRPA_OBDM_solve(Params,Orb,Orb_ph,Orb_1p1h,X,Y,Rho,C_HF_NAT,C_LHO_NAT)`.
- `nj_rpa_hamiltonian`:
  - `subroutine RPA_Hamiltonian(Params,Orb,Orb_NN,Orb_NNN,C,Rho,Sigma,T,V_NN_bare,V_NNN_bare,h_N,V_NN)` (C = C_LHO_NAT;
    Sigma, T, the bare interactions & Orb_NNN only read for ScH; Sigma only for TBCF; for `ScH = false` they may be
    unallocated, i.e. decision 6 can skip `V2b_read`/`V3b_no2b_read` there, Orb_NN then from `O2b_import` of H2b.bin).
  - `subroutine RPA_Hamiltonian_HF_reference(Params,Orb,Orb_NN,C,h_N,V_NN)`.
- `nj_rpa`:
  - `subroutine RPA_allocate(Params,Orb,Orb_1p1h,Orb_NN,h_N,V_NN,Rho,Sigma,TBCF,A,B,N)` (prints dev's two lines; TBCF:
    box diagrams + DEV `RPA_allocate_SE` exactly as dev).
  - `subroutine RPA_sigma_mix(Sigma,Sigma_old,gamma)` (in place).
  - existing: `RPA_diagonalize(Params,Orb,Orb_1p1h,A,B,N,TrOp,E,X,Y,Stability)`, `RPA_energy_QBA(Params,Orb,Orb_1p1h,E,X,Y)`,
    `RPA_energy(Params,Orb,Orb_1p1h,Orb_NN,E,X,Y,Rho,Rho_ref,h_N,Sigma,V_NN,Have_Sigma_AB,Sigma_AB,Density)` (pass
    `Density=.true.`; `Have_Sigma_AB = TBCF`, Sigma_AB = Sigma, an unallocated o2b_t otherwise),
    `RPA_sigma(Params,Orb,Orb_1p1h,Orb_NN,X,Y,beta,Report)`, `RPA_rM`; `nj_ph_common::ph_rB`;
    `nj_rpa_export::RPA_export(Params,Orb,Orb_1p1h,Orb_NN,Stability,E_corr,E,X,Y,rB,C,Rho_LHO,Sigma,h)`.
- Core: `O2b_transformation(Params,Orb,Orb_NN,O_NN,C,Report)` (in place), `O1b_transformation(O,C)` (= C^T O C),
  `O1b_transpose`, `O1b_unit_matrix(Params)`, `orbitals_update(Orb,Rho)`, `rho_chR2(Params,Orb,C,Rho)`,
  `Tr1b_initialize/Tr1b_transformation`, `V2b_read(Params,V_NN,Orb_NN)`, `V3b_no2b_read(Params,V_NNN,Orb_NNN)`,
  `T1b(Params,Orb)`. The `make_output_dirs`/`print_header` helpers of `nj_tda_rpa_driver` are private (part 2: make them
  public or copy the 10 lines).

Verification:
- Build: `make BUILD=build_s4 EXE=nujuliet_s4` clean (no warnings); debug build (`make debug`) clean, the new files also
  without `-warn all` warnings.
- Unit check (`unit_s4/run_unit.sh`, 4He bundled, common input = Julia initial X/Y), max abs deviations: oY4 3.5e-18,
  oY4 complex path vs Julia (generic column phases) 5.2e-18, complex vs real path (phases with D^2 = -1) 1.7e-18,
  OBDM_natural_orbitals C 1.7e-15, transform_amplitudes 1.2e-15 (complex 2.2e-16), HF_ScRPA_OBDM_solve (19 / 19 / 10
  converged inner iterations, printed deltas identical to Julia) Rho <= 1.0e-15, C_LHO_NAT / C_HF_NAT <= 2.3e-15 (default, OBDM_diagonal,
  OBDM_alpha = 1), RPA_sigma 4.9e-17, RPA_Hamiltonian h_N / V_NN <= 3.6e-14 (ScH, ScH+TBCF with V3N.Sigma, HF
  reference), after one ScRPA step E_QBA 8.1e-16 rel., E_sigma 1.9e-15, E_rho 1.1e-15, E_dens 2.7e-15 rel.
  Same with 1 thread; debug build (-O0, -check all) worst 1.1e-13.
- Regression: bundled O16 steps 01-04 (153 files: 141 ok, 12 SIGN, 0 DIFF) and rpa_options (337: 308 ok, 29 SIGN) as in
  stage 3; rpa_options output byte-identical to the stage-3 binary (`diff -r`), i.e. the `RPA_allocate` refactoring
  changes nothing. Legacy code untouched (HF_RRPA chain not re-run).
- Synthetic 208Pb 13-shell size (orbitals_make occupations, random real X/Y, sum N_ph 12770, largest block 689):
  `RPA_OBDM_oY4` 0.09-0.10 s, `RPA_OBDM_transform_amplitudes` 0.04-0.05 s per call (16/32 threads), one complex block
  0.16 s, so an inner OBDM iteration costs ~0.15 s (at most 150 per outer iteration).

Open items (part 2): `nj_scrpa` (HF_ScRPA, HF_ScRPA_solver, scrpa_beta, ScRPA_terminate_summary), `&run` action
`HF_ScRPA`, harness steps (bundled 08-10, ca40_5sh 03 + the options of the task list), decision 6 (skip the bare
reads for ScH = false). UPSTREAM_BUGS: dev's `M = -1/2 X^T W X^*` is phase dependent for complex amplitudes.

Not ported (dead in dev): `RPA_OBDM_sigma`, 9-argument `H1b(...,Sigma)` (V_NN-Sigma), the `Delta_E` bookkeeping of
the solver (part 2).

## 2026-11-20 — stage 4 part 2, HF-ScRPA solver & action (spec `scrpa.md` sections 1.1, 1.2, 2.3)

Done (deviations in `docs/PORTING_NOTES.md`, section "stage 4 part 2"):
- New `src/rpa/nj_scrpa.f90`: `HF_ScRPA` (RPA.jl:72-141: header, reference check, B2 guard, output tree replaced) and
  `HF_ScRPA_solver` (RPA_solver.jl:94-418): reference import, frozen `Orb_1p1h` / `Orb_ph`, initial RPA, outer loop
  (TBCF: `RPA_sigma` with ramp `scrpa_beta` + `RPA_sigma_mix`, `HF_ScRPA_OBDM_solve`, Sigma rotated with
  `C_step = C_HF_NAT_old^T C_HF_NAT`, `RPA_Hamiltonian`, `Tr1b_transformation`, `RPA_allocate` (DEV inside, as dev),
  `RPA_diagonalize`, E_QBA criterion in 2 consecutive iterations with beta = 1 for TBCF), Stability / E_QBA history,
  `ScRPA_terminate_summary` (non-convergence: only `RPA_summary.dat`), final step (`RPA_sigma` always, `rho_chR2`,
  `RPA_energy(Density=.true., Sigma_AB = mixed Sigma)`, `RPA_rM`, `ph_rB`, `RPA_export` with `C_LHO_NAT`, `Rho_LHO`,
  `h_N`; `RPA_Orb` & `sigma.bin` handled by `RPA_export`).
- `&run` action `HF_ScRPA` in `nj_dispatch_rpa`; `files.mk` (after `nj_tda_rpa_driver`); `print_header`,
  `print_threads`, `make_output_dirs` of `nj_tda_rpa_driver` public. `&rpa` already had all dev fields (stage 3).
- Decision 6: `ScH = false` skips `V2b_read` / `V3b_no2b_read` (Orb_NN from the `H2b.bin` import; identical table).
- Legacy `HF_RRPA` / `&rrpa` / `nj_hf_rrpa*` unchanged and reachable (`! LEGACY` markers).
- Harness: new case `scrpa_options` (`julia/scrpa_options.jl`, `fortran/scrpa_options/01..12`, Julia reference in
  `work/julia_t8/scrpa_options`): ScH = false (non-converging, summary only), ScH = false + TBCF with TBCF_alpha 0.5 /
  TBCF_beta 3, TBCF_alpha 1 + TBCF_beta 2 + OBDM_alpha 0.8 + RPA_Orb + Exact, TBCF + DEV, Imax = 3 non-convergence
  summary with TBCF, RPA externals on a ScRPA folder, HF-MBPT reference (16O), fractional open-shell reference (20O
  ESH = false, Tol_par 0.9 / Tol_hol 0.5). `compare.py` now compares numbers in parentheses (`(relative 0.029)`)
  numerically; counts of the stage-3 cases unchanged. README sections updated.

Verification (`BUILD=build_s4 EXE=nujuliet_s4`, clean, no warnings; `make debug` clean, `-warn all` no warning in
`nj_scrpa`; debug binary (`-check all`) runs all of scrpa_options without a runtime error, results equal to 1e-10):
- bundled 08-10 (He4 HF-ScRPA, -TBCF, -OBDMdiag) vs Julia: iteration histories (outer + inner, 390 lines) identical;
  28/29/28 files ok except `RPA_X/Y.bin` SIGN (eigenvector columns) and `sigma.bin` SIGN; gauge_check (C.bin column
  signs of the HF proton orbitals propagate to NAT): sigma <= 4.3e-15, rho <= 1.8e-14. Whole bundled case: 0 DIFF.
- ca40_5sh 03 (Ca40 HF-ScRPA) vs Julia: histories identical (227 lines), 26 ok + RPA_X/Y SIGN, rho gauge 3.8e-14.
- scrpa_options vs Julia (12 folders): 0 DIFF; SIGN only `RPA_X/Y.bin`, `sigma.bin`, `H1b.bin` (RPA_Orb), one 1-
  analysis table (eigenvector sign); gauge_check sigma <= 1.1e-14, rho <= 1.3e-14, H1b 6.7e-15. Histories identical
  except 4 last-digit lines (1e-9 MeV) late in the chaotic non-converged ScH = false run (100 iterations); its summary
  agrees to 1e-10.
- Threads: 1 vs 8 equal at 1e-12 (`--strict-sign`) for bundled He4 ScRPA and all of scrpa_options (1e-11 incl. the
  chaotic run); ca40_5sh HF-ScRPA: `RPA_X.bin` 3.2e-12 with 19/320 eigenvector columns of opposite sign (MKL `dgeevx`
  phase with different thread counts, stage-3 `RPA_diagonalize`; oY4/sigma are bilinear in the columns), all other
  files at 1e-12.
- Legacy: `runs/test_Ca40_9sh/stage4_legacy` (Ca40.nml: HF -> HF_LEGACY, HF_RPA, HF_RRPA, 16 threads, 2:59) byte-identical
  to `stage2_mbpt_legacy/io` (`diff -rq`).
- Timing (8 threads, HF_ScRPA step incl. reads): ca40_5sh Fortran 0.82 s vs Julia 33.8 s (1 thread 2.7 s); He4 bundled
  0.31-0.53 s vs 1.3-8.3 s; scrpa_options steps 0.09-0.86 s vs 1.5-16.8 s. 40Ca 9 shells (`runs/test_Ca40_9sh/
  stage4_scrpa`, default RPA parameters on the stage-3 dev HF, 16 threads): converged in 26 outer / 325 inner
  iterations (no cycled or diverging inner loop), 18.6 s total incl. the bare NN/NNN reads, peak RSS 0.9 GB; per
  outer iteration RPA_Hamiltonian 0.42 s, RPA_diagonalize 0.09 s, OBDM solve 0.03 s. The legacy HF_RRPA needed 166 s
  for the same nucleus. 208Pb 13 shells not run (the node is busy with the 64-thread production run): an outer
  iteration there costs one H1b+H2b build (as one dev-HF iteration) + at most 150 x ~0.15 s inner OBDM iterations
  (part 1 synthetic timing) + the RPA diagonalization.

Open items:
- 208Pb 13-shell HF_ScRPA production run (needs the dev HF reference `io/<Path>/HF`; the v2 run uses the legacy
  HF_LEGACY output layout).
- Cleanup stage: remove `HF_RRPA`, `&rrpa`, `nj_hf_rrpa*`, `rrpa_params_t`, `obdm_iteration_pairs_t`, the
  `HF_RRPA_transition_*` actions (spec 2.3, decision 2).
- Not ported (dead in dev): `RPA_OBDM_sigma`, 9-argument `H1b(...,Sigma)`, `Delta_E` bookkeeping of the solver.
- Kept on purpose: B3 (particle/hole test inconsistency), dev's `X^T W X^*` in oY4 (UPSTREAM_BUGS).

## 2026-11-20 — stage 4, independent verification (HF-ScRPA)

- Clean build (`build_s4`, `nujuliet_s4`, `build/`, `./nujuliet` removed; plain `make`, 105 s): no errors, no warnings.
- All cases, 8 threads, vs `work/julia_t8` (compare.py 1e-8), full io trees, **0 MISSING / 0 EXTRA everywhere**
  (directory listings also equal apart from the harness `parse/` dir and the stage-3 `.TDA_output` markers):
  bundled 411 ok + 36 SIGN; ca40_5sh 92 ok + 8 SIGN; o16_open 106 ok + 7 SIGN; hf_options 14 ok; mbpt_options 23 ok +
  6 SIGN; rpa_options 308 ok + 29 SIGN; qp_options 150 ok + 14 SIGN + 3 DIFF (the known ill-conditioned O20_HBCS HFB:
  summary dZ 4.5e-5 vs 5.5e-5, qpH2b 1.3e-8, pnQTDA_X 5.0e-8, as in stage 3); scrpa_options 197 ok + 22 SIGN.
- SIGN files of the ScRPA folders: `RPA_X/Y.bin` with a row (NAT orbital pair) + column (eigenvector) sign fit per
  block <= 2.5e-14 (He4), 1.6e-11 (Ca40 HF-ScRPA), 8.4e-13 (O16 ScRPA-MBPT), 9.5e-13 (O20 ScRPA-frac); `sigma.bin` /
  `rho.bin` / `H1b.bin` (RPA_Orb) with gauge_check.py (C.bin column signs, proton NAT columns 1,2,4,5,7 flipped in He4)
  <= 9.9e-15 / 1.8e-14 / 6.7e-15; one `levels/` analysis table (eigenvector sign). Non-ScRPA SIGN files as in stage 3.
- Iteration histories (outer `HF-ScRPA iteration` + inner `OBDM iteration` lines): ca40_5sh 258 and bundled 435 lines
  identical to Julia; scrpa_options 3776 lines, 4 differ in the last printed digit (iterations 62-91 of the chaotic
  non-converging `ScRPA-noScH` run).
- bundled Fortran -t 1 vs Julia -t 1: 428 ok + 19 SIGN, 0 DIFF. Fortran 1 vs 8 threads (`--strict-sign`, 1e-12): all
  cases identical except ca40_5sh HF-ScRPA `RPA_X.bin` 3.2e-12 (19/320 eigenvector columns with opposite sign, dgeevx)
  and the known rpa_options RPA-DEV-frac `RPA_E`/`rho` 1.5e-12. Two 8-thread runs are not bitwise equal: besides the
  known HF-MBPT `H1b.bin` (stage 2) the TBCF ScRPA folders (He4 HF-ScRPA-TBCF, ScRPA-DEV) differ at <= 7e-14, because
  the V3N.Sigma part of the core `H1b_H2b` (`nj_hamiltonian`, stage 1) sums into `ph, nh` with `reduction(+)` and
  `schedule(dynamic,4)` (not new stage-4 code; within the 1e-12 target; open item: fixed-order sum).
- Legacy: `Ca40.nml` (HF -> HF_LEGACY, HF_RPA, HF_RRPA, 16 threads, 2:43; again with the final binary after the import fix, 2:36) byte-identical to `stage2_mbpt_legacy/io`.
- Edge cases: `Name == Reference` fatal (`ERROR: &rpa Name must differ from Reference`, exit 1, reference untouched);
  `Imax = 2` -> only `RPA_summary.dat` (termination summary with history; empty `bin/` etc. subfolders as Julia);
  missing `Reference` folder -> dev's "No precomputed reference solutions ..." message, return, no output folder,
  exit 0 (dev behaviour); `Reference` = a folder without `orb.bin` -> `ERROR: orbitals_import: cannot open ...`;
  `Reference` = an RPA folder (no `H1b.bin`) aborted with an ifx runtime error + traceback -> **fixed**: `O1b_import` /
  `O2b_import` now stop with `ERROR: ... cannot open <path>` (PORTING_NOTES stage 4 part 2 item 7; outputs unchanged,
  all cases re-run at 8 threads with the same compare counts).
- Performance, 40Ca 9 shells (`runs/test_Ca40_9sh/stage4_scrpa/Ca40_scrpa.nml`: HF (BMF) -> HF_ScRPA, default `&rpa`,
  32 threads, node shared with the 64-thread Pb run): 18.4 s wall in total (HF + export 3.6 s, HF_ScRPA 14.6 s incl.
  1.5 s for the bare NN/NNN reads), peak RSS 0.92 GB; converged in 26 outer iterations, every inner OBDM loop converged
  (325 inner iterations, ~12 per outer; no cycling / divergence), E0 -144.8447 MeV. Per outer iteration 0.43-0.71 s:
  `RPA_Hamiltonian` 0.34 s (H1b + H2b with V3N on the NAT OBDM + O2b_transformation), `RPA_diagonalize` 0.12 s, OBDM
  solve 0.03-0.06 s, `RPA_allocate` 0.005 s. Legacy HF_RRPA: 154 s (52 iterations, 3 s each, stalled at eta 2e-6 with a
  cycled degenerate solution) -> ~10x faster in total. Thread scaling of the HF_ScRPA step: 24.7 / 19.0 / 14.9 s with
  8 / 16 / 32 threads (`RPA_Hamiltonian` 0.70 / 0.47 / 0.34 s; diagonalization flat at ~0.1 s).
  With `TBCF = .true.`: 22 outer iterations, 46 s, ~2 s per outer: `RPA_allocate` (box diagrams, `ph_allocate_tbcf`)
  0.7-1.2 s, `RPA_Hamiltonian` (H1b_H2b) 0.6 s, RPA_sigma + Sigma rotation ~0.2 s.
  Profile (perf, 32 threads, non-TBCF): 41 % OpenMP spin-wait (serial sections / imbalance at this small size),
  H2b V3N loop 14 % + V3b_no2b lookups 7 %, O2b_transformation loops 9 % + o2b_pp/pn/nn/orb_precomp index helpers 10 %,
  MKL (diagonalization) ~4 %. Dominant routine: `RPA_Hamiltonian` (H2b V3N contraction, then O2b_transformation).
- 208Pb 13 shells estimate (32 threads; scaling with the legacy Pb/Ca ratios of the same-kind steps in
  `runs/Pb208_13sh_hw12_v2` (read only): residual-interaction update 32.5 s @64 vs 2.1 s @16 (~60x CPU), RPA
  diagonalization 21.5 s vs 0.22 s (~100x)): `RPA_Hamiltonian` ~20-25 s, `RPA_diagonalize` ~10-15 s, OBDM solve
  ~2-3 s (12-20 inner iterations x 0.15 s, synthetic part-1 timing; at most 150 x 0.15 s = 22 s) -> ~35-50 s per outer
  iteration (legacy at 64 threads: ~1110 s, of which the legacy OBDM 1055 s), i.e. ~15-25 min for ~25 outer
  iterations. TBCF adds the box diagrams (~1 s here, ~60 s estimated) and RPA_sigma. Hot spots to watch: the V3N loop of
  `H2b` / `H1b_H2b`, `O2b_transformation` (index helpers dominate over the BLAS part at this size), the block loop of
  `RPA_diagonalize` not scaling with threads, `ph_allocate_tbcf` for TBCF.

## 2026-11-21 — stage 4, review follow-up (deterministic V3N.Sigma sum)

Review findings: no blocker / major / minor findings. Actionable open item fixed:
- `H1b_H2b` V3N.Sigma term (`src/core/nj_hamiltonian.f90`, stage-1 code): the OpenMP `reduction(+:ph,nh)` with dynamic
  scheduling is replaced by a fixed-order sum over 64 interleaved chunks of triples with per-chunk accumulators
  (PORTING_NOTES stage 4 part 2 item 8). Not re-read line by line against Hamiltonian.jl:1226-1336 in this pass; the
  formula is unchanged and the outputs agree with Julia as before (below).

Tests (plain `make`, no warnings; 8 threads unless noted):
- All cases vs `work/julia_t8` (compare.py 1e-8), same counts as the stage-4 verification: bundled 411 ok + 36 SIGN;
  ca40_5sh 92 + 8; o16_open 106 + 7; hf_options 14; mbpt_options 23 + 6; rpa_options 308 + 29; qp_options 150 + 14 +
  3 DIFF (known O20_HBCS HFB); scrpa_options 197 + 22. 0 MISSING / EXTRA.
- Fortran 1 vs 8 threads (`--strict-sign`, 1e-12): bundled 447/447 ok, scrpa_options 219/219 ok.
- Two 8-thread runs: bundled, scrpa_options, mbpt_options, rpa_options **bitwise identical** (`diff -rq`: 0 files). This
  also removes the stage-2 run-to-run variation of HF-MBPT `H1b.bin` (it came from this term) and the 7e-14 variation of
  the TBCF ScRPA folders.
- New vs previous binary (1e-13, strict): bundled and mbpt_options all ok; scrpa_options `ScRPA-MBPT` `RPA_X.bin`
  6.2e-13; rpa_options `RPA-DEV-frac` `RPA_E`/`rho` 1.4e-12 / 1.7e-12 (the known sensitive case).
- Legacy: `runs/test_Ca40_9sh/stage4_final_legacy` (Ca40.nml, 16 threads, 2:17) byte-identical to
  `stage2_mbpt_legacy/io` (`diff -rq`).

Still open (owner decisions / not run):
- Missing `Reference` folder: dev's message + return, exit status 0 (kept; make it fatal?).
- 208Pb 13-shell HF_ScRPA run (needs a dev HF reference; node busy). Timing estimate unchanged (35-50 s per outer
  iteration at 32 threads without TBCF).
- Edge cases checked only by reading: RPA 1- with `Ortho = .true.` and N_ph = 1; transition currents with M = 0 for
  E1/E3; TDA externals with `Reference_TBCF`. Stage-1 `OBDM_natural_orbitals` and the other routines listed by the
  review as not compared line by line are covered only by the matching test cases.

## 2026-11-21 — stage 5 part 1, EMPM on the dev HF reference (spec `consumers.md` sections 2, 3; owner decisions 1-8)

Done (details & deviations in `docs/PORTING_NOTES.md`, section "stage 5 part 1"):
- Shared helpers (section 2; the rest existed since stages 1-4): `nj_params::ref_path`; `nj_orbitals::read_h1b_diag`
  (diagonal check 1e-10) and the reference specific of the generic `orbitals_ph_make` (Orb from `orb.bin`, H1b path,
  holes/particles with 1e-8, fatal on fractional occupations); `nj_linalg::spurious_ortho` (core, moved from
  `nj_ph_common::ph_spurious_ortho`, used by the dev TDA/RPA and the EMPM). `&empm Reference` (default 'HF') in
  `nj_params` / `nj_input` (`read_empm` now assigns component-wise).
- EMPM (section 3.1-3.4): all seven steps read `io/<Path>/<Reference>/bin/{orb,C,H1b,H2b}.bin`; reference / BMF /
  diagonal / integer-occupation checks; dev spurious convention (E(1) = <spur|A|spur>, cut 1e-14, core `spurious_ortho`,
  `eig_sym`); output layout `io/<Path>/EMPM/{bin,spectra,transitions}` + `TDA_summary.dat` (dev text), other names
  unchanged; `States_Dir` default `EMPM/bin`; local `EMPM_phonon_count` and summary helpers, so `src/empm` no longer
  uses `nj_hf_rpa` / `nj_hf_rpa_export` (`files.mk` comment updated). `nujuliet.f90`: `EMPM*` no longer route `'HF'`
  to `HF_LEGACY` (`NCSM` still does; part 2).
- New test case `empm_ncsm` (no Julia counterpart): legacy baseline inputs `tests/dev_port/fortran/empm_ncsm_legacy/`
  (frozen binary `nujuliet-snap-d753f29`, `HF_LEGACY` + EMPM_* + NCSM), dev inputs `fortran/empm_ncsm/` (dev HF +
  EMPM_* + dev TDA; the NCSM twins are `*.nml.todo` for part 2), error cases `fortran/empm_ncsm/errors/`,
  `run_empm_ncsm.sh baseline|dev|compare`, `empm_compare.py` (legacy→dev name map, max deviation per file, eigenvector
  gauges, `--tda` mode). README section "Stage-5 EMPM / NCSM checks".

Verification (binary `nujuliet_s5`, build dir `build_s5`; baseline `work/empm_baseline`, dev `work/empm_dev_t16`):
- Baseline (frozen binary, 16 threads): 18 files, 45 s wall in all (02_He4c_EMPM 38.8 s, every other file <= 1 s).
  Dev chain: 12 files, 45 s (02_He4c_EMPM 38.1 s).
- dev HF vs legacy HF (Tol 1e-10): `C/H1b/H2b/rho.bin` identical for 4He, <= 3.6e-15 for 16O / 40Ca (C, H2b column
  signs for 40Ca); same level order (closed shells, Aufbau = naive filling here).
- Dev vs baseline, `empm_compare.py` (1714 files, tol 1e-10): energies / spectra (875 files: `*_E*.bin`, `Spectrum.dat`,
  `Full_*.dat`, `TDA.dat`) max 7.0e-11 relative (O16w 2PH 0+, ill-conditioned spurious D block; all others <= 9.4e-12),
  `*_Ncm_*.bin` (329) <= 3.0e-13, transitions & SF (39: TDA_E<L>, Full<n>_<lab>_E1*, SF_E1_*) 0 beyond the printed
  digits, `Levels.dat`, `TDA_Structure.dat` (|X|), `Basis_*.dat` (without the Ind column), `States*.dat`,
  `Summary*.dat` (without timings) equal; He4c (complete 4-phonon chain, Ortho false) and He4t (TRANS, SF) every file
  exact. Spurious E(1) (12 files `ok*`/`SIGN*`): legacy 0, dev 1.935505 (He4o), 0.518022 (O16), 0.557930 (O16w),
  -3.957443 (Ca40) MeV, all other values equal. GAUGE (20 files, O16/O16w/Ca40/He4o 2PH): another choice of
  independent rows among dependent ones in a few rows (`Ind`, `Kc`), energies equal. DIFF (47): 45 He4o amplitude
  files (2PH C, S2/S3, FULL Z: the 1- phonons of the Ortho TDA come out with other signs/rotations, a basis gauge;
  energies, N_cm, weights, B(E1) of He4o equal) and O16w `EMPM_FULL_{0+,1-}_Z.bin`. MISSING 21: NCSM (part 2).
- Proof that the remaining differences are only the spurious convention + 1e-16 HF input: scratch build with the legacy
  spurious algorithm (`Spur_Ortho`, `eig_general`, E(1) = 0, cut 1e-6): He4o all files ok (0 DIFF); the same build on the
  legacy HF files (symlinked as `HF/bin/*`) reproduces the frozen binary's O16 and Ca40 EMPM bit for bit (47 files, max 0).
- EMPM_TDA (full space) vs the dev TDA action (Reference 'HF'; He4c Ortho false, He4o, O16, Ca40 Ortho true):
  `spectra/TDA.dat` and `transitions/TDA_E{0,1,2,3}.dat` identical to the printed digits, `EMPM_TDA_E.bin` vs
  `TDA_E.bin` <= 1.6e-13 MeV (incl. the spurious E(1)).
- Error paths (`errors/`): fractional occupations (20O, ESH false), `BMF = false` (no H2b.bin), `Reference = 'HF-MBPT'`
  (non-diagonal H1b), missing reference (message + skip, exit 0), EMPM_2PH without EMPM_TDA (message + skip),
  `Reference` with '/': all give the expected message.
- Threads: 8 vs 16 threads 191 of 1805 files differ at 1e-12, 1 vs 16 threads 198 (energies <= 2.2e-11, SF Lanczos
  coefficients 9e-8, otherwise eigenvector signs / degenerate rotations / pivot choices); the frozen binary shows the
  same (8 vs 16: 194 of 1732 files, the same files and sizes): threaded MKL in EMPM_TDA (dsyevr signs) .. FULL, not a
  race of the stage-5 code. Wall time of the dev chain: 16 threads 44 s, 8 threads 70 s, 1 thread 439 s (He4c 409 s).
  Debug build (`-check all -warn all`, `build_s5_dbg`, removed afterwards): all dev steps except He4c run
  without runtime errors (array-temporary warnings of the existing EMPM code only), no compiler warnings.
- Regression of the dev TDA/RPA after the `spurious_ortho` move: bundled O16 steps 01-04 vs Julia 141 ok, 12 SIGN,
  0 DIFF (as stage 3).

Open items (part 2): NCSM on the reference (`window_make` already reads the EMPM reference H1b through the default of
`EMPM_space_make`, but its `Orb` still comes from `orbitals_make`), the `NCSM` routing in `nujuliet.f90`, the NCSM
dev inputs (`*.nml.todo`), the app layer and README/CLAUDE.md (consumers.md sections 4, 5, 7).

## 2026-11-21 — stage 5 part 2, NCSM on the dev HF reference, app layer, README (spec `consumers.md` sections 4, 5, 7)

Done (details & deviations in `docs/PORTING_NOTES.md`, section "stage 5 part 2"):
- NCSM (section 4): orbitals from `<&empm Reference>/bin/orb.bin` for `Hamiltonian = 'EMPM'` or `Basis = 'HF'`, else
  `orbitals_make(Params)` (Calc.Nmax size, was LEGACY Int.Nmax); `'EMPM'` Hamiltonian reads `H2b.bin` (missing: fatal),
  `H1b.bin` via `read_h1b_diag` (non-diagonal: fatal), `C.bin`, `rho.bin` (diagonal 0/1 within 1e-8, else fatal, values
  used as read), E_HF from `<Reference>/HF_summary.dat`; `trans_init` reads `<Reference>/bin/C.bin`; `window_make` gets
  the reference orbitals (NCSM `Window` works again); `V2b_read(Params,V,Orb_NN)` signature was already the dev one.
  Outputs unchanged in `io/<Path>/NCSM/`.
- `nujuliet.f90`: `NCSM` no longer routes `'HF'` to `HF_LEGACY` (only `HF_RPA*`, `HF_RRPA*` do); usage text names
  `io/`; dispatcher comments (incl. the stubs). Removed-action fatals not added (owner decision, stage 6);
  `&transition_currents M` default stays 0 (deviation, same result as dev's -1); `&parse` + `V2b_parse` /
  `V3b_no2b_parse` existed.
- Tests: `fortran/empm_ncsm/*NCSM*.nml.todo` renamed to `.nml` (06, 07, 08, 13, 16); new error cases
  `errors/e7_ncsm_fractional` (20O ESH false, NCSM EMPM/HO), `e8_ncsm_hf_mbpt` (Reference 'HF-MBPT', NCSM EMPM/HF),
  `e9_ncsm_fractional_window` (20O, NCSM `Window`). `tests/dev_port/README.md` updated.
- README.md (section 7: actions table, namelists, output layout, EMPM, NCSM, compatibility, verification);
  CLAUDE.md "Runs" solver mapping and "Original Julia code" only.

Verification (binary `nujuliet_s5`, `BUILD=build_s5`; clean full build: 0 warnings, identical binary to the incremental one):
- `run_empm_ncsm.sh dev 16` (work/empm_dev_t16): all 17 files run (NCSM steps 0.07-0.14 s each); all 9 error cases give
  the expected message (e7/e9: integer occupations, e8: "need an HF reference with diagonal H1b").
- `run_empm_ncsm.sh compare 16`: DIFF 47, GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614, ok* 12, **0 MISSING** (the
  EMPM rows are those of part 1). The 21 NCSM files: 17 ok (He4t 14 files incl. NN / EMPM-HO / EMPM-HF npnh with
  transitions **byte-identical** to the baseline; O16 spectra; O16w `Window` 2.4e-15), 4 SIGN (O16 `_E0/_E1/_E2/_E3.dat`:
  |d| <= 9.8e-13, the sign of the M_p, M_n amplitude columns = eigenvector phase). Cause: the 16O dev HF differs from the
  legacy HF at <= 3.6e-15; the frozen binary run on the dev HF files (copied to the legacy names) gives the same kind
  of roundoff/sign differences, and on identical HF files (4He) the dev NCSM is bit-identical. A first version that
  snapped the occupations to exact 0/1 flipped signs in 8 He4t files; removed (values used as read).
- Regression vs `work/julia_t8` (8 threads, compare.py 1e-8): bundled 411 ok + 36 SIGN; ca40_5sh 92 + 8; o16_open
  106 + 7; hf_options 14; mbpt_options 23 + 6; rpa_options 308 + 29; qp_options 150 + 14 + 3 DIFF (the known O20_HBCS
  HFB files); scrpa_options 197 + 22. 0 MISSING / EXTRA / SIZE. Same counts as stage 4.
- Legacy: `runs/test_Ca40_9sh/stage5_final_legacy` (Ca40.nml: HF, HF_RPA, HF_RRPA, 16 threads, 2:17, 0.87 GB)
  byte-identical to `stage2_mbpt_legacy/io` (`diff -rq`: no output).

Open: `tests/empm_fci_check.py` still needs the HRF tables of `HF_LEGACY`; NCSM with a missing reference is fatal
(EMPM: message + skip); stage 6 (cleanup): removed-action fatals, `&hf_efa`/`&rrpa`/legacy chain removal.

## 2026-11-21 — stage 5 independent verification (final tree, clean `make`, `./nujuliet`)

- Build: `build_s5`/`nujuliet_s5`/`build`/`nujuliet` removed, clean `make`: 0 warnings. No source changes by the verifier.
- Fresh legacy baseline (frozen `nujuliet-snap-d753f29`, 16 threads, own work dir) vs fresh dev run of `./nujuliet`
  (`empm_compare.py`, tol 1e-10): DIFF 47, GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614, ok* 12, 0 MISSING; no `ok` file
  above 1e-10. Every legacy EMPM/NCSM/HF-reference file is mapped (only `HF/Densities` unmapped; dev-only files are
  `orb.bin`, HF densities and the dev TDA). Independent checks of the non-ok classes:
  He4o (HF inputs identical): `EMPM_2PH_*_C`, `EMPM_S2_*` (X and C) and `EMPM_FULL1_1-_Z_*` agree to <= 5e-13 once the
  sign of the two flipped 1- phonons (i = 10, 11; spurious-projection algorithm) is applied to the 2-phonon / 1-phonon
  basis; TDA amplitude sign changes are whole-phonon flips in He4o/O16/O16w, and in Ca40 a consistent product of phonon
  signs and HF orbital signs (GF(2) solve, 0 contradictions; HF `C.bin` SIGN); Ca40/O16 `Basis_*` GAUGE = swapped
  pivots |a b> vs |b a>; O16w `EMPM_FULL_{0+,1-}_Z` sit on the 0+/1- 2PH bases whose Ind flags differ (GAUGE);
  O16 NCSM `_E*.dat` SIGN = whole-row phase flips of the M_p/M_n columns, other columns <= 1e-24 absolute.
  EMPM_TDA vs dev TDA (`--tda`): spectra/transitions 0 beyond the printed digits, `EMPM_TDA_E` vs `TDA_E` <= 1.6e-13 MeV;
  `EMPM_TDA_X.bin` and `TDA_X.bin` have the same sorted |X| (<= 1.1e-13), stored in another component order.
- vs `work/julia_t8` (8 threads, compare.py 1e-8): bundled 411+36 SIGN, ca40_5sh 92+8, o16_open 106+7, hf_options 14,
  mbpt_options 23+6, rpa_options 308+29, qp_options 150+14+3 DIFF (known O20_HBCS), scrpa_options 197+22;
  0 SIZE/MISSING/EXTRA. Bundled 1 vs 8 threads (`--strict-sign`, 1e-12): 447 ok.
- Legacy Ca40.nml (fresh dir, 16 threads, 2:29, 0.86 GB): `diff -rq` vs `stage2_mbpt_legacy/io` empty.
- Error cases e1..e9: all give their EXPECT text. `Solvers = 'HF', 'EMPM_TDA', 'NCSM'` in one `&run` (16O): dev HF
  (`HF/bin` identical to the separate run), no `Bin/` produced, NCSM energies equal.
- grep `src/empm`, `src/ncsm`: no `Bin/`, `_HF.bin`, `IO/`; `nj_hf_rpa*` only in comments, no such `use`.
- Review follow-up (same day): `empm_compare.py` allows SIGN only for sign-gauge files (`HF/bin/{C,H1b,H2b}.bin`,
  `EMPM_TDA_X.bin`, `EMPM_F.bin`, `TDA_Structure.dat`, NCSM `_E<L>.dat` M_p/M_n columns per row), compares text numbers
  beyond their own printed unit (was: the coarsest unit of the file, up to 1e-1 in the NCSM tables), checks the energies
  of GAUGE `EMPM_2PH_*_C` files, exits 1 on DIFF/SIZE/MISSING; `--tda` reports missing files as MISSING; dead
  `spurious_masks()` removed. EMPM_2PH..SF check the reference first; `'EMPM'` stops after a skipped EMPM_TDA (PORTING_NOTES
  stage 5 part 2 item 7). Clean `make` 0 warnings; `run_empm_ncsm.sh dev 16` vs `work/empm_baseline`: DIFF 47 (He4o 45,
  O16w 2: the known amplitude-basis files), GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614, ok* 12, 0 MISSING (classes
  per file identical to before); `--tda` worst 1.5e-14, 0 MISSING; error cases e1..e11 ok; bundled vs `julia_t8`
  (1e-8): 411 ok + 36 SIGN, 0 DIFF.

## 2026-11-21 — stage 6, cleanup: removal of the LEGACY code (owner approved; details in `docs/PORTING_NOTES.md`, stage-6 section)

Removed (as dev did upstream):
- Files (11, 5809 lines): `src/hf/nj_hf_legacy.f90`, `nj_hf_legacy_driver.f90`, `nj_hf_mbpt_legacy.f90`, `nj_hf_efa.f90`,
  `src/rpa/nj_hf_rpa.f90`, `nj_hf_rpa_driver.f90`, `nj_hf_rpa_export.f90`, `nj_hf_rpa_external.f90`, `nj_hf_rrpa.f90`,
  `nj_hf_rrpa_driver.f90`, `nj_hf_rrpa_export.f90`; Makefile `SRC_HF`, `src/rpa/files.mk`.
- Actions `HF_LEGACY`, `HF_RPA`, `HF_RPA_transition_{densities,currents}`, `HF_RRPA`, `HF_RRPA_transition_{densities,currents}`
  and the auto-routing of `HF` to `HF_LEGACY`: now fatal with a hint, checked for the whole `&run` before any action of
  the file runs (`'HF_RPA' was removed upstream, use 'TDA', 'RPA'`; `'HF_RRPA' was removed upstream, use 'HF_ScRPA'`; ...).
- Parameters/namelists: `hf_efa_params_t`/`&hf_efa`, `rrpa_params_t`/`&rrpa` (now: WARNING, group ignored), the `&hfb`
  dummies `Broy_Qmax`/`Broy_Bmax` (now: fatal with the hint to delete them); text pre-scan of the input file in `nj_input`.
- Types/routines/aliases used only by the legacy code: `obdm_iteration_pairs_t`, `pnarray4_t`; `nj_orbitals` LEGACY
  `orbitals_ph_make` specific, `orbitals_ph_list`, `orbitals_one_phonon_make`, `read_h_HF_diag`, HRF `orbitals_export`
  (generics reduced to the dev/reference forms); `nj_o2b` 3-block Pandya form, `O2b_export_HRF`; `nj_obdm` `OBDM_*`
  aliases; `nj_linalg` `eig_general_complex`, `isposdef`, `isposdef_c`, `zmhm` (+4 LAPACK interfaces); `nj_utils`
  `Verbose` (only HF-RRPA switched it off). All `LEGACY` markers / stale comments in `src` updated.
- `&transition_currents M` default -1 as in dev (the 0 was kept only for the LEGACY `HF_RPA_transition_currents`).
- Line count: `src` 41539 -> 35257 (-6282 net; 130 insertions, 6423 deletions incl. Makefile / files.mk).
- Kept: Fortran-only EMPM, NCSM, `CM_Lawson`, `PARSE`; `HRF` dead inputs (as dev); `Use_Int` of the interaction readers
  (dev `Int = true`).
- Docs: README (actions table, "Removed actions and input" migration table, no legacy namelists, `empm_fci_check.py`
  note), CLAUDE.md ("Runs", "Original Julia code", "Comparing with Julia"), PORTING_NOTES (obsolete items marked
  "(removed in stage 6)", new stage-6 section), `tests/dev_port/README.md` (legacy baseline inputs need the frozen binary,
  stage-6 regression recipe). New error cases `fortran/empm_ncsm/errors/e12_removed_hf_rpa`, `e13_removed_broy`.
  `fortran/empm_ncsm_legacy/` stays unchanged (uses `HF_LEGACY`: frozen binary only).

Verification (frozen binary `nujuliet-snap-ec14d69/nujuliet` vs the final `./nujuliet`, clean `make`, 0 warnings):
- All dev_port cases, 8 threads (`work/s6_frozen_t8`, `work/s6_final_t8`): `diff -rq` of the io trees **byte-identical**
  (1433 files; only `run.log`/`time.log` differ: timing lines and the last printed digit of 42-44 full-precision energy
  log lines, the same kind as frozen vs frozen run twice, `work/s6_frozen2_t8`: OpenMP reduction order).
- `empm_ncsm` (`run_empm_ncsm.sh dev 16`, TAG `_s6frozen`, `_s6frozen2`, `_s6final`): io trees byte-identical apart from the
  seconds column of `EMPM/{3,4}Phonon/Summary.dat` (also differs frozen vs frozen); logs differ in timings, RSS, the
  free-memory limit of the NCSM stored H and last-digit E_HF prints (noise floor identical). No MKL-thread differences at
  equal thread count in these runs. All 13 error cases (e1..e13) give their EXPECT text. `compare` vs `work/empm_baseline`:
  DIFF 47, GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614, ok* 12 (as stage 5).
- vs `work/julia_t8` (compare.py 1e-8): bundled 411+36 SIGN, ca40_5sh 92+8, o16_open 106+7, hf_options 14, mbpt_options
  23+6, rpa_options 308+29, qp_options 150+14+3 DIFF (the known O20_HBCS files), scrpa_options 197+22; 0 SIZE/MISSING/EXTRA.
- Removed actions / input: all 7 removed actions give their hint (exit 1, nothing written); `&rrpa`/`&hf_efa` warn,
  `Broy_Qmax` fatal; text in comments ignored. Default `M` (-1) of `TDA_/RPA_transition_currents`: 16 files identical
  to the `M = 0` run.
- Builds: intel clean (0 warnings), stub chain-A build (`SRC_RPA=src/rpa/stub/...`) ok. `make COMPILER=gnu` fails in
  `nj_ph_external.f90` (local `hc, m_n, m_p` parameters clash with `nj_kinds`; gfortran rejects, ifx accepts) — already
  broken at ec14d69 (stage 3 file, unchanged here). gfortran `-Wall` on the edited core/app modules: no new warnings.

Open: `tests/empm_fci_check.py` still needs the HRF tables of a pre-stage-6
binary; `runs/test_Ca40_9sh` legacy regression and `fortran/empm_ncsm_legacy/` baseline only with frozen binaries.

Independent verification, 2026-11-21 (fresh work dirs, clean serial `make`, final tree):
- Build: intel 0 warnings; no references to the removed modules/types/actions/namelists in `src`, Makefile, `files.mk`
  outside the hint code. **gfortran build fixed**: `nj_ph_external.f90` uses `m_n`, `m_p` of `nj_kinds` (same values,
  `hc = HbarC`); the link then failed on the MKL-only `dsygvd_64`/`dsyevd_64`/`dgemm_64`/`dgemv_64` of `nj_empm_eqm2`
  (stage 5): new `src/core/nj_lapack64_ref.f90` (wrappers of the default-integer routines, gnu only, Makefile
  `BLASTHR_SRC`); `V2b_parse`/`V3b_no2b_parse` opened the input inside the parallel region through a deferred-length
  `In_File` that gfortran passes as garbage (stage 1): now `Params%Int%NN_File`/`NNN_File`. `make COMPILER=gnu` 0 warnings;
  gnu vs Julia (1e-8), all 8 cases: 0 DIFF/SIZE/MISSING/EXTRA apart from the 3 known O20_HBCS files (SIGN counts differ
  from intel: LAPACK gauge). README "Building" updated.
- Intel `./nujuliet` vs frozen ec14d69, 8 threads, 8 dev_port cases: io trees byte-identical (1433 files). vs Julia: bundled
  411/36/0, ca40_5sh 92/8/0, o16_open 106/7/0, hf_options 14/0/0, mbpt_options 23/6/0, rpa_options 308/29/0, qp_options
  150/14/3 (O20_HBCS), scrpa_options 197/22/0.
- `empm_ncsm` (dev 16, both binaries): 1830 io files, identical apart from the seconds column of 4 `Summary.dat`; all 13
  error cases ok (frozen fails e12 as expected); compare vs baseline DIFF 47, GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614,
  ok* 12. Removed actions (7): hint, exit 1, nothing written; `&rrpa`/`&hf_efa` warn, io identical to the plain run;
  `Broy_Bmax` fatal, `Broy_Qmax` in a comment ignored; copy of `runs/test_Ca40_9sh/Ca40.nml`: `&rrpa` warning + HF_RPA
  hint, HF not run.
- Review follow-ups (coordinator): the removed-namelist scan blanks quoted strings (a value like 'out/broy_qmax' no longer
  matches); all input files of a call are read and checked before any action runs (a removed action in b.nml stops
  `nujuliet a.nml b.nml` before a.nml computes anything; tested); README `&transition_currents` sample M=-1.
  bundled 411/36/0 and scrpa_options 197/22/0 vs Julia after the change.

## 2026-11-21 — upstream dev a69ecf0 "HF-MBPT hotfix" ported (new parity target)

Done:
- Julia: `NuJuliet_dev` moved to origin/dev a69ecf0 (only `src/HF/HF_MBPT.jl` changed vs f3f7a33). `port-ref` rebased onto
  a69ecf0: F1 dropped (identical upstream change: nn TBCF blocks use `nO`), F2 3d41af6, F3 f3999d2, F4 841b9a6, F5 47f840d
  (diff vs a69ecf0: 4 files, +11/-5, as before minus F1). Old f3f7a33-based commits: F1 204b5be .. F5 55c5880.
- Fortran (`src/hf/nj_hf_mbpt_common.f90`, `nj_hf_mbpt_tbcf.f90`, `nj_hf_mbpt.f90`): `HF_MBPT3_V2b_Pandya` pn exchange
  orderings (-1,+1,+1,-1), (+1,-1,-1,+1) with `-phase(...)`; `HF_MBPT_TBCF` `Phase_ab/cd = -phase(j+j-J)` for every isospin;
  `HF_MBPT_OBDM` proton pp block `P_hr = mod(l_h + l_r,2)+1`; log "HF-MBPT(2+3) energy  ...   E^(tot)", summary
  "LO HF-MBPT(2) ... per nucleon" and "E_0^(3) / A ... per nucleon", export print "Constructing HF-MBPT 1-body Hamiltonian
  (correlated OBDM) & HF residual 2-body interaction in the NAT basis ..."; F1 comments now say "upstream".
- Numerical effect: none. P_hr: the loop runs only for `l_p == l_h` (the 1p1h block requires equal l, j), so
  `mod(l_p+l_r,2) == mod(l_h+l_r,2)`. Pandya/TBCF phases: the flipped elements connect pi-pi with nu-nu ph pairs; the
  E^(3)_ph ring crosses an even number of times and in the TBCF ph diagram the Pandya flip cancels the new Phase_ab/cd of
  pn elements (PORTING_NOTES). Julia old vs new refs: all binary files byte-identical, only 2 text lines of each MBPT
  `HF_summary.dat` differ. E^(3) unchanged, e.g. Ca40 (ca40_5sh) E^(3) = -5.29683304 MeV (pp -1.57694752, ph -0.152934047,
  hh -3.56695147), E^(3)/A = -0.13242083; O16_MBPT E^(3)/A = -0.18169379; O28_TBCF -0.05700789; He4 -0.12404847 (old =
  new in Julia and in the port).
- Julia refs regenerated (port-ref 47f840d): `work/julia_t8/{bundled,ca40_5sh,mbpt_options,rpa_options,scrpa_options}`
  (MBPT summaries changed) and, as a check, `o16_open`, `hf_options` (no MBPT: identical to the old refs); `work/julia_t1/
  bundled`. `qp_options` (no HF-MBPT) kept. Wall (4 lanes x 8 threads concurrently): bundled 3:30, ca40_5sh 3:06,
  scrpa_options 2:25, rpa_options 2:09, o16_open 1:26, mbpt_options 0:57, hf_options 0:38, bundled t1 4:46.
- Tests: intel clean build 0 warnings; `make COMPILER=gnu BUILD=build_gnu EXE=nujuliet_gnu` 0 warnings (deleted). vs new Julia
  refs (8 threads, 1e-8): bundled 411/36/0, ca40_5sh 92/8/0, o16_open 106/7/0, hf_options 14/0/0, mbpt_options 23/6/0,
  rpa_options 308/29/0, qp_options 150/14/3 (the known O20_HBCS), scrpa_options 197/22/0 (ok/SIGN/DIFF; 0 SIZE/MISSING/
  EXTRA). gauge_check (O28_TBCF, O16_ScH, O16_MBPT HF-MBPT/bin, O28_TBCF HF/bin): <= 1.3e-13. bundled Fortran 1 vs 8
  threads `--strict-sign --tol 1e-12`: 447 ok; Fortran t1 vs Julia t1: 428/19/0. `empm_ncsm` (dev 16): 13 error cases ok;
  compare vs `work/empm_baseline`: DIFF 47, GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614, ok* 12 (unchanged).
- Docs: CLAUDE.md, README.md, tests/dev_port/README.md, PORTING_NOTES (target, F1 upstream, hotfix note), UPSTREAM_BUGS
  ("per nucleon" fixed upstream, ScH=false export documented as intended), hf.md bug list annotated.

## 2026-11-21 — NCSM: missing reference non-fatal (owner decision)

- `NCSM` with `Hamiltonian = 'EMPM'` or `Basis = 'HF'` calls `EMPM_reference_check` first: a missing
  `io/<Path>/<Reference>/bin` prints "Error! ... No precomputed reference ..." and skips the NCSM (exit 0, no
  `io/<Path>/NCSM` created), as TDA/RPA/EMPM and dev. A reference without `H2b.bin` stays fatal (as EMPM). NCSM with
  `NN`/HO needs no reference and is unchanged. Closes the open item of the stage-5 entry.
- New error case `fortran/empm_ncsm/errors/e14_ncsm_no_reference`. `run_empm_ncsm.sh dev 16`: all 14 error cases give
  their EXPECT text; `compare` vs `work/empm_baseline`: DIFF 47, GAUGE 20, SIGN 11, SIGN* 4, SKIP 6, ok 1614, ok* 12
  (unchanged). Builds: intel and `COMPILER=gnu` 0 warnings.

## 2026-11-21 — fix F2 revised (RPA_sigma hpph uses Y in every branch)

- Finding: the first F2 (pn branch `X` in `RPA_sigma_hpph`, as pp/nn) was wrong. Only the Ket <= Bra triangle of Sigma
  is computed, so the hpph element `<h p|Sigma|p' h'>` must equal its transpose, the phhp element (Y Y*). Full-square
  check on 12C (bundled, RPA TBCF): pn with Y 1.7e-18, pn with X 0.49; pp with X (dev) 7.3e-2, pp with Y 3.5e-18. The
  real dev bug is the pp/nn `X` (RPA_sigma.jl:153-160, 273-280). The term enters the computed triangle only when a
  particle orbital precedes a hole orbital of the same isospin (Z or N = 6, 14, 28, 50, 82, 126, or fractional
  occupations); 4He/16O/40Ca never reach it, so no existing step saw the difference.
- Fortran `nj_rpa::RPA_sigma` `s_hpph` uses Y (all branches). `port-ref`: old F2 dropped, F3 fff2142, F4 0faf549,
  F5 76a58f5, new F2 ddcd86d (pp/nn pass Y).
- 12C RPA TBCF E_corr (density): dev -30.712444, old F2 -6.742027, new F2 -38.331296 MeV (quasi-boson -30.898290 in all).
- Regenerated every Julia reference (julia_t8 all cases, julia_t1 bundled) with the new port-ref: identical to the
  old ones at 1e-12 `--strict-sign` except `scrpa_options` O20/ScRPA-frac `RPA_summary.dat` (fractional occupations:
  density E_corr 22.108791 -> 8.290909 MeV) and the new 12C files.
- New `rpa_options` steps 13-15 (`io/C12`: HF, RPA TBCF, HF_ScRPA; with TBCF = true the 12C HF_ScRPA diverges in dev
  itself, unstable initial RPA, so the plain HF_ScRPA tests the final Sigma via E_corr density).
- Fortran vs Julia (8 threads, 1e-8): bundled 411/36/0, ca40_5sh 92/8/0, o16_open 106/7/0, hf_options 14/0/0,
  mbpt_options 23/6/0, rpa_options 368/35/0, qp_options 150/14/3 (known O20_HBCS), scrpa_options 197/22/0. 12C
  sigma.bin 6.5e-15 up to the HF-orbital sign; E_corr equal to all printed digits. intel and gnu builds 0 warnings.
- Affected earlier results: E_corr (density) and sigma.bin of HF_ScRPA / RPA TBCF runs of nuclei with such shell
  closures, e.g. the Pb208 HF-ScRPA run (no TBCF, so its iteration, E_QBA and the TDA/RPA results are unaffected; its
  density E_corr -690.923000 MeV was computed with the old F2 and must be recomputed).
- Upstream reports for F2-F5: `docs/dev_port/upstream_reports/` (not filed).
