
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
