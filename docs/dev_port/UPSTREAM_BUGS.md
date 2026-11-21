# Possible upstream bugs (NuJuliet dev f3f7a33, re-checked for a69ecf0) ported faithfully

Ambiguous behaviour of the Julia dev code that the Fortran port reproduces on purpose (bit parity with the
patched reference `port-ref`) and only reports here. Clear bugs are fixed instead (F1..F5, see `PORTING_NOTES.md`; F1 is upstream since dev a69ecf0).

## Group common A (core)

* `orbitals_make` (`src/common/orb.jl:36-79`): after the last (partially) filled shell the remaining orbitals get
  `Occ = min(1.0, (A - Z - N_occ)/j_hat)` where `A - Z - N_occ` is a rounding residual of about ±1e-16 instead of 0
  (only for open shells with a fractional occupation, e.g. 2/6). The tiny (possibly negative) occupations enter every
  density built from `Orb`. Ported bit for bit; a clean version would set `Occ = 0` once the filling is complete.
* `orbitals_1p1h_make` (`orb.jl:348-442`): with fractional occupations an open-shell orbital whose occupation lies in
  `[Tol_par, Tol_hol]` (e.g. 1/3 with the defaults 0.1 / 0.7) is neither a particle nor a hole, so it is silently dropped
  from the TDA/RPA 1p1h space. Ported as is.
* `OBDM_natural_orbitals` (`orb.jl:551-650`): the returned occupations keep the deterministic splitting `+1e-10*a`
  added to the diagonal before the diagonalization. Ported as is (harmless at the 1e-10 level).
* `orbitals_make` with `A - Z` or `Z` larger than the capacity of the basis loops forever (`while N_fill`); the port
  stops with an error (robustness guard, no Julia patch).

## Stage 2, HF (dev `src/HF/HF_orbitals.jl`, `HF_energy.jl`)

* `HF_orbitals_make` with `ESH = true` (`HF_orbitals.jl:76-88`, neutrons 105-117): the "closed-subshell" error is
  thrown as soon as the slot `a_max` is *examined* (its SPE is below the best eligible candidate found so far) and is
  too large, even when a smaller eligible slot exists later in the energy order. Observed: 20O (bundled interaction,
  Nmax 3) with `ESH = true` stops in the first iteration in both codes (0f7/2 = slot a_max lies below the
  remaining eligible 1p1/2), while 18O works. Ported as is (same trigger, `fatal` with the same text).
* `HF_orbitals_make`: when no unfilled slot is eligible, `a_last = 1` is used and the occupation of slot 1 is
  overwritten (only with invalid A/Z or ESH corner cases). The port stops with an error instead (robustness guard).
* `HF_Parameters.HRF` is accepted but never read (the HRF export was removed); kept as an ignored namelist field.
* `HF_energy` summary text and `HF_SPS_summary` table are unchanged from base; the `T.p` used for the neutron kinetic
  term (`HF_energy.jl`, kinetic line) is harmless because `T.p == T.n`.

## Stage 2, HF-MBPT (dev `src/HF/HF_MBPT.jl`)

* `HF_MBPT_export` with `ScH = false` (`HF_MBPT.jl:2288-2295`): the comment says "take HF Hamiltonian", but
  `h_N_NAT = H1b(...,Rho_NAT)` is built from the **correlated** (HF-MBPT natural-orbital) density, while `V` is the HF
  residual interaction `V_NN_res` (HF density). Possibly intended: `h = C_NAT^T (C h_N C^T) C_NAT` (HF mean field in the
  NAT basis). Ported as is (owner decision). **Resolved in a69ecf0**: the comment and printout now state the correlated
  OBDM 1-body Hamiltonian as intended.
* ~~`HF_MBPT_energy` summary (`HF_MBPT.jl:744`): the line `E_0^(3) / A ... correction to ground-state energy` lacks
  "per nucleon"~~ (fixed upstream in a69ecf0, ported).
* `HF_MBPT3_V2b_Pandya` (`HF_MBPT.jl:100-114`): for `t_a + t_b = -2/+2` the isospins of c, d are not checked (all
  callers conserve charge, so no wrong value is produced). Ported as is.
* MBPT denominators would vanish for an orbital that is both particle and hole (fractional occupation); this is
  prevented by `HF_MBPT_occupation_check` (non-idempotent HF density => no MBPT, e.g. 20O in `o16_open`). No action.

## Stage 3, quasiparticle chain (dev `src/BCS`, `src/HFB`, `src/QTDA`, `src/QRPA`, `src/pnQTDA`, `src/common/qpO2b*.jl`)

* BCS reference HF (`BCS_solver.jl:96-133`): `HF_solve(Params,Orb_HF,...)` refills the occupations with the **target**
  (A,Z) every iteration (dev `HF_orbitals_make`), so the "reference closed-shell" HF uses the target's fractional
  filling; `Params_Ref` only seeds the first density and enters `HF_BCS_initialize_chemical_potential` (whose `== 1`/
  `== 0` Float64 tests skip fractionally occupied levels). Ported as is (open question to the owner).
* `HFB_diagonalize_BCS` (`HFB_solver.jl:409`) calls `HF_orbitals_make`, which throws for `HF.ESH = true` and an open-shell
  target (HFB with BCS pairing is meant for open shells) and mutates HFB's `Orb` occupations (harmless). Ported (fatal).
* pnQTDA pair set (`pnQTDA_orbitals.jl`) is not truncated by N2max, so `qpO2b_1111_pn` throws a `KeyError` for
  `N2max < 2*Nmax`; already `qpO2b` itself needs `N2max >= 2*Nmax` (direct O2b terms over all pairs up to 2*Nmax).
  The port stops with a clear message in both places.
* `qpO2b_export` loops `J = 0:2Nmax+1` over arrays sized `N2max+1` (pre-existing; port note B.3); the new import size
  check uses the port's channel range.
* `BCS()` and `HFB()` both `rm -r io/<Path>`, so `Reference = "BCS"` and `"HFB"` results cannot coexist under one Path.
* QTDA/QRPA plot export writes only N of the N+1 sorted rows (drops the highest level); pnQTDA writes all N rows.
  `QTDA_plot.dat` vs `QRPA_Plot.dat` case. `HFB_Broyden.jl` header says `Broy_Amax` 0.9, the struct default is 0.8.
* QRPA with complex eigenvectors of (A-B)(A+B) stops with an `InexactError` (e.g. 22Ne, bundled interaction, Nmax 3);
  the port stops with an explicit error.
* HFB with `Pairing = "BCS"` for 20O (closed proton shell): the converged solution is sensitive at the 1e-8 level
  (Tol = 1e-7; the proton number dispersion `dZ = sqrt(|...|)` of a vanishing quantity prints 4.5e-5 in Julia, 5.5e-5
  in the port, 5.1e-5 in a debug build). Not a bug in either code; noted for the test tolerance.

## Stage 3, TDA & RPA part 1 (dev `src/TDA/*`, `src/RPA/*`; track RPA)

* `RPA_OBDM_oY2` (`RPA_rho.jl:1-73`) sums `|Y|^2` over all states, including unstable (complex-energy) ones, while
  `RPA_energy_QBA` excludes them (the comment claims the same set). Ported as is.
* `RPA_sigma` classifies particles / holes with inclusive tests `O <= Tol_par`, `O >= Tol_hol`, while
  `orbitals_1p1h_make` uses strict ones (`< Tol_par`, `> Tol_hol`). Ported as is.
* A spurious 1- RPA state with a real, small `E_s` is a normalized solution that enters the OBDM, the QBA energy and
  Sigma (only the E1 moments skip `nu = 1`); `|Y|^2 ~ 1/E_s` can be large. Ported as dev (owner decision).
* `RPA.jl:34-38` deletes `io/<Path>/<Name>` unconditionally and has the dead `if !isdir(io/<Path>) mkdir; return`.
  The port refuses `Name == Reference` in the reader, otherwise follows dev.
* `TDA.jl:49` checks for a marker `.TDA_output` that is never written (robustness fix in the port, see PORTING_NOTES).

## Stage 3, TDA & RPA part 2 (dev `TDA_external.jl`, `RPA_external.jl`, `RPA_allocate_SE`; track RPA)

* Transition current maps, TDA vs RPA (`TDA_external.jl:1223-1283`, `RPA_external.jl:894-942`): the TDA pair maps use
  `(-1)^((j_b - m_b)/2)`, the RPA ones `(-1)^((j_b + m_b)/2)` (opposite sign, `m_b` odd), and the level amplitudes are
  `g (-1)^J Gamma` (TDA) vs `-Re(g Gamma)/hat(J)` (RPA), so the TDA and RPA maps of the same multipole differ by a
  factor `(-1)^J hat(J)` beyond the physics (spec 7.4-6). Ported as is.
* Current map grid: a map point is snapped to the nearest `(r, theta)` node with `round` (ties to even); on the
  `x = 0` column `theta = pi` and `pi/dtheta = 127.5` is an exact tie, so the chosen node (and the sign of the printed
  `x = r sin(theta)`) depends on the last bit of `x_min + ind_x*dx` and of `atan`. Fragile but deterministic in Julia;
  the port reproduces it (PORTING_NOTES).
* `TDA_transition_operator_densities` allocates two `a_max x a_max x 8193` Float64 buffers (0.5 GB each for 13 shells)
  that only hold `<a|r^X Y_J|b>` products; `TDA/RPA_transitions_recalculate` create `transitions/E0..E3` only when
  `transitions/` is missing (non-recursive `mkdir`). Performance / robustness only.
* `RPA_allocate_SE`: the keyword arguments `Mode` (`"QP"`, `"Fixed"`), `Regularization` (`"Imaginary"`, `"None"`),
  `Gamma`, `Omega`, `dE_bin` and `Chunk` cannot be reached (both callers fix `Mode = "OnShell"`); `RPA_allocate_SE_Old`
  is unused. The OnShell evaluation builds a full `Sigma(w_i)` per 1p-1h row (N_ph^3 N_bin work, N_ph^2 N_bin memory
  per channel) where the needed entries follow from one product per block of configurations (see PORTING_NOTES).
  With `TBCF = true` the self-energy is computed twice (A and the dressing energy), with identical results.

## Stage 4 (HF-ScRPA building blocks)

* `RPA_OBDM_oY4` (`RPA_rho.jl:356`): `M = -0.5 * transpose(X) * (W .* conj(X))` is `-1/2 X^T W X^*`, not the
  Hermitian `X^+ W X`. For real amplitudes both agree; for complex amplitudes (unstable RPA blocks) the result depends
  on the arbitrary phase of each eigenvector: with `X -> X D`, `Y -> Y D` (D diagonal phases) `G` changes by
  `Y (D^2 M D^*2 - M) Y^+`. Verified with the port (generic phases change the OBDM at the 1e-3..1e-2 level on 4He;
  phases with `D^2 = -1` leave it invariant). Ported as is (bit parity with dev); the Fortran and Julia complex
  eigenvectors come from the same LAPACK `geev` convention, but parity in unstable cases relies on that.
* Spec B1..B9 of `scrpa.md` section 7 stand; B1 is fixed by F2 (stage 3), B2 is a fatal in the input reader.
