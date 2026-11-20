# Possible upstream bugs (NuJuliet dev f3f7a33) ported faithfully

Ambiguous behaviour of the Julia dev code that the Fortran port reproduces on purpose (bit parity with the
patched reference `port-ref`) and only reports here. Clear bugs are fixed instead (F1..F5, see `PORTING_NOTES.md`).

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
  NAT basis). Ported as is (owner decision).
* `HF_MBPT_energy` summary (`HF_MBPT.jl:744`): the line `E_0^(3) / A ... correction to ground-state energy` lacks
  "per nucleon" (cosmetic, kept for file parity).
* `HF_MBPT3_V2b_Pandya` (`HF_MBPT.jl:100-114`): for `t_a + t_b = -2/+2` the isospins of c, d are not checked (all
  callers conserve charge, so no wrong value is produced). Ported as is.
* MBPT denominators would vanish for an orbital that is both particle and hole (fractional occupation); this is
  prevented by `HF_MBPT_occupation_check` (non-idempotent HF density => no MBPT, e.g. 20O in `o16_open`). No action.
