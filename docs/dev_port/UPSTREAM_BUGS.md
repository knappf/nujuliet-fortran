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
