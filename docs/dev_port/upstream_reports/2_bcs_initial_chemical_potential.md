# BCS: initial chemical potential is 0.5*(Eps_hole - Eps_particle) instead of the mid-gap value

**Version:** dev a69ecf0, `src/BCS/BCS_chemical_potential.jl:90`, `HF_BCS_initialize_chemical_potential`

## Summary

The function finds the highest occupied single-particle energy, `Eps_hole`, and the lowest empty one,
`Eps_particle`, of the HF reference. It then sets

```julia
Lambda = pnFloat(0.5 * (pEps_hole - pEps_particle), 0.5 * (nEps_hole - nEps_particle))
```

That expression is minus half of the shell gap, not a point in the gap. Its value depends only on the gap size,
not on where the Fermi level is. The midpoint is

```julia
Lambda = pnFloat(0.5 * (pEps_hole + pEps_particle), 0.5 * (nEps_hole + nEps_particle))
```

## Why the starting value changes the result

For a closed-shell species the pairing gap goes to zero, and the particle number is then exact for any Lambda
inside the shell gap. So `dZ` (or `dN`) is zero from the first iteration, `BCS_chemical_potential` skips the
update (`abs(dA.p) > epsilon` is false), and the starting Lambda becomes the final one.

With Delta = 0 the quasiparticle energies are `|e_a - Lambda|`. Two-quasiparticle states of the particle-hole
type do not depend on Lambda. Pair states, with two quasiparticles above or two below Lambda, shift by
-/+ 2 delta_Lambda. This moves QTDA, QRPA and pnQTDA levels and reorders the spectra.

If the wrong value lands outside the gap (small gap, deep Fermi level), the first iteration also starts with the
wrong particle number.

## Reproduction

Bundled interaction (hw = 16, Nmax = 3, N2max = 6, N3max = 9). 20O:

```julia
Calculation_Parameters(A = 20, Z = 8, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9,
    BCS  = BCS_Parameters(ScBCS = true, A0 = 16, Z0 = 8),
    QTDA = QTDA_Parameters(Reference = "BCS"), QRPA = QRPA_Parameters(Ortho = false, Reference = "BCS"))
```

Run `BCS` and then `QTDA`/`QRPA`. The protons (Z = 8) form a closed shell.

| | current | mid-gap fix |
|---|---|---|
| pLambda (`BCS_summary.dat`) | -7.583 MeV | -10.101 MeV |
| proton qp energy of 0d5/2 | 5.134 MeV | 7.651 MeV |
| QTDA 0+ levels (MeV) | 10.020, 14.239, 19.511, 24.628 | 14.239, 14.477, 15.055, 24.628 |

The particle-hole-type levels (14.239, 24.628) are the same in both. The pair levels move by
2 x 2.518 = 5.035 MeV (15.055 to 10.020, and 14.477 to 19.511). The QRPA spectrum and transition tables change in
the same way.

The default 18O BCS case changes only slightly: pLambda -6.3547 vs -6.3020 MeV, E_tot differs by 3e-7 MeV.

## Proposed fix

Replace `-` with `+` in both components on line 90.
