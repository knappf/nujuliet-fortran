# NuJuliet dev: summary of bugs found during the Fortran port

**Version checked:** dev `a69ecf0` ("HF-MBPT hotfix"). Each bug has its own report in this folder, with the
reproduction and the proposed fix.

| # | Where | Bug | Effect on results | Severity |
|---|---|---|---|---|
| 1 | `src/RPA/RPA_sigma.jl:153-160, 273-280` | pp and nn branches pass `X` to `RPA_sigma_hpph`; should be `Y`, as the pn branch does | Σ is not symmetric; wrong `sigma.bin`, density E_corr and all TBCF quantities for some nuclei | **high** |
| 2 | `src/BCS/BCS_chemical_potential.jl:90` | initial chemical potential is `0.5*(e_hole - e_particle)` instead of `0.5*(e_hole + e_particle)` | wrong λ for closed-shell species; shifted quasiparticle energies and QTDA/QRPA spectra | **high** |
| 3 | `src/common/T1b.jl:3-4` | operator sized by `Params.Int.Nmax`, but the solvers pass `Calc.Nmax` orbitals | out-of-bounds read under `@inbounds`; results correct today, but crashes with bounds checking (`Pkg.test`) | medium |
| 4 | `src/common/orb.jl:623-631` | neutron natural orbitals take l, j from the proton eigenvectors `pC` | none in practice (the block structure keeps p/n consistent), but the assumption is unchecked | low |

## 1. `RPA_sigma`: X instead of Y in the pp/nn h-p-p-h terms

**Report:** [1_rpa_sigma_pp_nn_hpph.md](1_rpa_sigma_pp_nn_hpph.md)

- **Why it is wrong:** only the `Ket <= Bra` triangle of Σ is computed. The element ⟨h p|Σ|p′ h′⟩ (hpph) must
  therefore equal its transpose ⟨p′ h′|Σ|h p⟩, which is computed by phhp with Y. Both use the same 1p1h pairs, so
  they agree only with Y.
- **Numerical check, 12C:**
  - pp branch with X: asymmetry 7e-2.
  - pp branch with Y: asymmetry 1e-18.
  - pn branch (Y, current code): symmetric.
- **When it matters:** the term enters the computed triangle only when a particle orbital comes before a hole
  orbital of the same isospin. That happens for Z or N = 6, 14, 28, 50, 82, 126 (12C, 28Si, 56Ni, Sn, 208Pb) and
  with fractional occupations. 4He, 16O and 40Ca are unaffected.
- **Effect:** 12C RPA (TBCF), E_corr (density) is -30.71 MeV now and -38.33 MeV with the fix. The quasi-boson
  E_corr is unchanged.
- **Fix:** pass `Y` at the 8 call sites.

## 2. BCS: wrong initial chemical potential

**Report:** [2_bcs_initial_chemical_potential.md](2_bcs_initial_chemical_potential.md)

- **Why it is wrong:** `0.5*(e_hole - e_particle)` is minus half of the shell gap, not the middle of it.
- **Why it changes results:** for a closed-shell species the pairing gap is zero and the particle number is exact
  for any λ inside the gap. λ is therefore never updated, and the wrong starting value becomes the final one.
- **Effect:** 20O (ScBCS, A0 = 16): λ_p is -7.58 MeV instead of -10.10 MeV. Pair (pp/hh) 2qp states shift by
  5.04 MeV and the QTDA/QRPA spectra reorder. Particle-hole states are unchanged.
- **Fix:** change `-` to `+`.

## 3. `T1b`: out-of-bounds read

**Report:** [3_t1b_out_of_bounds.md](3_t1b_out_of_bounds.md)

- **The bug:** `a_max` is taken from the interaction Nmax, while the solvers pass the shorter list of
  calculation orbitals.
- **Effect:**
  - The leading block of T is correct, and it is the only part used, so current results are fine.
  - The rest of T is built from memory past the end of the array.
  - `--check-bounds=yes` gives `BoundsError: attempt to access 6-element Vector{Orb1B} at index [7]`.
- **Fix:** `a_max = length(Orb)`.

## 4. `OBDM_natural_orbitals`: neutron l, j taken from `pC`

**Report:** [4_obdm_natural_orbitals_neutron_lj.md](4_obdm_natural_orbitals_neutron_lj.md)

- **The bug:** l and j are averaged over the proton eigenvectors only.
- **Why it is harmless now:** after the block cleanup and the Hungarian reordering, column `a` of both `pC` and
  `nC` lies in the (l, j) block of `Orb[a]`.
- **Fix:** take l and j from `nC` as well and raise an error on a mismatch, or simply use `Orb[a].l` and `Orb[a].j`.
