# RPA_sigma: pp and nn branches pass X to RPA_sigma_hpph, should be Y (Sigma not symmetric)

**Version:** dev a69ecf0, `src/RPA/RPA_sigma.jl`

## Summary

`RPA_sigma` builds Sigma only for `Ket <= Bra` and stores that packed triangle, so every element also stands
for its transpose. For that to work, the kernel used for a given configuration must give the same value as the
kernel used for the transposed configuration.

The h-p-p-h configuration `<h p|Sigma|p' h'>` is computed with `RPA_sigma_hpph`. Its transpose
`<p' h'|Sigma|h p>` is computed with `RPA_sigma_phhp`, which is called with `Y`. Both kernels look up the same
1p1h pairs `(p',h)` and `(p,h')`. They therefore agree only when `RPA_sigma_hpph` is also called with `Y`.

- The **pn** branch passes `Y` (line 209), which is correct.
- The **pp** branch (lines 153-154 and 159-160) and the **nn** branch (lines 273-274 and 279-280) pass `X`, which is wrong.

In quasi-boson terms both orderings reduce to `<B+_{ph'} B_{p'h}>`, which is proportional to
sum_nu Y_{ph',nu} Y*_{p'h,nu}, so `Y` is the right amplitude.

## When it matters

The h-p-p-h case lands in the computed triangle only when a particle orbital comes before a hole orbital of the
same isospin in the orbital list. The list is ordered by N, then l, then j, so this happens when the Fermi surface
splits a major shell. Examples:

- Z or N = 6, 14, 28, 50, 82, 126 (12C, 28Si, 56Ni, 100Sn/132Sn, 208Pb, ...).
- Fractional occupations, where an orbital can count as both particle and hole.

16O and 40Ca are not affected.

Affected quantities:

- RPA with `TBCF = true`: `sigma.bin`, `E_corr (density)`, and every consumer of the RPA Sigma (`Reference_TBCF`, exact transitions).
- HF_ScRPA: the final Sigma and the density energy. With `TBCF = true` the iteration is affected too.

## Reproduction

Bundled interaction (`io/NN.bin`, `io/NNN.bin`, hw = 16, Nmax = 3, N2max = 6, N3max = 9); 12C (A = 12, Z = 6),
`HF` then `RPA` with `RPA_Parameters(Name = "RPA-TBCF", TBCF = true)`.

Symmetry check: evaluate each branch body for `(Bra,Ket)` and for `(Ket,Bra)` over the full square, using the
kernels with the RPA amplitudes from `RPA_import_binary`, and take `max |S(ab,cd) - S(cd,ab)|`:

| branch | amplitude passed to `RPA_sigma_hpph` | max asymmetry | max \|Sigma\| |
|---|---|---|---|
| pp | X (current) | 7.3e-2 | 0.48 |
| pp | Y | 3.5e-18 | 0.98 |
| pn | Y (current) | 1.7e-18 | 0.16 |
| pn | X | 4.9e-1 | 0.50 |

For 16O both variants give the same result (asymmetry 9e-19), because the h-p-p-h case never lands in the
computed triangle there.

Effect on 12C `RPA_summary.dat`:

| | E_corr (quasi-boson) | E_corr (density) |
|---|---|---|
| current (pp/nn X) | -30.898290 MeV | -30.712444 MeV |
| proposed fix (all Y) | -30.898290 MeV | **-38.331296 MeV** |

## Proposed fix

Pass `Y` instead of `X` to `RPA_sigma_hpph` at lines 153, 154, 159, 160, 273, 274, 279 and 280, as the pn branch
already does at line 209. Renaming the kernel argument from `X` to a neutral name would avoid the same confusion
later.
