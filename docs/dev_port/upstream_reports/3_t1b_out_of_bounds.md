# T1b reads past the end of Orb when Int.Nmax > Calc.Nmax

**Version:** dev a69ecf0, `src/common/T1b.jl:3-4`

## Summary

`T1b` sizes the kinetic operator from the interaction truncation, `Params.Int.Nmax`, and then reads `Orb[k]` and
`Orb[l]` for every `k, l <= a_max`:

```julia
N_max, hw = Params.Int.Nmax, Params.Calc.hw
a_max = div((N_max + 1)*(N_max + 2),2)
```

`HF_solver`, `BCS_solver`, `HFB_solver` and `RPA_solver` pass the calculation orbitals:
`orbitals_make(Params)` or the imported `orb.bin`, both of length `a_max(Calc.Nmax)`. Whenever the interaction
file has a larger Nmax than the calculation, which is the normal case, the loops read past the end of `Orb`.
Because of `@inbounds`, that is an unchecked out-of-bounds memory read: undefined behaviour, not an error.

Only `V2b_read` passes the `Int = true` orbitals that match this size.

## Current effect

Today's results are not wrong. The leading `a_max(Calc) x a_max(Calc)` block of `T` uses only valid orbitals, and
the callers use only that block (see the comment in `HF_radial_potential`). The rest of `T` is built from garbage
`n, l, j`. Problems:

- Any run with bounds checking aborts. This includes `--check-bounds=yes` and `Pkg.test`, which enables it by default.
- A large enough overrun could segfault.

## Reproduction

```julia
P = Parameters(Interaction_Parameters(NN_File = "x", NNN_File = "x", hw = 16.0, Nmax = 3, N2max = 6, N3max = 9),
               Calculation_Parameters(A = 16, Z = 8, hw = 16.0, Nmax = 2, N2max = 4, N3max = 6))
Orb = NuJuliet.orbitals_make(P)     # length 6
T = NuJuliet.T1b(P, Orb)            # returns a 10 x 10 operator
```

With `julia --check-bounds=yes`:

```
ERROR: BoundsError: attempt to access 6-element Vector{Orb1B} at index [7]
 [2] T1b(Params::Parameters, Orb::Vector{Orb1B})
```

## Proposed fix

Size the operator from the orbitals that are passed in:

```julia
hw = Params.Calc.hw
a_max = length(Orb)
```

With this change `V2b_read` still gets the `Int.Nmax`-sized `T`, because it passes the `Int = true` orbitals. The
solvers get a `Calc.Nmax`-sized `T`, and `T.p[1:a_max,1:a_max]` in `HF_radial_potential` keeps working unchanged.
