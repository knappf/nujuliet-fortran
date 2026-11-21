# OBDM_natural_orbitals: neutron natural orbitals get their l, j from the proton eigenvectors (minor)

**Version:** dev a69ecf0, `src/common/orb.jl:623-631`

## Summary

`l` and `j` of each natural orbital come from the weighted average over the **proton** eigenvectors only:

```julia
lSum += abs2(pC[b,a]) * l_b
jSum += abs2(pC[b,a]) * j_b
```

The same `(l, j)` is then used for the neutron natural orbital `nC[:,a]`, which is never checked.

## Severity

Low. With the current code this cannot go wrong in practice. Rho is cleaned to be block-diagonal in `(l, j)`, and
after the Hungarian reordering column `a` of both `pC` and `nC` lies in the `(l, j)` block of `Orb[a]`, so the
proton and neutron values agree.

The code depends on that silently, though. If the proton and neutron natural orbitals ever ended up in different
columns, the neutron orbitals would be mislabelled and nothing would report it.

## Proposed fix

Either:

- compute the neutron `(l, j)` from `nC` the same way and raise an error if it differs from the proton value
  (one `Orb1B` holds a single `(l, j)`), or
- simply use `Orb[a].l` and `Orb[a].j`, which is what the reordering already guarantees.

```julia
nlSum += abs2(nC[b,a]) * l_b
njSum += abs2(nC[b,a]) * j_b
...
if round(Int64,nlSum) != l_a || round(Int64,njSum) != j_a
    error("OBDM_natural_orbitals: proton & neutron natural orbital " * string(a) * " have different l/j")
end
```
