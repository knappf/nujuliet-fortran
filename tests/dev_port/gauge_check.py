"""Element-wise comparison of O1b/O2b files of a rotated basis (HF or HF-MBPT natural orbitals) with the column-sign
gauge of C.bin applied: O1b element (a,b) gets s_a s_b, O2b element <ab|O|cd> gets s_a s_b s_c s_d (proton signs for
protons, neutron signs for neutrons; pn blocks: lower triangle only, see compare.py).

    python3 -I gauge_check.py REF_BIN NEW_BIN NMAX N2MAX [files...]      (default files H1b.bin H2b.bin sigma.bin rho.bin)

REF_BIN / NEW_BIN are bin/ directories with C.bin, e.g. work/julia_t8/bundled/io/O16_MBPT/HF-MBPT/bin. Prints the
max scaled error |ref - s new| / max(1,|ref|) per file (and per pp/pn/nn block)."""
import sys, os
import numpy as np

ref, new, nmax, n2max = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
files = sys.argv[5:] or ['H1b.bin', 'H2b.bin', 'sigma.bin', 'rho.bin']

orb = []
for N in range(nmax + 1):
    for l in range(N + 1):
        n = (N - l) // 2
        if 2 * n + l == N:
            for j in range(abs(2 * l - 1), 2 * l + 2, 2):
                orb.append((n, l, j))
amax = len(orb)
lj_pairs = [(a, b) for a in range(amax) for b in range(amax) if orb[a][1:] == orb[b][1:]]


def read_o1b(path):
    v = np.fromfile(path, dtype='<f8')
    h = len(v) // 2
    P = np.zeros((amax, amax)); Nn = np.zeros((amax, amax))
    for k, (a, b) in enumerate(lj_pairs):
        P[a, b] = v[k]; Nn[a, b] = v[h + k]
    return P, Nn


Cr = read_o1b(os.path.join(ref, 'C.bin'))
Cn = read_o1b(os.path.join(new, 'C.bin'))
signs = []
for t in range(2):
    s = np.ones(amax)
    for b in range(amax):
        k = np.argmax(np.abs(Cr[t][:, b]))
        s[b] = np.sign(Cr[t][k, b]) * np.sign(Cn[t][k, b])
    signs.append(s)
print('flipped columns p:', np.where(signs[0] < 0)[0] + 1, ' n:', np.where(signs[1] < 0)[0] + 1)

jmax = 2 * nmax + 1
T1 = {}; T0 = {}
for J in range(n2max + 2):
    for P in range(2):
        T1[(J, P)] = []; T0[(J, P)] = []
for a in range(amax):
    for b in range(amax):
        na, la, ja = orb[a]; nb, lb, jb = orb[b]
        if 2 * (na + nb) + la + lb <= n2max:
            P = (la + lb) % 2
            for J in range(abs(ja - jb) // 2, (ja + jb) // 2 + 1):
                if b <= a:
                    T1[(J, P)].append((a, b))
                T0[(J, P)].append((a, b))


def scale(v):
    return np.maximum(1.0, np.abs(v))


for f in files:
    pr, pn_ = os.path.join(ref, f), os.path.join(new, f)
    if not os.path.exists(pr):
        continue
    if f in ('H1b.bin', 'rho.bin', 'C.bin'):
        R = read_o1b(pr); Nw = read_o1b(pn_)
        err = 0.0
        for t in range(2):
            g = Nw[t] * np.outer(signs[t], signs[t])
            err = max(err, np.max(np.abs(R[t] - g) / np.maximum(1.0, np.abs(R[t]))))
        print(f'{f:10s} gauge max err {err:.2e}')
        continue
    vr = np.fromfile(pr, dtype='<f8'); vn = np.fromfile(pn_, dtype='<f8')
    assert len(vr) == len(vn)
    off = 0; err = {'pp': 0.0, 'pn': 0.0, 'nn': 0.0}
    for kind in ('pp', 'pn', 'nn'):
        for J in range(jmax + 1):
            for P in range(2):
                if kind == 'pn':
                    L = T0.get((J, P), []); N = len(L)
                    for bra in range(N):
                        for ket in range(N):
                            if ket <= bra:
                                a, b = L[bra]; c, d = L[ket]
                                s = signs[0][a] * signs[1][b] * signs[0][c] * signs[1][d]
                                x, y = vr[off], vn[off] * s
                                err['pn'] = max(err['pn'], abs(x - y) / max(1.0, abs(x)))
                            off += 1
                else:
                    t = 0 if kind == 'pp' else 1
                    L = T1.get((J, P), []); N = len(L)
                    for bra in range(N):
                        for ket in range(bra + 1):
                            a, b = L[bra]; c, d = L[ket]
                            s = signs[t][a] * signs[t][b] * signs[t][c] * signs[t][d]
                            x, y = vr[off], vn[off] * s
                            err[kind] = max(err[kind], abs(x - y) / max(1.0, abs(x)))
                            off += 1
    assert off == len(vr), (off, len(vr))
    print(f'{f:10s} gauge max err pp {err["pp"]:.2e} pn {err["pn"]:.2e} nn {err["nn"]:.2e}  (max |ref| {np.max(np.abs(vr)):.3g})')
