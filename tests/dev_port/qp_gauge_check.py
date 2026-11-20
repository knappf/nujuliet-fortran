"""Element-wise comparison of quasiparticle-basis files (BCS/HFB qpH2b.bin, QTDA_X.bin, QRPA_X.bin + QRPA_Y.bin,
pnQTDA_X.bin) with the sign gauge of the quasiparticle states applied.

The quasiparticle state k (per isospin) has a free sign s_k: it follows the canonical basis C (BCS: s_k = sign of C
column k) and the HFB eigenvectors (U, V columns). s_k is found from the reference & new C.bin, U.bin, V.bin of the
reference solver:  s_k = sign( sum_i sC_i (U_ref U_new + V_ref V_new)[i,k] ),  sC = column signs of C.

  * qpH2b.bin: element <ab|H|cd> gets s_a s_b s_c s_d (isospins per block: pp, nn; 40 pn & 22 pn2002: a,b protons,
    c,d neutrons; 22 pn1111: a,c protons, b,d neutrons).
  * X files (column-major blocks J = 0..2Nmax+1, P = +,-): row (a,b) gets s_a s_b, then each column (phonon) keeps a
    free global sign (QRPA_X and QRPA_Y share it; QRPA files are ComplexF64).

    python3 -I qp_gauge_check.py REF_IO NEW_IO NMAX N2MAX REFERENCE      (REFERENCE = BCS or HFB)

REF_IO / NEW_IO are io/<Path> directories, e.g. work/julia_t8/bundled/io/O18_HFB. Prints the max scaled error
|ref - gauged new| / max(1,|ref|) per file."""
import sys, os
import numpy as np

ref, new, nmax, n2max, R = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]

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


def rd(base, f):
    return read_o1b(os.path.join(base, R, 'bin', f))


Cr, Cn = rd(ref, 'C.bin'), rd(new, 'C.bin')
Ur, Un = rd(ref, 'U.bin'), rd(new, 'U.bin')
Vr, Vn = rd(ref, 'V.bin'), rd(new, 'V.bin')
sq = []
for t in range(2):
    sC = np.ones(amax)
    for b in range(amax):
        k = np.argmax(np.abs(Cr[t][:, b]))
        sC[b] = np.sign(Cr[t][k, b]) * np.sign(Cn[t][k, b])
    s = np.ones(amax)
    for k in range(amax):
        ov = np.sum(sC * (Ur[t][:, k] * Un[t][:, k] + Vr[t][:, k] * Vn[t][:, k]))
        s[k] = 1.0 if ov >= 0 else -1.0
    sq.append(s)
print('flipped qp states p:', np.where(sq[0] < 0)[0] + 1, ' n:', np.where(sq[1] < 0)[0] + 1)
sp, sn = sq

# --- qpH2b.bin -------------------------------------------------------------------------------------------------
jmax = 2 * nmax + 1
je = min(jmax, n2max + 1)
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
H_orb = [2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2]
H_sym = [1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1]
# isospin signs of the 4 indices per block (0 = proton, 1 = neutron)
H_iso = [(0, 0, 0, 0), (0, 0, 1, 1), (1, 1, 1, 1), None, None, None, None, (0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1),
         (1, 1, 1, 1)]


def check_h2b(pr, pn_):
    vr = np.fromfile(pr, dtype='<f8'); vn = np.fromfile(pn_, dtype='<f8')
    assert len(vr) == len(vn)
    off = 0; err = 0.0
    for H in range(11):
        for J in range(je + 1):
            for P in range(2):
                L = (T1 if H_orb[H] == 2 else T0)[(J, P)]; N = len(L)
                for bra in range(N):
                    for ket in range(bra + 1 if H_sym[H] else N):
                        if H_iso[H] is not None:
                            (a, b), (c, d) = L[bra], L[ket]
                            i1, i2, i3, i4 = H_iso[H]
                            s = sq[i1][a] * sq[i2][b] * sq[i3][c] * sq[i4][d]
                        else:
                            s = 1.0
                        err = max(err, abs(vr[off] - s * vn[off]) / max(1.0, abs(vr[off])))
                        off += 1
    assert off == len(vr), (off, len(vr))
    return err


# --- X files ---------------------------------------------------------------------------------------------------
def pairs_2qp(J, P):
    L = []
    for t in (0, 1):
        for a in range(amax):
            for b in range(a + 1):
                ja, jb = orb[a][2], orb[b][2]
                if abs(ja - jb) // 2 <= J <= (ja + jb) // 2 and (a != b or J % 2 == 0) and (orb[a][1] + orb[b][1]) % 2 == P:
                    L.append((a, b, t))
    return [(sq[t][a] * sq[t][b]) for a, b, t in L]


def pairs_pn(J, P):
    L = []
    for a in range(amax):
        for b in range(amax):
            ja, jb = orb[a][2], orb[b][2]
            if abs(ja - jb) // 2 <= J <= (ja + jb) // 2 and (orb[a][1] + orb[b][1]) % 2 == P:
                L.append(sp[a] * sn[b])
    return L


def check_x(files, rowsign, cplx):
    vr = [np.fromfile(os.path.join(ref, f), dtype='<f8') for f in files]
    vn = [np.fromfile(os.path.join(new, f), dtype='<f8') for f in files]
    if cplx:
        vr = [v[0::2] + 1j * v[1::2] for v in vr]; vn = [v[0::2] + 1j * v[1::2] for v in vn]
    off = 0; err = 0.0
    for J in range(jmax + 1):
        for P in range(2):
            rs = np.array(rowsign(J, P)); N = len(rs)
            if N == 0:
                continue
            Br = [v[off:off + N * N].reshape(N, N, order='F') for v in vr]
            Bn = [v[off:off + N * N].reshape(N, N, order='F') * rs[:, None] for v in vn]
            for k in range(N):
                i = np.argmax(np.abs(Br[0][:, k]))
                c = np.sign(np.real(Br[0][i, k] * np.conj(Bn[0][i, k]))) or 1.0
                for x, y in zip(Br, Bn):
                    err = max(err, np.max(np.abs(x[:, k] - c * y[:, k]) / np.maximum(1.0, np.abs(x[:, k]))))
            off += N * N
    assert off == len(vr[0]), (off, len(vr[0]))
    return err


h2b = os.path.join(R, 'bin', 'qpH2b.bin')
print(f'{h2b:28s} gauge max err {check_h2b(os.path.join(ref, h2b), os.path.join(new, h2b)):.2e}')
for files, rs, cx in ((['QTDA/bin/QTDA_X.bin'], pairs_2qp, False),
                      (['QRPA/bin/QRPA_X.bin', 'QRPA/bin/QRPA_Y.bin'], pairs_2qp, True),
                      (['pnQTDA/bin/pnQTDA_X.bin'], pairs_pn, False)):
    if all(os.path.exists(os.path.join(ref, f)) and os.path.exists(os.path.join(new, f)) for f in files):
        print(f'{" + ".join(files):28s} gauge max err {check_x(files, rs, cx):.2e}')
