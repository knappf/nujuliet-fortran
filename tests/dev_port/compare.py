#!/usr/bin/env python3
"""Compare two NuJuliet dev-layout output trees (io/<Path>/...) file by file.

    python3 -I compare.py REF NEW [--tol 1e-8] [--strict-sign] [--quiet] [--nmax N --n2max M]

REF and NEW are io/ directories (or any matching subdirectories, e.g. io/O16). Every file of REF is compared
with the file at the same relative path in NEW; files only in REF are MISSING, files only in NEW are EXTRA.

Numbers a, b agree when |a - b| <= tol * max(1, |a|, |b|)  (absolute below 1, relative above).

File formats (dev f3f7a33):
  orb.bin                     records of 4 Int64 (a, n, l, j) + 2 Float64 (pO, nO): ints exact, floats with tol.
  *_parsed.bin                V2b_parse / V3b_no2b_parse output: Float32 stream.
  other *.bin                 little-endian Float64 stream (ComplexF64 files are re/im pairs).
  *.dat                       token-wise: numeric tokens with tol (plus one unit in the last printed digit, the
                              files print 6-8 decimals), text tokens exactly.

Gauge freedom (eigenvector signs are an arbitrary LAPACK choice, Julia/OpenBLAS vs Fortran/MKL):
  * per-column sign freedom, applied before the element-wise test:
      C.bin                   O1b format (only elements with l_a == l_b, j_a == j_b, rows a, columns b), p & n
      U.bin + V.bin           O1b format, joint column sign (BCS / HFB, canonical basis)
      TDA_X, RPA_X + RPA_Y    blocks (J = 0..J_max, P = +,-) of N x N, written row-major (eigenvector = column of
                              X, strided in the file); RPA_X & RPA_Y share the sign; RPA/QRPA are ComplexF64
      QTDA_X, QRPA_X + QRPA_Y, pnQTDA_X   blocks written column-major (eigenvector contiguous)
    Block sizes are taken from <Name>/spectra/<Solver>.dat (one line per state with J, P) and validated against
    the E.bin length and the X.bin length; if that fails the file falls back to the |.| test below.
  * 2-body files H2b.bin, sigma.bin (O2b_export): the pn block of every (J,P) is written as an N x N square with
    the packed lower-triangle index, so Julia's upper triangle holds aliased entries while the port writes the
    symmetric matrix: only the lower triangle (Ket <= Bra) of the pn blocks is compared. The block sizes follow
    from (Nmax, N2max), found from the file length (or --nmax/--n2max).
  * phonon-dependent .dat files (name contains 'Phonon': transition densities, currents, operator densities,
    transition analysis) may differ by a global sign of each numeric column (the phonon phase).
  * Anything that is equal only in absolute value after the above (e.g. H2b/H1b in the HF basis when an HF
    orbital came out with the opposite sign, which also changes the sign of X_ph amplitudes) is reported as
    SIGN: it passes by default and fails with --strict-sign.

Exit status: 0 if every file is ok (or SIGN without --strict-sign), 1 otherwise.
"""
import argparse
import os
import re
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore', category=RuntimeWarning)   # NaN/Inf entries (e.g. in *_parsed.bin)

NUM = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$')


# ----------------------------------------------------------------------------------------------------------------
# helpers
def close(a, b, tol):
    """element-wise closeness, returns max scaled error"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size == 0:
        return 0.0
    s = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
    d = np.abs(a - b) / s
    d[~np.isfinite(a) & ~np.isfinite(b) & (np.isnan(a) == np.isnan(b))] = 0.0
    return float(np.nanmax(np.where(np.isnan(d), np.inf, d)))


def ho_orbitals(nmax):
    """(n, l, j) of the LHO orbitals in the order of dev orbitals_make"""
    orb = []
    for N in range(nmax + 1):
        for l in range(N + 1):
            n = (N - l) // 2
            if 2 * n + l == N:
                for j in range(abs(2 * l - 1), 2 * l + 2, 2):
                    orb.append((n, l, j))
    return orb


def o1b_len(nmax):
    orb = ho_orbitals(nmax)
    return 2 * sum(1 for a in orb for b in orb if a[1:] == b[1:])


def find_nmax_o1b(n):
    for nmax in range(0, 40):
        if o1b_len(nmax) == n:
            return nmax
    return None


def o2b_blocks(nmax, n2max):
    """sizes [(N_T1, N_T0)] of the (J,P) blocks in O2b_export order (J = 0..2Nmax+1, P = 1,2)"""
    orb = ho_orbitals(nmax)
    jmax = 2 * nmax + 1
    n1 = np.zeros((2, n2max + 2), dtype=np.int64)
    n0 = np.zeros((2, n2max + 2), dtype=np.int64)
    for ia, (na, la, ja) in enumerate(orb):
        for ib, (nb, lb, jb) in enumerate(orb):
            if 2 * (na + nb) + la + lb <= n2max:
                P = (la + lb) % 2
                for J in range(abs(ja - jb) // 2, (ja + jb) // 2 + 1):
                    if ib <= ia:
                        n1[P, J] += 1
                    n0[P, J] += 1
    out = []
    for J in range(jmax + 1):
        for P in range(2):
            t1 = int(n1[P, J]) if J < n1.shape[1] else 0
            t0 = int(n0[P, J]) if J < n0.shape[1] else 0
            out.append((t1, t0))
    return out


def o2b_len(blocks):
    return 2 * sum(t1 * (t1 + 1) // 2 for t1, _ in blocks) + sum(t0 * t0 for _, t0 in blocks)


def o2b_pn_mask(blocks):
    """boolean mask over the O2b file: False on the aliased upper triangle of the pn blocks"""
    pp = sum(t1 * (t1 + 1) // 2 for t1, _ in blocks)
    pn = sum(t0 * t0 for _, t0 in blocks)
    mask = np.ones(2 * pp + pn, dtype=bool)
    off = pp
    for _, t0 in blocks:
        if t0:
            bra, ket = np.meshgrid(np.arange(t0), np.arange(t0), indexing='ij')   # file order: Bra outer, Ket inner
            mask[off:off + t0 * t0] = (ket <= bra).ravel()
            off += t0 * t0
    return mask


_o2b_cache = {}


def find_o2b_geometry(n, args):
    if args.nmax is not None and args.n2max is not None:
        cands = [(args.nmax, args.n2max)]
    else:
        cands = [(nm, n2) for nm in range(0, 21) for n2 in range(nm, 2 * nm + 1)]
    for nm, n2 in cands:
        key = (nm, n2)
        if key not in _o2b_cache:
            _o2b_cache[key] = o2b_blocks(nm, n2)
        if o2b_len(_o2b_cache[key]) == n:
            return _o2b_cache[key]
    return None


def read_spectra_blocks(spectra_file):
    """[(J, P, count)] in file order from a spectra .dat (lines 'J P E ...')"""
    counts = {}
    order = []
    with open(spectra_file, errors='replace') as f:
        for line in f:
            t = line.split()
            if len(t) >= 3 and t[0].isdigit() and t[1] in '+-' and len(t[1]) == 1:
                key = (int(t[0]), 0 if t[1] == '+' else 1)
                if key not in counts:
                    counts[key] = 0
                    order.append(key)
                counts[key] += 1
    return [counts[k] for k in sorted(order)]


# ----------------------------------------------------------------------------------------------------------------
# column-sign gauge
def column_signs(ref_cols, new_cols):
    """ref_cols/new_cols: lists of 1-D arrays (one per column, possibly concatenated parts) -> +-1 per column"""
    s = []
    for r, n in zip(ref_cols, new_cols):
        d = float(np.dot(r, n))
        s.append(-1.0 if d < 0.0 else 1.0)
    return s


def o1b_columns(vec, nmax):
    """index arrays (into the O1b vector) of each column b, for p and n: list of (isospin, b, idx)"""
    orb = ho_orbitals(nmax)
    amax = len(orb)
    pos = {}
    k = 0
    for t in range(2):
        for a in range(amax):
            for b in range(amax):
                if orb[a][1:] == orb[b][1:]:
                    pos.setdefault((t, b), []).append(k)
                    k += 1
    return [np.array(pos[key]) for key in sorted(pos)]


def block_columns(nblocks, ncomplex, rowmajor, total):
    """index arrays of each eigenvector column for square blocks n x n; ncomplex = 2 for ComplexF64"""
    cols = []
    off = 0
    for n in nblocks:
        for k in range(n):
            if rowmajor:      # element [i,k] at i*n + k
                el = np.arange(n) * n + k
            else:             # element [qp,nu] at nu*n + qp
                el = k * n + np.arange(n)
            if ncomplex == 2:
                el = np.concatenate([2 * el, 2 * el + 1])
            cols.append(off + el)
        off += n * n * ncomplex
    return cols if off == total else None


# ----------------------------------------------------------------------------------------------------------------
class Result:
    def __init__(self, status, msg):
        self.status, self.msg = status, msg   # status in ok, SIGN, DIFF, SIZE


def judge(a, b, tol, label='', mask=None):
    if mask is not None:
        a, b = a[mask], b[mask]
    e = close(a, b, tol)
    if e <= tol:
        return Result('ok', 'max err %.2e%s' % (e, label))
    ea = close(np.abs(a), np.abs(b), tol)
    if ea <= tol:
        return Result('SIGN', 'max err %.2e, |.| %.2e%s' % (e, ea, label))
    return Result('DIFF', 'max err %.2e, |.| %.2e%s' % (e, ea, label))


def eig_partner(name):
    """X/Y and U/V pairs that share the eigenvector sign"""
    for x, y in (('_X.bin', '_Y.bin'), ('U.bin', 'V.bin')):
        if name.endswith(x) and (x != 'U.bin' or name == 'U.bin'):
            return name[:-len(x)] + y
        if name.endswith(y) and (y != 'V.bin' or name == 'V.bin'):
            return name[:-len(y)] + x
    return None


def eig_columns(ref_path, name, n):
    """column index arrays for eigenvector-like files, or None (unknown structure), or [] (not eigenvector-like)"""
    d = os.path.dirname(ref_path)
    if name in ('C.bin', 'U.bin', 'V.bin'):
        nm = find_nmax_o1b(n)
        return None if nm is None else o1b_columns(n, nm)
    m = re.match(r'^(TDA|RPA|QTDA|QRPA|pnQTDA)_(X|Y)\.bin$', name)
    if not m:
        return []
    solver = m.group(1)
    spec = os.path.join(os.path.dirname(d), 'spectra', solver + '.dat')
    efile = os.path.join(d, solver + '_E.bin')
    if not (os.path.isfile(spec) and os.path.isfile(efile)):
        return None
    nb = read_spectra_blocks(spec)
    ne = os.path.getsize(efile) // 8
    if sum(nb) == ne:
        ncomplex = 1
    elif 2 * sum(nb) == ne:
        ncomplex = 2
    else:
        return None
    return block_columns(nb, ncomplex, rowmajor=solver in ('TDA', 'RPA'), total=n)


def cmp_eig(rp, np_, name, a, b, tol, rel_ref, rel_new):
    cols = eig_columns(rp, name, a.size)
    if cols is None:
        return judge(a, b, tol, '  (column structure unknown)')
    if cols == []:
        return judge(a, b, tol)
    # joint sign with the partner file (X/Y, U/V) when present in both trees
    partner = eig_partner(name)
    pa = pb = None
    if partner:
        pr, pn = os.path.join(os.path.dirname(rp), partner), os.path.join(os.path.dirname(np_), partner)
        if os.path.isfile(pr) and os.path.isfile(pn) and os.path.getsize(pr) == os.path.getsize(pn) == a.size * 8:
            pa, pb = np.fromfile(pr, '<f8'), np.fromfile(pn, '<f8')
    # the sign is taken from the X (U) part first, the partner only breaks ties of a zero column
    first_is_main = partner is None or name.endswith('_X.bin') or name == 'U.bin'
    rc, nc = [], []
    for c in cols:
        if pa is None:
            rc.append(a[c]); nc.append(b[c])
        elif first_is_main:
            rc.append(np.concatenate([a[c], 1e-3 * pa[c]])); nc.append(np.concatenate([b[c], 1e-3 * pb[c]]))
        else:
            rc.append(np.concatenate([pa[c], 1e-3 * a[c]])); nc.append(np.concatenate([pb[c], 1e-3 * b[c]]))
    s = column_signs(rc, nc)
    bb = b.copy()
    for c, sg in zip(cols, s):
        bb[c] *= sg
    nflip = sum(1 for x in s if x < 0)
    return judge(a, bb, tol, '  (%d of %d columns sign-flipped)' % (nflip, len(cols)) if nflip else '')


def cmp_bin(rp, np_, name, args):
    if name == 'orb.bin':
        dt = np.dtype([('i', '<i8', 4), ('o', '<f8', 2)])
        if os.path.getsize(rp) != os.path.getsize(np_) or os.path.getsize(rp) % dt.itemsize:
            return Result('SIZE', 'size %d vs %d' % (os.path.getsize(rp), os.path.getsize(np_)))
        a, b = np.fromfile(rp, dt), np.fromfile(np_, dt)
        if not np.array_equal(a['i'], b['i']):
            bad = np.argwhere(a['i'] != b['i'])[0]
            return Result('DIFF', 'integer (a,n,l,j) records differ, first at orbital %d' % (bad[0] + 1))
        return judge(a['o'].ravel(), b['o'].ravel(), args.tol)
    dt = '<f4' if name.endswith('_parsed.bin') else '<f8'
    a, b = np.fromfile(rp, dt), np.fromfile(np_, dt)
    if a.size != b.size:
        return Result('SIZE', '%d vs %d elements' % (a.size, b.size))
    if name in ('H2b.bin', 'sigma.bin'):
        geo = find_o2b_geometry(a.size, args)
        if geo is None:
            return judge(a, b, args.tol, '  (O2b geometry unknown, full compare)')
        return judge(a, b, args.tol, '  (pn: lower triangle)', mask=o2b_pn_mask(geo))
    return cmp_eig(rp, np_, name, a, b, args.tol, None, None)


def dat_tokens(path):
    out = []
    with open(path, errors='replace') as f:
        for li, line in enumerate(f):
            for t in line.split():
                out.append((li, t))
    return out


def last_digit_unit(tok):
    m = re.match(r'^[-+]?(\d*)\.?(\d*)(?:[eEdD]([-+]?\d+))?$', tok)
    if not m:
        return 0.0
    dec = len(m.group(2)) if '.' in tok else 0
    ex = int(m.group(3)) if m.group(3) else 0
    return 10.0 ** (ex - dec)


def cmp_dat(rp, np_, name, args):
    ta, tb = dat_tokens(rp), dat_tokens(np_)
    if len(ta) != len(tb):
        return Result('SIZE', '%d vs %d tokens' % (len(ta), len(tb)))
    sign_free = 'Phonon' in name
    va, vb, ul, cols, text = [], [], [], [], 0
    first_text = None
    colpos = {}
    lastline, col = -1, 0
    for (la, x), (_, y) in zip(ta, tb):
        if la != lastline:
            lastline, col = la, 0
        if NUM.match(x) and NUM.match(y):
            va.append(float(x.replace('d', 'e').replace('D', 'e')))
            vb.append(float(y.replace('d', 'e').replace('D', 'e')))
            ul.append(max(last_digit_unit(x), last_digit_unit(y)))
            cols.append(col)
        elif x != y:
            text += 1
            if first_text is None:
                first_text = 'line %d: %r vs %r' % (la + 1, x, y)
        col += 1
    if text:
        return Result('DIFF', '%d text tokens differ (first %s)' % (text, first_text))
    va, vb, ul, cols = np.array(va), np.array(vb), np.array(ul), np.array(cols, dtype=np.int64)
    if va.size == 0:
        return Result('ok', 'text only')

    def err(a, b):
        s = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
        d = np.maximum(np.abs(a - b) - 1.01 * ul, 0.0) / s
        return float(d.max())

    e = err(va, vb)
    if e <= args.tol:
        return Result('ok', 'max err %.2e' % e)
    if sign_free:
        bb = vb.copy()
        for c in np.unique(cols):
            m = cols == c
            if np.dot(va[m], vb[m]) < 0.0:
                bb[m] = -bb[m]
        e2 = err(va, bb)
        if e2 <= args.tol:
            return Result('ok', 'max err %.2e  (phonon phase: column signs)' % e2)
    ea = err(np.abs(va), np.abs(vb))
    if ea <= args.tol:
        return Result('SIGN', 'max err %.2e, |.| %.2e' % (e, ea))
    i = int(np.argmax(np.abs(va - vb) / np.maximum(1.0, np.maximum(np.abs(va), np.abs(vb)))))
    return Result('DIFF', 'max err %.2e, |.| %.2e (e.g. %r vs %r)' % (e, ea, float(va[i]), float(vb[i])))


# ----------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description='Compare two NuJuliet dev-layout io/ trees.')
    ap.add_argument('ref')
    ap.add_argument('new')
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--strict-sign', action='store_true', help='SIGN (equal only in |.|) counts as a difference')
    ap.add_argument('--quiet', action='store_true', help='print only files that are not ok')
    ap.add_argument('--nmax', type=int, default=None, help='Calc Nmax for the O2b geometry (default: from size)')
    ap.add_argument('--n2max', type=int, default=None, help='Calc N2max for the O2b geometry')
    args = ap.parse_args()

    for d in (args.ref, args.new):
        if not os.path.isdir(d):
            print('not a directory: %s' % d)
            sys.exit(2)

    ref_files = set()
    for root, _, files in os.walk(args.ref):
        for f in files:
            ref_files.add(os.path.relpath(os.path.join(root, f), args.ref))
    new_files = set()
    for root, _, files in os.walk(args.new):
        for f in files:
            new_files.add(os.path.relpath(os.path.join(root, f), args.new))

    count = {'ok': 0, 'SIGN': 0, 'DIFF': 0, 'SIZE': 0, 'MISSING': 0, 'EXTRA': 0}
    for rel in sorted(ref_files):
        rp, np_ = os.path.join(args.ref, rel), os.path.join(args.new, rel)
        if rel not in new_files:
            count['MISSING'] += 1
            print('MISSING  %s' % rel)
            continue
        name = os.path.basename(rel)
        try:
            if name.endswith('.bin'):
                r = cmp_bin(rp, np_, name, args)
            else:
                r = cmp_dat(rp, np_, name, args)
        except Exception as ex:          # unreadable / malformed file
            r = Result('DIFF', 'error: %s' % ex)
        count[r.status] += 1
        if not (args.quiet and r.status == 'ok'):
            print('%-8s %-90s %s' % (r.status, rel, r.msg))
    for rel in sorted(new_files - ref_files):
        count['EXTRA'] += 1
        print('EXTRA    %s' % rel)

    print('\n%d files: %d ok, %d SIGN (equal up to signs), %d DIFF, %d SIZE, %d MISSING, %d EXTRA  (tol %.1e)'
          % (len(ref_files), count['ok'], count['SIGN'], count['DIFF'], count['SIZE'], count['MISSING'],
             count['EXTRA'], args.tol))
    bad = count['DIFF'] + count['SIZE'] + count['MISSING'] + count['EXTRA']
    if args.strict_sign:
        bad += count['SIGN']
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
