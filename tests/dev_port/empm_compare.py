#!/usr/bin/env python3
"""EMPM / NCSM consumer checks (stage 5): legacy baseline vs dev-form chain, and EMPM_TDA vs the dev TDA.

  python3 -I empm_compare.py BASE_IO DEV_IO [--tol 1e-10] [--quiet]
      BASE_IO = work/empm_baseline/io (frozen stage-4 binary: legacy HF, io/<Path>/Bin/*, EMPM/{Spectra,Transitions}, ...)
      DEV_IO  = work/empm_dev_t<N>/io (dev HF io/<Path>/HF/bin/*, EMPM/{bin,spectra,transitions}, ...)
      Every legacy file of io/<Path>/{Bin/EMPM_*, EMPM/**, NCSM/**} and the HF reference files is mapped to its dev name
      (Bin/EMPM_* -> EMPM/bin/EMPM_*, EMPM/Spectra -> EMPM/spectra, EMPM/Transitions -> EMPM/transitions,
      EMPM/TDA_Summary.dat -> EMPM/TDA_summary.dat, Bin/{C_HF,h_HF,V2B_HF,Rho_HF}.bin -> HF/bin/{C,H1b,H2b,rho}.bin,
      HF/HF_Summary.dat -> HF/HF_summary.dat) and compared; the max deviation is printed per file:
        d     = max |a - b| / max(1, |a|, |b|)     (signed)
        d_abs = the same for |a| and |b|           (eigenvector / orbital sign gauges removed)
      Text files: every number is compared up to its own printed unit (one unit of the last printed digit, i.e. a
      rounding flip; integers exactly), d / d_abs are the deviations beyond that unit, the measure used for pass/fail.
      Status: ok (d <= tol), SIGN (only d_abs <= tol; ONLY for files that carry a sign gauge, see below), DIFF, SIZE,
      MISSING (no dev file), TEXT (text tokens differ; numbers are still compared), SKIP (TDA_summary.dat: the dev text
      replaced the HF_RPA-style header; only its numbers compared), GAUGE (see below).
      Sign-gauge files (a sign-only change is SIGN, not a failure): HF/bin/{C,H1b,H2b}.bin (HF orbital signs),
      EMPM/bin/EMPM_TDA_X.bin & EMPM_F.bin (phonon / orbital signs), EMPM/TDA_Structure.dat (X components); the NCSM
      transition tables NCSM_*_E<L>.dat: only the amplitude columns M_p, M_n, with one sign per row (state phase), the
      other columns signed. Every other file (energies *_E.bin, *_E_<lab>.bin, Ncm, Spectrum.dat, TDA.dat, B values,
      NCSM energy tables, ...) is compared signed: a sign-only change there is DIFF.
      Binary state files EMPM_S<n>_<JP>.bin are parsed (int32 header & index arrays compared exactly, E signed, C and X
      per column up to a sign, by projectors in degenerate clusters, R by singular values); EMPM_2PH_<JP>_C.bin,
      EMPM_FULL*_Z*.bin likewise (energies from the matching _E file of both runs, compared signed). GAUGE: the
      independent rows chosen by the pivoted Cholesky differ (Ind of Basis_<JP>.dat, Kc of the state files), so C on
      those rows is not comparable; the energies (and X) are compared and must agree, else DIFF.
      Other .bin files: Float64 streams.
      Known, documented differences (consumers.md 3.2.3): the spurious E(1) of EMPM_TDA (legacy 0, dev <spur|A|spur>) in
      EMPM/bin/EMPM_TDA_E.bin, EMPM/spectra/TDA.dat, EMPM/transitions/TDA_E1.dat, EMPM/TDA_Structure.dat; they are
      reported separately with the spurious value excluded (status 'ok*'; 'SIGN*' only for TDA_Structure.dat).

  python3 -I empm_compare.py --tda DEV_IO [--tol 1e-10]
      For every io/<Path> with both EMPM/ and TDA/: EMPM/spectra/TDA.dat vs TDA/spectra/TDA.dat (per J^P channel),
      EMPM/bin/EMPM_TDA_E.bin vs TDA/bin/TDA_E.bin, and EMPM/transitions/TDA_E<L>.dat vs TDA/transitions/E<L>/TDA_E<L>.dat
      (L = 0..3; E, B_ph, B_is, B_iv beyond the printed unit 1e-8); a missing file is MISSING (failure).

Exit status 1 on DIFF/SIZE/MISSING (legacy mode), on a deviation above tol or a MISSING file (TDA mode).
"""
import os
import re
import sys

import warnings

import numpy as np

warnings.filterwarnings('ignore', category=UserWarning)


def rel(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0:
        return 0.0
    if a.shape != b.shape:
        return float('inf')
    both_nan = np.isnan(a) & np.isnan(b)
    a = np.where(both_nan, 0.0, a)
    b = np.where(both_nan, 0.0, b)
    with np.errstate(invalid='ignore'):
        r = np.abs(a - b) / np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
    return float('inf') if np.any(np.isnan(r)) else float(np.max(r))


def clusters(E, tol=1e-6):
    """index groups of (numerically) degenerate eigenvalues, E in the stored order"""
    order = np.argsort(E, kind='stable')
    groups, cur = [], [order[0]] if len(order) else []
    for i, j in zip(order[:-1], order[1:]):
        if abs(E[j] - E[i]) < tol * max(1.0, abs(E[i])):
            cur.append(j)
        else:
            groups.append(cur)
            cur = [j]
    if cur:
        groups.append(cur)
    return groups


def vec_dev(E, A, B):
    """max deviation of eigenvector columns A vs B (same E), invariant under signs & rotations in degenerate clusters:
    singletons sign-aligned, clusters by the projector A_g A_g^T"""
    if A.shape != B.shape:
        return float('inf')
    d = 0.0
    for g in clusters(E):
        if len(g) == 1:
            a, b = A[:, g[0]], B[:, g[0]]
            sgn = 1.0 if np.dot(a, b) >= 0 else -1.0
            d = max(d, rel(a, sgn * b))
        else:
            Ag, Bg = A[:, g], B[:, g]
            d = max(d, rel(Ag @ Ag.T, Bg @ Bg.T))
    return d


NUM = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$')
NUMP = re.compile(r'^[+-]?\d*(?:\.(\d*))?(?:[eEdD]([+-]?\d+))?$')


def unit_of(t):
    """one unit of the last printed digit of the number token t (0 for an integer: compared exactly)"""
    m = NUMP.match(t)
    dec, ex = m.group(1), m.group(2)
    if dec is None and ex is None:
        return 0.0
    return 10.0 ** ((int(ex) if ex else 0) - len(dec or ''))


def num_tokens(tl):
    """numbers, their printed units & the text tokens of a token list"""
    nums, units, text = [], [], []
    for t in tl:
        if NUM.match(t):
            nums.append(float(t.replace('d', 'e').replace('D', 'e')))
            units.append(unit_of(t))
        else:
            text.append(t)
    return nums, units, text


def file_lines(path):
    """lines of a text file; TDA_Structure.dat: the component lines of every phonon sorted by (|X| descending, label),
    since the order among components with equal printed |X| (e.g. p/n partners) is arbitrary"""
    with open(path, errors='replace') as f:
        lines = f.read().splitlines()
    if os.path.basename(path) != 'TDA_Structure.dat':
        return lines
    out, block = [], []

    def flush():
        block.sort(key=lambda l: (-abs(float(l.split()[-1])), l.split()[:-1]))
        out.extend(block)
        block.clear()
    for line in lines:
        t = line.split()
        if len(t) == 5 and t[0] in ('p', 'n') and NUM.match(t[-1]):
            block.append(line)
        else:
            flush()
            out.append(line)
    flush()
    return out


def tokens(path, drop_seconds=False):
    """numbers, their printed units & the text tokens ... drop_seconds: the last column of the tables headed by a
    'seconds' column (timings)"""
    nums, units, text = [], [], []
    ncol = 0
    for line in file_lines(path):
        tl = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').split()
        if drop_seconds and tl and tl[-1] == 'seconds':
            ncol = len(tl)
        elif drop_seconds and ncol and len(tl) == ncol and NUM.match(tl[-1]):
            tl = tl[:-1]
        elif not tl:
            ncol = 0
        n, u, t = num_tokens(tl)
        nums += n
        units += u
        text += t
    return np.array(nums), np.array(units), text


def dev_units(a, b, u):
    """max (|a - b| - u)_+ / max(1, |a|, |b|): the deviation beyond the printed unit u (a rounding flip of the last
    printed digit is allowed), the measure used for pass/fail of text files"""
    a, b, u = (np.asarray(x, dtype=float) for x in (a, b, u))
    if a.size == 0:
        return 0.0
    if a.shape != b.shape:
        return float('inf')
    r = np.maximum(0.0, np.abs(a - b) - u * (1.0 + 1e-6)) / np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
    return float(np.max(r))


# ---------------------------------------------------------------------------------------------------------------------
def map_legacy(rel_path):
    """legacy path below io/ -> dev path, or None if not compared"""
    p = rel_path.split('/')
    if len(p) < 3:
        return None
    path, rest = p[0], p[1:]
    if rest[0] == 'Bin':
        name = rest[1]
        if name.startswith('EMPM_'):
            return '/'.join([path, 'EMPM', 'bin', name])
        m = {'C_HF.bin': 'C.bin', 'h_HF.bin': 'H1b.bin', 'V2B_HF.bin': 'H2b.bin', 'Rho_HF.bin': 'rho.bin'}
        if name in m:
            return '/'.join([path, 'HF', 'bin', m[name]])
        return None
    if rest[0] == 'EMPM':
        r = rest[1:]
        if r[0] == 'Spectra':
            r[0] = 'spectra'
        elif r[0] == 'Transitions':
            r[0] = 'transitions'
        elif r == ['TDA_Summary.dat']:
            r = ['TDA_summary.dat']
        return '/'.join([path, 'EMPM'] + r)
    if rest[0] == 'NCSM':
        return '/'.join([path] + rest)
    if rest == ['HF', 'HF_Summary.dat']:
        return '/'.join([path, 'HF', 'HF_summary.dat'])
    return None


def read_states(path):
    with open(path, 'rb') as f:
        buf = f.read()
    o = 0

    def ints(n):
        nonlocal o
        v = np.frombuffer(buf, dtype='<i4', count=n, offset=o)
        o += 4 * n
        return v

    def dbl(n):
        nonlocal o
        v = np.frombuffer(buf, dtype='<f8', count=n, offset=o)
        o += 8 * n
        return v

    mk = ints(1)[0]
    if mk not in (-3, -2):
        return None
    N_full, n, N_ind = ints(3)
    N_phys = ints(1)[0] if mk == -3 else n
    d = {'hdr': np.array([mk, N_full, n, N_ind, N_phys]), 'Al': ints(N_full), 'La': ints(N_full), 'E': dbl(n),
         'Kc': ints(N_ind)}
    d['C'] = dbl(N_ind * n).reshape((n, N_ind)).T
    d['X'] = dbl(N_full * n).reshape((n, N_full)).T
    d['R'] = dbl(N_phys * (n - N_phys)) if mk == -3 else np.zeros(0)
    return d


def cmp_states(a_path, b_path):
    a, b = read_states(a_path), read_states(b_path)
    if a is None or b is None:
        return None
    for k in ('hdr', 'Al', 'La'):
        if a[k].shape != b[k].shape or np.any(a[k] != b[k]):
            return ('DIFF', float('inf'), float('inf'), 'integer part %s differs' % k)
    d = rel(a['E'], b['E'])
    d = max(d, vec_dev(a['E'], a['X'], b['X']))
    if np.array_equal(a['Kc'], b['Kc']):
        d = max(d, vec_dev(a['E'], a['C'], b['C']))
        note = 'state file: E, X & C (degenerate clusters by projectors), R by singular values'
        st = None
    else:
        note = ('state file: independent rows Kc differ in %d of %d (pivot choice among dependent rows, a gauge): E, X '
                'compared, C not comparable; R by singular values' % (int(np.sum(a['Kc'] != b['Kc'])), a['Kc'].size))
        st = 'GAUGE'
    if a['R'].size:
        n, Np = a['hdr'][2], a['hdr'][4]
        sa = np.linalg.svd(a['R'].reshape((n - Np, Np)).T, compute_uv=False)
        sb = np.linalg.svd(b['R'].reshape((n - Np, Np)).T, compute_uv=False)
        d = max(d, rel(sa, sb))
    if st == 'GAUGE':
        st = 'GAUGE' if d <= TOL[0] else 'DIFF'
    return (st, d, d, note)


def basis_ind(path):
    """Ind column (4th) of a 2Phonon/Basis_<JP>.dat file"""
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) == 8 and t[0].isdigit():
                rows.append(int(t[3]))
    return np.array(rows)


def cmp_eigfile(rp, dp, relp, base_dir, dev_dir):
    """eigenvector files with their energy file: EMPM_2PH_<JP>_C.bin (+ _E.bin), EMPM_FULL_<JP>_Z.bin (+ _E.bin),
    EMPM_FULL<n>_<JP>_Z_<lab>.bin (5 int32 Nb, then Z; + _E_<lab>.bin) ... None if not such a file"""
    name = os.path.basename(rp)
    m2 = re.match(r'EMPM_2PH_(\d+[+-])_C\.bin$', name)
    mf = re.match(r'EMPM_FULL_(\d+[+-])_Z\.bin$', name)
    mn = re.match(r'(EMPM_FULL\d_\d+[+-]_)Z_(\w+)\.bin$', name)
    if not (m2 or mf or mn):
        return None
    ef = name.replace('_C.bin', '_E.bin') if m2 else (name.replace('_Z.bin', '_E.bin') if mf else
                                                      mn.group(1) + 'E_' + mn.group(2) + '.bin')
    E = np.fromfile(os.path.join(os.path.dirname(rp), ef), dtype='<f8')
    fe_dev = os.path.join(os.path.dirname(dp), ef)
    if not os.path.exists(fe_dev):
        return ('MISSING', float('nan'), float('nan'), 'energy file %s of the dev run missing' % ef)
    dE = rel(E, np.fromfile(fe_dev, dtype='<f8'))    # energies (signed), inf if the sizes differ
    if mn:
        def rd(p):
            buf = open(p, 'rb').read()
            Nb = np.frombuffer(buf, '<i4', 5)
            return Nb, np.frombuffer(buf, '<f8', offset=20)
        Nba, a = rd(rp)
        Nbb, b = rd(dp)
        if not np.array_equal(Nba, Nbb):
            return ('DIFF', float('inf'), float('inf'), 'Nb differs')
    else:
        a, b = np.fromfile(rp, dtype='<f8'), np.fromfile(dp, dtype='<f8')
    n = E.size
    if n == 0 and a.size == b.size:
        d = rel(a, b)
        d = max(d, dE)
        return (status(d, d, TOL[0]), d, d, 'no states')
    if n == 0 or a.size != b.size or a.size % n:
        return ('SIZE', float('inf'), float('inf'), 'sizes')
    A, B = a.reshape((n, a.size // n)).T, b.reshape((n, b.size // n)).T
    if m2:
        # C on the independent rows ... comparable only if the Ind flags of Basis_<JP>.dat agree ...
        def basis_of(p):
            up = os.path.dirname(os.path.dirname(p))
            for c in (os.path.join(up, '2Phonon'), os.path.join(up, 'EMPM', '2Phonon')):
                if os.path.isdir(c):
                    return os.path.join(c, 'Basis_%s.dat' % m2.group(1))
            return ''
        base, dev = basis_of(rp), basis_of(dp)
        ia, ib = basis_ind(base), basis_ind(dev)
        if not np.array_equal(ia, ib):
            return (('GAUGE' if dE <= TOL[0] else 'DIFF'), dE, dE,
                    'independent rows differ in %d of %d basis states (pivot choice among dependent rows): C not '
                    'comparable; d = energies of %s' % (int(np.sum(ia != ib)), ia.size, ef))
    d = max(vec_dev(E, A, B), dE)
    return (status(d, d, TOL[0]), d, d, 'eigenvectors (degenerate clusters by projectors) & energies of %s' % ef)


def cmp_ncsm_trans(rp, dp, tol):
    """NCSM transition table NCSM_*_E<L>.dat: '#' header lines & the columns k, E, Ex, J, B_ph, B_is, B_iv signed; the
    amplitude columns M_p, M_n (the last two) up to one sign per row (phase of the state k)"""
    def parse(p):
        hn, hu, rn, ru, text = [], [], [], [], []
        with open(p, errors='replace') as f:
            for line in f:
                n, u, t = num_tokens(line.split())
                text += t
                if line.lstrip().startswith('#'):
                    hn += n
                    hu += u
                elif n:
                    rn.append(n)
                    ru.append(u)
        return (np.array(hn), np.array(hu), np.array(rn, dtype=float).reshape((len(rn), -1)),
                np.array(ru, dtype=float).reshape((len(ru), -1)), text)
    ha, hua, Ra, Ua, ta = parse(rp)
    hb, hub, Rb, Ub, tb = parse(dp)
    if ha.shape != hb.shape or Ra.shape != Rb.shape:
        return 'SIZE', float('inf'), float('inf'), 'header numbers / table shape differ'
    dh = dev_units(ha, hb, np.maximum(hua, hub))
    U = np.maximum(Ua, Ub)
    d = max(dh, dev_units(Ra, Rb, U))
    nflip = 0
    da = d
    if Ra.shape[0] and Ra.shape[1] >= 9:
        Ma, Mb = Ra[:, -2:], Rb[:, -2:]
        sg = np.where(np.sum(Ma * Mb, axis=1) >= 0.0, 1.0, -1.0)
        nflip = int(np.sum(sg < 0))
        da = max(dh, dev_units(Ra[:, :-2], Rb[:, :-2], U[:, :-2]), dev_units(Ma, sg[:, None] * Mb, U[:, -2:]))
    st = 'ok' if d <= tol else ('SIGN' if da <= tol else 'DIFF')
    if ta != tb and st == 'ok':
        st = 'TEXT'
    note = ('M_p, M_n up to a sign per row (%d of %d rows flipped), other columns signed, per-number printed unit'
            % (nflip, Ra.shape[0]))
    return st, d, da, note


GAUGE_FILES = re.compile(r'(/HF/bin/(C|H1b|H2b)\.bin|/EMPM/bin/EMPM_(TDA_X|F)\.bin|/EMPM/TDA_Structure\.dat)$')


def has_gauge(dev_rel):
    """files whose numbers carry a sign gauge (eigenvector / orbital / phonon signs): a sign-only change is SIGN"""
    return GAUGE_FILES.search('/' + dev_rel) is not None


def cmp_files(rp, dp, tol, dev_rel):
    name = os.path.basename(rp)
    note = ''
    gauge = has_gauge(dev_rel)
    if name.endswith('.bin'):
        if name.startswith('EMPM_S') and re.match(r'EMPM_S\d_', name):
            r = cmp_states(rp, dp)
            if r is not None:
                st, d, da, note = r
                if st:
                    return st, d, da, note
                return status(d, da, tol), d, da, note
        r = cmp_eigfile(rp, dp, dev_rel, None, None)
        if r is not None:
            return r
        a = np.fromfile(rp, dtype='<f8')
        b = np.fromfile(dp, dtype='<f8')
        if a.size != b.size:
            return 'SIZE', float('inf'), float('inf'), '%d vs %d values' % (a.size, b.size)
        d, da = rel(a, b), rel(np.abs(a), np.abs(b))
        return status(d, da, tol, gauge), d, da, ('sign gauge' if gauge else '')
    if re.match(r'NCSM_.*_E\d+\.dat$', name):
        return cmp_ncsm_trans(rp, dp, tol)
    if name.lower() == 'hf_summary.dat':
        # dev HF summary (MBPT block, other text) vs the legacy HF_Summary.dat ... the E_HF line only ...
        def ehf(path):
            return np.array([float(x) for x in re.findall(r'E_HF\s*=\s*([-+0-9.eE]+)', open(path).read())[:1]])
        ea, eb = ehf(rp), ehf(dp)
        d = rel(ea, eb) if ea.size and eb.size else float('inf')
        return status(d, d, tol), d, d, 'E_HF only (dev summary has other text & the MBPT block)'
    if name.startswith('Basis_'):
        ia, ib = basis_ind(rp), basis_ind(dp)
        A, B = np.loadtxt(rp, skiprows=3, ndmin=2), np.loadtxt(dp, skiprows=3, ndmin=2)
        if A.shape != B.shape:
            return 'SIZE', float('inf'), float('inf'), 'rows'
        A[:, 3] = 0
        B[:, 3] = 0
        d = rel(A, B)
        d = max(0.0, d - 1e-5) if d <= 1e-5 * 1.0000001 + tol else d
        if np.array_equal(ia, ib):
            return status(d, d, tol), d, d, 'printed unit 1e-05'
        return ('GAUGE' if d <= tol else 'DIFF'), d, d, ('Ind flags differ in %d of %d rows (pivot choice among dependent '
                                                         'rows), other columns compared' % (int(np.sum(ia != ib)), ia.size))
    a, ua, ta = tokens(rp, drop_seconds=True)
    b, ub, tb = tokens(dp, drop_seconds=True)
    ta = [re.sub(r'^io/([^/]+)/Bin/', r'io/\1/EMPM/bin/', t) for t in ta]
    if name in ('TDA_Summary.dat', 'TDA_summary.dat'):
        # only the dimension table is common ... compare its numbers (from 'Information on dimensions' on) ...
        def dims(path):
            s = open(path).read()
            i = s.find('Information on dimensions')
            return np.array([float(t) for t in re.findall(r'(?<![\w.])\d+(?![\w.])', s[i:])])
        da_ = dims(rp)
        db_ = dims(dp)
        if da_.shape != db_.shape:
            return 'DIFF', float('inf'), float('inf'), 'dimension table differs'
        d = rel(da_, db_)
        return ('SKIP' if d == 0 else 'DIFF'), d, d, 'dev TDA_summary text; dimension table compared'
    if a.size != b.size:
        return 'SIZE', float('inf'), float('inf'), '%d vs %d numbers' % (a.size, b.size)
    u = np.maximum(ua, ub)
    d, da = dev_units(a, b, u), dev_units(np.abs(a), np.abs(b), u)
    if ta != tb:
        note = 'text tokens differ'
    st = status(d, da, tol, gauge)
    if ta != tb and st == 'ok':
        st = 'TEXT'
    if gauge:
        note = (note + ', ' if note else '') + 'sign gauge'
    return st, d, da, note


TOL = [1e-10]


def status(d, da, tol, gauge=False):
    """ok, SIGN (equal in |.|; only for a file with a sign gauge, else DIFF), DIFF"""
    if d <= tol:
        return 'ok'
    if gauge and da <= tol:
        return 'SIGN'
    return 'DIFF'


SPUR_FILES = ('EMPM/bin/EMPM_TDA_E.bin', 'EMPM/spectra/TDA.dat', 'EMPM/transitions/TDA_E1.dat', 'EMPM/TDA_Structure.dat')


def spur_cmp(rp, dp, relp, tol):
    """compare with the spurious E(1) excluded ... returns (d, d_abs, E_spur_dev) or None"""
    if relp.endswith('EMPM_TDA_E.bin'):
        a = np.fromfile(rp, dtype='<f8')
        b = np.fromfile(dp, dtype='<f8')
        if a.size != b.size:
            return None
        k = np.nonzero((a == 0.0) & (b != 0.0))[0]
        if k.size != 1:
            return None
        m = np.ones(a.size, bool)
        m[k] = False
        return rel(a[m], b[m]), rel(np.abs(a[m]), np.abs(b[m])), b[k[0]]
    a, ua, ta = tokens(rp)
    b, ub, tb = tokens(dp)
    if a.size != b.size or ta != tb:
        return None
    k = np.nonzero((a == 0.0) & (np.abs(b) > 0.0) & (np.abs(a - b) > 1e-3))[0]
    if k.size != 1:
        return None
    m = np.ones(a.size, bool)
    m[k] = False
    u = np.maximum(ua, ub)[m]
    return dev_units(a[m], b[m], u), dev_units(np.abs(a[m]), np.abs(b[m]), u), b[k[0]]


def legacy_mode(base, dev, tol, quiet):
    files = []
    for root, _, fs in os.walk(base):
        for f in fs:
            files.append(os.path.relpath(os.path.join(root, f), base))
    files.sort()
    counts = {}
    bad = 0
    print('%-62s %-6s %10s %10s  %s' % ('file (dev name)', 'status', 'd', 'd_abs', 'note'))
    for f in files:
        m = map_legacy(f)
        if m is None:
            continue
        rp, dp = os.path.join(base, f), os.path.join(dev, m)
        if not os.path.exists(dp):
            st, d, da, note = 'MISSING', float('nan'), float('nan'), f
        else:
            st, d, da, note = cmp_files(rp, dp, tol, m)
            if st in ('DIFF', 'SIZE') and any(m.endswith(s) for s in SPUR_FILES):
                r = spur_cmp(rp, dp, m, tol)
                if r is not None and (r[0] <= tol or (has_gauge(m) and r[1] <= tol)):
                    st = 'ok*' if r[0] <= tol else 'SIGN*'
                    d, da = r[0], r[1]
                    note = 'spurious E(1) excluded (legacy 0, dev %.6f MeV)' % r[2]
        counts[st] = counts.get(st, 0) + 1
        if st in ('DIFF', 'SIZE', 'MISSING'):
            bad += 1
        if not quiet or st not in ('ok',):
            print('%-62s %-6s %10.2e %10.2e  %s' % (m, st, d, da, note))
    print('summary: ' + ', '.join('%s %d' % (k, v) for k, v in sorted(counts.items())))
    return 1 if bad else 0


# ---------------------------------------------------------------------------------------------------------------------
def spectra_channels(path):
    ch = {}
    with open(path) as f:
        next(f)
        for line in f:
            t = line.split()
            if len(t) >= 3:
                ch.setdefault((int(t[0]), t[1]), []).append(float(t[2]))
    return ch


def tda_mode(dev, tol):
    worst = 0.0
    missing = 0

    def have(*fs):
        nonlocal missing
        ok = True
        for f in fs:
            if not os.path.exists(f):
                print('%-8s MISSING %s' % (path, os.path.relpath(f, dev)))
                missing += 1
                ok = False
        return ok

    for path in sorted(os.listdir(dev)):
        e, t = os.path.join(dev, path, 'EMPM'), os.path.join(dev, path, 'TDA')
        if not (os.path.isdir(e) and os.path.isdir(t)):
            continue
        fa, fb = os.path.join(e, 'spectra', 'TDA.dat'), os.path.join(t, 'spectra', 'TDA.dat')
        if have(fa, fb):
            a, b = spectra_channels(fa), spectra_channels(fb)
            if set(a) != set(b):
                print('%-8s spectra/TDA.dat channel sets differ: %s' % (path, sorted(set(a) ^ set(b))))
                worst = float('inf')
            else:
                d = max((rel(a[k], b[k]) if len(a[k]) == len(b[k]) else float('inf') for k in a), default=0.0)
                n = sum(len(v) for v in a.values())
                print('%-8s spectra/TDA.dat          %4d levels, %2d channels, max dE / max(1,|E|) = %.2e'
                      % (path, n, len(a), d))
                worst = max(worst, d)
        # full-precision energies: EMPM/bin/EMPM_TDA_E.bin & TDA/bin/TDA_E.bin (both J outer, P inner, empty channels skipped)
        fa, fb = os.path.join(e, 'bin', 'EMPM_TDA_E.bin'), os.path.join(t, 'bin', 'TDA_E.bin')
        if have(fa, fb):
            Ea, Eb = np.fromfile(fa, dtype='<f8'), np.fromfile(fb, dtype='<f8')
            dE = rel(Ea, Eb)
            print('%-8s bin EMPM_TDA_E vs TDA_E   %4d values, max dE / max(1,|E|) = %.2e' % (path, Ea.size, dE))
            worst = max(worst, dE)
        for L in range(4):
            fe = os.path.join(e, 'transitions', 'TDA_E%d.dat' % L)
            ft = os.path.join(t, 'transitions', 'E%d' % L, 'TDA_E%d.dat' % L)
            if not have(fe, ft):
                continue
            A = np.loadtxt(fe, skiprows=1, ndmin=2)
            B = np.loadtxt(ft, skiprows=1, ndmin=2)
            if A.shape != B.shape:
                print('%-8s transitions/TDA_E%d.dat shapes differ %s %s' % (path, L, A.shape, B.shape))
                worst = float('inf')
                continue
            dd = dev_units(A, B, np.full(A.shape, 1e-8))   # 8 printed decimals
            print('%-8s transitions/TDA_E%d.dat  %4d rows, max dev (E, B_ph, B_is, B_iv) beyond the printed unit 1e-08 = '
                  '%.2e' % (path, L, A.shape[0], dd))
            worst = max(worst, dd)
    print('EMPM_TDA vs TDA worst deviation = %.2e (tol %.0e), %d MISSING' % (worst, tol, missing))
    return 1 if (worst > tol or missing) else 0


def main():
    args = sys.argv[1:]
    tol = 1e-10
    quiet = '--quiet' in args
    if '--tol' in args:
        i = args.index('--tol')
        tol = float(args[i + 1])
        del args[i:i + 2]
    args = [a for a in args if a != '--quiet']
    TOL[0] = tol
    if args and args[0] == '--tda':
        sys.exit(tda_mode(args[1], tol))
    if len(args) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(legacy_mode(args[0], args[1], tol, quiet))


if __name__ == '__main__':
    main()
