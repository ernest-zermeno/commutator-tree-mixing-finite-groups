"""Intersection-six stratum of exterior octahedral pairs in PGL2(F_q): exact checks.

For prime q = 5, 11, 13, 19, 29, 37 (q = +-3 mod 8) this script verifies, with
exact arithmetic in the projective Hamilton quaternions of octahedral_pair_check.py:

  * the census of intersection orders |A cap gAg^{-1}| over all double cosets
    A\\G/A, the bound <= 6, and the uniqueness of the double coset of order six
    (it is AhA, h = i+j+k), hence exactly four subgroups B with |A cap B| = 6;
  * N_G(V) = A, N_G(A) = A, C_G(A) = 1, C_G(K) = {1, h}, N_G(K) = K u Kh, by
    brute force over G;
  * the adjoint (SO3) dictionary used in the manuscript: A acts by the signed
    permutation matrices of determinant one, V diagonally, K by sgn(pi) P_pi,
    h by (2/3)J - I; conjugation has determinant one on all of G for q <= 13;
  * the certified natural value c_q: the double-coset quotient of the Hecke graph
    has an integer characteristic polynomial; SymPy factors it over Z, isolates
    its real roots in rational intervals and refines the extremal one. Every
    distinct eigenvalue of the Hecke graph occurs in this quotient (the
    spectral projector of a vertex-transitive graph has constant positive
    diagonal, so it maps the root delta to a nonzero A-fixed eigenvector);
  * the all-pairs maxima of the manuscript table, computed by a single numerical
    procedure for every prime: for each nontrivial double coset D the coset
    operator T_D on l2(G/A) commutes with the left G-action, so its singular
    values are those of its restriction to A-invariant functions, i.e. of the
    orbit-weighted quotient matrix N_D(i,j) sqrt(w_i/w_j)/d. One copy of the
    singular value one (the constants) is deflated. These values are floating
    point; for q <= 19 they are cross-checked against dense SVDs on l2(G/A).

Usage: python octahedral_cap6_check.py [q ...] [--dense-all] [--write PATH]
Exit status is nonzero if any assertion fails.
"""
from collections import Counter
from fractions import Fraction
import argparse
import json
import math
import os
import sys
import time

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from octahedral_pair_check import field_group, closure  # noqa: E402

EXPECTED_NATURAL = {
    5: ('1/4', 'z + 1'),
    11: ('1/4 + sqrt(3)/4', 'z**2 + 2*z - 2'),
    13: ('3/4', 'z + 3'),
    19: ('1/4 + sqrt(5)/4', 'z**2 + 2*z - 4'),
}
EXPECTED_20_DIGITS = {
    29: '0.85436118917716762327',
    37: '0.84565136442071691397',
}
EXPECTED_TABLE = {11: '0.683013', 13: '0.750000', 19: '0.809017',
                  29: '0.854361', 37: '0.845651'}


def qmul_raw(a, b, q):
    w, x, y, z = a
    r, s, t, u = b
    return ((w*r - x*s - y*t - z*u) % q, (w*s + x*r + y*u - z*t) % q,
            (w*t - x*u + y*r + z*s) % q, (w*u + x*t - y*s + z*r) % q)


def adjoint_matrix(alpha, q):
    """Matrix of x -> alpha x alpha^{-1} on pure quaternions, basis (i, j, k)."""
    N = sum(c * c for c in alpha) % q
    Ninv = pow(N, -1, q)
    conj = (alpha[0], -alpha[1] % q, -alpha[2] % q, -alpha[3] % q)
    cols = []
    for x in ((0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)):
        y = qmul_raw(qmul_raw(alpha, x, q), conj, q)
        assert y[0] == 0
        cols.append(tuple(c * Ninv % q for c in y[1:]))
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def det3(M, q):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    return (a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)) % q


def signed_perm_form(M, q):
    """(perm, signs) if M is a signed permutation matrix with entries in {1, -1}."""
    perm, signs = [], []
    for r in range(3):
        nz = [c for c in range(3) if M[r][c]]
        if len(nz) != 1:
            return None
        v = M[r][nz[0]]
        s = 1 if v == 1 else (-1 if v == q - 1 else None)
        if s is None:
            return None
        perm.append(nz[0])
        signs.append(s)
    return (tuple(perm), tuple(signs)) if sorted(perm) == [0, 1, 2] else None


def perm_sign(p):
    return 1 if sum(1 for i in range(3) for j in range(i+1, 3) if p[i] > p[j]) % 2 == 0 else -1


def deflated_top_singular(M, u):
    """Largest singular value of (I-uu^T) M (I-uu^T) for a unit vector u."""
    P = np.eye(len(u)) - np.outer(u, u)
    return float(np.linalg.svd(P @ M @ P, compute_uv=False)[0])


def check(q, dense_all=False):
    t0 = time.time()
    norm, mul, inv, G = field_group(q)
    E = (1, 0, 0, 0)
    out = {'q': q, 'group_order': len(G)}
    i_, j_, k_ = norm((0, 1, 0, 0)), norm((0, 0, 1, 0)), norm((0, 0, 0, 1))
    A = closure(mul, [i_, norm((1, 1, 0, 0)), norm((1, 1, 1, 1))])
    assert len(A) == 24
    h = norm((0, 1, 1, 1))
    S = {norm((0, 1, 1, 1)), norm((0, 1, 1, -1)), norm((0, 1, -1, 1)), norm((0, 1, -1, -1))}
    assert len(S) == 4 and h in S and all(mul(s, s) == E for s in S)
    assert all(mul(mul(a, s), inv(a)) in S for a in A for s in S)
    V = {E, i_, j_, k_}
    assert V <= A
    K = {a for a in A if mul(mul(a, h), inv(a)) == h}
    hAh = {mul(mul(h, a), h) for a in A}
    assert len(K) == 6 and A & hAh == K and h not in A
    out['K_order'] = 6

    # --- adjoint dictionary ---
    adjA = {a: adjoint_matrix(a, q) for a in A}
    assert all(signed_perm_form(M, q) is not None and det3(M, q) == 1 for M in adjA.values())
    assert len({tuple(map(tuple, M)) for M in adjA.values()}) == 24
    assert all(all(adjA[v][r][c] == 0 for r in range(3) for c in range(3) if r != c) for v in V)
    # K consists of sgn(pi) P_pi: one global sign per matrix, equal to the parity of pi
    assert all((lambda f: f is not None and f[1] == (perm_sign(f[0]),) * 3)(signed_perm_form(adjA[k], q)) for k in K)
    inv3 = pow(3, -1, q)
    Hm = adjoint_matrix(h, q)
    assert all(Hm[r][c] == ((2 * inv3 - (1 if r == c else 0)) % q) for r in range(3) for c in range(3))
    if q <= 13:
        assert all(det3(adjoint_matrix(g, q), q) == 1 for g in G)
        out['adjoint_det_one_checked_on_all_of_G'] = True
    out['adjoint_dictionary'] = 'A -> signed permutations det 1; V diagonal; K = sgn(pi)P_pi; h -> (2/3)J - I'

    # --- cosets, double cosets, intersection census ---
    cid, reps = {}, []
    for g in G:
        if g in cid:
            continue
        j = len(reps)
        reps.append(g)
        for a in A:
            cid[mul(g, a)] = j
    n = len(reps)
    assert n * 24 == len(G)
    root = cid[E]
    orbit_of, orbits = {}, []
    for idx, x in enumerate(reps):
        if idx in orbit_of:
            continue
        o = sorted({cid[mul(a, x)] for a in A})
        for j in o:
            orbit_of[j] = len(orbits)
        orbits.append(o)
    m = len(orbits)
    caps = []
    for o in orbits:
        x = reps[o[0]]
        xi = inv(x)
        caps.append(len(A & {mul(mul(x, a), xi) for a in A}))
        assert caps[-1] * len(o) == 24
    natural = sorted({cid[mul(a, h)] for a in A})
    nat_idx = orbits.index(natural)
    census = Counter(c for idx, c in enumerate(caps) if idx != orbit_of[root])
    assert max(census) <= 6 and census[6] == 1 and caps[nat_idx] == 6
    out['cosets'] = n
    out['double_cosets_nontrivial'] = m - 1
    out['intersection_census'] = dict(sorted(census.items()))
    out['unique_order_six_double_coset_is_AhA'] = True
    Bs = set()
    for j in orbits[nat_idx]:
        x = reps[j]
        xi = inv(x)
        Bs.add(frozenset(mul(mul(x, a), xi) for a in A))
    assert len(Bs) == 4 and Bs == {frozenset(mul(mul(a, b), inv(a)) for b in hAh) for a in A}
    out['subgroups_B_with_intersection_six'] = 4

    # --- normalisers and centralisers ---
    NV = NA = CA = CK = NK = 0
    CKset = set()
    Kh = {mul(k, h) for k in K}
    for g in G:
        gi = inv(g)
        cj = lambda T: {mul(mul(g, t), gi) for t in T}
        if cj(V) == V:
            NV += 1
        if cj(A) == A:
            NA += 1
        if all(mul(g, a) == mul(a, g) for a in A):
            CA += 1
        if all(mul(g, k) == mul(k, g) for k in K):
            CK += 1
            CKset.add(g)
        if cj(K) == K:
            NK += 1
    assert (NV, NA, CA, CK, NK) == (24, 24, 1, 2, 12) and CKset == {E, h}
    out['N_G_V'] = 24
    out['N_G_A'] = 24
    out['C_G_A'] = 1
    out['C_G_K'] = '{1, h}'
    out['N_G_K'] = 12
    assert {g for g in G if {mul(mul(g, k), inv(g)) for k in K} == K} == K | Kh

    # --- Hecke graph (right multiplication by S) and its double-coset quotient ---
    Sl = sorted(S)
    adj = [sorted({cid[mul(x, s)] for s in Sl}) for x in reps]
    assert all(len(a) == 4 and idx not in a for idx, a in enumerate(adj))
    assert all(idx in adj[j] for idx, a in enumerate(adj) for j in a)
    seen, todo = {root}, [root]
    while todo:
        for j in adj[todo.pop()]:
            if j not in seen:
                seen.add(j)
                todo.append(j)
    assert len(seen) == n
    w = [len(o) for o in orbits]
    Q = [[0] * m for _ in range(m)]
    for i, o in enumerate(orbits):
        c0 = Counter(orbit_of[j] for j in adj[o[0]])
        for v in o:
            assert Counter(orbit_of[j] for j in adj[v]) == c0
        for j, c in c0.items():
            Q[i][j] = c
    assert all(w[i] * Q[i][j] == w[j] * Q[j][i] for i in range(m) for j in range(m))
    z = sp.Symbol('z')
    cp = sp.Matrix(Q).charpoly(z).as_poly(z)
    nt = cp.exquo(sp.Poly(z - 4, z))
    assert nt.eval(4) != 0
    assert nt.count_roots() == sp.Poly(sp.sqf_part(nt.as_expr()), z).degree()
    factors = sp.factor_list(nt.as_expr())[1]
    intervals = nt.intervals(eps=sp.Rational(1, 10**15))
    best = None
    for (a, b), mult in intervals:
        lo = Fraction(0) if a <= 0 <= b else min(abs(Fraction(a)), abs(Fraction(b)))
        hi = max(abs(Fraction(a)), abs(Fraction(b)))
        if best is None or lo > best[0]:
            best = (lo, hi, a, b, mult)
    lo, hi, a, b, mult = best
    for (a2, b2), _ in intervals:
        if (a2, b2) == (a, b):
            continue
        assert max(abs(Fraction(a2)), abs(Fraction(b2))) < lo, 'extremal root not separated'
    hits = [(f, e) for f, e in factors if sp.Poly(f, z).count_roots(a, b) > 0]
    assert len(hits) == 1
    fpoly, fmult = hits[0]
    ra, rb = sp.Poly(fpoly, z).refine_root(a, b, eps=sp.Rational(1, 10**30))
    lo_c, hi_c = sorted((abs(sp.Rational(ra)) / 4, abs(sp.Rational(rb)) / 4))
    # c_q lies in [lo_c, hi_c], an interval of width < 10^-30; both ends round
    # to the same 20 significant digits, which are therefore certified.
    digits = str(sp.N(lo_c, 20))
    assert str(sp.N(hi_c, 20)) == digits
    out['c_q_interval'] = [str(lo_c), str(hi_c)]
    out['natural_quotient_size'] = m
    out['natural_charpoly_factors'] = [(str(f), int(e)) for f, e in factors]
    out['natural_extremal_root_sign'] = 'negative' if b < 0 else 'positive'
    out['natural_extremal_root_isolating_interval'] = [str(a), str(b)]
    out['natural_extremal_root_minimal_polynomial'] = str(fpoly)
    out['natural_extremal_root_minimal_polynomial_degree'] = int(sp.Poly(fpoly, z).degree())
    out['c_q_certified_20_digits_rounded'] = digits
    if q in EXPECTED_NATURAL:
        value, poly = EXPECTED_NATURAL[q]
        assert sp.Poly(fpoly, z) == sp.Poly(sp.sympify(poly), z)
        val = sp.sympify(value)
        assert lo_c <= val <= hi_c
        out['c_q_exact'] = value
    if q in EXPECTED_20_DIGITS:
        assert digits == EXPECTED_20_DIGITS[q]
    ev = np.linalg.eigvalsh(np.array([[1.0 if j in adj[i] else 0.0 for j in range(n)] for i in range(n)]))
    ev = sorted(ev, reverse=True)
    assert abs(ev[0] - 4) < 1e-9 and ev[1] < 4 - 1e-6
    numeric = max(abs(v) for v in ev[1:]) / 4
    assert abs(numeric - float(lo_c)) < 1e-9
    out['c_q_dense_numeric'] = round(numeric, 12)
    evq = sorted(np.linalg.eigvals(np.array(Q, dtype=float)).real)
    assert {round(v, 6) for v in evq} == {round(v, 6) for v in ev}
    out['quotient_contains_all_distinct_eigenvalues'] = True

    # --- all-pairs numerical maxima (orbit-weighted quotients) ---
    sw = np.sqrt(np.array(w, dtype=float))
    u = sw / np.linalg.norm(sw)
    per_coset = []
    for didx, D in enumerate(orbits):
        if didx == orbit_of[root]:
            continue
        d = len(D)
        ys = [reps[j] for j in D]
        N = np.zeros((m, m))
        for i, o in enumerate(orbits):
            x = reps[o[0]]
            for y in ys:
                N[i, orbit_of[cid[mul(x, y)]]] += 1
        assert np.allclose(N.sum(axis=1), d)
        M = N * (sw[:, None] / sw[None, :]) / d
        assert np.allclose(M @ u, u) and np.allclose(M.T @ u, u)
        s2 = deflated_top_singular(M, u)
        row = {'intersection_order': caps[didx], 'coset_degree': d, 'norm': round(s2, 12)}
        if dense_all or q <= 19:
            T = np.zeros((n, n))
            for i, x in enumerate(reps):
                for y in ys:
                    T[i, cid[mul(x, y)]] += 1.0 / d
            one = np.ones(n) / math.sqrt(n)
            s2d = deflated_top_singular(T, one)
            assert abs(s2d - s2) < 1e-8, (q, didx, s2, s2d)
            row['dense_norm'] = round(s2d, 12)
        per_coset.append((didx, row))
    best_idx, best_row = max(per_coset, key=lambda t: t[1]['norm'])
    assert best_idx == nat_idx
    assert all(r['norm'] < 1 - 1e-9 for _, r in per_coset)
    per_cap = {}
    for _, r in per_coset:
        per_cap[r['intersection_order']] = max(per_cap.get(r['intersection_order'], 0.0), r['norm'])
    out['allpairs_max_norm'] = round(best_row['norm'], 6)
    out['allpairs_max_at_intersection_order'] = 6
    out['allpairs_max_by_intersection_order'] = {str(k): round(v, 6) for k, v in sorted(per_cap.items())}
    out['allpairs_dense_cross_check'] = bool(dense_all or q <= 19)
    assert abs(best_row['norm'] - float(lo_c)) < 1e-9
    if q in EXPECTED_TABLE:
        assert '%.6f' % best_row['norm'] == EXPECTED_TABLE[q]
    out['seconds'] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('q', nargs='*', type=int, default=[5, 11, 13, 19, 29, 37])
    ap.add_argument('--dense-all', action='store_true', help='dense SVD cross-check for every prime')
    ap.add_argument('--write', help='write the JSON report to this path')
    args = ap.parse_args()
    rows = [check(q, args.dense_all) for q in args.q]
    report = {'python': sys.version.split()[0], 'numpy': np.__version__, 'sympy': sp.__version__,
              'results': rows}
    text = json.dumps(report, indent=1)
    print(text)
    print('table:', ' '.join('%d:%.6f' % (r['q'], r['allpairs_max_norm']) for r in rows if r['q'] in EXPECTED_TABLE))
    if args.write:
        with open(args.write, 'w') as f:
            f.write(text + '\n')


if __name__ == '__main__':
    main()
