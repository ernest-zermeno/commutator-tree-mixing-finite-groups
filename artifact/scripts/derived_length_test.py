"""Exploratory augmentation-representation sweeps, with involution leaves.
Level-2 counts use reduced gate representatives, not every ordered tuple.
Level-3 sweeps may be capped. These numerical extrema and first observed levels
are cross-checks, not proofs of an exact all-configuration contraction level.
"""
import sys
import numpy as np
from itertools import combinations

# ---------- permutations ----------
def pmul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def pinv(p):
    r = [0]*len(p)
    for i, j in enumerate(p): r[j] = i
    return tuple(r)
def is_id(p): return all(p[i] == i for i in range(len(p)))
def comm(a, b): return pmul(pmul(a, b), pmul(pinv(a), pinv(b)))
def n_fixed(p): return sum(1 for i in range(len(p)) if p[i] == i)

def transpositions(n):
    out = []
    for i in range(n):
        for j in range(i+1, n):
            p = list(range(n)); p[i], p[j] = p[j], p[i]; out.append(tuple(p))
    return out

def pgl2(q):
    """PGL_2(F_q) as permutations of P^1(F_q) = {0..q-1} U {inf=index q}."""
    pts = list(range(q)) + ['inf']; idx = {pt: i for i, pt in enumerate(pts)}
    def ap(a, b, c, d, pt):
        if pt == 'inf':
            return 'inf' if c % q == 0 else (a * pow(c % q, -1, q)) % q
        num = (a*pt + b) % q; den = (c*pt + d) % q
        return 'inf' if den == 0 else (num * pow(den, -1, q)) % q
    S = set()
    for a in range(q):
        for b in range(q):
            for c in range(q):
                for d in range(q):
                    if (a*d - b*c) % q == 0: continue
                    S.add(tuple(idx[ap(a, b, c, d, pt)] for pt in pts))
    return list(S), q + 1

def involutions(perms):
    return [p for p in perms if (not is_id(p)) and is_id(pmul(p, p))]

# ---------- augmentation representation ----------
def std_basis(m):
    B = np.zeros((m-1, m))
    for k in range(1, m):
        B[k-1, :k] = 1.0
        B[k-1, k] = -k
        B[k-1] /= np.linalg.norm(B[k-1])
    return B

def perm_matrix(p):
    m = len(p); P = np.zeros((m, m))
    for i in range(m): P[p[i], i] = 1.0
    return P

def opnorm(M): return float(np.linalg.svd(M, compute_uv=False)[0])

# ---------- experiment ----------
def run(name, leaves, m, predict_dl, do_level3=False, l3_cap=3_000_000):
    B = std_basis(m)
    Pc = {}
    for p in leaves:
        M = B @ perm_matrix(p) @ B.T
        Pc[p] = 0.5 * (M + np.eye(m-1))
    # sanity: projection
    assert np.allclose(Pc[leaves[0]] @ Pc[leaves[0]], Pc[leaves[0]], atol=1e-9)
    fixed_profile = sorted(set(n_fixed(p) for p in leaves))

    pairs = [(s, t) for s, t in combinations(leaves, 2) if not is_id(comm(s, t))]
    D1 = {}; l1_one = 0
    for (s, t) in pairs:
        D = np.linalg.matrix_power(Pc[s] @ Pc[t], 2); D1[(s, t)] = D
        if opnorm(D) > 1 - 1e-9: l1_one += 1

    ptgt = {pr: comm(pr[0], pr[1]) for pr in pairs}
    quads = []
    for i in range(len(pairs)):
        for j in range(i+1, len(pairs)):
            if not is_id(comm(ptgt[pairs[i]], ptgt[pairs[j]])):
                quads.append((i, j))
    D2 = {}; qtgt = {}; l2_one = 0; max2 = 0.0
    for (i, j) in quads:
        (s1, t1) = pairs[i]; (s2, t2) = pairs[j]
        DA = np.linalg.matrix_power(Pc[s1] @ Pc[t1], 2)
        DB = np.linalg.matrix_power(Pc[s2] @ Pc[t2], 2)
        DAi = np.linalg.matrix_power(Pc[t1] @ Pc[s1], 2)
        DBi = np.linalg.matrix_power(Pc[t2] @ Pc[s2], 2)
        D = DA @ DB @ DAi @ DBi; D2[(i, j)] = D
        qtgt[(i, j)] = comm(ptgt[pairs[i]], ptgt[pairs[j]])
        o = opnorm(D); max2 = max(max2, o)
        if o > 1 - 1e-9: l2_one += 1

    print(f"\n### {name}")
    print(f"  leaves={len(leaves)} (fixed-pts {fixed_profile})  predict dl(stab)={predict_dl}")
    print(f"  L1: {len(pairs)} pairs, op-norm=1 in {l1_one}/{len(pairs)}")
    print(f"  L2: {len(quads)} valid quads, op-norm=1 in {l2_one}/{len(quads)}, max op-norm={max2:.6f}")
    first = None
    if l1_one == 0 and len(pairs) > 0: first = 1
    elif l2_one == 0 and len(quads) > 0: first = 2

    if do_level3 or (l2_one > 0 and len(quads) > 0):
        l3 = 0; l3_one = 0; max3 = 0.0
        qlist = list(quads)
        for a in range(len(qlist)):
            for b in range(a+1, len(qlist)):
                if is_id(comm(qtgt[qlist[a]], qtgt[qlist[b]])): continue
                l3 += 1
                if l3 > l3_cap: break
                D = D2[qlist[a]] @ D2[qlist[b]] @ D2[qlist[a]].T @ D2[qlist[b]].T
                o = opnorm(D); max3 = max(max3, o)
                if o > 1 - 1e-9: l3_one += 1
            if l3 > l3_cap: break
        capped = " (CAPPED)" if l3 > l3_cap else ""
        print(f"  L3: {min(l3,l3_cap)} valid configs{capped}, op-norm=1 in {l3_one}, max op-norm={max3:.6f}")
        if first is None and l3_one == 0 and l3 > 0: first = 3
    print(f"  => first fully-contracting level = {first}   (predicted {predict_dl})   "
          f"{'MATCH' if first == predict_dl else 'CHECK'}")
    return first


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "small"

    if which in ("small", "s5_5"):
        # S_5 on 5 points, transposition leaves. Stabilizer S_4, dl=3. (paper baseline)
        run("S_5 on 5 points (transpositions) [paper baseline]",
            transpositions(5), 5, predict_dl=3, do_level3=True)

    if which in ("small", "s5_6"):
        perms5, m5 = pgl2(5)
        assert len(perms5) == 120, len(perms5)   # |PGL_2(F_5)| = 120 ~ S_5
        run("S_5 on 6 points = PGL_2(F_5) on P^1 (involutions) [stab F_20]",
            involutions(perms5), m5, predict_dl=2)

    if which in ("small", "s6_6"):
        # S_6 on 6 points, transpositions. Stabilizer S_5 (non-solvable). No contraction.
        run("S_6 on 6 points (transpositions) [stab S_5 non-solvable]",
            transpositions(6), 6, predict_dl=None, do_level3=True, l3_cap=200000)

    if which in ("pgl7",):
        perms7, m7 = pgl2(7)
        assert len(perms7) == 336, len(perms7)   # |PGL_2(F_7)| = 336
        run("PGL_2(F_7) on 8 points (involutions) [stab Borel, dl=2]",
            involutions(perms7), m7, predict_dl=2, l3_cap=300000)
