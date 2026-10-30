"""
GL_3(F_2) ~= PSL_2(7) (order 168, SIMPLE) acting on the 7 nonzero vectors of F_2^3
(= 7 points of PG(2,2)). Point-stabilizer = AGL_2(F_2) ~= S_4 (derived length 3).
Prediction/evidence: standard rep contracts at level 3; the regular-rep bottleneck
also drops in the sampled/capped level-3 run.
Regular rep is 168x168: the level-2 bottleneck is exhaustive only over REDUCED gate representatives; level 3 is CAPPED at 60k
configs (evidence, not a certificate).
"""
import itertools
import numpy as np
from derived_length_test import (run, pmul, is_id, comm, involutions,
                                  std_basis, perm_matrix, opnorm)
from itertools import combinations

def gl3f2():
    vecs = [v for v in itertools.product([0, 1], repeat=3) if any(v)]   # 7 nonzero
    vidx = {v: i for i, v in enumerate(vecs)}
    out = set()
    for e in itertools.product([0, 1], repeat=9):
        rows = (e[0:3], e[3:6], e[6:9]); per = []; ok = True
        for v in vecs:
            w = tuple(sum(rows[r][c]*v[c] for c in range(3)) % 2 for r in range(3))
            if not any(w): ok = False; break
            per.append(vidx[w])
        if ok and len(set(per)) == 7:
            out.add(tuple(per))
    return list(out), 7

perms, m = gl3f2()
assert len(perms) == 168, len(perms)
invs = involutions(perms)
def n_fixed(p): return sum(1 for i in range(len(p)) if p[i] == i)
prof = {}
for p in invs: prof[n_fixed(p)] = prof.get(n_fixed(p), 0) + 1
print(f"GL_3(F_2)=PSL_2(7): |G|={len(perms)}, m={m}, involutions={len(invs)} (by #fixed: {prof})")

# ---- standard rep (6-dim) contraction level ----
run("GL_3(F_2) on 7 points (involutions) [stab S_4, dl=3]",
    invs, m, predict_dl=3, do_level3=True, l3_cap=2_000_000)

# ---- regular-rep bottleneck (168x168): L2 EXHAUSTIVE, L3 CAPPED at 60k ----
N = len(perms); gidx = {g: k for k, g in enumerate(perms)}
def reg_proj(s):
    rp = [gidx[pmul(s, g)] for g in perms]; P = np.zeros((N, N))
    for col, row in enumerate(rp): P[row, col] = 1.0
    return 0.5 * (P + np.eye(N))
Preg = {s: reg_proj(s) for s in invs}
def blk(s, t): return np.linalg.matrix_power(Preg[s] @ Preg[t], 2)
def sv2(M): sv = np.linalg.svd(M, compute_uv=False); return float(sv[1])

pairs = [(s, t) for s, t in combinations(invs, 2) if not is_id(comm(s, t))]
ptgt = {pr: comm(pr[0], pr[1]) for pr in pairs}
max1 = max(sv2(blk(s, t)) for (s, t) in pairs)

quads = []
for i in range(len(pairs)):
    for j in range(i+1, len(pairs)):
        if not is_id(comm(ptgt[pairs[i]], ptgt[pairs[j]])):
            quads.append((i, j))
D2 = {}; qtgt = {}; max2 = 0.0
for (i, j) in quads:
    (s1, t1), (s2, t2) = pairs[i], pairs[j]
    D = blk(s1, t1) @ blk(s2, t2) @ blk(t1, s1) @ blk(t2, s2)
    D2[(i, j)] = D; qtgt[(i, j)] = comm(ptgt[pairs[i]], ptgt[pairs[j]])
    max2 = max(max2, sv2(D))

max3 = 0.0; n3 = 0
ql = list(quads)
for a in range(len(ql)):
    for b in range(a+1, len(ql)):
        if is_id(comm(qtgt[ql[a]], qtgt[ql[b]])): continue
        n3 += 1
        if n3 > 60000: break
        max3 = max(max3, sv2(D2[ql[a]] @ D2[ql[b]] @ D2[ql[a]].T @ D2[ql[b]].T))
    if n3 > 60000: break

print("\n" + "=" * 60)
print("  FULL regular-rep bottleneck (2nd singular value = max op-norm over irreps):")
print(f"   L1={max1:.6f}   L2={max2:.6f} ({len(quads)} quads, EXHAUSTIVE)   "
      f"L3={max3:.6f} ({min(n3,60000)} configs, CAPPED at 60k -> evidence, not a certificate)")
first = next((L for L, v in [(1, max1), (2, max2), (3, max3)] if v < 1-1e-6), None)
print(f"  => GL_3(F_2) mixing level = {first}   (L1/L2 exhaustive; L3 CAPPED)   (predicted dl(S_4)=3)")
print("=" * 60)
