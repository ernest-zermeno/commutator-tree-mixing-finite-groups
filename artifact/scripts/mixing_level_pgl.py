"""Sampled regular-representation norms for EXTERIOR involutions of PGL2(F_p).
Only odd primes p>=5 are implemented here. Exterior means determinant non-square;
split type depends on p modulo 4 and is not the selector. --class-only skips norms.
"""
import sys
import numpy as np
from itertools import combinations
from derived_length_test import (pgl2, involutions, pmul, pinv, is_id, comm,
                                  std_basis, perm_matrix, opnorm)

rng = np.random.default_rng(0)
q = int(sys.argv[1]) if len(sys.argv) > 1 else 7
if q < 5 or q % 2 == 0 or any(q % d == 0 for d in range(3, int(q**0.5)+1, 2)):
    raise SystemExit("requires an odd prime q >= 5; extension fields use pgl2_classification.py")
perms, m = pgl2(q)
N = len(perms)

def n_fixed(p): return sum(1 for i in range(len(p)) if p[i] == i)
def parity(p):
    seen = [False]*len(p); cyc = 0
    for i in range(len(p)):
        if not seen[i]:
            cyc += 1; j = i
            while not seen[j]: seen[j] = True; j = p[j]
    return (-1) ** (len(p) - cyc)

invs = involutions(perms)
# classify: confirm fixed-points <-> parity/sign
byfix = {}
for p in invs:
    byfix.setdefault((n_fixed(p), parity(p)), 0)
    byfix[(n_fixed(p), parity(p))] += 1
print(f"PGL_2(F_{q}): |G|={N}, m={m}, involutions={len(invs)}")
print(f"  involution classes (n_fixed, parity/sign): {dict(byfix)}")

from pgl2_classification import make_field, pgl2 as matrix_pgl2
info = matrix_pgl2(make_field(q))["info"]
split = [p for p in invs if not info[p][0]]  # determinant non-square
assert all(parity(p) == -1 for p in split)
print(f"  using EXTERIOR (non-square determinant) leaves: {len(split)}")

if "--class-only" in sys.argv:
    raise SystemExit(0)

# ---- regular rep projections for the split leaves ----
gidx = {g: k for k, g in enumerate(perms)}
def reg_proj(s):
    rp = [gidx[pmul(s, g)] for g in perms]
    P = np.zeros((N, N))
    for col, row in enumerate(rp): P[row, col] = 1.0
    return 0.5 * (P + np.eye(N))
Preg = {s: reg_proj(s) for s in split}

def top2(M):
    sv = np.linalg.svd(M, compute_uv=False); return float(sv[0]), float(sv[1])
def block(s, t):  # (P_s P_t)^2 in regular rep
    return np.linalg.matrix_power(Preg[s] @ Preg[t], 2)

# ---- level 1 (all non-commuting pairs) ----
pairs = [(s, t) for s, t in combinations(split, 2) if not is_id(comm(s, t))]
ptgt = {pr: comm(pr[0], pr[1]) for pr in pairs}
max1 = 0.0
for (s, t) in pairs:
    _, sv1 = top2(block(s, t)); max1 = max(max1, sv1)
print(f"\n  L1: {len(pairs)} pairs, max bottleneck op-norm (2nd sv) = {max1:.6f}")

# ---- level 2 (sample of valid quads) ----
quads = []
for i in range(len(pairs)):
    for j in range(i+1, len(pairs)):
        if not is_id(comm(ptgt[pairs[i]], ptgt[pairs[j]])):
            quads.append((i, j))
qtgt = {}
sample2 = list(rng.choice(len(quads), size=min(800, len(quads)), replace=False))
D2cache = {}
max2 = 0.0
for idx in sample2:
    (i, j) = quads[idx]; (s1, t1), (s2, t2) = pairs[i], pairs[j]
    DA = block(s1, t1); DB = block(s2, t2)
    DAi = block(t1, s1); DBi = block(t2, s2)
    D2 = DA @ DB @ DAi @ DBi
    D2cache[idx] = D2
    qtgt[idx] = comm(ptgt[pairs[i]], ptgt[pairs[j]])
    _, sv1 = top2(D2); max2 = max(max2, sv1)
print(f"  L2: {len(quads)} valid quads ({len(sample2)} sampled), max bottleneck op-norm = {max2:.6f}")

# ---- level 3 (valid pairs among cached level-2 ops) ----
keys = list(D2cache.keys())
max3 = 0.0; n3 = 0
for a in range(len(keys)):
    for b in range(a+1, len(keys)):
        ka, kb = keys[a], keys[b]
        if is_id(comm(qtgt[ka], qtgt[kb])): continue
        n3 += 1
        if n3 > 4000: break
        D3 = D2cache[ka] @ D2cache[kb] @ D2cache[ka].T @ D2cache[kb].T
        _, sv1 = top2(D3); max3 = max(max3, sv1)
    if n3 > 4000: break
print(f"  L3: {min(n3,4000)} valid configs, max bottleneck op-norm = {max3:.6f}")

print("\n" + "=" * 60)
levels = [("1", max1), ("2", max2), ("3", max3)]
first = next((L for L, v in levels if v < 1 - 1e-6), None)
for L, v in levels:
    print(f"  level {L}: bottleneck op-norm = {v:.6f}  {'< 1  (all SAMPLED irreps contract)' if v < 1-1e-6 else '= 1  (some irrep preserved)'}")
print(f"  NOTE: SAMPLED (~800 quads, L3 capped 4000) -> strong evidence, not a certificate.")
print(f"  => PGL_2(F_{q}) full-mixing level (exterior leaves, sampled) = {first}")
print("=" * 60)
