"""
Bottleneck check for PGL_2(F_7) at level 2.

In the REGULAR representation, the level-2 commutator operator is block-diagonal
over irreps (each irrep rho appears d_rho times), so its singular values are the
union over irreps of the singular values of D_hat_2(rho). The trivial rep gives the
single largest singular value = 1 (the uniform/stationary direction); therefore

    second-largest singular value  =  max over NON-TRIVIAL irreps of op-norm.

If sv2 < 1 for a config, EVERY non-trivial irrep contracts there ⇒ the Diaconis-
Shahshahani bound is < (something)·1 ⇒ that config mixes toward uniform. We test this
over a sample of valid level-2 quads, including the worst ones in the standard rep.
"""
import numpy as np
from itertools import combinations
from derived_length_test import (pgl2, involutions, pmul, pinv, is_id, comm,
                                  std_basis, perm_matrix, opnorm)

rng = np.random.default_rng(0)

perms, m = pgl2(7)
assert len(perms) == 336
invs = involutions(perms)
N = len(perms)
print(f"PGL_2(F_7): |G|={N}, involutions={len(invs)}, P^1 size m={m}")

# ---- standard rep projections + enumerate valid level-2 quads ----
B = std_basis(m)
Pstd = {p: 0.5 * (B @ perm_matrix(p) @ B.T + np.eye(m - 1)) for p in invs}
pairs = [(s, t) for s, t in combinations(invs, 2) if not is_id(comm(s, t))]
ptgt = {pr: comm(pr[0], pr[1]) for pr in pairs}
quads = []
for i in range(len(pairs)):
    for j in range(i + 1, len(pairs)):
        if not is_id(comm(ptgt[pairs[i]], ptgt[pairs[j]])):
            quads.append((i, j))
print(f"valid level-2 quads: {len(quads)}")

def std_op(q):
    (s1, t1), (s2, t2) = pairs[q[0]], pairs[q[1]]
    DA = np.linalg.matrix_power(Pstd[s1] @ Pstd[t1], 2)
    DB = np.linalg.matrix_power(Pstd[s2] @ Pstd[t2], 2)
    DAi = np.linalg.matrix_power(Pstd[t1] @ Pstd[s1], 2)
    DBi = np.linalg.matrix_power(Pstd[t2] @ Pstd[s2], 2)
    return opnorm(DA @ DB @ DAi @ DBi)

# ---- regular rep projections (left multiplication), one per involution ----
gidx = {g: k for k, g in enumerate(perms)}
def reg_proj(sigma):
    rp = [gidx[pmul(sigma, g)] for g in perms]      # left mult: g -> sigma·g
    P = np.zeros((N, N))
    for col, row in enumerate(rp):
        P[row, col] = 1.0
    return 0.5 * (P + np.eye(N))
print("building regular-rep projections...")
Preg = {s: reg_proj(s) for s in invs}

def reg_top2(q):
    (s1, t1), (s2, t2) = pairs[q[0]], pairs[q[1]]
    DA = np.linalg.matrix_power(Preg[s1] @ Preg[t1], 2)
    DB = np.linalg.matrix_power(Preg[s2] @ Preg[t2], 2)
    DAi = np.linalg.matrix_power(Preg[t1] @ Preg[s1], 2)
    DBi = np.linalg.matrix_power(Preg[t2] @ Preg[s2], 2)
    sv = np.linalg.svd(DA @ DB @ DAi @ DBi, compute_uv=False)
    return float(sv[0]), float(sv[1])

# ---- sample: worst-200-by-standard-op-norm + 1500 uniform random ----
print("scoring standard-rep op-norms over all quads...")
std_scores = np.array([std_op(q) for q in quads])
order = np.argsort(-std_scores)
worst = list(order[:200])
rand = list(rng.choice(len(quads), size=min(1500, len(quads)), replace=False))
sample = sorted(set(worst) | set(rand))
print(f"regular-rep bottleneck over {len(sample)} sampled quads "
      f"(worst-200 + {len(rand)} random)...")

max_sv1 = 0.0; max_sv0 = 0.0; bad = 0
worst_sv1_overall = 0.0
for idx in sample:
    sv0, sv1 = reg_top2(quads[idx])
    max_sv0 = max(max_sv0, sv0)
    max_sv1 = max(max_sv1, sv1)
    if sv1 > 1 - 1e-9:
        bad += 1

print("\n" + "=" * 64)
print(f"  max  top  singular value (should be ~1, uniform dir): {max_sv0:.8f}")
print(f"  MAX  2nd  singular value (bottleneck op-norm)        : {max_sv1:.8f}")
print(f"  configs with bottleneck op-norm >= 1                 : {bad}/{len(sample)}")
print(f"  worst standard-rep op-norm in sample                 : {std_scores[order[0]]:.8f}")
print("=" * 64)
if bad == 0:
    print("  ⇒ Every sampled config: ALL non-trivial irreps contract at level 2.")
    print("    Evidence that PGL_2(F_7) MIXES to uniform at level 2 (not just std rep).")
else:
    print("  ⇒ Some config has a non-contracting irrep at level 2 — standard rep is NOT the bottleneck.")
