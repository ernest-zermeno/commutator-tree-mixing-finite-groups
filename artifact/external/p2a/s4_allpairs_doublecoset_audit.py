#!/usr/bin/env python3
"""
Audit the proposed "all exterior S4 pairs" sharp bound.

For a fixed exterior octahedral subgroup A < PGL2(F_q), the true global
quantity for a conjugate B = g A g^{-1} is

    sup_{rho != 1} || P_A rho(g) P_A ||,

which is the second singular value of the normalized double-coset operator
P_A R_g P_A on l2(G/A).  This script enumerates A\G/A for small prime q,
bucketizes by |A cap gAg^{-1}|, and estimates the non-trivial singular norm
by power iteration on the sparse coset operator.

It is a numerical audit, not a proof.
"""

import math
import random
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, "scripts")
from s4_family_exact_sup import build, find_S4, closure  # noqa: E402


SQRT3_OVER_2 = math.sqrt(3.0) / 2.0


def cosets_right(B, A):
    mm, norm = B["mm"], B["norm"]
    cos_of = {}
    reps = []
    for g in B["G"]:
        if g in cos_of:
            continue
        C = frozenset(norm(mm(g, a)) for a in A)
        r = min(C)
        reps.append(r)
        for x in C:
            cos_of[x] = r
    return reps, {r: i for i, r in enumerate(reps)}, cos_of


def double_coset(B, A, g):
    mm, norm = B["mm"], B["norm"]
    return frozenset(norm(mm(mm(a, g), b)) for a in A for b in A)


def conjugate(B, A, g):
    mm, norm, minv = B["mm"], B["norm"], B["minv"]
    gi = minv(g)
    return frozenset(norm(mm(mm(g, a), gi)) for a in A)


def generates_full(B, A, Bg):
    # For prime q in this audit, exterior S4 is maximal; this check is defensive.
    return closure(B, list(A | Bg), cap=200) is None


def source_targets(B, A, reps, idx, cos_of, g):
    mm, norm = B["mm"], B["norm"]
    targets = []
    for r in reps:
        targets.append([idx[cos_of[norm(mm(mm(r, a), g))]] for a in A])
    return targets


def second_singular_power(targets, iters=80, restarts=4, seed=12345):
    """Power iteration for ||M|| on constants^\perp, where M averages over targets."""
    targets = np.asarray(targets, dtype=np.int64)
    n, deg = targets.shape
    rng = np.random.default_rng(seed)

    def proj(x):
        return x - np.mean(x)

    def apply_m(v):
        weights = np.repeat(v / deg, deg)
        return np.bincount(targets.ravel(), weights=weights, minlength=n)

    def apply_mt(w):
        return np.mean(w[targets], axis=1)

    best = 0.0
    for _ in range(restarts):
        v = proj(rng.normal(size=n))
        v /= np.linalg.norm(v)
        last = 0.0
        for _ in range(iters):
            mv = proj(apply_m(v))
            av = proj(apply_mt(mv))
            na = np.linalg.norm(av)
            if na < 1e-14:
                last = 0.0
                break
            v = av / na
            last = float(np.linalg.norm(proj(apply_m(v))))
        best = max(best, last)
    return best


def audit_q(q):
    B = build(q)
    A, _ = find_S4(B, seed=1)
    if A is None:
        raise RuntimeError(f"q={q}: no exterior S4 found")

    reps, idx, cos_of = cosets_right(B, A)
    seen = set()
    rows = []
    for g in B["G"]:
        if g in seen:
            continue
        D = double_coset(B, A, g)
        seen |= D
        Bg = conjugate(B, A, g)
        cap = len(A & Bg)
        if cap == 24:
            continue
        if not generates_full(B, A, Bg):
            continue
        targets = source_targets(B, A, reps, idx, cos_of, g)
        nrm = second_singular_power(targets)
        rows.append((cap, 24 // cap, nrm, g))

    buckets = defaultdict(list)
    for cap, deg, nrm, g in rows:
        buckets[cap].append((nrm, deg, g))
    return len(reps), buckets


def main():
    qs = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [11, 13, 19, 29]
    print("=== true all-pairs double-coset audit on l2(G/A) ===")
    print(f"threshold sqrt(3)/2 = {SQRT3_OVER_2:.6f}\n")
    for q in qs:
        ncos, buckets = audit_q(q)
        print(f"q={q:>3} |G/A|={ncos}")
        global_best = (0.0, None, None)
        for cap in sorted(buckets):
            vals = buckets[cap]
            best = max(vals, key=lambda x: x[0])
            if best[0] > global_best[0]:
                global_best = (best[0], cap, best[1])
            flag = "  **OVER**" if best[0] > SQRT3_OVER_2 + 5e-4 else ""
            print(
                f"  cap={cap:>2} deg={best[1]:>2} "
                f"#double-cosets={len(vals):>3} max_nontriv_sv={best[0]:.6f}{flag}"
            )
        print(
            f"  global max={global_best[0]:.6f} "
            f"(cap={global_best[1]}, deg={global_best[2]})\n"
        )


if __name__ == "__main__":
    main()
