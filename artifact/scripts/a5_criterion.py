"""
Exact group/count certificate and numerical SVD cross-check for cor:A5.
The SVD values below are float64 estimates, not interval certificates.

For G = A_5 (order 60) with C = the involution class (the 15 double transpositions):
  (1) closure lattice (EXHAUSTIVE over NONEMPTY subsets): the closures <S> for
      nonempty S subset C are exactly 37 subgroups = 36 proper nontrivial ones
      (Counter: C_2 x15, V_4 x5, S_3 x10, D_10 x6; derived lengths C_2=1, V_4=1,
      S_3=2, D_10=2) plus A_5 itself; max derived length over proper closures = 2.
      The trivial subgroup arises only as the closure of the EMPTY subset (dl = 0)
      and is deliberately omitted -- it cannot affect the max. In particular A_4
      is NOT C-generated -- its involutions close to V_4 -- so it does not enter
      the max. Hence D = 2 by thm:dichotomy(2).
  (2) L1 (EXHAUSTIVE, regular rep): ORDERED universe = 180 valid ordered pairs,
      REDUCED = 90 representatives (free orbits of the swap); ALL 180 have 2nd
      singular value 1 (two involutions generate a dihedral, hence proper,
      subgroup) -- level 1 never contracts.
  (3) L2 (EXHAUSTIVE, regular rep): ORDERED universe = 29,640 valid ordered
      quadruples, REDUCED = 3,705 representatives (free orbits of the order-8
      symmetry group: swap inside each gate, swap of gates); ALL 29,640 generate
      A_5, 0 with 2nd-sv = 1; bottleneck (max 2nd-sv) = BOTTLENECK_REF < 1
      =>  exact level 2.

Certificate discipline: no sampling anywhere; every count is asserted for the
ORDERED universe and for the REDUCED representatives separately, orbit-freeness
is verified (|ordered| = orbit_size * |reduced| with per-representative preimage
counts), and the script exits nonzero if any invariant fails.
"""
import sys
import numpy as np
from itertools import combinations, permutations

# Reference bottleneck (max 2nd singular value over the L2 ordered universe) and
# tolerance. Justification of TOL: entries of the 60x60 operators are exact
# dyadic rationals; float64 SVD backward error is ~ n*eps*||A|| ~ 60*2.2e-16
# ~ 1.3e-14, and cross-platform LAPACK variation is orders of magnitude below
# 1e-9. TOL = 1e-9 is therefore generous, while the certificate's decision
# margins are ~0.3464 (bottleneck to 0) and ~0.6536 (bottleneck to 1) -- both
# many orders larger than the SVD error.
BOTTLENECK_REF = 0.346367621209230   # derived 2026-07-11 by this script (derivation mode)
TOL = 1e-9

def pmul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def pinv(p):
    r = [0]*len(p)
    for i, j in enumerate(p): r[j] = i
    return tuple(r)
def is_id(p): return all(p[i] == i for i in range(len(p)))
def comm(a, b): return pmul(pmul(a, b), pmul(pinv(a), pinv(b)))
def involutions(perms): return [p for p in perms if not is_id(p) and is_id(pmul(p, p))]

def a5():
    def sign(p):
        s = 1
        for i, j in combinations(range(len(p)), 2):
            if p[i] > p[j]: s = -s
        return s
    return [p for p in permutations(range(5)) if sign(p) == 1]

def closure(gens, n=5):
    ID = tuple(range(n))
    seen = {ID}; fr = [ID]
    while fr:
        x = fr.pop()
        for s in gens:
            y = pmul(x, s)
            if y not in seen: seen.add(y); fr.append(y)
    return frozenset(seen)

def derived_length(H):
    """dl of a subgroup given as a frozenset; None if the series stalls (perfect)."""
    dl = 0
    while len(H) > 1:
        K = closure([comm(a, b) for a in H for b in H if not is_id(comm(a, b))] or [tuple(range(5))])
        if K == H: return None   # perfect, non-solvable
        H = K; dl += 1
    return dl

def iso_name(H):
    n = len(H); ab = all(pmul(a, b) == pmul(b, a) for a in H for b in H)
    return {2: "C2", 4: "V4" if ab else "?4", 6: "C6" if ab else "S3",
            10: "C10" if ab else "D10", 12: "A4", 60: "A5"}.get(n, f"?{n}")

def cgen_subgroups(invs):
    """EXHAUSTIVE lattice of C-generated subgroups: <S> for every nonempty S subset C.
    BFS adding one involution at a time reaches every such closure (any <S> is a
    chain <s1> < <s1,s2> < ... of one-at-a-time closures)."""
    seen = set(); frontier = {closure([s]) for s in invs}
    while frontier:
        nxt = set()
        for H in frontier:
            if H in seen: continue
            seen.add(H)
            for t in invs:
                if t not in H:
                    K = closure(list(H) + [t])
                    if K not in seen: nxt.add(K)
        frontier = nxt
    return seen

def analyze(perms):
    N = len(perms); gidx = {g: i for i, g in enumerate(perms)}
    invs = involutions(perms)
    print(f"A_5: |G|={N}, involutions={len(invs)}", flush=True)

    # ---- (1) exhaustive C-generated closure lattice ----
    subs = cgen_subgroups(invs)
    proper = sorted((H for H in subs if len(H) < N), key=len)
    from collections import Counter
    types = Counter(iso_name(H) for H in proper)
    dls = {iso_name(H): derived_length(H) for H in proper}
    maxdl = max(dls.values())
    print(f"\n  (1) C-generated closures (EXHAUSTIVE over NONEMPTY subsets): "
          f"{len(subs)} total", flush=True)
    print(f"      proper nontrivial closures: {dict(sorted(types.items(), key=lambda kv: kv[0]))}  "
          f"(dl: {dict(sorted(dls.items()))})", flush=True)
    a4 = frozenset(p for p in perms if p[4] == 4)
    a4cl = closure(involutions(list(a4)))
    print(f"      A_4 (point stabiliser, order {len(a4)}): involution closure = "
          f"{iso_name(a4cl)} (order {len(a4cl)})  => A_4 NOT C-generated", flush=True)
    # exact multiplicities and derived lengths, pinned (not just the name set):
    assert dict(types) == {"C2": 15, "V4": 5, "S3": 10, "D10": 6}, \
        f"proper nontrivial closure multiplicities changed: {dict(types)}"
    assert dls == {"C2": 1, "V4": 1, "S3": 2, "D10": 2}, \
        f"derived lengths changed: {dls}"
    # 37 closures of nonempty subsets = 36 proper nontrivial + A_5 itself;
    # the trivial subgroup (empty subset, dl 0) is outside this enumeration.
    assert len(subs) == 37 and len(proper) == 36, (len(subs), len(proper))
    assert sum(1 for H in subs if len(H) == N) == 1, "A_5 itself must appear exactly once"
    assert all(len(H) > 1 for H in subs), "trivial subgroup must not arise from nonempty subsets"
    assert maxdl == 2
    assert len(a4cl) == 4 and iso_name(a4cl) == "V4"
    print(f"      => proper nontrivial C-generated = {{C2 x15, V4 x5, S3 x10, D10 x6}}, "
          f"max dl = {maxdl}  => D = 2 by thm:dichotomy(2)", flush=True)

    # ---- regular-rep operators (2nd-sv < 1 <=> <leaves> = G, cf. generates_G.py) ----
    def regproj(s):
        P = np.zeros((N, N))
        for c, g in enumerate(perms): P[gidx[pmul(s, g)], c] = 1.0
        return 0.5*(P + np.eye(N))
    Preg = {s: regproj(s) for s in invs}
    def sv2(M): return float(np.linalg.svd(M, compute_uv=False)[1])
    def blk(s, t): return np.linalg.matrix_power(Preg[s] @ Preg[t], 2)

    # ---- (2) level 1 ----
    # ORDERED universe: (s,t), s != t, [s,t] != e.  REDUCED: unordered {s,t}
    # (the swap (s,t)->(t,s) preserves validity since [t,s]=[s,t]^{-1}, and is
    # free: a fixed point would need s=t, excluded).
    opairs = [(s, t) for s in invs for t in invs if s != t and not is_id(comm(s, t))]
    pairs = [(s, t) for s, t in combinations(invs, 2) if not is_id(comm(s, t))]
    assert len(opairs) == 180, f"L1 ORDERED = {len(opairs)} != 180"
    assert len(pairs) == 90, f"L1 REDUCED = {len(pairs)} != 90"
    # orbit-freeness, per representative: each unordered pair has exactly 2 ordered forms
    from collections import Counter as _C
    canon1 = _C(frozenset(pr) for pr in opairs)
    assert set(canon1) == {frozenset(pr) for pr in pairs} and set(canon1.values()) == {2}, \
        "L1 orbits are not free of size 2"
    D1o = {pr: blk(*pr) for pr in opairs}          # all 180 ordered blocks
    sv1 = {pr: sv2(D1o[pr]) for pr in opairs}
    assert max(sv1.values()) <= 1 + TOL, "singular value above 1: not a contraction"
    l1_one_o = sum(1 for pr in opairs if sv1[pr] > 1 - TOL)
    l1_one = sum(1 for pr in pairs if sv1[pr] > 1 - TOL)
    print(f"\n  (2) L1 (EXHAUSTIVE): ORDERED universe {len(opairs)} valid pairs, "
          f"REDUCED {len(pairs)} representatives (free swap orbits, 180 = 2 x 90)", flush=True)
    print(f"      2nd-sv=1 in {l1_one_o}/{len(opairs)} ORDERED ({l1_one}/{len(pairs)} REDUCED)"
          f"  => level 1 never contracts", flush=True)
    assert l1_one_o == 180 and l1_one == 90

    # ---- (3) level 2 ----
    # ORDERED universe: (p,q) with p,q valid ORDERED pairs and [c_p,c_q] != e.
    # REDUCED: orbits of the order-8 group <swap p, swap q, swap (p,q)>; validity
    # is invariant (an element commutes with c iff with c^{-1}; [c2,c1]=[c1,c2]^{-1})
    # and orbits are free (any fixed point forces a gate with alpha = e, excluded).
    ptgt = {pr: comm(*pr) for pr in opairs}
    oquads = [(p, q) for p in opairs for q in opairs if not is_id(comm(ptgt[p], ptgt[q]))]
    quads = [(i, j) for i in range(len(pairs)) for j in range(i+1, len(pairs))
             if not is_id(comm(ptgt[pairs[i]], ptgt[pairs[j]]))]
    assert len(oquads) == 29640, f"L2 ORDERED = {len(oquads)} != 29640"
    assert len(quads) == 3705, f"L2 REDUCED = {len(quads)} != 3705"
    canon2 = _C(frozenset((frozenset(p), frozenset(q))) for (p, q) in oquads)
    assert set(canon2) == {frozenset((frozenset(pairs[i]), frozenset(pairs[j])))
                           for (i, j) in quads} and set(canon2.values()) == {8}, \
        "L2 orbits are not free of size 8 (29640 = 8 x 3705)"
    l2_one_o = n_genG_o = 0; max2 = 0.0
    for (p, q) in oquads:
        (s1, t1), (s2, t2) = p, q
        D = D1o[p] @ D1o[q] @ D1o[(t1, s1)] @ D1o[(t2, s2)]
        o = sv2(D); max2 = max(max2, o)
        assert o <= 1 + TOL, "singular value above 1: not a contraction"
        if o > 1 - TOL: l2_one_o += 1
        if len(closure([s1, t1, s2, t2])) == N: n_genG_o += 1
    print(f"  (3) L2 (EXHAUSTIVE): ORDERED universe {len(oquads)} valid quads, "
          f"REDUCED {len(quads)} representatives (free order-8 orbits, 29640 = 8 x 3705)", flush=True)
    print(f"      generate A_5 in {n_genG_o}/{len(oquads)} ORDERED, 2nd-sv=1 in {l2_one_o}, "
          f"bottleneck (max 2nd-sv over ORDERED universe) = {max2:.15f}", flush=True)
    assert n_genG_o == 29640 and l2_one_o == 0 and max2 < 1
    if BOTTLENECK_REF is None:
        print(f"      [derivation mode] pin BOTTLENECK_REF = {max2:.15f}", flush=True)
    else:
        assert abs(max2 - BOTTLENECK_REF) < TOL, \
            f"bottleneck {max2:.15f} != pinned {BOTTLENECK_REF:.15f} (+-{TOL})"

    ok = (maxdl == 2 and l1_one_o == len(opairs) and l2_one_o == 0)
    print(f"\n  => exact level = 2 (L1 all fail, L2 all contract, all irreps; "
          f"EXHAUSTIVE over the ORDERED universes)   "
          f"(predicted D = max dl = 2)   {'MATCH' if ok else 'CHECK'}", flush=True)
    assert ok
    print("\nA5 CERTIFICATE: PASS", flush=True)

if __name__ == "__main__":
    # Any failed assert raises AssertionError -> nonzero exit. Explicit 0 on success.
    G = a5()
    assert len(G) == 60, len(G)
    analyze(G)
    sys.exit(0)
