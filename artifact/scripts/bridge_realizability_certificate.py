#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 4A certificate: rigid realizability over the transposition class of S5.

Backs (never replaces) the primary mathematical proofs in
notes/rigid_realizability.md:

  PART 1  Target-generation dichotomy at level 2 (balanced valid configs over
          S5 transpositions): census 3,600 = 1,080 + 2,160 + 360, generating
          <=> 3-cycle target, trapped <=> (2,2) target with H = S4 (order 24);
          analytic identities 20*6*9 / 20*12*9 / 20*2*9; the 8 canonical
          orientation cases.
  PART 2  Exact-level-3 self-similar family: seed = the s4xs4 configuration
          (two valid S4-trapped children, generating root with 5-cycle
          target rho) and FIXED conjugators (g,h) with [g rho g^-1,
          h rho h^-1] = rho; tree-level verification at levels 3..6.
  PART 3  Threshold-2 family: G1 witness (S3), generating level-2 seed with
          3-cycle target delta, and the DIRECT self-similar 3-cycle loop
          [gd delta gd^-1, hd delta hd^-1] = delta from level 2 (primary
          construction; tree levels 2..6), plus an independent gate pair.
          The two-step promotion (3)->(2,2)->(5)->(5) is kept as an OPTIONAL
          ALTERNATIVE construction (tree levels 2..6) - it is NOT necessary.
  PART 4  Level 0/1 leaf-subgroup facts for the {2,3} classification (every
          valid level-<=1 configuration has proper leaf subgroup).

Pure python, exhaustive where stated, deterministic (no sampling, no seeds;
conjugator constants below were derived by lexicographically-first search in
derivation mode and then pinned). Exit 0 = all invariants hold; exit 1 = any
failure. Scope: tau-level (group shadows) only; no support/C/D claims.
"""
import sys
from itertools import combinations, permutations, product

N = 5
E = tuple(range(N))

def mul(a, b):
    return tuple(a[b[i]] for i in range(N))

def inv(a):
    r = [0] * N
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)

def comm(a, b):                     # P1 convention [a,b] = a b a^-1 b^-1
    return mul(mul(a, b), mul(inv(a), inv(b)))

def conj(g, x):
    return mul(mul(g, x), inv(g))

def tp(i, j):
    p = list(range(N))
    p[i], p[j] = p[j], p[i]
    return tuple(p)

def closure(gens):
    S = {E}
    fr = [E]
    while fr:
        x = fr.pop()
        for s in gens:
            y = mul(x, s)
            if y not in S:
                S.add(y)
                fr.append(y)
    return S

def ctype(p):
    seen = set()
    out = []
    for i in range(N):
        if i in seen:
            continue
        c = 0
        j = i
        while j not in seen:
            seen.add(j)
            j = p[j]
            c += 1
        if c > 1:
            out.append(c)
    return tuple(sorted(out, reverse=True))

def cyc(*pts):                      # 1-indexed cycle -> permutation tuple
    r = list(range(N))
    for i in range(len(pts)):
        r[pts[i] - 1] = pts[(i + 1) % len(pts)] - 1
    return tuple(r)

FAILURES = []

def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (tag, name, ("  " + str(detail)) if detail else ""), flush=True)
    if not cond:
        FAILURES.append(name)

# --------------------------------------------------------------------------
# Pinned constants (derived by lexicographically-first search, then pinned).
# Set to None to re-derive.
# --------------------------------------------------------------------------
PIN = {
    # PART 2: [g rho g^-1, h rho h^-1] = rho for the s4xs4 root 5-cycle rho
    "rho": (2, 3, 1, 4, 0),
    "g": (0, 1, 3, 4, 2), "h": (0, 2, 3, 1, 4),
    # PART 3 primary: DIRECT 3-cycle loop [gd d gd^-1, hd d hd^-1] = delta
    "delta": (0, 4, 1, 3, 2),
    "gd": (0, 1, 2, 4, 3), "hd": (1, 0, 4, 3, 2),   # = the reviewer's pair
    # PART 3 optional alternative: (3) -> (2,2) -> (5) -> stable
    "g2": (0, 1, 2, 3, 4), "h2": (0, 1, 2, 4, 3), "eps": (0, 2, 1, 4, 3),
    "g3": (0, 1, 2, 3, 4), "h3": (1, 0, 3, 2, 4), "rho2": (1, 3, 0, 4, 2),
    "g4": (0, 2, 3, 1, 4), "h4": (0, 1, 4, 2, 3),
}

def first_pair(target_test, base):
    """Lexicographically-first (g,h) in S5 x S5 with target_test([g b g^-1, h b h^-1])."""
    perms = sorted(permutations(range(N)))
    for g in perms:
        cg = conj(g, base)
        for h in perms:
            w = comm(cg, conj(h, base))
            if target_test(w):
                return g, h, w
    return None

def main():
    T = [tp(i, j) for i, j in combinations(range(N), 2)]
    S5_order = 120

    # ---------------- PART 1: dichotomy ----------------
    print("PART 1 - target-generation dichotomy at level 2 (EXHAUSTIVE)", flush=True)
    gates = [(s, t) for s in T for t in T if s != t and comm(s, t) != E]
    check("valid ordered L1 gates = 60", len(gates) == 60, len(gates))
    tcyc = {g: comm(*g) for g in gates}
    check("every valid gate value is a 3-cycle",
          all(ctype(c) == (3,) for c in tcyc.values()))
    three_cycles = {c for c in tcyc.values()}
    check("20 distinct 3-cycles, 3 ordered realisations each",
          len(three_cycles) == 20 and
          all(sum(1 for c in tcyc.values() if c == x) == 3 for x in three_cycles))

    # literal traversal of the FULL 10^4 leaf universe (no factorisation)
    stats = {"gen": 0, "trap": 0, "invalid": 0, "gate_invalid": 0}
    ok_gen = ok_trap = ok_traporder = True
    gen_orbits = {}
    trap_orbits = {}
    for a, b, c, d in product(T, repeat=4):
        c1, c2 = comm(a, b), comm(c, d)
        if c1 == E or c2 == E:
            stats["gate_invalid"] += 1
            continue
        root = comm(c1, c2)
        if root == E:
            stats["invalid"] += 1
            continue
        H = closure([a, b, c, d])
        canon = min(((p, q) if p <= q else (q, p))
                    for p in ((a, b), (b, a)) for q in ((c, d), (d, c)))
        if len(H) == S5_order:
            stats["gen"] += 1
            ok_gen &= (ctype(root) == (3,))
            gen_orbits[canon] = gen_orbits.get(canon, 0) + 1
        else:
            stats["trap"] += 1
            ok_trap &= (ctype(root) == (2, 2))
            ok_traporder &= (len(H) == 24)
            trap_orbits[canon] = trap_orbits.get(canon, 0) + 1
    check("10^4 universe accounting: 10000 = 3600 valid-gates + 6400 discarded",
          stats["gate_invalid"] == 6400 and
          sum(v for k, v in stats.items() if k != "gate_invalid") == 3600, stats)
    check("census 3600 = 1080 gen + 2160 trapped + 360 root-invalid",
          (stats["gen"], stats["trap"], stats["invalid"]) == (1080, 2160, 360))
    check("ALL generating have 3-cycle target", ok_gen)
    check("ALL trapped have (2,2) target", ok_trap)
    check("ALL trapped leaf subgroups are full S4 (order 24)", ok_traporder)
    check("reduced census: 135 generating orbits, all free of size 8",
          len(gen_orbits) == 135 and set(gen_orbits.values()) == {8},
          len(gen_orbits))
    check("reduced census: 270 trapped orbits, all free of size 8",
          len(trap_orbits) == 270 and set(trap_orbits.values()) == {8},
          len(trap_orbits))
    check("analytic identities 20*6*9 / 20*12*9 / 20*2*9",
          (20 * 6 * 9, 20 * 12 * 9, 20 * 2 * 9) == (1080, 2160, 360))

    # canonical orientation cases (support intersection 1 and 2)
    a, b = cyc(1, 2, 3), cyc(1, 4, 5)
    case1 = [ctype(comm(x, y)) for x in (a, inv(a)) for y in (b, inv(b))]
    check("canonical cases |supp cap|=1 (4 cases, 8 total) -> all 3-cycles",
          set(case1) == {(3,)}, case1)
    b2 = cyc(1, 2, 4)
    case2 = [ctype(comm(x, y)) for x in (a, inv(a)) for y in (b2, inv(b2))]
    check("canonical cases |supp cap|=2 (4 cases, 8 total) -> all (2,2)",
          set(case2) == {(2, 2)}, case2)

    # ---------------- tree machinery ----------------
    def alpha_of(leaves):
        lv = list(leaves)
        while len(lv) > 1:
            lv = [comm(lv[i], lv[i + 1]) for i in range(0, len(lv), 2)]
        return lv[0]

    def step(leaves, g, h):
        return ([conj(g, l) for l in leaves] + [conj(h, l) for l in leaves])

    def family_check(tag, leaves_seed, k_seed, target_seq, conj_seq, kmax=6):
        """target_seq[i] = expected alpha at level k_seed+i;
        conj_seq[i] = (g,h) used to go from level k_seed+i to k_seed+i+1."""
        lv = list(leaves_seed)
        k = k_seed
        i = 0
        while True:
            al = alpha_of(lv)
            H = closure(lv)
            check("%s: level %d alpha == expected, type %s" % (tag, k, ctype(al)),
                  al == target_seq[i])
            check("%s: level %d generates S5, leaves are transpositions, count 2^%d"
                  % (tag, k, k),
                  len(H) == S5_order and all(ctype(l) == (2,) for l in lv)
                  and len(lv) == 2 ** k)
            if k == kmax:
                break
            g, h = conj_seq[i]
            lv = step(lv, g, h)
            k += 1
            i += 1

    # ---------------- PART 2: exact-3 self-similar family ----------------
    print("PART 2 - exact-level-3 self-similar family", flush=True)
    A = [tp(0, 1), tp(1, 2), tp(2, 3), tp(0, 3)]     # cherry pair on {1..4}
    B = [tp(1, 2), tp(2, 4), tp(3, 4), tp(1, 3)]     # cherry pair on {2..5}
    aA, aB = alpha_of(A), alpha_of(B)
    HA, HB = closure(A), closure(B)
    check("child tau_A: valid, H = S4 (order 24), (2,2) target",
          aA != E and len(HA) == 24 and ctype(aA) == (2, 2))
    check("child tau_B: valid, H = S4 (order 24), (2,2) target",
          aB != E and len(HB) == 24 and ctype(aB) == (2, 2))
    seed3 = A + B
    rho = alpha_of(seed3)
    check("seed tau_3: root target is a 5-cycle, H = S5",
          ctype(rho) == (5,) and len(closure(seed3)) == S5_order, rho)
    if PIN["rho"] is not None:
        check("pinned rho matches", rho == PIN["rho"], rho)
    found = first_pair(lambda w: w == rho, rho)
    if PIN["g"] is None:
        g, h, _ = found
        print("  [derivation] g = %s  h = %s" % (g, h), flush=True)
    else:
        g, h = PIN["g"], PIN["h"]
        check("pinned (g,h) IS the lexicographically-first pair for rho",
              found == (g, h, rho), found[:2])
    check("fixed conjugators: [g rho g^-1, h rho h^-1] == rho",
          comm(conj(g, rho), conj(h, rho)) == rho, (g, h))
    family_check("exact3", seed3, 3, [rho] * 4, [(g, h)] * 3, kmax=6)

    # ---------------- PART 3: threshold-2 family ----------------
    print("PART 3 - threshold-2 family: direct 3-cycle loop (primary)", flush=True)
    w1 = [tp(0, 1), tp(1, 2)]
    check("G1 witness: valid, H = S3 (order 6) proper",
          alpha_of(w1) != E and len(closure(w1)) == 6)
    seed2 = [tp(0, 1), tp(1, 2), tp(2, 3), tp(3, 4)]
    delta = alpha_of(seed2)
    check("level-2 seed: generates S5, 3-cycle target delta",
          len(closure(seed2)) == S5_order and ctype(delta) == (3,), delta)
    if PIN["delta"] is not None:
        check("pinned delta matches", delta == PIN["delta"], delta)
    # PRIMARY: direct self-similar 3-cycle loop from level 2
    found = first_pair(lambda w: w == delta, delta)
    if PIN["gd"] is None:
        gd, hd, _ = found
        print("  [derivation] gd = %s  hd = %s" % (gd, hd), flush=True)
    else:
        gd, hd = PIN["gd"], PIN["hd"]
        check("pinned (gd,hd) IS the lexicographically-first pair for delta",
              found == (gd, hd, delta), found[:2])
    check("DIRECT loop: [gd delta gd^-1, hd delta hd^-1] == delta",
          comm(conj(gd, delta), conj(hd, delta)) == delta, (gd, hd))
    family_check("thr2-direct", seed2, 2, [delta] * 5, [(gd, hd)] * 4, kmax=6)

    print("PART 3b - optional ALTERNATIVE promotion (3)->(2,2)->(5) "
          "(not necessary)", flush=True)
    found = first_pair(lambda w: ctype(w) == (2, 2), delta)
    if PIN["g2"] is None:
        g2, h2, eps = found
        print("  [derivation] g2 = %s  h2 = %s  eps = %s" % (g2, h2, eps), flush=True)
    else:
        g2, h2, eps = PIN["g2"], PIN["h2"], PIN["eps"]
        check("pinned (g2,h2,eps) IS the lexicographically-first (2,2) promotion",
              found == (g2, h2, eps))
    check("promotion (3)->(2,2): [g2 d g2^-1, h2 d h2^-1] = eps of type (2,2)",
          comm(conj(g2, delta), conj(h2, delta)) == eps and ctype(eps) == (2, 2))
    found = first_pair(lambda w: ctype(w) == (5,), eps)
    if PIN["g3"] is None:
        g3, h3, rho2 = found
        print("  [derivation] g3 = %s  h3 = %s  rho2 = %s" % (g3, h3, rho2), flush=True)
    else:
        g3, h3, rho2 = PIN["g3"], PIN["h3"], PIN["rho2"]
        check("pinned (g3,h3,rho2) IS the lexicographically-first (5,) promotion",
              found == (g3, h3, rho2))
    check("promotion (2,2)->(5): [g3 e g3^-1, h3 e h3^-1] = rho2 of type (5,)",
          comm(conj(g3, eps), conj(h3, eps)) == rho2 and ctype(rho2) == (5,))
    found = first_pair(lambda w: w == rho2, rho2)
    if PIN["g4"] is None:
        g4, h4, _ = found
        print("  [derivation] g4 = %s  h4 = %s" % (g4, h4), flush=True)
    else:
        g4, h4 = PIN["g4"], PIN["h4"]
        check("pinned (g4,h4) IS the lexicographically-first pair for rho2",
              found == (g4, h4, rho2))
    check("stable rule: [g4 rho2 g4^-1, h4 rho2 h4^-1] == rho2",
          comm(conj(g4, rho2), conj(h4, rho2)) == rho2, (g4, h4))
    family_check("thr2-alt", seed2, 2, [delta, eps, rho2, rho2, rho2],
                 [(g2, h2), (g3, h3), (g4, h4), (g4, h4)], kmax=6)

    # ---------------- PART 4: level 0/1 leaf-subgroup facts ----------------
    print("PART 4 - level 0/1 leaf-subgroup facts for the {2,3} classification",
          flush=True)
    check("level 0: every leaf is valid with proper H (order 2)",
          all(s != E and len(closure([s])) == 2 for s in T))
    check("level 1: all 60 valid pairs have H = S3 (order 6) proper",
          all(len(closure(list(p))) == 6 for p in gates))

    print()
    if FAILURES:
        print("RIGID REALIZABILITY CERTIFICATE: FAIL (%d)" % len(FAILURES), flush=True)
        return 1
    print("RIGID REALIZABILITY CERTIFICATE: PASS", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
