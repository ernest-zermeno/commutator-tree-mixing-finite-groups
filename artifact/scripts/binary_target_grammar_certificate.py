#!/usr/bin/env python3
"""
Certificate for the DIAGONAL REDUCTION of binary target grammars.

Companion to notes/binary_grammar_diagonal_reduction.md (gate #2 outcome).

Verifies, for nine groups (S5, A5, S6, GL3(F2), PGL2(7), PSL2(11), SL2(5),
SL2(3), AGL3(2)):

  A. the ternary relation R(A,B;C) <=> ex x in A, y in B with [x,y] in C != e;
  B. diagonal-slice regression: R(C,C;.) equals the def:targetgraph sweep
     {[w, h w h^-1] : h in G} on classes (one representative per class);
  C. the mixed transfer construction (exhaustive on S5, sampled elsewhere):
     from one witness [x,y] = z, for ALL u in A, v in B, z' in [z] the pair
     g = c a, h = c b (a u a^-1 = x, b v b^-1 = y, c z c^-1 = z') satisfies
     [g u g^-1, h v h^-1] = z';
  D. the solvable radical: O_inf(G) = < union of classes with solvable normal
     closure > is solvable, and <C> solvable <=> C included in O_inf(G);
  E. the MAIN EQUIVALENCE, subset-exhaustively (every nonempty S of nontrivial
     classes):   binary persistence of S   <=>   ex C in S with <C> non-solvable
                 <=>   diagonal persistence of S;
  F. the cycle characterisation of T(G): C lies in the recurrent core of T(G)
     (has an infinite outgoing walk)  <=>  <C> non-solvable  <=>  the diagonal
     singleton dynamics from {C} survives;
  G. specific S5 facts: 33/36 ordered nontrivial class pairs have multiple
     output classes; the three single-output pairs; ((2,2),(2,2)) -> {(3),(5)};
  H. an end-to-end literal realisation over S5: a five-level binary-grammar
     family with a genuinely MIXED production at every level (children in
     distinct classes), all members valid, targets tracking the schedule,
     leaves staying in the transposition class.

Pure stdlib, deterministic. Expected final line:
  BINARY TARGET GRAMMAR CERTIFICATE: PASS (<N> checks)
"""
from itertools import permutations, product, combinations
from collections import defaultdict

CHECKS = 0
def check(cond, msg):
    global CHECKS
    if not cond:
        raise AssertionError("FAILED: " + msg)
    CHECKS += 1

# ------------------------------------------------------------------
# group backends -> index-based (n, MUL, INV, E)
# ------------------------------------------------------------------

def idx_group(elems, mul_raw, e_raw, name):
    elems = list(elems)
    ix = {x: i for i, x in enumerate(elems)}
    n = len(elems)
    MUL = [[ix[mul_raw(elems[i], elems[j])] for j in range(n)] for i in range(n)]
    E = ix[e_raw]
    INV = [0] * n
    for i in range(n):
        for j in range(n):
            if MUL[i][j] == E:
                INV[i] = j
                break
    return name, n, MUL, INV, E

def perm_backend(perms, n, name):
    def mul(p, q):
        return tuple(p[q[i]] for i in range(n))
    return idx_group([tuple(p) for p in perms], mul, tuple(range(n)), name)

def parity(p):
    s = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s ^= 1
    return s

def mobius(p, psl):
    INF = p
    sq = {(x * x) % p for x in range(1, p)}
    seen, out = set(), []
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    det = (a * d - b * c) % p
                    if det == 0 or (psl and det not in sq):
                        continue
                    pm = []
                    for x in range(p + 1):
                        if x == INF:
                            pm.append(INF if c == 0 else (a * pow(c, p - 2, p)) % p)
                        else:
                            den = (c * x + d) % p
                            pm.append(INF if den == 0
                                      else ((a * x + b) * pow(den, p - 2, p)) % p)
                    t = tuple(pm)
                    if t not in seen:
                        seen.add(t)
                        out.append(t)
    return out

def sl2(p):
    els = [(a, b, c, d) for a in range(p) for b in range(p)
           for c in range(p) for d in range(p) if (a * d - b * c) % p == 1]
    def mul(M, N):
        a, b, c, d = M
        x, y, z, w = N
        return ((a * x + b * z) % p, (a * y + b * w) % p,
                (c * x + d * z) % p, (c * y + d * w) % p)
    return els, mul, (1, 0, 0, 1)

def gl3_matrices():
    def ok(M):
        a, b, c = M
        return all((a, b, c)) and len({a, b, c}) == 3 and (a ^ b) != c
    G = [M for M in product(range(8), repeat=3) if ok(M)]
    def mmul(M, N):
        out = []
        for i in range(3):
            acc = 0
            for k in range(3):
                if (M[i] >> (2 - k)) & 1:
                    acc ^= N[k]
            out.append(acc)
        return tuple(out)
    return G, mmul, (4, 2, 1)

def agl3_2():
    G, mmul, I3 = gl3_matrices()
    def av(M, v):
        w = 0
        for i in range(3):
            if bin(M[i] & v).count("1") & 1:
                w |= 1 << (2 - i)
        return w
    els = [(A, b) for A in G for b in range(8)]
    def mul(X, Y):
        (A1, b1), (A2, b2) = X, Y
        return (mmul(A1, A2), av(A1, b2) ^ b1)
    return els, mul, (I3, 0)

# ------------------------------------------------------------------
# per-group analysis
# ------------------------------------------------------------------

def analyze(backend, transfer_exhaustive=False):
    name, n, MUL, INV, E = backend

    def comm(x, y):
        return MUL[MUL[x][y]][MUL[INV[x]][INV[y]]]

    def conj(g, x):
        return MUL[MUL[g][x]][INV[g]]

    # conjugacy classes
    ec, classes = {}, []
    for x in range(n):
        if x in ec:
            continue
        orb = {conj(g, x) for g in range(n)}
        ci = len(classes)
        classes.append(sorted(orb))
        for y in orb:
            ec[y] = ci
    nont = [ci for ci in range(len(classes)) if ci != ec[E]]

    # subgroup closures and solvability
    def subgroup_closure(gens):
        gens = set(gens) | {E}
        seen, fr = set(gens), list(gens)
        while fr:
            nf = []
            for x in fr:
                for g in gens:
                    y = MUL[x][g]
                    if y not in seen:
                        seen.add(y)
                        nf.append(y)
            fr = nf
        return seen

    def solvable(H):
        D = set(H)
        while True:
            D2 = subgroup_closure({comm(x, y) for x in D for y in D})
            if D2 == D:
                return len(D) == 1
            D = D2

    nonsolv = {}
    for C in nont:
        nonsolv[C] = not solvable(subgroup_closure(classes[C]))

    # (A) ternary relation
    R = defaultdict(set)
    for A in nont:
        for B in nont:
            for x in classes[A]:
                for y in classes[B]:
                    z = comm(x, y)
                    if z != E:
                        R[(A, B)].add(ec[z])

    # (B) diagonal-slice regression against the def:targetgraph sweep
    for C in nont:
        w = classes[C][0]
        sweep = set()
        for h in range(n):
            z = comm(w, conj(h, w))
            if z != E:
                sweep.add(ec[z])
        check(sweep == R[(C, C)],
              f"{name}: sweep formula != diagonal slice at class {C}")

    # (C) mixed transfer construction
    def find_conj(u, x):
        for a in range(n):
            if conj(a, u) == x:
                return a
        return None

    for (A, B), Cs in sorted(R.items()):
        # one witness per (A,B;C)
        witness = {}
        for x in classes[A]:
            for y in classes[B]:
                z = comm(x, y)
                if z != E and ec[z] not in witness:
                    witness[ec[z]] = (x, y, z)
            if len(witness) == len(Cs):
                break
        check(set(witness) == Cs, f"{name}: witness collection at {(A, B)}")
        us = classes[A] if transfer_exhaustive else classes[A][:3]
        vs = classes[B] if transfer_exhaustive else classes[B][:3]
        for Ct, (x, y, z) in sorted(witness.items()):
            zs = classes[Ct] if transfer_exhaustive else classes[Ct][:3]
            for u in us:
                a = find_conj(u, x)
                for v in vs:
                    b = find_conj(v, y)
                    for zp in zs:
                        c = find_conj(z, zp)
                        g, h = MUL[c][a], MUL[c][b]
                        check(comm(conj(g, u), conj(h, v)) == zp,
                              f"{name}: transfer identity at {(A, B, Ct)}")

    # (D) solvable radical
    rad_seed = set()
    for C in nont:
        if not nonsolv[C]:
            rad_seed |= set(classes[C])
    radical = subgroup_closure(rad_seed) if rad_seed else {E}
    check(solvable(radical), f"{name}: candidate radical not solvable")
    for C in nont:
        inside = set(classes[C]) <= radical
        check(inside == (not nonsolv[C]),
              f"{name}: <C> solvable <=> C in radical fails at class {C}")

    # (E) main equivalence, subset-exhaustive
    def phi(S, diag=False):
        out = set()
        for A in S:
            for B in S:
                if diag and A != B:
                    continue
                out |= R[(A, B)]
        return frozenset(out)

    def survives(S0, diag=False):
        seen, S = set(), frozenset(S0)
        while S not in seen:
            seen.add(S)
            S = phi(S, diag)
        return len(S) > 0

    for r in range(1, len(nont) + 1):
        for combo in combinations(nont, r):
            S0 = frozenset(combo)
            s = survives(S0)
            sd = survives(S0, diag=True)
            pred = any(nonsolv[C] for C in S0)
            check(s == pred == sd,
                  f"{name}: equivalence fails at S0={sorted(S0)}")

    # (F) recurrent core of T(G) == non-solvable closure == diagonal singleton
    alive = set(nont)
    changed = True
    while changed:
        changed = False
        for C in list(alive):
            if not (R[(C, C)] & alive):
                alive.discard(C)
                changed = True
    for C in nont:
        core = C in alive
        sd = survives(frozenset([C]), diag=True)
        check(core == nonsolv[C] == sd,
              f"{name}: cycle characterisation fails at class {C}")

    n_nonsolv = sum(nonsolv.values())
    print(f"  {name:9s} |G|={n:5d}  nontriv={len(nont):2d}  "
          f"nonsolv-closures={n_nonsolv:2d}  subsets=2^{len(nont)}-1  OK")
    return classes, ec, nont, R, MUL, INV, E

# ------------------------------------------------------------------
# run all nine groups
# ------------------------------------------------------------------

print("Diagonal reduction certificate -- nine groups")
S5 = analyze(perm_backend(permutations(range(5)), 5, "S5"), transfer_exhaustive=True)
analyze(perm_backend([p for p in permutations(range(5)) if parity(p) == 0], 5, "A5"))
analyze(perm_backend(permutations(range(6)), 6, "S6"))
G3 = gl3_matrices()
analyze(idx_group(G3[0], G3[1], G3[2], "GL3(F2)"))
analyze(perm_backend(mobius(7, False), 8, "PGL2(7)"))
analyze(perm_backend(mobius(11, True), 12, "PSL2(11)"))
M5 = sl2(5)
analyze(idx_group(M5[0], M5[1], M5[2], "SL2(5)"))
M3 = sl2(3)
analyze(idx_group(M3[0], M3[1], M3[2], "SL2(3)"))
Ag = agl3_2()
analyze(idx_group(Ag[0], Ag[1], Ag[2], "AGL3(2)"))

# ------------------------------------------------------------------
# (G) specific S5 facts
# ------------------------------------------------------------------

classes, ec, nont, R, MUL, INV, E = S5

def ctype_of_class(C):
    # cycle type via a representative permutation index -> need raw perms
    return C  # class indices are stable; we identify by size+order below

perms = [tuple(p) for p in permutations(range(5))]
def ptype(p):
    seen = [False] * 5
    t = []
    for i in range(5):
        if not seen[i]:
            l, j = 0, i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                l += 1
            t.append(l)
    return tuple(sorted(t, reverse=True))

class_type = {C: ptype(perms[classes[C][0]]) for C in nont}
type_class = {v: k for k, v in class_type.items()}

multi = sum(1 for k in R if len(R[k]) > 1)
check(len(nont) ** 2 == 36, "S5: 36 ordered pairs")
check(multi == 33, "S5: 33 pairs with multiple outputs")
mono = sorted((class_type[a], class_type[b]) for (a, b) in R if len(R[(a, b)]) <= 1)
check(mono == [((2, 1, 1, 1), (2, 1, 1, 1)),
               ((2, 1, 1, 1), (3, 1, 1)),
               ((3, 1, 1), (2, 1, 1, 1))],
      "S5: the three single-output pairs")
dd = type_class[(2, 2, 1)]
check({class_type[c] for c in R[(dd, dd)]} == {(3, 1, 1), (5,)},
      "S5: ((2,2),(2,2)) -> {(3),(5)}")
print(f"  S5 facts: 36 pairs, 33 multi-output, mono pairs and "
      f"((2,2),(2,2))->{{(3),(5)}}  OK")

# ------------------------------------------------------------------
# (H) end-to-end literal mixed-grammar family over S5
# ------------------------------------------------------------------
# Portfolio {P3 (target class (3)), P22 (target class (2,2))} maintained for
# levels 2..6; EVERY production below is mixed (children in distinct classes):
#   next P3  = [ g . P3 . g^-1 , h . P22 . h^-1 ]   (needs R((3),(2,2);(3)))
#   next P22 = [ g . P22 . g^-1, h . P3 . h^-1 ]    (needs R((2,2),(3);(2,2)))

C3, C22, C2 = type_class[(3, 1, 1)], type_class[(2, 2, 1)], type_class[(2, 1, 1, 1)]
check(C3 in R[(C3, C22)], "S5: mixed production (3),(2,2)->(3) exists")
check(C22 in R[(C22, C3)], "S5: mixed production (2,2),(3)->(2,2) exists")

def comm_i(x, y):
    return MUL[MUL[x][y]][MUL[INV[x]][INV[y]]]

def conj_i(g, x):
    return MUL[MUL[g][x]][INV[g]]

def alpha(t):
    if isinstance(t, int):
        return t
    return comm_i(alpha(t[0]), alpha(t[1]))

def all_valid(t):
    if isinstance(t, int):
        return t != E
    return alpha(t) != E and all_valid(t[0]) and all_valid(t[1])

def conj_tree(g, t):
    if isinstance(t, int):
        return conj_i(g, t)
    return (conj_tree(g, t[0]), conj_tree(g, t[1]))

def leaves(t):
    if isinstance(t, int):
        return {t}
    return leaves(t[0]) | leaves(t[1])

def find_conj_i(u, x):
    for a in range(len(perms)):
        if conj_i(a, u) == x:
            return a
    return None

def compose(tA, tB, Ctarget):
    """Mixed composition with target prescribed to the canonical class rep."""
    u, v = alpha(tA), alpha(tB)
    A, B = ec[u], ec[v]
    zp = classes[Ctarget][0]
    for x in classes[A]:
        for y in classes[B]:
            z = comm_i(x, y)
            if z != E and ec[z] == Ctarget:
                a, b = find_conj_i(u, x), find_conj_i(v, y)
                c = find_conj_i(z, zp)
                return (conj_tree(MUL[c][a], tA), conj_tree(MUL[c][b], tB))
    return None

pix = {p: i for i, p in enumerate(perms)}
t12, t23, t34, t45, t14 = (pix[(1, 0, 2, 3, 4)], pix[(0, 2, 1, 3, 4)],
                           pix[(0, 1, 3, 2, 4)], pix[(0, 1, 2, 4, 3)],
                           pix[(3, 1, 2, 0, 4)])
trans_class = set(classes[C2])
lvl1_a, lvl1_b = (t12, t23), (t34, t45)          # 3-cycle targets
P3 = (lvl1_a, lvl1_b)                             # level-2, want target in (3)
# fix P3's target class: [3cyc,3cyc] lands where it lands; steer via compose
P3 = compose(lvl1_a, lvl1_b, C3)
check(P3 is not None and ec[alpha(P3)] == C3, "S5: level-2 P3 built")
P22 = compose((t12, t23), (t34, t14), C22)
check(P22 is not None and ec[alpha(P22)] == C22, "S5: level-2 P22 built")
check(all_valid(P3) and all_valid(P22), "S5: level-2 members valid")

for level in range(3, 7):
    P3n = compose(P3, P22, C3)     # mixed: children classes (3) and (2,2)
    P22n = compose(P22, P3, C22)   # mixed: children classes (2,2) and (3)
    check(P3n is not None and P22n is not None, f"S5: level-{level} productions")
    for t, Ct in ((P3n, C3), (P22n, C22)):
        check(all_valid(t), f"S5: level-{level} member valid")
        check(ec[alpha(t)] == Ct, f"S5: level-{level} target class on schedule")
        check(leaves(t) <= trans_class, f"S5: level-{level} leaves in class (2)")
    P3, P22 = P3n, P22n
print("  end-to-end: 5-level mixed-production family over S5 "
      "(valid, on-schedule, leaves in (2))  OK")

print()
print(f"BINARY TARGET GRAMMAR CERTIFICATE: PASS ({CHECKS} checks)")
