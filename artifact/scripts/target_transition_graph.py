#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 4A.2 certificate: the target transition graph T(G).

For a finite group G, the SUCCESSOR relation on nontrivial conjugacy
classes is
    [w] -> [w']   iff   w' = [u, v] != e  for some  u, v  in the class [w]
(equivalently, w' = [g w g^-1, h w h^-1] for some g, h in G). The relation
is class-well-defined, every literal element of [w'] is realisable
(transfer lemma), and a class SELF-LOOP is equivalent to the existence of a
FIXED pair (g, h) with [g w g^-1, h w h^-1] = w literally (loop lemma) -
the engine behind the self-similar families of Phase 4A.

Computed instances (each class swept by [w, h w h^-1] over all h in G,
which covers all successor classes since [g w g^-1, h w h^-1] is conjugate
to [w, (g^-1 h) w (g^-1 h)^-1]):

  PART 1  S5 (6 nontrivial classes, labelled by cycle type) - pinned
          adjacency, pinned loop set, and the Phase-4A re-derivations:
          successors of (2) = {(3)} (the gate lemma), successors of (3) =
          {(3),(2,2)} (in particular NO direct (3)->(5) edge), and the
          reachability facts behind the {2,3} threshold classification.
  PART 2  GL3(F2) (order 168) - small finite Lie-type instance; classes
          labelled (element order, class size); transvection class
          highlighted.
  PART 3  PGL2(7) (order 336) - small finite Lie-type instance via the
          shared core builder.

Pure python + bridge_sweep_core (imported read-only, never modified).
Deterministic, no sampling, no seeds; pinned constants derived by
lexicographically-first search and re-verified on every run.
Exit 0 = all invariants hold; exit 1 = any failure. tau-level target
combinatorics ONLY: no mixing/contraction claims, no support/C/D claims.
"""
import sys
from itertools import combinations, permutations

import bridge_sweep_core as core

FAILURES = []

def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (tag, name, ("  " + str(detail)) if detail else ""),
          flush=True)
    if not cond:
        FAILURES.append(name)

# --------------------------------------------------------------------------
# Generic machinery over core.Group (table-backed; works for any small G)
# --------------------------------------------------------------------------

def conj_idx(G, g, x):
    return int(G.tbl[G.tbl[g, x], G.inv[g]])

def comm_idx(G, a, b):
    return int(G.tbl[G.tbl[G.tbl[a, b], G.inv[a]], G.inv[b]])

def element_order(G, x):
    o, y = 1, x
    while y != G.e:
        y = int(G.tbl[y, x])
        o += 1
    return o

def rank_g_minus_I_f2(e):
    """rank of (g - I) over F2, for a 9-bit row-major 3x3 matrix tuple."""
    rows = [list(e[0:3]), list(e[3:6]), list(e[6:9])]
    for r in range(3):
        rows[r][r] ^= 1
    rk = 0
    for c in range(3):
        piv = next((r for r in range(rk, 3) if rows[r][c]), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for r in range(3):
            if r != rk and rows[r][c]:
                rows[r] = [a ^ b for a, b in zip(rows[r], rows[rk])]
        rk += 1
    return rk

def conjugacy_classes(G):
    """All conjugacy classes (sorted index lists), identity class first."""
    seen = set()
    classes = []
    for i in range(G.n):
        if i in seen:
            continue
        cls = {conj_idx(G, g, i) for g in range(G.n)}
        classes.append(sorted(cls))
        seen |= cls
    classes.sort(key=lambda c: (element_order(G, c[0]), len(c), c[0]))
    return classes

def successor_classes(G, w, class_of):
    """Classes of [w, h w h^-1] over all h (e excluded)."""
    out = set()
    for h in range(G.n):
        w2 = comm_idx(G, w, conj_idx(G, h, w))
        if w2 != G.e:
            out.add(class_of[w2])
    return out

def first_loop_pair(G, w):
    """Lexicographically-first (g,h) with [g w g^-1, h w h^-1] == w, or None."""
    for g in range(G.n):
        wg = conj_idx(G, g, w)
        for h in range(G.n):
            if comm_idx(G, wg, conj_idx(G, h, w)) == w:
                return (g, h)
    return None

def build_graph(G, label_of):
    classes = conjugacy_classes(G)
    class_of = {}
    for ci, cls in enumerate(classes):
        for i in cls:
            class_of[i] = ci
    labels = [label_of(G, cls) for cls in classes]
    # disambiguate duplicate labels (e.g. two order-7 classes in GL3(F2)):
    # suffix a/b/... in order of appearance, so dict keys never collide
    from collections import Counter as _C
    dup = {k for k, v in _C(labels).items() if v > 1}
    seen = _C()
    for i, l in enumerate(labels):
        if l in dup:
            labels[i] = l + "abcdefgh"[seen[l]]
            seen[l] += 1
    adj = {}
    loops = {}
    for ci, cls in enumerate(classes):
        if cls[0] == G.e:
            continue
        succ = successor_classes(G, cls[0], class_of)
        adj[labels[ci]] = sorted(labels[cj] for cj in succ)
        if ci in succ:
            loops[labels[ci]] = first_loop_pair(G, cls[0])
    return classes, labels, adj, loops

# --------------------------------------------------------------------------
# Instances
# --------------------------------------------------------------------------

def s5_label(G, cls):
    p = G.elems[cls[0]]
    seen, out = set(), []
    for i in range(5):
        if i in seen:
            continue
        c, j = 0, i
        while j not in seen:
            seen.add(j)
            j = p[j]
            c += 1
        if c > 1:
            out.append(c)
    return str(tuple(sorted(out, reverse=True))) if out else "e"

def order_size_label(G, cls):
    return "o%d_s%d" % (element_order(G, cls[0]), len(cls))

def gl3f2_group():
    """GL3(F2), order 168, elements as 9-bit row-major tuples."""
    def mmul(a, b):
        return tuple((a[3 * r] & b[3 * 0 + c]) ^ (a[3 * r + 1] & b[3 * 1 + c])
                     ^ (a[3 * r + 2] & b[3 * 2 + c])
                     for r in range(3) for c in range(3))
    ident = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    gens = [(0, 1, 0, 0, 0, 1, 1, 0, 0),      # cyclic rotation
            (1, 1, 0, 0, 1, 0, 0, 0, 1)]      # transvection
    elems = {ident}
    frontier = [ident]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = mmul(x, g)
            if y not in elems:
                elems.add(y)
                frontier.append(y)
    return core.Group("GL3_2", sorted(elems), mmul, ident)

# Pinned expectations (derived by this script in derivation mode, then
# pinned; None => derivation mode prints them instead of checking).
PIN = {
    "S5_adj": {'(2,)': ['(3,)'], '(2, 2)': ['(3,)', '(5,)'],
               '(3,)': ['(2, 2)', '(3,)'],
               '(4,)': ['(2, 2)', '(3,)', '(5,)'],
               '(5,)': ['(3,)', '(5,)'], '(3, 2)': ['(3,)', '(5,)']},
    "S5_loops": {'(3,)': (6, 43), '(5,)': (7, 4)},
    "GL3_2_adj": {'o2_s21': ['o2_s21', 'o3_s56'],
                  'o3_s56': ['o2_s21', 'o3_s56', 'o7_s24a', 'o7_s24b'],
                  'o4_s42': ['o3_s56', 'o4_s42', 'o7_s24a', 'o7_s24b'],
                  'o7_s24a': ['o3_s56', 'o4_s42'],
                  'o7_s24b': ['o3_s56', 'o4_s42']},
    "GL3_2_loops": {'o2_s21': (2, 4), 'o3_s56': (11, 50), 'o4_s42': (11, 23)},
    "PGL2_7_adj": {'o2_s21': ['o2_s21', 'o3_s56'],
                   'o2_s28': ['o2_s21', 'o3_s56', 'o7_s48'],
                   'o3_s56': ['o2_s21', 'o3_s56', 'o7_s48'],
                   'o4_s42': ['o3_s56', 'o4_s42', 'o7_s48'],
                   'o6_s56': ['o2_s21', 'o3_s56', 'o4_s42', 'o7_s48'],
                   'o7_s48': ['o3_s56', 'o4_s42'],
                   'o8_s42a': ['o3_s56', 'o4_s42', 'o7_s48'],
                   'o8_s42b': ['o2_s21', 'o3_s56', 'o4_s42']},
    "PGL2_7_loops": {'o2_s21': (1, 6), 'o3_s56': (1, 5), 'o4_s42': (0, 26)},
}

def run_instance(name, G, label_of, extra=None):
    print("PART - %s (|G| = %d)" % (name, G.n), flush=True)
    classes, labels, adj, loops = build_graph(G, label_of)
    print("  classes: %s" % ", ".join("%s(x%d)" % (l, len(c))
                                      for l, c in zip(labels, classes)), flush=True)
    for k in sorted(adj):
        print("    %s -> %s%s" % (k, adj[k],
              "   LOOP pair %s" % (loops[k],) if k in loops else ""), flush=True)
    if PIN[name + "_adj"] is None:
        print("  [derivation] %s_adj = %r" % (name, adj), flush=True)
        print("  [derivation] %s_loops = %r" % (name, {k: v for k, v in loops.items()}),
              flush=True)
    else:
        check("%s adjacency matches pins" % name, adj == PIN[name + "_adj"])
        check("%s loop set + lexicographically-first pairs match pins" % name,
              loops == PIN[name + "_loops"])
    if extra:
        extra(G, classes, labels, adj, loops)
    return adj, loops

def s5_extra(G, classes, labels, adj, loops):
    # Phase-4A re-derivations through the graph
    check("gate lemma: successors of (2) = {(3,)}", adj["(2,)"] == ["(3,)"])
    check("successors of (3) = {(2,2),(3)}: NO direct (3)->(5) edge",
          adj["(3,)"] == ["(2, 2)", "(3,)"])
    check("(3) has a self-loop (threshold-2 direct family engine)",
          "(3,)" in loops)
    check("(5) has a self-loop (exact-3 family engine)", "(5,)" in loops)
    check("(5) reachable from (3) only through (2,2)",
          "(5,)" not in adj["(3,)"] and "(5,)" in adj["(2, 2)"])
    # every target class lies in A5 (successors always land in G')
    a5types = {"(3,)", "(2, 2)", "(5,)"}
    check("all successor classes lie in A5-classes",
          all(set(v) <= a5types for v in adj.values()))

# GL3(F2) transvection seed (literal 3x3 matrices over F2, row-major bits).
# Found independently by the 4A.2 gate reviewer and re-verified here.
GL_SEED = {
    "t1": (0, 0, 1, 0, 1, 0, 1, 0, 0),   # [[0,0,1],[0,1,0],[1,0,0]]
    "t2": (0, 1, 0, 1, 0, 0, 1, 1, 1),   # [[0,1,0],[1,0,0],[1,1,1]]
    "t3": (0, 1, 1, 1, 0, 1, 0, 0, 1),   # [[0,1,1],[1,0,1],[0,0,1]]
    "g":  (0, 0, 1, 0, 1, 0, 1, 1, 0),   # loop conjugator (order 4)
    "h":  (0, 0, 1, 0, 1, 1, 1, 0, 0),   # loop conjugator (order 4)
}

def gl3f2_extra(G, classes, labels, adj, loops):
    # exact transvection census: 21 involutions, every one has rank(g-I) = 1
    invs = G.involutions
    check("GL3(F2): exactly 21 involutions, all with rank(g-I) = 1",
          len(invs) == 21 and
          all(rank_g_minus_I_f2(G.elems[i]) == 1 for i in invs))
    # level-1 mirror of the S5 gate fact: exact census of valid ordered pairs
    census = {}
    for a in invs:
        for b in invs:
            if a != b and comm_idx(G, a, b) != G.e:
                o = len(G.closure([a, b]))
                census[o] = census.get(o, 0) + 1
    check("GL3(F2): level-1 census = 336 valid ordered pairs, "
          "168 close to order 6 and 168 to order 8 (max |H| == 8, all PROPER)",
          census == {6: 168, 8: 168} and max(census) == 8, census)
    # generating transvection seed at level 2, with literal matrices
    t1, t2, t3 = (G.idx[GL_SEED[k]] for k in ("t1", "t2", "t3"))
    g, h = G.idx[GL_SEED["g"]], G.idx[GL_SEED["h"]]
    for name in ("t1", "t2", "t3"):
        e = GL_SEED[name]
        check("seed leaf %s is a transvection (order 2, rank(g-I)=1)" % name,
              element_order(G, G.idx[e]) == 2 and rank_g_minus_I_f2(e) == 1)
    c1, c2 = comm_idx(G, t1, t2), comm_idx(G, t1, t3)
    root = comm_idx(G, c1, c2)
    check("seed gates valid; root target is again the transvection t1 (literal)",
          c1 != G.e and c2 != G.e and root == t1)
    check("seed leaves generate the FULL GL3(F2) (order 168)",
          len(G.closure([t1, t2, t3])) == G.n)
    check("pinned loop pair (g,h) fixes the literal root: [g r g^-1, h r h^-1] = r",
          comm_idx(G, conj_idx(G, g, root), conj_idx(G, h, root)) == root)
    check("(g,h) coincides with the pinned lexicographically-first loop pair "
          "of the transvection class", (g, h) == loops["o2_s21"], (g, h))
    check("both loop conjugators have order 4",
          element_order(G, g) == 4 and element_order(G, h) == 4)
    # self-similar family from level 2: leaves (t1,t2,t1,t3), doubling by (g,h)
    lv = [t1, t2, t1, t3]
    k = 2
    while True:
        al = lv[0]
        cur = list(lv)
        while len(cur) > 1:
            cur = [comm_idx(G, cur[i], cur[i + 1]) for i in range(0, len(cur), 2)]
        al = cur[0]
        check("GL3(F2) family: level %d target = t1, leaves transvections, "
              "count 2^%d, H full" % (k, k),
              al == t1 and len(lv) == 2 ** k and
              all(element_order(G, x) == 2 and
                  rank_g_minus_I_f2(G.elems[x]) == 1 for x in lv) and
              len(G.closure(sorted(set(lv)))) == G.n)
        if k == 5:
            break
        lv = ([conj_idx(G, g, x) for x in lv] + [conj_idx(G, h, x) for x in lv])
        k += 1

def pgl2_7_extra(G, classes, labels, adj, loops):
    # interior/exterior certification: det mod squares is PGL-invariant
    squares = {1, 2, 4}
    def detmod(i):
        a, b, c, d = G.elems[i]
        return (a * d - b * c) % 7
    for cls, lab in zip(classes, labels):
        if element_order(G, cls[0]) != 2:
            continue
        dets = {detmod(i) for i in cls}
        if len(cls) == 21:
            check("PGL2(7) class o2_s21 is INTERIOR (all dets are squares)",
                  dets <= squares, sorted(dets))
        else:
            check("PGL2(7) class o2_s28 is EXTERIOR (all dets are non-squares)",
                  dets.isdisjoint(squares), sorted(dets))

def main():
    S5 = core.sym_group(5, "S5")
    run_instance("S5", S5, s5_label, extra=s5_extra)
    GL = gl3f2_group()
    check("GL3(F2) has order 168", GL.n == 168)
    run_instance("GL3_2", GL, order_size_label, extra=gl3f2_extra)
    P7 = core.pgl2(7, "PGL2_7")
    run_instance("PGL2_7", P7, order_size_label, extra=pgl2_7_extra)
    print()
    if FAILURES:
        print("TARGET TRANSITION GRAPH CERTIFICATE: FAIL (%d)" % len(FAILURES),
              flush=True)
        return 1
    print("TARGET TRANSITION GRAPH CERTIFICATE: PASS", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
