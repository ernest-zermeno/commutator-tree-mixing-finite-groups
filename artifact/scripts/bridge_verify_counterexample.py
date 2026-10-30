#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Deterministic counterexample verifier for the bridge hypotheses (Phase 1).

Each case ASSERTS (not merely prints) every invariant; any failure exits 1.

  --case s5-mixed        S5, depth 2, leaves (1 2)(3 4), (1 5) | (1 2)(4 5), (3 5):
                         tower (T) holds, saturation (S) fails (3 vs 60).
                         Kills "per-step + (T) => (S)" with mixed involution types.
  --case a5-order3       A5, depth 1, sigma=(1 2 3), tau=(3 4 5): SAME ambient
                         class (order 3), generate A5, (T) holds at level 1,
                         (S) fails (3 vs 60). Kills any single-class statement
                         that allows non-involution leaves.
  --case a7-singleclass  A7, depth 2, four leaves in ONE ambient class of A7
                         ((2,2)-involutions), |H|=120, |HL|=|HR|=12, (T) holds,
                         (S) fails (3 vs 60). Kills the ambient single-class
                         INVOLUTION conjecture as stated.
  --case all             run the three cases.

Uses only the pure permutation helpers of bridge_sweep_core (no tables).
Exit codes: 0 all requested cases pass; 1 any assertion fails; 2 usage error.
"""
import argparse
import sys

from bridge_sweep_core import (pmul, pinv, pcomm, pidentity, pparity,
                               pclosure, pcomm_subgroup, pderived,
                               cycles_to_perm, an_conjugator)

def V_gate(s, t, n):
    """Per-step value set of a depth-1 gate with leaves s, t:
    V = {e,s}.{e,t}.{e,s}^-1.{e,t}^-1 (element inverses)."""
    Vl = {pidentity(n), s}
    Vr = {pidentity(n), t}
    return {pmul(pmul(pmul(x, y), pinv(xp)), pinv(yp))
            for x in Vl for y in Vr for xp in Vl for yp in Vr}

def gate_TS(WL_V, WR_V, HL, HR, Hv, k, n):
    """(T) and (S) data at a gate of level k given children value sets and
    leaf subgroups. Returns (T_ok, S_lhs, S_rhs)."""
    HLd = pderived(HL, n, k - 1)
    HRd = pderived(HR, n, k - 1)
    Hvk = pderived(Hv, n, k)
    T_lhs = pcomm_subgroup(HLd, HRd, n)
    T_ok = (T_lhs == Hvk)
    WL = WL_V & HLd
    WR = WR_V & HRd
    S_lhs = pclosure([pcomm(a, b) for a in WL for b in WR], n)
    S_rhs = pcomm_subgroup(pclosure(WL, n), pclosure(WR, n), n)
    return T_ok, S_lhs, S_rhs

def case_s5_mixed():
    n = 5
    a1 = cycles_to_perm([(1, 2), (3, 4)], n)
    a2 = cycles_to_perm([(1, 5)], n)
    b1 = cycles_to_perm([(1, 2), (4, 5)], n)
    b2 = cycles_to_perm([(3, 5)], n)
    for s in (a1, a2, b1, b2):
        assert pmul(s, s) == pidentity(n), "leaf is not an involution"
    HL = pclosure([a1, a2], n)
    HR = pclosure([b1, b2], n)
    Hv = pclosure([a1, a2, b1, b2], n)
    assert len(HL) == 12 and len(HR) == 12, (len(HL), len(HR))
    assert len(Hv) == 120, len(Hv)
    assert len(pderived(HL, n, 1)) == 3 and len(pderived(HR, n, 1)) == 3
    assert len(pderived(Hv, n, 2)) == 60
    T_ok, S_lhs, S_rhs = gate_TS(V_gate(a1, a2, n), V_gate(b1, b2, n),
                                 HL, HR, Hv, 2, n)
    assert T_ok, "(T) must hold"
    assert len(S_lhs) == 3 and len(S_rhs) == 60, (len(S_lhs), len(S_rhs))
    assert S_lhs != S_rhs, "(S) must fail"
    root = pcomm(pcomm(a1, a2), pcomm(b1, b2))
    assert root != pidentity(n), "root target must be nontrivial"
    print("  s5-mixed: PASS  (|HL|=|HR|=12, |H|=120, |HL'|=|HR'|=3, |H''|=60; "
          "(T) holds; (S) fails 3 vs 60; root target nontrivial)")

def case_a5_order3():
    n = 5
    sg = cycles_to_perm([(1, 2, 3)], n)
    ta = cycles_to_perm([(3, 4, 5)], n)
    # order 3, even, same A5 conjugacy class (explicit conjugator)
    assert pmul(pmul(sg, sg), sg) == pidentity(n)
    assert pmul(pmul(ta, ta), ta) == pidentity(n)
    assert pparity(sg) == 0 and pparity(ta) == 0
    assert an_conjugator(sg, ta, n) is not None, \
        "sigma and tau must be A5-conjugate"
    Hv = pclosure([sg, ta], n)
    assert len(Hv) == 60, "leaves must generate A5 (got %d)" % len(Hv)
    HL = pclosure([sg], n)
    HR = pclosure([ta], n)
    # level-1 gate: children are the leaves themselves; V(leaf) = {e, sigma}
    T_ok, S_lhs, S_rhs = gate_TS({pidentity(n), sg}, {pidentity(n), ta},
                                 HL, HR, Hv, 1, n)
    assert T_ok, "(T) must hold at level 1"
    assert len(S_lhs) == 3 and len(S_rhs) == 60, (len(S_lhs), len(S_rhs))
    assert S_lhs != S_rhs, "(S) must fail"
    assert pcomm(sg, ta) != pidentity(n), "root target must be nontrivial"
    print("  a5-order3: PASS  (single A5 class of 3-cycles; <sigma,tau>=A5; "
          "(T) holds at level 1; (S) fails 3 vs 60)")

def case_a7_singleclass():
    n = 7
    a1 = cycles_to_perm([(1, 2), (3, 4)], n)
    a2 = cycles_to_perm([(1, 5), (6, 7)], n)
    b1 = cycles_to_perm([(1, 2), (4, 5)], n)
    b2 = cycles_to_perm([(3, 5), (6, 7)], n)
    for s in (a1, a2, b1, b2):
        assert pmul(s, s) == pidentity(n) and pparity(s) == 0
    # one single ambient class of A7: explicit conjugators from a1
    for other in (a2, b1, b2):
        assert an_conjugator(a1, other, n) is not None, \
            "leaves must lie in one A7 class"
    HL = pclosure([a1, a2], n)
    HR = pclosure([b1, b2], n)
    Hv = pclosure([a1, a2, b1, b2], n)
    assert len(HL) == 12 and len(HR) == 12, (len(HL), len(HR))
    assert len(Hv) == 120, len(Hv)
    assert len(pderived(Hv, n, 2)) == 60
    T_ok, S_lhs, S_rhs = gate_TS(V_gate(a1, a2, n), V_gate(b1, b2, n),
                                 HL, HR, Hv, 2, n)
    assert T_ok, "(T) must hold"
    assert len(S_lhs) == 3 and len(S_rhs) == 60, (len(S_lhs), len(S_rhs))
    assert S_lhs != S_rhs, "(S) must fail"
    root = pcomm(pcomm(a1, a2), pcomm(b1, b2))
    assert root != pidentity(n), "root target must be nontrivial"
    print("  a7-singleclass: PASS  (four leaves in ONE A7 class; |H|=120, "
          "|HL|=|HR|=12; (T) holds; (S) fails 3 vs 60; root nontrivial)")

CASES = {"s5-mixed": case_s5_mixed,
         "a5-order3": case_a5_order3,
         "a7-singleclass": case_a7_singleclass}

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Deterministic bridge counterexample verifier (Phase 1)")
    ap.add_argument("--case", required=True, choices=sorted(CASES) + ["all"])
    try:
        a = ap.parse_args(argv)
    except SystemExit as exc:
        # preserve argparse's own contract: --help -> 0, usage errors -> 2
        return int(exc.code)
    todo = sorted(CASES) if a.case == "all" else [a.case]
    try:
        for c in todo:
            CASES[c]()
    except AssertionError as e:
        print("WITNESS FAILURE (%s): %s" % (c, e), file=sys.stderr)
        return 1
    print("BRIDGE WITNESSES: PASS (%s)" % ", ".join(todo))
    return 0

if __name__ == "__main__":
    sys.exit(main())
