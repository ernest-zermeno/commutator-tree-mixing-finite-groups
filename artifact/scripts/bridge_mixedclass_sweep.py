#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mixed-class bridge sweep certificate (Phase 1).

Suites:
  legacy-regression : reproduces the original 12-sweep run (seed 20260711)
                      and pins its counterexample census.
  s5-depth2-exact   : exhaustive 25^4 = 390,625 quadruples over all S5
                      involutions; pins 24,480 (T)&~(S) failures, all of
                      them mixing ambient involution classes.
  smoke-s4          : tiny exhaustive suite (S4, depth 2) used by the unit
                      tests for byte-identity / overwrite / tamper checks.

Exit codes: 0 = all pinned invariants hold; 1 = mathematical failure or
count regression; 2 = usage error (bad arguments, wrong seed for a pinned
suite, or existing outputs without --overwrite).

Outputs (per suite, inside --output-dir): <suite>_summary.json (canonical,
sorted keys, no timestamps, no absolute paths) and
<suite>_counterexamples.jsonl. Repeat runs are byte-identical.
"""
import argparse
import hashlib
import json
import os
import platform
import sys

import numpy as np

import bridge_sweep_core as core

SCHEMA_VERSION = "1.0"
LEGACY_SEED = 20260711

S = ("templates", "gates_k2plus", "T_and_S", "T_and_not_S", "not_T",
     "a_and_not_T", "T_and_not_a")

def _st(*vals):
    return dict(zip(S, vals))

# Pinned expectations. Derivation history: the 12 legacy sweep rows were
# derived on 2026-07-11 by the frozen scripts (bridge_mixedclass_results.txt);
# the census metrics (records / unique leafsets / unique legacy keys) and the
# s5-depth2-exact numbers were reproduced independently by the Phase-1 gate
# before being pinned here. A mismatch aborts with exit 1.
EXPECTED = {
    "legacy-regression": {
        "sweeps": [
            ("S4", "all-inv", 2, None, _st(6561, 6561, 3609, 0, 2952, 0, 2025)),
            ("A5", "all-inv", 2, None, _st(50625, 50625, 32385, 0, 18240, 2160, 4905)),
            ("S5", "all-inv", 2, 30000, _st(30000, 30000, 11302, 1878, 16820, 6262, 2680)),
            ("S6", "transpositions", 2, None, _st(50625, 50625, 25785, 0, 24840, 0, 12825)),
            ("S6", "all-inv", 2, 30000, _st(30000, 30000, 6231, 2592, 21177, 11398, 1381)),
            ("PGL2_7", "all-inv", 2, 30000, _st(30000, 30000, 8750, 5771, 15479, 7175, 1160)),
            ("S4", "all-inv", 3, 8000, _st(8000, 24000, 16666, 0, 7334, 0, 12899)),
            ("A5", "all-inv", 3, 8000, _st(8000, 24000, 16887, 0, 7113, 613, 5454)),
            ("S5", "all-inv", 3, 6000, _st(6000, 18000, 9461, 785, 7754, 2456, 4292)),
            ("S6", "transpositions", 3, 3000, _st(3000, 9000, 4390, 0, 4610, 17, 2648)),
            ("S6", "all-inv", 3, 3000, _st(3000, 9000, 4012, 497, 4491, 2235, 1849)),
            ("PGL2_7", "all-inv", 3, 3000, _st(3000, 9000, 4522, 1142, 3336, 1425, 1638)),
        ],
        "records_total": 12665,
        "unique_gate_leafsets": 10151,
        "unique_legacy_key": 12575,
        # new metric (reported alongside, does not replace the legacy ones):
        # distinct under (group, depth, k, full_template_leaves, gate_path)
        "unique_with_gate_path": 12575,
        "all_k2": True,
        "all_mixed_ambient": True,
        "mono_class_failures": 0,
        "requires_seed": LEGACY_SEED,
    },
    "s5-depth2-exact": {
        "sweeps": [
            # full stats derived 2026-07-11 by this suite and pinned after derivation
            ("S5", "all-inv", 2, None, _st(390625, 390625, 146425, 24480,
                                           219720, 82320, 33865)),
        ],
        "templates_total": 390625,
        "records_total": 24480,
        "all_k2": True,
        "all_mixed_ambient": True,
        "mono_class_failures": 0,
    },
    "smoke-s4": {
        "sweeps": [
            ("S4", "all-inv", 2, None, _st(6561, 6561, 3609, 0, 2952, 0, 2025)),
        ],
        "records_total": 0,
        "all_k2": True,
        "all_mixed_ambient": True,
        "mono_class_failures": 0,
    },
}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()

def canonical_dumps(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)

def build_groups(names):
    out = {}
    for n in names:
        if n == "S4":
            out[n] = core.sym_group(4, "S4")
        elif n == "S5":
            out[n] = core.sym_group(5, "S5")
        elif n == "S6":
            out[n] = core.sym_group(6, "S6")
        elif n == "A5":
            out[n] = core.alt_group(5, "A5")
        elif n == "PGL2_7":
            out[n] = core.pgl2(7, "PGL2_7")
        else:
            raise ValueError(n)
    return out

def pool_for(G, tag):
    if tag == "all-inv":
        return G.involutions
    if tag == "transpositions":
        assert G.name in ("S5", "S6")
        npts = 5 if G.name == "S5" else 6
        return [G.idx[p] for p in G.elems
                if sum(p[i] != i for i in range(npts)) == 2]
    raise ValueError(tag)

def run_suite(suite, seed, outdir, overwrite, expected=None):
    if suite not in EXPECTED:
        print("unknown suite: %s" % suite, file=sys.stderr)
        return 2
    exp = EXPECTED[suite] if expected is None else expected
    req = exp.get("requires_seed")
    if req is not None and seed != req:
        print("suite %s pins a fixed-seed regression; it requires --seed %d"
              % (suite, req), file=sys.stderr)
        return 2
    os.makedirs(outdir, exist_ok=True)
    p_sum = os.path.join(outdir, "%s_summary.json" % suite)
    p_cex = os.path.join(outdir, "%s_counterexamples.jsonl" % suite)
    if (os.path.exists(p_sum) or os.path.exists(p_cex)) and not overwrite:
        print("outputs already exist in %s; pass --overwrite" % outdir,
              file=sys.stderr)
        return 2

    plan = exp["sweeps"]
    groups = build_groups({g for (g, _, _, _, _) in plan})
    classes = {gn: core.involution_classes(G) for gn, G in groups.items()}
    cls_of = {gn: core.ambient_class_index(G, classes[gn])
              for gn, G in groups.items()}

    # sanity gate before any sweep (same spirit as the frozen scripts)
    if "S5" in groups:
        core.s5_p0_l4_sanity()
    core.core_selftest()

    records = []
    sweep_rows = []
    mismatches = []
    for (gname, ptag, depth, n_samples, pin) in plan:
        G = groups[gname]
        pool = pool_for(G, ptag)
        stats = core.sweep(G, depth, pool, seed=seed, n_samples=n_samples,
                           tag=ptag, record_sink=records.append)
        row = {"group": gname, "pool": ptag, "pool_size": len(pool),
               "depth": depth, "mode": stats["mode"],
               "n_samples": n_samples,
               "stats": {k: stats[k] for k in S}}
        sweep_rows.append(row)
        print("  == %s %s d=%d %s: %s" % (gname, ptag, depth, stats["mode"],
                                          row["stats"]), flush=True)
        # pins may be partial: only the pinned fields are compared
        if pin is not None and any(row["stats"][k] != v for k, v in pin.items()):
            mismatches.append({"sweep": [gname, ptag, depth],
                               "expected": pin, "got": row["stats"]})

    # derived census over all (T)&~(S) records
    derived = {
        "records_total": len(records),
        "all_k2": all(r["k"] == 2 for r in records),
        "unique_gate_leafsets": len({tuple(r["gate_leafset"]) for r in records}),
        "unique_legacy_key": len({(r["group"], r["depth"], r["k"],
                                   tuple(r["leaves"]), tuple(r["gate_leafset"]))
                                  for r in records}),
        "unique_with_gate_path": len({(r["group"], r["depth"], r["k"],
                                       tuple(r["leaves"]), r["gate_path"])
                                      for r in records}),
    }
    mono_fail = 0
    for r in records:
        cs = {cls_of[r["group"]][i] for i in r["gate_leaf_idx"]}
        if len(cs) < 2:
            mono_fail += 1
    derived["mono_class_failures"] = mono_fail
    derived["all_mixed_ambient"] = (mono_fail == 0)
    derived["templates_total"] = sum(row["stats"]["templates"]
                                     for row in sweep_rows)

    for key in ("records_total", "unique_gate_leafsets", "unique_legacy_key",
                "unique_with_gate_path", "all_k2", "all_mixed_ambient",
                "mono_class_failures", "templates_total"):
        if key in exp and derived[key] != exp[key]:
            mismatches.append({"census": key, "expected": exp[key],
                               "got": derived[key]})

    verdict = "PASS" if not mismatches else "FAIL"
    here = os.path.dirname(os.path.abspath(__file__))
    summary = {
        "schema_version": SCHEMA_VERSION,
        "suite": suite,
        "seed": seed,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "script_hashes": {
            "bridge_sweep_core.py":
                sha256_file(os.path.join(here, "bridge_sweep_core.py")),
            "bridge_mixedclass_sweep.py": sha256_file(os.path.abspath(__file__)),
        },
        "sweeps": sweep_rows,
        "derived": derived,
        "mismatches": mismatches,
        "verdict": verdict,
    }
    with open(p_sum, "w") as f:            # truncate, never append
        f.write(canonical_dumps(summary) + "\n")
    with open(p_cex, "w") as f:            # truncate, never append
        for r in records:
            f.write(canonical_dumps(r) + "\n")
    print("suite %s: %s (%d records) -> %s" % (suite, verdict, len(records),
                                               os.path.basename(p_sum)),
          flush=True)
    if mismatches:
        print("MISMATCHES: %s" % json.dumps(mismatches, default=str),
              file=sys.stderr)
        return 1
    return 0

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Mixed-class bridge sweep certificate (Phase 1)")
    ap.add_argument("--suite", required=True, choices=sorted(EXPECTED))
    ap.add_argument("--seed", type=int, default=LEGACY_SEED)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--overwrite", action="store_true")
    try:
        a = ap.parse_args(argv)
    except SystemExit as exc:
        # preserve argparse's own contract: --help -> 0, usage errors -> 2
        return int(exc.code)
    return run_suite(a.suite, a.seed, a.output_dir, a.overwrite)

if __name__ == "__main__":
    sys.exit(main())
