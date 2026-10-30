#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-class bridge sweep certificate (Phase 1).

Sweeps every conjugacy class of involutions SEPARATELY (leaves restricted to
one ambient class) and pins the observed counts. These finite zeros are data,
NOT a proof of any general conjecture — witnesses outside the swept range
exist (see bridge_verify_counterexample.py: a5-order3, a7-singleclass).

Suites:
  depth2-exact    : exhaustive depth-2 sweeps for S4 (both classes), A5,
                    S5 (both), S6 (transpositions and (2,2,2)), and BOTH
                    involution classes of PGL2(7) (614,656 and 194,481
                    templates, exhaustive). The S6 (2,2) class (45^4 = 4.1M)
                    remains SAMPLED at 60,000 and is labelled as such.
  depth3-evidence : reproduces exactly the original depth-3 samples
                    (4,000 templates per class, seed 20260711).

Exit codes: 0 = all pinned invariants hold; 1 = mathematical failure or
count regression; 2 = usage error / existing outputs without --overwrite.

Outputs are canonical and byte-identical across repeat runs (see the
mixed-class script for the format).
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

# plan rows: (group, class_index, depth, n_samples, pinned_stats or None)
# Pins for the previously-exhaustive sweeps were derived 2026-07-11 by the
# frozen scripts (singleclass_results.txt); the exhaustive PGL2(7) rows and
# the A5 row were derived by this script and reproduced by the Phase-1 gate
# before pinning. T_and_not_S is pinned to 0 in EVERY row.
EXPECTED = {
    "depth2-exact": {
        "sweeps": [
            ("S4", 0, 2, None, _st(1296, 1296, 792, 0, 504, 0, 360)),
            ("S4", 1, 2, None, _st(81, 81, 81, 0, 0, 0, 81)),
            ("A5", 0, 2, None, _st(50625, 50625, 32385, 0, 18240, 2160, 4905)),
            ("S5", 0, 2, None, _st(10000, 10000, 5560, 0, 4440, 0, 2320)),
            ("S5", 1, 2, None, _st(50625, 50625, 32385, 0, 18240, 2160, 4905)),
            ("S6", 0, 2, None, _st(50625, 50625, 25785, 0, 24840, 0, 12825)),
            ("S6", 1, 2, 60000, _st(60000, 60000, 26789, 0, 33211, 21039, 1594)),
            ("S6", 2, 2, None, _st(50625, 50625, 25785, 0, 24840, 0, 12825)),
            # exhaustive PGL2(7): full stats derived 2026-07-11 by this suite
            # (first exhaustive run) and pinned after derivation.
            ("PGL2_7", 0, 2, None, _st(614656, 614656, 351904, 0, 262752, 108192, 30352)),
            ("PGL2_7", 1, 2, None, _st(194481, 194481, 60417, 0, 134064, 57792, 13041)),
        ],
        "records_total": 0,
        "requires_seed": LEGACY_SEED,   # for the single sampled row
    },
    "depth3-evidence": {
        "sweeps": [
            ("S4", 0, 3, 4000, _st(4000, 12000, 8925, 0, 3075, 0, 6221)),
            ("S4", 1, 3, 4000, _st(4000, 12000, 12000, 0, 0, 0, 12000)),
            ("S5", 0, 3, 4000, _st(4000, 12000, 7026, 0, 4974, 0, 3928)),
            ("S5", 1, 3, 4000, _st(4000, 12000, 8427, 0, 3573, 311, 2741)),
            ("S6", 0, 3, 4000, _st(4000, 12000, 5851, 0, 6149, 21, 3580)),
            ("S6", 1, 3, 4000, _st(4000, 12000, 7321, 0, 4679, 2803, 1693)),
            ("S6", 2, 3, 4000, _st(4000, 12000, 5882, 0, 6118, 24, 3590)),
            ("PGL2_7", 0, 3, 4000, _st(4000, 12000, 8201, 0, 3799, 1414, 2093)),
            ("PGL2_7", 1, 3, 4000, _st(4000, 12000, 5721, 0, 6279, 2467, 2656)),
        ],
        "records_total": 0,
        "requires_seed": LEGACY_SEED,
    },
}

EXPECTED_CLASS_SIZES = {"S4": [6, 3], "A5": [15], "S5": [10, 15],
                        "S6": [15, 45, 15], "PGL2_7": [28, 21]}

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

def run_suite(suite, seed, outdir, overwrite, expected=None):
    if suite not in EXPECTED:
        print("unknown suite: %s" % suite, file=sys.stderr)
        return 2
    exp = EXPECTED[suite] if expected is None else expected
    req = exp.get("requires_seed")
    if req is not None and seed != req:
        print("suite %s pins fixed-seed samples; it requires --seed %d"
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
    for gn in groups:
        got = [len(c) for c in classes[gn]]
        assert got == EXPECTED_CLASS_SIZES[gn], (gn, got)
    core.core_selftest()

    records = []
    sweep_rows = []
    mismatches = []
    for (gname, ci, depth, n_samples, pin) in plan:
        G = groups[gname]
        pool = classes[gname][ci]
        stats = core.sweep(G, depth, pool, seed=seed, n_samples=n_samples,
                           tag="class%d" % ci, record_sink=records.append)
        row = {"group": gname, "class_index": ci, "class_size": len(pool),
               "depth": depth, "mode": stats["mode"], "n_samples": n_samples,
               "stats": {k: stats[k] for k in S}}
        sweep_rows.append(row)
        print("  == %s class%d(sz%d) d=%d %s: %s"
              % (gname, ci, len(pool), depth, stats["mode"], row["stats"]),
              flush=True)
        # pins may be partial: only the pinned fields are compared
        if pin is not None and any(row["stats"][k] != v for k, v in pin.items()):
            mismatches.append({"sweep": [gname, ci, depth],
                               "expected": pin, "got": row["stats"]})

    derived = {
        "records_total": len(records),
        "templates_total": sum(r["stats"]["templates"] for r in sweep_rows),
        "sampled_rows": [i for i, r in enumerate(sweep_rows)
                         if r["mode"] == "sampled"],
    }
    if derived["records_total"] != exp["records_total"]:
        mismatches.append({"census": "records_total",
                           "expected": exp["records_total"],
                           "got": derived["records_total"]})

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
            "bridge_singleclass_sweep.py": sha256_file(os.path.abspath(__file__)),
        },
        "sweeps": sweep_rows,
        "derived": derived,
        "mismatches": mismatches,
        "verdict": verdict,
        "note": ("finite zeros over the swept classes; NOT a proof of any "
                 "general single-class statement — see "
                 "bridge_verify_counterexample.py cases a5-order3 and "
                 "a7-singleclass"),
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
        description="Single-class bridge sweep certificate (Phase 1)")
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
