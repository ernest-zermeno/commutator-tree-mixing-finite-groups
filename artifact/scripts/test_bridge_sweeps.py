#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the Phase-1 bridge sweep infrastructure (standard unittest).

Heavy enumerations are NOT run here; suites used by these tests are the tiny
smoke-s4 (exhaustive, ~2 s) and the witness verifier (~3 s).

Run:  python3 -m unittest test_bridge_sweeps -v      (from scripts/)
"""
import hashlib
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

sys.path.insert(0, HERE)
import bridge_sweep_core as core            # noqa: E402
import bridge_mixedclass_sweep as mcs       # noqa: E402

def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()

def run(args, cwd=None):
    return subprocess.run([PY] + args, cwd=cwd or HERE,
                          capture_output=True, text=True)

class TestCoreImportSideEffects(unittest.TestCase):
    def test_import_core_writes_no_files(self):
        with tempfile.TemporaryDirectory() as td:
            r = subprocess.run(
                [PY, "-c", "import sys; sys.path.insert(0, %r); "
                           "import bridge_sweep_core" % HERE],
                cwd=td, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            leftover = [f for f in os.listdir(td) if f != "__pycache__"]
            self.assertEqual(leftover, [], "import wrote files: %s" % leftover)

class TestCoreSanity(unittest.TestCase):
    def test_selftest(self):
        self.assertTrue(core.core_selftest())

    def test_s5_p0_l4_sanity(self):
        # <V(root) & A5> for the P0 l=4 template has order 60
        self.assertTrue(core.s5_p0_l4_sanity())

    def test_involution_class_sizes(self):
        self.assertEqual([len(c) for c in
                          core.involution_classes(core.sym_group(4, "S4"))], [6, 3])
        self.assertEqual([len(c) for c in
                          core.involution_classes(core.alt_group(5, "A5"))], [15])
        self.assertEqual([len(c) for c in
                          core.involution_classes(core.sym_group(5, "S5"))], [10, 15])
        self.assertEqual([len(c) for c in
                          core.involution_classes(core.pgl2(7, "PGL2_7"))], [28, 21])

class TestCLI(unittest.TestCase):
    def test_help_exits_zero(self):
        for script in ("bridge_mixedclass_sweep.py", "bridge_singleclass_sweep.py",
                       "bridge_verify_counterexample.py"):
            r = run([os.path.join(HERE, script), "--help"])
            self.assertEqual(r.returncode, 0, (script, r.stderr))
            self.assertIn("usage:", r.stdout)

    def test_invalid_cli_returns_2(self):
        for script in ("bridge_mixedclass_sweep.py", "bridge_singleclass_sweep.py"):
            r = run([os.path.join(HERE, script), "--suite", "no-such-suite",
                     "--output-dir", "/tmp/x"])
            self.assertEqual(r.returncode, 2, (script, r.stderr))
            r = run([os.path.join(HERE, script)])
            self.assertEqual(r.returncode, 2, (script, r.stderr))
        r = run([os.path.join(HERE, "bridge_verify_counterexample.py"),
                 "--case", "no-such-case"])
        self.assertEqual(r.returncode, 2, r.stderr)

    def test_wrong_seed_for_pinned_suite_returns_2(self):
        with tempfile.TemporaryDirectory() as td:
            r = run([os.path.join(HERE, "bridge_mixedclass_sweep.py"),
                     "--suite", "legacy-regression", "--seed", "1",
                     "--output-dir", td])
            self.assertEqual(r.returncode, 2, r.stderr)

    def test_overwrite_protection_and_byte_identity(self):
        with tempfile.TemporaryDirectory() as td:
            a = run([os.path.join(HERE, "bridge_mixedclass_sweep.py"),
                     "--suite", "smoke-s4", "--output-dir", td])
            self.assertEqual(a.returncode, 0, a.stderr)
            p_sum = os.path.join(td, "smoke-s4_summary.json")
            p_cex = os.path.join(td, "smoke-s4_counterexamples.jsonl")
            h1 = (sha256(p_sum), sha256(p_cex))
            # second run without --overwrite must abort with 2, files untouched
            b = run([os.path.join(HERE, "bridge_mixedclass_sweep.py"),
                     "--suite", "smoke-s4", "--output-dir", td])
            self.assertEqual(b.returncode, 2, b.stderr)
            self.assertEqual((sha256(p_sum), sha256(p_cex)), h1)
            # third run WITH --overwrite must be byte-identical
            c = run([os.path.join(HERE, "bridge_mixedclass_sweep.py"),
                     "--suite", "smoke-s4", "--output-dir", td, "--overwrite"])
            self.assertEqual(c.returncode, 0, c.stderr)
            self.assertEqual((sha256(p_sum), sha256(p_cex)), h1)

class TestWitnesses(unittest.TestCase):
    def test_all_witnesses_pass(self):
        r = run([os.path.join(HERE, "bridge_verify_counterexample.py"),
                 "--case", "all"])
        self.assertEqual(r.returncode, 0, r.stderr)
        for tag in ("s5-mixed", "a5-order3", "a7-singleclass"):
            self.assertIn(tag + ": PASS", r.stdout)

class TestTamperedExpectationFails(unittest.TestCase):
    def test_tampered_expectation_exits_1(self):
        tampered = {
            "sweeps": [("S4", "all-inv", 2, None,
                        dict(mcs._st(6561, 6561, 3609, 0, 2952, 0, 2025),
                             T_and_not_S=7))],   # deliberate lie
            "records_total": 7,
            "all_k2": True,
            "all_mixed_ambient": True,
            "mono_class_failures": 0,
        }
        with tempfile.TemporaryDirectory() as td:
            rc = mcs.run_suite("smoke-s4", mcs.LEGACY_SEED, td, True,
                               expected=tampered)
            self.assertEqual(rc, 1)

if __name__ == "__main__":
    unittest.main()
