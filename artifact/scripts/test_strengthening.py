"""Regression tests for the new statements and previously silent failures.

Group computations are exact. Spectral tests are float64 cross-checks of the
analytic inequalities, not interval certificates or substitutes for their proofs.
"""
import ast
import contextlib
import io
import itertools
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

import numpy as np
from pgl2_classification import make_field, FIELDS, predicted_levels
from verification_groups import Group, mul, parity, pgl

ROOT = Path(__file__).resolve().parent


def definitions(name):
    """Load definitions without executing a script's expensive top-level sweep."""
    tree = ast.parse((ROOT / name).read_text(), filename=name)
    tree.body = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom,
                                                     ast.FunctionDef, ast.ClassDef))]
    ns = {'__name__': 'regression_loaded'}
    exec(compile(tree, name, 'exec'), ns)
    return ns


class StrengtheningTests(unittest.TestCase):
    def test_extension_fields_and_subfield_parity(self):
        for q in (9, 25, 27, 125):
            p, redux = FIELDS[q]
            F = make_field(p, redux)
            for a in range(1, q):
                self.assertEqual(F['mul'](a, F['inv'](a)), 1)
            # The embedded prime-field determinant 2 is nonsquare in F5.
            if p == 5:
                self.assertEqual(2 in F['squares'], q == 25)
        # A PGL2(5) subgroup inside PGL2(125) is nonsolvable and remains exterior-
        # generated; no enumeration of the ambient 1,953,000 elements is needed.
        G = Group(pgl(5))
        C = [i for i, g in enumerate(G.elements)
             if i != G.e and G.m[i][i] == G.e and parity(g)]
        H = G.close(C)
        self.assertEqual(len(H), 120)
        self.assertIsNone(G.derived_sizes(H)[1])
        self.assertLess(len(H), 125 * (125**2 - 1))

    def test_prediction_boundaries(self):
        # Tests formula encoding only; Dickson's classification supplies its proof.
        expected = {5:(3,2), 7:(2,3), 9:(2,None), 11:(3,None),
                    13:(3,2), 17:(2,3), 25:(2,None), 27:(3,2),
                    125:(None,None)}
        for q, (ext, inner) in expected.items():
            p, redux = FIELDS[q]
            row = predicted_levels(p, len(redux) if redux else 1)
            self.assertEqual(row, {'all_pgl':None, 'exterior_pgl':ext,
                                   'involutions_psl':inner})

    def test_affine_extension_group_not_dihedral(self):
        F = make_field(*FIELDS[27])
        add, neg = F['add'], F['neg']
        translations = [tuple(add(x,b) for x in range(27)) for b in range(27)]
        reflections = [tuple(add(neg(x),b) for x in range(27)) for b in range(27)]
        G = Group(translations + reflections)
        C = [G.index[g] for g in reflections]
        self.assertEqual(len(G.close(C)), 54)
        self.assertEqual(G.derived_sizes(frozenset(range(G.n))), ([54,27,1],2))
        e = tuple(range(27))
        self.assertTrue(all(mul(mul(t,t),t)==e for t in translations))

    def test_exterior_selector_and_composite_rejection(self):
        for q, count in ((5,10), (7,28)):
            r = subprocess.run([sys.executable, str(ROOT/'mixing_level_pgl.py'),
                                str(q), '--class-only'], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn('EXTERIOR (non-square determinant) leaves: '+str(count), r.stdout)
        r = subprocess.run([sys.executable, str(ROOT/'mixing_level_pgl.py'),
                            '9', '--class-only'], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('odd prime', r.stderr)

    def test_cayley_comparison(self):
        G = Group(itertools.permutations(range(3)))
        C = [i for i in range(G.n) if i!=G.e and G.m[i][i]==G.e]
        # Work on the regular representation minus constants.
        B = np.zeros((G.n, G.n-1))
        for k in range(1,G.n):
            B[:k,k-1]=1; B[k,k-1]=-k
            B[:,k-1] /= np.linalg.norm(B[:,k-1])
        Ps = {}
        for s in C:
            R = np.zeros((G.n,G.n))
            for j in range(G.n): R[G.m[s][j],j]=1
            Ps[s] = B.T @ ((np.eye(G.n)+R)/2) @ B
        for N in range(1,6):
            for seq in itertools.product(C,repeat=N):
                A = sum(Ps[s] for s in seq)/N
                delta = max(0.0,1-np.linalg.eigvalsh(A)[-1])
                T = np.eye(G.n-1)
                for s in seq: T=Ps[s]@T
                loss=1-np.linalg.norm(T,2)**2
                self.assertGreaterEqual(loss+2e-12,2*delta/(N+1))
                self.assertLessEqual(loss,min(1,4*N*delta)+2e-12)
                self.assertEqual(delta>1e-10,len(G.close(seq))==G.n)

    def test_product_boundary(self):
        A = [p for p in itertools.permutations(range(5)) if not parity(p)]
        G = Group([p+(5+b,5+(1-b)) for p in A for b in (0,1)])
        C = [i for i,g in enumerate(G.elements) if g[5]==6
             and g[:5]!=tuple(range(5)) and G.m[i][i]==G.e]
        # Complete lattice, including the trivial subgroup.
        E=frozenset([G.e]); seen={E:()}; queue=[E]
        for H in queue:
            for c in C:
                if c not in H:
                    gens=seen[H]+(c,); K=G.close(gens)
                    if K not in seen: seen[K]=gens; queue.append(K)
        self.assertEqual(len(G.close(C)),120)
        self.assertEqual(max(G.derived_sizes(H)[1] for H in seen if len(H)<120),2)
        # For A5 x V4, conjugates of (t,b) with a fixed nonidentity b
        # generate A5 x <b>, a proper nonsolvable subgroup of index two.
        self.assertIsNone(G.derived_sizes(G.close(C))[1])
        self.assertEqual(60*4,2*len(G.close(C)))

    def test_gl2_z4_derived_series_independently(self):
        vectors=list(itertools.product(range(4), repeat=2))
        ix={v:i for i,v in enumerate(vectors)}
        els=[]
        for a,b,c,d in itertools.product(range(4),repeat=4):
            if (a*d-b*c)%4 in (1,3):
                els.append(tuple(ix[((a*x+b*y)%4,(c*x+d*y)%4)] for x,y in vectors))
        G=Group(els)
        self.assertEqual(G.derived_sizes(frozenset(range(G.n))),([96,24,4,1],3))

    def test_mismatch_exits_are_not_silent(self):
        # Force a wrong norm while keeping exact group computations unchanged.
        for filename, fn in (('geometric_lemma.py','check'),('generates_G.py','analyze')):
            ns=definitions(filename)
            with contextlib.redirect_stdout(io.StringIO()):
                if fn=='check':
                    leaves=ns['transpositions'](5)
                    ns['opnorm']=lambda _: 0.0
                    with self.assertRaisesRegex(AssertionError,'discrepancies'):
                        ns[fn]('injected failure',leaves,5)
                else:
                    with patch.object(ns['np'].linalg,'svd',return_value=np.zeros(6)):
                        with self.assertRaisesRegex(AssertionError,'discrepancies'):
                            ns[fn]('injected failure',[p+(4+b,4+(1-b)) for p in itertools.permutations(range(4)) for b in (0,1)],6)


if __name__ == '__main__':
    unittest.main(verbosity=2)
