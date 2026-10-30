"""Small-group checks of arbitrary-order stability; analytic proofs are in P1.

The group operations/diameters are exact; NumPy norm checks are corroboration,
not certified upper bounds for a family of groups.
"""
import itertools
import unittest
import numpy as np
from verification_groups import Group


def augmentation(n):
    b = np.zeros((n, n-1))
    for k in range(1, n):
        b[:k, k-1] = 1
        b[k, k-1] = -k
        b[:, k-1] /= np.linalg.norm(b[:, k-1])
    return b


def diameter(group, generators):
    steps = set(generators)
    steps.update(next(t for t in range(group.n) if group.m[s][t] == group.e)
                 for s in generators)
    distances = {group.e: 0}
    frontier = [group.e]
    for x in frontier:
        for s in steps:
            y = group.m[s][x]
            if y not in distances:
                distances[y] = distances[x]+1
                frontier.append(y)
    return distances, max(distances.values())


class QuantitativeTests(unittest.TestCase):
    def test_arbitrary_order_stability_and_fixed_spaces(self):
        cases = [([tuple((i+j)%n for i in range(n)) for j in range(n)], None)
                 for n in (3, 4)]
        cases.append((list(itertools.permutations(range(3))), None))
        rng = np.random.default_rng(20260923)
        tested_orders = set()
        for elements, _ in cases:
            g = Group(elements)
            b = augmentation(g.n)
            matrices = []
            regular = []
            for s in range(g.n):
                r = np.zeros((g.n,g.n))
                for j in range(g.n): r[g.m[s][j],j] = 1
                u = b.T @ r @ b
                regular.append(u)
                matrices.append((np.eye(g.n-1)+u)/2)
                x, order = s, 1
                while x != g.e: x = g.m[s][x]; order += 1
                tested_orders.add(order)
            for length in (1, 2, 3):
                for seq in itertools.product(range(g.n), repeat=length):
                    if all(s == g.e for s in seq): continue
                    h, diam = diameter(g, seq)
                    fixed = sum(regular[s] for s in h)/len(h)
                    t = np.eye(g.n-1)
                    for s in seq: t = matrices[s] @ t
                    v = rng.normal(size=g.n-1)+1j*rng.normal(size=g.n-1)
                    v /= np.linalg.norm(v)
                    energy = max(0., 1-np.linalg.norm(t@v)**2)
                    self.assertLessEqual(np.linalg.norm(v-fixed@v),
                                         2*diam*np.sqrt(length*energy)+2e-12)
                    if len(h) == g.n:
                        self.assertLessEqual(np.linalg.norm(t,2)**2,
                                             1-1/(4*length*diam**2)+2e-12)
                    else:
                        self.assertAlmostEqual(np.linalg.norm(t,2),1.,places=11)
        self.assertTrue({2,3,4}.issubset(tested_orders))

    def test_support_obstruction_for_an_actual_distribution(self):
        g = Group(itertools.permutations(range(3)))
        leaf = next(s for s in range(g.n) if s != g.e and g.m[s][s] == g.e)
        p = {g.e: 1.}
        for _ in range(4):
            nxt = {}
            for x, weight in p.items():
                for y in (x, g.m[x][leaf]): nxt[y] = nxt.get(y,0.)+weight/2
            p = nxt
        tv = sum(abs(p.get(x,0.)-1/g.n) for x in range(g.n))/2
        self.assertAlmostEqual(tv, 1-len(p)/g.n)


if __name__ == '__main__':
    unittest.main()
