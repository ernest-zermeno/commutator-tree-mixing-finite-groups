"""Independent integer group computations; imports no repository verification code.

Checks complete ordered level-2 universes, class-generated subgroup lattices,
derived lengths, and split/external involution distinction. NumPy is used only
for a vectorised integer lookup, not for group or solvability decisions.
"""
from collections import Counter, deque
from itertools import permutations, product
from pathlib import Path
import json
import numpy as np



def mul(p, q):
    return tuple(p[x] for x in q)


def parity(p):
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p))) % 2


def pgl(q):
    answer = set()
    for a,b,c,d in product(range(q), repeat=4):
        if (a*d-b*c) % q == 0:
            continue
        p = []
        for x in range(q+1):
            num,den = (a,c) if x == q else ((a*x+b)%q,(c*x+d)%q)
            p.append(q if den == 0 else num*pow(den,-1,q)%q)
        answer.add(tuple(p))
    return answer


def gl32():
    answer = set()
    for cols in product(range(1,8), repeat=3):
        images = [0]
        for x in range(1,8):
            y = 0
            for i in range(3):
                if x & (1 << i):
                    y ^= cols[i]
            images.append(y)
        if len(set(images)) == 8:
            answer.add(tuple(x-1 for x in images[1:]))
    return answer


class Group:
    def __init__(self, elements):
        self.elements = sorted(elements)
        self.degree = len(self.elements[0])
        self.n = len(self.elements)
        self.index = {p:i for i,p in enumerate(self.elements)}
        self.e = self.index[tuple(range(self.degree))]
        self.m = [[self.index[mul(p,q)] for q in self.elements] for p in self.elements]
        self.inv = [row.index(self.e) for row in self.m]
        self.comm = [[self.m[self.m[x][y]][self.m[self.inv[x]][self.inv[y]]]
                      for y in range(self.n)] for x in range(self.n)]
        self.closures = {}

    def close(self, gens):
        key = tuple(sorted(set(gens)))
        if key in self.closures:
            return self.closures[key]
        seen, stack = {self.e}, [self.e]
        while stack:
            x = stack.pop()
            for s in key:
                y = self.m[x][s]
                if y not in seen:
                    seen.add(y); stack.append(y)
        result = frozenset(seen)
        self.closures[key] = result
        return result

    def derived_sizes(self, H):
        sizes = [len(H)]
        while len(H) > 1:
            K = self.close({self.comm[x][y] for x in H for y in H})
            sizes.append(len(K))
            if K == H:
                return sizes, None
            H = K
        return sizes, len(sizes)-1

    def orbit_count(self, H):
        unseen = set(range(self.degree)); count = 0
        while unseen:
            x = min(unseen)
            unseen.difference_update(self.elements[g][x] for g in H)
            count += 1
        return count

    def analyze(self, name, C):
        assert len(self.close(C)) == self.n
        # Enumerate every C-generated subgroup by adding a leaf generator.
        E = frozenset([self.e]); lattice = {E:()}; queue = deque([E])
        while queue:
            H = queue.popleft()
            for c in C:
                if c in H:
                    continue
                gens = lattice[H] + (c,)
                K = self.close(gens)
                if K not in lattice:
                    lattice[K] = gens; queue.append(K)
        profiles = Counter()
        for H in lattice:
            if len(H) < self.n:
                sizes, dl = self.derived_sizes(H)
                profiles[(len(H), dl)] += 1
        # Full ordered pairs, full ordered pairs-of-pairs (not combinations).
        pairs = [(s,t,self.comm[s][t]) for s in C for t in C
                 if self.comm[s][t] != self.e]
        counts = Counter(); bad_orbits = Counter(); cache = {}
        comm = np.asarray(self.comm, dtype=np.int32)
        targets = np.asarray([p[2] for p in pairs],dtype=np.int32)
        for s,t,a in pairs:
            for j in np.flatnonzero(comm[a,targets] != self.e):
                u,v,_ = pairs[j]
                key = tuple(sorted(set((s,t,u,v))))
                if key not in cache:
                    H = self.close(key)
                    cache[key] = (len(H), self.orbit_count(H))
                order,orbits = cache[key]
                counts[order] += 1
                if orbits > 1:
                    bad_orbits[order] += 1
        row = dict(name=name, order=self.n, leaf_count=len(C),
                   class_generated_subgroups_including_trivial=len(lattice),
                   proper_lattice=[dict(order=o, derived_length=d, count=c)
                                   for (o,d),c in sorted(profiles.items(), key=lambda x:(x[0][0],str(x[0][1])))],
                   valid_ordered_L1=len(pairs),valid_ordered_L2=sum(counts.values()),
                   ordered_L2_leaf_subgroup_orders=dict(sorted(counts.items())),
                   ordered_L2_intransitive=dict(sorted(bad_orbits.items())))
        print(json.dumps(row),flush=True)
        return row


if __name__ == '__main__':
    cases = [
        ('S5_transpositions',list(permutations(range(5))),'transpositions',3240,2160),
        ('A5_involutions',[p for p in permutations(range(5)) if not parity(p)],'all',29640,0),
        ('GL3_F2_involutions',gl32(),'all',105168,22176),
        ('PGL2_5_all_involutions',pgl(5),'all',219360,37560),
        ('PGL2_7_external_involutions',pgl(7),'external',429744,0)]
    for name, els, kind, total, proper in cases:
        G=Group(els)
        C=[i for i,p in enumerate(G.elements) if i!=G.e and G.m[i][i]==G.e
           and (kind!='transpositions' or sum(j!=x for j,x in enumerate(p))==2)
           and (kind!='external' or parity(p)==1)]
        row=G.analyze(name,C)
        assert row['valid_ordered_L2']==total, row
        assert sum(v for k,v in row['ordered_L2_leaf_subgroup_orders'].items() if k<G.n)==proper, row
    print('ORDERED GROUP CHECKS: PASS')
