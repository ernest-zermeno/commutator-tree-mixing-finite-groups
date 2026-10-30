"""
Verify the geometric lemma:  op-norm(D_hat_2) = 1  <=>  <leaves> is intransitive on Omega.
If this holds with ZERO mismatches, the spectral half of the derived-length principle
is proved (it is a general fact about products of orthogonal projections + the
permutation-module fixed space).
"""
import itertools
import numpy as np
from itertools import combinations
from derived_length_test import (pgl2, involutions, pmul, is_id, comm,
                                  std_basis, perm_matrix, opnorm, transpositions)

def gl3f2():
    vecs = [v for v in itertools.product([0,1],repeat=3) if any(v)]
    vidx = {v:i for i,v in enumerate(vecs)}; out=set()
    for e in itertools.product([0,1],repeat=9):
        rows=(e[0:3],e[3:6],e[6:9]); per=[]; ok=True
        for v in vecs:
            w=tuple(sum(rows[r][c]*v[c] for c in range(3))%2 for r in range(3))
            if not any(w): ok=False;break
            per.append(vidx[w])
        if ok and len(set(per))==7: out.add(tuple(per))
    return list(out),7

def orbits(leaf_perms, m):
    parent=list(range(m))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for p in leaf_perms:
        for i in range(m):
            a,b=find(i),find(p[i])
            if a!=b: parent[a]=b
    return len(set(find(i) for i in range(m)))

def n_fixed(p): return sum(1 for i in range(len(p)) if p[i]==i)

def check(name, leaves, m):
    B=std_basis(m); Pc={p:0.5*(B@perm_matrix(p)@B.T+np.eye(m-1)) for p in leaves}
    pairs=[(s,t) for s,t in combinations(leaves,2) if not is_id(comm(s,t))]
    ptgt={pr:comm(pr[0],pr[1]) for pr in pairs}
    quads=[(i,j) for i in range(len(pairs)) for j in range(i+1,len(pairs))
           if not is_id(comm(ptgt[pairs[i]],ptgt[pairs[j]]))]
    mism=0; n_op1=0; n_intr=0; n_op1_fixespoint=0
    for (i,j) in quads:
        (s1,t1),(s2,t2)=pairs[i],pairs[j]
        D=np.linalg.matrix_power(Pc[s1]@Pc[t1],2)@np.linalg.matrix_power(Pc[s2]@Pc[t2],2)@\
          np.linalg.matrix_power(Pc[t1]@Pc[s1],2)@np.linalg.matrix_power(Pc[t2]@Pc[s2],2)
        op1 = opnorm(D) > 1-1e-9
        intr = orbits([s1,t1,s2,t2], m) > 1
        if op1 != intr: mism += 1
        if op1:
            n_op1 += 1
            # does <leaves> fix a COMMON point? (orbit of size 1 exists)
            # count common fixed points of all 4 leaves:
            cf = sum(1 for x in range(m) if s1[x]==x and t1[x]==x and s2[x]==x and t2[x]==x)
            if cf >= 1: n_op1_fixespoint += 1
        if intr: n_intr += 1
    print(f"  {name}: {len(quads)} reduced quads | op-norm=1: {n_op1} | intransitive: {n_intr} | "
          f"MISMATCHES: {mism}  ({'LEMMA HOLDS' if mism==0 else 'FAIL'})")
    print(f"      of the {n_op1} op-norm-1 configs, {n_op1_fixespoint} fix a COMMON POINT "
          f"({'all' if n_op1_fixespoint==n_op1 else 'NOT all -> intransitive is strictly weaker'})")

    if mism:
        raise AssertionError(f"criterion discrepancies: {mism}")

print("Geometric lemma:  op-norm(D2)=1  <=>  <leaves> intransitive on Omega")
print("="*70)
p5,_=pgl2(5)
p7,_=pgl2(7)
g3,_=gl3f2()
check("S_5 / 5pt (transp.)", transpositions(5), 5)
check("S_5 / 6pt = PGL2(F5) (invol.)", involutions(p5), 6)
check("GL_3(F2) / 7pt (invol.)", involutions(g3), 7)
check("PGL2(F7) / 8pt (split invol.)", [p for p in involutions(p7) if n_fixed(p)==2], 8)
