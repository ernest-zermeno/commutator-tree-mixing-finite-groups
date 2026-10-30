"""Exact natural-pair algebra and Hecke graph checks; spectra are numerical.

Uses projective Hamilton quaternions over F_q, isomorphic to PGL2(F_q).
No all-pairs theorem is tested or asserted. NumPy eigenvalues corroborate
the external Ramanujan theorem; they are not interval certificates.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
import json
import math
import sys
import numpy as np


def field_group(q):
    if q < 5 or q % 8 not in (3,5) or any(q % p == 0 for p in range(2,math.isqrt(q)+1)):
        raise ValueError('Expected a prime q >= 5 congruent to +/-3 mod 8')
    def norm(a):
        a = tuple(x % q for x in a)
        k = next(x for x in a if x)
        return tuple(x * pow(k,-1,q) % q for x in a)
    def mul(a,b):
        w,x,y,z = a; r,s,t,u = b
        return norm((w*r-x*s-y*t-z*u,w*s+x*r+y*u-z*t,
                     w*t-x*u+y*r+z*s,w*u+x*t-y*s+z*r))
    def inv(a): return norm((a[0],-a[1],-a[2],-a[3]))
    elements = []
    for first in range(4):
        for tail in product(range(q), repeat=3-first):
            a = (0,)*first+(1,)+tail
            if sum(x*x for x in a) % q: elements.append(a)
    return norm,mul,inv,elements


def closure(mul, generators):
    seen = {(1,0,0,0)}
    queue = list(seen)
    for x in queue:
        for s in generators:
            y = mul(x,s)
            if y not in seen: seen.add(y); queue.append(y)
    return seen


def rank_exact(a):
    a = [[Fraction(x) for x in row] for row in a]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[rank],a[pivot] = a[pivot],a[rank]
        d = a[rank][col]; a[rank] = [x/d for x in a[rank]]
        for i in range(rank+1,len(a)):
            d = a[i][col]
            a[i] = [x-d*y for x,y in zip(a[i],a[rank])]
        rank += 1
    return rank


def check(q):
    normalize,mul,inv,elements = field_group(q)
    one = (1,0,0,0)
    a = closure(mul,[normalize(x) for x in ((0,1,0,0),(1,1,0,0),(1,1,1,1))])
    h = normalize((0,1,1,1))
    s = {normalize((0,1,b,c)) for b,c in product((-1,1),repeat=2)}
    b = {mul(mul(h,x),inv(h)) for x in a}
    assert len(elements) == q*(q*q-1) and len(a) == 24 and len(s) == 4
    assert mul(h,h) == one and len(a & b) == 6 and h not in a
    assert len(closure(mul,a | b)) == len(elements)
    assert all({mul(mul(x,t),inv(x)) for t in s} == s for x in a)
    # Reduced norm is determinant under a splitting. The sign is nontrivial on A.
    character = lambda x: pow(sum(v*v for v in x)%q,(q-1)//2,q)
    assert Counter(character(x) for x in a) == {1:12,q-1:12}
    epsilon = 1 if pow(3,(q-1)//2,q)==1 else -1
    generated = closure(mul,s)
    assert len(generated) == len(elements)//(2 if epsilon==1 else 1)
    # Equality of rational group-algebra coefficients, in every representation.
    sa = Counter(mul(t,x) for t in s for x in a)
    ass = Counter(mul(x,t) for t in s for x in a)
    aha = Counter(mul(mul(x,h),y) for x in a for y in a)
    assert sa == ass
    assert aha == Counter({x:6*n for x,n in sa.items()})
    assert len(sa) == 96 and set(sa.values()) == {1}
    # Left cosets xA; right-convolution is well-defined because S is A-invariant.
    coset = {}; reps = []
    for x in elements:
        if x in coset: continue
        idx = len(reps); reps.append(x)
        for y in a: coset[mul(x,y)] = idx
    n = len(reps)
    adjacency = np.zeros((n,n),dtype=np.int64)
    for i,x in enumerate(reps):
        targets = {coset[mul(x,t)] for t in s}
        assert len(targets) == 4 and i not in targets
        adjacency[i,list(targets)] = 1
        for y in a:
            assert {coset[mul(mul(x,y),t)] for t in s} == targets
    assert np.array_equal(adjacency,adjacency.T)
    seen = {0}; queue = [0]
    for i in queue:
        for j in np.flatnonzero(adjacency[i]):
            if int(j) not in seen: seen.add(int(j)); queue.append(int(j))
    assert len(seen) == n
    eig = np.linalg.eigvalsh(adjacency.astype(float)/4)
    assert abs(eig[-1]-1)<1e-12
    c_q = float(max(abs(eig[:-1])))
    assert c_q <= math.sqrt(3)/2+1e-11  # numerical corroboration only
    if q == 5:
        assert np.array_equal(adjacency,np.ones((5,5),dtype=np.int64)-np.eye(5,dtype=np.int64))
        left = np.zeros((5,5),dtype=np.int64)
        for j,x in enumerate(reps):
            for t in s: left[coset[mul(t,x)],j] += 1
        # Exact nullities certify the two distinct full spectra at q=5.
        for val,multiplicity in ((4,1),(3,3),(-1,1)):
            assert 5-rank_exact((left-val*np.eye(5,dtype=np.int64)).tolist()) == multiplicity
        for val,multiplicity in ((4,1),(-1,4)):
            assert 5-rank_exact((adjacency-val*np.eye(5,dtype=np.int64)).tolist()) == multiplicity
        assert not np.array_equal(left,adjacency)
    return dict(q=q,group_order=len(elements),cosets=n,intersection=6,
                legendre_3=epsilon,hecke_norm_numerical=c_q,
                algebra='exact',spectrum='float64 corroboration')


def check_gluing():
    for angle in (0.0001,0.1,0.7,math.pi/2):
        u = np.array([math.cos(angle),math.sin(angle)])
        p = np.diag([1.,0.]); qq = np.outer(u,u); ident = np.eye(2)
        c = abs(math.cos(angle))
        for lam in (0.,0.085,0.8,0.999):
            x = p+lam*1j*(ident-p)
            y = qq-lam*(ident-qq)
            t = x@y@x.conj().T@y.conj().T
            assert np.linalg.norm(t,2)**2 <= 1-(1-lam)*(1-c)/(8*(1+lam))+2e-12


def main():
    check_gluing()
    rows = [check(int(q)) for q in (sys.argv[1:] or [5,11,13,19])]
    print(json.dumps(dict(results=rows,gluing='numerical stress checks passed'),indent=2))


if __name__ == '__main__': main()
