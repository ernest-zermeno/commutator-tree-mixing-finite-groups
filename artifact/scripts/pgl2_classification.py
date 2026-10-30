"""
pgl2_classification.py  —  Phase 1: classify which involution leaf classes of
PGL_2(F_q) / PSL_2(F_q) mix and which are obstructed, via the solvability
dichotomy (P1A/P2A, thm:dichotomy) + Dickson's subgroup classification.

For an involution leaf class C (conjugation-closed) of an ambient G with <C>=G:
  * OBSTRUCTION  iff there is a proper non-solvable C-generated subgroup.
  * else MIXES at exact level  D = max{ derived_length(H) : H<G proper, C-generated }.

Three rows per q:
  (1) G=PGL_2(q), C = ALL involutions          -> interior ones generate PSL_2 -> obstruction.
  (2) G=PGL_2(q), C = EXTERIOR involutions      -> det a non-square; <C>=G.
  (3) G=PSL_2(q), C = involutions -> proper A5 or non-solvable subfield obstruction.

Interior = det a square (lies in PSL_2);  Exterior = det a non-square (PGL_2 minus PSL_2).
A_5's involutions are all even => interior => A_5 is NEVER exterior-generated.
S_4 is exterior-generated iff S_4 not in PSL_2, i.e. q = +-3 mod 8 (Dickson).
Subfield PGL_2(q0): its exterior involutions stay exterior in PGL_2(q) iff [F_q:F_q0] is odd.

SAMPLED witness search; non-detection never certifies absence. Exact levels are ANALYTIC predictions.
Migrated from P2A; maintained canonically in P1. Pure Python finite-field arithmetic.
"""
import itertools, random

# ----------------------------------------------------------------------
# finite field F_q, q = p^k.  Elements are ints 0..q-1 = base-p digits = poly coeffs.
# redux (length k) encodes x^k = sum_i redux[i] x^i ; prime field -> redux=None.
# ----------------------------------------------------------------------
def make_field(p, redux=None):
    if redux is None:
        k, q = 1, p
        add = lambda a, b: (a + b) % p
        mul = lambda a, b: (a * b) % p
        neg = lambda a: (-a) % p
    else:
        k = len(redux); q = p ** k
        def digs(n):
            d = []
            for _ in range(k):
                d.append(n % p); n //= p
            return d
        def undig(d):
            s = 0
            for i in range(k - 1, -1, -1):
                s = s * p + (d[i] % p)
            return s
        def add(a, b):
            da, db = digs(a), digs(b)
            return undig([(da[i] + db[i]) % p for i in range(k)])
        def neg(a):
            return undig([(-x) % p for x in digs(a)])
        def mul(a, b):
            da, db = digs(a), digs(b)
            res = [0] * (2 * k)
            for i in range(k):
                if da[i]:
                    for j in range(k):
                        res[i + j] = (res[i + j] + da[i] * db[j]) % p
            for d in range(2 * k - 1, k - 1, -1):
                if res[d]:
                    c = res[d]; res[d] = 0
                    for i in range(k):
                        res[d - k + i] = (res[d - k + i] + c * redux[i]) % p
            return undig(res[:k])
    elems = list(range(q))
    invtab = {}
    for x in elems:
        if x:
            for y in elems:
                if mul(x, y) == 1:
                    invtab[x] = y; break
    inv = lambda a: invtab[a]
    squares = {mul(x, x) for x in elems if x}
    # multiplicative generator (primitive element)
    def order(g):
        o, c = 1, g
        while c != 1:
            c = mul(c, g); o += 1
        return o
    prim = next(g for g in elems if g and order(g) == q - 1)
    return dict(p=p, k=k, q=q, add=add, mul=mul, neg=neg, inv=inv,
                elems=elems, squares=squares, prim=prim)

# ----------------------------------------------------------------------
# PGL_2(q) on P^1(F_q): points 0..q-1 and q (= infinity). element = canonical matrix (a,b,c,d).
# ----------------------------------------------------------------------
def pgl2(F):
    p, q = F['p'], F['q']; mul, add, inv, neg = F['mul'], F['add'], F['inv'], F['neg']
    INF = q
    def det(M):
        a, b, c, d = M
        return add(mul(a, d), neg(mul(b, c)))
    def normalize(M):
        for x in M:
            if x:
                li = inv(x)
                return tuple(mul(li, y) for y in M)
        return M
    def matmul(M, N):
        a, b, c, d = M; e, f, g, h = N
        return (add(mul(a, e), mul(b, g)), add(mul(a, f), mul(b, h)),
                add(mul(c, e), mul(d, g)), add(mul(c, f), mul(d, h)))
    def act(M, x):
        a, b, c, d = M
        if x == INF:
            return INF if c == 0 else mul(a, inv(c))
        den = add(mul(c, x), d)
        if den == 0:
            return INF
        return mul(add(mul(a, x), b), inv(den))
    def perm(M):
        return tuple(act(M, x) for x in range(q + 1))
    # generators: translation, inversion, scaling by primitive (exterior)
    g = F['prim']
    gens = [normalize((1, 1, 0, 1)), normalize((0, neg(1), 1, 0)), normalize((g, 0, 0, 1))]
    I = normalize((1, 0, 0, 1))
    seen = {I}; frontier = [I]
    while frontier:
        M = frontier.pop()
        for G in gens:
            N = normalize(matmul(M, G))
            if N not in seen:
                seen.add(N); frontier.append(N)
    elems = list(seen)
    info = {}  # perm -> (in_PSL, is_involution)
    permmap = {}
    for M in elems:
        pm = perm(M)
        in_psl = det(M) in F['squares']
        is_inv = (pm != tuple(range(q + 1))) and all(pm[pm[i]] == i for i in range(q + 1))
        info[pm] = (in_psl, is_inv)
        permmap[pm] = M
    return dict(F=F, n=q + 1, order=len(elems), perms=list(info.keys()), info=info)

# ----------------------------------------------------------------------
# permutation group helpers
# ----------------------------------------------------------------------
def pmul(a, b): return tuple(a[b[i]] for i in range(len(b)))
def pinv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)
def comm(a, b): return pmul(pmul(a, b), pmul(pinv(a), pinv(b)))
def closure_capped(gens, cap):
    n = len(gens[0]); e = tuple(range(n)); S = {e}; fr = [e]
    while fr:
        x = fr.pop()
        for g in gens:
            y = pmul(g, x)
            if y not in S:
                S.add(y)
                if len(S) > cap: return None
                fr.append(y)
    return S
def derived_length(H):
    cur = set(H)
    for d in range(8):
        if len(cur) == 1: return d
        cl = list(cur); C = {comm(a, b) for a in cl for b in cl}
        S = set(C); fr = list(C)
        while fr:
            x = fr.pop()
            for g in C:
                y = pmul(g, x)
                if y not in S: S.add(y); fr.append(y)
        if S == cur: return float('inf')   # perfect -> non-solvable
        cur = S
    return float('inf')

# ----------------------------------------------------------------------
# classify (G, C):  search proper C-generated subgroups.
# ----------------------------------------------------------------------
def classify(G_perms, C, name, tries=6000, cap=130, seed=0):
    """Sampled witness search; absence of a witness is never a mixing certificate."""
    rng = random.Random(seed)
    Gset = set(G_perms)
    seen = {}; orders = {}; deepest = 0; nonsolv_order = None
    for _ in range(tries):
        gens = [rng.choice(C) for _ in range(rng.choice([2, 3, 4]))]
        H = closure_capped(gens, cap)
        if H is None or len(H) <= 2 or len(H) == len(Gset):
            continue
        key = frozenset(H)
        if key not in seen:
            seen[key] = derived_length(H)
        dl = seen[key]
        orders.setdefault(len(H), set()).add(dl)
        if dl == float('inf'):
            nonsolv_order = len(H)
        else:
            deepest = max(deepest, dl)
    print(name, 'SAMPLED', {'tries': tries, 'closure_cap': cap, 'seed': seed})
    print('proper witnessed orders:derived_lengths', {o: sorted(v) for o,v in sorted(orders.items())})
    print('OBSTRUCTION WITNESS order '+str(nonsolv_order) if nonsolv_order else
          'NO OBSTRUCTION FOUND; absence and the exact level are not certified by sampling')
    return dict(obstruction=nonsolv_order is not None, witnessed_depth=deepest,
                nonsolv_order=nonsolv_order)

# ----------------------------------------------------------------------
# arithmetic predictions (Dickson)
# ----------------------------------------------------------------------
def A5_in_PSL2(q, p, k):
    return (q % 5 in (1, 4)) or (p == 5)         # q = +-1 mod 5, or power of 5
def S4_exterior(q):  # S4 not in PSL2  <=>  q = +-3 mod 8
    return q % 8 in (3, 5)
def subfield_obstruction(p, k):
    # PGL_2(p^j), j|k proper, p^j>=5 non-solvable, exterior stays exterior iff (k/j) odd
    for j in range(1, k):
        if k % j == 0:
            q0 = p ** j
            if q0 >= 5 and (k // j) % 2 == 1:
                return q0
    return None

# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------
FIELDS = {
    5: (5, None), 7: (7, None), 11: (11, None), 13: (13, None),
    17: (17, None), 19: (19, None), 23: (23, None),
    9: (3, [2, 0]), 25: (5, [2, 0]), 27: (3, [1, 1, 0]),
    125: (5, [4, 4, 0]),  # x^3+x+1 has no root over F5
}

def predicted_levels(p, k):
    """Theorem-derived values, not computational subgroup certificates. None = obstructed."""
    q = p**k
    if q < 5 or p == 2:
        raise ValueError('the theorem requires odd q >= 5')
    ext = None if subfield_obstruction(p,k) else (3 if S4_exterior(q) else 2)
    proper_a5 = q > 5 and A5_in_PSL2(q,p,k)
    proper_nonsolvable_subfield = any(k % j == 0 and p**j >= 5 for j in range(1,k))
    inner = None if proper_a5 or proper_nonsolvable_subfield else (3 if q % 8 in (1,7) else 2)
    return {'all_pgl': None, 'exterior_pgl': ext, 'involutions_psl': inner}

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', type=int, nargs='*', default=[5,7,9,11,13,17])
    parser.add_argument('--tries', type=int, default=6000)
    parser.add_argument('--cap', type=int, default=130)
    args = parser.parse_args()
    for q in args.q:
        if q not in FIELDS or q == 125:
            parser.error('supported sampled group sizes: '+str(sorted(set(FIELDS)-{125})))
        p, redux = FIELDS[q]; F=make_field(p,redux); G=pgl2(F)
        assert G['order'] == q*(q*q-1)
        invol=[g for g in G['perms'] if G['info'][g][1]]
        ext=[g for g in invol if not G['info'][g][0]]
        inn=[g for g in invol if G['info'][g][0]]
        PSL=[g for g in G['perms'] if G['info'][g][0]]
        print('q=',q, 'ANALYTIC theorem predictions:',predicted_levels(p,F['k']))
        assert len(closure_capped(inn, G['order'])) == len(PSL)
        classify(G['perms'],ext,'exterior',args.tries,args.cap)
        classify(PSL,inn,'PSL involutions',args.tries,args.cap)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
