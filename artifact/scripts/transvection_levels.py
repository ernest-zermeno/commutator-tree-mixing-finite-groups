"""Exact transvection-generated subgroup census in SL3(3)=PSL3(3).

Enumerates the closure lattice up to conjugacy and computes derived series in
integer arithmetic modulo 3. Cross-checks the structural proof in P1: level 4,
attained at order 216, and no proper non-solvable transvection-generated subgroup.
Adapted from the independent audit's sl33_transvections.py (26 September 2026);
the original and hashes are recorded in notes/open-problems-implementation-2026-09-26.
No third-party dependencies. Runtime is about half a minute.
"""
import time
p=3
def mm(a,b): return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3))%p for i in range(3) for j in range(3))
e=(1,0,0,0,1,0,0,0,1)
def inv(a):
    c=[[(a[3*((j+1)%3)+((i+1)%3)]*a[3*((j+2)%3)+((i+2)%3)]-a[3*((j+1)%3)+((i+2)%3)]*a[3*((j+2)%3)+((i+1)%3)])%p for j in range(3)] for i in range(3)]
    r=tuple(c[i][j] for i in range(3) for j in range(3)); assert mm(a,r)==e; return r
gens=[]
for i in range(3):
    for j in range(3):
        if i!=j:
            m=list(e); m[3*i+j]=1; gens.append(tuple(m))
def closure(gs):
    S={e}; fr=[e]
    while fr:
        nx=[]
        for x in fr:
            for g in gs:
                y=mm(x,g)
                if y not in S: S.add(y); nx.append(y)
        fr=nx
    return S
G=closure(gens); assert len(G)==5616
def is_transvection(g):
    # g != e, (g-I)^2 = 0, rank(g-I)=1  <=> g-I has rank 1 and (g-I)^2=0
    n=tuple((g[k]-e[k])%p for k in range(9))
    if g==e: return False
    n2=tuple(sum(n[3*i+k]*n[3*k+j] for k in range(3))%p for i in range(3) for j in range(3))
    if any(n2): return False
    rows=[n[0:3],n[3:6],n[6:9]]; nz=[r for r in rows if any(r)]
    # rank 1: all nonzero rows proportional
    r0=nz[0]
    for r in nz[1:]:
        ok=False
        for c in range(1,p):
            if tuple((c*x)%p for x in r0)==r: ok=True
        if not ok: return False
    return True
T=[g for g in G if is_transvection(g)]
print("transvections:",len(T),"(expected (27-1)(9-1)/2=104)")
def small_gens(H):
    gs=[]; S={e}
    for x in H:
        if x not in S:
            gs.append(x); S=closure(gs)
            if len(S)==len(H): break
    return gs
def derived_series(H):
    cur=set(H); orders=[len(cur)]
    while len(cur)>1:
        gs=small_gens(cur); comm=set()
        for x in gs:
            xi=inv(x)
            for y in cur: comm.add(mm(mm(x,y),mm(xi,inv(y))))
        D=closure(comm)
        changed=True
        while changed:
            changed=False; new=set()
            for g in gs:
                gi=inv(g)
                for x in list(D):
                    y=mm(mm(g,x),gi)
                    if y not in D: new.add(y)
            if new: D=closure(small_gens(D)+list(new)); changed=True
        if len(D)==len(cur): return orders+['nonsolvable']
        cur=D; orders.append(len(cur))
    return orders
t0=time.time(); seen=set(); queue=[]; results=[]
def orbit_mark(H):
    Hf=frozenset(H)
    if Hf in seen: return False
    orb={Hf}; fr=[Hf]
    while fr:
        nx=[]
        for X in fr:
            for g in gens:
                gi=inv(g); Y=frozenset(mm(mm(g,x),gi) for x in X)
                if Y not in orb: orb.add(Y); nx.append(Y)
        fr=nx
    seen.update(orb); return True
for t in T:
    H=closure([t])
    if orbit_mark(H): queue.append((H,[t]))
while queue:
    H,gs=queue.pop()
    if len(H)==5616: results.append((5616,'G')); continue
    results.append((len(H),derived_series(H)))
    for t in T:
        if t in H: continue
        H2=closure(gs+[t])
        if orbit_mark(H2): queue.append((H2,gs+[t]))
proper=[r for r in results if r[0]<5616]
print("classes:",len(results),"proper non-solvable:",[r for r in proper if r[1][-1]=='nonsolvable'])
mx=max(len(r[1])-1 for r in proper if r[1][-1]!='nonsolvable')
print("max derived length of proper transvection-generated subgroups:",mx,"attained by orders",sorted({r[0] for r in proper if len(r[1])-1==mx}))
print("all proper classes (order, series):",sorted(set((r[0],tuple(r[1])) for r in proper)))
print(f"{time.time()-t0:.0f}s")

assert len(T)==104 and all(g in T for g in gens)
assert len(results)==8 and all(x[1][-1]!="nonsolvable" for x in proper)
assert mx==4
assert sorted(set((x[0],tuple(x[1])) for x in proper)) == [(3,(3,1)),(9,(9,1)),(24,(24,8,2,1)),(27,(27,3,1)),(216,(216,72,18,9,1))]
print("TRANSVECTION LEVEL CERTIFICATE: PASS")
