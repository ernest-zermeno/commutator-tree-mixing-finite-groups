"""
(b) Robustness to NON-involution (general matrix) leaves + the determinant layer.

The geometric lemma ( op-norm=1 <=> <leaves> intransitive ) does not use that leaves
are involutions, so the derived-length principle is LEAF-AGNOSTIC. We confirm it with
general (order>2) leaves in the standard rep of GL_3(F_2) (simple -> no abelian layer),
and make the determinant condition (3b) concrete in GL_2(F_5).

General leaf in a real orthogonal rep: P_s = 1/2 (I + rho(s)); since rho(s)^T = rho(s^-1),
D_hat_1 = P_s P_t P_{s^-1} P_{t^-1} and D_hat_2 = D_A D_B D_A^T D_B^T (the involution
case is s^-1 = s).
"""
import itertools
import numpy as np
from itertools import combinations

# ---------- (b1) determinant layer in GL_2(F_5) (cheap; run first) ----------
def det5(M): a, b, c, d = M; return (a*d - b*c) % 5
def mm5(A, B):
    a, b, c, d = A; e, f, g, h = B
    return ((a*e+b*g) % 5, (a*f+b*h) % 5, (c*e+d*g) % 5, (c*f+d*h) % 5)
I5 = (1, 0, 0, 1)
G5 = [M for M in itertools.product(range(5), repeat=4) if det5(M) != 0]
invol5 = [M for M in G5 if M != I5 and mm5(M, M) == I5]
def gen_Fstar(S):
    seen = {1}; fr = [1]
    while fr:
        x = fr.pop()
        for s in S:
            y = (x*s) % 5
            if y not in seen: seen.add(y); fr.append(y)
    return sorted(seen)
det_invol = sorted(set(det5(M) for M in invol5))
print("(b1) GL_2(F_5): F_5^* = {1,2,3,4} cyclic order 4; squares = {1,4}.", flush=True)
print(f"     det(involutions) = {det_invol}  ->  <det> = {gen_Fstar(det_invol)} (= squares, index 2)",
      flush=True)
print(f"     => involution leaves trap det in {gen_Fstar(det_invol)} != F_5^*  -> NOT uniform on GL_2(F_5).",
      flush=True)
print(f"     adding a leaf with det=2: <det> = {gen_Fstar(det_invol+[2])} = F_5^*  -> det mixes.",
      flush=True)
print("     (Involutions have det in {1,4}; full GL_n(F_q) mixing for q-1>2 needs NON-involution leaves.)\n",
      flush=True)

# ---------- (b2) non-involution leaves, GL_3(F_2), STANDARD rep (6-dim) ----------
def gl3f2():
    vecs = [v for v in itertools.product([0,1],repeat=3) if any(v)]; vidx={v:i for i,v in enumerate(vecs)}
    out=set()
    for e in itertools.product([0,1],repeat=9):
        rows=(e[0:3],e[3:6],e[6:9]); per=[]; ok=True
        for v in vecs:
            w=tuple(sum(rows[r][c]*v[c] for c in range(3))%2 for r in range(3))
            if not any(w): ok=False; break
            per.append(vidx[w])
        if ok and len(set(per))==7: out.add(tuple(per))
    return list(out), 7
def pmul(a,b): return tuple(a[b[i]] for i in range(len(a)))
def pinv(p):
    r=[0]*len(p)
    for i,j in enumerate(p): r[j]=i
    return tuple(r)
def is_id(p): return all(p[i]==i for i in range(len(p)))
def comm(a,b): return pmul(pmul(a,b),pmul(pinv(a),pinv(b)))
def order(p):
    x=p; k=1
    while not is_id(x): x=pmul(x,p); k+=1
    return k
def n_fixed(p): return sum(1 for i in range(len(p)) if p[i]==i)
def std_basis(m):
    B=np.zeros((m-1,m))
    for k in range(1,m):
        B[k-1,:k]=1.0; B[k-1,k]=-k; B[k-1]/=np.linalg.norm(B[k-1])
    return B
def pm(p):
    M=np.zeros((len(p),len(p)))
    for i in range(len(p)): M[p[i],i]=1.0
    return M
def opn(M): return float(np.linalg.svd(M,compute_uv=False)[0])

perms,m=gl3f2(); B=std_basis(m)
cand=[p for p in perms if order(p) in (3,4) and n_fixed(p)>=1]
leaves=cand[:14]
print(f"(b2) GL_3(F_2) standard rep (dim {m-1}); NON-involution point-fixing leaves: "
      f"{len(leaves)} of {len(cand)} (orders {sorted(set(order(p) for p in leaves))})", flush=True)
need=set(leaves)|{pinv(p) for p in leaves}
P={s:0.5*(B@pm(s)@B.T+np.eye(m-1)) for s in need}
def D1(s,t): return P[s]@P[t]@P[pinv(s)]@P[pinv(t)]
pairs=[(s,t) for s,t in combinations(leaves,2) if not is_id(comm(s,t))]
ptgt={pr:comm(pr[0],pr[1]) for pr in pairs}
o1=sum(1 for (s,t) in pairs if opn(D1(s,t))>1-1e-9); max1=max((opn(D1(s,t)) for (s,t) in pairs),default=0)
quads=[(i,j) for i in range(len(pairs)) for j in range(i+1,len(pairs))
       if not is_id(comm(ptgt[pairs[i]],ptgt[pairs[j]]))]
D2={}; qtgt={}; o2=0; max2=0.0
for (i,j) in quads:
    (s1,t1),(s2,t2)=pairs[i],pairs[j]
    D=D1(s1,t1)@D1(s2,t2)@D1(s1,t1).T@D1(s2,t2).T
    D2[(i,j)]=D; qtgt[(i,j)]=comm(ptgt[pairs[i]],ptgt[pairs[j]])
    op=opn(D); max2=max(max2,op); o2+=op>1-1e-9
max3=0.0; n3=0; ql=list(quads)
for a in range(len(ql)):
    for b in range(a+1,len(ql)):
        if is_id(comm(qtgt[ql[a]],qtgt[ql[b]])): continue
        n3+=1
        if n3>20000: break
        max3=max(max3,opn(D2[ql[a]]@D2[ql[b]]@D2[ql[a]].T@D2[ql[b]].T))
    if n3>20000: break
print(f"     L1: {len(pairs)} pairs, op-norm=1 in {o1}, max {max1:.4f}", flush=True)
print(f"     L2: {len(quads)} quads, op-norm=1 in {o2}, max {max2:.4f}", flush=True)
print(f"     L3: {min(n3,20000)} configs, max op-norm {max3:.4f}", flush=True)
first=next((L for L,v in [(1,max1),(2,max2),(3,max3)] if v<1-1e-6),None)
print(f"     => GL_3(F_2), NON-involution leaves: first sampled STANDARD-rep contractive level = {first}  "
      f"(std rep only; L3 capped 20k; generation not checked; selected 14 of 98 candidates, not asserted conjugation-closed) (predicted dl(S_4)=3)",
      flush=True)
