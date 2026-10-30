"""
Closes finding #1: the full-mixing criterion is <leaves> = G (across ALL irreps),
of which the standard-rep "intransitive" test is only the projection to one irrep.

Per valid level-2 config we check:
    regular-rep 2nd singular value < 1   <=>   <4 leaves> = G          (the real lemma)
and compare to the standard-rep op-norm=1 (= <leaves> intransitive). The gap
(proper TRANSITIVE leaf-subgroups: not G, yet transitive) is exactly what the
standard-rep lemma misses — e.g. leaves trapped in proper transitive subgroups inside
S_5 (including A_5 and S_4-type closures).
"""
import itertools
import numpy as np
from itertools import combinations

def pmul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def pinv(p):
    r = [0]*len(p)
    for i, j in enumerate(p): r[j] = i
    return tuple(r)
def is_id(p): return all(p[i] == i for i in range(len(p)))
def comm(a, b): return pmul(pmul(a, b), pmul(pinv(a), pinv(b)))
def involutions(perms): return [p for p in perms if not is_id(p) and is_id(pmul(p, p))]

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
def pgl2(q):
    pts=list(range(q))+['inf']; idx={pt:i for i,pt in enumerate(pts)}
    def ap(a,b,c,d,pt):
        if pt=='inf': return 'inf' if c%q==0 else (a*pow(c%q,-1,q))%q
        num=(a*pt+b)%q; den=(c*pt+d)%q
        return 'inf' if den==0 else (num*pow(den,-1,q))%q
    S=set()
    for a in range(q):
        for b in range(q):
            for c in range(q):
                for d in range(q):
                    if (a*d-b*c)%q==0: continue
                    S.add(tuple(idx[ap(a,b,c,d,pt)] for pt in pts))
    return list(S), q+1

def stdbasis(m):
    B=np.zeros((m-1,m))
    for k in range(1,m): B[k-1,:k]=1; B[k-1,k]=-k; B[k-1]/=np.linalg.norm(B[k-1])
    return B
def pm(p,m):
    M=np.zeros((m,m))
    for i in range(m): M[p[i],i]=1.0
    return M

def analyze(name, perms, m):
    N=len(perms); gidx={g:i for i,g in enumerate(perms)}; ID=tuple(range(m))
    invs=involutions(perms)
    def closure(gens):
        seen={ID}; fr=[ID]
        while fr:
            x=fr.pop()
            for s in gens:
                y=pmul(x,s)
                if y not in seen: seen.add(y); fr.append(y)
        return len(seen)
    def regproj(s):
        P=np.zeros((N,N))
        for c,g in enumerate(perms): P[gidx[pmul(s,g)],c]=1.0
        return 0.5*(P+np.eye(N))
    Preg={s:regproj(s) for s in invs}
    B=stdbasis(m); Pstd={s:0.5*(B@pm(s,m)@B.T+np.eye(m-1)) for s in invs}
    def blk(Pc,s,t): return np.linalg.matrix_power(Pc[s]@Pc[t],2)
    def D2(Pc,s1,t1,s2,t2): return blk(Pc,s1,t1)@blk(Pc,s2,t2)@blk(Pc,t1,s1)@blk(Pc,t2,s2)
    pairs=[(s,t) for s,t in combinations(invs,2) if not is_id(comm(s,t))]
    ptgt={pr:comm(*pr) for pr in pairs}
    quads=[(i,j) for i in range(len(pairs)) for j in range(i+1,len(pairs))
           if not is_id(comm(ptgt[pairs[i]],ptgt[pairs[j]]))]
    n_intr=n_notG=n_genG=mism=n_tp=0
    for (i,j) in quads:
        (s1,t1)=pairs[i]; (s2,t2)=pairs[j]
        genG = closure([s1,t1,s2,t2])==N
        sv2 = float(np.linalg.svd(D2(Preg,s1,t1,s2,t2),compute_uv=False)[1])
        contracts_all = sv2 < 1-1e-9
        intr = float(np.linalg.svd(D2(Pstd,s1,t1,s2,t2),compute_uv=False)[0]) > 1-1e-9
        n_intr += intr; n_notG += (not contracts_all); n_genG += genG
        mism += (contracts_all != genG); n_tp += (not genG) and (not intr)
    print(f"\n{name}: |G|={N}, {len(quads)} valid reduced level-2 configs", flush=True)
    print(f"  standard-rep op-norm=1 (intransitive):           {n_intr}", flush=True)
    print(f"  regular-rep 2nd-sv=1  (<leaves> != G):           {n_notG}", flush=True)
    print(f"  <leaves> generate G:                             {n_genG}", flush=True)
    print(f"  proper-TRANSITIVE leaf-groups (std rep MISSES):  {n_tp}", flush=True)
    print(f"  mismatches  (2nd-sv<1) == (<leaves>=G):          {mism}  "
          f"({'LEMMA HOLDS' if mism==0 else 'FAIL'})", flush=True)

    if mism:
        raise AssertionError(f"criterion discrepancies: {mism}")

if __name__ == "__main__":
    g3,_=gl3f2(); analyze("GL_3(F_2) [maximal subgrps S_4 (intr) + 7:3 (transitive, no involutions)]", g3, 7)
    p5,_=pgl2(5); analyze("S_5 on 6pts = PGL_2(F_5) [proper transitive subgroups have involution leaves]", p5, 6)
