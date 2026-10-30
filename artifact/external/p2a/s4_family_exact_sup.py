r"""
s4_family_exact_sup.py
======================
EXACT per-family TRACE of the S4-compression for the rank-one frontier (ii), with the
multiplicity dim rho^{S4} and the resulting norm LOWER bound, by representation family
(principal / Steinberg / cuspidal), sup over h in N_G(S3)\A (verified <A,B>=G) and over
all regular characters.

WHAT THIS COMPUTES (and what it does NOT):
    T_rho(h) := (1/24) sum_{a in S4} chi_rho(a h) = tr( P_A^rho rho(h) )   [TRACE of the compression]
This is NOT the spherical function and NOT ||P_A rho(h) P_A|| unless dim rho^{S4} = 1.
For k = dim rho^{S4} >= 1 the only clean inequality is
    ||P_A rho(h) P_A||  >=  |T_rho(h)| / k          (NOT ||.|| >= |T|).
So |T| can EXCEED the norm (extreme: rho(h)|_{rho^{S4}} = I gives |T| = k, norm = 1).
[Correction over a previous version that mislabeled |T| as the norm / a lower bound for it.]

Multiplicities via character average over the subgroup A = S4 (identity handled separately,
chi_rho(e) = dim rho):
    dim rho^{S4} = (1/24) sum_{a in A} chi_rho(a)  =  (1/24)( dim rho + sum_{a in A, a!=e} chi_rho(a) ).
NOTE: this is the correct route; #S4-orbits on P^1 only gives dim(1 + St)^{S4}, not principal/cuspidal.

PGL2(q) character values (q odd) on the relevant classes:
    family      identity   unipotent   split(ratio r)        elliptic(param xi)
    principal   q+1        1           chi(r)+chi(r^-1)       0
    Steinberg   q          0           1                     -1
    cuspidal    q-1        -1          0                      -(theta(xi)+theta(xi^-1))

Regular characters: principal 2m !≡ 0 (mod q-1); cuspidal 2m !≡ 0 (mod q+1).

CONCLUSION THIS SUPPORTS: for large q, dim rho^{S4} ~ dim rho / 24 ~ q/24, so multiplicity > 1
is GENERIC. The scalar Bessel/Soto-Andrade reduction (rem:bessel) resolves ONLY multiplicity-one
blocks; the general frontier of (ii) is the operator NORM of a growing-dimension matrix
compression (as rem:spherical already states). This trace scan is a STRUCTURE/multiplicity
diagnostic + a norm lower bound, NOT the norm itself.

Self-contained (pure arithmetic + math.cos). Reuses the exact arithmetic of s4_eigenvalue_params.py.
"""
import math
import random
from collections import Counter


def build(q):
    inv = [0] * q
    for a in range(1, q):
        inv[a] = pow(a, q - 2, q)

    def mm(M, N):
        a, b, c, d = M
        e, f, g_, h = N
        return ((a * e + b * g_) % q, (a * f + b * h) % q,
                (c * e + d * g_) % q, (c * f + d * h) % q)

    def det(M):
        a, b, c, d = M
        return (a * d - b * c) % q

    def norm(M):
        for v in M:
            if v % q != 0:
                s = inv[v % q]
                return tuple((x * s) % q for x in M)
        return M

    def minv(M):
        a, b, c, d = M
        di = inv[det(M)]
        return norm(((d * di) % q, (-b * di) % q, (-c * di) % q, (a * di) % q))

    I = norm((1, 0, 0, 1))
    sq = set((x * x) % q for x in range(1, q))
    g_ns = next(c for c in range(2, q) if c not in sq)

    gens = [norm((1, 1, 0, 1)), norm((0, q - 1, 1, 0)), norm((g_ns, 0, 0, 1))]
    G = {I}
    fr = [I]
    while fr:
        x = fr.pop()
        for gg in gens:
            y = norm(mm(gg, x))
            if y not in G:
                G.add(y)
                fr.append(y)
    assert len(G) == q ** 3 - q, f"q={q}: |G|={len(G)} != q^3-q"

    def fmul(X, Y):
        a, b = X
        c, d = Y
        return ((a * c + b * d * g_ns) % q, (a * d + b * c) % q)

    def forder(X):
        o, Z = 1, X
        while Z != (1, 0):
            Z = fmul(Z, X)
            o += 1
        return o

    fgen = None
    for a in range(q):
        for b in range(q):
            if (a, b) != (0, 0) and forder((a, b)) == q * q - 1:
                fgen = (a, b)
                break
        if fgen:
            break
    dl = {}
    Z = (1, 0)
    for k in range(q * q - 1):
        dl[Z] = k
        Z = fmul(Z, fgen)

    gq = next(c for c in range(2, q) if len(set(pow(c, k, q) for k in range(q - 1))) == q - 1)
    dlq = {pow(gq, k, q): k for k in range(q - 1)}
    sqrt = {(x * x) % q: x for x in range(q)}

    return dict(q=q, I=I, mm=mm, det=det, norm=norm, minv=minv, sq=sq, g=g_ns,
                G=list(G), Gset=G, dl=dl, dlq=dlq, sqrt=sqrt)


def comm(B, x, y):
    mm, norm, minv = B["mm"], B["norm"], B["minv"]
    return norm(mm(mm(x, y), mm(minv(x), minv(y))))


def closure(B, gens, cap=None):
    mm, norm, I = B["mm"], B["norm"], B["I"]
    S = {I}
    fr = [I]
    while fr:
        x = fr.pop()
        for g in gens:
            y = norm(mm(g, x))
            if y not in S:
                S.add(y)
                if cap is not None and len(S) > cap:
                    return None
                fr.append(y)
    return S


def order(B, a):
    mm, norm, I = B["mm"], B["norm"], B["I"]
    o, c = 1, a
    while c != I:
        c = norm(mm(c, a))
        o += 1
    return o


def derived_length(B, H, cap=6):
    cur = set(H)
    for d in range(cap):
        if len(cur) == 1:
            return d
        cl = list(cur)
        C = {comm(B, a, b) for a in cl for b in cl}
        nxt = closure(B, C)
        if nxt == cur:
            return float("inf")
        cur = nxt
    return float("inf")


def find_S4(B, seed=1, tries=200000):
    I, mm, norm, sq = B["I"], B["mm"], B["norm"], B["sq"]
    ext = [M for M in B["G"] if M != I and norm(mm(M, M)) == I and B["det"](M) not in sq]
    random.seed(seed)
    for _ in range(tries):
        H = closure(B, [random.choice(ext) for _ in range(random.choice([3, 4]))], cap=25)
        if (H and len(H) == 24
                and sum(1 for x in H if x != I and norm(mm(x, x)) == I) == 9
                and derived_length(B, H) == 3):
            return frozenset(H), set(ext)
    return None, set(ext)


def classify(B, elems):
    """Classify matrices by PGL2 conjugacy type. Returns (n_id, u, splits, ells).
       n_id = #identity, u = #unipotent, splits = [dl(r) mod q-1], ells = [dl(xi) mod q+1]."""
    q, mm, norm, sq, I = B["q"], B["mm"], B["norm"], B["sq"], B["I"]
    sqrt, inv_g, dl, dlq = B["sqrt"], B["g"], B["dl"], B["dlq"]
    n_id = u = 0
    ells, spls = [], []
    for M in elems:
        M = norm(M)
        if M == I:
            n_id += 1
            continue
        aa, bb, cc, dd = M
        t = (aa + dd) % q
        dt = (aa * dd - bb * cc) % q
        D = (t * t - 4 * dt) % q
        if D == 0:
            u += 1
        elif D in sq:
            sd = sqrt[D]
            den = (t - sd) % q
            if den:
                r = ((t + sd) * pow(den, q - 2, q)) % q
                if r:
                    spls.append(dlq[r] % (q - 1))
        else:
            s = sqrt[(D * pow(inv_g, q - 2, q)) % q]
            x = ((t * pow(2, q - 2, q)) % q, (s * pow(2, q - 2, q)) % q)
            ells.append(dl[x] % (q + 1))
    return n_id, u, spls, ells


def fam_value(fam, q, m, n_id, u, spls, ells):
    """(1/24) sum chi_rho over the classified element set, for family `fam` and character index m.
       With n_id=0 (coset A*h) this is the TRACE; with the identity present (A) it gives 24*dim."""
    if fam == "cusp":
        n1 = q + 1
        s = n_id * (q - 1) - u
        s += -sum(2.0 * math.cos(2 * math.pi * m * xi / n1) for xi in ells)
        return s / 24.0
    if fam == "princ":
        n2 = q - 1
        s = n_id * (q + 1) + u
        s += sum(2.0 * math.cos(2 * math.pi * m * r / n2) for r in spls)
        return s / 24.0
    # Steinberg (no character index)
    s = n_id * q + len(spls) - len(ells)
    return s / 24.0


def family_report(B, A, h):
    q = B["q"]
    n0, u0, sp0, el0 = classify(B, A)                 # over A (identity present): dimensions
    nh, uh, sph, elh = classify(B, [B["mm"](a, h) for a in A])  # over A*h: traces
    out = {}
    # cuspidal
    n1 = q + 1
    best_tr = best_lb = 0.0
    for m in range(n1):
        if (2 * m) % n1 == 0:
            continue
        tr = abs(fam_value("cusp", q, m, 0, uh, sph, elh))
        dim = round(fam_value("cusp", q, m, n0, u0, sp0, el0))
        best_tr = max(best_tr, tr)
        if dim >= 1:
            best_lb = max(best_lb, tr / dim)
    out["cusp"] = (best_tr, best_lb)
    # principal
    n2 = q - 1
    best_tr = best_lb = 0.0
    for m in range(n2):
        if (2 * m) % n2 == 0:
            continue
        tr = abs(fam_value("princ", q, m, 0, uh, sph, elh))
        dim = round(fam_value("princ", q, m, n0, u0, sp0, el0))
        best_tr = max(best_tr, tr)
        if dim >= 1:
            best_lb = max(best_lb, tr / dim)
    out["princ"] = (best_tr, best_lb)
    # Steinberg
    tr = abs(fam_value("stein", q, 0, 0, uh, sph, elh))
    dim = round(fam_value("stein", q, 0, n0, u0, sp0, el0))
    out["stein"] = (tr, (tr / dim if dim >= 1 else 0.0))
    return out


def scan_q(q, seed=1):
    B = build(q)
    I, mm, norm, minv, Gset = B["I"], B["mm"], B["norm"], B["minv"], B["Gset"]
    A, extset = find_S4(B, seed=seed)
    if A is None:
        return None
    extA = [x for x in A if x in extset]
    S3 = None
    for s in extA:
        for t in extA:
            if s != t and order(B, norm(mm(s, t))) == 3:
                S3 = frozenset(closure(B, [s, t], cap=6))
                break
        if S3:
            break

    agg = {"cusp": (0.0, 0.0), "princ": (0.0, 0.0), "stein": (0.0, 0.0)}
    seen = set()
    n_cos = 0
    all_genG = True
    for h in B["G"]:
        if h in A:
            continue
        hi = minv(h)
        if frozenset(norm(mm(mm(h, x), hi)) for x in S3) != S3:
            continue
        Bgrp = frozenset(norm(mm(mm(h, a), hi)) for a in A)
        if Bgrp == A:
            continue
        rep = min(norm(mm(a, h)) for a in A)
        if rep in seen:
            continue
        seen.add(rep)
        # verify <A,B> = G
        gen_clo = closure(B, list(A | Bgrp))
        genG = (len(gen_clo) == len(Gset))
        all_genG = all_genG and genG
        if not genG:
            continue
        n_cos += 1
        rep_fam = family_report(B, A, h)
        for fam in agg:
            agg[fam] = (max(agg[fam][0], rep_fam[fam][0]),
                        max(agg[fam][1], rep_fam[fam][1]))

    def small_d(n):
        for d in range(3, 13):
            if n % d == 0:
                return d
        return None

    norm_lb = max(agg[f][1] for f in agg)
    return dict(q=q, n_cos=n_cos, agg=agg, norm_lb=norm_lb,
                d_cusp=small_d(q + 1), d_prin=small_d(q - 1), all_genG=all_genG)


def main():
    qs = [11, 13, 19, 29, 37, 43, 53, 59]
    print("=== TRACE of S4-compression by family (NOT the norm); norm LOWER bound |tr|/dim ===")
    print("    T_rho(h)=(1/24)sum_a chi_rho(ah)=tr(compression);  ||P_A rho(h) P_A|| >= |T|/dim")
    print("    sup over h in N(S3)\\A with <A,B>=G, and over all regular characters.\n")
    print(f"{'q':>4} {'mod8':>4} {'#h':>3} "
          f"{'cusp_tr':>8} {'cusp_lb':>8} {'prin_tr':>8} {'prin_lb':>8} "
          f"{'stein_tr':>8} {'NORM_lb':>8} {'<A,B>=G':>8} {'d|q+1':>6} {'d|q-1':>6}")
    for q in qs:
        r = scan_q(q)
        if r is None:
            print(f"{q:>4} {q % 8:>4}  (no exterior S4)")
            continue
        a = r["agg"]
        print(f"{q:>4} {q % 8:>4} {r['n_cos']:>3} "
              f"{a['cusp'][0]:>8.4f} {a['cusp'][1]:>8.4f} {a['princ'][0]:>8.4f} {a['princ'][1]:>8.4f} "
              f"{a['stein'][0]:>8.4f} {r['norm_lb']:>8.4f} {str(r['all_genG']):>8} "
              f"{str(r['d_cusp']):>6} {str(r['d_prin']):>6}")
    print("\nNORM_lb = max over families/characters of |trace|/dim  <=  ||P_A P_B|| (true norm, "
          "from power iteration).  trace can exceed 1 (=> dim rho^{S4} >= 2); the norm cannot.")


if __name__ == "__main__":
    main()
