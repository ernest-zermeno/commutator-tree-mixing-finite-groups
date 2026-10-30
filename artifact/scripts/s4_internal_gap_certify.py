"""
Exact certificate for the S_4 internal gap.

The exploratory script s4_internal_gap.py computes

    max ||Dhat_2^pi|| = 0.084351...

with floating-point SVDs.  This script gives a rational certificate for the
weaker bound

    ||Dhat_2^pi|| <= 17/200 = 0.085

for every non-trivial irreducible pi of S_4 and every valid level-2
configuration of transposition leaves generating S_4.

Representations are written in rational bases.  Their invariant inner products
have Gram matrices G, so the operator-norm certificate is the exact Loewner
inequality

    D^T G D <= (17/200)^2 G.

Since dimensions are at most 3, positive semidefiniteness is checked by all
principal minors, exactly over fractions.
"""

from fractions import Fraction
from itertools import combinations, permutations, product


PTS = (0, 1, 2, 3)
E = (0, 1, 2, 3)
BOUND = Fraction(17, 200)
BOUND2 = BOUND * BOUND


def mult(a, b):
    return tuple(a[b[i]] for i in range(4))


def inv(a):
    r = [0] * 4
    for i in range(4):
        r[a[i]] = i
    return tuple(r)


def comm(a, b):
    return mult(mult(a, b), mult(inv(a), inv(b)))


def sgn(a):
    s = 1
    for i in range(4):
        for j in range(i + 1, 4):
            if a[i] > a[j]:
                s = -s
    return s


S4 = [tuple(p) for p in permutations(PTS)]
TRANSPOSITIONS = [
    tuple({i: j, j: i}.get(k, k) for k in PTS)
    for i in range(4)
    for j in range(i + 1, 4)
]


def closure(gens):
    S = {E}
    fr = [E]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mult(g, x)
            if y not in S:
                S.add(y)
                fr.append(y)
    return S


def mat_zero(n, m):
    return [[Fraction(0) for _ in range(m)] for _ in range(n)]


def mat_eye(n):
    M = mat_zero(n, n)
    for i in range(n):
        M[i][i] = Fraction(1)
    return M


def mat_add(A, B):
    return [
        [A[i][j] + B[i][j] for j in range(len(A[0]))]
        for i in range(len(A))
    ]


def mat_sub(A, B):
    return [
        [A[i][j] - B[i][j] for j in range(len(A[0]))]
        for i in range(len(A))
    ]


def mat_scalar(c, A):
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def mat_mul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    C = mat_zero(n, m)
    for i in range(n):
        for t in range(k):
            if A[i][t]:
                for j in range(m):
                    C[i][j] += A[i][t] * B[t][j]
    return C


def mat_transpose(A):
    return [list(row) for row in zip(*A)]


def det(A):
    n = len(A)
    if n == 0:
        return Fraction(1)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    total = Fraction(0)
    for j in range(n):
        minor = [
            [A[i][k] for k in range(n) if k != j]
            for i in range(1, n)
        ]
        total += ((-1) ** j) * A[0][j] * det(minor)
    return total


def principal_submatrix(A, idx):
    return [[A[i][j] for j in idx] for i in idx]


def is_psd_exact(A):
    n = len(A)
    min_minor = None
    for r in range(1, n + 1):
        for idx in combinations(range(n), r):
            d = det(principal_submatrix(A, idx))
            min_minor = d if min_minor is None else min(min_minor, d)
            if d < 0:
                return False, min_minor, idx, d
    return True, min_minor, None, None


def std_matrix_from_action(n, action):
    """Standard sum-zero representation for a permutation action on n points.

    Basis: b_i=e_i-e_{n-1}, i=0,...,n-2.  The invariant Gram matrix is I+J.
    """
    m = n - 1
    M = mat_zero(m, m)
    last = action(n - 1)
    for col in range(m):
        p, q = action(col), last
        if p < m:
            M[p][col] += 1
        if q < m:
            M[q][col] -= 1
    return M


PAIRINGS = [
    frozenset({frozenset({0, 1}), frozenset({2, 3})}),
    frozenset({frozenset({0, 2}), frozenset({1, 3})}),
    frozenset({frozenset({0, 3}), frozenset({1, 2})}),
]


def pairing_action(a, i):
    image = frozenset(frozenset({a[x] for x in pair}) for pair in PAIRINGS[i])
    return PAIRINGS.index(image)


def gram_sum_zero(n):
    m = n - 1
    G = mat_eye(m)
    for i in range(m):
        for j in range(m):
            G[i][j] += 1
    return G


def gram_sum_zero_inv(n):
    m = n - 1
    Ginv = mat_eye(m)
    for i in range(m):
        for j in range(m):
            Ginv[i][j] -= Fraction(1, n)
    return Ginv


def rep_sign(a):
    return [[Fraction(sgn(a))]]


def rep_std(a):
    return std_matrix_from_action(4, lambda i: a[i])


def rep_2dim(a):
    return std_matrix_from_action(3, lambda i: pairing_action(a, i))


def rep_std_twist(a):
    return mat_scalar(Fraction(sgn(a)), rep_std(a))


IRREPS = {
    "sign(1)": (rep_sign, [[Fraction(1)]], [[Fraction(1)]]),
    "std(3)": (rep_std, gram_sum_zero(4), gram_sum_zero_inv(4)),
    "2dim(2)": (rep_2dim, gram_sum_zero(3), gram_sum_zero_inv(3)),
    "std'(3)": (rep_std_twist, gram_sum_zero(4), gram_sum_zero_inv(4)),
}


def adjoint(A, G, Ginv):
    return mat_mul(mat_mul(Ginv, mat_transpose(A)), G)


def projection(sig, rep):
    M = rep(sig)
    return mat_scalar(Fraction(1, 2), mat_add(mat_eye(len(M)), M))


def block(s, t, rep):
    M = mat_mul(projection(s, rep), projection(t, rep))
    return mat_mul(M, M)


def level2(quad, rep, G, Ginv):
    s1, t1, s2, t2 = quad
    A = block(s1, t1, rep)
    B = block(s2, t2, rep)
    return mat_mul(mat_mul(mat_mul(A, B), adjoint(A, G, Ginv)), adjoint(B, G, Ginv))


def valid_configs():
    for quad in product(TRANSPOSITIONS, repeat=4):
        if len(closure(quad)) != 24:
            continue
        s1, t1, s2, t2 = quad
        if comm(comm(s1, t1), comm(s2, t2)) == E:
            continue
        yield quad


def certify_operator_bound(D, G):
    lhs = mat_mul(mat_mul(mat_transpose(D), G), D)
    rhs = mat_scalar(BOUND2, G)
    return is_psd_exact(mat_sub(rhs, lhs))


def main():
    configs = list(valid_configs())
    print(f"valid level-2 S_4 configs: {len(configs)}")
    print(f"certifying ||Dhat_2^pi|| <= {BOUND} = {float(BOUND):.6f}")

    global_min_minor = None
    for name, (rep, G, Ginv) in IRREPS.items():
        min_minor = None
        worst_cfg = None
        for quad in configs:
            D = level2(quad, rep, G, Ginv)
            ok, mm, idx, bad = certify_operator_bound(D, G)
            min_minor = mm if min_minor is None else min(min_minor, mm)
            if not ok:
                print(f"FAILED for {name}, principal minor {idx}: {bad}")
                print(f"config={quad}")
                return 1
            worst_cfg = quad
        global_min_minor = min_minor if global_min_minor is None else min(global_min_minor, min_minor)
        print(f"  {name:9s}: certified; smallest principal minor margin = {min_minor}")
    print()
    print(f"ALL CERTIFIED: lambda_2^(S_4) <= {BOUND} = {float(BOUND):.6f}")
    print(f"UNIFORM internal gap: g_0 >= 1 - {BOUND} = {1 - float(BOUND):.6f}")
    print(f"global smallest principal-minor margin: {global_min_minor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
