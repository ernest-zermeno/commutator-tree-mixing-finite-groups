#!/usr/bin/env python3
"""Exact audit of the Hamilton S4 intersection after reduction modulo q.

In the trace-zero model, A is the 24 determinant-one signed permutation
matrices and Ad(h_0) = M/3.  A new element in A cap hAh^{-1} modulo a prime q
would give M R M = 9 R' modulo q for a pair R,R' in A that is not an equality
over the integers.  The gcd census below proves that this can happen only in
characteristic 2 or 3, both excluded in the paper.
"""

import itertools
import math
from functools import reduce


def matmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def determinant(a):
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def signed_permutation_group():
    group = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            matrix = [[0] * 3 for _ in range(3)]
            for row, column in enumerate(permutation):
                matrix[row][column] = signs[row]
            matrix = tuple(tuple(row) for row in matrix)
            if determinant(matrix) == 1:
                group.append(matrix)
    return group


def main():
    group = signed_permutation_group()
    m = ((-1, 2, 2), (2, -1, 2), (2, 2, -1))
    exact = 0
    gcd_census = {}

    for r in group:
        mrm = matmul(matmul(m, r), m)
        for r_prime in group:
            differences = [
                mrm[i][j] - 9 * r_prime[i][j]
                for i in range(3)
                for j in range(3)
            ]
            common_gcd = reduce(math.gcd, map(abs, differences), 0)
            if common_gcd == 0:
                exact += 1
            else:
                gcd_census[common_gcd] = gcd_census.get(common_gcd, 0) + 1

    assert len(group) == 24
    assert exact == 6
    assert gcd_census == {1: 360, 2: 54, 4: 18, 9: 120, 18: 18}
    for common_gcd in gcd_census:
        residue = common_gcd
        for prime in (2, 3):
            while residue % prime == 0:
                residue //= prime
        assert residue == 1

    print("|A| = 24")
    print("exact intersection size = 6")
    print("nonzero gcd census =", sorted(gcd_census.items()))
    print("No reduction modulo a prime q != 2,3 enlarges the intersection: PASS")


if __name__ == "__main__":
    main()
