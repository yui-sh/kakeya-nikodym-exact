#!/usr/bin/env python3
"""Exact replay of the hand certificate K(3,5) >= 46 (kakeya-fq WRITEUP4 section 2), plain Python.

Setting: K = union of one witness line per direction in AG(3,5) (31 directions, 155 planes),
n = |K|; for a plane pi let r = number of witness lines inside pi and k = |K ∩ pi|.
Identities (counted here from the geometry where they are incidence constants):
  (I2) sum_pi 1 = 155, sum_pi r = 31*6 = 186 (each line lies in 6 planes), sum_pi k = 31 n;
  (I4) sum_pi k^2 = 6 n^2 + 25 n   (a pair of distinct points lies in 6 planes, a point in 31);
  (I5) sum_pi k r = 31 n + 775     (for p in K: sum over planes through p of r = 31 + 5 d_p,
                                    and sum_p d_p = 155).
  (P)  k >= f(r) with f = 0,5,9,12,14,15,17 for r = 0..6: r lines with distinct directions in a plane
       pairwise meet, so their union has >= 5r - C(r,2) points for r <= 5; for r = 6 the union is a
       planar Kakeya set of AG(2,5), of size >= 17 (Blokhuis-Mazzocca: K(2,q) = q(q+1)/2 + (q-1)/2
       for odd q; this is the only external input besides Bukh-Chao's n >= 39 for n <= 37).
Certificate: c(r,k) = 10k^2 - (150+17r)k + 560 + 143r >= 0 on all admissible (r,k) with
  sum_pi c = 60 n^2 - 4927 n + 100223 < 0 for 38 <= n <= 44, and = 8 for n = 45 while the values of c
  on admissible pairs are {0, 6, 7, 10, ...} (8 is not a sum of them).  Hence n >= 46.
"""
import itertools
from fractions import Fraction as Fr
from ff import GF

f = {0: 0, 1: 5, 2: 9, 3: 12, 4: 14, 5: 15, 6: 17}


def incidence_constants():
    F = GF(5)
    pts = list(itertools.product(range(5), repeat=3))
    dirs = [v for v in pts if any(v) and v[next(i for i, c in enumerate(v) if c)] == 1]
    dot = lambda u, v: sum(F.mul[a][b] for a, b in zip(u, v)) % 5
    planes = [frozenset(p for p in pts if dot(nv, p) == b) for nv in dirs for b in range(5)]
    line = frozenset(tuple(F.mul[t][1] * 0 + (t if i == 0 else 0) for i in range(3)) for t in range(5))  # x-axis
    return dict(n_planes=len(planes), planes_per_point=sum(1 for P in planes if (0, 0, 0) in P),
                planes_per_pair=sum(1 for P in planes if {(0, 0, 0), (1, 0, 0)} <= P),
                planes_per_line=sum(1 for P in planes if line <= P), n_dirs=len(dirs))


def c(r, k):
    return 10 * k * k - (150 + 17 * r) * k + 560 + 143 * r


def main():
    inc = incidence_constants()
    assert inc == dict(n_planes=155, planes_per_point=31, planes_per_pair=6, planes_per_line=6, n_dirs=31), inc
    # (P) for r <= 5 from pairwise intersections: 5r - C(r,2)
    for r in range(6):
        assert f[r] == 5 * r - r * (r - 1) // 2 or r == 0
    vals = set()
    for r in range(7):
        for k in range(max(f[r], 1), 26):
            assert c(r, k) >= 0, (r, k)
            vals.add(c(r, k))
        if r <= 2:
            assert all(c(r, k) >= 0 for k in range(0, 200))
    # sum over planes, symbolic in n: 10*(6n^2+25n) - 150*31n - 17*(31n+775) + 560*155 + 143*186
    A, B, C0 = 10 * 6, 10 * 25 - 150 * 31 - 17 * 31, -17 * 775 + 560 * 155 + 143 * 186
    assert (A, B, C0) == (60, -4927, 100223)
    total = lambda n: A * n * n + B * n + C0
    assert all(total(n) < 0 for n in range(38, 45)) and total(45) == 8
    assert 8 not in vals and all(v not in vals for v in range(1, 6))
    print("incidence constants recounted:", inc)
    print("c(r,k) >= 0 on admissible pairs; value set starts", sorted(vals)[:6])
    print("sum_pi c(r,k) = 60 n^2 - 4927 n + 100223:", {n: total(n) for n in range(38, 47)})
    print("K35 LOWER BOUND 46 CERTIFICATE VERIFIED (external inputs: Blokhuis-Mazzocca K(2,5)=17; Bukh-Chao n>=39 for n<=37)")


if __name__ == "__main__":
    main()
