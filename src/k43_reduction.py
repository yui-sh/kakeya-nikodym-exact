#!/usr/bin/env python3
"""Exact (integer/Fraction) replay of the hand reduction for K(4,3) >= 27.

Claim: a Kakeya set S in AG(4,3) with |S| = 26 has a hyperplane direction whose three parallel
slices have sizes (a,b,c) in the list TYPES below (9 types).  Together with the SAT refutation of
every (type, orbit representative of the smallest slice) case this gives K(4,3) >= 27; since
adding points keeps a set Kakeya, excluding size exactly 26 excludes every smaller size.

Ingredients (all re-derived here; the incidence constants are recounted from AG(4,3) itself):
 (C)  capacity table f(k) = max number of directions of F_3^3 having a full line inside some
      k-point set (certified separately by the SAT instances cap_k_r: 'k points, >= f(k)+1
      directions' is UNSAT; k >= 13 uses the trivial bound 13).
 (S)  slice constraint: the 13 directions parallel to the hyperplane must each have a witness
      line inside one of the 3 slices, so f(a)+f(b)+f(c) >= 13.
 (P)  the 27 transverse directions: a transverse witness line meets every slice once and is
      determined by its points in two slices, so ab, ac, bc >= 27.
 (E)  additive energy: for slices A,B (|A|=a,|B|=b, t=min(a,b), m=ab) the difference
      multiplicities r(d) = #{(x,y) in A x B : y-x = d} satisfy 1 <= r(d) <= t (every transverse
      direction occurs) and sum r(d) = m, so by convexity E(A,B) = sum r(d)^2 <= 27 + qq(t^2-1) + rr(rr+2)
      where m-27 = qq(t-1)+rr, 0 <= rr < t-1.  E(A,B) also equals sum_s #{(x,y): x+y=s}^2; in
      characteristic 3 the three points of a line sum to 0, so the 27 transverse witness pairs
      (x,y) (distinct, one per direction) have x+y in -C; Cauchy-Schwarz gives E >= 27^2/c.
 (M)  moment identities over the 40 hyperplane directions h, with Q_h = a^2+b^2+c^2 and
      T_h = C(a,3)+C(b,3)+C(c,3):  sum_h Q_h = 13 N^2 + 27 N  and  sum_h T_h = 4 C(N,3) + 9 L,
      where L >= 40 is the number of full lines in S (a point lies in 40 hyperplanes, a pair in 13,
      a line in 13, a non-collinear triple in 4; the script recounts these).
 (R)  for every admissible type with c <= 11:  4 T <= 13 Q - 2018  (checked by enumeration).
      Summing over h: if all 40 types had c <= 11 then 4 sum T <= 13 sum Q - 80720 = 42650,
      but 4 sum T >= 16 C(26,3) + 1440 = 43040.  So some hyperplane direction has c >= 12, and
      its type is an admissible type with c >= 12: exactly the 9 types listed.
usage: k43_reduction.py [--json out.json]
"""
import itertools, json, sys
from fractions import Fraction as Fr
from math import comb

F_TABLE = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 2, 6: 3, 7: 4, 8: 4, 9: 6, 10: 7, 11: 9, 12: 10}
N = 26


def f(k):
    return F_TABLE[k] if k <= 12 else 13


def energy_ok(a, b, c):
    """True if the energy test does NOT exclude ordering (A,B,C)"""
    t, m = min(a, b), a * b
    if t <= 1:
        upper = 27  # all multiplicities equal 1 (m must be 27)
        if m != 27:
            return False
    else:
        qq, rr = divmod(m - 27, t - 1)
        assert 0 <= qq <= 27
        upper = 27 + qq * (t * t - 1) + rr * (rr + 2)
    return Fr(upper) >= Fr(729, c)


def admissible(a, b, c):
    if f(a) + f(b) + f(c) < 13:
        return False, "slice-capacity"
    if min(a * b, a * c, b * c) < 27:
        return False, "pair-product"
    for A, B, C in itertools.permutations((a, b, c)):
        if not energy_ok(A, B, C):
            return False, "energy"
    return True, ""


def recount_incidences():
    """count hyperplanes of AG(4,3) through a point / a pair / a line / a non-collinear triple"""
    pts = list(itertools.product(range(3), repeat=4))
    normals = [v for v in pts if any(v) and v[next(i for i, c in enumerate(v) if c)] == 1]
    assert len(normals) == 40
    hyper = [frozenset(p for p in pts if sum(x * y for x, y in zip(nv, p)) % 3 == b) for nv in normals for b in range(3)]
    assert len(hyper) == 120
    p0, p1, p2 = (0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0)  # collinear
    p3 = (0, 1, 0, 0)
    cnt = lambda S: sum(1 for H in hyper if S <= H)
    return dict(point=cnt({p0}), pair=cnt({p0, p1}), line=cnt({p0, p1, p2}), triangle=cnt({p0, p1, p3}))


def main():
    inc = recount_incidences()
    assert inc == dict(point=40, pair=13, line=13, triangle=4), inc
    sumQ = inc["pair"] * N * (N - 1) + inc["point"] * N            # 13 N^2 + 27 N
    assert sumQ == 13 * N * N + 27 * N
    Lmin = 40
    sumT_min = inc["triangle"] * comb(N, 3) + (inc["line"] - inc["triangle"]) * Lmin  # 4 C(N,3) + 9 L
    types = [(a, b, c) for a in range(N + 1) for b in range(a, N + 1) for c in range(b, N + 1) if a + b + c == N]
    adm, excluded = [], {}
    for t in types:
        ok, why = admissible(*t)
        if ok:
            adm.append(t)
        else:
            excluded.setdefault(why, []).append(t)
    # (R): the separating inequality on every admissible type with c <= 11
    K = 2018
    for a, b, c in adm:
        if c <= 11:
            Q = a * a + b * b + c * c
            T = comb(a, 3) + comb(b, 3) + comb(c, 3)
            assert 4 * T <= 13 * Q - K, (a, b, c, 4 * T, 13 * Q - K)
    upper = 13 * sumQ - 40 * K      # bound on 4 sum T if every direction had c <= 11
    lower = 4 * sumT_min
    assert upper == 42650 and lower == 43040 and upper < lower
    final = [t for t in adm if t[2] >= 12]
    expected = [(4, 9, 13), (4, 10, 12), (5, 7, 14), (5, 8, 13), (5, 9, 12), (6, 6, 14), (6, 7, 13), (6, 8, 12), (7, 7, 12)]
    print("incidence constants recounted:", inc)
    print(f"admissible slice types (f/pair/energy): {len(adm)}; of these with c<=11: {sum(1 for t in adm if t[2] <= 11)} "
          f"(all satisfy 4T <= 13Q - {K})")
    print(f"energy-excluded types with c>=12: {[t for t in excluded.get('energy', []) if t[2] >= 12]}")
    print(f"if all 40 directions had c<=11: 4*sum T <= {upper} < {lower} <= 4*sum T  -> contradiction")
    print(f"remaining types (c>=12): {final}")
    assert final == expected, final
    print("K43 REDUCTION VERIFIED: 9 slice types")
    if "--json" in sys.argv:
        json.dump({"types": final, "f_table": F_TABLE}, open(sys.argv[sys.argv.index("--json") + 1], "w"))


if __name__ == "__main__":
    main()
