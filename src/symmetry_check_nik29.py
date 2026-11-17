#!/usr/bin/env python3
"""Explicit check of every WLOG step of the max-secant split in enc_nik29.py (planar q).

For each step the script exhibits the affine / semilinear map as a point permutation, verifies that it
is a collineation (maps every line onto a line, hence preserves 'Nikodym complement', |C| and the
maximum secant size), and verifies the claimed effect.  Nothing relies on a group being a full
stabiliser; only coverage is needed.
  (1) for every line L: an affine map sending the x-axis onto L           (transitivity on lines)
  (2) the maps (x,y) -> (a x^sigma + b, y^sigma) are collineations fixing the axis setwise; the orbits of
      s-subsets of the axis under them are exactly the orbits enumerated by the encoder; every
      representative contains 0 and 1; the orbit sizes sum to C(q,s)      (axis normalisation)
  (3) for every line T0 != axis through (0,0): a shear fixing the axis pointwise, T0 -> {x = 0}
  (4) for every line T1 != axis through (1,0): a map (x,y) -> (x, d y) fixing the axis pointwise and
      {x = 0} setwise, T1 -> {x = 1} or {x = 1 + y}
usage: symmetry_check_nik29.py [q ...]   (default 9; also run for 7, 8 as a cross-check)
"""
import itertools, sys
from ag import AG
from enc_nik29 import axis_orbit_reps


def perm_of(g, M, b, frob=0):
    F = g.F
    perm = []
    for p in g.points:
        x, y = p
        for _ in range(frob):
            x, y = F.frob[x], F.frob[y]
        img = g.add((g.dot(M[0], (x, y)), g.dot(M[1], (x, y))), b)
        perm.append(g.pidx[img])
    assert sorted(perm) == list(range(len(g.points))), "not a bijection"
    for _, pts in g.lines:
        assert frozenset(perm[p] for p in pts) in g.line_of, "not a collineation"
    return perm


def line_img(g, perm, l):
    return g.line_of[frozenset(perm[p] for p in g.lines[l][1])]


def main(q):
    g = AG(2, q)
    F = g.F
    axis = g.line_through((0, 0), (1, 0))
    x0 = g.line_through((0, 0), (0, 1))
    x1 = g.line_through((1, 0), (0, 1))
    diag1 = g.line_through((1, 0), (1, 1))
    # (1) transitivity on lines
    for l, (di, pts) in enumerate(g.lines):
        d = g.dirs[di]
        p = g.points[pts[0]]
        # columns of M: d and any vector completing it to a basis
        for v in g.points:
            M = [[d[0], v[0]], [d[1], v[1]]]
            try:
                perm = perm_of(g, M, p)
            except AssertionError:
                continue
            break
        assert line_img(g, perm, axis) == l
    print(f"q={q}: (1) every one of the {len(g.lines)} lines is the image of the axis under an explicit affine map")
    # (2) axis maps
    axis_pts = set(g.lines[axis][1])
    maps = []
    for a in range(1, q):
        for b in range(q):
            for e in range(F.k):
                perm = perm_of(g, [[a, 0], [0, 1]], (b, 0), e)
                assert {perm[p] for p in axis_pts} == axis_pts
                maps.append(perm)
    reps = axis_orbit_reps(F)
    for s in range(2, q + 1):
        all_subsets = set(itertools.combinations(range(q), s))
        covered = set()
        for rep, size in reps[s]:
            assert 0 in rep and 1 in rep
            orb = set()
            for perm in maps:
                img = tuple(sorted(g.points[perm[g.pidx[(t, 0)]]][0] for t in rep))
                orb.add(img)
            assert len(orb) == size
            assert not (orb & covered), "orbits overlap"
            covered |= orb
        assert covered == all_subsets, f"s={s}: orbits do not cover all subsets"
    print(f"q={q}: (2) {len(maps)} axis-preserving semilinear collineations; orbit representatives cover all s-subsets of the axis disjointly, s=2..{q}: "
          + ", ".join(f"s={s}:{len(reps[s])}" for s in range(2, q + 1)))
    # (3) shear
    for l in g.lines_through[g.pidx[(0, 0)]]:
        if l == axis:
            continue
        d = g.dirs[g.lines[l][0]]  # (b, 1) after normalisation? dirs are normalised with first nonzero = 1
        # line through origin with direction d: find b with d ~ (b, 1)
        b = F.mul[d[0]][F.inv[d[1]]]
        perm = perm_of(g, [[1, F.neg[b]], [0, 1]], (0, 0))
        assert all(perm[p] == p for p in axis_pts)
        assert line_img(g, perm, l) == x0
    print(f"q={q}: (3) every non-axis line through (0,0) is mapped to x = 0 by a shear fixing the axis pointwise")
    # (4) scaling
    for l in g.lines_through[g.pidx[(1, 0)]]:
        if l == axis:
            continue
        d = g.dirs[g.lines[l][0]]
        b = F.mul[d[0]][F.inv[d[1]]]  # line x = 1 + b y
        dd = 1 if b == 0 else b
        perm = perm_of(g, [[1, 0], [0, dd]], (0, 0))
        assert all(perm[p] == p for p in axis_pts)
        assert line_img(g, perm, x0) == x0
        assert line_img(g, perm, l) == (x1 if b == 0 else diag1)
    print(f"q={q}: (4) every non-axis line through (1,0) is mapped to x = 1 or x = 1 + y by a y-scaling fixing the axis pointwise and x = 0 setwise")


if __name__ == "__main__":
    qs = [int(t) for t in sys.argv[1:]] or [9]
    for q in qs:
        main(q)
    print("NIK29 SYMMETRY CHECKS PASSED")
