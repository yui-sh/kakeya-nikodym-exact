#!/usr/bin/env python3
"""Computational check of every symmetry-breaking (WLOG) step used by the Nikodym encodings.

For each step we (1) list the group G used, as explicit affine maps x -> M x + b, (2) verify every
element is an affine bijection mapping lines onto lines (so it preserves 'Nikodym complement'
and |C|), (3) verify the claimed transitivity / orbit structure by explicit orbit computation.
Nothing here relies on a group being the full stabiliser.
usage: symmetry_check.py [q ...]     (default: planar q = 2,3,4,5,7,8 and n=3 q=3,4)
"""
import itertools, sys
from ag import AG


def affine_map(g, M, b):
    n = g.n
    perm = [g.pidx[g.add(tuple(g.dot(M[r], p) for r in range(n)), b)] for p in g.points]
    if sorted(perm) != list(range(len(g.points))):
        return None
    for _, pts in g.lines:  # collineation check
        assert frozenset(perm[p] for p in pts) in g.line_of
    return perm


def line_image(g, perm, l):
    return g.line_of[frozenset(perm[p] for p in g.lines[l][1])]


def check_flag_transitive(g):
    """AGL(n,q) is transitive on flags: every flag (p, L through p) is the image of (O, x1-axis)."""
    n, q = g.n, g.q
    O = tuple([0] * n)
    e1 = tuple(int(i == 0) for i in range(n))
    base = g.line_through(O, e1)
    count = 0
    for p in g.points:
        for l in g.lines_through[g.pidx[p]]:
            d = g.dirs[g.lines[l][0]]
            # columns of M: d, then vectors completing d to a basis (found by search)
            for rest in itertools.product(g.points, repeat=n - 1):
                cols = [d] + list(rest)
                M = [[cols[c][r] for c in range(n)] for r in range(n)]
                perm = affine_map(g, M, p)
                if perm is not None:
                    break
            assert perm is not None
            assert perm[g.pidx[O]] == g.pidx[p] and line_image(g, perm, base) == l
            count += 1
    assert count == len(g.points) * len(g.dirs)
    return count


def check_second_point(g):
    """the stabiliser of the base flag contains maps fixing the x1-axis pointwise and sending e2 to any
    point off the x1-axis (so a second point of C may be taken to be e2)."""
    n, q = g.n, g.q
    O = tuple([0] * n)
    e1 = tuple(int(i == 0) for i in range(n))
    e2 = tuple(int(i == 1) for i in range(n))
    axis = set(g.lines[g.line_through(O, e1)][1])
    count = 0
    for r in g.points:
        if g.pidx[r] in axis:
            continue
        found = False
        for rest in itertools.product(g.points, repeat=n - 2):
            cols = [e1, r] + list(rest)   # M e1 = e1, M e2 = r
            M = [[cols[c][rr] for c in range(n)] for rr in range(n)]
            perm = affine_map(g, M, O)
            if perm is not None:
                found = True
                break
        assert found and all(perm[p] == p for p in axis) and perm[g.pidx[e2]] == g.pidx[r]
        count += 1
    assert count == len(g.points) - q
    return count


def residual_orbits(g):
    """orbits of the residual group {M : M e1 = a e1, M e2 = e2, M e_j = any (j>=3), M invertible}
    (for n = 2 this is just (x,y) -> (ax, y))
    (it fixes O, the x1-axis setwise and e2) on the directions of possible tangents at e2
    (all directions except e2's own, whose line contains O in C)."""
    n, q = g.n, g.q
    O = tuple([0] * n)
    e1 = tuple(int(i == 0) for i in range(n))
    e2 = tuple(int(i == 1) for i in range(n))
    group = []
    for a in range(1, q):
        for rest in itertools.product(g.points, repeat=n - 2):
            cols = [g.scale(a, e1), e2] + list(rest)
            M = [[cols[c][r] for c in range(n)] for r in range(n)]
            perm = affine_map(g, M, O)
            if perm is not None:
                group.append(perm)
    for perm in group:
        assert perm[g.pidx[O]] == g.pidx[O] and perm[g.pidx[e2]] == g.pidx[e2]
        assert line_image(g, perm, g.line_through(O, e1)) == g.line_through(O, e1)
    dirs = [d for d in range(len(g.dirs)) if g.dirs[d] != e2]
    orbits = []
    seen = set()
    for d in dirs:
        if d in seen:
            continue
        orb = {d}
        stack = [d]
        while stack:
            x = stack.pop()
            for perm in group:
                # image direction of the tangent line at e2 with direction x
                l = g.line_through(e2, g.dirs[x])
                y = g.lines[line_image(g, perm, l)][0]
                if y not in orb:
                    orb.add(y)
                    stack.append(y)
        seen |= orb
        orbits.append(sorted(orb))
    return [[g.dirs[d] for d in orb] for orb in orbits]


if __name__ == "__main__":
    specs = [(2, q) for q in (2, 3, 4, 5, 7, 8)] + [(3, 3), (3, 4)]
    if len(sys.argv) > 1:
        specs = [(2, int(q)) for q in sys.argv[1:]]
    for n, q in specs:
        g = AG(n, q)
        nf = check_flag_transitive(g)
        ns = check_second_point(g)
        orbs = residual_orbits(g)
        reps = [o[0] for o in orbs]
        e1 = tuple(int(i == 0) for i in range(n))
        if n == 2:
            assert len(orbs) == 2 and reps == [e1, (1, 1)], orbs
        else:
            assert len(orbs) == 3 and set(reps) == {e1, (1, 1, 0), (0, 0, 1)}, orbs
        print(f"AG({n},{q}): flag-transitivity checked on {nf} flags; second point: {ns} images of e2; "
              f"tangent-direction orbits at e2: {[(o[0], len(o)) for o in orbs]}  OK")
    print("SYMMETRY CHECKS PASSED")
