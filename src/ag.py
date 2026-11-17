"""Affine geometry AG(n,q) = F_q^n: points, directions, lines, (for n=3) planes.

Everything is derived from field arithmetic in ff.py.  A direction is a nonzero vector
normalised so that its first nonzero coordinate is 1 (one representative per point of
PG(n-1,q)).  The line through p with direction d is {p + t d : t in F_q}.
"""
import itertools
from ff import GF


class AG:
    def __init__(self, n, q):
        self.n, self.q = n, q
        self.F = F = GF(q)
        self.points = list(itertools.product(range(q), repeat=n))
        self.pidx = {p: i for i, p in enumerate(self.points)}
        self.dirs = [v for v in self.points if any(v) and v[next(i for i, c in enumerate(v) if c)] == 1]
        assert len(self.dirs) == (q ** n - 1) // (q - 1)
        self.didx = {d: i for i, d in enumerate(self.dirs)}
        # lines: one per (direction, coset); stored as (dir index, sorted tuple of point indices)
        self.lines = []
        self.lines_by_dir = [[] for _ in self.dirs]
        self.line_of = {}  # frozenset(point idxs) -> line idx
        self.lines_through = [[] for _ in self.points]
        for di, d in enumerate(self.dirs):
            seen = set()
            for p in self.points:
                pi = self.pidx[p]
                if pi in seen:
                    continue
                pts = tuple(sorted(self.pidx[self.add(p, self.scale(t, d))] for t in range(q)))
                seen.update(pts)
                li = len(self.lines)
                self.lines.append((di, pts))
                self.lines_by_dir[di].append(li)
                self.line_of[frozenset(pts)] = li
                for x in pts:
                    self.lines_through[x].append(li)
        assert len(self.lines) == len(self.dirs) * q ** (n - 1)
        if n == 3:
            # planes: {x : <normal, x> = b}; one per (normal direction, b)
            self.planes = []
            for ni, nv in enumerate(self.dirs):
                for b in range(q):
                    pts = frozenset(i for i, p in enumerate(self.points) if self.dot(nv, p) == b)
                    assert len(pts) == q * q
                    lines_in = [li for li, (_, lp) in enumerate(self.lines) if set(lp) <= pts]
                    assert len(lines_in) == q * (q + 1)
                    self.planes.append((ni, b, pts, tuple(lines_in)))
            assert len(self.planes) == len(self.dirs) * q

    # vector helpers (tuples of field elements)
    def add(self, u, v):
        return tuple(self.F.add[a][b] for a, b in zip(u, v))

    def sub(self, u, v):
        return tuple(self.F.sub[a][b] for a, b in zip(u, v))

    def scale(self, t, v):
        return tuple(self.F.mul[t][a] for a in v)

    def dot(self, u, v):
        s = 0
        for a, b in zip(u, v):
            s = self.F.add[s][self.F.mul[a][b]]
        return s

    def normalize(self, v):
        i = next(i for i, c in enumerate(v) if c)
        return self.scale(self.F.inv[v[i]], v)

    def direction(self, p, r):
        """direction index of the line through two distinct points"""
        return self.didx[self.normalize(self.sub(r, p))]

    def line_through(self, p, d):
        """line index of the line through point p (tuple) with direction d (tuple)"""
        pts = frozenset(self.pidx[self.add(p, self.scale(t, d))] for t in range(self.q))
        return self.line_of[pts]

    def line_points(self, li):
        return self.lines[li][1]


if __name__ == "__main__":
    for n, q in ((2, 7), (3, 3), (3, 4), (4, 3)):
        g = AG(n, q)
        print(f"AG({n},{q}): {len(g.points)} points, {len(g.dirs)} directions, {len(g.lines)} lines"
              + (f", {len(g.planes)} planes" if n == 3 else ""))
