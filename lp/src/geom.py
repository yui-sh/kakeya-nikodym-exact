"""Affine geometry AG(n,q): points, directions (projective points), lines."""
import itertools
from gf import GF

class AG:
    def __init__(self, n, q):
        self.n, self.q = n, q
        self.F = F = GF(q)
        self.points = list(itertools.product(range(q), repeat=n))
        self.pidx = {p: i for i, p in enumerate(self.points)}
        # directions: nonzero vectors normalized so first nonzero coordinate is 1
        dirs = []
        for v in self.points:
            if any(v):
                i = next(i for i, c in enumerate(v) if c)
                if v[i] == 1:
                    dirs.append(v)
        self.dirs = dirs
        self.didx = {d: i for i, d in enumerate(dirs)}
        # lines per direction: cosets of span(d); line as sorted tuple of point indices
        self.lines = []       # list of (dir index, frozenset of point idx)
        self.lines_by_dir = [[] for _ in dirs]
        seen = set()
        for di, d in enumerate(dirs):
            for p in self.points:
                pts = frozenset(self.pidx[self.add(p, self.scale(t, d))] for t in range(q))
                if pts in seen:
                    continue
                seen.add(pts)
                self.lines_by_dir[di].append(len(self.lines))
                self.lines.append((di, pts))
        assert all(len(L) == q ** (n - 1) for L in self.lines_by_dir)
        self.line_of = {pts: i for i, (di, pts) in enumerate(self.lines)}

    def add(self, u, v):
        return tuple(self.F.add[a][b] for a, b in zip(u, v))
    def sub(self, u, v):
        return tuple(self.F.add[a][self.F.neg[b]] for a, b in zip(u, v))
    def scale(self, t, v):
        return tuple(self.F.mul[t][a] for a in v)
    def normalize_dir(self, v):
        i = next(i for i, c in enumerate(v) if c)
        return self.scale(self.F.inv[v[i]], v)
    def line_through(self, p, d):
        """index of the line through point p (tuple) with direction d (tuple)."""
        pts = frozenset(self.pidx[self.add(p, self.scale(t, d))] for t in range(self.q))
        return self.line_of[pts]
    def unit(self, i):
        v = [0] * self.n; v[i] = 1; return tuple(v)

def is_kakeya(ag, S):
    """S: iterable of point indices. Returns (ok, missing directions)."""
    S = set(S)
    missing = []
    for di in range(len(ag.dirs)):
        if not any(pts <= S for li in ag.lines_by_dir[di] for (_, pts) in [ag.lines[li]]):
            missing.append(ag.dirs[di])
    return (not missing), missing

if __name__ == "__main__":
    for n, q in [(2, 3), (2, 4), (3, 2), (3, 3), (3, 4)]:
        ag = AG(n, q)
        print(n, q, len(ag.points), "pts", len(ag.dirs), "dirs", len(ag.lines), "lines")
