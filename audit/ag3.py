"""AUDIT helper (fresh code, no imports from lp/src): AG(3,q) for q in {3,4,5}, witness systems, plane statistics."""
import random, itertools, math
from collections import Counter, defaultdict

class Field:
    def __init__(self, q):
        self.q = q
        if q in (2, 3, 5):
            self.add = lambda a, b: (a + b) % q; self.mul = lambda a, b: (a * b) % q; self.neg = lambda a: (-a) % q
        elif q == 4:   # GF(4) = {0,1,a,a+1} as 0,1,2,3; a^2 = a+1
            M = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
            self.add = lambda a, b: a ^ b; self.mul = lambda a, b: M[a][b]; self.neg = lambda a: a
        else: raise ValueError(q)
        self.inv = {a: next(b for b in range(1, q) if self.mul(a, b) == 1) for a in range(1, q)}
        # axioms spot-check
        els = range(q)
        assert all(self.mul(a, self.add(b, c)) == self.add(self.mul(a, b), self.mul(a, c)) for a in els for b in els for c in els)
        assert all(self.mul(a, self.inv[a]) == 1 for a in range(1, q))

class AG3:
    def __init__(self, q):
        F = Field(q); self.F = F; self.q = q
        self.points = list(itertools.product(range(q), repeat=3)); self.pidx = {p: i for i, p in enumerate(self.points)}
        def norm(v):   # scale so the first nonzero coordinate is 1
            i = next(i for i in range(3) if v[i]); s = F.inv[v[i]]; return tuple(F.mul(s, x) for x in v)
        self.dirs = sorted({norm(v) for v in self.points if any(v)}); self.didx = {d: i for i, d in enumerate(self.dirs)}
        assert len(self.dirs) == q * q + q + 1
        add = lambda p, v: tuple(F.add(a, b) for a, b in zip(p, v)); smul = lambda s, v: tuple(F.mul(s, x) for x in v)
        self.add, self.smul = add, smul
        self.lines = []; self.line_of = {}; self.lines_by_dir = defaultdict(list)
        for di, d in enumerate(self.dirs):
            for p in self.points:
                pts = frozenset(self.pidx[add(p, smul(t, d))] for t in range(q))
                if pts not in self.line_of:
                    self.line_of[pts] = len(self.lines); self.lines.append((di, pts)); self.lines_by_dir[di].append(len(self.lines) - 1)
        assert len(self.lines) == q * q * len(self.dirs)
        self.lines_through = defaultdict(list)
        for li, (di, pts) in enumerate(self.lines):
            for p in pts: self.lines_through[p].append(li)
        # planes: normal a (normalized) and constant c; points x with a.x = c
        dot = lambda a, x: F.add(F.add(F.mul(a[0], x[0]), F.mul(a[1], x[1])), F.mul(a[2], x[2]))
        self.dot = dot
        self.planes = []   # (class index, c, frozenset of point indices, set of direction indices at infinity, list of line indices inside)
        self.planes_of_line = defaultdict(list); self.planes_through = defaultdict(list)
        for ci, a in enumerate(self.dirs):
            atinf = frozenset(self.didx[d] for d in self.dirs if dot(a, d) == 0); assert len(atinf) == q + 1
            for c in range(q):
                pts = frozenset(self.pidx[x] for x in self.points if dot(a, x) == c); assert len(pts) == q * q
                inside = [li for li, (di, lp) in enumerate(self.lines) if di in atinf and lp <= pts]; assert len(inside) == q * (q + 1)
                pi = len(self.planes); self.planes.append((ci, c, pts, atinf, inside))
                for li in inside: self.planes_of_line[li].append(pi)
                for p in pts: self.planes_through[p].append(pi)
        assert all(len(v) == q + 1 for v in self.planes_of_line.values()) and all(len(v) == q * q + q + 1 for v in self.planes_through.values())
        # collinear triples of directions in PG(2,q): d1,d2,d3 collinear iff det = 0
        self.det0 = {}
    def collinear(self, d1, d2, d3):
        F = self.F; a, b, c = self.dirs[d1], self.dirs[d2], self.dirs[d3]
        m = lambda x, y: F.mul(x, y); s = F.add; n = F.neg
        det = s(s(m(a[0], s(m(b[1], c[2]), n(m(b[2], c[1])))), n(m(a[1], s(m(b[0], c[2]), n(m(b[2], c[0])))))), m(a[2], s(m(b[0], c[1]), n(m(b[1], c[0])))))
        return det == 0

def random_witness(ag, rng, mode='random', forced_planes=0):
    """one line per direction. mode: random | greedy (choose line maximising overlap with current union) | forced (first
    forced_planes planes get 5 (q=5) inside lines in general position, rest random)"""
    q = ag.q; W = {}
    if mode == 'forced':
        used = set()
        for _ in range(forced_planes):
            pi = rng.randrange(len(ag.planes)); ci, c, pts, atinf, inside = ag.planes[pi]
            ds = [d for d in atinf if d not in used]
            if len(ds) < 5: continue
            ds = rng.sample(ds, 5)
            for _try in range(200):
                choice = {d: rng.choice([li for li in inside if ag.lines[li][0] == d]) for d in ds}
                U = set().union(*(ag.lines[li][1] for li in choice.values()))
                if len(U) == 5 * q - 10: break
            for d, li in choice.items(): W[d] = li; used.add(d)
    K = set()
    for li in W.values(): K |= ag.lines[li][1]
    order = list(range(len(ag.dirs))); rng.shuffle(order)
    for d in order:
        if d in W: continue
        if mode == 'greedy':
            best = max(rng.sample(ag.lines_by_dir[d], len(ag.lines_by_dir[d])), key=lambda li: len(ag.lines[li][1] & K))
            W[d] = best
        else: W[d] = rng.choice(ag.lines_by_dir[d])
        K |= ag.lines[W[d]][1]
    return W, K

def stats(ag, W, K):
    """per-point degrees, per-plane (r, k, h, U, b, t, e, inside witness lines) etc."""
    q = ag.q
    deg = Counter(p for li in W.values() for p in ag.lines[li][1]); D = max(deg.values())
    Wset = set(W.values())
    P = []
    for pi, (ci, c, pts, atinf, inside) in enumerate(ag.planes):
        Kp = pts & K; ins = [li for li in inside if li in Wset]
        U = set().union(*(ag.lines[li][1] for li in ins)) if ins else set()
        b = {p: sum(p in ag.lines[li][1] for li in ins) for p in Kp}
        t = [0] * (q + 1)
        for li in inside: t[len(ag.lines[li][1] & Kp)] += 1
        marked = {ag.lines[li][0] for li in ins}
        e = {d: tuple(sorted(len(ag.lines[li][1] & Kp) for li in inside if ag.lines[li][0] == d)) for d in atinf if d not in marked}
        h = tuple(sum(deg[p] == j for p in Kp) for j in range(1, D + 1))
        P.append(dict(r=len(ins), k=len(Kp), h=h, U=U, u=len(U), b=b, t=tuple(t[1:]), e=e, ins=ins, Kp=Kp, marked=marked, c=tuple(sum(v == i for v in b.values()) for i in range(1, len(ins) + 1))))
    return deg, D, P
