"""Structure of a Kakeya set S in AG(3,q): contained lines, point degrees, plane sections,
and its stabilizer in AGL(3,q) (exact, by enumerating images of an affine frame in S with pruning).
Usage: python structure.py file.json [--nostab]"""
import json, sys, itertools, time
from geom import AG

def affine_map(ag, src, dst):
    """Affine map sending frame src=(p0,p1,p2,p3) to dst (tuples). Returns perm on point idx or None if dst dependent."""
    F, q = ag.F, ag.q
    # M sends e_i -> dst_i - dst_0 ; need M * (src_i - src_0) = dst_i - dst_0, i.e. M = D * Sinv
    Sm = [ag.sub(src[i], src[0]) for i in (1, 2, 3)]   # columns
    Dm = [ag.sub(dst[i], dst[0]) for i in (1, 2, 3)]
    Sinv = inv3(F, [[Sm[j][i] for j in range(3)] for i in range(3)])
    if Sinv is None: return None
    Dmat = [[Dm[j][i] for j in range(3)] for i in range(3)]
    M = matmul(F, Dmat, Sinv)
    if det3(F, M) == 0: return None
    def ap(p):
        v = ag.sub(p, src[0])
        w = tuple(_dot(F, M[i], v) for i in range(3))
        return ag.add(w, dst[0])
    return [ag.pidx[ap(p)] for p in ag.points], M

def _dot(F, r, v):
    s = 0
    for a, b in zip(r, v): s = F.add[s][F.mul[a][b]]
    return s
def matmul(F, A, B):
    return [[_dot(F, A[i], [B[k][j] for k in range(3)]) for j in range(3)] for i in range(3)]
def det3(F, M):
    a = F.add; m = F.mul; n = F.neg
    def t(i, j, k): return m[M[0][i]][m[M[1][j]][M[2][k]]]
    s = 0
    for (i, j, k, sg) in [(0,1,2,1),(1,2,0,1),(2,0,1,1),(0,2,1,-1),(1,0,2,-1),(2,1,0,-1)]:
        v = t(i, j, k); s = a[s][v if sg > 0 else n[v]]
    return s
def inv3(F, M):
    d = det3(F, M)
    if d == 0: return None
    di = F.inv[d]; m = F.mul; a = F.add; n = F.neg
    def cof(i, j):
        r = [x for x in range(3) if x != i]; c = [x for x in range(3) if x != j]
        v = a[m[M[r[0]][c[0]]][M[r[1]][c[1]]]][n[m[M[r[0]][c[1]]][M[r[1]][c[0]]]]]
        return v if (i + j) % 2 == 0 else n[v]
    return [[m[di][cof(j, i)] for j in range(3)] for i in range(3)]

def planes_of(ag):
    F = ag.F
    out = []
    for nrm in ag.dirs:
        cl = {}
        for p in ag.points:
            cl.setdefault(_dot(F, nrm, p), set()).add(ag.pidx[p])
        for c, pts in cl.items():
            out.append((nrm, c, frozenset(pts)))
    return out

def analyze(ag, S):
    S = frozenset(S)
    inside = [l for l, (di, pts) in enumerate(ag.lines) if pts <= S]
    deg = {p: 0 for p in S}
    for l in inside:
        for p in ag.lines[l][1]: deg[p] += 1
    degdist = {}
    for v in deg.values(): degdist[v] = degdist.get(v, 0) + 1
    dirs_cov = {}
    for l in inside: dirs_cov[ag.lines[l][0]] = dirs_cov.get(ag.lines[l][0], 0) + 1
    multi = {k: v for k, v in dirs_cov.items() if v > 1}
    secs = {}
    for nrm, c, pts in planes_of(ag):
        k = len(pts & S); secs[k] = secs.get(k, 0) + 1
    linesec = {}
    for di, pts in ag.lines:
        k = len(pts & S); linesec[k] = linesec.get(k, 0) + 1
    return dict(size=len(S), lines_inside=len(inside), dirs_with_multiple_lines=len(multi),
                degree_dist=dict(sorted(degdist.items())), plane_section_dist=dict(sorted(secs.items())),
                line_section_dist=dict(sorted(linesec.items())), deg=deg)

def stabilizer(ag, S, deg):
    """All affine maps g with g(S)=S. Returns list of (perm, M, b)."""
    S = frozenset(S); Sl = sorted(S)
    P = ag.points
    # profile per point: (degree, sorted line-section sizes through it)
    lines_through = {p: [] for p in range(len(P))}
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts: lines_through[p].append(l)
    sec = [len(pts & S) for di, pts in ag.lines]
    prof = {p: (deg[p], tuple(sorted(sec[l] for l in lines_through[p]))) for p in Sl}
    def line_sec(p, r):
        return sec[ag.line_through(P[p], ag.sub(P[r], P[p]))]
    # choose a source frame of affinely independent points of S with rarest profiles
    cnt = {}
    for p in Sl: cnt[prof[p]] = cnt.get(prof[p], 0) + 1
    order = sorted(Sl, key=lambda p: cnt[prof[p]])
    src = None
    for p0 in order:
        for p1 in order:
            if p1 == p0: continue
            for p2 in order:
                if p2 in (p0, p1): continue
                for p3 in order:
                    if p3 in (p0, p1, p2): continue
                    r = affine_map(ag, [P[p0], P[p1], P[p2], P[p3]], [P[p0], P[p1], P[p2], P[p3]])
                    if r is not None:
                        src = (p0, p1, p2, p3); break
                if src: break
            if src: break
        if src: break
    p0, p1, p2, p3 = src
    cand = {i: [r for r in Sl if prof[r] == prof[i]] for i in src}
    # plane-section invariant for triples
    def plane_sec(a, b, c):
        u, v = ag.sub(P[b], P[a]), ag.sub(P[c], P[a])
        pts = {ag.pidx[ag.add(P[a], ag.add(ag.scale(s, u), ag.scale(t, v)))] for s in range(ag.q) for t in range(ag.q)}
        return len(pts & S) if len(pts) == ag.q ** 2 else -1
    ls01, ls02, ls12, ls03, ls13, ls23 = (line_sec(p0, p1), line_sec(p0, p2), line_sec(p1, p2),
                                          line_sec(p0, p3), line_sec(p1, p3), line_sec(p2, p3))
    ps012 = plane_sec(p0, p1, p2)
    stab = []
    for r0 in cand[p0]:
        for r1 in cand[p1]:
            if r1 == r0 or line_sec(r0, r1) != ls01: continue
            for r2 in cand[p2]:
                if r2 in (r0, r1) or line_sec(r0, r2) != ls02 or line_sec(r1, r2) != ls12: continue
                if plane_sec(r0, r1, r2) != ps012: continue
                for r3 in cand[p3]:
                    if r3 in (r0, r1, r2) or line_sec(r0, r3) != ls03 or line_sec(r1, r3) != ls13 or line_sec(r2, r3) != ls23: continue
                    res = affine_map(ag, [P[p0], P[p1], P[p2], P[p3]], [P[r0], P[r1], P[r2], P[r3]])
                    if res is None: continue
                    perm, M = res
                    if all(perm[p] in S for p in Sl):
                        b = P[perm[ag.pidx[(0, 0, 0)]]]
                        stab.append((perm, M, b))
    return stab

def perm_order(perm):
    n = len(perm); seen = [False] * n; from math import lcm
    o = 1
    for i in range(n):
        if not seen[i]:
            j, c = i, 0
            while not seen[j]:
                seen[j] = True; j = perm[j]; c += 1
            o = lcm(o, c)
    return o

if __name__ == "__main__":
    fn = sys.argv[1]
    d = json.load(open(fn))
    ag = AG(d["n"], d["q"])
    S = [ag.pidx[tuple(p)] for p in d["points"]]
    t0 = time.time()
    a = analyze(ag, S)
    deg = a.pop("deg")
    print(json.dumps(a))
    if "--nostab" not in sys.argv:
        st = stabilizer(ag, S, deg)
        orders = {}
        for perm, M, b in st:
            o = perm_order(perm); orders[o] = orders.get(o, 0) + 1
        ntrans = sum(1 for perm, M, b in st if M == [[1,0,0],[0,1,0],[0,0,1]])
        print(json.dumps({"stab_order": len(st), "element_orders": dict(sorted(orders.items())),
                          "translations": ntrans, "time": round(time.time() - t0, 1)}))
        out = fn.replace(".json", "_stab.json").replace("out/", "out/").replace("solutions/", "out/")
        json.dump([{"M": M, "b": list(b)} for perm, M, b in st], open(out, "w"))
        print("wrote", out)
