"""Alternative WLOG for small n (session IV): a plane with the maximum number r_max of witness lines.
Lemma: sum_pi C(r_pi,2) = sum_p C(d_p,2) >= B(n) (balanced degree sequence with sum 155 on n points) and
sum_pi C(r_pi,2) <= 186 (r_max - 1)/2, so r_max >= 1 + 2B(n)/186 (n = 46: B = 189 -> r_max >= 4).
WLOG (AGL(3,5) is transitive on planes; the stabiliser of a plane induces AGL(2,5) on it) the plane z = 0 contains a
fixed representative C of the AGL(2,5)-orbit of the r_max-line configuration, and every plane has r_pi <= r_max.
Subcases: (r_max, C) for r_max from the lemma bound to 6. Model: as maxdeg_ident (identities, cuts), forced lines C,
cap r <= r_max on every plane, optional degree cap D (then the subcase is (r_max, C, D) with some point of degree D).
Usage: rmax_model.py --q 5 --n 46 --list ;  rmax_model.py --q 5 --n 46 --r 4 --idx 0 [--D 4] [--time 1800]"""
import argparse, json, math, time, itertools, sys
import numpy as np
from ortools.sat.python import cp_model
from geom import AG, is_kakeya
from structure import planes_of
from plane_f import get

def balanced_min_pairs(n, total):
    a, b = divmod(total, n)   # n-b points of degree a, b points of degree a+1
    return (n - b) * math.comb(a, 2) + b * math.comb(a + 1, 2)

def rmax_lower_bound(n, q=5):
    nd = q * q + q + 1; B = balanced_min_pairs(n, nd * q); S = (q + 1) * nd
    r = 1
    while S * (r - 1) / 2 < B: r += 1
    return r, B

def plane_configs(ag, r):
    """AGL(2,q)-orbit representatives of r-sets of lines with distinct directions in the plane z=0.
    Lines of z=0 as point-index sets; the group = all affine maps of the plane (enumerated: GL(2,q) x translations)."""
    q = ag.q; F = ag.F
    pl_pts = [p for p in ag.points if p[2] == 0]; pidx = {p: i for i, p in enumerate(pl_pts)}
    lines = [(di, pts) for di, pts in ag.lines if all(ag.points[p][2] == 0 for p in pts)]
    assert len(lines) == q * (q + 1)
    line_mask = {}
    for li, (di, pts) in enumerate(lines):
        line_mask[li] = sum(1 << pidx[ag.points[p]] for p in pts)
    # group elements as permutations of the 25 plane points
    perms = []
    for a_, b_, c_, d_ in itertools.product(range(q), repeat=4):
        if F.add[F.mul[a_][d_]][F.neg[F.mul[b_][c_]]] == 0: continue
        for tx, ty in itertools.product(range(q), repeat=2):
            perm = []
            for (x, y, _) in pl_pts:
                nx = F.add[F.add[F.mul[a_][x]][F.mul[b_][y]]][tx]; ny = F.add[F.add[F.mul[c_][x]][F.mul[d_][y]]][ty]
                perm.append(pidx[(nx, ny, 0)])
            perms.append(perm)
    # act on the 30 lines of the plane (canonical form on LINE sets: a union of lines may contain extra full lines)
    Gp = np.array(perms, dtype=np.int64)
    lmask_to_idx = {line_mask[li]: li for li in range(len(lines))}
    Gl = np.zeros((len(perms), len(lines)), dtype=np.int64)
    pts_of = {li: [i for i in range(len(pl_pts)) if (line_mask[li] >> i) & 1] for li in range(len(lines))}
    pow2p = (np.int64(1) << np.arange(len(pl_pts), dtype=np.int64))
    for li in range(len(lines)):
        im = pow2p[Gp[:, pts_of[li]]].sum(axis=1)
        Gl[:, li] = [lmask_to_idx[int(v)] for v in im]
    pow2l = (np.int64(1) << np.arange(len(lines), dtype=np.int64))
    def canon(lset):
        im = pow2l[Gl[:, list(lset)]].sum(axis=1)
        mask = sum(1 << li for li in lset)
        return int(im.min()), int((im == mask).sum())
    bydir = {}
    for li, (di, pts) in enumerate(lines): bydir.setdefault(di, []).append(li)
    dirs = sorted(bydir)
    seen = {}; reps = []
    for dset in itertools.combinations(dirs, r):
        for choice in itertools.product(*[bydir[d] for d in dset]):
            c, stab = canon(choice)
            if c in seen: continue
            seen[c] = len(reps)
            umask = 0
            for li in choice: umask |= line_mask[li]
            reps.append(dict(lines=[sorted(ag.points[p] for p in lines[li][1]) for li in choice], line_idx=[ag.line_of[lines[li][1]] for li in choice],
                             union=bin(umask).count("1"), orbit_size=len(perms) // stab))
    tot = sum(x["orbit_size"] for x in reps); assert tot == math.comb(q + 1, r) * q ** r, (tot, math.comb(q + 1, r) * q ** r)
    return reps

def build(ag, forced_lines, rmax, n, D=None, ident=True, planes=True, proj=True, split=None, mvec=None):
    q = ag.q; nd = len(ag.dirs); f = get(q)
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(ag.points))]
    y = [m.NewBoolVar(f"y{i}") for i in range(len(ag.lines))]
    for di in range(nd): m.AddExactlyOne([y[l] for l in ag.lines_by_dir[di]])
    thru = [[] for _ in ag.points]
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts: m.AddImplication(y[l], x[p]); thru[p].append(l)
    degs = []
    for p in range(len(ag.points)):
        m.AddBoolOr([y[l] for l in thru[p]]).OnlyEnforceIf(x[p])
        dp = m.NewIntVar(0, min(nd, D) if D else nd, ""); m.Add(dp == sum(y[l] for l in thru[p])); degs.append(dp)
    if D:
        bs = []
        for p in range(len(ag.points)):
            b = m.NewBoolVar(""); m.Add(degs[p] >= D).OnlyEnforceIf(b); bs.append(b)
        m.AddBoolOr(bs)
    Dcap = D if D else nd
    if mvec is not None:   # session V: fixed degree vector (m_1..m_D): exactly m_j points of degree j
        assert D and len(mvec) == D and sum(mvec) == n
        isdeg = [[m.NewBoolVar("") for j in range(D + 1)] for p in range(len(ag.points))]
        for p in range(len(ag.points)):
            m.AddExactlyOne(isdeg[p]); m.Add(degs[p] == sum(j * isdeg[p][j] for j in range(D + 1)))
        for j in range(1, D + 1): m.Add(sum(isdeg[p][j] for p in range(len(ag.points))) == mvec[j - 1])
    for l in forced_lines: m.Add(y[l] == 1)
    m.Add(sum(x) == n)
    if split is not None:
        e, k = split; m.Add(y[ag.lines_by_dir[e][k]] == 1)
    PL = planes_of(ag)
    line_in_plane = [[l for l, (di, lp) in enumerate(ag.lines) if lp <= pts] for nrm, c, pts in PL]
    ks, rs = [], []
    sq = [i * i for i in range(q * q + 1)]; c3 = [math.comb(i, 3) for i in range(q * q + 1)]; c2 = [math.comb(i, 2) for i in range(q * q + 2)]
    for (nrm, c, pts), lines in zip(PL, line_in_plane):
        k = m.NewIntVar(1, q * q, ""); m.Add(k == sum(x[p] for p in pts))
        r = m.NewIntVar(0, rmax, ""); m.Add(r == sum(y[l] for l in lines))
        ks.append(k); rs.append(r)
        if planes:
            t = m.NewIntVar(0, q * q, ""); m.AddElement(r, [f[i] for i in range(q + 2)], t); m.Add(k >= t)
        m.Add(k <= (nd - q - 1) + q * r)
        if D: m.Add((nd - q - 1) + q * r <= D * k)
    if proj:
        bm = f[q + 1]
        for di in range(nd):
            zs = []
            for l in ag.lines_by_dir[di]:
                z = m.NewBoolVar(""); pts = ag.lines[l][1]
                m.AddBoolOr([x[p] for p in pts]).OnlyEnforceIf(z)
                for p in pts: m.AddImplication(x[p], z)
                zs.append(z)
            m.Add(sum(zs) >= bm)
    if ident:
        ksq, kr, kc3, rc2 = [], [], [], []
        for k, r in zip(ks, rs):
            a = m.NewIntVar(0, q ** 4, ""); m.AddElement(k, sq, a); ksq.append(a)
            b = m.NewIntVar(0, q * q * (q + 1), ""); m.AddMultiplicationEquality(b, [k, r]); kr.append(b)
            e = m.NewIntVar(0, c3[-1], ""); m.AddElement(k, c3, e); kc3.append(e)
            g = m.NewIntVar(0, c2[q + 1], ""); m.AddElement(r, c2[:q + 2], g); rc2.append(g)
        m.Add(sum(ksq) == (q + 1) * n * n + (nd - q - 1) * n)
        m.Add(sum(kr) == nd * n + q * q * nd)
        dc2 = []
        for p in range(len(ag.points)):
            g = m.NewIntVar(0, math.comb(Dcap, 2), ""); m.AddElement(degs[p], [math.comb(i, 2) for i in range(Dcap + 1)], g); dc2.append(g)
        m.Add(sum(rc2) == sum(dc2))
        lc3, lc2 = [], []
        for l, (di, pts) in enumerate(ag.lines):
            lam = m.NewIntVar(0, q, ""); m.Add(lam == sum(x[p] for p in pts))
            e = m.NewIntVar(0, math.comb(q, 3), ""); m.AddElement(lam, [math.comb(i, 3) for i in range(q + 1)], e); lc3.append(e)
            g = m.NewIntVar(0, math.comb(q, 2), ""); m.AddElement(lam, [math.comb(i, 2) for i in range(q + 1)], g); lc2.append(g)
        m.Add(sum(kc3) == math.comb(n, 3) + q * sum(lc3))
        m.Add(sum(lc2) == math.comb(n, 2))
    return m, x, y

def validate(ag, fn, tl=300):
    """model must accept a known set: force the witness lines of one of its max-r planes, cap r by that r_max, hint the set."""
    d = json.load(open(fn)); S = {ag.pidx[tuple(p)] for p in d["points"]}
    W = [next(l for l in ag.lines_by_dir[di] if ag.lines[l][1] <= S) for di in range(len(ag.dirs))]
    best = None
    for nrm, c, pts in planes_of(ag):
        ins = [l for l in W if ag.lines[l][1] <= pts]
        if best is None or len(ins) > len(best): best = ins
    rmax = len(best)
    m, x, y = build(ag, best, rmax, len(S))
    for i in range(len(ag.points)): m.AddHint(x[i], int(i in S))
    for l in range(len(ag.lines)): m.AddHint(y[l], int(l in W))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 1
    t0 = time.time(); st = s.Solve(m)
    print(json.dumps(dict(validate=fn, n=len(S), rmax=rmax, status=s.StatusName(st), time=round(time.time() - t0, 1))))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
    ap.add_argument("--list", action="store_true"); ap.add_argument("--r", type=int); ap.add_argument("--idx", type=int)
    ap.add_argument("--D", type=int, default=None); ap.add_argument("--time", type=float, default=1800); ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--validate", default=None, help="solution json: the model must accept it (its own max-r plane forced)")
    ap.add_argument("--mvec", default=None, help="fixed degree vector m_1,...,m_D (session V)")
    a = ap.parse_args(); ag = AG(3, a.q)
    if a.validate:
        validate(ag, a.validate, a.time); sys.exit()
    rlo, B = rmax_lower_bound(a.n, a.q)
    if a.list:
        print(f"n={a.n}: sum_p C(d_p,2) >= {B}, so r_max >= {rlo}")
        for r in range(rlo, a.q + 2):
            reps = plane_configs(ag, r)
            print(f"r_max={r}: {len(reps)} configurations (orbit sizes sum to C(6,{r})*5^{r})")
            for i, c in enumerate(reps): print(f"   idx {i}: union {c['union']} pts, orbit {c['orbit_size']}, lines {c['lines']}")
        sys.exit()
    reps = plane_configs(ag, a.r); cfg = reps[a.idx]
    mvec = [int(v) for v in a.mvec.split(",")] if a.mvec else None
    m, x, y = build(ag, cfg["line_idx"], a.r, a.n, a.D, mvec=mvec)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = a.time; s.parameters.num_workers = a.workers
    t0 = time.time(); st = s.Solve(m); dt = time.time() - t0
    res = dict(q=a.q, n=a.n, rmax=a.r, idx=a.idx, D=a.D, mvec=mvec, union=cfg["union"], orbit_size=cfg["orbit_size"], status=s.StatusName(st), time=round(dt, 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = [i for i in range(len(ag.points)) if s.Value(x[i])]; ok, miss = is_kakeya(ag, S); assert ok
        res["size"] = len(S); res["points"] = [ag.points[i] for i in S]
    print(json.dumps({k: v for k, v in res.items() if k != "points"}), flush=True)

