"""CP-SAT subcase model (max-degree WLOG, as maxdeg.py) + the session-IV counting identities as redundant constraints:
  k_pi = |K ∩ pi|, r_pi = # witness lines in pi, lam_L = |K ∩ L| (all 775 lines), d_p = degree.
  (I2) sum_pi k = 31 n, sum_pi r = 186;  (I3) sum_pi C(r,2) = sum_p C(d_p,2);  (I4) sum_pi k^2 = 6n^2 + 25n;
  (I5) sum_pi k r = 31 n + 775;  (I6) sum_pi C(k,3) = C(n,3) + 5 sum_L C(lam_L,3);  sum_L C(lam,2) = C(n,2);
  planar cut k >= f(r); projection cut; exact size sum x = n.
Usage: maxdeg_ident.py --q 5 --n 46 --D 4 --rep 0 [--time 600] [--no-ident] [--workers 1]"""
import argparse, json, math, time
from ortools.sat.python import cp_model
from geom import AG, is_kakeya
from orbits_fast import load_or_compute
from maxdeg import describe
from structure import planes_of
from plane_f import get

def build(ag, D, R, n, ident=True, planes=True, proj=True, split=None):
    """split=(direction index e, line index k in lines_by_dir[e]): also force that witness line (25-way cube)."""
    q = ag.q; nd = len(ag.dirs); f = get(q)
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(ag.points))]
    y = [m.NewBoolVar(f"y{i}") for i in range(len(ag.lines))]
    for di in range(nd):
        m.AddExactlyOne([y[l] for l in ag.lines_by_dir[di]])
    thru = [[] for _ in ag.points]
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts:
            m.AddImplication(y[l], x[p]); thru[p].append(l)
    for p in range(len(ag.points)):
        m.AddBoolOr([y[l] for l in thru[p]]).OnlyEnforceIf(x[p])
        m.Add(sum(y[l] for l in thru[p]) <= D)
    O = (0, 0, 0)
    for d in R: m.Add(y[ag.line_through(O, ag.dirs[d])] == 1)
    m.Add(sum(x) == n)
    if split is not None:
        e, k = split; m.Add(y[ag.lines_by_dir[e][k]] == 1)
    PL = planes_of(ag)   # (normal, c, frozenset pts)
    # lines inside each plane
    line_in_plane = []
    for nrm, c, pts in PL:
        line_in_plane.append([l for l, (di, lp) in enumerate(ag.lines) if lp <= pts])
    ks, rs = [], []
    sq = [i * i for i in range(q * q + 1)]; c3 = [math.comb(i, 3) for i in range(q * q + 1)]; c2 = [math.comb(i, 2) for i in range(q * q + 2)]
    for (nrm, c, pts), lines in zip(PL, line_in_plane):
        k = m.NewIntVar(1, q * q, ""); m.Add(k == sum(x[p] for p in pts))
        r = m.NewIntVar(0, q + 1, ""); m.Add(r == sum(y[l] for l in lines))
        ks.append(k); rs.append(r)
        if planes:
            t = m.NewIntVar(0, q * q, ""); m.AddElement(r, [f[i] for i in range(q + 2)], t); m.Add(k >= t)
        # degree-sum identity per plane (linear, implied): k <= 25 + 5r <= D k
        m.Add(k <= (nd - q - 1) + q * r); m.Add((nd - q - 1) + q * r <= D * k)
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
            dp = m.NewIntVar(0, D, ""); m.Add(dp == sum(y[l] for l in thru[p]))
            g = m.NewIntVar(0, math.comb(D, 2), ""); m.AddElement(dp, [math.comb(i, 2) for i in range(D + 1)], g); dc2.append(g)
        m.Add(sum(rc2) == sum(dc2))
        lc3, lc2 = [], []
        for l, (di, pts) in enumerate(ag.lines):
            lam = m.NewIntVar(0, q, ""); m.Add(lam == sum(x[p] for p in pts))
            e = m.NewIntVar(0, math.comb(q, 3), ""); m.AddElement(lam, [math.comb(i, 3) for i in range(q + 1)], e); lc3.append(e)
            g = m.NewIntVar(0, math.comb(q, 2), ""); m.AddElement(lam, [math.comb(i, 2) for i in range(q + 1)], g); lc2.append(g)
        m.Add(sum(kc3) == math.comb(n, 3) + q * sum(lc3))
        m.Add(sum(lc2) == math.comb(n, 2))
    return m, x, y

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
    ap.add_argument("--D", type=int, required=True); ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--time", type=float, default=600); ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--no-ident", action="store_true"); ap.add_argument("--log", action="store_true"); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    ag = AG(3, a.q)
    lv = load_or_compute(a.q, a.D, verbose=False)
    rep, sz, stab = lv[a.D][a.rep]
    m, x, y = build(ag, a.D, rep, a.n, ident=not a.no_ident)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = a.time; s.parameters.num_workers = a.workers
    s.parameters.log_search_progress = a.log
    t0 = time.time(); st = s.Solve(m); dt = time.time() - t0
    res = dict(q=a.q, n=a.n, D=a.D, rep=a.rep, dirs=[ag.dirs[d] for d in rep], orbit_size=sz, ident=not a.no_ident,
               status=s.StatusName(st), time=round(dt, 1), workers=a.workers)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = [i for i in range(len(ag.points)) if s.Value(x[i])]
        ok, miss = is_kakeya(ag, S); assert ok, miss
        res["size"] = len(S); res["points"] = [ag.points[i] for i in S]
    print(json.dumps({k: v for k, v in res.items() if k != "points"}), flush=True)
    if a.out: json.dump(res, open(a.out, "w"))
