"""Exact minimum Kakeya set in AG(n,q) with OR-tools CP-SAT.
Variables: x_p (point chosen), y_l (line l chosen as witness for its direction).
Constraints: each direction has >= 1 chosen line; chosen line => its points chosen.
Symmetry breaking (all WLOG under AGL(n,q), see WRITEUP.md):
  --fix 1 : the line through origin in direction e_1 is a witness.
  --fix 2 : additionally the witness of e_2 is the line through origin ('meet') OR
            the line {t e_2 + e_3} ('skew', n>=3). For n=2 only 'meet' exists.
  --fix 3 : (n>=3) additionally the witness of e_3 is {t e_3 + (a,b,0)} with
            (a,b) in {00,10,01,11} (--sub). The residual group is diag(l,m,1)-ish.
  --case both / --sub all : run all cases, take min; proof only if every case OPTIMAL.
Optional strengthening:
  --planes : for every affine plane pi, sum_{p in pi} x_p >= f(#witness lines inside pi),
             f = exact planar table from plane_f.py (valid, redundant).
Writes out/sol_n{n}_q{q}.json (or --out)."""
import argparse, json, time, sys, itertools
from ortools.sat.python import cp_model
from geom import AG, is_kakeya

def planes_of(ag):
    """All affine planes of AG(3,q) as (frozenset of point idx, list of line idx inside)."""
    assert ag.n == 3
    q = ag.q
    # plane direction = 2-dim subspace; enumerate via its normal (projective point) using dot product
    F = ag.F
    def dot(u, v):
        s = 0
        for a, b in zip(u, v):
            s = F.add[s][F.mul[a][b]]
        return s
    planes = []
    for nrm in ag.dirs:
        classes = {}
        for p in ag.points:
            classes.setdefault(dot(nrm, p), set()).add(ag.pidx[p])
        for c, pts in classes.items():
            pts = frozenset(pts)
            lines = [l for l, (di, lp) in enumerate(ag.lines) if lp <= pts]
            assert len(pts) == q * q and len(lines) == q * (q + 1)
            planes.append((pts, lines))
    assert len(planes) == q * (q * q + q + 1)
    return planes

def build(ag, fix, case, sub, ub=None, lb_hint=None, planes=False, exact_one=False):
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(ag.points))]
    y = [m.NewBoolVar(f"y{i}") for i in range(len(ag.lines))]
    for di in range(len(ag.dirs)):
        if exact_one:
            m.AddExactlyOne([y[l] for l in ag.lines_by_dir[di]])
        else:
            m.AddBoolOr([y[l] for l in ag.lines_by_dir[di]])
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts:
            m.AddImplication(y[l], x[p])
    for p in range(len(ag.points)):
        m.AddBoolOr([y[l] for l, (di, pts) in enumerate(ag.lines) if p in pts]).OnlyEnforceIf(x[p])
    n = ag.n
    origin = tuple([0] * n)
    fixed = []
    if fix >= 1:
        fixed.append(ag.line_through(origin, ag.unit(0)))
    if fix >= 2:
        if case == "meet":
            fixed.append(ag.line_through(origin, ag.unit(1)))
        elif case == "skew":
            assert n >= 3
            fixed.append(ag.line_through(ag.unit(2), ag.unit(1)))
        else:
            raise ValueError(case)
    if fix >= 3:
        assert n >= 3
        a, b = int(sub[0]), int(sub[1])
        fixed.append(ag.line_through((a, b) + (0,) * (n - 2), ag.unit(2)))
    for l in fixed:
        m.Add(y[l] == 1)
    if planes:
        from plane_f import get
        f = get(ag.q)
        table = [f[r] for r in range(ag.q + 2)]
        for pts, lines in planes_of(ag):
            z = m.NewIntVar(0, ag.q + 1, "")
            m.Add(z == sum(y[l] for l in lines))
            t = m.NewIntVar(0, ag.q * ag.q, "")
            m.AddElement(z, table, t)
            m.Add(sum(x[p] for p in pts) >= t)
    size = sum(x)
    if ub is not None:
        m.Add(size <= ub)
    if lb_hint is not None:
        m.Add(size >= lb_hint)
    m.Minimize(size)
    return m, x, y

def solve(n, q, fix=1, case="meet", sub="00", time_limit=600, workers=4, ub=None, lb_hint=None,
          log=False, planes=False, core=False, exact_one=False, hint=None):
    ag = AG(n, q)
    m, x, y = build(ag, fix, case, sub, ub, lb_hint, planes, exact_one)
    if hint:
        hs = set(map(tuple, hint))
        for i, p in enumerate(ag.points):
            m.AddHint(x[i], int(p in hs))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = workers
    s.parameters.log_search_progress = log
    s.parameters.symmetry_level = 2
    if core:
        s.parameters.optimize_with_core = True
    t0 = time.time()
    st = s.Solve(m)
    dt = time.time() - t0
    name = s.StatusName(st)
    res = {"n": n, "q": q, "fix": fix, "case": case, "sub": sub if fix >= 3 else None,
           "planes": planes, "core": core, "status": name, "time": round(dt, 1)}
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = [i for i in range(len(ag.points)) if s.Value(x[i])]
        ok, missing = is_kakeya(ag, S)
        assert ok, missing
        res.update(size=len(S), bound=int(round(s.BestObjectiveBound())), points=[ag.points[i] for i in S],
                   lines=[sorted(ag.lines[l][1]) for l in range(len(ag.lines)) if s.Value(y[l])])
    else:
        res.update(size=None, bound=int(round(s.BestObjectiveBound())))
    return res

def strip(r):
    return {k: v for k, v in r.items() if k not in ("points", "lines")}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--fix", type=int, default=1)
    ap.add_argument("--case", default="meet")
    ap.add_argument("--sub", default="00")
    ap.add_argument("--time", type=float, default=600)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--ub", type=int, default=None)
    ap.add_argument("--lb", type=int, default=None)
    ap.add_argument("--planes", action="store_true")
    ap.add_argument("--core", action="store_true")
    ap.add_argument("--exact-one", action="store_true")
    ap.add_argument("--hint", default=None, help="json solution file to use as hint")
    ap.add_argument("--log", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    cases = ["meet", "skew"] if a.case == "both" else [a.case]
    subs = ["00", "10", "01", "11"] if (a.fix >= 3 and a.sub == "all") else [a.sub]
    hint = json.load(open(a.hint))["points"] if a.hint else None
    results = []
    for c in cases:
        for sb in subs:
            r = solve(a.n, a.q, a.fix, c, sb, a.time, a.workers, a.ub, a.lb, a.log, a.planes, a.core,
                      a.exact_one, hint)
            print(json.dumps(strip(r)), flush=True)
            results.append(r)
    feas = [r for r in results if r.get("size") is not None]
    best = dict(min(feas, key=lambda r: r["size"])) if feas else dict(results[0])
    if len(results) > 1:
        best["case"] = a.case; best["sub"] = a.sub
        best["status"] = "OPTIMAL" if all(r["status"] == "OPTIMAL" for r in results) else "FEASIBLE"
        best["bound"] = min(r["bound"] for r in results)
        best["time"] = sum(r["time"] for r in results)
        best["cases"] = [strip(r) for r in results]
        print(json.dumps(strip(best)), flush=True)
    out = a.out or f"out/sol_n{a.n}_q{a.q}.json"
    json.dump(best, open(out, "w"))
    print("wrote", out)
