"""Session V. Exact planar tables for AG(2,q):
  F(q, a)        = min |union of a lines| over unions that contain >= 1 line of every direction, 6 <= a <= q(q+1)
  Ft(q, tuple)   = min |union| when the (sorted) numbers of lines per direction are the given (q+1)-tuple (entries 1..q)
  G(n, lam, q)   = min sum_x C(j_x,2) over weights j_x in 1..q on lam points with sum n and at least one weight = q
Computed by CP-SAT (tiny) and cached in out/union_tables_q{q}.json."""
import json, os, sys, itertools, math
from ortools.sat.python import cp_model
from geom import AG
TIME_LIMIT = 120

def _solve_union(ag, per_dir=None, total=None):
    m = cp_model.CpModel()
    x = [m.NewBoolVar("") for _ in ag.points]; y = [m.NewBoolVar("") for _ in ag.lines]
    counts = []
    for di in range(len(ag.dirs)):
        count = m.NewIntVar(1, ag.q, "")
        m.Add(count == sum(y[l] for l in ag.lines_by_dir[di]))
        counts.append(count)
    if per_dir is not None:
        # Ft is indexed by a SORTED multiset, so minimize over every placement on directions.
        # PGL(2,q) need not act as the full symmetric group on q+1 directions.
        assert len(per_dir) == len(ag.dirs)
        for value in sorted(set(per_dir)):
            indicators = [m.NewBoolVar("") for _ in counts]
            for count, indicator in zip(counts, indicators):
                m.Add(count == value).OnlyEnforceIf(indicator)
                m.Add(count != value).OnlyEnforceIf(indicator.Not())
            m.Add(sum(indicators) == per_dir.count(value))
    if total is not None: m.Add(sum(y) == total)
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts: m.AddImplication(y[l], x[p])
    m.Minimize(sum(x))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = TIME_LIMIT
    st = s.Solve(m); assert st in (cp_model.OPTIMAL, cp_model.FEASIBLE), s.StatusName(st)
    if st == cp_model.OPTIMAL: return int(s.ObjectiveValue())
    return int(math.ceil(s.BestObjectiveBound() - 1e-9))   # proven lower bound only (valid, possibly weaker); q=7 large tuples

def tables(q):
    fn = f"out/union_tables_q{q}.json"
    if os.path.exists(fn):
        d = json.load(open(fn))
        if q <= 4 or d.get("unordered_counts_version") == 2:
            return {int(k): v for k, v in d["F"].items()}, {tuple(map(int, k.split(","))): v for k, v in d["Ft"].items()}
        # Older q>=5 tables fixed one ordering of a sorted tuple; regenerate safely.
    ag = AG(2, q); nd = q + 1
    F = {a: _solve_union(ag, total=a) for a in range(nd, q * nd + 1)}
    Ft = {}
    for t in itertools.combinations_with_replacement(range(1, q + 1), nd):
        Ft[t] = _solve_union(ag, per_dir=list(t))
    json.dump({"unordered_counts_version": 2, "F": F, "Ft": {",".join(map(str, k)): v for k, v in Ft.items()}}, open(fn, "w"))
    return F, Ft

def G(n, lam, q):
    """min sum C(j,2) for lam weights in [1,q] summing to n with one weight forced = q; None if impossible."""
    rest_n, rest_l = n - q, lam - 1
    if rest_l < 0 or rest_n < rest_l or rest_n > q * rest_l: return None
    base, extra = divmod(rest_n, rest_l) if rest_l else (0, 0)
    return math.comb(q, 2) + extra * math.comb(base + 1, 2) + (rest_l - extra) * math.comb(base, 2)

if __name__ == "__main__":
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    F, Ft = tables(q)
    print("F(a):", F)
    byF = {}
    for t, v in Ft.items(): byF.setdefault(v, []).append(t)
    for v in sorted(byF): print(f"  Ft = {v}: {len(byF[v])} tuples, e.g. {byF[v][:4]}")
    # sanity: F(a) == min over tuples with sum a
    for a in F: assert F[a] == min(v for t, v in Ft.items() if sum(t) == a), a
    print("F consistent with Ft; G(46,lam,5):", {lam: G(46, lam, 5) for lam in range(17, 26)})
