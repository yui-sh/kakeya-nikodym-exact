"""f(q, r) = minimum size of the union of r lines with pairwise distinct directions in AG(2,q).
Computed exactly with CP-SAT (tiny instances); cached in out/plane_f.json."""
import json, os, sys
from ortools.sat.python import cp_model
from geom import AG

def compute(q):
    ag = AG(2, q)
    f = {0: 0}
    for r in range(1, q + 2):
        m = cp_model.CpModel()
        x = [m.NewBoolVar("") for _ in ag.points]
        y = [m.NewBoolVar("") for _ in ag.lines]
        for di in range(len(ag.dirs)):
            m.Add(sum(y[l] for l in ag.lines_by_dir[di]) <= 1)
        m.Add(sum(y) == r)
        for l, (di, pts) in enumerate(ag.lines):
            for p in pts:
                m.AddImplication(y[l], x[p])
        m.Add(y[ag.line_through((0, 0), (1, 0))] == 1)  # WLOG
        m.Minimize(sum(x))
        s = cp_model.CpSolver(); s.parameters.num_workers = 4; s.parameters.max_time_in_seconds = 120
        st = s.Solve(m)
        assert st == cp_model.OPTIMAL, s.StatusName(st)
        f[r] = int(s.ObjectiveValue())
    return f

def get(q):
    fn = "out/plane_f.json"
    cache = json.load(open(fn)) if os.path.exists(fn) else {}
    if str(q) not in cache:
        cache[str(q)] = compute(q)
        json.dump(cache, open(fn, "w"))
    return {int(k): v for k, v in cache[str(q)].items()}

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        print(q, get(q))
