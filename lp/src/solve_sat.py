"""Independent exact method: incremental SAT (CaDiCaL via python-sat) with a
totalizer cardinality encoding. Decreases k until UNSAT: proves min = last SAT k.
Usage: python solve_sat.py --n 2 --q 5 [--fix 1] [--start K]"""
import argparse, json, time
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from geom import AG, is_kakeya

def run(n, q, fix=1, case="meet", start=None, solver="cadical153", time_limit=None):
    ag = AG(n, q)
    pool = IDPool()
    x = [pool.id(("x", i)) for i in range(len(ag.points))]
    y = [pool.id(("y", i)) for i in range(len(ag.lines))]
    cls = []
    for di in range(len(ag.dirs)):
        cls.append([y[l] for l in ag.lines_by_dir[di]])
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts:
            cls.append([-y[l], x[p]])
    origin = tuple([0] * n)
    if fix >= 1:
        cls.append([y[ag.line_through(origin, ag.unit(0))]])
    if fix >= 2:
        if case == "meet":
            cls.append([y[ag.line_through(origin, ag.unit(1))]])
        else:
            cls.append([y[ag.line_through(ag.unit(2), ag.unit(1))]])
    k = start if start is not None else len(ag.points)
    card = CardEnc.atmost(lits=x, bound=k, vpool=pool, encoding=EncType.totalizer)
    # totalizer: atmost(k) clauses with outputs; we rebuild per k for simplicity (cheap instances)
    best = None
    t0 = time.time()
    while True:
        with Solver(name=solver, bootstrap_with=cls) as s:
            card = CardEnc.atmost(lits=x, bound=k, vpool=pool, encoding=EncType.totalizer)
            for c in card.clauses:
                s.add_clause(c)
            sat = s.solve()
            if sat:
                model = set(v for v in s.get_model() if v > 0)
                S = [i for i in range(len(ag.points)) if x[i] in model]
                ok, miss = is_kakeya(ag, S); assert ok
                best = S
                k = len(S) - 1
                print(f"  SAT with {len(S)} points ({time.time()-t0:.1f}s); trying {k}")
            else:
                print(f"  UNSAT at {k} ({time.time()-t0:.1f}s)")
                break
    return {"n": n, "q": q, "fix": fix, "case": case, "min": (len(best) if best else None), "unsat_at": k,
            "points": ([ag.points[i] for i in best] if best else None), "time": round(time.time() - t0, 1)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--fix", type=int, default=1)
    ap.add_argument("--case", default="meet")
    ap.add_argument("--start", type=int, default=None)
    a = ap.parse_args()
    r = run(a.n, a.q, a.fix, a.case, a.start)
    print(json.dumps({k: v for k, v in r.items() if k != "points"}))
    json.dump(r, open(f"out/sat_n{a.n}_q{a.q}_{a.case}.json", "w"))
