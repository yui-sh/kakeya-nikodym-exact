"""Large-neighbourhood search for small Kakeya sets in AG(3,q) with CP-SAT.
State: a Kakeya set S with one witness line per direction. Each round: keep a random fraction of the
witness lines fixed (y=1), require |S| <= |S|-1, solve with a short time limit; accept if found.
Also tries 'remove-and-repair': forbid a random subset of current points.
Usage: python lns.py --q 5 --seed out/sol_n3_q5.json --rounds 100 --time 20 --keep 0.5 --out solutions/lns_q5.json"""
import argparse, json, random, time, sys
from ortools.sat.python import cp_model
from geom import AG, is_kakeya

def witness_lines(ag, S):
    S = set(S); W = []
    for di in range(len(ag.dirs)):
        ls = [l for l in ag.lines_by_dir[di] if ag.lines[l][1] <= S]
        W.append(ls)
    return W

def solve_round(ag, S, fixed_lines, forbid_pts, ub, tl, workers):
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(ag.points))]
    y = [m.NewBoolVar(f"y{i}") for i in range(len(ag.lines))]
    for di in range(len(ag.dirs)):
        m.AddBoolOr([y[l] for l in ag.lines_by_dir[di]])
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts: m.AddImplication(y[l], x[p])
    for p in range(len(ag.points)):
        m.AddBoolOr([y[l] for l, (di, pts) in enumerate(ag.lines) if p in pts]).OnlyEnforceIf(x[p])
    for l in fixed_lines: m.Add(y[l] == 1)
    for p in forbid_pts: m.Add(x[p] == 0)
    m.Add(sum(x) <= ub)
    m.Minimize(sum(x))
    Sset = set(S)
    for i in range(len(ag.points)): m.AddHint(x[i], int(i in Sset))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return [i for i in range(len(ag.points)) if s.Value(x[i])]
    return None

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--rounds", type=int, default=50)
    ap.add_argument("--time", type=float, default=20)
    ap.add_argument("--keep", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rseed", type=int, default=1)
    a = ap.parse_args()
    random.seed(a.rseed)
    ag = AG(3, a.q)
    S = [ag.pidx[tuple(p)] for p in json.load(open(a.seed))["points"]]
    assert is_kakeya(ag, S)[0]
    t0 = time.time()
    for r in range(a.rounds):
        W = witness_lines(ag, S)
        mode = r % 3
        if mode < 2:
            keep = [random.choice(ls) for ls in W if random.random() < a.keep]
            forbid = []
        else:
            keep = []
            forbid = random.sample(S, max(1, len(S) // 8))
        T = solve_round(ag, S, keep, forbid, len(S) - 1, a.time, a.workers)
        msg = f"round {r} mode {mode} keep {len(keep)} forbid {len(forbid)} |S|={len(S)} t={time.time()-t0:.0f}s -> "
        if T is not None:
            ok, miss = is_kakeya(ag, T); assert ok
            S = T
            json.dump({"n": 3, "q": a.q, "points": [list(ag.points[i]) for i in S], "method": "lns", "status": "heuristic"}, open(a.out, "w"))
            msg += f"IMPROVED to {len(S)} (saved {a.out})"
        else:
            msg += "no"
        print(msg, flush=True)
