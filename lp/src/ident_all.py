"""Campaign driver for maxdeg_ident.py: decide 'exists witness system with |K| = n in AG(3,q)' over all (D, orbit rep)
subcases (D from ceil(31*5/n) to floor((n-1)/4)), CP-SAT with the counting identities, 1 worker each, resumable JSONL.
Usage: ident_all.py --q 5 --n 46 [--time 1800] [--part 0/2] [--order asc|desc]
Validation: ident_all.py --q 5 --validate solutions/best_n3_q5_53.json   (model must accept the known set in its own subcase)"""
import argparse, json, math, os, time, sys
from ortools.sat.python import cp_model
from geom import AG, is_kakeya
from orbits_fast import load_or_compute
from maxdeg_ident import build
from maxdeg import describe

def solve(ag, D, R, n, tl, workers=1, hint=None, split=None):
    m, x, y = build(ag, D, R, n, split=split)
    if hint is not None:
        S, W = hint
        for i in range(len(ag.points)): m.AddHint(x[i], int(i in S))
        for l in range(len(ag.lines)): m.AddHint(y[l], int(l in W))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    t0 = time.time(); st = s.Solve(m); dt = time.time() - t0
    res = dict(status=s.StatusName(st), time=round(dt, 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = [i for i in range(len(ag.points)) if s.Value(x[i])]
        ok, miss = is_kakeya(ag, S); assert ok, miss
        res["size"] = len(S); res["points"] = [ag.points[i] for i in S]
    return res

def validate(ag, fn, tl):
    d = json.load(open(fn)); S = {ag.pidx[tuple(p)] for p in d["points"]}
    W = []
    for di in range(len(ag.dirs)):
        W.append(next(l for l in ag.lines_by_dir[di] if ag.lines[l][1] <= S))
    deg = {p: sum(1 for l in W if p in ag.lines[l][1]) for p in S}
    O = max(S, key=lambda p: deg[p]); D = deg[O]
    # translate so O -> origin
    Ov = ag.points[O]
    T = [ag.pidx[ag.sub(ag.points[p], Ov)] for p in range(len(ag.points))]
    S2 = {T[p] for p in S}; W2 = [ag.line_of[frozenset(T[p] for p in ag.lines[l][1])] for l in W]
    R = sorted(ag.lines[l][0] for l in W2 if ag.pidx[(0, 0, 0)] in ag.lines[l][1])
    assert len(R) == D
    res = solve(ag, D, R, len(S), tl, hint=(S2, set(W2)))
    print(json.dumps(dict(validate=fn, n=len(S), D=D, R=[ag.dirs[r] for r in R], status=res["status"], time=res["time"], size=res.get("size"))))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
    ap.add_argument("--time", type=float, default=1800); ap.add_argument("--part", default="0/1"); ap.add_argument("--order", default="asc")
    ap.add_argument("--validate", default=None); ap.add_argument("--tag", default="")
    ap.add_argument("--split", default=None, help="D:rep:e -> 25 sub-subcases forcing each witness line of direction e")
    a = ap.parse_args()
    ag = AG(3, a.q)
    if a.validate:
        validate(ag, a.validate, a.time); sys.exit()
    if a.split:
        D, i, e = map(int, a.split.split(":")); q = a.q; n = a.n
        lv = load_or_compute(q, D, verbose=False); rep, sz, st = lv[D][i]
        out = f"out/ident_q{q}_n{n}_split_D{D}_r{i}_e{e}.jsonl"
        done = set()
        if os.path.exists(out):
            for line in open(out): done.add(json.loads(line)["k"])
        pi, pk = map(int, a.part.split("/"))
        for k in range(len(ag.lines_by_dir[e])):
            if k % pk != pi or k in done: continue
            res = solve(ag, D, rep, n, a.time, split=(e, k))
            rec = dict(q=q, n=n, D=D, rep=i, e=e, k=k, line=sorted(ag.points[p] for p in ag.lines[ag.lines_by_dir[e][k]][1]), **res)
            print(json.dumps({kk: v for kk, v in rec.items() if kk != "points"}), flush=True)
            with open(out, "a") as fo: fo.write(json.dumps(rec) + "\n")
        print("SPLIT_DONE", a.split, a.part); sys.exit()
    q = a.q; nd = len(ag.dirs); n = a.n
    Dmin = math.ceil(nd * q / n); Dmax = (n - 1) // (q - 1)
    lv = load_or_compute(q, Dmax, verbose=False)
    cases = [(D, i, rep, sz) for D in range(Dmin, Dmax + 1) for i, (rep, sz, st) in enumerate(lv[D])]
    if a.order == "desc": cases = cases[::-1]
    pi, pk = map(int, a.part.split("/"))
    out = f"out/ident_q{q}_n{n}{a.tag}.jsonl"
    done = set()
    if os.path.exists(out):
        for line in open(out):
            r = json.loads(line); done.add((r["D"], r["rep"]))
    print(f"q={q} n={n}: D in [{Dmin},{Dmax}], {len(cases)} subcases, {len(done)} done, part {a.part} -> {out}", flush=True)
    for idx, (D, i, rep, sz) in enumerate(cases):
        if idx % pk != pi or (D, i) in done: continue
        res = solve(ag, D, rep, n, a.time)
        rec = dict(q=q, n=n, D=D, rep=i, dirs=[ag.dirs[d] for d in rep], orbit_size=sz, max_coplanar=describe(ag, rep), **res)
        print(json.dumps({k: v for k, v in rec.items() if k != "points"}), flush=True)
        with open(out, "a") as fo: fo.write(json.dumps(rec) + "\n")
    print("PART_DONE", a.part)
