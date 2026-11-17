"""Campaign driver for rmax_model.py: all (r_max, configuration) subcases for |K| = n, resumable JSONL, optional part i/k.
Usage: rmax_all.py --q 5 --n 46 [--time 1800] [--part 0/2] [--D]  (--D: additionally split each configuration by the max degree D)"""
import argparse, json, math, os, sys, time
from ortools.sat.python import cp_model
from geom import AG, is_kakeya
from rmax_model import rmax_lower_bound, plane_configs, build

ap = argparse.ArgumentParser()
ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
ap.add_argument("--time", type=float, default=1800); ap.add_argument("--part", default="0/1"); ap.add_argument("--D", action="store_true")
ap.add_argument("--tag", default="")
a = ap.parse_args(); ag = AG(3, a.q); q = a.q; n = a.n; nd = len(ag.dirs)
rlo, B = rmax_lower_bound(n, q)
cases = []
for r in range(rlo, q + 2):
    for i, cfg in enumerate(plane_configs(ag, r)):
        if a.D:
            for D in range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1): cases.append((r, i, D, cfg))
        else: cases.append((r, i, None, cfg))
out = f"out/rmax_q{q}_n{n}{'_D' if a.D else ''}{a.tag}.jsonl"
done = set()
if os.path.exists(out):
    for line in open(out):
        rec = json.loads(line); done.add((rec["rmax"], rec["idx"], rec["D"]))
pi, pk = map(int, a.part.split("/"))
print(f"q={q} n={n}: r_max >= {rlo} (sum C(d,2) >= {B}); {len(cases)} subcases, {len(done)} done, part {a.part} -> {out}", flush=True)
for k, (r, i, D, cfg) in enumerate(cases):
    if k % pk != pi or (r, i, D) in done: continue
    m, x, y = build(ag, cfg["line_idx"], r, n, D)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = a.time; s.parameters.num_workers = 1
    t0 = time.time(); st = s.Solve(m); dt = time.time() - t0
    rec = dict(q=q, n=n, rmax=r, idx=i, D=D, union=cfg["union"], orbit_size=cfg["orbit_size"], lines=cfg["lines"], status=s.StatusName(st), time=round(dt, 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = [j for j in range(len(ag.points)) if s.Value(x[j])]; ok, miss = is_kakeya(ag, S); assert ok
        rec["size"] = len(S); rec["points"] = [ag.points[j] for j in S]
    print(json.dumps({kk: v for kk, v in rec.items() if kk not in ("points", "lines")}), flush=True)
    with open(out, "a") as fo: fo.write(json.dumps(rec) + "\n")
print("PART_DONE", a.part)
