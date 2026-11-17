"""Session V: certified integer branch-and-bound on the merged LP for one (config, D, degree vector) leaf.
All variables are counts (integers). Depth-first: solve the node LP; if infeasible store an exact Farkas certificate
(planecount5.farkas_phase1) for the node (base rows + branching rows); if feasible and some variable is fractional,
branch on the most fractional one (x <= floor | x >= ceil); if feasible and integral -> an integer aggregated solution
(reported, the case is NOT excluded). Output: out/bbint_{cfg}_D{D}_{mvec}.json with the pruned-node list (paths +
certificates); cert5_check.py --bbint re-verifies certificates and that the pruned nodes form a complete binary tree.
Usage: bb_int.py --D 4 --config 5,15 --mvec 0,4,21,21 [--max-nodes 5000] [--prefer-aggregates]"""
import argparse, json, time, math, sys, copy
import numpy as np
import planecount5 as p5, planecount6 as p6

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
    ap.add_argument("--D", type=int, required=True); ap.add_argument("--config", required=True); ap.add_argument("--mvec", required=True)
    ap.add_argument("--max-nodes", type=int, default=5000); ap.add_argument("--tag", default="")
    a = ap.parse_args(); q, n, D = a.q, a.n, a.D; r, u = map(int, a.config.split(",")); mvec = [int(v) for v in a.mvec.split(",")]
    lp0, info, types = p6.build(q, n, D, u_for_r=r)
    idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]; lp0.add([(lp0.col[f"N{i}"], -1) for i in idx], '<', -1, f"config ({r},{u}) present")
    idx = [i for i, t in enumerate(types) if t[0] > r]
    if idx: lp0.add([(lp0.col[f"N{i}"], 1) for i in idx], '<', 0, f"rmax<={r}")
    for j, v in enumerate(mvec, start=1): lp0.add([(lp0.col[f"m{j}"], 1)], '=', v, f"fix m{j} = {v}")
    nv = len(lp0.names)
    # branching priority: aggregates first is emulated by variable order: m, L, M, then C (classes), A, P, V, W, then N types
    def prio(name):
        return {"m": 0, "L": 1, "M": 2, "C": 3, "A": 4, "P": 5, "V": 6, "W": 7, "N": 8}[name[0]]
    pr = np.array([prio(nm) for nm in lp0.names])
    pruned = []; integral = []; nodes = 0; t0 = time.time(); stop = False
    def rec(path):
        nonlocal nodes, stop
        if stop: return
        nodes += 1
        if nodes > a.max_nodes: stop = True; return
        lp = copy.copy(lp0); lp.rows = list(lp0.rows)
        for (v, s, c) in path:
            if s == "<=": lp.add([(v, 1)], '<', c, f"branch {lp0.names[v]} <= {c}")
            else: lp.add([(v, -1)], '<', -c, f"branch {lp0.names[v]} >= {c}")
        A, b, sense = lp.matrices(); eq = [i for i, s in enumerate(sense) if s == '=']; ub = [i for i, s in enumerate(sense) if s == '<']
        yr, msg, x = p5.farkas_phase1(lp, A, b, sense, eq, ub)
        if x is None:
            if yr is None: pruned.append(dict(path=path, status="UNDECIDED: " + msg)); return
            pruned.append(dict(path=path, cert={lp.rows[i][3]: str(yr[i]) for i in range(len(yr)) if yr[i] != 0})); return
        # the phase-1 x is a feasible point but not necessarily a vertex; re-solve to get a vertex solution for branching
        res, *_ = p5.solve(lp); x = res.x
        frac = np.abs(x - np.round(x)); cand = np.where(frac > 1e-6)[0]
        if len(cand) == 0:
            integral.append(dict(path=path, solution={lp0.names[i]: int(round(x[i])) for i in range(nv) if x[i] > 0.5})); return
        # most fractional within the highest-priority group that has a fractional variable
        g = pr[cand].min(); cand = cand[pr[cand] == g]; v = int(cand[np.argmax(frac[cand])]); c = math.floor(x[v])
        rec(path + [(v, "<=", c)]); rec(path + [(v, ">=", c + 1)])
    rec([])
    und = [p for p in pruned if "cert" not in p]
    out = dict(config=a.config, n=n, D=D, mvec=mvec, nodes=nodes, pruned=len(pruned), undecided=len(und), integral_solutions=integral,
               stopped=stop, pruned_nodes=[dict(path=[(lp0.names[v], s, c) for v, s, c in p["path"]], **{k: p[k] for k in p if k != "path"}) for p in pruned], time=round(time.time() - t0, 1))
    fn = f"out/bbint_{a.config.replace(',', '_')}_D{D}_{a.mvec.replace(',', '_')}{a.tag}.json"
    json.dump(out, open(fn, "w"))
    print(json.dumps({k: v for k, v in out.items() if k not in ("pruned_nodes",)}), flush=True)

if __name__ == "__main__": main()
