"""Session V: probe structural lemmas with the merged LP (planecount6.build) by adding side constraints.
Usage: lp_probe.py --n 46 --D 4 5 --rmax 4 | --forbid 5,15 2,9 | --min-r5 k | --kmax 14"""
import argparse, math, json, time, sys
import numpy as np
import planecount5 as p5, planecount6 as p6

def probe(q, n, D, rmax=None, forbid=(), max_r5=None, families=(True, True, True, True), cert=False, show=False, extra=None, config=None):
    lp, info, types = p6.build(q, n, D, *families, u_for_r=(config[0] if config else None))
    if config is not None:   # r_max = r and some plane has an inside configuration of r lines with union size u
        r, u = config; rmax = r
        idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]
        lp.add([(lp.col[f"N{i}"], -1) for i in idx], '<', -1, f"config ({r},{u}) present")
    if rmax is not None:
        idx = [i for i, t in enumerate(types) if t[0] > rmax]
        if idx: lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', 0, f"rmax<={rmax}")
    for (r, k) in forbid:
        idx = [i for i, t in enumerate(types) if (t[0], t[1]) == (r, k)]
        if idx: lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', 0, f"forbid ({r},{k})")
    if max_r5 is not None:
        idx = [i for i, t in enumerate(types) if t[0] >= 5]
        lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', max_r5, f"#planes r>=5 <= {max_r5}")
    if extra: extra(lp, types)
    A, b, sense = lp.matrices(); eq = [i for i, s in enumerate(sense) if s == '=']; ub = [i for i, s in enumerate(sense) if s == '<']
    if cert:
        yr, msg, x = p5.farkas_phase1(lp, A, b, sense, eq, ub)
        st = "feasible" if x is not None else ("INFEASIBLE" if yr is not None else "UNDECIDED: " + msg)
        out = dict(status=st, cert_msg=msg)
        if yr is not None: out["cert"] = {lp.rows[i][3]: str(yr[i]) for i in range(len(yr)) if yr[i] != 0}
        class _R: pass
        res = _R(); res.status = 0 if x is not None else 2; res.x = x
    else:
        res, A, b, sense, eq, ub = p5.solve(lp)
        st = {0: "feasible", 2: "INFEASIBLE"}.get(res.status, res.message); out = dict(status=st)
    if res.status == 0 and show:
        sol = p5.show_solution(lp, res.x)
        out["planes"] = [(types[int(k[1:])][:2], types[int(k[1:])][2], v) for k, v in sol.items() if k[0] == 'N']
        out["m"] = {k: v for k, v in sol.items() if k[0] == 'm'}
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, nargs="+", default=[46])
    ap.add_argument("--D", type=int, nargs="+", default=None); ap.add_argument("--rmax", type=int, default=None)
    ap.add_argument("--forbid", nargs="*", default=[]); ap.add_argument("--max-r5", type=int, default=None)
    ap.add_argument("--cert", action="store_true"); ap.add_argument("--show", action="store_true"); ap.add_argument("--tag", default="probe")
    ap.add_argument("--config", default=None, help="r,u")
    for fl in ("classes", "lines", "dir", "hist"): ap.add_argument(f"--no-{fl}", action="store_true")
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    forbid = [tuple(map(int, s.split(","))) for s in a.forbid]
    fam = (not a.no_classes, not a.no_lines, not a.no_dir, not a.no_hist)
    import os
    outfn = f"out/lp_probe_{a.tag}.json"
    res = json.load(open(outfn)) if os.path.exists(outfn) else {}
    for n in a.n:
        Ds = a.D or list(range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1))
        for D in Ds:
            t0 = time.time()
            cfg = tuple(map(int, a.config.split(","))) if a.config else None
            o = probe(q, n, D, a.rmax, forbid, a.max_r5, fam, a.cert, a.show, config=cfg)
            res[f"{n},{D}"] = o; json.dump(res, open(outfn, "w"), indent=1)
            print(f"n={n} D={D} config={a.config} rmax={a.rmax} forbid={forbid} max_r5={a.max_r5}: {o['status']} ({time.time()-t0:.1f}s) {o.get('cert_msg','')}", flush=True)
            if a.show and o["status"] == "feasible": print("   ", o["planes"], o["m"])
    json.dump(res, open(outfn, "w"), indent=1)
