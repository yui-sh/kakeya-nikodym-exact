"""Session V: independent exact check of the Farkas certificates produced by lp_probe.py / planecount6.py.
No LP solver is used: the constraint rows are regenerated (planecount6.build + the probe's side constraints, the same
generator that pc6_validate.py checks against the known Kakeya sets), the stored integer multipliers y are applied
with exact Fractions, and we verify  (a) y_i >= 0 on every '<=' row, (b) sum_i y_i * row_i >= 0 on every column
(all variables are counts, hence >= 0), (c) sum_i y_i * rhs_i < 0.  Together: no nonnegative solution exists.
Usage: cert5_check.py out/lp_probe_rmax4.json [--rmax 4] [--forbid 5,15 ...]"""
import sys, json, argparse
from fractions import Fraction
import planecount6 as p6, planecount5 as p5

ap = argparse.ArgumentParser(); ap.add_argument("file"); ap.add_argument("--q", type=int, default=5)
ap.add_argument("--rmax", type=int, default=None); ap.add_argument("--forbid", nargs="*", default=[]); ap.add_argument("--max-r5", type=int, default=None); ap.add_argument("--config", default=None)
ap.add_argument("--bb", action="store_true", help="file is a degvec_bb.py output: verify every pruned node's certificate and the tree cover")
a = ap.parse_args(); forbid = [tuple(map(int, s.split(","))) for s in a.forbid]
data = json.load(open(a.file)); allok = True
if a.bb:
    import copy
    from degvec_bb import children
    n, D = data["n"], data["D"]; r, u = map(int, data["config"].split(",")); nd = a.q * a.q + a.q + 1; tot = nd * a.q
    lp0, info, types = p6.build(a.q, n, D, u_for_r=r)
    idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]; lp0.add([(lp0.col[f"N{i}"], -1) for i in idx], '<', -1, f"config ({r},{u}) present")
    idx = [i for i, t in enumerate(types) if t[0] > r]
    if idx: lp0.add([(lp0.col[f"N{i}"], 1) for i in idx], '<', 0, f"rmax<={r}")
    mcol = {j: lp0.col[f"m{j}"] for j in range(1, D + 1)}
    pruned = {tuple(p["prefix"]): p for p in data["pruned_nodes"]}; surv = {tuple(s) for s in data["survivors"]}
    nver = 0; nbad = 0; ncover = 0
    def rec(prefix):
        global nver, nbad, ncover
        t = tuple(prefix)
        if t in pruned:
            p = pruned[t]
            if "cert" not in p: nbad += 1; print("  node", prefix, "pruned without certificate:", p.get("status")); return
            lp = copy.copy(lp0); lp.rows = list(lp0.rows)
            for j, v in enumerate(prefix, start=1): lp.add([(mcol[j], 1)], '=', v, f"fix m{j} = {v}")
            names = [row[3] for row in lp.rows]; y = [Fraction(p["cert"].get(nm, "0")) for nm in names]
            ok, rhs = p5.verify(lp, y); extra = set(p["cert"]) - set(names)
            if ok and not extra: nver += 1
            else: nbad += 1; print("  node", prefix, "certificate FAILED")
            return
        if len(prefix) == D:
            if t in surv: ncover += 1
            else: nbad += 1; print("  leaf", prefix, "neither pruned nor listed as survivor")
            return
        for h in children(prefix, n, tot, D): rec(prefix + [h])
    rec([])
    print(f"bb {data['config']} n={n} D={D}: pruned nodes verified={nver}, survivors covered={ncover}/{len(surv)}, problems={nbad} -> {'VERIFIED' if nbad == 0 else 'FAILED'}")
    sys.exit(0 if nbad == 0 else 1)
for key, rec in data.items():
    n, D = map(int, key.split(","))
    if "cert" not in rec: print(key, rec["status"], "- no certificate stored"); allok = allok and rec["status"] != "INFEASIBLE"; continue
    lp, info, types = p6.build(a.q, n, D, u_for_r=(int(a.config.split(",")[0]) if a.config else None))
    rmax = a.rmax
    if a.config:
        r, u = map(int, a.config.split(",")); rmax = r
        idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]
        lp.add([(lp.col[f"N{i}"], -1) for i in idx], '<', -1, f"config ({r},{u}) present")
    if rmax is not None:
        idx = [i for i, t in enumerate(types) if t[0] > rmax]
        if idx: lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', 0, f"rmax<={rmax}")
    for (r, k) in forbid:
        idx = [i for i, t in enumerate(types) if (t[0], t[1]) == (r, k)]
        if idx: lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', 0, f"forbid ({r},{k})")
    if a.max_r5 is not None:
        idx = [i for i, t in enumerate(types) if t[0] >= 5]
        lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', a.max_r5, f"#planes r>=5 <= {a.max_r5}")
    names = [row[3] for row in lp.rows]; assert len(set(names)) == len(names), "row names not unique"
    y = [Fraction(rec["cert"].get(nm, "0")) for nm in names]
    missing = set(rec["cert"]) - set(names)
    ok, rhs = p5.verify(lp, y)
    print(f"n={n} D={D}: rows={len(lp.rows)} vars={len(lp.names)} nonzero multipliers={sum(1 for v in y if v)} missing-rows={len(missing)} y^T b={rhs} -> {'VERIFIED' if ok and not missing else 'FAILED'}", flush=True)
    allok = allok and ok and not missing
print("ALL CERTIFICATES VERIFIED" if allok else "SOME CERTIFICATE FAILED")
