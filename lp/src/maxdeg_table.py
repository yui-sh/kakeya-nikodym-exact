"""Markdown table from out/maxdeg_q*_ub*.jsonl (+ optional CP-SAT json cross-check files)."""
import json, sys, glob, os
fn = sys.argv[1]
rows = [json.loads(l) for l in open(fn)]
cp = {}; split = {}
for f in sorted(set(glob.glob("out/md_D*_r*.json")) | set(glob.glob("out/md_D*_r*_cpsat.json"))):
    d = json.load(open(f))
    key = (d["D"], d["rep"])
    tag = "core" if d.get("core") else "feas"
    if d.get("split"):
        split.setdefault(key, []).append(d)
    else:
        cp.setdefault(key, []).append(f"{d['status']} {d['time']}s ({tag})")
for key, ds in split.items():
    n_inf = sum(1 for d in ds if d["status"] == "INFEASIBLE")
    cp.setdefault(key, []).append(f"split 16: {n_inf}/{len(ds)} INFEASIBLE, max {max(d['time'] for d in ds)}s, total {sum(d['time'] for d in ds):.0f}s (core)")
print("| D | rep | directions through O | orbit size | max coplanar | vars | clauses | kissat | time (s) | DRAT verified (s, MB) | CP-SAT cross-check (planes+proj) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
tot = 0
for r in rows:
    v = f"yes ({r['verify_time']}, {r['drat_bytes']/1e6:.0f})" if r.get("drat_verified") else ("-" if r["status"] != "UNSAT" else "NOT VERIFIED")
    extra = f" size {r['size']} ok={r['kakeya_ok']}" if r["status"] == "SAT" else ""
    dirs = " ".join("".join(map(str, d)) for d in r["dirs"])
    print(f"| {r['D']} | {r['rep']} | {dirs} | {r['orbit_size']} | {r['max_coplanar']} | {r['nvars']} | {r['nclauses']} | {r['status']}{extra} | {r['solve_time']} | {v} | {'; '.join(cp.get((r['D'], r['rep']), ['-']))} |")
    tot += r["solve_time"] + r.get("verify_time", 0)
print(f"\n{len(rows)} subcases, statuses: { {s: sum(1 for r in rows if r['status']==s) for s in set(r['status'] for r in rows)} }, total solve+verify time {tot:.0f} s")
