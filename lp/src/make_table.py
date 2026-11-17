"""Collect out/*.json into the markdown results table."""
import json, glob, os
from math import comb
from constructions import plane_min, dvir_bound
rows = []
def load(fn):
    return json.load(open(fn)) if os.path.exists(fn) else None
print("| n | q | min size | status | Dvir LB C(q+n-1,n) | Dvir/BM construction | notes |")
print("|---|---|---|---|---|---|---|")
for q in [2, 3, 4, 5, 7, 8, 9]:
    r = load(f"out/sol_n2_q{q}.json")
    sat = load(f"out/sat_n2_q{q}_meet.json")
    st = r["status"] if r else "-"
    if sat and sat.get("unsat_at") == (r["size"] - 1 if r else None):
        st = "OPTIMAL (CP-SAT) " if st == "OPTIMAL" else "CP-SAT feasible only; "
        st += "+ SAT UNSAT@%d" % sat["unsat_at"]
    print(f"| 2 | {q} | {r['size'] if r else '?'} | {st} | {dvir_bound(2,q)} | {plane_min(q)} (Blokhuis-Mazzocca) | matches formula: {r and r['size']==plane_min(q)} |")
for q in [2, 3, 4, 5]:
    r = load(f"out/sol_n3_q{q}.json")
    dv = load(f"out/dvir_n3_q{q}.json")
    if r:
        st = r["status"]; sz = r["size"]; bd = max(r["bound"], dvir_bound(3, q))
        cell = f"{sz}" if st == "OPTIMAL" else f"{bd} <= min <= {sz}"
    else:
        cell = "?"; st = "-"
    print(f"| 3 | {q} | {cell} | {st} | {dvir_bound(3,q)} | {dv['size'] if dv else '-'} (Dvir) | |")
