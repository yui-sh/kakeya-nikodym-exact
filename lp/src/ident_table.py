"""Summarise out/ident_q{q}_n{n}.jsonl: per-D counts by status, times, list of non-INFEASIBLE subcases.
Usage: ident_table.py out/ident_q5_n46.jsonl [--md]"""
import json, sys
from collections import defaultdict
fn = sys.argv[1]; md = "--md" in sys.argv
recs = [json.loads(l) for l in open(fn) if l.strip()]
byD = defaultdict(list)
for r in recs: byD[r["D"]].append(r)
tot = len(recs); st = defaultdict(int)
for r in recs: st[r["status"]] += 1
print(f"{fn}: {tot} subcases done; statuses {dict(st)}; total solve time {sum(r['time'] for r in recs)/3600:.2f} h")
if md: print("| D | subcases done | INFEASIBLE | other | min / median / max time (s) |\n|---|---|---|---|---|")
for D in sorted(byD):
    rs = byD[D]; ts = sorted(r["time"] for r in rs); inf = sum(1 for r in rs if r["status"] == "INFEASIBLE")
    other = [(r["rep"], r["status"], r["time"]) for r in rs if r["status"] != "INFEASIBLE"]
    if md: print(f"| {D} | {len(rs)} | {inf} | {len(other)} | {ts[0]} / {ts[len(ts)//2]} / {ts[-1]} |")
    else: print(f"D={D}: {len(rs)} done, {inf} INFEASIBLE, other: {other}; time min/med/max {ts[0]}/{ts[len(ts)//2]}/{ts[-1]}")
