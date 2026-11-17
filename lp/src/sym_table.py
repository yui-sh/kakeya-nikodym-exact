"""Summarize sym.py jsonl results into a markdown table (one row per group), plus a per-order summary."""
import json, sys, glob
from collections import defaultdict
rows = []
for fn in sys.argv[1:]:
    for l in open(fn):
        r = json.loads(l); r["file"] = fn.split("/")[-1]; rows.append(r)
def key(r): return (r["q"], r["size"] is None, r["size"] or 10**9, -r["order"])
rows.sort(key=key)
print("| q | group | order | pt-orbits | status | G-inv min (size) | lower bound in G | time s | file |")
print("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['q']} | `{r['name']}` | {r['order']} | {r['npt_orbits']} | {r['status']} | {r['size']} | {r['bound']} | {r['time']} | {r['file']} |")
print()
print("Per order (best size found / #groups / #OPTIMAL / min proven G-invariant value among OPTIMAL):")
agg = defaultdict(list)
for r in rows: agg[(r["q"], r["order"])].append(r)
print("| q | order | #groups | best size found | #OPTIMAL | smallest proven G-min |")
print("|---|---|---|---|---|---|")
for (q, o), rs in sorted(agg.items()):
    sizes = [r["size"] for r in rs if r["size"] is not None]
    opt = [r["size"] for r in rs if r["status"] == "OPTIMAL"]
    print(f"| {q} | {o} | {len(rs)} | {min(sizes) if sizes else '-'} | {len(opt)} | {min(opt) if opt else '-'} |")
