"""Session VII diagnostic: for an integer aggregated solution (milp_probe JSON), test Lemma D per pencil type using only the
solution's own plane histograms. Usage: f8_diag.py out/milpf7_solution_6_17_D5.json"""
import json, sys, itertools
from collections import Counter
d = json.load(open(sys.argv[1])); D = 5
m = [d['scalars'].get(f'm{j}', 0) for j in range(1, D+1)]
H = {}  # (r,k) -> set of hist (with r>=1)
for (rk, h, u, c) in d['planes']:
    if rk[0] >= 1: H.setdefault(tuple(rk), set()).add(tuple(h))
print('m =', m)
for rk in sorted(H): print(rk, sorted(H[rk]))
dms = [tuple(x[0]) for x in d['wline_degrees']]
for pen, cnt in d['pencils']:
    items = [(r, k) for r, lo, hi in pen for k in [lo] ]
    s = sum(r for r, k in items) - 1
    res = []
    for dm in dms:
        if sum(dm) != s: continue
        dv = [dm.count(j) for j in range(1, D+1)]
        target = tuple(m[j] + 5*dv[j] for j in range(D))
        # each plane must contain dm
        opts = [[h for h in H.get(it, ()) if all(h[j] >= dv[j] for j in range(D))] for it in items]
        if any(not o for o in opts): res.append((dm, 'no-containing-hist')); continue
        # subset-sum over vectors
        reach = {tuple([0]*D)}
        for o in opts:
            reach = {tuple(a+b for a, b in zip(x, h)) for x in reach for h in o}
            reach = {x for x in reach if all(x[j] <= target[j] for j in range(D))}
        res.append((dm, 'OK' if target in reach else f'sum-fail'))
    print(cnt, 'x', items, 's=', s, res)
