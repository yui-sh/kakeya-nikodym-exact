"""Session V: for an LP-surviving (r_max configuration, D) case, scan every degree vector (m_1..m_D) with the merged LP
(planecount6.build + config rows + m fixed) and report which vectors survive; certificates are not stored here (the
scan is a splitting device: survivors go to CP-SAT with the degree vector fixed).
Usage: degvec_scan.py --n 46 --D 4 --config 5,15 [--limit N]"""
import argparse, json, time, math, sys
import numpy as np
import planecount5 as p5, planecount6 as p6
from planecount3_lp import degree_vectors

ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
ap.add_argument("--D", type=int, default=4); ap.add_argument("--config", required=True); ap.add_argument("--limit", type=int, default=None)
a = ap.parse_args(); q, n, D = a.q, a.n, a.D; r, u = map(int, a.config.split(","))
lp, info, types = p6.build(q, n, D, u_for_r=r)
idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]; lp.add([(lp.col[f"N{i}"], -1) for i in idx], '<', -1, f"config ({r},{u}) present")
idx = [i for i, t in enumerate(types) if t[0] > r]
if idx: lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', 0, f"rmax<={r}")
base_rows = len(lp.rows)
A0, b0, sense0 = lp.matrices()
vecs = degree_vectors(q, n, D)
if a.limit: vecs = vecs[:a.limit]
surv = []; t0 = time.time()
mcols = [lp.col[f"m{j}"] for j in range(1, D + 1)]
import scipy.sparse as sps
for mv in vecs:
    # append rows m_j = mv_j
    rows = sps.lil_matrix((D, A0.shape[1]))
    for j in range(D): rows[j, mcols[j]] = 1
    A = sps.vstack([A0, rows.tocsr()]).tocsr(); b = np.concatenate([b0, np.array(mv, float)]); sense = sense0 + ['='] * D
    eq = [i for i, s in enumerate(sense) if s == '=']; ub = [i for i, s in enumerate(sense) if s == '<']
    from scipy.optimize import linprog
    res = linprog(np.zeros(A.shape[1]), A_ub=A[ub], b_ub=b[ub], A_eq=A[eq], b_eq=b[eq], bounds=(0, None), method="highs")
    if res.status == 0: surv.append(list(mv))
print(json.dumps(dict(config=a.config, n=n, D=D, nvec=len(vecs), survivors=surv, time=round(time.time() - t0, 1))), flush=True)
json.dump(dict(config=a.config, n=n, D=D, nvec=len(vecs), survivors=surv), open(f"out/degvec_{a.config.replace(',', '_')}_D{D}.json", "w"))
