"""Second-solver cross-check of planecount.py: same aggregated integer model solved by HiGHS MILP (scipy.optimize.milp).
Usage: planecount_milp.py --q 5 --nmin 45 --nmax 46"""
import argparse, math, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from planecount_lp import build
ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--nmin", type=int, default=45); ap.add_argument("--nmax", type=int, default=46)
a = ap.parse_args(); q = a.q; nd = q * q + q + 1
for n in range(a.nmin, a.nmax + 1):
    row = []
    for D in range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1):
        Aeq, beq, eqn, Aub, bub, ubn, names = build(q, n, D)
        cons = [LinearConstraint(Aeq, beq, beq), LinearConstraint(Aub, -np.inf, bub)]
        res = milp(np.zeros(len(names)), constraints=cons, integrality=np.ones(len(names)), bounds=Bounds(0, np.inf), options=dict(time_limit=600))
        st = {0: "ok", 2: "X"}.get(res.status, res.message)
        row.append(f"D={D}:{st}")
    print(f"n={n}: " + "  ".join(row), flush=True)
