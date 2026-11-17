"""Session V: HiGHS MILP (all variables integer) on the merged model for one (config, D, degree vector) survivor.
Verdict is solver-certified only (no exportable proof); used to gauge whether integer branching closes the case."""
import sys, time, json, numpy as np, scipy.sparse as sps
from scipy.optimize import milp, LinearConstraint, Bounds
import planecount5 as p5, planecount6 as p6
q, n = 5, 46; D = int(sys.argv[1]); r, u = map(int, sys.argv[2].split(",")); mvec = [int(v) for v in sys.argv[3].split(",")] if sys.argv[3] != "-" else None; tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
lp, info, types = p6.build(q, n, D, u_for_r=r)
idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]; lp.add([(lp.col[f"N{i}"], -1) for i in idx], '<', -1, "cfg")
idx = [i for i, t in enumerate(types) if t[0] > r]; lp.add([(lp.col[f"N{i}"], 1) for i in idx], '<', 0, "rmax")
for j, v in enumerate(mvec or [], start=1): lp.add([(lp.col[f"m{j}"], 1)], '=', v, f"m{j}")
A, b, sense = lp.matrices(); lo = np.array([b[i] if s == '=' else -np.inf for i, s in enumerate(sense)]); hi = b
t0 = time.time()
res = milp(np.zeros(A.shape[1]), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(A.shape[1]), bounds=Bounds(0, np.inf), options=dict(time_limit=tl, disp=False))
print(json.dumps(dict(D=D, config=sys.argv[2], mvec=mvec, status=int(res.status), message=res.message, time=round(time.time() - t0, 1))))
if res.status == 0 and len(sys.argv) > 5:   # save the integer aggregated solution (clean-negative witness)
    sol = {lp.names[i]: int(round(v)) for i, v in enumerate(res.x) if v > 0.5}
    dec = {}
    for k, v in sol.items():
        i = int(k[1:]) if k[1:].isdigit() else None
        if k[0] == 'N': dec.setdefault("planes", []).append([list(types[i][:2]), list(types[i][2]), types[i][3], v])
        elif k[0] == 'C': dec.setdefault("classes", []).append([info["_ctypes"][i], v])
        elif k[0] == 'A': dec.setdefault("directions", []).append([info["_utypes"][i], v])
        elif k[0] == 'P': dec.setdefault("line_vectors", []).append([info["_wvecs"][i], v])
        elif k[0] == 'W': dec.setdefault("pencils", []).append([[info["_items"][ii][:3] for ii in info["_wt"][i]], v])
        elif k[0] == 'V': dec.setdefault("wline_degrees", []).append([info["_dms"][i], v])
        else: dec.setdefault("scalars", {})[k] = v
    json.dump(dec, open(sys.argv[5], "w"), indent=0, default=str); print("solution ->", sys.argv[5])
