"""Aggregated plane-counting feasibility model (session IV, hand-provable identities; see WRITEUP4 §2).
K = union of one witness line per direction in AG(3,q), |K| = n, degrees d_p (# witness lines through p),
D = max degree. For a plane pi: r = # witness lines inside, k = |K ∩ pi|. Identities:
 (I1) sum_{p in pi} d_p = (q^2+q+1 - (q+1)) + q r   [lines not in pi and not parallel to pi meet it once]
 (I2) sum_pi 1 = q(q^2+q+1); sum_pi r = (q+1)(q^2+q+1); sum_pi k = (q^2+q+1) n
 (I3) sum_pi C(r,2) = sum_p C(d_p,2)   [two witness lines meet iff coplanar; distinct directions]
 (I4) sum_pi k^2 = (q+1) n^2 + (q^2 - q) n  ... general: sum_{pi∋p} k = (q^2+q+1) + (n-1)(q+1)
 (I5) sum_pi k r = (q^2+q+1) n + q*(q+1)(q^2+q+1)/... general: sum_{pi∋p} r = (q^2+q+1) + q d_p
 (I6) line sections: sum_j L_j = Lambda >= (q^2+q+1)*BM(q) (projection), sum_j j L_j = (q^2+q+1) n,
      sum_j C(j,2) L_j = C(n,2), L_q >= q^2+q+1; sum_pi C(k,3) = C(n,3) + q sum_j C(j,3) L_j
 plus per plane: f(q,r) <= k (planar table), k <= sum d <= D k, and k <= q^2.
Variables: counts N[r,k] of planes, m[j] of points of degree j, L[j] of lines with j points of K.
Usage: planecount.py --q 5 --nmin 39 --nmax 52 [--planes-through-O]"""
import argparse, json, math, sys
from ortools.sat.python import cp_model
from plane_f import get

def feasible(q, n, D, extra=True, time_limit=60, verbose=False):
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]
    m = cp_model.CpModel()
    pairs = [(r, k) for r in range(q + 2) for k in range(1, q * q + 1)
             if k >= f[r] and k <= (nd - (q + 1)) + q * r <= D * k]
    N = {rk: m.NewIntVar(0, npl, f"N{rk}") for rk in pairs}
    m.Add(sum(N.values()) == npl)
    m.Add(sum(r * v for (r, k), v in N.items()) == (q + 1) * nd)
    m.Add(sum(k * v for (r, k), v in N.items()) == nd * n)
    deg = {j: m.NewIntVar(0, n, f"m{j}") for j in range(1, D + 1)}
    m.Add(sum(deg.values()) == n); m.Add(sum(j * v for j, v in deg.items()) == nd * q); m.Add(deg[D] >= 1)
    m.Add(sum(math.comb(r, 2) * v for (r, k), v in N.items()) == sum(math.comb(j, 2) * v for j, v in deg.items()))
    if extra:
        m.Add(sum(k * k * v for (r, k), v in N.items()) == (q + 1) * n * n + (nd - q - 1) * n)
        m.Add(sum(k * r * v for (r, k), v in N.items()) == nd * n + q * nd * q)   # sum_p (nd + q d_p) = nd n + q*(q nd)
        L = {j: m.NewIntVar(0, nd * q * q, f"L{j}") for j in range(1, q + 1)}
        m.Add(sum(L.values()) >= nd * bm)
        m.Add(sum(j * v for j, v in L.items()) == nd * n)
        m.Add(sum(math.comb(j, 2) * v for j, v in L.items()) == math.comb(n, 2))
        m.Add(L[q] >= nd)
        m.Add(sum(math.comb(k, 3) * v for (r, k), v in N.items()) == math.comb(n, 3) + q * sum(math.comb(j, 3) * v for j, v in L.items()))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = 1
    st = s.Solve(m)
    res = s.StatusName(st)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = dict(N={f"{r},{k}": s.Value(v) for (r, k), v in N.items() if s.Value(v)}, deg={j: s.Value(v) for j, v in deg.items() if s.Value(v)})
    return res, sol

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5); ap.add_argument("--nmin", type=int, default=39); ap.add_argument("--nmax", type=int, default=52)
    ap.add_argument("--basic", action="store_true", help="only I1-I3 + planar table")
    a = ap.parse_args()
    q = a.q; nd = q * q + q + 1
    table = {}
    for n in range(a.nmin, a.nmax + 1):
        Dmin = math.ceil(nd * q / n); Dmax = (n - 1) // (q - 1)
        row = []
        for D in range(Dmin, Dmax + 1):
            st, sol = feasible(q, n, D, extra=not a.basic)
            row.append((D, st))
            if st not in ("INFEASIBLE",) and sol:
                table[f"{n},{D}"] = sol
        print(f"n={n}: " + "  ".join(f"D={D}:{'X' if st=='INFEASIBLE' else ('ok' if st in ('OPTIMAL','FEASIBLE') else st)}" for D, st in row), flush=True)
    json.dump(table, open(f"out/planecount_q{q}.json", "w"), indent=0)
