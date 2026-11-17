"""Degree-refined plane counting (session IV). Plane type = (r, k, h_1..h_D): h_j = # points of degree j in the plane,
sum h = k, sum j h_j = 25 + 5r (identity I1), k >= f(r). Point classes A_j = {p : d_p = j}, |A_j| = m_j.
Identities (all from 'two points lie in 6 common planes, one point in 31, a line not through p spans one plane with p'):
  sum_pi h_j = 31 m_j ;  sum_pi h_j k = m_j (6n+25) ;  sum_pi h_j r = m_j (31 + 5j) ;
  sum_pi h_j^2 = 6 m_j^2 + 25 m_j ;  sum_pi h_i h_j = 6 m_i m_j (i != j) ;
plus planecount.py's (I2)-(I6) and the line-section variables. CP-SAT; products m_i m_j via AddMultiplicationEquality.
Usage: planecount3.py --q 5 --n 46 --D 4 [--D 5 ...] [--time 600]"""
import argparse, math, json, time, itertools, sys
from ortools.sat.python import cp_model
from plane_f import get

def compositions(k, total, D):
    """all (h_1..h_D) with h_j >= 0, sum h_j = k, sum j h_j = total."""
    out = []
    def rec(j, left_k, left_t, cur):
        if j == D:
            if left_t == j * left_k:  # remaining all degree D
                out.append(tuple(cur + [left_k]))
            return
        for h in range(left_k + 1):
            rem_k = left_k - h; rem_t = left_t - j * h
            if rem_t < (j + 1) * rem_k or rem_t > D * rem_k: continue
            rec(j + 1, rem_k, rem_t, cur + [h])
    rec(1, k, total, [])
    return out

def feasible(q, n, D, time_limit=600, workers=1, verbose=False):
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]; base = nd - q - 1
    types = []
    for r in range(q + 2):
        for k in range(max(f[r], 1), q * q + 1):
            tot = base + q * r
            if tot < k or tot > D * k: continue
            for h in compositions(k, tot, D):
                types.append((r, k, h))
    m = cp_model.CpModel()
    N = [m.NewIntVar(0, npl, f"N{i}") for i in range(len(types))]
    deg = {j: m.NewIntVar(0, n, f"m{j}") for j in range(1, D + 1)}
    m.Add(sum(deg.values()) == n); m.Add(sum(j * v for j, v in deg.items()) == nd * q); m.Add(deg[D] >= 1)
    m.Add(sum(N) == npl)
    m.Add(sum(t[0] * v for t, v in zip(types, N)) == (q + 1) * nd)
    m.Add(sum(t[1] * v for t, v in zip(types, N)) == nd * n)
    m.Add(sum(math.comb(t[0], 2) * v for t, v in zip(types, N)) == sum(math.comb(j, 2) * v for j, v in deg.items()))
    m.Add(sum(t[1] ** 2 * v for t, v in zip(types, N)) == (q + 1) * n * n + base * n)
    m.Add(sum(t[1] * t[0] * v for t, v in zip(types, N)) == nd * n + q * q * nd)
    L = {j: m.NewIntVar(0, nd * q * q, f"L{j}") for j in range(1, q + 1)}
    m.Add(sum(L.values()) >= nd * bm); m.Add(sum(j * v for j, v in L.items()) == nd * n)
    m.Add(sum(math.comb(j, 2) * v for j, v in L.items()) == math.comb(n, 2)); m.Add(L[q] >= nd)
    m.Add(sum(math.comb(t[1], 3) * v for t, v in zip(types, N)) == math.comb(n, 3) + q * sum(math.comb(j, 3) * v for j, v in L.items()))
    # degree-class identities
    prod = {}
    for i in range(1, D + 1):
        for j in range(i, D + 1):
            p = m.NewIntVar(0, n * n, f"p{i}{j}"); m.AddMultiplicationEquality(p, [deg[i], deg[j]]); prod[(i, j)] = p
    for j in range(1, D + 1):
        hj = lambda t: t[2][j - 1]
        m.Add(sum(hj(t) * v for t, v in zip(types, N)) == nd * deg[j])
        m.Add(sum(hj(t) * t[1] * v for t, v in zip(types, N)) == (q + 1) * n * deg[j] + base * deg[j])
        m.Add(sum(hj(t) * t[0] * v for t, v in zip(types, N)) == nd * deg[j] + q * j * deg[j])
        m.Add(sum(hj(t) ** 2 * v for t, v in zip(types, N)) == (q + 1) * prod[(j, j)] + base * deg[j])
        for i in range(1, j):
            m.Add(sum(t[2][i - 1] * hj(t) * v for t, v in zip(types, N)) == (q + 1) * prod[(i, j)])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    s.parameters.log_search_progress = verbose
    t0 = time.time(); st = s.Solve(m); dt = time.time() - t0
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = dict(deg={j: s.Value(v) for j, v in deg.items()}, L={j: s.Value(v) for j, v in L.items()},
                   planes=[(t, s.Value(v)) for t, v in zip(types, N) if s.Value(v)])
    return s.StatusName(st), dt, len(types), sol

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, nargs="+", default=[46])
    ap.add_argument("--D", type=int, nargs="+", default=None); ap.add_argument("--time", type=float, default=600); ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    res = {}
    for n in a.n:
        Ds = a.D or list(range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1))
        for D in Ds:
            st, dt, nt, sol = feasible(q, n, D, a.time, verbose=a.verbose)
            print(f"n={n} D={D}: {st} ({nt} plane types, {dt:.1f}s)" + (f"  deg={sol['deg']}" if sol else ""), flush=True)
            res[f"{n},{D}"] = dict(status=st, types=nt, time=round(dt, 1), sol=sol)
    json.dump(res, open(f"out/planecount3_q{q}_n{'_'.join(map(str, a.n))}.json", "w"), default=str)
