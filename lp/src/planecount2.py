"""Strengthened aggregated model (session IV): planecount.py identities plus
 (A) parallel classes: the 31 classes of 5 parallel planes each have sum r = q+1 = 6 and sum k = n
     (planes of a class partition the points; the 6 witness lines whose directions lie in the class's line at
     infinity lie in exactly one plane of the class each).
 (C) witness lines: the 6 planes through a witness line l have sum k = n + (q-1)*q = n + 20 ... general: 
     sum_{pi ⊃ l} k = q(q+1) + (n - q) = n + q^2, and sum_{pi ⊃ l} r = (q+1) + m_l where m_l = # witness lines meeting l
     = sum_{p in l} (d_p - 1); so sum_{pi ⊃ l} r = 1 + sum_{p in l} d_p, and each plane through l has r >= 1.
     The multiset of degrees on l: 5 values in [1, D] with sum s_l = sum_{pi⊃l} r - 1.
Columns: class types (5-multisets of (r,k) pairs with sum r = 6, sum k = n) and witness-line types
(6-multisets of (r,k) pairs with r >= 1, sum k = n + 25, 6 <= sum r <= 1 + 5D) linked to plane counts N[r,k]:
  N[r,k] = sum_types count_type(r,k) * C_type ;  r * N[r,k] = sum_ltypes count(r,k) * W_ltype ; sum C = 31; sum W = 31.
Also sum over witness lines of (sum r - 1) = sum_p d_p^2 = sum_j j^2 m_j.
Usage: planecount2.py --q 5 --nmin 46 --nmax 50"""
import argparse, math, sys, json, itertools
from ortools.sat.python import cp_model
from plane_f import get

def allowed_pairs(q, D):
    nd = q * q + q + 1; f = get(q)
    return [(r, k) for r in range(q + 2) for k in range(1, q * q + 1) if k >= f[r] and k <= (nd - (q + 1)) + q * r <= D * k]

def multisets(pairs, size, sum_r, sum_k_target, rmin=0, sum_r_range=None):
    """all multisets of `size` pairs (r,k) with r >= rmin, sum k = sum_k_target and sum r == sum_r (or in range)."""
    pairs = [p for p in pairs if p[0] >= rmin]
    out = []
    def rec(start, left, sr, sk, cur):
        if left == 0:
            if sk == sum_k_target and ((sum_r_range is None and sr == sum_r) or (sum_r_range and sum_r_range[0] <= sr <= sum_r_range[1])):
                out.append(tuple(cur))
            return
        for i in range(start, len(pairs)):
            r, k = pairs[i]
            if sk + k > sum_k_target: continue
            if sum_r_range is None and sr + r > sum_r: continue
            if sum_r_range and sr + r > sum_r_range[1]: continue
            # remaining feasibility: min k of the remaining is >= smallest allowed k
            cur.append((r, k)); rec(i, left - 1, sr + r, sk + k, cur); cur.pop()
    rec(0, size, 0, 0, [])
    return out

def feasible(q, n, D, use_classes=True, use_lines=True, time_limit=300, workers=1):
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]
    pairs = allowed_pairs(q, D)
    m = cp_model.CpModel()
    N = {rk: m.NewIntVar(0, npl, f"N{rk}") for rk in pairs}
    m.Add(sum(N.values()) == npl)
    m.Add(sum(r * v for (r, k), v in N.items()) == (q + 1) * nd)
    m.Add(sum(k * v for (r, k), v in N.items()) == nd * n)
    deg = {j: m.NewIntVar(0, n, f"m{j}") for j in range(1, D + 1)}
    m.Add(sum(deg.values()) == n); m.Add(sum(j * v for j, v in deg.items()) == nd * q); m.Add(deg[D] >= 1)
    m.Add(sum(math.comb(r, 2) * v for (r, k), v in N.items()) == sum(math.comb(j, 2) * v for j, v in deg.items()))
    m.Add(sum(k * k * v for (r, k), v in N.items()) == (q + 1) * n * n + (nd - q - 1) * n)
    m.Add(sum(k * r * v for (r, k), v in N.items()) == nd * n + q * nd * q)
    L = {j: m.NewIntVar(0, nd * q * q, f"L{j}") for j in range(1, q + 1)}
    m.Add(sum(L.values()) >= nd * bm)
    m.Add(sum(j * v for j, v in L.items()) == nd * n)
    m.Add(sum(math.comb(j, 2) * v for j, v in L.items()) == math.comb(n, 2))
    m.Add(L[q] >= nd)
    m.Add(sum(math.comb(k, 3) * v for (r, k), v in N.items()) == math.comb(n, 3) + q * sum(math.comb(j, 3) * v for j, v in L.items()))
    info = {}
    if use_classes:
        ctypes = multisets(pairs, q, q + 1, n)
        info["ntypes_class"] = len(ctypes)
        C = [m.NewIntVar(0, nd, f"C{i}") for i in range(len(ctypes))]
        m.Add(sum(C) == nd)
        for rk in pairs:
            m.Add(sum(t.count(rk) * C[i] for i, t in enumerate(ctypes) if rk in t) == N[rk])
    if use_lines:
        ltypes = multisets(pairs, q + 1, None, n + q * q, rmin=1, sum_r_range=(q + 1, 1 + q * D))
        info["ntypes_line"] = len(ltypes)
        W = [m.NewIntVar(0, nd, f"W{i}") for i in range(len(ltypes))]
        m.Add(sum(W) == nd)
        for rk in pairs:
            m.Add(sum(t.count(rk) * W[i] for i, t in enumerate(ltypes) if rk in t) == rk[0] * N[rk])
        # sum over witness lines of (sum_{pi ⊃ l} r - 1) = sum_{l} sum_{p in l} d_p = sum_p d_p^2
        m.Add(sum((sum(r for r, k in t) - 1) * W[i] for i, t in enumerate(ltypes)) == sum(j * j * v for j, v in deg.items()))
        # degrees on a witness line: 5 values in [1,D] summing to s = sum r - 1 -> need 5 <= s <= 5D (already by range)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    st = s.Solve(m)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = dict(N={f"{r},{k}": s.Value(v) for (r, k), v in N.items() if s.Value(v)}, deg={j: s.Value(v) for j, v in deg.items() if s.Value(v)},
                   L={j: s.Value(v) for j, v in L.items()})
        if use_classes: sol["classes"] = [(ctypes[i], s.Value(C[i])) for i in range(len(ctypes)) if s.Value(C[i])]
        if use_lines: sol["wlines"] = [(ltypes[i], s.Value(W[i])) for i in range(len(ltypes)) if s.Value(W[i])]
    return s.StatusName(st), sol, info

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--nmin", type=int, default=46); ap.add_argument("--nmax", type=int, default=50)
    ap.add_argument("--no-classes", action="store_true"); ap.add_argument("--no-lines", action="store_true"); ap.add_argument("--time", type=float, default=300)
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    allsol = {}
    for n in range(a.nmin, a.nmax + 1):
        row = []
        for D in range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1):
            st, sol, info = feasible(q, n, D, not a.no_classes, not a.no_lines, a.time)
            row.append(f"D={D}:{'X' if st == 'INFEASIBLE' else ('ok' if sol else st)}")
            if sol: allsol[f"{n},{D}"] = sol
            print(f"  n={n} D={D}: {st} {info}", flush=True)
        print(f"n={n}: " + "  ".join(row), flush=True)
    json.dump(allsol, open(f"out/planecount2_q{q}.json", "w"), indent=0, default=str)
