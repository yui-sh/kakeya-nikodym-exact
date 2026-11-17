"""LP version of planecount3.py: enumerate every degree vector (m_1..m_D) with sum m = n, sum j m_j = 155, m_D >= 1,
and solve the degree-refined plane-counting LP (HiGHS) for each; (n, D) is excluded iff every LP is infeasible.
Usage: planecount3_lp.py --q 5 --n 46 --D 4 5 [--cert]"""
import argparse, math, json, time, sys
import numpy as np
from scipy.optimize import linprog
from plane_f import get
from planecount3 import compositions

def degree_vectors(q, n, D):
    nd = q * q + q + 1; tot = nd * q
    out = []
    def rec(j, left_n, left_t, cur):
        if j == D:
            if left_t == D * left_n and left_n >= 1: out.append(tuple(cur + [left_n]))
            return
        for h in range(left_n + 1):
            rn, rt = left_n - h, left_t - j * h
            if rt < (j + 1) * rn or rt > D * rn: continue
            rec(j + 1, rn, rt, cur + [h])
    rec(1, n, tot, [])
    return out

def types_for(q, D):
    f = get(q); nd = q * q + q + 1; base = nd - q - 1
    types = []
    for r in range(q + 2):
        for k in range(max(f[r], 1), q * q + 1):
            tot = base + q * r
            if tot < k or tot > D * k: continue
            for h in compositions(k, tot, D): types.append((r, k, h))
    return types

def lp(q, n, D, mvec, types, want_cert=False):
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]; base = nd - q - 1
    T = len(types); nL = q
    nv = T + nL
    rows, rhs, names = [], [], []
    def eq(coef_types, coef_L, b, name):
        row = np.zeros(nv); row[:T] = coef_types; row[T:] = coef_L; rows.append(row); rhs.append(b); names.append(name)
    r_ = np.array([t[0] for t in types], float); k_ = np.array([t[1] for t in types], float)
    H = np.array([t[2] for t in types], float)  # T x D
    z = np.zeros(nL); jj = np.arange(1, nL + 1, dtype=float)
    eq(np.ones(T), z, npl, "sum 1"); eq(r_, z, (q + 1) * nd, "sum r"); eq(k_, z, nd * n, "sum k")
    eq(k_ ** 2, z, (q + 1) * n * n + base * n, "sum k^2"); eq(k_ * r_, z, nd * n + q * q * nd, "sum kr")
    eq(np.array([math.comb(int(r), 2) for r in r_]), z, sum(math.comb(j, 2) * mvec[j - 1] for j in range(1, D + 1)), "sum C(r,2) = sum C(d,2)")
    eq(np.zeros(T), jj, nd * n, "sum j L_j"); eq(np.zeros(T), np.array([math.comb(j, 2) for j in range(1, nL + 1)], float), math.comb(n, 2), "sum C(j,2) L_j")
    eq(np.array([math.comb(int(k), 3) for k in k_]), -q * np.array([math.comb(j, 3) for j in range(1, nL + 1)], float), math.comb(n, 3), "sum C(k,3) - q sum C(j,3)L_j")
    for j in range(1, D + 1):
        hj = H[:, j - 1]; mj = mvec[j - 1]
        eq(hj, z, nd * mj, f"sum h{j}"); eq(hj * k_, z, ((q + 1) * n + base) * mj, f"sum h{j} k")
        eq(hj * r_, z, (nd + q * j) * mj, f"sum h{j} r"); eq(hj ** 2, z, (q + 1) * mj * mj + base * mj, f"sum h{j}^2")
        for i in range(1, j):
            eq(H[:, i - 1] * hj, z, (q + 1) * mvec[i - 1] * mj, f"sum h{i} h{j}")
    Aeq = np.array(rows); beq = np.array(rhs, float)
    # inequalities: sum L_j >= nd*bm ; L_q >= nd
    Aub = np.zeros((2, nv)); Aub[0, T:] = -1; Aub[1, T + nL - 1] = -1; bub = np.array([-nd * bm, -nd], float)
    res = linprog(np.zeros(nv), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * nv, method="highs")
    cert = None
    if res.status == 2 and want_cert:
        me, mu = len(rows), 2
        c = np.ones(2 * me + mu)
        A1 = np.hstack([-Aeq.T, Aeq.T, -Aub.T]); A2 = np.hstack([beq, -beq, bub]).reshape(1, -1)
        r2 = linprog(c, A_ub=A1, b_ub=np.zeros(nv), A_eq=A2, b_eq=[-1.0], bounds=[(0, None)] * (2 * me + mu), method="highs")
        if r2.status == 0:
            y = r2.x[:me] - r2.x[me:2 * me]; zz = r2.x[2 * me:]
            cert = dict(eq={names[i]: round(float(y[i]), 6) for i in range(me) if abs(y[i]) > 1e-9}, ineq={["sum L_j >= nd*BM", "L_q >= nd"][i]: round(float(zz[i]), 6) for i in range(mu) if zz[i] > 1e-9})
    return res.status, cert

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, nargs="+", default=[46])
    ap.add_argument("--D", type=int, nargs="+", default=None); ap.add_argument("--cert", action="store_true")
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    summary = {}
    for n in a.n:
        Ds = a.D or list(range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1))
        for D in Ds:
            t0 = time.time(); types = types_for(q, D); mv = degree_vectors(q, n, D)
            feas = []; certs = {}
            for mvec in mv:
                st, cert = lp(q, n, D, mvec, types, a.cert)
                if st == 0: feas.append(mvec)
                elif st == 2 and cert: certs[str(mvec)] = cert
                elif st not in (0, 2): feas.append(("?", mvec, st))
            verdict = "X" if not feas else "ok"
            print(f"n={n} D={D}: {verdict}  ({len(mv)} degree vectors, {len(types)} plane types, {len(feas)} LP-feasible, {time.time()-t0:.1f}s)" + (f" e.g. {feas[:3]}" if feas else ""), flush=True)
            summary[f"{n},{D}"] = dict(verdict=verdict, nvec=len(mv), ntypes=len(types), feasible=[list(x) if isinstance(x, tuple) else str(x) for x in feas[:50]], certs=certs if a.cert else None)
    json.dump(summary, open(f"out/planecount3_lp_q{q}_n{'_'.join(map(str, a.n))}.json", "w"), default=str)
