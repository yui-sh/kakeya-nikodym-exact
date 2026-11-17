"""LP relaxation of planecount.py's aggregated model, solved with HiGHS (scipy). If the LP is infeasible we extract a
Farkas certificate (dual ray) so the (n, D) exclusion becomes a hand-checkable linear combination of the identities.
Usage: planecount_lp.py --q 5 --nmin 39 --nmax 47"""
import argparse, math, json
import numpy as np
from scipy.optimize import linprog
from plane_f import get

def build(q, n, D):
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]
    pairs = [(r, k) for r in range(q + 2) for k in range(1, q * q + 1) if k >= f[r] and k <= (nd - (q + 1)) + q * r <= D * k]
    degs = list(range(1, D + 1)); Ls = list(range(1, q + 1))
    names = [f"N[{r},{k}]" for r, k in pairs] + [f"m[{j}]" for j in degs] + [f"L[{j}]" for j in Ls]
    nv = len(names)
    iN = {rk: i for i, rk in enumerate(pairs)}; im = {j: len(pairs) + i for i, j in enumerate(degs)}; iL = {j: len(pairs) + len(degs) + i for i, j in enumerate(Ls)}
    Aeq, beq, eqn = [], [], []
    def eq(coefs, rhs, name):
        row = np.zeros(nv)
        for i, c in coefs: row[i] += c
        Aeq.append(row); beq.append(rhs); eqn.append(name)
    eq([(iN[rk], 1) for rk in pairs], npl, "sum_pi 1 = q(q^2+q+1)")
    eq([(iN[(r, k)], r) for r, k in pairs], (q + 1) * nd, "sum_pi r = (q+1)(q^2+q+1)")
    eq([(iN[(r, k)], k) for r, k in pairs], nd * n, "sum_pi k = (q^2+q+1) n")
    eq([(im[j], 1) for j in degs], n, "sum_j m_j = n")
    eq([(im[j], j) for j in degs], nd * q, "sum_j j m_j = q(q^2+q+1)")
    eq([(iN[(r, k)], math.comb(r, 2)) for r, k in pairs] + [(im[j], -math.comb(j, 2)) for j in degs], 0, "sum_pi C(r,2) = sum_p C(d,2)")
    eq([(iN[(r, k)], k * k) for r, k in pairs], (q + 1) * n * n + (nd - q - 1) * n, "sum_pi k^2 = (q+1)n^2 + (q^2-q)n")
    eq([(iN[(r, k)], k * r) for r, k in pairs], nd * n + q * q * nd, "sum_pi k r = (q^2+q+1) n + q^2 (q^2+q+1)")
    eq([(iL[j], j) for j in Ls], nd * n, "sum_j j L_j = (q^2+q+1) n")
    eq([(iL[j], math.comb(j, 2)) for j in Ls], math.comb(n, 2), "sum_j C(j,2) L_j = C(n,2)")
    eq([(iN[(r, k)], math.comb(k, 3)) for r, k in pairs] + [(iL[j], -q * math.comb(j, 3)) for j in Ls], math.comb(n, 3), "sum_pi C(k,3) - q sum_j C(j,3) L_j = C(n,3)")
    # inequalities as A_ub x <= b_ub
    Aub, bub, ubn = [], [], []
    def le(coefs, rhs, name):
        row = np.zeros(nv)
        for i, c in coefs: row[i] += c
        Aub.append(row); bub.append(rhs); ubn.append(name)
    le([(im[D], -1)], -1, "m_D >= 1")
    le([(iL[j], -1) for j in Ls], -nd * bm, "sum_j L_j >= (q^2+q+1) BM(q)")
    le([(iL[q], -1)], -nd, "L_q >= q^2+q+1")
    return np.array(Aeq), np.array(beq, float), eqn, np.array(Aub), np.array(bub, float), ubn, names

def farkas(Aeq, beq, Aub, bub):
    """find y (free), z >= 0 with y^T Aeq + z^T Aub >= 0 (componentwise, since x >= 0) and y^T beq + z^T bub < 0.
    Normalise y^T beq + z^T bub = -1 and minimise sum|y| + sum z for a sparse certificate."""
    me, mu, nv = Aeq.shape[0], Aub.shape[0], Aeq.shape[1]
    # variables: y+ (me), y- (me), z (mu)
    c = np.ones(2 * me + mu)
    A1 = np.hstack([-Aeq.T, Aeq.T, -Aub.T])          # -(y^T Aeq + z^T Aub) <= 0
    b1 = np.zeros(nv)
    A2 = np.hstack([beq, -beq, bub]).reshape(1, -1)  # y^T beq + z^T bub = -1
    res = linprog(c, A_ub=A1, b_ub=b1, A_eq=A2, b_eq=[-1.0], bounds=[(0, None)] * (2 * me + mu), method="highs")
    if res.status != 0: return None
    y = res.x[:me] - res.x[me:2 * me]; z = res.x[2 * me:]
    return y, z

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--nmin", type=int, default=39); ap.add_argument("--nmax", type=int, default=47)
    ap.add_argument("--cert", action="store_true")
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    out = {}
    for n in range(a.nmin, a.nmax + 1):
        row = []
        for D in range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1):
            Aeq, beq, eqn, Aub, bub, ubn, names = build(q, n, D)
            res = linprog(np.zeros(len(names)), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * len(names), method="highs")
            st = {0: "LPfeas", 2: "LPinfeas"}.get(res.status, res.message)
            row.append(f"D={D}:{st}")
            if res.status == 2 and a.cert:
                fk = farkas(Aeq, beq, Aub, bub)
                if fk is not None:
                    y, z = fk
                    out[f"{n},{D}"] = dict(eq={eqn[i]: round(float(y[i]), 6) for i in range(len(eqn)) if abs(y[i]) > 1e-9},
                                           ineq={ubn[i]: round(float(z[i]), 6) for i in range(len(ubn)) if z[i] > 1e-9})
        print(f"n={n}: " + "  ".join(row), flush=True)
    if a.cert:
        json.dump(out, open(f"out/planecount_farkas_q{q}.json", "w"), indent=1)
        print("certificates ->", f"out/planecount_farkas_q{q}.json")
