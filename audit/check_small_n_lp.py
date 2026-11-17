"""AUDIT: can the counting identities alone (no Bukh-Chao citation) exclude witness systems with n <= 38 in AG(3,5)?
Plain plane-counting LP: N[r,k] >= 0 over admissible (r,k) (k >= f(r), 1 <= k <= 25, 25+5r <= D k for some D <= 11 i.e. k >= ceil((25+5r)/11)),
rows I2a,b,c, I4, I5, plus degree rows I3 with m_j. Infeasibility -> exact Farkas certificate verified with Fractions."""
import math, sys, json
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
f = [0, 5, 9, 12, 14, 15, 17]
def lp(n, D):
    pairs = [(r, k) for r in range(7) for k in range(1, 26) if k >= f[r] and 25 + 5 * r <= D * k]
    rows = []  # (dict var->coef, rhs, sense)
    var = {}
    for p in pairs: var[('N',) + p] = len(var)
    for j in range(1, D + 1): var[('m', j)] = len(var)
    def add(coefs, s, rhs): rows.append(({var[v]: c for v, c in coefs if c}, rhs, s))
    add([(('N', r, k), 1) for r, k in pairs], '=', 155); add([(('N', r, k), r) for r, k in pairs], '=', 186); add([(('N', r, k), k) for r, k in pairs], '=', 31 * n)
    add([(('N', r, k), k * k) for r, k in pairs], '=', 6 * n * n + 25 * n); add([(('N', r, k), k * r) for r, k in pairs], '=', 31 * n + 775)
    add([(('m', j), 1) for j in range(1, D + 1)], '=', n); add([(('m', j), j) for j in range(1, D + 1)], '=', 155); add([(('m', D), -1)], '<', -1)
    add([(('N', r, k), math.comb(r, 2)) for r, k in pairs] + [(('m', j), -math.comb(j, 2)) for j in range(1, D + 1)], '=', 0)
    return var, rows
def farkas(var, rows):
    nv = len(var); eq = [i for i, r in enumerate(rows) if r[2] == '=']; ub = [i for i, r in enumerate(rows) if r[2] == '<']
    A = np.zeros((len(rows), nv)); b = np.zeros(len(rows))
    for i, (d, rhs, s) in enumerate(rows):
        for v, c in d.items(): A[i, v] = c
        b[i] = rhs
    res = linprog(np.zeros(nv), A_ub=A[ub], b_ub=b[ub], A_eq=A[eq], b_eq=b[eq], bounds=(0, None), method='highs')
    if res.status != 2: return None, 'feasible' if res.status == 0 else res.message
    # Farkas: find y with y_ub >= 0, y^T A >= 0, y^T b < 0 (LP in y)
    me, mu = len(eq), len(ub)
    A1 = np.hstack([-A[eq].T, A[eq].T, -A[ub].T]); c = np.ones(2 * me + mu)
    r2 = linprog(c, A_ub=A1, b_ub=np.zeros(nv), A_eq=np.hstack([b[eq], -b[eq], b[ub]]).reshape(1, -1), b_eq=[-1.0], bounds=(0, None), method='highs')
    if r2.status != 0: return None, 'farkas LP failed'
    y = np.zeros(len(rows)); y[eq] = r2.x[:me] - r2.x[me:2 * me]; y[ub] = r2.x[2 * me:]
    def combo(yi):
        col = [Fraction(0)] * nv; rhs = Fraction(0)
        for yv, (d, bb, s) in zip(yi, rows):
            for v, c in d.items(): col[v] += yv * c
            rhs += yv * bb
        return col, rhs
    bump = [Fraction(0)] * len(rows); bump[0] = Fraction(1); bump[5] = Fraction(1)   # rows 'sum N = 155' and 'sum m = n': +1 on every column
    bcol, brhs = combo(bump); assert all(c > 0 for c in bcol) and brhs > 0
    for sc in (10**3, 10**6, 10**9):
        yi = [Fraction(int(round(v * sc))) for v in y]
        for i in ub: yi[i] = max(yi[i], Fraction(0))
        col, rhs = combo(yi)
        T = max([Fraction(0)] + [-col[j] / bcol[j] for j in range(nv) if col[j] < 0]); T = Fraction(math.ceil(T))
        yr = [yi[i] + T * bump[i] for i in range(len(yi))]; col, rhs = combo(yr)
        if all(c >= 0 for c in col) and rhs < 0: return yr, f'exact certificate (scale {sc}, bump {T}), y^T b = {rhs}'
    return None, 'rounding failed'
out = {}; certs = {}
for n in range(20, 46):
    verdicts = {}
    for D in range(math.ceil(155 / n), (n - 1) // 4 + 1):
        var, rows = lp(n, D); y, msg = farkas(var, rows); verdicts[D] = 'INFEASIBLE(cert)' if y else msg
        if y: certs[f'{n},{D}'] = {str(i): str(v) for i, v in enumerate(y) if v}
    out[n] = verdicts; print(n, verdicts)
print('n excluded by the plain counting LP (every D infeasible with an exact certificate):', [n for n, v in out.items() if all(x.startswith('INFEASIBLE') for x in v.values())])
json.dump({'verdicts': out, 'certificates (row index -> multiplier; rows as built by lp())': certs}, open('out/small_n_lp_result.json', 'w'), indent=1)
