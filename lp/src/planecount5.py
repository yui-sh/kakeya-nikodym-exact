"""Session V: unified aggregated LP for witness systems in AG(3,q) with the new per-direction projection family,
exact (rational) Farkas certificates when infeasible, explicit fractional solution when feasible.

Families (all hand-provable; WRITEUP5 §2):
 N[r,k]  plane types (session IV identities I1-I6, planar table f)
 C[t]    parallel-class types: 5-multisets of (r,k), sum r = q+1, sum k = n   (link N = sum cnt*C, sum C = nd)
 W[t]    witness-line types: 6-multisets of (r,k), r >= 1, sum k = n + q^2, sum r = 1 + s, q <= s <= qD
         (link r N[r,k] = sum cnt*W, sum W = nd, sum (s) W = sum j^2 m_j)
 A[u]    direction types: sorted (q+1)-tuple u of a_{d,L} = # planes of class L through d containing a witness line
         other than l_d (1..q); linked to classes: sum_u cnt_u(v) A[u] = sum_t phi_v(t) C[t]
 P[w]    per-direction line vectors w = (L_1..L_q), sum j L_j = n, L_q >= 1, lambda(w) = sum L_j = |projection|
         dominance link: for each theta, sum_{w: lambda >= theta} P[w] >= sum_{u: Ft(u) >= theta} A[u]
         and global L_j = sum_w L_j(w) P[w] in the 3rd-moment identity.
Usage: planecount5.py --q 5 --n 46 [--D 4 5 ...] [--no-classes] [--no-lines] [--no-dir] [--cert]"""
import argparse, math, json, sys, time, itertools
from fractions import Fraction
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog
from plane_f import get
from union_tables import tables

def allowed_pairs(q, D, n):
    nd = q * q + q + 1; f = get(q)
    return [(r, k) for r in range(q + 2) for k in range(1, min(q * q, n) + 1) if k >= f[r] and k <= (nd - (q + 1)) + q * r <= D * k]

def multisets(pairs, size, sum_k_target, sum_r_exact=None, sum_r_range=None, rmin=0):
    pairs = [p for p in pairs if p[0] >= rmin]
    out = []
    def rec(start, left, sr, sk, cur):
        if left == 0:
            if sk == sum_k_target and ((sum_r_exact is not None and sr == sum_r_exact) or (sum_r_range is not None and sum_r_range[0] <= sr <= sum_r_range[1])):
                out.append(tuple(cur))
            return
        for i in range(start, len(pairs)):
            r, k = pairs[i]
            if sk + k > sum_k_target: continue
            if sum_r_exact is not None and sr + r > sum_r_exact: continue
            if sum_r_range is not None and sr + r > sum_r_range[1]: continue
            cur.append((r, k)); rec(i, left - 1, sr + r, sk + k, cur); cur.pop()
    rec(0, size, 0, 0, [])
    return out

def line_vectors(q, n):
    """all (L_1..L_q) with sum j L_j = n, L_q >= 1, sum L_j <= q^2"""
    out = []
    def rec(j, left, cur):
        if j == 0:
            if left == 0 and sum(cur) <= q * q: out.append(tuple(reversed(cur)))
            return
        lo = 1 if j == q else 0
        for c in range(lo, left // j + 1): rec(j - 1, left - j * c, cur + [c])
    rec(q, n, [])
    return out

class LP:
    def __init__(self):
        self.names = []; self.rows = []  # rows: (dict col->coef, rhs, sense, name); sense in '=', '<'
        self.col = {}
    def var(self, name):
        if name not in self.col: self.col[name] = len(self.names); self.names.append(name)
        return self.col[name]
    def add(self, coefs, sense, rhs, name):
        d = {}
        for v, c in coefs:
            if c: d[v] = d.get(v, 0) + c
        d = {v: c for v, c in d.items() if c}
        self.rows.append((d, rhs, sense, name))
    def matrices(self):
        nv = len(self.names); A = sps.lil_matrix((len(self.rows), nv)); b = np.zeros(len(self.rows)); sense = []
        for i, (d, rhs, s, nm) in enumerate(self.rows):
            for v, c in d.items(): A[i, v] = c
            b[i] = rhs; sense.append(s)
        return A.tocsr(), b, sense

def build(q, n, D, use_classes=True, use_lines=True, use_dir=True, lines_third=True):
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]; base = nd - q - 1
    pairs = allowed_pairs(q, D, n)
    lp = LP(); info = {}
    N = {rk: lp.var(f"N{rk}") for rk in pairs}
    m = {j: lp.var(f"m{j}") for j in range(1, D + 1)}
    L = {j: lp.var(f"L{j}") for j in range(1, q + 1)}
    lp.add([(N[rk], 1) for rk in pairs], '=', npl, "I2a: sum_pi 1")
    lp.add([(N[(r, k)], r) for r, k in pairs], '=', (q + 1) * nd, "I2b: sum_pi r")
    lp.add([(N[(r, k)], k) for r, k in pairs], '=', nd * n, "I2c: sum_pi k")
    lp.add([(m[j], 1) for j in m], '=', n, "deg: sum m_j")
    lp.add([(m[j], j) for j in m], '=', nd * q, "deg: sum j m_j")
    lp.add([(m[D], -1)], '<', -1, "m_D >= 1")
    lp.add([(N[(r, k)], math.comb(r, 2)) for r, k in pairs] + [(m[j], -math.comb(j, 2)) for j in m], '=', 0, "I3: sum C(r,2) = sum C(d,2)")
    lp.add([(N[(r, k)], k * k) for r, k in pairs], '=', (q + 1) * n * n + base * n, "I4: sum k^2")
    lp.add([(N[(r, k)], k * r) for r, k in pairs], '=', nd * n + q * q * nd, "I5: sum k r")
    lp.add([(L[j], j) for j in L], '=', nd * n, "I6a: sum j L_j")
    lp.add([(L[j], math.comb(j, 2)) for j in L], '=', math.comb(n, 2), "I6b: sum C(j,2) L_j")
    lp.add([(L[j], -1) for j in L], '<', -nd * bm, "I6c: sum L_j >= nd*BM")
    lp.add([(L[q], -1)], '<', -nd, "I6d: L_q >= nd")
    if lines_third:
        lp.add([(N[(r, k)], math.comb(k, 3)) for r, k in pairs] + [(L[j], -q * math.comb(j, 3)) for j in L], '=', math.comb(n, 3), "I6e: sum C(k,3) - q sum C(j,3) L_j")
    ctypes = []
    if use_classes:
        ctypes = multisets(pairs, q, n, sum_r_exact=q + 1); info["class_types"] = len(ctypes)
        C = [lp.var(f"C{i}") for i in range(len(ctypes))]
        lp.add([(c, 1) for c in C], '=', nd, "classes: sum C = nd")
        for rk in pairs:
            lp.add([(C[i], t.count(rk)) for i, t in enumerate(ctypes) if rk in t] + [(N[rk], -1)], '=', 0, f"class link N{rk}")
    if use_lines:
        ltypes = multisets(pairs, q + 1, n + q * q, sum_r_range=(q + 1, 1 + q * D), rmin=1); info["line_types"] = len(ltypes)
        W = [lp.var(f"W{i}") for i in range(len(ltypes))]
        lp.add([(w, 1) for w in W], '=', nd, "wlines: sum W = nd")
        for rk in pairs:
            lp.add([(W[i], t.count(rk)) for i, t in enumerate(ltypes) if rk in t] + [(N[rk], -rk[0])], '=', 0, f"wline link r*N{rk}")
        lp.add([(W[i], sum(r for r, k in t) - 1) for i, t in enumerate(ltypes)] + [(m[j], -j * j) for j in m], '=', 0, "wlines: sum s_l = sum d_p^2")
    if use_dir:
        assert use_classes
        F, Ft = tables(q)
        utypes = sorted(Ft); info["dir_types"] = len(utypes)
        A = [lp.var(f"A{i}") for i in range(len(utypes))]
        lp.add([(a, 1) for a in A], '=', nd, "dirs: sum A = nd")
        def phi(t, v):  # number of directions d in the class's line at infinity with a_{d,L} = v
            rho = sum(1 for r, k in t if r >= 1); n1 = sum(1 for r, k in t if r == 1)
            return (q + 1 - n1) * (rho == v) + n1 * (rho - 1 == v)
        for v in range(1, q + 1):
            lp.add([(A[i], u.count(v)) for i, u in enumerate(utypes)] + [(C[i], -phi(t, v)) for i, t in enumerate(ctypes) if phi(t, v)], '=', 0, f"dir flags a=={v}")
        wvecs = line_vectors(q, n); info["line_vectors"] = len(wvecs)
        P = [lp.var(f"P{i}") for i in range(len(wvecs))]
        lp.add([(p, 1) for p in P], '=', nd, "dirs: sum P = nd")
        for j in range(1, q + 1):
            lp.add([(P[i], w[j - 1]) for i, w in enumerate(wvecs) if w[j - 1]] + [(L[j], -1)], '=', 0, f"dir link L{j}")
        lam = [sum(w) for w in wvecs]
        for theta in range(bm + 1, q * q + 1):
            lp.add([(P[i], -1) for i in range(len(wvecs)) if lam[i] >= theta] + [(A[i], 1) for i, u in enumerate(utypes) if Ft[u] >= theta], '<', 0, f"dominance lambda>={theta}")
    return lp, info

def solve(lp):
    A, b, sense = lp.matrices()
    eq = [i for i, s in enumerate(sense) if s == '=']; ub = [i for i, s in enumerate(sense) if s == '<']
    res = linprog(np.zeros(A.shape[1]), A_ub=A[ub] if ub else None, b_ub=b[ub] if ub else None, A_eq=A[eq] if eq else None, b_eq=b[eq] if eq else None,
                  bounds=(0, None), method="highs")
    return res, A, b, sense, eq, ub

def farkas_exact(lp, A, b, sense, eq, ub, scale=10**6):
    """numerical Farkas via HiGHS, then an exactly verified integer certificate:
    y_int = round(scale*y); then add T * (bump) where bump is an integer combination of equality rows that is
    strictly positive on every column (the 'sum of group' rows), T minimal so that y^T A >= 0 exactly; finally check y^T b < 0."""
    me, mu, nv = len(eq), len(ub), A.shape[1]
    Aeq, Aub = A[eq], A[ub]
    c = np.ones(2 * me + mu)
    A1 = sps.hstack([-Aeq.T, Aeq.T, -Aub.T]).tocsr()
    A2 = np.hstack([b[eq], -b[eq], b[ub]]).reshape(1, -1)
    res = linprog(c, A_ub=A1, b_ub=np.zeros(nv), A_eq=A2, b_eq=[-1.0], bounds=(0, None), method="highs")
    if res.status != 0: return None, "farkas LP failed: " + res.message
    y = np.zeros(len(sense)); y[eq] = res.x[:me] - res.x[me:2 * me]; y[ub] = res.x[2 * me:]
    bump = bump_rows(lp)
    if bump is None: return None, "no positive bump combination"
    for sc in (scale, scale * 100, scale * 10**4):
        yi = [Fraction(int(round(v * sc))) for v in y]
        for i in ub:
            if yi[i] < 0: yi[i] = Fraction(0)
        col, rhs = combo(lp, yi)
        bcol, brhs = combo(lp, bump)
        assert all(v > 0 for v in bcol) and brhs > 0
        T = max([Fraction(0)] + [-col[j] / bcol[j] for j in range(nv) if col[j] < 0])
        T = Fraction(math.ceil(T))
        yr = [yi[i] + T * bump[i] for i in range(len(yi))]
        ok, r = verify(lp, yr)
        if ok:
            g = 0
            for v in yr: g = math.gcd(g, int(v))
            yr = [v / g for v in yr]
            return yr, f"exact integer certificate (scale {sc}, bump T={T}, gcd {g}); y^T b = {r / g}"
    return None, "rationalisation failed"

def farkas_phase1(lp, A, b, sense, eq, ub, scale=10**6):
    """one LP: minimise the total artificial slack; optimum 0 <=> feasible. If > 0 the dual marginals are a Farkas
    vector (up to sign); it is then rounded to integers and repaired exactly as in farkas_exact."""
    me, mu, nv = len(eq), len(ub), A.shape[1]
    Aeq, Aub = A[eq], A[ub]
    I = sps.identity(me, format="csr"); Iu = sps.identity(mu, format="csr")
    Ae = sps.hstack([Aeq, I, -I, sps.csr_matrix((me, mu))]).tocsr() if me else None
    Au = sps.hstack([Aub, sps.csr_matrix((mu, 2 * me)), -Iu]).tocsr() if mu else None
    c = np.concatenate([np.zeros(nv), np.ones(2 * me + mu)])
    res = linprog(c, A_ub=Au, b_ub=b[ub] if mu else None, A_eq=Ae, b_eq=b[eq] if me else None, bounds=(0, None), method="highs")
    if res.status != 0: return None, "phase-1 LP failed: " + res.message, None
    if res.fun < 1e-7: return None, "feasible (phase-1 optimum 0)", res.x[:nv]
    y = np.zeros(len(sense)); y[eq] = res.eqlin.marginals; y[ub] = res.ineqlin.marginals
    bump = bump_rows(lp)
    for sgn in (-1, 1):
        ys = sgn * y
        if any(ys[i] < -1e-9 for i in ub): continue
        for sc in (scale, scale * 100, scale * 10**4):
            yi = [Fraction(int(round(v * sc))) for v in ys]
            for i in ub:
                if yi[i] < 0: yi[i] = Fraction(0)
            col, rhs = combo(lp, yi); bcol, brhs = combo(lp, bump)
            T = max([Fraction(0)] + [-col[j] / bcol[j] for j in range(nv) if col[j] < 0]); T = Fraction(math.ceil(T))
            yr = [yi[i] + T * bump[i] for i in range(len(yi))]
            ok, r = verify(lp, yr)
            if ok:
                g = 0
                for v in yr: g = math.gcd(g, int(v))
                yr = [v / g for v in yr]
                return yr, f"exact integer certificate from phase-1 dual (scale {sc}, bump T={T}, gcd {g}); y^T b = {r / g}", None
    return None, "phase-1 dual rationalisation failed", None

def combo(lp, y):
    nv = len(lp.names); col = [Fraction(0)] * nv; rhs = Fraction(0)
    for yi, (d, bi, s, nm) in zip(y, lp.rows):
        if yi == 0: continue
        for v, c in d.items(): col[v] += yi * c
        rhs += yi * bi
    return col, rhs

def bump_rows(lp):
    """integer multipliers on equality rows whose combination is > 0 on every column: start with all 'sum' rows
    (coefficient pattern >= 0), then repair the remaining columns with known equality rows: degree pairs M via the
    H6 rows (+ the degree-sum row), session-VI flags Y via the F6 link rows (negated), X via the F6 R1 rows, and the
    plane-type columns N (made negative by the R1 rows) via the plane-sum row I2a. Rows are equalities, so any
    integer multipliers are allowed; only strict positivity on every column matters."""
    nv = len(lp.names); y = [Fraction(0)] * len(lp.rows)
    for i, (d, bi, s, nm) in enumerate(lp.rows):
        if s == '=' and all(c > 0 for c in d.values()): y[i] = Fraction(1)
    def rows_with(prefix): return [i for i, r in enumerate(lp.rows) if r[3].startswith(prefix)]
    def worst(prefix):
        col, rhs = combo(lp, y)
        vals = [col[j] for j in range(nv) if lp.names[j].startswith(prefix)]
        return min(vals) if vals else Fraction(1)
    h6, degrow, link, r1, i2a = rows_with("H6"), rows_with("deg: sum m_j"), rows_with("F6 link"), rows_with("F6 R1"), rows_with("I2a")
    for _ in range(4):
        col, rhs = combo(lp, y)
        if all(c > 0 for c in col): return y
        if worst("M") <= 0 and h6 and degrow:
            for i in h6: y[i] += 1
            w = worst("m")
            if w <= 0: y[degrow[0]] += Fraction(1) - w
        if worst("Y") <= 0 and link:
            w = worst("Y")
            for i in link: y[i] -= Fraction(1) - w
        if worst("X") <= 0 and r1:
            w = worst("X")
            for i in r1: y[i] += Fraction(1) - w
        if worst("N") <= 0 and i2a:
            w = worst("N")
            y[i2a[0]] += Fraction(1) - w
    col, rhs = combo(lp, y)
    return y if all(c > 0 for c in col) else None

def verify(lp, yr):
    """exact check with Fractions on lp.rows. returns (ok, y^T b)."""
    nv = len(lp.names); col = [Fraction(0)] * nv; rhs = Fraction(0)
    for yi, (d, bi, s, nm) in zip(yr, lp.rows):
        if yi == 0: continue
        if s == '<' and yi < 0: return False, None
        for v, c in d.items(): col[v] += yi * c
        rhs += yi * bi
    return (all(c >= 0 for c in col) and rhs < 0), rhs

def show_solution(lp, x, tol=1e-7):
    out = {}
    for i, v in enumerate(x):
        if v > tol: out[lp.names[i]] = round(float(v), 6)
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, nargs="+", default=[46])
    ap.add_argument("--D", type=int, nargs="+", default=None); ap.add_argument("--no-classes", action="store_true"); ap.add_argument("--no-lines", action="store_true")
    ap.add_argument("--no-dir", action="store_true"); ap.add_argument("--cert", action="store_true"); ap.add_argument("--show", action="store_true"); ap.add_argument("--tag", default="")
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    results = {}
    for n in a.n:
        Ds = a.D or list(range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1))
        row = []
        for D in Ds:
            t0 = time.time()
            lp, info = build(q, n, D, not a.no_classes, not a.no_lines, not a.no_dir)
            res, A, b, sense, eq, ub = solve(lp)
            st = {0: "feasible", 2: "INFEASIBLE"}.get(res.status, res.message)
            rec = dict(status=st, info=info, nvars=len(lp.names), nrows=len(lp.rows), time=round(time.time() - t0, 1))
            if res.status == 2 and a.cert:
                yr, msg = farkas_exact(lp, A, b, sense, eq, ub); rec["cert_msg"] = msg
                if yr is not None:
                    rec["cert"] = {lp.rows[i][3]: str(yr[i]) for i in range(len(yr)) if yr[i] != 0}
            if res.status == 0 and a.show: rec["solution"] = show_solution(lp, res.x)
            results[f"{n},{D}"] = rec
            print(f"n={n} D={D}: {st}  vars={rec['nvars']} rows={rec['nrows']} {info} {rec['time']}s" + (f"  cert: {rec.get('cert_msg','')}" if a.cert else ""), flush=True)
            if a.show and res.status == 0:
                sol = rec["solution"]
                print("   N:", {k: v for k, v in sol.items() if k.startswith("N")}); print("   m:", {k: v for k, v in sol.items() if k.startswith("m")}, " L:", {k: v for k, v in sol.items() if k.startswith("L")})
                print("   A:", {k: v for k, v in sol.items() if k.startswith("A")}); print("   P:", {k: v for k, v in sol.items() if k.startswith("P")})
            row.append(f"D={D}:{'X' if st == 'INFEASIBLE' else 'ok'}")
        print(f"n={n}: " + "  ".join(row), flush=True)
    tag = a.tag or ("full" if not (a.no_classes or a.no_lines or a.no_dir) else "".join(s for s, fl in (("c", a.no_classes), ("l", a.no_lines), ("d", a.no_dir)) if not fl))
    json.dump(results, open(f"out/planecount5_q{q}_{tag}.json", "w"), indent=1)
