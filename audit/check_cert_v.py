"""AUDIT: independent replay of the 76 session-V reduction certificates (r_max<=4 for D=4..11; the ten (r,u) configuration
probes for D=4..11) with a FROM-SCRATCH row generator for the session-V merged LP (WRITEUP5 §2-3; no sections/F6/F7).
Only the row NAMES are taken from the research code (the saved certificates are keyed by row name); every coefficient,
every column universe and the exact Farkas check are re-derived here. Needs out/planar_raw.txt (fresh Ft table)."""
import json, math, itertools, sys
from fractions import Fraction
from collections import defaultdict
C2 = math.comb
q, n = 5, 46; nd = 31; npl = 155; base = 25; f = [0, 5, 9, 12, 14, 15, 17]; bm = 17
# fresh tables
prof = {}; inside = defaultdict(set)
for ln in open('out/planar_raw.txt'):
    t = ln.split()
    if t[0] == 'PROF': prof[tuple(map(int, t[1:7]))] = int(t[7])
    elif t[0] == 'SEC':
        a = list(map(int, ln[4:].split('|')[0].split())); k, r, u = a[:3]; inside[r].add((u,) + tuple(a[3:3 + r]))
Ft = {a: min(v for p, v in prof.items() if all(x >= y for x, y in zip(p, a))) for a in itertools.combinations_with_replacement(range(1, 6), 6)}
def types_for(D):
    out = []
    for r in range(7):
        S = base + q * r
        for k in range(max(f[r], 1), 26):
            if S < k or S > D * k: continue
            profs = [p for p in inside[r] if p[0] <= k]
            def rec(j, lk, lS, cur):
                if j == D:
                    if lS == D * lk: yield tuple(cur + [lk])
                    return
                for h in range(lk + 1):
                    rk, rS = lk - h, lS - j * h
                    if rS < (j + 1) * rk or rS > D * rk: continue
                    yield from rec(j + 1, rk, rS, cur + [h])
            for h in rec(1, k, S, []):
                tail = [sum(h[j - 1:]) for j in range(1, D + 2)]
                us = {p[0] for p in profs if all(sum(p[1:][i - 1:]) <= (tail[i - 1] if i <= D else 0) for i in range(1, r + 1))}
                for u in sorted(us): out.append((r, k, h, u))
    return out
def build(D, rmax=None, config=None):
    rows = []   # (name, sense, rhs, {col: coef})
    def add(name, sense, rhs, coefs):
        d = defaultdict(int)
        for c, v in coefs: d[c] += v
        rows.append((name, sense, rhs, {c: v for c, v in d.items() if v}))
    T = types_for(D); N = [('N',) + t for t in T]
    m = {j: ('m', j) for j in range(1, D + 1)}; L = {j: ('L', j) for j in range(1, 6)}
    M = {(i, j): ('M', i, j) for i in range(1, D + 1) for j in range(i, D + 1)}
    add("I2a: sum_pi 1", '=', npl, [(c, 1) for c in N]); add("I2b: sum_pi r", '=', 6 * nd, [(c, c[1]) for c in N]); add("I2c: sum_pi k", '=', nd * n, [(c, c[2]) for c in N])
    add("deg: sum m_j", '=', n, [(m[j], 1) for j in m]); add("deg: sum j m_j", '=', nd * q, [(m[j], j) for j in m]); add("m_D >= 1", '<', -1, [(m[D], -1)])
    add("I3: sum C(r,2) = sum C(d,2)", '=', 0, [(c, C2(c[1], 2)) for c in N] + [(m[j], -C2(j, 2)) for j in m])
    add("I4: sum k^2", '=', 6 * n * n + 25 * n, [(c, c[2] ** 2) for c in N]); add("I5: sum k r", '=', nd * n + 25 * nd, [(c, c[1] * c[2]) for c in N])
    add("I6a: sum j L_j", '=', nd * n, [(L[j], j) for j in L]); add("I6b: sum C(j,2) L_j", '=', C2(n, 2), [(L[j], C2(j, 2)) for j in L])
    add("I6c: sum L_j >= nd*BM", '<', -nd * bm, [(L[j], -1) for j in L]); add("I6d: L_q >= nd", '<', -nd, [(L[5], -1)])
    add("I6e: sum C(k,3) - q sum C(j,3) L_j", '=', C2(n, 3), [(c, C2(c[2], 3)) for c in N] + [(L[j], -q * C2(j, 3)) for j in L])
    for j in range(1, D + 1):
        hj = lambda c: c[3][j - 1]
        add(f"H1: sum_pi h{j} = nd m{j}", '=', 0, [(c, hj(c)) for c in N] + [(m[j], -nd)])
        add(f"H2: sum_pi h{j} k = ((q+1)n+q^2-q) m{j}", '=', 0, [(c, hj(c) * c[2]) for c in N] + [(m[j], -(6 * n + 25))])
        add(f"H3: sum_pi h{j} r = (nd+q j) m{j}", '=', 0, [(c, hj(c) * c[1]) for c in N] + [(m[j], -(nd + q * j))])
        add(f"H4: sum_pi C(h{j},2) = (q+1) M{j}{j}", '=', 0, [(c, C2(hj(c), 2)) for c in N] + [(M[(j, j)], -6)])
        for i in range(1, j): add(f"H5: sum_pi h{i} h{j} = (q+1) M{i}{j}", '=', 0, [(c, c[3][i - 1] * hj(c)) for c in N] + [(M[(i, j)], -6)])
        add(f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}", '=', 0, [(M[(min(i, j), max(i, j))], 1) for i in range(1, D + 1) if i != j] + [(M[(j, j)], 2), (m[j], -(n - 1))])
        add(f"H7: 2 M{j}{j} <= (n-1) m{j}", '<', 0, [(M[(j, j)], 2), (m[j], -(n - 1))])
        for a in range(1, n): add(f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)", '<', C2(a + 1, 2), [(M[(j, j)], -1), (m[j], a)])
    pairs = sorted({(t[0], t[1]) for t in T}); rk_cols = defaultdict(list)
    for c in N: rk_cols[(c[1], c[2])].append(c)
    # classes: 5-multisets of (r,k) with sum r = 6, sum k = n
    ctypes = []
    def crec(start, left, sr, sk, cur):
        if left == 0:
            if sk == n and sr == 6: ctypes.append(tuple(cur))
            return
        for i in range(start, len(pairs)):
            r, k = pairs[i]
            if sk + k > n or sr + r > 6: continue
            cur.append((r, k)); crec(i, left - 1, sr + r, sk + k, cur); cur.pop()
    crec(0, 5, 0, 0, []); Cc = [('C', t) for t in ctypes]
    add("classes: sum C = nd", '=', nd, [(c, 1) for c in Cc])
    for rk in pairs: add(f"class link N{rk}", '=', 0, [(c, c[1].count(rk)) for c in Cc if rk in c[1]] + [(c, -1) for c in rk_cols[rk]])
    # directions
    utypes = sorted(Ft); A = [('A', u) for u in utypes]
    add("dirs: sum A = nd", '=', nd, [(c, 1) for c in A])
    def phi(t, v):
        rho = sum(1 for r, k in t if r >= 1); n1 = sum(1 for r, k in t if r == 1); return (6 - n1) * (rho == v) + n1 * (rho - 1 == v)
    for v in range(1, 6): add(f"dir flags a=={v}", '=', 0, [(c, c[1].count(v)) for c in A] + [(c, -phi(c[1], v)) for c in Cc if phi(c[1], v)])
    wvecs = []
    def lrec(j, left, cur):
        if j == 0:
            if left == 0 and sum(cur) <= 25: wvecs.append(tuple(reversed(cur)))
            return
        for cnt in range(1 if j == 5 else 0, left // j + 1): lrec(j - 1, left - j * cnt, cur + [cnt])
    lrec(5, n, []); P = [('P', w) for w in wvecs]
    add("dirs: sum P = nd", '=', nd, [(c, 1) for c in P])
    for j in range(1, 6): add(f"dir link L{j}", '=', 0, [(c, c[1][j - 1]) for c in P if c[1][j - 1]] + [(L[j], -1)])
    for theta in range(bm + 1, 26): add(f"dominance lambda>={theta}", '<', 0, [(c, -1) for c in P if sum(c[1]) >= theta] + [(c, 1) for c in A if Ft[c[1]] >= theta])
    # witness-line pencils (bucketed items) and degree multisets
    items = []
    for r in range(1, 7):
        ks = sorted(k for (rr, k) in pairs if rr == r)
        if not ks: continue
        lo = ks[0]
        for (a, b) in [(lo, lo), (lo + 1, lo + 2), (lo + 3, ks[-1])]:
            kk = [k for k in ks if a <= k <= b]
            if kk: items.append((r, min(kk), max(kk), tuple(kk)))
    wt = []
    def wrec(start, left, sr, slo, shi, cur):
        if left == 0:
            if slo <= n + 25 <= shi and 6 <= sr <= 1 + q * D: wt.append(tuple(cur))
            return
        for i in range(start, len(items)):
            r, lo, hi, kk = items[i]
            if slo + lo > n + 25 or sr + r > 1 + q * D: continue
            cur.append(i); wrec(i, left - 1, sr + r, slo + lo, shi + hi, cur); cur.pop()
    wrec(0, 6, 0, 0, 0, []); W = [('W', t) for t in wt]; s_of = {t: sum(items[i][0] for i in t) - 1 for t in wt}
    add("wlines: sum W = nd", '=', nd, [(c, 1) for c in W])
    for ii, (r, lo, hi, kk) in enumerate(items):
        add(f"wline link r*N[{r},{lo}..{hi}]", '=', 0, [(c, c[1].count(ii)) for c in W if ii in c[1]] + [(c, -r) for k in kk for c in rk_cols[(r, k)]])
    add("wlines: sum s_l = sum d_p^2", '=', 0, [(c, s_of[c[1]]) for c in W] + [(m[j], -j * j) for j in m])
    dms = list(itertools.combinations_with_replacement(range(1, D + 1), 5)); V = [('V', dm) for dm in dms]
    add("wdeg: sum V = nd", '=', nd, [(c, 1) for c in V])
    for j in range(1, D + 1): add(f"wdeg: sum_l #(deg {j} on l) = {j} m{j}", '=', 0, [(c, c[1].count(j)) for c in V if j in c[1]] + [(m[j], -j)])
    for s in range(5, 5 * D + 1): add(f"wdeg/wline link s={s}", '=', 0, [(c, 1) for c in V if sum(c[1]) == s] + [(c, -1) for c in W if s_of[c[1]] == s])
    for i in range(1, D + 1):
        for j in range(i, D + 1):
            coef = [(c, (C2(c[1].count(i), 2) if i == j else c[1].count(i) * c[1].count(j))) for c in V]
            coef = [(c, v) for c, v in coef if v]
            if coef: add(f"wdeg: witness pairs ({i},{j}) <= M{i}{j}", '<', 0, coef + [(M[(i, j)], -1)])
    if config is not None:
        r, u = config; rmax = r
        add(f"config ({r},{u}) present", '<', -1, [(c, -1) for c in N if c[1] == r and c[4] == u])
    if rmax is not None:
        idx = [c for c in N if c[1] > rmax]
        if idx: add(f"rmax<={rmax}", '<', 0, [(c, 1) for c in idx])
    cols = N + list(m.values()) + list(L.values()) + list(M.values()) + Cc + A + P + W + V
    return rows, cols
def verify(rows, cols, cert):
    names = {r[0] for r in rows}; missing = [k for k in cert if k not in names]
    if missing: return False, f'certificate rows not in my model: {missing[:5]}'
    y = {k: Fraction(v) for k, v in cert.items()}
    comb = defaultdict(Fraction); rhs = Fraction(0)
    for name, sense, b, d in rows:
        yi = y.get(name, 0)
        if yi == 0: continue
        if sense == '<' and yi < 0: return False, f'negative multiplier on <= row {name}'
        for c, v in d.items(): comb[c] += yi * v
        rhs += yi * b
    bad = [c for c in cols if comb[c] < 0]
    return (not bad and rhs < 0), f'rhs={rhs} min_col={min(comb[c] for c in cols)} bad_cols={len(bad)} e.g. {bad[:2]}'
if __name__ == '__main__':
    CE = '../lp/certs/'; results = {}
    jobs = [('lp_probe_rmax4.json', dict(rmax=4))] + [(f'lp_probe_cfg_{r}_{u}.json', dict(config=(r, u))) for r, u in [(5, 15), (5, 16), (5, 17), (5, 18), (5, 21), (6, 17), (6, 18), (6, 19), (6, 21), (6, 25)]]
    which = sys.argv[1:] 
    for fn, kw in jobs:
        if which and fn not in which: continue
        d = json.load(open(CE + fn))
        for key, rec in sorted(d.items(), key=lambda x: int(x[0].split(',')[1])):
            if 'cert' not in rec: results[fn + ':' + key] = rec['status']; continue
            D = int(key.split(',')[1]); rows, cols = build(D, **kw)
            ok, msg = verify(rows, cols, rec['cert'])
            results[fn + ':' + key] = ('VERIFIED' if ok else 'FAILED') + ' ' + msg
            print(fn, key, 'rows', len(rows), 'cols', len(cols), results[fn + ':' + key], flush=True)
    json.dump(results, open('out/cert_v_result.json', 'w'), indent=1)
    print('verified', sum(v.startswith('VERIFIED') for v in results.values()), 'failed', sum(v.startswith('FAILED') for v in results.values()))
