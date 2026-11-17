"""AUDIT: independent replay of the six session-VI certificates (D=4 for (5,15),(5,16),(6,17); (5,16) at D=5,6,7) with a
FROM-SCRATCH generator of the session-VI model (WRITEUP6): section-profile plane types (F5, from out/planar_raw.txt),
exact pencil items (EXACT_K=4) with the Lemma-C filter (F7), sixth-direction flags (F6 rows R1-R4, Lemma B / R4 capacity),
plus the session-V families. Nothing is imported from lp/src; only the saved certificates (keyed by row NAME) are read.
Column-wise verification: for every column of my universe the signed combination of the certificate multipliers must be
>= 0, multipliers on <= rows >= 0, and y^T b < 0 (exact integers)."""
import json, math, itertools, sys, time
from collections import defaultdict
C2 = math.comb
q, n = 5, 46; nd = 31; npl = 155; base = 25; f = [0, 5, 9, 12, 14, 15, 17]; bm = 17; RMIN = 5; EXACT_K = 4
prof = {}; secs = []
for ln in open('out/planar_raw.txt'):
    t = ln.split()
    if t[0] == 'PROF': prof[tuple(map(int, t[1:7]))] = int(t[7])
    elif t[0] == 'SEC':
        parts = ln[4:].split('|'); a = list(map(int, parts[0].split())); tt = tuple(map(int, parts[1].split()))
        k, r, u = a[:3]; c = tuple(a[3:3 + r]); e = ((tuple(map(int, parts[2].split())),) if len(parts) > 2 else None)
        secs.append((k, r, u, c, tt, e))
secs.sort(key=lambda x: x[:5])
Ft = {a: min(v for p, v in prof.items() if all(x >= y for x, y in zip(p, a))) for a in itertools.combinations_with_replacement(range(1, 6), 6)}
def hists(k, S, D):
    def rec(j, lk, lS, cur):
        if j == D:
            if lS == D * lk: yield tuple(cur + [lk])
            return
        for h in range(lk + 1):
            rk, rS = lk - h, lS - j * h
            if rS < (j + 1) * rk or rS > D * rk: continue
            yield from rec(j + 1, rk, rS, cur + [h])
    return list(rec(1, k, S, []))
def types_for(D):
    out = []; hc = {}
    for (k, r, u, c, tt, e) in secs:
        if k > 25 or r > 6 or k < max(f[r], 1): continue
        S = base + q * r
        if S < k or S > D * k: continue
        if (k, S) not in hc: hc[(k, S)] = hists(k, S, D)
        for h in hc[(k, S)]:
            tail = [sum(h[j - 1:]) for j in range(1, D + 2)]
            if all(sum(c[i - 1:]) <= (tail[i - 1] if i <= D else 0) for i in range(1, r + 1)): out.append((r, k, h, u, tt, e))
    return out
def cap3(item):
    r, lo, hi, kk = item
    if r == 5 and hi == 15: return 2
    if r == 1: return min(4, (hi - 5) // 3)
    if r == 2: return min(4, (hi - 9) // 2)
    if r == 3: return min(4, hi - 12)
    return 4
class Model:
    def __init__(self, D, rmax=None, config=None):
        t0 = time.time(); self.D = D; self.rows = {}   # name -> (sense, rhs)
        T = self.T = types_for(D)
        def row(name, sense, rhs): assert name not in self.rows, name; self.rows[name] = (sense, rhs)
        row("I2a: sum_pi 1", '=', npl); row("I2b: sum_pi r", '=', 6 * nd); row("I2c: sum_pi k", '=', nd * n)
        row("deg: sum m_j", '=', n); row("deg: sum j m_j", '=', nd * q); row("m_D >= 1", '<', -1)
        row("I3: sum C(r,2) = sum C(d,2)", '=', 0); row("I4: sum k^2", '=', 6 * n * n + 25 * n); row("I5: sum k r", '=', nd * n + 25 * nd)
        row("I6a: sum j L_j", '=', nd * n); row("I6b: sum C(j,2) L_j", '=', C2(n, 2)); row("I6c: sum L_j >= nd*BM", '<', -nd * bm); row("I6d: L_q >= nd", '<', -nd)
        row("I6e: sum C(k,3) - q sum C(j,3) L_j", '=', C2(n, 3))
        for j in range(1, 6): row(f"F5: sum_pi t{j} = (q+1) L{j}", '=', 0)
        for j in range(1, D + 1):
            row(f"H1: sum_pi h{j} = nd m{j}", '=', 0); row(f"H2: sum_pi h{j} k = ((q+1)n+q^2-q) m{j}", '=', 0); row(f"H3: sum_pi h{j} r = (nd+q j) m{j}", '=', 0)
            row(f"H4: sum_pi C(h{j},2) = (q+1) M{j}{j}", '=', 0)
            for i in range(1, j): row(f"H5: sum_pi h{i} h{j} = (q+1) M{i}{j}", '=', 0)
            row(f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}", '=', 0); row(f"H7: 2 M{j}{j} <= (n-1) m{j}", '<', 0)
            for a in range(1, n): row(f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)", '<', C2(a + 1, 2))
        pairs = self.pairs = sorted({(t[0], t[1]) for t in T})
        # classes
        ctypes = []
        def crec(start, left, sr, sk, cur):
            if left == 0:
                if sk == n and sr == 6: ctypes.append(tuple(cur))
                return
            for i in range(start, len(pairs)):
                r, k = pairs[i]
                if sk + k > n or sr + r > 6: continue
                cur.append((r, k)); crec(i, left - 1, sr + r, sk + k, cur); cur.pop()
        crec(0, 5, 0, 0, []); self.ctypes = ctypes
        row("classes: sum C = nd", '=', nd)
        for rk in pairs: row(f"class link N{rk}", '=', 0)
        self.utypes = sorted(Ft); row("dirs: sum A = nd", '=', nd)
        for v in range(1, 6): row(f"dir flags a=={v}", '=', 0)
        wvecs = []
        def lrec(j, left, cur):
            if j == 0:
                if left == 0 and sum(cur) <= 25: wvecs.append(tuple(reversed(cur)))
                return
            for cnt in range(1 if j == 5 else 0, left // j + 1): lrec(j - 1, left - j * cnt, cur + [cnt])
        lrec(5, n, []); self.wvecs = wvecs; row("dirs: sum P = nd", '=', nd)
        for j in range(1, 6): row(f"dir link L{j}", '=', 0)
        for theta in range(bm + 1, 26): row(f"dominance lambda>={theta}", '<', 0)
        # exact pencil items (EXACT_K = 4) and pencil types with the Lemma C (F7) filter
        items = []
        for r in range(1, 7):
            ks = sorted(k for (rr, k) in pairs if rr == r)
            if not ks: continue
            lo = ks[0]
            for (a, b) in [(lo + i, lo + i) for i in range(EXACT_K)] + [(lo + EXACT_K, ks[-1])]:
                kk = [k for k in ks if a <= k <= b]
                if kk: items.append((r, min(kk), max(kk), tuple(kk)))
        self.items = items; wt = []
        def wrec(start, left, sr, slo, shi, cur):
            if left == 0:
                if slo <= n + 25 <= shi and 6 <= sr <= 1 + q * D:
                    s = sr - 1
                    if all(s >= base + q * items[i][0] - D * (items[i][2] - q) for i in cur): wt.append(tuple(cur))
                return
            for i in range(start, len(items)):
                r, lo, hi, kk = items[i]
                if slo + lo > n + 25 or sr + r > 1 + q * D: continue
                cur.append(i); wrec(i, left - 1, sr + r, slo + lo, shi + hi, cur); cur.pop()
        wrec(0, 6, 0, 0, 0, []); self.wt = wt; self.s_of = [sum(items[i][0] for i in t) - 1 for t in wt]
        row("wlines: sum W = nd", '=', nd)
        for (r, lo, hi, kk) in items: row(f"wline link r*N[{r},{lo}..{hi}]", '=', 0)
        row("wlines: sum s_l = sum d_p^2", '=', 0)
        self.dms = list(itertools.combinations_with_replacement(range(1, D + 1), 5)); row("wdeg: sum V = nd", '=', nd)
        for j in range(1, D + 1): row(f"wdeg: sum_l #(deg {j} on l) = {j} m{j}", '=', 0)
        for s in range(5, 5 * D + 1): row(f"wdeg/wline link s={s}", '=', 0)
        for i in range(1, D + 1):
            for j in range(i, D + 1):
                if any((C2(dm.count(i), 2) if i == j else dm.count(i) * dm.count(j)) for dm in self.dms): row(f"wdeg: witness pairs ({i},{j}) <= M{i}{j}", '<', 0)
        # F6 flags (rmin = 5)
        prof_of = defaultdict(set)
        for t in T:
            if t[0] >= RMIN and t[5]:
                for sp in t[5]: prof_of[(t[0], t[1])].add(sp)
        kinds = set(); self.gco = defaultdict(dict)   # (rk, rk2) -> {class index: -g}
        for ci, ct in enumerate(ctypes):
            cnt = defaultdict(int)
            for rk in ct: cnt[rk] += 1
            for rk, m1 in cnt.items():
                if rk[0] < RMIN or rk not in prof_of: continue
                for rk2, m2 in cnt.items():
                    if rk2[0] < 1: continue
                    g = m1 * (m2 - (1 if rk2 == rk else 0)) * rk2[0]
                    if g: self.gco[(rk, rk2)][ci] = -g; kinds.update((rk[0], rk[1], sp, rk2) for sp in prof_of[rk])
        self.kinds = sorted(kinds); self.srk2 = sorted({(kd[2], kd[3]) for kd in kinds}); srk2set = set(self.srk2)
        def ok_item(ii, sv): r, lo, hi, kk = items[ii]; return hi >= f[r] + max(0, sv - r + 1)
        cache = {}
        def matchable(rest, e):
            key = (rest, e)
            if key in cache: return cache[key]
            masks = [sum(1 << j for j in range(5) if ok_item(ii, e[j])) for ii in rest]
            reach = {0}
            for mk in masks:
                reach = {rm | (1 << j) for rm in reach for j in range(5) if (mk >> j) & 1 and not (rm >> j) & 1}
                if not reach: break
            cache[key] = bool(reach); return cache[key]
        sps = {kd[2] for kd in kinds}; self.X2 = []   # ((sp, rk2), i)
        for i, w in enumerate(wt):
            for pos in {w[p]: p for p in range(len(w))}.values():
                r2, lo, hi, kk = items[w[pos]]
                if r2 < 1: continue
                rest = tuple(sorted(w[:pos] + w[pos + 1:]))
                for sp in sps:
                    if matchable(rest, sp):
                        for k2 in kk:
                            if (sp, (r2, k2)) in srk2set: self.X2.append(((sp, (r2, k2)), i))
        self.X2 = list(dict.fromkeys(self.X2))
        self.types_rks = defaultdict(list)
        for j, t in enumerate(T):
            if t[0] >= RMIN and t[5]:
                for sp in set(t[5]): self.types_rks[(t[0], t[1], sp)].append((j, -t[5].count(sp)))
        for key in self.types_rks: row(f"F6 R1 flags (r={key[0]},k={key[1]},s={key[2]})", '=', 0)
        for (rk, rk2) in self.gco: row(f"F6 R2 flags ({rk} hosted by {rk2})", '=', 0)
        by_srk2 = {(kd[2], kd[3]) for kd in kinds} | {k for k, i in self.X2}
        for key in sorted(by_srk2): row(f"F6 link flags (s={key[0]},mate={key[1]})", '=', 0)
        by_w = defaultdict(list)
        for (key, i) in self.X2: by_w[i].append(key)
        self.cap = {}; self.xmax = {}
        full3 = (3,) * 5
        for i in by_w:
            cap = sum((6 - items[ii][0]) // RMIN for ii in wt[i] if items[ii][0] >= 1); self.cap[i] = cap
            row(f"F6 R3 flags on pencil {i} <= {cap} W{i}", '<', 0)
        for i in by_w:
            if not any(key[0] == full3 for key in by_w[i]): continue
            tot = sum(cap3(items[ii]) for ii in wt[i]); xmax = 0
            while xmax + 1 <= 6 and 5 * (xmax + 1) - C2(xmax + 1, 2) <= tot: xmax += 1
            self.xmax[i] = xmax; row(f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {xmax} W{i}", '<', 0)
        self.config = config; self.rmax = rmax
        if config is not None: row(f"config ({config[0]},{config[1]}) present", '<', -1); self.rmax = config[0]
        if self.rmax is not None and any(t[0] > self.rmax for t in T): row(f"rmax<={self.rmax}", '<', 0)
        self.info = dict(types=len(T), classes=len(ctypes), pencils=len(wt), items=len(items), X1=len(self.kinds), X2=len(self.X2), rows=len(self.rows), build_s=round(time.time() - t0, 1))
    def verify(self, cert):
        D = self.D; y = {k: int(v) for k, v in cert.items() if int(v)}
        missing = [k for k in y if k not in self.rows]
        if missing: return False, f'certificate rows not in my model: {missing[:5]} (+{len(missing)-5 if len(missing)>5 else 0})'
        neg = [k for k in y if self.rows[k][0] == '<' and y[k] < 0]
        if neg: return False, f'negative multiplier on <= rows {neg[:3]}'
        rhs = sum(y[k] * self.rows[k][1] for k in y)
        Y = lambda name: y.get(name, 0)
        bad = []; mincol = None; ncol = 0
        def check(v, col):
            nonlocal mincol, ncol
            ncol += 1
            if mincol is None or v < mincol: mincol = v
            if v < 0 and len(bad) < 5: bad.append((col, v))
        items = self.items; pairs_items = defaultdict(list)
        for ii, (r, lo, hi, kk) in enumerate(items):
            for k in kk: pairs_items[(r, k)].append(ii)
        yH = {j: (Y(f"H1: sum_pi h{j} = nd m{j}"), Y(f"H2: sum_pi h{j} k = ((q+1)n+q^2-q) m{j}"), Y(f"H3: sum_pi h{j} r = (nd+q j) m{j}"), Y(f"H4: sum_pi C(h{j},2) = (q+1) M{j}{j}")) for j in range(1, D + 1)}
        yH5 = {(i, j): Y(f"H5: sum_pi h{i} h{j} = (q+1) M{i}{j}") for j in range(1, D + 1) for i in range(1, j)}
        yF5 = [Y(f"F5: sum_pi t{j} = (q+1) L{j}") for j in range(1, 6)]
        ycl = {rk: Y(f"class link N{rk}") for rk in self.pairs}
        ywl = [Y(f"wline link r*N[{r},{lo}..{hi}]") for (r, lo, hi, kk) in items]
        yR1 = {key: Y(f"F6 R1 flags (r={key[0]},k={key[1]},s={key[2]})") for key in self.types_rks}
        yI2a, yI2b, yI2c, yI3, yI4, yI5, yI6e = (Y(s) for s in ("I2a: sum_pi 1", "I2b: sum_pi r", "I2c: sum_pi k", "I3: sum C(r,2) = sum C(d,2)", "I4: sum k^2", "I5: sum k r", "I6e: sum C(k,3) - q sum C(j,3) L_j"))
        ycfg = Y(f"config ({self.config[0]},{self.config[1]}) present") if self.config else 0
        yrmax = Y(f"rmax<={self.rmax}") if self.rmax is not None else 0
        for t in self.T:
            r, k, h, u, tt, e = t
            v = yI2a + yI2b * r + yI2c * k + yI3 * C2(r, 2) + yI4 * k * k + yI5 * k * r + yI6e * C2(k, 3) - ycl[(r, k)]
            v += sum(yF5[j] * tt[j] for j in range(5))
            for j in range(1, D + 1):
                hj = h[j - 1]
                if hj:
                    a, b, c, d = yH[j]; v += hj * (a + b * k + c * r) + d * C2(hj, 2)
                    for i in range(1, j): v += yH5[(i, j)] * h[i - 1] * hj
            for ii in pairs_items[(r, k)]: v -= r * ywl[ii]
            if r >= RMIN and e:
                for sp in set(e): v -= e.count(sp) * yR1[(r, k, sp)]
            if self.config and r == self.config[0] and u == self.config[1]: v -= ycfg
            if self.rmax is not None and r > self.rmax: v += yrmax
            check(v, ('N', t))
        for j in range(1, D + 1):
            v = Y("deg: sum m_j") + j * Y("deg: sum j m_j") - (Y("m_D >= 1") if j == D else 0) - C2(j, 2) * yI3 - nd * yH[j][0] - (6 * n + 25) * yH[j][1] - (nd + 5 * j) * yH[j][2]
            v += -(n - 1) * Y(f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}") - (n - 1) * Y(f"H7: 2 M{j}{j} <= (n-1) m{j}") + sum(a * Y(f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)") for a in range(1, n))
            v += -j * j * Y("wlines: sum s_l = sum d_p^2") - j * Y(f"wdeg: sum_l #(deg {j} on l) = {j} m{j}")
            check(v, ('m', j))
        for j in range(1, 6):
            v = j * Y("I6a: sum j L_j") + C2(j, 2) * Y("I6b: sum C(j,2) L_j") - Y("I6c: sum L_j >= nd*BM") - (Y("I6d: L_q >= nd") if j == 5 else 0) - 5 * C2(j, 3) * yI6e - 6 * yF5[j - 1] - Y(f"dir link L{j}")
            check(v, ('L', j))
        for i in range(1, D + 1):
            for j in range(i, D + 1):
                if i == j:
                    v = -6 * yH[j][3] + 2 * Y(f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}") + 2 * Y(f"H7: 2 M{j}{j} <= (n-1) m{j}") - sum(Y(f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)") for a in range(1, n))
                else: v = -6 * yH5[(i, j)] + Y(f"H6: sum_i M_i{i} + 2 M{i}{i} = (n-1) m{i}") + Y(f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}")
                v -= Y(f"wdeg: witness pairs ({i},{j}) <= M{i}{j}")
                check(v, ('M', i, j))
        def phi(t, v):
            rho = sum(1 for r, k in t if r >= 1); n1 = sum(1 for r, k in t if r == 1); return (6 - n1) * (rho == v) + n1 * (rho - 1 == v)
        ydf = [Y(f"dir flags a=={v}") for v in range(1, 6)]
        yR2 = {key: Y(f"F6 R2 flags ({key[0]} hosted by {key[1]})") for key in self.gco}
        r2_by_class = defaultdict(list)
        for key, d in self.gco.items():
            for ci, g in d.items(): r2_by_class[ci].append((key, g))
        ycs = Y("classes: sum C = nd")
        for ci, ct in enumerate(self.ctypes):
            v = ycs + sum(ycl[rk] for rk in ct) - sum(ydf[vv - 1] * phi(ct, vv) for vv in range(1, 6)) + sum(g * yR2[key] for key, g in r2_by_class[ci])
            check(v, ('C', ct))
        ydom = {theta: Y(f"dominance lambda>={theta}") for theta in range(bm + 1, 26)}
        for u in self.utypes: check(Y("dirs: sum A = nd") + sum(ydf[vv - 1] * u.count(vv) for vv in range(1, 6)) + sum(ydom[th] for th in ydom if Ft[u] >= th), ('A', u))
        ydl = [Y(f"dir link L{j}") for j in range(1, 6)]
        for w in self.wvecs: check(Y("dirs: sum P = nd") + sum(ydl[j] * w[j] for j in range(5)) - sum(ydom[th] for th in ydom if sum(w) >= th), ('P', w))
        yws, ysl = Y("wlines: sum W = nd"), Y("wlines: sum s_l = sum d_p^2"); ylink = {s: Y(f"wdeg/wline link s={s}") for s in range(5, 5 * D + 1)}
        for i, t in enumerate(self.wt):
            v = yws + sum(ywl[ii] for ii in t) + ysl * self.s_of[i] - ylink[self.s_of[i]]
            if i in self.cap: v -= self.cap[i] * Y(f"F6 R3 flags on pencil {i} <= {self.cap[i]} W{i}")
            if i in self.xmax: v -= self.xmax[i] * Y(f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {self.xmax[i]} W{i}")
            check(v, ('W', i))
        yvs = Y("wdeg: sum V = nd"); yvd = [Y(f"wdeg: sum_l #(deg {j} on l) = {j} m{j}") for j in range(1, D + 1)]
        ywp = {(i, j): Y(f"wdeg: witness pairs ({i},{j}) <= M{i}{j}") for i in range(1, D + 1) for j in range(i, D + 1)}
        for dm in self.dms:
            v = yvs + sum(yvd[j - 1] * dm.count(j) for j in set(dm)) + ylink[sum(dm)]
            for i in range(1, D + 1):
                for j in range(i, D + 1):
                    c = C2(dm.count(i), 2) if i == j else dm.count(i) * dm.count(j)
                    if c: v += c * ywp[(i, j)]
            check(v, ('V', dm))
        ylk = {key: Y(f"F6 link flags (s={key[0]},mate={key[1]})") for key in self.srk2}
        for kd in self.kinds: check(yR1[(kd[0], kd[1], kd[2])] + yR2.get(((kd[0], kd[1]), kd[3]), 0) + ylk[(kd[2], kd[3])], ('X1', kd))
        full3 = (3,) * 5
        for (key, i) in self.X2:
            v = -ylk[key] + Y(f"F6 R3 flags on pencil {i} <= {self.cap[i]} W{i}")
            if key[0] == full3 and i in self.xmax: v += Y(f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {self.xmax[i]} W{i}")
            check(v, ('X2', key, i))
        ok = not bad and rhs < 0
        return ok, f'rhs={rhs} min_col={mincol} columns={ncol} bad={bad[:3]}'
if __name__ == '__main__':
    CE = '../lp/certs/'; results = {}
    jobs = [('lp_probe_f7_cfg_5_15.json', (5, 15)), ('lp_probe_f7_cfg_6_17.json', (6, 17)), ('lp_probe_f7_cfg_5_16.json', (5, 16))]
    want = sys.argv[1:]
    for fn, cfg in jobs:
        d = json.load(open(CE + fn))
        for key, rec in sorted(d.items(), key=lambda x: int(x[0].split(',')[1])):
            if 'cert' not in rec: continue
            D = int(key.split(',')[1])
            if want and f'{cfg[0]},{cfg[1]},{D}' not in want: continue
            M = Model(D, config=cfg); print(fn, key, M.info, flush=True)
            ok, msg = M.verify(rec['cert']); results[fn + ':' + key] = ('VERIFIED ' if ok else 'FAILED ') + msg
            print('  ->', results[fn + ':' + key], ' their y^T b:', rec['cert_msg'].split('y^T b =')[-1].strip(), flush=True)
            json.dump(results, open('out/cert_vi_result.json', 'w'), indent=1)
    print('verified', sum(v.startswith('VERIFIED') for v in results.values()), 'of', len(results))
