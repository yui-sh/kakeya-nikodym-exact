"""AUDIT: independent replay of the EIGHT ROOT certificates of the K(3,5) >= 47 chain (D = 5: (5,15), (6,17) roots; D = 6: both
roots; D = 7: (5,15) root, 617:0, 617:27 and the m1 >= 1 domain certificate) with a FROM-SCRATCH generator of the full root models:
the session-VI model of check_cert_vi.py (my own planar tables, F5, exact pencil items + Lemma C filter F7, F6 flags R1-R4) plus
  * prefix caps (plane types with h_j > max_j over the node's raw degree vectors or r > r_max are dropped: WRITEUP11), the
    maximal-only filter (AUDIT §3.13) and the plane-budget filter c(r,k) <= 541 (§3.16) for the session-XIV models,
  * F8 = Lemma D(a): incidence columns G[(r,k,h),dm] (dm subset of h) and Z[dm,pencil] with rows Y1, Y2, Z1, Z2 (WRITEUP7 §2),
  * F9 = Lemma C on the mate plane: deficit rows on the X2 flag columns (WRITEUP8 §1),
  * the lifts and cut families, ALL PRICED HERE from my own planar/geometry catalogues (no price table is read for the
    verification; the research tables are only compared afterwards, informationally): coordinate lift XI, pairwise lift XIII,
    XIV additive subset lift, within-line pair/triple costs (exhaustive colourings), r = 6 extra-point costs, point crossing /
    coplanar-triple / triples-minus-quadruples / rich-secant / mixed / degree-subset costs (allocation enumeration with Hall
    thresholds; global PG(2,5) bounds from my own check_projective.py output), star bounds, rounded plane budget, degree hull,
    pair-bound envelopes, parallel-class degree rounding, the compact-objective rule,
  * the node rows (config present, r_max, prefix bounds on m_j, domain row).
Nothing is imported from lp/src. The certificates are POSITIONAL multiplier lists, so the research row-name list dumped once by
dump_research_rownames.py / dump_research_rownames14.py (the only research-code helpers) is used as the position -> row-name
map; my rows are keyed by the same names. Before the Farkas check the two row lists are compared (names, sense, rhs, and the
number of nonzeros per row counted over MY columns) and the column counts per family are compared with the research model.
Farkas check: every multiplier on a <= row is >= 0, the combination y^T A is >= 0 on EVERY column of my universe, and
y^T b < 0, all in exact integer arithmetic.
Usage (from audit/):  python3 check_root.py 515 | 617 | D6_515 | D6_617 | D7_515 | D7_617_27 | D7_617_0 | D7_617_dom
   (output: out/root_<tag>.log, out/root_result.json; needs out/research_rownames_<tag>.json from the dump helpers)"""
import json, math, itertools, sys, time, os
from collections import defaultdict
from functools import lru_cache
import check_cert_vi as vi                      # my own session-VI tables/enumerations (planar_raw.txt, Ft, hists, cap3)
C2 = math.comb
q, n, nd, npl, base, bm, RMIN, EXACT_K = vi.q, vi.n, vi.nd, vi.npl, vi.base, vi.bm, vi.RMIN, vi.EXACT_K
f = vi.f; Ft = vi.Ft; secs = vi.secs
CORE = ((5, 15), (6, 17))
D = 5   # set by RootModel; module-level default

def raw_vectors(D):
    """all degree vectors m with sum m = 46, sum j m_j = 155, m_D >= 1 (first D-2 coordinates free, last two solved)."""
    out = []
    def rec(prefix, left):
        if len(prefix) == D - 2:
            s = 155 - sum((i + 1) * v for i, v in enumerate(prefix)); last = s - (D - 1) * left; pen = left - last
            if last >= 1 and pen >= 0: out.append(prefix + (pen, last))
            return
        for v in range(left + 1): rec(prefix + (v,), left - v)
    rec((), n); return out
def types_dedup(D):
    """my session-VI type universe (r,k,h,u,t,e), one column per distinct tuple (check_cert_vi.types_for lists a tuple once per
    crossing profile c that admits it; the research model keys columns by the tuple only)."""
    return list(dict.fromkeys(vi.types_for(D)))
def plane_cost(r, k): return 10 * k * k - (150 + 17 * r) * k + 560 + 143 * r   # WRITEUP14 budget polynomial, sum over planes = 541

# ---------------- exact pricer for additive (b, d)-costs (WRITEUP11/13: min over degree assignments with d >= max(1,b)) ----------
profiles_of = defaultdict(set)   # (r,k,u,t,e) -> set of multiplicity profiles (count of points with b = 0..r)
for (k, r, u, c, tt, e) in secs: profiles_of[(r, k, u, tt, e)].add((k - u,) + tuple(c))
@lru_cache(None)
def assign_min(profile, h, cost):
    D = len(h)
    """min over assignments of the degree multiset h (h[d-1] points of degree d) to the points (profile[b] points of inside
    multiplicity b, each needing d >= max(1,b)) of sum cost[b][d-1]; None if no assignment exists. Exhaustive group DP."""
    groups = [(b, cnt) for b, cnt in enumerate(profile) if cnt]
    @lru_cache(None)
    def rec(gi, rem):
        if gi == len(groups): return 0 if not any(rem) else None
        b, cnt = groups[gi]; ds = [d for d in range(max(1, b), D + 1) if rem[d - 1]]; best = None
        def comp(i, left, rem2, acc):
            nonlocal best
            if i == len(ds):
                if left == 0:
                    v = rec(gi + 1, rem2)
                    if v is not None and (best is None or acc + v < best): best = acc + v
                return
            d = ds[i]
            for x in range(min(left, rem2[d - 1]) + 1):
                r2 = list(rem2); r2[d - 1] -= x; comp(i + 1, left - x, tuple(r2), acc + x * cost[b][d - 1])
        comp(0, cnt, rem, 0); return best
    return rec(0, tuple(h))
def price(t, S, sign):
    """mu_a for the full type a = t and the objective sign * sum_p b_p [d(p) in S]: minimum over every multiplicity profile
    realisable by (r,k,u,t,e) and every admissible degree assignment with histogram h."""
    r, k, h, u, tt, e = t; key = (r, k, u, tt, e); DD = len(h)
    cost = tuple(tuple(sign * b * (d in S) for d in range(1, DD + 1)) for b in range(r + 1))
    vals = [v for p in profiles_of[key] if (v := assign_min(p, tuple(h), cost)) is not None]
    assert vals, ('no admissible assignment for a type in the universe', t)
    return min(vals)
OBJ_XI = [((j,), s) for j in range(1, 6) for s in (1, -1)]
OBJ_XIII = [((i, j), s) for i, j in itertools.combinations(range(1, 6), 2) for s in (1, -1)]
def obj_subsets(D):
    """WRITEUP14 additive subset family: singletons, then all subsets of size 2..D//2, one of each complementary pair at size D/2."""
    out = [(i,) for i in range(1, D + 1)]
    for sz in range(2, D // 2 + 1):
        for S in itertools.combinations(range(1, D + 1), sz):
            comp = tuple(i for i in range(1, D + 1) if i not in S)
            if 2 * sz == D and S > comp: continue
            out.append(S)
    return [(S, s) for S in out for s in (1, -1)]
@lru_cache(maxsize=300000)
def alloc_extrema(profile, h, S, costs):
    """Additive costs that depend on a point only through its inside multiplicity b and on whether its degree lies in S:
    min and max of sum_b costs[b] x_b over allocations x_b (= # points of multiplicity b given a degree in S), 0 <= x_b <= profile[b],
    sum x_b = #S-degrees, such that the S-points can be matched to the S-degrees and the rest to the other degrees (Hall threshold
    inequalities, AUDIT §3.11/§3.12). None if the histogram is not realisable at all."""
    r = len(profile) - 1; DD = len(h)
    nS = sum(h[d - 1] for d in S)
    degS = tuple(sum(h[d - 1] for d in S if d >= b) for b in range(r + 2)); degO = tuple(sum(h[d - 1] for d in range(1, DD + 1) if d not in S and d >= b) for b in range(r + 2))
    tot = tuple(sum(profile[b:]) for b in range(r + 2))
    best = [None, None]
    def rec(b, left, selected_tail, acc):
        # b descends from r to 0; selected_tail = sum of x_{b'} for b' > b
        if b < 0:
            if left == 0:
                if best[0] is None or acc < best[0]: best[0] = acc
                if best[1] is None or acc > best[1]: best[1] = acc
            return
        for x in range(min(profile[b], left) + 1):
            st = selected_tail + x
            if b >= 1 and (st > degS[b] or tot[b] - st > degO[b]): continue   # Hall at threshold b for both sides
            rec(b - 1, left - x, st, acc + costs[b] * x)
    rec(r, nS, 0, 0)
    return None if best[0] is None else tuple(best)
def price_alloc(t, S, sign, costs):
    """mu_a = min over the realisable multiplicity profiles of the full type and all admissible assignments of
    sign * sum_p costs[b_p] [d(p) in S]  (costs[b] >= 0 nondecreasing: b, C(b,2) or C(b,3))."""
    r, k, h, u, tt, e = t; key = (r, k, u, tt, e); vals = []
    for p in profiles_of[key]:
        ex = alloc_extrema(p, tuple(h), tuple(S), costs)
        if ex is not None: vals.append(ex[0] if sign == 1 else -ex[1])
    assert vals, ('no admissible assignment for a type in the universe', t)
    return min(vals)
def obj_withinline(D): return [(S, s) for sz in range(1, D) for S in itertools.combinations(range(1, D + 1), sz) for s in (1, -1)]
def jkey(a): return json.dumps(a, separators=(',', ':'))
# ---------------- within-line pair costs (WRITEUP14 "Stronger within-line pair costs"; AUDIT §3.12): exhaustive colourings ----------
geo_masks = defaultdict(list)   # (r,k,u,t,e) -> list of mask tuples (one per geometry of my catalogue, extra points = mask 0)
for ln in open('out/geometry_raw.txt'):
    a, t, e, mm = ln.split('|'); k, r, u = map(int, a.split()); t = tuple(map(int, t.split())); e = tuple(map(int, e.split())) or None
    geo_masks[(r, k, u, t, ((e,) if e else None))].append(tuple(map(int, mm.split())))
@lru_cache(None)
def colour_extrema(masks, r, power=2):
    """all 2^k colourings X of the points of one geometry (Gray code); for every profile (x_1..x_r) = #selected points of inside
    multiplicity b, the min and max of sum_lines C(|X ∩ line|, power)."""
    n_ = len(masks); bc = [bin(m).count('1') for m in masks]; line = [0] * r; cnt = [0] * (r + 1); cost = 0; prev = 0; out = {}
    for num in range(1 << n_):
        g = num ^ (num >> 1)
        if num:
            ch = g ^ prev; p = ch.bit_length() - 1; add = bool(g & ch); m = masks[p]; cnt[bc[p]] += 1 if add else -1
            for i in range(r):
                if m >> i & 1:
                    if add: cost += C2(line[i], power - 1); line[i] += 1
                    else: line[i] -= 1; cost -= C2(line[i], power - 1)
        prev = g; pr = tuple(cnt[1:])
        if pr in out:
            lo, hi = out[pr]
            if cost < lo or cost > hi: out[pr] = (min(lo, cost), max(hi, cost))
        else: out[pr] = (cost, cost)
    return out
def hall_ok(sel, degs_by_b):
    """a set of points with sel[b] points of multiplicity b (b >= 1; b = 0 points need only d >= 1) can be matched to a degree
    multiset iff for every threshold b the points with multiplicity >= b are at most the degrees >= b (Hall, nested sets)."""
    return all(sum(sel[b - 1:]) <= degs_by_b[b] for b in range(1, len(sel) + 1))
@lru_cache(None)
def price_withinline(key, h, S, sign, power=2):
    r = key[0]; DD = len(h); extra = key[1] - key[2]
    nS = sum(h[d - 1] for d in S); nO = sum(h) - nS
    degS = {b: sum(h[d - 1] for d in S if d >= b) for b in range(1, r + 1)}; degO = {b: sum(h[d - 1] for d in range(1, DD + 1) if d not in S and d >= b) for b in range(1, r + 1)}
    best = None
    for mm in geo_masks[key]:
        tot = [0] * (r + 1)
        for m in mm: tot[bin(m).count('1')] += 1
        mm_in = tuple(m for m in mm if m)          # extra points (mask 0) lie on no inside line: cost-free, handled by x0
        for pr, (lo, hi) in colour_extrema(mm_in, r, power).items():
            x0 = nS - sum(pr)                      # selected extra points
            if not 0 <= x0 <= extra: continue
            oth = tuple(tot[b] - pr[b - 1] for b in range(1, r + 1))
            if nO - sum(oth) != extra - x0: continue
            if hall_ok(pr, degS) and hall_ok(oth, degO):
                v = lo if sign == 1 else -hi
                if best is None or v < best: best = v
    return best

RAW_COUNT = {5: 2599, 6: 25622, 7: 137209}
class RootModel:
    def __init__(self, spec):
        global D
        F = self.F = spec; D = spec['D']; cfg = spec['cfg']; prefix = spec.get('prefix', ())
        t0 = time.time(); self.cfg = cfg; rr, uu = (5, 15) if cfg == '515' else (6, 17); self.rr, self.uu = rr, uu
        self.rmax = rr; raw = raw_vectors(D); assert len(raw) == RAW_COUNT[D]
        self.xiv = (D >= 6)                       # session-XIV model (D >= 6): 'XIV' node rows, degree hull
        self.pair = F.get('xiii', False)          # session-XIII pairwise lift (D = 5, (6,17) root only)
        cand = [m for m in raw if m[:len(prefix)] == tuple(prefix) and m[0] >= F.get('domain_m1', 0)]; assert cand; self.cand = cand
        cap = self.cap = tuple(max(m[j] for m in cand) for j in range(D)); self.floor = tuple(min(m[j] for m in cand) for j in range(D))
        T = [t for t in types_dedup(D) if t[0] <= rr and all(hj <= cj for hj, cj in zip(t[2], cap))]
        if F.get('maximal'): T = [t for t in T if t[0] != rr or t[3] == uu]      # maximal-only (AUDIT §3.13): every r_max-plane has union uu
        if F.get('budget'): T = [t for t in T if plane_cost(t[0], t[1]) <= 541]  # plane budget (AUDIT §3.16): every actual plane has cost <= 541
        self.T = T
        self.rows = {}; self.order = []
        def row(name, sense, rhs):
            assert name not in self.rows, name; self.rows[name] = (sense, rhs); self.order.append(name)
        # ---- session-V/VI rows (same statements as check_cert_vi.Model; see AUDIT.md §3) ----
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
        items = []
        for r in range(1, 7):
            ks = sorted(k for (r2, k) in pairs if r2 == r)
            if not ks: continue
            lo = ks[0]
            for (a, b) in [(lo + i, lo + i) for i in range(EXACT_K)] + [(lo + EXACT_K, ks[-1])]:
                kk = [k for k in ks if a <= k <= b]
                if kk: items.append((r, min(kk), max(kk), tuple(kk)))
        self.items = items; self.item_of = {}
        for ii, (r, lo, hi, kk) in enumerate(items):
            for k in kk: assert (r, k) not in self.item_of; self.item_of[(r, k)] = ii
        wt = []
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
        self.wpair_rows = set()
        for i in range(1, D + 1):
            for j in range(i, D + 1):
                if any((C2(dm.count(i), 2) if i == j else dm.count(i) * dm.count(j)) for dm in self.dms):
                    row(f"wdeg: witness pairs ({i},{j}) <= M{i}{j}", '<', 0); self.wpair_rows.add((i, j))
        # ---- F6 flags (Lemma B, R1-R4) ----
        prof_of = defaultdict(set)
        for t in T:
            if t[0] >= RMIN and t[5]:
                for sp in t[5]: prof_of[(t[0], t[1])].add(sp)
        kinds = set(); self.gco = defaultdict(dict)
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
        sps = {kd[2] for kd in kinds}; X2 = []
        for i, w in enumerate(wt):
            for pos in {w[p]: p for p in range(len(w))}.values():
                r2, lo, hi, kk = items[w[pos]]
                if r2 < 1: continue
                rest = tuple(sorted(w[:pos] + w[pos + 1:]))
                for sp in sps:
                    if matchable(rest, sp):
                        for k2 in kk:
                            if (sp, (r2, k2)) in srk2set: X2.append(((sp, (r2, k2)), i))
        self.X2 = list(dict.fromkeys(X2))
        self.types_rks = defaultdict(list)
        for j, t in enumerate(T):
            if t[0] >= RMIN and t[5]:
                for sp in set(t[5]): self.types_rks[(t[0], t[1], sp)].append(j)
        for key in self.types_rks: row(f"F6 R1 flags (r={key[0]},k={key[1]},s={key[2]})", '=', 0)
        for (rk, rk2) in self.gco: row(f"F6 R2 flags ({rk} hosted by {rk2})", '=', 0)
        by_srk2 = {(kd[2], kd[3]) for kd in kinds} | {k for k, i in self.X2}
        assert by_srk2 <= srk2set
        for key in sorted(by_srk2): row(f"F6 link flags (s={key[0]},mate={key[1]})", '=', 0)
        # F9 (Lemma C on the mate plane pi' of the flag, WRITEUP8 §1): s(w) >= 25 + 5 r' - D (k' - 5); positive deficits
        # force the flag column to 0. One row per (profile, mate) key.
        self.deficit = {}
        for (key, i) in self.X2:
            d = max(0, base + q * key[1][0] - D * (key[1][1] - q) - self.s_of[i])
            if d: self.deficit[(key, i)] = d
        for key in self.srk2: row(f"F9 unmarked C (s={key[0]},mate={key[1]})", '<', 0)
        by_w = defaultdict(list)
        for (key, i) in self.X2: by_w[i].append(key)
        self.cap3cap = {}; self.xmax = {}; full3 = (3,) * 5
        for i in by_w:
            cap_ = sum((6 - items[ii][0]) // RMIN for ii in wt[i] if items[ii][0] >= 1); self.cap3cap[i] = cap_
            row(f"F6 R3 flags on pencil {i} <= {cap_} W{i}", '<', 0)
        for i in by_w:
            if not any(key[0] == full3 for key in by_w[i]): continue
            tot = sum(vi.cap3(items[ii]) for ii in wt[i]); xmax = 0
            while xmax + 1 <= 6 and 5 * (xmax + 1) - C2(xmax + 1, 2) <= tot: xmax += 1
            self.xmax[i] = xmax; row(f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {xmax} W{i}", '<', 0)
        # ---- F8 (Lemma D(a)): G[(r,k,h), dm] for r >= 1 and dm subset of h (as multisets); Z[dm, w] for sum(dm) = s(w) and every
        # item of w realisable with a histogram containing dm. Rows: Y1 (each plane of key (r,k,h) carries r inside witness
        # lines), Y2 (incidences counted plane-side = pencil-side, per item and dm), Z1, Z2 (marginals to W and V). ----
        rkh = {}
        for i, t in enumerate(T):
            if t[0] >= 1: rkh.setdefault((t[0], t[1], tuple(t[2])), []).append(i)
        self.rkh = rkh; self.key_of_type = {}
        for key, idxs in rkh.items():
            for i in idxs: self.key_of_type[i] = key
        dmv = self.dmv = [tuple(dm.count(j) for j in range(1, D + 1)) for dm in self.dms]
        self.G = [(key, di) for key in rkh for di, dv in enumerate(dmv) if all(key[2][j] >= dv[j] for j in range(D))]
        rk_dm_ok = {(key[0], key[1], di) for key, di in self.G}
        self.Z = []
        for i, t in enumerate(wt):
            for di, dm in enumerate(self.dms):
                if sum(dm) != self.s_of[i]: continue
                if all(any((items[ii][0], k, di) in rk_dm_ok for k in items[ii][3]) for ii in set(t)): self.Z.append((di, i))
        for key in rkh: row(f"F8 Y1 {key}", '=', 0)
        for i in range(len(wt)): row(f"F8 Z1 w{i}", '=', 0)
        for di in range(len(self.dms)): row(f"F8 Z2 dm{self.dms[di]}", '=', 0)
        y2 = set()
        for key, di in self.G: y2.add((self.item_of[(key[0], key[1])], di))
        for di, i in self.Z:
            for ii in set(wt[i]): y2.add((ii, di))
        self.y2 = y2
        for ii in range(len(items)):
            for di in range(len(self.dms)):
                if (ii, di) in y2: row(f"F8 Y2 item{ii} dm{self.dms[di]}", '=', 0)
        # ---- lifts on histogram keys (r,k,h) ----
        self.groups = {}       # CORE keys (r,k) in {(5,15),(6,17)}
        self.groups_all = {}   # every key with r >= 5 (D = 7 'all' scope: non-core keys get the 2D coordinate objectives only)
        for i, t in enumerate(T):
            if t[0] >= 5:
                a = (t[0], t[1], tuple(t[2])); self.groups_all.setdefault(a, []).append(i)
                if (t[0], t[1]) in CORE: self.groups.setdefault(a, []).append(i)
        self.priceXI = {}; self.priceXIII = {}; self.priceSub = {}; self.priceNL = {}; self.priceTr = {}; self.pricePt = {}; self.pricePt3 = {}
        self.subrows = {}; self.nlrows = {}; self.trrows = {}   # key -> list of (objective index, (S, sign)) rows actually present (compact objectives)
        if not self.xiv:    # D = 5: coordinate lift (WRITEUP11) and, for (6,17), pairwise lift (WRITEUP13)
            for a, inds in self.groups.items():
                for (S, s) in OBJ_XI: row(f"XI incidence {jkey(a)} {S[0]} {s}", '<', 0)
                for i in inds: self.priceXI[i] = [price(T[i], S, s) for (S, s) in OBJ_XI]
            if self.pair:
                for a, inds in self.groups.items():
                    for (S, s) in OBJ_XIII: row(f"XIII pair {jkey(a)} {S[0]} {S[1]} {s}", '<', 0)
                    for i in inds: self.priceXIII[i] = [price(T[i], S, s) for (S, s) in OBJ_XIII]
        else:               # D >= 6: additive subset lift, within-line pair (and D = 7: triple) costs, WRITEUP14
            self.OBJ_SUB = obj_subsets(D); self.OBJ_NL = obj_withinline(D); compact = F.get('compact', False)
            def select(a, objs, scope):
                """scope 'pair': subsets of size <= 2 only. compact objectives: drop an objective whose subset misses the support of
                h (cost identically 0) and repeat signatures (S ∩ support, sign) (identical cost function on every column of the key)."""
                objs = [(z, o) for z, o in enumerate(objs) if (scope != 'pair' or len(o[0]) <= 2) and (scope != 'coordinate' or len(o[0]) <= 1)]
                if not compact: return objs
                support = {j for j, hj in enumerate(a[2], 1) if hj}; seen = set(); out = []
                for z, (S, sg) in objs:
                    eff = tuple(d for d in S if d in support)
                    if not eff or (eff, sg) in seen: continue
                    seen.add((eff, sg)); out.append((z, (S, sg)))
                return out
            cb = tuple(b for b in range(8))
            for a, inds in (self.groups_all if F['sub'] == 'all' else self.groups).items():
                objs = self.OBJ_SUB if (a[0], a[1]) in CORE else self.OBJ_SUB[:2 * D]
                sel = select(a, objs, F['sub']); self.subrows[a] = sel
                for z, (S, sg) in sel: row(f"XIV subset {jkey(a)} {jkey(S)} {sg}", '<', 0)
                for i in inds: self.priceSub[i] = {z: (price_alloc(T[i], S, sg, cb) if D == 7 else price(T[i], S, sg)) for z, (S, sg) in sel}
            for a, inds in (self.groups.items() if F.get('nl') else []):
                sel = select(a, self.OBJ_NL, F['nl']); self.nlrows[a] = sel
                for z, (S, sg) in sel: row(f"XIV nonlinear {jkey(a)} {jkey(S)} {sg}", '<', 0)
                for i in inds:
                    t = T[i]; key = (t[0], t[1], t[3], t[4], t[5]); assert key in geo_masks, ('no geometry for type', t)
                    pv = {z: price_withinline(key, tuple(t[2]), S, sg, 2) for z, (S, sg) in sel}
                    assert all(v is not None for v in pv.values()), ('within-line price missing', t); self.priceNL[i] = pv
            if F.get('tr'):
                for a, inds in self.groups.items():
                    sel = select(a, self.OBJ_NL, F['tr']); self.trrows[a] = sel
                    for z, (S, sg) in sel: row(f"XIV triple {jkey(a)} {jkey(S)} {sg}", '<', 0)
                    for i in inds:
                        t = T[i]; key = (t[0], t[1], t[3], t[4], t[5])
                        pv = {z: price_withinline(key, tuple(t[2]), S, sg, 3) for z, (S, sg) in sel}
                        assert all(v is not None for v in pv.values()), ('triple price missing', t); self.priceTr[i] = pv
            if F.get('point'):
                # point families (AUDIT §3.15): sum_{pi ∋ p} C(b_pi(p),2) = C(d_p,2); sum_pi C(b_pi(p),3) = T(p) in [lo_j, hi_j] (PG(2,5) bounds)
                self.OBJ_PT = [(j, sg) for j in range(2, D + 1) for sg in (1, -1)]; self.OBJ_PT3 = [(j, sg) for j in range(3, D + 1) for sg in (1, -1)]
                pj = json.load(open('out/projective_result.json'))   # my own exhaustive PG(2,5) enumeration (check_projective.py)
                self.tri_bounds = {j: tuple(pj[str(j)]['T']) for j in range(3, D + 1)}
                assert self.tri_bounds[7] == (2, 20) and all(self.tri_bounds[j] == (0, C2(j, 3)) for j in range(3, 7))
                for (j, sg) in self.OBJ_PT: row(f"XIV point crossing degree {j} sign {sg}", '<', 0)
                for (j, sg) in self.OBJ_PT3: row(f"XIV point coplanar triples degree {j} sign {sg}", '<', 0)
                c2 = tuple(C2(b, 2) for b in range(8)); c3 = tuple(C2(b, 3) for b in range(8))
                for i, t in enumerate(T):
                    if t[0] >= 2: self.pricePt[i] = [price_alloc(t, (j,), sg, c2) for (j, sg) in self.OBJ_PT]
                    if t[0] >= 3: self.pricePt3[i] = [price_alloc(t, (j,), sg, c3) for (j, sg) in self.OBJ_PT3]
            if F.get('star'):       # AUDIT §3.15: b <= 4 in the retained cases, T(p) <= star bound 4,5,8 for j = 5,6,7 (my enumeration)
                assert F.get('point') and F.get('maximal') and F.get('budget')
                self.star = {j: pj[str(j)]['star_bound_cap4'] for j in range(5, D + 1)}; assert [self.star[j] for j in (5, 6, 7)] == [4, 5, 8]
                for j in range(5, D + 1): row(f"XIV star-four triples degree {j}", '<', 0)
            if F.get('nonfano'):    # AUDIT §3.15: sum_pi (C(b,3) - C(b,4)) at a point = g(p) in [lo_j, hi_j] (my enumeration)
                self.OBJ_NF = [(j, sg) for j in range(3, D + 1) for sg in (1, -1)]; self.nf_bounds = {j: tuple(pj[str(j)]['g']) for j in range(3, D + 1)}
                assert self.nf_bounds == {3: (0, 1), 4: (0, 3), 5: (0, 5), 6: (0, 5), 7: (2, 6)}
                for (j, sg) in self.OBJ_NF: row(f"XIV point triples-minus-quadruples degree {j} sign {sg}", '<', 0)
                c34 = tuple(C2(b, 3) - C2(b, 4) for b in range(8)); self.priceNF = {}
                for i, t in enumerate(T):
                    if t[0] >= 3: self.priceNF[i] = [price_alloc(t, (j,), sg, c34) for (j, sg) in self.OBJ_NF]
            if F.get('rbudget'):    # AUDIT §3.16: sum_pi floor(c_pi / t) <= floor(541 / t); one row per distinct floor(541/t)
                assert F.get('budget'); self.dens = []; last = None
                for t_ in range(2, 542):
                    if 541 // t_ != last: self.dens.append(t_); last = 541 // t_
                for t_ in self.dens: row(f"XIV rounded plane budget denominator {t_}", '<', 541 // t_)
            if F.get('rich'):       # AUDIT §3.15: #planes through p with b >= 3 = #rich secants (>= 3 selected points) of p's direction set, in [lo, hi] under cap 4
                assert F.get('star'); self.rich_bounds = {j: tuple(pj[str(j)]['rich_cap4']) for j in range(3, D + 1)}
                self.OBJ_RICH = [(j, sg) for j in range(3, D + 1) for sg in (1, -1)]; crich = tuple(int(b >= 3) for b in range(8)); self.priceRich = {}
                for (j, sg) in self.OBJ_RICH: row(f"XIV point rich-secants degree {j} sign {sg}", '<', 0)
                for i, t in enumerate(T):
                    if t[0] >= 3: self.priceRich[i] = [price_alloc(t, (j,), sg, crich) for (j, sg) in self.OBJ_RICH]
            if F.get('mix'):        # joint objective lam * C(b,2) + [b >= 3]: lam C(j,2) + rich(p) in [lo + lam C(j,2), hi + lam C(j,2)]
                assert F.get('rich'); self.LAM = (-2, -1, 1, 2); self.OBJ_MIX = [(j, lam, sg) for j in range(3, D + 1) for lam in self.LAM for sg in (1, -1)]; self.priceMix = {}
                for (j, lam, sg) in self.OBJ_MIX: row(f"XIV point mixed degree {j} lambda {lam} sign {sg}", '<', 0)
                cm = {lam: tuple(lam * C2(b, 2) + int(b >= 3) for b in range(8)) for lam in self.LAM}
                for i, t in enumerate(T):
                    if t[0] >= 2: self.priceMix[i] = [price_alloc(t, (j,), sg, cm[lam]) for (j, lam, sg) in self.OBJ_MIX]
            if F.get('psubset'):    # crossing pairs summed over a SET of degrees: sum_{d in S} C(d,2) m_d
                degs = tuple(range(2, D + 1)); self.PSUB = [S for sz in range(2, len(degs) // 2 + 1) for S in itertools.combinations(degs, sz) if 2 * sz < len(degs) or S < tuple(d for d in degs if d not in S)]
                self.OBJ_PS = [(S, sg) for S in self.PSUB for sg in (1, -1)]; self.pricePS = {}; c2 = tuple(C2(b, 2) for b in range(8))
                for (S, sg) in self.OBJ_PS: row(f"XIV point subset crossing degrees {S} sign {sg}", '<', 0)
                for i, t in enumerate(T):
                    if t[0] >= 2: self.pricePS[i] = [price_alloc(t, S, sg, c2) for (S, sg) in self.OBJ_PS]
            if F.get('extra'):      # within-line pair/triple costs on (6, k > 17) planes (u = 17 + extra points), WRITEUP14 extra14
                assert F.get('maximal') and cfg == '617'
                self.OBJ_EX = [(S, sg, pw) for sz in (1, 2) for S in itertools.combinations(range(1, D + 1), sz) for pw in (2, 3) for sg in (1, -1)]
                self.groups_ex = {}; self.exrows = {}; self.priceEx = {}
                for i, t in enumerate(T):
                    if t[0] == 6 and t[1] > 17: self.groups_ex.setdefault((t[0], t[1], tuple(t[2])), []).append(i)
                for a, inds in self.groups_ex.items():
                    support = {j for j, hj in enumerate(a[2], 1) if hj}; seen = set(); sel = []
                    for z, (S, sg, pw) in enumerate(self.OBJ_EX):
                        eff = tuple(d for d in S if d in support)
                        if not eff or (compact and (eff, sg, pw) in seen): continue
                        seen.add((eff, sg, pw)); sel.append((z, (S, sg, pw)))
                    self.exrows[a] = sel
                    for z, (S, sg, pw) in sel: row(f"XIV extra {jkey(a)} {jkey(S)} {sg} {pw}", '<', 0)
                    for i in inds:
                        t = T[i]; key = (t[0], t[1], t[3], t[4], t[5]); assert key in geo_masks, ('no geometry', t)
                        pv = {z: price_withinline(key, tuple(t[2]), S, sg, pw) for z, (S, sg, pw) in sel}
                        assert all(v is not None for v in pv.values()), ('extra price missing', t); self.priceEx[i] = pv
            if F.get('rdeg'):       # AUDIT §3.16 parallel-class rounding: per degree subset S and t, sum_classes sum_i floor(h_i(S)/t) in [31*lower(M,t), 31*floor(M/t)], M = sum_{j in S} m_j
                assert F.get('hull'); self.RDS = sorted({(j,) for j in range(1, D + 1)} | {tuple(range(1, j + 1)) for j in range(2, D)})
                self.RD = []; self.rd_M_of = {}
                for S in self.RDS:
                    Ms = sorted({sum(v[j - 1] for j in S) for v in cand}); self.rd_M_of[S] = Ms
                    row(f"XIV rounded degree normalization {jkey(S)}", '=', 1)
                    for M in Ms: row(f"XIV rounded degree marginal {jkey(S)} {M}", '=', 0); self.RD.append((S, M))
                    for t_ in range(2, 26): row(f"XIV rounded degree upper {jkey(S)} {t_}", '<', 0); row(f"XIV rounded degree lower {jkey(S)} {t_}", '<', 0)
            if F.get('domain_m1'): row("XIV domain m1 lower", '<', -F['domain_m1'])
            if F.get('pairb'):
                # pair bounds (AUDIT §3.14): box L <= m <= U from the raw vectors, McCormick products, diagonal chord
                for j in range(1, D + 1): row(f"PB14 upper {j}", '<', cap[j - 1]); row(f"PB14 lower {j}", '<', -self.floor[j - 1]); row(f"PB14 diagonal chord {j}", '<', -self.floor[j - 1] * cap[j - 1])
                for i, j in itertools.combinations(range(1, D + 1), 2):
                    Li, Ui, Lj, Uj = self.floor[i - 1], cap[i - 1], self.floor[j - 1], cap[j - 1]
                    row(f"PB14 product lower LL {i},{j}", '<', Li * Lj); row(f"PB14 product lower UU {i},{j}", '<', Ui * Uj)
                    row(f"PB14 product upper UL {i},{j}", '<', -Ui * Lj); row(f"PB14 product upper LU {i},{j}", '<', -Li * Uj)
            # degree hull (AUDIT §3.14): Q_v >= 0 per raw vector, sum Q = 1, m = sum v Q, M = sum moment(v) Q
            row("XIV degree hull normalization", '=', 1)
            for j in range(1, D + 1): row(f"XIV degree hull m{j}", '=', 0)
            for i in range(1, D + 1):
                for j in range(i, D + 1): row(f"XIV degree hull M{i},{j}", '=', 0)
        # ---- node rows ----
        self.pfx = 'XIV' if self.xiv else 'XI'
        row(f"{self.pfx} config{'uration' if self.xiv else ''}", '<', -1); row(f"{self.pfx} rmax", '<', 0)
        for j in range(1, D + 1):
            lo = hi = prefix[j - 1] if j <= len(prefix) else None
            if j > len(prefix): lo, hi = 0, n
            row(f"{self.pfx} upper m{j}", '<', hi); row(f"{self.pfx} lower m{j}", '<', -lo)
        nQ = len(cand) if F.get('hull') else 0; nRD = len(self.RD) if F.get('rdeg') else 0
        self.info = dict(D=D, types=len(T), classes=len(ctypes), dirtypes=len(self.utypes), linevecs=len(wvecs), pencils=len(wt), items=len(items),
                         X1=len(self.kinds), X2=len(self.X2), G=len(self.G), Z=len(self.Z), Q=nQ, lift_keys=len(self.groups), lift_types=sum(map(len, self.groups.values())),
                         columns=len(T) + D + 5 + D * (D + 1) // 2 + len(ctypes) + len(self.utypes) + len(wvecs) + len(wt) + len(self.dms) + len(self.kinds) + len(self.X2) + len(self.G) + len(self.Z) + nQ + nRD, RD=nRD,
                         rows=len(self.rows), cap=cap, build_s=round(time.time() - t0, 1))

    def columns(self):
        """yield (column id, [(row name, coefficient), ...]) for every column of the universe (zero coefficients allowed; dropped by verify)."""
        T, items, wt, dms = self.T, self.items, self.wt, self.dms; full3 = (3,) * 5; pfx = self.pfx; xiv = self.xiv; F = self.F
        for i, t in enumerate(T):
            r, k, h, u, tt, e = t
            c = [("I2a: sum_pi 1", 1), ("I2b: sum_pi r", r), ("I2c: sum_pi k", k), ("I3: sum C(r,2) = sum C(d,2)", C2(r, 2)), ("I4: sum k^2", k * k),
                 ("I5: sum k r", k * r), ("I6e: sum C(k,3) - q sum C(j,3) L_j", C2(k, 3)), (f"class link N{(r, k)}", -1)]
            c += [(f"F5: sum_pi t{j} = (q+1) L{j}", tt[j - 1]) for j in range(1, 6)]
            for j in range(1, D + 1):
                hj = h[j - 1]
                if hj:
                    c += [(f"H1: sum_pi h{j} = nd m{j}", hj), (f"H2: sum_pi h{j} k = ((q+1)n+q^2-q) m{j}", hj * k), (f"H3: sum_pi h{j} r = (nd+q j) m{j}", hj * r),
                          (f"H4: sum_pi C(h{j},2) = (q+1) M{j}{j}", C2(hj, 2))]
                    c += [(f"H5: sum_pi h{i2} h{j} = (q+1) M{i2}{j}", h[i2 - 1] * hj) for i2 in range(1, j)]
            if r >= 1:
                ii = self.item_of[(r, k)]; c.append((f"wline link r*N[{r},{items[ii][1]}..{items[ii][2]}]", -r))
                key = self.key_of_type[i]; c.append((f"F8 Y1 {key}", -r))
            if r >= RMIN and e:
                for sp in set(e): c.append((f"F6 R1 flags (r={r},k={k},s={sp})", -e.count(sp)))
            a = (r, k, tuple(h))
            if i in self.priceXI:
                c += [(f"XI incidence {jkey(a)} {S[0]} {s}", p) for (S, s), p in zip(OBJ_XI, self.priceXI[i])]
                if self.pair: c += [(f"XIII pair {jkey(a)} {S[0]} {S[1]} {s}", p) for (S, s), p in zip(OBJ_XIII, self.priceXIII[i])]
            if i in self.priceSub:
                objs = self.OBJ_SUB if (r, k) in CORE else self.OBJ_SUB[:2 * D]
                c += [(f"XIV subset {jkey(a)} {jkey(objs[z][0])} {objs[z][1]}", p) for z, p in self.priceSub[i].items()]
            if i in self.priceNL: c += [(f"XIV nonlinear {jkey(a)} {jkey(self.OBJ_NL[z][0])} {self.OBJ_NL[z][1]}", p) for z, p in self.priceNL[i].items()]
            if i in self.priceTr: c += [(f"XIV triple {jkey(a)} {jkey(self.OBJ_NL[z][0])} {self.OBJ_NL[z][1]}", p) for z, p in self.priceTr[i].items()]
            if F.get('point'):
                if i in self.pricePt: c += [(f"XIV point crossing degree {j} sign {sg}", p) for (j, sg), p in zip(self.OBJ_PT, self.pricePt[i])]
                if i in self.pricePt3: c += [(f"XIV point coplanar triples degree {j} sign {sg}", p) for (j, sg), p in zip(self.OBJ_PT3, self.pricePt3[i])]
            if F.get('star') and i in self.pricePt3: c += [(f"XIV star-four triples degree {j}", self.pricePt3[i][self.OBJ_PT3.index((j, 1))]) for j in range(5, D + 1)]
            if F.get('nonfano') and i in self.priceNF: c += [(f"XIV point triples-minus-quadruples degree {j} sign {sg}", p) for (j, sg), p in zip(self.OBJ_NF, self.priceNF[i])]
            if F.get('rbudget'):
                cost = plane_cost(r, k); c += [(f"XIV rounded plane budget denominator {t_}", cost // t_) for t_ in self.dens if cost >= t_]
            if F.get('rich') and i in self.priceRich: c += [(f"XIV point rich-secants degree {j} sign {sg}", p) for (j, sg), p in zip(self.OBJ_RICH, self.priceRich[i])]
            if F.get('mix') and i in self.priceMix: c += [(f"XIV point mixed degree {j} lambda {lam} sign {sg}", p) for (j, lam, sg), p in zip(self.OBJ_MIX, self.priceMix[i])]
            if F.get('psubset') and i in self.pricePS: c += [(f"XIV point subset crossing degrees {S} sign {sg}", p) for (S, sg), p in zip(self.OBJ_PS, self.pricePS[i])]
            if F.get('extra') and i in self.priceEx: c += [(f"XIV extra {jkey(a)} {jkey(self.OBJ_EX[z][0])} {self.OBJ_EX[z][1]} {self.OBJ_EX[z][2]}", p) for z, p in self.priceEx[i].items()]
            if F.get('rdeg'):
                for S in self.RDS:
                    hS = sum(h[j - 1] for j in S)
                    for t_ in range(2, 26):
                        if hS >= t_: c += [(f"XIV rounded degree upper {jkey(S)} {t_}", hS // t_), (f"XIV rounded degree lower {jkey(S)} {t_}", -(hS // t_))]
            if r == self.rr and u == self.uu: c.append((f"{pfx} config{'uration' if xiv else ''}", -1))
            if r > self.rr: c.append((f"{pfx} rmax", 1))
            yield ('N', t), c
        for j in range(1, D + 1):
            c = [("deg: sum m_j", 1), ("deg: sum j m_j", j), ("I3: sum C(r,2) = sum C(d,2)", -C2(j, 2)), (f"H1: sum_pi h{j} = nd m{j}", -nd),
                 (f"H2: sum_pi h{j} k = ((q+1)n+q^2-q) m{j}", -(6 * n + 25)), (f"H3: sum_pi h{j} r = (nd+q j) m{j}", -(nd + 5 * j)),
                 (f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}", -(n - 1)), (f"H7: 2 M{j}{j} <= (n-1) m{j}", -(n - 1)),
                 ("wlines: sum s_l = sum d_p^2", -j * j), (f"wdeg: sum_l #(deg {j} on l) = {j} m{j}", -j), (f"{pfx} upper m{j}", 1), (f"{pfx} lower m{j}", -1)]
            if j == D: c.append(("m_D >= 1", -1))
            c += [(f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)", a) for a in range(1, n)]
            if xiv: c.append((f"XIV degree hull m{j}", 1))
            if F.get('point'):
                c += [(f"XIV point crossing degree {j} sign {sg}", -sg * C2(j, 2)) for sg in (1, -1) if j >= 2]
                c += [(f"XIV point coplanar triples degree {j} sign {sg}", -self.tri_bounds[j][1] if sg == 1 else self.tri_bounds[j][0]) for sg in (1, -1) if j >= 3]
            if F.get('star') and j >= 5: c.append((f"XIV star-four triples degree {j}", -self.star[j]))
            if F.get('rich') and j >= 3: c += [(f"XIV point rich-secants degree {j} sign {sg}", -self.rich_bounds[j][1] if sg == 1 else self.rich_bounds[j][0]) for sg in (1, -1)]
            if F.get('mix') and j >= 3:
                for lam in self.LAM:
                    lo, hi = self.rich_bounds[j][0] + lam * C2(j, 2), self.rich_bounds[j][1] + lam * C2(j, 2)
                    c += [(f"XIV point mixed degree {j} lambda {lam} sign {sg}", -hi if sg == 1 else lo) for sg in (1, -1)]
            if F.get('psubset') and j >= 2: c += [(f"XIV point subset crossing degrees {S} sign {sg}", -sg * C2(j, 2)) for (S, sg) in self.OBJ_PS if j in S]
            if F.get('domain_m1') and j == 1: c.append(("XIV domain m1 lower", -1))
            if F.get('nonfano') and j >= 3: c += [(f"XIV point triples-minus-quadruples degree {j} sign {sg}", -self.nf_bounds[j][1] if sg == 1 else self.nf_bounds[j][0]) for sg in (1, -1)]
            if F.get('pairb'):
                c += [(f"PB14 upper {j}", 1), (f"PB14 lower {j}", -1), (f"PB14 diagonal chord {j}", -(self.floor[j - 1] + self.cap[j - 1] - 1))]
                for i2 in range(1, D + 1):
                    if i2 == j: continue
                    a_, b_ = min(i2, j), max(i2, j); Lo, Uo = self.floor[i2 - 1], self.cap[i2 - 1]   # bounds of the OTHER index
                    c += [(f"PB14 product lower LL {a_},{b_}", Lo), (f"PB14 product lower UU {a_},{b_}", Uo), (f"PB14 product upper UL {a_},{b_}", -Uo if i2 < j else -Lo), (f"PB14 product upper LU {a_},{b_}", -Lo if i2 < j else -Uo)]
            yield ('m', j), c
        for j in range(1, 6):
            c = [("I6a: sum j L_j", j), ("I6b: sum C(j,2) L_j", C2(j, 2)), ("I6c: sum L_j >= nd*BM", -1), ("I6e: sum C(k,3) - q sum C(j,3) L_j", -5 * C2(j, 3)),
                 (f"F5: sum_pi t{j} = (q+1) L{j}", -6), (f"dir link L{j}", -1)]
            if j == 5: c.append(("I6d: L_q >= nd", -1))
            yield ('L', j), c
        for i in range(1, D + 1):
            for j in range(i, D + 1):
                if i == j:
                    c = [(f"H4: sum_pi C(h{j},2) = (q+1) M{j}{j}", -6), (f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}", 2), (f"H7: 2 M{j}{j} <= (n-1) m{j}", 2)]
                    c += [(f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)", -1) for a in range(1, n)]
                else:
                    c = [(f"H5: sum_pi h{i} h{j} = (q+1) M{i}{j}", -6), (f"H6: sum_i M_i{i} + 2 M{i}{i} = (n-1) m{i}", 1), (f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}", 1)]
                if (i, j) in self.wpair_rows: c.append((f"wdeg: witness pairs ({i},{j}) <= M{i}{j}", -1))
                if xiv: c.append((f"XIV degree hull M{i},{j}", 1))
                if F.get('pairb'):
                    if i == j: c.append((f"PB14 diagonal chord {j}", 2))
                    else: c += [(f"PB14 product lower LL {i},{j}", -1), (f"PB14 product lower UU {i},{j}", -1), (f"PB14 product upper UL {i},{j}", 1), (f"PB14 product upper LU {i},{j}", 1)]
                yield ('M', i, j), c
        def phi(t, v):
            rho = sum(1 for r, k in t if r >= 1); n1 = sum(1 for r, k in t if r == 1); return (6 - n1) * (rho == v) + n1 * (rho - 1 == v)
        r2_by_class = defaultdict(list)
        for key, d in self.gco.items():
            for ci, g in d.items(): r2_by_class[ci].append((key, g))
        for ci, ct in enumerate(self.ctypes):
            c = [("classes: sum C = nd", 1)]
            cnt = defaultdict(int)
            for rk in ct: cnt[rk] += 1
            c += [(f"class link N{rk}", v) for rk, v in cnt.items()]
            c += [(f"dir flags a=={v}", -phi(ct, v)) for v in range(1, 6) if phi(ct, v)]
            c += [(f"F6 R2 flags ({key[0]} hosted by {key[1]})", g) for key, g in r2_by_class[ci]]
            yield ('C', ct), c
        for u in self.utypes:
            c = [("dirs: sum A = nd", 1)] + [(f"dir flags a=={v}", u.count(v)) for v in range(1, 6) if u.count(v)]
            c += [(f"dominance lambda>={th}", 1) for th in range(bm + 1, 26) if Ft[u] >= th]
            yield ('A', u), c
        for w in self.wvecs:
            c = [("dirs: sum P = nd", 1)] + [(f"dir link L{j}", w[j - 1]) for j in range(1, 6) if w[j - 1]]
            c += [(f"dominance lambda>={th}", -1) for th in range(bm + 1, 26) if sum(w) >= th]
            yield ('P', w), c
        for i, t in enumerate(wt):
            c = [("wlines: sum W = nd", 1), ("wlines: sum s_l = sum d_p^2", self.s_of[i]), (f"wdeg/wline link s={self.s_of[i]}", -1), (f"F8 Z1 w{i}", -1)]
            cnt = defaultdict(int)
            for ii in t: cnt[ii] += 1
            c += [(f"wline link r*N[{items[ii][0]},{items[ii][1]}..{items[ii][2]}]", v) for ii, v in cnt.items()]
            if i in self.cap3cap: c.append((f"F6 R3 flags on pencil {i} <= {self.cap3cap[i]} W{i}", -self.cap3cap[i]))
            if i in self.xmax: c.append((f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {self.xmax[i]} W{i}", -self.xmax[i]))
            yield ('W', i), c
        for di, dm in enumerate(dms):
            c = [("wdeg: sum V = nd", 1), (f"wdeg/wline link s={sum(dm)}", 1), (f"F8 Z2 dm{dm}", -1)]
            c += [(f"wdeg: sum_l #(deg {j} on l) = {j} m{j}", dm.count(j)) for j in sorted(set(dm))]
            for i in range(1, D + 1):
                for j in range(i, D + 1):
                    v = C2(dm.count(i), 2) if i == j else dm.count(i) * dm.count(j)
                    if v: c.append((f"wdeg: witness pairs ({i},{j}) <= M{i}{j}", v))
            yield ('V', dm), c
        for kd in self.kinds:
            yield ('X1', kd), [(f"F6 R1 flags (r={kd[0]},k={kd[1]},s={kd[2]})", 1), (f"F6 R2 flags ({(kd[0], kd[1])} hosted by {kd[3]})", 1), (f"F6 link flags (s={kd[2]},mate={kd[3]})", 1)]
        for (key, i) in self.X2:
            c = [(f"F6 link flags (s={key[0]},mate={key[1]})", -1), (f"F6 R3 flags on pencil {i} <= {self.cap3cap[i]} W{i}", 1)]
            if key[0] == full3 and i in self.xmax: c.append((f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {self.xmax[i]} W{i}", 1))
            if (key, i) in self.deficit: c.append((f"F9 unmarked C (s={key[0]},mate={key[1]})", self.deficit[(key, i)]))
            yield ('X2', key, i), c
        for (key, di) in self.G:
            dm = dms[di]; ii = self.item_of[(key[0], key[1])]
            c = [(f"F8 Y1 {key}", 1), (f"F8 Y2 item{ii} dm{dm}", -1)]
            if key in self.groups and not xiv:
                if True:
                    c += [(f"XI incidence {jkey(key)} {S[0]} {s}", -s * dm.count(S[0])) for (S, s) in OBJ_XI if dm.count(S[0])]
                    if self.pair: c += [(f"XIII pair {jkey(key)} {S[0]} {S[1]} {s}", -s * (dm.count(S[0]) + dm.count(S[1]))) for (S, s) in OBJ_XIII if dm.count(S[0]) + dm.count(S[1])]
            if xiv and key in self.subrows:
                c += [(f"XIV subset {jkey(key)} {jkey(S)} {s}", -s * sum(d in S for d in dm)) for z, (S, s) in self.subrows[key] if any(d in S for d in dm)]
                if key in self.nlrows: c += [(f"XIV nonlinear {jkey(key)} {jkey(S)} {s}", -s * C2(sum(d in S for d in dm), 2)) for z, (S, s) in self.nlrows[key] if sum(d in S for d in dm) >= 2]
                if key in self.trrows: c += [(f"XIV triple {jkey(key)} {jkey(S)} {s}", -s * C2(sum(d in S for d in dm), 3)) for z, (S, s) in self.trrows[key] if sum(d in S for d in dm) >= 3]
            if F.get('extra') and key in self.exrows:
                c += [(f"XIV extra {jkey(key)} {jkey(S)} {s} {pw}", -s * C2(sum(d in S for d in dm), pw)) for z, (S, s, pw) in self.exrows[key] if sum(d in S for d in dm) >= pw]
            yield ('G', key, di), c
        for (di, i) in self.Z:
            dm = dms[di]; cnt = defaultdict(int)
            for ii in wt[i]: cnt[ii] += 1
            yield ('Z', di, i), [(f"F8 Z1 w{i}", 1), (f"F8 Z2 dm{dm}", 1)] + [(f"F8 Y2 item{ii} dm{dm}", v) for ii, v in cnt.items()]
        if F.get('hull'):
            for v in self.cand:
                c = [("XIV degree hull normalization", 1)] + [(f"XIV degree hull m{j}", -v[j - 1]) for j in range(1, D + 1) if v[j - 1]]
                for i in range(1, D + 1):
                    for j in range(i, D + 1):
                        mo = C2(v[i - 1], 2) if i == j else v[i - 1] * v[j - 1]
                        if mo: c.append((f"XIV degree hull M{i},{j}", -mo))
                if F.get('rdeg'): c += [(f"XIV rounded degree marginal {jkey(S)} {sum(v[j - 1] for j in S)}", -1) for S in self.RDS]
                yield ('Q', v), c
        if F.get('rdeg'):
            def rd_lower(M, t_): return max(0, -(-(M - 5 * (t_ - 1)) // t_))
            for (S, M) in self.RD:
                c = [(f"XIV rounded degree normalization {jkey(S)}", 1), (f"XIV rounded degree marginal {jkey(S)} {M}", 1)]
                for t_ in range(2, 26):
                    c += [(f"XIV rounded degree upper {jkey(S)} {t_}", -31 * (M // t_)), (f"XIV rounded degree lower {jkey(S)} {t_}", 31 * rd_lower(M, t_))]
                yield ('RD', S, M), c

    def verify(self, y_by_name, research_nnz=None):
        rows = self.rows; y = {k: v for k, v in y_by_name.items() if v}
        neg = [k for k, v in y.items() if rows[k][0] == '<' and v < 0]
        rhs = sum(v * rows[k][1] for k, v in y.items())
        nnz = defaultdict(int); mincol = None; ncol = 0; bad = []; fam = defaultdict(int); zero_cols = 0
        for cid, entries in self.columns():
            ncol += 1; fam[cid[0]] += 1; v = 0
            for name, coef in entries:
                if not coef: continue
                assert name in rows, (cid, name); nnz[name] += 1; v += coef * y.get(name, 0)
            if mincol is None or v < mincol: mincol = v
            if v == 0: zero_cols += 1
            if v < 0 and len(bad) < 5: bad.append((cid, v))
        struct = None
        if research_nnz is not None:
            diff = [(k, nnz.get(k, 0), research_nnz[k]) for k in rows if nnz.get(k, 0) != research_nnz.get(k, 0)]
            struct = dict(rows_with_different_support=len(diff), examples=diff[:5])
        ok = not bad and not neg and rhs < 0
        return ok, dict(rhs=rhs, min_column=mincol, columns=ncol, columns_by_family=dict(fam), tight_columns=zero_cols, negative_multipliers=neg[:5], bad=bad, structure=struct)

# per-certificate model specification (families as declared in the certificate's own metadata; names are the research flags)
XIV6 = dict(maximal=True, sub='core', nl='all', hull=True)
XIV7 = dict(maximal=True, budget=True, compact=True, point=True, pairb=True, hull=True)
SPECS = {'515': dict(D=5, cfg='515', sub='XI', cert='../lp/certs/coverage11_cert_515_root.json', dump='research_rownames_515.json'),
         '617': dict(D=5, cfg='617', sub='XI', xiii=True, cert='../lp/certs/coverage13_cert_pair_root_ipm.json', dump='research_rownames_617.json'),
         'D6_515': dict(D=6, cfg='515', cert='../lp/certs/cert14_D6_515_root.json', dump='research_rownames_D6_515.json', **XIV6),
         'D6_617': dict(D=6, cfg='617', cert='../lp/certs/cert14_D6_617_root.json', dump='research_rownames_D6_617.json', **XIV6),
         'D7_515': dict(D=7, cfg='515', sub='all', nl='all', tr='all', cert='../lp/certs/cert14_D7_515_root_strongroot.json', dump='research_rownames_D7_515.json', **XIV7),
         'D7_617_27': dict(D=7, cfg='617', prefix=(27,), sub='pair', nl='pair', star=True, rbudget=True, nonfano=True,
                           cert='../lp/certs/cert14_D7_617_27_d7high1_27.json', dump='research_rownames_D7_617_27.json', **XIV7),
         'D7_617_0': dict(D=7, cfg='617', prefix=(0,), sub='all', nl='all', tr='all', star=True, rbudget=True, nonfano=True, rich=True, mix=True, extra=True, rdeg=True,
                          cert='../lp/certs/cert14_D7_617_0_scaled0.json', dump='research_rownames_D7_617_0.json', **XIV7),
         'D7_617_dom': dict(D=7, cfg='617', prefix=(), domain_m1=1, sub='coordinate', star=True, rbudget=True, nonfano=True, rich=True, mix=True, psubset=True, rdeg=True,
                            cert='../lp/certs/cert14_D7_617_domain_m1ge1_subset_bounded.json', dump='research_rownames_D7_617_dom.json', **XIV7)}
def tkey(t): return jkey(list(t[:2]) + [list(t[2]), t[3], list(t[4]), [list(x) for x in t[5]] if t[5] else None])
def main():
    tag = sys.argv[1]; spec = SPECS[tag]; cert_file, dump, DD, cfg = spec['cert'], spec['dump'], spec['D'], spec['cfg']; t0 = time.time()
    log = open(f'out/root_{tag}.log', 'w')
    def say(*a):
        s = ' '.join(map(str, a)); print(s, flush=True); log.write(s + '\n'); log.flush()
    M = RootModel(spec); say('MODEL', tag, json.dumps(M.info))
    res = json.load(open('out/' + dump)); rr = res['rows']
    names = [r[0] for r in rr]; assert len(set(names)) == len(names), 'research row names not unique'
    mine, theirs = set(M.rows), set(names)
    say('ROWS mine', len(mine), 'research', len(theirs), 'only mine', len(mine - theirs), 'only research', len(theirs - mine), sorted(mine - theirs)[:3], sorted(theirs - mine)[:3])
    sr = {r[0]: (r[1], r[2]) for r in rr}
    srdiff = [(k, M.rows[k], sr[k]) for k in mine & theirs if M.rows[k] != sr[k]]
    say('SENSE/RHS differences', len(srdiff), srdiff[:5])
    say('row order identical to research:', names == M.order)
    say('COLUMNS mine', M.info['columns'], 'research', res['ncols'], res['col_prefix_counts'], 'types mine', M.info['types'], 'research', res['ntypes'])
    cert = json.load(open(cert_file)); mult = cert['multipliers']; assert all(type(v) == int for v in mult)
    assert len(mult) == len(names), (len(mult), len(names))
    y = {nm: v for nm, v in zip(names, mult) if v}
    unknown = [k for k in y if k not in M.rows]
    if unknown: say('FAILED: certificate uses rows absent from my model', unknown[:5]); return
    research_nnz = {r[0]: r[3] for r in rr}
    ok, det = M.verify(y, research_nnz)
    verdict = 'VERIFIED' if ok and not (mine ^ theirs) and not srdiff and det['structure']['rows_with_different_support'] == 0 else ('VERIFIED-but-model-differs' if ok else 'FAILED')
    say('RESULT', tag, verdict, json.dumps(det), 'research rhs', cert['detail']['rhs'], 'research min_column', cert['detail']['min_column'])
    # informational: compare my prices with the research tables (NOT used above)
    cmp = {}
    try:
        if DD == 5:
            tab = json.load(open('../big/lp_inputs/lift11_prices.json'))['prices']
            cmp['lift11'] = dict(types=len(M.priceXI), mismatches=sum(tab.get(tkey(M.T[i])) != p for i, p in M.priceXI.items()))
            if M.pair:
                tab = json.load(open('../lp/certs/lift13_pair_prices.json'))['prices']
                cmp['lift13'] = dict(types=len(M.priceXIII), mismatches=sum(tab.get(tkey(M.T[i])) != p for i, p in M.priceXIII.items()))
        else:
            def cmpdict(fn, pr):
                tab = json.load(open(fn))['prices']; mism = cnt = 0
                for i, p in pr.items():
                    th = tab.get(tkey(M.T[i]))
                    for z, v in p.items():
                        cnt += 1; mism += (th is None or z >= len(th) or th[z] != v)
                return dict(types=len(pr), prices=cnt, mismatches=mism)
            cmp['lift14'] = cmpdict(f'../big/lp_inputs/lift14_D{DD}_prices.json', M.priceSub)
            cmp['nonlinear14'] = cmpdict(f'../big/lp_inputs/nonlinear14_D{DD}_prices.json', M.priceNL)
            if M.F.get('tr'):
                cmp['triple14'] = cmpdict(f'../big/lp_inputs/triple14_D{DD}_prices.json', M.priceTr)
            if M.F.get('point'):
                cmp['point14'] = cmpdict(f'../big/lp_inputs/point14_D{DD}_prices.json', {i: dict(enumerate(p)) for i, p in M.pricePt.items()})
                cmp['pointtriple14'] = cmpdict(f'../big/lp_inputs/pointtriple14_D{DD}_prices.json', {i: dict(enumerate(p)) for i, p in M.pricePt3.items()})
    except Exception as ex: cmp['error'] = repr(ex)
    say('PRICE TABLE COMPARISON (informational)', json.dumps(cmp))
    out = {}
    if os.path.exists('out/root_result.json'): out = json.load(open('out/root_result.json'))
    out[tag] = dict(verdict=verdict, model=M.info, research_columns=res['ncols'], research_column_families=res['col_prefix_counts'], research_types=res['ntypes'], detail=det, research_rhs=cert['detail']['rhs'],
                    row_set_equal=not (mine ^ theirs), sense_rhs_diffs=len(srdiff), row_order_identical=(names == M.order), price_comparison=cmp, wall_s=round(time.time() - t0, 1))
    json.dump(out, open('out/root_result.json', 'w'), indent=1); say('wall', round(time.time() - t0, 1), 's')
if __name__ == '__main__': main()
