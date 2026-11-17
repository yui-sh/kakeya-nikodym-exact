"""AUDIT: independent re-pricing of the m-free cut families (WRITEUP11/13/14) from a FRESH geometry catalogue
(out/geometry_raw.txt from check_geometry.c; out/planar_raw.txt from check_planar.c), compared with the price tables the
certificates were built from (kakeya-fq/out/*_prices.json, read-only). Also compares the catalogue itself with
lp/certs/geometry14_D6.json. Exact minimisers: group DP for additive (b,d)-costs, exhaustive colourings for within-line costs."""
import json, math, random, sys, itertools, os, time
from functools import lru_cache
from collections import defaultdict
C = math.comb; O = '../../kakeya-fq/out/'; rng = random.Random(7)
NS = int(os.environ.get('NSAMPLE', '250'))
def freeze(x): return tuple(map(freeze, x)) if isinstance(x, list) else x
# ---------- fresh catalogues ----------
geo = defaultdict(set)   # key (r,k,u,t,e) -> set of sorted mask multisets (masks over the r inside lines)
for ln in open('out/geometry_raw.txt'):
    a, t, e, mm = ln.split('|'); k, r, u = map(int, a.split()); t = tuple(map(int, t.split())); e = tuple(map(int, e.split())) or None
    e = (e,) if e else None
    geo[(r, k, u, t, e)].add(tuple(map(int, mm.split())))
prof = defaultdict(set)  # key (r,k,u,t,e) -> set of multiplicity profiles (count of points with b = 0..r), ALL r
for ln in open('out/planar_raw.txt'):
    f = ln.split()
    if f[0] != 'SEC': continue
    parts = ln[4:].split('|'); a = list(map(int, parts[0].split())); t = tuple(map(int, parts[1].split()))
    k, r, u = a[:3]; c = tuple(a[3:3 + r]); e = ((tuple(map(int, parts[2].split())),) if len(parts) > 2 else None)
    prof[(r, k, u, t, e)].add((k - u,) + c)
print('fresh catalogue keys (r>=5):', len(geo), ' profile keys (all r):', len(prof))
# compare with their geometry14 catalogue
theirs = defaultdict(set)
for g in json.load(open('../lp/certs/geometry14_D6.json')):
    s = freeze(g['signature']); theirs[s[:2] + s[3:]].add(tuple(sorted(g['masks'])))
# masks in their catalogue are bit positions in THEIR line order; compare multiplicity profiles (popcounts) per key and the
# multiset of line-incidence structures up to relabelling via a canonical form
def canon(mm, r):
    best = None
    for perm in itertools.permutations(range(r)):
        cm = tuple(sorted(sum(((m >> i) & 1) << perm[i] for i in range(r)) for m in mm))
        if best is None or cm < best: best = cm
    return best
same_keys = set(geo) == set(theirs); print('same key set as geometry14_D6:', same_keys, len(geo), len(theirs))
mism = 0
for key in geo:
    A = {canon(mm, key[0]) for mm in geo[key]}; B = {canon(mm, key[0]) for mm in theirs[key]}
    if A != B: mism += 1; print('  catalogue mismatch at', key, len(A), len(B)) if mism < 5 else None
print('keys whose geometry sets differ (up to line relabelling):', mism)
# ---------- exact pricers ----------
@lru_cache(None)
def additive_min(profile, h, costf_id, costrows):
    """min over assignments of degrees (histogram h, degrees 1..D) to points (profile[b] points of multiplicity b,
    constraint d >= max(1,b)) of sum cost(b,d); costrows[b][d-1] = cost. Group DP over b descending."""
    D = len(h); groups = [(b, n) for b, n in enumerate(profile) if n]
    @lru_cache(None)
    def rec(gi, rem):
        if gi == len(groups): return 0 if not any(rem) else None
        b, n = groups[gi]; best = None
        ds = [d for d in range(max(1, b), D + 1) if rem[d - 1]]
        def comp(i, left, rem2, acc):
            nonlocal best
            if i == len(ds):
                if left == 0:
                    v = rec(gi + 1, rem2)
                    if v is not None and (best is None or acc + v < best): best = acc + v
                return
            d = ds[i]
            for x in range(min(left, rem2[d - 1]) + 1):
                r2 = list(rem2); r2[d - 1] -= x; comp(i + 1, left - x, tuple(r2), acc + x * costrows[b][d - 1])
        comp(0, n, rem, 0); return best
    return rec(0, tuple(h))
def as_profile(x, r):
    if isinstance(x, tuple) and x and max(x) > r: pass
    return x
geo_prof = {key: {tuple(sum(bin(m).count('1') == b for m in mm) for b in range(key[0] + 1)) for mm in mms} for key, mms in geo.items()}
def price_additive(key, h, costrows, catalogue):
    cat = geo_prof if catalogue is geo else catalogue
    vals = [v for p in cat[key] if (v := additive_min(p, tuple(h), id(costrows), costrows)) is not None]
    return min(vals) if vals else None
def hall_ok(sel_counts, degs):   # sel_counts[b] points needing degree >= max(1,b); degs: multiset tuple of degrees available
    if sum(sel_counts) != len(degs): return False
    return all(sum(sel_counts[b:]) <= sum(1 for d in degs if d >= b) for b in range(1, len(sel_counts)))
@lru_cache(None)
def colour_extrema(mm, r, power):
    """for a geometry (sorted masks incl. zeros = extra points): colour inside points; return dict profile -> (min,max) of
    sum_lines C(#selected on line, power); profile = counts of selected points per b >= 1 (b=0 handled separately)."""
    inside = [m for m in mm if m]; out = {}
    n = len(inside); bc = [bin(m).count('1') for m in inside]
    line_cnt = [0] * r; cnt = [0] * (r + 1); cost = 0; prev = 0
    for num in range(1 << n):
        g = num ^ (num >> 1)
        if num:
            ch = g ^ prev; p = ch.bit_length() - 1; add = bool(g & ch); m = inside[p]
            cnt[bc[p]] += 1 if add else -1
            for i in range(r):
                if m >> i & 1:
                    if add: cost += C(line_cnt[i] + 1, power) - C(line_cnt[i], power); line_cnt[i] += 1
                    else: cost -= C(line_cnt[i], power) - C(line_cnt[i] - 1, power); line_cnt[i] -= 1
        prev = g
        pr = tuple(cnt[1:])
        if pr in out:
            lo, hi = out[pr]; out[pr] = (min(lo, cost), max(hi, cost))
        else: out[pr] = (cost, cost)
    return out
def price_withinline(key, h, S, sign, power):
    r = key[0]; extra = key[1] - key[2]; best = None
    sel = tuple(x for d, x in enumerate(h, 1) if d in S); seldegs = tuple(d for d, x in enumerate(h, 1) if d in S for _ in range(x))
    othdegs = tuple(d for d, x in enumerate(h, 1) if d not in S for _ in range(x)); s = len(seldegs)
    for mm in geo[key]:
        tot = [0] * (r + 1)
        for m in mm: tot[bin(m).count('1')] += 1
        for pr, (lo, hi) in colour_extrema(mm, r, power).items():
            v0 = s - sum(pr)
            if not 0 <= v0 <= extra: continue
            selc = (v0,) + pr; othc = tuple(t - c for t, c in zip(tot, selc))
            if hall_ok(selc, seldegs) and hall_ok(othc, othdegs):
                v = lo if sign == 1 else -hi
                best = v if best is None else min(best, v)
    return best
report = {}
def sample_types(prices, n): return rng.sample(list(prices), min(n, len(prices)))
def run_additive(name, fn, catalogue, costfun, objs_key='objectives', keyfilter=None, keytuple=lambda t: (t[0], t[1], t[3], t[4], t[5])):
    t0 = time.time(); tab = json.load(open(O + fn)); objs = tab[objs_key]; prices = tab['prices']
    types = sample_types(prices, NS); mism = 0; checked = 0; missing = 0
    for tk in types:
        t = freeze(json.loads(tk)); key = keytuple(t); h = t[2]
        if keyfilter and not keyfilter(t): continue
        if key not in catalogue: missing += 1; print('  MISSING KEY', name, key); continue
        theirs = prices[tk]
        objlist = objs if len(objs) == len(theirs) else objs[:len(theirs)]
        for z, ob in enumerate(objlist):
            costrows = costfun(ob, len(h))
            mine = price_additive(key, h, costrows, catalogue); checked += 1
            if mine != theirs[z]:
                mism += 1
                if mism <= 5: print('  MISMATCH', name, t, ob, 'mine', mine, 'theirs', theirs[z])
    report[name] = dict(types=len(types), prices_checked=checked, mismatches=mism, missing_keys=missing, seconds=round(time.time() - t0, 1))
    print(name, report[name], flush=True)
def rows(f, D, r=7): return tuple(tuple(f(b, d) for d in range(1, D + 1)) for b in range(r + 1))
# lift14 D6/D7: additive subset costs sign*b*[d in S]; objectives (subset, sign); non-core types have only the first 2D objectives
for D in (6, 7):
    run_additive(f'lift14_D{D} (additive subsets)', f'lift14_D{D}_prices.json', geo, lambda ob, DD: rows(lambda b, d: ob[1] * b * (d in ob[0]), DD))
# lift11 (D5 coordinate costs), lift13 (pairs) : objectives implicit
tab = json.load(open(O + 'lift11_prices.json')); tab['objectives'] = [((j,), s) for j in range(1, 6) for s in (1, -1)]; json.dump(tab, open('out/_lift11_tmp.json', 'w'))
run_additive('lift11 D5 (coordinate costs)', '../../kakeya-public/audit/out/_lift11_tmp.json', geo, lambda ob, DD: rows(lambda b, d: ob[1] * b * (d in ob[0]), DD))
tab = json.load(open(O + 'lift13_pair_prices.json')); tab['objectives'] = [((i, j), s) for i, j in itertools.combinations(range(1, 6), 2) for s in (1, -1)]; json.dump(tab, open('out/_lift13_tmp.json', 'w'))
run_additive('lift13 D5 (pair costs)', '../../kakeya-public/audit/out/_lift13_tmp.json', geo, lambda ob, DD: rows(lambda b, d: ob[1] * b * (d in ob[0]), DD))
# point families (all r, catalogue = profiles): objectives (j, sign) -> cost sign*[d==j]*w(b)
run_additive('point14 D7 (C(b,2) by degree)', 'point14_D7_prices.json', prof, lambda ob, DD: rows(lambda b, d: ob[1] * (d == ob[0]) * C(b, 2), DD))
run_additive('pointtriple14 D7 (C(b,3) by degree)', 'pointtriple14_D7_prices.json', prof, lambda ob, DD: rows(lambda b, d: ob[1] * (d == ob[0]) * C(b, 3), DD))
run_additive('pointnonfano14 D7 (C(b,3)-C(b,4))', 'pointnonfano14_D7_prices.json', prof, lambda ob, DD: rows(lambda b, d: ob[1] * (d == ob[0]) * (C(b, 3) - C(b, 4)), DD))
run_additive('pointrich14 D7 ([b>=3])', 'pointrich14_D7_prices.json', prof, lambda ob, DD: rows(lambda b, d: ob[1] * (d == ob[0]) * (b >= 3), DD))
run_additive('pointmix14 D7 (lam C(b,2)+[b>=3])', 'pointmix14_D7_prices.json', prof, lambda ob, DD: rows(lambda b, d: ob[2] * (d == ob[0]) * (ob[1] * C(b, 2) + (b >= 3)), DD))
run_additive('pointsubset14 D7 (C(b,2), degree subsets)', 'pointsubset14_D7_prices.json', prof, lambda ob, DD: rows(lambda b, d: ob[1] * (d in ob[0]) * C(b, 2), DD))
# within-line families: nonlinear14 (pairs), triple14 (triples), extra14 (r=6,u=17 incl. extra points; objectives (S,sign,power))
def run_withinline(name, fn, nsample, objmap):
    t0 = time.time(); tab = json.load(open(O + fn)); objs = tab['objectives']; prices = tab['prices']
    types = sample_types(prices, nsample); mism = checked = 0
    for tk in types:
        t = freeze(json.loads(tk)); key = (t[0], t[1], t[3], t[4], t[5]); h = t[2]; theirs = prices[tk]
        for z, ob in enumerate(objs):
            S, sign, power = objmap(ob)
            mine = price_withinline(key, h, set(S), sign, power); checked += 1
            if mine != theirs[z]:
                mism += 1
                if mism <= 5: print('  MISMATCH', name, t, ob, 'mine', mine, 'theirs', theirs[z])
    report[name] = dict(types=len(types), prices_checked=checked, mismatches=mism, seconds=round(time.time() - t0, 1)); print(name, report[name], flush=True)
run_withinline('nonlinear14 D6 (within-line pairs)', 'nonlinear14_D6_prices.json', 12, lambda ob: (ob[0], ob[1], 2))
run_withinline('triple14 D6 (within-line triples)', 'triple14_D6_prices.json', 8, lambda ob: (ob[0], ob[1], 3))
run_withinline('extra14 D7 (r=6,u=17 with extra points)', 'extra14_D7_prices.json', 12, lambda ob: (ob[0], ob[1], ob[2]))
json.dump(report, open('out/prices_result.json', 'w'), indent=1)
print('TOTAL mismatches:', sum(v['mismatches'] for v in report.values()))
