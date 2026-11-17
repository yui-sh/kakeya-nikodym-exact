"""AUDIT: exhaustive facts about direction sets in PG(2,5) used by the D=7 cuts (WRITEUP14):
  T = # collinear triples among j directions; g-sum = sum_L (C(|L∩S|,3) - C(|L∩S|,4)); rich secants (|L∩S|>=3);
  star bounds T <= C(a,3)+C(j-a,2)+C(j-a,3) / floor(C(j,2)/3) when the largest secant a <= 4.
  Fresh code: PG(2,5) built from scratch."""
import itertools, math, json
q = 5
pts = sorted({tuple(v) for v in itertools.product(range(q), repeat=3) if any(v) and v[next(i for i in range(3) if v[i])] == 1})
assert len(pts) == 31
def dot(a, b): return sum(x * y for x, y in zip(a, b)) % q
lines = [frozenset(i for i, p in enumerate(pts) if dot(a, p) == 0) for a in pts]
assert all(len(l) == 6 for l in lines) and all(len(l1 & l2) == 1 for l1, l2 in itertools.combinations(lines, 2))
C = math.comb
res = {}
for j in range(3, 8):
    Tmin = Tmax = gmin = gmax = richmin_cap = richmax_cap = Tcapmax = None
    star = {}
    cnt = 0
    for S in itertools.combinations(range(31), j):
        Sset = set(S); secs = [len(l & Sset) for l in lines]
        T = sum(C(s, 3) for s in secs); g = sum(C(s, 3) - C(s, 4) for s in secs); rich = sum(s >= 3 for s in secs); a = max(secs)
        cnt += 1
        Tmin = T if Tmin is None else min(Tmin, T); Tmax = T if Tmax is None else max(Tmax, T)
        gmin = g if gmin is None else min(gmin, g); gmax = g if gmax is None else max(gmax, g)
        if a <= 4:
            richmin_cap = rich if richmin_cap is None else min(richmin_cap, rich); richmax_cap = rich if richmax_cap is None else max(richmax_cap, rich)
            Tcapmax = T if Tcapmax is None else max(Tcapmax, T)
    # star bound formula
    bounds = [0]
    if j >= 3: bounds.append(C(j, 2) // 3)
    for a in range(4, min(j, 4) + 1): bounds.append(C(a, 3) + C(j - a, 2) + C(j - a, 3))
    star_bound = max(bounds)
    res[j] = dict(subsets=cnt, T=(Tmin, Tmax), g=(gmin, gmax), rich_cap4=(richmin_cap, richmax_cap), Tmax_cap4=Tcapmax, star_bound_cap4=star_bound, star_ok=(Tcapmax <= star_bound))
    print(j, res[j])
assert res[7]['T'] == (2, 20) and res[7]['g'] == (2, 6) and res[7]['rich_cap4'][0] == 2
assert all(res[j]['star_ok'] for j in (5, 6, 7)) and [res[j]['star_bound_cap4'] for j in (5, 6, 7)] == [4, 5, 8]
print('PG(2,5) claims for degree-7 points: T in [2,20], g in [2,6], >=2 rich secants under cap 4, star bounds 4,5,8 (j=5,6,7): ALL CONFIRMED')
json.dump({str(k): v for k, v in res.items()}, open('out/projective_result.json', 'w'))
