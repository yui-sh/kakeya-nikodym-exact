"""AUDIT: compare the fresh 2^25 enumeration (out/planar_raw.txt from check_planar.c) with the tables the LP chain trusts:
   lp/certs/{inside_profiles_q5,plane_f,sections3_q5_r5,union_tables_q5}.json, and check Lemma A, the planar table (P),
   the plane-size filter caps (WRITEUP14) and the b<=4 claim."""
import json, itertools, sys, math
C = '../lp/certs/'
raw = open('out/planar_raw.txt').read().splitlines()
sec = set(); prof = {}
for ln in raw:
    f = ln.split()
    if f[0] == 'SEC':
        parts = ln[4:].split('|'); a = list(map(int, parts[0].split())); t = tuple(map(int, parts[1].split()))
        k, r, u = a[:3]; c = tuple(a[3:3 + r]); e = tuple(map(int, parts[2].split())) if len(parts) > 2 else None
        sec.add((k, r, u, c, t, e))
    elif f[0] == 'PROF': prof[tuple(map(int, f[1:7]))] = int(f[7])
print('fresh SEC tuples:', len(sec))
# --- sections table ---
theirs = set()
for k, r, u, c, t, e in json.load(open(C + 'sections3_q5_r5.json')):
    ee = None if e is None else tuple(e[0])   # r==5: one unmarked direction -> tuple of one sorted 5-tuple
    theirs.add((k, r, u, tuple(c), tuple(t), ee))
print('their sections3_q5_r5 tuples:', len(theirs), ' equal:', theirs == sec, ' only-theirs:', len(theirs - sec), ' only-mine:', len(sec - theirs))
# --- inside profiles and planar table f ---
ip = {int(r): {tuple(p) for p in v} for r, v in json.load(open(C + 'inside_profiles_q5.json')).items()}
mine_ip = {}
for k, r, u, c, t, e in sec: mine_ip.setdefault(r, set()).add((u,) + c)
print('inside_profiles equal:', all(mine_ip[r] == ip[r] for r in range(7)), {r: sorted(mine_ip[r]) for r in range(7)})
f = {int(r): v for r, v in json.load(open(C + 'plane_f.json'))['5'].items()}
mine_f = {r: min(u for (u, *c) in mine_ip[r]) for r in range(7)}
print('planar table f(5,r):', mine_f, ' theirs:', f, ' equal:', mine_f == f)
# --- Lemma A: every (k=15,r=5,u=15) tuple has t=(0,10,15,0,5) and e=(3,3,3,3,3); and every r=5,u=15 plane (any k) has
#     e-sections >= 3 each and the parallel sections of a marked direction are >= 2,2,3,3 (checked via t only loosely) ---
la = [x for x in sec if x[1] == 5 and x[2] == 15 and x[0] == 15]
print('Lemma A tuples (k=15,r=5,u=15):', la)
print('Lemma A(i) for all r=5,u=15 (any k): min e-section >= 3:', all(min(x[5]) >= 3 for x in sec if x[1] == 5 and x[2] == 15))
# --- Ft table (min |S| over S with sorted full-line profile dominating a) and F(a) ---
Ft_theirs = {tuple(map(int, k.split(','))): v for k, v in json.load(open(C + 'union_tables_q5.json'))['Ft'].items()}
F_theirs = {int(k): v for k, v in json.load(open(C + 'union_tables_q5.json'))['F'].items()}
def dominates(p, a): return all(pi >= ai for pi, ai in zip(p, a))   # both sorted ascending: sufficient and necessary for a matching
Ft_mine = {}
for a in itertools.combinations_with_replacement(range(1, 6), 6):
    Ft_mine[a] = min(v for p, v in prof.items() if dominates(p, a))
bad = {a: (Ft_mine[a], Ft_theirs[a]) for a in Ft_mine if Ft_mine[a] != Ft_theirs.get(a)}
print('Ft table: entries', len(Ft_mine), ' mismatches vs corrected table:', bad)
print('Ft(1,1,2,2,3,3) mine =', Ft_mine[(1, 1, 2, 2, 3, 3)])
orig = {tuple(map(int, k.split(','))): v for k, v in json.load(open(C + 'union14_original_table.json'))['Ft'].items()}
print('mismatches vs ORIGINAL (pre-correction) table:', {a: (Ft_mine[a], orig[a]) for a in Ft_mine if Ft_mine[a] != orig.get(a)})
unsafe = {a: (Ft_mine[a], Ft_theirs[a]) for a in Ft_mine if Ft_theirs[a] > Ft_mine[a]}
print('UNSAFE entries (their bound > true min):', unsafe)
F_mine = {a: min(v for t, v in Ft_mine.items() if sum(t) == a) for a in range(6, 31)}
print('F(a) equal:', F_mine == F_theirs, F_mine)
# --- plane-size filter from the budget: c(r,k) <= 541 ---
def c(r, k): return 10 * k * k - (150 + 17 * r) * k + 560 + 143 * r
caps = {r: max(k for k in range(mine_f[r], 26) if c(r, k) <= 541) for r in range(7)}
print('budget caps (max k with c<=541):', caps, ' claimed 14,15,16,17,18,19,21; r=0 k=0 excluded:', c(0, 0) > 541)
print('c(r,k) >= 0 on all attainable (r,k):', all(c(x[1], x[0]) >= 0 for x in sec), ' min over sec:', min(c(x[1], x[0]) for x in sec))
# --- b<=4 claim: among retained section profiles (r<=rmax, u==15/17 at r==rmax, c<=541) the max inside multiplicity is 4 ---
for rmax, u0 in ((5, 15), (6, 17)):
    adm = [x for x in sec if x[1] <= rmax and (x[1] < rmax or x[2] == u0) and c(x[1], x[0]) <= 541]
    B = max(max((b for b, n in enumerate(x[3], 1) if n), default=0) for x in adm)
    print(f'config rmax={rmax}: retained profiles {len(adm)}, max inside multiplicity {B}')
# --- hand check of c(5,21)=750 and f(6)=17 (Blokhuis-Mazzocca) ---
print('c(5,21) =', c(5, 21), '; f(5,6) =', mine_f[6])
json.dump({'sections_equal': theirs == sec, 'Ft_mismatch': {','.join(map(str, a)): v for a, v in bad.items()}, 'unsafe': len(unsafe), 'f': mine_f, 'caps': caps}, open('out/planar_tables_result.json', 'w'))
