"""AUDIT: every identity / lemma / cut the LP chain trusts, checked on many random, greedy and structured witness systems
in AG(3,q), q = 3,4,5 (they must hold for EVERY witness system, not only Kakeya-minimal ones). Fresh code (ag3.py)."""
import sys, json, math, random, itertools
from collections import Counter, defaultdict
from ag3 import AG3, random_witness, stats
C2 = math.comb
qs = [int(a) for a in sys.argv[1:]] or [5, 4, 3]
NSYS = int(__import__('os').environ.get('NSYS', '60'))
# planar tables from the fresh enumeration (q=5 only)
sec5 = set(); prof5 = {}
for ln in open('out/planar_raw.txt'):
    f = ln.split()
    if f[0] == 'SEC':
        parts = ln[4:].split('|'); a = list(map(int, parts[0].split())); t = tuple(map(int, parts[1].split()))
        k, r, u = a[:3]; c = tuple(a[3:3 + r]); e = tuple(map(int, parts[2].split())) if len(parts) > 2 else None
        sec5.add((k, r, u, c, t, e))
    elif f[0] == 'PROF': prof5[tuple(map(int, f[1:7]))] = int(f[7])
def Ft5(a):
    a = tuple(sorted(a)); return min(v for p, v in prof5.items() if all(x >= y for x, y in zip(p, a)))
fq = {3: [0, 3, 5, 6, 7], 4: [0, 4, 7, 9, 10, 10], 5: [0, 5, 9, 12, 14, 15, 17]}   # planar table f(q,r) (q=5 verified fresh; 3,4 from plane_f, only used for Lemma B(c) check)
def cap3(r, k):
    if r == 5 and k == 15: return 2
    if r == 1: return min(4, (k - 5) // 3)
    if r == 2: return min(4, (k - 9) // 2)
    if r == 3: return min(4, k - 12)
    return 4
fails = Counter(); checks = Counter()
def chk(name, cond, info=''):
    checks[name] += 1
    if not cond:
        fails[name] += 1
        if fails[name] <= 3: print('FAIL', name, info)
for q in qs:
    ag = AG3(q); nd = q * q + q + 1; rng = random.Random(1000 + q)
    f = fq[q]
    modes = [('random', 0)] * (NSYS // 3) + [('greedy', 0)] * (NSYS // 3) + [('forced', 3)] * (NSYS - 2 * (NSYS // 3))
    for mode, fp in modes:
        W, K = random_witness(ag, rng, mode, fp); n = len(K)
        deg, D, P = stats(ag, W, K); m = Counter(deg.values())
        Wl = {li: d for d, li in W.items()}
        # ---- I1..I6 ----
        for pl in P:
            chk('I1', sum(deg[p] for p in pl['Kp']) == q * q + q * pl['r'])
            chk('P: k >= f(r)', pl['k'] >= f[pl['r']] and pl['u'] >= f[pl['r']])
        chk('I2', len(P) == q * nd and sum(pl['r'] for pl in P) == (q + 1) * nd and sum(pl['k'] for pl in P) == nd * n)
        chk('I3', sum(C2(pl['r'], 2) for pl in P) == sum(C2(d, 2) for d in deg.values()))
        chk('I4', sum(pl['k'] ** 2 for pl in P) == (q + 1) * n * n + q * q * n)
        chk('I5', sum(pl['k'] * pl['r'] for pl in P) == nd * n + q * q * nd)
        L = Counter(len(pts & K) for di, pts in ag.lines)
        chk('I6a', sum(j * L[j] for j in L) == nd * n); chk('I6b', sum(C2(j, 2) * L[j] for j in L) == C2(n, 2))
        chk('I6c: L_q >= nd', L[q] >= nd)
        chk('I6e', sum(C2(pl['k'], 3) for pl in P) == C2(n, 3) + q * sum(C2(j, 3) * L[j] for j in L))
        if q == 5:
            chk('quadratic identity', sum(10 * pl['k'] ** 2 - (150 + 17 * pl['r']) * pl['k'] + 560 + 143 * pl['r'] for pl in P) == 60 * n * n - 4927 * n + 100223)
        # F5: sum_pi t_j = (q+1) L_j
        for j in range(1, q + 1): chk('F5', sum(pl['t'][j - 1] for pl in P) == (q + 1) * L[j])
        # ---- F1 histogram identities ----
        M = Counter(tuple(sorted((deg[a], deg[b]))) for a, b in itertools.combinations(sorted(K), 2))
        for j in range(1, D + 1):
            chk('H1', sum(pl['h'][j - 1] for pl in P) == nd * m[j])
            chk('H2', sum(pl['h'][j - 1] * pl['k'] for pl in P) == ((q + 1) * n + q * q) * m[j])
            chk('H3', sum(pl['h'][j - 1] * pl['r'] for pl in P) == (nd + q * j) * m[j])
            chk('H4', sum(C2(pl['h'][j - 1], 2) for pl in P) == (q + 1) * M[(j, j)])
            for i in range(1, j): chk('H5', sum(pl['h'][i - 1] * pl['h'][j - 1] for pl in P) == (q + 1) * M[(i, j)])
            chk('H6', sum(M[(min(i, j), max(i, j))] for i in range(1, D + 1) if i != j) + 2 * M[(j, j)] == (n - 1) * m[j])
            chk('M_jj = C(m_j,2)', M[(j, j)] == C2(m[j], 2))
            for i in range(1, j): chk('M_ij = m_i m_j', M[(i, j)] == m[i] * m[j])
        # crossing filter: #points on >= i inside lines <= #points of degree >= i
        for pl in P:
            for i in range(1, pl['r'] + 1): chk('crossing Hall filter', sum(pl['c'][i - 1:]) <= sum(pl['h'][i - 1:]))
        # ---- sections table membership (q=5): every real plane's (k,r,u,c,t,e) is in the fresh 2^25 table ----
        if q == 5:
            for pl in P:
                e = tuple(next(iter(pl['e'].values()))) if pl['r'] == 5 else None
                chk('plane type in sections table', (pl['k'], pl['r'], pl['u'], pl['c'], pl['t'], e) in sec5, (pl['k'], pl['r'], pl['u'], pl['c'], pl['t'], e))
        # ---- F3 pencils, Lemma C, Lemma D, wline degree pairs ----
        for li, d in Wl.items():
            pen = ag.planes_of_line[li]; s = sum(deg[p] for p in ag.lines[li][1]); dm = Counter(deg[p] for p in ag.lines[li][1])
            chk('pencil sum k', sum(P[pi]['k'] for pi in pen) == n + q * q)
            chk('pencil sum r = 1+s', sum(P[pi]['r'] for pi in pen) == 1 + s)
            chk('pencil r>=1', all(P[pi]['r'] >= 1 for pi in pen))
            for pi in pen:
                chk('Lemma C', s >= q * q + q * P[pi]['r'] - D * (P[pi]['k'] - q), (s, P[pi]['r'], P[pi]['k'], D))
                chk('Lemma D(a)', all(P[pi]['h'][j - 1] >= dm[j] for j in dm))
            for j in range(1, D + 1): chk('Lemma D(b)', sum(P[pi]['h'][j - 1] for pi in pen) == m[j] + q * dm[j])
        chk('sum s_l = sum d^2', sum(sum(deg[p] for p in ag.lines[li][1]) for li in Wl) == sum(d * d for d in deg.values()))
        pairs_on_wl = Counter()
        for li in Wl:
            for a, b in itertools.combinations(sorted(ag.lines[li][1]), 2): pairs_on_wl[tuple(sorted((deg[a], deg[b])))] += 1
        chk('witness pairs <= M', all(pairs_on_wl[k] <= M[k] for k in pairs_on_wl))
        # ---- F2 classes and F4 directions/projections ----
        classes = defaultdict(list)
        for pi, pl in enumerate(P): classes[ag.planes[pi][0]].append(pi)
        for ci, pis in classes.items():
            chk('class sum r = q+1', sum(P[pi]['r'] for pi in pis) == q + 1); chk('class sum k = n', sum(P[pi]['k'] for pi in pis) == n)
            rho = sum(P[pi]['r'] >= 1 for pi in pis)
            for d in ag.planes[pis[0]][3]:
                # a_{d,L} = # planes of the class containing a witness line other than l_d
                a = sum(any(Wl[x] != d for x in P[pi]['ins']) for pi in pis)
                alone = any(P[pi]['ins'] == [W[d]] for pi in pis)
                chk('a_{d,L} class formula', a == (rho - 1 if alone else rho))
        for d in range(nd):
            # projection along d: lines parallel to d meeting K; full lines in quotient direction L correspond to planes of class L containing a witness line != l_d
            proj = [li for li in ag.lines_by_dir[d] if ag.lines[li][1] & K]
            chk('|P_d| = sum L_j^{(d)}', len(proj) == sum(1 for li in ag.lines_by_dir[d] if ag.lines[li][1] & K))
            avec = []
            for ci, pis in classes.items():
                if d in ag.planes[pis[0]][3]: avec.append(sum(any(Wl[x] != d for x in P[pi]['ins']) for pi in pis))
            chk('a in [1,q]', len(avec) == q + 1 and all(1 <= a <= q for a in avec))
            if q == 5: chk('|P_d| >= Ft(a)', len(proj) >= Ft5(avec), (len(proj), sorted(avec)))
            # projection is Kakeya in the quotient plane: for every class L through d, some plane of L has all q lines parallel to d meeting K
            for ci, pis in classes.items():
                if d in ag.planes[pis[0]][3]:
                    chk('P_d Kakeya', any(all(ag.lines[li][1] & K for li in ag.planes[pi][4] if ag.lines[li][0] == d) for pi in pis))
        # ---- Lemma A / B / F6 / R3 / R4 / F9 ----
        for pi, pl in enumerate(P):
            if pl['r'] == 5 and pl['u'] == 5 * q - 10 and q == 5:
                chk('Lemma A(i) e-sections', all(v == 3 for v in next(iter(pl['e'].values()))) if pl['k'] == 15 else all(v >= 3 for v in next(iter(pl['e'].values()))))
                if pl['k'] == 15:
                    for li in pl['ins']:
                        par = sorted(len(ag.lines[x][1] & pl['Kp']) for x in ag.planes[pi][4] if ag.lines[x][0] == ag.lines[li][0] and x != li)
                        chk('Lemma A(ii) parallel sections 2,2,3,3', par == [2, 2, 3, 3], par)
            if pl['r'] >= 1:
                for e, sprof in pl['e'].items():
                    le = W[e]; mates = [pj for pj in classes[ag.planes[pi][0]] if le in ag.planes[pj][4]]
                    chk('Lemma B(a) unique mate', len(mates) == 1 and mates[0] != pi)
                    mate = mates[0]; chk('mate r>=1', P[mate]['r'] >= 1)
                    if pl['r'] == 5: chk('Lemma B(a) r=5: mate r=1, others 0', P[mate]['r'] == 1 and all(P[pj]['r'] == 0 for pj in classes[ag.planes[pi][0]] if pj not in (pi, mate)))
                    pen = [pj for pj in ag.planes_of_line[le] if pj != mate]; chk('pencil minus mate = q planes', len(pen) == q)
                    elines = [x for x in ag.planes[pi][4] if ag.lines[x][0] == e]
                    inter = {pj: [x for x in elines if x in ag.planes[pj][4]] for pj in pen}
                    chk('Lemma B(b) bijection', all(len(v) == 1 for v in inter.values()) and len({v[0] for v in inter.values()}) == q)
                    for pj, (x,) in inter.items():
                        sj = len(ag.lines[x][1] & K); rs, ks = P[pj]['r'], P[pj]['k']
                        chk('Lemma B(c) exact', ks >= len(P[pj]['U']) + max(0, sj - (rs - 1)))
                        chk('Lemma B(c) with f', ks >= f[rs] + max(0, sj - rs + 1))
                    # F9: unmarked Lemma C on the mate
                    s_e = sum(deg[p] for p in ag.lines[le][1])
                    chk('F9 unmarked Lemma C', s_e >= q * q + q * P[mate]['r'] - D * (P[mate]['k'] - q))
        if q == 5:
            for e in range(nd):
                le = W[e]; pen = ag.planes_of_line[le]
                X = [pi for pi, pl in enumerate(P) if pl['r'] == 5 and pl['k'] == 15 and pl['u'] == 15 and e in pl['e']]
                x = len(X)
                hosted = set()
                for pi in X:
                    for li in ag.planes[pi][4]:
                        if ag.lines[li][0] == e and len(ag.lines[li][1] & K) == 3: hosted.add(li)
                chk('R4 lower: distinct 3-point e-lines >= 5x-C(x,2)', len(hosted) >= 5 * x - C2(x, 2), (x, len(hosted)))
                tot = sum(cap3(P[pj]['r'], P[pj]['k']) for pj in pen)
                chk('R4 capacity: hosted <= sum cap3', len(hosted) <= tot, (x, len(hosted), tot))
                for pj in pen:   # per-plane capacity
                    hp = [li for li in ag.planes[pj][4] if ag.lines[li][0] == e and li != le and len(ag.lines[li][1] & K) == 3]
                    chk('cap3 per plane', len(hp) <= cap3(P[pj]['r'], P[pj]['k']), (P[pj]['r'], P[pj]['k'], len(hp)))
                # R3 with rmin = 5: flags (pi,e) with r_pi >= 5 on the pencil of l_e <= sum_{pj in pencil} floor((6-r_pj)/5)
                flags = sum(1 for pi, pl in enumerate(P) if pl['r'] >= 5 and e in pl['e'])
                chk('R3 (rmin=5)', flags <= sum((q + 1 - P[pj]['r']) // 5 for pj in pen if P[pj]['r'] >= 1))
                chk('R3 (rmin=2)', sum(1 for pi, pl in enumerate(P) if pl['r'] >= 2 and e in pl['e']) <= sum((q + 1 - P[pj]['r']) // 2 for pj in pen if P[pj]['r'] >= 1))
        # ---- point-level identities (WRITEUP14 point families) ----
        for p in K:
            j = deg[p]; ds = [Wl[li] for li in ag.lines_through[p] if li in Wl]
            chk('sum_pi C(b,2) = C(j,2)', sum(C2(P[pi]['b'][p], 2) for pi in ag.planes_through[p]) == C2(j, 2))
            T = sum(1 for a, b, c in itertools.combinations(ds, 3) if ag.collinear(a, b, c))
            chk('sum_pi C(b,3) = collinear triples', sum(C2(P[pi]['b'][p], 3) for pi in ag.planes_through[p]) == T)
            chk('b <= j', all(P[pi]['b'][p] <= j for pi in ag.planes_through[p]))
            if q == 5 and j == 7: chk('T in [2,20] at degree 7', 2 <= T <= 20)
            # parallel-class partition: sum over the q planes of a class of h = m (Lemma D-like) -- trivially true; skip
        # ---- prefix caps: h_j <= m_j for all planes ----
        for pl in P: chk('h_j <= m_j', all(pl['h'][j - 1] <= m[j] for j in range(1, D + 1)))
        # ---- degree range ----
        chk('D >= ceil(nd q / n) and 1+(q-1)D <= n', D >= math.ceil(nd * q / n) and 1 + (q - 1) * D <= n)
    print(f'q={q}: {len(modes)} witness systems checked')
print('CHECKS:', dict(checks))
print('FAILS:', dict(fails) if fails else 'none')
json.dump({'checks': dict(checks), 'fails': dict(fails)}, open('out/identities_result.json', 'w'))
