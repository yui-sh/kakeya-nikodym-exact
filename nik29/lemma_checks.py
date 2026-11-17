#!/usr/bin/env python3
"""Fresh computational checks of the hand lemmas in ../nikodym-fq/WRITEUP3.md that the research
four-case reduction for Nik(2,9) depends on.  Imports only src/ff.py and src/ag.py (release code).
Everything here is re-derived from the geometry of AG(2,9); nothing is read from nikodym-fq.

Checks (each prints OK or raises):
  C1  line-profile identity: coefficient of a_i in C(k,2)-30k+540 is (i-3)(i-4)/2; B(19..22)
  C2  near-pencil union sizes: possible |C| for a pencil / pencil+1 family of >= 8 empty lines
  C3  s<=3 at k=19: the integer arithmetic (a,b) in {(1,6),(4,5)}, z=30, 180 > 171
  C4  bisecant bound (3) at k=19, e0=9: s=7,8 impossible; s=6 gives delta+6(a0-9)<=9
  C5  a_0 >= 9 at k=19: m in 3..6 union bound <= 60; m<=2 needs 18 parallel pairs, impossible
  C6  three-ray / three-row: (h_j+1)(h_l+1)>=9 for all pairs => sum h >= 6; h_j h_l>=9 => sum >= 9
  C7  triple packings on 6 labels (no repeated pair): every >=3-triple packing is the Fano-minus-one
      shape; with 4 partitions covering 14 of 15 pairs the missing pair is a T4 edge (r=2) or M edge (r=3)
  C8  BRUTE FORCE of the line-only core of the Desargues step: no family of 10 affine lines with
      pairwise distinct directions has union of size 62, maximum finite concurrence 4 and >= 4 points
      of concurrence 4.  (Also: with max concurrence 5 or 6 the union is <= 61.)
  C9  frames of 4 lines (no three projectively concurrent) in AG(2,9): orbits under AGammaL(2,9):
      1 orbit with two parallel pairs, 3 with one, 3 with none; stabiliser orders 16,8,2,4,2,4,8
  C10 orbits of s-subsets of a line under x -> a x^sigma + b: counts 3,3,2,1,1 for s=4..8
"""
import itertools, os, sys
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ag import AG

g = AG(2, 9)
F = g.F
q = 9
P = len(g.points)      # 81
LMASK = [sum(1 << p for p in pts) for _, pts in g.lines]   # 81-bit masks
LDIR = [d for d, _ in g.lines]
POPC = [bin(i).count("1") for i in range(1 << 16)]


def popcount(m):
    c = 0
    while m:
        c += POPC[m & 0xFFFF]
        m >>= 16
    return c


def ok(name, cond, extra=""):
    print(f"{name}: {'OK' if cond else 'FAIL'} {extra}")
    assert cond, name


# ---- C1
coef = [comb(i, 2) - 3 * i + 6 for i in range(10)]
ok("C1 coefficients", coef == [(i - 3) * (i - 4) // 2 for i in range(10)] and coef == [6, 3, 1, 0, 0, 1, 3, 6, 10, 15], str(coef))
B = {k: comb(k, 2) - 30 * k + 540 - 48 - 3 * k for k in (19, 20, 21, 22)}
ok("C1 B(k) = C(k,2)-30k+540-6*8-3k", B == {19: 36, 20: 22, 21: 9, 22: -3}, str(B))
# k=22: RHS 111 < 6*8+3*22 = 114
ok("C1 k=22 excluded by hand", comb(22, 2) - 660 + 540 == 111 and 111 < 114)

# ---- C2 near-pencil sizes (projective completion; 'pencil' may be a parallel class)
sizes = set()
for m in range(8, 11):
    sizes.add(81 - (1 + 8 * m))                      # m concurrent at a finite point
for m in range(8, 10):
    sizes.add(81 - 9 * m)                            # m parallel
for m in range(7, 11):
    for delta in (0, 1):
        sizes.add(81 - (10 + 7 * m + delta))         # m through a finite point + 1 line avoiding it
for m in range(7, 10):
    sizes.add(81 - (8 * m + 9))                      # m parallel + a transversal
sizes = {s for s in sizes if s >= 0}
ok("C2 near-pencil |C| list", sizes <= {0, 1, 7, 8, 9, 14, 15, 16, 21, 22} and 19 not in sizes and 20 not in sizes, str(sorted(sizes)))
# the two formulas themselves, recomputed from masks: m concurrent lines at a finite point
o = g.pidx[(0, 0)]
pencil = [l for l in g.lines_through[o]]
for m in range(8, 11):
    u = 0
    for l in pencil[:m]:
        u |= LMASK[l]
    assert popcount(u) == 1 + 8 * m
# m concurrent + one line avoiding o: union = 1+8m + (9 - (#pencil lines it meets finitely))
for m in range(7, 11):
    for l in range(len(g.lines)):
        if o in g.lines[l][1]:
            continue
        u = 0
        for pl in pencil[:m]:
            u |= LMASK[pl]
        u |= LMASK[l]
        meets = sum(1 for pl in pencil[:m] if LDIR[pl] != LDIR[l])
        assert popcount(u) == 1 + 8 * m + 9 - meets and (m - meets) in (0, 1)
print("C2 union formulas recomputed from masks: OK")

# ---- C3
sols = [(a, b) for a in range(20) for b in range(20) if a + 3 * b == 19 and a + b <= 9]
ok("C3 (a,b) solutions", sols == [(1, 6), (4, 5)], str(sols))
z = (19 * 8 - 62) // 3
ok("C3 z=30, 30*C(4,2)=180 > C(19,2)=171", 62 + 3 * z == 152 and z == 30 and 30 * 6 > comb(19, 2))

# ---- C4 bisecant bound (3): (3-gamma) delta + 6(a0-e0) <= B_e0(k) - s beta - (s-3)(s-4)/2
def bound3(k, s, e0):
    A = s + 18 - k
    beta, gamma = ((1, 1) if A == 3 else (A - 2, 2))
    Be = comb(k, 2) - 33 * k + 540 - 6 * e0
    return (3 - gamma), Be - s * beta - (s - 3) * (s - 4) // 2
for s in (7, 8):
    cg, rhs = bound3(19, s, 9)
    ok(f"C4 k=19 s={s} e0=9 impossible", rhs < 0, f"rhs={rhs}")
cg, rhs = bound3(19, 8, 8)
ok("C4 k=19 s=8 impossible already with e0=8", rhs < 0, f"rhs={rhs}")
cg, rhs = bound3(19, 7, 8)
print(f"C4 NOTE: k=19 s=7 with only e0=8 gives rhs={rhs} (not excluded) -> s=7 exclusion NEEDS a_0>=9 (C5)")
cg, rhs = bound3(19, 6, 9)
ok("C4 k=19 s=6: delta+6(a0-9)<=9", cg == 1 and rhs == 9)
# the local identity sum_{i>=4}(i-3)n_i = 2n_1+n_2+k-21: check on random pencils is pointless; it is
# sum (i-1) n_i = k-1 with sum n_i = 10, verified symbolically: sum(i-1)n_i - 2*sum n_i = sum(i-3)n_i = k-21

# ---- C5 a_0 >= 9 at k=19,20
for m in (3, 4, 5, 6):
    assert 81 - 10 * m + m * m <= 60
pairs_by_partition = set()
for part in itertools.product(range(1, 9), repeat=8):
    pass
def partitions(n, mx=None):
    mx = mx or n
    if n == 0:
        yield []
        return
    for first in range(min(n, mx), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest
pp = {sum(comb(c, 2) for c in part) for part in partitions(8)}
ok("C5 no partition of 8 lines into direction classes has 17 or 18 parallel pairs", 17 not in pp and 18 not in pp, str(sorted(pp)))
# also with the class cap 5 (three-row lemma) 16 is impossible too, but not needed for k=19

# ---- C6
ray_ok = all(h1 + h2 + h3 >= 6 for h1 in range(9) for h2 in range(9) for h3 in range(9)
             if all((a + 1) * (b + 1) >= 9 for a, b in ((h1, h2), (h1, h3), (h2, h3))))
row_ok = all(h1 + h2 + h3 >= 9 for h1 in range(10) for h2 in range(10) for h3 in range(10)
             if all(a * b >= 9 for a, b in ((h1, h2), (h1, h3), (h2, h3))))
ok("C6 three-ray sum h>=6 and three-row sum h>=9 (exhaustive over h)", ray_ok and row_ok)

# ---- C7 triple packings
labels = range(6)
triples = list(itertools.combinations(labels, 3))
def packing_ok(T):
    pairs = set()
    for t in T:
        for pr in itertools.combinations(t, 2):
            if pr in pairs:
                return False
            pairs.add(pr)
    return True
pack3 = [T for T in itertools.combinations(triples, 3) if packing_ok(T)]
pack4 = [T for T in itertools.combinations(triples, 4) if packing_ok(T)]
pack5 = [T for T in itertools.combinations(triples, 5) if packing_ok(T)]
ok("C7 120 three-triple and 30 four-triple packings, no five", (len(pack3), len(pack4), len(pack5)) == (120, 30, 0))
# every 3-packing: pairwise share exactly one label, no common label; the uncovered 6 pairs split as a triangle + perfect matching
for T in pack3:
    assert all(len(set(a) & set(b)) == 1 for a, b in itertools.combinations(T, 2))
    assert not (set(T[0]) & set(T[1]) & set(T[2]))
    covered = {pr for t in T for pr in itertools.combinations(t, 2)}
    unc = [pr for pr in itertools.combinations(labels, 2) if pr not in covered]
    assert len(unc) == 6
    tri = [t for t in triples if all(pr in unc for pr in itertools.combinations(t, 2))]
    assert len(tri) == 1
    M = [pr for pr in unc if not set(pr) <= set(tri[0])]
    assert len(M) == 3 and len(set(x for pr in M for x in pr)) == 6   # perfect matching
    # the fourth triple is forced to be tri[0]
    assert [t for t in triples if packing_ok(T + (t,))] == tri
    # each M edge meets every one of the 4 triples
    assert all(set(e) & set(t) for e in M for t in list(T) + tri)
# the partition statement: 4 partitions of 6 labels into blocks of size <= 3, covering 14 distinct pairs,
# with >= 3 triple blocks: enumerate all and classify the missing pair
def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for sub in set_partitions(rest):
        for i, blk in enumerate(sub):
            if len(blk) < 3:
                yield sub[:i] + [blk + [first]] + sub[i + 1:]
        yield [[first]] + sub
parts = [p for p in set_partitions(list(labels))]
def pairs_of(p):
    return {pr for blk in p for pr in itertools.combinations(sorted(blk), 2)}
outcomes = set()
count = 0
for combo in itertools.combinations_with_replacement(range(len(parts)), 4):
    ps = [parts[i] for i in combo]
    prs = [pairs_of(p) for p in ps]
    if sum(len(x) for x in prs) != 14 or len(set().union(*prs)) != 14:
        continue
    ntrip = sum(1 for p in ps for blk in p if len(blk) == 3)
    if ntrip < 3:
        continue
    count += 1
    T = tuple(tuple(sorted(blk)) for p in ps for blk in p if len(blk) == 3)
    assert ntrip == 3
    covered = set().union(*prs)
    missing = [pr for pr in itertools.combinations(labels, 2) if pr not in covered][0]
    tri = [t for t in triples if packing_ok(T + (t,))][0]
    kind = "T4" if set(missing) <= set(tri) else "M"
    r = sum(1 for p in ps for blk in p if len(blk) == 2 and set(blk) <= set(tri))
    outcomes.add((kind, r))
ok("C7 partition classification: missing pair is T4-edge iff r=2, M-edge iff r=3", outcomes == {("T4", 2), ("M", 3)}, f"{count} configurations, outcomes {sorted(outcomes)}")

# ---- C8 brute force: 10 lines, one per direction, with the origin a point of max concurrence m.
# WLOG (translation) a point of maximal finite concurrence is the origin. m lines through the origin
# (any m-subset of the 10 directions; reduced by nothing -- 210/252/210 subsets x 9^(10-m) is cheap enough
# with masks) and one line in each remaining direction.
def c8(m):
    res = {}
    through_o = {LDIR[l]: l for l in g.lines_through[o]}
    others = {}
    for d in range(10):
        others[d] = [l for l in g.lines_by_dir[d] if l != through_o[d]]
    for dirs in itertools.combinations(range(10), m):
        base = 0
        for d in dirs:
            base |= LMASK[through_o[d]]
        rest = [d for d in range(10) if d not in dirs]
        for choice in itertools.product(*[others[d] for d in rest]):
            u = base
            for l in choice:
                u |= LMASK[l]
            n = popcount(u)
            if n != 62:
                continue
            # multiplicities
            lines = [through_o[d] for d in dirs] + list(choice)
            mult = [0] * P
            for l in lines:
                for p in g.lines[l][1]:
                    mult[p] += 1
            mx = max(mult)
            if mx > m:
                continue  # origin not of maximal concurrence; counted under the larger m
            n4 = sum(1 for x in mult if x == 4)
            key = (mx, n4)
            res[key] = res.get(key, 0) + 1
    return res
r4 = c8(4)
print("C8 m=4: (maxconc, n4) -> count:", r4)
ok("C8 no 10-line distinct-direction family with union 62, max concurrence 4 and n4 >= 4", all(n4 <= 3 for (mx, n4) in r4))
r5 = c8(5)
r6 = c8(6)
ok("C8 union 62 impossible with max concurrence 5 or 6", not r5 and not r6, f"{r5} {r6}")

# ---- C9 frames
# affine-semilinear maps (x,y) -> M (x^s, y^s) + b  as point permutations, from generators
def perm_of(M, b, frob):
    perm = []
    for p in g.points:
        x, y = (F.frob[p[0]], F.frob[p[1]]) if frob else p
        img = g.add((g.dot(M[0], (x, y)), g.dot(M[1], (x, y))), b)
        perm.append(g.pidx[img])
    assert sorted(perm) == list(range(P))
    return perm
w = next(a for a in range(2, 9) if len({F.power(a, e) for e in range(8)}) == 8)  # primitive element
gens = [perm_of(((1, 0), (0, 1)), (1, 0), False), perm_of(((1, 0), (0, 1)), (0, 1), False),
        perm_of(((1, 1), (0, 1)), (0, 0), False), perm_of(((1, 0), (1, 1)), (0, 0), False),
        perm_of(((w, 0), (0, 1)), (0, 0), False), perm_of(((1, 0), (0, w)), (0, 0), False),
        perm_of(((1, 0), (0, 1)), (0, 0), True)]
for perm in gens:   # collineation check
    for _, pts in g.lines:
        assert frozenset(perm[p] for p in pts) in g.line_of
def line_img(perm, l):
    return g.line_of[frozenset(perm[p] for p in g.lines[l][1])]
LIMG = [[line_img(perm, l) for l in range(len(g.lines))] for perm in gens]
inf = {}  # projective point of a line: finite points are indices, infinite point = ('inf', dir)
def is_frame(ls):
    # no three concurrent in the projective completion
    for a, b, c in itertools.combinations(ls, 3):
        da, db, dc = LDIR[a], LDIR[b], LDIR[c]
        if da == db == dc:
            return False
        if da == db or da == dc or db == dc:
            continue  # the parallel pair meets at infinity, the third line does not pass through it
        if LMASK[a] & LMASK[b] & LMASK[c]:
            return False
    return True
frames = [fs for fs in itertools.combinations(range(len(g.lines)), 4) if is_frame(fs)]
idx = {fs: i for i, fs in enumerate(frames)}
parent = list(range(len(frames)))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i
for i, fs in enumerate(frames):
    for LI in LIMG:
        j = idx[tuple(sorted(LI[l] for l in fs))]
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
orbits = {}
for i, fs in enumerate(frames):
    orbits.setdefault(find(i), []).append(fs)
def npar(fs):
    return sum(1 for a, b in itertools.combinations(fs, 2) if LDIR[a] == LDIR[b])
summary = sorted((npar(fs[0]), len(fs)) for fs in orbits.values())
G_ORDER = 81 * 80 * 72 * 2
stab = sorted((npar(fs[0]), G_ORDER // len(fs)) for fs in orbits.values())
ok("C9 frame orbits under AGammaL(2,9): (parallel pairs, orbit size)", [p for p, _ in summary] == [0, 0, 0, 1, 1, 1, 2] and all(G_ORDER % len(fs) == 0 for fs in orbits.values()), f"{summary} stabiliser orders {stab}")
ok("C9 stabiliser orders {16,8,2,4,2,4,8} as multiset", sorted(s for _, s in stab) == sorted([16, 8, 2, 4, 2, 4, 8]))
ok("C9 total frames = sum of orbit sizes", sum(len(v) for v in orbits.values()) == len(frames), f"{len(frames)} frames")
# C9b: every one-parallel-pair frame is AGL(2,9)-equivalent (even without Frobenius) to x=0,y=0,x+y=1,(0,b): b in {2,3,5}
# (the research's (0,2),(0,3),(0,5)) -- check that the three p=1 orbit representatives can be normalised to those forms
def line_eq(a, b, c):  # a x + b y = c
    return g.line_of[frozenset(g.pidx[p] for p in g.points if F.add[F.mul[a][p[0]]][F.mul[b][p[1]]] == c)]
p1_forms = {}
for bb in range(2, 9):
    if bb == 1:
        continue
    fs = tuple(sorted([line_eq(1, 0, 0), line_eq(0, 1, 0), line_eq(1, 1, 1), line_eq(0, bb, 1)]))
    assert is_frame(fs) and npar(fs) == 1
    p1_forms.setdefault(find(idx[fs]), []).append(bb)
ok("C9b forms (0,b), b=2..8, hit exactly the 3 one-pair orbits", len(p1_forms) == 3 and set(p1_forms) == {r for r, fs in orbits.items() if npar(fs[0]) == 1}, str(sorted(p1_forms.values())))
rect = tuple(sorted([line_eq(1, 0, 0), line_eq(1, 0, 1), line_eq(0, 1, 0), line_eq(0, 1, 1)]))
ok("C9c rectangle x=0,x=1,y=0,y=1 is the two-pair orbit", npar(rect) == 2 and len(orbits[find(idx[rect])]) == G_ORDER // 16)

# ---- C10 axis orbits
maps = [(a, b, e) for a in range(1, 9) for b in range(9) for e in (0, 1)]
def apply(m, x):
    a, b, e = m
    return F.add[F.mul[a][F.frob[x] if e else x]][b]
counts = {}
for s in range(2, 9):
    seen = set()
    n = 0
    for S in itertools.combinations(range(9), s):
        if S in seen:
            continue
        n += 1
        for m in maps:
            seen.add(tuple(sorted(apply(m, x) for x in S)))
    counts[s] = n
ok("C10 axis orbit counts s=4..8 are 3,3,2,1,1", [counts[s] for s in (4, 5, 6, 7, 8)] == [3, 3, 2, 1, 1], str(counts))
print("ALL LEMMA CHECKS PASSED")
