#!/usr/bin/env python3
"""Fresh encoder and case generator for Nik(3,4) = 22 (no Nikodym set of size <= 21 in AG(3,4)).

Definition: N is Nikodym iff every point x has a line L through x with L \\ {x} inside N.
Consequences used (all elementary, re-derived in README):
 (F) every point of N lies on a FULL line (4 points of N): take its own witness line.
 (T) every point x outside N has a witness line with exactly 3 N-points; distinct x give distinct
     lines; so t >= 64 - n where t = number of 3-point lines, n = |N|.
 (L) u + 3t + 6f = C(n,2) (u, t, f = lines with 2, 3, 4 N-points), so 6f <= C(n,2) - 3(64-n).
 (E) exterior-pair lemma for a plane P with m = |N ∩ P|, r = n - m: the 48-r complement points
     outside P have distinct witness lines each using a pair of N-points outside P; the full lines
     not inside P cover N \\ P (a with 3 points outside, b with 4; 3a+4b >= r) and use 3a+6b further
     pairs; the b_C complement points of P without a 3-point witness line inside P use 3 more pairs
     each; every N-point of P on no full line of P needs a full line leaving P (so a >= b_N).  Hence
        C(r,2) >= (48-r) + 3 b_C + min{3a+6b : a >= b_N, 3a+4b >= r}        (*)
     and in particular r >= 11, i.e. |N ∩ P| <= n - 11 for every plane.
 (S) if the full lines are pairwise disjoint they partition N, so 4 | n: n = 20 with 5 full lines.
Case split for n <= 21 (every case is a CNF; UNSAT of all cases excludes n <= 21):
  intersect : two full lines meet.  WLOG (AGL(3,4) is transitive on ordered pairs of intersecting
              lines, checked by exhibiting explicit maps) they are the x- and y-axes in P = {z=0}.
              The plane section S = N ∩ P is a 16-bit mask admitted by (*) (enumerated here over all
              2^16 masks) and containing both axes, taken up to the group G36 of maps
              (x,y,z) -> (a x^s, b y^s, z^s) and (b y^s, a x^s, z^s), a,b in F4*, s in {id, Frobenius}
              (semilinear collineations fixing P and the pair of axes; orbits computed here).
              For each representative S: if some point p of S lies on no full line inside P, a full
              line through p leaves P and a shear (x,y,z) -> (x+Az, y+Bz, cz), which fixes P pointwise,
              makes it the vertical line through p ('anchor').  Otherwise N is not inside P (a point
              off P needs 3 N-points on a line through it), so some full line leaves P: either through
              a point p of S (shear -> vertical line through p; p up to the stabiliser of S in G36) or
              parallel to P (shear + z-scaling -> a line in {z=1} through (0,0,1), direction up to the
              stabiliser of S).  Every map used fixes P pointwise or stabilises S, so S is unchanged.
  parallel  : full lines pairwise disjoint and two of them parallel: WLOG the x-axis and {(t,1,0)}
              (AGL transitive on ordered pairs of distinct parallel lines, checked).
  skew      : full lines pairwise disjoint, none parallel: WLOG the x-axis and {(0,t,1)} (checked).
CNF (variables n_p, f_L <-> all 4 points of L in N, w_{x,L} -> points of L\\{x} in N; OR_L w_{x,L};
n_p -> OR f_L; sum n_p <= k; cuts (E) per plane, (L) on f; for k = 20 the two facts verified on the
admitted masks: no plane contains two parallel full lines or three full lines).
usage: enc_nik34.py --spec DIR  |  --case NAME --k K --out f.cnf  |  --list  |  --model m.txt --case NAME --k K
"""
import argparse, itertools, json, os, sys
from ag import AG
import card


def plane_geometry():
    """AG(2,4) in the mask convention bit = 4x+y; returns (lines as masks, directions as line sets)"""
    g2 = AG(2, 4)
    idx = lambda p: 4 * p[0] + p[1]
    lines = [sum(1 << idx(g2.points[p]) for p in pts) for _, pts in g2.lines]
    return g2, lines


def admitted_masks(n, lines):
    """all 16-bit masks S = N ∩ P passing (*) for |N| = n"""
    mincost = {}
    for r in range(n + 1):
        for bn in range(17):
            mincost[(r, bn)] = min(3 * a + 6 * b for a in range(0, 22) for b in range(0, 22) if a >= bn and 3 * a + 4 * b >= r)
    ok = []
    for S in range(1 << 16):
        m = bin(S).count("1")
        r = n - m
        if r < 11:
            continue
        full = [L for L in lines if L & S == L]
        covered = 0
        for L in full:
            covered |= L
        bN = bin(S & ~covered).count("1")
        witnessed = 0  # complement points of P with a 3-point witness line inside P
        for L in lines:
            inside = L & S
            if bin(inside).count("1") == 3:
                witnessed |= L & ~S
        bC = bin((0xFFFF & ~S) & ~witnessed).count("1")
        if r * (r - 1) // 2 >= (48 - r) + 3 * bC + mincost[(r, bN)]:
            ok.append(S)
    return ok


def group36(g2):
    """the 36 maps (x,y) -> (a x^s, b y^s) / swapped, as permutations of the 16 plane bits"""
    F = g2.F
    perms = []
    for a in range(1, 4):
        for b in range(1, 4):
            for swap in (0, 1):
                for s in (0, 1):
                    perm = []
                    for x in range(4):
                        for y in range(4):
                            xx, yy = (F.frob[x], F.frob[y]) if s else (x, y)
                            if swap:
                                xx, yy = yy, xx
                            perm.append(4 * F.mul[a][xx] + F.mul[b][yy])
                    assert sorted(perm) == list(range(16))
                    perms.append(perm)
    assert len(perms) == 36
    return perms


def act(mask, perm):
    return sum(1 << perm[i] for i in range(16) if mask >> i & 1)


def orbit_reps(masks, perms):
    todo = set(masks)
    reps = []
    while todo:
        m0 = min(todo)
        orb = {act(m0, p) for p in perms}
        assert orb <= todo, "group does not preserve the admitted set"
        todo -= orb
        reps.append((m0, len(orb)))
    return reps


AXES = sum(1 << (4 * x) for x in range(4)) | sum(1 << y for y in range(4))  # x-axis bits 4x+0, y-axis bits 0+y


def cases_for(n, g2, lines, perms):
    adm = admitted_masks(n, lines)
    with_axes = [S for S in adm if S & AXES == AXES]
    reps = orbit_reps(with_axes, perms)
    plane_dirs = sorted({L for L in lines if L & 1})  # the 5 lines through the origin = directions
    cases = []
    for S, orb in reps:
        full = [L for L in lines if L & S == L]
        covered = 0
        for L in full:
            covered |= L
        unsupported = [i for i in range(16) if S >> i & 1 and not covered >> i & 1]
        if unsupported:
            p = unsupported[0]
            cases.append(dict(name=f"n{n}_m{S}_anchor", n=n, mask=S, kind="vertical", point=[p // 4, p % 4], anchor=True))
            continue
        stab = [p for p in perms if act(S, p) == S]
        # orbits of the stabiliser on the points of S
        pts = [i for i in range(16) if S >> i & 1]
        seen = set()
        for p in pts:
            if p in seen:
                continue
            orb = {q[p] for q in stab}
            seen |= orb
            cases.append(dict(name=f"n{n}_m{S}_v{p}", n=n, mask=S, kind="vertical", point=[p // 4, p % 4], anchor=False))
        seen = set()
        for d in plane_dirs:
            if d in seen:
                continue
            orb = {act(d, q) for q in stab}
            seen |= orb
            cases.append(dict(name=f"n{n}_m{S}_d{d}", n=n, mask=S, kind="parallel", dirmask=d))
    return adm, reps, cases


def check_transitivity(g):
    """AGL(3,4) moves any ordered pair of intersecting / parallel / skew lines to the base pair"""
    O = (0, 0, 0)
    base = {"intersect": ((O, (1, 0, 0)), (O, (0, 1, 0))), "parallel": ((O, (1, 0, 0)), ((0, 1, 0), (1, 0, 0))), "skew": ((O, (1, 0, 0)), ((0, 0, 1), (0, 1, 0)))}
    counts = {}
    for kind, ((p1, d1), (p2, d2)) in base.items():
        L1 = g.line_through(p1, d1)
        L2 = g.line_through(p2, d2)
        cnt = 0
        for la in range(len(g.lines)):
            for lb in range(len(g.lines)):
                if la == lb:
                    continue
                A, B = set(g.lines[la][1]), set(g.lines[lb][1])
                da, db = g.dirs[g.lines[la][0]], g.dirs[g.lines[lb][0]]
                if A & B:
                    k = "intersect"
                elif da == db:
                    k = "parallel"
                else:
                    k = "skew"
                if k != kind:
                    continue
                # construct an affine map sending (L1, L2) to (la, lb): M e1 = da, M e2 chosen per kind
                pa = g.points[min(A)]
                if k == "intersect":
                    meet = g.points[next(iter(A & B))]
                    cols = [da, db]
                    b = meet
                elif k == "parallel":
                    pb = g.points[min(B)]
                    cols = [da, g.sub(pb, pa)]
                    b = pa
                else:
                    pb = g.points[min(B)]
                    cols = [da, db]
                    b = pa
                found = None
                for c3 in g.points:
                    M = [[cols[0][r], cols[1][r], c3[r]] for r in range(3)]
                    if k == "skew":
                        # third column must send e3 to pb - pa  (then (0,t,1) -> pa + t db + (pb-pa) = pb + t db)
                        c3 = g.sub(pb, pa)
                        M = [[cols[0][r], cols[1][r], c3[r]] for r in range(3)]
                    perm = [g.pidx[g.add(tuple(g.dot(M[r], p) for r in range(3)), b)] for p in g.points]
                    if sorted(perm) == list(range(64)):
                        found = perm
                        break
                assert found is not None, (kind, la, lb)
                img = lambda l: g.line_of[frozenset(found[p] for p in g.lines[l][1])]
                assert img(L1) == la and img(L2) == lb, (kind, la, lb)
                cnt += 1
        counts[kind] = cnt
    return counts


def build_cnf(g, k, case):
    P, NL = 64, len(g.lines)
    nvar = lambda p: p + 1
    f = lambda l: P + l + 1
    top = P + NL
    cls = []
    through = [[] for _ in range(P)]
    for l, (_, pts) in enumerate(g.lines):
        for p in pts:
            cls.append([-f(l), nvar(p)])
            through[p].append(l)
        cls.append([f(l)] + [-nvar(p) for p in pts])
    for x in range(P):
        wl = []
        for l in through[x]:
            top += 1
            wl.append(top)
            for p in g.lines[l][1]:
                if p != x:
                    cls.append([-top, nvar(p)])
        cls.append(wl)
        cls.append([-nvar(x)] + [f(l) for l in through[x]])  # (F)
    c, top = card.atmost([nvar(p) for p in range(P)], k, top)
    cls += c
    fmax = (k * (k - 1) // 2 - 3 * (64 - k)) // 6  # (L)
    c, top = card.atmost([f(l) for l in range(NL)], fmax, top)
    cls += c
    for (_, _, pts, lines_in) in g.planes:  # (E): at most k-11 N-points per plane
        c, top = card.atmost([nvar(p) for p in sorted(pts)], k - 11, top)
        cls += c
        if k == 20:
            c, top = card.atmost([f(l) for l in lines_in], 2, top)
            cls += c
    if k == 20:
        for ls in g.lines_by_dir:
            c, top = card.atmost([f(l) for l in ls], 1, top)
            cls += c
    kind = case["kind"]
    xaxis = [(t, 0, 0) for t in range(4)]
    if kind in ("vertical", "parallel"):
        for p in xaxis + [(0, t, 0) for t in range(4)]:
            cls.append([nvar(g.pidx[p])])
        S = case["mask"]
        for x in range(4):
            for y in range(4):
                v = nvar(g.pidx[(x, y, 0)])
                cls.append([v] if S >> (4 * x + y) & 1 else [-v])
        if kind == "vertical":
            x, y = case["point"]
            for z in range(4):
                cls.append([nvar(g.pidx[(x, y, z)])])
        else:
            d = case["dirmask"]
            bits = [i for i in range(16) if d >> i & 1 and i != 0]
            dx, dy = bits[0] // 4, bits[0] % 4
            for t in range(4):
                cls.append([nvar(g.pidx[(g.F.mul[t][dx], g.F.mul[t][dy], 1)])])
    else:  # parallel / skew: all full lines pairwise disjoint
        second = [(t, 1, 0) for t in range(4)] if kind == "parallel_pair" else [(0, t, 1) for t in range(4)]
        for p in xaxis + second:
            cls.append([nvar(g.pidx[p])])
        for p in range(P):
            for a, b in itertools.combinations(through[p], 2):
                cls.append([-f(a), -f(b)])
        if kind == "skew_pair":
            for ls in g.lines_by_dir:
                c, top = card.atmost([f(l) for l in ls], 1, top)
                cls += c
    return top, cls


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--case")
    ap.add_argument("--k", type=int)
    ap.add_argument("--out")
    ap.add_argument("--model")
    a = ap.parse_args()
    g2, lines = plane_geometry()
    perms = group36(g2)
    allcases = {}
    summary = {}
    for n in (20, 21):
        adm, reps, cases = cases_for(n, g2, lines, perms)
        if n == 20:  # facts used as cuts for k = 20, verified on every admitted mask
            par = [(A, B) for A, B in itertools.combinations(lines, 2) if not A & B]
            for S in adm:
                full = [L for L in lines if L & S == L]
                assert len(full) <= 2 and not any(A in full and B in full for A, B in par)
        adm20 = set(admitted_masks(20, lines)) if n == 21 else None
        if n == 21:
            assert adm20 <= set(adm), "n=20 admitted masks must be admitted for n=21"
        summary[n] = dict(admitted=len(adm), with_axes=sum(1 for S in adm if S & AXES == AXES), reps=reps, n_cases=len(cases))
        for c in cases:
            allcases[c["name"]] = c
    # positive control: the 22-point Baer cone {0} u {t v : v in {0,1}^3 \ 0, t in F4*} has full lines
    # through the origin, so it lies in the 'intersect' case with plane section mask 38207 (axes + diagonal)
    # and a vertical full line through the origin; its orbit representative under G36 is used.
    adm22, reps22, _ = cases_for(22, g2, lines, perms)
    rep22 = min(act(38207, p) for p in perms)
    assert rep22 in dict(reps22), "Baer cone plane section must be admitted at n=22"
    allcases["control22_cone"] = dict(name="control22_cone", n=22, mask=rep22, kind="vertical", point=[0, 0], anchor=False)
    for n in (20, 21):
        for kind in ("parallel_pair", "skew_pair"):
            name = f"n{n}_{kind}"
            allcases[name] = dict(name=name, n=n, kind=kind)
    if a.list:
        for n in (20, 21):
            s = summary[n]
            print(f"n={n}: {s['admitted']} admitted plane masks, {s['with_axes']} containing both axes, {len(s['reps'])} orbits under G36 "
                  f"{[(m, sz) for m, sz in s['reps']]}; {s['n_cases']} intersect-cases + 2 disjoint-cases")
        for c in allcases.values():
            print("  ", c["name"], {k: v for k, v in c.items() if k != "name"})
    elif a.spec:
        g = AG(3, 4)
        print("transitivity checks:", check_transitivity(g), file=sys.stderr)
        spec = [{"name": c["name"], "gen": ["python3", "enc_nik34.py", "--case", c["name"], "--k", str(c["n"]), "--out", f"{a.spec}/{c['name']}.cnf"]} for c in allcases.values()]
        print(json.dumps(spec, indent=0))
    else:
        g = AG(3, 4)
        case = allcases[a.case]
        nv, cls = build_cnf(g, a.k, case)
        if a.model:
            lits = set()
            for line in open(a.model):
                if line.startswith("v"):
                    lits.update(int(t) for t in line.split()[1:])
            N = [list(g.points[p]) for p in range(64) if (p + 1) in lits]
            print(json.dumps({"n": 3, "q": 4, "N": N, "case": a.case, "k": a.k}))
        else:
            card.write_dimacs(a.out, nv, cls, [f"Nik(3,4) case {a.case} |N|<={a.k}: {json.dumps(case)}"])
            print(json.dumps({"case": a.case, "k": a.k, "nvars": nv, "nclauses": len(cls)}))
