#!/usr/bin/env python3
"""Fresh encoder for the Nikodym bound in AG(2,9) (and any planar q): max-secant split.

Definition (every-point variant, as in enc_nikodym.py): N is Nikodym iff every point x lies on a line L
with L \\ {x} contained in N.  C = AG(2,q) \\ N.  Formula 'C is a Nikodym complement with |C| = K':
   variables c_p (p in C), w_{x,L} (L is the free line at x);
   w_{x,L} -> not c_p (p in L, p != x);   OR_L w_{x,L} for every x;   sum c_p = K (sequential counters).
Exactly K (not >= K) is sound by heredity: a subset of a Nikodym complement is a Nikodym complement
(the free line at x avoids C \\ {x}, hence avoids any subset), so if some complement has >= K points
one of its K-subsets is a complement with exactly K points.

Case split (every case is one CNF; UNSAT of all cases => no complement with |C| = K):
  Let C be a complement with |C| = K >= 2 and let s = max_L |L cap C| (2 <= s <= q, since two points
  of C span a line).  Fix a line A attaining s.
  (1) AGL(2,q) is transitive on lines: WLOG A is the x-axis y = 0.
  (2) The maps (x,y) -> (a x^sigma + b, y^sigma), a != 0, fix the axis setwise and act on it as
      AGammaL(1,q); WLOG A cap C is the orbit representative S (chosen with 0,1 in S, possible since
      AGL(1,q) is 2-transitive).  Encoded as units: c(t,0) = [t in S].
  (3) (0,0) in C has a free line (tangent) T0 != A, of the form x = b y; the shear (x,y) -> (x - b y, y)
      fixes A pointwise and maps T0 to x = 0.  Encoded: c(0,t) = 0 for t != 0.
  (4) (1,0) in C has a tangent T1 != A, of the form x = 1 + b y; the scaling (x,y) -> (x, d y) fixes A
      pointwise and x = 0 setwise and maps T1 to x = 1 + (b/d) y, so WLOG b = 0 (tan=h: c(1,t) = 0,
      t != 0) or b = 1 (tan=d: c(1+t,t) = 0, t != 0).
  (5) every line has at most s points of C (at-most-s counters on all lines; preserved by every step,
      because the maps are collineations).
  All transitivity/orbit claims are checked by explicit maps in symmetry_check_nik29.py.
  Cases: (s, orbit representative, tan) for s = 2..q.  Cubes (--cube) split the witness clause of
  chosen points exactly as enc_nikodym.py does (plain clause split, --cover certifies exhaustiveness).
usage: enc_nik29.py --q 9 --K 19 --s 4 --orbit 0 --tan h --out f.cnf
       enc_nik29.py --list                       (print the orbit representatives per s)
       enc_nik29.py ... --model model.txt         (decode a SAT model to JSON for check_nikodym.py)
       enc_nik29.py ... --cube D --cube-points i,j,... --outdir DIR [--cover DIR/cover.cnf]
"""
import argparse, itertools, json, os, sys
from ag import AG
import card


def axis_orbit_reps(F):
    """orbit representatives (containing 0 and 1) of s-subsets of F_q under x -> a x^sigma + b"""
    q = F.q
    maps = [(a, b, e) for a in range(1, q) for b in range(q) for e in range(F.k)]

    def frob(x, e):
        for _ in range(e):
            x = F.frob[x]
        return x

    reps = {}
    for s in range(2, q + 1):
        seen = set()
        reps[s] = []
        for S in itertools.combinations(range(q), s):
            if S in seen:
                continue
            orbit = {tuple(sorted(F.add[F.mul[a][frob(x, e)]][b] for x in S)) for (a, b, e) in maps}
            seen |= orbit
            rep = min(T for T in orbit if 0 in T and 1 in T)
            reps[s].append((rep, len(orbit)))
    return reps


def base_formula(g, K, s):
    P = len(g.points)
    c = lambda p: p + 1
    top = P
    w = {}
    cls = []
    witness_clauses = []
    for x in range(P):
        lits = []
        for l in g.lines_through[x]:
            top += 1
            w[(x, l)] = top
            lits.append(top)
            for p in g.lines[l][1]:
                if p != x:
                    cls.append([-top, -c(p)])
        cls.append(lits)
        witness_clauses.append(lits)
    cc, top = card.equals([c(p) for p in range(P)], K, top)
    cls += cc
    for _, pts in g.lines:
        cc, top = card.atmost([c(p) for p in pts], s, top)
        cls += cc
    return top, cls, w, witness_clauses


def case_units(g, S, tan):
    q = g.q
    pidx = g.pidx
    units = []
    for t in range(q):  # axis y = 0 meets C exactly in S
        units.append((1 if t in S else -1) * (pidx[(t, 0)] + 1))
    for t in range(1, q):  # tangent at (0,0) is x = 0
        units.append(-(pidx[(0, t)] + 1))
    for t in range(1, q):  # tangent at (1,0) is x = 1 (h) or x = 1 + y (d)
        p = (1, t) if tan == "h" else (g.F.add[1][t], t)
        units.append(-(pidx[p] + 1))
    return units


def cube_list(g, w, order):
    choices = [[w[(x, l)] for l in g.lines_through[x]] for x in order]
    return [list(cube) for cube in itertools.product(*choices)]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=9)
    ap.add_argument("--K", type=int, default=19)
    ap.add_argument("--s", type=int)
    ap.add_argument("--orbit", type=int)
    ap.add_argument("--tan", choices=["h", "d"])
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--out")
    ap.add_argument("--model")
    ap.add_argument("--cube", type=int, default=0)
    ap.add_argument("--cube-points", default=None)
    ap.add_argument("--outdir")
    ap.add_argument("--cover")
    a = ap.parse_args()
    g = AG(2, a.q)
    reps = axis_orbit_reps(g.F)
    if a.list:
        for s in reps:
            print(s, reps[s], "sum of orbit sizes", sum(n for _, n in reps[s]), "=", len(list(itertools.combinations(range(a.q), s))))
        sys.exit(0)
    S = reps[a.s][a.orbit][0]
    nv, cls, w, wcl = base_formula(g, a.K, a.s)
    units = case_units(g, S, a.tan)
    cls = cls + [[u] for u in units]
    header = f"nikodym complement n=2 q={a.q} |C|={a.K} max-secant s={a.s} axis S={list(S)} tan={a.tan}"
    if a.model:
        lits = set()
        for line in open(a.model):
            if line.startswith("v"):
                lits.update(int(t) for t in line.split()[1:])
        C = [list(g.points[p]) for p in range(len(g.points)) if (p + 1) in lits]
        print(json.dumps({"n": 2, "q": a.q, "C": C, "s": a.s, "S": list(S), "tan": a.tan, "K": a.K}))
    elif a.cube:
        order = [int(t) for t in a.cube_points.split(",")]
        assert len(order) == a.cube
        cubes = cube_list(g, w, order)
        os.makedirs(a.outdir, exist_ok=True)
        names = []
        for i, cube in enumerate(cubes):
            name = f"cube{i}"
            card.write_dimacs(os.path.join(a.outdir, name + ".cnf"), nv, cls + [[u] for u in cube], [header + f" cube {cube}"])
            names.append({"name": name, "cube": cube})
        if a.cover:
            card.write_dimacs(a.cover, nv, wcl + [[-u for u in cube] for cube in cubes], [header + " COVER CHECK: witness clauses and the negation of every cube"])
        json.dump(names, open(os.path.join(a.outdir, "cubes.json"), "w"))
        print(json.dumps({"cubes": len(cubes), "nvars": nv, "nclauses": len(cls)}))
    else:
        card.write_dimacs(a.out, nv, cls, [header])
        print(json.dumps({"q": a.q, "K": a.K, "s": a.s, "S": list(S), "tan": a.tan, "nvars": nv, "nclauses": len(cls)}))
