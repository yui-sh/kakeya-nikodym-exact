#!/usr/bin/env python3
"""Fresh encoder for Nikodym lower bounds in AG(n,q) (planar q <= 8 and AG(3,3)).

Definition (every-point variant): N is Nikodym iff every point x (in N or not) lies on a line L with
L \\ {x} contained in N.  In terms of the complement C = AG \\ N: every point x has a line L through x
with (L \\ {x}) disjoint from C ('free line': a tangent if x in C, a passant otherwise).

Formula  nik(n,q,K,case):  'C is a Nikodym complement with |C| >= K' plus the case constraints.
Variables c_p (p in C) and w_{x,L} for every incident point-line pair:
   w_{x,L} -> not c_p  (p in L, p != x);    OR_{L through x} w_{x,L}  for every point x;
   sum_p c_p >= K (sequential counter).
UNSAT for every case of an exhaustive split  =>  max |C| <= K-1  =>  Nik(n,q) >= q^n - K + 1.

Cases (symmetry breaking; the group is AGL(n,q), whose elements map lines to lines, preserve
|C| and the Nikodym property, so every C may be replaced by an image g(C)):
  none   : no symmetry breaking (exhaustive by itself).
  flag_h, flag_d  (n = 2; exhaustive together for |C| >= 2):
     AGL(2,q) is transitive on flags (point, line through it), so WLOG the origin O is in C and its
     free line (tangent) is the x-axis: c(t,0) = 0 for t != 0.  Any second point of C is off the
     tangent; the flag stabiliser {(x,y) -> (ax+by, dy)} is transitive on points off the x-axis, so
     WLOG e2 = (0,1) in C.  The tangent at e2 has direction (1,v); the residual group
     {(x,y) -> (ax+by, y)} sends (1,v) to (1, v/(a+bv)), so v = 0 (flag_h: the line y = 1 is free)
     or v != 0 WLOG v = 1 (flag_d: the line {(t, 1+t)} is free) are the two orbits.
  tan2_e1, tan2_e1e2, tan2_e3  (n = 3; exhaustive together for |C| >= 2): same first two steps
     (O in C, tangent x1-axis, e2 in C); the residual group {M : M e1 = a e1, M e2 = e2, M e3 = b}
     has three orbits on the possible directions d of the tangent at e2: d = e1, d = (1,1,0), d = e3.
  The transitivity claims are checked computationally by symmetry_check.py (explicit affine maps).
Cube splitting (--cube DEPTH): the witness clause of a point x is OR_L w_{x,L}; replacing the formula
by the (q+1) formulas 'formula + unit w_{x,L}' is a plain clause split, so the cubes are exhaustive
by construction; --cover writes the CNF '(witness clauses) and (negation of every cube)' whose UNSAT
proof certifies the coverage mechanically.
usage: enc_nikodym.py --n 2 --q 7 --K 14 --case flag_h --out f.cnf
       enc_nikodym.py ... --model model.txt   (decode: print C as JSON for check_nikodym.py)
       enc_nikodym.py ... --cube 2 --outdir DIR [--cover DIR/cover.cnf]   (write all cubes)
"""
import argparse, itertools, json, os, sys
from ag import AG
import card


def base_formula(g, K):
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
    cc, top = card.atleast([c(p) for p in range(P)], K, top)
    cls += cc
    return top, cls, w, witness_clauses


def case_units(g, case):
    n, q = g.n, g.q
    pidx = g.pidx
    zero = tuple([0] * n)
    e1 = tuple(int(i == 0) for i in range(n))
    e2 = tuple(int(i == 1) for i in range(n))
    units = []
    if case == "none":
        return units
    units.append(pidx[zero] + 1)
    for t in range(1, q):  # x1-axis is the tangent at the origin
        units.append(-(pidx[g.scale(t, e1)] + 1))
    units.append(pidx[e2] + 1)
    d = {"flag_h": e1, "flag_d": tuple([1, 1] + [0] * (n - 2)), "tan2_e1": e1, "tan2_e1e2": (1, 1, 0), "tan2_e3": (0, 0, 1)}[case]
    for t in range(1, q):  # tangent at e2 with direction d
        units.append(-(pidx[g.add(e2, g.scale(t, d))] + 1))
    return units


def cube_list(g, w, depth, order=None):
    """cubes = choices of one witness line for each of the first `depth` points in `order`"""
    order = order or list(range(len(g.points)))[:depth]
    choices = [[w[(x, l)] for l in g.lines_through[x]] for x in order]
    return [list(cube) for cube in itertools.product(*choices)]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--case", default="none")
    ap.add_argument("--out")
    ap.add_argument("--model")
    ap.add_argument("--cube", type=int, default=0)
    ap.add_argument("--cube-points", default=None, help="comma-separated point indices to split on")
    ap.add_argument("--outdir")
    ap.add_argument("--cover")
    a = ap.parse_args()
    g = AG(a.n, a.q)
    nv, cls, w, wcl = base_formula(g, a.K)
    units = case_units(g, a.case)
    cls = cls + [[u] for u in units]
    header = f"nikodym complement n={a.n} q={a.q} |C|>={a.K} case={a.case}"
    if a.model:
        lits = set()
        for line in open(a.model):
            if line.startswith("v"):
                lits.update(int(t) for t in line.split()[1:])
        C = [list(g.points[p]) for p in range(len(g.points)) if (p + 1) in lits]
        print(json.dumps({"n": a.n, "q": a.q, "C": C, "case": a.case, "K": a.K}))
    elif a.cube:
        order = [int(t) for t in a.cube_points.split(",")] if a.cube_points else None
        cubes = cube_list(g, w, a.cube, order)
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
        print(json.dumps({"n": a.n, "q": a.q, "K": a.K, "case": a.case, "nvars": nv, "nclauses": len(cls)}))
