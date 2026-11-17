"""Session-I WLOG (3 fixed witness lines, 8 subcases meet/skew x {00,10,01,11}; proof in WRITEUP.md §3) as plain CNF,
without the degree reduction: independent second proof path for K(3,4) >= 28.
Usage: wlog3_cnf.py --q 4 --ub 27 --case meet --sub 00 --cnf out/w3_meet_00.cnf"""
import argparse, json
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from geom import AG

ap = argparse.ArgumentParser()
ap.add_argument("--q", type=int, default=4); ap.add_argument("--ub", type=int, default=27)
ap.add_argument("--case", default="meet"); ap.add_argument("--sub", default="00"); ap.add_argument("--cnf", required=True)
a = ap.parse_args()
ag = AG(3, a.q); pool = IDPool()
x = [pool.id(("x", i)) for i in range(len(ag.points))]
y = [pool.id(("y", i)) for i in range(len(ag.lines))]
cls = []; thru = [[] for _ in ag.points]
for di in range(len(ag.dirs)):
    lits = [y[l] for l in ag.lines_by_dir[di]]
    cls.append(lits); cls += CardEnc.atmost(lits=lits, bound=1, vpool=pool, encoding=EncType.pairwise).clauses
for l, (di, pts) in enumerate(ag.lines):
    for p in pts: cls.append([-y[l], x[p]]); thru[p].append(l)
for p in range(len(ag.points)): cls.append([-x[p]] + [y[l] for l in thru[p]])
O = (0, 0, 0)
fixed = [ag.line_through(O, ag.unit(0))]
fixed.append(ag.line_through(O, ag.unit(1)) if a.case == "meet" else ag.line_through(ag.unit(2), ag.unit(1)))
fixed.append(ag.line_through((int(a.sub[0]), int(a.sub[1]), 0), ag.unit(2)))
for l in fixed: cls.append([y[l]])
cls += CardEnc.atmost(lits=x, bound=a.ub, vpool=pool, encoding=EncType.totalizer).clauses
with open(a.cnf, "w") as fo:
    fo.write(f"p cnf {pool.top} {len(cls)}\n")
    for c in cls: fo.write(" ".join(map(str, c)) + " 0\n")
print(json.dumps(dict(case=a.case, sub=a.sub, ub=a.ub, nvars=pool.top, nclauses=len(cls), fixed_lines=[sorted(ag.points[p] for p in ag.lines[l][1]) for l in fixed])))
