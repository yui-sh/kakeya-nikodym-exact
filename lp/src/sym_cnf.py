"""G-invariant Kakeya sets in AG(3,q) of size <= ub as plain CNF for kissat (session IV).
WLOG a minimal G-invariant Kakeya set is a union of G-orbits of lines (close one witness line per direction under G).
Variables: Y_k per line orbit, X_i per point orbit. Clauses: every direction has a chosen line orbit meeting it;
Y_k -> X_i for point orbits on the orbit's lines; X_i -> OR of covering Y_k; sum |orbit_i| X_i <= ub as a totalizer
over the multiset (X_i repeated |orbit_i| times).
Usage: sym_cnf.py --q 5 --ub 52 --gens out/sym_q5_followup.jsonl --idx 0 --cnf out/symc_0.cnf"""
import argparse, json
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from geom import AG
from sym import perm_of, orbits_of

ap = argparse.ArgumentParser()
ap.add_argument("--q", type=int, default=5); ap.add_argument("--ub", type=int, default=52)
ap.add_argument("--gens", required=True); ap.add_argument("--idx", type=int, required=True); ap.add_argument("--cnf", required=True)
a = ap.parse_args()
recs = [json.loads(l) for l in open(a.gens) if l.strip()]
rec = recs[a.idx]
ag = AG(3, a.q)
gens = [perm_of(ag, g["M"], g["b"]) for g in rec["gens"]]
plab, porbs = orbits_of(len(ag.points), gens)
line_perms = [[ag.line_of[frozenset(g[p] for p in pts)] for (di, pts) in ag.lines] for g in gens]
llab, lorbs = orbits_of(len(ag.lines), line_perms)
pool = IDPool()
X = [pool.id(("X", i)) for i in range(len(porbs))]
Y = [pool.id(("Y", k)) for k in range(len(lorbs))]
cls = []
for di in range(len(ag.dirs)):
    cls.append(sorted({Y[llab[l]] for l in ag.lines_by_dir[di]}))
cover = [set() for _ in porbs]
for k, lo in enumerate(lorbs):
    pos = {plab[p] for p in ag.lines[lo[0]][1]}
    for po in pos:
        cls.append([-Y[k], X[po]]); cover[po].add(k)
for i in range(len(porbs)):
    cls.append([-X[i]] + [Y[k] for k in sorted(cover[i])])
lits = []
for i, po in enumerate(porbs): lits += [X[i]] * len(po)
assert len(lits) == len(ag.points)
cls += CardEnc.atmost(lits=lits, bound=a.ub, vpool=pool, encoding=EncType.totalizer).clauses
with open(a.cnf, "w") as fo:
    fo.write(f"p cnf {pool.top} {len(cls)}\n")
    for c in cls: fo.write(" ".join(map(str, c)) + " 0\n")
json.dump(dict(name=rec["name"], order=rec["order"], gens=rec["gens"], npt_orbits=len(porbs), nline_orbits=len(lorbs),
               X=X, porbs=porbs, nvars=pool.top, nclauses=len(cls)), open(a.cnf.replace(".cnf", "_map.json"), "w"))
print(json.dumps(dict(name=rec["name"], order=rec["order"], npt_orbits=len(porbs), nline_orbits=len(lorbs), nvars=pool.top, nclauses=len(cls))))
