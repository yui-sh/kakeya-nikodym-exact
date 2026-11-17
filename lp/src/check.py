"""Independent checker: verify a point set in AG(n,q) is a Kakeya set.
Usage: python check.py out/sol_n3_q3.json
Re-derives lines directly from field arithmetic (does not use geom.AG.lines)."""
import json, sys, itertools
from gf import GF

def check(n, q, S):
    F = GF(q)
    S = set(tuple(p) for p in S)
    dirs = [v for v in itertools.product(range(q), repeat=n) if any(v) and v[next(i for i, c in enumerate(v) if c)] == 1]
    assert len(dirs) == (q**n - 1) // (q - 1)
    missing = []
    for d in dirs:
        found = False
        for p in S:
            if all(tuple(F.add[a][F.mul[t][b]] for a, b in zip(p, d)) in S for t in range(q)):
                found = True; break
        if not found:
            missing.append(d)
    return missing

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        d = json.load(open(fn))
        missing = check(d["n"], d["q"], d["points"])
        print(f"{fn}: n={d['n']} q={d['q']} |S|={len(d['points'])} -> {'KAKEYA OK' if not missing else 'MISSING ' + str(missing)}")
        if missing: sys.exit(1)
