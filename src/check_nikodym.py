#!/usr/bin/env python3
"""Independent Nikodym witness checker.  No solver, no imports beyond ff.py (own field arithmetic).

Definition used throughout this release (the 'strong' / every-point variant of
Lund-Saraf-Wolf, Tao, Georgiev et al.):  N in AG(n,q) is a NIKODYM SET if for EVERY point x of
AG(n,q) -- including the points of N itself -- there is a line L through x with L \\ {x} contained in N.
usage: check_nikodym.py file.json [...]   (JSON: {"n","q", "N":[[..]..]} or {"n","q","C":[[..]..]} with C the complement)
Exit status 0 iff every file is a Nikodym set.  Prints |N| so the size claim can be read off.
"""
import itertools, json, sys
from ff import GF


def is_nikodym(n, q, N):
    F = GF(q)
    N = set(map(tuple, N))
    assert all(len(p) == n and all(0 <= c < q for c in p) for p in N)
    pts = list(itertools.product(range(q), repeat=n))
    dirs = [v for v in pts if any(v) and v[next(i for i, c in enumerate(v) if c)] == 1]
    bad = []
    for x in pts:
        if not any(all(tuple(F.add[a][F.mul[t][b]] for a, b in zip(x, d)) in N for t in range(1, q)) for d in dirs):
            bad.append(x)
    return bad, len(N)


if __name__ == "__main__":
    rc = 0
    for fn in sys.argv[1:]:
        d = json.load(open(fn))
        n, q = d["n"], d["q"]
        if "N" in d:
            N = [tuple(p) for p in d["N"]]
        else:
            C = set(map(tuple, d["C"]))
            N = [p for p in itertools.product(range(q), repeat=n) if p not in C]
        bad, size = is_nikodym(n, q, N)
        ok = not bad
        print(f"{fn}: AG({n},{q}) |N|={size} |C|={q**n-size} -> {'NIKODYM OK' if ok else 'NOT NIKODYM, uncovered points ' + str(bad[:4])}")
        rc |= (not ok)
    sys.exit(rc)
