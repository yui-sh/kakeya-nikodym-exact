#!/usr/bin/env python3
"""Independent Kakeya witness checker.  No solver, no imports beyond ff.py (own field arithmetic).

A set S in AG(n,q) is a KAKEYA SET if for every direction d (nonzero vector up to scalars)
there is a point p with the whole line {p + t d : t in F_q} contained in S.
usage: check_kakeya.py file.json [...]   (JSON: {"n":..,"q":..,"points":[[..],..]})
Exit status 0 iff every file is a Kakeya set.  Prints |S| so the size claim can be read off.
"""
import itertools, json, sys
from ff import GF


def is_kakeya(n, q, S):
    F = GF(q)
    S = set(map(tuple, S))
    assert all(len(p) == n and all(0 <= c < q for c in p) for p in S)
    dirs = [v for v in itertools.product(range(q), repeat=n) if any(v) and v[next(i for i, c in enumerate(v) if c)] == 1]
    assert len(dirs) == (q ** n - 1) // (q - 1)
    missing = []
    for d in dirs:
        if not any(all(tuple(F.add[a][F.mul[t][b]] for a, b in zip(p, d)) in S for t in range(q)) for p in S):
            missing.append(d)
    return missing, len(S)


if __name__ == "__main__":
    rc = 0
    for fn in sys.argv[1:]:
        d = json.load(open(fn))
        missing, size = is_kakeya(d["n"], d["q"], d["points"])
        ok = not missing
        print(f"{fn}: AG({d['n']},{d['q']}) |S|={size} -> {'KAKEYA OK' if ok else 'NOT KAKEYA, missing directions ' + str(missing[:4])}")
        rc |= (not ok)
    sys.exit(rc)
