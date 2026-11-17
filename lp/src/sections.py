"""Session VI: exact per-plane line-section profiles. For AG(2,q), all tuples (k, r, u, c, t) attainable by a point set
S with |S| = k: t_j = # lines of the plane meeting S in exactly j points (j = 1..q), and (u; c_1..c_r) the crossing
profile of some choice of r lines contained in S with pairwise distinct directions. Brute force over all 2^(q^2)
subsets (bin/sections_enum, C; line masks supplied from geom.AG so any prime power works). Cached in out/sections_q{q}.json."""
import json, os, subprocess, sys
from geom import AG

def sections(q, rmin=5):
    fn = f"out/sections3_q{q}_r{rmin}.json"
    def conv(e):
        if not e: return None
        return tuple(e) if not isinstance(e[0], list) else tuple(tuple(x) for x in e)
    if os.path.exists(fn): return [(k, r, u, tuple(c), tuple(t), conv(e)) for k, r, u, c, t, e in json.load(open(fn))]
    ag = AG(2, q)
    inp = [f"{q*q} {q*(q+1)}"] + [f"{di} {sum(1 << p for p in pts)}" for di, pts in ag.lines]
    res = subprocess.run(["bin/sections_enum", str(q), str(rmin)], input="\n".join(inp) + "\n", capture_output=True, text=True, check=True)
    out = []
    for line in res.stdout.splitlines():
        parts = line.split("|"); a = list(map(int, parts[0].split())); t = tuple(map(int, parts[1].split()))
        profs = tuple(tuple(map(int, pp.split())) for pp in parts[2:])
        e = profs if profs else None   # tuple of sorted q-tuples, one per unmarked direction (r >= rmin), else None
        k, r, u = a[:3]; c = tuple(a[3:3 + r]); out.append((k, r, u, c, t, e))
    out.sort(key=lambda x: x[:5]); json.dump([list(x[:3]) + [list(x[3]), list(x[4]), (list(x[5]) if x[5] else None)] for x in out], open(fn, "w"))
    return out

if __name__ == "__main__":
    q = int(sys.argv[1]); rmin = int(sys.argv[2]) if len(sys.argv) > 2 else 5; S = sections(q, rmin); print(len(S), "tuples")
    from inside_profiles import profiles
    P = profiles(q)
    # consistency: every (u,c) crossing profile with r lines appears in inside_profiles and vice versa
    A = {(r, (u,) + c) for (k, r, u, c, t, e) in S}; B = {(r, p) for r, v in P.items() for p in v}
    print("crossing profiles match inside_profiles:", A == B, len(A), len(B))
    for (k, r, u, c, t, e) in S:
        assert (e is None) == (r < rmin or r == q + 1) and (e is None or all(sum(pr) == k for pr in e))
        assert sum(t) <= q * (q + 1) and sum(j * t[j - 1] for j in range(1, q + 1)) == (q + 1) * k
        assert sum(j * (j - 1) // 2 * t[j - 1] for j in range(1, q + 1)) == k * (k - 1) // 2 and t[q - 1] >= r
    print("sanity identities ok")
