"""Session V. For AG(2,q): all achievable profiles of a union of r lines with pairwise distinct directions,
profile = (u, c_1, c_2, ..., c_r) with u = |union| and c_i = # points lying on exactly i of the r lines.
Brute force over all choices (one line per chosen direction); cached in out/inside_profiles_q{q}.json."""
import json, os, sys, itertools
from geom import AG

def profiles(q):
    fn = f"out/inside_profiles_q{q}.json"
    if os.path.exists(fn):
        d = json.load(open(fn)); return {int(r): [tuple(p) for p in v] for r, v in d.items()}
    ag = AG(2, q); nd = q + 1
    out = {0: [(0,) + (0,) * 0]}
    for r in range(1, nd + 1):
        seen = set()
        for dirs in itertools.combinations(range(nd), r):
            for choice in itertools.product(*[ag.lines_by_dir[d] for d in dirs]):
                cnt = {}
                for l in choice:
                    for p in ag.lines[l][1]: cnt[p] = cnt.get(p, 0) + 1
                c = [0] * (r + 1)
                for v in cnt.values(): c[v] += 1
                seen.add((len(cnt),) + tuple(c[1:]))
        out[r] = sorted(seen)
    json.dump({r: [list(p) for p in v] for r, v in out.items()}, open(fn, "w"))
    return out

if __name__ == "__main__":
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    P = profiles(q)
    for r, v in P.items(): print(r, len(v), v[:12])
