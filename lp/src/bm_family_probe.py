"""Probe the two witness shapes from bm_overlap_fast at larger q (restricted families, NOT exhaustive).
A: M = [[1,k],[0,k]]  (x,y) -> (x+ky, ky) + v   (q = 3 mod 4 witnesses, k nonsquare)
B: M = [[1,k],[0,l]]  for k,l in a small set, all v (generic upper-triangular slice)
Overlap = |image(P_c) ∩ G| maximised over v and line c, same count as bm_overlap_fast.run.
Usage: python3 bm_family_probe.py 47 53 59 ..."""
import sys, json
import numpy as np


def overlap_for(q, r1, r2, base, nb, G, vx, vy):
    MG0 = (r1[0] * G[:, 0] + r1[1] * G[:, 1]) % q
    MG1 = (r2[0] * G[:, 0] + r2[1] * G[:, 1]) % q
    X = (vx[:, None] + MG0[None]) % q; Y = (vy[:, None] + MG1[None]) % q
    b = base[X, Y].sum(1)
    nbv = nb[X, Y]
    extra_col = np.zeros((q * q, q), np.int32)
    np.add.at(extra_col, (np.repeat(np.arange(q * q), X.shape[1]), X.ravel()), nbv.ravel())
    return b, extra_col.max(1)


def run(q, ks):
    S = sorted({a * a % q for a in range(q)}); sq = np.zeros(q, bool); sq[S] = True
    G = np.array([(a, b) for a in S for b in S])
    x, y = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    base = sq[(x * x - y) % q].astype(np.int32); nb = 1 - base
    vx, vy = [a.ravel() for a in np.meshgrid(np.arange(q), np.arange(q), indexing="ij")]
    res = {}
    for k, l in ks:
        # map is x' = x + k y, y' = l y ; in bm_overlap_fast the reps act on G (inverse direction) — same
        # parametrisation as the stored witnesses M=[[1,k],[0,l]].
        b, e = overlap_for(q, (1, k), (0, l), base, nb, G, vx, vy)
        cnt = b + e; i = int(cnt.argmax())
        res[(k, l)] = (int(cnt[i]), int(b[i]), int(e[i]))
    return res


if __name__ == "__main__":
    out = []
    for q in map(int, sys.argv[1:]):
        sqs = {a * a % q for a in range(1, q)}
        ns = min(a for a in range(2, q) if a not in sqs)
        ks = [(1, 1), (ns, ns), (2, 1), (1, ns), (ns, 1)]
        r = run(q, ks)
        best = max(v[0] for v in r.values())
        g = ((q + 1) // 2) ** 2
        rec = dict(q=q, mod8=q % 8, nonsq=ns, half_G_plus_q=g / 2 + q,
                   fam={f"{k},{l}": v for (k, l), v in r.items()}, best=best)
        print(rec, flush=True); out.append(rec)
    json.dump(out, open("out/bm_family_probe.json", "w"), indent=1)
