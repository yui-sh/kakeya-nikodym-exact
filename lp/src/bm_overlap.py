"""Can a MINIMUM planar Kakeya set contain the grid Sq x Sq?  (odd prime q)

Minimum planar Kakeya size (Blokhuis-Mazzocca 2008) is q(q+1)/2 + (q-1)/2, attained by
the dual-conic example: union of the q affine tangents of the parabola y = x^2
(= {(x,y): x^2 - y is a square or 0}, size q(q+1)/2, covers every non-vertical slope)
plus one vertical line x = c ((q-1)/2 new points).

We enumerate EVERY affine image of that example (all M in GL(2,q), all translations v,
all c) and compute max |image ∩ G|, G = S x S, S = squares incl. 0.  If the max equals
|G| = ((q+1)/2)^2, a minimum planar Kakeya set contains the grid and the completion
extra equals BM_min - |G| exactly.  Otherwise this gives an upper bound
extra <= BM_min - max_overlap restricted to this family (CP-SAT may beat it with larger sets).

Caveat: "every minimum set is an affine image of this example" is the BM characterisation
claim as I remember it -- NOT re-verified here.  The enumeration itself is exhaustive over
affine images of the one example only.

Own mod-p arithmetic; numpy only.  Usage: python3 bm_overlap.py 3 5 7 11 13
"""
import itertools, sys, json
import numpy as np


def run(q):
    S = sorted({a * a % q for a in range(q)})
    G = np.array([(a, b) for a in S for b in S])            # grid points
    # base example P_c, as q x q indicator per c
    sqset = set(S)
    base = np.zeros((q, q), dtype=np.int16)
    for x in range(q):
        for y in range(q):
            if (x * x - y) % q in sqset:
                base[x, y] = 1
    Ps = []
    for c in range(q):
        P = base.copy(); P[c, :] = 1; Ps.append(P)
    Ps = np.stack(Ps)                                        # (q, q, q): c, x, y
    bm_min = q * (q + 1) // 2 + (q - 1) // 2
    assert all(P.sum() == bm_min for P in Ps)
    # |(M^-1)(P) - w ∩ G| over all images  <=>  |P ∩ (M G + v)|; enumerate M in GL2, v, c
    vx, vy = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    vx = vx.ravel(); vy = vy.ravel()                         # q^2 translations
    best = -1; best_arg = None; hist = {}
    for a, b, cc, d in itertools.product(range(q), repeat=4):
        if (a * d - b * cc) % q == 0:
            continue
        MG = np.stack([(a * G[:, 0] + b * G[:, 1]) % q, (cc * G[:, 0] + d * G[:, 1]) % q], 1)
        X = (vx[:, None] + MG[None, :, 0]) % q               # (q^2, |G|)
        Y = (vy[:, None] + MG[None, :, 1]) % q
        cnt = Ps[:, X, Y].sum(axis=2)                        # (q, q^2): c, v
        m = int(cnt.max())
        if m > best:
            ci, vi = np.unravel_index(int(cnt.argmax()), cnt.shape)
            best = m; best_arg = dict(M=[[a, b], [cc, d]], v=[int(vx[vi]), int(vy[vi])], c=int(ci))
        hist[m] = hist.get(m, 0) + 1
    # verify the witness independently: build image explicitly, check Kakeya + overlap
    M = best_arg["M"]; v = best_arg["v"]; c = best_arg["c"]
    P = {(x, y) for x in range(q) for y in range(q) if Ps[c, x, y]}
    # point p in P corresponds to grid point g with M g + v = p  ->  image of P under inverse map
    det = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % q; di = pow(det, q - 2, q)
    Mi = [[M[1][1] * di % q, -M[0][1] * di % q], [-M[1][0] * di % q, M[0][0] * di % q]]
    img = {((Mi[0][0] * (x - v[0]) + Mi[0][1] * (y - v[1])) % q,
            (Mi[1][0] * (x - v[0]) + Mi[1][1] * (y - v[1])) % q) for (x, y) in P}
    Gs = {tuple(g) for g in G.tolist()}
    assert len(img & Gs) == best, (len(img & Gs), best)
    assert is_kakeya2(q, img)
    return dict(q=q, grid=len(G), bm_min=bm_min, max_overlap=best, contains_grid=best == len(G),
                extra_ub_family=bm_min - best, extra_lb_trivial=bm_min - len(G), witness=best_arg,
                hist={int(k): v for k, v in sorted(hist.items())[-4:]})


def is_kakeya2(q, P):
    dirs = [(1, m) for m in range(q)] + [(0, 1)]
    return all(any(all(((x + t * dx) % q, (y + t * dy) % q) in P for t in range(q)) for (x, y) in P)
               for dx, dy in dirs)


if __name__ == "__main__":
    out = []
    for q in map(int, sys.argv[1:]):
        r = run(q); print(r, flush=True); out.append(r)
    json.dump(out, open("out/bm_overlap.json", "w"), indent=1)
