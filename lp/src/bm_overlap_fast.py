"""Same as bm_overlap.py, but M ranges over coset reps of H\\GL(2,q), H = linear parts
[[a,0],[2ab,a^2]] of the stabiliser of the parabola family {P_c} (q^2-1 reps instead of |GL2|).
Cross-checked against bm_overlap.py on q<=13.  Usage: python3 bm_overlap_fast.py 3 5 ... 31"""
import sys, json
import numpy as np

def run(q):
    S = sorted({a * a % q for a in range(q)}); sq = np.zeros(q, bool); sq[S] = True
    G = np.array([(a, b) for a in S for b in S])
    x, y = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    base = sq[(x * x - y) % q].astype(np.int32)
    # count = |P_c ∩ T| = |base ∩ T| + |{t in T: x=c, not base}|  -> no need to materialise q copies
    nb = 1 - base
    reps = [((1, t), (0, l)) for t in range(q) for l in range(1, q)] + [((0, 1), (l, 0)) for l in range(1, q)]
    vx, vy = [a.ravel() for a in np.meshgrid(np.arange(q), np.arange(q), indexing="ij")]
    best = -1; arg = None; hist = {}
    for r1, r2 in reps:
        MG0 = (r1[0] * G[:, 0] + r1[1] * G[:, 1]) % q
        MG1 = (r2[0] * G[:, 0] + r2[1] * G[:, 1]) % q
        X = (vx[:, None] + MG0[None]) % q; Y = (vy[:, None] + MG1[None]) % q   # (q^2,|G|)
        b = base[X, Y].sum(1)                                              # (q^2,)
        nbv = nb[X, Y]                                                     # (q^2,|G|)
        # for each translation, best vertical line = column x with most non-base points
        extra_col = np.zeros((q * q, q), np.int32)
        np.add.at(extra_col, (np.repeat(np.arange(q * q), X.shape[1]), X.ravel()), nbv.ravel())
        cnt = b + extra_col.max(1)
        m = int(cnt.max())
        for k, n in zip(*np.unique(cnt, return_counts=True)): hist[int(k)] = hist.get(int(k), 0) + int(n)
        if m > best: best = m; arg = dict(M=[list(r1), list(r2)], v=[int(vx[cnt.argmax()]), int(vy[cnt.argmax()])])
    bm_min = q * (q + 1) // 2 + (q - 1) // 2
    return dict(q=q, qmod4=q % 4, grid=len(G), bm_min=bm_min, max_overlap=best, extra_ub_family=bm_min - best,
                extra_lb_trivial=bm_min - len(G), witness=arg, top_hist=dict(sorted(hist.items())[-3:]))

if __name__ == "__main__":
    out = []
    for q in map(int, sys.argv[1:]):
        r = run(q); print(r, flush=True); out.append(r)
    json.dump(out, open("out/bm_overlap_fast.json", "w"), indent=1)
