"""Rebuild the 3D set (Dvir almost-Kakeya K' ∪ best BM image in slice t=0) from bm_fast.log witnesses
and verify with bc_complete.is_kakeya3 (independent of the counting code).  Usage: bm_verify3d.py 17 19 23"""
import sys, ast
from bc_complete import dvir_almost, is_kakeya3
W = {r["q"]: r for r in map(ast.literal_eval, open("out/bm_fast.log"))}
for q in map(int, sys.argv[1:]):
    r = W[q]; (a, b), (c_, d) = r["M"] if "M" in r else r["witness"]["M"]; v = r["witness"]["v"]
    S = {t * t % q for t in range(q)}
    Kp = dvir_almost(q); grid = {(x, y) for (x, y, t) in Kp if t == 0}
    det = (a * d - b * c_) % q; di = pow(det, q - 2, q)
    inv = lambda x, y: ((d * (x - v[0]) - b * (y - v[1])) * di % q, (-c_ * (x - v[0]) + a * (y - v[1])) * di % q)
    best = None
    for c in range(q):  # vertical line not stored in the log: try all, keep the smallest
        P = {inv(x, y) for x in range(q) for y in range(q) if (x * x - y) % q in S or x == c}
        if best is None or len(P - grid) < len(best - grid): best = P
    K = Kp | {(x, y, 0) for (x, y) in best}
    print(dict(q=q, size=len(K), claimed=len(Kp) + r["extra_ub_family"], kakeya=is_kakeya3(q, K)), flush=True)
