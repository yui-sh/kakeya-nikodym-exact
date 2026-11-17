"""Validate the plane-counting identities of planecount.py on actual Kakeya sets (witness systems) and
run the aggregated model on q=2,3,4 at the known optimum (must be feasible) and below."""
import json, math, sys, itertools
from geom import AG
from structure import planes_of
from plane_f import get
from planecount import feasible

def stats(ag, S, W):
    """S: point idx set, W: witness line idx list (one per direction). Returns all aggregated statistics."""
    q = ag.q; nd = len(ag.dirs); n = len(S)
    deg = {p: 0 for p in S}
    for l in W:
        for p in ag.lines[l][1]: deg[p] += 1
    Wset = set(W)
    N = {}
    for nrm, c, pts in planes_of(ag):
        r = sum(1 for l in W if ag.lines[l][1] <= pts); k = len(pts & S)
        N[(r, k)] = N.get((r, k), 0) + 1
        assert sum(deg[p] for p in pts & S) == (nd - q - 1) + q * r, "I1"
    L = {}
    for di, pts in ag.lines:
        j = len(pts & S)
        if j: L[j] = L.get(j, 0) + 1
    D = max(deg.values())
    assert sum(N.values()) == q * nd and sum(r * v for (r, k), v in N.items()) == (q + 1) * nd and sum(k * v for (r, k), v in N.items()) == nd * n, "I2"
    assert sum(math.comb(r, 2) * v for (r, k), v in N.items()) == sum(math.comb(d, 2) for d in deg.values()), "I3"
    assert sum(k * k * v for (r, k), v in N.items()) == (q + 1) * n * n + (nd - q - 1) * n, "I4"
    assert sum(k * r * v for (r, k), v in N.items()) == nd * n + q * nd * q, "I5"
    f = get(q); bm = f[q + 1]
    assert sum(L.values()) >= nd * bm and sum(j * v for j, v in L.items()) == nd * n and sum(math.comb(j, 2) * v for j, v in L.items()) == math.comb(n, 2) and L.get(q, 0) >= nd, "I6"
    assert sum(math.comb(k, 3) * v for (r, k), v in N.items()) == math.comb(n, 3) + q * sum(math.comb(j, 3) * v for j, v in L.items()), "I6b"
    for (r, k) in N: assert f[r] <= k <= (nd - q - 1) + q * r <= D * k, ("plane", r, k)
    return dict(n=n, D=D, deg=sorted(deg.values()), N=sorted(N.items()), L=sorted(L.items()))

def witness(ag, S):
    S = set(S); W = []
    for di in range(len(ag.dirs)):
        ls = [l for l in ag.lines_by_dir[di] if ag.lines[l][1] <= S]
        assert ls; W.append(ls[0])
    return W

if __name__ == "__main__":
    for fn in ["solutions/best_n3_q3_13.json", "solutions/best_n3_q4_28_dvir.json", "solutions/best_n3_q5_53.json",
               "solutions/bc_n3_q5_53.json", "solutions/sym_n3_q5_53_ord4_1.json", "solutions/sym_n3_q5_53_ord2_3.json", "solutions/bc_n3_q7_128.json"]:
        d = json.load(open(fn)); ag = AG(3, d["q"])
        S = [ag.pidx[tuple(p)] for p in d["points"]]
        # Dvir q=4 has 36 lines: try several witness choices
        st = stats(ag, set(S), witness(ag, S))
        print(fn, "identities OK; n", st["n"], "D", st["D"], "planes (r,k):", st["N"], "lines:", st["L"])
    for q, nopt in [(2, 5), (3, 13), (4, 28)]:
        nd = q * q + q + 1
        for n in range(max(nd * q // 20, 3), nopt + 1):
            Dmin = math.ceil(nd * q / n); Dmax = (n - 1) // (q - 1)
            row = []
            for D in range(Dmin, Dmax + 1):
                stt, sol = feasible(q, n, D)
                row.append(f"D={D}:{'X' if stt == 'INFEASIBLE' else 'ok'}")
            print(f"q={q} n={n}: " + " ".join(row))
