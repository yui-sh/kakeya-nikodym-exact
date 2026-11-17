"""Dvir almost-Kakeya set + optimally completed horizontal directions (n=3, odd PRIME q).

Bukh-Chao (Discrete Analysis 2021:26, Prop. 11) note that Dvir's odd-q construction
K' = {(a1,a2,t): a_i + t^2 is a square (incl. 0)}, |K'| = q((q+1)/2)^2,
handles every non-horizontal direction, and that the horizontal directions can be handled by a
2-D Kakeya set in ONE horizontal plane instead of the whole plane (Dvir/Saraf-Sudan add the
plane: |K| = q((q+1)/2)^2 + q^2 - ((q+1)/2)^2).  They estimate a RANDOM translate.

Here we OPTIMISE it exactly: for a horizontal plane t=c, the slice of K' is the grid
G_c = {(a1,a2): a_i + c^2 square}; we need a planar Kakeya set P in AG(2,q) minimising |P \\ G_c|.
All slices are translates of Sq x Sq, so one CP-SAT model per q.  Result size = |K'| + min|P\\G|.

Independent of gf.py / check.py on purpose (own mod-p arithmetic, own checker).
Usage: .venv/bin/python bc_complete.py 3 5 7 11 [--time 120]
"""
import itertools, json, sys
from ortools.sat.python import cp_model


def dvir_almost(q):
    sq = {a * a % q for a in range(q)}
    return {(a1, a2, t) for a1 in range(q) for a2 in range(q) for t in range(q)
            if (a1 + t * t) % q in sq and (a2 + t * t) % q in sq}


def lines2(q):
    """All lines of AG(2,q) grouped by direction; each line a frozenset of points."""
    dirs = [(1, m) for m in range(q)] + [(0, 1)]
    out = {}
    for d in dirs:
        seen = set()
        for x in range(q):
            for y in range(q):
                L = frozenset(((x + t * d[0]) % q, (y + t * d[1]) % q) for t in range(q))
                seen.add(L)
        out[d] = sorted(seen, key=sorted)
    return out


def complete(q, grid, tlimit):
    m = cp_model.CpModel()
    pts = [(x, y) for x in range(q) for y in range(q)]
    X = {p: m.NewBoolVar(f"x{p}") for p in pts}
    for p in grid:
        m.Add(X[p] == 1)  # free: already in K'
    for d, Ls in lines2(q).items():
        Y = [m.NewBoolVar("") for _ in Ls]
        m.AddBoolOr(Y)
        for y, L in zip(Y, Ls):
            for p in L:
                m.AddImplication(y, X[p])
    m.Minimize(sum(X[p] for p in pts if p not in grid))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tlimit; s.parameters.num_workers = 4
    st = s.Solve(m)
    P = {p for p in pts if s.Value(X[p])}
    return s.StatusName(st), P, int(s.BestObjectiveBound())


def is_kakeya3(q, S):
    S = {tuple(p) for p in S}
    dirs = [d for d in itertools.product(range(q), repeat=3) if any(d) and d[next(i for i in range(3) if d[i])] == 1]
    for d in dirs:
        if not any(all(tuple((p[i] + t * d[i]) % q for i in range(3)) in S for t in range(q)) for p in S):
            return False
    return True


if __name__ == "__main__":
    args = sys.argv[1:]
    tl = 120.0
    if "--time" in args:
        i = args.index("--time"); tl = float(args[i + 1]); del args[i:i + 2]
    res = []
    for q in map(int, args):
        Kp = dvir_almost(q)
        c = 0
        grid = {(a1, a2) for (a1, a2, t) in Kp if t == c}
        st, P, lb = complete(q, grid, tl)
        extra = len(P - grid)
        K = Kp | {(x, y, c) for (x, y) in P}
        ok = is_kakeya3(q, K)
        h = (q + 1) // 2
        row = dict(q=q, almost=len(Kp), grid=len(grid), planar_extra=extra, status=st, extra_lb=lb,
                   size=len(K), kakeya=ok, dvir_plane=q * h * h + q * q - h * h,
                   bc_random=round(len(Kp) + (q * (q + 1) // 2 + (q - 1) // 2) * (1 - len(Kp) / q**3), 1),
                   bc_lb=-(-q**3 * q * q // (2 * q - 1) ** 2))
        print(row, flush=True)
        res.append(row)
        json.dump({"n": 3, "q": q, "points": sorted(K), "note": "Dvir almost-Kakeya + optimal planar completion"},
                  open(f"out/bc_n3_q{q}.json", "w"))
    json.dump(res, open("out/bc_complete.json", "w"), indent=1)
