"""Are the known best sets affine images of 'Dvir almost-Kakeya + optimal planar completion'?

1. Enumerate ALL optimal completions P (|P \\ grid| = min) for the t=0 slice (all slices equivalent:
   (a,t) -> (a - 2ct - c^2, t + c) maps K' to itself).
2. For each, A = K' u P; test affine equivalence with each target set B (frame mapping with
   profile pruning, like structure.stabilizer but source != target).
Usage: .venv/bin/python bc_iso.py 5 solutions/best_n3_q5_53.json solutions/sym_*.json
"""
import json, sys, itertools
from ortools.sat.python import cp_model
from geom import AG
from structure import affine_map, analyze
from bc_complete import dvir_almost, lines2


def all_completions(q, grid):
    pts = [(x, y) for x in range(q) for y in range(q)]
    # first get optimum
    def model(fix_obj=None):
        m = cp_model.CpModel()
        X = {p: m.NewBoolVar("") for p in pts}
        for p in grid: m.Add(X[p] == 1)
        for d, Ls in lines2(q).items():
            Y = [m.NewBoolVar("") for _ in Ls]
            m.AddBoolOr(Y)
            for y, L in zip(Y, Ls):
                for p in L: m.AddImplication(y, X[p])
        obj = sum(X[p] for p in pts if p not in grid)
        if fix_obj is None: m.Minimize(obj)
        else: m.Add(obj == fix_obj)
        return m, X
    m, X = model(); s = cp_model.CpSolver(); s.Solve(m); opt = int(s.ObjectiveValue())
    m, X = model(opt)
    sols = []
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            sols.append(frozenset(p for p in pts if self.Value(X[p])))
    s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True
    s.Solve(m, CB())
    return opt, sols


def profiles(ag, S):
    S = frozenset(S)
    a = analyze(ag, S); deg = a["deg"]
    sec = [len(pts & S) for di, pts in ag.lines]
    lt = {p: [] for p in range(len(ag.points))}
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts: lt[p].append(l)
    return {p: (deg[p], tuple(sorted(sec[l] for l in lt[p]))) for p in S}, sec


def equivalent(ag, A, B):
    A, B = frozenset(A), frozenset(B)
    if len(A) != len(B): return None
    pa, seca = profiles(ag, A); pb, secb = profiles(ag, B)
    if sorted(pa.values()) != sorted(pb.values()): return None
    P = ag.points
    def ls(sec, x, y): return sec[ag.line_through(P[x], ag.sub(P[y], P[x]))]
    cnt = {}
    for p in A: cnt[pa[p]] = cnt.get(pa[p], 0) + 1
    order = sorted(A, key=lambda p: cnt[pa[p]])
    src = None
    for fr in itertools.permutations(order[:12], 4):
        if affine_map(ag, [P[i] for i in fr], [P[i] for i in fr]) is not None:
            src = fr; break
    cand = {i: [r for r in B if pb[r] == pa[i]] for i in src}
    p0, p1, p2, p3 = src
    for r0 in cand[p0]:
        for r1 in cand[p1]:
            if r1 == r0 or ls(secb, r0, r1) != ls(seca, p0, p1): continue
            for r2 in cand[p2]:
                if r2 in (r0, r1) or ls(secb, r0, r2) != ls(seca, p0, p2) or ls(secb, r1, r2) != ls(seca, p1, p2): continue
                for r3 in cand[p3]:
                    if r3 in (r0, r1, r2): continue
                    if any(ls(secb, x, r3) != ls(seca, y, p3) for x, y in ((r0, p0), (r1, p1), (r2, p2))): continue
                    res = affine_map(ag, [P[i] for i in src], [P[i] for i in (r0, r1, r2, r3)])
                    if res is None: continue
                    perm, M = res
                    if all(perm[p] in B for p in A):
                        return M, P[perm[ag.pidx[(0, 0, 0)]]]
    return None


if __name__ == "__main__":
    q = int(sys.argv[1]); targets = sys.argv[2:]
    ag = AG(3, q)
    Kp = dvir_almost(q)
    grid = {(a, b) for (a, b, t) in Kp if t == 0}
    opt, sols = all_completions(q, grid)
    print(f"q={q}: |K'|={len(Kp)}, optimal extra={opt}, #optimal completions={len(sols)}", flush=True)
    As = [frozenset(ag.pidx[p] for p in Kp | {(x, y, 0) for (x, y) in P}) for P in sols]
    for fn in targets:
        d = json.load(open(fn)); B = frozenset(ag.pidx[tuple(p)] for p in d["points"])
        hit = None
        for k, A in enumerate(As):
            r = equivalent(ag, A, B)
            if r: hit = (k, r); break
        print(f"  {fn}: {'EQUIVALENT to completion #%d, M=%s b=%s' % (hit[0], hit[1][0], hit[1][1]) if hit else 'not equivalent to any'}", flush=True)
