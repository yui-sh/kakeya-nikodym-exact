"""Long CP-SAT run on the planar completion for one q, hinted with the best BM-family image.
Usage: .venv/bin/python bc_long.py 11 --time 2400 --workers 2"""
import sys, json, time
from ortools.sat.python import cp_model
from bc_complete import dvir_almost, lines2
args = sys.argv[1:]; q = int(args[0])
tl = float(args[args.index("--time") + 1]) if "--time" in args else 600
nw = int(args[args.index("--workers") + 1]) if "--workers" in args else 2
grid = {(a1, a2) for (a1, a2, t) in dvir_almost(q) if t == 0}
bm = {r["q"]: r for r in json.load(open("out/bm_overlap.json"))}[q]["witness"]
# rebuild hint image (inverse affine map of P_c)
M, v, c = bm["M"], bm["v"], bm["c"]
S = {a * a % q for a in range(q)}
P = {(x, y) for x in range(q) for y in range(q) if (x * x - y) % q in S or x == c}
det = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % q; di = pow(det, q - 2, q)
Mi = [[M[1][1] * di % q, -M[0][1] * di % q], [-M[1][0] * di % q, M[0][0] * di % q]]
hint = {((Mi[0][0] * (x - v[0]) + Mi[0][1] * (y - v[1])) % q, (Mi[1][0] * (x - v[0]) + Mi[1][1] * (y - v[1])) % q) for x, y in P}
m = cp_model.CpModel(); pts = [(x, y) for x in range(q) for y in range(q)]
X = {p: m.NewBoolVar("") for p in pts}
for p in grid: m.Add(X[p] == 1)
for p in pts: m.AddHint(X[p], int(p in hint))
for d, Ls in lines2(q).items():
    Y = [m.NewBoolVar("") for _ in Ls]; m.AddBoolOr(Y)
    for y, L in zip(Y, Ls):
        for p in L: m.AddImplication(y, X[p])
m.Add(sum(X[p] for p in pts) >= q * (q + 1) // 2 + (q - 1) // 2)  # Blokhuis-Mazzocca planar minimum (valid cut)
m.Minimize(sum(X[p] for p in pts if p not in grid))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = nw
t0 = time.time(); st = s.Solve(m)
print(dict(q=q, status=s.StatusName(st), extra=int(s.ObjectiveValue()), lb=int(s.BestObjectiveBound()),
           hint_extra=len(hint - grid), secs=round(time.time() - t0)), flush=True)
