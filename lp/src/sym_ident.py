"""G-invariant Kakeya sets of exact size n in AG(3,q): orbit-collapsed CP-SAT model (as sym.py) + the session-IV counting
identities and planar/projection cuts (k_pi, r_pi are linear in the orbit variables). Exact within the class G.
Usage: sym_ident.py --q 5 --gens out/sym_q5_followup.jsonl --idx 4 --n 46 47 ... [--time 900]"""
import argparse, json, math, time
from ortools.sat.python import cp_model
from geom import AG, is_kakeya
from sym import perm_of, orbits_of
from structure import planes_of
from plane_f import get

def build(ag, gens, n, ident=True, planes=True, proj=True):
    q = ag.q; nd = len(ag.dirs); f = get(q)
    plab, porbs = orbits_of(len(ag.points), gens)
    line_perms = [[ag.line_of[frozenset(g[p] for p in pts)] for (di, pts) in ag.lines] for g in gens]
    llab, lorbs = orbits_of(len(ag.lines), line_perms)
    m = cp_model.CpModel()
    X = [m.NewBoolVar(f"X{i}") for i in range(len(porbs))]
    Y = [m.NewBoolVar(f"Y{k}") for k in range(len(lorbs))]
    x = [X[plab[p]] for p in range(len(ag.points))]      # point indicator (shared per orbit)
    y = [Y[llab[l]] for l in range(len(ag.lines))]
    for di in range(nd):
        m.AddBoolOr(list({Y[llab[l]] for l in ag.lines_by_dir[di]}))
    cover = [set() for _ in porbs]
    for k, lo in enumerate(lorbs):
        for po in {plab[p] for p in ag.lines[lo[0]][1]}:
            m.AddImplication(Y[k], X[po]); cover[po].add(k)
    for i in range(len(porbs)):
        m.AddBoolOr([Y[k] for k in cover[i]]).OnlyEnforceIf(X[i])
    m.Add(sum(len(po) * X[i] for i, po in enumerate(porbs)) == n)
    # NOTE: a G-invariant set is a union of line ORBITS; a direction may get several lines. The identities below assume
    # one witness line per direction, so we add explicit witness choice: w_l per line, exactly one per direction, w_l -> y_l.
    w = [m.NewBoolVar(f"w{l}") for l in range(len(ag.lines))]
    for di in range(nd): m.AddExactlyOne([w[l] for l in ag.lines_by_dir[di]])
    for l in range(len(ag.lines)): m.AddImplication(w[l], y[l])
    # the set K' = union of witness lines has size n' <= n; identities hold for K' with its own n'. Use K' variables:
    xw = [m.NewBoolVar(f"xw{p}") for p in range(len(ag.points))]
    thru = [[] for _ in ag.points]
    for l, (di, pts) in enumerate(ag.lines):
        for p in pts: m.AddImplication(w[l], xw[p]); thru[p].append(l)
    for p in range(len(ag.points)):
        m.AddBoolOr([w[l] for l in thru[p]]).OnlyEnforceIf(xw[p]); m.AddImplication(xw[p], x[p])
    nn = m.NewIntVar(1, n, "nn"); m.Add(nn == sum(xw))
    PL = planes_of(ag)
    line_in_plane = [[l for l, (di, lp) in enumerate(ag.lines) if lp <= pts] for nrm, c, pts in PL]
    ks, rs = [], []
    sq = [i * i for i in range(q * q + 1)]; c3 = [math.comb(i, 3) for i in range(q * q + 1)]; c2 = [math.comb(i, 2) for i in range(q * q + 2)]
    for (nrm, c, pts), lines in zip(PL, line_in_plane):
        k = m.NewIntVar(1, q * q, ""); m.Add(k == sum(xw[p] for p in pts))
        r = m.NewIntVar(0, q + 1, ""); m.Add(r == sum(w[l] for l in lines))
        ks.append(k); rs.append(r)
        if planes:
            t = m.NewIntVar(0, q * q, ""); m.AddElement(r, [f[i] for i in range(q + 2)], t); m.Add(k >= t)
        m.Add(k <= (nd - q - 1) + q * r)
    if proj:
        bm = f[q + 1]
        for di in range(nd):
            zs = []
            for l in ag.lines_by_dir[di]:
                z = m.NewBoolVar(""); pts = ag.lines[l][1]
                m.AddBoolOr([xw[p] for p in pts]).OnlyEnforceIf(z)
                for p in pts: m.AddImplication(xw[p], z)
                zs.append(z)
            m.Add(sum(zs) >= bm)
    if ident:
        nsq = m.NewIntVar(0, n * n, ""); m.AddMultiplicationEquality(nsq, [nn, nn])
        nc3 = m.NewIntVar(0, math.comb(n, 3), ""); m.AddElement(nn, [math.comb(i, 3) for i in range(n + 1)], nc3)
        nc2 = m.NewIntVar(0, math.comb(n, 2), ""); m.AddElement(nn, [math.comb(i, 2) for i in range(n + 1)], nc2)
        ksq, kr, kc3, rc2 = [], [], [], []
        for k, r in zip(ks, rs):
            a = m.NewIntVar(0, q ** 4, ""); m.AddElement(k, sq, a); ksq.append(a)
            b = m.NewIntVar(0, q * q * (q + 1), ""); m.AddMultiplicationEquality(b, [k, r]); kr.append(b)
            e = m.NewIntVar(0, c3[-1], ""); m.AddElement(k, c3, e); kc3.append(e)
            g = m.NewIntVar(0, c2[q + 1], ""); m.AddElement(r, c2[:q + 2], g); rc2.append(g)
        m.Add(sum(ksq) == (q + 1) * nsq + (nd - q - 1) * nn)
        m.Add(sum(kr) == nd * nn + q * q * nd)
        dc2 = []
        for p in range(len(ag.points)):
            dp = m.NewIntVar(0, nd, ""); m.Add(dp == sum(w[l] for l in thru[p]))
            g = m.NewIntVar(0, math.comb(nd, 2), ""); m.AddElement(dp, [math.comb(i, 2) for i in range(nd + 1)], g); dc2.append(g)
        m.Add(sum(rc2) == sum(dc2))
        lc3, lc2 = [], []
        for l, (di, pts) in enumerate(ag.lines):
            lam = m.NewIntVar(0, q, ""); m.Add(lam == sum(xw[p] for p in pts))
            e = m.NewIntVar(0, math.comb(q, 3), ""); m.AddElement(lam, [math.comb(i, 3) for i in range(q + 1)], e); lc3.append(e)
            g = m.NewIntVar(0, math.comb(q, 2), ""); m.AddElement(lam, [math.comb(i, 2) for i in range(q + 1)], g); lc2.append(g)
        m.Add(sum(kc3) == nc3 + q * sum(lc3))
        m.Add(sum(lc2) == nc2)
    return m, X, porbs

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5); ap.add_argument("--gens", required=True); ap.add_argument("--idx", type=int, required=True)
    ap.add_argument("--n", type=int, nargs="+", required=True); ap.add_argument("--time", type=float, default=900); ap.add_argument("--workers", type=int, default=1); ap.add_argument("--hint", default=None)
    a = ap.parse_args(); ag = AG(3, a.q)
    rec = [json.loads(l) for l in open(a.gens) if l.strip()][a.idx]
    gens = [perm_of(ag, g["M"], g["b"]) for g in rec["gens"]]
    for n in a.n:
        m, X, porbs = build(ag, gens, n)
        if a.hint:
            S = {ag.pidx[tuple(p)] for p in json.load(open(a.hint))["points"]}
            inv = all(g[p] in S for g in gens for p in S)
            print("hint set invariant under the class generators:", inv)
            for i, po in enumerate(porbs): m.AddHint(X[i], int(po[0] in S))
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = a.time; s.parameters.num_workers = a.workers
        t0 = time.time(); st = s.Solve(m); dt = time.time() - t0
        res = dict(name=rec["name"], order=rec["order"], n=n, status=s.StatusName(st), time=round(dt, 1))
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            S = [p for i, po in enumerate(porbs) if s.Value(X[i]) for p in po]; ok, miss = is_kakeya(ag, S); assert ok
            res["size"] = len(S); res["points"] = [ag.points[p] for p in S]
            json.dump(dict(n=3, q=a.q, points=res["points"], group=rec["name"], status="G-invariant, CP-SAT"), open(f"out/sym_ident_{a.idx}_n{n}.json", "w"))
        print(json.dumps({k: v for k, v in res.items() if k != "points"}), flush=True)
