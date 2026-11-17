"""Session V: branch-and-bound over the degree vector (m_1, ..., m_D) for one LP-surviving (r_max configuration, D)
case. Node = prefix (m_1..m_j) fixed by equality rows; the merged LP (planecount6.build + configuration rows) is solved
at every node; infeasible nodes are pruned WITH an exact integer Farkas certificate (planecount5.farkas_phase1);
feasible leaves (all m fixed) are the survivors. Output JSON lists pruned nodes (prefix + certificate) and survivors;
cert5_check.py --bb re-verifies every certificate and that pruned nodes + survivors tile the whole tree.
Usage: degvec_bb.py --n 46 --D 5 --config 5,16 [--no-cert]"""
import argparse, json, time, math, sys, copy
import numpy as np
import planecount5 as p5, planecount6 as p6

def children(prefix, n, tot, D):
    """possible values of m_{j+1} given prefix m_1..m_j (sum m = n, sum i m_i = tot, m_D >= 1)."""
    j = len(prefix); left_n = n - sum(prefix); left_t = tot - sum((i + 1) * v for i, v in enumerate(prefix))
    i = j + 1
    if i == D: return [left_n] if (left_t == D * left_n and left_n >= 1) else []
    out = []
    for h in range(left_n + 1):
        rn, rt = left_n - h, left_t - i * h
        if rt < (i + 1) * rn or rt > D * rn: continue
        out.append(h)
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, default=46)
    ap.add_argument("--D", type=int, required=True); ap.add_argument("--config", required=True); ap.add_argument("--no-cert", action="store_true")
    a = ap.parse_args(); q, n, D = a.q, a.n, a.D; r, u = map(int, a.config.split(",")); nd = q * q + q + 1; tot = nd * q
    lp0, info, types = p6.build(q, n, D, u_for_r=r)
    idx = [i for i, t in enumerate(types) if t[0] == r and t[3] == u]; lp0.add([(lp0.col[f"N{i}"], -1) for i in idx], '<', -1, f"config ({r},{u}) present")
    idx = [i for i, t in enumerate(types) if t[0] > r]
    if idx: lp0.add([(lp0.col[f"N{i}"], 1) for i in idx], '<', 0, f"rmax<={r}")
    mcol = {j: lp0.col[f"m{j}"] for j in range(1, D + 1)}
    pruned, survivors, nodes = [], [], 0; t0 = time.time()
    def solve_node(prefix):
        lp = copy.copy(lp0); lp.rows = list(lp0.rows)
        for j, v in enumerate(prefix, start=1): lp.add([(mcol[j], 1)], '=', v, f"fix m{j} = {v}")
        A, b, sense = lp.matrices(); eq = [i for i, s in enumerate(sense) if s == '=']; ub = [i for i, s in enumerate(sense) if s == '<']
        if a.no_cert:
            res, *_ = p5.solve(lp); return ("feasible" if res.status == 0 else "INFEASIBLE"), None
        yr, msg, x = p5.farkas_phase1(lp, A, b, sense, eq, ub)
        if x is not None: return "feasible", None
        if yr is None: return "UNDECIDED: " + msg, None
        return "INFEASIBLE", {lp.rows[i][3]: str(yr[i]) for i in range(len(yr)) if yr[i] != 0}
    def rec(prefix):
        nonlocal nodes
        nodes += 1
        st, cert = solve_node(prefix)
        if st == "INFEASIBLE": pruned.append(dict(prefix=prefix, cert=cert)); return
        if st != "feasible": pruned.append(dict(prefix=prefix, status=st)); return
        if len(prefix) == D: survivors.append(prefix); return
        for h in children(prefix, n, tot, D): rec(prefix + [h])
    rec([])
    und = [p for p in pruned if "cert" not in p]
    out = dict(config=a.config, n=n, D=D, nodes=nodes, pruned=len(pruned), undecided=len(und), survivors=survivors, pruned_nodes=pruned, time=round(time.time() - t0, 1))
    json.dump(out, open(f"out/degvec_bb_{a.config.replace(',', '_')}_D{D}.json", "w"))
    print(json.dumps({k: v for k, v in out.items() if k != "pruned_nodes"}), flush=True)

if __name__ == "__main__": main()
