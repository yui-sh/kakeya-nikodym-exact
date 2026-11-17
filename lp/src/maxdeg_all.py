"""Driver: decide 'exists Kakeya set of size <= ub in AG(3,q)' by the max-degree decomposition (maxdeg.py).
D ranges over ceil(#dirs*q/ub) .. floor((ub-1)/(q-1)); every PGammaL(3,q)-orbit of D-subsets of PG(2,q) is a subcase.
Each subcase -> DIMACS CNF -> kissat (DRAT proof) -> drat-trim verification (UNSAT) or check.py-style
verification of the model (SAT). Results appended to out/maxdeg_q{q}_ub{ub}.jsonl.
Usage: maxdeg_all.py --q 4 --ub 27 [--time 7200] [--only D:rep,...] [--no-verify] [--cuts]"""
import argparse, json, subprocess, time, os, sys, math, re
from geom import AG, is_kakeya
from maxdeg import gammaL_gens, dir_perms, subset_orbits, describe, write_cnf
from check import check as indep_check

def run(cmd, timeout=None):
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr, time.time() - t0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=4)
    ap.add_argument("--ub", type=int, default=27)
    ap.add_argument("--time", type=int, default=7200)
    ap.add_argument("--only", default=None)
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--cuts", action="store_true", help="add planar + projection cuts to the CNF")
    ap.add_argument("--solver", default="kissat", choices=["kissat", "cadical"])
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    ag = AG(3, a.q)
    q = a.q; nd = len(ag.dirs)
    Dmin = math.ceil(nd * q / a.ub); Dmax = (a.ub - 1) // (q - 1)
    perms = dir_perms(ag, gammaL_gens(ag))
    only = None
    if a.only:
        only = {tuple(map(int, s.split(":"))) for s in a.only.split(",")}
    out = f"out/maxdeg_q{q}_ub{a.ub}{a.tag}.jsonl"
    print(f"q={q} ub={a.ub}: D in [{Dmin},{Dmax}] (D>=ceil({nd}*{q}/{a.ub}); 1+(q-1)D<=ub); results -> {out}", flush=True)
    summary = []
    for D in range(Dmin, Dmax + 1):
        orbs = subset_orbits(nd, perms, D)
        assert sum(sz for _, sz in orbs) == math.comb(nd, D)
        for i, (rep, sz) in enumerate(orbs):
            if only and (D, i) not in only: continue
            base = f"out/md_q{q}_ub{a.ub}_D{D}_r{i}{a.tag}"
            info = write_cnf(ag, D, rep, a.ub, base + ".cnf", planes=a.cuts, proj=a.cuts)
            if a.solver == "kissat":
                cmd = ["nice", "-n", "10", "bin/kissat", "-q", f"--time={a.time}", base + ".cnf", base + ".drat"]
            else:
                cmd = ["nice", "-n", "10", "bin/cadical", "-q", "-t", str(a.time), base + ".cnf", base + ".drat"]
            rc, so, se, dt = run(cmd)
            status = "UNSAT" if "s UNSATISFIABLE" in so else ("SAT" if "s SATISFIABLE" in so else "UNKNOWN")
            rec = dict(q=q, ub=a.ub, D=D, rep=i, dirs=[ag.dirs[d] for d in rep], orbit_size=sz, max_coplanar=describe(ag, rep),
                       nvars=info["nvars"], nclauses=info["nclauses"], solver=a.solver, cuts=a.cuts,
                       status=status, solve_time=round(dt, 1), cnf=base + ".cnf")
            if status == "SAT":
                lits = set()
                for line in so.splitlines():
                    if line.startswith("v "):
                        lits.update(int(t) for t in line[2:].split())
                S = [p for p in range(len(ag.points)) if info["x"][p] in lits]
                ok, miss = is_kakeya(ag, S)
                miss2 = indep_check(3, q, [ag.points[p] for p in S])
                rec.update(size=len(S), kakeya_ok=bool(ok and not miss2), points=[ag.points[p] for p in S])
                maxdeg = max(sum(1 for l in range(len(ag.lines)) if info["y"][l] in lits and p in ag.lines[l][1]) for p in S)
                rec["model_maxdeg"] = maxdeg
            elif status == "UNSAT" and not a.no_verify:
                rc2, so2, se2, dt2 = run(["nice", "-n", "10", "bin/drat-trim", base + ".cnf", base + ".drat"])
                rec.update(drat_verified=("s VERIFIED" in so2), verify_time=round(dt2, 1),
                           drat_bytes=os.path.getsize(base + ".drat"))
                if rec["drat_verified"]:
                    os.remove(base + ".drat")
            print(json.dumps({k: v for k, v in rec.items() if k != "points"}), flush=True)
            with open(out, "a") as fo: fo.write(json.dumps(rec) + "\n")
            summary.append(rec)
    st = {r["status"] for r in summary}
    print("SUMMARY:", {s: sum(1 for r in summary if r["status"] == s) for s in st},
          "all verified:", all(r.get("drat_verified", False) for r in summary if r["status"] == "UNSAT"))

if __name__ == "__main__":
    main()
