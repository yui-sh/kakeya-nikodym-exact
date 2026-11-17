#!/usr/bin/env python3
"""Write specs/nik29*.json for the max-secant split (enc_nik29.py).  Run from anywhere."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from ff import GF
from enc_nik29 import axis_orbit_reps

def cases(q, K, outdir, expect, s_list=None, cube=None):
    reps = axis_orbit_reps(GF(q))
    out = []
    for s in (s_list or range(2, q + 1)):
        for o in range(len(reps[s])):
            for tan in ("h", "d"):
                name = f"q{q}K{K}_s{s}o{o}{tan}"
                gen = ["python3", "enc_nik29.py", "--q", str(q), "--K", str(K), "--s", str(s), "--orbit", str(o), "--tan", tan]
                if cube:
                    d = f"{outdir}/{name}"
                    out.append({"name": name + "_cover", "gen": gen + ["--cube", str(len(cube)), "--cube-points", ",".join(map(str, cube)), "--outdir", d, "--cover", f"{outdir}/{name}_cover.cnf"], "expect": "UNSAT"})
                    ncubes = (q + 1) ** len(cube)
                    for i in range(ncubes):
                        out.append({"name": f"{name}_cube{i}", "gen": ["cp", f"{d}/cube{i}.cnf", f"{outdir}/{name}_cube{i}.cnf"], "expect": expect})
                else:
                    out.append({"name": name, "gen": gen + ["--out", f"{outdir}/{name}.cnf"], "expect": expect})
    return out

specs = os.path.join(HERE, "..", "specs")
# q = 7 control of the encoder: Nik(2,7) = 36, max |C| = 13: all cases UNSAT at K = 14, some case SAT at K = 13
json.dump(cases(7, 14, "../big/nik29/ctrl7", "UNSAT") + cases(7, 13, "../big/nik29/ctrl7", "SAT-or-UNSAT"), open(os.path.join(specs, "nik29ctrl7.json"), "w"), indent=0)
# q = 9: UNSAT at K = 19 for every case; SAT control at K = 18 (s = 6 cases; the 3-row witness has s = 6)
json.dump(cases(9, 19, "../big/nik29", "UNSAT"), open(os.path.join(specs, "nik29.json"), "w"), indent=0)
json.dump(cases(9, 18, "../big/nik29", "SAT-or-UNSAT", s_list=[6]), open(os.path.join(specs, "nik29sat.json"), "w"), indent=0)
print("written")
