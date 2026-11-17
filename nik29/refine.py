#!/usr/bin/env python3
"""Depth-3 refinement spec for cubes that timed out (result UNKNOWN) in the depth-2 cube batches.

For a timed-out cube i of case NAME (cube points p1,p2) the ten depth-3 cubes with points p1,p2,EXTRA and
indices 10*i .. 10*i+9 are exactly 'cube i + one witness choice for EXTRA' (checked here by comparing the
literal lists in the two cubes.json files once the depth-3 cubes are generated), so UNSAT of the ten
children, together with the depth-3 cover CNF, closes cube i.  Writes specs/nik29cubes_refine.json.
usage: refine.py [EXTRA_POINT_INDEX=10]      then  nik29/run_lanes_staged.sh specs/nik29cubes_refine.json nik29/cubes_refine 1800
"""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIG = os.path.join(ROOT, "big", "nik29")
extra = int(sys.argv[1]) if len(sys.argv) > 1 else 10
out = []
for sp in ("nik29cubes_s4.json", "nik29cubes_s5.json"):
    spec = json.load(open(os.path.join(ROOT, "specs", sp)))
    covers = {c["name"].rsplit("_", 1)[0]: c for c in spec if c["name"].endswith("_cover")}
    timed = {}
    for c in spec:
        p = os.path.join(BIG, c["name"] + ".json")
        if c["name"].endswith("_cover") or not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("verified") and r.get("result") != "SAT":
            case, cube = c["name"].rsplit("_cube", 1)
            timed.setdefault(case, []).append(int(cube))
    for case, cubes in sorted(timed.items()):
        gen = list(covers[case]["gen"])
        i = gen.index("--cube-points")
        pts = gen[i + 1] + f",{extra}"
        gen[gen.index("--cube") + 1] = "3"
        gen[i + 1] = pts
        d3 = f"../big/nik29/cubes3/{case}"
        gen[gen.index("--outdir") + 1] = d3
        gen[gen.index("--cover") + 1] = f"../big/nik29/{case}_cover3.cnf"
        # generate now so the parent/child literal containment can be checked
        subprocess.run(gen, cwd=os.path.join(ROOT, "src"), check=True, stdout=subprocess.DEVNULL)
        parent = {c["name"]: c["cube"] for c in json.load(open(os.path.join(BIG, "cubes", case, "cubes.json")))}
        child = {c["name"]: c["cube"] for c in json.load(open(os.path.join(BIG, "cubes3", case, "cubes.json")))}
        out.append({"name": f"{case}_cover3", "gen": gen, "expect": "UNSAT"})
        for i in cubes:
            for j in range(10):
                k = 10 * i + j
                assert child[f"cube{k}"][:2] == parent[f"cube{i}"], (case, i, k)
                out.append({"name": f"{case}_cube3_{k}", "gen": ["cp", f"{d3}/cube{k}.cnf", f"../big/nik29/{case}_cube3_{k}.cnf"], "expect": "UNSAT"})
        print(case, "timed-out cubes:", cubes, "-> refined to", 10 * len(cubes), "depth-3 cubes")
json.dump(out, open(os.path.join(ROOT, "specs", "nik29cubes_refine.json"), "w"), indent=0)
print(len(out), "entries written to specs/nik29cubes_refine.json")
