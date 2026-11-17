#!/usr/bin/env python3
"""Progress of the cube batches: per case, cubes verified / total, timeouts, CPU so far.
usage: progress.py [specs/nik29cubes_s4.json specs/nik29cubes_s5.json]"""
import glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIG = os.path.join(ROOT, "big", "nik29")
specs = sys.argv[1:] or [os.path.join(ROOT, "specs", "nik29cubes_s4.json"), os.path.join(ROOT, "specs", "nik29cubes_s5.json")]
grand = {"verified": 0, "total": 0, "cadical": 0.0, "trim": 0.0, "cake": 0.0, "timeouts": [], "unverified": []}
for sp in specs:
    cases = {}
    for c in json.load(open(sp)):
        cases.setdefault(c["name"].rsplit("_", 1)[0], []).append(c["name"])
    print(os.path.basename(sp))
    for case, names in sorted(cases.items()):
        v = t = 0
        cad = trim = cake = 0.0
        mx = 0.0
        mxl = 0
        for n in names:
            p = os.path.join(BIG, n + ".json")
            if not os.path.exists(p):
                continue
            r = json.load(open(p))
            if r.get("verified"):
                v += 1
                cad += r["cadical_seconds"]; trim += r.get("trim_seconds", 0); cake += r["cake_seconds"]
                mx = max(mx, r["cadical_seconds"]); mxl = max(mxl, r["lrat_bytes"])
            elif r.get("result") == "UNKNOWN":
                t += 1; grand["timeouts"].append(n)
            else:
                grand["unverified"].append(n)
        grand["verified"] += v; grand["total"] += len(names); grand["cadical"] += cad; grand["trim"] += trim; grand["cake"] += cake
        print(f"  {case}: {v}/{len(names)} verified, {t} timeouts; cadical {cad:.0f} s (max {mx:.0f}), trim {trim:.0f} s, cake {cake:.0f} s, max LRAT {mxl/1e9:.2f} GB")
print(json.dumps(grand))
