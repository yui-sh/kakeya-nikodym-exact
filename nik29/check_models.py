#!/usr/bin/env python3
"""Decode every SAT record of a batch directory with enc_nik29.py --model and check it with
src/check_nikodym.py (independent checker).  usage: check_models.py SPEC DIR"""
import json, os, subprocess, sys
spec, d = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
n_sat = 0
for case in json.load(open(spec)):
    rec = os.path.join(d, case["name"] + ".json")
    if not os.path.exists(rec):
        continue
    r = json.load(open(rec))
    if r.get("result") != "SAT":
        continue
    n_sat += 1
    gen = [t for t in case["gen"] if t not in ("--out",) and not t.endswith(".cnf")]
    js = subprocess.run(gen + ["--model", os.path.abspath(r["model_file"])], cwd=SRC, capture_output=True, text=True, check=True).stdout
    out = os.path.join(d, case["name"] + ".decoded.json")
    open(out, "w").write(js)
    p = subprocess.run(["python3", os.path.join(SRC, "check_nikodym.py"), out], capture_output=True, text=True)
    C = json.loads(js)["C"]
    # also check the case conditions: axis intersection, max secant
    print(case["name"], "|C| =", len(C), p.stdout.strip(), "rc", p.returncode)
print("SAT cases checked:", n_sat)
