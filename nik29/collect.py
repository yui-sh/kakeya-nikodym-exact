#!/usr/bin/env python3
"""Collect the Nik(2,9) records (big/nik29/*.json) per spec into nik29/records/<spec>/ and write
nik29/cases.md + nik29/summary.json.  Requires cake_lpr verification for every UNSAT; SAT models are
decoded and re-checked with src/check_nikodym.py.  usage: collect.py specs/nik29*.json ..."""
import json, os, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
BIG = os.path.join(ROOT, "big", "nik29")
summary = {}
md = ["# Nik(2,9) cases (max-secant split, enc_nik29.py)\n"]
bad = 0
for spec_path in sys.argv[1:]:
    spec_name = os.path.basename(spec_path)[:-5]
    spec = json.load(open(spec_path))
    outdir = os.path.join(HERE, "records", spec_name)
    os.makedirs(outdir, exist_ok=True)
    bigdir = BIG if "ctrl7" not in spec_name else os.path.join(BIG, "ctrl7")
    rows, tot = [], {"cases": len(spec), "verified": 0, "sat": 0, "missing": 0, "other": 0, "cadical": 0.0, "trim": 0.0, "cake": 0.0, "max_cadical": 0.0, "max_lrat": 0}
    for case in spec:
        rec_path = os.path.join(bigdir, case["name"] + ".json")
        if not os.path.exists(rec_path):
            tot["missing"] += 1
            rows.append((case["name"], case["expect"], "MISSING", "", "", "", "", "not run yet"))
            continue
        r = json.load(open(rec_path))
        shutil.copy(rec_path, outdir)
        note = ""
        if r.get("result") == "UNSAT" and r.get("verified"):
            tot["verified"] += 1
            note = "cake_lpr VERIFIED"
            tot["cadical"] += r["cadical_seconds"]; tot["trim"] += r.get("trim_seconds", 0); tot["cake"] += r["cake_seconds"]
            tot["max_cadical"] = max(tot["max_cadical"], r["cadical_seconds"]); tot["max_lrat"] = max(tot["max_lrat"], r["lrat_bytes"])
            if case["expect"] != "UNSAT" and case["expect"] != "SAT-or-UNSAT":
                bad += 1; note += " UNEXPECTED"
        elif r.get("result") == "SAT":
            tot["sat"] += 1
            gen = []
            skip = False
            for t in case["gen"]:
                if skip:
                    skip = False; continue
                if t == "--out":
                    skip = True; continue
                gen.append(t)
            out = subprocess.run(gen + ["--model", r["model_file"]], cwd=SRC, capture_output=True, text=True, check=True).stdout
            d = json.loads(out)
            dec = os.path.join(outdir, case["name"] + ".decoded.json")
            open(dec, "w").write(out)
            chk = subprocess.run(["python3", "check_nikodym.py", dec], cwd=SRC, capture_output=True, text=True)
            ok = chk.returncode == 0 and len(d["C"]) == d["K"]
            note = f"SAT, model decoded: {chk.stdout.strip().split(': ', 1)[-1]}" + ("" if ok else " CHECK FAILED")
            if not ok or case["expect"] == "UNSAT":
                bad += 1; note += " UNEXPECTED"
        else:
            tot["other"] += 1
            note = f"NOT VERIFIED: {r.get('result')}"
            bad += 1
        rows.append((case["name"], case["expect"], r.get("result"), r.get("cadical_seconds", ""), r.get("trim_seconds", ""), r.get("cake_seconds", ""), r.get("lrat_bytes", ""), note))
    summary[spec_name] = tot
    md.append(f"\n## {spec_name}\n\n{json.dumps(tot)}\n\n| case | expected | result | cadical s | lrat-trim s | cake_lpr s | LRAT bytes | note |\n|---|---|---|---|---|---|---|---|")
    for row in rows:
        md.append("| " + " | ".join(str(x) for x in row) + " |")
open(os.path.join(HERE, "cases.md"), "w").write("\n".join(md) + "\n")
json.dump(summary, open(os.path.join(HERE, "summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
print("DISCREPANCIES:", bad)
sys.exit(1 if bad else 0)
