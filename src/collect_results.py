#!/usr/bin/env python3
"""Collect the per-case records of every batch into results/ and write results/RESULTS.md.

For each spec in specs/*.json: compare the recorded verdict with the expected one, require
cake_lpr verification for every UNSAT, decode every SAT model and re-check it with the independent
checker (check_kakeya.py / check_nikodym.py) including the size claim.  Copies the small JSON
records (no proofs) to results/<spec>/.  Exit 1 on any discrepancy.
"""
import glob, json, os, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")


def decode_and_check(case, rec, bigdir):
    gen = case["gen"]
    model = rec.get("model_file") or os.path.join(bigdir, case["name"] + ".cnf.model")
    if not os.path.exists(model):
        return "no model file"
    cmd = []
    skip = False
    for t in gen:
        if skip:
            skip = False
            continue
        if t == "--out":
            skip = True
            continue
        cmd.append(t)
    cmd += ["--model", model]
    out = subprocess.run(cmd, cwd=SRC, capture_output=True, text=True)
    if out.returncode != 0:
        return "decode failed: " + out.stderr[-200:]
    d = json.loads(out.stdout)
    dec = os.path.join(bigdir, case["name"] + ".decoded.json")
    json.dump(d, open(dec, "w"))
    if "points" in d:  # Kakeya
        r = subprocess.run(["python3", "check_kakeya.py", dec], cwd=SRC, capture_output=True, text=True)
        ok = r.returncode == 0
        n = len(d["points"])
        if "ub" in d:
            ok = ok and n <= d["ub"]
            return f"KAKEYA OK |S|={n}<={d['ub']}" if ok else f"CHECK FAILED: {r.stdout.strip()} |S|={n}"
        return f"KAKEYA OK |S|={n}" if ok else f"CHECK FAILED: {r.stdout.strip()}"
    if "X" in d:  # capacity positive control
        return f"capacity model: {d['size']} points, {d['directions_with_full_line']} directions"
    r = subprocess.run(["python3", "check_nikodym.py", dec], cwd=SRC, capture_output=True, text=True)
    ok = r.returncode == 0
    if "C" in d:
        size = len(d["C"])
        ok = ok and size >= d["K"]
        return f"NIKODYM OK |C|={size}>={d['K']}" if ok else f"CHECK FAILED: {r.stdout.strip()}"
    size = len(d["N"])
    ok = ok and size <= d["k"]
    return f"NIKODYM OK |N|={size}<={d['k']}" if ok else f"CHECK FAILED: {r.stdout.strip()}"


def main():
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    rows = []
    problems = 0
    for spec_path in sorted(glob.glob(os.path.join(ROOT, "specs", "*.json"))):
        spec_name = os.path.basename(spec_path)[:-5]
        if spec_name.startswith("nik29"):  # Nik(2,9): records live in big/nik29 (+ctrl7) and are collected by nik29/collect.py
            continue
        spec = json.load(open(spec_path))
        bigdir = os.path.join(ROOT, "big", spec_name)
        resdir = os.path.join(ROOT, "results", spec_name)
        os.makedirs(resdir, exist_ok=True)
        n_unsat = n_sat = n_missing = 0
        cad = trim = cake = 0.0
        maxlrat = 0
        details = []
        for case in spec:
            rp = os.path.join(bigdir, case["name"] + ".json")
            if not os.path.exists(rp):
                n_missing += 1
                details.append(f"| {case['name']} | {case.get('expect','')} | MISSING | | | | |")
                continue
            rec = json.load(open(rp))
            shutil.copy(rp, resdir)
            exp = case.get("expect")
            res = rec.get("result")
            note = ""
            if res == "UNSAT":
                if rec.get("verified"):
                    n_unsat += 1
                    note = "cake_lpr VERIFIED"
                else:
                    note = "NOT VERIFIED"
                    problems += 1
                maxlrat = max(maxlrat, rec.get("lrat_bytes", 0))
            elif res == "SAT":
                n_sat += 1
                note = decode_and_check(case, rec, bigdir)
                if "OK" not in note and "capacity model" not in note:
                    problems += 1
            else:
                note = str(res)
                problems += 1
            if exp and res != exp:
                note += f" **UNEXPECTED (expected {exp})**"
                problems += 1
            cad += rec.get("cadical_seconds", 0)
            trim += rec.get("trim_seconds", 0)
            cake += rec.get("cake_seconds", 0)
            details.append(f"| {case['name']} | {exp} | {res} | {rec.get('cadical_seconds','')} | {rec.get('trim_seconds','')} | {rec.get('cake_seconds','')} | {rec.get('lrat_bytes','')} | {note} |")
        rows.append((spec_name, len(spec), n_unsat, n_sat, n_missing, cad, trim, cake, maxlrat))
        with open(os.path.join(resdir, "cases.md"), "w") as f:
            f.write(f"# {spec_name}\n\n| case | expected | result | cadical s | lrat-trim s | cake_lpr s | LRAT bytes | note |\n|---|---|---|---|---|---|---|---|\n")
            f.write("\n".join(details) + "\n")
    with open(os.path.join(ROOT, "results", "RESULTS.md"), "w") as f:
        f.write("# Batch results (measured on this box, nice -n 19, one solver at a time)\n\n")
        f.write("| batch | cases | UNSAT verified by cake_lpr | SAT (model re-checked) | missing | cadical s | lrat-trim s | cake_lpr s | largest LRAT (bytes) |\n|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]:.0f} | {r[6]:.0f} | {r[7]:.0f} | {r[8]} |\n")
        f.write(f"\nproblems (unexpected verdicts, unverified UNSAT, failed model checks): {problems}\n")
        f.write("\nNik(2,9) = 63 (specs/nik29*.json, 1128 UNSAT + 4 SAT controls + 32 q = 7 controls) is collected separately by `nik29/collect.py` -> `nik29/cases.md`, `nik29/summary.json`; audited by `nik29/final_audit.py` -> `nik29/FINAL_AUDIT.md`.\n")
    print(open(os.path.join(ROOT, "results", "RESULTS.md")).read())
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
