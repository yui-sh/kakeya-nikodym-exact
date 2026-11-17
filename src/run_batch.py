#!/usr/bin/env python3
"""Sequential, resumable driver: generate each CNF of a case list and run prove.py on it.

spec.json: [{"name": "...", "gen": ["python3", "enc_x.py", ...args writing <dir>/<name>.cnf]}, ...]
Records go to <dir>/<name>.json; cases with an existing record whose 'verified' is true
(or result SAT) are skipped, so the batch can be resumed.  Before every case the driver
waits while the 1-minute load average exceeds --max-load (shared box) and, if
--quota-script is given, stops when the Claude usage quota is above the configured limits.
usage: run_batch.py spec.json DIR [--timeout S] [--max-load 6] [--quota-script PATH]
"""
import argparse, json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))


def quota_ok(script):
    try:
        out = subprocess.run(["python3", script], capture_output=True, text=True, timeout=60).stdout
    except Exception as e:  # noqa
        print("quota check failed:", e)
        return True
    m = re.search(r"claude.*?wk ([0-9.]+)% used.*?5h ([0-9.]+)%", out)
    if not m:
        print("quota line not parsed:", out[:200])
        return True
    wk, h5 = float(m.group(1)), float(m.group(2))
    print(f"quota: weekly {wk}% 5h {h5}%")
    return h5 < 65 and wk < 60


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("dir")
    ap.add_argument("--timeout", type=float, default=None)
    ap.add_argument("--max-load", type=float, default=6.0)
    ap.add_argument("--quota-script", default=None)
    ap.add_argument("--summary", default=None)
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    spec = json.load(open(a.spec))
    results = []
    for case in spec:
        name = case["name"]
        rec_path = os.path.join(a.dir, name + ".json")
        cnf = os.path.join(a.dir, name + ".cnf")
        if os.path.exists(rec_path):
            rec = json.load(open(rec_path))
            if rec.get("verified") or rec.get("result") == "SAT":
                results.append(rec)
                print("skip (done):", name, rec.get("result"))
                continue
        while os.getloadavg()[0] > a.max_load:
            print(f"load {os.getloadavg()[0]:.2f} > {a.max_load}, waiting 60 s", flush=True)
            time.sleep(60)
        if a.quota_script and not quota_ok(a.quota_script):
            print("QUOTA LIMIT REACHED - stopping batch before", name, flush=True)
            break
        t = time.time()
        subprocess.run(case["gen"] + [], check=True, cwd=HERE, stdout=subprocess.DEVNULL)
        print(f"[{time.strftime('%H:%M:%S')}] {name}: generated in {time.time()-t:.1f}s", flush=True)
        cmd = ["python3", os.path.join(HERE, "prove.py"), cnf, rec_path]
        if a.timeout:
            cmd += ["--timeout", str(a.timeout)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        print(p.stdout.strip()[-400:], flush=True)
        if p.returncode not in (0, 10):
            print("  FAILED/UNKNOWN:", name, p.stderr[-500:], flush=True)
        rec = json.load(open(rec_path))
        results.append(rec)
        if os.path.exists(cnf) and rec.get("result") == "UNSAT" and rec.get("verified"):
            os.remove(cnf)
    if a.summary:
        json.dump(results, open(a.summary, "w"), indent=1)
    nver = sum(1 for r in results if r.get("verified"))
    nsat = sum(1 for r in results if r.get("result") == "SAT")
    print(f"BATCH DONE: {len(results)}/{len(spec)} cases run, {nver} UNSAT verified by cake_lpr, {nsat} SAT, "
          f"total cadical {sum(r.get('cadical_seconds',0) for r in results):.0f}s trim {sum(r.get('trim_seconds',0) for r in results):.0f}s "
          f"cake {sum(r.get('cake_seconds',0) for r in results):.0f}s")


if __name__ == "__main__":
    main()
