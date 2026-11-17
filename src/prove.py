#!/usr/bin/env python3
"""Run one CNF through the certified pipeline and record the result.

  cadical --lrat  ->  lrat-trim  ->  cake_lpr (formally verified LRAT checker, CakeML)

Writes a JSON record {cnf_sha256, result (UNSAT/SAT/UNKNOWN), cadical_seconds, trim_seconds,
cake_seconds, cake_verdict, lrat_bytes, model_file?} and deletes the proof files afterwards
(they are large; their sha256 is kept).  Everything runs under nice -n 19, single process.
usage: prove.py case.cnf record.json [--timeout S] [--tmp DIR] [--keep]
Exit 0 if cake_lpr printed 's VERIFIED UNSAT', 10 if SAT, 1 otherwise.
"""
import argparse, hashlib, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "repos", "bin")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, timeout=None, stdout=None):
    t = time.time()
    try:
        p = subprocess.run(["nice", "-n", "19"] + cmd, stdout=stdout or subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
        return p.returncode, (p.stdout.decode(errors="replace") if stdout is None else ""), time.time() - t
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT", time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cnf")
    ap.add_argument("record")
    ap.add_argument("--timeout", type=float, default=None)
    ap.add_argument("--tmp", default=None)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    tmp = a.tmp or os.path.dirname(os.path.abspath(a.cnf))
    base = os.path.join(tmp, os.path.basename(a.cnf))
    lrat, trimmed, model = base + ".lrat", base + ".trim.lrat", base + ".model"
    rec = {"cnf": os.path.basename(a.cnf), "cnf_sha256": sha256(a.cnf), "tools": {k: sha256(os.path.join(BIN, k)) for k in ("cadical", "lrat-trim", "cake_lpr")}}
    rc, out, t = run([os.path.join(BIN, "cadical"), "-q", "--lrat", "--no-binary", a.cnf, lrat], a.timeout)
    rec["cadical_seconds"] = round(t, 2)
    if rc == 10:
        rec["result"] = "SAT"
        with open(model, "w") as f:
            f.write(out)
        rec["model_file"] = model
        for p in (lrat,):
            if os.path.exists(p) and not a.keep:
                os.remove(p)
        json.dump(rec, open(a.record, "w"), indent=1)
        print(json.dumps(rec))
        sys.exit(10)
    if rc != 20:
        rec["result"] = "UNKNOWN" if rc is None else f"cadical_rc={rc}"
        rec["cadical_output"] = out[-2000:]
        if os.path.exists(lrat) and not a.keep:
            os.remove(lrat)
        json.dump(rec, open(a.record, "w"), indent=1)
        print(json.dumps(rec))
        sys.exit(1)
    rec["result"] = "UNSAT"
    rec["lrat_bytes"] = os.path.getsize(lrat)
    rec["lrat_sha256"] = sha256(lrat)
    rc, out, t = run([os.path.join(BIN, "lrat-trim"), lrat, trimmed])
    rec["trim_seconds"] = round(t, 2)
    rec["trim_rc"] = rc
    if rc == 0:
        rec["trim_lrat_bytes"] = os.path.getsize(trimmed)
        rec["trim_lrat_sha256"] = sha256(trimmed)
        proof = trimmed
    else:  # fall back to checking the untrimmed proof
        rec["trim_output"] = out[-1000:]
        proof = lrat
    rc, out, t = run([os.path.join(BIN, "cake_lpr"), a.cnf, proof])
    rec["cake_seconds"] = round(t, 2)
    rec["cake_verdict"] = out.strip()[-200:]
    rec["verified"] = (rc == 0 and "s VERIFIED UNSAT" in out)
    if not a.keep:
        for p in (lrat, trimmed):
            if os.path.exists(p):
                os.remove(p)
    json.dump(rec, open(a.record, "w"), indent=1)
    print(json.dumps({k: rec[k] for k in ("cnf", "result", "cadical_seconds", "trim_seconds", "cake_seconds", "verified", "lrat_bytes")}))
    sys.exit(0 if rec["verified"] else 1)


if __name__ == "__main__":
    main()
