"""Session V: run the HiGHS MILP (all variables integer) on every degree-vector survivor listed in out/degvec_bb_*.json;
resumable JSONL out/milp_sweep.jsonl. Solver-certified verdicts (no exportable proof).  Usage: milp_sweep.py [--time 900]"""
import glob, json, os, sys, time, argparse, subprocess
ap = argparse.ArgumentParser(); ap.add_argument("--time", type=float, default=900); ap.add_argument("--only", default=None)
ap.add_argument("--D", type=int, nargs="*", default=None, help="restrict to these D"); ap.add_argument("--wait-for", nargs="*", default=[], help="log files whose DONE marker must appear before giving up")
a = ap.parse_args()
out = "out/milp_sweep.jsonl"; done = set()
if os.path.exists(out):
    for line in open(out):
        d = json.loads(line); done.add((d["config"], d["D"], tuple(d["mvec"])))
while True:
    todo = []
    todo_sorted = True
    for fn in sorted(glob.glob("out/degvec_bb_*_D*.json")):
        d = json.load(open(fn))
        for mv in d["survivors"]:
            key = (d["config"], d["D"], tuple(mv))
            if key not in done and (a.only is None or d["config"] == a.only) and (a.D is None or d["D"] in a.D): todo.append(key)
    if not todo:
        pending = [f for f in a.wait_for if not (os.path.exists(f) and "DONE" in open(f).read())]
        if not pending: break
        time.sleep(60); continue
    cfg, D, mv = todo[0]
    r = subprocess.run([".venv/bin/python", "milp_probe.py", str(D), cfg, ",".join(map(str, mv)), str(a.time)], capture_output=True, text=True)
    line = [l for l in r.stdout.splitlines() if l.startswith("{")]
    rec = json.loads(line[-1]) if line else dict(D=D, config=cfg, mvec=list(mv), status=-1, message="no output: " + r.stderr[-300:])
    with open(out, "a") as fo: fo.write(json.dumps(rec) + "\n")
    done.add((cfg, D, tuple(mv))); print(json.dumps(rec), flush=True)
print("MILP_SWEEP_DONE", flush=True)
