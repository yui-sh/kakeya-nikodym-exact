#!/usr/bin/env python3
"""Run a vendored kakeya-fq script while logging every file it opens (audit hook), so the exact
input closure of the LP replay is recorded.  usage: trace_open.py LOG script.py [args...]"""
import runpy, sys, os
log = open(sys.argv[1], "a")
seen = set()


def hook(event, args):
    if event == "open":
        p = str(args[0])
        if p not in seen and not p.startswith("/usr") and ".venv" not in p and "/proc" not in p:
            seen.add(p)
            log.write(f"{args[1]}\t{p}\n")
            log.flush()


sys.addaudithook(hook)
sys.argv = sys.argv[2:]
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
runpy.run_path(sys.argv[0], run_name="__main__")
