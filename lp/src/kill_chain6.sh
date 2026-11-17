#!/bin/bash
for p in $(pgrep -f "run_sec6""_cfg.sh"; pgrep -f "lp_pro""be.py --n 46"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo chain-killed
