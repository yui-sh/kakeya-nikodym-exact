#!/bin/bash
for p in $(pgrep -f "run_f7_""D.sh 6,17"; pgrep -f "lp_pro""be.py --n 46 --D 7 --config 6,17"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo 617-chain-killed
