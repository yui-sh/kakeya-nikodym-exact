#!/bin/bash
for p in $(pgrep -f "run_f7_""D.sh 5,16"; pgrep -f "lp_pro""be.py --n 46 --D 7 --config 5,15"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo 515-chain-killed
