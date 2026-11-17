#!/bin/bash
for p in $(pgrep -f "run_r""4.sh"; pgrep -f "run_r""4c.sh"; pgrep -f "milp_pro""be.py 4 6,17"; pgrep -f "pc6_valid""ate.py"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo r4-jobs-killed
