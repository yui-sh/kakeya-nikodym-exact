#!/bin/bash
for p in $(pgrep -f "run_exact_""D.sh"; pgrep -f "run_r""4b.sh"; pgrep -f "lp_pro""be.py --n 46 --D 5"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo d567-killed
