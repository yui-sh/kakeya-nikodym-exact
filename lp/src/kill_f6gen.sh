#!/bin/bash
for p in $(pgrep -f "run_f6g""en.sh"; pgrep -f "milp_pr""obe.py 4"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo f6gen-killed
