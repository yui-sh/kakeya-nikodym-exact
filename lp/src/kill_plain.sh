#!/bin/bash
for p in $(pgrep -f "run_f7_pl""ain.sh"; pgrep -f "planecount""6.py --q 5 --n 46"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo plain-killed
