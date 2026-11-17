#!/bin/bash
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python milp_probe.py 4 6,17 - 2400 out/milpr4_solution_6_17_D4.json > out/milpr4_6_17_D4.log 2>&1
echo MILP_DONE >> out/milpr4_6_17_D4.log
