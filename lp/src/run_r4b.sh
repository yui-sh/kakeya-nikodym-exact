#!/bin/bash
# after the D=5..7 exact probes: MILP (6,17) D=4 with the R4 capacity cut
cd "$(dirname "$0")"
until grep -q EXACT_D_DONE out/exact_cfg_D567.log 2>/dev/null; do sleep 20; done
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python milp_probe.py 4 6,17 - 2400 out/milpr4_solution_6_17_D4.json > out/milpr4_6_17_D4.log 2>&1
echo MILP_DONE >> out/milpr4_6_17_D4.log
