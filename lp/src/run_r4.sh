#!/bin/bash
# Session VI: after CP-SAT frees its core: MILP (5,15) D=4 with the R4 capacity cut, then validation of the final model
cd "$(dirname "$0")"
until grep -q CPSAT6_DONE out/cpsat6_515_D4_3_1_18_24.log 2>/dev/null; do sleep 20; done
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python milp_probe.py 4 5,15 - 2400 out/milpr4_solution_5_15_D4.json > out/milpr4_5_15_D4.log 2>&1
echo MILP_DONE >> out/milpr4_5_15_D4.log
nice -n 5 .venv/bin/python pc6_validate.py solutions/best_n3_q4_28_dvir.json solutions/best_n3_q3_13.json solutions/best_n3_q5_53.json solutions/bc_n3_q5_53.json > out/r4_validate.log 2>&1
echo VALIDATE_DONE >> out/r4_validate.log
