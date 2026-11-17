#!/bin/bash
# wait for the CP-SAT leaf run to finish, then validate the generalised F6 model on the known sets (1 core)
cd "$(dirname "$0")"
until grep -q CPSAT6_DONE out/cpsat6_515_D4_3_1_18_24.log 2>/dev/null; do sleep 20; done
PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 nice -n 5 .venv/bin/python pc6_validate.py solutions/best_n3_q4_28_dvir.json solutions/best_n3_q3_13.json solutions/best_n3_q5_53.json solutions/bc_n3_q5_53.json > out/f6gen_validate.log 2>&1
echo VALIDATE_DONE >> out/f6gen_validate.log
