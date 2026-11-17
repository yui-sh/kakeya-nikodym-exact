#!/bin/bash
# Session VI: generalised F6 (r>=2 flags): validate, LP probe (5,15) D=4, MILP (5,15) D=4 and (6,17) D=4
cd "$(dirname "$0")"
export PC6_SECTIONS=1
nice -n 5 .venv/bin/python pc6_validate.py solutions/best_n3_q4_28_dvir.json solutions/best_n3_q3_13.json solutions/best_n3_q5_53.json solutions/bc_n3_q5_53.json > out/f6gen_validate.log 2>&1
echo VALIDATE_DONE >> out/f6gen_validate.log
nice -n 5 .venv/bin/python lp_probe.py --n 46 --D 4 --config 5,15 --cert --show --tag f6gen_cfg_5_15 > out/f6gen_cfg.log 2>&1
echo LP_DONE >> out/f6gen_cfg.log
nice -n 5 .venv/bin/python milp_probe.py 4 5,15 - 2400 out/milp6g_solution_5_15_D4.json > out/milp6g_5_15_D4.log 2>&1
echo MILP_DONE >> out/milp6g_5_15_D4.log
nice -n 5 .venv/bin/python milp_probe.py 4 6,17 - 2400 out/milp6g_solution_6_17_D4.json > out/milp6g_6_17_D4.log 2>&1
echo MILP_DONE >> out/milp6g_6_17_D4.log
echo F6GEN_DONE >> out/f6gen_cfg.log
