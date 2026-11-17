#!/bin/bash
# Session VI: HiGHS MILP (all integer) on whole configurations at D=4 with F5+F6; $1 = config r,u ; $2 = time limit
cd "$(dirname "$0")"
PC6_SECTIONS=1 nice -n 5 .venv/bin/python milp_probe.py 4 $1 - $2 out/milp6_solution_${1/,/_}_D4.json > out/milp6_${1/,/_}_D4.log 2>&1
echo MILP6_DONE >> out/milp6_${1/,/_}_D4.log
