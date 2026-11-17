#!/bin/bash
# final model: whole-configuration MILP at D=5 for config $1 (time limit $2)
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python milp_probe.py 5 $1 - $2 out/milpf7_solution_${1/,/_}_D5.json > out/milpf7_${1/,/_}_D5.log 2>&1
echo MILP_DONE >> out/milpf7_${1/,/_}_D5.log
