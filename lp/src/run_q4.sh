#!/bin/bash
# n=3 q=4: 8 symmetry subcases, core-based CP-SAT, 1 worker each, 3 in parallel.
cd "$(dirname "$0")"
T=${T:-1500}
for c in meet skew; do for s in 00 10 01 11; do echo "$c $s"; done; done | \
xargs -P 3 -n 2 sh -c '.venv/bin/python solve_cpsat.py --n 3 --q 4 --fix 3 --case $0 --sub $1 --planes --core --ub 28 --time '$T' --workers 1 --out out/q4_$0_$1.json > out/q4_$0_$1.log 2>&1'
