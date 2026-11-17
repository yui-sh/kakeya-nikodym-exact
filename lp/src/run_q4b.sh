#!/bin/bash
# remaining 5 subcases of the q=4 campaign, shorter budget (default 840 s each, 3 in parallel)
cd "$(dirname "$0")"
T=${T:-840}
printf "meet 11\nskew 00\nskew 10\nskew 01\nskew 11\n" | \
xargs -P 3 -n 2 sh -c '.venv/bin/python solve_cpsat.py --n 3 --q 4 --fix 3 --case $0 --sub $1 --planes --core --ub 28 --time '$T' --workers 1 --out out/q4_$0_$1.json > out/q4_$0_$1.log 2>&1'
