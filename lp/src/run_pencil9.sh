#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX
# Run the whole driver with nice -n 10 (including compilation and child processes).
g++ -O3 -std=c++17 -Wall -Wextra pencil9_price.cpp -o out/pencil9_price
.venv/bin/python -u pencil9_validate.py > out/pencil9_validate.log 2>&1
.venv/bin/python -u pencil9.py --config 6,17 --m 0,3,27,12,4 --tag pencil9_617 > out/pencil9_617.log 2>&1
.venv/bin/python -u pencil9.py --config 5,15 --m 0,4,24,15,3 --tag pencil9_515 > out/pencil9_515.log 2>&1
.venv/bin/python -u pencil9_check.py pencil9_617 pencil9_515 > out/pencil9_check.log 2>&1
.venv/bin/python -u pencil9_check.py --rational-only pencil9_617 pencil9_515 > out/pencil9_rational_check.log 2>&1
echo PENCIL9_DONE
