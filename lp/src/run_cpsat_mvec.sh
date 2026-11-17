#!/bin/bash
# Session V: CP-SAT with configuration + D + fixed degree vector on MILP-feasible leaves (gauge).
cd "$(dirname "$0")"
nice -n 10 .venv/bin/python rmax_model.py --q 5 --n 46 --r 6 --idx 5 --D 4 --mvec 0,5,19,22 --time 1800 --workers 1 >> out/cpsat_mvec.log 2>&1
nice -n 10 .venv/bin/python rmax_model.py --q 5 --n 46 --r 5 --idx 5 --D 4 --mvec 0,7,15,24 --time 1800 --workers 1 >> out/cpsat_mvec.log 2>&1
echo CPSAT_MVEC_DONE >> out/cpsat_mvec.log
