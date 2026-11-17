#!/bin/bash
# Session V: MILP sweep, D=4 leaves first, then everything else; keeps polling while the tree chains run.
cd "$(dirname "$0")"
nice -n 10 .venv/bin/python milp_sweep.py --time 900 --D 4 >> out/milp_sweep.log 2>&1
nice -n 10 .venv/bin/python milp_sweep.py --time 900 --wait-for out/degvec_bb_Y.log out/degvec_bb_X.log >> out/milp_sweep.log 2>&1
echo SWEEP2_DONE >> out/milp_sweep.log
