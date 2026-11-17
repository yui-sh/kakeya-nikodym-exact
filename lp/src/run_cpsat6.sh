#!/bin/bash
# Session VI: CP-SAT on the (5,15) configuration, D=4, degree vector of the F5+F6 integer solution (construction attempt / leaf test)
cd "$(dirname "$0")"
nice -n 10 .venv/bin/python rmax_model.py --q 5 --n 46 --r 5 --idx 5 --D 4 --mvec 3,1,18,24 --time 3000 --workers 1 > out/cpsat6_515_D4_3_1_18_24.log 2>&1
echo CPSAT6_DONE >> out/cpsat6_515_D4_3_1_18_24.log
