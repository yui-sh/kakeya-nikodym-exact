#!/bin/bash
# plain merged LP (no WLOG) at n=46 for every D with certificates, final model; then n=47.. if infeasible
cd "$(dirname "$0")"
until grep -q LP_DONE out/f7_cfg_5_16.log 2>/dev/null; do sleep 10; done
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python planecount6.py --q 5 --n 46 --D 4 5 6 7 8 9 10 11 --cert --tag f7 > out/f7_plain_n46.log 2>&1
echo PLAIN_DONE >> out/f7_plain_n46.log
