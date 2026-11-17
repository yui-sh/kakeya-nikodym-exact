#!/bin/bash
# Session V: CP-SAT fallback for the LP-surviving configuration r_max=5, idx 4 (union 16), per degree cap D; starts when chain A is done.
cd "$(dirname "$0")"
until grep -q CHAIN_A_DONE out/lp_probe_cfgA.log; do sleep 30; done
for D in 4 5 6 7; do
  nice -n 10 .venv/bin/python rmax_model.py --q 5 --n 46 --r 5 --idx 4 --D $D --time 1800 --workers 1 >> out/rmax516_n46.log 2>&1
done
echo RMAX516_DONE >> out/rmax516_n46.log
