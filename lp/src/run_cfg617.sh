#!/bin/bash
# Session V: CP-SAT fallback for r_max=6, idx 5 (union 17, Blokhuis-Mazzocca plane), D=4..7; starts after the (5,15) run.
cd "$(dirname "$0")"
until grep -q RMAX515_DONE out/rmax515_n46.log; do sleep 30; done
for D in 4 5 6 7; do
  nice -n 10 .venv/bin/python rmax_model.py --q 5 --n 46 --r 6 --idx 5 --D $D --time 1800 --workers 1 >> out/rmax617_n46.log 2>&1
done
echo RMAX617_DONE >> out/rmax617_n46.log
