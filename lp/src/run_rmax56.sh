#!/bin/bash
# Session V: the 12 r_max in {5,6} configurations at n=46 (session IV left them unrun), 1 core, 900 s each.
cd "$(dirname "$0")"
for r in 5 6; do for i in 0 1 2 3 4 5; do
  nice -n 10 .venv/bin/python rmax_model.py --q 5 --n 46 --r $r --idx $i --time 900 --workers 1
done; done > out/rmax56_n46.log 2>&1
echo RMAX56_DONE >> out/rmax56_n46.log
