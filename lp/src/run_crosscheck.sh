#!/bin/bash
# session III CP-SAT cross-check (second encoding: native model + planar + projection cuts, core mode, 1 worker):
#   D=4 subcases split into 16 sub-subcases each by the witness line of direction index 20 = (1,3,3); D=5 subcases whole.
# (pure-feasibility CP-SAT without the split was UNKNOWN after 3600 s on D=4 rep 0: out/md_D4_r0_cpsat.json)
cd "$(dirname "$0")"
for r in 0 1 2; do for k in $(seq 0 15); do
  nice -n 10 .venv/bin/python maxdeg.py --q 4 --ub 27 --D 4 --rep $r --time 1200 --workers 1 --planes --proj --core --split 20:$k --out out/md_D4_r${r}_s${k}_cpsat.json
done; done
for r in 0 1 2 3 4; do
  nice -n 10 .venv/bin/python maxdeg.py --q 4 --ub 27 --D 5 --rep $r --time 1200 --workers 1 --planes --proj --core --out out/md_D5_r${r}_cpsat.json
done
echo CPSAT_DONE
