#!/bin/bash
# Session VI: exact pencil items (k exact up to f(r)+3), F6 flags for r=5 planes: LP probes with certificates, then MILP
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
for cfg in 5,15 6,17 5,16; do
  nice -n 5 .venv/bin/python lp_probe.py --n 46 --D 4 --config $cfg --cert --show --tag exact_cfg_${cfg/,/_} >> out/exact_cfg.log 2>&1
done
echo LP_DONE >> out/exact_cfg.log
for cfg in 5,15 6,17; do
  nice -n 5 .venv/bin/python milp_probe.py 4 $cfg - 2400 out/milpx_solution_${cfg/,/_}_D4.json > out/milpx_${cfg/,/_}_D4.log 2>&1
  echo MILP_DONE >> out/milpx_${cfg/,/_}_D4.log
done
echo EXACT_DONE >> out/exact_cfg.log
