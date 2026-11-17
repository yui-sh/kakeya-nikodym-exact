#!/bin/bash
# after the exact MILP chain: LP probes (with certificates) for D=5,6,7 on the three configurations with the exact model
cd "$(dirname "$0")"
until grep -q EXACT_DONE out/exact_cfg.log 2>/dev/null; do sleep 20; done
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
for D in 5 6 7; do for cfg in 5,16 5,15 6,17; do
  nice -n 5 .venv/bin/python lp_probe.py --n 46 --D $D --config $cfg --cert --tag exact_cfg_${cfg/,/_} >> out/exact_cfg_D567.log 2>&1
done; done
echo EXACT_D_DONE >> out/exact_cfg_D567.log
