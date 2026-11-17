#!/bin/bash
# F7 final model: configuration probes with certificates at D = 5, 6, 7 for configs given as args
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
for D in 5 6 7; do for cfg in "$@"; do
  nice -n 5 .venv/bin/python lp_probe.py --n 46 --D $D --config $cfg --cert --tag f7_cfg_${cfg/,/_} >> out/f7_cfg_${cfg/,/_}.log 2>&1
done; done
echo D567_DONE >> out/f7_cfg_$1.log
