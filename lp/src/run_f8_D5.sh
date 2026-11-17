#!/bin/bash
# session VII: F8 LP probes (with certificates) at D = 5 for configs given as args
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1
for cfg in "$@"; do
  nice -n 5 .venv/bin/python lp_probe.py --n 46 --D 5 --config $cfg --cert --tag f8_cfg_${cfg/,/_} >> out/f8_cfg_${cfg/,/_}.log 2>&1
done
echo F8_D5_DONE >> out/f8_cfg_$1.log
