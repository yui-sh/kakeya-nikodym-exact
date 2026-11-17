#!/bin/bash
# Session VI: F7 (pencil degree cut) + exact items + F5/F6/R4: LP probes with certificates at D=4, configs $@
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
for cfg in "$@"; do
  nice -n 5 .venv/bin/python lp_probe.py --n 46 --D 4 --config $cfg --cert --show --tag f7_cfg_${cfg/,/_} >> out/f7_cfg_${cfg/,/_}.log 2>&1
  echo LP_DONE >> out/f7_cfg_${cfg/,/_}.log
done
