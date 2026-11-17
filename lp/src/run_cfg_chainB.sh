#!/bin/bash
cd "$(dirname "$0")"
for cfg in 6,25 6,21 6,19 6,18 6,17; do
  nice -n 10 .venv/bin/python lp_probe.py --n 46 --config $cfg --cert --tag cfg_${cfg/,/_} >> out/lp_probe_cfgB.log 2>&1
done
echo CHAIN_B_DONE >> out/lp_probe_cfgB.log
