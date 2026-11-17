#!/bin/bash
cd "$(dirname "$0")"
for cfg in 5,21 5,18 5,17 5,16 5,15; do
  nice -n 10 .venv/bin/python lp_probe.py --n 46 --config $cfg --cert --tag cfg_${cfg/,/_} >> out/lp_probe_cfgA.log 2>&1
done
echo CHAIN_A_DONE >> out/lp_probe_cfgA.log
