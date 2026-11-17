#!/bin/bash
# Session VI: q=5, n=46 configuration probes with F5+F6, exact certificates
cd "$(dirname "$0")"
for cfg in 5,15 6,17 5,16; do
  for D in 4 5 6 7; do
    PC6_SECTIONS=1 nice -n 5 .venv/bin/python lp_probe.py --n 46 --D $D --config $cfg --cert --show --tag sec6_cfg_${cfg/,/_} >> out/sec6_cfg.log 2>&1
  done
done
echo SEC6_CFG_DONE >> out/sec6_cfg.log
