#!/bin/bash
# Session VI: q=5, n=46 configuration probes with F5 (sections), exact certificates
cd "$(dirname "$0")"
for cfg in 5,15 6,17 5,16; do
  PC6_SECTIONS=1 nice -n 5 .venv/bin/python lp_probe.py --n 46 --D 4 5 6 7 --config $cfg --cert --show --tag sec_cfg_${cfg/,/_} >> out/sec_cfg.log 2>&1
done
echo SEC_CFG_DONE >> out/sec_cfg.log
