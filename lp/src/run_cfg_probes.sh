#!/bin/bash
# Session V: LP probes "r_max = r and some plane carries an r-line configuration with union u", n=46, all D, exact certs.
cd "$(dirname "$0")"
for cfg in 5,21 5,18 5,17 5,16 5,15 6,25 6,21 6,19 6,18 6,17; do
  nice -n 10 .venv/bin/python lp_probe.py --n 46 --config $cfg --cert --tag cfg_${cfg/,/_} >> out/lp_probe_cfg.log 2>&1
done
echo CFG_PROBES_DONE >> out/lp_probe_cfg.log
