#!/bin/bash
# Session V: standalone exact re-verification of every stored certificate; starts when chain B is done.
cd "$(dirname "$0")"
until grep -q CHAIN_B_DONE out/lp_probe_cfgB.log; do sleep 30; done
until grep -q CHAIN_A_DONE out/lp_probe_cfgA.log; do sleep 30; done
{ echo "== rmax<=4 lemma"; nice -n 10 .venv/bin/python cert5_check.py out/lp_probe_rmax4.json --rmax 4;
  for cfg in 5,21 5,18 5,17 5,16 5,15 6,25 6,21 6,19 6,18 6,17; do echo "== config $cfg"; nice -n 10 .venv/bin/python cert5_check.py out/lp_probe_cfg_${cfg/,/_}.json --config $cfg; done; echo CERTCHECK_DONE; } > out/cert5_check.log 2>&1
