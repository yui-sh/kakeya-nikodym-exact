#!/bin/bash
# Session V: standalone exact re-verification of every stored LP certificate (run after the campaign is stopped).
cd "$(dirname "$0")"
{ echo "== rmax<=4 lemma"; nice -n 10 .venv/bin/python cert5_check.py out/lp_probe_rmax4.json --rmax 4;
  for cfg in 5,21 5,18 5,17 5,16 5,15 6,25 6,21 6,19 6,18 6,17; do echo "== config $cfg"; nice -n 10 .venv/bin/python cert5_check.py out/lp_probe_cfg_${cfg/,/_}.json --config $cfg; done;
  for f in out/degvec_bb_*_D*.json; do echo "== tree $f"; c=$(basename $f .json); c=${c#degvec_bb_}; cfg=${c%_D*}; cfg=${cfg/_/,}; nice -n 10 .venv/bin/python cert5_check.py $f --bb --config $cfg; done; echo CERTCHECK_DONE; } > out/cert5_check.log 2>&1
