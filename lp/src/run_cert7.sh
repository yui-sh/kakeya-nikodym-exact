#!/bin/bash
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python cert5_check.py out/lp_probe_f7_cfg_5_16.json --config 5,16 > out/cert6_check_516_all.log 2>&1
echo CERT7_DONE >> out/cert6_check_516_all.log
