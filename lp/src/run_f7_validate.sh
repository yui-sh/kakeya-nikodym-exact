#!/bin/bash
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5
nice -n 5 .venv/bin/python pc6_validate.py solutions/best_n3_q3_13.json solutions/best_n3_q4_28_dvir.json solutions/best_n3_q5_53.json solutions/bc_n3_q5_53.json > out/f7_validate.log 2>&1
echo VALIDATE_DONE >> out/f7_validate.log
