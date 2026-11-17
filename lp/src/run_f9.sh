#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX
nice -n 10 .venv/bin/python -u f9_validate.py solutions/best_n3_q3_13.json solutions/best_n3_q4_28_dvir.json solutions/best_n3_q5_53.json solutions/bc_n3_q5_53.json > out/f9_validate.log 2>&1
nice -n 10 .venv/bin/python -u f9_probe.py
