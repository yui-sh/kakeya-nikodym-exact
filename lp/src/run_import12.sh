#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX PLANE10_FULL
remaining=$(( $(date -u -d '2026-11-04 11:29:40' +%s) - $(date -u +%s) ))
if (( remaining <= 0 )); then exit 1; fi
timeout --signal=TERM --kill-after=10s "$remaining" .venv/bin/python -u import12.py
