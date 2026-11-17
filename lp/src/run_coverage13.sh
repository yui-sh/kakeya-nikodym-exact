#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX PLANE10_FULL
export COVERAGE_SESSION=13
.venv/bin/python -u lift11.py > out/lift13_base_validation.log 2>&1
remaining=$(( $(date -u -d '2026-11-04 11:14:30' +%s) - $(date -u +%s) ))
if (( remaining <= 0 )); then echo HARD_STOP; exit 1; fi
if [[ "${1:-}" == "--joint" ]]; then
    shift
    timeout --signal=TERM --kill-after=10s "$remaining" .venv/bin/python -u joint11.py "$@"
else
    timeout --signal=TERM --kill-after=10s "$remaining" .venv/bin/python -u coverage11.py "$@"
fi
