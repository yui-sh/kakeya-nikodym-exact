#!/bin/bash
# Verification-only entry point for later sessions; launches no search solver.
# The XIV campaign itself used run14.sh with its fixed overnight hard stop.
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX PLANE10_FULL
case "${1:-all}" in
  5) exec nice -n 10 .venv/bin/python -u replay14_D5.py --fresh-prices ;;
  6|7) exec nice -n 10 .venv/bin/python -u coverage14.py --D "$1" --replay --fresh-prices ;;
  all)
    nice -n 10 .venv/bin/python -u replay14_D5.py --fresh-prices
    nice -n 10 .venv/bin/python -u coverage14.py --D 6 --replay --fresh-prices
    nice -n 10 .venv/bin/python -u coverage14.py --D 7 --replay --fresh-prices
    nice -n 10 .venv/bin/python -u final14.py
    ;;
  manifest) exec nice -n 10 .venv/bin/python -u final14.py ;;
  *) echo 'Usage: bash replay14.sh [all|5|6|7|manifest]' >&2; exit 2 ;;
esac
