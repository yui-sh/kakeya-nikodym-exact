#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX PLANE10_FULL
g++ -O3 -std=c++17 -Wall -Wextra plane10_enum.cpp -o out/plane10_enum
.venv/bin/python -u plane10_validate.py > out/plane10_validate.log 2>&1
bash run_plane10_catalog.sh > out/plane10_catalog_driver.log 2>&1
.venv/bin/python -u plane10_catalog_check.py > out/plane10_catalog_check.log 2>&1
for leaf in 617 515; do
    .venv/bin/python -u plane10_solve.py "$leaf" > "out/plane10_${leaf}.log" 2>&1
    .venv/bin/python -u - "$leaf" > "out/plane10_${leaf}_final_check.log" 2>&1 <<'PY'
import json,sys
from plane10_cert import run
from plane10_check import check
short=sys.argv[1]
r=json.load(open(f'out/plane10_{short}.json'))
if r['status']=='LP-certified infeasible':run(short,replay=True)
elif r['status']=='computed_feasible_pending_exact':
    check(short);check(short,replay=True)
else:raise RuntimeError('No certified verdict: '+r['status'])
PY
done
.venv/bin/python -u plane10_diag.py > out/plane10_diag.log 2>&1
echo PLANE10_DONE
