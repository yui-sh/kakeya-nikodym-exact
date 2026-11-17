#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_UNMARKED=1
unset PC6_MFIX
g++ -O3 -std=c++17 -Wall -Wextra plane10_enum.cpp -o out/plane10_enum
g++ -O3 -std=c++17 -Wall -Wextra plane10_catalog.cpp -o out/plane10_catalog
.venv/bin/python -u - <<'PY'
from plane10 import build,PlaneOracle
for short,config,m in [('617',(6,17),(0,3,27,12,4)),('515',(5,15),(0,4,24,15,3))]:
    lp,inf,ts,gr,zr,nr,pgr=build(config,m)
    PlaneOracle(ts,inf,nr,pgr,'plane10_catalog_'+short).price([0]*len(lp.rows),first=True)
PY
out/plane10_catalog out/plane10_catalog_617_plane.in out/plane10_catalog_617.bin 2> out/plane10_catalog_617.log
out/plane10_catalog out/plane10_catalog_515_plane.in out/plane10_catalog_515.bin 2> out/plane10_catalog_515.log
echo CATALOGS_DONE
