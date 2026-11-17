#!/bin/bash
# session III independent path: (1) cadical195 (pysat) on the D=4,5 subcase CNFs; (2) session-I 3-line WLOG, 8 subcases, kissat + drat-trim
cd "$(dirname "$0")"
nice -n 10 .venv/bin/python - <<'PY' > out/cadical_D45.log 2>&1
import time, json
from pysat.formula import CNF
from pysat.solvers import Solver
for D, n in [(4, 3), (5, 5)]:
    for r in range(n):
        fn = f"out/md_q4_ub27_D{D}_r{r}.cnf"
        f = CNF(from_file=fn); t0 = time.time()
        with Solver(name="cadical195", bootstrap_with=f.clauses) as s:
            res = s.solve()
        print(json.dumps(dict(cnf=fn, solver="cadical195 (pysat)", sat=res, time=round(time.time() - t0, 1))), flush=True)
PY
for c in meet skew; do for s in 00 10 01 11; do
  .venv/bin/python wlog3_cnf.py --q 4 --ub 27 --case $c --sub $s --cnf out/w3_${c}_${s}.cnf
  /usr/bin/time -f "kissat wall %e s maxrss %M KB" nice -n 10 bin/kissat -q --time=5400 out/w3_${c}_${s}.cnf out/w3_${c}_${s}.drat 2>&1 | grep -E "^s |wall"
  if grep -q . out/w3_${c}_${s}.drat 2>/dev/null; then /usr/bin/time -f "drat-trim wall %e s" nice -n 10 bin/drat-trim out/w3_${c}_${s}.cnf out/w3_${c}_${s}.drat 2>&1 | grep -E "^s |wall"; fi
done; done > out/wlog3_kissat.log 2>&1
