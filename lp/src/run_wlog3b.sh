#!/bin/bash
# session III, independent proof path (no degree reduction): session-I 3-line WLOG, 8 subcases,
# kissat (DRAT) + drat-trim; one DRAT file on disk at a time (they reach GBs). Resumable: skips cases with a .res verdict.
cd "$(dirname "$0")"
while pgrep -f "w3_meet_00.cnf" >/dev/null; do sleep 20; done     # let an already-running meet/00 kissat finish
for c in meet skew; do for s in 00 10 01 11; do
  f=out/w3_${c}_${s}
  if [ ! -s $f.cnf ]; then .venv/bin/python wlog3_cnf.py --q 4 --ub 27 --case $c --sub $s --cnf $f.cnf; fi
  if ! grep -q "^s " $f.res 2>/dev/null; then
    /usr/bin/time -f "kissat wall %e s maxrss %M KB" nice -n 10 bin/kissat -q --time=5400 $f.cnf $f.drat > $f.res 2>&1
  fi
  echo "== $c $s"; grep -E "^s |wall" $f.res
  if grep -q "^s UNSATISFIABLE" $f.res && [ -s $f.drat ]; then
    ls -la $f.drat | awk '{print "drat bytes", $5}'
    /usr/bin/time -f "drat-trim wall %e s" nice -n 10 bin/drat-trim $f.cnf $f.drat 2>&1 | grep -E "^s |wall|verification time"
  fi
  rm -f $f.drat
done; done
echo WLOG3_DONE
