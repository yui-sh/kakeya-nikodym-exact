#!/bin/bash
# session IV: the four unresolved symmetry classes at ub=52 (WRITEUP2 §2), kissat + drat-trim, 1 core
cd "$(dirname "$0")"
for i in 0 1 2 4; do
  f=out/symc_q5_ub52_$i
  /usr/bin/time -f "kissat wall %e s maxrss %M KB" nice -n 10 bin/kissat -q --time=3600 $f.cnf $f.drat > $f.res 2>&1
  echo "== class $i: $(grep -E '^s ' $f.res) $(grep wall $f.res)"
  if grep -q "^s UNSATISFIABLE" $f.res; then
    /usr/bin/time -f "drat-trim wall %e s" nice -n 10 bin/drat-trim $f.cnf $f.drat 2>&1 | grep -E "^s |wall"
  fi
  rm -f $f.drat
done
echo SYMC_DONE
