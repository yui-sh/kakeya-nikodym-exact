#!/bin/bash
# session IV continuation: verify class 2 DRAT (kissat already UNSAT in out/symc_q5_ub52_2.res), run class 4, then start campaign part 1/2
cd "$(dirname "$0")"
f=out/symc_q5_ub52_2
echo "== class 2: $(grep -E '^s ' $f.res) $(grep wall $f.res)"
ls -la $f.drat | awk '{print "drat bytes", $5}'
/usr/bin/time -f "drat-trim wall %e s" nice -n 10 bin/drat-trim $f.cnf $f.drat 2>&1 | tr '\r' '\n' | grep -E "VERIFIED|wall"
rm -f $f.drat
f=out/symc_q5_ub52_4
/usr/bin/time -f "kissat wall %e s maxrss %M KB" nice -n 10 bin/kissat -q --time=3600 $f.cnf $f.drat > $f.res 2>&1
echo "== class 4: $(grep -E '^s ' $f.res) $(grep wall $f.res)"
if grep -q "^s UNSATISFIABLE" $f.res; then
  ls -la $f.drat | awk '{print "drat bytes", $5}'
  /usr/bin/time -f "drat-trim wall %e s" nice -n 10 bin/drat-trim $f.cnf $f.drat 2>&1 | tr '\r' '\n' | grep -E "VERIFIED|wall"
fi
rm -f $f.drat
echo SYMC_DONE
nice -n 10 .venv/bin/python ident_all.py --q 5 --n 46 --time 1800 --part 1/2 > out/ident_q5_n46_p1.log 2>&1
