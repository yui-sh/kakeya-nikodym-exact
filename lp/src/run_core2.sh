#!/bin/bash
# core-2 queue after the symmetric classes: replace the auto-started campaign part 1 (asc) by: LNS at 52 (3 seeds, ~1 h), then part 1 in descending-D order
cd "$(dirname "$0")"
until grep -q SYMC_DONE out/symc_q5_ub52.log; do sleep 30; done
sleep 5; for p in $(pgrep -f "[p]art 1/2"); do kill $p; done
for s in sym_n3_q5_53_ord4_1 best_n3_q5_53 bc_n3_q5_53; do
  nice -n 10 .venv/bin/python lns.py --q 5 --seed solutions/$s.json --rounds 40 --time 30 --workers 1 --keep 0.5 --rseed 7 --out out/lns52_$s.json > out/lns52_$s.log 2>&1
  echo "LNS $s: $(grep -c IMPROVED out/lns52_$s.log) improvements; last: $(tail -1 out/lns52_$s.log)"
done
echo LNS_DONE
nice -n 10 .venv/bin/python ident_all.py --q 5 --n 46 --time 1800 --part 1/2 --order desc > out/ident_q5_n46_p1.log 2>&1
