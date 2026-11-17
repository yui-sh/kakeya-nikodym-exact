#!/bin/bash
# after LNS: replace the queued part-1 (1800 s cap) by: class 4 with the G-invariant+identities model at n=46..52, then part 1 desc with 600 s cap
cd "$(dirname "$0")"
until grep -q LNS_DONE out/core2.log; do sleep 30; done
sleep 10; for p in $(pgrep -f "[p]art 1/2"); do kill $p; done
nice -n 10 .venv/bin/python sym_ident.py --q 5 --gens out/sym_q5_followup.jsonl --idx 4 --n 46 47 48 49 50 51 52 --time 900 > out/sym_ident_4.log 2>&1
echo SYMIDENT_DONE
nice -n 10 .venv/bin/python ident_all.py --q 5 --n 46 --time 600 --part 1/2 --order desc > out/ident_q5_n46_p1desc.log 2>&1
