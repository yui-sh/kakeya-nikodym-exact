#!/bin/bash
cd "$(dirname "$0")"
for spec in "5,16 5" "5,16 6" "5,16 7" "5,15 5" "5,15 6" "5,15 7"; do set -- $spec
  nice -n 10 .venv/bin/python degvec_bb.py --n 46 --D $2 --config $1 >> out/degvec_bb_X.log 2>&1
done
echo BB_X_DONE >> out/degvec_bb_X.log
