#!/bin/bash
cd "$(dirname "$0")"
for spec in "5,15 4" "6,17 4" "6,17 5" "6,17 6" "6,17 7"; do set -- $spec
  nice -n 10 .venv/bin/python degvec_bb.py --n 46 --D $2 --config $1 >> out/degvec_bb_Y.log 2>&1
done
echo BB_Y_DONE >> out/degvec_bb_Y.log
