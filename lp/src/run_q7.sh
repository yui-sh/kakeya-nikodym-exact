#!/bin/bash
# Session V: q=7 with the merged LP without degree histograms (classes + directions + bucketed pencils).
cd "$(dirname "$0")"
nice -n 10 .venv/bin/python union_tables.py 7 > out/union_tables_q7.log 2>&1
nice -n 10 .venv/bin/python planecount6.py --q 7 --n 106 107 108 --no-hist --cert --tag q7nohist > out/planecount6_q7.log 2>&1
echo Q7_DONE >> out/planecount6_q7.log
