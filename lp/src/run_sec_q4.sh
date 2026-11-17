#!/bin/bash
# Session VI: q=4 calibration of the F5 (sections) family
cd "$(dirname "$0")"
PC6_SECTIONS=1 nice -n 5 .venv/bin/python planecount6.py --q 4 --n 24 25 26 27 28 --tag sec > out/sec_q4.log 2>&1
echo SEC_Q4_DONE >> out/sec_q4.log
