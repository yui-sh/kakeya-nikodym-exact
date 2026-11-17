#!/bin/bash
pkill -f run_bb_X.sh; pkill -f "degvec_bb.py --n 46 --D 6 --config 5,16"; pkill -f run_milp_sweep.sh
