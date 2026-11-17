#!/bin/bash
# Session V: MILP sweep over all LP-surviving degree vectors; starts when B&B chain Y is done.
cd "$(dirname "$0")"
until grep -q BB_Y_DONE out/degvec_bb_Y.log; do sleep 30; done
nice -n 10 .venv/bin/python milp_sweep.py --time 900 >> out/milp_sweep.log 2>&1
