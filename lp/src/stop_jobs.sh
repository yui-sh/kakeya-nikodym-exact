#!/bin/bash
pkill -f run_rmax56.sh; pkill -f "rmax_model.py --q 5 --n 46"; pkill -f run_cfg_probes.sh
sleep 1; ps -o pid,args -C python,bash | grep -E "rmax|cfg_probes|lp_probe" | grep -v grep
