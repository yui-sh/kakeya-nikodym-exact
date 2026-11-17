#!/bin/bash
# Session V: final cleanup - stop every campaign process and log tail started this session.
for pat in run_q7.sh union_tables.py "planecount6.py --q 7" rmax_model.py milp_probe.py milp_sweep.py degvec_bb.py "lp_probe.py" bb_int.py run_cfg run_bb_ run_sweep run_certcheck run_cpsat_mvec "tail -q -n 0 -F out/" "tail -n 0 -F /DATA" "tail -n +1 -F /DATA"; do pkill -f "$pat"; done
sleep 1; ps -eo pid,etime,args | grep -E "kakeya-fq|rmax_model|planecount|degvec|milp|lp_probe|union_tables" | grep -v -E "grep|cleanup" | sed "s/.*python //" | cut -c1-100
echo CLEANUP_DONE
