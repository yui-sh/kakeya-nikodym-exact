#!/bin/bash
# Session V: (re)launch the two configuration-probe chains; safe to call from the tool shell (patterns live here, not in the caller's command line).
cd "$(dirname "$0")"
pkill -f "run_cfg_chain"; pkill -f "lp_probe.py --n 46 --config"; sleep 2
[ -s out/lp_probe_cfgA.log ] && mv out/lp_probe_cfgA.log out/lp_probe_cfgA_old2.log
[ -s out/lp_probe_cfgB.log ] && mv out/lp_probe_cfgB.log out/lp_probe_cfgB_old2.log
printf '%s\n' '#!/bin/bash' 'cd "$(dirname "$0")"' 'for cfg in 5,21 5,18 5,17 5,16 5,15; do' '  nice -n 10 .venv/bin/python lp_probe.py --n 46 --config $cfg --cert --tag cfg_${cfg/,/_} >> out/lp_probe_cfgA.log 2>&1' 'done' 'echo CHAIN_A_DONE >> out/lp_probe_cfgA.log' > run_cfg_chainA.sh
printf '%s\n' '#!/bin/bash' 'cd "$(dirname "$0")"' 'for cfg in 6,25 6,21 6,19 6,18 6,17; do' '  nice -n 10 .venv/bin/python lp_probe.py --n 46 --config $cfg --cert --tag cfg_${cfg/,/_} >> out/lp_probe_cfgB.log 2>&1' 'done' 'echo CHAIN_B_DONE >> out/lp_probe_cfgB.log' > run_cfg_chainB.sh
chmod +x run_cfg_chainA.sh run_cfg_chainB.sh; touch out/lp_probe_cfgA.log out/lp_probe_cfgB.log
rm -f out/lp_probe_cfg_5_21.json out/lp_probe_cfg_6_25.json
nohup ./run_cfg_chainA.sh >/dev/null 2>&1 &
nohup ./run_cfg_chainB.sh >/dev/null 2>&1 &
sleep 3; ps -eo pid,etime,args | grep -E "lp_probe|run_cfg_chain" | grep -v grep | sed 's/\.venv\/bin\/python //'
