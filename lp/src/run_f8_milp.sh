#!/bin/bash
# session VII: whole-configuration MILP with F8 at D=5: $1 config ; $2 m vector or '-' ; $3 time limit
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1
[ "$2" != "-" ] && export PC6_MFIX=$2
tag=milpf8_${1/,/_}_D5_${2//,/_}
nice -n 5 .venv/bin/python milp_probe.py 5 $1 $2 $3 out/${tag}_solution.json > out/$tag.log 2>&1
echo MILP_DONE >> out/$tag.log
