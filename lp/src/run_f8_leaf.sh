#!/bin/bash
# session VII: F8 + Y3 (degree vector fixed) LP probe at D=5: $1 = config r,u ; $2 = m vector
cd "$(dirname "$0")"
export PC6_SECTIONS=1 PC6_EXACT_K=4 PC6_F6_RMIN=5 PC6_DMH=1 PC6_MFIX=$2
tag=f8leaf_${1/,/_}_${2//,/_}
nice -n 5 .venv/bin/python lp_probe.py --n 46 --D 5 --config $1 --cert --tag $tag >> out/$tag.log 2>&1
echo LEAF_DONE >> out/$tag.log
