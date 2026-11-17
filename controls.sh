#!/bin/bash
# Controls: (1) SAT models at n = K are decoded and checked by the independent checkers (positive controls,
# done inside collect_results.py for every SAT record with a model); (2) deliberately broken encodings
# must give WRONG verdicts: 'drop-cap' (size cap omitted) is SAT at ub = 12 for K(3,3) where UNSAT is
# expected (its model is a Kakeya set, but larger than 12); 'extra-cap' (cap one too small) is UNSAT at
# ub = 13 where the honest encoding is SAT.  Exit 0 iff both broken encodings misbehave as expected.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p big/controls && cd big/controls
python3 ../../src/enc_kakeya_maxdeg.py --q 3 --ub 12 --D 4 --rep 0 --broken drop-cap --out broken_dropcap.cnf >/dev/null
set +e; python3 ../../src/prove.py broken_dropcap.cnf broken_dropcap.json >/dev/null; rc1=$?; set -e
python3 ../../src/enc_kakeya_maxdeg.py --q 3 --ub 12 --D 4 --rep 0 --model broken_dropcap.cnf.model > broken_dropcap.points.json
sz=$(python3 -c "import json; print(len(json.load(open('broken_dropcap.points.json'))['points']))")
python3 ../../src/enc_kakeya_maxdeg.py --q 3 --ub 13 --D 3 --rep 1 --broken extra-cap --out broken_extracap.cnf >/dev/null
set +e; python3 ../../src/prove.py broken_extracap.cnf broken_extracap.json >/dev/null; rc2=$?; set -e
echo "broken drop-cap: result=$(python3 -c "import json;print(json.load(open('broken_dropcap.json'))['result'])") (expected SAT, honest encoding UNSAT), model size $sz > 12"
echo "broken extra-cap: result=$(python3 -c "import json;print(json.load(open('broken_extracap.json'))['result'])") (expected UNSAT, honest encoding SAT)"
[ $rc1 -eq 10 ] && [ "$sz" -gt 12 ] && [ $rc2 -eq 0 ] && echo "CONTROLS OK: both broken encodings detected" || { echo "CONTROLS FAILED"; exit 1; }
