#!/bin/bash
# One-command replay of the release.   usage: bash replay.sh quick | full | nik29 | lp
#   quick : witness checks, symmetry/orbit checks, exact reductions, K(3,3) (18 subcases), tiny Nikodym
#           cases (q <= 5, Nik(3,3)), controls (positive + deliberately broken)  -- about 2 CPU-minutes
#   full  : quick + all certified SAT batches: K(3,4) 39 subcases, K(4,3) 77 instances, Nik(3,4) 43
#           instances, Nik(2,7) 4 instances, Nik(2,8) 2 SAT controls + 162 cubes + 2 cover CNFs (CPU budget: see README)
#   nik29 : Nik(2,9) = 63: q=7 controls, |C|=18 SAT controls, 17 whole cases, 11 x (cover + 100 cubes); ~50 CPU-h CaDiCaL + 4.5 h cake_lpr
#   lp    : the K(3,5) >= 47 LP/Farkas chain (lp/replay_lp.sh all; needs the local .venv and omitted LP inputs, see README)
# Every UNSAT verdict is checked by cake_lpr (formally verified LRAT checker); proofs are deleted after
# checking, their sha256 and sizes are kept in results/.  Everything runs under nice -n 19, one solver at a time.
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)
mkdir -p results tmp big
stage=${1:-quick}
echo "== toolchain"; for b in cadical lrat-trim cake_lpr; do [ -x repos/bin/$b ] || { echo "missing repos/bin/$b: run build_tools.sh"; exit 1; }; done
sha256sum repos/bin/cadical repos/bin/lrat-trim repos/bin/cake_lpr
echo "== self-tests of the fresh library"; (cd src && python3 ff.py && python3 ag.py >/dev/null && python3 card.py)
echo "== witnesses (independent checkers)"; (cd src && python3 check_kakeya.py ../witnesses/kakeya_*.json && python3 check_nikodym.py ../witnesses/nikodym_*.json)
echo "== symmetry checks (Nikodym WLOG steps)"; (cd src && python3 symmetry_check.py | tail -1)
echo "== exact reductions"; (cd src && python3 k43_reduction.py | tail -1 && python3 k35_lower46_check.py | tail -1)
echo "== orbit enumerations (collineation + coverage asserts inside)"; (cd src && python3 enc_kakeya_maxdeg.py --q 4 --ub 27 --list | tail -1 && python3 enc_k43.py --orbits 4 | python3 -c "import json,sys; print('K(4,3) size-4 slice orbits:', json.load(sys.stdin)['n_orbits'])" && python3 enc_nik34.py --list | head -2)
run_batch() { python3 src/run_batch.py "$1" "$2" --timeout "$3" --max-load 8 --summary "$2/summary.json" | tail -1; }
echo "== K(3,3) = 13 (positive control of the max-degree pipeline)"; run_batch specs/k33.json big/k33 600
echo "== K(3,4) ub=28 positive controls (SAT, 28-sets checked)"; run_batch specs/k34ctrl.json big/k34ctrl 900
echo "== K(4,3) capacity table positive controls (k points, f(k) directions: SAT)"; run_batch specs/k43ctrl.json big/k43ctrl 600
echo "== tiny Nikodym cases"; run_batch specs/nik_small.json big/nik_small 600
echo "== controls"; bash controls.sh
if [ "$stage" = full ]; then
  echo "== K(3,4) = 28: 39 subcases";          run_batch specs/k34.json big/k34 3600
  echo "== K(4,3) = 27: 10 capacity + 67 cases"; run_batch specs/k43.json big/k43 3600
  echo "== Nik(3,4) = 22: 43 instances";       run_batch specs/nik34.json big/nik34 1800
  echo "== Nik(2,7) = 36 / Nik(2,8) = 48";     run_batch specs/nik27.json big/nik27 1200; run_batch specs/nik28.json big/nik28 1200; run_batch specs/nik28cubes.json big/nik28cubes 2400
fi
if [ "$stage" = nik29 ]; then
  echo "== Nik(2,9) = 63 (max-secant split; see nik29/REPORT.md)"; (cd src && python3 symmetry_check_nik29.py 7 8 9 | tail -1)
  run_batch specs/nik29ctrl7.json big/nik29/ctrl7 600; run_batch specs/nik29sat.json big/nik29 1200; run_batch specs/nik29.json big/nik29 3600
  run_batch specs/nik29cubes_s4.json big/nik29 1800; run_batch specs/nik29cubes_s5.json big/nik29 1800
  python3 nik29/collect.py specs/nik29ctrl7.json specs/nik29sat.json specs/nik29.json specs/nik29cubes_s4.json specs/nik29cubes_s5.json | tail -1
  nice -n 19 python3 nik29/final_audit.py | tail -3
fi
if [ "$stage" = lp ]; then bash lp/replay_lp.sh all; fi
python3 src/collect_results.py
echo "== manifest"; python3 src/manifest.py --check MANIFEST.sha256 || echo "(manifest mismatch: regenerate with 'python3 src/manifest.py MANIFEST.sha256 src witnesses specs results nik29 lp/certs lp/src')"
echo "REPLAY $stage DONE"
