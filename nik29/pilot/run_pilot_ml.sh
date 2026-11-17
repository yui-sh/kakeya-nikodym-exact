#!/bin/bash
cd "$(dirname "$0")/../.."
for n in "$@"; do
  ( s=$(date +%s); timeout 1200 nice -n 19 repos/bin/cadical -q --lrat --no-binary big/nik29/pilot_ml/$n.cnf big/nik29/pilot_ml/$n.lrat > big/nik29/pilot_ml/$n.out 2>&1; rc=$?; e=$(date +%s); echo "$n rc=$rc $((e-s))s $(stat -c %s big/nik29/pilot_ml/$n.lrat 2>/dev/null)" >> nik29/pilot/pilot_ml.log ) &
done
wait
