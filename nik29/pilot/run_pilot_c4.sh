#!/bin/bash
cd "$(dirname "$0")/../.."
until [ -f nik29/k18/lanes.done ]; do sleep 5; done
for n in "$@"; do
  ( s=$(date +%s); timeout 900 nice -n 19 repos/bin/cadical -q --lrat --no-binary big/nik29/pilot_c4/cube$n.cnf big/nik29/pilot_c4/cube$n.lrat > /dev/null 2>&1; rc=$?; e=$(date +%s); echo "c4 cube$n rc=$rc $((e-s))s $(stat -c %s big/nik29/pilot_c4/cube$n.lrat 2>/dev/null)" >> nik29/pilot/pilot_ml.log; rm -f big/nik29/pilot_c4/cube$n.lrat ) &
done
wait
