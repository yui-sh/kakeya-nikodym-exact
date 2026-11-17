#!/bin/bash
cd "$(dirname "$0")/../.."
for n in $(seq 0 99); do
  s=$(date +%s); timeout 300 nice -n 19 repos/bin/cadical -q --lrat --no-binary big/nik29/pilot_c4/cube$n.cnf big/nik29/pilot_c4/cube$n.lrat > /dev/null 2>&1; rc=$?; e=$(date +%s)
  echo "c4 cube$n rc=$rc $((e-s))s $(stat -c %s big/nik29/pilot_c4/cube$n.lrat 2>/dev/null)" >> nik29/pilot/pilot_c4.log; rm -f big/nik29/pilot_c4/cube$n.lrat
done
echo DONE >> nik29/pilot/pilot_c4.log
