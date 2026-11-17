#!/bin/bash
# pilot timing: cadical --lrat on selected cubes, 3 in parallel, timeout 900 s
cd "$(dirname "$0")/../.."
BIN=repos/bin
for i in "$@"; do
  ( s=$(date +%s); timeout 900 nice -n 19 $BIN/cadical -q --lrat --no-binary big/nik29/pilot_h/cube$i.cnf big/nik29/pilot_h/cube$i.lrat > big/nik29/pilot_h/cube$i.out 2>&1; rc=$?; e=$(date +%s); echo "cube$i rc=$rc $((e-s))s $(stat -c %s big/nik29/pilot_h/cube$i.lrat 2>/dev/null)" >> nik29/pilot/pilot_h.log ) &
done
wait
