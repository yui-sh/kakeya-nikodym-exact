#!/bin/bash
# lane0 starts now; lanes 1,2 start when nik29/k19/lanes.done exists (whole-case lanes finished).
# usage: run_lanes_staged.sh SPEC LOGDIR TIMEOUT
SPEC=$1; LOG=$2; TO=$3
cd "$(dirname "$0")/.."
mkdir -p "$LOG" big/nik29/cubes
python3 - "$SPEC" "$LOG" <<'PY'
import json, sys
spec, log = json.load(open(sys.argv[1])), sys.argv[2]
# keep each case's cover + cubes together in one lane (the cover entry's generator writes the cube files)
groups = {}
for c in spec:
    groups.setdefault(c["name"].rsplit("_", 1)[0], []).append(c)
keys = sorted(groups)
for i in range(3):
    json.dump([c for k in keys[i::3] for c in groups[k]], open(f"{log}/lane{i}.json", "w"))
PY
nice -n 19 python3 src/run_batch.py "$LOG/lane0.json" big/nik29 --timeout "$TO" --max-load 8 > "$LOG/lane0.log" 2>&1 &
until [ -f nik29/k19/lanes.done ]; do sleep 20; done
for i in 1 2; do
  nice -n 19 python3 src/run_batch.py "$LOG/lane$i.json" big/nik29 --timeout "$TO" --max-load 8 > "$LOG/lane$i.log" 2>&1 &
done
wait
echo "ALL LANES DONE" >> "$LOG/lanes.done"
