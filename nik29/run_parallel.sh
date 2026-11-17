#!/bin/bash
# Run the cases of a spec in up to N parallel lanes (each lane = one sequential run_batch over its own spec slice).
# usage: run_parallel.sh SPEC OUTDIR LOGDIR N TIMEOUT
SPEC=$1; OUT=$2; LOG=$3; N=$4; TO=$5
cd "$(dirname "$0")/.."
mkdir -p "$OUT" "$LOG"
python3 - "$SPEC" "$LOG" "$N" <<'PY'
import json, sys
spec, log, n = json.load(open(sys.argv[1])), sys.argv[2], int(sys.argv[3])
for i in range(n):
    json.dump(spec[i::n], open(f"{log}/lane{i}.json", "w"))
PY
for i in $(seq 0 $((N-1))); do
  nice -n 19 python3 src/run_batch.py "$LOG/lane$i.json" "$OUT" --timeout "$TO" --max-load 8 > "$LOG/lane$i.log" 2>&1 &
done
wait
echo "ALL LANES DONE" >> "$LOG/lanes.done"
