#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
bash run14.sh replay14_D5.py --fresh-prices > out/replay14_D5_fresh_final.log 2>&1
bash run14.sh coverage14.py --D 6 --replay --fresh-prices > out/replay14_D6_fresh_final.log 2>&1
