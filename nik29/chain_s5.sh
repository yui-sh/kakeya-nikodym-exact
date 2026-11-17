#!/bin/bash
# after ALL s=4 cube lanes finish (including the relaunched lane0), run the s=5 cube spec (3 grouped lanes, 1200 s per cube)
cd "$(dirname "$0")/.."
until [ -f nik29/cubes_s4/lanes.done ]; do sleep 60; done
while pgrep -f "cubes_s4/lane0.json" > /dev/null; do sleep 60; done
nik29/run_lanes_staged.sh specs/nik29cubes_s5.json nik29/cubes_s5 1200
