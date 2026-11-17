#!/bin/bash
# kill stale pc6_validate runs; exclude this script, its parent shell and the parent's parent
for p in $(pgrep -f "pc6_valid""ate.py"); do
  [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null
done
echo killed-validate
