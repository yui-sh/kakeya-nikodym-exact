#!/bin/bash
for p in $(pgrep -f "run_validate_af""ter.sh"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo waiter-killed
