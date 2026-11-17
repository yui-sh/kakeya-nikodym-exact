#!/bin/bash
# kill_pat.sh PATTERN : kill processes whose command line matches PATTERN (run from a file to avoid self-match)
pkill -f "$1"; sleep 1; pgrep -af "$1" | grep -v kill_pat || echo "none left matching $1"
