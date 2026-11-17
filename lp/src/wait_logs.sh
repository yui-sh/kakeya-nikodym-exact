#!/bin/bash
# usage: wait_logs.sh seconds logfile...   (polls until all logs contain DONE or timeout)
t=$1; shift
for ((i=0;i<t;i+=15)); do ok=1; for f in "$@"; do grep -q DONE "$f" 2>/dev/null || ok=0; done; [ $ok = 1 ] && break; sleep 15; done
for f in "$@"; do echo "== $f"; tail -n 12 "$f"; done
