#!/bin/bash
# usage: wait_task.sh seconds file pattern  -- wait until file contains pattern (or timeout), then print it
t=$1; f=$2; p=$3
for ((i=0;i<t;i+=10)); do grep -q "$p" "$f" 2>/dev/null && break; sleep 10; done
cat "$f"
