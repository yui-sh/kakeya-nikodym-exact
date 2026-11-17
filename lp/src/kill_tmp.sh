#!/bin/bash
for p in $(pgrep -f "sixth-dir""ection pencil"); do [ "$p" = "$$" ] || [ "$p" = "$PPID" ] || kill "$p" 2>/dev/null; done; echo done
