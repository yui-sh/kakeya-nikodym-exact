#!/bin/bash
# Pinned toolchain for the certified pipeline (builds into repos/, git-ignored).  Network required.
#   cadical   https://github.com/arminbiere/cadical    commit c60730422e758ef1cebe7aeddf2dda31c996bf04 (3.0.1)
#   lrat-trim https://github.com/arminbiere/lrat-trim  commit b30f400f4ee5c32b77ee566a7c006081b521534f (0.2.0)
#   cake_lpr  https://github.com/tanyongkiam/cake_lpr  commit 033aad53d1733bd4b8ca41cf6f1b5fb5d84f8de2
#             (verified CakeML binary shipped as cake_lpr.S / cake_lpr_arm8.S + basis_ffi.c, compiled with gcc)
# On this box (aarch64, gcc 13.3.0) the resulting binaries have the sha256 listed in TOOLCHAIN.md.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p repos/bin
if [ ! -x repos/bin/cadical ]; then
  git clone -q https://github.com/arminbiere/cadical repos/cadical && (cd repos/cadical && git checkout -q c60730422e758ef1cebe7aeddf2dda31c996bf04 && ./configure -q && nice -n 19 make -j2 >/dev/null && cp build/cadical ../bin/)
fi
if [ ! -x repos/bin/lrat-trim ]; then
  git clone -q https://github.com/arminbiere/lrat-trim repos/lrat-trim && (cd repos/lrat-trim && git checkout -q b30f400f4ee5c32b77ee566a7c006081b521534f && ./configure && nice -n 19 make >/dev/null && cp lrat-trim ../bin/)
fi
if [ ! -x repos/bin/cake_lpr ]; then
  C=033aad53d1733bd4b8ca41cf6f1b5fb5d84f8de2; mkdir -p repos/cake_lpr-src; cd repos/cake_lpr-src
  for f in cake_lpr.S cake_lpr_arm8.S basis_ffi.c Makefile cake_lpr.sha256 LICENSE README.md example.cnf example.lpr; do curl -sL -o $f https://raw.githubusercontent.com/tanyongkiam/cake_lpr/$C/$f; done
  sha256sum -c cake_lpr.sha256 --ignore-missing
  nice -n 19 make cake_lpr && cp cake_lpr ../bin/ && ./cake_lpr example.cnf example.lpr; cd ../..
fi
sha256sum repos/bin/*
