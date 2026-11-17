# Pinned toolchain (measured on this box: Linux 7.0 aarch64, 4 cores, gcc 13.3.0, Python 3.12.3)

| tool | role | source / commit | binary sha256 (this box) |
|---|---|---|---|
| CaDiCaL 3.0.1 | SAT solver, emits LRAT proofs (`--lrat --no-binary`) | github.com/arminbiere/cadical @ c60730422e758ef1cebe7aeddf2dda31c996bf04 | fafdc017e81c48e0418986460cc7738b28c5efb823764c6a7b875c1a321756e2 |
| lrat-trim 0.2.0 | trims LRAT proofs (optional; cake_lpr also accepts the untrimmed proof) | github.com/arminbiere/lrat-trim @ b30f400f4ee5c32b77ee566a7c006081b521534f | 1c4a4ec77c953dabc6c8b47835c15d94d7bc48fff75fbbcde8b4ff6fb5626129 |
| **cake_lpr** | **formally verified LRAT checker** (CakeML; Tan, Heule, Myreen). Its correctness theorem: if it prints `s VERIFIED UNSAT`, the input CNF is unsatisfiable. | github.com/tanyongkiam/cake_lpr @ 033aad53d1733bd4b8ca41cf6f1b5fb5d84f8de2, built from `cake_lpr_arm8.S` (sha256 50b63e18…4852b, matches the repo's `cake_lpr.sha256`) + `basis_ffi.c` with gcc -O2 | 130816693a5392bbba02cf0c382ac3569542dcbb482ad0b409e8a1136ab4126d |
| Python 3.12.3 | encoders, checkers, drivers (standard library only) | Ubuntu 24.04 | — |
| kissat 4.0.4 / drat-trim | NOT used by the release pipeline (the research sessions used them); historical local copies in ignored repos/bin for comparison only; not distributed | — | 354c1891… / 683fca42… |

Trust base of an UNSAT verdict in this release: the CNF file (produced by the fresh encoder in `src/`, whose
meaning is argued in README.md) and cake_lpr. CaDiCaL and lrat-trim are untrusted proof producers.
The LRAT proofs are deleted after checking (the largest was 3.9 GB, Nik(2,9) whole case s5o0h; 1.1 GB outside Nik(2,9)); their sha256 and byte size are
recorded in `results/*/<case>.json` and `nik29/records/*/<case>.json`, and `replay.sh` regenerates and re-checks them.

`build_tools.sh` rebuilds all three from the pinned commits.  `repos/` is git-ignored.

No third-party tool source or binaries are included in the tracked release tree. `build_tools.sh`
downloads upstream sources into ignored `repos/`; their upstream licences apply independently.
Alternatively, create `repos/bin/` and copy or symlink existing executables named `cadical`,
`lrat-trim`, and `cake_lpr` there. These are the documented external tool locations; the replay
has no dependency on a private quota script. LP replay additionally uses the root `.venv/bin/python`
(with numpy, scipy, sympy) or an absolute `PY` override; standard shell utilities and `/usr/bin/time`
are required. The historical binary hashes above are machine-specific.
