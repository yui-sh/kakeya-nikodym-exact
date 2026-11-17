# FINAL AUDIT of Nik(2,9) = 63 (max-secant split, enc_nik29.py) — 2026-11-13 00:27 UTC

Produced by `nik29/final_audit.py` (this file is its verbatim output).  The audit is independent of the lane scripts: it enumerates the cases with its own orbit code, parses the regenerated DIMACS files and checks the cube/cover structure clause by clause, and applies its own acceptance predicate to every record.

## (a) the 28 cases and how each is closed

own orbit enumeration (union-find over all images, no encoder code): representatives per s = s=2:1, s=3:2, s=4:3, s=5:3, s=6:2, s=7:1, s=8:1, s=9:1; 14 representatives x 2 tangent flags = 28 cases; identical to the encoder's list (same representatives, same orbit sizes; sizes sum to C(9,s) for every s).

| case | s | axis set S | tangent | closed by | superseded whole-case entry |
|---|---|---|---|---|---|
| q9K19_s2o0h | 2 | [0, 1] | h | whole-case record |  |
| q9K19_s2o0d | 2 | [0, 1] | d | whole-case record |  |
| q9K19_s3o0h | 3 | [0, 1, 2] | h | whole-case record |  |
| q9K19_s3o0d | 3 | [0, 1, 2] | d | whole-case record |  |
| q9K19_s3o1h | 3 | [0, 1, 3] | h | whole-case record |  |
| q9K19_s3o1d | 3 | [0, 1, 3] | d | whole-case record |  |
| q9K19_s4o0h | 4 | [0, 1, 2, 3] | h | cover + 100 cubes | whole case never run |
| q9K19_s4o0d | 4 | [0, 1, 2, 3] | d | cover + 100 cubes | whole case never run |
| q9K19_s4o1h | 4 | [0, 1, 3, 4] | h | cover + 100 cubes | whole case never run |
| q9K19_s4o1d | 4 | [0, 1, 3, 4] | d | cover + 100 cubes | whole case never run |
| q9K19_s4o2h | 4 | [0, 1, 3, 5] | h | cover + 100 cubes | whole case never run |
| q9K19_s4o2d | 4 | [0, 1, 3, 5] | d | cover + 100 cubes | whole case never run |
| q9K19_s5o0h | 5 | [0, 1, 2, 3, 4] | h | whole-case record |  |
| q9K19_s5o0d | 5 | [0, 1, 2, 3, 4] | d | cover + 100 cubes | whole-case attempt: UNKNOWN after 1800.11 s |
| q9K19_s5o1h | 5 | [0, 1, 2, 3, 6] | h | cover + 100 cubes | whole case never run |
| q9K19_s5o1d | 5 | [0, 1, 2, 3, 6] | d | cover + 100 cubes | whole case never run |
| q9K19_s5o2h | 5 | [0, 1, 2, 3, 7] | h | cover + 100 cubes | whole case never run |
| q9K19_s5o2d | 5 | [0, 1, 2, 3, 7] | d | cover + 100 cubes | whole case never run |
| q9K19_s6o0h | 6 | [0, 1, 2, 3, 4, 5] | h | whole-case record |  |
| q9K19_s6o0d | 6 | [0, 1, 2, 3, 4, 5] | d | whole-case record |  |
| q9K19_s6o1h | 6 | [0, 1, 2, 3, 4, 6] | h | whole-case record |  |
| q9K19_s6o1d | 6 | [0, 1, 2, 3, 4, 6] | d | whole-case record |  |
| q9K19_s7o0h | 7 | [0, 1, 2, 3, 4, 5, 6] | h | whole-case record |  |
| q9K19_s7o0d | 7 | [0, 1, 2, 3, 4, 5, 6] | d | whole-case record |  |
| q9K19_s8o0h | 8 | [0, 1, 2, 3, 4, 5, 6, 7] | h | whole-case record |  |
| q9K19_s8o0d | 8 | [0, 1, 2, 3, 4, 5, 6, 7] | d | whole-case record |  |
| q9K19_s9o0h | 9 | [0, 1, 2, 3, 4, 5, 6, 7, 8] | h | whole-case record |  |
| q9K19_s9o0d | 9 | [0, 1, 2, 3, 4, 5, 6, 7, 8] | d | whole-case record |  |

17 cases closed by a verified whole-case record, 11 by a verified cover + 100 cubes, total 28 / 28.

## (b) structure of covers and cubes, (c) byte-identical regeneration

q9K19_s4o0h: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s4o0d: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s4o1h: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [27, 36] = [(3, 0), (4, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s4o1d: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [27, 36] = [(3, 0), (4, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s4o2h: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [27, 45] = [(3, 0), (5, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s4o2d: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [27, 45] = [(3, 0), (5, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s5o0d: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s5o1h: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s5o1d: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s5o2h: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).
q9K19_s5o2d: cover = 81 witness clauses (all clauses of the whole CNF) + 100 negated cubes; 100 cubes = whole CNF + 2 witness units of points [18, 27] = [(2, 0), (3, 0)], full 10 x 10 product; 101/101 regenerated CNFs byte-identical to the records (sha256).

Regenerated and hashed 1165 CNFs (28 whole cases incl. the superseded ones that have a record, 11 covers, 1100 cubes, 4 SAT controls, 32 q=7 controls): 1165 byte-identical to the recorded cnf_sha256, 0 mismatches.

Hardest 10 cubes by CaDiCaL time (all among the byte-identical regenerations):

| cube | CaDiCaL s | LRAT bytes |
|---|---|---|
| q9K19_s4o0d_cube86 | 667.12 | 1798629627 |
| q9K19_s4o0d_cube84 | 664.62 | 1755146277 |
| q9K19_s4o1d_cube84 | 634.99 | 1803481094 |
| q9K19_s4o2d_cube68 | 629.45 | 1890669432 |
| q9K19_s4o2d_cube89 | 622.77 | 1799292683 |
| q9K19_s4o2d_cube48 | 622.72 | 1940921419 |
| q9K19_s4o0d_cube83 | 622.08 | 1671050922 |
| q9K19_s4o1d_cube83 | 618.24 | 1627709476 |
| q9K19_s4o1d_cube59 | 607.98 | 1662489304 |
| q9K19_s4o1d_cube68 | 604.63 | 1590161890 |

## (d) symmetry checks and SAT-control models

```
q=7: (1) every one of the 56 lines is the image of the axis under an explicit affine map
q=7: (2) 42 axis-preserving semilinear collineations; orbit representatives cover all s-subsets of the axis disjointly, s=2..7: s=2:1, s=3:2, s=4:2, s=5:1, s=6:1, s=7:1
q=7: (3) every non-axis line through (0,0) is mapped to x = 0 by a shear fixing the axis pointwise
q=7: (4) every non-axis line through (1,0) is mapped to x = 1 or x = 1 + y by a y-scaling fixing the axis pointwise and x = 0 setwise
q=8: (1) every one of the 72 lines is the image of the axis under an explicit affine map
q=8: (2) 168 axis-preserving semilinear collineations; orbit representatives cover all s-subsets of the axis disjointly, s=2..8: s=2:1, s=3:1, s=4:2, s=5:1, s=6:1, s=7:1, s=8:1
q=8: (3) every non-axis line through (0,0) is mapped to x = 0 by a shear fixing the axis pointwise
q=8: (4) every non-axis line through (1,0) is mapped to x = 1 or x = 1 + y by a y-scaling fixing the axis pointwise and x = 0 setwise
q=9: (1) every one of the 90 lines is the image of the axis under an explicit affine map
q=9: (2) 144 axis-preserving semilinear collineations; orbit representatives cover all s-subsets of the axis disjointly, s=2..9: s=2:1, s=3:2, s=4:3, s=5:3, s=6:2, s=7:1, s=8:1, s=9:1
q=9: (3) every non-axis line through (0,0) is mapped to x = 0 by a shear fixing the axis pointwise
q=9: (4) every non-axis line through (1,0) is mapped to x = 1 or x = 1 + y by a y-scaling fixing the axis pointwise and x = 0 setwise
NIK29 SYMMETRY CHECKS PASSED
```
```
q9K18_s6o0h |C| = 18 big/nik29/q9K18_s6o0h.decoded.json: AG(2,9) |N|=63 |C|=18 -> NIKODYM OK rc 0
q9K18_s6o0d |C| = 18 big/nik29/q9K18_s6o0d.decoded.json: AG(2,9) |N|=63 |C|=18 -> NIKODYM OK rc 0
SAT cases checked: 2
```
```
q7K13_s4o1h |C| = 13 big/nik29/ctrl7/q7K13_s4o1h.decoded.json: AG(2,7) |N|=36 |C|=13 -> NIKODYM OK rc 0
q7K13_s4o1d |C| = 13 big/nik29/ctrl7/q7K13_s4o1d.decoded.json: AG(2,7) |N|=36 |C|=13 -> NIKODYM OK rc 0
SAT cases checked: 2
```
q9K18_s6o0h: decoded model has |C| = 18, axis section = S = [0, 1, 2, 3, 4, 5], max secant = 6, Nikodym-checked (|N| = 63).
q9K18_s6o0d: decoded model has |C| = 18, axis section = S = [0, 1, 2, 3, 4, 5], max secant = 6, Nikodym-checked (|N| = 63).

## (e) record audit

1160 UNSAT records accepted (each: result UNSAT, verified true, cake_verdict 's VERIFIED UNSAT', LRAT sha256 + byte size, CNF sha256, pinned tool hashes); 4 SAT records (models re-checked in (d)); 1164/1164 tracked copies in nik29/records byte-identical to big/nik29; repos/bin binaries match TOOLCHAIN.md.  No instance is counted on solver status alone: the acceptance predicate requires the cake_lpr verdict string.

| spec | instances | UNSAT verified | SAT | CaDiCaL s | lrat-trim s | cake_lpr s | max CaDiCaL s | max LRAT | trim fallbacks | record dates (UTC) |
|---|---|---|---|---|---|---|---|---|---|---|
| nik29ctrl7 | 32 | 30 | 2 | 5 | 0 | 2 | 1 | 0.00 GB | 0 | 2026-11-07 09:32 - 2026-11-07 09:33 |
| nik29sat | 4 | 2 | 2 | 128 | 5 | 23 | 70 | 0.38 GB | 0 | 2026-11-07 09:51 - 2026-11-07 09:54 |
| nik29 | 17 | 17 | 0 | 1955 | 46 | 167 | 1460 | 3.90 GB | 0 | 2026-11-07 09:33 - 2026-11-07 10:05 |
| nik29cubes_s4 | 606 | 606 | 0 | 155898 | 3982 | 13619 | 667 | 1.94 GB | 0 | 2026-11-07 09:59 - 2026-11-08 02:53 |
| nik29cubes_s5 | 505 | 505 | 0 | 20111 | 761 | 2544 | 113 | 0.42 GB | 0 | 2026-11-08 02:54 - 2026-11-08 05:50 |

**Proof of |C| = 19 impossible (the 17 whole cases + 1111 cube instances = 1128 UNSAT instances): CaDiCaL 177964 s (49.4 h), lrat-trim 4790 s, cake_lpr 16330 s (4.5 h); max CaDiCaL 1460 s; largest LRAT 3896123395 bytes (3.90 GB).**

## Verdict

**ALL CHECKS PASSED.**  Every one of the 28 cases is closed by cake_lpr-verified instances whose CNFs are regenerated byte-identically by the shipped encoder; the cube splits are exhaustive by construction (full product) and by the verified cover CNFs; the SAT controls behave as required.  Nik(2,9) = 63 is READY at the Nik(2,8) standard.

(audit wall time 70 s)
