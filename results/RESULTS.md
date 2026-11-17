# Batch results (measured on this box, nice -n 19, one solver at a time)

| batch | cases | UNSAT verified by cake_lpr | SAT (model re-checked) | missing | cadical s | lrat-trim s | cake_lpr s | largest LRAT (bytes) |
|---|---|---|---|---|---|---|---|---|
| k33 | 18 | 17 | 1 | 0 | 0 | 0 | 1 | 253518 |
| k34 | 39 | 39 | 0 | 0 | 1171 | 49 | 175 | 2217264791 |
| k34ctrl | 2 | 0 | 2 | 0 | 1 | 0 | 0 | 0 |
| k43 | 77 | 77 | 0 | 0 | 562 | 22 | 66 | 313051636 |
| k43ctrl | 10 | 0 | 10 | 0 | 0 | 0 | 0 | 0 |
| nik27 | 4 | 2 | 2 | 0 | 45 | 3 | 8 | 118084706 |
| nik28 | 2 | 0 | 2 | 0 | 1 | 0 | 0 | 0 |
| nik28cubes | 164 | 164 | 0 | 0 | 11916 | 412 | 1516 | 1018860434 |
| nik34 | 43 | 42 | 1 | 0 | 281 | 12 | 39 | 113285236 |
| nik_small | 12 | 6 | 6 | 0 | 4 | 0 | 1 | 21686342 |

problems (unexpected verdicts, unverified UNSAT, failed model checks): 0

Nik(2,9) = 63 (specs/nik29*.json, 1128 UNSAT + 4 SAT controls + 32 q = 7 controls) is collected separately by `nik29/collect.py` -> `nik29/cases.md`, `nik29/summary.json`; audited by `nik29/final_audit.py` -> `nik29/FINAL_AUDIT.md`.
