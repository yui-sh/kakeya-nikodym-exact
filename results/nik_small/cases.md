# nik_small

| case | expected | result | cadical s | lrat-trim s | cake_lpr s | LRAT bytes | note |
|---|---|---|---|---|---|---|---|
| nik_n2_q2_K3_none | UNSAT | UNSAT | 0.01 | 0.0 | 0.0 | 309 | cake_lpr VERIFIED |
| nik_n2_q2_K2_none | SAT | SAT | 0.01 |  |  |  | NIKODYM OK |C|=2>=2 |
| nik_n2_q3_K5_none | UNSAT | UNSAT | 0.01 | 0.0 | 0.01 | 4671 | cake_lpr VERIFIED |
| nik_n2_q3_K4_none | SAT | SAT | 0.0 |  |  |  | NIKODYM OK |C|=4>=4 |
| nik_n2_q4_K7_none | UNSAT | UNSAT | 0.02 | 0.01 | 0.03 | 120190 | cake_lpr VERIFIED |
| nik_n2_q4_K6_none | SAT | SAT | 0.04 |  |  |  | NIKODYM OK |C|=6>=6 |
| nik_n2_q5_K9_none | UNSAT | UNSAT | 0.73 | 0.05 | 0.15 | 5425022 | cake_lpr VERIFIED |
| nik_n2_q5_K8_none | SAT | SAT | 0.01 |  |  |  | NIKODYM OK |C|=8>=8 |
| nik_n3_q2_K7_none | UNSAT | UNSAT | 0.0 | 0.0 | 0.0 | 1391 | cake_lpr VERIFIED |
| nik_n3_q2_K6_none | SAT | SAT | 0.04 |  |  |  | NIKODYM OK |C|=6>=6 |
| nik_n3_q3_K19_none | UNSAT | UNSAT | 3.2 | 0.19 | 0.61 | 21686342 | cake_lpr VERIFIED |
| nik_n3_q3_K18_none | SAT | SAT | 0.01 |  |  |  | NIKODYM OK |C|=18>=18 |
