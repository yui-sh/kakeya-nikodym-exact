# nik34

| case | expected | result | cadical s | lrat-trim s | cake_lpr s | LRAT bytes | note |
|---|---|---|---|---|---|---|---|
| n20_m4383_v0 | UNSAT | UNSAT | 4.4 | 0.21 | 0.65 | 26422113 | cake_lpr VERIFIED |
| n20_m4383_v1 | UNSAT | UNSAT | 8.92 | 0.36 | 1.49 | 44641582 | cake_lpr VERIFIED |
| n20_m4383_d15 | UNSAT | UNSAT | 0.04 | 0.03 | 0.08 | 156975 | cake_lpr VERIFIED |
| n20_m4383_d10305 | UNSAT | UNSAT | 7.95 | 0.67 | 1.04 | 36842822 | cake_lpr VERIFIED |
| n20_m4415_anchor | UNSAT | UNSAT | 13.56 | 0.42 | 1.68 | 41299205 | cake_lpr VERIFIED |
| n20_m4479_anchor | UNSAT | UNSAT | 5.92 | 0.16 | 0.63 | 21426448 | cake_lpr VERIFIED |
| n20_m5023_anchor | UNSAT | UNSAT | 6.77 | 0.3 | 1.06 | 27314170 | cake_lpr VERIFIED |
| n21_m4383_v0 | UNSAT | UNSAT | 20.06 | 0.93 | 2.81 | 113285236 | cake_lpr VERIFIED |
| n21_m4383_v1 | UNSAT | UNSAT | 16.85 | 0.71 | 3.45 | 102996694 | cake_lpr VERIFIED |
| n21_m4383_d15 | UNSAT | UNSAT | 8.54 | 0.34 | 1.02 | 45574913 | cake_lpr VERIFIED |
| n21_m4383_d10305 | UNSAT | UNSAT | 12.83 | 0.53 | 1.5 | 67889916 | cake_lpr VERIFIED |
| n21_m4415_anchor | UNSAT | UNSAT | 13.03 | 0.48 | 1.59 | 65125596 | cake_lpr VERIFIED |
| n21_m4479_anchor | UNSAT | UNSAT | 6.09 | 0.3 | 0.73 | 31561698 | cake_lpr VERIFIED |
| n21_m4607_v0 | UNSAT | UNSAT | 3.59 | 0.16 | 0.5 | 23726807 | cake_lpr VERIFIED |
| n21_m4607_v1 | UNSAT | UNSAT | 3.12 | 0.14 | 0.47 | 20137934 | cake_lpr VERIFIED |
| n21_m4607_v4 | UNSAT | UNSAT | 3.82 | 0.27 | 0.92 | 21532857 | cake_lpr VERIFIED |
| n21_m4607_v5 | UNSAT | UNSAT | 6.43 | 0.42 | 0.97 | 31639642 | cake_lpr VERIFIED |
| n21_m4607_v8 | UNSAT | UNSAT | 6.06 | 0.23 | 0.76 | 31341488 | cake_lpr VERIFIED |
| n21_m4607_d15 | UNSAT | UNSAT | 2.18 | 0.1 | 0.36 | 12868202 | cake_lpr VERIFIED |
| n21_m4607_d4369 | UNSAT | UNSAT | 3.3 | 0.09 | 0.3 | 14319333 | cake_lpr VERIFIED |
| n21_m4607_d10305 | UNSAT | UNSAT | 2.86 | 0.14 | 0.41 | 17613221 | cake_lpr VERIFIED |
| n21_m4959_v0 | UNSAT | UNSAT | 5.89 | 0.21 | 0.63 | 29132352 | cake_lpr VERIFIED |
| n21_m4959_v1 | UNSAT | UNSAT | 4.58 | 0.2 | 0.64 | 28244384 | cake_lpr VERIFIED |
| n21_m4959_v3 | UNSAT | UNSAT | 6.75 | 0.25 | 0.83 | 37412384 | cake_lpr VERIFIED |
| n21_m4959_v6 | UNSAT | UNSAT | 5.51 | 0.25 | 0.9 | 33130374 | cake_lpr VERIFIED |
| n21_m4959_d15 | UNSAT | UNSAT | 2.35 | 0.11 | 0.34 | 15002694 | cake_lpr VERIFIED |
| n21_m4959_d10305 | UNSAT | UNSAT | 3.65 | 0.16 | 0.43 | 15177824 | cake_lpr VERIFIED |
| n21_m4959_d33825 | UNSAT | UNSAT | 4.05 | 0.16 | 0.45 | 20729912 | cake_lpr VERIFIED |
| n21_m4991_anchor | UNSAT | UNSAT | 4.71 | 0.16 | 0.53 | 23120753 | cake_lpr VERIFIED |
| n21_m5023_anchor | UNSAT | UNSAT | 11.11 | 0.43 | 1.31 | 57821121 | cake_lpr VERIFIED |
| n21_m5055_anchor | UNSAT | UNSAT | 3.0 | 0.14 | 0.56 | 18535502 | cake_lpr VERIFIED |
| n21_m5087_anchor | UNSAT | UNSAT | 3.77 | 0.18 | 0.57 | 23717156 | cake_lpr VERIFIED |
| n21_m14687_v0 | UNSAT | UNSAT | 5.41 | 0.24 | 1.5 | 33828031 | cake_lpr VERIFIED |
| n21_m14687_v1 | UNSAT | UNSAT | 4.1 | 0.18 | 0.57 | 25985479 | cake_lpr VERIFIED |
| n21_m14687_v6 | UNSAT | UNSAT | 6.94 | 0.26 | 0.72 | 37918134 | cake_lpr VERIFIED |
| n21_m14687_d15 | UNSAT | UNSAT | 1.92 | 0.08 | 0.29 | 12356763 | cake_lpr VERIFIED |
| n21_m14687_d10305 | UNSAT | UNSAT | 1.44 | 0.08 | 0.23 | 10184471 | cake_lpr VERIFIED |
| n21_m14687_d17025 | UNSAT | UNSAT | 3.43 | 0.11 | 0.35 | 14398604 | cake_lpr VERIFIED |
| control22_cone | SAT | SAT | 2.73 |  |  |  | NIKODYM OK |N|=22<=22 |
| n20_parallel_pair | UNSAT | UNSAT | 0.04 | 0.01 | 0.09 | 93767 | cake_lpr VERIFIED |
| n20_skew_pair | UNSAT | UNSAT | 17.23 | 0.81 | 2.15 | 79628117 | cake_lpr VERIFIED |
| n21_parallel_pair | UNSAT | UNSAT | 12.69 | 0.57 | 1.92 | 73510809 | cake_lpr VERIFIED |
| n21_skew_pair | UNSAT | UNSAT | 13.26 | 0.75 | 2.0 | 72891945 | cake_lpr VERIFIED |
