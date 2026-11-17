# audit/ — fresh-code audit of the K(3,5) >= 47 hand lemmas (see AUDIT.md)

Nothing here imports from lp/src or ../kakeya-fq; their JSON tables/certificates are only read for comparison.

| script | what it reproduces | runtime |
|---|---|---|
| `check_planar.c` → `out/planar_raw.txt` | all 2^25 subsets of AG(2,5): section tuples (k,r,u,c,t,e), full-line profile minima | 31 s |
| `check_planar_tables.py` | sections table, inside profiles, f(5,r), Ft (projection table), F(a), Lemma A, budget caps, b<=4 | 1 s |
| `ag3.py`, `check_identities.py [q ...]` | every identity/lemma on random, greedy and forced-(5,15)-plane witness systems in AG(3,q) (NSYS env) | 2 min |
| `check_projective.py`, `check_projective_tables.py` | PG(2,5) direction subsets j<=7: T, g, rich secants, star bounds; vs their JSON | 3 min |
| `check_hand46.py` | hand certificate K(3,5)>=46 (sympy, knapsack) | 1 s |
| `check_small_n_lp.py` | plain counting LP excludes every n<=44 with exact Farkas certificates (needs scipy; run with ../../kakeya-fq/.venv/bin/python) | 1 min |
| `check_coverage.py` | raw degree-vector universes, ledger tree walks, side rows of the 82 reduction certificates | 10 s |
| `check_geometry.c` → `out/geometry_raw.txt` | all r=5,6 line configurations x extra points with line-incidence masks | 11 s |
| `check_prices.py` (NSAMPLE env) | catalogue vs geometry14_D6.json; exact re-pricing of samples of every m-free price table | 15 min |
| `check_cert_vi.py [r,u,D ...]` | from-scratch session-VI generator (sections F5, exact items + F7, F6 flags R1–R4); exact replay of the 6 session-VI certificates | 1 min (D=4) / 10 min (D=5–7) |
| `check_cert_v.py [cert files]` | from-scratch session-V LP row generator; exact replay of the 76 session-V reduction certificates | 3 min (rmax4) / ~30 min (configs) |
| `dump_research_rownames.py`, `dump_research_rownames14.py` | **the only files here that import lp/src** (bookkeeping only): dump the research model's row list (name, sense, rhs, nnz) and column counts for a positional root certificate; run from `lp/src` with its venv | 15 s – 4 min each |
| `check_root.py <tag>` (515, 617, D6_515, D6_617, D7_515, D7_617_27, D7_617_0, D7_617_dom) | from-scratch generator of the full D = 5/6/7 root models (session-VI model + prefix caps, F8, F9, every lift/point/budget/hull family, own prices from own tables); exact replay of the 8 root certificates, column universe and row supports compared with the research model (AUDIT.md §4.4) | 15 s – 12 min (domain cert: ~2 h) |
