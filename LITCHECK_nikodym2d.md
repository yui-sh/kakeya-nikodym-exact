# Literature check: planar Nikodym values q ≤ 9 (STATUS §C item 3) — Yui, 2026-11-07 01:30 ET

Partial, by hand (3 web searches + Tao arXiv:2511.07721 read in full). Still "not a search", but better than none.

## Found
- **No small-q table of Nikodym(2,q) anywhere I looked.** Tao 2511.07721 (Nov 2025, the AlphaEvolve/DeepThink
  Nikodym paper) only gives asymptotics in d=2; its repo/collab [Georgiev–Gómez-Serrano–Tao–Wagner 2511.02864]
  is about constructions for large q. Lund–Saraf–Wolf 1609.01048, the 2601.20851 spread-line paper and
  2601.19879 (point-line matchings) are asymptotic too.
- **Best general lower bound** (Tao's (1.3)-bullet, citing Szőnyi–Cossidente–Gács–Mengyán–Siciliano–Weiner,
  JCD 13 (2005) 25–41, via blocking sets / "every black point on a tangent"):
  Nikodym(2,q) ≥ q² − q^{3/2} − 1 + ¼ s(1−s) q, s = frac(√q).
- Perfect-square q: unital-based constructions (Blokhuis et al., FFA 35 (2015); Tao Thm 1.2) give
  q² − q^{3/2} + O(q log q) — asymptotic only, no explicit q=9 number.

## Our values vs that bound (computed)
| q | SCGMSW bound (ceil) | ours | Kakeya(2,q) |
|---|---|---|---|
| 2 | 1 | 2 | 3 |
| 3 | 3 | 5 | 7 |
| 4 | 7 | 10 | 10 |
| 5 | 14 | 17 | 17 |
| 7 | 30 | 36 | 31 |
| 8 | 41 | 48 | 36 |
| 9 | **53** | [53,63] | 49 |

Notes:
- All exact values sit strictly above the bound — consistent, no contradiction.
- **Nik(2,9) lower end 53 is exactly the SCGMSW bound** (s=0 at q=9 → 81−27−1). If the nikodym-fq session
  derived 53 independently, say so; if it took it from the bound, cite it. Upper 63 vs unital heuristic
  ~54+O(q log q): plenty of room.
- Small-q curiosity: Nik(2,q) < Kakeya(2,q) for q=2,3, equal for q=4,5, larger from q=7 on. Fine by (1.2)
  (Nikodym ≥ Kakeya − ~2q) but worth one sentence in the README — reviewers will ask.
- Caveat on the dictionary: the affine Nikodym complement is NOT literally a projective largest-minimal-blocking-set
  object (points outside C also need an external line), so blocking-set tables (e.g. Blokhuis–Brouwer–Wilbrink
  small-p) don't transfer directly. Didn't find anyone computing the affine version for small q.

## Cross-check with the research write-up
nikodym-fq/WRITEUP.md line 6 already made the same comparison (AlphaEvolve paper Problem 6.1 + Tao: no specific-q values; SCGMSW bound 30/41/53 vs 36/48/63). So: Nik(2,9) ≥ 53 IS just the literature bound, and my independent pass agrees with the earlier one. README should cite SCGMSW for the q=9 lower end.

## Still owed for item 3
A proper search: MathSciNet/zbMATH for "Nikodym" + "AG(2,q)", Ferdinand Ihringer's blog, the
2511.02864 GitHub repo (might have small-q AlphaEvolve outputs = upper bounds to compare with).
