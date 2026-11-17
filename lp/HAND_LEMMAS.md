# Hand lemmas the K(3,5) >= 47 LP chain depends on (NOT re-derived in this release)

Links are to the research writeups in `../kakeya-fq/` (vendored copies of the code are in `lp/src/`).
"Validated" means the row family was checked numerically on the known Kakeya sets (q=3 13-set, q=4
Dvir 28-set, two q=5 53-sets); that is evidence, not proof. Replayed in this release: the hand
certificate for n >= 46 (`src/k35_lower46_check.py`, fresh code). Everything below is trusted as
written by the research sessions.

| lemma / row family | where stated | used for |
|---|---|---|
| Identities (I1)-(I6), planar table (P) with Blokhuis-Mazzocca K(2,5)=17; Bukh-Chao n>=39 | WRITEUP4 §2 | base model rows; K(3,5)>=46 |
| Counting families of session V ("all hand-provable; validated") | WRITEUP5 §2 | merged LP |
| r_max >= 5 (LP-certified, `lp_probe_rmax4.json`) | WRITEUP5 §4 | reduction step 1 |
| twelve r_max configurations, classes (5,21),(5,18),(5,17),(5,16),(5,15),(6,25),(6,21),(6,19),(6,18),(6,17); D<=7 | WRITEUP5 §5 | reduction step 2 |
| Lemma A (sections of a (5,15)-plane) | WRITEUP6 §1 | F5 rows |
| Lemma B (sixth-direction pencil) | WRITEUP6 §2 | F6 rows |
| Lemma C (pencil degree cut: s_l >= 25 + 5r - D(k-5)) | WRITEUP6 §6 | F7 rows; D >= 5 |
| (5,16) excluded at D=5,6,7 (LP-certified) | WRITEUP6 §8 | reduction step 4 |
| Lemma D | WRITEUP7 §1 | session VII rows |
| Lemma C for unmarked directions | WRITEUP8 §1 | session VIII rows |
| pairwise m-free lift (validity proof) | WRITEUP13 "The new pairwise m-free lift — proven validity" | D=5 root certificates |
| projection-table correction ((1,1,2,2,3,3) minimum is 21 not 22; exhaustive 2^25 enumeration, two implementations) | WRITEUP14 "Inherited projection-table correction" | corrected `union_tables_q5.json` |
| additive cuts / within-line pair costs (exact Hungarian pricing); "using the reduction on every maximal-r plane"; degree-pair products and degree-vector hull | WRITEUP14 "Strongest current additive cuts", "Stronger within-line pair costs", "Using the reduction on every maximal-r plane", "Global degree-pair products..." | D=6,7 root certificates |
| coplanar triples at a degree-seven point; plane-size filter; rounded plane-budget inequalities; triples-minus-quadruples point cost | WRITEUP14 sections marked "proven" | D=7 certificates |

Degree range: D >= ceil(155/46) = 4 and 1+4D <= 46 so D <= 11 (hand, WRITEUP14 "Re-derivation of the inherited reduction").
Final coverage: D=5: 5198, D=6: 51244, D=7: 274418 cases, eight root certificates (WRITEUP14 "Final coverage and fresh-price replay").
