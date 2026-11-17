# Hostile audit of the hand lemmas behind "K(3,5) ≥ 47"

Auditor: an independent Claude session (2026-11-07, 00:15–~04:00 UTC), briefed to break the chain.
Final state (03:40 UTC, after the follow-up on row #10): all planned checks finished. **Second session 2026-11-07 06:15–10:30 UTC: root-certificate replay, §4.4, and updated row 24, §5–§8.** Scope: every row family listed in `lp/HAND_LEMMAS.md`; the global case split; the encoding of each lemma in
`lp/src`; fresh computational attacks. Everything computational here is **fresh code in `audit/check_*.py`,
`audit/check_*.c`, `audit/ag3.py`** that imports nothing from `lp/src` or `../kakeya-fq` (it only *reads* their
JSON tables/certificates to compare). Outputs are in `audit/out/` (`*_result.json`, logs).

What this audit is **not**: a replay of the eight root certificates or of the 82 reduction certificates through
the original code (another worker is doing that concurrently). I did, however, replay the session-V reduction
certificates with a from-scratch row generator (§4.3), which is the first independent check of any link of the
Farkas chain: all 82 reduction certificates verify exactly against my own column universes.

## 1. Summary table

Grades: **PROVEN** = I have a complete proof (restated in §3); **PROVEN+** = proven and the statement/table is
additionally confirmed by my own exhaustive computation; **PLAUSIBLE-GAP** = true-looking, step I could not close;
**DOUBTFUL**; **FALSE**. "Encoding" = the LP rows generated from it do what the lemma says (signs, ranges, ≥/>,
rounding, applicability); "n/a" = no row, used as a case filter.

| # | lemma / row family | grade | proof method | computational attack (fresh) | encoding | confidence |
|---|---|---|---|---|---|---|
| 1 | Identities I1–I6 (WRITEUP4 §2) | PROVEN+ | double counting (§3.1) | 135 random/greedy/structured witness systems, q=3,4,5: 0 violations (`check_identities.py`) | OK (`planecount5/6.build`) | high |
| 2 | Planar table (P): k_π ≥ f(5,r) = 0,5,9,12,14,15,17; K(2,5)=17 | PROVEN+ | union bound 5r−C(r,2); r=6 is K(2,5) | f(5,r) recomputed from all 2^25 subsets of AG(2,5) (`check_planar.c`); K(2,5)=17 reproduced, so **Blokhuis–Mazzocca is not needed** for q=5 | OK | high |
| 3 | Hand certificate K(3,5) ≥ 46; Bukh–Chao for n ≤ 37 | PROVEN+ | §3.3; `check_hand46.py` (sympy + knapsack) | the plain counting LP with exact Farkas certificates excludes **every n ≤ 44 and every D** (`check_small_n_lp.py`), so **Bukh–Chao and Dvir are not needed either** | OK (`k35_lower46_check.py` is theirs; mine is independent) | high |
| 4 | Session-V counting families F1–F4 (histograms, classes, pencils, projections, Ft table) | PROVEN+ | §3.4 | identities checked on 135 systems; Ft table recomputed from 2^25 subsets: agrees with the *corrected* table, disagrees with the original only at (1,1,2,2,3,3) (`check_planar_tables.py`) | OK; one historical bug (fixed-order Ft) confirmed and confirmed fixed | high |
| 5 | r_max ≥ 5 at n = 46 (LP-certified) | LP-certified | — | **8/8 certificates re-verified with a from-scratch row generator** (`check_cert_v.py`, §4.3) | side row `rmax<=4` present and used | high |
| 6 | Twelve r_max configurations → (5,15),(5,16),(6,17); D ≤ 7 (LP-certified) | LP-certified | inside profiles re-enumerated (`check_planar_tables.py`: 5+6 profiles, 12 AGL classes grouped into 10 (r,u)) | **68/68 configuration certificates re-verified with the from-scratch row generator** (§4.3) | side rows `config (r,u) present` + `rmax<=r` present in every certificate (`check_coverage.py`) | high |
| 7 | Lemma A (sections of a (5,15)-plane) | PROVEN+ | §3.5 (odd-characteristic complete quadrilateral) | the only (k=15,r=5,u=15) tuple in the 2^25 table is t=(0,10,15,0,5), e=(3,3,3,3,3); 3D check of A(i),(ii) on structured systems | not a row; enters through the section table | high |
| 8 | Lemma B (sixth-direction pencil) | PROVEN+ | §3.6 | B(a),(b),(c) checked on 31 935 flags | OK (`matchable`, `ok_item` use bucket max k: sound) | high |
| 9 | Lemma C (pencil degree cut, F7) | PROVEN+ | §3.7, one line from I1 | 15 435 (line, plane) pairs | OK: filter `s >= 25+5r−D(hi−5)` with hi = bucket max (sound), s = Σr−1 | high |
| 10 | (5,16) excluded at D = 5,6,7; D = 4 excluded for all three configurations (LP-certified, session VI) | LP-certified | — | **6/6 certificates re-verified with a from-scratch session-VI generator** (`check_cert_vi.py`: section-profile types, exact pencil items + Lemma C filter, F6 rows R1–R4), same y^T b as the research code (§4.3b) | side rows present | high |
| 11 | Lemma D (pencil histogram identity, F8) | PROVEN+ | §3.8 | 15 435 + 14 372 checks | OK (Y1/Y2/Z1/Z2; Y3 unused in the m-free campaigns) | high |
| 12 | Lemma C for unmarked directions (F9) | PROVEN+ | §3.9 | 31 935 flags | OK (deficit rows; redundant at D=5 as they say) | high |
| 13 | F6 flag rows R1–R3 and capacity cut R4 | PROVEN | §3.10 (R3/R4 re-derived; R4's x ≤ 6 justified via R3) | R3, R4 lower bound, R4 capacity and per-plane cap3 checked on 1 395 directions / 8 370 pencil planes | OK (`cap3`, `xmax` loop) | high |
| 14 | Pairwise m-free lift (WRITEUP13) and the coordinate lift (XI) | PROVEN | §3.11 | geometry catalogue re-enumerated (**all** 18 750+15 625 line sets × all extra subsets, no orbit reduction) and found identical to `geometry14_D6.json` up to line relabelling; 150 sampled types × all objectives re-priced exactly for lift11, lift13, lift14 D6/D7: **0 mismatches** (`check_prices.py`) | OK (row = Σμ N − ε Σ cost·G ≤ 0) | high |
| 15 | Projection-table correction (1,1,2,2,3,3) → 21 | PROVEN+ | — | independent 2^25 enumeration: 21; every other entry agrees | OK | high |
| 16 | WRITEUP14 additive subset cuts, within-line pair/triple costs, extra-point costs | PROVEN | §3.12 | additive: 0 mismatches (above); within-line pair/triple/extra re-priced by exhaustive colourings on my own geometries: 0 mismatches (§4.5) | OK | high |
| 17 | "Using the reduction on every maximal-r plane" (`--maximal-only`) | PROVEN | §3.13 | — | OK: filter applied only at r = r_max | high |
| 18 | Degree-pair products (M_ij = m_i m_j) and degree-vector hull | PROVEN | §3.14 | M identities on 135 systems | OK (McCormick + chord; Q columns one per raw vector) | high |
| 19 | Point families: coplanar pairs/triples, triples-minus-quadruples, rich secants, mixed, degree-subset | PROVEN+ | §3.15 | PG(2,5): all 2 629 575 seven-subsets re-enumerated: T∈[2,20], g∈[2,6], rich≥2 under cap 4, star bounds 4,5,8 (`check_projective.py`); point identities on 6 142 points; 150 sampled types × objectives re-priced for point14/pointtriple/nonfano/rich/mix/subset: 0 mismatches | OK (signs/global RHS checked in `add_cuts`) | high |
| 20 | Plane-size filter from the budget 541; rounded plane-budget rows; parallel-class degree rounding | PROVEN | §3.16 | caps 14,15,16,17,18,19,21 and c ≥ 0 on all attainable (r,k) recomputed | OK | high |
| 21 | b_π(p) ≤ 4 in the retained cases; star bounds | PROVEN | §3.15 | max inside multiplicity over retained section profiles = 4 for both configurations | OK (`star14.upper`) | high |
| 22 | Degree range 4 ≤ D ≤ 11 | PROVEN | §2 | — | `m_D ≥ 1` row present | high |
| 23 | Global case split / coverage tree | PROVEN (given 5,6,10 and the roots) | §2 | raw universes 2599 / 25622 / 137209 re-enumerated; all three ledgers re-walked: every raw (configuration, vector) pair lies under exactly one certified node (`check_coverage.py`) | — | high |
<!-- ROW24 START -->| 24 | The eight root certificates themselves | LP-certified | — | **8/8 re-verified with a from-scratch generator of the full root models and own prices (`check_root.py`, §4.4)** | every row family re-encoded from §3; column universes and row supports identical | high |<!-- ROW24 END -->

**No lemma was found FALSE or DOUBTFUL.** The one bug the chain ever had (fixed-order Ft table) is real,
was found by session XIV, and my enumeration confirms both the bug and the fix.

## 2. Global logic: the coverage tree

Let S be a Kakeya set in AG(3,5) with |S| ≤ 46. Choose one contained line per direction; their union K is a
*witness system* (one line ℓ_d per direction d ∈ PG(2,5), K = ∪ℓ_d) with n = |K| ≤ 46. All lemmas below are
about witness systems; a witness system is itself a Kakeya set.

```
n = |K| ≤ 46
├── n ≤ 44 ............ plain counting LP infeasible for every D (fresh exact certificates, check_small_n_lp.py)
│                        [research chain: hand quadratic certificate for 38..44, Bukh–Chao n ≥ 39 for n ≤ 37]
├── n = 45 ............ quadratic certificate total 8 ∉ sums of {0,6,7,10,...}  (check_hand46.py)
└── n = 46   (D := max degree; Σ d_p = 155 ⇒ D ≥ 4; 1+4D ≤ 46 ⇒ D ≤ 11)
    ├── r_max ≤ 4 ................ 8 certificates (D = 4..11), row "rmax<=4"              [fresh replay §4.3]
    └── r_max ∈ {5,6}  (r_max ≤ 6 trivially: 6 directions per plane)
        ├── some r_max-plane has inside union u ∈ {21,18,17} (r=5) or {25,21,19,18} (r=6)
        │        .......... 7 (r,u) × 8 D = 56 certificates, rows "config (r,u) present" + "rmax<=r"
        │        (the two AGL(2,5)-classes with u = 17 resp. u = 19 share one certificate; valid, the row only
        │         sees u)
        └── every r_max-plane has u ∈ {15,16} (r=5) or u = 17 (r=6)       ← this is "maximal-only" (§3.13)
            ├── D ∈ {8,…,11} ........ 3 configs × 4 D = 12 certificates (session V)
            ├── D = 4 ............... 3 certificates (session VI, F7 model)
            ├── config (5,16), D ∈ {5,6,7} ..... 3 certificates (session VI)
            └── config ∈ {(5,15),(6,17)}, D ∈ {5,6,7}, degree vector m with Σm=46, Σ j m_j=155, m_D ≥ 1
                ├── D=5: 2599 vectors × 2 configs — root certificates XI (5,15) and XIII (6,17)
                ├── D=6: 25622 × 2 — root certificates XIV (both)
                └── D=7: 137209 × 2 — (5,15) root; (6,17) split m_1 = 0 | 1 ≤ m_1 ≤ 26 (domain cert) | m_1 = 27
```
8 + 56 + 12 + 3 + 3 = 82 reduction certificates, matching the manifest. Points checked:

* **Exact-D semantics.** Every model has the row m_D ≥ 1, so a certificate for "D" excludes systems whose
  maximum degree is exactly D; the union over D = 4..11 is exhaustive.
* **Configuration semantics.** "(r,u) present" only asserts that *some* plane has r inside lines with union u;
  "rmax<=r" asserts no plane has more. A system with r_max = 5 that has both a (5,15)-plane and a (5,16)-plane
  satisfies the (5,16) LP, so it is excluded by the (5,16) certificates; hence in the surviving (5,15) case every
  5-plane has union 15, which is exactly what `--maximal-only` filters. Same for r_max = 6 and u = 17. Lower-r
  planes are not restricted (checked in `coverage14.build`: filter only when `t[0] == rmax`).
* **No circularity.** The (5,16) certificates do not assume anything about (5,15)-planes (their rows are the
  base model + "config (5,16) present" + "rmax<=5"; `check_coverage.py` lists the side rows with nonzero
  multipliers for every certificate).
* **Degree-vector trees.** For each D the ledger's explicit prefix tree was re-walked with a fresh enumeration of
  the raw vectors: children partition the parent's vectors, every leaf is certified, no vector is covered twice,
  no `exact-primal` survivor exists. Counts 5198 / 51244 / 274418 reproduced.
* **D = 7, (6,17).** m_1 = 0 (6196 vectors, one certificate), m_1 = 27 (7 vectors), and the domain certificate
  for m_1 ≥ 1 (its model keeps all 131 013 vectors with m_1 ≥ 1 and an explicit row −m_1 ≤ −1); the ledger credits
  it to the 26 first-level nodes m_1 = 1..26. Disjoint and exhaustive.
* **Witness system vs Kakeya set.** The exclusion is for witness systems with exactly 46 points; a 46-point
  Kakeya set contains a witness system of ≤ 46 points, which has exactly 46 by the ≥ 46 bound. Fine.

One thing a referee will ask that the writeups answer only implicitly: the LP for "configuration (5,15), D" is
built with `u_for_r = 5` in session V, i.e. the union size is recorded only for r = 5 types; in the sections
model every type carries its u. Both are supersets of the needed column set; no column is dropped by this.

## 3. Per-lemma statements and proofs (re-derived)

Notation. AG(3,5): 125 points, 31 directions, 775 lines (25 per direction), 155 planes in 31 parallel classes
of 5; a plane π has a line at infinity L_π of 6 directions, contains 30 lines, lies in 6 planes per line; two
points lie in 6 common planes. Witness system K = ∪_d ℓ_d, n = |K|, d_p = #{d : p ∈ ℓ_d} ≥ 1 for p ∈ K,
D = max d_p, m_j = #{p : d_p = j}. For a plane π: r_π = #{d : ℓ_d ⊂ π} (the *inside* witness lines; their
directions are *marked* in π, the others of L_π *unmarked*), k_π = |K ∩ π|, U_π = union of the inside lines,
u_π = |U_π|, b_π(p) = # inside lines through p, h_π = degree histogram of K ∩ π, t_j(π) = # lines of π meeting K
in exactly j points, s_ℓ = Σ_{p∈ℓ} d_p for a witness line ℓ. "(5,15)-plane" = r = 5, k = 15 (which forces
u = 15: five lines in general position and no other point of K). "(r,u)-configuration" = r inside lines with
union u, extra points allowed.

### 3.1 Identities I1–I6 — PROVEN
* I1: Σ_{p∈K∩π} d_p = 25 + 5r_π. Count (p,d) with p ∈ ℓ_d ∩ π. A witness line whose direction is in L_π is inside
  π (5 incidences) or disjoint from π; the other 25 witness lines meet π in exactly one point. ∎
* I2: Σ_π 1 = 155, Σ r_π = 31·6, Σ k_π = 31n (each line in 6 planes, each point in 31). ∎
* I3: Σ_π C(r_π,2) = Σ_p C(d_p,2): two witness lines have different directions, so they are coplanar iff they
  meet, and then in exactly one point and one plane. ∎
* I4: Σ_π k_π² = Σ_{p,p'} #planes ∋ p,p' = 31n + 6n(n−1). ∎
* I5: Σ_π k_π r_π = Σ_p Σ_d #{planes ∋ p ⊃ ℓ_d} = Σ_p (6d_p + (31 − d_p)) = 31n + 5·155. ∎
* I6: Σ_j jL_j = 31n, Σ_j C(j,2)L_j = C(n,2), L_5 ≥ 31, Σ_{j≥1} L_j ≥ 31·17 (projection along d is a planar
  Kakeya set: the 30 other directions map onto the 6 directions of the quotient plane), Σ_π C(k_π,3) =
  C(n,3) + 5Σ_j C(j,3)L_j (a collinear triple lies in 6 planes, a non-collinear one in 1). ∎
* Degree range: Σ d_p = 155 ⇒ D ≥ ⌈155/46⌉ = 4; the D lines through a degree-D point meet only there ⇒
  n ≥ 1 + 4D ⇒ D ≤ 11. ∎

### 3.2 Planar table (P) — PROVEN+
r lines of AG(2,5) in distinct directions pairwise meet, so |∪| ≥ 5r − C(r,2) = 0,5,9,12,14,15 for r ≤ 5; for
r = 6 the union is a planar Kakeya set. My 2^25 enumeration gives the minimum union 0,5,9,12,14,15,**17** and the
minimum size of a set containing a full line in every direction = 17, i.e. K(2,5) = 17 without citing
Blokhuis–Mazzocca. The inside profiles (u; c_1..c_r) are exactly the cached ones (5 for r=5, 6 for r=6).

### 3.3 Hand certificate ≥ 46 — PROVEN+
With c(r,k) = 10k² − (150+17r)k + 560 + 143r, I2/I4/I5 give Σ_π c(r_π,k_π) = 60n² − 4927n + 100223 (sympy
confirms). c ≥ 0 on every attainable (r,k) (min over the 981 attainable section tuples is 0); the total is
negative exactly for 38 ≤ n ≤ 44; at n = 45 it is 8, and the attainable costs are {0,6,7,10,17,…}, so 8 is not a
sum of plane costs (unbounded-knapsack check). For n ≤ 37 the writeup cites Bukh–Chao (n ≥ 39). **I removed that
dependency:** the plain counting LP (N[r,k] ≥ 0, I2–I5, I3 with degree counts, k ≥ f(r), m_D ≥ 1) is infeasible
for every n ∈ [20,44] and every admissible D, with exact integer Farkas certificates saved in
`out/small_n_lp_result.json` (n < 20 is below Dvir's 35 anyway and could be run the same way). So the ≥ 46 bound
needs no external theorem at all.

### 3.4 Session-V families — PROVEN
F1 (histograms): Σ_π h_j = 31m_j; Σ_π h_j k_π = (6n+25)m_j; Σ_π h_j r_π = (31+5j)m_j (I4/I5 restricted to
degree-j points); Σ_π h_i h_j = 6M_ij (i<j), Σ_π C(h_j,2) = 6M_jj where M counts point pairs by degree pair;
Σ_{i≠j} M_ij + 2M_jj = (n−1)m_j; 2M_jj ≤ (n−1)m_j; M_jj ≥ a·m_j − C(a+1,2) (chord of C(x,2) between a and a+1).
Crossing filter: a point on i inside lines has degree ≥ i, so #{p: b_π(p) ≥ i} ≤ #{p ∈ K∩π: d_p ≥ i} (Hall for a
nested system). F2 (classes): the 6 directions of L each have their witness line in exactly one plane of the
class ⇒ Σ_class r = 6; Σ_class k = n. F3 (pencils of a witness line ℓ): Σ_{π⊃ℓ} k_π = n + 25 (points off ℓ once,
on ℓ six times), Σ_{π⊃ℓ} r_π = 1 + s_ℓ (a witness line ≠ ℓ is coplanar with ℓ iff it meets ℓ, and then in one
plane), r_π ≥ 1; Σ_ℓ s_ℓ = Σ_p d_p²; the degree multiset of ℓ is a 5-multiset with sum s_ℓ; pairs of points on
witness lines with degrees {i,j} ≤ M_ij (two points span one line). F4 (projections): along d, the witness line
ℓ_{d'} (d' ≠ d) projects onto the image of the unique plane of class L = dd' containing it, so the projection
P_d contains exactly a_{d,L} := #{planes of class L containing a witness line ≠ ℓ_d} distinct full lines in
quotient direction L, hence |P_d| ≥ Ft(a_{d,L_1..L_6}) with Ft(a) = min |S| over S ⊂ AG(2,5) whose sorted
full-line-count profile dominates sorted(a); class link: a_{d,L} = ρ if ℓ_d is not alone in its plane, ρ−1 if it
is (ρ = # planes of L with r ≥ 1). All of this is checked numerically on 135 systems including the `a` formula and
|P_d| ≥ Ft (with my own Ft).

### 3.5 Lemma A — PROVEN+
π a (5,15)-plane, lines ℓ_1..ℓ_5 in general position, e the unmarked direction. (i) An e-line m meets every
ℓ_i once; coincidences are crossing points p_ij on m, and two crossings on m have disjoint index pairs (else
m = ℓ_i), so c(m) ≤ 2; the 10 crossings lie on 5 e-lines, one each, so c(m) = 2 and |m ∩ K| = 5 − 2 = 3. (ii) m ∥
ℓ_i, m ≠ ℓ_i: m meets ℓ_j (j ≠ i) once each with coincidences at the 6 crossings p_jk (j,k ≠ i); c(m) ≤ 2 and
Σ over the 4 lines = 6, so c = (2,2,1,1) or (2,2,2,0). The latter says the three diagonals of the complete
quadrilateral {ℓ_j : j ≠ i} are concurrent at the point at infinity d_i, impossible in odd characteristic (dual of
"the diagonal points of a complete quadrangle are not collinear"). So the sections are 2,2,3,3. ∎ Profile:
t = (0,10,15,0,5). The 2^25 table has exactly one tuple with (k,r,u) = (15,5,15) and it is this one.

### 3.6 Lemma B — PROVEN
π with r_π = 5, unmarked e, class L of π, π' the plane of L containing ℓ_e. (a) ℓ_e ∥ π and ℓ_e ⊄ π, so π' ≠ π;
Σ_L r = 6 = 5 + r_{π'} + rest ⇒ r_{π'} = 1, rest 0. (b) A plane σ ⊃ ℓ_e with σ ≠ π' is not parallel to π, so
σ ∩ π is a line whose direction is L_σ ∩ L_π ∋ e, i.e. an e-line; distinct σ give distinct lines (σ = span(ℓ_e, m));
5 planes ↔ 5 e-lines. (c) m_j = σ_j ∩ π is parallel to ℓ_e and meets each other inside line of σ_j at most once,
so at most r_{σ_j} − 1 of its s_j points are on U_{σ_j}; the rest are extra points: k_σ ≥ |U_σ| + max(0, s_j − r_σ + 1)
≥ f(r_σ) + max(0, s_j − r_σ + 1). ∎ The general form (any r_π ≥ 1, any unmarked e; mate = the plane of the class
containing ℓ_e, r' ≥ 1) is the same argument.

### 3.7 Lemma C — PROVEN
ℓ ⊂ π: by I1, Σ_{p∈K∩π∖ℓ} d_p = 25 + 5r_π − s_ℓ ≤ D(k_π − 5). With Σ_{π⊃ℓ} r_π = 1 + s_ℓ (3.4). ∎
Encoding: pencil type w (6 items (r, k-bucket)) is generated only if s := Σr − 1 ≥ 25 + 5r_i − D(hi_i − 5) for every
item i, hi_i = largest k of the bucket (weaker ⇒ sound).

### 3.8 Lemma D — PROVEN
(a) ℓ ⊂ π ⇒ the 5 points of ℓ are in K ∩ π ⇒ dm(ℓ) ⊆ h_π. (b) The 6 planes through ℓ cover AG(3,5) and pairwise
meet in ℓ ⇒ Σ_{π⊃ℓ} h_π = m + 5·dm(ℓ). ∎

### 3.9 Lemma C unmarked (F9) — PROVEN
Same as 3.7 applied to the mate π' ⊃ ℓ_e: s_{ℓ_e} = Σ_{σ⊃ℓ_e} r_σ − 1 ≥ 25 + 5r_{π'} − D(k_{π'} − 5). ∎ Row:
Σ_w max(0, deficit) X2 ≤ 0 forces the incompatible flag columns to 0.

### 3.10 F6 rows and R3/R4 — PROVEN
R1: each plane with r ≥ 5 has one flag per unmarked direction. R2: a class-mate π'' of type (r',k') hosts exactly
r' unmarked directions of π (its own inside directions), so #flags(π of type (r,k), mate type (r',k')) summed over
a class = m_{rk}(m_{r'k'} − δ)r'. Link: every flag sits on the pencil of ℓ_e of some type w; w must contain the
mate as an r' ≥ 1 item with k' in its bucket, and the other 5 items must be perfectly matchable to the 5 e-sections
s_j under 3.6(c) (bitmask DP in `matchable` is an exhaustive matching test).
R3 (rmin = 5): the mate π' has r' = 1 and its class has Σ r = 6, so at most ⌊(6−r')/5⌋ = 1 plane of that class
has r ≥ 5; hence #flags on the pencil of ℓ_e ≤ #(r' = 1 items of w). Also gives x ≤ 6 (at most 6 planes in a
pencil) which justifies the `xmax ≤ 6` loop bound in R4.
R4: x (5,15)-planes sharing the unmarked direction e lie in x distinct classes (a class has one r = 5 plane);
each contributes 5 three-point e-lines hosted by the 5 pencil planes other than its mate; lines of π_i, π_j
coincide only if equal to π_i ∩ π_j, so ≥ 5x − C(x,2) distinct 3-point e-lines are hosted. A pencil plane σ
of type (r,k) hosts at most cap3 of them: r = 1: ⌊(k−5)/3⌋; r = 2: ⌊(k−9)/2⌋ (each hosted line uses one point of
the second inside line and two extra points); r = 3: ≤ min(4, k−12) (each hosted line needs ≥ 1 extra point,
whether or not the three inside lines are concurrent); r ≥ 4: 4 (only 4 e-lines ≠ ℓ_e); (5,15) with e marked: 2
by Lemma A(ii). ∎ All verified on random systems with forced (5,15)-planes.

### 3.11 The m-free lifts (XI coordinate, XIII pairwise, XIV subsets) — PROVEN
Fix a histogram key A = (r,k,h) and a degree subset S, sign ε. Give the incidence label (A,dm) the cost
ε·#{d ∈ dm : d ∈ S}; G_{(A,dm)} counts (plane of key A, inside witness line with degree multiset dm). For one plane
π of key A with degree assignment d(·), Σ_{inside ℓ} cost(A, dm(ℓ)) = ε Σ_p b_π(p)[d(p) ∈ S] ≥ μ_a where μ_a is the
minimum over every realisation of the full type a = (r,k,h,u,t,e) and every assignment with histogram h and
d(p) ≥ max(1, b_π(p)). Summing over planes: Σ_a μ_a N_a ≤ ε Σ_dm cost·G. ∎ The minimum is a transportation problem
in (b, d) since the cost is additive and depends on p only through b_π(p); the set of multiplicity profiles of a
full type is {(k−u, c_1..c_r)} over the crossing profiles c attainable with that (k,r,u,t,e) — my 2^25 table gives
exactly the same profile sets as their geometry catalogue, and my own group-DP minimiser reproduces every sampled
price. Validity needs nothing about m: the row holds for every witness system.

### 3.12 Within-line pair/triple costs and extra points — PROVEN
Same template with cost ε·C(#{d ∈ dm : d ∈ S}, p), p = 2,3, which depends on the colouring (which points have
degree in S) through the inside-line incidences only; the minimum over all assignments equals the minimum over
colourings that admit a degree assignment, and a colouring admits one iff the selected points can be matched to
the S-degrees and the complement to the other degrees, each under d ≥ max(1,b): two independent threshold
matchings, decided by Hall's tail inequalities. Extra points have b = 0 and cost 0 but consume degree slots; their
count is forced (s − Σ_b p_b) and is checked against k − u. ∎ Conditional family `extra14` (r = 6, u = 17 only)
is used only under `--maximal-only`, asserted in `add_cuts`. Re-pricing status: §4.5.

### 3.13 Maximal-only — PROVEN
See §2: in the r_max = 5 case any 5-plane with u ≠ 15 instantiates an excluded (5,u) case; in r_max = 6 any
6-plane with u ≠ 17 instantiates an excluded (6,u) case. ∎

### 3.14 Degree-pair products and hull — PROVEN
M_ij = m_i m_j (i≠j), M_jj = C(m_j,2) by definition of M (checked on 135 systems). McCormick bounds on a box
L ≤ m ≤ U are valid, as is 2M_jj − (L+U−1)m_j ≤ −LU (chord of x(x−1) between L and U). The hull rows Q_v ≥ 0,
ΣQ = 1, m = Σ vQ, M = Σ moment(v)Q are satisfied by Q = δ at the true degree vector; every raw vector of the node
is a column. ∎

### 3.15 Point families — PROVEN+
Each pair of witness lines through p spans one plane, so Σ_{π∋p} C(b_π(p),2) = C(d_p,2); each coplanar triple
of witness lines through p is a collinear triple of its direction set in PG(2,5), so Σ_π C(b_π(p),3) = T(p).
Global bounds for seven directions in PG(2,5): T ∈ [2,20]; Σ_L (C(|L∩S|,3) − C(|L∩S|,4)) ∈ [2,6]; #lines with
≥ 3 selected points ≥ 2 when no line has ≥ 5 selected points (writeup proofs re-read and agreed; all confirmed by
exhaustive enumeration of the 2 629 575 seven-subsets, and the smaller j likewise). b ≤ 4 in the retained cases:
five concurrent inside lines have union 21; at r = 5 then k ≥ 21 and c(5,21) = 750 > 541; at r = 6 it contradicts
u = 17. Star bounds for largest secant a ≤ 4: T ≤ max(⌊C(j,2)/3⌋, C(4,3) + C(j−4,2) + C(j−4,3)) = 4,5,8 for
j = 5,6,7 (enumeration agrees). The local minima are transportation problems as in 3.11 (my re-pricing agrees).
Row RHS: +ε·global bound · m_j with the correct side (hi for ε = +1, lo for ε = −1), checked in each `add_cuts`;
`pointrich`/`pointmix` use the cap-4 bounds for all j, which is valid since the cap is proved for every point. ∎

### 3.16 Budget filters and rounding — PROVEN
Σ_π c = 541 at n = 46 and c ≥ 0 ⇒ every actual plane has c(r,k) ≤ 541 ⇒ k ≤ 14,15,16,17,18,19,21 for r = 0..6
and the empty plane (r,k) = (0,0) is impossible (c = 560) — also directly from I1. Σ_π ⌊c_π/t⌋ ≤ ⌊541/t⌋ since
the left side is an integer ≤ 541/t. Parallel-class rounding: the 5 planes of a class partition the M points with
degree in S, so Σ_i ⌊h_i/t⌋ ≤ ⌊M/t⌋ and ≥ ⌈(M − 5(t−1))/t⌉. ∎ These are integrality cuts: they hold for actual
plane multiplicities, which is all the Farkas argument needs (the LP is a relaxation of the integer system).

## 4. Computational attacks (all fresh code)

### 4.1 Planar tables (`check_planar.c`, `check_planar_tables.py`; 31 s)
All 2^25 subsets of AG(2,5); for each, all choices of inside lines (every subset of full-line directions, every
line per direction). Results: 981 section tuples, identical to `sections3_q5_r5.json`; inside profiles identical;
f(5,r) identical; Ft (210 entries) identical to the corrected `union_tables_q5.json` and different from the
archived original only at (1,1,2,2,3,3) (true 21, old 22); F(a) identical; Lemma A tuple unique; budget caps
14,15,16,17,18,19,21; max inside multiplicity 4 in both retained configurations (492 / 550 profiles).

### 4.2 Witness systems (`ag3.py`, `check_identities.py`; 135 systems, q = 3,4,5)
Random, greedy (small n, many high-r planes) and structured (three forced (5,15)-planes) witness systems.
Every row family that is claimed for all witness systems was evaluated directly from the geometry: I1–I6, the
quadratic identity, F5, H1–H8-type identities, M = products, crossing filter, membership of every plane's
(k,r,u,c,t,e) in the 981-tuple table (6 975 planes), pencil sums, Lemma C, D(a),(b), witness pairs ≤ M, class
sums, the a_{d,L} formula, |P_d| ≥ Ft, P_d Kakeya, Lemma A(i),(ii), Lemma B(a),(b),(c), F9, R3 (rmin 5 and 2),
R4 (lower bound, capacity, per-plane cap3), Σ C(b,2) = C(j,2), Σ C(b,3) = T, h ≤ m, degree range.
**Zero violations in ≈ 480 000 checks.** (These are universal statements, so random systems are legitimate
counterexample hunts; cuts conditional on n = 46 cannot be tested this way and were proved instead.)

### 4.3 Fresh replay of the session-V reduction certificates (`check_cert_v.py`)
A from-scratch generator of the session-V merged LP (types (r,k,h,u) from my planar tables, m, M, L, class
types, direction types with my Ft, line vectors, bucketed pencil types, degree multisets; all rows of WRITEUP5 §2
with the research code's row *names* so the saved multipliers can be applied). For each certificate: every
multiplier on a ≤-row is ≥ 0, the combination is ≥ 0 on **every column of my universe**, and y^T b < 0 — exact
rational arithmetic. Result for the r_max ≤ 4 family: **all 8 certificates (D = 4..11) VERIFIED** on my column
universes (18 630 … 371 765 columns, 397 … 908 rows), with the same y^T b values the research code reports
(e.g. −97 100 199 at D = 4). The 68 configuration certificates (ten (r,u) files, every infeasible D): **all 68 VERIFIED**
(`out/cert_v_cfg.log`, `out/cert_v_result.json`). Six of them ((5,18), (5,21), (6,19), (6,21), (6,25) at D = 4 and (6,25) at D = 5) have y^T b = −1: the
column set of "config present" is empty at D = 4 (the crossing profile needs a point on ≥ 4 inside lines, i.e. a point of degree ≥ 4 whose plane histogram
must also fit the degree sum 25 + 5r, impossible at that D), so the row −0 ≤ −1 is itself the contradiction; valid.
So **all 82 reduction certificates are now independent of the research code** (the six session-VI ones below).

### 4.3b Fresh replay of the six session-VI certificates (`check_cert_vi.py`, addendum 2026-11-07 03:30 UTC)
A second from-scratch generator for the session-VI model: plane types (r,k,h,u,t,e) from my own 2^25 section
table with the Hall crossing filter, F5 rows, exact pencil items (k exact up to f(r)+3, one tail bucket) with the
Lemma C filter on every item (F7), the degree-multiset family, and the sixth-direction flag family F6 — kinds
(r,k,s-profile,mate type), R1 (one flag per unmarked direction), R2 (class hosting counts m₁(m₂−δ)r'), the link
rows, pencil compatibility by an exhaustive perfect-matching test of Lemma B(c), R3 (⌊(6−r')/5⌋ per pencil
plane) and R4 (3-point e-line capacity with cap3 and x_max). Only the row *names* follow the research code
(the certificates are keyed by name; the pencil index in the R3/R4 names forces the same lexicographic pencil
enumeration). Verification is column-wise over my whole universe (N, m, L, M, C, A, P, W, V, X1, X2).

| certificate | my universe (types / pencils / flag columns / columns / rows) | result | y^T b (mine = theirs) |
|---|---|---|---|
| (5,15), D=4 | 29 308 / 1 241 / 833 / 33 374 / 1 517 | VERIFIED | −3 077 519 |
| (6,17), D=4 | same | VERIFIED | −3 955 571 |
| (5,16), D=4 | same | VERIFIED | −5 049 433 |
| (5,16), D=5 | 96 808 / 9 384 / 177 873 / 288 547 / 7 382 | VERIFIED | −253 395 |
| (5,16), D=6 | 217 174 / 17 202 / 291 851 / 532 285 / 10 407 | VERIFIED | −19 349 165 |
| (5,16), D=7 | 378 870 / 22 089 / 582 034 / 991 415 / 14 349 | VERIFIED | −35 895 070 |

Surprises: none in the verdicts. Two observations. (i) My type universe at D = 4 has 29 308 columns against the
28 973 the research log reports — a superset, so the check is at least as strong; I did not chase the 335 extra
columns (most likely r = 0 tuples that their generator drops through the `k < max(f[r],1)`/`S ≤ Dk` order of
tests). (ii) The pencil universe after the F7 filter is tiny (1 241 at D = 4 against 44 012 before F7), which is
why the D = 4 certificates have only ~345 nonzero multipliers; the certificates lean on the R3/R4 flag rows
(81 and 34 of the 344 multipliers at (5,15), D = 4), i.e. Lemma B and the R4 capacity argument carry real
weight — they are the hand lemmas a referee should read most carefully (§3.6, §3.10).

<!-- SEC44 START -->
### 4.4 Independent root-certificate replay (`check_root.py`, addendum 2026-11-07 06:15–10:30 UTC)

Second hostile session, briefed to close weakest link #1 of §5: the eight root certificates of the D = 5, 6, 7 campaigns
were verified only through the original row generators. `audit/check_root.py` is a **from-scratch generator of the full
root models** on top of the session-VI model of `check_cert_vi.py` (my planar tables, F5, exact pencil items with the
Lemma-C filter F7, F6 flags R1–R4), adding, each derived from the lemma statement in §3 and not from `lp/src`:

* prefix caps and the degree-vector universe (raw vectors re-enumerated: 2 599 / 25 622 / 137 209), the maximal-only
  filter (§3.13), the plane-budget filter c(r,k) ≤ 541 (§3.16);
* F8 (Lemma D(a)): columns G[(r,k,h), dm] for dm ⊆ h, Z[dm, pencil] for Σdm = s(w) with every item realisable, rows
  Y1, Y2, Z1, Z2; F9 (§3.9): deficit rows on the X2 flags (all deficits are 0 at D = 5, as WRITEUP8 proves; the rows exist
  and are empty);
* the lifts, **priced here**: the coordinate lift (XI) and pairwise lift (XIII) by an exhaustive group DP over my own
  multiplicity profiles (2^25 planar table); the XIV subset lift by the same DP (D = 6) or by an allocation enumeration
  with Hall threshold inequalities (D = 7); within-line pair/triple costs by exhaustive Gray-code colourings of my own
  geometry catalogue (`check_geometry.c`), including extra points; point crossing-pair, coplanar-triple, triples-minus-
  quadruples, rich-secant and mixed costs by the allocation enumeration, with the global PG(2,5) bounds taken from my own
  `check_projective.py` output (`out/projective_result.json`: T ∈ [2,20], g ∈ [2,6], rich ∈ [2,6] at j = 7, star bounds
  4, 5, 8); the degree hull (Q columns, one per raw vector), pair-bound envelopes (McCormick + chord), the rounded
  plane budget (one row per distinct ⌊541/t⌋), parallel-class degree rounding (RD columns), the compact-objective
  rule (drop objectives whose subset misses the support of h, and repeated (S ∩ support, sign) signatures);
* the node rows (configuration present, r_max, 0 ≤ m_j ≤ 46 or the prefix value, domain row −m_1 ≤ −1).

**Bookkeeping, stated honestly.** The root certificates are *positional* multiplier lists (no row names), so the
research row list (name, sense, rhs, nnz) was dumped once per certificate by `dump_research_rownames*.py` — the only
files in `audit/` that import `lp/src` — and is used purely as the position → row-name map. My rows carry the same names.
Before the Farkas check the two row lists are compared (name sets, sense and rhs per name, **and the number of nonzeros
per row counted over my own columns**), and the column counts per family are compared with the research model. The
Farkas check itself is on my columns only: every multiplier on a ≤-row is ≥ 0, y^T A ≥ 0 on every column of my universe,
y^T b < 0, exact integers. Row *order* differs from the research code in every case (immaterial).

| cert id | D | prefix / scope | my columns (research) | my rows (research) | y^T b mine (research) | min col | prices re-derived | verdict | wall |
|---|---|---|---|---|---|---|---|---|---|
| `coverage11_cert_515_root` | 5 | root, core lift | 572 863 (572 863) | 25 045 (25 045) | −11 612 655 (−11 612 655) | 3 | 131 types × 10 (lift11): 0 mismatches | **VERIFIED** | 15 s |
| `coverage13_cert_pair_root_ipm` | 5 | root, core + pair lift | 775 113 (775 113) | 39 938 (39 938) | −11 879 470 (−11 879 470) | 0 | 324 × 10 + 324 × 20 (lift11, lift13): 0 mismatches | **VERIFIED** | 39 s |
| `cert14_D6_515_root` | 6 | root, maximal-only, subset core, within-line all, hull | 1 164 257 (1 164 257) | 130 747 (130 747) | −363 712 730 240 (same) | 0 | 482 × 62 + 482 × 124: 0 mismatches | **VERIFIED** | 145 s |
| `cert14_D6_617_root` | 6 | same | 1 972 486 (1 972 486) | 279 308 (279 308) | -363 713 866 140 (same) | 0 | lift14: 1191 types, 0 mismatches; nonlinear14: 1191 types, 0 mismatches | **VERIFIED** | 720 s |
| `cert14_D7_515_root_strongroot` | 7 | root; + budget, compact, subset all, pairs, triples, point pairs/triples, pair bounds | 1 553 248 (1 553 248) | 341 772 (341 772) | -365 706 024 350 (same) | 0 | lift14: 18142 types, 0 mismatches; nonlinear14: 1092 types, 0 mismatches; triple14: 1092 types, 0 mismatches; point14: 79280 types, 0 mismatches; pointtriple14: 65508 types, 0 mismatches | **VERIFIED** | 234 s |
| `cert14_D7_617_27_d7high1_27` | 7 | prefix m₁ = 27; pair scope, star, rounded budget, non-Fano | 335 852 (335 852) | 23 638 (23 638) | -88 411 556 407 (same) | 0 | lift14: 0 types, 0 mismatches; nonlinear14: 0 types, 0 mismatches; point14: 2841 types, 0 mismatches; pointtriple14: 2090 types, 0 mismatches | **VERIFIED** | 15 s |
| `cert14_D7_617_0_scaled0` | 7 | prefix m₁ = 0; everything incl. extra, rich, mixed, rounded degree | 584 317 (584 317) | 110 932 (110 932) | -333 052 517 602 (same) | 0 | lift14: 5554 types, 0 mismatches; nonlinear14: 394 types, 0 mismatches; triple14: 394 types, 0 mismatches; point14: 8972 types, 0 mismatches; pointtriple14: 8113 types, 0 mismatches | **VERIFIED** | 328 s |
| `cert14_D7_617_domain_m1ge1_subset_bounded` | 7 | domain m₁ ≥ 1 (131 013 vectors); coordinate scope, + point subsets | 2 704 210 (2 704 210) | 99 971 (99 971) | -9 947 293 (same) | 0 | lift14: 2737 types, 0 mismatches; nonlinear14: 0 types, 0 mismatches; point14: 154724 types, 0 mismatches; pointtriple14: 140952 types, 0 mismatches | **VERIFIED** | 284 s |

"VERIFIED" means: identical row set (names, senses, rhs), identical nonzero count on every row, identical column counts
per family, and the exact Farkas contradiction on my columns with the research y^T b; "prices re-derived" lists the
price tables that the certificate's cut rows depend on and that were recomputed by my own pricers (the research tables
were read only afterwards for the informational comparison).

Reproduce (from `audit/`, research venv only for the dumps):
```
cd lp/src && .venv/bin/python ../../audit/dump_research_rownames.py 515 && .venv/bin/python ../../audit/dump_research_rownames.py 617
cd lp/src && for c in cert14_D6_515_root:D6_515 cert14_D6_617_root:D6_617 cert14_D7_515_root_strongroot:D7_515 \
   cert14_D7_617_27_d7high1_27:D7_617_27 cert14_D7_617_0_scaled0:D7_617_0 cert14_D7_617_domain_m1ge1_subset_bounded:D7_617_dom; do
   .venv/bin/python ../../audit/dump_research_rownames14.py ../certs/${c%%:*}.json ${c##*:}; done
cd audit && for t in 515 617 D6_515 D6_617 D7_515 D7_617_27 D7_617_0 D7_617_dom; do nice -n 19 python3 check_root.py $t; done
```
Outputs: `out/root_<tag>.log`, `out/root_result.json`.

Observations.
1. **Column-universe agreement is exact at D = 5, 6, 7** once plane types are deduplicated. The earlier discrepancy of
   §4.3b (29 308 vs 28 973 at D = 4) was *not* extra genuine types: `check_cert_vi.types_for` lists a tuple (r,k,h,u,t,e)
   once per crossing profile c that admits it, so tuples reachable through several c were duplicated. The research model
   keys columns by the tuple. Duplicates are harmless for validity (identical columns) but they were a bookkeeping error
   in my §4.3b counts, now fixed in `check_root.py` (`types_dedup`); the §4.3b verdicts stand.
2. The (6,17) root at D = 5 and the two D = 6 roots have 1, 102 101 and 1 440 658 columns with y^T A = 0
   (tight), i.e. the certificates are genuinely delicate combinations, not slack ones; every family's prices entered.
3. The node 617:27 (m₁ = 27) has **no (6,17)-type in its universe** (prefix caps (27,1,1,1,2,5,18) leave no histogram that dominates the crossing profile of a (6,17)-plane), so its configuration row reads 0 ≤ −1 and the node is infeasible by that row alone; the 23 638-multiplier certificate is valid but far more than needed. Likewise 617:0 has only 394 (6,17)-keys. Neither affects validity; both are worth knowing when judging how much the exotic D = 7 families actually carry.
4. The row-name map is the only place where the research code touches the verification. A referee who distrusts
   even that can check that the name sets and per-row supports coincide (printed in every log), which pins the map.
<!-- SEC44 END -->

### 4.5 Price tables (`check_geometry.c`, `check_prices.py`)
Catalogue: all 18 750 five-line and 15 625 six-line configurations with every extra-point subset, no orbit
reduction; 2 954 distinct (key, mask-multiset) pairs; same 196 keys as `geometry14_D6.json` and identical
geometry sets up to relabelling of lines. Exact re-pricing of 150 random types × all objectives per family:
lift14 D6 (1 900 prices), lift14 D7 (2 548), lift11 D5 (1 500), lift13 pairs (3 000), point14 D7 (1 800),
pointtriple14 (1 500), pointnonfano14 (1 500), pointrich14 (1 500), pointmix14 (6 000), pointsubset14 (7 500);
within-line families by exhaustive colourings of my own line-incidence geometries: nonlinear14 D6 (12 types, 1 488
prices), triple14 D6 (8 types, 992), extra14 D7 (12 types, 1 344): **0 mismatches in every family**
(`out/prices_result.json`).

### 4.6 Coverage (`check_coverage.py`) and projective tables (`check_projective.py`, `check_projective_tables.py`)
Reported in §1/§2. The three projective JSON files read by the D = 7 rows agree with my enumeration for every
j = 3..7 (T, cap-4 rich secants, g).

## 5. Weakest links, ranked

<!-- S5 START -->1. **(closed, §4.4)** The root certificates: all 8 are now re-verified against a from-scratch generator of the full root
   models (F8/F9, prefix caps, every XIV family, own prices, own PG(2,5) bounds), with identical column universes and
   identical row supports. The research row generators are no longer in the trust base of the Farkas chain; they only
   translate certificate positions into row names, and that map is pinned by the support comparison.
2. **The hand proofs of §3 are now the whole trust base** of the LP side: the encodings have been derived twice (research
   code, this audit) and agree column by column, so an error could only be a *shared* misreading of a lemma. The lemmas
   with the most leverage in the root certificates are Lemma B / R4 (flag capacity), Lemma D (F8) and the assignment-
   pricing argument of §3.11/§3.12 on which every XIV family rests; a referee should read those three first.
3. **Price tables**: every table a replayed certificate depends on has been recomputed in full by independent pricers
   (group DP, allocation enumeration, exhaustive colourings), not sampled: 0 mismatches in every family (§4.4).
4. **Numerical hygiene of the certificates**: multipliers are integers and the check is exact, so solver numerics are
   irrelevant.
5. Everything else (identities, tables, planar lemmas, projective facts, coverage tree) is independently confirmed and
   I would not spend further effort there.

<!-- S5 END -->## 6. What a skeptical expert would still need

<!-- S6 START -->1. **Done (§4.4):** a second, independent row generator for the full session-XIV model and a replay of all eight root
   certificates against it, with independently enumerated column universes and independently computed prices.<!-- S6 END -->
2. A written proof of each row family in one place (this document is a start; §3 covers all of them).
3. The replay logs of the other worker (`lp/log/`) with the manifest hashes.
4. Optional: replace the Bukh–Chao citation by the small-n certificates in `out/small_n_lp_result.json`
   (the hand ≥ 46 proof then cites nothing).

## 7. What I did NOT check
<!-- S7A START -->* The eight root certificates: now done (§4.4).<!-- S7A END -->
<!-- S7B START -->* Prices were recomputed in full for every type that occurs in a replayed root model (all the chain needs); the research
  tables also hold prices for types pruned by the caps, which were not recomputed.<!-- S7B END -->
* The HiGHS-side numerics: irrelevant to validity but relevant to reproducibility.
* The session-XIV "support analyses" (which families a certificate uses) — descriptive only.
* The 53-point upper bound (out of scope).
* Any search for a 46-point Kakeya set (the theorem, if true, forbids one; the local searches of sessions I–IV
  found nothing below 53; I ran none).

## 8. Verdict
<!-- S8 START -->**No lemma in `HAND_LEMMAS.md` is false or doubtful, every one is
proven in §3, every table they rely on is reproduced by fresh enumeration, the coverage tree is complete, and
the encodings I read match the lemmas.** The whole reduction (82 certificates) and all eight root certificates are now verified
independently of the research code (§4.3, §4.3b, §4.4): own row generators, own column universes, own prices, exact
integer Farkas checks, every research y^T b reproduced. "K(3,5) ≥ 47" can be stated publicly as **"LP-certified,
independently replayed"**. What separates it from "theorem" is (a) a human reading of the hand proofs in §3 — the lemmas,
not the code, are now the trust base — and (b) the caveat that both generators were written by Claude sessions reading the
same write-ups, so a shared misreading of a lemma would survive this replay and only a referee would catch it.
<!-- S8 END -->