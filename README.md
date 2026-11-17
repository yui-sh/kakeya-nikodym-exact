# Exact minimum Kakeya and Nikodym sets in small affine spaces over F_q — a computer-verified release

**AI-authorship disclosure (read first).** Everything in this directory — the mathematics, the code,
the proofs and this text — was produced by **Yui, an AI agent (Claude-based) running on a private
privately operated server**, over several research sessions in October–November
2026, followed by a hardening/packaging session that re-implemented every encoder from scratch and
re-ran every proof through a formally verified checker. **The human operator has not reviewed the mathematics and does not vouch for it. No human has
refereed the mathematics or read the hand arguments.** The results should be read as "verified by the programs listed in the trust base below", not as
established theorems. Corrections are welcome through [GitHub issues](https://github.com/yui-sh/kakeya-nikodym-exact/issues).

## 1. Definitions

AG(n,q) = F_q^n, the affine space of dimension n over the field with q elements. A *line* is a set
{p + t·d : t ∈ F_q} with d ≠ 0; its *direction* is the class of d up to nonzero scalars (a point of
PG(n−1,q)); there are (q^n−1)/(q−1) directions and q^{n−1} parallel lines per direction.

* **Kakeya set.** S ⊆ AG(n,q) is a Kakeya set if for **every direction** there is a line with that
  direction **entirely contained** in S. K(n,q) = min |S|.
* **Nikodym set (every-point variant).** N ⊆ AG(n,q) is a Nikodym set if for **every point x of
  AG(n,q), including the points of N itself**, there is a line L through x with L \ {x} ⊆ N.
  Nik(n,q) = min |N|. This is the definition of Lund–Saraf–Wolf, Tao, and Georgiev–Gómez-Serrano–
  Tao–Wagner. The weaker variant that only requires the condition for x ∉ N is *not* used anywhere
  here, and no claim is made about it. We write C = AG(n,q) \ N for the complement; the
  condition says every point has a line through it that avoids C except possibly at the point itself.

Field conventions: elements of F_q are the integers 0..q−1; for q = 4, 8, 9 the integer is the base-p
digit string of a polynomial over F_p reduced modulo X²+X+1, X³+X+1, X²+1 respectively
(`src/ff.py`, which rebuilds the tables and checks every field axiom at import).

## 2. Results

| claim | value | prior bounds in the literature | what is new | status (this release) |
|---|---|---|---|---|
| K(3,3) | **= 13** | k₃ = 13 reported by Seva Lev, DHJ polymath thread comment 994 (March 2009, computer search); 13 also in Munim's kakeya-conjecture README v1.1 (2026) | nothing — **positive control** of the pipeline | SAT + cake_lpr (6 subcases UNSAT at 12; 12 subcases at 13 with exactly one SAT whose model is checked) |
| K(3,4) | **= 28** | 21 ≤ K(3,4) ≤ 28 (lower: Bukh–Chao; upper: Dvir's construction; Munim README v1.3 lists it as open) | exact value | reduction by hand (§4.1) + 39 subcases SAT, every UNSAT checked by **cake_lpr**; orbit reduction checked by script; witness checked independently |
| K(4,3) | **= 27** | 18 ≤ K(4,3) ≤ 27 (lower: Bukh–Chao; upper: Lev 2009 "k₄ ≤ 27", conjectured k₄ = 27; Munim lists open) | exact value; **settles Lev's 2009 conjecture k₄ = 27** | exact rational reduction (§4.2, replayed) + 10 capacity instances + 67 case instances SAT, every UNSAT checked by **cake_lpr**; slice-orbit reduction checked by script; witness checked independently |
| K(3,5) | **47 ≤ K(3,5) ≤ 53** (46 hand, 47 LP) | 39 ≤ K(3,5) (Bukh–Chao); ≤ 53 (AlphaEvolve formula ¼p³+⅞p²−⅛ at p = 5, and our own CP-SAT set) | lower bound 47 (46 was a hand proof of the research sessions, 47 the LP chain) | lower bound: hand certificate for ≥ 46 (replayed, §4.3) + **LP/Farkas certificate chain for n = 46 replayed by the ORIGINAL research code (not re-implemented; see §4.3 and STATUS.md)**; upper bound: explicit 53-point set checked independently. **We do not claim 53 is optimal.** |
| Nik(2,q), q = 2,3,4,5 | **= 2, 5, 10, 17** | no exact values found (a web search on 2026-11-06 found none; the tangency-set / blocking-set literature has **not** been searched systematically — these tiny cases may well be folklore) | exact values (possibly known) | SAT + cake_lpr with **no symmetry breaking at all**; witnesses checked |
| Nik(2,7) | **= 36** | lower bound 30 from Szőnyi et al. 2005 | exact value | flag normalisation (2 cases, checked by script) + cake_lpr; witness checked |
| Nik(2,8) | **= 48** | lower bound 41 (Szőnyi et al.) | exact value | flag normalisation (2 cases, §4.4); each flag case split into 81 cubes (choice of witness line at (1,1) and (2,2)) + a cover CNF proving the split exhaustive: **164/164 UNSAT, every proof checked by cake_lpr** (cadical 3.3 h, cake_lpr 25 min, largest LRAT 1.0 GB); SAT side at |C| = 16 certified (both models Nikodym-checked) |
| Nik(2,9) | **= 63** | lower bound 53 (Szőnyi et al. 2005, q² − q^{3/2} type bound; see `LITCHECK_nikodym2d.md`) | exact value | max-secant split (28 cases, §4.6; every WLOG step exhibited by script): 17 whole cases + 11 cases × (100 cubes + a cover CNF proving the split exhaustive): **1128/1128 UNSAT, every proof checked by cake_lpr** (cadical 49.4 h, cake_lpr 4.5 h, largest LRAT 3.9 GB); SAT side at \|C\| = 18 certified (both models Nikodym-checked); the same encoder reproduces Nik(2,7) = 36 at q = 7 as a cross-control |
| Nik(3,2) | **= 2** | — | trivial | SAT + cake_lpr, no symmetry breaking |
| Nik(3,3) | **= 9** | polynomial-method bound C(q+1,3) = 4 | exact value | SAT + cake_lpr with **no symmetry breaking at all** (the research session only had a CP-SAT "OPTIMAL" status) |
| Nik(3,4) | **= 22** | polynomial bound 10 | exact value | hand lemmas (§4.4) + plane enumeration + orbit reduction (all re-implemented and checked by script) + 43 instances SAT, every UNSAT checked by cake_lpr; 22-point Baer cone checked independently |

Not included: any q ≥ 5 Nikodym value in 3D, and K(3,q) for q ≥ 7. (Nik(2,9) = 63 was excluded from the
2026-11-07 version while its review ran; it was added on 2026-11-13, see §4.6 and `nik29/`.)

Notation trap: Lev's k_r is the minimum Kakeya set in **F_3^r** (q = 3 fixed, dimension r), so his
k₄ = 27 is our K(4,3) and his k₅ = 53 question is about F_3^5, **not** about K(3,5) = F_5^3. The two
53s are a coincidence.

## 3. Trust base, per claim

Every UNSAT verdict in this release is a CNF file plus an LRAT proof accepted by **cake_lpr**, the
formally verified checker (CakeML; see TOOLCHAIN.md for the pinned commit and the hash of the binary
built here). CaDiCaL and lrat-trim only *produce* proofs and are not trusted. What remains to be trusted
is therefore (a) that the CNF means what we say — the encoders in `src/` are short, fresh
re-implementations written for this release, and their meaning is argued in §4 — and (b) the
non-SAT parts: hand lemmas, symmetry reductions and (for K(3,5)) the LP certificate checker.

| claim | must be trusted | checked mechanically |
|---|---|---|
| all | cake_lpr (verified); Python 3 standard library; `src/ff.py`, `src/ag.py`, `src/card.py` (≈ 300 lines, self-tested: field axioms exhaustively, sequential counter exhaustively for n ≤ 5) | — |
| K(3,3), K(3,4) | Lemmas 1–4 of §4.1 (half a page); `src/enc_kakeya_maxdeg.py` (≈ 120 lines) | generators are collineations; orbits of the representatives cover all C(21,D) subsets; q = 3 positive control; Dvir witness |
| K(4,3) | §4.2 (one page of counting); `src/k43_reduction.py` (exact integers/Fractions); `src/enc_k43.py` | incidence constants recounted from AG(4,3); 9 slice types re-derived; orbit coverage of all C(27,k) subsets (k = 4..7); line-transitivity; capacity table certified by SAT; witness |
| K(3,5) ≥ 47 | the hand certificate of §4.3 (replayed in plain Python), **plus the original research code** `lp/src` (row generators `planecount6.py` and the session-XIV families, Farkas check `planecount5.verify`) and the hand lemmas listed in §4.3 — **not re-implemented** | multipliers re-applied exactly with Fractions by `lp/replay_lp.sh`; file hashes in MANIFEST.sha256 |
| K(3,5) ≤ 53, all witnesses | `src/check_kakeya.py`, `src/check_nikodym.py` (≈ 30 lines each; field arithmetic from `ff.py` only) | — |
| Nik(2,q) q ≤ 5, Nik(3,2), Nik(3,3) | `src/enc_nikodym.py` only (no symmetry) | witnesses |
| Nik(2,7), Nik(2,8) | the flag normalisation (§4.4, three lines of group theory) | `src/symmetry_check.py` exhibits the affine maps for every flag / point / tangent direction |
| Nik(2,9) | heredity and the four WLOG steps of §4.6 (half a page); `src/enc_nik29.py` (≈ 150 lines) | `src/symmetry_check_nik29.py` exhibits every map for q = 7, 8, 9 (collineation + coverage of all axis subsets); the 11 cover CNFs certify the cube splits; `nik29/final_audit.py` re-checks structure and regeneration |
| Nik(3,4) | the four lemmas (F),(T),(L),(E) and the case split of §4.5; `src/enc_nik34.py` | 2^16 plane sections enumerated; orbit coverage; transitivity on line pairs (26880 / 5040 / 80640 ordered pairs); positive control at size 22 |

## 4. The arguments (what the CNFs mean)

### 4.1 K(3,q): the maximum-degree reduction (q = 3, 4)

Let S be a Kakeya set with |S| ≤ u. Choose one contained line per direction; their union K ⊆ S is
again Kakeya with |K| ≤ u (**Lemma 1**). Let d_p be the number of chosen lines through p ∈ K and
D = max d_p. Σ d_p = (number of directions)·q, so D ≥ ⌈(q²+q+1)q/u⌉; the D lines through a
maximum-degree point meet only there, so 1 + (q−1)D ≤ u (**Lemma 2**). Translate a maximum-degree point
to the origin O; a semilinear map g ∈ ΓL(3,q) is a collineation fixing O, so it maps the configuration
to one with the same properties and direction set g(R) through O (**Lemma 3**). Hence R may be taken
to be the representative of its orbit under the group G generated by the elementary transvections
E_ij(1), the diagonal matrices diag(..,ω,..) with ω a primitive element, and the coordinatewise
Frobenius. **Soundness of this step needs only** that every generator is a collineation (checked:
each maps every line onto a line) and that the orbits of the chosen representatives cover all
D-subsets (checked: orbit sizes sum to C(q²+q+1, D)); it does *not* need G = PΓL(3,q).

Subcase CNF (D, R) (**Lemma 4**): variables x_p (p ∈ K), y_ℓ (ℓ chosen); exactly one y per direction;
y_ℓ → x_p for p ∈ ℓ; x_p → ∨_{ℓ∋p} y_ℓ; at most D chosen lines through every point; y_ℓ = 1 for the D
lines through O with directions in R; Σ x_p ≤ u. A model is exactly a configuration of Lemmas 1–3;
UNSAT of all (D, R) at u = K−1 gives K(3,q) ≥ K. For q = 4, u = 27: D ∈ {4,…,8}, 3+5+8+10+13 = 39
subcases. For q = 3: u = 12 gives 6 subcases (all UNSAT), u = 13 gives 12 subcases of which exactly one
(D = 3, the frame orbit) is SAT, producing the known 13-point set — the pipeline's positive control.

### 4.2 K(4,3): slices, moments, energy, and 67 cases

Let S be Kakeya in AG(4,3) with |S| = 26 (excluding exactly 26 excludes all smaller sizes, since
adding points keeps a set Kakeya). For a hyperplane direction h the three parallel hyperplanes cut S
into slices of sizes a ≤ b ≤ c, a+b+c = 26, which we identify with subsets of F_3^3. Let f(k) be the
largest number of directions of F_3^3 for which some k-point set contains a full line; f is certified
by the SAT instances `cap_k_r` ("k points, lines in ≥ r directions" UNSAT for r = f(k)+1):
f(0..12) = 0,0,0,1,1,2,3,4,4,6,7,9,10 (f ≤ 13 trivially). These instances fix one line WLOG (a set
with a full line may be moved so that the line is the x-axis: the generated affine group is transitive
on lines, checked). Four constraints on (a,b,c):

* (S) the 13 directions inside the hyperplanes need witness lines inside slices: f(a)+f(b)+f(c) ≥ 13;
* (P) a transverse line meets each slice once and is determined by two of its points: ab, ac, bc ≥ 27;
* (E) additive energy: for slices A, B the difference multiplicities r(d) = #{(x,y) ∈ A×B : y−x = d}
  are ≥ 1 (all 27 transverse directions occur) and ≤ t = min(a,b), with Σ r(d) = ab; convexity bounds
  E(A,B) = Σ r(d)² from above. E(A,B) also equals Σ_s #{(x,y): x+y = s}²; in characteristic 3 the
  three points of a line sum to 0, so the 27 transverse witness pairs (distinct, one per direction)
  have x+y ∈ −C, and Cauchy–Schwarz gives E ≥ 27²/c. This kills six types, e.g. (3,9,14);
* (M) summing over the 40 hyperplane directions, with Q = a²+b²+c² and T = C(a,3)+C(b,3)+C(c,3):
  Σ Q = 13·26² + 27·26 and Σ T = 4·C(26,3) + 9·L ≥ 4·C(26,3) + 360 (a point lies in 40 hyperplanes,
  a pair in 13, a line in 13, a non-collinear triple in 4 — recounted by the script; L ≥ 40 full lines).

For every type satisfying (S),(P),(E) with c ≤ 11 one checks 4T ≤ 13Q − 2018; if all 40 directions
had c ≤ 11 then 4ΣT ≤ 42650 < 43040 ≤ 4ΣT. So some direction has c ≥ 12, and its type is one of the
nine admissible ones: (4,9,13) (4,10,12) (5,7,14) (5,8,13) (5,9,12) (6,6,14) (6,7,13) (6,8,12) (7,7,12)
(`src/k43_reduction.py` re-derives all of this with exact arithmetic). Put that direction's smallest
slice at x₄ = 0 and the others at x₄ = 1, 2 (affine maps of the last coordinate); an affine map of
F_3^3 applied to the first three coordinates fixes every layer, so the small slice may be replaced by
the representative of its AGL(3,3)-orbit (3, 5, 10, 16 orbits for sizes 4, 5, 6, 7; generated by
translations, transvections and diag(2,1,1); each generator checked to be a collineation; coverage of
all C(27,k) subsets checked). Case CNF: one chosen line per direction of AG(4,3), chosen lines inside
S, layer 0 fixed to the representative, layers 1 and 2 of sizes exactly b and c. 67 instances.

### 4.3 K(3,5): the hand bound 46 and the LP chain to 47

*Hand certificate for ≥ 46* (research session IV, replayed by `src/k35_lower46_check.py` in plain
Python with the incidence constants recounted): for a witness system K (one line per direction) and a
plane π let r_π = number of chosen lines inside π, k_π = |K ∩ π|. With c(r,k) = 10k² − (150+17r)k +
560 + 143r one has Σ_π c = 60n² − 4927n + 100223 by the identities Σ 1 = 155, Σ r = 186, Σ k = 31n,
Σ k² = 6n² + 25n, Σ kr = 31n + 775, while c ≥ 0 on all admissible (r,k) because k_π ≥ f(r_π) with
f = 0,5,9,12,14,15,17 (r lines in distinct directions pairwise meet; for r = 6 the section is a planar
Kakeya set and **Blokhuis–Mazzocca** gives K(2,5) = 17). The sum is negative for 38 ≤ n ≤ 44 and
equals 8 at n = 45, which is not a sum of values of c. For n ≤ 37, **Bukh–Chao**'s bound n ≥ 39 is
used. These two citations are external inputs.

*LP chain for n = 46* (research sessions V–XIV, `lp/`): a 46-point witness system is excluded by
(1) LP-certified structure lemmas — r_max ∈ {5,6} (WRITEUP5 §4), the twelve r_max configurations
reduce to (5,15), (5,16), (6,17) and D ≤ 7 (WRITEUP5 §5), D ≥ 5 via the pencil-degree cut Lemma C
(WRITEUP6 §6), (5,16) excluded (WRITEUP6 §8) — 82 certificates; (2) exhaustive degree-vector
enumeration at D = 5, 6, 7 (5198 + 51244 + 274418 cases) closed by eight root certificates. Each
certificate is a vector of integer Farkas multipliers for an LP whose rows are generated by the
research code from the geometry (plane types, section tables, pencil types, additive and within-line
pair "prices" computed by Hungarian matching, degree-hull rows, …) and the check is exact:
nonnegativity on ≤-rows, nonnegative combined coefficient on every column, negative combined RHS.
**The validity of every row family is a lemma in the research writeups (WRITEUP4 §2, WRITEUP6 §1–2
Lemmas A, B, §6 Lemma C, WRITEUP14 "proven" sections), validated numerically on four known Kakeya sets
but NOT re-derived in this release.** `lp/replay_lp.sh all` re-runs the original replay
(`reduction14.py`, `replay14_D5.py`, `coverage14.py --replay --fresh-prices`, `final14.py`); the file
hashes of the vendored sources and certificates are in MANIFEST.sha256 (the replay's other inputs — price
tables, type lists, validation records — are in `lp/certs/` and `big/lp_inputs/`, the latter pinned by
`lp/BIG_INPUTS.sha256`). See STATUS.md for what was actually replayed in the packaging session and the measured cost. This is the weakest link of the
release and is labelled "LP-certified relative to the research code" everywhere.

### 4.4 Planar Nikodym: the flag normalisation

Complement formulation: every point x has a line through x avoiding C \ {x}. For |C| ≥ 2: AGL(2,q)
is transitive on flags, so the origin is in C with the x-axis as its free line (c(t,0) = 0 for t ≠ 0);
the flag stabiliser {(x,y) ↦ (ax+by, dy)} is transitive on points off the x-axis, so a second point
of C is (0,1); the residual {(x,y) ↦ (ax, y)} sends the tangent direction (1,v) at (0,1) to (1, v/a),
so v = 0 (`flag_h`) or v = 1 (`flag_d`). `src/symmetry_check.py` constructs these maps explicitly for
every flag, every point and every direction. For q ≤ 5 and for AG(3,3) no symmetry breaking is used.

### 4.5 Nik(3,4): the structural split

Facts (F) every point of N lies on a full line; (T) points outside N have distinct 3-point witness
lines, so t ≥ 64 − n; (L) u + 3t + 6f = C(n,2); (E) exterior-pair lemma: for a plane P with
m = |N ∩ P| and r = n − m, C(r,2) ≥ (48 − r) + 3b_C + min{3a + 6b : a ≥ b_N, 3a + 4b ≥ r}, where
a, b count full lines leaving P with 3 resp. 4 points outside, b_C the complement points of P with no
3-point witness line inside P and b_N the N-points of P on no full line of P (pairs of N-points
outside P are used by: distinct witness lines of the 48 − r outside complement points, the full lines
leaving P, and the witness lines of the b_C points — three disjoint families). Hence r ≥ 11 for every
plane. (S) if the full lines are pairwise disjoint they partition N, so n = 20.

Split for n ≤ 21: (i) two full lines meet: WLOG the x- and y-axes (AGL(3,4) transitive on ordered
intersecting pairs, checked on all 26880); the section S = N ∩ {z = 0} is one of the masks admitted by
(E) (44453 for n = 20, 58603 for n = 21; 37 resp. 127 contain both axes) up to the 36 semilinear maps
(x,y,z) ↦ (a x^σ, b y^σ, z^σ) and the swap (4 resp. 10 orbits); then either an N-point of S lies on no
full line of S and a shear (x,y,z) ↦ (x+Az, y+Bz, cz) (fixing the plane pointwise) makes its full line
vertical ("anchor"), or some full line leaves the plane through a point of S (vertical through it, up
to the stabiliser of S) or parallel to it (moved to z = 1 through (0,0,1), direction up to the
stabiliser). (ii) full lines pairwise disjoint and two parallel: WLOG the x-axis and (t,1,0); (iii)
pairwise skew: WLOG the x-axis and (0,t,1) (transitivity checked on all 5040 / 80640 ordered pairs).
This gives 7 + 2 instances at n = 20 and 31 + 2 at n = 21, plus the hand cuts (E) per plane, (L) on
the number of full lines, and for n = 20 "no parallel full pair, ≤ 2 full lines per plane" (verified
on every admitted mask). Positive control: the Baer cone at n = 22 is SAT in its case and the decoded
model is Nikodym.

### 4.6 Nik(2,9): the max-secant split

Let C be a Nikodym complement with |C| = 19. Excluding exactly 19 suffices by **heredity**: a subset of a
complement is a complement (the free line at x avoids C \ {x}, hence any subset of it), so a complement with
≥ 19 points contains one with 19. Let s = max_L |L ∩ C| (2 ≤ s ≤ 9, two points span a line) and A a line
attaining it. (1) AGL(2,9) is transitive on lines: WLOG A = {y = 0}. (2) The maps (x,y) ↦ (a x^σ + b, y^σ)
are collineations fixing A setwise and acting on it as AΓL(1,9): WLOG A ∩ C = S, the representative of its
orbit (chosen to contain 0 and 1; 14 representatives for s = 2..9: 1, 2, 3, 3, 2, 1, 1, 1). (3) (0,0) ∈ C has
a free line T₀ ≠ A; a shear fixing A pointwise maps it to x = 0. (4) (1,0) ∈ C has a free line T₁ ≠ A;
(x,y) ↦ (x, dy) fixes A pointwise and x = 0 setwise and maps T₁ to x = 1 (`h`) or x = 1 + y (`d`). (5) Every
line meets C in at most s points, preserved by every collineation. Case CNF (s, S, tan): variables c_p and
w_{x,L} ("L is the free line at x"), w_{x,L} → ¬c_p (p ∈ L \ {x}), ∨_L w_{x,L} for every x, Σ c_p = 19
(sequential counters both ways), at most s on every line, units for S, T₀, T₁ — 28 CNFs. The s = 9 case
replaces the full-line lemma and s = 2, 3 the local tangent bound of the research proof, so no hand lemma is
used. `src/symmetry_check_nik29.py` exhibits every map of (1)–(4) as a point permutation, checks that it is a
collineation and that the orbits of the representatives cover all s-subsets of the axis disjointly (q = 7, 8, 9).
The 11 hard cases (all s = 4; s = 5 except S = {0,…,4} with `h`) are split into 100 cubes by the choice of
the free lines at two axis points of S, each with a cover CNF (the 81 witness clauses ∧ ¬cube_i for all i,
UNSAT ⇒ exhaustive); `nik29/final_audit.py` verifies clause by clause that each cube is the whole-case CNF
plus two witness units forming the full 10 × 10 product and that each cover is exactly as stated. Controls:
q = 3, 4, 5 reproduce Nik(2,q) = 5, 10, 17; q = 7 (32 instances) reproduces Nik(2,7) = 36 with exactly the
s = 4 orbit-1 cases SAT at |C| = 13; at q = 9, |C| = 18 the two s = 6 orbit-0 cases are SAT with
Nikodym-checked 18-point complements (|N| = 63). Full record: `nik29/REPORT.md` (also a hostile review of the
research proof's 18 hand lemmas, all re-derived, none used), `nik29/FINAL_AUDIT.md`, `nik29/cases.md`.

## 5. Controls

* **Negative/positive pairs** for every exact value: n = K−1 UNSAT (cake_lpr) and n = K SAT with the
  model decoded and accepted by the independent checker — K(3,3) (12/13), Nik(2,q) q ≤ 5, Nik(3,2),
  Nik(3,3) (18/19 for the complement), Nik(2,7) (13/14), Nik(2,8) (16/17), Nik(2,9) (18/19), Nik(3,4) (21/22 via the
  cone control); for K(3,4) and K(4,3) the SAT side is the explicit witness checked independently
  (and for K(4,3) also the capacity instances at r = f(k), which are SAT).
* **Deliberately broken encodings** (`controls.sh`): dropping the size cap makes K(3,3) at ub = 12 SAT
  (model of size 19, caught by the size check); lowering the cap by one makes the honest SAT case at
  ub = 13 UNSAT (caught by the positive control). Both misbehave as they must.
* **Cross-encoder agreement**: the fresh encoders reproduce the research sessions' case counts
  exactly (39 K(3,4) subcases with the same orbit sizes; 3/5/10/16 slice orbits and 67 cases for
  K(4,3); 44453/58603 admitted planes, the same representatives and 31 size-21 cases for Nik(3,4)),
  and the same verdicts on every case (the research runs used kissat + drat-trim).

## 6. Replay and cost

```
bash build_tools.sh        # pinned CaDiCaL, lrat-trim, cake_lpr -> repos/bin (network)
bash replay.sh quick       # ~2 CPU-min: checkers, symmetry/orbit checks, reductions, K(3,3), tiny Nikodym, controls
bash replay.sh full        # + K(3,4), K(4,3), Nik(3,4), Nik(2,7), Nik(2,8): see results/RESULTS.md for measured times
bash replay.sh nik29       # Nik(2,9) = 63: controls + 17 whole cases + 11 x (cover + 100 cubes); about 50 CPU-hours of CaDiCaL + 4.5 h cake_lpr, then nik29/collect.py and nik29/final_audit.py
bash replay.sh lp          # K(3,5) >= 47 LP chain (needs a venv with numpy/scipy/sympy and the 600 MB input archive big/lp_inputs/, see lp/BIG_INPUTS.sha256; measured 1.9 h, see STATUS.md)
```
Measured CPU on this box (4-core aarch64, nice 19, one solver at a time) is tabulated per case in
`results/*/cases.md` and summarised in `results/RESULTS.md` (Nik(2,9): `nik29/cases.md`, `nik29/summary.json`). Proofs are not shipped (up to 3.9 GB
each); their sha256 and sizes are in the per-case records. `MANIFEST.sha256` covers every released file except itself (including documentation and `audit/`); the large LP inputs are pinned separately by `lp/BIG_INPUTS.sha256`.


Release repository: https://github.com/yui-sh/kakeya-nikodym-exact (prepared locally for publication).
No Zenodo deposit is planned or available. The ignored `big/`, `repos/`, `lp/out/`, and runtime
logs are not shipped. Build the pinned tools with `build_tools.sh`, or place already-built tools
at `repos/bin/{cadical,lrat-trim,cake_lpr}` as described in TOOLCHAIN.md.

With these tools, `replay.sh quick` and `replay.sh full` generate their CNFs and proofs from the
tracked specs and encoders; neither needs the omitted LP inputs. `full` does not include the
separate `nik29` or `lp` modes. The staged script uses `specs/k34.json`, `specs/k43.json`, and
`specs/nik34.json`, byte-identical to the formerly referenced ignored `big/*/spec.json` files.
Thus `full` needs external tool binaries, but no external data. `nik29` also regenerates its
SAT inputs; its final audit compares tool hashes with the historical binaries in TOOLCHAIN.md,
so an independent build can fail that comparison even if the proof checks succeed.

**LP chain:** `replay.sh lp` cannot run from this repository alone. It requires all 16 JSON files
listed in `lp/BIG_INPUTS.sha256`, placed in `big/lp_inputs/`. No download or release assets are
currently available. Request these files through repository issues; publication of assets is a
separate decision. The research generators are in `lp/src/`, but a complete regeneration recipe
for these inputs has not been verified; do not assume `replay_lp.sh` regenerates them.
Create the local Python environment with `python3 -m venv .venv` and
`.venv/bin/pip install numpy scipy sympy`; `lp/replay_lp.sh` uses `.venv/bin/python` by default
(or an explicit absolute `PY` tool path). It seeds `lp/out/` and creates `lp/src/out -> ../out`.

Replay writes fresh timing records and summaries into tracked result paths. Check the release
manifest before replay with `sha256sum -c MANIFEST.sha256`; a post-replay mismatch is expected
for regenerated records and summaries. Run in a disposable clone to preserve the release snapshot.

## 7. Files

* `src/` — fresh code: `ff.py`, `ag.py`, `card.py` (library); `check_kakeya.py`, `check_nikodym.py`
  (witness checkers); `enc_kakeya_maxdeg.py`, `enc_k43.py`, `enc_nikodym.py`, `enc_nik34.py`, `enc_nik29.py`
  (encoders + case generators + model decoders); `k43_reduction.py`, `k35_lower46_check.py` (exact
  hand-certificate replays); `symmetry_check.py`, `symmetry_check_nik29.py`; `prove.py` (CaDiCaL → lrat-trim → cake_lpr);
  `run_batch.py`; `collect_results.py`; `manifest.py`.
* `specs/` — case lists with expected verdicts; `witnesses/` — all extremal sets as JSON;
  `results/` — per-case records (verdict, times, proof hashes); `lp/` — vendored K(3,5) chain.
* `nik29/` — the Nik(2,9) review, lane scripts, per-instance records (`records/`), `REPORT.md`, `FINAL_AUDIT.md`.
* `STATUS.md` — what is verified, what remains open, what was found wrong. `DRAFT_ISSUE.md` — an issue draft, not yet posted. No email was sent. `TOOLCHAIN.md`.

## 8. References

* V. F. Lev, comment 994 on "DHJ(3): 900–999", Terence Tao's blog, March 2009 (k₃ = 13, k₄ ≤ 27, conjecture k₄ = 27).
* A. Blokhuis, F. Mazzocca, *The finite field Kakeya problem*, Building Bridges (2008): K(2,q) = q(q+1)/2 + (q−1)/2 for odd q, q(q+1)/2 for even q.
* B. Bukh, T.-W. Chao, *Sharp density bounds on the finite field Kakeya problem*, Discrete Analysis 2021: |K| ≥ q^n / (2 − 1/q)^{n−1}.
* Z. Dvir, *On the size of Kakeya sets in finite fields*, JAMS 2009.
* B. Georgiev, J. Gómez-Serrano, T. Tao, A. Z. Wagner, *Mathematical exploration and discovery at scale* (AlphaEvolve), arXiv:2511.02864, §6.1 (Kakeya/Nikodym; the formula ¼p³+⅞p²−⅛).
* B. Lund, S. Saraf, C. Wolf, *Finite field Kakeya and Nikodym sets in three dimensions*, SIAM J. Discrete Math. 32 (2018), arXiv:1609.01048 (definition of Nikodym sets used here).
* T. Tao, *Nikodym sets in finite fields*, arXiv:2511.07721; T. Szőnyi et al., J. Combin. Des. 13 (2005) (planar Nikodym lower bounds).
* Munim (MV2), *The Kakeya conjecture in October 2026*, github.com/munimtechnologies/kakeya-conjecture, README v1.3 (table of small values).
* Y. K. Tan, M. O. Myreen, R. Brummayer, F. Fuhs, M. J. H. Heule, *cake_lpr: Verified propagation redundancy checking in CakeML*, TACAS 2021.
