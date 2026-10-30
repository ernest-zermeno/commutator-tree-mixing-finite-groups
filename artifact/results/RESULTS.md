# Collected results

The original validation environment was Python 3.9.6, NumPy 2.0.2; those logs are
recorded in `notes/strengthening-2026-09-23/`. The 26 September additions are
recorded in `notes/open-problems-implementation-2026-09-26/`, including the
environment and the selected checks from a fresh artifact extraction. The 16
October additions (intersection-six stratum, `octahedral_cap6_check.py`, SymPy
1.14.0 added to the requirements) are recorded in
`notes/octahedral-cap6-implementation-2026-10-16-1950/`, including the full
26-command suite run from a fresh artifact extraction with every command
exiting 0 (Python 3.9.6, NumPy 2.0.2, SymPy 1.14.0).
"op-norm" = operator norm of the level-k commutator-tree Fourier operator;
"bottleneck" = 2nd singular value of the regular-rep operator = max op-norm over all
non-trivial irreducibles (so < 1 ⟺ strict contraction in all non-trivial irreps).

**Certification scope:** exact group counts and derived series use integer arithmetic.
SVD norms use float64 and tolerances, not interval arithmetic. Unless explicitly labelled
ORDERED (A5 and `verification_groups.py`), level-2 sweeps use reduced representatives;
validity and leaf closure survive the reductions, but worst norms need not. Level-3
regular-representation sweeps are capped/sampled. The external P0 certificate separately
covers all 9,447,840 oriented S5 level-3 configurations in 39,591 orbits.

## A. The derived-length principle (`derived_length_test.py`, `gl3f2.py`, `mixing_level_pgl.py`)

**This table's "level" column is the augmentation-rep contraction level.** Full contraction (all
irreps) is the stronger ⟨leaves⟩=G criterion of §A2 — for S₅/6pt it is NOT reached at
level 2 with these leaves. Level-3 regular-rep figures are sampled/capped (§A2, C).

| action | stabilizer | dl | level-2 op-norm=1 | aug-rep level | matches dl? |
|---|---|:--:|---|:--:|:--:|
| S₅ / 5 pts (transpositions) | S₄ | 3 | 270/405 | **3** (safe full-oriented bound 0.0355) | ✓ |
| S₅ / 6 pts = PGL₂(𝔽₅) (involutions) | F₂₀ | 2 | 0/27420 | **2** *(augmentation rep only — NOT full; see §A2)* | ✓ |
| S₆ / 6 pts (transpositions) | S₅ (non-solv.) | ∞ | 1620/1620 | **never** (λ=1.0) | ✓ |
| GL₃(𝔽₂) / 7 pts (involutions) | S₄ | 3 | 2772/13146 | **3** (L1/L2 all-irreps exhaustive; L3 **capped 60k**, =0.108) | ✓ |
| PGL₂(𝔽₇) / 8 pts (split, sign−1 leaves) | Borel | 2 | 0 | **2** (all irreps, **SAMPLED**: L2 0.772, L3 0.092) | ✓ |

S₅/5pt reproduces the reduced enumeration to the digit (30 unordered pairs, 405 quads,
270 op-norm-1, 73810 reduced level-3 configs, λ = 0.0353047532). The corrected full oriented
enumeration has max λ≈0.03540023208379892, so quote 0.0355 as a safe bound.

## A2. The real criterion: `2nd-sv < 1 ⟺ ⟨leaves⟩ = G` (`generates_G.py`, exhaustive L2)

| group | configs | std op-norm=1 (intransitive) | regular 2nd-sv=1 (⟨leaves⟩≠G) | generate G | std MISSES (proper-transitive) | mismatches |
|---|--:|--:|--:|--:|--:|:--:|
| GL₃(𝔽₂) | 13146 | 2772 | 2772 | 10374 | 0 | **0** |
| S₅ on 6 pts = PGL₂(𝔽₅) | 27420 | **0** | **4695** | 22725 | **4695** | **0** |

This closes finding #1: the full criterion is ⟨leaves⟩=G (0 mismatches with regular-rep
2nd-sv<1). Transitivity in the ambient transitive action is necessary, not sufficient, for full contraction: in S₅/6pt
it reports 0 obstructions at level 2, but **4695 configs have leaves in a proper
*transitive* subgroup** — verified split: **3705 close to A₅ (order 60) and 990 to an
S₄-type subgroup (order 24)** (6 distinct proper closures = 1 A₅ + 5 conjugate S₄;
all transitive on the 6 points, so the augmentation rep contracts but other irreps — the
sign rep for the A₅ ones — do not). In GL₃(𝔽₂) the std and full criteria coincide
because its only transitive maximal subgroup (7:3) has no involutions.

## A3. A₅ — the exhaustive check behind the paper's cor:A5 (`a5_criterion.py`)

| group | proper **nontrivial** C-generated closures (EXHAUSTIVE lattice, nonempty subsets) | max dl | L1 universe (ORDERED / REDUCED) | L1: 2nd-sv=1 | L2 universe (ORDERED / REDUCED) | generate A₅ | L2: 2nd-sv=1 | L2 bottleneck | exact level |
|---|---|:--:|:--:|--:|:--:|--:|--:|--:|:--:|
| A₅ (involutions) | C₂ (15) · V₄ (5) · S₃ (10) · D₁₀ (6) | **2** | 180 / 90 | **180/180** | 29,640 / 3,705 | 29,640/29,640 | **0** | 0.346367621209230 ± 1e-9 | **2** |

**All universes EXHAUSTIVE; SVD values NUMERICAL** (regular rep, so all irreps at once; the closure lattice covers
⟨S⟩ for *every nonempty* subset S of the 15 involutions — 37 closures = the 36 proper
nontrivial ones above + A₅ itself; the trivial closure of the *empty* subset is omitted
because dl(1) = 0 and it cannot affect the max D). Counts are stated for both universes:
**ORDERED** = ordered tuples of leaves with every gate valid (α ≠ e); **REDUCED** =
representatives of the free orbits of the gate symmetries (swap inside a gate, swap of
gates: orbit size 2 at L1, size 8 at L2; freeness asserted *per representative* by the
certificate, so 180 = 2×90 and 29,640 = 8×3,705). The operator sweep (2nd singular value,
generation) runs over the full ORDERED universe — reduced representatives are never
presented as the universe. A₄ is **not** C-generated — its involutions close to V₄ — so it
does not enter the max, and every proper closure is solvable with dl ≤ 2.
This certifies cor:A5: D = 2 by the dichotomy, level 1 never contracts (180/180 at 2nd-sv 1),
level 2 always does (bottleneck 0.3464 < 1, value pinned with justified tolerance) ⇒ exact
level 2. `a5_criterion.py` exits nonzero if any of these invariants fails.

## B. Geometric Lemma — `op-norm(D̂₂)=1 ⟺ ⟨leaves⟩ intransitive` (augmentation rep) (`geometric_lemma.py`)

| action | quads | op-norm=1 | intransitive | mismatches |
|---|--:|--:|--:|:--:|
| S₅ / 5 pts | 405 | 270 | 270 | **0** |
| S₅ / 6 pts = PGL₂(𝔽₅) | 27420 | 0 | 0 | **0** |
| GL₃(𝔽₂) / 7 pts | 13146 | 2772 | 2772 | **0** |
| PGL₂(𝔽₇) / 8 pts (split) | 53718 | 0 | 0 | **0** |

Refinement: of GL₃(𝔽₂)'s 2772 op-norm-1 configs, only **1386 fix a common point** — the
rest fix a *block*. "Intransitive" is strictly weaker than "fixes a point"; the S₅
simplex was the special case.

## C. Sign / determinant subtlety (`bottleneck_pgl7.py`, `mixing_level_pgl.py`)

PGL₂(𝔽₇) involution classes: {fix 0 pts → even → sign +1 : 21}, {fix 2 pts → odd →
sign −1 : 28}. With arbitrary involution leaves a 1-dim (sign) rep stays at op-norm 1
(44/1700 sampled configs) — the walk is trapped in PSL₂. With the correct **sign −1
(split) leaves** the sign rep vanishes and the regular-rep bottleneck is
L1=1.0, **L2=0.772 (<1)**, L3=0.092. Exact level 2 follows from the analytic rank-one theorem; these norms are **SAMPLED** (~800 quads, L3
capped 4000): strong evidence, not a certificate.

## D. Composite modulus (`gl2z4.py`)

`test_strengthening.py` independently checks the full derived series; `gl2z4.py` alone only sweeps commutator trees. GL₂(ℤ₄): |G|=96, derived series **96→24→4→1 (solvable, dl 3)**; G″ order 4, abelian ⇒
level-2 targets commute ⇒ **0 valid level-3 configs**. q=4 is degenerate: base field 𝔽₂
has solvable GL₂(𝔽₂)=S₃ ⇒ solvable ambient ⇒ Barrington's commutator tree collapses.
(Same for GL₂(ℤ₉).) First non-degenerate composite case GL₂(ℤ₂₅), order 300000 — out of
direct reach.

## E. Robustness to non-involution leaves + determinant layer (`noninvolution_leaves.py`)

**(b1) Determinant (GL₂(𝔽₅)).** det(involutions) = {1,4} = squares; ⟨{1,4}⟩ = {1,4} (index
2 in 𝔽₅^*) ⇒ involution leaves trap the output determinant in the squares ⇒ NOT uniform
on GL₂(𝔽₅). Adding any leaf with det = 2: ⟨det⟩ = {1,2,3,4} = 𝔽₅^* ⇒ determinant mixes.
Since every involution has det ∈ {±1}, full GL_n(𝔽_q) mixing for q−1 > 2 **requires
non-involution leaves**.

**(b2) Non-involution leaves:** only the first 14 of 98 candidates in GL₃(𝔽₂) are used. The set is not asserted conjugation-closed, and generation is not checked. Norms concern the augmentation representation.
L1: 22/91 op-norm=1 (max 1.0); L2: 21/3994 op-norm=1 (max 1.0); **L3: max op-norm 0.011**
over the sampled/capped level-3 configurations. This does not establish an exact level for this selected set, nor an all-irreducibles bound. The analytic fixed-space criterion separately permits non-involution leaves.

## F. Rigid templates — realizability and target-transition certificates (paper §5.1)

**Class tables and finite censuses in this section are exact; mixed-transfer identities are exhaustive on S5 and checked on the first three conjugates per class elsewhere.**
pinned constants are re-verified as lexicographically-first on every run. The Phase-1
sweep infrastructure (`bridge_*_sweep.py`, pinned censuses, `test_bridge_sweeps.py`)
remains DATA — reproducible counts with counterexample records — and is **not** to be
reinterpreted as general theorems.

### F1. `bridge_realizability_certificate.py` — {2,3} realizability over S₅ (60 checks, exit 0)

| item | value | status |
|---|---|---|
| level-2 leaf universe (literal traversal) | 10⁴ = 3,600 valid-gate + 6,400 gate-invalid | EXHAUSTIVE |
| valid census | 3,600 = **1,080 generating** + **2,160 trapped** + **360 root-invalid** | EXHAUSTIVE |
| target dichotomy | all generating → 3-cycle target; all trapped → (2,2), H = full S₄ | CERTIFIED |
| reduced orbits (order-8 gate symmetry, all free) | **135** generating / **270** trapped | CERTIFIED |
| exact non-vacuous level 2 family | seed γ₂, direct 3-cycle loop, fixed pair (4 5), (1 2)(3 5); levels 2–6 | CERTIFIED |
| exact non-vacuous level 3 family | seed τ₃ = [τ_A,τ_B] (S₄-trapped children), 5-cycle loop, fixed pair (3 4 5), (2 3 4); levels 3–6 | CERTIFIED |

### F2. `target_transition_graph.py` — T(G) for S₅, GL₃(𝔽₂), PGL₂(7) (29 checks, exit 0)

| item | value | status |
|---|---|---|
| full class digraphs | S₅ (6 classes), GL₃(𝔽₂) (5), PGL₂(7) (8) — adjacency + loop pairs pinned | EXHAUSTIVE |
| S₅ facts | succ(2)={(3)}; succ(3)={(3),(2,2)} (no direct (3)→(5)); loops exactly at (3),(5); (2,2) has NO loop | CERTIFIED |
| GL₃(𝔽₂) seed | transvection seed (t1,t2,t1,t3): gates valid, root = t1, ⟨leaves⟩ = 168; loop pair of orders (4,4) | CERTIFIED |
| GL₃(𝔽₂) level-1 census | 336 valid ordered transvection pairs = 168 (order 6) + 168 (order 8); all proper | EXHAUSTIVE |
| PGL₂(7) parities | class of size 21: square dets (interior); size 28: non-square dets (exterior) | CERTIFIED |
| scope | small-group data; **not** a uniform theorem over large families | — |

### F3. `binary_target_grammar_certificate.py` — diagonal reduction over nine groups (694,181 checks, exit 0)

Verifies exhaustively the class-level dynamics, the solvable-radical characterisation, the
transfer identities and an end-to-end example — the ingredients of the paper's
diagonal-reduction block (`lem:mixedtransfer`, `thm:diagred`, `cor:binaryred`). The exact-D
corollary `cor:binarylevel` is established analytically; this script supports and
cross-checks its ingredients. Two audited runs, approximately 8 s.

| item | value | status |
|---|---|---|
| groups | S₅, A₅, S₆, GL₃(𝔽₂), PGL₂(7), PSL₂(11), SL₂(5), SL₂(3) (solvable control), AGL₃(2) | — |
| diagonal-slice regression | R(C,C;·) = def:targetgraph sweep {[w,hwh⁻¹]}, every class, all 9 groups | CERTIFIED |
| mixed transfer identities | constructed (g,h) verified; EXHAUSTIVE over (u,v,z′) on S₅, first three conjugates per class elsewhere | CERTIFIED |
| solvable radical | O∞(G) = ⟨∪{C : ⟨C⟩ solvable}⟩ solvable; ⟨C⟩ solvable ⟺ C ⊆ O∞(G), every class | CERTIFIED |
| main equivalence | binary persistence ⟺ some ⟨C⟩ non-solvable ⟺ diagonal persistence, EVERY non-empty subset of Ω(G) (up to 2¹⁰−1) | EXHAUSTIVE |
| cycle characterisation | recurrent core of T(G) ⟺ non-solvable closure ⟺ diagonal singleton dynamics (three formulations, shared tables) | CERTIFIED |
| survival census | S₅ 63/63 · A₅ 15/15 · S₆ 1023/1023 · GL₃(𝔽₂) 31/31 · PGL₂(7) 255/255 · PSL₂(11) 127/127 · SL₂(5) 254/255 (central class dies) · AGL₃(2) 1022/1023 (translation class dies) · SL₂(3) 0/63 | EXHAUSTIVE |
| S₅ pair facts | 36 ordered nontrivial class pairs; 33 multi-output; mono pairs ((2),(2)),((2),(3)),((3),(2)); ((2,2),(2,2))→{(3),(5)} | EXHAUSTIVE |
| end-to-end family | 5-level S₅ family, genuinely mixed production at every level, all members valid, targets on schedule, leaves in (2) | CERTIFIED |
| essential mixing | NONE in any group — every surviving portfolio is diagonally sustainable (the theorem explains why: impossible in every finite group) | CERTIFIED |

## G. Independent ordered census and new checks

`verification_groups.py` enumerates all ordered level-2 tuples and all C-generated
subgroups independently of the numerical scripts, with pinned exact counts:

| Group / leaves | Ordered valid level-2 | Proper leaf subgroup |
|---|---:|---:|
| S5 / transpositions | 3,240 | 2,160 |
| A5 / involutions | 29,640 | 0 |
| GL3(F2) / involutions | 105,168 | 22,176 |
| PGL2(5) / all involutions | 219,360 | 37,560 |
| PGL2(7) / exterior | 429,744 | 0 |

`test_strengthening.py` checks extension-field arithmetic at q=9,25,27,125,
subfield square-character parity, the q=5 exception, the non-dihedral affine
3³:C2 of order 54, determinant-based class selection, rejection of composite
moduli by the prime-only numerical script, both spectral comparison bounds,
the A5×C2 versus A5×V4 boundary, GL2(Z4)'s derived series, and nonzero failure
behaviour under injected norm discrepancies. Analytic formulas are distinguished
from independent group checks in its comments.

`psl33_maximals_check.py` checks |PSL3(3)|=5616, 117 involutions generating it,
a parabolic of order 432 and derived length 5 generated by involutions, and the
maximal S4 of order 24 and derived length 3. Completeness of the maximal list is
an external ATLAS input, not inferred from enumeration of these examples.

`pgl2_classification.py` uses exact finite-field arithmetic to search sampled
subgroup witnesses. A found non-solvable proper subgroup certifies an obstruction;
not finding one does not certify its absence. Printed exact levels are predictions
of the analytic theorem. No ambient PGL2(125) sweep is claimed.

`external/p0/20_level3_exact_certificate.py` is copied unchanged from P0, with
source commit and SHA-256 in its provenance. It certifies the bound 71/2000=0.0355
(and a tighter rational bound) for the P0 S5 calculation. P1 cites it as an input.

## Uniform family rates (23 September 2026)

The analytic diameter-to-mixing theorem and support lower bound give optimal
order Theta(log log |G|) of fixed-error depth for the families specified in
the manuscript. The generalized stability estimate allows non-involutions.
`python3 scripts/test_quantitative_families.py` checks small groups with leaf
orders 2, 3 and 4, and fixed-space obstructions. Its numerical operator checks
corroborate the proof; they do not certify an infinite-family statement.

## Natural octahedral pair (23 September 2026)

`octahedral_pair_check.py` verifies exact finite algebra for q=5,11,13,19:
A order 24, intersection order 6, generation of G, the group-algebra compression
identity, and well-defined connected simple four-regular Hecke graphs.
At q=5 it certifies the distinction between fixed-left and Hecke spectra
(4,3,3,3,−1 versus 4,−1,−1,−1,−1) by rational ranks.
The nonconstant normalized Hecke radii at q=5,11,13,19 are numerically
0.25, 0.6830127019, 0.75, 0.8090169944. These are not certificates for the
infinite-family bound. Generic almost-projection gluing stress checks are numerical.

`s4_internal_gap_certify.py` proves the exact rational bound 17/200 for all
432 valid ordered S4-generating level-2 transposition configurations, on all
four non-trivial irreducibles, using invariant Gram matrices and all principal
minors. `s4_hamilton_reduction_check.py` gives a supplementary exact integer
certificate; it is no longer needed for the analytic intersection proof.
The two scripts were moved unchanged from P2A; provenance is in the artifact README.

## Transvection levels and all-pairs evidence (26 September 2026)

`transvection_levels.py` enumerates the transvection-generated subgroup lattice
of SL3(3), up to conjugacy, with exact integer arithmetic. It finds 104
transvections and eight non-trivial subgroup conjugacy classes including G. There are no
proper non-solvable cases. The distinct proper orders and derived series are:

| Order | Derived orders |
|---:|---|
| 3 | 3→1 |
| 9 | 9→1 |
| 24 | 24→8→2→1 |
| 27 | 27→3→1 |
| 216 | 216→72→18→9→1 |

The largest derived length is 4. This independently cross-checks the short
structural proof using the cited maximal-subgroup list and 3²:SL2(3).

`external/p2a/s4_allpairs_doublecoset_audit.py 11 13 19 29 37` reproduces:

| q | Largest numerical estimate | Intersection order at that estimate |
|---:|---:|---:|
| 11 | 0.683013 | 6 |
| 13 | 0.750000 | 6 |
| 19 | 0.809017 | 6 |
| 29 | 0.854361 | 6 |
| 37 | 0.845645 | 6 |

Double cosets are enumerated exactly. The norms are estimated by 80 power
iterations and four starts; they are not certified upper bounds. The generator
test uses maximality of the exterior S4 at these primes, with a capped closure.
Both required scripts were copied unchanged from P2A; source hashes and origin
are in `external/p2a/PROVENANCE.json`.

The congruence-lifting proposition is analytic and uses Serre's lemma; the
exceptional-configuration corollary uses Becker–Breuillard's dimension-zero
theorem and the Cayley comparison. Finite experiments do not certify either
infinite-family statement.

## Intersection-six stratum and recomputed all-pairs table (16 October 2026)

`octahedral_cap6_check.py` (q = 5, 11, 13, 19, 29, 37; SymPy required) checks the
ingredients of the intersection-six proposition with exact arithmetic, and
recomputes the manuscript's all-pairs table by one numerical procedure.

**Exact (integer/rational arithmetic).** For every double coset AgA the order
|A ∩ gAg⁻¹| is computed; the census never exceeds 6, the double coset of order
six is unique and equals AhA, and the four subgroups B with |A ∩ B| = 6 are the
A-conjugates of hAh⁻¹. Brute force over G gives N(V)=A, N(A)=A, C(A)=1,
C(K)={1,h}, |N(K)|=12. The adjoint action is checked to be the stated signed
permutation / (2/3)J−I dictionary, with determinant one on all of G for q ≤ 13.

| q | cosets | double cosets (nontrivial) | intersection census (order: count) |
|---:|---:|---:|---|
| 5 | 5 | 1 | 6:1 |
| 11 | 55 | 5 | 1:1, 2:1, 3:1, 4:1, 6:1 |
| 13 | 91 | 7 | 1:2, 2:2, 3:1, 4:1, 6:1 |
| 19 | 285 | 17 | 1:9, 2:3, 3:2, 4:2, 6:1 |
| 29 | 1015 | 51 | 1:37, 2:6, 3:4, 4:3, 6:1 |
| 37 | 2109 | 99 | 1:81, 2:8, 3:5, 4:4, 6:1 |

**Certified natural values.** The Hecke graph on G/A (right multiplication by
the four norm-3 classes) is vertex-transitive, and the quotient by the A-orbits
(double cosets) carries every distinct eigenvalue. SymPy factors the integer
characteristic polynomial of that quotient over Z, isolates its real roots in
rational intervals and refines the extremal one to width 10⁻³⁰.

| q | c_q | extremal eigenvalue (4c_q) | minimal polynomial |
|---:|---|---|---|
| 5 | 1/4 | −1 | z + 1 |
| 11 | (1+√3)/4 | −(1+√3) | z² + 2z − 2 |
| 13 | 3/4 | −3 | z + 3 |
| 19 | (1+√5)/4 | −(1+√5) | z² + 2z − 4 |
| 29 | 0.85436118917716762327 (rounded to 20 digits) | positive, degree 6 | z⁶ − 2z⁵ − 13z⁴ + 17z³ + 52z² − 32z − 64 |
| 37 | 0.84565136442071691397 (rounded to 20 digits) | negative, degree 9 | z⁹ + 7z⁸ + 7z⁷ − 45z⁶ − 100z⁵ + 13z⁴ + 170z³ + 139z² + 36z + 1 |

c_q is a maximum of moduli: at q = 11, 19, 37 the extremal eigenvalue is
negative. The 20-digit values are rounded; the certified rational enclosures of
width below 10⁻³⁰, the full factorisations and the isolating intervals are in
`results/octahedral_cap6_certificates.json`.

**Numerical (float64).** For every nontrivial double coset D, the coset operator
T_D on ℓ²(G/A) commutes with the left G-action, so its singular values are those
of its restriction to A-invariant functions: the orbit-weighted quotient matrix
N_D(i,j)·√(w_i/w_j)/d. One copy of the trivial singular value (the constants) is
deflated by projection. For q ≤ 19 every value is cross-checked against a dense
SVD on ℓ²(G/A); the same cross-check was run once for q = 29, 37
(`--dense-all`, agreement to 10⁻⁸, logged in
`notes/octahedral-cap6-implementation-2026-10-16-1950/`).

| q | all-pairs maximum | attained at intersection order | maxima by intersection order |
|---:|---:|---:|---|
| 11 | 0.683013 | 6 | see JSON |
| 13 | 0.750000 | 6 | see JSON |
| 19 | 0.809017 | 6 | see JSON |
| 29 | 0.854361 | 6 | see JSON |
| 37 | 0.845651 | 6 | see JSON |

These maxima are numerical. They agree, to six decimals, with the separately
certified natural values; they do not certify that the natural pair is extremal
for other primes. The historical table produced by
`external/p2a/s4_allpairs_doublecoset_audit.py` (80 power iterations) gave the
same values except 0.845645 at q = 37, an underestimate of the certified value by
6·10⁻⁶; the historical script and its provenance remain bundled unchanged.
