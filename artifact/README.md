# Standalone verification artifact for P1

Paper: *The solvability dichotomy for commutator-tree mixing in finite groups*.
Python 3.9+, NumPy 1.24+ (below major version 3) and SymPy 1.12+ (below major
version 2) are required; SymPy is used only by `octahedral_cap6_check.py`, for
exact certificates (integer characteristic polynomials, factorisation over Z and
rational root isolation). The scripts use no network services or sibling
repositories. Install `requirements.txt`, then:

```sh
python3 scripts/run_validation.py --output validation-output --jobs 2
```

The runner records every command, exit status, duration and log hash. It limits
BLAS threads for reproducibility. Individual commands can be run from `scripts/`.
The full suite includes 26 commands; regular-representation and PSL3(3) sweeps
can take several minutes. `SHA256SUMS` records every distributed file except itself.

## What is checked

- `bridge_realizability_certificate.py`, `target_transition_graph.py`,
  `binary_target_grammar_certificate.py`, and `test_bridge_sweeps.py`: finite
  target/rigid-family data. Mixed-transfer checks are exhaustive on S5 and use
  the first three conjugates per class elsewhere.
- `a5_criterion.py`: complete closure lattice and ordered level-2 census;
  operator norms use numerical SVD, not interval arithmetic.
- `verification_groups.py`: independent exact ordered level-2 censuses and
  C-generated subgroup lattices for five small-group cases.
- `test_strengthening.py`: field/class/obstruction regressions, both spectral
  inequalities (numerical checks), product boundary, exact GL2(Z4) derived series,
  and intentionally injected failures in the criterion scripts.
- `psl33_maximals_check.py`: exact finite data for level 5. Completeness of the
  maximal-subgroup list is a cited ATLAS input, not a computational conclusion.
- `transvection_levels.py`: exact closure lattice up to conjugacy for the 104
  transvections of SL3(3); eight non-trivial subgroup classes including the ambient group,
  largest proper derived length 4 attained at order 216. This cross-checks the
  manuscript's structural proof.
- `octahedral_cap6_check.py`: the manuscript's all-pairs table (q = 11, 13, 19,
  29, 37) by one floating-point procedure (orbit-weighted double-coset quotients,
  one copy of the trivial singular value deflated, dense SVD cross-check for
  q <= 19), together with the exact checks listed under "Natural octahedral pair".
- `external/p2a/s4_allpairs_doublecoset_audit.py 11 13 19 29 37`: historical
  numerical all-pairs estimates (80 power iterations, four starts) with an exact
  double-coset enumeration. Two unchanged P2A scripts and provenance are bundled
  together; its estimates are not certified upper bounds and its q = 37 value
  (0.845645) underestimates the certified natural value 0.845651 by 6e-6.
- `pgl2_classification.py`: finite-field subgroup witness sampling. Its level
  predictions come from the analytic theorem; failure to find a witness proves nothing.
- `generates_G.py`, `geometric_lemma.py`, `derived_length_test.py`, `gl3f2.py`,
  `gl2z4.py`: reduced numerical sweeps, with scope recorded in `results/RESULTS.md`.
- `mixing_level_pgl.py 7`, `bottleneck_pgl7.py`: sampled norms. The former selects
  exterior involutions by determinant and rejects non-prime q.
- `noninvolution_leaves.py`: 14 selected leaves out of 98 candidates; capped
  numerical experiment, no generation check or conjugation-closure claim.
- `bridge_verify_counterexample.py --case all`: exact explicit witnesses.
- `external/p0/20_level3_exact_certificate.py`: unchanged P0 certificate for the
  S5 bound 71/2000. Provenance, hash and attribution are beside the script.

The analytic dichotomy, rank-one classification and gap comparison are proved in
the manuscript. The code distinguishes exact counts from numerical estimates,
ordered universes from reduced sweeps, and sampled evidence from certificates.
See `results/RESULTS.md`. No q-uniform gap is inferred from a finite mixing level.

P1 code is distributed under `LICENSE`. The external P0 component retains its
separate attribution and source provenance. This is a local artifact prepared for
release; it has not been deposited by this build command.

## Arbitrary-order stability

Run `python3 scripts/test_quantitative_families.py` for small-group stability
and support checks. Norm calculations are numerical; the rate and depth
theorems in the manuscript are analytic and use a cited diameter theorem.

## Natural octahedral pair

- `octahedral_pair_check.py`: exact finite group and compression identities for
  q=5,11,13,19; q=5 exact spectral distinction; other spectra and generic gluing
  tests are numerical. These checks do not prove the infinite-family theorem.
- `s4_hamilton_reduction_check.py`: supplementary integer reduction certificate;
  the manuscript now has a short analytic intersection proof.
- `s4_internal_gap_certify.py`: exact rational λ≤17/200 certificate for all 432
  valid ordered transposition configurations and four non-trivial S4 irreducibles.
- `octahedral_cap6_check.py` (q = 5, 11, 13, 19, 29, 37): exact census of the
  intersection orders |A ∩ gAg⁻¹| over all double cosets (never above 6), the
  unique double coset of order six (AhA) and its four subgroups B, the brute-force
  computations N(V)=A, N(A)=A, C(A)=1, C(K)={1,h}, |N(K)|=12, and the orthogonal
  dictionary (A as signed permutation matrices of determinant one, h as (2/3)J−I;
  determinant one of the adjoint action checked on all of G for q ≤ 13). The
  natural value c_q is certified exactly from the integer characteristic
  polynomial of the double-coset quotient of the Hecke graph, which carries every
  distinct eigenvalue: closed forms 1/4, (1+√3)/4, 3/4, (1+√5)/4 at q = 5, 11,
  13, 19, and rational enclosures of width below 1e-30 at q = 29, 37, rounding
  to 0.85436118917716762327 and 0.84565136442071691397 (extremal eigenvalues of
  minimal polynomial degrees 6 and 9). Polynomials, factorisations, isolating
  intervals and enclosures are in `results/octahedral_cap6_certificates.json`.
  Its all-pairs norms are numerical.

The two P2A scripts were transferred unchanged from P2A on 2026-09-23.
Their SHA-256 values are respectively
`224a3ba8906c24c1b5962355d68e184aada376fff24d94bca7c957baa3a2a0f3` and
`98a68f2bc5c5f794f45e002b743d0db0f220adaa431766a6bf0c7813ec3a07d5`.
The angle bound for every exterior octahedral pair with intersection of order
six reduces, by the simultaneous-conjugacy proposition, to the natural pair; it
uses DSV (2003), Theorem 4.2.2 at p=3, and asymptotic sharpness uses
Alon–Boppana. Neither is inferred from finite sampling.
