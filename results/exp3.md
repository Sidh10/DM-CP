# Experiment 3: Algebraic verification

- run date: 2026-10-31
- command: `python -m experiments.exp3_algebraic_verification`
- Python 3.13.11 on Windows

## Question
Do the strand-orientation maps {I, RC} and {I, R, C, RC} satisfy the group axioms when verified exhaustively, do Burnside's orbit counts agree with direct enumeration, and why do substitution operators fail to form a group?

## Input
All words over {A,C,G,T} of length 0..6 (function equality) and of each length 0..8 (fixed points, orbits); Dataset A haplotypes; Dataset B haplotypes.

## Mathematical model
Maps I, R, C, RC on Σ^n under composition (a∘b)(x) = a(b(x)); group axioms; Burnside #orbits = (1/|G|)·Σ|Fix(g)|; orbit–stabiliser |Gx|·|G_x| = |G|; the monoid of substitution operators f_p ('set position i to b').

## Procedure
Build each Cayley-table entry by composing the maps and identifying the result with an element, comparing on all 5,461 words of length ≤ 6. Check closure, associativity (all |G|³ triples), identity, inverses and commutativity from the table. Count fixed points by formula and by enumeration; count orbits by Burnside and by enumeration; generate the substitution monoid by closure under composition and test the same axioms. Apply the groups to Datasets A and B.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | Cayley table of {I, RC} | {I: [I, RC], RC: [RC, I]} | docs/unit3_expected.md §1 |
| 2 | Cayley table of {I, R, C, RC} | {I: [I, R, C, RC], R: [R, I, RC, C], C: [C, RC, I, R], RC: [RC, C, R, I]} | docs/unit3_expected.md §1 |
| 3 | axioms for {I, RC} | {closure: True, associativity: True, identity: True, inverses: True, commutative: True} | docs/unit3_expected.md §1 |
| 4 | axioms for {I, R, C, RC} | {closure: True, associativity: True, identity: True, inverses: True, commutative: True} | docs/unit3_expected.md §1 |
| 5 | element orders in V4 | {I: 1, R: 2, C: 2, RC: 2} | docs/unit3_expected.md §1 |
| 6 | V4 cyclic? | False | docs/unit3_expected.md §1 |
| 7 | isomorphic to Z2×Z2; homomorphic bijections onto Z4 | (True, 0) | docs/unit3_expected.md §1 |
| 8 | subgroups of V4 | [[I], [I, R], [I, C], [I, RC], [I, R, C, RC]] | docs/unit3_expected.md §2 |
| 9 | left cosets of {I,RC} in V4 | [[I, RC], [R, C]] | docs/unit3_expected.md §2 |
| 10 | which maps coincide on Σ^n, n = 0,1,2 | {0: [[I, R, C, RC]], 1: [[I, R], [C, RC]], 2: [[I], [R], [C], [RC]]} | docs/unit3_expected.md §1 |
| 11 | \|Fix(g)\| by formula, n = 0..8 | {I: [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536], R: [1, 4, 4, 16, 16, 64, 64, 256, 256], C: [1, 0, 0, 0, 0, 0, 0, 0, 0], RC: [1, 0, 4, 0, 16, 0, 64, 0, 256]} | docs/unit3_expected.md §3 |
| 12 | \|Fix(g)\| by enumeration, n = 0..8 | {I: [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536], R: [1, 4, 4, 16, 16, 64, 64, 256, 256], C: [1, 0, 0, 0, 0, 0, 0, 0, 0], RC: [1, 0, 4, 0, 16, 0, 64, 0, 256]} | docs/unit3_expected.md §3 |
| 13 | {I,RC} orbits by Burnside, n = 0..8 | [1, 2, 10, 32, 136, 512, 2080, 8192, 32896] | docs/unit3_expected.md §4 |
| 14 | V4 orbits by Burnside, n = 0..8 | [1, 2, 6, 20, 72, 272, 1056, 4160, 16512] | docs/unit3_expected.md §4 |
| 15 | {I,RC} orbits by enumeration, n = 0..8 | [1, 2, 10, 32, 136, 512, 2080, 8192, 32896] | docs/unit3_expected.md §4 |
| 16 | V4 orbits by enumeration, n = 0..8 | [1, 2, 6, 20, 72, 272, 1056, 4160, 16512] | docs/unit3_expected.md §4 |
| 17 | orbit counts at n = 12 (Burnside) | {Z2: 8390656, V4: 4196352} | docs/unit3_expected.md §4 |
| 18 | V4 orbit-size distribution, n = 1..4 | {1: {2: 2}, 2: {2: 4, 4: 2}, 3: {2: 8, 4: 12}, 4: {2: 16, 4: 56}} | docs/unit3_expected.md §4 |
| 19 | orbit–stabiliser holds for every orbit, n = 1..4 | True | docs/unit3_expected.md §4 — *consistency check* |
| 20 | complement as a permutation of Σ | {cycles: [(A, T), (C, G)], order: 2, parity: even} | docs/unit3_expected.md §5 |
| 21 | reversal of positions 1..12 | {transpositions: 6, order: 2, parity: even} | docs/unit3_expected.md §5 |
| 22 | substitution monoid: size and units, n = 1..3 | {1: [5, [·]], 2: [25, [··]], 3: [125, [···]]} | docs/unit3_expected.md §6 |
| 23 | substitution monoid is a group / is commutative | (False, False) | docs/unit3_expected.md §6 |
| 24 | non-injectivity witness (n = 2) | *(none: reported as observation)* |  |
| 25 | concatenation laws, words of length ≤ 3 | {C(xy) = C(x)C(y): True, R(xy) = R(y)R(x): True, RC(xy) = RC(y)RC(x): True} | docs/unit3_expected.md §7 |
| 26 | A: {I,RC} classes among 9 haplotypes, and the pair | (8, [syn03, syn08]) | docs/unit3_expected.md §8 |
| 27 | A: V4 orbit of syn01 | {size: 2, stabiliser: [I, RC]} | docs/unit3_expected.md §8 |
| 28 | A: V4 orbit sizes of syn03, syn09 | [4, 4] | docs/unit3_expected.md §8 |
| 29 | A: canonical representative of syn03 and syn08 | [GATTACGCAAGC, GATTACGCAAGC] | docs/unit3_expected.md §8 |
| 30 | A: R-fixed haplotypes | [] | docs/unit3_expected.md §8 |
| 31 | B: RC pairs, RC-fixed, {I,RC} classes | {RC_pairs: 0, RC_fixed: 0, Z2_classes: 30} | docs/unit3_expected.md §8 |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | Cayley table of {I, RC} | {I: [I, RC], RC: [RC, I]} |
| 2 | Cayley table of {I, R, C, RC} | {I: [I, R, C, RC], R: [R, I, RC, C], C: [C, RC, I, R], RC: [RC, C, R, I]} |
| 3 | axioms for {I, RC} | {closure: True, associativity: True, identity: True, inverses: True, commutative: True} |
| 4 | axioms for {I, R, C, RC} | {closure: True, associativity: True, identity: True, inverses: True, commutative: True} |
| 5 | element orders in V4 | {I: 1, R: 2, C: 2, RC: 2} |
| 6 | V4 cyclic? | False |
| 7 | isomorphic to Z2×Z2; homomorphic bijections onto Z4 | (True, 0) |
| 8 | subgroups of V4 | [[I], [I, R], [I, C], [I, RC], [I, R, C, RC]] |
| 9 | left cosets of {I,RC} in V4 | [[I, RC], [R, C]] |
| 10 | which maps coincide on Σ^n, n = 0,1,2 | {0: [[I, R, C, RC]], 1: [[I, R], [C, RC]], 2: [[I], [R], [C], [RC]]} |
| 11 | \|Fix(g)\| by formula, n = 0..8 | {I: [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536], R: [1, 4, 4, 16, 16, 64, 64, 256, 256], C: [1, 0, 0, 0, 0, 0, 0, 0, 0], RC: [1, 0, 4, 0, 16, 0, 64, 0, 256]} |
| 12 | \|Fix(g)\| by enumeration, n = 0..8 | {I: [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536], R: [1, 4, 4, 16, 16, 64, 64, 256, 256], C: [1, 0, 0, 0, 0, 0, 0, 0, 0], RC: [1, 0, 4, 0, 16, 0, 64, 0, 256]} |
| 13 | {I,RC} orbits by Burnside, n = 0..8 | [1, 2, 10, 32, 136, 512, 2080, 8192, 32896] |
| 14 | V4 orbits by Burnside, n = 0..8 | [1, 2, 6, 20, 72, 272, 1056, 4160, 16512] |
| 15 | {I,RC} orbits by enumeration, n = 0..8 | [1, 2, 10, 32, 136, 512, 2080, 8192, 32896] |
| 16 | V4 orbits by enumeration, n = 0..8 | [1, 2, 6, 20, 72, 272, 1056, 4160, 16512] |
| 17 | orbit counts at n = 12 (Burnside) | {Z2: 8390656, V4: 4196352} |
| 18 | V4 orbit-size distribution, n = 1..4 | {1: {2: 2}, 2: {2: 4, 4: 2}, 3: {2: 8, 4: 12}, 4: {2: 16, 4: 56}} |
| 19 | orbit–stabiliser holds for every orbit, n = 1..4 | True |
| 20 | complement as a permutation of Σ | {cycles: [(A, T), (C, G)], order: 2, parity: even} |
| 21 | reversal of positions 1..12 | {transpositions: 6, order: 2, parity: even} |
| 22 | substitution monoid: size and units, n = 1..3 | {1: [5, [·]], 2: [25, [··]], 3: [125, [···]]} |
| 23 | substitution monoid is a group / is commutative | (False, False) |
| 24 | non-injectivity witness (n = 2) | (AA, AA, AC, AA) |
| 25 | concatenation laws, words of length ≤ 3 | {C(xy) = C(x)C(y): True, R(xy) = R(y)R(x): True, RC(xy) = RC(y)RC(x): True} |
| 26 | A: {I,RC} classes among 9 haplotypes, and the pair | (8, [syn03, syn08]) |
| 27 | A: V4 orbit of syn01 | {size: 2, stabiliser: [I, RC]} |
| 28 | A: V4 orbit sizes of syn03, syn09 | [4, 4] |
| 29 | A: canonical representative of syn03 and syn08 | [GATTACGCAAGC, GATTACGCAAGC] |
| 30 | A: R-fixed haplotypes | [] |
| 31 | B: RC pairs, RC-fixed, {I,RC} classes | {RC_pairs: 0, RC_fixed: 0, Z2_classes: 30} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | Cayley table of {I, RC} | yes |
| 2 | Cayley table of {I, R, C, RC} | yes |
| 3 | axioms for {I, RC} | yes |
| 4 | axioms for {I, R, C, RC} | yes |
| 5 | element orders in V4 | yes |
| 6 | V4 cyclic? | yes |
| 7 | isomorphic to Z2×Z2; homomorphic bijections onto Z4 | yes |
| 8 | subgroups of V4 | yes |
| 9 | left cosets of {I,RC} in V4 | yes |
| 10 | which maps coincide on Σ^n, n = 0,1,2 | yes |
| 11 | \|Fix(g)\| by formula, n = 0..8 | yes |
| 12 | \|Fix(g)\| by enumeration, n = 0..8 | yes |
| 13 | {I,RC} orbits by Burnside, n = 0..8 | yes |
| 14 | V4 orbits by Burnside, n = 0..8 | yes |
| 15 | {I,RC} orbits by enumeration, n = 0..8 | yes |
| 16 | V4 orbits by enumeration, n = 0..8 | yes |
| 17 | orbit counts at n = 12 (Burnside) | yes |
| 18 | V4 orbit-size distribution, n = 1..4 | yes |
| 19 | orbit–stabiliser holds for every orbit, n = 1..4 | yes (consistency) |
| 20 | complement as a permutation of Σ | yes |
| 21 | reversal of positions 1..12 | yes |
| 22 | substitution monoid: size and units, n = 1..3 | yes |
| 23 | substitution monoid is a group / is commutative | yes |
| 24 | non-injectivity witness (n = 2) | n/a |
| 25 | concatenation laws, words of length ≤ 3 | yes |
| 26 | A: {I,RC} classes among 9 haplotypes, and the pair | yes |
| 27 | A: V4 orbit of syn01 | yes |
| 28 | A: V4 orbit sizes of syn03, syn09 | yes |
| 29 | A: canonical representative of syn03 and syn08 | yes |
| 30 | A: R-fixed haplotypes | yes |
| 31 | B: RC pairs, RC-fixed, {I,RC} classes | yes |

**29/29 independent expected values matched; 1/1 consistency checks passed; 1 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

## Interpretation
{I, RC} is a group of order 2 and {I, R, C, RC} is the Klein four-group: abelian, every element self-inverse, not cyclic, isomorphic to Z2×Z2 and to no Z4 labelling. Burnside's formula matches enumeration at every length tested, including n = 0 and n = 1, where the four maps are not distinct as functions (the action is not faithful there, and Burnside still holds). Substitution operators close under composition into a monoid with 5^n elements whose only invertible element is the identity, so the statement 'mutations do not form a group' is backed by a computed witness, not just asserted. On Dataset A the RC pair syn03/syn08 is merged into one class. On Dataset B no sequence's reverse complement is present, as predicted from the start/stop codon argument, so Unit III does no analytical work on full annotated coding sequences.

## Limitations
Function equality is checked exhaustively only up to length 6; general lengths rest on the proof in `docs/unit3_expected.md` §1. R and C on their own are not natural strand operations, so the Klein group is a mathematical extension (pedagogical). RC equivalence concerns strand orientation, not mutation. Canonical representatives here are lexicographic; Unit IV restates them as min(enc(S), enc(RC(S))).
