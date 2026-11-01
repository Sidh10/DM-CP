# Experiment 1: Mutation-set analysis

- run date: 2026-10-31
- command: `python -m experiments.exp1_mutation_sets`
- Python 3.13.11 on Windows

## Question
Do mutation sets, set operations, Hamming distances and the threshold relation R_k computed by `genomath` reproduce the hand-derived values for Dataset A, and do the stated identities hold on real data?

## Input
Dataset A (`data/synthetic/dataset_a.fasta`, reference syn01) plus the two records in `dataset_a_invalid.fasta`; Dataset B (`data/public_genomic/h3n2_ha_cds.fasta`, medoid reference).

## Mathematical model
M(R,S) = {i : R_i ≠ S_i}; M*(R,S) = {(i, S_i) : R_i ≠ S_i}; d_H(S,T) = |M(R,S) Δ M(R,T)| + |{i ∈ M(R,S) ∩ M(R,T) : S_i ≠ T_i}|; S R_k T ⇔ d_H(S,T) ≤ k; components of G_k = classes of the transitive closure of R_k (Warshall).

## Procedure
Load and preprocess each dataset (uppercase, drop non-ACGT, require one length, collapse haplotypes, choose reference). Compute every quantity in the Expected table with `genomath.set_logic` and compare against values fixed beforehand. Dataset B has no hand values: its expected entries come from the independent QC script or from theorems, and items without an independent value are reported as observations.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | A: haplotype ids | [syn01, syn02, syn03, syn04, syn05, syn06, syn07, syn08, syn09] | docs/synthetic_expected.md §3 |
| 2 | A: identical group | [syn03, syn10] | docs/synthetic_expected.md §3 |
| 3 | A + bad01: dropped record and reason | [(bad01, non-ACGT characters Nx1)] | data/synthetic/dataset_a_invalid.fasta |
| 4 | A + bad02 (length 11): outcome | rejected: sequences have different lengths {11: 1, 12: 10} | final_spec.md §4 rule 3 |
| 5 | A: M*(R,S) for all 9 haplotypes | {syn01: [], syn02: [(2, C)], syn03: [(2, C), (5, G)], syn04: [(2, C), (5, G), (9, T)], syn05: [(11, A)], syn06: [(2, C), (11, A)], syn07: [(5, C), (7, A)], syn08: [(8, C), (11, G)], syn09: [(1, C), (3, G), (4, C), (6, A), (10, G), (12, A)]} | docs/synthetic_expected.md §1 |
| 6 | A: M3 vs M7 set operations | {intersection: [5], union: [2, 5, 7], a_minus_b: [2], b_minus_a: [7], symmetric_difference: [2, 7]} | docs/synthetic_expected.md §2 |
| 7 | A: M*3 ∩ M*7 | [] | docs/synthetic_expected.md §2 |
| 8 | A: chain M2 ∩ M4, M4 \ M3 | {M2 ∩ M4: [2], M4 \ M3: [9]} | docs/synthetic_expected.md §2 |
| 9 | A: d_H(syn03,syn07) decomposition | {sym_diff: [2, 7], shared_differing: [5], d_H: 3} | docs/synthetic_expected.md §2 |
| 10 | A: d_H matrix entries equal to hand matrix | 81 | docs/synthetic_expected.md §3 |
| 11 | A: pairs where the d_H identity holds | 36 | math_audit.md Unit I (theorem) — *consistency check* |
| 12 | A: \|R_1\| (ordered pairs) | 21 | docs/synthetic_expected.md §3 |
| 13 | A: R_1 reflexive / symmetric / transitive | {reflexive: True, symmetric: True, transitive: False} | docs/synthetic_expected.md §3 |
| 14 | A: R_1 transitivity counterexample | (syn01, syn02, syn03) | docs/synthetic_expected.md §3 |
| 15 | A: components of closure(R_k) for k in {0,1,2,5,6} | {0: 9, 1: 4, 2: 2, 5: 2, 6: 1} | docs/synthetic_expected.md §4 |
| 16 | A: predicate filters | {P ∧ Q: [syn02, syn03, syn06], P ∧ ¬Q: [syn01, syn05, syn07, syn08], ¬P: [syn04, syn09]} | docs/synthetic_expected.md §2 |
| 17 | A: orientation class of syn03 | [syn03, syn08, syn10] | docs/synthetic_expected.md §6 |
| 18 | A: RC-fixed records | [syn01] | docs/synthetic_expected.md §6 |
| 19 | B: dropped records | [MH246583.1, MK181701.1] | docs/dataset_b_qc.md (independent script) |
| 20 | B: haplotypes | 30 | docs/dataset_b_qc.md (independent script) |
| 21 | B: pairs where the d_H identity holds | 435 | math_audit.md Unit I (theorem) — *consistency check* |
| 22 | B: \|R_1\| = 30 + 2 × (d_H = 1 pairs) | 46 | docs/dataset_b_qc.md (independent script) |
| 23 | B: R_1 transitive? (witness) | *(none: reported as observation)* |  |
| 24 | B: medoid reference | *(none: reported as observation)* |  |
| 25 | B: \|M(R,S)\| distribution (size: count) | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | A: haplotype ids | [syn01, syn02, syn03, syn04, syn05, syn06, syn07, syn08, syn09] |
| 2 | A: identical group | [syn03, syn10] |
| 3 | A + bad01: dropped record and reason | [(bad01, non-ACGT characters Nx1)] |
| 4 | A + bad02 (length 11): outcome | rejected: sequences have different lengths {11: 1, 12: 10} |
| 5 | A: M*(R,S) for all 9 haplotypes | {syn01: [], syn02: [(2, C)], syn03: [(2, C), (5, G)], syn04: [(2, C), (5, G), (9, T)], syn05: [(11, A)], syn06: [(2, C), (11, A)], syn07: [(5, C), (7, A)], syn08: [(8, C), (11, G)], syn09: [(1, C), (3, G), (4, C), (6, A), (10, G), (12, A)]} |
| 6 | A: M3 vs M7 set operations | {intersection: [5], union: [2, 5, 7], a_minus_b: [2], b_minus_a: [7], symmetric_difference: [2, 7]} |
| 7 | A: M*3 ∩ M*7 | [] |
| 8 | A: chain M2 ∩ M4, M4 \ M3 | {M2 ∩ M4: [2], M4 \ M3: [9]} |
| 9 | A: d_H(syn03,syn07) decomposition | {sym_diff: [2, 7], shared_differing: [5], d_H: 3} |
| 10 | A: d_H matrix entries equal to hand matrix | 81 |
| 11 | A: pairs where the d_H identity holds | 36 |
| 12 | A: \|R_1\| (ordered pairs) | 21 |
| 13 | A: R_1 reflexive / symmetric / transitive | {reflexive: True, symmetric: True, transitive: False} |
| 14 | A: R_1 transitivity counterexample | (syn01, syn02, syn03) |
| 15 | A: components of closure(R_k) for k in {0,1,2,5,6} | {0: 9, 1: 4, 2: 2, 5: 2, 6: 1} |
| 16 | A: predicate filters | {P ∧ Q: [syn02, syn03, syn06], P ∧ ¬Q: [syn01, syn05, syn07, syn08], ¬P: [syn04, syn09]} |
| 17 | A: orientation class of syn03 | [syn03, syn08, syn10] |
| 18 | A: RC-fixed records | [syn01] |
| 19 | B: dropped records | [MH246583.1, MK181701.1] |
| 20 | B: haplotypes | 30 |
| 21 | B: pairs where the d_H identity holds | 435 |
| 22 | B: \|R_1\| = 30 + 2 × (d_H = 1 pairs) | 46 |
| 23 | B: R_1 transitive? (witness) | (False, (CY249987.1, MH637643.1, MH082676.1)) |
| 24 | B: medoid reference | CY249987.1 |
| 25 | B: \|M(R,S)\| distribution (size: count) | {0: 1, 1: 6, 2: 3, 3: 5, 4: 1, 5: 3, 6: 2, 8: 1, 27: 1, 29: 2, 30: 1, 31: 1, 32: 1, 48: 1, 49: 1} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | A: haplotype ids | yes |
| 2 | A: identical group | yes |
| 3 | A + bad01: dropped record and reason | yes |
| 4 | A + bad02 (length 11): outcome | yes |
| 5 | A: M*(R,S) for all 9 haplotypes | yes |
| 6 | A: M3 vs M7 set operations | yes |
| 7 | A: M*3 ∩ M*7 | yes |
| 8 | A: chain M2 ∩ M4, M4 \ M3 | yes |
| 9 | A: d_H(syn03,syn07) decomposition | yes |
| 10 | A: d_H matrix entries equal to hand matrix | yes |
| 11 | A: pairs where the d_H identity holds | yes (consistency) |
| 12 | A: \|R_1\| (ordered pairs) | yes |
| 13 | A: R_1 reflexive / symmetric / transitive | yes |
| 14 | A: R_1 transitivity counterexample | yes |
| 15 | A: components of closure(R_k) for k in {0,1,2,5,6} | yes |
| 16 | A: predicate filters | yes |
| 17 | A: orientation class of syn03 | yes |
| 18 | A: RC-fixed records | yes |
| 19 | B: dropped records | yes |
| 20 | B: haplotypes | yes |
| 21 | B: pairs where the d_H identity holds | yes (consistency) |
| 22 | B: \|R_1\| = 30 + 2 × (d_H = 1 pairs) | yes |
| 23 | B: R_1 transitive? (witness) | n/a |
| 24 | B: medoid reference | n/a |
| 25 | B: \|M(R,S)\| distribution (size: count) | n/a |

**20/20 independent expected values matched; 2/2 consistency checks passed; 3 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

## Interpretation
Dataset A: the shared-position example (syn03 vs syn07) shows why position-level sets are not enough: they share position 5 but no mutation, and that position still contributes 1 to d_H. R_1 is reflexive and symmetric but not transitive, so it is a similarity relation, not an equivalence; its transitive closure is an equivalence whose classes are the components of the threshold graph.

Dataset B (observation, no biological claim): relative to the medoid CY249987.1, mutation-set sizes are distributed as {0: 1, 1: 6, 2: 3, 3: 5, 4: 1, 5: 3, 6: 2, 8: 1, 27: 1, 29: 2, 30: 1, 31: 1, 32: 1, 48: 1, 49: 1}. R_1 is not transitive on real data either: CY249987.1 R_1 MH637643.1 and MH637643.1 R_1 MH082676.1 but not CY249987.1 R_1 MH082676.1.

## Limitations
Substitutions only; positional comparison assumes equal-length, colinear sequences (Dataset B: equal CDS length and ≤ 3.4% divergence, not an alignment). Mutation sets depend on the chosen reference; pairwise distances do not. Two Dataset B records with ambiguity codes were dropped, not resolved.
