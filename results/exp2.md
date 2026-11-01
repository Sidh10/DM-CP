# Experiment 2: Counting versus enumeration

- run date: 2026-10-31
- command: `python -m experiments.exp2_counting_enumeration`
- Python 3.13.11 on Windows

## Question
Does recursive enumeration of all sequences at Hamming distance exactly k from a reference produce exactly C(n,k)·3^k distinct strings, matching the recurrence M(n,k) = M(n−1,k) + 3M(n−1,k−1)? Where does enumeration stop being feasible, and how sparse are real observations within the possible space?

## Input
Reference prefixes of GATTACGT for n = 0..8; Dataset A (syn01, syn04) for orderings; Dataset B (34 clean sequences, 30 haplotypes) for pigeonhole and sparsity.

## Mathematical model
N(n) = 4^n; M(n,k) = C(n,k)·3^k; M(n,k) = M(n−1,k) + 3M(n−1,k−1), M(0,0) = 1; |B(n,r)| = Σ_{k≤r} M(n,k); inclusion–exclusion U(n) = 4^n − 4·3^n + 6·2^n − 4; pigeonhole L − k + 1 > 4^k; k! shortest substitution orderings.

## Procedure
For each (n,k) with n ≤ 8, generate all mutants recursively (choose k positions in increasing order, then 3 bases each), check that they are distinct and at distance exactly k, and count them. Compare with the hand table, the closed form and the recurrence table. Exercise the 10^5 enumeration cap. Count all-four-base strings by inclusion–exclusion and by enumeration. Apply the pigeonhole bound to every clean Dataset B sequence. Compare observed neighbour counts in Dataset B with the Hamming-ball sizes.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | M(n,k) by recursive enumeration, n ≤ 8 | [[1], [1, 3], [1, 6, 9], [1, 9, 27, 27], [1, 12, 54, 108, 81], [1, 15, 90, 270, 405, 243], [1, 18, 135, 540, 1215, 1458, 729], [1, 21, 189, 945, 2835, 5103, 5103, 2187], [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561]] | docs/unit2_expected.md §2 |
| 2 | M(n,k) by closed form C(n,k)·3^k, n ≤ 8 | [[1], [1, 3], [1, 6, 9], [1, 9, 27, 27], [1, 12, 54, 108, 81], [1, 15, 90, 270, 405, 243], [1, 18, 135, 540, 1215, 1458, 729], [1, 21, 189, 945, 2835, 5103, 5103, 2187], [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561]] | docs/unit2_expected.md §2 |
| 3 | M(n,k) by recurrence table, n ≤ 8 | [[1], [1, 3], [1, 6, 9], [1, 9, 27, 27], [1, 12, 54, 108, 81], [1, 15, 90, 270, 405, 243], [1, 18, 135, 540, 1215, 1458, 729], [1, 21, 189, 945, 2835, 5103, 5103, 2187], [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561]] | docs/unit2_expected.md §2 |
| 4 | cells where all strings are distinct and at distance exactly k | 45 | docs/unit2_expected.md §2 — *consistency check* |
| 5 | row sums Σ_k M(n,k), n ≤ 8 | [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536] | docs/unit2_expected.md §1 |
| 6 | M(12,k) for k = 1,2,3 | {1: 36, 2: 594, 3: 5940} | docs/synthetic_expected.md §8 |
| 7 | P(12,3) | 1320 | docs/unit2_expected.md §3 |
| 8 | Hamming-ball sizes \|B(n,r)\| | {B(8, 2): 277, B(12, 2): 631, B(1701, 1): 5104} | docs/unit2_expected.md §2 |
| 9 | enumerated ball B(8,8) = whole Σ^8 | 65536 | docs/unit2_expected.md §2 |
| 10 | decimal digits of 4^1701 | 1025 | docs/unit2_expected.md §1 |
| 11 | guard outcome (cap 10^5) and count | {(12, 6): (refused, 673596), (1701, 2): (refused, 13012650), (1701, 1): (allowed, 5103)} | docs/unit2_expected.md §2 |
| 12 | strings using all four bases, by inclusion–exclusion | {4: 24, 5: 240, 6: 1560, 7: 8400, 8: 40824, 12: 14676024} | docs/unit2_expected.md §4 |
| 13 | strings using all four bases, by enumeration (n = 4..8) | {4: 24, 5: 240, 6: 1560, 7: 8400, 8: 40824} | docs/unit2_expected.md §4 |
| 14 | syn01 → syn04: shortest orderings, fully observed ones | {count: 6, fully_observed: 1} | docs/unit2_expected.md §3 |
| 15 | pigeonhole on L = 1701 | {max_forced_k: 5, forced_k5: True, forced_k6: False, clean_with_repeat_5mer: 34} | docs/unit2_expected.md §5 |
| 16 | clean B sequences with a repeated 6-mer (not forced) | *(none: reported as observation)* |  |
| 17 | B: Σ observed neighbours at d = 1 (= 2 × pairs) | 16 | docs/dataset_b_qc.md (independent script) |
| 18 | B: Σ observed neighbours at d ≤ 2 (= 2 × pairs) | 62 | docs/dataset_b_qc.md (independent script) |
| 19 | B: max observed vs possible neighbours | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | M(n,k) by recursive enumeration, n ≤ 8 | [[1], [1, 3], [1, 6, 9], [1, 9, 27, 27], [1, 12, 54, 108, 81], [1, 15, 90, 270, 405, 243], [1, 18, 135, 540, 1215, 1458, 729], [1, 21, 189, 945, 2835, 5103, 5103, 2187], [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561]] |
| 2 | M(n,k) by closed form C(n,k)·3^k, n ≤ 8 | [[1], [1, 3], [1, 6, 9], [1, 9, 27, 27], [1, 12, 54, 108, 81], [1, 15, 90, 270, 405, 243], [1, 18, 135, 540, 1215, 1458, 729], [1, 21, 189, 945, 2835, 5103, 5103, 2187], [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561]] |
| 3 | M(n,k) by recurrence table, n ≤ 8 | [[1], [1, 3], [1, 6, 9], [1, 9, 27, 27], [1, 12, 54, 108, 81], [1, 15, 90, 270, 405, 243], [1, 18, 135, 540, 1215, 1458, 729], [1, 21, 189, 945, 2835, 5103, 5103, 2187], [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561]] |
| 4 | cells where all strings are distinct and at distance exactly k | 45 |
| 5 | row sums Σ_k M(n,k), n ≤ 8 | [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536] |
| 6 | M(12,k) for k = 1,2,3 | {1: 36, 2: 594, 3: 5940} |
| 7 | P(12,3) | 1320 |
| 8 | Hamming-ball sizes \|B(n,r)\| | {B(8, 2): 277, B(12, 2): 631, B(1701, 1): 5104} |
| 9 | enumerated ball B(8,8) = whole Σ^8 | 65536 |
| 10 | decimal digits of 4^1701 | 1025 |
| 11 | guard outcome (cap 10^5) and count | {(12, 6): (refused, 673596), (1701, 2): (refused, 13012650), (1701, 1): (allowed, 5103)} |
| 12 | strings using all four bases, by inclusion–exclusion | {4: 24, 5: 240, 6: 1560, 7: 8400, 8: 40824, 12: 14676024} |
| 13 | strings using all four bases, by enumeration (n = 4..8) | {4: 24, 5: 240, 6: 1560, 7: 8400, 8: 40824} |
| 14 | syn01 → syn04: shortest orderings, fully observed ones | {count: 6, fully_observed: 1} |
| 15 | pigeonhole on L = 1701 | {max_forced_k: 5, forced_k5: True, forced_k6: False, clean_with_repeat_5mer: 34} |
| 16 | clean B sequences with a repeated 6-mer (not forced) | 34/34 |
| 17 | B: Σ observed neighbours at d = 1 (= 2 × pairs) | 16 |
| 18 | B: Σ observed neighbours at d ≤ 2 (= 2 × pairs) | 62 |
| 19 | B: max observed vs possible neighbours | {d=1: 6 of 5103, d≤2: 9 of 13017753} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | M(n,k) by recursive enumeration, n ≤ 8 | yes |
| 2 | M(n,k) by closed form C(n,k)·3^k, n ≤ 8 | yes |
| 3 | M(n,k) by recurrence table, n ≤ 8 | yes |
| 4 | cells where all strings are distinct and at distance exactly k | yes (consistency) |
| 5 | row sums Σ_k M(n,k), n ≤ 8 | yes |
| 6 | M(12,k) for k = 1,2,3 | yes |
| 7 | P(12,3) | yes |
| 8 | Hamming-ball sizes \|B(n,r)\| | yes |
| 9 | enumerated ball B(8,8) = whole Σ^8 | yes |
| 10 | decimal digits of 4^1701 | yes |
| 11 | guard outcome (cap 10^5) and count | yes |
| 12 | strings using all four bases, by inclusion–exclusion | yes |
| 13 | strings using all four bases, by enumeration (n = 4..8) | yes |
| 14 | syn01 → syn04: shortest orderings, fully observed ones | yes |
| 15 | pigeonhole on L = 1701 | yes |
| 16 | clean B sequences with a repeated 6-mer (not forced) | n/a |
| 17 | B: Σ observed neighbours at d = 1 (= 2 × pairs) | yes |
| 18 | B: Σ observed neighbours at d ≤ 2 (= 2 × pairs) | yes |
| 19 | B: max observed vs possible neighbours | n/a |

**16/16 independent expected values matched; 1/1 consistency checks passed; 2 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

### Enumeration cost at n = 8 (observation; wall-clock time depends on the machine)

| k | strings M(8,k) | recursive calls | time (ms) |
|---|---|---|---|
| 0 | 1 | 1 | 0.0 |
| 1 | 24 | 25 | 0.0 |
| 2 | 252 | 274 | 0.1 |
| 3 | 1512 | 1720 | 0.5 |
| 4 | 5670 | 6766 | 1.9 |
| 5 | 13608 | 17086 | 4.8 |
| 6 | 20412 | 27064 | 7.4 |
| 7 | 17496 | 24604 | 6.7 |
| 8 | 6561 | 9841 | 2.7 |

## Interpretation
Enumeration, closed form and recurrence agree on every cell for n ≤ 8; the induction proof in `docs/unit2_expected.md` §2 extends this to all n. The cost of enumeration is proportional to the number of strings, so the cap is reached quickly: M(12,6) is already 673,596 and M(1701,2) is 13,012,650. In Dataset B, each haplotype has at most 6 observed one-substitution neighbours among 5103 possible ones. The single-substitution graph built in Unit V is therefore an extremely sparse sample of the Hamming graph. This describes the sample, not the virus population. Of the 6 shortest orderings from syn01 to syn04, only one passes entirely through observed sequences.

## Limitations
Counts are combinatorial possibilities, not probabilities: they ignore transition/transversion bias, codon structure and viability. Substitutions only. Timings are machine-dependent observations. The finite checks for n ≤ 8 are evidence; the general statement rests on the induction proof.
