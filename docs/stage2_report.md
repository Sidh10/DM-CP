# Stage 2 report: Unit II (counting and recursion), Experiment 2

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 3 has not started.

## Decision applied provisionally
- Open decision 7: the enumeration cap is **10^5 strings** (`DEFAULT_CAP` in `counting_recursion.py`; every generator
  also takes a `cap=` argument). Change it if you prefer another value.
- The Stage 1 provisional decisions (trace shape, M*, Warshall closure, medoid reference) are unchanged.

## What exists
| Path | Contents |
|---|---|
| `docs/unit2_expected.md` | hand-derived M(n,k) table for n ≤ 8 (each row checked against 4^n); recurrence derivation; induction proof; ball sizes; guard numbers; k! orderings; inclusion–exclusion table with a Stirling-number cross-check; pigeonhole bound |
| `genomath/counting_recursion.py` | `factorial`, `binomial` (multiplicative, traced steps), `permutations`, `sequence_space`, `growth_table`, `mutation_count` (closed form), `recurrence_table` (bottom-up DP), `mutation_count_recursive` (top-down, explicit memo, reports calls), `verify_recurrence`, `ball_size`, `generate_mutants` / `generate_within` (recursive, capped), `inclusion_exclusion_all_bases`, `forced_repeat`, `max_forced_k`, `find_repeated_kmer`, `mutation_orderings` (k! shortest substitution orders, hand-written permutation recursion) |
| `experiments/expected_unit2.py` | expected values transcribed from `unit2_expected.md` |
| `experiments/exp2_counting_enumeration.py` → `results/exp2.md` | Experiment 2 |
| `tests/test_counting.py` (21 tests); `tests/test_experiments.py` now covers Exp 1 and Exp 2 | |
| `scripts/mutation_check.py` | 8 Unit II mutants added (16 in total) |

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp2_counting_enumeration
.venv/Scripts/python scripts/mutation_check.py
```
Raw output from my last run:
```
88 passed in 4.37s
wrote results\exp1.md: 22/22 expected values matched
wrote results\exp2.md: 17/17 expected values matched
16/16 mutants caught
genomath: no banned imports
```
Oracles in `test_counting.py`:
- `math.comb`, `math.perm` and `math.factorial` for n ≤ 60;
- brute-force `itertools.product` enumeration of Σ^n, comparing the generated *sets*, not just counts, for 7 references and every k;
- brute-force inclusion–exclusion for n ≤ 7;
- NetworkX `all_shortest_paths` in the full Hamming graph H(4,4) (256 vertices) for the k! orderings;
- the hand table.

## Experiment 2 summary (from `results/exp2.md`)
17/17 expected values matched; 2 items are observation-only.
- **Agreement on all 45 cells.** For n ≤ 8, recursive enumeration, the closed form and the recurrence table all equal
  the hand table. All 45 cells contain distinct strings at exactly distance k, and row sums equal 4^n.
- **Guard.** M(12,6) = 673,596 and M(1701,2) = 13,012,650 are refused; M(1701,1) = 5103 is generated.
- **Inclusion–exclusion.** Agrees with enumeration for n = 4..8 (24, 240, 1560, 8400, 40824); n = 12 gives 14,676,024.
- **Orderings.** syn01 → syn04 has 3! = 6 shortest orderings; exactly 1 passes only through observed haplotypes.
- **Pigeonhole.** For L = 1701 the largest forced k is 5, and all 34 clean Dataset B sequences contain a repeated
  5-mer. Observation: all 34 also contain a repeated 6-mer, although the bound does not force it.
- **Sparsity** (observation). Dataset B degree sums are 16 at d = 1 and 62 at d ≤ 2, matching the QC script. The most
  connected haplotype has 6 observed one-substitution neighbours out of 5103 possible.
- **Enumeration cost at n = 8** (observation): recursive calls track the output size (e.g. k = 6: 20,412 strings,
  27,064 calls). Wall-clock times are listed but machine-dependent.

## Complexity of the Stage 2 algorithms
| Algorithm | Time | Space |
|---|---|---|
| `binomial` (multiplicative) | O(min(k, n−k)) big-integer steps | O(1) numbers |
| `recurrence_table(n)` | Θ(n²) | Θ(n²) |
| `mutation_count_recursive(n,k)` | O(n·k) distinct subproblems; recursion depth n (refused above n = 900 because of Python's recursion limit) | O(n·k) memo |
| `generate_mutants(ref,k)` | Θ(n·M(n,k)) character writes; calls ≈ output size | recursion depth k; Θ(n·M(n,k)) output |
| `generate_within(ref,r)` | Θ(n·\|B(n,r)\|) | Θ(n·\|B(n,r)\|) |
| `inclusion_exclusion_all_bases(n)` | 5 terms, each a big-integer power | O(1) |
| `find_repeated_kmer(seq,k)` | O(L·k) (slicing plus hashing each window) | O(L·k) |
| `mutation_orderings(s,t)` | counting O(n + k); listing Θ(k!·k) (only for k ≤ 8) | Θ(k!·n) when listed |

The exponential pieces are generation (guarded by the cap) and ordering enumeration (listed only for k ≤ 8; above that
the count k! is returned without a list).

## Problems found and fixed during the stage
- **One surviving mutant.** Changing the pigeonhole test from `>` to `>=` was not caught at first: no test used the tight
  case where windows = holes. I added `forced_repeat(4,1)` and `forced_repeat(17,2)` must be False ("ACGT" has 4
  windows, 4 holes and no repeat), and a consistency check of `max_forced_k` against `forced_repeat` for L = 1..2000.
  Now 16/16 mutants are caught.
- Experiment 2 had an unused variable and a redundant guard in one comprehension; both were removed before the final run.

## Open issues and doubts
- The **generator recurses on position choice (depth k)**, so it mirrors the closed form C(n,k)·3^k rather than the
  recurrence. A version that recursed on the last position would mirror the recurrence exactly, but would need depth n,
  which overflows Python's recursion limit for n = 1701. The recurrence is demonstrated by `recurrence_table` and
  `mutation_count_recursive` instead. Tell me if you want a recurrence-shaped generator for small n as a teaching aid.
- `math.comb` and the others are used only in tests. `genomath` computes them by hand.
- The **Stirling numbers** S(n,4) in `unit2_expected.md` §4 are standard table values that I recalled rather than looked
  up. They agree with the inclusion–exclusion results, and the inclusion–exclusion results agree with brute-force
  enumeration, so they are cross-checked but not sourced.
- **The "all 34 sequences contain a repeated 6-mer" result is an observation only.** It is not forced by the pigeonhole
  bound, and I make no biological reading of it.
- Still open from Stage 0, for later stages: k-mer analysis on Dataset B (Stage 4), the UPGMA tie-break (Stage 6) and the
  modulus set (Stage 4).
