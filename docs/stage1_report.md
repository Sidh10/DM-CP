# Stage 1 report: Unit I (sequence processing, set logic), trace type, Experiment 1

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 2 has not started.

## Decisions applied provisionally
You said "go stage 1" without ruling on the open Stage 0 decisions. I applied my recommendations **only for the ones Stage 1
needs**. Each is a local change if you want it different:
1. Trace shape: `Trace` (a list) of `{"kind", "text", "data"}`, with `kind` restricted to
   `observation | calculation | algorithm | interpretation | limitation` (`genomath/trace.py`).
2. Allele-level sets M* are added alongside position sets M.
3. Warshall transitive closure of R_k is added as the Unit I ↔ V bridge.
4. Dataset B reference = medoid haplotype (minimum Σ d_H; ties broken by smallest id). Result: **CY249987.1**.

Not done: the optional (M*, ⊆) Hasse diagram. The partial-order *property check* on (M*, ⊆) is tested
(`test_subset_order_on_mutation_sets_is_partial_order`), but nothing draws the diagram.
Still open for later stages: k-mer analysis on B, UPGMA tie-break, enumeration cap, modulus set.

## What exists
| Path | Contents |
|---|---|
| `genomath/trace.py` | `Trace.add(kind, text, **data)`; rejects unknown kinds |
| `genomath/_core.py` | untraced primitives: `mismatch_positions`, `mismatch_count`, `complement`, `reverse_complement`, `SequenceError` (carries the trace so far) |
| `genomath/sequence_processing.py` | `parse_fasta`, `read_fasta`, `validate`, `collapse_haplotypes`, `hamming`, `position_map`, `choose_reference`, `build_dataset`, `load_dataset` |
| `genomath/set_logic.py` | `mutation_set`, `allele_mutation_set`, `reference_table`, `set_operations`, `distance_decomposition`; `Predicate` with ∧ ∨ ¬ and `within_distance`, `base_at`, `carries_mutation`, `shares_mutation_with`, `filter_sequences`; `threshold_relation`, `identity_relation`, `orientation_relation`, `relation_properties` (with witnesses), `transitive_closure` (Warshall), `equivalence_classes`; `function_properties` (injective/surjective with witnesses) |
| `experiments/expected_dataset_a.py` | expected values transcribed from `synthetic_expected.md`, with the section noted on each |
| `experiments/_report.py` | writes `results/expN.md` with separate Expected / Actual / Comparison tables |
| `experiments/exp1_mutation_sets.py` → `results/exp1.md` | Experiment 1 |
| `tests/conftest.py`, `test_sequences.py`, `test_set_logic.py`, `test_experiments.py`; `test_skeleton.py` updated | 66 tests |
| `scripts/mutation_check.py` | injects 8 deliberate bugs one at a time and confirms the suite catches each |

Two modules were added beyond the §7 names (`trace.py`, `_core.py`), plus `experiments/__init__.py` (needed for
`python -m experiments.exp1_mutation_sets`).

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp1_mutation_sets
.venv/Scripts/python scripts/mutation_check.py
```
Raw output from my last run:
```
66 passed in 2.22s
wrote results\exp1.md: 22/22 expected values matched
baseline: 66 passed in 2.36s
CAUGHT   Warshall skips the last pivot                    2 failed, 64 passed
CAUGHT   transitivity tests (c,a) instead of (a,c)        3 failed, 63 passed
CAUGHT   d_H identity drops the shared-differing term     4 failed, 62 passed
CAUGHT   B-only elements counted as A-only                3 failed, 63 passed
CAUGHT   shared mutation ignores the base                 1 failed, 65 passed
CAUGHT   no uppercasing                                   1 failed, 65 passed
CAUGHT   medoid tie broken by sequence, not id            1 failed, 65 passed
CAUGHT   mismatch count skips positions where S has A     12 failed, 54 passed
8/8 mutants caught
```
Independent oracles used in the tests: Biopython (FASTA parsing, reverse complement), NetworkX (transitive closure,
connected components), Python's built-in set operators, brute-force definitions over all 512 relations on 3 elements,
the `synthetic_expected.md` hand values, and the genomath-free QC script for Dataset B counts.

## Experiment 1 summary (from `results/exp1.md`)
22/22 expected values matched; 3 items are observation-only.
- Dataset A: all 9 allele-level mutation sets; set operations; the d_H decomposition for syn03/syn07 (2 + 1 = 3);
  all 81 matrix entries; R_1 reflexive and symmetric but not transitive, with counterexample (syn01, syn02, syn03);
  components 9/4/2/2/1 for k = 0/1/2/5/6; predicate filters; orientation class {syn03, syn08, syn10}; bad01 dropped;
  bad02 rejected for length.
- Dataset B: 2 records dropped, 30 haplotypes, the d_H identity holds on 435/435 pairs, |R_1| = 46.
  Observations without an expected value:
  - R_1 is also non-transitive on real data (CY249987.1, MH637643.1, MH082676.1).
  - The medoid is CY249987.1, a haplotype shared by 3 records.
  - Mutation-set sizes relative to the medoid range from 0 to 49.

## Complexity of the Stage 1 algorithms (N sequences, length n)
| Algorithm | Time | Space |
|---|---|---|
| parse / validate | O(total characters) | O(total characters) |
| mismatch positions, d_H, M, M* | O(n) | O(|M|) |
| collapse to haplotypes | O(N·n) expected (string hashing) | O(N·n) |
| medoid reference | O(N²·n) | O(N) |
| set operations | O(m log m), m = \|A\| + \|B\| (sort for a deterministic trace) | O(m) |
| threshold relation R_k | O(N²·n) | O(N²) |
| relation properties | O(N³) worst case (transitivity over successor lists) | O(N²) |
| Warshall closure | Θ(N³) | O(N²) |
| equivalence classes | O(N³) (property check) + O(N²) | O(N²) |

## Problems found and fixed during the stage
- **Two weak tests.** On the first mutation run, two of four bugs survived (Warshall skipping the last pivot;
  transitivity testing (c,a)). Cause: every relation I had tested was symmetric, and on symmetric relations those bugs
  are invisible. I added exhaustive tests over all 512 relations on 3 elements, non-symmetric relations (divisibility,
  a chain), and a closure case that needs the last pivot.
- **My oracle was misused.** `nx.transitive_closure(reflexive=None)` never adds self-loops, so it disagreed with the
  correct closure on 252 of 512 relations (e.g. {(0,1),(1,0)} must yield (0,0)). The correct option is `reflexive=False`,
  which agreed on all 512.
- **A bug in the Experiment 1 script** (not in `genomath`). It counted every single-member orientation class as RC-fixed,
  but a singleton class can just mean the RC partner is absent. The Expected/Actual comparison caught it (21/22).
  It now tests RC(S) = S directly.
- **One equivalent mutant.** "Ignore N in the mismatch count" survived because validation guarantees N never reaches
  that function. I replaced it with a mutant that changes behaviour.
- File handles left open by `SeqIO.parse(path)` in a test and in the QC script (ResourceWarning under `-W error`) were
  closed. The QC output was confirmed byte-identical afterwards.

## Open issues and doubts
- **Exception to "(result, trace)":** the primitives in `_core.py` return plain values because they run inside O(N²)
  loops. Every public function in `sequence_processing` and `set_logic` returns `(result, trace)`.
  Please confirm this exception is acceptable.
- `Predicate` objects are callables returning `bool`, not `(result, trace)`. `filter_sequences` produces the traced
  truth table.
- Set operations iterate explicitly, but membership tests use Python's `frozenset` (`x in A`). I consider membership the
  definition, not a hidden algorithm. Say if you want sorted-list merging instead.
- Trace size: `mutation_set` emits one step per differing position (at most 49 on Dataset B), and `set_operations` one
  per element. This is fine for N = 30. A UI might need to collapse long traces.
- `reverse_complement` lives in `_core` for the orientation relation. Stage 3 will build the traced group structure on
  top of it rather than duplicate it.
- The interpretation text in `results/exp1.md` is partly fixed prose about Dataset A's design. The parts about Dataset B
  are generated from the run's values.
