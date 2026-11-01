# Stage 6 report: Unit VI (trees), Experiment 6 completed, Experiment 7

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 7 (real-data evaluation, Exp 8, plus D9) has not
started.

Decisions followed:
- your Stage 0 rules: MST and UPGMA kept strictly separate; UPGMA's clock assumption stated; MST rooting marked
  arbitrary; optional counting link n^(n−2) and (2n−3)!!;
- **D8:** UPGMA ties resolved by the lexicographically smallest pair of sorted member tuples. A `reverse` option exists
  only to show sensitivity to the rule.

## What exists
| Path | Contents |
|---|---|
| `docs/unit6_expected.md` | hand-derived Kruskal trace (29 examined / 8 accepted / 21 rejected), Prim trace, validity, bottleneck = k\*, MST count 8 via weight levels and Kirchhoff, rooted statistics and traversals for roots 01 and 04, the full 8-step UPGMA with exact fractions and tie counts, Newick, in-order traversal, the ultrametric witness, the cophenetic fit, MST-vs-UPGMA bipartitions, and the counting link |
| `genomath/tree_analysis.py` | `merge_sort` (hand-written, stable), `UnionFind` (path compression, union by rank), `kruskal`, `prim` (O(V²)), `validate_spanning_tree`, `kirchhoff_count` (exact rational determinant), `count_minimum_spanning_trees`, `bottleneck`, `cayley_count`, `double_factorial`, `rooted_binary_topologies`, `RootedTree`, `root_tree`, `tree_statistics`, `preorder`, `postorder`, `level_order`, `inorder` (binary only), `edge_splits` |
| `genomath/phylogenetic_analysis.py` | `upgma` (exact `Fraction` arithmetic, D8 tie-break, merge table with tie counts, molecular-clock limitation in the trace), `branch_lengths`, `newick`, `cophenetic`, `is_ultrametric`, `fit_to_input`, `upgma_splits`, `compare_splits` |
| `experiments/expected_unit6.py`; `exp6_paths_and_trees.py` (MST half added) → `results/exp6.md`; `exp7_tree_analysis.py` → `results/exp7.md` | |
| `tests/test_trees.py` (14 tests); `test_experiments.py` covers Exp 1–7 | |
| `scripts/mutation_check.py` | 12 Unit VI mutants (57 in total) |

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp6_paths_and_trees
.venv/Scripts/python -m experiments.exp7_tree_analysis
.venv/Scripts/python scripts/mutation_check.py        # ~20 minutes; run it in the background
```
Raw output from my last run:
```
157 passed in 14.98s
wrote results\exp1.md: 20/20 independent expected values matched; 2/2 consistency checks passed; 3 observation-only items
wrote results\exp2.md: 16/16 independent expected values matched; 1/1 consistency checks passed; 2 observation-only items
wrote results\exp3.md: 29/29 independent expected values matched; 1/1 consistency checks passed; 1 observation-only items
wrote results\exp4.md: 23/23 independent expected values matched; 5/5 consistency checks passed; 4 observation-only items
wrote results\exp5.md: 14/14 independent expected values matched; 2/2 consistency checks passed; 2 observation-only items
wrote results\exp6.md: 14/14 independent expected values matched; 4/4 consistency checks passed; 4 observation-only items
wrote results\exp7.md: 11/11 independent expected values matched; 2/2 consistency checks passed; 5 observation-only items
57/57 mutants caught   (baseline: 157 passed)
```
Oracles in `test_trees.py`:
- NetworkX `minimum_spanning_tree` weights and `is_tree` (Datasets A and B);
- NetworkX `SpanningTreeIterator` to count MSTs on 25 random tie-heavy graphs;
- SymPy exact determinants (Laplacians and general matrices) and Cayley's formula for K_n, n ≤ 11;
- Python `sorted` for merge-sort order and stability;
- an **independent naive UPGMA** in the test that recomputes every cluster mean from scratch, compared step by step on
  25 random sequence sets;
- the theorems "UPGMA output is ultrametric", "UPGMA reproduces ultrametric input exactly", "every leaf-to-root path sums
  to the root height", and "largest MST edge = k\*" on A, B and 25 random sets.

## Experiment 6 (completed) summary
14/14 independent, 4/4 consistency, 4 observation-only.
- **Dataset A.** The Kruskal and Prim edge lists and orders match the hand derivation exactly; both totals are 15; both
  outputs are valid spanning trees; there are 8 MSTs; max MST edge 6 = k\*; Kirchhoff on K_9 gives 9^7.
- **Dataset B** (observation). MST weight **171**. The largest edge is **47 = k\***, so the bottleneck property holds on
  real data. There are **4** MSTs, and Prim and Kruskal return different ones, which ties allow.

## Experiment 7 summary
11/11 independent, 2/2 consistency, 5 observation-only.
- **Dataset A.**
  - Rooted statistics and traversals match for both roots: height 3 vs 4, and 6 vs 5 leaves on the same 8 edges.
  - All 8 UPGMA merges (distances and heights as exact fractions, tie counts), the Newick string, the 4 tie steps, the
    in-order traversal and the non-ultrametric witness (01, 02, 03) match.
  - UPGMA reproduces exactly 3/36 distances, with largest deviation 12/7 at 04–08. The MST and UPGMA share 1 of their
    2 vs 6 non-trivial bipartitions.
- **Dataset B** (observation).
  - The distances are not ultrametric, and UPGMA reproduces only 9 of 435 exactly.
  - 11 merge steps had ties, and the reverse tie rule changes the topology.
  - The MST and UPGMA share 5 of their 10 vs 27 non-trivial bipartitions.

## Complexity
| Algorithm | Time | Space |
|---|---|---|
| merge sort | Θ(E log E) | Θ(E) |
| Kruskal | Θ(E log E) sort + near-linear union–find | Θ(V + E) |
| Prim (array) | Θ(V²) | Θ(V) |
| validate tree | O(V + E·α(V)) | Θ(V) |
| Kirchhoff (exact) | O(V³) rational operations | Θ(V²) |
| MST count | O(E log E + Σ_levels V_level³) | Θ(V²) |
| root / traversals / statistics | Θ(V) | Θ(V) |
| UPGMA (naive) | Θ(N³) with exact fractions | Θ(N²) |
| cophenetic / fit | Θ(N²) | Θ(N²) |
| ultrametric test | Θ(N³) triples | O(1) |
| edge splits | O(V·E) | Θ(V) |

## Problems found and fixed during the stage
- **One of my hand values was wrong.** I first wrote the largest cophenetic deviation as 13/8 at 01–09. The
  implementation said 12/7 at 04–08. Re-deriving by hand showed I had missed d(04,08) = 5 at merge U7. I corrected the doc
  before the experiment used it, and recorded the correction there (`unit6_expected.md` §6).
- **Two surviving mutants exposed real gaps.**
  1. UPGMA heights were stored twice (`height` dict and `merges[...]["height"]`), and only the untested copy fed the
     Newick branch lengths. It is now a single source, with an exact `branch_lengths` function, hand values for Dataset A
     (e.g. U6 → U7 = 5/84, U7 → U8 = 243/112) and a leaf-to-root sum test.
  2. The determinant's row-swap branch never runs on Laplacian minors. They are positive semidefinite, so a zero pivot
     means its whole column below is zero. The mutant was therefore equivalent where the helper is used. A direct test
     on general matrices now covers the branch anyway.

  Both mutants are now caught.
- **Wording corrections.** "Not ultrametric, so UPGMA reproduces only 3" overstated causation and is now two facts.
  Newick branch lengths were mislabelled "in units of d_H/2"; they are in d_H units.
- The usual heredoc escaping slip broke two mutation entries once; I fixed them by direct edit.
- NetworkX's `number_of_spanning_trees` needs SciPy, which isn't installed. I used SymPy determinants instead rather than
  add a dependency.

## Open issues and doubts
- **Mutation-check time has grown to about 20 minutes** (57 mutants; some mutants make the experiment tests run
  longer before failing). Run it in the background. I could add a "quick" mode that runs only the tests for the mutated
  module, if you want.
- The Dataset B UPGMA tree has 11 tied merge steps, and the reverse rule changes the topology. Any Dataset B tree
  figure therefore depends on D8 and must say so.
- "Robinson–Foulds" is the usual name for the bipartition symmetric difference. I report the count without that name,
  because I have not located the original reference.
- Next, Stage 7: the Exp 8 real-data run across all units, including D9 (transition/transversion counts via the
  base-level Klein group).
