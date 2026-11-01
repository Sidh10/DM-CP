# Experiment 6: Shortest paths and spanning trees

- run date: 2026-10-31
- command: `python -m experiments.exp6_paths_and_trees`
- Python 3.13.11 on Windows

## Question
On the sparse threshold graph (edge iff d_H ≤ k, weight d_H), what does Dijkstra's shortest path cost between observed haplotypes, and when does it exceed the direct Hamming distance? Do Prim and Kruskal reproduce the hand-derived minimum spanning tree, how many minimum spanning trees exist, and is the largest MST edge equal to the smallest connecting threshold?

## Input
Dataset A, G_1, G_2 and K_9; Dataset B, G_1, G_k* with k* = 47 (the smallest threshold at which the graph is connected) and K_30. All edge weights are d_H.

## Mathematical model
Dijkstra's algorithm, O(V²) array version, non-negative weights d_H ≥ 1; all optimal paths via the predecessor DAG; triangle inequality cost ≥ d_H. Kruskal (merge sort + union–find with path compression and union by rank) and Prim (O(V²) array). MST count by weight-level contraction with Kirchhoff's matrix-tree theorem. Bottleneck (minimax) property: max MST edge = k*.

## Procedure
Compare the hand-derived Dataset A rows (cost, direct distance, number of optimal paths, full distance table from 04, unreachable pair in G_1). On Dataset B, run Dijkstra from every haplotype on G_k* and check cost ≥ d_H for all ordered pairs; report the detours from the medoid as observations. Run Kruskal and Prim on the complete distance graph of each dataset, compare Dataset A edge by edge with the hand derivation, validate both outputs as spanning trees, count the MSTs and test the bottleneck property.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | A, G_2: (cost, direct d_H, #optimal paths) | {04→08: [5, 5, 3], 07→03: [4, 3, 2], 07→08: [4, 4, 1]} | docs/synthetic_expected.md §4 |
| 2 | A, G_2: optimal paths 04 → 08 | [[04, 02, 01, 08], [04, 03, 01, 08], [04, 03, 02, 01, 08]] | docs/synthetic_expected.md §4 |
| 3 | A, G_2: Dijkstra distances from 04 | {01: 3, 02: 2, 03: 1, 04: 0, 05: 4, 06: 3, 07: 5, 08: 5, 09: None} | docs/unit5_expected.md §5 |
| 4 | A, G_2: targets where cost > direct d_H (detour) | {07: 1} | docs/unit5_expected.md §5 |
| 5 | A, G_1: cost 03 → 07 (unreachable) | None | docs/synthetic_expected.md §4 |
| 6 | B, G_k* (k* = 47): cost ≥ direct d_H for all ordered pairs | True | triangle inequality (math_audit.md Unit V) — *consistency check* |
| 7 | B, G_k*: detours from the medoid CY249987.1 (target: extra cost) | *(none: reported as observation)* |  |
| 8 | B, G_k*: largest detour from the medoid | *(none: reported as observation)* |  |
| 9 | B, G_1: haplotypes reachable from the medoid | *(none: reported as observation)* |  |
| 10 | A, K_9: Kruskal accepted edges | [(01, 02, 1), (01, 05, 1), (02, 03, 1), (02, 06, 1), (03, 04, 1), (01, 07, 2), (01, 08, 2), (01, 09, 6)] | docs/unit6_expected.md §1 |
| 11 | A, K_9: Kruskal examined / accepted / rejected | {examined: 29, accepted: 8, rejected: 21} | docs/unit6_expected.md §1 |
| 12 | A, K_9: Prim order from 01 | [01, 02, 03, 04, 05, 06, 07, 08, 09] | docs/unit6_expected.md §2 |
| 13 | A, K_9: Prim edges in order added | [(01, 02, 1), (02, 03, 1), (03, 04, 1), (01, 05, 1), (02, 06, 1), (01, 07, 2), (01, 08, 2), (01, 09, 6)] | docs/unit6_expected.md §2 |
| 14 | A, K_9: total weights | {prim: 15, kruskal: 15} | docs/unit6_expected.md §1–2 |
| 15 | A, K_9: both are spanning trees | {prim: True, kruskal: True} | docs/unit6_expected.md §3 |
| 16 | A, K_9: number of minimum spanning trees | 8 | docs/unit6_expected.md §4 |
| 17 | A: bottleneck = smallest connecting threshold | {max MST edge: 6, k*: 6} | docs/unit6_expected.md §3 |
| 18 | A: Kirchhoff spanning-tree count of K_9 (Cayley 9^7) | 4782969 | docs/unit6_expected.md §4 |
| 19 | B, K_30: Prim total = Kruskal total | True | MST weight is unique — *consistency check* |
| 20 | B, K_30: both outputs are spanning trees | True | definition — *consistency check* |
| 21 | B: max MST edge = k* (bottleneck property) | True | docs/unit6_expected.md §3 — *consistency check* |
| 22 | B, K_30: MST weight, number of MSTs, largest edge, Prim/Kruskal same edges | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | A, G_2: (cost, direct d_H, #optimal paths) | {04→08: [5, 5, 3], 07→03: [4, 3, 2], 07→08: [4, 4, 1]} |
| 2 | A, G_2: optimal paths 04 → 08 | [[04, 02, 01, 08], [04, 03, 01, 08], [04, 03, 02, 01, 08]] |
| 3 | A, G_2: Dijkstra distances from 04 | {01: 3, 02: 2, 03: 1, 04: 0, 05: 4, 06: 3, 07: 5, 08: 5, 09: None} |
| 4 | A, G_2: targets where cost > direct d_H (detour) | {07: 1} |
| 5 | A, G_1: cost 03 → 07 (unreachable) | None |
| 6 | B, G_k* (k* = 47): cost ≥ direct d_H for all ordered pairs | True |
| 7 | B, G_k*: detours from the medoid CY249987.1 (target: extra cost) | {MH885199.1: 2} |
| 8 | B, G_k*: largest detour from the medoid | {target: MH885199.1, cost: 51, direct: 49, hops: 3, optimal_paths: 1} |
| 9 | B, G_1: haplotypes reachable from the medoid | 8 |
| 10 | A, K_9: Kruskal accepted edges | [(01, 02, 1), (01, 05, 1), (02, 03, 1), (02, 06, 1), (03, 04, 1), (01, 07, 2), (01, 08, 2), (01, 09, 6)] |
| 11 | A, K_9: Kruskal examined / accepted / rejected | {examined: 29, accepted: 8, rejected: 21} |
| 12 | A, K_9: Prim order from 01 | [01, 02, 03, 04, 05, 06, 07, 08, 09] |
| 13 | A, K_9: Prim edges in order added | [(01, 02, 1), (02, 03, 1), (03, 04, 1), (01, 05, 1), (02, 06, 1), (01, 07, 2), (01, 08, 2), (01, 09, 6)] |
| 14 | A, K_9: total weights | {prim: 15, kruskal: 15} |
| 15 | A, K_9: both are spanning trees | {prim: True, kruskal: True} |
| 16 | A, K_9: number of minimum spanning trees | 8 |
| 17 | A: bottleneck = smallest connecting threshold | {max MST edge: 6, k*: 6} |
| 18 | A: Kirchhoff spanning-tree count of K_9 (Cayley 9^7) | 4782969 |
| 19 | B, K_30: Prim total = Kruskal total | True |
| 20 | B, K_30: both outputs are spanning trees | True |
| 21 | B: max MST edge = k* (bottleneck property) | True |
| 22 | B, K_30: MST weight, number of MSTs, largest edge, Prim/Kruskal same edges | {total: 171, count: 4, max_edge: 47, same_edges: False} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | A, G_2: (cost, direct d_H, #optimal paths) | yes |
| 2 | A, G_2: optimal paths 04 → 08 | yes |
| 3 | A, G_2: Dijkstra distances from 04 | yes |
| 4 | A, G_2: targets where cost > direct d_H (detour) | yes |
| 5 | A, G_1: cost 03 → 07 (unreachable) | yes |
| 6 | B, G_k* (k* = 47): cost ≥ direct d_H for all ordered pairs | yes (consistency) |
| 7 | B, G_k*: detours from the medoid CY249987.1 (target: extra cost) | n/a |
| 8 | B, G_k*: largest detour from the medoid | n/a |
| 9 | B, G_1: haplotypes reachable from the medoid | n/a |
| 10 | A, K_9: Kruskal accepted edges | yes |
| 11 | A, K_9: Kruskal examined / accepted / rejected | yes |
| 12 | A, K_9: Prim order from 01 | yes |
| 13 | A, K_9: Prim edges in order added | yes |
| 14 | A, K_9: total weights | yes |
| 15 | A, K_9: both are spanning trees | yes |
| 16 | A, K_9: number of minimum spanning trees | yes |
| 17 | A: bottleneck = smallest connecting threshold | yes |
| 18 | A: Kirchhoff spanning-tree count of K_9 (Cayley 9^7) | yes |
| 19 | B, K_30: Prim total = Kruskal total | yes (consistency) |
| 20 | B, K_30: both outputs are spanning trees | yes (consistency) |
| 21 | B: max MST edge = k* (bottleneck property) | yes (consistency) |
| 22 | B, K_30: MST weight, number of MSTs, largest edge, Prim/Kruskal same edges | n/a |

**14/14 independent expected values matched; 4/4 consistency checks passed; 4 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

## Interpretation
On Dataset A a shortest path through observed haplotypes costs more than the direct distance only for 04 → 07, where the route must pass through syn01 and pay for position 5 twice. Ties are real: 04 → 08 has three optimal paths. On Dataset B (observation) 1 of the other 29 haplotypes needs a detour from the medoid CY249987.1 at k* = 47; the largest is to MH885199.1 (cost 51 vs d_H 49, 3 hops, 1 optimal path(s)). In G_1 only 8 haplotypes are reachable from the medoid at all. Prim and Kruskal both reach weight 15 on Dataset A; 8 different minimum spanning trees exist, so a test can require equal weight and validity but not equal edge sets. The largest MST edge equals k* on both datasets: on Dataset B the MST weighs 171, its largest edge is 47, and 4 minimum spanning tree(s) exist (observation).

## Limitations
Paths run through observed haplotypes only and are not inferred evolutionary trajectories. The weight d_H counts substitutions without a substitution model. k* is a data-driven choice of threshold; other thresholds give other paths. Dijkstra here returns one optimal path, plus all optimal paths when asked. A minimum spanning tree minimises total edge weight among trees on the observed sequences; it is not a phylogeny and does not establish ancestry (compared with UPGMA in Experiment 7).
