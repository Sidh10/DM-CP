# Experiment 5: Mutation-network structure

- run date: 2026-10-31
- command: `python -m experiments.exp5_network_structure`
- Python 3.13.11 on Windows

## Question
How do edges, connected components, cycles and connectivity of the mutation-relationship graph change as the distance threshold k grows, and do the hand-computed graph quantities for Dataset A match?

## Input
Dataset A haplotypes (9 vertices); Dataset B haplotypes (30 vertices, n = 1701), thresholds [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 57].

## Mathematical model
G_=1: edge iff d_H = 1. G_k: edge iff 1 ≤ d_H ≤ k, weight d_H. Degrees and the handshake lemma; adjacency and incidence matrices; BFS / DFS; components; cyclomatic number |E| − |V| + c; bipartiteness; Whitney κ ≤ λ ≤ δ; smallest connecting threshold k*.

## Procedure
Build G_1, G_2 and the full sweep for Dataset A and compare every quantity with the hand derivations. Compute κ and λ by exhaustive search on three small graphs. For Dataset B, sweep the listed thresholds, compare edge counts with the genomath-free QC script, and report the remaining statistics as observations.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | A: edges of G_1 | [(01, 02), (01, 05), (02, 03), (02, 06), (03, 04), (05, 06)] | docs/unit5_expected.md §1 |
| 2 | A: adjacency matrix of G_1 | [[0, 1, 0, 0, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]] | docs/unit5_expected.md §1 |
| 3 | A: incidence matrix of G_1 | [[1, 1, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0], [0, 0, 1, 0, 1, 0], [0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1], [0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]] | docs/unit5_expected.md §1 |
| 4 | A: degrees in G_1 and G_2 | [[2, 3, 2, 1, 2, 2, 0, 0, 0], [6, 5, 4, 2, 4, 4, 1, 2, 0]] | docs/unit5_expected.md §1–2 |
| 5 | A: handshake Σ deg = 2\|E\| (G_1, G_2) | True | docs/unit5_expected.md §1–2 — *consistency check* |
| 6 | A: BFS from syn01 in G_1 | {order: [01, 02, 05, 03, 06, 04], hops: {01: 0, 02: 1, 05: 1, 03: 2, 06: 2, 04: 3}} | docs/unit5_expected.md §1 |
| 7 | A: DFS from syn01 in G_1 | {pre: [01, 02, 03, 04, 06, 05], post: [04, 03, 05, 06, 02, 01], back_edges: [(05, 01)]} | docs/unit5_expected.md §1 |
| 8 | A: cycle found in G_1 | [01, 02, 06, 05] | docs/unit5_expected.md §1 |
| 9 | A: cyclomatic number \|E\| − \|V\| + c | {1: 1, 2: 7} | docs/unit5_expected.md §1–2 |
| 10 | A: bipartite? | {1: True, 2: False} | docs/unit5_expected.md §1–2 |
| 11 | A: triangles in G_1 | 0 | docs/unit5_expected.md §1 |
| 12 | A: sweep (edges, components, largest, isolated, cyclomatic) | {0: [0, 9, 1, 9, 0], 1: [6, 4, 6, 3, 1], 2: [14, 2, 8, 1, 7], 3: [22, 2, 8, 1, 15], 4: [27, 2, 8, 1, 20], 5: [28, 2, 8, 1, 21], 6: [29, 1, 9, 0, 21], 9: [36, 1, 9, 0, 28]} | docs/unit5_expected.md §3 |
| 13 | A: smallest connecting threshold k* | 6 | docs/unit5_expected.md §3 |
| 14 | A: (κ, λ, δ) | {G_1: [0, 0, 0], square: [2, 2, 2], G_6: [1, 1, 1]} | docs/unit5_expected.md §4 |
| 15 | B: edges of G_k | {1: 8, 2: 31, 3: 57, 4: 90, 5: 113, 6: 148, 7: 178, 8: 202, 9: 223, 10: 227, 57: 435} | docs/dataset_b_qc.md (independent script) |
| 16 | B: components non-increasing in k | True | math_audit.md Unit V (G_j ⊆ G_k) — *consistency check* |
| 17 | B: smallest connecting threshold k* | *(none: reported as observation)* |  |
| 18 | B: G_=1 structure | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | A: edges of G_1 | [(01, 02), (01, 05), (02, 03), (02, 06), (03, 04), (05, 06)] |
| 2 | A: adjacency matrix of G_1 | [[0, 1, 0, 0, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]] |
| 3 | A: incidence matrix of G_1 | [[1, 1, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0], [0, 0, 1, 0, 1, 0], [0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1], [0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]] |
| 4 | A: degrees in G_1 and G_2 | [[2, 3, 2, 1, 2, 2, 0, 0, 0], [6, 5, 4, 2, 4, 4, 1, 2, 0]] |
| 5 | A: handshake Σ deg = 2\|E\| (G_1, G_2) | True |
| 6 | A: BFS from syn01 in G_1 | {order: [01, 02, 05, 03, 06, 04], hops: {01: 0, 02: 1, 05: 1, 03: 2, 06: 2, 04: 3}} |
| 7 | A: DFS from syn01 in G_1 | {pre: [01, 02, 03, 04, 06, 05], post: [04, 03, 05, 06, 02, 01], back_edges: [(05, 01)]} |
| 8 | A: cycle found in G_1 | [01, 02, 06, 05] |
| 9 | A: cyclomatic number \|E\| − \|V\| + c | {1: 1, 2: 7} |
| 10 | A: bipartite? | {1: True, 2: False} |
| 11 | A: triangles in G_1 | 0 |
| 12 | A: sweep (edges, components, largest, isolated, cyclomatic) | {0: [0, 9, 1, 9, 0], 1: [6, 4, 6, 3, 1], 2: [14, 2, 8, 1, 7], 3: [22, 2, 8, 1, 15], 4: [27, 2, 8, 1, 20], 5: [28, 2, 8, 1, 21], 6: [29, 1, 9, 0, 21], 9: [36, 1, 9, 0, 28]} |
| 13 | A: smallest connecting threshold k* | 6 |
| 14 | A: (κ, λ, δ) | {G_1: [0, 0, 0], square: [2, 2, 2], G_6: [1, 1, 1]} |
| 15 | B: edges of G_k | {1: 8, 2: 31, 3: 57, 4: 90, 5: 113, 6: 148, 7: 178, 8: 202, 9: 223, 10: 227, 57: 435} |
| 16 | B: components non-increasing in k | True |
| 17 | B: smallest connecting threshold k* | 47 |
| 18 | B: G_=1 structure | {components: 22, sizes: [9, 1, 1, 1, 1], max_degree: 6, cyclomatic: 0, triangles: 0, bipartite: True} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | A: edges of G_1 | yes |
| 2 | A: adjacency matrix of G_1 | yes |
| 3 | A: incidence matrix of G_1 | yes |
| 4 | A: degrees in G_1 and G_2 | yes |
| 5 | A: handshake Σ deg = 2\|E\| (G_1, G_2) | yes (consistency) |
| 6 | A: BFS from syn01 in G_1 | yes |
| 7 | A: DFS from syn01 in G_1 | yes |
| 8 | A: cycle found in G_1 | yes |
| 9 | A: cyclomatic number \|E\| − \|V\| + c | yes |
| 10 | A: bipartite? | yes |
| 11 | A: triangles in G_1 | yes |
| 12 | A: sweep (edges, components, largest, isolated, cyclomatic) | yes |
| 13 | A: smallest connecting threshold k* | yes |
| 14 | A: (κ, λ, δ) | yes |
| 15 | B: edges of G_k | yes |
| 16 | B: components non-increasing in k | yes (consistency) |
| 17 | B: smallest connecting threshold k* | n/a |
| 18 | B: G_=1 structure | n/a |

**14/14 independent expected values matched; 2/2 consistency checks passed; 2 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

### Dataset B threshold sweep (observation)

| k | edges | components | largest | isolated | max degree | cyclomatic |
|---|---|---|---|---|---|---|
| 0 | 0 | 30 | 1 | 30 | 0 | 0 |
| 1 | 8 | 22 | 9 | 21 | 6 | 0 |
| 2 | 31 | 16 | 13 | 13 | 9 | 17 |
| 3 | 57 | 12 | 15 | 8 | 14 | 39 |
| 4 | 90 | 10 | 17 | 6 | 15 | 70 |
| 5 | 113 | 8 | 19 | 4 | 18 | 91 |
| 6 | 148 | 6 | 21 | 2 | 20 | 124 |
| 7 | 178 | 6 | 21 | 2 | 20 | 154 |
| 8 | 202 | 5 | 22 | 1 | 21 | 177 |
| 9 | 223 | 5 | 22 | 1 | 21 | 198 |
| 10 | 227 | 5 | 22 | 1 | 21 | 202 |
| 15 | 238 | 4 | 22 | 0 | 21 | 212 |
| 20 | 244 | 3 | 22 | 0 | 21 | 217 |
| 30 | 284 | 2 | 28 | 0 | 25 | 256 |
| 40 | 379 | 2 | 28 | 0 | 27 | 351 |
| 50 | 397 | 1 | 30 | 0 | 29 | 368 |
| 57 | 435 | 1 | 30 | 0 | 29 | 406 |

## Interpretation
Dataset A behaves as designed: one 4-cycle at k = 1 (cyclomatic number 1), odd cycles appearing at k = 2 (the triangle 01–02–03), and the outlier syn09 isolated until k = 6. The non-transitivity of R_k from Unit I is exactly why components merge in stages rather than all at once. Dataset B (observation, no biological claim): G_=1 has 22 components (sizes, largest first: [9, 1, 1, 1, 1], …), maximum degree 6, cyclomatic number 0 and 0 triangles, and it is bipartite. The graph first becomes connected at k* = 47. Unit VI will check that this equals the largest edge of every minimum spanning tree.

## Limitations
Vertices are observed haplotypes only, so a missing edge can mean an unsampled intermediate. Any single k is arbitrary, which is why results are reported as a sweep. κ and λ are computed by exhaustive search and are restricted to graphs with at most 16 vertices. Graph structure describes this sample of 30 haplotypes, not the virus population.
