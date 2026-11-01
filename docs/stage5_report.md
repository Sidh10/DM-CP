# Stage 5 report: Unit V (graph theory), Experiment 5 and the Dijkstra half of Experiment 6

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 6 has not started.
No new decisions were needed. This stage follows your Stage 0 rule: Dijkstra runs on the sparse threshold graph (edge iff
d_H ≤ k, weight d_H) and reports path cost vs direct d_H. It also follows D3: components = classes of the transitive
closure of R_k (tested since Stage 1).

## What exists
| Path | Contents |
|---|---|
| `docs/unit5_expected.md` | hand-derived G_1 adjacency and incidence matrices, degrees, BFS/DFS orders, back edge and cycle, cyclomatic numbers, bipartiteness, G_2 degrees, the full threshold sweep, k\* = 6, three κ/λ/δ examples, and all Dijkstra distances from 04 in G_2 |
| `genomath/mutation_graph.py` | `Graph`; `distance_matrix`, `build_graph` (modes `threshold` and `single`, with the edge rule written into the trace), `graph_from_matrix`, `degrees` (handshake check), `adjacency_matrix`, `incidence_matrix`, `bfs` (list + head-pointer queue), `dfs` (pre/post order, back edges), `connected_components`, `find_cycle`, `cyclomatic_number`, `is_bipartite`, `count_triangles`, `induced_subgraph`, `complement`, `edge_set_operation` (union / intersection / difference), `vertex_connectivity` and `edge_connectivity` (exhaustive, refused above 16 vertices), `dijkstra` (O(V²) array version), `path_from`, `all_shortest_paths`, `path_versus_direct`, `threshold_sweep`, `smallest_connecting_threshold` |
| `experiments/expected_unit5.py`, `experiments/exp5_network_structure.py` → `results/exp5.md` | Experiment 5 |
| `experiments/exp6_paths_and_trees.py` → `results/exp6.md` | Experiment 6, Dijkstra half. The file states that Prim/Kruskal are added in Stage 6 |
| `tests/test_graphs.py` (12 tests); `test_experiments.py` covers Exp 1–6 | |
| `scripts/mutation_check.py` | 10 Unit V mutants (45 in total). Mutants now run with `pytest -x`; the baseline still runs the full suite |

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp5_network_structure
.venv/Scripts/python -m experiments.exp6_paths_and_trees
.venv/Scripts/python scripts/mutation_check.py        # ~6 minutes; run it in the background
```
Raw output from my last run:
```
142 passed in 17.27s
wrote results\exp1.md: 20/20 independent expected values matched; 2/2 consistency checks passed; 3 observation-only items
wrote results\exp2.md: 16/16 independent expected values matched; 1/1 consistency checks passed; 2 observation-only items
wrote results\exp3.md: 29/29 independent expected values matched; 1/1 consistency checks passed; 1 observation-only items
wrote results\exp4.md: 23/23 independent expected values matched; 5/5 consistency checks passed; 4 observation-only items
wrote results\exp5.md: 14/14 independent expected values matched; 2/2 consistency checks passed; 2 observation-only items
wrote results\exp6.md: 5/5 independent expected values matched; 1/1 consistency checks passed; 3 observation-only items
45/45 mutants caught   (baseline: 142 passed)
```
NetworkX oracles in `test_graphs.py`:
- **On Dataset A for k = 0..9 and on Dataset B for k ∈ {1, 2, 5, 9, 20, 57}:** components, degrees, bipartiteness,
  triangles, cycle basis size vs cyclomatic number, BFS hop distances, Dijkstra distances.
- **On 60 random 7-vertex weighted graphs:** κ, λ, forest detection, the cycle returned, back edges = cyclomatic number,
  all shortest paths, and Whitney.
- **On 20 relabelled barbell graphs:** λ = κ = 1 < δ = 3.

Dataset B edge counts are also checked against the genomath-free QC script.

## Experiment 5 summary
14/14 independent, 2/2 consistency, 2 observation-only.
- **Dataset A.** All hand values match: matrices, degrees, BFS order 01, 02, 05, 03, 06, 04; DFS back edge 05–01 and cycle
  01–02–06–05; cyclomatic numbers 1 and 7; G_1 bipartite, G_2 not (triangle 01–02–03); the sweep; k\* = 6; and (κ, λ, δ)
  = (0,0,0), (2,2,2), (1,1,1).
- **Dataset B** (observation). G_=1 has 8 edges in 22 components, and **it is a forest** (cyclomatic number 0, no
  triangles, bipartite). The largest component has 9 haplotypes, around one haplotype of degree 6. The graph becomes
  connected only at **k\* = 47**. At k = 10 there are 5 components, and at k = 30 there are 2. The full sweep table is in
  `results/exp5.md`.

## Experiment 6 (Dijkstra half) summary
5/5 independent, 1/1 consistency, 3 observation-only.
- **Dataset A.** Rows 04→08 (5 = 5, 3 optimal paths), 07→03 (4 > 3) and 07→08 (4 = 4) match. All distances from 04 match,
  and the only detour target is 07. 03→07 is unreachable in G_1.
- **Dataset B** (observation). At k\* = 47, cost ≥ d_H holds for all 900 ordered pairs. From the medoid CY249987.1,
  exactly one haplotype needs a detour: MH885199.1, with cost 51 vs d_H 49 over 3 hops. In G_1 only 8 haplotypes are
  reachable from the medoid.

## Complexity
| Algorithm | Time | Space |
|---|---|---|
| distance matrix | Θ(N²·n) | Θ(N²) |
| graph from matrix | Θ(N²) | O(N + E) |
| degrees, matrices | O(N + E), Θ(N²), Θ(N·E) | same |
| BFS, DFS, components, find_cycle, bipartite | O(N + E) | O(N) |
| triangles | O(E·Δ) | O(1) |
| κ / λ exhaustive | O(2^N·(N + E)) / O(C(E, δ)·(N + E)) | O(N) — refused above 16 vertices |
| Dijkstra (array) | O(N² + E) | O(N) |
| all shortest paths | O(N + E) to build the predecessor DAG, plus output size (capped at 100 paths) | |
| threshold sweep over t thresholds | Θ(N²·n) once, then O(t·N²) | |

## Problems found and fixed during the stage
- **A DFS bug, before any test ran.** My first back-edge condition compared vertex indices, which is wrong because an
  ancestor can have a larger index. The fix tests "visited, not the parent, not yet finished". The random-graph test
  "back edges = cyclomatic number" now guards it, and the old form is a mutant.
- **Wasteful recomputation.** `threshold_sweep` and `smallest_connecting_threshold` recomputed the 1701-nt distance matrix
  for every threshold, which made the graph tests take 30 s. `graph_from_matrix` now reuses one matrix.
- **One random 8-vertex dense graph** made the exhaustive λ search dominate the test time (16 s). With 7 vertices the
  test runs in 0.5 s and keeps the same checks.
- **A surviving mutant revealed a test gap.** Breaking the edge-orientation lookup in the edge-removal check went
  unnoticed. When no smaller cut is found, λ falls back to δ, and the random graphs almost never had λ < δ with the cut
  edge in the "wrong" orientation. The 20 relabelled barbell graphs now catch it.
- Three candidate mutants were dropped as **equivalent**:
  - re-relaxing settled vertices in Dijkstra (impossible to improve with non-negative weights);
  - `<=` instead of `==` for shortest-path predecessors (dist[u] + w ≥ dist[v] always);
  - shortening the λ loop (the function returns δ anyway).
- Two wording slips in the experiment interpretations were corrected.

## Open issues and doubts
- **k\* = 47 for Dataset B is large compared with typical within-cluster distances**, because the distance distribution is
  multimodal (gaps around 22–27 and 37–47 nt). This is reported as an observation. I make no clade or lineage claim,
  because I have no located source for clade assignments.
- The suite takes about 17–31 s depending on machine load. The mutation check takes about 6 minutes with `-x`.
- κ and λ are exhaustive and limited to 16 vertices, so they are demonstrated on Dataset A graphs only, as planned in the audit.
- The DFS is recursive (depth ≤ N). That is fine for N = 30; a 1000-vertex graph would hit Python's recursion limit.
  `find_cycle` is iterative.
- Experiment 6 is half-complete by design. Stage 6 adds Prim/Kruskal and the bottleneck check k\* = max MST edge
  (A: 6; B: expected 47).
- Still open: nothing new. D8 (UPGMA tie-break) and D9 (transitions/transversions in Stage 7) are scheduled.
