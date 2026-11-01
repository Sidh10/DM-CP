# Unit V: hand-derived graph results for Dataset A

Status: **expected values only**, worked by hand from the distance matrix in `synthetic_expected.md` §3; nothing here
was produced by `genomath`. Vertices are the 9 haplotypes, written 01..09 for syn01..syn09.
Neighbours are always visited in increasing id order.

**Edge rules.**
- G_=1 (single-substitution graph): edge {S,T} iff d_H(S,T) = 1.
- G_k (threshold graph): edge iff 1 ≤ d_H(S,T) ≤ k, with weight d_H.

On haplotypes d_H ≥ 1 for distinct vertices, so G_=1 = G_1.

## 1. G_1
Edges, in lexicographic order: e1 = 01–02, e2 = 01–05, e3 = 02–03, e4 = 02–06, e5 = 03–04, e6 = 05–06.

**Adjacency matrix** (rows and columns 01..09):
```
     01 02 03 04 05 06 07 08 09
01    0  1  0  0  1  0  0  0  0
02    1  0  1  0  0  1  0  0  0
03    0  1  0  1  0  0  0  0  0
04    0  0  1  0  0  0  0  0  0
05    1  0  0  0  0  1  0  0  0
06    0  1  0  0  1  0  0  0  0
07    0  0  0  0  0  0  0  0  0
08    0  0  0  0  0  0  0  0  0
09    0  0  0  0  0  0  0  0  0
```
**Incidence matrix** (rows 01..09, columns e1..e6):
```
     e1 e2 e3 e4 e5 e6
01    1  1  0  0  0  0
02    1  0  1  1  0  0
03    0  0  1  0  1  0
04    0  0  0  0  1  0
05    0  1  0  0  0  1
06    0  0  0  1  0  1
07..09 all 0
```
Every column has exactly two 1s. The row sums are the degrees: 2, 3, 2, 1, 2, 2, 0, 0, 0, which sum to 12 = 2|E| (handshake).

**BFS from 01:**
- dequeue 01, enqueue 02 and 05;
- dequeue 02, enqueue 03 and 06;
- dequeue 05 (06 is already seen);
- dequeue 03, enqueue 04;
- dequeue 06;
- dequeue 04.

Visit order **01, 02, 05, 03, 06, 04**. Hop distances: 01:0, 02:1, 05:1, 03:2, 06:2, 04:3.
Parents: 02←01, 05←01, 03←02, 06←02, 04←03. Unreached: 07, 08, 09.

**DFS from 01** (recursive):
- 01 → 02 → 03 → 04 (dead end; back up to 03, then to 02);
- 02 → 06 → 05. From 05, neighbour 01 is visited and is not 05's parent, so 05–01 is a **back edge**.

Preorder **01, 02, 03, 04, 06, 05**. Postorder **04, 03, 05, 06, 02, 01**.
Tree edges: 01–02, 02–03, 03–04, 02–06, 06–05. The back edge closes the cycle **01–02–06–05–01** (length 4).

**Components:** {01,02,03,04,05,06}, {07}, {08}, {09}.
Cyclomatic number |E| − |V| + c = 6 − 9 + 4 = **1** (exactly one independent cycle).
**Bipartite:** yes. The only cycle has even length; the 2-colouring is {01,03,06} vs {02,04,05}.
Triangles: 0. A triangle in G_=1 would need three haplotypes with three different bases at one position and identical
elsewhere, and there are none.

## 2. G_2
Edges with d ≤ 2: the six above, plus 01–03, 01–06, 01–07, 01–08, 02–04, 02–05, 03–06 and 05–08, giving **14**.
Degrees: 01:6, 02:5, 03:4, 04:2, 05:4, 06:4, 07:1, 08:2, 09:0. They sum to 28 = 2·14.
Components: {01..08} and {09}. Cyclomatic number 14 − 9 + 2 = **7**.
**Not bipartite:** 01–02–03 is a triangle (d = 1, 1, 2), an odd cycle.

## 3. Threshold sweep
| k | \|E\| | components | largest component | isolated vertices | cyclomatic \|E\|−\|V\|+c |
|---|---|---|---|---|---|
| 0 | 0 | 9 | 1 | 9 | 0 |
| 1 | 6 | 4 | 6 | 3 (07, 08, 09) | 1 |
| 2 | 14 | 2 | 8 | 1 (09) | 7 |
| 3 | 22 | 2 | 8 | 1 | 15 |
| 4 | 27 | 2 | 8 | 1 | 20 |
| 5 | 28 | 2 | 8 | 1 | 21 |
| 6 | 29 | 1 | 9 | 0 | 21 |
| 9 | 36 | 1 | 9 | 0 | 28 |

The **smallest connecting threshold** is **k\* = 6** (syn09's nearest haplotype is at distance 6).
The number of components is non-increasing in k because G_j ⊆ G_k for j ≤ k.

## 4. Connectivity: κ (vertex), λ (edge), δ (minimum degree)
- G_1 (disconnected): κ = λ = 0, δ = 0.
- The square, G_1 induced on {01,02,05,06}, is a 4-cycle: κ = 2, λ = 2, δ = 2.
  Removing one vertex leaves a path; removing 02 and 05 disconnects 01 from 06.
- G_6 (connected): syn09's only neighbour is 01 (distance 6; every other distance to 09 is ≥ 7). So δ = 1, λ = 1
  (01–09 is a bridge), and κ = 1 (01 is a cut vertex).

All three satisfy Whitney's κ ≤ λ ≤ δ.

## 5. Dijkstra on G_2 (sparse threshold graph, weight d_H)
Single-pair rows are in `synthetic_expected.md` §4:
- 04→08: cost 5, direct 5, 3 optimal paths.
- 07→03: cost 4, direct 3, 2 optimal paths.
- 07→08: cost 4, direct 4, 1 optimal path.
- In G_1, 07 is unreachable from 03.

**All distances from 04 in G_2**, settling vertices in order of distance:
- 04: 0
- 03: 1 (edge 04–03)
- 02: 2 (04–02 directly, or 04–03–02 = 1 + 1)
- 01: 3 (04–03–02–01, 04–02–01 or 04–03–01)
- 06: 3 (04–03–02–06, 04–02–06 or 04–03–06)
- 05: 4 (via 01: 3 + 1; via 06: 3 + 1; via 02: 2 + 2)
- 07: 5 (only via 01: 3 + 2)
- 08: 5 (via 01: 3 + 2 = 5, beating via 05: 4 + 2 = 6)
- 09: unreachable

Direct d_H from 04: 01:3, 02:2, 03:1, 05:4, 06:3, 07:4, 08:5.
Path cost exceeds direct d_H for **07 only** (5 vs 4, detour 1). Every other reachable vertex has cost = d_H.

**Why cost ≥ d_H always:** d_H is a metric, so the triangle inequality makes any path's total weight at least d_H
between its endpoints. Equality needs a path that never changes the same position twice.
