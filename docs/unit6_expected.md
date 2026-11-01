# Unit VI: hand-derived tree results for Dataset A

Status: **expected values only**, worked by hand from the distance matrix in `synthetic_expected.md` §3; nothing here
was produced by `genomath`. Vertices are written 01..09 for syn01..syn09. The input graph for both MST algorithms is the
complete graph K_9 with weights d_H. The MST and UPGMA are kept strictly separate (your Stage 0 decision).

## 1. Kruskal (edges sorted by (weight, first endpoint, second endpoint); union–find; stop at 8 accepted)
- **Weight 1:** accept 01–02, 01–05, 02–03, 02–06, 03–04; **reject 05–06**, because 05 and 06 are already joined via 01–02–06.
- **Weight 2:** reject 01–03 and 01–06 (same component); **accept 01–07** and **accept 01–08**; reject 02–04, 02–05,
  03–06 and 05–08.
- **Weights 3, 4, 5:** all 8 + 5 + 1 = 14 edges are rejected, because only 09 is outside the tree and it is ≥ 6 away from everything.
- **Weight 6:** **accept 01–09**. That is the 8th edge, so the algorithm stops.

Edges examined: 6 + 8 + 8 + 5 + 1 + 1 = **29**; accepted **8**; rejected **21**.
Kruskal tree: {01–02, 01–05, 02–03, 02–06, 03–04, 01–07, 01–08, 01–09}, total weight 5·1 + 2·2 + 6 = **15**.

## 2. Prim from 01 (O(V²): key[v] = cheapest edge to the tree; pick the unvisited v with the smallest (key, index); update only on strict improvement)
1. Start 01. Keys: 02:1, 03:2, 04:3, 05:1, 06:2, 07:2, 08:2, 09:6 (all via 01).
2. Add **02** (key 1, via 01). Improvements: 03 → 1 (via 02), 04 → 2 (via 02), 06 → 1 (via 02).
3. Add **03** (key 1, via 02; it beats 05 and 06 on index). Improvement: 04 → 1 (via 03).
4. Add **04** (1, via 03). No improvements.
5. Add **05** (1, via 01). 06 stays at 1 (d(05,06) = 1 is not *strictly* smaller), and 08 stays at 2.
6. Add **06** (1, via 02).
7. Add **07** (2, via 01; it beats 08 on index).
8. Add **08** (2, via 01).
9. Add **09** (6, via 01).

Prim order: 01, 02, 03, 04, 05, 06, 07, 08, 09. Edges in the order added: 01–02, 02–03, 03–04, 01–05, 02–06, 01–07,
01–08, 01–09. Total **15**: the same edge set as Kruskal, by coincidence of the tie rules (8 MSTs exist, §4).

## 3. Validity and the bottleneck
Both trees have 8 = 9 − 1 edges, are acyclic and span all 9 vertices.
The largest MST edge is **6**, which equals the smallest connecting threshold k\* = 6 from `unit5_expected.md` §3
(bottleneck property: every MST contains a minimax path between any two vertices).

## 4. Number of MSTs (weight-level contraction; spanning trees per level by Kirchhoff)
- **Weight 1:** the unit edges form one component {01..06}, whose only cycle is the square 01–02–06–05. The spanning
  trees of the unit graph number **4**, by Kirchhoff on a 4-cycle (or directly: drop one of its 4 edges).
- **Weight 2:** the super-vertices are C = {01..06}, {07}, {08}, {09}. The crossing weight-2 edges are C–07 (01–07),
  and C–08 twice (01–08, 05–08). Spanning trees: 1 × 2 = **2**.
- **Weight 6:** C′–09 via 01–09 only, giving **1**.

Number of MSTs = 4 · 2 · 1 = **8** (matches the brute-force count in `synthetic_expected.md` §5).
Kirchhoff on K_9 (all spanning trees, ignoring weights) gives 9^7 = **4,782,969**, Cayley's formula.

## 5. The MST as a rooted tree (the root is arbitrary)
**Rooted at 01:**
- Children: 01 → [02, 05, 07, 08, 09]; 02 → [03, 06]; 03 → [04]. All other vertices are leaves.
- Depth: 01:0; 02, 05, 07, 08, 09:1; 03, 06:2; 04:3. **Height 3.**
- **Leaves (6):** 04, 05, 06, 07, 08, 09. **Internal (3):** 01, 02, 03.
- **Preorder:** 01, 02, 03, 04, 06, 05, 07, 08, 09.
- **Postorder:** 04, 03, 06, 02, 05, 07, 08, 09, 01.
- **Level order:** 01, 02, 05, 07, 08, 09, 03, 06, 04.

**Rooted at 04** (same tree, different root):
- Children: 04 → [03]; 03 → [02]; 02 → [01, 06]; 01 → [05, 07, 08, 09].
- Depth: 04:0, 03:1, 02:2, 01:3, 06:3, 05/07/08/09:4. **Height 4.**
- **Leaves (5):** 06, 05, 07, 08, 09. **Internal (4):** 04, 03, 02, 01.
- **Preorder:** 04, 03, 02, 01, 05, 07, 08, 09, 06.
- **Postorder:** 05, 07, 08, 09, 01, 06, 02, 03, 04.

Rooting changes the depth, height and leaf statistics but not the tree (the same 8 edges). This is why "rooting the
MST is arbitrary".

## 6. UPGMA (average linkage; tie-break D8: the tied pair whose sorted member tuples come first lexicographically)
Distances between clusters are averages over all member pairs. Heights are distance / 2. Exact fractions are used throughout.

| step | merge | distance | height | tied pairs at this minimum |
|---|---|---|---|---|
| 1 | U1 = (01, 02) | 1 | 1/2 | 6 |
| 2 | U2 = (03, 04) | 1 | 1/2 | 2 ((03,04), (05,06)) |
| 3 | U3 = (05, 06) | 1 | 1/2 | 1 |
| 4 | U4 = (U1, U3) = {01,02,05,06} | 3/2 | 3/4 | 1 |
| 5 | U5 = (U4, U2) = {01..06} | 5/2 | 5/4 | 2 ((U4,U2), (U4,08)) |
| 6 | U6 = (U5, 07) | 19/6 | 19/12 | 2 ((U5,07), (U5,08)) |
| 7 | U7 = (U6, 08) | 23/7 | 23/14 | 1 |
| 8 | U8 = (U7, 09) | 61/8 | 61/16 | 1 |

Working for the non-trivial distances (sums of the original d_H over member pairs):
- d(U1, U3) = (1 + 2 + 2 + 1)/4 = 3/2.
- d(U4, U2) = 20/8 = 5/2, and d(U4, 08) = 10/4 = 5/2. These tie; D8 picks (U4, U2) because the second member tuple
  (03,04) < (08).
- d(U5, 07) = (2+3+3+4+3+4)/6 = 19/6 and d(U5, 08) = (2+3+4+5+2+3)/6 = 19/6. These tie; (U5, 07) comes first.
- d(U6, 08) = (19 + 4)/7 = 23/7.
- d(U7, 09) = (6+7+8+9+7+8+8+8)/8 = 61/8.

Steps with more than one tied pair: **4** (steps 1, 2, 5, 6). At steps 5 and 6 the alternative choice changes the
topology, so the tie rule matters.
**Newick:** `((((((01,02),(05,06)),(03,04)),07),08),09);`
**In-order traversal** (left, node, right):
01, U1, 02, U4, 05, U3, 06, U5, 03, U2, 04, U6, 07, U7, 08, U8, 09.
Heights are non-decreasing up the tree, as UPGMA guarantees.

**Ultrametric?** The input distances are not ultrametric. Witness (01, 02, 03): d = 1, 2, 1. The two largest values
(2 and 1) differ.
**Cophenetic distances** (2 × the height of the lowest common ancestor) reproduce d_H exactly for **3 of 36** pairs
(01–02, 03–04, 05–06). Largest deviation per merge level:
- U4: 1/2 (1 vs 3/2);
- U5: 3/2 (01–04 = 3 vs 5/2 is 1/2; 02–03 = 1 and 05–04 = 4 are both 3/2);
- U6: 7/6 (01–07 = 2 vs 19/6);
- U7: **12/7 (04–08 = 5 vs 23/7)**;
- U8: 13/8 (01–09 = 6 vs 61/8).

The largest overall is **12/7, at 04–08**.
*(Corrected during Stage 6: my first hand pass checked only the d = 2 pairs at U7, missed d(04,08) = 5, and reported
13/8 at 01–09. The implementation disagreed, and re-deriving by hand confirmed 12/7.)*
The difference is the cost of forcing a clock-like (ultrametric) tree onto non-ultrametric data.

## 7. MST vs UPGMA: non-trivial bipartitions (both sides have ≥ 2 vertices)
- **MST** (cut each edge): 01–02 separates {02,03,04,06}, and 02–03 separates {03,04}. Every other edge cuts off a single
  vertex. Total: **2** splits.
- **UPGMA** (each cluster vs its complement): {01,02}, {03,04}, {05,06}, {01,02,05,06}, {01..06} (complement
  {07,08,09}), {01..07} (complement {08,09}). Total: **6** splits.

**Shared: 1** ({03,04}). Symmetric difference: 2 + 6 − 2 = **6**.

## 8. Counting link (optional, your Stage 0 decision)
For n = 9 labelled vertices: n^(n−2) = **4,782,969** labelled trees.
(2n − 3)!! = 15!! = **2,027,025** rooted binary topologies with 9 labelled leaves.
