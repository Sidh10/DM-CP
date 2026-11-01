# Experiment 7: Tree analysis: MST versus UPGMA

- run date: 2026-10-31
- command: `python -m experiments.exp7_tree_analysis`
- Python 3.13.11 on Windows

## Question
How do a minimum spanning tree and a UPGMA tree built from the same Hamming distances differ in structure, assumptions and what they can be used to claim?

## Input
Dataset A (9 haplotypes) and Dataset B (30 haplotypes); complete distance graphs and distance matrices with entries d_H.

## Mathematical model
MST (Kruskal) on K_N, rooted at an arbitrary vertex for traversal; rooted-tree statistics; pre-, post- and level-order. UPGMA: average linkage, heights = distance / 2, tie-break D8 (lexicographic); in-order traversal of the binary tree; ultrametric test; cophenetic distances; non-trivial bipartitions of both trees.

## Procedure
Root the Dataset A MST at 01 and at 04 and compare statistics and traversals with the hand derivation. Run UPGMA and compare every merge, the Newick string, the in-order traversal, the ultrametric test and the cophenetic fit. Count shared bipartitions. Repeat on Dataset B and report the results as observations.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | A: MST rooted at 01 | {height: 3, leaves: [04, 05, 06, 07, 08, 09], internal: [01, 02, 03], preorder: [01, 02, 03, 04, 06, 05, 07, 08, 09], postorder: [04, 03, 06, 02, 05, 07, 08, 09, 01], level: [01, 02, 05, 07, 08, 09, 03, 06, 04]} | docs/unit6_expected.md §5 |
| 2 | A: MST rooted at 04 | {height: 4, leaves: [05, 06, 07, 08, 09], internal: [01, 02, 03, 04], preorder: [04, 03, 02, 01, 05, 07, 08, 09, 06], postorder: [05, 07, 08, 09, 01, 06, 02, 03, 04]} | docs/unit6_expected.md §5 |
| 3 | A: UPGMA merges (node, left, right, distance, height, tied pairs) | [[U1, 01, 02, 1, 1/2, 6], [U2, 03, 04, 1, 1/2, 2], [U3, 05, 06, 1, 1/2, 1], [U4, U1, U3, 3/2, 3/4, 1], [U5, U4, U2, 5/2, 5/4, 2], [U6, U5, 07, 19/6, 19/12, 2], [U7, U6, 08, 23/7, 23/14, 1], [U8, U7, 09, 61/8, 61/16, 1]] | docs/unit6_expected.md §6 |
| 4 | A: UPGMA Newick | ((((((01,02),(05,06)),(03,04)),07),08),09); | docs/unit6_expected.md §6 |
| 5 | A: merge steps with ties | 4 | docs/unit6_expected.md §6 |
| 6 | A: reverse tie rule changes the topology | True | docs/unit6_expected.md §6 |
| 7 | A: UPGMA in-order traversal | [01, U1, 02, U4, 05, U3, 06, U5, 03, U2, 04, U6, 07, U7, 08, U8, 09] | docs/unit6_expected.md §6 |
| 8 | A: input distances ultrametric? (witness) | (False, (01, 02, 03)) | docs/unit6_expected.md §6 |
| 9 | A: cophenetic fit to d_H | {exact: 3, pairs: 36, max_deviation: 12/7, worst_pair: (04, 08)} | docs/unit6_expected.md §6 |
| 10 | A: MST vs UPGMA non-trivial bipartitions | {mst: 2, upgma: 6, shared: 1, symmetric_difference: 6} | docs/unit6_expected.md §7 |
| 11 | A: counting link for n = 9 | {cayley n^(n-2): 4782969, (2n-3)!!: 2027025} | docs/unit6_expected.md §8 |
| 12 | B: UPGMA heights non-decreasing | True | UPGMA property — *consistency check* |
| 13 | B: UPGMA cophenetic matrix is ultrametric | True | UPGMA property — *consistency check* |
| 14 | B: input distances ultrametric? (witness) | *(none: reported as observation)* |  |
| 15 | B: cophenetic fit to d_H | *(none: reported as observation)* |  |
| 16 | B: UPGMA merge steps with ties; reverse rule changes topology | *(none: reported as observation)* |  |
| 17 | B: MST vs UPGMA non-trivial bipartitions | *(none: reported as observation)* |  |
| 18 | B: MST rooted at the medoid (height, #leaves, #internal) | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | A: MST rooted at 01 | {height: 3, leaves: [04, 05, 06, 07, 08, 09], internal: [01, 02, 03], preorder: [01, 02, 03, 04, 06, 05, 07, 08, 09], postorder: [04, 03, 06, 02, 05, 07, 08, 09, 01], level: [01, 02, 05, 07, 08, 09, 03, 06, 04]} |
| 2 | A: MST rooted at 04 | {height: 4, leaves: [05, 06, 07, 08, 09], internal: [01, 02, 03, 04], preorder: [04, 03, 02, 01, 05, 07, 08, 09, 06], postorder: [05, 07, 08, 09, 01, 06, 02, 03, 04]} |
| 3 | A: UPGMA merges (node, left, right, distance, height, tied pairs) | [[U1, 01, 02, 1, 1/2, 6], [U2, 03, 04, 1, 1/2, 2], [U3, 05, 06, 1, 1/2, 1], [U4, U1, U3, 3/2, 3/4, 1], [U5, U4, U2, 5/2, 5/4, 2], [U6, U5, 07, 19/6, 19/12, 2], [U7, U6, 08, 23/7, 23/14, 1], [U8, U7, 09, 61/8, 61/16, 1]] |
| 4 | A: UPGMA Newick | ((((((01,02),(05,06)),(03,04)),07),08),09); |
| 5 | A: merge steps with ties | 4 |
| 6 | A: reverse tie rule changes the topology | True |
| 7 | A: UPGMA in-order traversal | [01, U1, 02, U4, 05, U3, 06, U5, 03, U2, 04, U6, 07, U7, 08, U8, 09] |
| 8 | A: input distances ultrametric? (witness) | (False, (01, 02, 03)) |
| 9 | A: cophenetic fit to d_H | {exact: 3, pairs: 36, max_deviation: 12/7, worst_pair: (04, 08)} |
| 10 | A: MST vs UPGMA non-trivial bipartitions | {mst: 2, upgma: 6, shared: 1, symmetric_difference: 6} |
| 11 | A: counting link for n = 9 | {cayley n^(n-2): 4782969, (2n-3)!!: 2027025} |
| 12 | B: UPGMA heights non-decreasing | True |
| 13 | B: UPGMA cophenetic matrix is ultrametric | True |
| 14 | B: input distances ultrametric? (witness) | (False, (CY249987.1, CY257981.1, CY260146.1)) |
| 15 | B: cophenetic fit to d_H | {exact: 9, pairs: 435, max_deviation: 73/14, worst_pair: (MG974757.1, MH586436.1)} |
| 16 | B: UPGMA merge steps with ties; reverse rule changes topology | {tie_steps: 11, reverse_differs: True} |
| 17 | B: MST vs UPGMA non-trivial bipartitions | {mst: 10, upgma: 27, shared: 5, symmetric_difference: 27} |
| 18 | B: MST rooted at the medoid (height, #leaves, #internal) | {root: CY249987.1, height: 4, leaves: 19, internal: 11} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | A: MST rooted at 01 | yes |
| 2 | A: MST rooted at 04 | yes |
| 3 | A: UPGMA merges (node, left, right, distance, height, tied pairs) | yes |
| 4 | A: UPGMA Newick | yes |
| 5 | A: merge steps with ties | yes |
| 6 | A: reverse tie rule changes the topology | yes |
| 7 | A: UPGMA in-order traversal | yes |
| 8 | A: input distances ultrametric? (witness) | yes |
| 9 | A: cophenetic fit to d_H | yes |
| 10 | A: MST vs UPGMA non-trivial bipartitions | yes |
| 11 | A: counting link for n = 9 | yes |
| 12 | B: UPGMA heights non-decreasing | yes (consistency) |
| 13 | B: UPGMA cophenetic matrix is ultrametric | yes (consistency) |
| 14 | B: input distances ultrametric? (witness) | n/a |
| 15 | B: cophenetic fit to d_H | n/a |
| 16 | B: UPGMA merge steps with ties; reverse rule changes topology | n/a |
| 17 | B: MST vs UPGMA non-trivial bipartitions | n/a |
| 18 | B: MST rooted at the medoid (height, #leaves, #internal) | n/a |

**11/11 independent expected values matched; 2/2 consistency checks passed; 5 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

### Dataset B UPGMA tree (Newick; branch length = difference of UPGMA heights, same units as d_H; the path between two leaves sums to their cophenetic distance; observation)

`(((((((((((((((((((CY249987.1:0.5,CY257981.1:0.5):0.25,MG914163.1:0.75):0.08333,MH126465.1:0.8333):0.04167,MH294937.1:0.875):0.025,MK401406.1:0.9):0.5167,((MH082676.1:0.5,MH637643.1:0.5):0.25,MN538542.1:0.75):0.6667):0.1944,MH233375.1:1.611):0.03889,MN538725.1:1.65):0.2591,CY264957.1:1.909):0.2576,MH898008.1:2.167):0.141,MH127265.1:2.308):0.5495,MH540789.1:2.857):0.009524,(MH701052.1:2,MK503016.1:2):0.8667):0.5157,MH792626.1:3.382):0.03431,MN538624.1:3.417):0.5307,(CY264115.1:1,MN538442.1:1):2.947):1.1,MK241188.1:5.048):10.94,(((CY260146.1:3,CY261017.1:3):3,MN568911.1:6):3.667,(MG974757.1:1.5,(MH778984.1:1,MK382008.1:1):0.5):8.167):6.318):9.908,(MH586436.1:1.5,MH885199.1:1.5):24.39);`

## Interpretation
The two trees answer different questions. The MST connects the observed sequences directly with the least total substitution count; its internal nodes are observed haplotypes, and its root is arbitrary (on Dataset A moving the root from 01 to 04 changes the height from 3 to 4 and the number of leaves from 6 to 5 without changing a single edge). UPGMA builds hypothetical ancestors (U1, U2, …) and places all sequences at the leaves under a molecular-clock assumption. Dataset A's distances are not ultrametric, and UPGMA reproduces only 3 of 36 distances exactly (an ultrametric input would be reproduced in full), and the two trees share 1 of their 2 and 6 non-trivial bipartitions. Dataset B (observation): the distances are not ultrametric; the UPGMA tree reproduces 9 of 435 distances exactly (largest deviation 73/14 at ('MG974757.1', 'MH586436.1')); 11 merge steps had ties, and the reverse tie rule changes the topology; the MST and UPGMA share 5 of their 10 and 27 non-trivial bipartitions.

## Limitations
Neither tree is evidence of evolutionary history. The MST is a graph-theoretic object: it minimises total edge weight and nothing else, and its rooting is arbitrary. UPGMA assumes equal evolutionary rates (a molecular clock); where distances are not ultrametric its branch heights are averages, not times. Raw Hamming distance has no multiple-substitution correction. Ties are resolved by a fixed rule (D8), and a different rule can give a different tree. Established phylogenetic tools are out of scope.
