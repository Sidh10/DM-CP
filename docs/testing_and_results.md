# Deliverable E — Testing and results

## The validation method

The guiding rule is that **the code is checked against mathematics, not against itself**. Three mechanisms enforce it.

**1. Expected values are derived by hand, before implementation.** For each unit, a document in `docs/` works the
answers out on paper — `synthetic_expected.md` for Dataset A and `unit2_expected.md` … `unit6_expected.md` for the
rest. Each experiment imports those values from an `experiments/expected_*.py` file that only transcribes them, with
the source section named per entry. A results file prints **Expected** and **Actual** as separate tables and then a
comparison, so a disagreement is visible rather than absorbed.

**2. Every check is labelled by what it can actually detect.**
- **Independent** — the actual value is compared with one fixed outside `genomath`: a hand derivation, the numeric
  prediction of a hand-proved theorem, or the quality-control script, which shares no code with `genomath`.
- **Consistency** — the check can only fail if `genomath` disagrees with itself: a round trip, two of its own routes
  agreeing, or a bound that holds for any data. Useful, but weaker, so it is counted separately.
- **Observation** — a measurement with no independent expected value, reported and never scored.

**3. Independent library oracles in the test suite.** `genomath` imports no third-party library, so the tests can use
them as referees: NetworkX for graphs and trees, SymPy for group theory and number theory, Biopython for parsing and
reverse complement, and Python's own `math`, `pow` and `int(…, 4)`.

| Oracle | Used for |
|---|---|
| NetworkX | components, degrees, bipartiteness, triangles, cycle basis, BFS distances, Dijkstra, all shortest paths, κ and λ, MST weight, MST count via `SpanningTreeIterator` |
| SymPy | `PermutationGroup` order/abelian/cyclic/orbits, `Permutation` cycle form and parity, `factorint`, `totient`, `n_order`, `crt`, exact determinants |
| Biopython | FASTA parsing, reverse complement, translation in the QC script |
| `math`, `pow`, `int(s, 4)` | `comb`, `perm`, `factorial`, `gcd`, modular inverse, base-4 parsing |
| Brute force | all 512 relations on 3 points, all 24 permutations of ACGT, enumeration of Σ^n for n ≤ 8, every linear congruence for m < 25, the bucket criterion for every m < 70 |

### Mutation checking
Passing tests prove nothing about tests that cannot fail. `scripts/mutation_check.py` injects **67** deliberate bugs
one at a time — a wrong recurrence factor, Warshall skipping its last pivot, Dijkstra settling the farthest vertex,
UPGMA's tie-break reversed, an inverted layout force — and requires the suite to catch each. Files are always restored.
Current status: **67 of 67 caught**.

Three candidate mutants were dropped as **equivalent**, with the reason recorded in the script: swapping the
composition order in an abelian group's Cayley table, re-relaxing settled vertices in Dijkstra, and removing CRT's
explicit coprimality check (the inverse step fails anyway). An equivalent mutant cannot change behaviour, so requiring
a test to catch it would be requiring a test to detect nothing.

## Test suite

179 tests, run with warnings promoted to errors (`pytest -W error`).

| File | Tests | Covers |
|---|---|---|
| `test_sequences.py` | 24 | parsing, validation, haplotypes, medoid, Hamming distance, trace capping |
| `test_set_logic.py` | 23 | mutation sets, set operations, the distance identity, relations, Warshall, predicates |
| `test_counting.py` | 21 | binomials, the recurrence, generation and its guard, inclusion–exclusion, pigeonhole, k! orderings |
| `test_algebra.py` | 24 | Cayley tables, axioms, orbits, Burnside, the substitution monoid, the base-level group |
| `test_number_theory.py` | 16 | encoding, Euclid, inverses, congruences, CRT, orders, the bucket criterion |
| `test_graphs.py` | 12 | construction, matrices, traversals, components, connectivity, Dijkstra |
| `test_trees.py` | 14 | Prim, Kruskal, MST counting, rooted trees, traversals, UPGMA |
| `test_app.py` | 16 | the pure presentation layer and the import boundaries |
| `test_experiments.py` | 9 | every experiment still matches its expected values; every figure quoted in the documentation still matches the code |
| `test_skeleton.py` | 20 | modules import; no UI or oracle imports leak into `genomath` |

### Required edge cases (specification §10)
| Case | Where |
|---|---|
| invalid input | non-ACGT dropped with a reason; malformed FASTA; duplicate ids |
| empty and single-sequence | `validate([])` raises; a one-sequence dataset is its own reference |
| repeated sequences | syn03 = syn10 collapse; two triples collapse in Dataset B |
| different sequence lengths | rejected with `SequenceError` carrying the trace |
| ambiguous symbols | the two Dataset B records with R and Y, dropped and logged |
| counting boundaries | k = 0, k = n, k > n, n = 0; the 10⁵ generation guard |
| number-theory boundaries | gcd(0,0), non-invertible inverses, non-coprime CRT moduli, unsolvable congruences |
| disconnected graphs | syn09 isolated for k ≤ 5; Dataset B G_=1 has 22 components; unreachable Dijkstra targets |
| multiple optimal solutions | 8 MSTs on Dataset A, 4 on Dataset B, 3 equally short paths syn04 → syn08, 11 tied UPGMA merges |

## Experiment results

All eight experiments pass. Totals: **147 independent expected values matched, 28 consistency checks passed, 35
observations**.

| Experiment | Independent | Consistency | Observations |
|---|---|---|---|
| 1 Mutation-set analysis | 20/20 | 2/2 | 3 |
| 2 Counting versus enumeration | 16/16 | 1/1 | 2 |
| 3 Algebraic verification | 29/29 | 1/1 | 1 |
| 4 Encoding and modular indexing | 23/23 | 5/5 | 4 |
| 5 Mutation-network structure | 14/14 | 2/2 | 2 |
| 6 Shortest paths and spanning trees | 14/14 | 4/4 | 4 |
| 7 Tree analysis: MST versus UPGMA | 11/11 | 2/2 | 5 |
| 8 Real-data evaluation | 20/20 | 11/11 | 14 |

### Dataset A — hand-verified throughout
Every value below was derived on paper first and then reproduced by the implementation.

| Quantity | Value |
|---|---|
| haplotypes from 10 records | 9 (syn03 = syn10) |
| d_H matrix | all 81 entries match |
| distance identity (Theorem 1) | holds on all 36 pairs |
| R_1 | reflexive, symmetric, **not transitive**; witness (syn01, syn02, syn03) |
| components at k = 0, 1, 2, 5, 6 | 9, 4, 2, 2, 1 |
| G_2 | 14 edges, degrees 6 5 4 2 4 4 1 2 0, Σ deg = 28 = 2\|E\|, cyclomatic 7, not bipartite |
| Dijkstra on G_2 | syn04→syn08 cost 5 = d_H with 3 optimal paths; syn07→syn03 cost 4 > d_H = 3 |
| MST | weight 15 by both Prim and Kruskal; **8 distinct MSTs**; largest edge 6 = k\* |
| rooted MST | at syn01: height 3, 6 leaves; at syn04: height 4, 5 leaves — same 8 edges |
| UPGMA | 4 tied merge steps; heights 1/2, 3/4, 5/4, 19/12, 23/14, 61/16; reproduces 3 of 36 distances |
| orbits under {I, RC} | 8 classes from 9 haplotypes; syn03/syn08 merge; syn01 is RC-fixed |
| encodings | all 9 match, with enc(comp S) = 4¹² − 1 − enc(S) |
| buckets | m = 16 collides on 4 of the 6 distance-1 pairs; every odd modulus on none |

### Dataset B — observations, with no biological claim
| Quantity | Value |
|---|---|
| records read / kept / haplotypes | 36 / 34 / 30 |
| dropped | MH246583.1 (R), MK181701.1 (R, Y) |
| length | 1701 nt, all records |
| pairwise distances | 435 pairs; min 1, median 9, max 57 (3.4% of length) |
| pairs at distance 1 | 8 |
| reference (medoid) | CY249987.1 |
| sparsity | best-connected haplotype has 6 observed neighbours out of 5,103 possible |
| G_=1 | 8 edges, 22 components, **a forest**: no cycles, max degree 6 |
| connects at | k\* = 47 |
| MST | weight 171; 4 distinct MSTs; largest edge 47 = k\*; Prim and Kruskal return **different** trees |
| UPGMA | not ultrametric; reproduces 9 of 435 distances; 11 tied merge steps |
| MST vs UPGMA | 10 and 27 non-trivial bipartitions, 5 shared |
| buckets, m = 16 and 64 | **all 30 haplotypes in one bucket** (every sequence ends TGA) |
| buckets, m = 61 | 24 buckets occupied, max load 2 |
| distance-1 pairs colliding | 8/8 under m = 16 and 64; 0/8 under 7, 13, 17, 31, 61 |
| canonical k-mers | k = 4: 248 distinct → 136 canonical (all of the whole space); k = 8: 2,563 → 2,488 |
| substitutions (direction-free) | vs reference: 250 transitions, 91 transversions · along MST edges: 131 and 40 |

## Complexity summary

| Algorithm | Time | Note |
|---|---|---|
| mutation comparison, d_H | Θ(n) per pair, Θ(N²·n) for all | computed once and reused |
| recursive mutation generation | Θ(n · M(n,k)) | exponential; guarded at 10⁵ |
| sequence encoding / decoding | Θ(n) big-integer steps | |
| graph construction | Θ(N²) after the distance matrix | |
| BFS / DFS / components | O(V + E) | |
| Dijkstra | O(V² + E) | array version, no library heap |
| Prim | Θ(V²) | array version |
| Kruskal | Θ(E log E) + near-linear union–find | hand-written merge sort |
| MST counting | O(E log E + Σ V_level³) | Kirchhoff per weight level |
| tree traversals | Θ(V) | |
| UPGMA | Θ(N³) | exact `Fraction` arithmetic |
| κ, λ | exponential | refused above 16 vertices |
| Burnside by formula | O(\|G\|) big-integer powers | any n |
| orbit enumeration | Θ(4^n·\|G\|·n) | refused above n = 10 |

**On the exponential algorithms.** Three are genuinely exponential and all three are guarded rather than left to hang:
mutation generation (refused above 10⁵ strings, with the predicted count shown), exhaustive orbit enumeration
(refused above n = 10, where 4^10 ≈ 10⁶), and exhaustive κ/λ search (refused above 16 vertices). The guard computes the
size first from a closed form, so refusal is instant.

## Errors found during development

The method's value is best shown by what it caught. All were found before the work was reported as finished.

| # | Error | How it surfaced |
|---|---|---|
| 1 | DFS classified back edges by comparing vertex indices; an ancestor can have a larger index | writing the test "back edges = cyclomatic number" on random graphs |
| 2 | My NetworkX oracle was called with the wrong option, disagreeing with the correct transitive closure on 252 of 512 relations | exhaustive check over all relations on 3 points |
| 3 | The dendrogram was built from the last merge down, so a parent needed its children's positions first | a layout test written before looking at the page |
| 4 | `fixed_point_free_involutions` counted each matching twice, returning 6 | the hand derivation said 3!! = 3 |
| 5 | My hand value for the largest UPGMA deviation was wrong (13/8 at 01–09; it is 12/7 at 04–08) | the implementation disagreed; re-deriving by hand confirmed the code |
| 6 | My first statement of the bucket criterion was wrong for δ = ±2 | writing out the proof before coding |
| 7 | Experiment 1 treated "alone in its orientation class" as RC-fixed | the expected/actual comparison, at 21/22 |
| 8 | matplotlib figures were theme-dependent: white background, black text | checking rcParams rather than trusting a screenshot |
| 9 | Arrow serialisation failed on every dashboard table, hidden behind a silent fallback | reading the server log, not the page |
| 10 | Tests too weak to catch three injected bugs (Warshall's last pivot, transitivity direction, the λ edge lookup) | the mutation check |

Items 4, 5 and 6 are the ones that justify the method: two were caught **because** the expected values existed on
paper first, and one was a mistake in the hand derivation itself that the code exposed.
