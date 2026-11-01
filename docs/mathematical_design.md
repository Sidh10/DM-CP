# Deliverable B — Mathematical design document

Notation used throughout: Σ = {A, C, G, T}; [n] = {1, …, n}; a sequence is S ∈ Σ^n, equivalently a function
S: [n] → Σ with S_i the base at position i; R is the chosen reference; V is the set of haplotypes and N = |V|.
Full hand derivations are in `docs/synthetic_expected.md` and `docs/unit2_expected.md` … `docs/unit6_expected.md`;
this document states the mathematics and points to them for the long working.

---

## Unit I — Set theory and logic

### Objects and definitions
- **Alphabet and sequences.** Σ^n is the set of length-n strings. |Σ^n| = 4^n.
- **Mutation-position set.** M(R,S) = { i ∈ [n] : R_i ≠ S_i } ⊆ [n].
- **Allele-level mutation set.** M\*(R,S) = { (i, S_i) : R_i ≠ S_i } ⊆ [n] × Σ.
- **Hamming distance.** d_H(S,T) = |{ i : S_i ≠ T_i }|. Hence d_H(R,S) = |M(R,S)|.
- **Threshold relation.** S R_k T ⇔ d_H(S,T) ≤ k.
- **Predicates.** Boolean-valued functions on Σ^n, closed under ∧, ∨, ¬.

(The project specification writes the mutation set as "{i ≠ S_i}"; the intended definition is the one above.)

### Theorem 1 (distance from mutation sets)
For any reference R and any S, T ∈ Σ^n,

> d_H(S,T) = |M(R,S) Δ M(R,T)| + |{ i ∈ M(R,S) ∩ M(R,T) : S_i ≠ T_i }|.

*Proof.* Fix a position i and compare the contribution of each side.
- i in neither set: S_i = R_i = T_i, so i contributes 0 to d_H, 0 to the symmetric difference and 0 to the second term.
- i in exactly one set, say M(R,S): then S_i ≠ R_i = T_i, so i contributes 1 to d_H and 1 to the symmetric difference.
- i in both: S_i ≠ R_i and T_i ≠ R_i, so i contributes 1 to d_H exactly when S_i ≠ T_i, which is exactly when it is
  counted by the second term, and 0 to the symmetric difference.
Summing over positions gives the identity. ∎

The second term vanishes over a binary alphabet; it exists here because a position has three alternatives. This is why
M alone is insufficient and M\* is also computed.

### Theorem 2 (R_k is not an equivalence relation)
R_k is reflexive (d_H(S,S) = 0 ≤ k) and symmetric (d_H is symmetric), but **not transitive** for k ≥ 1. A single
counterexample suffices, and Dataset A provides one: d(syn01, syn02) = 1 and d(syn02, syn03) = 1, but
d(syn01, syn03) = 2 > 1. By contrast R_0 restricted to haplotypes is equality, which *is* an equivalence relation, and
its classes are the haplotypes.

The transitive closure of R_k (computed by Warshall's algorithm) **is** an equivalence relation, and its classes are
exactly the connected components of the threshold graph G_k of Unit V. This is the bridge between Units I and V.

### Worked example (Dataset A, R = syn01)
M(R, syn03) = {2, 5} and M(R, syn07) = {5, 7}, so M ∩ M' = {5}: they share a *position*.
But M\*(R, syn03) = {(2,C), (5,G)} and M\*(R, syn07) = {(5,C), (7,A)} have empty intersection: they share no *mutation*.
Theorem 1: d_H(syn03, syn07) = |{2,7}| + |{5}| = 2 + 1 = 3, which matches the direct count.

### Algorithms and complexity
| Operation | Algorithm | Time | Space |
|---|---|---|---|
| M, M\*, d_H | one linear scan of the pair | Θ(n) | O(\|M\|) |
| set operations | explicit membership loop, sorted for a stable trace | O(m log m), m = \|A\|+\|B\| | O(m) |
| relation properties | successor lists, early exit with a witness | O(N³) worst case | O(N²) |
| transitive closure | Warshall, three nested loops | Θ(N³) | O(N²) |
| equivalence classes | closure, then grouping | O(N³) | O(N²) |

### Contribution
Unit I supplies the vocabulary everything else uses: the haplotype quotient that makes the graph simple, the
distance that weights its edges, and the failure of transitivity that explains why components merge gradually as the
threshold grows rather than all at once.

---

## Unit II — Counting and recursion

### Definitions and results
- **Sequence space.** N(n) = 4^n, by the product rule over n independent positions.
- **Mutant count.** M(n,k) = C(n,k)·3^k: choose which k positions change, then one of 3 alternative bases at each.
  Substitutions only — insertions and deletions change the length and leave Σ^n.
- **Recurrence.** M(n,k) = M(n−1,k) + 3·M(n−1,k−1), with M(0,0) = 1 and M(n,k) = 0 outside 0 ≤ k ≤ n.
  *Derivation:* split on position n. Either it is unchanged, and the first n−1 positions carry all k changes; or it is
  changed, in 3 ways, and the first n−1 positions carry k−1. The cases are disjoint and exhaustive.
- **Hamming ball.** |B(n,r)| = Σ_{k ≤ r} M(n,k).
- **Inclusion–exclusion.** Strings using all four bases: U(n) = Σ_j (−1)^j C(4,j) (4−j)^n = 4^n − 4·3^n + 6·2^n − 4.
- **Pigeonhole.** A string of length L has L − k + 1 windows of length k and there are 4^k possible k-mers, so a
  repeated k-mer is forced whenever L − k + 1 > 4^k.
- **Permutations.** Between sequences at distance k there are exactly **k!** shortest substitution orderings: each
  differing position must change exactly once, directly to its target base.

### Theorem 3 (the closed form satisfies the recurrence)
By induction on n. Base: f(0,0) = C(0,0)·3^0 = 1 and f(0,k) = 0 for k ≠ 0. Step: assuming the claim for n−1,

  f(n−1,k) + 3·f(n−1,k−1) = C(n−1,k)·3^k + 3·C(n−1,k−1)·3^(k−1) = 3^k·[C(n−1,k) + C(n−1,k−1)] = 3^k·C(n,k),

using Pascal's rule. ∎ A free consistency check follows from the binomial theorem: Σ_k M(n,k) = (3+1)^n = 4^n.

### Worked example
n = 12: 4^12 = 16,777,216; M(12,1) = 12·3 = 36; M(12,2) = 66·9 = 594; |B(12,2)| = 1 + 36 + 594 = 631.
Recurrence check: M(12,2) = M(11,2) + 3·M(11,1) = 495 + 99 = 594.
For Dataset B, n = 1701: 4^1701 has ⌊1701·log₁₀4⌋ + 1 = 1025 digits, and each sequence has 3n = 5,103
single-substitution neighbours.

### Algorithms and complexity
| Operation | Algorithm | Time |
|---|---|---|
| C(n,k) | multiplicative, exact integer division at each step | O(min(k, n−k)) |
| M(n,k) table | bottom-up dynamic programme over the recurrence | Θ(n²) |
| M(n,k) recursive | top-down with an explicit memo | O(n·k) subproblems |
| generation | recursion over which k positions change, 3 bases each | Θ(n · M(n,k)) |
| repeated k-mer | sliding window into a hash set | O(L·k) |
| k! orderings | counted always; listed only for k ≤ 8 | Θ(k!·k) when listed |

Generation is exponential, so it is guarded: the count is computed from the closed form first and the run is refused
above 10⁵ strings (decision D5). M(12,6) = 673,596 and M(1701,2) = 13,012,650 are both refused; M(1701,1) = 5,103 runs.

### Contribution
Unit II sizes the problem and sets the enumeration limits. Its most load-bearing result in this project is the
sparsity comparison: 6 observed neighbours against 5,103 possible is what licenses calling the single-substitution
graph a thin sample rather than a picture of the mutational landscape.

---

## Unit III — Algebraic structures

### Two different groups acting on two different sets
**(a) On strings.** I (identity), R (reverse), C (complement, A↔T and C↔G letterwise), and RC = R∘C act on Σ^n.
**(b) On bases (decision D9).** The three fixed-point-free involutions of Σ — (A T)(C G), (A G)(C T), (A C)(G T) —
together with the identity act on Σ itself.

Both are Klein four-groups, abstractly isomorphic to Z2 × Z2, but they act on different sets and must not be conflated.
They are linked: the string map C is the letterwise extension of the base permutation (A T)(C G).

### Theorem 4 (both are Klein four-groups)
R and C are commuting involutions: R²= I because reversing twice restores the order; C² = I because complementing a
base twice restores it; and RC = CR because C changes letters without moving them while R moves letters without
changing them. Hence the four maps are closed under composition, each is its own inverse, composition of functions is
associative, and no non-identity element has order 4. So the group is Klein, not cyclic. ∎
For the base-level group, σ∘τ is computed directly: A → G → C, C → T → A, G → A → T, T → C → G, which is (A C)(G T).

**Subtlety (faithfulness).** On Σ⁰ all four string maps coincide; on Σ¹, R = I and RC = C; they are distinct only for
n ≥ 2. Burnside's lemma does not require a faithful action, so the orbit counts remain valid at n = 0 and n = 1, but
"the four maps are distinct" is false there. Function equality is therefore decided exhaustively over all 5,461 words
of length ≤ 6, never on a sample.

### Theorem 5 (fixed points and Burnside)
|Fix(I)| = 4^n. |Fix(R)| = 4^⌈n/2⌉ (the first half is free, the rest is forced). |Fix(C)| = 0 for n ≥ 1, since no base
is its own complement. |Fix(RC)| = 4^(n/2) for even n and 0 for odd n, because an odd middle base would have to equal
its own complement. Burnside's lemma, #orbits = (1/|G|)·Σ_g |Fix(g)|, then gives

> {I, RC}: (4^n + |Fix(RC)|)/2 — ten orbits at n = 2.
> {I, R, C, RC}: (4^n + 4^⌈n/2⌉ + 0 + |Fix(RC)|)/4 — six orbits at n = 2.

### Theorem 6 (mutations do not form a group)
Let σ_{i,b} be the operator "set position i to base b". Compositions of such operators are described by patterns
p ∈ {·, A, C, G, T}^n, where a later write wins; so they are closed, contain the identity and are associative — a
**monoid** with 5^n elements, every one of them idempotent. But any f_p with p_i = b ≠ · is not injective: two strings
differing only at position i have the same image, so it has no inverse. The group of units is {I}. ∎
This gives the specification's warning a proof: the {I, RC} group is legitimate because RC is a *bijection* describing
strand orientation, not a mutation.

### Theorem 7 (substitution classification, decision D9)
Each non-identity element of the base-level group moves all four bases, so each realises exactly 4 ordered pairs
(x, y) with y ≠ x; 3 × 4 = 12 covers every ordered pair, so every substitution comes from exactly one element. Because
each element is an involution, the same element also takes y to x, so the classification is well defined on **unordered**
pairs — which is the only honest option, since two observed sequences do not say which base came first.
τ = (A G)(C T) pairs purine with purine and pyrimidine with pyrimidine: these are the **transitions**. The other two
elements give the **transversions**. Of the 6 unordered base pairs, 2 are transitions and 4 are transversions, so
equal likelihood would give the ratio 1 : 2.

### Worked example
RC(syn03) = syn08 in Dataset A, so those two haplotypes form one orientation class and the dataset's 9 haplotypes fall
into 8 classes. syn01 is RC-fixed (an RC-palindrome). In Dataset B, no sequence's reverse complement is present, as
predicted: every sequence begins ATG, so its reverse complement would begin with the complement of a stop codon.

### Complexity
Cayley tables cost O(|G|³ · W · n) with W = 5,461 the number of test words; axiom checking is O(|G|³) table lookups;
fixed-point and Burnside counts are O(|G|) big-integer powers for any n; direct enumeration of orbits is Θ(4^n · |G| · n)
and is refused above n = 10.

### Contribution
Unit III supplies the orientation quotient used for canonical k-mers in Unit IV, the proof that bounds what "mutation
algebra" can claim, and the substitution classification applied to the real dataset.

---

## Unit IV — Number theory

### Definitions
Encoding A = 0, C = 1, G = 2, T = 3 and enc(S) = Σ_{i=1..n} s_i · 4^(n−i).

### Theorem 8 (encoding is a bijection)
enc maps Σ^n into {0, …, 4^n − 1}, and every integer in that range has exactly one n-digit base-4 representation
(leading zeros allowed), so enc is injective and surjective. Decoding applies the **division algorithm** n times:
x = 4q + r with 0 ≤ r < 4 yields the last digit r, and the process repeats on q. ∎

### Theorem 9 (complement identity)
Complement maps digit d ↦ 3 − d, so enc(comp S) = Σ (3 − s_i)4^(n−i) = 3·(4^n − 1)/3 − enc(S) = **4^n − 1 − enc(S)**. ∎

### Theorem 10 (exact single-substitution bucket criterion)
A substitution at position i changes the code by δ · 4^t with δ ∈ {±1, ±2, ±3} and t = n − i. The bucket x mod m is
unchanged iff m | δ·4^t. Write m = 2^a · q with q odd.
- Since gcd(q, 4^t) = 1, q | δ·4^t iff q | δ, which with |δ| ≤ 3 requires q ∈ {1, 3}.
- The power of 2 in δ·4^t is 2t + v, where v = 1 for δ = ±2 and v = 0 otherwise; we need 2t + v ≥ a.

Hence: **if the odd part of m is at least 5, no single substitution ever preserves the bucket, at any length.**
Otherwise one exists exactly when the sequence has a position with t ≥ t_min, where t_min = ⌈a/2⌉ for q = 3 (δ must be
±3) and ⌊a/2⌋ for q = 1 (δ = ±2 contributes its own factor of 2). ∎

This is sharper than the sufficient condition gcd(m,4) = 1: even moduli such as m = 10 also separate all single
substitutions. The criterion is verified by brute force for every m < 70 and n ≤ 5.

**Corollary (the bad hash).** If m | 4^j then every positional weight 4^t with t ≥ j is ≡ 0, so x mod m depends only on
the last j bases. m = 16 sees the last 2 bases, m = 64 the last 3.

### Theorem 11 (two substitutions can cancel)
For odd m the weights 4^t are periodic mod m with period ord_m(4), which divides φ(m) by Euler's theorem and divides
p − 1 for prime p by Fermat's little theorem. An opposite pair of changes (+δ at t₁, −δ at t₂) cancels exactly when
ord_m(4) divides |t₁ − t₂|. Worked example: changing syn01 at positions 6 and 12 shifts the code by 1 − 4⁶ = −4095 = −7·585,
so the result is at distance 2 yet shares a bucket under m = 7 (order 3) and m = 13 (order 6), but not under 17, 31 or 61.

### Further results used
- **Canonical k-mer.** c(S) = min(enc(S), enc(RC(S))), one representative per {I, RC}-orbit. With A < C < G < T, numeric
  order of equal-length codes agrees with lexicographic order, so this is the lexicographic minimum.
- **Sliding window.** The next k-mer code is x' = (4x mod 4^k) + d_new, which is O(1) per window.
- **Pigeonhole on buckets.** N keys in m buckets force a bucket with at least ⌈N/m⌉ keys.
- **CRT.** For pairwise coprime moduli with product P, an integer in [0, P) is uniquely recoverable from its residues.
  3-mer codes lie in [0, 64), so moduli with product ≥ 64 suffice; (3, 5) with product 15 recovers only 15 of 64.
- **Euclid / extended Euclid / modular inverses / linear congruences.** Standard; ax ≡ b (mod m) is solvable iff
  g = gcd(a,m) divides b, and then has exactly g solutions modulo m.

### Complexity
Euclid and extended Euclid take O(log min(a,b)) division steps; encoding, decoding and RC-on-codes are Θ(n)
big-integer steps; k-mer codes are Θ(L) for the whole sequence; CRT with t moduli is O(t²) coprimality checks plus
O(t log M); φ and factorisation are O(√n) by trial division.

### Contribution
Unit IV turns "represent DNA as integers" into provable statements about what such a representation preserves and
destroys. On real data it exhibits the destruction: under m = 16 or 64 all 30 haplotypes occupy a single bucket,
because they all end TGA.

---

## Unit V — Graph theory

### Definitions
- **Single-substitution graph** G_=1: vertices are haplotypes, edge {S,T} iff d_H(S,T) = 1.
- **Threshold graph** G_k: edge iff 1 ≤ d_H(S,T) ≤ k, weighted by d_H.
- Degrees, adjacency and incidence matrices, paths, cycles, components, vertex connectivity κ, edge connectivity λ.

### Results used
- **Handshake lemma.** Σ_v deg(v) = 2|E|.
- **Monotonicity.** G_j ⊆ G_k for j ≤ k, so edges never decrease and components never increase with k.
- **Components and closure.** The components of G_k are the classes of the transitive closure of R_k (Unit I).
- **Cyclomatic number.** |E| − |V| + c counts independent cycles; it is 0 exactly for a forest.
- **Whitney's inequality.** κ ≤ λ ≤ δ, where δ is the minimum degree.
- **Triangle inequality.** d_H is a metric, so any path's total weight is at least the direct distance between its
  endpoints: dist_{G_k}(u,v) ≥ d_H(u,v). Equality holds exactly when some shortest route never changes a position twice.
- **Bottleneck (minimax) property.** The smallest k for which G_k is connected equals the largest edge in any minimum
  spanning tree. Verified on both datasets: 6 on Dataset A, 47 on Dataset B.

Unlike the binary hypercube, the Hamming graph over a 4-letter alphabet contains triangles, since one position has
three alternatives; so G_=1 need not be bipartite.

### Worked example (Dataset A, G_2)
14 edges; degrees 6, 5, 4, 2, 4, 4, 1, 2, 0 summing to 28 = 2|E|; two components, with syn09 isolated; cyclomatic
number 14 − 9 + 2 = 7; not bipartite, because syn01–syn02–syn03 is a triangle at distances 1, 1, 2.
Dijkstra from syn04 reaches syn07 at cost 5 while d_H(syn04, syn07) = 4: the only route leaves syn07 through syn01,
which pays for position 5 twice. syn04 → syn08 has three equally short paths, so a test may compare costs but not edge
sets.

### Algorithms and complexity
| Algorithm | Time | Space |
|---|---|---|
| distance matrix | Θ(N²·n) | Θ(N²) |
| graph construction from the matrix | Θ(N²) | O(N+E) |
| BFS, DFS, components, cycle search, bipartiteness | O(V+E) | O(V) |
| triangles | O(E·Δ) | O(1) |
| Dijkstra, array version (no library heap) | O(V² + E) | O(V) |
| all shortest paths | O(V+E) for the predecessor DAG, plus output | — |
| κ, λ by exhaustive removal | exponential; refused above 16 vertices | O(V) |

### Contribution
Unit V is the project's main analytical component. The threshold sweep, rather than any single k, is what the
project reports, because a single threshold is an arbitrary choice.

---

## Unit VI — Trees

### Definitions and results
- **Spanning tree.** A connected acyclic subgraph on all N vertices; it has exactly N − 1 edges.
- **MST.** A spanning tree of minimum total weight. The minimum *weight* is unique; the *tree* need not be.
- **Kruskal.** Sort edges by weight and accept an edge iff it joins two different components (union–find).
- **Prim.** Grow one tree, repeatedly adding the cheapest edge leaving it.
- **Counting MSTs.** Process weight levels in increasing order; contract the components formed by lighter edges; at each
  level the number of ways is the number of spanning forests of the multigraph of crossing edges, computed by
  **Kirchhoff's matrix-tree theorem** (any cofactor of the Laplacian). Multiply across levels.
- **Cayley's formula.** n^(n−2) labelled trees on n vertices; (2n−3)!! rooted binary topologies on n labelled leaves.
- **UPGMA.** Repeatedly merge the two closest clusters; the distance between clusters is the mean over all member
  pairs, d(A∪B, W) = (|A|·d(A,W) + |B|·d(B,W)) / (|A|+|B|); a merge at distance d has height d/2.
- **Ultrametric.** A distance is ultrametric iff in every triple the two largest values are equal. UPGMA's output is
  always ultrametric and reproduces an ultrametric input exactly.
- **Cophenetic distance.** Twice the height of the lowest common ancestor; comparing it with d_H measures the cost of
  forcing a clock-like shape onto data that is not ultrametric.

### Worked example (Dataset A)
Kruskal examines 29 edges, accepts 8 and rejects 21, for total weight 5·1 + 2·2 + 6 = 15. Prim from syn01 gives the
same weight. The square syn01–syn02–syn06–syn05 of unit edges forces one rejection, and syn08 can be attached two
ways, so by the level-contraction argument there are 4 × 2 × 1 = **8 distinct MSTs**. Rooting the MST at syn01 gives
height 3 and 6 leaves; rooting the *same* tree at syn04 gives height 4 and 5 leaves — the clearest demonstration that
rooting an MST is a presentation choice. UPGMA on the same distances needs exact fractions (19/6, 23/7, 61/8), has 4
tied merge steps, and reproduces only 3 of 36 distances exactly, since the input is not ultrametric (witness: syn01,
syn02, syn03 at distances 1, 2, 1).

### Complexity
Kruskal Θ(E log E) for the sort plus near-linear union–find; Prim Θ(V²) in the array version; Kirchhoff O(V³) exact
rational operations; UPGMA Θ(N³) with exact fractions; traversals and tree statistics Θ(V); the ultrametric test Θ(N³).

### Contribution
Unit VI produces the project's two summary structures and, more importantly, keeps them apart: an MST minimises total
edge weight on observed sequences, while UPGMA builds hypothetical ancestors under a molecular-clock assumption.
Neither is evidence of evolutionary history.

---

## Syllabus → implementation mapping

Topic names for Units I–III are taken from the course slide titles; the folders for Units IV–VI were not available, so
those topic names follow the project specification §5. "Verified by" names the experiment and the pytest file.

| Unit | Syllabus topic | Implemented in | Verified by |
|---|---|---|---|
| I | Set theory | `set_logic.mutation_set`, `allele_mutation_set`, `set_operations` | Exp 1 · `test_set_logic.py` |
| I | Relations; composition of relations | `threshold_relation`, `identity_relation`, `orientation_relation` | Exp 1 · `test_set_logic.py` |
| I | Relation properties; equivalence | `relation_properties`, `equivalence_classes`, `transitive_closure` (Warshall) | Exp 1 · exhaustive over all 512 relations on 3 points |
| I | Partial order relation | `relation_properties` antisymmetry; the (M\*, ⊆) poset | `test_set_logic.py` |
| I | Functions | `function_properties` (injective, surjective, bijective with witnesses); `position_map` | `test_set_logic.py` |
| I | Mathematical logic | `Predicate` with ∧, ∨, ¬; `filter_sequences` truth tables | Exp 1 |
| II | Mathematical induction | Theorem 3 proof; `verify_recurrence` as the finite check | Exp 2 · `docs/unit2_expected.md` §2 |
| II | Basics of counting | `sequence_space`, `growth_table`, `mutation_count` | Exp 2 |
| II | Permutations and combinations | `binomial`, `permutations`, `factorial`, `mutation_orderings` (k!) | Exp 2 · vs `math.comb`, `math.perm` |
| II | Inclusion–exclusion | `inclusion_exclusion_all_bases` | Exp 2 · vs brute force n ≤ 7 |
| II | Pigeonhole | `forced_repeat`, `max_forced_k`, `find_repeated_kmer`; bucket bound in Unit IV | Exp 2, Exp 8 |
| II | Recursion / recurrence relations | `recurrence_table`, `mutation_count_recursive`, `generate_mutants` | Exp 2 |
| III | Algebraic structures; semigroup, monoid | `substitution_monoid` (Theorem 6) | Exp 3 · `test_algebra.py` |
| III | Groups; group axioms | `cayley_table`, `check_axioms`, `element_orders` | Exp 3 · vs SymPy `PermutationGroup` |
| III | Group examples | {I, RC} ≅ Z2; {I, R, C, RC} ≅ Z2×Z2; base-level Klein group | Exp 3, Exp 8 |
| III | Subgroups and cosets | `subgroups`, `left_cosets` (Lagrange index) | Exp 3 |
| III | Permutation groups | `permutation_cycles` (cycle form, order, parity) | Exp 3 · vs SymPy `Permutation` |
| III | Group actions, orbits, Burnside | `orbit`, `orbits_enumerated`, `burnside`, `canonical`, `dataset_orbits` | Exp 3 · formula vs enumeration |
| IV | Division algorithm; divisibility | `division`, `divides`, `decode` | Exp 4 |
| IV | GCD; Euclidean and extended Euclidean | `gcd`, `extended_gcd` | Exp 4 · vs `math.gcd`, Bézout identity |
| IV | Congruences; modular inverses | `mod_inverse`, `solve_linear_congruence` | Exp 4 · vs `pow(a,-1,m)` and brute force |
| IV | Chinese Remainder Theorem | `crt` | Exp 4 · vs SymPy, all 64 3-mer codes |
| IV | Euler's φ; Euler and Fermat theorems | `euler_phi`, `mod_pow`, `multiplicative_order` (labelled pedagogical) | Exp 4 · vs SymPy `totient`, `n_order` |
| IV | Modular arithmetic / indexing | `encode`, `bucket_report`, `last_bases_seen`, `substitution_collision_t_min` (Theorem 10) | Exp 4, Exp 8 |
| V | Graphs; adjacency and incidence matrices | `build_graph`, `adjacency_matrix`, `incidence_matrix` | Exp 5 |
| V | Degree; handshake lemma | `degrees` | Exp 5 |
| V | Paths, cycles, connectivity, components | `bfs`, `dfs`, `find_cycle`, `connected_components`, `cyclomatic_number`, `is_bipartite` | Exp 5 · vs NetworkX |
| V | Vertex and edge connectivity | `vertex_connectivity`, `edge_connectivity` (Whitney) | Exp 5 · vs NetworkX on random graphs |
| V | Graph operations | `induced_subgraph`, `complement`, `edge_set_operation` | `test_graphs.py` |
| V | Shortest paths; Dijkstra | `dijkstra`, `all_shortest_paths`, `path_versus_direct` | Exp 6 · vs NetworkX |
| VI | Trees; spanning trees | `validate_spanning_tree`, `kirchhoff_count`, `cayley_count` | Exp 6 |
| VI | Minimum spanning trees; Prim, Kruskal | `prim`, `kruskal`, `merge_sort`, `UnionFind` | Exp 6 · vs NetworkX; MST count vs `SpanningTreeIterator` |
| VI | Rooted trees; depth, height, leaves | `root_tree`, `tree_statistics` | Exp 7 |
| VI | Tree traversals (pre-, in-, post-order) | `preorder`, `inorder`, `postorder`, `level_order` | Exp 7 |
| VI | Distance-based trees | `phylogenetic_analysis.upgma`, `cophenetic`, `is_ultrametric`, `fit_to_input` | Exp 7 · vs an independent naive UPGMA |

Topics deliberately **not** claimed: rings (Unit III slides include one; no honest genomic use was found, so none is
asserted), and graph colouring, planarity and Euler/Hamiltonian circuits (no meaningful interpretation on a mutation
graph was identified).
