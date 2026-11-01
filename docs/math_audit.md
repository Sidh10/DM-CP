# Phase 1: Mathematical feasibility audit

Scope: the six units of MM0903 against `docs/master_prompt.md` §4–§5, with your Stage 0 decisions applied.
Each demo is tagged **core** (needed to answer the research question on real sequences) or **pedagogical**
(mathematically correct and on-syllabus, but not needed for genomic analysis). Items marked **[flag]** are claims I am not
fully sure of, or that need your decision.

## Verdict
The project is suitable. Units I, II, V and VI connect naturally: aligned equal-length sequences *are* functions
[n] → Σ, differences *are* sets, and the "one substitution apart" structure *is* a graph. Unit III is genuine but narrow:
reverse complement is a real strand symmetry, but it only does analytical work on k-mers of unknown orientation, not on
annotated coding sequences. Unit IV is strongest as representation mathematics: encoding, decoding, the complement identity,
canonical k-mers and a provable hashing result. CRT and Euler/Fermat are pedagogical. The proposed fixes are in "Weak links" below.

## Spec correction
§5 Unit I B prints M(R,S) = {i ≠ S_i}. The intended definition is **M(R,S) = { i ∈ [n] : R_i ≠ S_i }**, with [n] = {1,…,n}.
The spec's example {3,6,10} ∩ {6,10,15} needs n ≥ 15, so it is illustrative only for Dataset A (n = 12).

---

## Unit I — Set Theory and Logic
- **Problem.** State exactly where sequences differ from a reference and from each other, and filter a collection by such properties.
- **Objects.** Σ = {A,C,G,T}; a sequence S ∈ Σ^n viewed as a function S: [n] → Σ; the collection as a list, and its haplotypes
  as a set; position-level mutation sets M(R,S) ⊆ [n]; allele-level sets M*(R,S) = {(i, S_i) : S_i ≠ R_i} ⊆ [n] × Σ;
  relations R_k ⊆ V × V.
- **Equations.** d_H(S,T) = |{i : S_i ≠ T_i}|, so d_H(R,S) = |M(R,S)|.
  d_H(S,T) = |M(R,S) Δ M(R,T)| + |{i ∈ M(R,S) ∩ M(R,T) : S_i ≠ T_i}| (case proof in `synthetic_expected.md` §2).
  S R_k T ⇔ d_H(S,T) ≤ k. Composition: R_j ∘ R_k ⊆ R_{j+k}, by the triangle inequality.
- **Algorithms.** One linear scan per pair gives M (O(n)). Set operations are written as explicit membership loops, so the trace
  can show each element's decision. Relation tests: reflexive and symmetric in O(|V|²); transitivity in O(|V|³), returning a
  counterexample triple. Predicates are Python callables combined with ∧, ∨, ¬ and evaluated per sequence into a truth table.
- **Output.** Mutation sets, cardinalities, ∩ ∪ \ Δ, predicate truth tables and filters, a relation-property report with witnesses,
  and equivalence classes.
- **Why appropriate.** For aligned, equal-length sequences the comparison literally *is* set-valued.
  Hamming distance counts substitutions exactly.
- **Limitations.** Requires positional homology (equal length, no indels). Substitutions only. M depends on the choice of R, but
  d_H(S,T) does not. A position-level set mixes different alleles at the same position, which is why M* exists.

| Demo | Tag |
|---|---|
| Mutation sets, cardinality, ∩ ∪ \ Δ, shared vs specific | core |
| Allele-level M* vs position-level M (syn03/syn07 example) | core |
| Predicates and Boolean filters | core |
| R_k: reflexive, symmetric, *not* transitive (counterexample) | core |
| Identity (d_H = 0) as an equivalence relation → haplotypes | core |
| Orientation equivalence S ~ RC(S) | core for k-mers; **pedagogical on Dataset B full CDS** (all are in annotated sense orientation) |
| Functions: S as a map [n]→Σ; enc as a bijection with inverse dec | core |
| *Proposed:* transitive closure of R_k (Warshall) = "same component of G_k" | core bridge to Unit V **[decision]** |
| *Proposed:* poset (M*, ⊆) with Hasse diagram (the syn02 ⊂ syn03 ⊂ syn04 chain) | pedagogical-plus **[decision]** |

## Unit II — Recursion and Counting
- **Problem.** How large is the space of sequences and of k-substitution mutants, and when is exhaustive enumeration feasible?
- **Objects.** Σ^n; the Hamming sphere S(X,k) and ball B(X,r); counts M(n,k).
- **Equations.** N(n) = 4^n. M(n,k) = C(n,k)·3^k (choose k positions, then 3 alternatives each; substitutions only, no indels).
  M(n,k) = M(n−1,k) + 3·M(n−1,k−1): either position n is unchanged, or it is one of the k changed positions with 3 choices.
  Base case M(0,0) = 1, and M = 0 outside 0 ≤ k ≤ n. Σ_k M(n,k) = (1+3)^n = 4^n (binomial theorem; a free consistency check).
  |B(X,r)| = Σ_{k≤r} M(n,k). Inclusion–exclusion for "uses all four bases": 4^n − 4·3^n + 6·2^n − 4 (24 at n = 4).
- **Algorithms.** Memoised recurrence table, O(n·k) time. Recursive generator over positions, which either keeps the base or
  substitutes one of 3: Θ(n·M(n,k)) time to emit all strings, O(n) stack. Before generating, a guard computes M(n,k) from the
  closed form and refuses above a fixed cap (proposed 10^5), with a warning that shows the count.
- **Output.** Growth table for 4^n, M(n,k) tables (closed form vs recurrence), generated lists with a duplicate check, complexity statement.
- **Why appropriate.** It quantifies how sparse observed data are. For n = 1701 a single sequence has 3n = 5,103
  one-substitution neighbours, yet Dataset B has only 8 observed d_H = 1 pairs.
- **Limitations.** Counts possibilities, not probabilities. It ignores transition/transversion bias and biological viability:
  most of Σ^n contains premature stop codons.

| Demo | Tag |
|---|---|
| 4^n growth; M(n,k); recurrence; recurrence = closed form | core |
| Bounded recursive generation + count match (Exp 2) | core |
| Hamming-ball sizes vs observed neighbours | core |
| Permutations: k! shortest mutation orders between sequences at distance k in the Hamming graph | core (links to Unit V) |
| Induction proof that the closed form satisfies the recurrence (Pascal's rule) | pedagogical |
| Inclusion–exclusion, all four bases | pedagogical |
| Pigeonhole: a sequence of length L contains a repeated k-mer once L − k + 1 > 4^k | pedagogical (true and checkable) |
| Codon count 4^3 = 64 | pedagogical context |

## Unit III — Algebraic Structures
- **Problem.** DNA is double-stranded, so the same molecule can be written as S or RC(S). Analyses of k-mers of unknown
  orientation must treat the two as one object.
- **Objects.** Maps I, R (reverse), C (complement), RC on Σ^n. The groups G2 = ({I,RC}, ∘) ≅ Z2 and
  V4 = ({I,R,C,RC}, ∘) ≅ Z2 × Z2 (Klein). Their action on Σ^n, with orbits, fixed points and stabilisers; cosets of G2 in V4
  ({I,RC} and {R,C}). Complement as the permutation (A T)(C G) ∈ S4.
- **Equations.** RC∘RC = I; R∘C = C∘R = RC; R² = C² = I. |Fix(RC)| = 4^(n/2) for even n and 0 for odd n.
  |Fix(R)| = 4^⌈n/2⌉; |Fix(C)| = 0. Burnside: #orbits = (1/|G|)·Σ_g |Fix(g)|, which gives (4^n + Fix_RC)/2 for G2
  (10 at n = 2) and (4^n + 4^⌈n/2⌉ + Fix_RC)/4 for V4 (6 at n = 2). Both are verified by brute force in `synthetic_expected.md` §6.
- **Algorithms.** Transformations in O(n). Cayley tables are built by composing the maps and identifying each result. Function
  equality is decided **on all of Σ^n for n ≤ 6** (exhaustive), plus a short general-n proof, never on a random sample.
  Axiom checks: closure, identity and inverses from the table; associativity over all |G|³ triples. This is guaranteed anyway,
  because function composition is associative; the check is evidence, not proof.
- **Output.** Composition and Cayley tables, an axiom report, orbits of the dataset, fixed-point and orbit counts with brute-force agreement.
- **Why appropriate.** RC is a real involution corresponding to reading the opposite strand 5'→3'.
- **Limitations.** R and C on their own are *not* natural strand operations, so V4 is a mathematical extension
  (**pedagogical**). Mutations do not form a group: the point-substitution map "set S_i := b" is idempotent and not injective,
  so it has no inverse. These maps generate a *monoid* of transformations, not a group. This gives the spec's "do not pretend
  mutations form a group" a proof instead of an assertion. No arbitrary position permutations, per your decision.

| Demo | Tag |
|---|---|
| RC, RC∘RC = I, G2 Cayley table and axioms | core |
| Orbits under G2 = orientation classes; canonical representative | core (k-mers) |
| RC-fixed count and Burnside orbit count for G2 | core (how many canonical k-mers exist) |
| Klein four-group V4, its Cayley table, cosets of G2 | pedagogical |
| Point substitutions form a monoid, not a group (proof) | core (a limitation stated as mathematics) |
| Free monoid (Σ*, concatenation); RC(xy) = RC(y)RC(x) (anti-automorphism) | pedagogical |

## Unit IV — Number Theory ("Genomic Encoding and Modular Indexing Laboratory")
- **Problem.** Represent fixed-length words (k-mers) as integers for indexing, and understand exactly when modular bucketing
  loses information.
- **Objects.** enc: Σ^n → {0,…,4^n−1}, digits A0 C1 G2 T3; dec as its inverse; buckets b_m(x) = x mod m in Z_m; CRT moduli.
- **Equations.**
  enc(S) = Σ s_i·4^(n−i). It is a bijection, by uniqueness of base-4 representation. dec repeatedly applies the division
  algorithm x = 4q + r, which is a genuine use of it.
  enc(comp(S)) = 4^n − 1 − enc(S). enc(RC(S)) = the reversed digits of enc(comp(S)).
  Canonical k-mer: c(S) = min(enc(S), enc(RC(S))), a representative of each G2-orbit.
  **Bad hash:** if m | 4^j, then b_m(S) depends only on the last j bases, because 4^(n−i) ≡ 0 (mod m) for n − i ≥ j.
  **Good-hash fact (proved here):** if gcd(m, 4) = 1 and m > 3, any single substitution changes the bucket. A substitution at
  position i shifts x by δ·4^(n−i) with δ ∈ {±1, ±2, ±3}. Since 4 is invertible mod m, m divides this only if m | δ, which is
  impossible because 0 < |δ| < m. Dataset A agrees for m = 7, 13 and 17 on all six d = 1 pairs, and fails for m = 16
  (syn01 and syn02 collide). The hand-checked pairs are in `synthetic_expected.md` §7.
  *Stage 4 update:* this condition is sufficient but not necessary. The exact criterion depends only on the odd part q
  of m: there is no collision at any length iff q ≥ 5 (so even m = 10 also works). It is proved in
  `unit4_expected.md` §3 and checked by brute force for m < 70, n ≤ 5.
  Pigeonhole: N distinct keys in m < N buckets force a collision.
  CRT: for pairwise coprime m_1..m_t with product P, x ∈ [0, P) is uniquely recovered from its residues. 3-mers need P ≥ 64.
  Extended Euclid computes the inverses used by CRT.
- **Algorithms.** Horner encode in O(n); dec in O(n); the sliding-window update x' = (4x mod 4^k) + d_new gives O(1) per k-mer.
  Euclid takes O(log min(a,b)) steps; CRT is O(t) inverse computations.
- **Output.** Encoder/decoder, bucket histograms and collision counts per modulus, theorem checks on real d = 1 pairs, a CRT
  reconstruction table, and Euclid traces.
- **Why appropriate.** Canonical k-mer counting is documented practice (source below). The hashing theorems explain
  mathematically why one modulus is bad and another good.
- **Limitations.** Sharing a bucket says nothing about similarity (syn07 and syn03 share bucket 13 mod 16 at d_H = 3).
  Modular indexing is not presented as a clinical or established genomics technique. CRT and Euler/Fermat have no genomic necessity.

| Demo | Tag |
|---|---|
| enc/dec bijection; division algorithm in dec | core |
| Complement identity 4^n − 1 − enc | core |
| Canonical k-mer min(enc, enc∘RC) | core |
| Buckets, collisions, pigeonhole bound; m \| 4^j theorem; gcd(m,4) = 1 theorem | core (as representation analysis) |
| GCD, Euclid, extended Euclid, modular inverse | core *inside* CRT; otherwise pedagogical |
| CRT reconstruction of small k-mer codes | pedagogical |
| Euler φ / Fermat: ord_p(4) divides p − 1, so positional weights 4^t mod p repeat with that period (e.g. ord_7(4) = 3, ord_17(4) = 4) | pedagogical |

**Source for canonical k-mers (located and read).** G. Marçais and C. Kingsford, *Jellyfish: A fast k-mer counter*,
user manual, version 0.9 (2010-10-01), https://www.cbcb.umd.edu/software/jellyfish/jellyfish-manual.html.
Under the `--both-strands` option of `count`, the manual defines the canonical representation as whichever of a k-mer and its
reverse complement comes first in lexicographic order. It recommends the option when input orientation is unknown, for
example sequencing reads. Our numeric minimum is the same thing: for equal-length strings with A<C<G<T, comparing encodings
compares the first differing digit, which is lexicographic order.
**[flag]** This is the 2010 v0.9 manual. I did not verify whether current releases (where the option is `-C`) word it the same
way, and I did **not** locate the Jellyfish journal paper, so it is not cited. I also have no located source for "2-bit
encoding is standard in k-mer tools", so the project will not make that claim.

## Unit V — Graph Theory (principal component)
- **Problem.** Which observed haplotypes are linked by few substitutions, how that structure changes with the threshold, and
  what the cheapest routes through *observed* haplotypes cost.
- **Objects.** G_=1 (edge iff d_H = 1) and G_k (edge iff 1 ≤ d_H ≤ k, weight d_H) on haplotypes; adjacency and incidence
  matrices; degrees; components; κ (vertex) and λ (edge) connectivity.
- **Equations.** Handshake: Σ deg = 2|E|. G_j ⊆ G_k for j ≤ k, so #components is non-increasing in k. Components of G_k are the
  classes of the transitive closure of R_k. Whitney: κ ≤ λ ≤ δ. Triangle inequality: dist_{G_k}(u,v) ≥ d_H(u,v), with equality
  iff some path never revisits a position (Dataset A, syn07 → syn03: 4 > 3). **Bottleneck property:** the smallest k with G_k
  connected equals the largest edge in any MST (Dataset A: 6 and 6).
  G_=1 is an induced subgraph of the Hamming graph H(n,4). Unlike the binary hypercube, H(n,4) contains triangles
  (three bases at one position), so G_=1 need not be bipartite.
- **Algorithms.** All pairwise distances: O(N²n). BFS/DFS: O(V+E). Components by BFS. Dijkstra on the sparse G_k, using the
  O(V²) array version (no library heap; adequate for V ≤ 40). κ and λ by exhaustive removal, guarded to small graphs.
  Cycle detection by DFS back-edges.
- **Output.** Graph drawings, matrices, degree tables, component sweeps over k, path results with cost vs direct d_H,
  and complexity notes.
- **Why appropriate.** A d_H = 1 edge is an exact, checkable statement: the two observed sequences differ by one substitution.
- **Limitations.** Only observed sequences are vertices. A missing edge may mean an unsampled intermediate. A path is not an
  evolutionary trajectory. Any single k is arbitrary, so results are always reported as a sweep.

| Demo | Tag |
|---|---|
| G_=1, G_k, degrees, matrices, BFS/DFS, components, k-sweep | core |
| Dijkstra on sparse G_k with cost vs direct d_H | core |
| Bottleneck property and transitive-closure bridge | core |
| κ, λ, Whitney inequality | pedagogical (small graphs only) |
| Incidence matrix | pedagogical |

## Unit VI — Trees
- **Problem.** (a) Find the cheapest structure that connects all haplotypes. (b) *Separately*, build a hierarchical
  distance-based hypothesis.
- **Objects.** K_N weighted by d_H; spanning trees; MST (non-unique under ties); a rooted version of the MST, rooted at the
  reference **arbitrarily**; the UPGMA rooted binary tree with node heights.
- **Equations.** |E(T)| = N − 1. Cut and cycle properties. Number of MSTs via weight levels (8 for Dataset A).
  Cayley: n^(n−2) labelled trees. (2n−3)!! rooted binary topologies with n labelled leaves.
  UPGMA merge distance: d(A∪B, W) = (|A|·d(A,W) + |B|·d(B,W)) / (|A| + |B|), with node height = d/2.
- **Algorithms.** Prim: O(V²) array. Kruskal: O(E log E), with a hand-written sort-and-union-find (path compression, union by
  rank). Pre/post-order on rooted trees, and in-order **only** on the binary UPGMA tree, all O(V). UPGMA: naive O(N³).
- **Output.** MST edge lists and weights from both algorithms, validity checks (N−1 edges, acyclic, spanning), rooted-tree
  statistics (depth, height, leaves, internal nodes), traversals, the UPGMA dendrogram, and a side-by-side comparison.
- **Why appropriate.** An MST is a genuine optimisation object on the distance data. UPGMA is a standard distance-based
  clustering method.
- **Limitations.** An MST is **not** a phylogeny: its internal nodes are observed sequences, it is not unique under ties, and
  its root is arbitrary. UPGMA assumes ultrametric distances, i.e. a molecular clock with equal rates, which is generally
  violated. Its output is a hypothesis, not history. Raw p-distance (no multiple-hit correction) is acceptable at ≤ 3.4%
  divergence but is not an evolutionary model. **[decision]** UPGMA needs a fixed tie-break rule; I propose merging the pair
  whose sorted label pair is lexicographically smallest.

| Demo | Tag |
|---|---|
| Prim and Kruskal, equal weights, validity checks, tie case | core |
| MST vs UPGMA comparison, with assumptions stated | core |
| Rooted-tree statistics and traversals | core (pedagogical on the MST, since the root is arbitrary) |
| Bottleneck link to Unit V | core |
| n^(n−2) and (2n−3)!! | pedagogical |

---

## Relevance ratings and evidence needed

| Unit | Mathematical relevance | Role in the research question | Evidence that would justify it |
|---|---|---|---|
| I Set theory & logic | **High** | foundation: the sets, relations and equivalences everything else uses | set results match the hand-derived tables; non-transitivity counterexample found automatically; d_H identity holds on all pairs of A and B |
| II Counting & recursion | **High** | sizes the space; justifies enumeration limits; k! path count | enumeration = C(n,k)·3^k for all n ≤ 8; recurrence = closed form; guard triggers above the cap |
| III Algebraic structures | **Medium** | real strand symmetry; operative on k-mers, not on annotated CDS | exhaustive Cayley table for n ≤ 6; Burnside = brute force for n ≤ 8; canonical k-mer counts on Dataset B |
| IV Number theory | **Medium** (high as representation maths, low as genomics) | k-mer representation and hashing behaviour | dec(enc(S)) = S exhaustively for n ≤ 6; both hashing theorems verified on every d = 1 pair of A and B; CRT round-trip on all 64 3-mers |
| V Graph theory | **High** | principal analysis | handshake, component sweep and Dijkstra match hand values on A; match NetworkX on B; bottleneck equality on A and B |
| VI Trees | **High** | summary structures; the MST-vs-phylogeny contrast | Prim weight = Kruskal weight = NetworkX; validity checks; UPGMA ultrametric check; documented contrast on A and B |

## Weak links and proposed alternatives
1. **Unit III does no work on Dataset B as full CDS.** All sequences are in annotated sense orientation.
   *Fix:* apply G2 to the **k-mers** of Dataset B. Report distinct k-mers vs distinct canonical k-mers and the number of
   RC-palindromic k-mers (even k), and compare against Burnside's whole-space count. Proposed as part of Exp 4 **[decision]**.
2. **CRT and Euler/Fermat are pedagogical.** Keep them small and labelled. CRT round-trips 3-mer codes, and Fermat explains
   positional-weight periodicity mod p.
3. **"Inorder traversal"** is defined only for binary trees. Use the UPGMA tree, not the MST.
4. **Vertex and edge connectivity** by exhaustive removal is exponential. Restrict it to Dataset A and small subgraphs and say so.
5. **Inclusion–exclusion and pigeonhole** are on-topic but not needed. Keep one worked example each.
6. **Spec's "paths and cycles"**: cycles in G_=1 need two mutations observed in both orders, like the Dataset A square. On real
   data there may be none, which is a valid finding rather than a failure.

## Items flagged for your decision
- [decision] Add M* (allele-level sets) alongside M. I recommend yes; it is needed to state "shared mutation" correctly.
- [decision] Add Warshall's transitive closure of R_k as the Unit I ↔ V bridge.
- [decision] Add the (M*, ⊆) poset / Hasse diagram (optional).
- [decision] Canonical-k-mer analysis of Dataset B (see weak link 1), and the k value(s) to use. I suggest k ∈ {4, 6, 8}.
- [decision] UPGMA tie-break rule (lexicographic, as proposed above).
- [decision] Dataset B reference sequence. I propose the **medoid** haplotype (minimum Σ d_H; ties broken by smallest
  accession), because it is an observed sequence and is defined mathematically.

## Claims I am not certain of
- The Jellyfish manual wording is from v0.9 (2010). See the flag above.
- I state no clade or lineage names for Dataset B. The distance histogram is multimodal (gaps between 22–27 and 37–47 nt),
  but naming clades would need a reference I have not located.
- "3n one-substitution neighbours" and the Hamming-graph triangle remark are elementary and were checked by hand. They are
  not sourced, because they are mathematics rather than citations.
