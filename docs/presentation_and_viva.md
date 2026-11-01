# Deliverable G — Presentation and viva preparation

## Presentation outline (12–15 minutes)

| # | Slide | Content | Time |
|---|---|---|---|
| 1 | Title | Project, course, name. One line: "All six units, applied to one problem, with every result traceable." | 0:30 |
| 2 | Why DNA is a discrete object | Σ = {A,C,G,T}; a sequence is a function [n] → Σ. Not an approximation — this *is* the data. | 1:00 |
| 3 | The question | Not "what does the tool output" but "which mathematical object is computed, from which definition, with what assumptions". | 1:00 |
| 4 | Two datasets | A: 10 synthetic sequences, every answer hand-derived first. B: 36 GenBank H3N2 HA coding sequences, 1701 nt. | 1:00 |
| 5 | Method | Derive by hand → implement → compare in separate columns. Independent vs consistency checks. 147 independent matches. | 1:30 |
| 6 | Units I–II | Mutation sets; the distance identity; R_k is **not** transitive; M(n,k) = C(n,k)·3^k and its recurrence. | 2:00 |
| 7 | **The sparsity result** | 5,103 possible neighbours, 6 observed. This single counting fact constrains every later interpretation. | 1:00 |
| 8 | Unit III | Two Klein four-groups; Burnside; and the proof that substitution operators form a **monoid, not a group**. | 1:30 |
| 9 | Unit IV + **the bad hash** | Encoding is a bijection; m = 16 sees only the last 2 bases; on Dataset B all 30 haplotypes land in **one bucket** (every sequence ends TGA). Theorem: odd part ≥ 5 ⟹ never. | 2:00 |
| 10 | Units V–VI | Threshold sweep; k\* = 47 = largest MST edge; 4 distinct MSTs on real data. | 1:30 |
| 11 | **MST ≠ phylogeny** | Same distances, 5 of 10 and 27 bipartitions shared. Rooting at syn01 vs syn04: height 3 → 4, same 8 edges. | 1:30 |
| 12 | Live demo | Dashboard: pick a result, open its trace, show the five category tags. | 1:30 |
| 13 | Validation | 179 tests, 67/67 mutants caught, library oracles. Ten errors found during development, two in my own hand derivations. | 1:00 |
| 14 | Limitations | Substitutions only; no alignment; convenience sample; not a diagnostic or phylogenetic tool. | 0:45 |
| 15 | Conclusion | Six units doing real work; two theorems proved; results reported without overclaiming. | 0:30 |

**If time runs short**, cut slides 6 and 10 to one each. Never cut 7, 9 or 11 — they are the project's three results
that are genuinely about mathematics deciding an interpretation.

---

## Demonstration script (about 4 minutes)

Start with the dashboard already running (`.venv/Scripts/python -m streamlit run app/main.py`), Dataset A selected.

1. **Section 1 (20s).** "Ten records, nine haplotypes — syn03 and syn10 are identical, so they collapse. That collapse
   is an equivalence relation: d_H = 0 is reflexive, symmetric and transitive."
2. **Section 2 (50s).** Pick syn03 and syn07. "They share position 5 — but not the same base there. The allele-level
   sets have empty intersection. And the distance identity: 2 from the symmetric difference, plus 1 from the shared
   position carrying different bases, equals 3." Then set k = 1 and point at the relation table: "reflexive yes,
   symmetric yes, transitive **no**, and here is the counterexample the program found."
3. **Section 4 (40s).** Show the Cayley table and the axiom table. Then scroll to the monoid: "setting both positions
   to A sends AA and AC to the same string, so it has no inverse. That is the proof that mutations do not form a group."
4. **Section 5 (50s).** Switch the dataset to **B**. Show the bucket table: "m = 16 — one bucket, thirty haplotypes.
   They all end TGA, and this modulus sees only the last two bases. Every odd modulus separates all eight
   distance-1 pairs, and the theorem says exactly why."
5. **Section 6 (40s).** "k\* = 47 — that is when the graph first connects." Show the sweep plot.
6. **Section 7 (60s).** "Prim and Kruskal both give 171. But they returned *different trees*, because there are four
   minimum trees. Largest edge 47 — the same number as k\*." Then the dendrogram: "UPGMA assumes a molecular clock.
   These distances are not ultrametric, so it reproduces only 9 of 435 exactly. The two trees share 5 bipartitions out
   of 10 and 27. They answer different questions."
7. **Close (20s).** Open any trace expander. "Every result on this page has this, and every step is tagged by what kind
   of claim it is."

---

## Viva questions and answers

### A. The MST / phylogenetic tree distinction

**Q1. Your MST connects sequences by mutational distance. Why is it not a phylogenetic tree?**
Three reasons, each structural rather than a matter of degree.
(i) *Its internal nodes are observed sequences.* An MST connects what was sampled; a phylogeny posits unobserved common
ancestors at its internal nodes. UPGMA on the same data creates nodes U1…U29 that are not in the dataset.
(ii) *It optimises the wrong thing.* An MST minimises the total of all edge weights. Evolution does not minimise a
global sum; and nothing about minimum total weight corresponds to a claim about descent.
(iii) *It is not unique and not rooted.* Dataset A has 8 minimum spanning trees of weight 15 and Dataset B has 4. Any
root may be chosen: rooting Dataset A's MST at syn01 gives height 3 and 6 leaves, at syn04 height 4 and 5 leaves — the
same eight edges. A tree whose root and shape are presentation choices cannot be asserting an ancestral order.

**Q2. Then why compute an MST at all?**
Because it answers a real and different question: what is the cheapest set of links connecting every observed
sequence? It also yields the bottleneck property — the largest MST edge equals the smallest threshold at which the
graph connects (6 on A, 47 on B) — which links Units V and VI exactly.

**Q3. Is the UPGMA tree the evolutionary history, then?**
No. It is a clustering hypothesis under a stated assumption: a molecular clock, meaning equal rates on all lineages,
which makes distances ultrametric. Dataset B's distances are *not* ultrametric, and the tree reproduces only 9 of 435
distances exactly. Its branch heights are averages, not times. It is also sensitive to a tie rule: 11 merge steps were
ties, and the reverse rule changes the topology.

**Q4. How different are the two trees in practice?**
On Dataset B they share 5 non-trivial bipartitions out of the MST's 10 and UPGMA's 27. A bipartition is what you get
by cutting one edge and recording the split. Agreeing on 5 of a possible 27 is weak agreement, which is what one
should expect from two structures built for different purposes.

**Q5. What would you need to make a defensible evolutionary claim?**
A substitution model correcting for multiple hits at one position, a method whose assumptions match the data (neighbour
joining drops the clock assumption; likelihood or Bayesian methods state a model explicitly), branch support such as
bootstrapping, an outgroup to root on, and sampling that is not a convenience sample. This project has none of those
and claims none of them.

### B. Group axioms and Unit III

**Q6. State the group axioms and verify them for {I, RC}.**
Closure, associativity, identity, inverses. RC∘RC = I because reverse-complementing twice restores both the order and
the bases, so the set is closed; composition of functions is always associative; I is the identity; and each element is
its own inverse. Order 2, so it is isomorphic to Z2.

**Q7. Why is {I, R, C, RC} the Klein four-group and not Z4?**
Both have order 4, but in the Klein group every non-identity element has order 2, while Z4 has an element of order 4.
Here R² = C² = (RC)² = I, so no element has order 4. The program checks this two ways: element orders are 1, 2, 2, 2,
and none of the 24 bijections onto Z4 is a homomorphism.

**Q8. Why do mutations not form a group?**
Model a substitution as the operator "set position i to base b". These compose — a later write wins — and the identity
is among them, so they form a monoid, with 5^n elements for length n. But such an operator is not injective: setting
both positions of a length-2 string to A sends AA and AC to the same string. A non-injective map has no inverse, so the
only invertible element is the identity. A monoid, not a group. The program generates the monoid and produces that
witness.

**Q9. There is a subtlety in your group action. What is it?**
Faithfulness. As maps on Σ^n the four elements are distinct only for n ≥ 2. On strings of length 1 reversal does
nothing, so R = I and RC = C; on the empty string all four coincide. Burnside's lemma does not require a faithful
action, so the orbit counts stay correct at n = 0 and 1 — but "four distinct maps" is false there. That is why the
program decides function equality over all 5,461 words of length ≤ 6 rather than on one length.

**Q10. State Burnside's lemma and use it.**
The number of orbits equals the average number of fixed points: (1/|G|)·Σ_g |Fix(g)|. For {I, RC} on Σ²: I fixes all
16 strings, and RC fixes 4 (choose the first base freely, the second is forced), so (16+4)/2 = 10 orbits. The program
agrees with direct enumeration for every n from 0 to 8.

**Q11. What does the second Klein group, the one on bases, buy you?**
It classifies substitutions. The three fixed-point-free involutions of a four-element set — there are 3!! = 3 — each
move all four bases, so each realises exactly 4 ordered pairs; 3 × 4 = 12 covers every ordered pair, meaning every
substitution comes from exactly one element. (A G)(C T) pairs purine with purine and pyrimidine with pyrimidine: the
transitions. And since each element is an involution, the classification works on *unordered* pairs — which matters,
because two observed sequences never say which base came first.

### C. Recurrence relations and counting

**Q12. Derive M(n,k) = C(n,k)·3^k and its recurrence.**
Choose which k of the n positions change — C(n,k) ways — then one of the 3 other bases at each, giving 3^k.
For the recurrence, split on the last position: either it is unchanged and the first n−1 positions carry all k changes,
giving M(n−1,k); or it is changed in 3 ways and the first n−1 carry k−1, giving 3·M(n−1,k−1). Disjoint and exhaustive.

**Q13. Prove the closed form satisfies the recurrence.**
Induction on n. Base: M(0,0) = 1 = C(0,0)·3⁰. Step: C(n−1,k)·3^k + 3·C(n−1,k−1)·3^(k−1) = 3^k[C(n−1,k) + C(n−1,k−1)]
= 3^k·C(n,k) by Pascal's rule. A free check: Σ_k C(n,k)3^k = (3+1)^n = 4^n, which every row of the program's table
satisfies.

**Q14. Where does the pigeonhole principle appear, twice?**
First on k-mers: a string of length L has L − k + 1 windows but only 4^k possible k-mers, so a repeat is forced when
L − k + 1 > 4^k. For L = 1701 the largest such k is 5, and all 34 clean sequences do contain a repeated 5-mer.
Second on hash buckets: N keys in m buckets force a bucket with at least ⌈N/m⌉.

**Q15. Why is your sparsity number important?**
Because it decides how every later picture should be read. A 1701-base sequence has 3n = 5,103 possible
single-substitution neighbours; the best-connected haplotype has 6. The graph is a very thin sample, so a *missing*
edge is weak evidence of anything. That is a counting argument constraining interpretation.

### D. Modular arithmetic and Unit IV

**Q16. Why is base-4 encoding a bijection, and how do you invert it?**
Every integer in [0, 4^n) has exactly one n-digit base-4 representation with leading zeros allowed, so the map is
injective and surjective. Inversion is the division algorithm applied n times: x = 4q + r gives the last digit r, then
repeat on q.

**Q17. Why is m = 16 a bad hash for sequences?**
Because 16 = 4², every positional weight 4^t with t ≥ 2 is ≡ 0 mod 16, so the bucket depends only on the last two
bases. On Dataset B this is dramatic: every sequence ends TGA, so all 30 haplotypes land in one bucket. All eight
distance-1 pairs collide.

**Q18. Give the exact criterion for when a single substitution preserves the bucket.**
A substitution at position i changes the code by δ·4^t with δ ∈ {±1, ±2, ±3} and t = n − i; the bucket survives iff
m | δ·4^t. Write m = 2^a·q with q odd. Since gcd(q, 4^t) = 1 we need q | δ, which with |δ| ≤ 3 forces q ∈ {1, 3}.
**So if the odd part of m is at least 5, it never happens, at any length.** Otherwise it happens exactly when a
position with t ≥ t_min exists, where t_min = ⌈a/2⌉ for q = 3 and ⌊a/2⌋ for q = 1, since δ = ±2 contributes its own
factor of 2. Note this is sharper than requiring gcd(m,4) = 1 — m = 10 also works.

**Q19. Can two substitutions hide from a good modulus?**
Yes, and Euler and Fermat say where. For odd m the weights 4^t repeat with period ord_m(4), which divides φ(m), and
divides p−1 for prime p. An opposite pair of changes cancels exactly when that order divides the gap between the
positions. Changing syn01 at positions 6 and 12 shifts the code by 1 − 4⁶ = −4095 = −7·585: distance 2, yet the same
bucket mod 7 (order 3) and mod 13 (order 6), but not mod 17, 31 or 61.

**Q20. Does sharing a bucket mean two sequences are similar?**
No, and that is the point of the experiment. The residue depends on the whole integer code. Under m = 16 on Dataset A,
syn07 shares a bucket with syn03 at distance 3, while syn01 and syn05 at distance 1 are in different buckets.

**Q21. State the Chinese Remainder Theorem and how you tested it.**
For pairwise coprime moduli with product P, any integer in [0, P) is uniquely recoverable from its residues. 3-mer
codes lie in [0, 64), so moduli with product ≥ 64 suffice: (7, 11) and (5, 13) recover all 64 codes. (3, 5) has product
15 and recovers only 15 of 64 — which is the uniqueness range doing exactly what the theorem says. Non-coprime moduli
such as (4, 6) are refused.

### E. Graph algorithms

**Q22. Why does Dijkstra apply here, and what is its complexity in your implementation?**
Weights are Hamming distances, so they are non-negative — Dijkstra's requirement. I use the array version: each round
scans for the closest unsettled vertex, giving O(V² + E). For 30 vertices that is preferable to a heap, and it keeps
the algorithm visible rather than hidden in a library.

**Q23. When does a shortest path cost more than the direct distance, and why can it never cost less?**
It can never cost less because d_H is a metric: the triangle inequality makes any path's total at least the direct
distance. It costs *more* when every available route changes some position twice. In Dataset A's G_2, syn04 → syn07
costs 5 against a direct distance of 4, because the only route out of syn07 passes through syn01 and pays for
position 5 twice.

**Q24. What is the bottleneck property, and did it hold?**
The smallest threshold k at which the graph becomes connected equals the largest edge of any minimum spanning tree.
It held on both datasets: 6 on Dataset A and 47 on Dataset B.

**Q25. Prim and Kruskal gave different trees on the real data. Is one wrong?**
Neither. The minimum *weight* is unique, but the minimum *tree* need not be when weights tie. Dataset B has 4 distinct
minimum spanning trees, all of weight 171, and the two algorithms break ties differently. This is why the tests compare
weights and validity rather than edge sets. I count the minimum trees by contracting weight levels and applying
Kirchhoff's matrix-tree theorem at each level.

**Q26. How do you know your graph code is right?**
Three ways. Hand-derived values for Dataset A — the adjacency matrix, degrees, BFS and DFS orders, the cycle, the
sweep, and Dijkstra's full distance table — were written down before the code. NetworkX acts as an independent referee
on both datasets and on 60 random graphs. And the mutation check injects bugs such as Warshall skipping its last pivot
or Dijkstra settling the farthest vertex, and requires the suite to catch them.

### F. Why each unit is included

**Q27. Justify each of the six units honestly, including the weak ones.**

| Unit | Why it is here | Honest status |
|---|---|---|
| I Sets and logic | Sequences *are* functions [n] → Σ and their differences *are* sets. Supplies the haplotype quotient, the distance identity, and the non-transitivity that explains gradual component merging. | **Core.** Nothing else works without it. |
| II Counting | Sizes the space, sets the enumeration guards, and supplies the sparsity argument that constrains every interpretation. | **Core.** |
| III Algebra | Reverse complement is a genuine strand symmetry; it gives canonical k-mers, and the monoid proof bounds what "mutation algebra" may claim. | **Medium.** Does no work on full annotated coding sequences — predicted and then measured. Real contribution is at the k-mer level and in the substitution classification. R and C alone are pedagogical. |
| IV Number theory | Encoding, the complement identity, canonical representatives, and two proved theorems about what a hash preserves. | **Medium.** Strong as representation mathematics, weak as genomics. CRT, Euler and Fermat are labelled pedagogical. |
| V Graphs | The principal analysis: the structure of "few substitutions apart" and how it changes with the threshold. | **Core.** |
| VI Trees | The summary structures, and the MST-versus-phylogeny distinction that the whole project turns on. | **Core.** |

**Q28. Which part is weakest, and what would you do about it?**
Unit III on full sequences. Every Dataset B sequence begins ATG, so its reverse complement would begin with the
complement of a stop codon and cannot be in the dataset — meaning the orientation group does nothing at that level. I
predicted this before running it and then measured it (0 reverse-complement pairs). I addressed it two ways: applying
the group to k-mers, where 2,563 distinct 8-mers collapse to 2,488 canonical classes, and adding the base-level group
that classifies every observed substitution. Both are honest uses; neither pretends the full-sequence application
works.

### G. Method and general

**Q29. How do you know your code is correct and not just self-consistent?**
Expected values were derived by hand *before* implementation, and each experiment prints expected and actual in
separate tables. I also distinguish *independent* checks — against a hand derivation, a hand-proved theorem, or a
quality-control script sharing no code with the implementation — from *consistency* checks that can only fail if the
code disagrees with itself, and I count them separately: 147 independent and 28 consistency. Library oracles referee
the tests. And the mutation check injects 67 deliberate bugs to confirm the tests can fail at all.

**Q30. Did the method actually catch anything?**
Ten errors. Two are the best illustration. My `fixed_point_free_involutions` returned 6 because it generated each
pairing twice — the hand derivation said 3!! = 3, so the disagreement was immediate. And my first statement of the
bucket criterion was wrong for δ = ±2; writing the proof before coding exposed it. It also worked in reverse once: a
hand-derived cophenetic value of 13/8 was wrong, the code said 12/7, and re-deriving confirmed the code.

**Q31. What is an equivalent mutant, and why does it matter?**
A deliberately injected change that cannot alter behaviour, so no test can catch it. I found three — for example,
swapping the composition order in the Cayley table of an *abelian* group leaves the table identical. I verified that
one by comparing both tables entry by entry rather than asserting it, then removed it with the reason recorded.
Requiring a test to catch an equivalent mutant would be requiring a test to detect nothing.

**Q32. Why is Dataset A artificial? Isn't real data better?**
They do different jobs. Dataset A is engineered so that every structural case — an identical pair, an RC-palindrome,
an MST tie, a disconnecting outlier — has an answer a human can check on paper. That is what makes it a *test*. Real
data cannot serve that purpose, because there is nothing independent to compare against. Dataset B then shows the
validated pipeline works at a realistic scale.

**Q33. What would you do differently?**
I would write the hand derivations even earlier. The two errors it caught were both found at the boundary where
derivation met implementation, and the one error it failed to prevent — the wrong cophenetic value — happened in a
section I derived in a hurry. I would also add rendering tests for the dashboard, which is the one place where a
regression could currently slip past the suite.
