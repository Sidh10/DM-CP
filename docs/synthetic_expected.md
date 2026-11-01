# Dataset A: hand-derived expected answers

Status: **expected values only**. Nothing here was produced by `genomath`. Every number was derived by hand
with the working shown below, then checked against an independent brute-force script that is not project code
(see "How these were checked" at the end). When the experiments run, they will report *actual* values next to these.

File: `data/synthetic/dataset_a.fasta` (10 sequences, all length n = 12).
Companion file `data/synthetic/dataset_a_invalid.fasta` holds two deliberately bad records for input-validation tests.

Conventions used throughout:
- Positions are numbered 1..12 from the left (5' end), matching the spec's examples such as {3, 6, 10}.
- Encoding A=0, C=1, G=2, T=3; the most significant digit is position 1.
- The reference is R = syn01.

## 1. The sequences and why each one is there

R = syn01 = `G A T T A C G T A A T C`

| id | sequence | allele-level mutation set M*(R,S) | role |
|---|---|---|---|
| syn01 | GATTACGTAATC | ∅ | reference; the only RC-palindrome |
| syn02 | GCTTACGTAATC | {(2,C)} | chain step 1; corner of the tie square |
| syn03 | GCTTGCGTAATC | {(2,C),(5,G)} | chain step 2 |
| syn04 | GCTTGCGTTATC | {(2,C),(5,G),(9,T)} | chain step 3 (end of chain) |
| syn05 | GATTACGTAAAC | {(11,A)} | corner of the tie square |
| syn06 | GCTTACGTAAAC | {(2,C),(11,A)} | corner of the tie square |
| syn07 | GATTCCATAATC | {(5,C),(7,A)} | same *position* 5 as syn03, different *base* |
| syn08 | GATTACGCAAGC | {(8,C),(11,G)} | equals RC(syn03); same {I,RC}-orbit as syn03 |
| syn09 | CAGCAAGTAGTA | {(1,C),(3,G),(4,C),(6,A),(10,G),(12,A)} | outlier: d_H ≥ 6 from everything |
| syn10 | GCTTGCGTAATC | {(2,C),(5,G)} | identical to syn03 |

Design checklist (each item is verified in later sections):
- two identical sequences: syn03 = syn10 (§3);
- one RC-palindrome: syn01 (§6);
- chain: syn01 – syn02 – syn03 – syn04, each step one substitution (§4);
- MST tie: the 4-cycle syn01–syn02–syn06–syn05 of weight-1 edges, plus a second tie for attaching syn08 (§5);
- outlier: syn09 is isolated for every threshold k ≤ 5 (§4).

## 2. Unit I: mutation sets and set operations

The position-level mutation set is M(R,S) = { i : S_i ≠ R_i }. Writing each changed position next to R:

```
pos    1 2 3 4 5 6 7 8 9 10 11 12
R      G A T T A C G T A A  T  C
syn03  G C T T G C G T A A  T  C   -> differs at 2 (A->C), 5 (A->G)          M = {2,5}
syn07  G A T T C C A T A A  T  C   -> differs at 5 (A->C), 7 (G->A)          M = {5,7}
```

Position-level operations on M3 = M(R,syn03) = {2,5} and M7 = M(R,syn07) = {5,7}:
- |M3| = 2, |M7| = 2
- M3 ∩ M7 = {5}
- M3 ∪ M7 = {2,5,7}
- M3 \ M7 = {2}, M7 \ M3 = {7}
- M3 Δ M7 = {2,7}

**Shared position does not mean shared mutation.** syn03 has G at position 5 and syn07 has C. Using the
allele-level set M*(R,S) = { (i, S_i) : S_i ≠ R_i }:
M*3 ∩ M*7 = {(2,C),(5,G)} ∩ {(5,C),(7,A)} = ∅.

Identity linking sets to Hamming distance (proved by checking the three cases for a position i:
in neither set ⇒ S_i = T_i = R_i, contributes 0; in exactly one ⇒ one equals R_i and the other does not, contributes 1;
in both ⇒ contributes 1 exactly when S_i ≠ T_i):

  d_H(S,T) = |M(R,S) Δ M(R,T)| + |{ i ∈ M(R,S) ∩ M(R,T) : S_i ≠ T_i }|

Example: d_H(syn03, syn07) = |{2,7}| + |{5}| = 2 + 1 = **3**.
The second term vanishes for a binary alphabet; it exists here because each position has three alternative bases.

Second example (nested sets, the chain): M(R,syn02) = {2} ⊂ M(R,syn03) = {2,5} ⊂ M(R,syn04) = {2,5,9}.
Shared with syn04: M(R,syn02) ∩ M(R,syn04) = {2}; syn04-specific relative to syn03: {9}.

### Predicates (Unit I logic)
On the 9 haplotypes, with P(S) ≡ d_H(R,S) ≤ 2 and Q(S) ≡ S_2 = C:
- P is true for syn01, syn02, syn03, syn05, syn06, syn07, syn08 (distances 0,1,2,1,2,2,2; read from row syn01 in §3).
- Q is true for syn02, syn03, syn04, syn06.
- P ∧ Q = {syn02, syn03, syn06}; P ∧ ¬Q = {syn01, syn05, syn07, syn08}; ¬P = {syn04, syn09}.

## 3. Distance matrix and haplotypes

syn10 is identical to syn03, so d_H(syn03, syn10) = 0. "Identical" (d_H = 0) is an equivalence relation
(reflexive, symmetric, transitive), and its classes are the haplotypes:
{syn03, syn10} and seven singletons, which gives **9 haplotypes**. Each haplotype is named after its first member.
From here on, the graph and tree work uses the 9 haplotypes.

Each entry is computed with the §2 identity from the mutation sets in §1.
Worked entries: d(syn05, syn08): positions {11} Δ {8,11} = {8}, plus shared position 11 with A vs G ⇒ 1 + 1 = 2.
d(syn06, syn08): {2,11} Δ {8,11} = {2,8}, plus 11 (A vs G) ⇒ 3. d(syn04, syn08): {2,5,9} Δ {8,11} = {2,5,8,9,11} ⇒ 5.
syn09's six positions are disjoint from every other set, so d(syn09, X) = 6 + |M(R,X)|.

```
       01 02 03 04 05 06 07 08 09
syn01   0  1  2  3  1  2  2  2  6
syn02   1  0  1  2  2  1  3  3  7
syn03   2  1  0  1  3  2  3  4  8
syn04   3  2  1  0  4  3  4  5  9
syn05   1  2  3  4  0  1  3  2  7
syn06   2  1  2  3  1  0  4  3  8
syn07   2  3  3  4  3  4  0  4  8
syn08   2  3  4  5  2  3  4  0  8
syn09   6  7  8  9  7  8  8  8  0
```

Distance histogram over the C(9,2) = 36 unordered pairs:

| d_H | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | total |
|---|---|---|---|---|---|---|---|---|---|---|
| pairs | 6 | 8 | 8 | 5 | 1 | 1 | 2 | 4 | 1 | 36 |

The six d_H = 1 pairs: 01–02, 01–05, 02–03, 02–06, 03–04, 05–06.

### Relation R_k: S R_k T ⇔ d_H(S,T) ≤ k
- Reflexive: d_H(S,S) = 0 ≤ k for every k ≥ 0.
- Symmetric: d_H is symmetric.
- **Not transitive for k = 1**: syn01 R_1 syn02 (d = 1) and syn02 R_1 syn03 (d = 1), but d(syn01, syn03) = 2, so ¬(syn01 R_1 syn03).
  (For k = 0 on haplotypes, R_0 is equality, which is transitive.)

## 4. Unit V: threshold graphs

G_k has the 9 haplotypes as vertices and an edge {S,T} iff 1 ≤ d_H(S,T) ≤ k, with weight d_H. Edge counts are cumulative sums of the histogram.

| k | new edges at weight k | |E| | components | component membership |
|---|---|---|---|---|
| 0 | – | 0 | 9 | all singletons |
| 1 | 6 | 6 | 4 | {01,02,03,04,05,06}, {07}, {08}, {09} |
| 2 | 8 | 14 | 2 | {01..08}, {09} |
| 3 | 8 | 22 | 2 | {01..08}, {09} |
| 4 | 5 | 27 | 2 | {01..08}, {09} |
| 5 | 1 | 28 | 2 | {01..08}, {09} |
| 6 | 1 | 29 | 1 | all |
| 9 | (all) | 36 | 1 | complete graph K9 |

Why k = 2 joins syn07 and syn08: d(syn01,syn07) = 2 and d(syn01,syn08) = 2, while both are at distance ≥ 2 from everything else.
Why syn09 stays isolated until k = 6: its minimum distance is d(syn01, syn09) = 6.

Degrees in G_1: syn01: 2 (02, 05); syn02: 3 (01, 03, 06); syn03: 2 (02, 04); syn04: 1; syn05: 2 (01, 06); syn06: 2 (02, 05); syn07, syn08, syn09: 0.
Sum = 12 = 2|E|, as required by the handshake lemma. G_1 contains exactly one cycle, syn01–syn02–syn06–syn05–syn01.

### Dijkstra on G_2 (sparse, weights d_H)
Adjacency of G_2 used below: syn07's only neighbour is syn01 (weight 2); syn08's neighbours are syn01 (2) and syn05 (2);
syn04's neighbours are syn03 (1) and syn02 (2).

| source → target | direct d_H | shortest-path cost in G_2 | optimal paths |
|---|---|---|---|
| syn04 → syn08 | 5 | **5** | 04–03–02–01–08 (1+1+1+2), 04–02–01–08 (2+1+2), 04–03–01–08 (1+2+2) |
| syn07 → syn03 | 3 | **4** | 07–01–03 (2+2), 07–01–02–03 (2+1+1) |
| syn07 → syn08 | 4 | **4** | 07–01–08 (2+2) |

Interpretation of the numbers: d_H is a metric, so by the triangle inequality every path cost is ≥ the direct d_H.
The cost is strictly larger (syn07 → syn03: 4 > 3) because the only route out of syn07 passes through syn01. That route
reverts position 5 from C to A and then mutates it to G, so position 5 is paid for twice. The syn04 → syn08 row has
three equally short paths, which is the "multiple optimal solutions" test case. In G_1, syn07 is unreachable from
syn03, so Dijkstra must return ∞ / no path.

## 5. Unit VI: minimum spanning tree (complete graph on the 9 haplotypes, weights d_H)

Kruskal, by weight level:
1. **Weight 1** (6 edges): 01–02, 01–05, 02–03, 02–06, 03–04, 05–06. The square 01–02–06–05 is a cycle of four
   weight-1 edges, so exactly one of them must be rejected. Accepted: 5 edges, which span {01,...,06}. → **4 choices** (drop any one square edge).
2. **Weight 2**: the remaining components are C = {01..06}, {07}, {08}, {09}. Weight-2 edges that cross components:
   01–07 (only way to reach syn07 at weight 2), 01–08 and 05–08 (two ways to reach syn08). Every other
   weight-2 edge lies inside C. Accept 01–07 and one of {01–08, 05–08}. → **2 choices**.
3. Weights 3, 4, 5 add no crossing edges ({09} is still at distance ≥ 6).
4. **Weight 6**: 01–09 is the unique cheapest edge to syn09. Accept it.

MST weight = 5·1 + 2·2 + 6 = **15**, with 8 edges (= 9 − 1).
Number of distinct MSTs = 4 × 2 = **8**. Prim and Kruskal must both return weight 15. They may legitimately
return *different* trees, so tests must compare weights and tree validity, not edge sets.

Counting link (optional): the number of labelled trees on 9 vertices is 9^(9−2) = 9^7 = **4,782,969** (Cayley).
The number of rooted binary tree topologies with 9 labelled leaves is (2·9−3)!! = 15!! = 15·13·11·9·7·5·3·1 = **2,027,025**.

UPGMA expected values are **not** derived yet. Dataset A has many tied minimum distances (six pairs at d = 1),
so the tree depends on the tie-breaking rule. That rule has to be fixed in the spec first.

## 6. Unit III: transformations

Complement A↔T, C↔G; RC(S) = reverse(complement(S)).
RC(syn01): complement of GATTACGTAATC = CTAATGCATTAG; reversed = GATTACGTAATC = syn01. So syn01 is an RC-palindrome.
RC(syn03): complement of GCTTGCGTAATC = CGAACGCATTAG; reversed = GATTACGCAAGC = syn08.
No other sequence in the set is RC-fixed.

Because R is RC-fixed, RC acts on mutations relative to R by (i, b) ↦ (13 − i, comp(b)):
{(2,C),(5,G)} ↦ {(11,G),(8,C)}, which is exactly M*(R,syn08).

Orbits of {I, RC} on the dataset: {syn01} (fixed point, size 1), {syn03, syn10, syn08} (syn10 is the same string as syn03;
as haplotypes the orbit is {syn03, syn08}). Every other haplotype is alone *within the dataset*, because its RC partner
is not present (e.g. RC(syn02) = GATTACGTAAGC).

Whole-space counts (Burnside); n = 12 is the dataset length:

| n | 4^n | RC-fixed | {I,RC} orbits = (4^n + fixed)/2 | Klein {I,R,C,RC} orbits = (4^n + 4^⌈n/2⌉ + 0 + RCfixed)/4 |
|---|---|---|---|---|
| 1 | 4 | 0 | 2 | 2 |
| 2 | 16 | 4 | **10** | 6 |
| 3 | 64 | 0 | 32 | 20 |
| 4 | 256 | 16 | 136 | 72 |
| 12 | 16,777,216 | 4^6 = 4,096 | 8,390,656 | 4,196,352 |

Fixed-point counts per element: I fixes all 4^n strings; R (plain reversal) fixes the 4^⌈n/2⌉ ordinary palindromes;
C fixes nothing because no base is its own complement; RC fixes 4^(n/2) strings for even n (choose the first half freely,
and the second half is forced) and 0 for odd n (the middle base would have to equal its own complement).

## 7. Unit IV: encodings

Weights 4^(12−i): i = 1..12 → 4194304, 1048576, 262144, 65536, 16384, 4096, 1024, 256, 64, 16, 4, 1.

enc(syn01) by Horner on digits 2 0 3 3 0 1 2 3 0 0 3 1:
0 → 2 → 8 → 35 → 143 → 572 → 2289 → 9158 → 36635 → 146540 → 586160 → 2344643 → **9,378,573**.

Every other encoding follows from enc(S) = enc(R) + Σ_{i∈M(R,S)} (s_i − r_i)·4^(12−i):

| id | correction to enc(R) | enc(S) | enc(comp(S)) = 4^12 − 1 − enc(S) | canonical = min(enc(S), enc(RC(S))) |
|---|---|---|---|---|
| syn01 | 0 | 9,378,573 | 7,398,642 | 9,378,573 (RC-fixed) |
| syn02 | +1·4^10 = +1,048,576 | 10,427,149 | 6,350,066 | 9,378,569 |
| syn03 | +1,048,576 + 2·4^7 (+32,768) | 10,459,917 | 6,317,298 | 9,378,057 |
| syn04 | syn03 + 3·4^3 (+192) | 10,460,109 | 6,317,106 | 9,181,449 |
| syn05 | −3·4 = −12 | 9,378,561 | 7,398,654 | 9,378,561 |
| syn06 | +1,048,576 − 12 | 10,427,137 | 6,350,078 | 10,427,137 |
| syn07 | +1·4^7 − 2·4^5 = +16,384 − 2,048 | 9,392,909 | 7,384,306 | 9,386,509 |
| syn08 | −2·4^4 − 1·4 = −512 − 4 | 9,378,057 | 7,399,158 | 9,378,057 |
| syn09 | −4^11 − 4^9 − 2·4^8 − 4^6 + 2·4^2 − 1 = −4,591,585 | 4,786,988 | 11,990,227 | 4,786,988 |

Complement identity: complement maps digit d ↦ 3 − d at every position, so
enc(comp(S)) = Σ (3 − s_i)·4^(n−i) = 3·(4^n − 1)/3 − enc(S) = 4^n − 1 − enc(S). Here 4^12 − 1 = 16,777,215.
Canonical k-mer check: syn03 and syn08 are an RC pair and share canonical value 9,378,057 = enc(syn08).
With A<C<G<T, numeric order of equal-length encodings is the same as lexicographic order, so min(enc) is the
lexicographic minimum used by Jellyfish's canonical representation.
enc(RC(S)) by hand: since R is RC-fixed, RC(S) = R with mutations (i, b) ↦ (13 − i, comp(b)). Apply the same correction rule:
- syn02: (2,C) ↦ (11,G): T→G at 11: −1·4 ⇒ 9,378,569; min(10,427,149, 9,378,569) = 9,378,569.
- syn04: ↦ (11,G),(8,C),(4,A): −4 − 2·256 − 3·65,536 = −197,124 ⇒ 9,181,449 (the minimum).
- syn05: (11,A) ↦ (2,T): +3·4^10 = +3,145,728 ⇒ 12,524,301; the minimum is enc(syn05) itself.
- syn06: ↦ (11,G),(2,T): −4 + 3,145,728 ⇒ 12,524,297; the minimum is enc(syn06) = 10,427,137.
- syn07: (5,C),(7,A) ↦ (8,G),(6,T): −256 + 2·4,096 = +7,936 ⇒ 9,386,509 (the minimum).
- syn09: ↦ (12,G),(10,C),(9,G),(7,T),(3,C),(1,T): +1 + 16 + 128 + 1,024 − 524,288 + 4,194,304 = +3,671,185 ⇒ 13,049,758; the minimum is enc(syn09).

### Modular buckets b(x) = x mod m

| m | syn01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | occupied buckets | colliding groups |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 = 4^2 | 13 | 13 | 13 | 13 | 1 | 1 | 13 | 9 | 12 | 4 | {01,02,03,04,07}, {05,06} |
| 7 | 1 | 5 | 6 | 2 | 3 | 0 | 1 | 3 | 3 | 6 | {01,07}, {05,08,09} |
| 13 | 9 | 5 | 0 | 10 | 10 | 6 | 6 | 0 | 11 | 6 | {03,08}, {04,05}, {06,07} |
| 17 | 13 | 12 | 4 | 9 | 1 | 0 | 1 | 7 | 9 | 7 | {05,07}, {04,09} |

Why m = 16 is a bad hash: x mod 4^j is the value of the last j base-4 digits, so it depends only on the last j bases.
Every sequence ending in "TC" (digits 3,1 → 3·4 + 1 = 13) lands in bucket 13. Check: 9,378,573 = 16·586,160 + 13.
Bucket-sharing ≠ similarity: under m = 16, syn07 shares a bucket with syn03 although d_H = 3, while syn05 and syn01
(d_H = 1) are in different buckets. Under m = 7, syn09 (the outlier) shares a bucket with syn05 and syn08.
Pigeonhole: 9 distinct haplotypes into m = 7 buckets ⇒ some bucket holds ≥ ⌈9/7⌉ = 2. This is confirmed above.

Single-substitution check (theorem in `math_audit.md` Unit IV: odd m > 3 separates every d_H = 1 pair). Read from the table:

| d = 1 pair | m = 16 | m = 7 | m = 13 | m = 17 |
|---|---|---|---|---|
| 01–02 | 13, 13 **collide** | 1, 5 | 9, 5 | 13, 12 |
| 01–05 | 13, 1 | 1, 3 | 9, 10 | 13, 1 |
| 02–03 | 13, 13 **collide** | 5, 6 | 5, 0 | 12, 4 |
| 02–06 | 13, 1 | 5, 0 | 5, 6 | 12, 0 |
| 03–04 | 13, 13 **collide** | 6, 2 | 0, 10 | 4, 9 |
| 05–06 | 1, 1 **collide** | 3, 0 | 10, 6 | 1, 0 |

All odd moduli separate all six pairs. m = 16 collides on four of them, because those substitutions are not in the last two positions.

### Number theory on a 3-mer (CRT)
The first 3-mer of R is GAT: enc = 2·16 + 0·4 + 3 = **35** (spec example: ACG = 0·16 + 1·4 + 2 = 6).
All 3-mers lie in [0, 63], so pairwise-coprime moduli with product ≥ 64 reconstruct uniquely. Take 7 and 11 (product 77):
35 ≡ 0 (mod 7), 35 ≡ 2 (mod 11).
M = 77; M1 = 11, 11 ≡ 4 (mod 7), 4·2 = 8 ≡ 1 ⇒ 11^(−1) ≡ 2 (mod 7); M2 = 7, 7·8 = 56 ≡ 1 (mod 11) ⇒ 7^(−1) ≡ 8 (mod 11).
x ≡ 0·11·2 + 2·7·8 = 112 ≡ 112 − 77 = **35** (mod 77). ✓

## 8. Unit II: counts for n = 12

- Sequence space: 4^12 = **16,777,216**.
- M(12,k) = C(12,k)·3^k: k=0: 1; k=1: 12·3 = **36**; k=2: 66·9 = **594**; k=3: 220·27 = **5,940**.
- Hamming ball of radius 2 around any sequence: 1 + 36 + 594 = **631**.
- Recurrence check, M(12,2) = M(11,2) + 3·M(11,1) = 55·9 + 3·33 = 495 + 99 = 594. ✓
- Strings that use all four bases (inclusion–exclusion over "base X missing"):
  4^n − C(4,1)3^n + C(4,2)2^n − C(4,3)1^n + C(4,4)0^n.
  n = 4: 256 − 324 + 96 − 4 + 0 = **24** (= 4!, as it must be, since each base appears exactly once).
  n = 12: 16,777,216 − 2,125,764 + 24,576 − 4 = **14,676,024**.

## How these were checked
`scripts/check_synthetic_expected.py` shares no code with `genomath` and uses only the standard library. It brute-forces: the full distance matrix, components for k = 0..9, the MST weight and the number of minimum spanning
trees (exhaustive search over all 8-edge subsets of the 36 edges), all simple paths in G_2 for the Dijkstra rows,
orbit counts for n ≤ 4 by enumerating every string, inclusion–exclusion for n = 4 by enumeration, and all encodings
and residues. Every value in this document matched. Each unit's stage (see `final_spec.md` §7) re-checks the relevant
numbers in pytest, using NetworkX/SymPy as independent oracles.
