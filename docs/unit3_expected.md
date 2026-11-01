# Unit III: hand-derived algebraic results

Status: **expected values and proofs only**, worked by hand; nothing here was produced by `genomath`.
Maps act on strings over Σ = {A,C,G,T}: I (identity), R (reverse), C (complement A↔T, C↔G letterwise), RC = R∘C.
**Composition convention:** (a∘b)(x) = a(b(x)): apply b first, then a. A Cayley-table entry in row a, column b is a∘b.

## 1. Why {I, R, C, RC} is a group (proof for every length)
- R∘R = I: reversing twice restores the order. C∘C = I: complementing a letter twice restores it.
- R∘C = C∘R: C changes letters without moving them and R moves letters without changing them, so the order does not matter.
  Write RC for this map.
- RC∘RC = R∘C∘R∘C = R∘R∘C∘C = I.
- Products of distinct non-identity elements: R∘C = RC, R∘RC = R∘R∘C = C, C∘RC = C∘C∘R = R.

So the set is closed, I is the identity, every element is its own inverse, and composition of functions is associative
(for any maps, ((a∘b)∘c)(x) = a(b(c(x))) = (a∘(b∘c))(x)). It is a group of order 4 in which every non-identity element
has order 2. Hence it is **not cyclic** (Z4 has an element of order 4), and it is the Klein four-group, isomorphic to
Z2 × Z2 via I↦(0,0), R↦(1,0), C↦(0,1), RC↦(1,1). {I, RC} is a subgroup of order 2, isomorphic to Z2.

### Cayley tables (row ∘ column)
```
 Z2 | I   RC          V4 | I   R   C   RC
 ---+--------         ---+----------------
 I  | I   RC          I  | I   R   C   RC
 RC | RC  I           R  | R   I   RC  C
                      C  | C   RC  I   R
                      RC | RC  C   R   I
```
Both tables are symmetric (abelian groups). Every row and column is a permutation of the elements (Latin square), and
the diagonal is all I.

### Faithfulness depends on length [subtle point]
As maps on Σ^n for a *single* n, the four elements need not be distinct:
- n = 0: all four maps send "" to "" (the action is trivial).
- n = 1: R = I (a one-letter word reversed is itself) and RC = C.
- n ≥ 2: all four are distinct. Witnesses: R(AC) = CA ≠ AC; C(AC) = TG; RC(AA) = TT ≠ AA; R(AC) = CA ≠ TG = C(AC);
  R(AA) = AA ≠ TT = RC(AA); C(AC) = TG ≠ GT = RC(AC).

Function equality is therefore decided on **all words of length 0..6 together**, where the four maps are distinct.
Burnside's lemma does not require a faithful action, so it still applies at n = 0 and n = 1.

## 2. Subgroups and cosets
Subgroups of V4 (a subset containing I, closed, with inverses): {I}, {I,R}, {I,C}, {I,RC}, V4, so **5** in all.
Left cosets of H = {I,RC}: I·H = {I,RC}, R·H = {R, R∘RC} = {R,C}, C·H = {C, R}, RC·H = {RC, I}.
There are two distinct cosets, {I,RC} and {R,C}, and index 2 = |V4|/|H| = 4/2 (Lagrange).

## 3. Fixed points |Fix(g)| on Σ^n
- I: all 4^n strings.
- R: s_i = s_{n+1−i}; the first ⌈n/2⌉ letters are free and the rest are forced, so **4^⌈n/2⌉**.
- C: needs s_i = comp(s_i), which is impossible because no base is its own complement. So **0** for n ≥ 1, and 1 for n = 0 (the empty word).
- RC: s_{n+1−i} = comp(s_i). For even n the first n/2 letters are free: **4^(n/2)**. For odd n the middle letter would
  have to equal its own complement: **0**. For n = 0: 1.

| n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| Fix(I) | 1 | 4 | 16 | 64 | 256 | 1024 | 4096 | 16384 | 65536 | 16,777,216 |
| Fix(R) | 1 | 4 | 4 | 16 | 16 | 64 | 64 | 256 | 256 | 4096 |
| Fix(C) | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Fix(RC) | 1 | 0 | 4 | 0 | 16 | 0 | 64 | 0 | 256 | 4096 |

## 4. Burnside orbit counts: #orbits = (1/|G|)·Σ_g |Fix(g)|
| n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| Z2 = {I,RC}: (Fix I + Fix RC)/2 | 1 | 2 | **10** | 32 | 136 | 512 | 2080 | 8192 | 32896 | 8,390,656 |
| V4: (Fix I + Fix R + Fix C + Fix RC)/4 | 1 | 2 | **6** | 20 | 72 | 272 | 1056 | 4160 | 16512 | 4,196,352 |

Worked examples: Z2, n = 8: (65536 + 256)/2 = 32896. V4, n = 7: (16384 + 256 + 0 + 0)/4 = 16640/4 = 4160.
V4, n = 6: (4096 + 64 + 0 + 64)/4 = 4224/4 = 1056.

### Orbit sizes (orbit–stabiliser: |orbit| · |stabiliser| = |G|)
- Z2: fixed strings form orbits of size 1, so there are Fix(RC) of them; all other strings pair up, giving (4^n − Fix RC)/2 orbits of size 2.
  n = 4: 16 of size 1 and 120 of size 2, total 136.
- V4, n ≥ 1: no string is fixed by C, so no stabiliser contains C. A stabiliser is therefore {I}, {I,R} or {I,RC}.
  {I,R} and {I,RC} cannot both lie in one stabiliser, because then C = R∘RC would fix the string.
  - Size-2 orbits: (Fix R + Fix RC)/2.
  - Size-4 orbits: (4^n − Fix R − Fix RC)/4.

| n | V4 size-2 orbits | V4 size-4 orbits | total |
|---|---|---|---|
| 1 | (4+0)/2 = 2 | 0 | 2 |
| 2 | (4+4)/2 = 4 | (16−8)/4 = 2 | 6 |
| 3 | (16+0)/2 = 8 | (64−16)/4 = 12 | 20 |
| 4 | (16+16)/2 = 16 | (256−32)/4 = 56 | 72 |

n = 2 written out: {AA,TT}, {CC,GG}, {AT,TA}, {CG,GC} (size 2) and {AC,CA,TG,GT}, {AG,GA,TC,CT} (size 4).

## 5. Permutations of the alphabet
The complement as a permutation of Σ is **(A T)(C G)**: two disjoint transpositions. Its order is lcm(2,2) = 2, and it is
even because it is a product of 2 transpositions.
Reversal as a permutation of positions 1..n is i ↦ n+1−i. For n = 12 that is (1 12)(2 11)(3 10)(4 9)(5 8)(6 7):
6 transpositions, order 2, even.
(Reversal is the *only* position permutation used. No arbitrary position permutations, by decision.)

## 6. Point substitutions form a monoid, not a group
Define σ_{i,b}: Σ^n → Σ^n by "set position i to b". A composition of such maps is described by a pattern
p ∈ {·, A, C, G, T}^n, where f_p keeps position i if p_i = · and otherwise writes p_i.
- Composition rule: (f_p ∘ f_q) has pattern r with r_i = p_i if p_i ≠ ·, else q_i (the later write wins).
  So the set of all f_p is closed and contains the identity (p = ··…·). It is associative (function composition).
  It is a **monoid** with **5^n** elements: 5, 25 and 125 for n = 1, 2, 3.
- Every f_p is idempotent (writing the same values twice changes nothing), so all 5^n elements are idempotents.
- Only the identity has an inverse. Any f_p with some p_i = b ≠ · is not injective: two strings that differ only at
  position i (e.g. A and C there) have the same image. A non-injective map has no inverse. So the group of units is **{I}**.
- It is not commutative: on n = 1, f_A ∘ f_C = f_A but f_C ∘ f_A = f_C.

**Conclusion:** "mutations" modelled as substitution operators do not form a group. The {I,RC} group is legitimate
because RC is a bijection (an involution) describing strand orientation, not a mutation.
(One could force a group by defining base "shifts" A→C→G→T→A mod 4 at each position, giving Z4^n. That structure is
an artefact of the arbitrary A=0..T=3 labelling and has no biochemical meaning, so it is not used.)

## 7. Concatenation laws (Σ* as a free monoid under concatenation, identity "")
C(xy) = C(x)C(y) (C is a monoid homomorphism). R(xy) = R(y)R(x) and RC(xy) = RC(y)RC(x) (anti-homomorphisms).
Example: RC(GA·TT) = RC(GATT) = AATC = RC(TT)·RC(GA) = AA·TC. ✓

## 8. Datasets
**Dataset A** (haplotypes; RC(syn01) = syn01 and RC(syn03) = syn08 from `synthetic_expected.md` §6):
- RC-fixed: syn01 only. R-fixed: none, because position 1 is G and position 12 is C in every haplotype except syn09,
  which has C and A, and G ≠ C, C ≠ A.
- Z2-orbits within the haplotypes: {syn03, syn08} plus 7 singletons, so **8 classes from 9 haplotypes**.
- V4 orbit of syn01: RC(syn01) = syn01 means R(syn01) = C(syn01) = CTAATGCATTAG. The orbit is {GATTACGTAATC, CTAATGCATTAG},
  of size **2**, with stabiliser {I,RC} (2·2 = 4). The V4 orbits of syn03 and of syn09 have size **4** (stabiliser {I}).
- Canonical representative (lexicographically smallest member of the Z2-orbit): syn03 and syn08 both map to
  GATTACGCAAGC = syn08 ("GA…" < "GC…").

## 9. A second Klein four-group, acting on bases (decision D9)
This is a **different group acting on a different set** from §1: §1's {I, R, C, RC} acts on strings Σ^n, while this one acts
on the four bases Σ. Both are abstractly Klein, and they are linked: the string map C is the letterwise extension of the
base permutation (A T)(C G).

**The three fixed-point-free involutions of Σ.** A fixed-point-free involution of a 4-element set is a pairing of all four
elements, and the number of perfect matchings on 4 labelled points is 3!! = 3·1 = **3**:
σ = (A T)(C G), τ = (A G)(C T), υ = (A C)(G T). With the identity they form K = {id, σ, τ, υ}.

**Closure** (apply the right factor first): σ∘τ sends A → τ(A) = G → σ(G) = C, so A → C; C → T → A; G → A → T; T → C → G.
That is (A C)(G T) = υ. By symmetry τ∘υ = σ and σ∘υ = τ, and each element is its own inverse. So K is closed, has the
identity, has inverses, and composition is associative: a Klein four-group, isomorphic to Z2 × Z2 and not to Z4
(every non-identity element has order 2).

**Cayley table** (row ∘ column):
```
     | id  σ   τ   υ
 ----+----------------
 id  | id  σ   τ   υ
 σ   | σ   id  υ   τ
 τ   | τ   υ   id  σ
 υ   | υ   τ   σ   id
```

**Every substitution comes from exactly one non-identity element.** Each of σ, τ, υ moves all four bases, so each
realises exactly 4 ordered pairs (x, y) with y ≠ x; 3 × 4 = **12**, which is every ordered pair. Hence the map
"(x, y) ↦ the element taking x to y" is a bijection from the 12 ordered pairs onto 3 elements × 4 bases.
Because each element is an **involution**, the same element also takes y to x, so the classification is well defined on
**unordered** pairs {x, y}. This matters here: two observed sequences differ at a position without telling us which base
came first, so a direction-free classification is the only honest one.

**Transitions and transversions.** τ = (A G)(C T) pairs purine with purine (A, G) and pyrimidine with pyrimidine (C, T);
these are the **transitions**. σ and υ each pair a purine with a pyrimidine; these are the **transversions**.
Of the 6 unordered pairs: **2 transitions** (A↔G, C↔T) and **4 transversions** (A↔T, C↔G from σ; A↔C, G↔T from υ).

**Combinatorial baseline.** If all 6 unordered pairs were equally likely, the expected ratio would be
2 : 4 = **1 : 2**, i.e. transitions would be 1/3 of substitutions. This is a counting statement about Σ only.
I have **not** located a source for what ratio real genomic data typically shows, so the experiment reports observed
counts against this 1:2 combinatorial baseline and makes no claim about what is biologically expected.

**Action on Σ.** K acts transitively: from A the four elements reach A, T, G, C respectively, so there is a single orbit
of size 4 with trivial stabiliser {id}. Orbit–stabiliser: 4 × 1 = 4 = |K|. The action is regular (simply transitive).
Burnside confirms the orbit count: (|Fix id| + |Fix σ| + |Fix τ| + |Fix υ|)/4 = (4 + 0 + 0 + 0)/4 = **1**.

**Dataset B** (34 clean CDS, each starting with ATG and ending in a stop codon; QC verified):
RC(s) starts with the reverse complement of s's stop codon (TAA→TTA, TAG→CTA, TGA→TCA), so it cannot start with ATG.
Therefore RC(s) is never in Dataset B, and no sequence is RC-fixed. Expected: **0 RC pairs, 0 RC-fixed**, and the
Z2-classes are exactly the 30 haplotypes. Unit III does no work on full annotated CDS (the weak link in `math_audit.md`).
