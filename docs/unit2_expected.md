# Unit II: hand-derived counts and derivations

Status: **expected values and proofs only**, worked by hand. Nothing here was produced by `genomath`. Where a number
could be checked by an independent identity (row sums = 4^n, Stirling numbers), that check is shown.

## 1. Sequence space (product rule)
A sequence in Σ^n is n independent choices from 4 bases, so N(n) = 4·4·…·4 = 4^n.
1, 4, 16, 64, 256, 1024, 4096, 16384, 65536 for n = 0..8. Codons: 4^3 = 64. Dataset A (n = 12): 16,777,216.
Dataset B (n = 1701): 4^1701 has ⌊1701·log10 4⌋ + 1 = ⌊1024.10…⌋ + 1 = **1025** decimal digits (log10 4 = 0.60206).

## 2. M(n,k): sequences at Hamming distance exactly k from a fixed reference
**Derivation.** Choose which k of the n positions change: C(n,k) ways. Each chosen position takes one of the 3 bases
different from the reference: 3^k ways (product rule). The choices are independent, so M(n,k) = C(n,k)·3^k.
Only substitutions are counted; insertions and deletions change the length and are outside Σ^n.

**Recurrence.** Split on position n:
- position n unchanged: the first n − 1 positions carry all k changes → M(n−1,k);
- position n changed (3 choices): the first n − 1 positions carry k − 1 changes → 3·M(n−1,k−1).

The two cases are disjoint and exhaustive, so M(n,k) = M(n−1,k) + 3·M(n−1,k−1). Base case: M(0,0) = 1 (the empty string
is the one string at distance 0 from the empty reference). M(n,k) = 0 when k < 0 or k > n.

**Induction proof that the closed form solves the recurrence.** Let f(n,k) = C(n,k)·3^k, with C(n,k) = 0 outside 0 ≤ k ≤ n.
- Base: f(0,0) = 1, and f(0,k) = 0 for k ≠ 0, matching the base cases.
- Step: assume f(n−1,·) = M(n−1,·). Then M(n,k) = f(n−1,k) + 3f(n−1,k−1) = C(n−1,k)3^k + 3·C(n−1,k−1)3^(k−1)
  = 3^k·(C(n−1,k) + C(n−1,k−1)) = 3^k·C(n,k) = f(n,k), by Pascal's rule.

By induction on n, M(n,k) = C(n,k)·3^k for all n ≥ 0 and all k.

**Row sums.** Σ_k C(n,k)·3^k·1^(n−k) = (3 + 1)^n = 4^n (binomial theorem). Every string is at exactly one distance
from the reference, which is the counting reason the identity must hold.

### Hand table M(n,k), n = 0..8 (built row by row with the recurrence; each row sum checked against 4^n)
| n \ k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | row sum | 4^n |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | | | | | | | | | 1 | 1 |
| 1 | 1 | 3 | | | | | | | | 4 | 4 |
| 2 | 1 | 6 | 9 | | | | | | | 16 | 16 |
| 3 | 1 | 9 | 27 | 27 | | | | | | 64 | 64 |
| 4 | 1 | 12 | 54 | 108 | 81 | | | | | 256 | 256 |
| 5 | 1 | 15 | 90 | 270 | 405 | 243 | | | | 1024 | 1024 |
| 6 | 1 | 18 | 135 | 540 | 1215 | 1458 | 729 | | | 4096 | 4096 |
| 7 | 1 | 21 | 189 | 945 | 2835 | 5103 | 5103 | 2187 | | 16384 | 16384 |
| 8 | 1 | 24 | 252 | 1512 | 5670 | 13608 | 20412 | 17496 | 6561 | 65536 | 65536 |

Sample recurrence steps: M(5,3) = M(4,3) + 3·M(4,2) = 108 + 162 = 270. M(8,6) = M(7,6) + 3·M(7,5) = 5103 + 15309 = 20412.
M(8,4) = M(7,4) + 3·M(7,3) = 2835 + 2835 = 5670.

Cells in the table: Σ_{n=0..8} (n + 1) = 45.

### Hamming balls B(n,r) = Σ_{k ≤ r} M(n,k)
B(8,2) = 1 + 24 + 252 = **277**. B(12,2) = 1 + 36 + 594 = **631**. B(1701,1) = 1 + 3·1701 = **5104**.
So a 1701-nt sequence has **5103** possible one-substitution neighbours.

### Enumeration guard (proposed cap 10^5 strings)
- M(12,6) = C(12,6)·3^6 = 924·729 = 646,800 + 26,796 = **673,596** > 10^5 → must be refused.
- M(1701,1) = 5103 ≤ 10^5 → allowed.
- M(1701,2) = C(1701,2)·9 = (1701·1700/2)·9 = 1,445,850·9 = **13,012,650** > 10^5 → must be refused.
- B(8,8) = 4^8 = 65,536 ≤ 10^5 → allowed (the whole space Σ^8 is the ball of radius 8).

## 3. Permutations: orderings of k substitutions
To turn S into T with d_H(S,T) = k using k single substitutions (a shortest path in the Hamming graph), each differing
position must be changed exactly once, directly to T's base. Changing a position twice, or changing a position where S and
T agree, would need more than k steps. A shortest path is therefore an ordering of the k differing positions:
**k!** paths. Example from Dataset A: syn01 → syn04 has M = {2,5,9}, k = 3, so 3! = 6 orderings.
Of these, only the order 2, 5, 9 passes through observed haplotypes only (syn02 = R+{2}, syn03 = R+{2,5}).
The other first steps R+{5:G} and R+{9:T}, and the intermediate R+{2,9}, are not in the dataset. Expected: **1 of 6**
orderings is fully observed.

Ordered selections: P(n,k) = n!/(n−k)! = C(n,k)·k!. For example, P(12,3) = 12·11·10 = 1320 = 220·6.

## 4. Inclusion–exclusion: strings using all four bases
Let A_x = {strings with no x}, for x ∈ Σ. Then |A_x| = 3^n, |A_x ∩ A_y| = 2^n, |A_x ∩ A_y ∩ A_z| = 1, and the
four-way intersection is empty (0^n = 0 for n ≥ 1). Subtracting |A_A ∪ A_C ∪ A_G ∪ A_T| from 4^n:
U(n) = 4^n − 4·3^n + 6·2^n − 4·1^n + 0.

| n | 4^n | 4·3^n | 6·2^n | U(n) | check: 4!·S(n,4) |
|---|---|---|---|---|---|
| 4 | 256 | 324 | 96 | **24** | 24·1 = 24 |
| 5 | 1024 | 972 | 192 | **240** | 24·10 = 240 |
| 6 | 4096 | 2916 | 384 | **1560** | 24·65 = 1560 |
| 7 | 16384 | 8748 | 768 | **8400** | 24·350 = 8400 |
| 8 | 65536 | 26244 | 1536 | **40824** | 24·1701 = 40824 |
| 12 | 16,777,216 | 2,125,764 | 24,576 | **14,676,024** | (not checked) |

U(n) also counts surjections [n] → Σ, which is why it equals 4!·S(n,4) (S = Stirling numbers of the second kind;
S(5,4) = 10, S(6,4) = 65, S(7,4) = 350 and S(8,4) = 1701 are standard table values).
For n = 1, 2, 3, U(n) = 0: fewer than four positions cannot use four bases (4 − 4·3 + 6·2 − 4 = 0 for n = 1, for example).

## 5. Pigeonhole: forced repeated k-mers
A sequence of length L has L − k + 1 windows of length k (the pigeons), and there are 4^k possible k-mers (the holes).
If L − k + 1 > 4^k, some k-mer must occur twice.
For L = 1701: k = 5 gives 1697 windows > 1024, so a repeat is **forced**. k = 6 gives 1696 ≤ 4096, so it is **not forced**.
The largest forced k is **5**. Every one of the 34 clean Dataset B sequences must therefore contain a repeated 5-mer.
The bound is tight for k = 1: "ACGT" has 4 windows and no repeated base, while any length-5 string must repeat one.
