# Unit IV: hand-derived number-theory results

Status: **expected values and proofs only**, worked by hand; nothing here was produced by `genomath`.
Dataset A encodings, complement codes, canonical codes and bucket tables are in `synthetic_expected.md` §7 and are not
repeated. Encoding: A=0, C=1, G=2, T=3, with the most significant digit at position 1.

## 1. Encoding is a bijection Σ^n → {0, …, 4^n − 1}
enc(S) = Σ_{i=1..n} s_i·4^(n−i) lies between 0 and 3·(4^n − 1)/3 = 4^n − 1. Every integer in that range has exactly one
base-4 representation with n digits (leading zeros allowed), so enc is injective and surjective.
**Decoding is the division algorithm applied n times:** x = 4q + r with 0 ≤ r < 4. Then r is the last digit, and the
process repeats on q.
Example: 35 = 4·8 + 3 (T), 8 = 4·2 + 0 (A), 2 = 4·0 + 2 (G), giving GAT. ✓
Complement: comp maps digit d ↦ 3 − d, so enc(comp S) = Σ(3 − s_i)4^(n−i) = (4^n − 1) − enc(S).
RC on codes: take digits from the least significant end, complement each, and push it onto a new number:
rc = 0; repeat n times {d = x mod 4; x = x div 4; rc = 4·rc + (3 − d)}.

**Sliding window:** the code of the next k-mer is x' = (4x mod 4^k) + d_new.
syn01 3-mers (hand): GAT 35, ATT 15, TTA 60, TAC 49, ACG 6, CGT 27, GTA 44, TAA 48, AAT 3, ATC 13.
Rolling check: from GAT, (4·35 mod 64) + 3 = 140 − 128 + 3 = 15 = ATT ✓. From ATT, (60 mod 64) + 0 = 60 = TTA ✓.
From AAT, (12 mod 64) + 1 = 13 = ATC ✓.

## 2. Canonical k-mers in syn01 (hand)
syn01 is RC-fixed, so its k-mer multiset is closed under RC.
- k = 3: 10 windows, all distinct. They pair up as GAT/ATC, ATT/AAT, TTA/TAA, TAC/GTA and ACG/CGT, giving
  **5 canonical classes and 0 palindromes** (odd k has no RC-palindromic k-mers).
- k = 4: 9 windows, all distinct: GATT, ATTA, TTAC, TACG, ACGT, CGTA, GTAA, TAAT, AATC. Pairs: GATT/AATC, ATTA/TAAT,
  TTAC/GTAA, TACG/CGTA, and ACGT is its own reverse complement. That gives **5 classes and 1 palindrome** (ACGT).

Whole-space numbers of canonical k-mers (Burnside, `unit3_expected.md` §4): k = 4: 136, k = 6: 2080, k = 8: 32896.
RC-palindromic k-mers: 4^(k/2) = 16, 64 and 256.
**Identity used as a check:** #distinct canonical = #distinct k-mers − #{unordered pairs {w, RC(w)}, w ≠ RC(w), both present}.

## 3. When can a single substitution keep the bucket? (exact criterion, proved here)
A substitution at position i changes the base digit by δ ∈ {±1, ±2, ±3}, so it changes x by δ·4^t, with t = n − i ∈ {0, …, n−1}.
The bucket is unchanged iff m | δ·4^t. Write m = 2^a·q with q odd.
- Since gcd(q, 4^t) = 1, q | δ·4^t iff q | δ. With |δ| ≤ 3, this needs q ∈ {1, 3}.
- The power of 2 in δ·4^t is 2t + v, where v = 1 for δ = ±2 and v = 0 for δ = ±1, ±3. We need 2t + v ≥ a.

**Theorem.** If the odd part q of m is ≥ 5, no single substitution ever keeps the bucket, at any length.
Otherwise a bucket-preserving single substitution exists exactly when the sequence has a position with t ≥ t_min:
- q = 3: δ must be ±3 (A↔T), so v = 0 and t_min = ⌈a/2⌉.
- q = 1: δ = ±2 (A↔G or C↔T) gives v = 1, so t_min = ⌊a/2⌋.

Checks: m = 16 = 2^4: t_min = 2. m = 64 = 2^6: 3. m = 8 = 2^3: 1 (δ = 2 at t = 1 gives 8). m = 12 = 2^2·3: 1 (3·4 = 12).
m = 6 = 2·3: 1. m = 3 and m = 1: 0.

Corollaries:
- Every odd prime p ≥ 5, and odd m ≥ 5 generally, separates all single substitutions. Even moduli with odd part ≥ 5,
  such as m = 10, also separate them; the audit's condition "gcd(m,4) = 1" is sufficient but not necessary.
- **Bad hash:** if m | 4^j (m a power of 2, m ≤ 4^j), then 4^t ≡ 0 (mod m) for t ≥ j. So b(x) depends only on the last
  j bases: m = 16 sees exactly the last 2 bases, and m = 64 the last 3.
  For m = 16 a substitution collides **iff** t ≥ 2, because at t = 1, 16 | 4δ would need 4 | δ.
  For m = 64 a substitution collides iff t ≥ 3, because at t = 2, 64 | 16δ would need 4 | δ.

Smallest j with m | 4^j: m = 16 → **2**; m = 64 → **3**; odd m → **none**.

## 4. Two substitutions can collide, and Fermat/Euler say where
For odd m, 4^t (mod m) is periodic with period ord_m(4), which divides φ(m) (Euler) and divides p − 1 for prime p (Fermat).
Orders by hand:
- mod 7: 4, 16 ≡ 2, 8 ≡ 1, so **ord = 3** (divides 6).
- mod 13: 4, 3, 12, 9, 10, 1, so **ord = 6** (divides 12).
- mod 17: 4, 16 ≡ −1, so 4^4 ≡ 1 and **ord = 4** (divides 16).
- mod 31: 4^5 = 1024 = 33·31 + 1, and 5 is prime, so **ord = 5** (divides 30).
- mod 61: 2^10 = 1024 ≡ 48, 2^12 ≡ 9, 2^20 ≡ 48² = 2304 ≡ 47, 2^30 ≡ 47·48 = 2256 ≡ 60 ≡ −1. The order of 2 divides 60
  but not 30, and 2^4, 2^12, 2^20 ≢ 1, so ord_61(2) = 60. Then ord_61(4) = 60/gcd(60,2) = **30** (divides 60).

**Cancellation example.** Take syn01 and apply (6: C→A, δ = −1, t = 6) and (12: C→G, δ = +1, t = 0).
The total change is 1 − 4^6 = −4095 = −7·585, so the new string S has d_H(syn01, S) = 2 but the **same bucket mod 7**.
In general the change is ≡ 0 (mod m) iff 4^6 ≡ 1 iff ord_m(4) | 6. Among {16, 64, 7, 13, 17, 31, 61}:
- 7 (ord 3) and 13 (ord 6) collide.
- 17 (ord 4), 31 (5) and 61 (30) do not.
- For 16 and 64, −4095 ≡ 1, so no collision (position 12 is inside the last 2 or 3 bases).

Expected collide set: **{7, 13}**. The modified string is GATTAAGTAATG.

φ values: φ(7) = 6, φ(13) = 12, φ(17) = 16, φ(31) = 30, φ(61) = 60, φ(77) = 6·10 = 60, φ(16) = 8, φ(64) = 32.
Euler: 4^φ(m) ≡ 1 (mod m) whenever gcd(4, m) = 1. For m = 16 it fails: 4^8 ≡ 0 (mod 16), since gcd(4, 16) = 4.

## 5. Euclid, inverses and congruences (hand)
- gcd(252, 198): 252 = 1·198 + 54; 198 = 3·54 + 36; 54 = 1·36 + 18; 36 = 2·18 + 0, so the gcd is **18**.
  Back-substitution: 18 = 54 − 36 = 54 − (198 − 3·54) = 4·54 − 198 = 4·(252 − 198) − 198 = **4·252 − 5·198**.
  Check: 1008 − 990 = 18. ✓
- gcd(enc(syn01), 4^12) = gcd(9,378,573, 16,777,216) = **1**, because the code is odd and 4^12 is a power of 2.
- Inverses: 11^(−1) mod 7 = **2** (11·2 = 22 = 3·7 + 1). 7^(−1) mod 11 = **8** (56 = 5·11 + 1). 4^(−1) mod 7 = **2**.
  4^(−1) mod 16 does **not exist** (gcd = 4).
- Congruences ax ≡ b (mod m), solvable iff g = gcd(a,m) divides b, with g solutions mod m:
  3x ≡ 6 (mod 9): g = 3 | 6, so x ≡ 2 (mod 3), giving solutions **{2, 5, 8}**.
  2x ≡ 1 (mod 4): g = 2 ∤ 1, so there is **no solution**.

## 6. Chinese Remainder Theorem on 3-mer codes (range 0..63)
- GAT = 35 has residues (0 mod 7, 2 mod 11) and reconstructs to **35** (worked in `synthetic_expected.md` §7).
- Moduli (7, 11), product 77 ≥ 64, and (5, 13), product 65 ≥ 64: all **64/64** codes are recovered.
- Moduli (3, 5), product 15 < 64: CRT returns x mod 15, which is correct only for codes 0..14, so **15/64** are recovered.
  For example, codes 0 and 15 (AAA and ATT = 0·16 + 3·4 + 3) share residues (0, 0).
- Moduli (4, 6) are not coprime (gcd 2), so the reconstruction must be **refused**.

## 7. Pigeonhole bounds on buckets
N keys in m buckets force some bucket to hold at least ⌈N/m⌉ keys.
- Dataset A (9 haplotypes), m = 7: bound 2; the observed maximum is 3 ({05, 08, 09}, `synthetic_expected.md` §7).
- Dataset B (30 haplotypes):

| m | 7 | 13 | 16 | 17 | 31 | 61 | 64 |
|---|---|---|---|---|---|---|---|
| ⌈30/m⌉ | 5 | 3 | 2 | 2 | 1 | 1 | 1 |
