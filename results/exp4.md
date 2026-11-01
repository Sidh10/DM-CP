# Experiment 4: Encoding and modular indexing

- run date: 2026-10-31
- command: `python -m experiments.exp4_modular_indexing`
- Python 3.13.11 on Windows

## Question
What does representing DNA words as base-4 integers make provable? In particular, how does the choice of modulus change bucket collisions, when can a single substitution leave the bucket unchanged, and how many {I,RC}-canonical k-mers does real data use?

## Input
Dataset A haplotypes (n = 12); all 3-mer codes 0..63; Dataset B haplotypes (n = 1701) and their k-mers for k = 4, 6, 8; moduli [16, 64, 7, 13, 17, 31, 61] (provisional decision 8).

## Mathematical model
enc(S) = Σ s_i 4^(n−i) (bijection; decoding by the division algorithm); enc(comp S) = 4^n − 1 − enc(S); canonical = min(enc S, enc RC S); b(x) = x mod m; pigeonhole ⌈N/m⌉; exact single-substitution criterion via m = 2^a·q; Euclid / extended Euclid / inverses / linear congruences; CRT; Euler φ and multiplicative order (pedagogical).

## Procedure
Encode every Dataset A haplotype and compare codes, complement codes, canonical codes and buckets with the hand tables. Test every one-substitution pair against the collision theorem on both datasets. Reconstruct all 64 3-mer codes by CRT under several moduli pairs. Compute orders and φ, and the two-substitution cancellation example. Bucket Dataset B full-sequence codes and distinct 8-mers. Take the {I,RC} k-mer census of Dataset B.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | A: enc(S) | {syn01: 9378573, syn02: 10427149, syn03: 10459917, syn04: 10460109, syn05: 9378561, syn06: 10427137, syn07: 9392909, syn08: 9378057, syn09: 4786988} | docs/synthetic_expected.md §7 |
| 2 | A: enc(comp S) = 4^12 − 1 − enc(S) | {syn01: 7398642, syn02: 6350066, syn03: 6317298, syn04: 6317106, syn05: 7398654, syn06: 6350078, syn07: 7384306, syn08: 7399158, syn09: 11990227} | docs/synthetic_expected.md §7 |
| 3 | A: canonical min(enc S, enc RC S) | {syn01: 9378573, syn02: 9378569, syn03: 9378057, syn04: 9181449, syn05: 9378561, syn06: 10427137, syn07: 9386509, syn08: 9378057, syn09: 4786988} | docs/synthetic_expected.md §7 |
| 4 | A: decode(enc(S)) = S | True | docs/unit4_expected.md §1 — *consistency check* |
| 5 | A: buckets for m = 16, 7, 13, 17 | {16: [13, 13, 13, 13, 1, 1, 13, 9, 12], 7: [1, 5, 6, 2, 3, 0, 1, 3, 3], 13: [9, 5, 0, 10, 10, 6, 6, 0, 11], 17: [13, 12, 4, 9, 1, 0, 1, 7, 9]} | docs/synthetic_expected.md §7 |
| 6 | A: colliding groups | {16: [[syn01, syn02, syn03, syn04, syn07], [syn05, syn06]], 7: [[syn01, syn07], [syn05, syn08, syn09]], 13: [[syn03, syn08], [syn04, syn05], [syn06, syn07]], 17: [[syn04, syn09], [syn05, syn07]]} | docs/synthetic_expected.md §7 |
| 7 | A: d_H = 1 pairs sharing a bucket | {16: 4, 7: 0, 13: 0, 17: 0} | docs/synthetic_expected.md §7 |
| 8 | A: pigeonhole, 9 keys in 7 buckets | {bound: 2, max_load: 3} | docs/unit4_expected.md §7 |
| 9 | last bases seen: smallest j with m \| 4^j | {16: 2, 64: 3, 7: None, 13: None, 17: None, 31: None, 61: None} | docs/unit4_expected.md §3 |
| 10 | t_min for a bucket-preserving single substitution | {16: 2, 64: 3, 8: 1, 12: 1, 6: 1, 3: 0, 1: 0, 7: None, 10: None} | docs/unit4_expected.md §3 |
| 11 | syn01 3-mer codes by sliding window | [35, 15, 60, 49, 6, 27, 44, 48, 3, 13] | docs/unit4_expected.md §1 |
| 12 | syn01 k-mer census | {3: {distinct: 10, canonical: 5, palindromic: 0}, 4: {distinct: 9, canonical: 5, palindromic: 1}} | docs/unit4_expected.md §2 |
| 13 | extended Euclid on (252, 198) | (18, 4, -5) | docs/unit4_expected.md §5 |
| 14 | modular inverses | {11 mod 7: 2, 7 mod 11: 8, 4 mod 7: 2, 4 mod 16: none} | docs/unit4_expected.md §5 |
| 15 | linear congruences | {3x≡6 (mod 9): [2, 5, 8], 2x≡1 (mod 4): []} | docs/unit4_expected.md §5 |
| 16 | CRT reconstruction of 3-mer codes | {GAT from (0 mod 7, 2 mod 11): 35, (7,11) round-trips: 64, (5,13) round-trips: 64, (3,5) round-trips: 15, (4,6): refused} | docs/unit4_expected.md §6 |
| 17 | ord_m(4) (pedagogical) | {7: 3, 13: 6, 17: 4, 31: 5, 61: 30} | docs/unit4_expected.md §4 |
| 18 | φ(m) | {7: 6, 13: 12, 17: 16, 31: 30, 61: 60, 77: 60, 16: 8, 64: 32} | docs/unit4_expected.md §4 |
| 19 | Euler: 4^φ(m) ≡ 1 (mod m) | {7: True, 13: True, 17: True, 31: True, 61: True, 77: True, 16: False} | docs/unit4_expected.md §4 |
| 20 | two-substitution cancellation (Fermat periodicity) | {string: GATTAAGTAATG, d_H: 2, collide: [7, 13]} | docs/unit4_expected.md §4 |
| 21 | B: pigeonhole bound ⌈30/m⌉ | {7: 5, 13: 3, 16: 2, 17: 2, 31: 1, 61: 1, 64: 1} | docs/unit4_expected.md §7 |
| 22 | B: max load ≥ bound for every m | True | docs/unit4_expected.md §7 — *consistency check* |
| 23 | B: (occupied, max load, collisions) per m | *(none: reported as observation)* |  |
| 24 | B: last 3 bases of every haplotype (count) | *(none: reported as observation)* |  |
| 25 | B: d_H = 1 pairs | 8 | docs/dataset_b_qc.md (independent script) |
| 26 | B: d_H = 1 collisions, observed = predicted by §3 theorem | {16: 8, 64: 8, 7: 0, 13: 0, 17: 0, 31: 0, 61: 0} | docs/unit4_expected.md §3 — *consistency check* |
| 27 | B: d_H = 1 collisions for odd primes (theorem: none) | {7: 0, 13: 0, 17: 0, 31: 0, 61: 0} | docs/unit4_expected.md §3 |
| 28 | B: canonical = distinct − present RC pairs (k = 4, 6, 8) | {4: True, 6: True, 8: True} | docs/unit4_expected.md §2 — *consistency check* |
| 29 | B: whole-space canonical classes (Burnside) | {4: 136, 6: 2080, 8: 32896} | docs/unit4_expected.md §2 |
| 30 | B: k-mer census (distinct, canonical, palindromic, RC pairs present) | *(none: reported as observation)* |  |
| 31 | B: distinct 8-mers bucketed, max load ≥ ⌈N/m⌉ | {64: True, 61: True} | docs/unit4_expected.md §7 — *consistency check* |
| 32 | B: 8-mer buckets, m = 64 (last 3 bases) vs m = 61 | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | A: enc(S) | {syn01: 9378573, syn02: 10427149, syn03: 10459917, syn04: 10460109, syn05: 9378561, syn06: 10427137, syn07: 9392909, syn08: 9378057, syn09: 4786988} |
| 2 | A: enc(comp S) = 4^12 − 1 − enc(S) | {syn01: 7398642, syn02: 6350066, syn03: 6317298, syn04: 6317106, syn05: 7398654, syn06: 6350078, syn07: 7384306, syn08: 7399158, syn09: 11990227} |
| 3 | A: canonical min(enc S, enc RC S) | {syn01: 9378573, syn02: 9378569, syn03: 9378057, syn04: 9181449, syn05: 9378561, syn06: 10427137, syn07: 9386509, syn08: 9378057, syn09: 4786988} |
| 4 | A: decode(enc(S)) = S | True |
| 5 | A: buckets for m = 16, 7, 13, 17 | {16: [13, 13, 13, 13, 1, 1, 13, 9, 12], 7: [1, 5, 6, 2, 3, 0, 1, 3, 3], 13: [9, 5, 0, 10, 10, 6, 6, 0, 11], 17: [13, 12, 4, 9, 1, 0, 1, 7, 9]} |
| 6 | A: colliding groups | {16: [[syn01, syn02, syn03, syn04, syn07], [syn05, syn06]], 7: [[syn01, syn07], [syn05, syn08, syn09]], 13: [[syn03, syn08], [syn04, syn05], [syn06, syn07]], 17: [[syn04, syn09], [syn05, syn07]]} |
| 7 | A: d_H = 1 pairs sharing a bucket | {16: 4, 7: 0, 13: 0, 17: 0} |
| 8 | A: pigeonhole, 9 keys in 7 buckets | {bound: 2, max_load: 3} |
| 9 | last bases seen: smallest j with m \| 4^j | {16: 2, 64: 3, 7: None, 13: None, 17: None, 31: None, 61: None} |
| 10 | t_min for a bucket-preserving single substitution | {16: 2, 64: 3, 8: 1, 12: 1, 6: 1, 3: 0, 1: 0, 7: None, 10: None} |
| 11 | syn01 3-mer codes by sliding window | [35, 15, 60, 49, 6, 27, 44, 48, 3, 13] |
| 12 | syn01 k-mer census | {3: {distinct: 10, canonical: 5, palindromic: 0}, 4: {distinct: 9, canonical: 5, palindromic: 1}} |
| 13 | extended Euclid on (252, 198) | (18, 4, -5) |
| 14 | modular inverses | {11 mod 7: 2, 7 mod 11: 8, 4 mod 7: 2, 4 mod 16: none} |
| 15 | linear congruences | {3x≡6 (mod 9): [2, 5, 8], 2x≡1 (mod 4): []} |
| 16 | CRT reconstruction of 3-mer codes | {GAT from (0 mod 7, 2 mod 11): 35, (7,11) round-trips: 64, (5,13) round-trips: 64, (3,5) round-trips: 15, (4,6): refused} |
| 17 | ord_m(4) (pedagogical) | {7: 3, 13: 6, 17: 4, 31: 5, 61: 30} |
| 18 | φ(m) | {7: 6, 13: 12, 17: 16, 31: 30, 61: 60, 77: 60, 16: 8, 64: 32} |
| 19 | Euler: 4^φ(m) ≡ 1 (mod m) | {7: True, 13: True, 17: True, 31: True, 61: True, 77: True, 16: False} |
| 20 | two-substitution cancellation (Fermat periodicity) | {string: GATTAAGTAATG, d_H: 2, collide: [7, 13]} |
| 21 | B: pigeonhole bound ⌈30/m⌉ | {7: 5, 13: 3, 16: 2, 17: 2, 31: 1, 61: 1, 64: 1} |
| 22 | B: max load ≥ bound for every m | True |
| 23 | B: (occupied, max load, collisions) per m | {16: (1, 30, 29), 64: (1, 30, 29), 7: (7, 8, 23), 13: (11, 8, 19), 17: (15, 3, 15), 31: (20, 4, 10), 61: (24, 2, 6)} |
| 24 | B: last 3 bases of every haplotype (count) | {TGA: 30} |
| 25 | B: d_H = 1 pairs | 8 |
| 26 | B: d_H = 1 collisions, observed = predicted by §3 theorem | {16: 8, 64: 8, 7: 0, 13: 0, 17: 0, 31: 0, 61: 0} |
| 27 | B: d_H = 1 collisions for odd primes (theorem: none) | {7: 0, 13: 0, 17: 0, 31: 0, 61: 0} |
| 28 | B: canonical = distinct − present RC pairs (k = 4, 6, 8) | {4: True, 6: True, 8: True} |
| 29 | B: whole-space canonical classes (Burnside) | {4: 136, 6: 2080, 8: 32896} |
| 30 | B: k-mer census (distinct, canonical, palindromic, RC pairs present) | {4: {distinct: 248, canonical: 136, palindromic: 16, rc_pairs_present: 112}, 6: {distinct: 1610, canonical: 1224, palindromic: 20, rc_pairs_present: 386}, 8: {distinct: 2563, canonical: 2488, palindromic: 7, rc_pairs_present: 75}} |
| 31 | B: distinct 8-mers bucketed, max load ≥ ⌈N/m⌉ | {64: True, 61: True} |
| 32 | B: 8-mer buckets, m = 64 (last 3 bases) vs m = 61 | {64: {keys: 2563, occupied: 64, max_load: 121, bound: 41}, 61: {keys: 2563, occupied: 61, max_load: 55, bound: 43}} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | A: enc(S) | yes |
| 2 | A: enc(comp S) = 4^12 − 1 − enc(S) | yes |
| 3 | A: canonical min(enc S, enc RC S) | yes |
| 4 | A: decode(enc(S)) = S | yes (consistency) |
| 5 | A: buckets for m = 16, 7, 13, 17 | yes |
| 6 | A: colliding groups | yes |
| 7 | A: d_H = 1 pairs sharing a bucket | yes |
| 8 | A: pigeonhole, 9 keys in 7 buckets | yes |
| 9 | last bases seen: smallest j with m \| 4^j | yes |
| 10 | t_min for a bucket-preserving single substitution | yes |
| 11 | syn01 3-mer codes by sliding window | yes |
| 12 | syn01 k-mer census | yes |
| 13 | extended Euclid on (252, 198) | yes |
| 14 | modular inverses | yes |
| 15 | linear congruences | yes |
| 16 | CRT reconstruction of 3-mer codes | yes |
| 17 | ord_m(4) (pedagogical) | yes |
| 18 | φ(m) | yes |
| 19 | Euler: 4^φ(m) ≡ 1 (mod m) | yes |
| 20 | two-substitution cancellation (Fermat periodicity) | yes |
| 21 | B: pigeonhole bound ⌈30/m⌉ | yes |
| 22 | B: max load ≥ bound for every m | yes (consistency) |
| 23 | B: (occupied, max load, collisions) per m | n/a |
| 24 | B: last 3 bases of every haplotype (count) | n/a |
| 25 | B: d_H = 1 pairs | yes |
| 26 | B: d_H = 1 collisions, observed = predicted by §3 theorem | yes (consistency) |
| 27 | B: d_H = 1 collisions for odd primes (theorem: none) | yes |
| 28 | B: canonical = distinct − present RC pairs (k = 4, 6, 8) | yes (consistency) |
| 29 | B: whole-space canonical classes (Burnside) | yes |
| 30 | B: k-mer census (distinct, canonical, palindromic, RC pairs present) | n/a |
| 31 | B: distinct 8-mers bucketed, max load ≥ ⌈N/m⌉ | yes (consistency) |
| 32 | B: 8-mer buckets, m = 64 (last 3 bases) vs m = 61 | n/a |

**23/23 independent expected values matched; 5/5 consistency checks passed; 4 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

## Interpretation
Encoding is exact and invertible, and the complement identity holds. m = 16 and m = 64 see only the last 2 or 3 bases, so they collide on every substitution further from the end. On Dataset B every haplotype ends in TGA, so m = 16 and m = 64 put all 30 haplotypes into 1 and 1 bucket(s): the stop codon is all they can see. Every odd modulus ≥ 5 separates all single substitutions, as proved, while two substitutions can still cancel: an opposite pair (+δ, −δ) cancels exactly when ord_m(4) divides the distance between the positions (mod 7 and 13 in the example). Bucket-sharing is therefore a property of the residue, not of biological similarity. On Dataset B the {I,RC} quotient does real work at the k-mer level: for k = 8, 2563 distinct 8-mers collapse to 2488 canonical classes, and for k = 4 all 136 canonical classes of the whole space occur. For 8-mers, m = 64 occupies 64 of 64 buckets with max load 121, while m = 61 occupies 61 with max load 55 (observation).

## Limitations
Modular indexing here is a representation experiment, not an established clinical or genomic technique. The encoding A=0..T=3 is a convention; any bijection would do. CRT, Euler and Fermat are pedagogical. The k-mer census describes 30 haplotypes from one season and implies nothing about other data. Canonical k-mer counting is documented practice for reads of unknown orientation (Jellyfish manual v0.9, `math_audit.md`); the CDS here have known orientation, so the census is a mathematical demonstration.
