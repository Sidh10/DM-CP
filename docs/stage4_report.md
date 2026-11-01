# Stage 4 report: Unit IV (number theory), Experiment 4

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 5 has not started.

## Decisions applied provisionally
- **Decision 5:** canonical k-mer census of Dataset B for k = 4, 6, 8 (`kmer_census`). This gives Unit III real-data work
  at the k-mer level.
- **Decision 8:** modulus set {16, 64, 7, 13, 17, 31, 61}.
- Not implemented: the Stage 3 proposal (a base-level Klein group classifying transitions and transversions). It still
  needs your decision.
- Stages 1–2 provisional decisions are unchanged.

## What exists
| Path | Contents |
|---|---|
| `docs/unit4_expected.md` | proof that encoding is a bijection, with decoding by the division algorithm; complement and RC on codes; sliding-window example; syn01 k-mer census; **exact single-substitution collision theorem**; Fermat/Euler orders and a two-substitution cancellation example; Euclid, inverses and congruences; CRT on 3-mers; pigeonhole bounds |
| `genomath/number_theory.py` | `division`, `divides`, `gcd`, `extended_gcd`, `mod_inverse`, `solve_linear_congruence`, `crt` (checks pairwise coprimality first), `factorize`, `euler_phi`, `mod_pow` (square-and-multiply), `multiplicative_order`, `encode` (Horner), `decode` (repeated division), `complement_code`, `rc_code` (RC computed arithmetically on the code), `canonical_code`, `kmer_codes` (sliding window), `kmer_census`, `bucket_report` (with pigeonhole bound), `last_bases_seen`, `substitution_collision_t_min`, `substitution_collisions` |
| `experiments/expected_unit4.py`, `experiments/exp4_modular_indexing.py` → `results/exp4.md` | Experiment 4 |
| `tests/test_number_theory.py` (16 tests); `test_experiments.py` now covers Exp 1–4 | |
| `scripts/mutation_check.py` | 10 Unit IV mutants added (35 in total) |
| `docs/math_audit.md` | Unit IV entry now points to the sharper theorem |

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp4_modular_indexing
.venv/Scripts/python scripts/mutation_check.py
```
Raw output from my last run:
```
126 passed in 14.91s
wrote results\exp1.md: 22/22 expected values matched
wrote results\exp2.md: 17/17 expected values matched
wrote results\exp3.md: 30/30 expected values matched
wrote results\exp4.md: 28/28 expected values matched
35/35 mutants caught
```
Oracles in `test_number_theory.py`:
- `math.gcd` and `pow(a, -1, m)`, over all a, b in −30..30 and all m < 60;
- SymPy `factorint`, `totient` (n < 400), `n_order` and `crt` (all 64 3-mer codes under four moduli sets);
- Python's `int(…, 4)` as an independent encoder;
- brute-force solutions of every linear congruence for m < 25;
- brute-force search over all single substitutions of all strings with n ≤ 5, for every m < 70, against the
  collision theorem.

## Experiment 4 summary (from `results/exp4.md`)
28/28 expected values matched; 4 items are observation-only. *(Corrected in the post-review section: originally
said 5. After the review, 5 of the 28 are reclassified as consistency checks, so it is 23/23 independent + 5/5 consistency.)*
- **Dataset A.** Codes, complement codes, canonical codes, buckets for m = 16/7/13/17 and colliding groups all equal the
  hand tables. m = 16 collides on 4 of the 6 d = 1 pairs, and every odd modulus on 0.
- **Theory.** `last_bases_seen`: 16 → 2, 64 → 3, odd → none. The t_min values match the theorem. Extended Euclid gives
  (18, 4, −5). Inverses: 2, 8, 2, and none for 4 mod 16. Congruences: {2,5,8} and ∅.
- **CRT.** 64/64 round-trips for (7,11) and (5,13); 15/64 for (3,5); (4,6) refused.
- **Orders and φ** (pedagogical). Orders 3, 6, 4, 5, 30; the φ values match; Euler's congruence holds exactly when gcd(4, m) = 1.
  The cancellation string GATTAAGTAATG (d_H = 2 from syn01) collides under exactly {7, 13}.
- **Dataset B.**
  - Under m = 16 and m = 64 **all 30 haplotypes fall in one bucket**, because every haplotype ends in TGA and those
    moduli see only the last 2–3 bases. This is the bad hash on real data, explained by a computed observation.
  - The 8 d_H = 1 pairs collide 8/8 under 16 and 64, as the theorem predicts, and 0/8 under every odd prime.
  - k-mer census: k = 4 has 248 distinct 4-mers covering all 136 canonical classes; k = 6 has 1610 → 1224; k = 8 has
    2563 → 2488. The identity canonical = distinct − (RC pairs present) holds for every k.
  - Distinct 8-mers: m = 64 has max load 121 against a pigeonhole bound of 41; m = 61 has max load 55 against 43.

## Complexity
| Algorithm | Time |
|---|---|
| `gcd`, `extended_gcd` | O(log min(a,b)) division steps (Lamé) |
| `mod_inverse`, `solve_linear_congruence` | O(log m) plus O(g) to list solutions |
| `crt` (t moduli) | O(t²) coprimality checks + O(t log M) |
| `factorize`, `euler_phi` | O(√n) trial division |
| `mod_pow` | O(log e) multiplications |
| `multiplicative_order` | O(ord) ≤ O(φ(m)) |
| `encode`, `decode`, `rc_code` | Θ(n) big-integer steps (each on numbers up to 4^n) |
| `kmer_codes` | Θ(L) after the first window |
| `kmer_census` | Θ(N·L) codes + O(D·k) for D distinct codes |
| `bucket_report` | Θ(N) |

## Problems found and fixed during the stage
- **My first statement of the collision theorem was wrong at one edge.** Writing the proof showed that δ = ±2 brings its
  own factor of 2, so for odd part q = 1 collisions start at t = ⌊a/2⌋, not ⌈a/2⌉ (e.g. m = 8 collides at t = 1).
  The doc was corrected before coding, the brute-force test confirms the corrected formula, and the original mistake is
  now a mutant that the suite catches ("t_min ignores the delta = 2 case").
- A hand example named the wrong 3-mer for code 15 (AGT instead of ATT); corrected.
- The brute-force collision test first took 20 s because it called the traced encoder about 2 million times. It now uses
  Python's own base-4 parser as the encoder, which is also an independent oracle. The test now takes 1.9 s.
- One interpretation sentence overstated when two substitutions cancel. It now says "an opposite pair (+δ, −δ)
  cancels exactly when ord_m(4) divides the distance between the positions".
- **One equivalent mutant was dropped:** removing CRT's explicit coprimality check. With non-coprime moduli, some
  M_i = M/m_i shares a factor with m_i, so the inverse step fails and CRT is still refused, only with a different message.

## Open issues and doubts
- "TGA is a stop codon" is standard genetic-code knowledge that I have not cited. The data itself shows only that every
  haplotype ends in TGA. The QC script separately confirmed a single terminal stop under the standard code (Biopython
  translation).
- The k-mer census is a **mathematical demonstration**: these CDS have known orientation, so canonicalisation is not
  needed for them biologically. This is stated in the experiment's limitations.
- Suite runtime is now about 15 s, mostly from re-running Experiments 3–4 inside `test_experiments.py`.
  The mutation check runs the suite 35 times (about 9 minutes).
- `encode` emits one trace step per base (1701 steps for a Dataset B sequence). That is fine for the experiments, but
  the UI will need to collapse it.
- Still open: the UPGMA tie-break (Stage 6) and the Stage 3 transitions/transversions proposal.

## Post-review changes (2026-10-30, after your Stage 4 review)
Your answers: confirm all provisional decisions; add the transitions/transversions group in Stage 7; UPGMA tie-break by
lexicographic pair; summarise long traces now. All are recorded in **`docs/decisions.md`** (D1–D10), which is now the
single place where decisions live.

**1. Capped trace steps (D10).** `Trace.add_steps(kind, texts)` keeps the first 10 per-element steps and adds one
"… N more steps" step whose `data["omitted"]` holds the rest. It is used by `mutation_set`, `set_operations`, `encode`,
`decode` and `binomial`.
- Example: `encode` of a 1701-nt sequence now has 12 trace steps instead of 1702 (10 shown, 1 omission step holding 1691
  texts, 1 result).
- Two tests were added. The suite became faster, since less trace text is built.

**2. Independent vs consistency checks.** Reviewing all four experiments (not just Exp 4) gave **9** checks that can
only fail if `genomath` disagrees with itself: a round trip, two genomath computations agreeing, or a bound that holds
for any data.
- Exp 1: d_H identity on A and on B.
- Exp 2: generated strings are distinct and at distance k.
- Exp 3: orbit–stabiliser.
- Exp 4: decode∘encode, max load ≥ pigeonhole bound (twice), observed = predicted collisions, and the canonical-count
  identity.

They are now tagged `consistency=True`. The results files mark them, and each headline counts them separately:
```
wrote results\exp1.md: 20/20 independent expected values matched; 2/2 consistency checks passed; 3 observation-only items
wrote results\exp2.md: 16/16 independent expected values matched; 1/1 consistency checks passed; 2 observation-only items
wrote results\exp3.md: 29/29 independent expected values matched; 1/1 consistency checks passed; 1 observation-only items
wrote results\exp4.md: 23/23 independent expected values matched; 5/5 consistency checks passed; 4 observation-only items
```
In my review message I said Exp 4 had 2 such checks. It has 5, and Exps 1–3 have 4 more.

**3. Correction.** This report originally said Exp 4 had 5 observation-only items. It has 4.

**Verification after the changes:** `128 passed`; all four experiments as above; the mutation check result is in the
next line.
```
35/35 mutants caught
```
(The mutation check now takes about 10 minutes because it runs the suite 36 times. Running it in the background is recommended.)
