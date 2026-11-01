# Stage 3 report: Unit III (algebraic structures), Experiment 3

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 4 has not started.
No new provisional decisions were needed. Scope follows your Stage 0 rules: {I, RC}, the Klein four-group {I, R, C, RC},
the RC-fixed count, Burnside, and **no arbitrary position permutations**. The only position permutation used is reversal.

## What exists
| Path | Contents |
|---|---|
| `docs/unit3_expected.md` | hand-derived Cayley tables and the general-length group proof; the faithfulness subtlety (n = 0, 1); subgroups and cosets; fixed-point formulas and table n = 0..8, 12; Burnside counts; orbit-size distributions; the complement as (A T)(C G); the proof that point substitutions form a monoid of size 5^n with only the identity as a unit; concatenation laws; expected dataset results with the start/stop-codon argument for Dataset B |
| `genomath/algebraic_structures.py` | `transform`, `identify` (exhaustive function equality on all words of length ≤ 6), `distinct_on`, `cayley_table`, `check_axioms` (closure, associativity over all triples, identity, inverses, commutativity, each with a witness; classifies the result as group, monoid or neither), `element_orders`, `check_isomorphism`, `subgroups`, `left_cosets`, `fixed_count` (formula), `fixed_count_enumerated`, `burnside`, `orbit` (with stabiliser and an orbit–stabiliser check), `canonical`, `orbits_enumerated`, `dataset_orbits`, `permutation_cycles`, `substitution_monoid`, `concatenation_laws` |
| `experiments/expected_unit3.py`, `experiments/exp3_algebraic_verification.py` → `results/exp3.md` | Experiment 3 |
| `tests/test_algebra.py` (20 tests); `test_experiments.py` now covers Exp 1–3 | |
| `scripts/mutation_check.py` | 9 Unit III mutants added (25 in total) |
| `genomath/trace.py` | `add(kind, text, /, **data)` is now positional-only, so a data field may be named `kind` |

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp3_algebraic_verification
.venv/Scripts/python scripts/mutation_check.py
```
Raw output from my last run:
```
109 passed in 7.28s
wrote results\exp1.md: 22/22 expected values matched
wrote results\exp2.md: 17/17 expected values matched
wrote results\exp3.md: 30/30 expected values matched
25/25 mutants caught
genomath: no banned imports
```
Oracles in `test_algebra.py`:
- **SymPy `PermutationGroup`** built from the four maps acting on Σ^n: group order (2 on Σ^1, 4 on Σ^2), abelian/cyclic, and orbit counts for n = 1..5 against Burnside.
- **SymPy `Permutation`**: cycle structure, order and parity for all 120 permutations of 5 points.
- Enumeration of fixed points and orbits for n ≤ 8, and the hand tables.
- An exhaustive check that none of the 24 bijections V4 → Z4 is a homomorphism.
- A deliberately non-associative operation, (a − b) mod 3, to show that `check_axioms` detects failures as well as successes.

## Experiment 3 summary (from `results/exp3.md`)
30/30 expected values matched; 1 item is observation-only (the non-injectivity witness).
- **Groups.** Both Cayley tables equal the hand tables, and every axiom holds for both groups. Element orders are 1, 2, 2, 2.
  V4 is not cyclic, is isomorphic to Z2×Z2, and has 0 homomorphic bijections onto Z4.
- **Subgroups and cosets.** There are 5 subgroups; the cosets of {I,RC} are {I,RC} and {R,C}.
- **Faithfulness.** On Σ^0 all four maps coincide; on Σ^1, R = I and RC = C; on Σ^2 all four are distinct.
- **Counts.** Fixed points agree by formula and by enumeration for n = 0..8. Burnside equals enumeration for both groups
  for n = 0..8. At n = 12, Burnside gives 8,390,656 orbits under {I,RC} and 4,196,352 under V4.
  Orbit-size distributions match for n = 1..4, and orbit–stabiliser holds for every orbit.
- **Substitution monoid.** For n = 1, 2, 3 it has 5, 25 and 125 elements, its only unit is the identity, and it is not
  a group and not commutative. Witness: f_AA(AA) = f_AA(AC) = AA.
- **Dataset A.** syn03 and syn08 merge, giving 8 classes from 9 haplotypes. syn01's V4 orbit has size 2 with
  stabiliser {I,RC}. The canonical representative of the syn03/syn08 pair is GATTACGCAAGC.
- **Dataset B.** 0 RC pairs, 0 RC-fixed, 30 classes, as predicted from the start/stop-codon argument.

## Complexity
| Algorithm | Time | Notes |
|---|---|---|
| `identify` / `cayley_table` | O(\|G\|² · \|G\| · W · n) with W = Σ_{k≤6} 4^k = 5,461 words | exhaustive by design |
| `check_axioms` | O(\|G\|³) table lookups | 64 triples for V4; 25³ = 15,625 for the n = 2 monoid |
| `fixed_count`, `burnside` | O(\|G\|) big-integer powers | any n |
| `fixed_count_enumerated`, `orbits_enumerated` | Θ(4^n · \|G\| · n) | refused above n = 10 |
| `subgroups` | O(2^\|G\| · \|G\|²) | fine for \|G\| ≤ 16 |
| `substitution_monoid(n)` | O(5^n · 4n · 4^n · n) closure, plus 5^(2n) · 4^n for the table | n ≤ 3 in practice |
| `permutation_cycles` | O(m) for m points | |

## Problems found and fixed during the stage
- `Trace.add` crashed when a data field was named `kind`, because it clashed with the parameter name. The fix makes the
  parameters positional-only; all earlier tests still pass.
- A mutation entry was mangled by shell escaping while I was adding it, so the checker failed to parse.
  I fixed it by a direct edit. This was a tooling slip, not a code defect.
- **Two candidate mutants were dropped as equivalent, not hidden.** Swapping the composition order changes nothing in an
  abelian group. Dropping "contains the identity" from the subgroup test changes nothing, because any non-empty closed
  subset of a finite group contains the identity.

## Open issues and doubts
- **The faithfulness subtlety deserves a sentence in the report and viva.** "The Klein four-group acts on Σ^n" is only
  faithful for n ≥ 2. Burnside's lemma is unaffected, but "the four maps are distinct" is false on Σ^1.
- **Confirmed weak link: Unit III does nothing on Dataset B as full CDS.** It is now a measured, predicted result rather
  than an assertion. Two ways to give Unit III real-data work, both needing your decision:
  1. *(already open, decision 5, planned for Stage 4)* Canonical k-mers of Dataset B under {I,RC}.
  2. *(new proposal, not implemented)* The three fixed-point-free involutions of Σ, together with the identity, form a
     second Klein four-group, acting on bases rather than strings: complement (A T)(C G), (A G)(C T) and (A C)(G T).
     Each substitution x → y with y ≠ x is produced by exactly one of the three. (A G)(C T) is the purine↔purine /
     pyrimidine↔pyrimidine swap, i.e. the transitions; the other two are transversions. On Dataset B this would
     classify every observed substitution, giving a real-data count of transitions vs transversions. I would report the
     counts only and make no claim about their typical ratio, because I have not located a source for that.
- Function equality is exhaustive only up to length 6; general n rests on the written proof.
- The test suite now takes about 7–15 s, mostly from enumerating orbits for n ≤ 8 inside Experiment 3. Acceptable, but
  noted in case you run it often.
