# Stage 7 report: real-data evaluation (Experiment 8) and decision D9

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** Stage 8 (the Streamlit dashboard) has not started.
This stage implements **D9**, the last open decision. `docs/decisions.md` now records D1–D10 as implemented.

## What exists
| Path | Contents |
|---|---|
| `docs/unit3_expected.md` §9 (new) | hand derivation of the base-level Klein group: the three fixed-point-free involutions of Σ (3!! perfect matchings), closure worked out element by element, the Cayley table, the proof that every ordered base pair comes from exactly one non-identity element, why involutions make the classification direction-free, the 2-transitions / 4-transversions split of the 6 unordered pairs, the 1:2 combinatorial baseline, and the regular action with its Burnside check |
| `genomath/algebraic_structures.py` (extended) | `BASE_PERMS`, `BASE_K4`, `TRANSITION_ELEMENT`, `PURINES`/`PYRIMIDINES`, `fixed_point_free_involutions`, `base_cayley_table`, `classify_substitution`, `ordered_pair_coverage`, `substitution_census` |
| `experiments/expected_dataset_b.py` | Dataset B expectations drawn only from the genomath-free QC script or from hand derivations, with the source named per entry |
| `experiments/exp8_real_data.py` → `results/exp8.md` | the integrated six-unit run on Dataset B |
| `tests/test_algebra.py` | 4 new D9 tests (24 in that file, 162 in total) |
| `scripts/mutation_check.py` | 6 D9 mutants (63 in total) |

The base-level group reuses the Stage 3 axiom checker (`check_axioms`, `element_orders`), so the same verified machinery
is applied to a second group.

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python -m experiments.exp8_real_data
.venv/Scripts/python scripts/mutation_check.py        # ~25 minutes; run it in the background
```
Raw output from my last run:
```
162 passed in 22.73s
wrote results\exp1.md: 20/20 independent expected values matched; 2/2 consistency checks passed; 3 observation-only items
wrote results\exp2.md: 16/16 independent expected values matched; 1/1 consistency checks passed; 2 observation-only items
wrote results\exp3.md: 29/29 independent expected values matched; 1/1 consistency checks passed; 1 observation-only items
wrote results\exp4.md: 23/23 independent expected values matched; 5/5 consistency checks passed; 4 observation-only items
wrote results\exp5.md: 14/14 independent expected values matched; 2/2 consistency checks passed; 2 observation-only items
wrote results\exp6.md: 14/14 independent expected values matched; 4/4 consistency checks passed; 4 observation-only items
wrote results\exp7.md: 11/11 independent expected values matched; 2/2 consistency checks passed; 5 observation-only items
wrote results\exp8.md: 20/20 independent expected values matched; 11/11 consistency checks passed; 14 observation-only items
genomath: no banned imports
62/62 mutants caught
```
All eight of the spec's experiments now run and are green. Oracles for D9 in `tests/test_algebra.py`:
- brute force over all 24 permutations of ACGT, keeping those that are involutions with no fixed point (exactly 3,
  matching the hand derivation);
- SymPy `PermutationGroup`: order 4, abelian, not cyclic, transitive, one orbit;
- an **independent definition** of a transition (the substitution keeps the purine/pyrimidine class) compared against the
  group-element classification for all 12 ordered pairs;
- an independent recount of the Dataset B census that classifies by class membership instead of by group element.

## Experiment 8 summary
20/20 independent, 11/11 consistency, 14 observation-only.

**Everything independently known is reproduced by the integrated pipeline:**
- QC figures: 36 records read, 2 dropped for ambiguity codes, 34 clean, 30 haplotypes, 1701 nt, 435 pairs, 8 pairs at
  distance 1, distances from 1 to 57, and both groups of identical sequences.
- Unit II hand values: 4^1701 has 1025 decimal digits; 5103 possible single-substitution neighbours; |B(n,1)| = 5104;
  a repeated 5-mer is forced, and all 34 sequences contain one.
- Unit III prediction: 0 reverse-complement pairs, 0 RC-fixed sequences, 30 orientation classes.
- Unit IV theorems: m = 16 sees 2 bases and m = 64 sees 3; those two moduli collide on all 8 distance-1 pairs, and
  every odd modulus on none.
- D9 hand values: 3 involutions, 12 ordered pairs each from exactly one element, 2 transitions vs 4 transversions,
  baseline 0.5.

**Observations about this sample, with no biological claim:**
- **Sparsity.** The best-connected haplotype has 6 observed single-substitution neighbours out of 5103 possible
  (1.18 × 10⁻³).
- **Bad hash on real data.** Under m = 16 and m = 64 all 30 haplotypes occupy a single bucket, because every haplotype
  ends in TGA and those moduli see only the last 2–3 bases. m = 61 spreads them over 24 buckets.
- **Structure.** G_=1 is a forest: 8 edges, 22 components, no cycles, maximum degree 6. The graph first connects at
  k\* = 47, which equals the largest MST edge (weight 171, 4 distinct MSTs).
- **Trees.** The distances are not ultrametric, so UPGMA reproduces only 9 of 435 distances exactly; 11 merge steps
  involved ties; the MST and UPGMA share 5 of their 10 and 27 non-trivial bipartitions.
- **Substitution classes (D9).** Direction-free counts: 250 transitions vs 91 transversions among differences from the
  reference (ratio 2.75), and 131 vs 40 along MST edges (ratio 3.28).

### On the transition/transversion numbers
Both ratios are far above the combinatorial baseline of 0.5, which is what equal likelihood over the six unordered base
pairs would give. The experiment reports this as a count and draws **no conclusion from it**. Explaining the gap would
need a substitution model and a source for what ratio real genomic data shows, and this project has neither. I am aware
of the general claim that transitions outnumber transversions, but I have not located a citable source, so the results
file does not assert it.

A cross-unit relation is now asserted rather than left as a coincidence: the substitution census along MST edges totals
171, which **must** equal the MST weight, because each edge contributes one substitution per differing position and the
weight is the sum of those d_H values. It is tagged as a consistency check.

## Complexity of the D9 additions
| Function | Time |
|---|---|
| `fixed_point_free_involutions`, `base_cayley_table`, `ordered_pair_coverage` | O(1): the alphabet and group are fixed size |
| `classify_substitution` | O(1) |
| `substitution_census` over P pairs of length n | Θ(P·n) |

## Problems found and fixed during the stage
- **The hand derivation caught a real bug.** My first `fixed_point_free_involutions` generated each pairing twice and
  reported 6, while the hand derivation said 3 (3!! perfect matchings). The loop paired every unordered pair with its
  complement, so (A C)+(G T) and (G T)+(A C) both appeared. Fixing it to pair the first base with each of the other
  three gives exactly 3, and the old form is now a mutant the suite catches. This is the second time (after the Stage 6
  cophenetic value) that writing the expected values by hand first caught something the code got wrong.
- **A test contained a no-op assertion** (`assert al.burnside.__name__`), left over from reorganising a comment. Removed,
  and replaced with the direct Burnside count for this action, since the string-level `burnside()` is a different action.
- **My shell heredocs mangled `\n` escapes for the third time** (Stages 3, 6, 7), each time inside a multi-line mutation
  entry. Fixed by direct edit, and both new mutants were rewritten to target single lines. For the rest of the project I
  will add mutation entries with the editor rather than through a shell heredoc.

## Open issues and doubts
- **Mutation-check runtime is now about 25 minutes** for 63 mutants, because each mutant runs the whole suite until the
  first failure and the suite re-runs all eight experiments. If this becomes a nuisance, I can add a mode that runs only
  the test file matching the mutated module; say the word and I will.
- The Dataset B tree figures depend on D8, since 11 merge steps had ties and the reverse rule changes the topology. Any
  figure must name the rule.
- "Robinson–Foulds" remains the usual name for the bipartition symmetric difference; I report the count without the name
  because I have not located the original reference.
- The 30 haplotypes are a convenience sample from one season, dominated by United States deposits. Every Dataset B
  statement in `results/exp8.md` is scoped to the sample.
- **Remaining stages:** 8, the Streamlit dashboard in `app/`, which renders the traces and tags each block with the five
  categories; and 9, the written deliverables (proposal, mathematical design document, system design, report, slides and
  viva questions).

## Mutation-check note
The run reported 62 of 63 caught. The survivor was **"base Cayley table composes in the other order"**, which is an
**equivalent mutant**: the base-level group is abelian, so a∘b = b∘a for every pair and the mutated table is literally
identical. I verified this directly (both tables compared equal entry by entry) rather than assuming it, then removed
the entry with the reason recorded in `scripts/mutation_check.py`. The suite stands at **62/62**.

This is the same equivalence I dropped for the string-level Klein group in Stage 3; I reintroduced it by accident for the
base-level group.
