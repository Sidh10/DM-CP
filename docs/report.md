# Genomic Mutation & Evolution Analyzer
### An integrated computational framework using discrete mathematics

**MM0903 Discrete Mathematics · B.Tech CSE (Artificial Intelligence), second year · Vishwakarma Institute of
Technology, Pune · 2026–27**

---

## 1. Introduction

DNA is a string over the four-letter alphabet Σ = {A, C, G, T}. That is not an approximation or a modelling choice: it
is what the data is. Discrete mathematics therefore describes genomic sequences exactly, and the usual apology for
applied mathematics — that the model is a simplification — does not apply to the statement "these two sequences differ
at exactly one position."

This project builds a working framework in which equal-length DNA sequences are analysed using all six units of the
MM0903 syllabus, implemented by hand in pure Python, and in which every output can be traced back to a definition, a
derivation or a named algorithm. It is validated first on a ten-sequence synthetic dataset whose every expected answer
was worked out on paper before any code existed, and then applied to 36 influenza A(H3N2) hemagglutinin coding
sequences downloaded from GenBank.

Across eight reproducible experiments, 147 values computed by the implementation matched values fixed independently
beforehand, with 28 further internal-consistency checks and 35 measurements reported as observations. The test suite
has 179 tests; a mutation check injects 67 deliberate bugs and confirms the suite catches all of them.

## 2. Problem definition

Sequence comparison is normally performed with tools whose internal mathematics is invisible. For a mathematics course
that is the wrong way round. The question is not *what does the tool output* but *which mathematical object is being
computed, from which definition, with what guarantees, and under what assumptions*.

The problem addressed is therefore: **build an analysis framework for equal-length DNA sequences in which that question
has an answer at every step, and in which the kind of each claim is explicit** — raw observation, mathematical
calculation, algorithmic output, biological interpretation, or limitation.

This is a stricter requirement than correctness. A program that prints the right number while concealing how it got
there fails this project's criterion, and so does one that presents a computed tree as if it were evolutionary history.

## 3. Motivation and related work

**Motivation.** Three things make this problem a good fit for discrete mathematics rather than a decorative application
of it.

First, the objects are discrete and the statements are exact. A substitution is a checkable claim about two strings.

Second, the mathematics decides how results should be read. A 1701-nucleotide sequence has 3n = 5,103 possible
single-substitution neighbours; the best-connected haplotype in the real dataset has 6 of them. That counting argument
is what licenses calling the observed graph a thin sample rather than a picture of the mutational landscape — a
conclusion no amount of plotting would establish.

Third, the failure modes are mathematical. This project demonstrates a hash that destroys all information on real data
and then proves exactly when that can happen, which is more informative than either a working hash or a vague warning.

**Related work and its role here.** Two external sources were located, read, and are cited; nothing else is.

The canonical k-mer convention — representing a k-mer and its reverse complement by whichever comes first
lexicographically — is documented in the Jellyfish manual, which recommends it when the orientation of the input
sequences is unknown, as in sequencing reads [1]. This project implements the same idea arithmetically as
min(enc(S), enc(RC(S))) and notes the honest caveat that its own coding sequences have *known* orientation, so here the
construction is a mathematical demonstration rather than a necessity.

The data source and its terms are documented by NCBI, which states that it places no restrictions on the use or
distribution of GenBank data while noting that submitters may claim rights NCBI cannot assess [2, 3].

Deliberately absent: no claim is made about what transition-to-transversion ratio real genomic data typically shows,
because no citable source for it was located. The observed counts are reported and no conclusion is drawn from them.
Established phylogenetic software is out of scope; UPGMA is implemented here as a syllabus exercise in distance-based
clustering, not as a competitor to such tools.

## 4. Mathematical foundations

Full definitions, theorem statements, proofs and worked examples are in `docs/mathematical_design.md`. The results the
rest of this report relies on are:

**Theorem 1 (Unit I).** d_H(S,T) = |M(R,S) Δ M(R,T)| + |{i ∈ M(R,S) ∩ M(R,T) : S_i ≠ T_i}|. The second term exists only
because a position has three alternatives; it vanishes over a binary alphabet. This is why the project computes
allele-level sets M\* as well as position sets M: two sequences can share a *position* without sharing a *mutation*.

**Theorem 2 (Unit I).** The threshold relation R_k is reflexive and symmetric but not transitive for k ≥ 1. Its
transitive closure is an equivalence relation whose classes are exactly the connected components of the threshold
graph — the bridge between Units I and V.

**Theorem 3 (Unit II).** M(n,k) = C(n,k)·3^k satisfies M(n,k) = M(n−1,k) + 3M(n−1,k−1), by induction using Pascal's
rule; and Σ_k M(n,k) = 4^n by the binomial theorem.

**Theorems 4–7 (Unit III).** The string maps {I, R, C, RC} and the base permutations {id, (A T)(C G), (A G)(C T),
(A C)(G T)} are both Klein four-groups acting on different sets. Burnside's lemma counts their orbits. Substitution
operators form a **monoid with 5^n elements whose only unit is the identity**, which turns the specification's warning
that mutations do not form a group into a proof. Every substitution is realised by exactly one non-identity base
permutation, and because these are involutions the resulting transition/transversion classification is well defined on
*unordered* pairs — the only honest option when no ancestor is known.

**Theorems 8–11 (Unit IV).** Base-4 encoding is a bijection, decoded by the division algorithm;
enc(comp S) = 4^n − 1 − enc(S). The **exact criterion** for a single substitution to preserve a bucket depends only on
the odd part q of the modulus: if q ≥ 5 it never happens at any length; otherwise it happens exactly when a position
with t ≥ t_min exists, where t_min depends on whether δ = ±2 can contribute its own factor of 2. Two substitutions can
cancel exactly when ord_m(4) divides the distance between their positions, which is where Euler and Fermat enter.

**Units V and VI.** d_H is a metric, so a path cost never falls below the direct distance. The smallest threshold at
which the graph connects equals the largest edge of any minimum spanning tree (the bottleneck property). MSTs are
counted by contracting weight levels and applying Kirchhoff's matrix-tree theorem. UPGMA's output is always
ultrametric and reproduces an ultrametric input exactly, which is precisely why its fit to non-ultrametric data
measures the cost of its molecular-clock assumption.

## 5. Proposed methodology

The project was built in stages, each ending with a report and an explicit stop.

**The core methodological commitment is that expected values are derived by hand before the code exists.** For each
unit, a document works the answers out on paper; each experiment imports those values from a transcription module and
prints **Expected** and **Actual** in separate tables. A disagreement is therefore visible rather than absorbed.

Checks are labelled by what they can detect: *independent* (compared with a hand derivation, a hand-proved theorem's
prediction, or the quality-control script that shares no code with the implementation), *consistency* (can only fail if
the implementation disagrees with itself), and *observation* (no independent expected value; reported, never scored).
Counting these separately matters: 147 independent matches is a much stronger claim than 175 undifferentiated ones.

Because `genomath` imports no third-party library, the test suite can use NetworkX, SymPy and Biopython as independent
referees. And because passing tests say nothing about tests that cannot fail, a mutation check injects 67 deliberate
bugs and requires the suite to catch each one.

## 6. System architecture

Described in full, with diagrams, in `docs/system_design.md`. The structure is a one-way dependency chain:

```
data/ → genomath/ (pure, hand-written) → experiments/ · app/ · tests/
                                       ↑
                      scripts/ (fetch, QC, independent checkers)
```

`genomath` contains the mathematics and imports nothing beyond the standard library; it has no user-interface code and
no oracle libraries. `app/presentation.py` imports only `math`, so the dashboard's data shaping is unit-testable
without Streamlit. Both boundaries are enforced by tests that parse import statements rather than by convention.

Every public function returns `(result, trace)`, where the trace is a list of steps each tagged with one of five kinds.
Those five kinds are exactly the five categories the specification requires the interface to distinguish, so the
dashboard's labelling is a direct rendering of what the mathematics layer already recorded rather than a presentation
afterthought.

## 7. Implementation

Eight modules, one per unit plus the input layer, with 144 public functions and about 2,400 lines in `genomath`.
Course algorithms are written by hand: Dijkstra uses an array scan rather than a library heap, Kruskal uses a
hand-written merge sort and union–find with path compression and union by rank, and the graph layout in the dashboard
is a small force-directed implementation, since NetworkX is confined to the tests.

Three design choices are worth recording.

**Exact arithmetic where exactness matters.** UPGMA uses Python `Fraction`s throughout, so merge heights such as 19/12
and 23/14 are exact and the tie structure is not an artefact of floating-point error. Kirchhoff determinants are
computed over the rationals for the same reason.

**Guards instead of hangs.** The three genuinely exponential routines compute their output size from a closed form
first and refuse above a limit, so a refusal is instant and reports the number that would have been produced.

**Traces are bounded.** Per-element steps are capped at ten shown, with the remainder kept in the step's data, so
encoding a 1701-nucleotide sequence produces a 12-step trace rather than a 1,702-step one.

The dashboard presents seven sections, one per specification section, and recomputes everything live from `genomath`
rather than reading the committed results files, so it cannot display a stale number.

## 8. Experiments and results

Eight experiments, each writing `results/expN.md`. Totals: **147 independent expected values matched, 28 consistency
checks passed, 35 observations**. Per-experiment figures and the full result tables are in
`docs/testing_and_results.md`.

**Dataset A** reproduces every hand-derived value: the full 81-entry distance matrix; R_1 reflexive and symmetric but
not transitive with witness (syn01, syn02, syn03); components 9, 4, 2, 2, 1 at thresholds 0, 1, 2, 5, 6; a minimum
spanning tree of weight 15 with **8 distinct minimum trees**; and a UPGMA tree with exact heights and 4 tied merges.

**Dataset B** results are observations about 30 haplotypes from one season:

- **Sparsity.** The best-connected haplotype has 6 observed single-substitution neighbours out of 5,103 possible.
- **The single-substitution graph is a forest**: 8 edges, 22 components, no cycles. It connects only at k\* = 47,
  which equals the largest edge of the minimum spanning tree (weight 171, 4 distinct trees).
- **The bad hash, on real data.** Under m = 16 or m = 64 all 30 haplotypes fall into a *single* bucket, because every
  sequence ends in TGA and those moduli see only the last two or three bases. All 8 distance-1 pairs collide under
  them, and none under any odd modulus tested — exactly as Theorem 10 predicts.
- **Tree disagreement.** The distances are not ultrametric, so UPGMA reproduces only 9 of 435 distances exactly; 11
  merge steps involved ties; and the MST and UPGMA share 5 of their 10 and 27 non-trivial bipartitions.
- **Substitution classes.** Direction-free counts give 250 transitions to 91 transversions relative to the reference,
  and 131 to 40 along MST edges.

## 9. Discussion

**What the results show.** The framework runs end to end on real data and reproduces every independently known figure,
so the implementation is doing the mathematics it claims. Three findings are worth drawing out.

*The sparsity result changes how everything else should be read.* Six observed neighbours out of 5,103 possible means
the mutation graph is a very thin sample. A missing edge is therefore weak evidence of anything, and this single
counting fact constrains the interpretation of every graph and tree in the project.

*The bad hash is the clearest case of mathematics explaining a failure.* It would be easy to observe that m = 16 works
badly and move on. Theorem 10 says precisely why — the bucket sees only the last two bases — and precisely when it
cannot happen: whenever the odd part of the modulus is at least 5. The real dataset then supplies a striking
demonstration, since a shared stop codon collapses all 30 haplotypes into one bucket.

*The MST-versus-UPGMA comparison is the project's most important distinction.* The two structures are built from the
same distance matrix and disagree on most of their bipartitions. They answer different questions: an MST minimises
total edge weight among observed sequences, while UPGMA builds hypothetical ancestors under an assumption of equal
rates. Rooting the Dataset A MST at syn01 rather than syn04 changes its height from 3 to 4 and its leaf count from 6 to
5 without changing a single edge, which demonstrates in one example that an MST's rooting carries no information.

**On the transition counts.** Both observed ratios are far above the combinatorial baseline of 0.5 that equal
likelihood over the six unordered base pairs would give. This project reports that and stops. Explaining it would
require a substitution model and a citable source for what real data shows, and neither is available here.

**On method.** The practice of deriving expected values first repaid itself concretely. It caught an involutions
function that returned 6 instead of 3, and a bucket criterion that was wrong for δ = ±2 — the latter before any code
was written. It also worked in reverse once: a hand-derived cophenetic value was wrong, and the implementation's
disagreement exposed it. Separately, the mutation check found three tests too weak to detect injected bugs, and one
surviving mutant revealed both a missing test and a badly chosen mutant on my part.

## 10. Limitations

**Mathematical and modelling.**
- Substitutions only. Insertions and deletions change the length and leave Σ^n.
- Positional comparison assumes homology inferred from equal length and low divergence (≤ 3.4%), not from an alignment.
  A compensating insertion-deletion pair could in principle pass both checks.
- Hamming distance has no multiple-substitution correction, so it is not an evolutionary model.
- Mutation sets depend on the chosen reference; pairwise distances do not.
- The number-theoretic subsystem is a representation experiment. Modular indexing is not presented as an established
  genomic technique, and CRT, Euler and Fermat are labelled pedagogical.
- Unit III does no analytical work on full annotated coding sequences — a predicted and then measured result. Its
  real-data contribution is at the k-mer level and in the substitution classification.

**Data.** The 30 haplotypes are a convenience sample from one season, with 31 of 36 records deposited from the United
States. Records lacking a collection month were excluded. Two records with ambiguity codes were dropped rather than
resolved. Nothing here supports a claim about influenza populations.

**Interpretation.** Neither the graph nor either tree is evidence of evolutionary history. A path through observed
haplotypes is a route in a graph, not a trajectory. Substitution counts are direction-free, because no ancestor is
known. No clade or lineage is named, since no source for clade assignment was located.

**Software.** κ and λ are computed exhaustively and restricted to at most 16 vertices. No automated test renders the
dashboard pages, so a regression confined to the drawing code would not fail the suite. The FASTA upload path shares
the tested code path but was not exercised through the browser.

**Scope.** Not a diagnostic system, not a pathogenicity predictor, not an epidemiological tool, and not a replacement
for established phylogenetic software.

## 11. Conclusion and future work

The project set out to make all six units of the syllabus do real work on a single problem while keeping every claim
traceable and every limitation explicit. The six units are not evenly weighted — Units I, II, V and VI are structural,
while Units III and IV earn their place through the k-mer quotient, the substitution classification and two proved
theorems about representation — and the audit in `docs/math_audit.md` says so plainly rather than overstating the
weaker connections.

The outcome is a framework in which 147 independently fixed values are reproduced, two new results are proved and
machine-checked, and the real-data analysis is reported without overclaiming. The method that produced it — derive by
hand first, separate expected from actual, label what each check can detect, and attack the tests themselves — found
ten errors during development, including two in my own hand derivations.

**Future work, in order of value.**
1. *Alignment*, which would remove the equal-length restriction and admit insertions and deletions, at the cost of
   leaving Σ^n and needing a scoring model.
2. *A substitution model.* The transition counts are currently uninterpretable by design; a model and a citable source
   would make them a finding rather than a measurement.
3. *Neighbour joining*, which drops UPGMA's molecular-clock assumption and would sharpen the tree comparison by adding
   a third structure with different assumptions.
4. *Larger datasets*, which would need the O(N³) routines replaced — a heap-based Dijkstra and a better UPGMA — and
   would make the exhaustive connectivity search infeasible rather than merely slow.
5. *Rendering tests for the dashboard*, to close the one gap the current suite leaves.

## 12. References

This list is short deliberately. It contains only sources that were located and read; no reference is cited from
memory, and claims for which no source was found are reported as unsupported in the text rather than attributed.

[1] G. Marçais and C. Kingsford, *Jellyfish: A fast k-mer counter*, user manual version 0.9, 2010-10-01.
<https://www.cbcb.umd.edu/software/jellyfish/jellyfish-manual.html> — the `--both-strands` option defines the canonical
representation of a k-mer as whichever of it and its reverse complement comes first in lexicographic order, and
recommends it when input orientation is unknown. *Caveat: this is the 2010 manual for version 0.9; current releases
spell the option `-C`, and the Jellyfish journal paper was not located, so it is not cited.*

[2] National Center for Biotechnology Information, *GenBank Overview — GenBank Data Usage*.
<https://www.ncbi.nlm.nih.gov/genbank/about/> — "NCBI places no restrictions on the use or distribution of the GenBank
data", with the qualification that some submitters may claim intellectual property rights that NCBI cannot assess.

[3] National Center for Biotechnology Information, *NCBI Website and Data Usage Policies and Disclaimers*.
<https://www.ncbi.nlm.nih.gov/home/about/policies/> — the same position on molecular data, and the scripted-access
guidance the download script follows: no more than three requests per second, and supply the `email` and `tool`
parameters.

[4] Course slide decks for Units I–III, MM0903 Discrete Mathematics, VIT Pune (`PPT-UNIT-1-*`, `PPT-UNIT-2-*`,
`PPT-UNIT-3-*`), used for the syllabus topic names in the mapping table of `docs/mathematical_design.md`. The slide
folders for Units IV–VI were not available, so those topic names follow the project specification instead.

[5] Project specification, `docs/master_prompt.md` (supplied with the assignment), which defines the required units,
experiments, deliverables and prohibitions.

**Primary data.** 36 GenBank nucleotide records, accessions listed in `data/public_genomic/accessions.csv`, retrieved
2026-10-30 via NCBI Entrez; the exact query and filtering steps are recorded in
`data/public_genomic/selection_log.md`.
