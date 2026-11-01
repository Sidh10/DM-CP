# Deliverable A — Project proposal

**Title.** Genomic Mutation & Evolution Analyzer: an integrated computational framework using discrete mathematics

**Course.** MM0903 Discrete Mathematics · B.Tech CSE (Artificial Intelligence), second year · Vishwakarma Institute of
Technology, Pune · 2026–27

## Abstract
DNA over the alphabet Σ = {A, C, G, T} is a discrete object, so discrete mathematics describes it exactly rather than
approximately. This project builds a framework in which equal-length DNA sequences are analysed by hand-written
implementations of the six MM0903 units, and in which every output can be traced back to a definition, a derivation or a
named algorithm. Sequences are modelled as functions [n] → Σ; their differences as sets; the "few substitutions apart"
relation as a graph; and their summary structures as spanning and clustering trees. The framework is validated on a
10-sequence synthetic dataset whose every expected answer was derived by hand before any code was written, and then
applied to 36 influenza A(H3N2) hemagglutinin coding sequences downloaded from GenBank. Across eight reproducible
experiments, **147 values computed by the implementation matched values fixed independently beforehand**, with 28
further internal-consistency checks and 35 measurements reported as observations. Two mathematical results were proved
in the course of the work: an exact criterion for when a single substitution preserves a modular bucket, and a proof
that substitution operators form a monoid rather than a group.

## Problem statement
Comparing genomic sequences is usually done with tools whose internal mathematics is invisible to the user. For a
mathematics course this is backwards: the question is not *what does the tool output* but *which mathematical object is
being computed, from which definition, with what guarantees and what assumptions*. The problem addressed here is to
build a working analysis framework in which that question has an answer at every step, and in which claims are separated
by kind — raw observation, mathematical calculation, algorithmic output, biological interpretation, and limitation.

## Motivation and real-world relevance
- **The data are genuinely discrete.** A substitution is an exact, checkable statement about two strings; no modelling
  approximation is involved in saying that two sequences differ at exactly one position.
- **The mathematics is load-bearing, not decorative.** The sparsity result is an example: a 1701-nucleotide sequence has
  3n = 5,103 possible single-substitution neighbours, and the best-connected haplotype in the real dataset has 6 of
  them. That is a counting argument deciding how a dataset should be interpreted.
- **The failure modes are mathematical too.** The project demonstrates a hash that destroys all information on real data
  (all 30 haplotypes land in one bucket under m = 16 or 64, because they all end in TGA), and proves exactly when that
  can and cannot happen.
- **Honest limits matter in biology.** Distinguishing "this is the cheapest tree on these distances" from "this is the
  evolutionary history" is a mathematical distinction with real consequences, and the project makes it explicit.

## Objectives
1. Implement all six syllabus units as hand-written pure-Python algorithms, each returning its result together with a
   step-by-step trace.
2. Derive every expected answer for a small synthetic dataset by hand, before implementation, so the code is checked
   against mathematics rather than against itself.
3. Validate each unit with an experiment that reports expected and actual values in separate fields.
4. Apply the validated pipeline to a real, openly licensed genomic dataset and report what it observes without
   overclaiming.
5. Present the results in a dashboard where every figure is accompanied by the calculation that produced it.

## Research questions
**Central.** How can discrete mathematical structures — sets, relations, counting, group actions, modular arithmetic,
graphs and trees — represent equal-length DNA sequences, analyse their substitution relationships, and yield graph and
tree summaries whose every output traces back to a stated definition or algorithm?

Supporting, each with the experiment that answers it:
1. Do mutation sets and the threshold relation reproduce all pairwise distances, and where does the relation fail to be
   transitive? (Experiment 1)
2. Does bounded enumeration match the closed-form count, and how sparse are observed neighbours? (Experiments 2, 8)
3. Do the strand-symmetry transformations satisfy the group axioms, and do Burnside's counts match enumeration?
   (Experiment 3)
4. How does the modulus change bucket collisions, and when can a single substitution be invisible? (Experiment 4)
5. How do edges and components change with the distance threshold, and when does a shortest path cost more than the
   direct distance? (Experiments 5, 6)
6. How do a minimum spanning tree and a UPGMA tree differ in structure and in what they license? (Experiments 6, 7)

## Scope
**In scope.** Equal-length sequences over {A, C, G, T}; substitutions only; the eight `genomath` modules; eight
reproducible experiments; a pytest suite with independent library oracles; a Streamlit dashboard; this documentation.

**Out of scope.** Sequence alignment; insertions and deletions; diagnosis, pathogenicity prediction or epidemiological
forecasting; clade or lineage naming; substitution models; maximum-likelihood or Bayesian phylogenetics; machine
learning.

## Expected contributions
1. A complete, tested implementation of all six units applied to one coherent problem, with 144 public functions that
   each return a result and a trace.
2. A validation method in which expected values are hand-derived first and recorded separately from results, so that a
   disagreement is visible rather than absorbed. This caught three errors during development — two in the code and one
   in a hand derivation.
3. Two proved results beyond the syllabus minimum: the exact single-substitution bucket criterion (depending only on the
   odd part of the modulus), and the monoid proof that makes "mutations do not form a group" a theorem rather than an
   assertion.
4. A dashboard that tags every result by kind, so a reader can tell a measurement from an interpretation.

## Limitations
Stated in full in `docs/final_spec.md` §2 and in each experiment's own limitations section. The principal ones: the
analysis assumes positional homology from equal length and low divergence rather than from an alignment; the real
dataset is a convenience sample from one season dominated by United States deposits; substitution counts are
direction-free because no ancestor is known; neither the graph nor either tree is evidence of evolutionary history; and
the number-theoretic subsystem is a representation experiment, not an established genomic technique.
