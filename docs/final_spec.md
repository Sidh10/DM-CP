# Phase 2: Final project specification

Title: **Genomic Mutation & Evolution Analyzer: an integrated computational framework using discrete mathematics**
(MM0903, VIT Pune, 2026–27). This document fixes *what* will be built and measured. `math_audit.md` justifies *why*.
Items marked **[decision]** are open until you confirm them (they are collected in `stage0_report.md`).

## 1. Research question
**Central.** How can discrete mathematical structures (sets, relations, counting, group actions, modular arithmetic,
graphs and trees) represent equal-length DNA sequences, analyse their substitution relationships, and yield graph and tree
summaries whose every output traces back to a stated definition or algorithm?

Supporting questions, each with the measurement that answers it:
1. Sets/relations: do mutation sets and R_k reproduce all pairwise distances (via the d_H identity), and where does R_k fail to be transitive? → Exp 1
2. Counting: does bounded enumeration match C(n,k)·3^k, and how sparse are observed neighbours relative to the 3n possible? → Exp 2, Exp 8
3. Algebra: do {I,RC} and V4 satisfy the group axioms exhaustively, and do Burnside counts match brute force? → Exp 3
4. Number theory: how does the modulus change bucket collisions, and do the two hashing theorems hold on real data? → Exp 4
5. Graphs: how do edges and components change with k, and when does a shortest path cost more than the direct d_H? → Exp 5, Exp 6
6. Trees: how do MSTs (Prim, Kruskal) and a UPGMA tree differ in structure and in what they can claim? → Exp 6, Exp 7

## 2. Scope
In scope: equal-length sequences over {A,C,G,T}; substitutions only; the eight `genomath` modules; eight experiments;
pytest validation; a Streamlit dashboard (Stage 8); the report and viva material (Stage 9).
Out of scope (spec §11): alignment, indels, diagnosis or pathogenicity, forecasting, clade naming, substitution models,
maximum-likelihood or Bayesian phylogenetics, and machine learning.

## 3. Engineering rules (from your Stage 0 decisions)
- Package `genomath/` contains hand-written pure-Python algorithms and **no Streamlit imports**. NetworkX, SymPy and
  `math.gcd` appear **only in tests** as oracles. Biopython appears only in `scripts/` for download and parsing.
- Every public algorithm returns `(result, trace)`. Proposed trace shape **[decision]**: a list of steps, each
  `{"kind": ..., "text": str, "data": dict}`, where `kind` is one of `observation | calculation | algorithm | interpretation | limitation`.
  These are your five UI tags. The UI renders the steps and tags each result block with the same vocabulary.
- `pathlib` everywhere; no OS-specific code. Short comments explaining design choices only.
- Each experiment `experiments/expN_*.py` writes `results/expN.md` with **separate fields** "Expected (source: …)" and
  "Actual (this run)", plus Match, Interpretation and Limitations. Expected values come only from `docs/synthetic_expected.md`
  (Dataset A) or from independently computed QC values (Dataset B). They are never written as results.

## 4. Data strategy
### Dataset A (synthetic, `data/synthetic/dataset_a.fasta`)
Ten sequences of length 12, chosen so that every structural case has a hand-checkable answer: identical pair, RC-palindrome,
4-step chain, MST ties (8 distinct MSTs), an RC pair, a shared-position/different-allele pair, and an outlier isolated for
k ≤ 5. All expected answers are in `docs/synthetic_expected.md`. `dataset_a_invalid.fasta` holds one N-containing and one
wrong-length record for validation tests.

### Dataset B (public, `data/public_genomic/`)
- **What.** Influenza A virus subtype H3N2, segment 4 hemagglutinin (HA) **coding sequences**, human host, collected in the
  2017–18 Northern-hemisphere season (Oct 2017 – Apr 2018). 36 records were selected; all have a 1701-nt CDS.
- **Source and terms.** GenBank (NCBI nucleotide), via Entrez E-utilities only. **No GISAID data.** All accessions are INSDC
  accessions retrieved from NCBI. NCBI's GenBank overview and its website policies state that NCBI places no restrictions on
  use or distribution of GenBank data, but that submitters may claim patent, copyright or other IP rights that NCBI cannot
  assess: https://www.ncbi.nlm.nih.gov/genbank/about/, https://www.ncbi.nlm.nih.gov/home/about/policies/.
  The policies page also limits scripted access to at most 3 requests per second and asks for the `email` and `tool`
  parameters. The script obeys both (0.4 s minimum interval, `--email` argument, fixed `tool` string, no keys).
- **Why appropriate.** Homologous, closely related, equal-length sequences, so positional comparison is meaningful without
  alignment. The QC check is that the maximum pairwise distance is 57/1701 = 3.4%, while a misaligned pair would sit near 75%.
  The set is small enough for O(N³) algorithms, and it has real d_H = 1 pairs (8).
- **Selection procedure (reproducible).** `fetch_public_data.py select` runs an esearch (query in `selection_log.md`) and
  collects all 11,887 hits by paging. It sorts the accessions, takes 600 evenly spaced ones, and fetches their GenBank records.
  It then keeps those with serotype H3N2, a human host, a complete (non-partial) HA CDS, and a collection month inside the
  season, and extracts the CDS (171 remain). It keeps the modal length, 1701 (all 171), and picks 36 evenly spaced over
  accession order. `fetch_public_data.py refetch` re-downloads exactly the accessions in `accessions.csv`.
- **Known selection biases.** 31 of 36 are from the USA (this reflects deposits, not prevalence). Records whose collection
  date has no month are excluded. Even spacing over accession order is deterministic but not a random sample.
- **Preprocessing rules** (implemented in Stage 1, with every decision logged in the trace and shown in the UI):
  1. Uppercase.
  2. Drop any record containing a character outside {A,C,G,T}, and log its id and the offending characters
     (Dataset B: MH246583.1 {R:1}, MK181701.1 {R:1, Y:1}).
  3. Require one common length, otherwise reject with an error (no alignment in scope).
  4. Collapse identical sequences into haplotypes, keeping the member lists (Dataset B: 34 → 30 haplotypes).
  5. Reference: syn01 for A. For B, **[decision]** the medoid haplotype (minimum Σ d_H; ties broken by smallest accession).

## 5. Mathematical models and measurable outputs per unit
Notation: V = haplotypes, N = |V|, n = sequence length, d = d_H.

| Unit | Model | Measurable outputs (all with traces) |
|---|---|---|
| I | S: [n]→Σ; M(R,S), M*(R,S) **[decision]**; R_k; predicates; identity and orientation equivalences | mutation sets and cardinalities; ∩ ∪ \ Δ; d-identity check on all pairs (count of pairs where it holds = N(N−1)/2); relation-property report with counterexample triple; predicate truth tables; haplotype and orientation classes; optional Warshall closure = components **[decision]** |
| II | 4^n; M(n,k) = C(n,k)·3^k; recurrence; ball sizes; IE; pigeonhole; k! geodesics | tables for n ≤ 12; enumeration count = formula for all n ≤ 8 with no duplicates; guard refusal above the cap (10^5); 3n vs observed d = 1 neighbours per haplotype |
| III | I, R, C, RC on Σ^n; G2 ≅ Z2; V4 ≅ Z2×Z2; action, orbits, Burnside; point-substitution monoid | Cayley tables; axiom report (closure, associativity over all triples, identity, inverses); exhaustive function equality n ≤ 6; Fix and orbit counts = brute force n ≤ 8; dataset orbits; Dataset B canonical k-mer counts **[decision]** |
| IV | enc/dec; comp identity; canonical k-mer; b_m; two hashing theorems; Euclid, extended Euclid, inverses; CRT; Fermat order (pedagogical) | dec∘enc = id exhaustive n ≤ 6; bucket histograms and collision counts for m ∈ {16, 64, 7, 13, 17, 31, 61} (proposed set); theorem checks on every d = 1 pair of A and B; CRT round-trip for all 64 3-mers; Euclid traces |
| V | G_=1, G_k (sparse, weight d); matrices; BFS/DFS; components; κ, λ (small only); Dijkstra O(V²) | |E| and #components for k = 0..max; degree table and handshake check; component membership; Dijkstra cost vs direct d with all optimal paths; bottleneck k* vs max MST edge |
| VI | K_N with weight d; Prim O(V²); Kruskal with union-find; rooted tree (root = reference, arbitrary); UPGMA | MST weight (Prim = Kruskal = NetworkX); validity (N−1 edges, acyclic, spanning); rooted-tree depth, height, leaves, internal nodes; pre/post-order (MST), in-order (UPGMA); UPGMA merge table and ultrametric check; MST-vs-UPGMA comparison |

Complexities to report (§10) are listed per unit in `math_audit.md`. The exponential pieces are recursive generation
(Θ(n·M(n,k)), guarded), exhaustive axiom checks (4^n, n ≤ 6), Burnside brute force (n ≤ 8) and κ/λ (subset removal, small graphs only).

## 6. Experiments
| Exp | Question | Data | Expected values from | Stage |
|---|---|---|---|---|
| 1 Mutation-set analysis | do sets/ops/distances match hand values; the d-identity on all pairs | A (+ B for the identity) | `synthetic_expected.md` §2–3 | 1 |
| 2 Counting vs enumeration | enumeration = C(n,k)·3^k; recurrence = closed form | generated, n ≤ 8 | closed form (§8 of expected) | 2 |
| 3 Algebraic verification | group axioms, Cayley tables, Burnside vs brute force | Σ^n, n ≤ 8; A | §6 | 3 |
| 4 Modular indexing | collisions vs modulus; both hashing theorems; canonical k-mers **[decision]** | A, B | §7 (A); theorem statements (B) | 4 |
| 5 Network structure | edges/components vs k; handshake; bottleneck | A, B | §4 (A); QC script's edge counts (B) | 5 |
| 6 Shortest paths and spanning trees | Dijkstra, Prim, Kruskal vs hand values and NetworkX | A, B | §4–5 (A); NetworkX at test time (B) | 5 (Dijkstra), 6 (MST) |
| 7 Tree analysis | MST vs UPGMA: structure, assumptions, interpretation | A, B | §5 (MST); UPGMA after tie rule **[decision]** | 6 |
| 8 Real-data evaluation | whole validated pipeline on B; observed patterns without unsupported claims | B | QC invariants: 34 clean, 30 haplotypes, 8 edges at d = 1, max d = 57 | 7 |

## 7. Stage plan (stage ↔ spec phase)
| Stage | Content | Spec phase |
|---|---|---|
| 0 | audit, final spec, skeleton, Dataset A, Dataset B feasibility (this stage) | 1, 2 (+ data check) |
| 1 | `sequence_processing`, `set_logic`; trace type; Exp 1 | 3–4 |
| 2 | `counting_recursion`; Exp 2 | 3–4 |
| 3 | `algebraic_structures`; Exp 3 | 3–4 |
| 4 | `number_theory`; Exp 4 | 3–4 |
| 5 | `mutation_graph`; Exp 5, Exp 6 (Dijkstra) | 3–4 |
| 6 | `tree_analysis`, `phylogenetic_analysis`; Exp 6 (MST), Exp 7 | 3–4 |
| 7 | real-data run on Dataset B; Exp 8 | 5 |
| 8 | Streamlit app in `app/` | 6 |
| 9 | proposal, design document, report, slides, viva questions | 7 |

Each stage ends with a report listing what exists, the raw test output, the commands to re-run, and open issues. Then I stop for your go/no-go.

## 8. Validation plan (spec §10)
Per module: hand-value tests from `synthetic_expected.md`; invalid input (non-ACGT, unequal lengths, empty list, single
sequence); repeated sequences; counting boundaries (k = 0, k = n, k > n, n = 0); number-theory boundaries (gcd(0, a),
non-invertible inverses, non-coprime CRT moduli); a disconnected graph (syn09); multiple optimal solutions (8 MSTs, 3 equal
Dijkstra paths); and NetworkX/SymPy oracles on Dataset B.
