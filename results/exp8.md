# Experiment 8: Real-data evaluation of the full pipeline

- run date: 2026-10-31
- command: `python -m experiments.exp8_real_data`
- Python 3.13.11 on Windows

## Question
Does the pipeline validated on Dataset A run end to end on real genomic data, reproducing the independently known quality-control figures and the hand-proved theorems, and what does it then observe about the 30 influenza A(H3N2) HA haplotypes?

## Input
`data/public_genomic/h3n2_ha_cds.fasta`: 36 GenBank hemagglutinin coding sequences, human host, 2017–18 northern-hemisphere season, each 1701 nt. Accessions in `data/public_genomic/accessions.csv`.

## Mathematical model
All six units in one pass: mutation sets and the threshold relation R_k (I); sequence-space and Hamming-ball counts and the pigeonhole bound (II); the {I,RC} action and the base-level Klein group {id, (A T)(C G), (A G)(C T), (A C)(G T)} (III, decision D9); base-4 encoding and modular bucketing (IV); threshold graphs, components and Dijkstra (V); Prim, Kruskal and UPGMA (VI).

## Procedure
Load and preprocess the FASTA (uppercase, drop non-ACGT, require one length, collapse haplotypes, choose the medoid reference). Then run each unit's algorithms in turn. Independent entries are compared with the genomath-free QC script or with hand derivations; entries that only check genomath against itself are tagged as consistency checks; everything else is reported as an observation.

## Expected
Values fixed before the run; never produced by `genomath`.

| # | item | expected | source |
|---|---|---|---|
| 1 | records read from FASTA | 36 | docs/dataset_b_qc.md (independent script) |
| 2 | records dropped for non-ACGT characters | [MH246583.1, MK181701.1] | docs/dataset_b_qc.md (independent script) |
| 3 | clean records / haplotypes / common length | [34, 30, 1701] | docs/dataset_b_qc.md (independent script) |
| 4 | groups of identical sequences | [[CY249987.1, MH082044.1, MH601290.1], [MH082676.1, MH294810.1, MH358590.1]] | docs/dataset_b_qc.md (independent script) |
| 5 | Unit I: haplotype pairs; pairs at d_H = 1 | [435, 8] | docs/dataset_b_qc.md (independent script) |
| 6 | Unit I: min and max pairwise d_H | {min: 1, max: 57} | docs/dataset_b_qc.md (independent script) |
| 7 | Unit I: d_H identity holds on every pair | 435 | math_audit.md Unit I (theorem) — *consistency check* |
| 8 | Unit I: R_1 properties and non-transitivity witness | *(none: reported as observation)* |  |
| 9 | Unit II: decimal digits of 4^n | 1025 | docs/unit2_expected.md §1 |
| 10 | Unit II: \|B(n,1)\| and 3n possible neighbours | [5104, 5103] | docs/unit2_expected.md §2 |
| 11 | Unit II: largest forced repeated k-mer length | 5 | docs/unit2_expected.md §5 |
| 12 | Unit II: clean sequences containing a repeated 5-mer | 34 | docs/unit2_expected.md §5 |
| 13 | Unit III: RC pairs present / RC-fixed / {I,RC} classes | [0, 0, 30] | docs/unit3_expected.md §8 |
| 14 | Unit III: canonical 4-mer classes in the whole space (Burnside) | 136 | docs/unit4_expected.md §2 |
| 15 | D9: fixed-point-free involutions of Σ | 3 | docs/unit3_expected.md §9 |
| 16 | D9: base-level group is an abelian group of order 4, every non-identity element an involution | [True, True, [1, 2, 2, 2]] | docs/unit3_expected.md §9 |
| 17 | D9: ordered base pairs covered, each by exactly one element | 12 | docs/unit3_expected.md §9 |
| 18 | D9: transitions / transversions among the 6 unordered pairs | [2, 4] | docs/unit3_expected.md §9 |
| 19 | D9: census by element sums to the substitution count | True | definition — *consistency check* |
| 20 | D9: combinatorial baseline ratio | 0.5 | docs/unit3_expected.md §9 |
| 21 | Unit IV: decode(encode(S)) = S for every haplotype | True | bijection (unit4 §1) — *consistency check* |
| 22 | Unit IV: smallest j with m \| 4^j (bases the bucket can see) | {16: 2, 64: 3, 7: None, 13: None, 17: None, 31: None, 61: None} | docs/unit4_expected.md §3 |
| 23 | Unit IV: d_H = 1 pairs sharing a bucket, by modulus | {16: 8, 64: 8, 7: 0, 13: 0, 17: 0, 31: 0, 61: 0} | docs/unit4_expected.md §3 |
| 24 | Unit IV: max bucket load ≥ pigeonhole bound for every modulus | True | pigeonhole — *consistency check* |
| 25 | Unit V: edges at k = 1 equal the QC count of d_H = 1 pairs | 8 | docs/dataset_b_qc.md (independent script) |
| 26 | Unit V: components never increase as k grows | True | G_j ⊆ G_k — *consistency check* |
| 27 | Unit V: path cost ≥ direct d_H for all 900 ordered pairs at k* = 47 | True | triangle inequality — *consistency check* |
| 28 | Unit VI: Prim total = Kruskal total | True | MST weight is unique — *consistency check* |
| 29 | Unit VI: both outputs are valid spanning trees | True | definition — *consistency check* |
| 30 | Unit VI: largest MST edge = k* (bottleneck property) | True | docs/unit6_expected.md §3 — *consistency check* |
| 31 | Unit VI: UPGMA cophenetic distances are ultrametric | True | UPGMA property — *consistency check* |
| 32 | D9 × Unit VI: substitutions along MST edges = MST weight | 171 | Σ d_H over tree edges — *consistency check* |
| 33 | D9: differences from the reference (29 pairs) | *(none: reported as observation)* |  |
| 34 | D9: substitutions along MST edges (29 edges) | *(none: reported as observation)* |  |
| 35 | observed: reference (medoid) | *(none: reported as observation)* |  |
| 36 | observed: mutation-set sizes vs the reference | *(none: reported as observation)* |  |
| 37 | observed: observed vs possible single-substitution neighbours | *(none: reported as observation)* |  |
| 38 | observed: distinct / canonical / RC-palindromic 4-mers | *(none: reported as observation)* |  |
| 39 | observed: bucket occupancy per modulus | *(none: reported as observation)* |  |
| 40 | observed: G_=1 structure | *(none: reported as observation)* |  |
| 41 | observed: smallest connecting threshold k* | *(none: reported as observation)* |  |
| 42 | observed: detours from the reference at k* | *(none: reported as observation)* |  |
| 43 | observed: MST weight / number of MSTs / largest edge | *(none: reported as observation)* |  |
| 44 | observed: UPGMA fit and ties | *(none: reported as observation)* |  |
| 45 | observed: MST vs UPGMA bipartitions | *(none: reported as observation)* |  |

## Actual
Values produced by this run of `genomath`.

| # | item | actual |
|---|---|---|
| 1 | records read from FASTA | 36 |
| 2 | records dropped for non-ACGT characters | [MH246583.1, MK181701.1] |
| 3 | clean records / haplotypes / common length | [34, 30, 1701] |
| 4 | groups of identical sequences | [[CY249987.1, MH082044.1, MH601290.1], [MH082676.1, MH294810.1, MH358590.1]] |
| 5 | Unit I: haplotype pairs; pairs at d_H = 1 | [435, 8] |
| 6 | Unit I: min and max pairwise d_H | {min: 1, max: 57} |
| 7 | Unit I: d_H identity holds on every pair | 435 |
| 8 | Unit I: R_1 properties and non-transitivity witness | {reflexive: True, symmetric: True, transitive: False, witness: (CY249987.1, MH637643.1, MH082676.1)} |
| 9 | Unit II: decimal digits of 4^n | 1025 |
| 10 | Unit II: \|B(n,1)\| and 3n possible neighbours | [5104, 5103] |
| 11 | Unit II: largest forced repeated k-mer length | 5 |
| 12 | Unit II: clean sequences containing a repeated 5-mer | 34 |
| 13 | Unit III: RC pairs present / RC-fixed / {I,RC} classes | [0, 0, 30] |
| 14 | Unit III: canonical 4-mer classes in the whole space (Burnside) | 136 |
| 15 | D9: fixed-point-free involutions of Σ | 3 |
| 16 | D9: base-level group is an abelian group of order 4, every non-identity element an involution | [True, True, [1, 2, 2, 2]] |
| 17 | D9: ordered base pairs covered, each by exactly one element | 12 |
| 18 | D9: transitions / transversions among the 6 unordered pairs | [2, 4] |
| 19 | D9: census by element sums to the substitution count | True |
| 20 | D9: combinatorial baseline ratio | 0.5 |
| 21 | Unit IV: decode(encode(S)) = S for every haplotype | True |
| 22 | Unit IV: smallest j with m \| 4^j (bases the bucket can see) | {16: 2, 64: 3, 7: None, 13: None, 17: None, 31: None, 61: None} |
| 23 | Unit IV: d_H = 1 pairs sharing a bucket, by modulus | {16: 8, 64: 8, 7: 0, 13: 0, 17: 0, 31: 0, 61: 0} |
| 24 | Unit IV: max bucket load ≥ pigeonhole bound for every modulus | True |
| 25 | Unit V: edges at k = 1 equal the QC count of d_H = 1 pairs | 8 |
| 26 | Unit V: components never increase as k grows | True |
| 27 | Unit V: path cost ≥ direct d_H for all 900 ordered pairs at k* = 47 | True |
| 28 | Unit VI: Prim total = Kruskal total | True |
| 29 | Unit VI: both outputs are valid spanning trees | True |
| 30 | Unit VI: largest MST edge = k* (bottleneck property) | True |
| 31 | Unit VI: UPGMA cophenetic distances are ultrametric | True |
| 32 | D9 × Unit VI: substitutions along MST edges = MST weight | 171 |
| 33 | D9: differences from the reference (29 pairs) | {substitutions: 341, transitions: 250, transversions: 91, by_base_pair: {AC: 46, AG: 154, AT: 31, CT: 96, GT: 14}} |
| 34 | D9: substitutions along MST edges (29 edges) | {substitutions: 171, transitions: 131, transversions: 40, by_base_pair: {AC: 20, AG: 77, AT: 14, CT: 54, GT: 6}} |
| 35 | observed: reference (medoid) | CY249987.1 |
| 36 | observed: mutation-set sizes vs the reference | {min: 0, median: 4, max: 49} |
| 37 | observed: observed vs possible single-substitution neighbours | {max_observed_neighbours: 6, possible: 5103, fraction: 1.18e-03} |
| 38 | observed: distinct / canonical / RC-palindromic 4-mers | {distinct: 248, canonical: 136, palindromic: 16} |
| 39 | observed: bucket occupancy per modulus | {16: {occupied: 1, max_load: 30, bound: 2}, 64: {occupied: 1, max_load: 30, bound: 1}, 7: {occupied: 7, max_load: 8, bound: 5}, 13: {occupied: 11, max_load: 8, bound: 3}, 17: {occupied: 15, max_load: 3, bound: 2}, 31: {occupied: 20, max_load: 4, bound: 1}, 61: {occupied: 24, max_load: 2, bound: 1}} |
| 40 | observed: G_=1 structure | {edges: 8, components: 22, cyclomatic: 0, max_degree: 6, is_forest: True} |
| 41 | observed: smallest connecting threshold k* | 47 |
| 42 | observed: detours from the reference at k* | {MH885199.1: 2} |
| 43 | observed: MST weight / number of MSTs / largest edge | {weight: 171, count: 4, max_edge: 47} |
| 44 | observed: UPGMA fit and ties | {ultrametric_input: False, exact_fit: 9, pairs: 435, max_deviation: 73/14, tie_steps: 11} |
| 45 | observed: MST vs UPGMA bipartitions | {mst: 10, upgma: 27, shared: 5, symmetric_difference: 27} |

## Comparison

| # | item | match |
|---|---|---|
| 1 | records read from FASTA | yes |
| 2 | records dropped for non-ACGT characters | yes |
| 3 | clean records / haplotypes / common length | yes |
| 4 | groups of identical sequences | yes |
| 5 | Unit I: haplotype pairs; pairs at d_H = 1 | yes |
| 6 | Unit I: min and max pairwise d_H | yes |
| 7 | Unit I: d_H identity holds on every pair | yes (consistency) |
| 8 | Unit I: R_1 properties and non-transitivity witness | n/a |
| 9 | Unit II: decimal digits of 4^n | yes |
| 10 | Unit II: \|B(n,1)\| and 3n possible neighbours | yes |
| 11 | Unit II: largest forced repeated k-mer length | yes |
| 12 | Unit II: clean sequences containing a repeated 5-mer | yes |
| 13 | Unit III: RC pairs present / RC-fixed / {I,RC} classes | yes |
| 14 | Unit III: canonical 4-mer classes in the whole space (Burnside) | yes |
| 15 | D9: fixed-point-free involutions of Σ | yes |
| 16 | D9: base-level group is an abelian group of order 4, every non-identity element an involution | yes |
| 17 | D9: ordered base pairs covered, each by exactly one element | yes |
| 18 | D9: transitions / transversions among the 6 unordered pairs | yes |
| 19 | D9: census by element sums to the substitution count | yes (consistency) |
| 20 | D9: combinatorial baseline ratio | yes |
| 21 | Unit IV: decode(encode(S)) = S for every haplotype | yes (consistency) |
| 22 | Unit IV: smallest j with m \| 4^j (bases the bucket can see) | yes |
| 23 | Unit IV: d_H = 1 pairs sharing a bucket, by modulus | yes |
| 24 | Unit IV: max bucket load ≥ pigeonhole bound for every modulus | yes (consistency) |
| 25 | Unit V: edges at k = 1 equal the QC count of d_H = 1 pairs | yes |
| 26 | Unit V: components never increase as k grows | yes (consistency) |
| 27 | Unit V: path cost ≥ direct d_H for all 900 ordered pairs at k* = 47 | yes (consistency) |
| 28 | Unit VI: Prim total = Kruskal total | yes (consistency) |
| 29 | Unit VI: both outputs are valid spanning trees | yes (consistency) |
| 30 | Unit VI: largest MST edge = k* (bottleneck property) | yes (consistency) |
| 31 | Unit VI: UPGMA cophenetic distances are ultrametric | yes (consistency) |
| 32 | D9 × Unit VI: substitutions along MST edges = MST weight | yes (consistency) |
| 33 | D9: differences from the reference (29 pairs) | n/a |
| 34 | D9: substitutions along MST edges (29 edges) | n/a |
| 35 | observed: reference (medoid) | n/a |
| 36 | observed: mutation-set sizes vs the reference | n/a |
| 37 | observed: observed vs possible single-substitution neighbours | n/a |
| 38 | observed: distinct / canonical / RC-palindromic 4-mers | n/a |
| 39 | observed: bucket occupancy per modulus | n/a |
| 40 | observed: G_=1 structure | n/a |
| 41 | observed: smallest connecting threshold k* | n/a |
| 42 | observed: detours from the reference at k* | n/a |
| 43 | observed: MST weight / number of MSTs / largest edge | n/a |
| 44 | observed: UPGMA fit and ties | n/a |
| 45 | observed: MST vs UPGMA bipartitions | n/a |

**20/20 independent expected values matched; 11/11 consistency checks passed; 14 observation-only items.**

Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees with itself.

### Substitutions by unordered base pair (observation, direction-free)

| view | substitutions | transitions | transversions | ratio | by pair |
|---|---|---|---|---|---|
| vs reference | 341 | 250 | 91 | 2.747 | {'AC': 46, 'AG': 154, 'AT': 31, 'CT': 96, 'GT': 14} |
| along MST edges | 171 | 131 | 40 | 3.275 | {'AC': 20, 'AG': 77, 'AT': 14, 'CT': 54, 'GT': 6} |

### Threshold sweep (observation)

| k | edges | components | largest | isolated | max degree |
|---|---|---|---|---|---|
| 0 | 0 | 30 | 1 | 30 | 0 |
| 1 | 8 | 22 | 9 | 21 | 6 |
| 2 | 31 | 16 | 13 | 13 | 9 |
| 3 | 57 | 12 | 15 | 8 | 14 |
| 4 | 90 | 10 | 17 | 6 | 15 |
| 5 | 113 | 8 | 19 | 4 | 18 |
| 6 | 148 | 6 | 21 | 2 | 20 |
| 7 | 178 | 6 | 21 | 2 | 20 |
| 8 | 202 | 5 | 22 | 1 | 21 |
| 9 | 223 | 5 | 22 | 1 | 21 |
| 10 | 227 | 5 | 22 | 1 | 21 |
| 15 | 238 | 4 | 22 | 0 | 21 |
| 20 | 244 | 3 | 22 | 0 | 21 |
| 30 | 284 | 2 | 28 | 0 | 25 |
| 40 | 379 | 2 | 28 | 0 | 27 |
| 47 | 380 | 1 | 30 | 0 | 28 |
| 57 | 435 | 1 | 30 | 0 | 29 |

## Interpretation
The pipeline runs end to end on real data and reproduces every independently known figure: the QC counts (36 read, 2 dropped for ambiguity codes, 34 clean, 30 haplotypes, 1701 nt, 8 pairs at distance 1, distances 1 to 57), the Unit II hand values (4^1701 has 1025 digits; 5103 possible single-substitution neighbours; a repeated 5-mer is forced in all 34 sequences), the Unit III prediction that no sequence's reverse complement is present, and both Unit IV hashing theorems (m = 16 and 64 collide on all 8 distance-1 pairs; every odd modulus separates all 8).

Observations about this sample, with no biological claim:
- **Sparsity.** The best-connected haplotype has 6 observed single-substitution neighbours out of 5103 possible (1.18e-03), so the single-substitution graph is an extremely thin sample of the Hamming graph.
- **Structure.** G_=1 has 8 edges in 22 components and is a forest; the graph first becomes connected at k* = 47, which equals the largest edge of the minimum spanning tree (weight 171, 4 distinct MSTs).
- **Tree comparison.** The distances are not ultrametric, so UPGMA reproduces only 9 of 435 distances exactly; 11 merge steps involved ties, and the MST and UPGMA share 5 of their 10 and 27 non-trivial bipartitions.
- **Substitution classes (D9).** Direction-free counts give 250 transitions and 91 transversions among differences from the reference (ratio 2.747), and 131 to 40 along MST edges (ratio 3.275). Both are far from the combinatorial baseline of 0.5, which is what equal likelihood over the 6 unordered base pairs would give. The counts are what they are; **explaining** the gap would need a substitution model and a source for what ratio real data shows, and this project has neither, so no conclusion is drawn from it.

## Limitations
All six units' limitations apply at once and are listed in Experiments 1–7. Specific to this run: the 30 haplotypes come from one season and are dominated by United States deposits, so they are a convenience sample and not a random sample of anything; positional comparison rests on equal CDS length and low divergence, not on an alignment; two records with ambiguity codes were dropped rather than resolved; substitution counts are direction-free, because no ancestor is known, and counting along MST edges counts one chosen tree's steps, not inferred events; the reference is the medoid, a mathematical choice with no biological standing; neither the graph nor either tree is evidence of evolutionary history; and nothing here supports any clinical, diagnostic or epidemiological claim.
