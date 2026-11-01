MASTER PROJECT PROMPT — DISCRETE MATHEMATICS COURSE PROJECT
Project Title
Genomic Mutation & Evolution Analyzer: An Integrated Computational Framework Using Discrete Mathematics
1. Your role
Act as a multidisciplinary team consisting of:

* A university professor specializing in Discrete Mathematics.
* A bioinformatics researcher specializing in genomic sequences, mutations and phylogenetic analysis.
* A computer science researcher specializing in algorithms, graph theory and computational modeling.
* A software architect responsible for designing a working, well-tested application.
* A strict academic evaluator who will assess whether the project genuinely uses mathematics rather than merely displaying mathematical formulas.

Your task is to help us design, implement, validate and document a complete undergraduate course project.
Do not simply agree with the proposed design. Critically evaluate it, identify weak mathematical connections, reject unjustified biological claims and improve the project wherever necessary.
2. Academic context
This is a second-year B.Tech Computer Science and Engineering (Artificial Intelligence) course project for Vishwakarma Institute of Technology, Pune.
Course: MM0903 — Discrete Mathematics
Academic year: 2026–27
The project must meaningfully incorporate all six units of the prescribed syllabus:

1. Unit I — Set Theory and Logic
2. Unit II — Recursion and Counting Methods
3. Unit III — Algebraic Structures
4. Unit IV — Number Theory
5. Unit V — Graph Theory
6. Unit VI — Trees

The primary objective is to demonstrate the practical application of discrete mathematical concepts to a real computational problem.
The project is not primarily an AI/ML project, a generic web application or a conventional bioinformatics tool. Its central contribution must be the mathematical modeling, algorithms, derivations, computational experiments and interpretation of results.
The software interface should be used to demonstrate the mathematics, not replace it.
3. Project overview
We propose a project called Genomic Mutation & Evolution Analyzer.
DNA sequences consist of discrete symbols drawn from a finite alphabet:
[
\Sigma = {A,C,G,T}
]
A DNA sequence can therefore be modeled as a finite string over this alphabet. Different sequences may differ at particular positions, share mutations, form similarity relationships and belong to connected groups of related sequences.
The proposed system will analyze genomic sequences through mathematical operations and construct useful representations of sequence differences, mutation relationships and hierarchical sequence groupings.
The system should be able to:

* Accept a collection of DNA sequences.
* Validate and represent the sequences mathematically.
* Identify substitutions between comparable sequences.
* Represent mutations as sets and relations.
* Calculate possible sequence configurations and mutation counts.
* Explore mathematically valid transformations and algebraic structures.
* Investigate numerical encodings and modular representations of sequences.
* Construct mutation-similarity graphs.
* Analyze graph connectivity, paths, clusters and other structural properties.
* Build tree representations of sequence relationships.
* Demonstrate relevant algorithms from the Discrete Mathematics syllabus.
* Display calculations and explain how each result was obtained.

The system must support both a small synthetic dataset for manually verifiable experiments and a publicly sourced genomic dataset for realistic evaluation.
Central research question
How can discrete mathematical structures be integrated into a computational framework to represent genomic sequences, analyze mutation relationships, and construct interpretable graph and tree representations?
Possible supporting questions include:

1. How can set theory, functions and relations represent genomic sequences and their differences?
2. How can counting and recurrence relations describe the mutation space?
3. How can algebraic transformations represent valid sequence symmetries?
4. What useful computational properties can be demonstrated through numerical encoding and modular arithmetic?
5. How can graph theory identify relationships and paths among observed sequences?
6. How can tree structures represent hierarchical relationships, and how do minimum spanning trees differ from evolutionary trees?

4. Critical mathematical requirement
Before implementation, conduct a mathematical feasibility audit.
For every unit, explain:

* The real problem being addressed.
* The mathematical objects being modeled.
* The relevant definitions and equations.
* The algorithm that implements the mathematics.
* The output produced.
* Why the mathematics is appropriate for the problem.
* What limitations or assumptions apply.

Do not claim that a mathematical concept has biological significance merely because it can be coded.
In particular:

* Do not pretend that arbitrary DNA mutations form a group.
* Do not describe modular indexing as an established clinical genomics technique without evidence.
* Do not equate a minimum spanning tree with a scientifically validated phylogenetic tree.
* Do not interpret a computational pattern as proof of evolutionary ancestry.
* Do not claim that the system diagnoses disease, predicts pathogenicity or discovers clinically meaningful mutations.
* Explicitly identify any mathematical demonstrations that are primarily educational rather than biologically necessary.

All six units must have a visible implementation or computational experiment. However, only include individual syllabus subtopics where their relevance can be explained honestly. If any unit is weakly integrated, explain why and propose better alternatives.
5. Mathematical design: all six units
UNIT I — Set Theory and Logic
The system should represent DNA sequences and their relationships through sets, functions, predicates and relations.
A. DNA sequence representation
Define:
[
\Sigma={A,C,G,T}
]
and represent a fixed-length sequence as:
[
S\in\Sigma^n
]
Explain how positions, nucleotides, sequences and sequence collections can be represented mathematically.
B. Mutation sets
Given a reference sequence (R) and a comparable sequence (S), define the mutation-position set:
[
M(R,S)={i\ne S_i}
]
Implement:

* Identification of mutation positions.
* Mutation-set cardinality.
* Intersection of mutation sets.
* Union of mutation sets.
* Set difference.
* Comparison of shared and sequence-specific mutations.

For example, if:
[
M(R,S_1)={3,6,10}
]
and
[
M(R,S_2)={6,10,15}
]
then the shared mutations are:
[
M(R,S_1)\cap M(R,S_2)={6,10}
]
Explain how these operations support sequence comparison.
C. Relations
Define a relation between sequences using a specified mutation-distance threshold:
[
S_iR_kS_j
\iff d_H(S_i,S_j)\le k
]
where (d_H) is Hamming distance for equal-length sequences.
Test and explain reflexivity, symmetry and transitivity separately. Do not incorrectly claim that every threshold-based similarity relation is transitive.
D. Logic
Represent sequence properties as predicates, such as:

* Sequence differs from the reference at no more than (k) positions.
* Sequence contains a specified nucleotide at a specified position.
* Sequence shares a specified mutation with another sequence.

Implement Boolean combinations of such predicates to filter sequences.
E. Functions and equivalence
Demonstrate functions from sequences to numerical encodings and define an appropriate equivalence relation for sequence orientation where applicable.
Explain equivalence classes and distinguish ordinary similarity from mathematical equivalence.
Required outputs: mutation sets, set-operation results, predicates, relation properties, sequence filters and mathematical explanations.
UNIT II — Recursion and Counting Methods
Use counting theory to describe the size of the sequence space and the number of possible mutations.
A. Sequence-space size
For DNA sequences of length (n), calculate:
[
N(n)=4^n
]
Explain how the search space grows as the sequence length increases.
B. Mutation-count formula
For fixed-length sequences where only substitutions are permitted, derive and implement:
[
M(n,k)=\binom nk 3^k
]
where:

* (n) is the sequence length.
* (k) is the exact number of substitutions.
* (\binom nk) selects the mutated positions.
* Each selected position has three alternative nucleotides.

Clearly state that this formula does not count insertions and deletions.
C. Recurrence relation
Derive and implement:
[
M(n,k)=M(n-1,k)+3M(n-1,k-1)
]
Explain the two cases responsible for this recurrence.
Include appropriate base cases, such as:
[
M(0,0)=1
]
and define values outside (0\le k\le n) as zero.
Verify the recurrence against the closed-form expression for small values.
D. Recursive mutation generation
Implement a recursive algorithm that generates sequences with a specified number of substitutions or within a specified mutation distance.
Explain its time and space complexity.
Prevent uncontrolled exponential generation for large inputs. Provide warnings and sensible input limits.
E. Additional counting concepts
Evaluate meaningful demonstrations of:

* Permutations and combinations.
* Binomial coefficients.
* Inclusion–Exclusion Principle.
* Pigeonhole Principle.
* Mathematical induction.

For each selected topic, provide a genomic-sequence-related example and derive the result instead of merely presenting a definition.
Required outputs: formulas, derivations, recurrence tables, recursive generation, complexity analysis and verification examples.
UNIT III — Algebraic Structures
This unit requires particularly careful mathematical justification.
The objective is to identify valid transformations of DNA sequences that genuinely satisfy algebraic properties.
A. Reverse-complement transformation
Define the DNA complement mapping:
[
A\leftrightarrow T,\qquad C\leftrightarrow G
]
Define the reverse-complement transformation:
[
RC(S)=\operatorname{reverse}(\operatorname{complement}(S))
]
Demonstrate the property:
[
RC(RC(S))=S
]
for valid DNA sequences.
B. Group structure
Define the set of transformations:
[
G={I,RC}
]
where (I) is the identity transformation, with the operation being composition.
Construct and explain the Cayley table.
Verify:

* Closure.
* Associativity.
* Identity.
* Inverses.

Determine the group and explain why the construction is mathematically valid.
Discuss the biological interpretation carefully: reverse-complement equivalence concerns DNA strand orientation and representation, not arbitrary mutation events.
C. Equivalence classes
Show how the transformations can induce equivalence classes of sequences.
Handle sequences fixed by the reverse-complement transformation as a special case.
D. Further algebraic exploration
Evaluate whether permutation groups, monoids, semigroups or other transformations provide a genuine extension.
If arbitrary position permutations are included, describe them as mathematical transformations and do not imply that they represent ordinary biological mutations.
Required outputs: transformation implementation, composition table, Cayley table, algebraic-property verification, equivalence classes and a clear statement of biological limitations.
UNIT IV — Number Theory
Develop a subsystem called Genomic Encoding and Modular Indexing Laboratory.
This subsystem should explore the mathematical consequences of representing fixed-length DNA words as integers.
It must be identified honestly as a computational representation experiment unless an independently validated biological application is established.
A. Base-4 encoding
Use the mapping:
[
A=0,\quad C=1,\quad G=2,\quad T=3
]
For a sequence of length (n), calculate its base-4 integer representation.
For example:
[
ACG\rightarrow 0\cdot4^2+1\cdot4+2=6
]
Demonstrate that fixed-length sequences have unique integer encodings under this mapping.
Implement both encoding and decoding.
B. Modular arithmetic
Calculate:
[
b(x)=x\bmod m
]
for a chosen modulus (m).
Use this to investigate indexing buckets for fixed-length sequence words or (k)-mers.
Measure collisions for different moduli and explain why sharing a bucket does not imply that two sequences are similar.
C. Number-theoretic algorithms
Implement relevant syllabus topics, including:

* Division algorithm.
* Divisibility.
* GCD.
* Euclidean algorithm.
* Extended Euclidean algorithm.
* Congruences.
* Modular inverses where defined.

Explain how each algorithm works and demonstrate its mathematical correctness on suitable examples.
D. Chinese Remainder Theorem
Provide an experiment showing how an encoded integer can be reconstructed from its residues under pairwise coprime moduli, subject to the appropriate conditions.
Explain the uniqueness range and verify the reconstruction.
E. Euler and Fermat results
Determine whether Euler's theorem, Euler's phi function or Fermat's Little Theorem can support a meaningful additional integer-representation experiment.
If included primarily to demonstrate the syllabus, label that experiment as pedagogical rather than claiming it is necessary for genomic analysis.
Required outputs: sequence encoder/decoder, modular-indexing experiments, collision analysis, number-theoretic calculations, CRT reconstruction and explicit limitations.
UNIT V — Graph Theory
Build the project's principal analytical component: the Mutation Relationship Graph.
A. Graph construction
Represent sequences as vertices:
[
G=(V,E)
]
Each vertex represents an observed sequence.
Provide at least two graph modes:

1. A single-substitution graph, where an edge connects sequences whose Hamming distance is exactly one.
2. A threshold-based similarity graph, where an edge connects sequences whose distance falls within a defined threshold.

Clearly state the edge rule used in each mode.
B. Mutation-distance calculation
Implement Hamming distance for equal-length sequences:
[
d_H(S_i,S_j)=\sum_{p=1}^{n}[S_{i,p}\ne S_{j,p}]
]
If the dataset includes unequal-length or unaligned sequences, do not apply Hamming distance directly. Align or preprocess comparable sequences appropriately before computing distances.
C. Graph analysis
Implement and demonstrate:

* Vertex degree.
* Paths and cycles.
* Connected components.
* Connectivity.
* Adjacency matrix.
* Incidence matrix.
* Graph traversal using BFS and DFS.
* Relevant graph operations.
* Vertex and edge connectivity where feasible.
* Shortest paths.
* Dijkstra's algorithm on a suitable non-negatively weighted graph.

Use a carefully chosen graph construction and weighting scheme so that shortest-path calculations are meaningful and not merely redundant with direct mutation-distance calculations.
D. Mutation clusters
Investigate how changing the distance threshold changes graph connectivity and cluster structure.
Display the results across several thresholds and explain the mathematical consequences.
E. Graph validation
Use small hand-verified examples to test graph construction, degrees, components, paths and shortest-path results.
Required outputs: interactive or static mutation graphs, matrices, connectivity reports, path results, cluster summaries and algorithmic complexity analysis.
UNIT VI — Trees
Use trees to represent hierarchical structures and investigate connections between graph theory and sequence relationships.
Keep distinct tree concepts separate.
A. Minimum spanning tree
Given a suitably defined weighted sequence graph, implement:

* Prim's algorithm.
* Kruskal's algorithm.

Compare the total edge weights and verify that the resulting spanning structures connect the appropriate vertices without cycles.
Explain what an MST means mathematically.
Do not automatically call an MST a phylogenetic tree. An MST minimizes total graph-edge weight; it does not, by itself, establish evolutionary ancestry.
B. Hierarchical sequence representation
Implement tree properties and operations, including:

* Rooted trees.
* Parent and child relationships.
* Internal nodes and leaves.
* Depth and height.
* Preorder traversal.
* Inorder traversal for suitable binary-tree representations.
* Postorder traversal.
* Tree construction and validation.

C. Evolutionary tree experiment
If the selected dataset supports it, construct a distance-based sequence-clustering or phylogenetic tree using an appropriate method such as UPGMA or neighbor joining.
Document the method's assumptions, its input distance matrix and the interpretation of branch relationships.
Distinguish the resulting evolutionary hypothesis from the MST produced by Prim's or Kruskal's algorithm.
Do not claim that an illustrative tree proves real evolutionary history.
Required outputs: tree visualizations, traversals, tree statistics, MST comparisons and a clear distinction between graph-theoretic and phylogenetic interpretations.
6. Data requirements
Use two kinds of data.
Dataset A — Synthetic validation dataset
Create a small collection of manually understandable DNA sequences.
Use it to verify:

* Sequence comparisons.
* Mutation positions.
* Set operations.
* Counting results.
* Algebraic transformations.
* Number-theoretic encoding.
* Graph structures.
* Tree algorithms.

For this dataset, manually verify several examples so the professor can reproduce the results.
Dataset B — Public genomic dataset
Identify an openly accessible dataset containing comparable genomic sequences.
Possible sources include:

* NCBI.
* Ensembl.
* EMBL-EBI.

Choose a dataset of homologous sequences or a collection of closely related sequences appropriate for the chosen research question.
Before selecting it:

1. Verify that the source and download page exist.
2. Identify the sequence type and organism or organisms.
3. Record accessions or sequence identifiers.
4. Check whether the sequences require alignment.
5. Check data availability, documentation and applicable terms.
6. Explain why the dataset is appropriate.

Prefer a modest initial dataset that can be processed reliably on an ordinary student laptop. Scale up only after the algorithms and experiments work correctly.
Do not invent accession numbers, dataset statistics or external references.
Include a reproducible data-loading and preprocessing procedure.
7. Software architecture
Propose a modular application in which every mathematical subsystem is independently testable.
Use Python as the preferred implementation language unless a specific reason justifies an alternative.
A reasonable starting stack is:

* Python 3.
* Biopython for sequence parsing and relevant bioinformatics support.
* NumPy for numerical arrays.
* Pandas for tabular data handling.
* NetworkX for graph visualization or independent result verification.
* Matplotlib or Plotly for visualizations.
* Streamlit for an interactive demonstration interface.
* pytest for testing.

Do not allow external libraries to hide all the mathematics that the course is supposed to demonstrate.
Implement key course algorithms explicitly where feasible, and use libraries for data parsing, visualization or independent cross-checking.
Propose a structure such as:

```
genomic-math-project/
│
├── app.py
├── data/
│   ├── synthetic/
│   └── public_genomic/
│
├── modules/
│   ├── sequence_processing.py
│   ├── set_logic.py
│   ├── counting_recursion.py
│   ├── algebraic_structures.py
│   ├── number_theory.py
│   ├── mutation_graph.py
│   ├── tree_analysis.py
│   └── phylogenetic_analysis.py
│
├── tests/
│   ├── test_sequences.py
│   ├── test_counting.py
│   ├── test_algebra.py
│   ├── test_number_theory.py
│   ├── test_graphs.py
│   └── test_trees.py
│
├── requirements.txt
└── README.md
```

Modify the structure as needed, but maintain clear separation between mathematical logic, input data, user interface and tests.
8. User interface requirements
Create a mathematical analysis dashboard with separate sections for each unit.
Section 1 — Dataset Explorer

* Load synthetic and public sequences.
* Show sequence identifiers and lengths.
* Report invalid characters and preprocessing decisions.
* Display the chosen reference sequence.

Section 2 — Set Theory and Logic

* Show mutation-position sets.
* Display intersections, unions and differences.
* Evaluate sequence predicates and relation properties.

Section 3 — Counting and Recursion

* Input sequence length and number of substitutions.
* Display (4^n), (\binom nk3^k) and recurrence results.
* Allow small recursive mutation-generation demonstrations.
* Explain when the search space becomes too large to enumerate.

Section 4 — Algebraic Structures

* Display the selected transformations.
* Show the transformation composition table.
* Show the Cayley table and algebraic-property checks.
* Explain the interpretation of equivalence classes.

Section 5 — Number Theory

* Encode and decode sequences.
* Show modular representations and bucket assignments.
* Display GCD, extended Euclidean algorithm and modular inverse calculations.
* Demonstrate CRT reconstruction.

Section 6 — Mutation Graph

* Visualize sequences as vertices.
* Show the chosen edge rule and distance threshold.
* Display degrees, connected components, matrices and selected paths.
* Allow the threshold to be changed for comparative analysis.

Section 7 — Tree Analysis

* Show a minimum spanning tree.
* Display Prim's and Kruskal's results.
* Display traversals and tree statistics.
* Show a separately constructed distance-based sequence tree if supported.

General UI requirements
Every result should have an explanation of the corresponding mathematical calculation.
Users should be able to distinguish:

* Raw observations.
* Mathematical calculations.
* Algorithmic outputs.
* Biological interpretations.
* Limitations and assumptions.

Do not prioritize decorative animations or excessive UI complexity over correctness.
9. Experiments and evaluation
Design a set of reproducible experiments that demonstrate more than successful program execution.
At minimum, include the following.
Experiment 1 — Mutation-set analysis
Compare selected sequences and verify their mutation sets, intersections, differences and distances manually.
Experiment 2 — Counting versus enumeration
Compare the theoretical number of sequences with a practical recursive enumeration for small values of (n) and (k).
Show that the theoretical count matches the enumerated count.
Experiment 3 — Algebraic verification
Verify the transformation group properties using computed results and provide a mathematical explanation.
Experiment 4 — Modular-indexing analysis
Measure how the choice of modulus changes bucket distribution and collisions.
Explain why numerical collisions do not prove sequence similarity.
Experiment 5 — Mutation-network structure
Study how the graph changes when the distance threshold changes.
Report the number of edges, connected components and other relevant statistics.
Experiment 6 — Shortest paths and spanning trees
Verify Dijkstra, Prim and Kruskal using hand-computable examples or independent checks.
Experiment 7 — Tree analysis
Compare the structures and interpretations of a minimum spanning tree and any optional phylogenetic tree.
Experiment 8 — Real-data evaluation
Apply the validated pipeline to the selected genomic dataset.
Summarize observed patterns without making biological claims unsupported by the data or analysis.
For each experiment, provide:

* Hypothesis or question.
* Input data.
* Mathematical model.
* Procedure.
* Expected result.
* Actual result.
* Interpretation.
* Limitations.

Never fabricate experimental results. Clearly separate expected results from results produced by the running implementation.
10. Correctness, complexity and validation
The project must include systematic software and mathematical validation.
Mathematical validation
For each mathematical subsystem:

* Define all variables and assumptions.
* Derive important formulas.
* Test boundary conditions.
* Provide hand-verified examples.
* Identify counterexamples where appropriate.
* Explain why the algorithm returns the claimed result.

Algorithm validation
Provide time and space complexity for the principal algorithms, including:

* Mutation comparison.
* Recursive mutation generation.
* Sequence encoding.
* Graph construction.
* BFS/DFS.
* Dijkstra.
* Prim/Kruskal.
* Tree traversal.
* Tree construction, if implemented.

For exponential algorithms, discuss their practical limitations.
Software validation
Include:

* Unit tests.
* Invalid-input tests.
* Empty and single-sequence cases.
* Repeated-sequence cases.
* Different sequence lengths.
* Ambiguous or unsupported nucleotide symbols.
* Boundary cases for counting and number theory.
* Disconnected graphs.
* Graphs with multiple equally optimal solutions.

Use independent library implementations as cross-checks when appropriate, while retaining explicit implementations for teaching.
11. Scope and non-goals
This project is an educational computational framework for discrete mathematical analysis of genomic sequences.
It is not:

* A medical diagnostic system.
* A mutation pathogenicity predictor.
* A validated epidemiological forecasting system.
* A replacement for established professional phylogenetic tools.
* A claim of discovering new genomic science.
* Proof that a computed graph or tree represents the true evolutionary history.

Do not add machine learning simply to make the project sound more advanced.
A machine-learning component is unnecessary unless the mathematical and bioinformatics objectives are already satisfied and a clearly justified extension is proposed.
The priority order is:

1. Mathematical correctness.
2. Meaningful integration of all six units.
3. Reproducible experiments.
4. Algorithmic implementation and testing.
5. Clear visual explanation.
6. Additional features only when justified.

12. Expected project deliverables
Develop the project with the following deliverables.
A. Proposal
Include:

* Final title.
* Abstract.
* Problem statement.
* Motivation and real-world relevance.
* Project objectives.
* Research questions.
* Scope.
* Expected contributions.
* Limitations.

B. Mathematical design document
For every unit, include:

* Definitions.
* Mathematical notation.
* Derivations.
* Worked examples.
* Algorithms.
* Complexity analysis.
* How the mathematics contributes to the project.

Include a complete mapping table between syllabus topics and implemented features.
C. System design
Include:

* Architecture diagram.
* Module descriptions.
* Data-flow diagram.
* Dataset-selection rationale.
* Data preprocessing steps.
* Algorithm flowcharts or pseudocode.

D. Working implementation
Provide complete, organized source code.
Include installation instructions, dependencies, commands to run the application and instructions to reproduce the experiments.
Do not provide disconnected code fragments when a working integration is required.
E. Testing and results
Provide:

* Synthetic dataset.
* Dataset-loading instructions for the public data.
* Test cases.
* Actual experiment results from executed code.
* Mathematical verification.
* Complexity analysis.
* Screenshots or plots where available.
* Limitations and interpretation.

F. Academic report
Prepare a report containing:

1. Introduction.
2. Problem definition.
3. Motivation and related work.
4. Mathematical foundations.
5. Proposed methodology.
6. System architecture.
7. Implementation.
8. Experiments and results.
9. Discussion.
10. Limitations.
11. Conclusion and future work.
12. References.

Use genuine, verifiable references for any external biological or computational claims.
G. Presentation and viva
Prepare:

* Presentation outline.
* A concise demonstration script.
* Likely viva questions.
* Detailed mathematical answers.
* Questions that test the distinction between an MST and a phylogenetic tree.
* Questions about group axioms, recurrence relations, modular arithmetic and graph algorithms.
* An explanation of why each of the six units is included.

13. Required working method
Do not attempt to implement everything simultaneously.
Proceed in these phases.
Phase 1 — Mathematical feasibility audit
Evaluate the project against all six units. Identify natural connections, weak connections, assumptions and potential alternatives.
Output a table rating mathematical relevance and identifying evidence needed to justify each unit.
Phase 2 — Final project specification
Define the central research question, exact scope, data strategy, mathematical models and measurable outputs.
Do not begin coding until the model is sufficiently clear.
Phase 3 — Data and mathematical prototype
Construct the synthetic sequences and implement the core formulas and algorithms in small, verifiable examples.
Phase 4 — Integrated backend
Implement sequence processing, the six mathematical modules, graph/tree algorithms and automated tests.
Phase 5 — Real-data evaluation
Select and document a public dataset, process it, run the validated pipeline and record actual results.
Phase 6 — User interface
Integrate the mathematical results into a clear dashboard.
Phase 7 — Final academic deliverables
Prepare the report, diagrams, experiment tables, presentation, viva questions and reproducibility instructions.
At the end of each phase, report what is complete, what remains, the tests passed and any unresolved mathematical problems.
Do not claim that a component works unless its implementation has been tested.
14. Final instruction
Begin with Phase 1: Mathematical Feasibility Audit.
Do not start by generating code or designing the interface.
First determine whether this project is genuinely suitable for a Discrete Mathematics course requiring all six syllabus units.
Pay special attention to Unit III (Algebraic Structures) and Unit IV (Number Theory). If the connections are weak, propose mathematically stronger alternatives within the genomic-analysis domain before committing to an implementation.
The final goal is a project in which a professor can trace every major output back to a mathematical definition, derivation or algorithm and understand why that mathematical operation is being performed.
The project succeeds only if the mathematics is central to the analysis, the algorithms are correct, the biological interpretation is appropriately cautious, and the complete system can be demonstrated and reproduced.
