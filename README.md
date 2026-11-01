# Genomic Mutation & Evolution Analyzer

A discrete-mathematics framework for analysing equal-length DNA sequences, built for **MM0903 Discrete Mathematics**
(B.Tech CSE-AI, second year, VIT Pune, 2026–27).

All six syllabus units are implemented as hand-written pure-Python algorithms. Every algorithm returns its result
together with a trace, so each output can be followed back to a definition, a derivation or a named algorithm.

## Status

| | |
|---|---|
| Tests | 179 passing (`pytest -W error`) |
| Mutation check | 67 of 67 deliberate bugs caught |
| Experiments | 8 of 8, with **147 independent expected values matched**, 28 consistency checks, 35 observations |
| Datasets | A: 10 synthetic sequences, every answer hand-derived · B: 36 GenBank influenza A(H3N2) HA coding sequences |

## Install

Requires Python 3.11 or newer (developed and tested on 3.13.11, Windows).

```
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
```
On Linux or macOS use `.venv/bin/python` throughout. `requirements-lock.txt` holds the full transitive pin set.

## Run the dashboard

```
.venv/Scripts/python -m streamlit run app/main.py
```
Then open <http://localhost:8501>. Choose the dataset, the reference sequence and the section in the sidebar, or upload
your own FASTA file. The seven sections follow the six units; each result is tagged as a raw observation, a
mathematical calculation, an algorithmic output, a biological interpretation, or a limitation.

## Reproduce the experiments

Each experiment writes `results/expN.md` with **expected** and **actual** values in separate tables:

```
.venv/Scripts/python -m experiments.exp1_mutation_sets
.venv/Scripts/python -m experiments.exp2_counting_enumeration
.venv/Scripts/python -m experiments.exp3_algebraic_verification
.venv/Scripts/python -m experiments.exp4_modular_indexing
.venv/Scripts/python -m experiments.exp5_network_structure
.venv/Scripts/python -m experiments.exp6_paths_and_trees
.venv/Scripts/python -m experiments.exp7_tree_analysis
.venv/Scripts/python -m experiments.exp8_real_data
```

## Verify

```
.venv/Scripts/python -m pytest -q -W error          # 179 tests, warnings treated as errors
.venv/Scripts/python scripts/mutation_check.py      # ~30 min; inject 67 bugs, confirm each is caught
.venv/Scripts/python -I scripts/check_synthetic_expected.py data/synthetic/dataset_a.fasta
.venv/Scripts/python scripts/qc_public_data.py      # regenerate the quality-control report for Dataset B
```

`check_synthetic_expected.py` and `qc_public_data.py` share no code with `genomath`, so the values they produce are
independent of the implementation they check.

## Re-download Dataset B

Needs network access and an email address, which NCBI requires with every Entrez request.

```
.venv/Scripts/python scripts/fetch_public_data.py --email YOU@example.com refetch
```
`refetch` re-downloads exactly the accessions already listed in `data/public_genomic/accessions.csv`, which is the
reproducible path. `select` runs the original GenBank search and may pick different records if the database has changed.

## Layout

```
genomath/     the six units, hand-written, no UI or oracle imports
app/          Streamlit dashboard (presentation.py is pure and unit-tested)
experiments/  expN_*.py -> results/expN.md, with expected values kept in expected_*.py
tests/        pytest suite; NetworkX, SymPy, Biopython and math used only here, as oracles
scripts/      data fetch, quality control, independent checkers, mutation check
docs/         derivations, design documents, the report, stage reports
data/         synthetic/ and public_genomic/
results/      experiment output
```

## Documentation

| Document | Contents |
|---|---|
| [docs/report.md](docs/report.md) | the academic report: the whole project in one narrative |
| [docs/proposal.md](docs/proposal.md) | title, abstract, objectives, research questions, scope |
| [docs/mathematical_design.md](docs/mathematical_design.md) | per unit: definitions, derivations, worked examples, algorithms, complexity, and the syllabus mapping table |
| [docs/system_design.md](docs/system_design.md) | architecture, data flow, module descriptions, pseudocode |
| [docs/testing_and_results.md](docs/testing_and_results.md) | the validation method and every result |
| [docs/presentation_and_viva.md](docs/presentation_and_viva.md) | presentation outline, demonstration script, viva questions with answers |
| [docs/decisions.md](docs/decisions.md) | every design decision, with when it was made |
| [docs/math_audit.md](docs/math_audit.md) | the feasibility audit: which connections are strong, which are weak |
| `docs/unit*_expected.md`, `docs/synthetic_expected.md` | hand derivations, written before the code |
| `docs/stage*_report.md` | what was built in each stage, including the mistakes found |

## Scope

An educational framework for the discrete-mathematical analysis of genomic sequences. It is **not** a diagnostic
system, a pathogenicity predictor, an epidemiological tool, a replacement for established phylogenetic software, or
evidence of evolutionary history.

## Data

Dataset B comes from GenBank (NCBI), which states that it places no restrictions on the use or distribution of GenBank
data, while noting that submitters may claim rights NCBI cannot assess
(<https://www.ncbi.nlm.nih.gov/genbank/about/>, <https://www.ncbi.nlm.nih.gov/home/about/policies/>). No GISAID data is
used. The download script obeys NCBI's scripted-access limit of three requests per second and sends the required
contact email. Accessions are listed in `data/public_genomic/accessions.csv`.
