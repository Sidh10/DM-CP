# Stage 0 report (spec Phases 1–2 + data feasibility)

Date: 2026-10-30. Status: **complete, waiting for your go/no-go.** No Stage 1 work has started, and there is no maths code in `genomath/`.

## What exists
| Path | What it is |
|---|---|
| `docs/master_prompt.md` | the spec, saved verbatim from your paste |
| `genomath/*.py` (8 modules + `__init__`) | one-line docstring stubs, no Streamlit |
| `app/`, `experiments/`, `results/` | empty (`.gitkeep`) |
| `tests/test_skeleton.py` | imports every module; asserts there is no `streamlit` in `genomath` (16 tests) |
| `pyproject.toml` | pytest config (`testpaths`, `pythonpath`) |
| `requirements.txt` / `requirements-lock.txt` | direct pins / full `pip freeze` (Python 3.13.11, Windows) |
| `.gitignore` | venv, caches, `data/public_genomic/raw/` |
| `docs/math_audit.md` | Phase 1: per-unit audit, core/pedagogical tags, relevance table, weak links, sources |
| `docs/final_spec.md` | Phase 2: research question, scope, data strategy, models, outputs, experiment → stage map |
| `data/synthetic/dataset_a.fasta`, `dataset_a_invalid.fasta` | Dataset A (10 × 12 nt) and 2 bad records |
| `docs/synthetic_expected.md` | hand derivations for Dataset A (expected values only) |
| `scripts/check_synthetic_expected.py` | stdlib-only brute-force check of that document |
| `scripts/fetch_public_data.py` | Entrez fetch: `select` (search + choose) and `refetch` (exact accessions) |
| `scripts/qc_public_data.py` | QC report generator → `docs/dataset_b_qc.md` |
| `data/public_genomic/accessions.csv`, `h3n2_ha_cds.fasta`, `selection_log.md`, `refetch_log.md` | Dataset B |

Additions beyond your list: `pyproject.toml` (holds the pytest config), `requirements-lock.txt`, the two extra scripts, and
`refetch_log.md`.

## How to verify (commands to re-run)
```
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt      # Linux/macOS: .venv/bin/python
.venv/Scripts/python -m pip check
.venv/Scripts/python -m pytest -q
.venv/Scripts/python -I scripts/check_synthetic_expected.py data/synthetic/dataset_a.fasta
.venv/Scripts/python scripts/qc_public_data.py
.venv/Scripts/python scripts/fetch_public_data.py --email YOU@example.com refetch
```
The last command needs network access. It rewrites `h3n2_ha_cds.fasta` and `accessions.csv` from NCBI; compare them against
your copies (they matched byte-for-byte when I ran it). Running `select` again may pick different records if GenBank has
changed since 2026-10-30, so `refetch` is the reproducible path.

Raw output from my last runs: pytest `16 passed in 0.06s`; `pip check` `No broken requirements found.`;
checker: `MST weight = 15  number of distinct MSTs = 8`, `k=1: |E|=6 components=4`, `n=2: ... Z2 orbits=10 (brute 10), Klein orbits=6 (brute 6)`.
The full Dataset B selection and QC outputs are in `selection_log.md` and `dataset_b_qc.md`.

## Dataset B in one paragraph (raw observations)
Influenza A H3N2 HA CDS from GenBank, human host, collected Oct 2017 – Apr 2018. Search: 11,887 hits; 600 evenly spaced
records examined; 171 passed the filters; 36 selected. All 36 are 1701 nt; all 34 ACGT-only records translate with ATG,
a single terminal stop and no internal stop. 2 records contain ambiguity codes (MH246583.1 R; MK181701.1 R, Y) → 34 clean →
30 haplotypes (two groups of 3 identical). Pairwise d_H: min 1, median 9, max 57 (3.4%). There are **8 d_H = 1 edges**, and
the distance distribution is multimodal. Countries: USA 31, Canada 4, Australia 1. GenBank only, no GISAID. Terms are cited in
`final_spec.md` §4.

## Things I got wrong and fixed during the stage
- My first selection filtered on serotype only. QC showed one haplotype ~225 nt from all others: MK462491.1,
  A/swine/Saskatchewan/SD0258/2017 (host: pig). The candidate pool also held canine and mallard H3N2. I added a human-host
  filter and re-selected. The QC divergence check (max d/n, "a misaligned pair would be near 0.75") remains as a guard.
- esearch returned only 10,000 of 11,887 ids. I added paging so the candidate pool no longer depends on NCBI's result order.
- `refetch` originally overwrote `selection_log.md`. It now writes `refetch_log.md`.

## Open decisions (answer these before or with "go stage 1")
1. Trace shape: `[{"kind", "text", "data"}]` with kinds `observation | calculation | algorithm | interpretation | limitation`.
2. Add allele-level mutation sets M* alongside position sets M (recommended).
3. Add Warshall transitive closure of R_k as the Unit I ↔ V bridge (recommended); the (M*, ⊆) Hasse diagram is optional.
4. Dataset B reference: medoid haplotype (minimum Σ d_H; tie → smallest accession).
5. Canonical k-mer analysis on Dataset B so that Unit III does real work (k ∈ {4, 6, 8}?).
6. UPGMA tie-break: merge the lexicographically smallest label pair. UPGMA expected values for Dataset A are **not** derived until this is fixed.
7. Enumeration cap for recursive generation: 10^5 strings.
8. Modulus set for Exp 4: {16, 64, 7, 13, 17, 31, 61}.

## Flags and doubts
- **Jellyfish source**: located and read the v0.9 manual (2010). The canonical definition is under `--both-strands`.
  I did not check current-release wording or locate the journal paper, so neither is cited.
- **E-utilities rate with an API key**: the NCBI policies page gives 3 requests/second (it does not mention keys). The
  E-utilities book page (NBK25497) returned a CAPTCHA, so the higher with-key rate is unverified. The script has no
  API-key support at all and always waits 0.4 s.
- **Positional homology on Dataset B** is inferred from equal CDS length plus low divergence (≤ 3.4%), not from an alignment.
  A compensating insertion/deletion pair could in principle pass both checks. I consider this unlikely here, but it is unproven.
- **No clade names** are given for the multimodal distances. Naming them would need a source I have not located.
- `fetch_public_data.py` has no unit tests. Its behaviour is evidenced only by the logs and the byte-identical `refetch`.
- `accessions.csv` adds a `collection_date` column beyond the five you listed; it is needed to show the season window.
- In `qc_public_data.py`, the collection-month line is sorted alphabetically, not chronologically (cosmetic).
