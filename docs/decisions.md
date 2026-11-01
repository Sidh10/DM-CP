# Design decisions

Single record of project decisions. Later documents defer to this file where they differ.
"Yours" = given in your Stage 0 brief. "Confirmed" = proposed by me and confirmed by you on 2026-10-30, after the Stage 4 review.

## Given in the Stage 0 brief (yours)
- Course algorithms are hand-written in pure Python. NetworkX, SymPy, `math` and Biopython appear only in tests as
  oracles (and Biopython in `scripts/` for download and parsing).
- Every algorithm returns `(result, trace)`. UI blocks are tagged raw observation / calculation / algorithmic output /
  biological interpretation / limitation.
- Package layout `genomath/` with the §7 module names; no Streamlit imports. Also `app/`, `tests/`, `experiments/`,
  `results/`, `docs/`, `scripts/`, `data/…`. pathlib only.
- `experiments/expN_*.py` write `results/expN.md` with expected and actual in separate fields.
- Dijkstra runs on the sparse threshold graph (edge iff d_H ≤ k, weight d_H); report path cost vs direct d_H.
- Identical sequences collapse into haplotypes; d_H ≤ k is tested as non-transitive; non-ACGT records are dropped and logged.
- Unit III: {I,RC} and the Klein four-group; RC-fixed count; Burnside; no arbitrary position permutations.
- Unit IV: complement identity; bad hash x mod 4^j vs prime moduli; pigeonhole; CRT on small k-mers; canonical k-mers;
  Euler/Fermat labelled pedagogical.
- Unit VI: MST and UPGMA kept separate; state UPGMA's clock assumption and that rooting the MST is arbitrary.
- No git commits or pushes by me.

## Confirmed 2026-10-30
| # | Decision | Introduced |
|---|---|---|
| D1 | Trace = list of `{"kind", "text", "data"}`; kind ∈ observation, calculation, algorithm, interpretation, limitation | Stage 1 |
| D2 | Allele-level mutation sets M* alongside position sets M | Stage 1 |
| D3 | Warshall transitive closure of R_k as the Unit I ↔ V bridge | Stage 1 |
| D4 | Dataset B reference = medoid haplotype (min Σ d_H; ties → smallest id); result CY249987.1 | Stage 1 |
| D5 | Recursive enumeration cap 10^5 strings | Stage 2 |
| D6 | Canonical k-mer census of Dataset B for k = 4, 6, 8 | Stage 4 |
| D7 | Modulus set {16, 64, 7, 13, 17, 31, 61} | Stage 4 |
| D8 | UPGMA tie-break: merge the tied pair whose sorted cluster labels come first lexicographically | implemented in Stage 6 |
| D9 | Base-level Klein group {id, (A T)(C G), (A G)(C T), (A C)(G T)} classifying substitutions as transitions/transversions, counts only, no claims about typical ratios | implemented in Stage 7 (Exp 8) |
| D10 | Long traces are summarised inside `genomath`: at most 10 per-element steps, then one "… N more" step whose data holds the omitted step texts | applied after the Stage 4 review |

## Conventions adopted along the way (no decision needed, recorded for transparency)
- Primitives in `genomath/_core.py` return plain values (they run in inner loops); every public function returns `(result, trace)`.
- Extra modules beyond §7: `genomath/trace.py`, `genomath/_core.py`.
- Experiment checks are either **independent** (actual compared with a value fixed beforehand by hand or by the
  genomath-free QC script) or **consistency** (an invariant or two genomath computations agreeing). Results report the
  two counts separately.
- Optional, not done: the Hasse diagram of (M*, ⊆).
