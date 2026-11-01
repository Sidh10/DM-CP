# Stage 8 report: the Streamlit dashboard

Date: 2026-10-31. Status: **complete, waiting for your go/no-go.** Stage 9 (written deliverables) has not started.
No new decisions were needed. The app renders decision D1's five trace kinds as the spec's five result categories.

## What exists
| Path | Contents |
|---|---|
| `app/presentation.py` | **no Streamlit import**: the five category labels, trace flattening (including the D10 "… N more" expansion), matrix / sweep / Cayley table shaping, a hand-written deterministic force-directed graph layout, dendrogram layout, edge styling, and exact-fraction formatting |
| `app/ui.py` | rendering only: category legend, colour-coded tags, trace expanders, tables, and four matplotlib figures (graph, dendrogram, bars, threshold sweep) |
| `app/sections.py` | the seven sections of spec §8 |
| `app/main.py` | entry point, sidebar (dataset, FASTA upload, reference override, section), cached loading |
| `.claude/launch.json` | runs the dashboard on port 8501 |
| `tests/test_app.py` | 13 tests of the presentation layer (175 in total) |
| `scripts/mutation_check.py` | 5 presentation mutants (67 in total) |

**Design decision: a pure presentation layer.** All data shaping lives in `app/presentation.py`, which imports only
`math`. A test asserts this by parsing its import statements, and a second test launches a clean interpreter and
confirms that importing it does not pull in Streamlit. So the layout and table logic is unit-tested like any other
module, and `app/ui.py` is a thin drawing layer.

**NetworkX is still test-only.** The graph layout in the app is a small Fruchterman–Reingold implementation written for
this project: a deterministic circular start, repulsion between all pairs, attraction along edges and a weak pull to
the centre so disconnected components stay in frame.

## How to run it
```
.venv/Scripts/python -m streamlit run app/main.py
```
Then open http://localhost:8501. The dataset, the reference sequence and the section are chosen in the sidebar, and a
FASTA file can be uploaded instead.

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python scripts/mutation_check.py        # ~30 minutes; run it in the background
```
Raw output from my last run:
```
178 passed in 26.11s
67/67 mutants caught
```
All eight experiments were re-run after the app was added and are unchanged (20/20, 16/16, 29/29, 23/23, 14/14,
14/14, 11/11 and 20/20 independent values, plus every consistency check).

## What I checked in the running app
I started the dashboard and drove it in a browser rather than only asserting that it imports. Every section was opened
on **both** datasets; no section produced a traceback and the browser console stayed clean.

Values read off the running page, against the hand derivations:
- **Section 1, Dataset A:** 10 records read, 10 kept, 9 haplotypes, n = 12; syn03 and syn10 collapsed; the d_H column
  reads 0, 1, 2, 3, 1, 2, 2, 2, 6, which is row syn01 of the matrix in `synthetic_expected.md` §3.
- **Section 3, Dataset B:** 4^1701 reported as a 1025-digit number, and the M(n,k) table for n ≤ 8 matches
  `unit2_expected.md` §2 cell for cell, with every row summing to 4^n.
- **Section 6, Dataset A at k = 2:** 14 edges, degrees 6, 5, 4, 2, 4, 4, 1, 2, 0, components 2, independent cycles 7,
  not bipartite, Σ deg = 28 = 2|E| — all matching `unit5_expected.md` §2.
- **Section 7, Dataset A:** MST weight 15 by both algorithms, 8 distinct MSTs, largest edge 6 = k\*, rooted height 3
  with 6 leaves, and the UPGMA merge table reproducing `unit6_expected.md` §6 exactly, including the exact fractions
  (19/6, 23/7, 61/8) and the tie counts 6, 2, 1, 1, 2, 2, 1, 1.
- **Section 7, Dataset B:** MST weight 171, 4 distinct MSTs, **Prim and Kruskal returned different trees** (which the
  ties permit), largest edge 47 = k\*.
- **Guards on Dataset B:** κ and λ are refused above 16 vertices with the reason shown, and the 30×30 adjacency and
  30×31 incidence matrices are omitted with an explanation instead of flooding the page.

## Problems found and fixed during the stage
- **A dendrogram bug, caught by a test I wrote before looking at the page.** I built the bars from the last merge
  downwards, so a parent's x-position was needed before its children had one (`KeyError: 'U7'`). Merge order is the
  correct order, because a node's children are always earlier merges. It is now a mutant.
- **A test that was too crude.** My first "no banned imports" test searched the raw text, and `Bio` matched inside
  "**Bio**logical interpretation". It now parses the import statements with `ast`, which is what I actually meant. I
  also added the reverse check: no `genomath` module imports the app.
- **Figures were theme-dependent.** matplotlib defaults to a white background with black text (I checked the rcParams
  rather than guessing from a screenshot), which would be unreadable under the dark theme. Every figure now has a
  transparent background and mid-grey furniture, so it reads under either theme.
- **Arrow serialisation failures hidden behind a fallback.** The page looked fine, but the server log showed
  `ArrowTypeError` on every table render: my display tables mix numbers with blanks and labels in one column. Cells are
  now rendered as text. I confirmed the fix by re-rendering the most table-heavy section and checking that no new error
  appeared in the log.
- A garbled dead conditional was left in the `table` helper while drafting; removed.
- **A surviving mutant exposed both a test gap and a weak mutant.** Inverting the layout's attraction force went
  unnoticed, because my tests checked only determinism, bounds and distinctness, never the one property a
  force-directed layout exists to provide. I added a test that neighbours end up closer than non-neighbours (on
  Dataset A at k = 1, 2, 3 the mean adjacent distance is 0.09, 0.20 and 0.31 against 0.60, 0.61 and 0.75 apart).
  The mutant *still* survived, because I had mutated only the x component and left attraction intact in y; inverting
  the force itself is the honest mutation, and that one is caught.

## Open issues and doubts
- **The upload path is only lightly exercised.** I verified the two bundled datasets thoroughly; the FASTA upload
  widget shares the same code path (`parse_fasta` → `build_dataset`), and the error path renders `SequenceError` with
  its partial trace, but I did not upload a file through the browser.
- **No automated test renders the Streamlit pages.** The presentation layer is unit-tested and I drove the app by hand,
  but a regression that only affects `ui.py` or `sections.py` drawing code would not fail the suite. A screenshot or
  Playwright test could close that gap if you want it.
- **Layout is not data.** Vertex positions carry no meaning, and the page says so under every graph. Worth repeating in
  the viva, since a reader may assume screen distance means genetic distance.
- Section 6 on Dataset B draws 30 vertices; at high thresholds the graph is dense and the picture is cluttered. The
  sweep table and the statistics remain the readable part, which is why both are shown.
- The dashboard reads the committed `results/*.md` for nothing; it recomputes everything live from `genomath`. That is
  deliberate, so the page can never show a stale number, but it means each section does its own computation on every
  rerun. Loading is cached; the per-section work is not.
