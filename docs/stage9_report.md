# Stage 9 report: written deliverables

Date: 2026-10-31. Status: **complete. This was the final stage of the plan.**
No new decisions were needed; no code in `genomath` was changed.

## What exists
All seven deliverables of specification §12:

| § | Deliverable | File | Lines |
|---|---|---|---|
| A | Proposal | `docs/proposal.md` | 90 |
| B | Mathematical design document | `docs/mathematical_design.md` | 364 |
| C | System design | `docs/system_design.md` | 203 |
| D | Working implementation — install, run, reproduce | `README.md` | 114 |
| E | Testing and results | `docs/testing_and_results.md` | 171 |
| F | Academic report (12 sections) | `docs/report.md` | 297 |
| G | Presentation and viva | `docs/presentation_and_viva.md` | 264 |

Also added: `scripts/check_report_figures.py`, which recomputes every figure the documents quote, and a test that runs
it (179 tests now).

**Deliverable B** states each unit's objects, definitions, theorems with proofs, a worked example with real numbers,
algorithms with complexity, and that unit's contribution — followed by the **complete syllabus-to-implementation
mapping table**, which names for every topic where it is implemented and where it is verified. Topic names for
Units I–III come from the actual course slide titles; the slide folders for Units IV–VI are not present, so those names
follow the specification, and the document says so.

**Deliverable G** contains a 15-slide outline with timings, a four-minute demonstration script keyed to the dashboard,
and **33 viva questions with full answers**, grouped as the specification asks: the MST-versus-phylogeny distinction
(5 questions), group axioms (6), recurrence relations and counting (4), modular arithmetic (6), graph algorithms (5),
why each unit is included (2), and method (5).

## Verify
```
.venv/Scripts/python -m pytest -q -W error
.venv/Scripts/python scripts/check_report_figures.py
.venv/Scripts/python scripts/mutation_check.py        # ~30 minutes; run it in the background
```
Raw output from my last run:
```
179 passed in 22.97s
78 figures verified against the live code
67/67 mutants caught
```

## Keeping the documents honest
A report full of numbers is a report full of opportunities to be wrong, and these documents quote well over a hundred.
Rather than trust my transcription, `scripts/check_report_figures.py` **recomputes each one from the live code and the
results files** and fails if any disagrees — the experiment totals, both datasets' structural results, and the specific
values quoted in the viva answers (CRT recovering 15 of 64 codes under moduli (3,5), the orders of 4 modulo 7 and 13,
the two-substitution cancellation, the monoid's size and units, and so on). A test runs it, so documentation drift now
fails the suite.

**It immediately caught an error.** I had written "33 observations" across four documents; the correct total is **35**.
I had mis-added the per-experiment counts by transcribing Experiment 6's 4 as 2. The independent total (147) and the
consistency total (28) were right. Corrected in `README.md`, `docs/proposal.md`, `docs/report.md` and
`docs/testing_and_results.md`.

This is the fourth time in the project that a check written to be independent of the thing it checks found a real
mistake, and the second time the mistake was mine rather than the code's.

## References: what is cited and what is not
The reference list in `docs/report.md` §12 has five entries, and that is deliberate. It contains only sources I located
and read:
- the Jellyfish manual for the canonical k-mer convention, cited with the caveat that it is the 2010 version 0.9
  manual and that I did not locate the journal paper;
- two NCBI pages for the GenBank terms and the scripted-access limits;
- the course slide titles, used only for syllabus topic names;
- the project specification.

Claims with no located source are reported as unsupported in the text rather than attributed to a plausible-sounding
reference. The clearest case is the transition/transversion ratio: the counts are given (250 to 91 relative to the
reference, 131 to 40 along MST edges) against a combinatorial baseline of 0.5, and the report states plainly that
explaining the gap would need a substitution model and a source, and so draws no conclusion. "Robinson–Foulds" is
likewise not used as a name for the bipartition symmetric difference, because I did not locate the original reference.

## Open issues and doubts
- **The deliverables are long** (about 1,500 lines in total, plus the existing derivation and stage documents). They
  are structured for reference rather than for reading end to end; `docs/report.md` is the one document meant to be
  read start to finish.
- **Two numbers in the documents are not machine-checked**: the line counts in this table, and "about 2,400 lines in
  `genomath`", which the figure check verifies only as a range (2,300–2,500) so it does not fail on a one-line edit.
- **The viva answers are mine, not a transcript.** They are written to be defensible from the project's own evidence,
  and every number in them is verified by the figure check — but an examiner may ask something not anticipated here.
- The presentation outline assumes 12–15 minutes and says which slides to cut if time is short. If your slot differs,
  the timings need adjusting.
- `docs/report.md` §3 discusses related work using only the two external sources located. If you can supply course
  references or a textbook for the discrete-mathematics results, that section would be the place to add them — I have
  not added any I could not verify.

## The project as a whole
All nine stages are complete: the audit and specification (0), the six units (1–6), the real-data evaluation (7), the
dashboard (8), and these documents (9). Across the project: **179 tests**, **67 of 67 injected bugs caught**, **eight
experiments** with 147 independent expected values matched, 28 consistency checks and 35 observations, and **two
theorems proved** that were not in the original plan — the exact single-substitution bucket criterion and the monoid
proof that bounds what "mutation algebra" may claim.

Ten errors were found during development and are listed in `docs/testing_and_results.md`. Three of them were in my own
hand derivations rather than in the code, which is the strongest evidence I can offer that the expected-values-first
method was doing real work rather than confirming what the code already said.
