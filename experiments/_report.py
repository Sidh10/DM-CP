"""Shared writer for results/expN.md: expected and actual are rendered as separate tables, then compared."""
import platform
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def fmt(v):
    if isinstance(v, (set, frozenset)):
        v = sorted(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(fmt(x) for x in v) + "]" if isinstance(v, list) else "(" + ", ".join(fmt(x) for x in v) + ")"
    if isinstance(v, dict):
        return "{" + ", ".join(f"{k}: {fmt(x)}" for k, x in v.items()) + "}"
    return str(v)


def normalise(v):
    if isinstance(v, (set, frozenset)):
        return sorted(normalise(x) for x in v)
    if isinstance(v, (list, tuple)):
        return [normalise(x) for x in v]
    if isinstance(v, dict):
        return {k: normalise(x) for k, x in v.items()}
    return v


def cell(v):
    # a literal | inside a cell (e.g. |R_1|) would split the markdown table column
    return (v if isinstance(v, str) else fmt(v)).replace("|", "\\|")


class Check:
    def __init__(self, item, expected, source, actual=None, actual_only=False, consistency=False):
        self.item, self.expected, self.source, self.actual = item, expected, source, actual
        self.actual_only = actual_only  # no independent expected value exists
        # consistency: can only fail if genomath disagrees with itself (round trip, two genomath paths, a data-free bound)
        self.consistency = consistency

    @property
    def match(self):
        return None if self.actual_only else normalise(self.expected) == normalise(self.actual)


def write(n, title, command, text_sections, checks, notes_after=()):
    lines = [f"# Experiment {n}: {title}", "",
             f"- run date: {date.today().isoformat()}",
             f"- command: `{command}`",
             f"- Python {platform.python_version()} on {platform.system()}", ""]
    for heading, body in text_sections[:4]:  # question, input, model, procedure
        lines += [f"## {heading}", body, ""]

    lines += ["## Expected", "Values fixed before the run; never produced by `genomath`.", "",
              "| # | item | expected | source |", "|---|---|---|---|"]
    for i, c in enumerate(checks, 1):
        exp = "*(none: reported as observation)*" if c.actual_only else fmt(c.expected)
        source = c.source + (" — *consistency check*" if c.consistency else "")
        lines.append(f"| {i} | {cell(c.item)} | {cell(exp)} | {source} |")

    lines += ["", "## Actual", "Values produced by this run of `genomath`.", "",
              "| # | item | actual |", "|---|---|---|"]
    for i, c in enumerate(checks, 1):
        lines.append(f"| {i} | {cell(c.item)} | {cell(c.actual)} |")

    independent = [c for c in checks if not c.actual_only and not c.consistency]
    consistency = [c for c in checks if c.consistency]
    ok_i = sum(bool(c.match) for c in independent)
    ok_c = sum(bool(c.match) for c in consistency)
    observed = len(checks) - len(independent) - len(consistency)
    summary = (f"{ok_i}/{len(independent)} independent expected values matched; "
               f"{ok_c}/{len(consistency)} consistency checks passed; {observed} observation-only items")
    lines += ["", "## Comparison", "", "| # | item | match |", "|---|---|---|"]
    for i, c in enumerate(checks, 1):
        verdict = "n/a" if c.actual_only else ("yes" if c.match else "**NO**")
        lines.append(f"| {i} | {cell(c.item)} | {verdict}{' (consistency)' if c.consistency else ''} |")
    lines += ["", f"**{summary}.**", "",
              "Independent: compared with a value fixed outside `genomath` (hand derivation, a hand-proved theorem's "
              "numeric prediction, or the genomath-free QC script). Consistency: can only fail if `genomath` disagrees "
              "with itself.", ""]
    for note in notes_after:
        lines += [note, ""]
    for heading, body in text_sections[4:]:  # interpretation, limitations
        lines += [f"## {heading}", body, ""]
    out = ROOT / "results" / f"exp{n}.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    return out, summary
