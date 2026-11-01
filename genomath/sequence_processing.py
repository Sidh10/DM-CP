"""Unit I input layer: parse, validate and filter DNA sequences over {A,C,G,T}; collapse identical sequences into haplotypes."""
from dataclasses import dataclass
from pathlib import Path

from ._core import ALPHABET, SequenceError, mismatch_count, mismatch_positions
from .trace import Trace


@dataclass(frozen=True)
class Record:
    id: str
    description: str
    seq: str


@dataclass(frozen=True)
class Haplotype:
    id: str  # named after its first member in input order
    seq: str
    members: tuple


@dataclass(frozen=True)
class Dataset:
    records: tuple
    dropped: tuple  # (id, reason) pairs
    n: int
    haplotypes: tuple
    reference: Haplotype

    def seqs(self):
        return {h.id: h.seq for h in self.haplotypes}


def parse_fasta(text):
    # hand-written so genomath has no third-party dependency; Biopython is used only to cross-check in tests
    trace = Trace()
    records, header, chunks = [], None, []

    def flush():
        if header is not None:
            rid, _, desc = header.partition(" ")
            records.append(Record(rid, desc.strip(), "".join(chunks)))

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header, chunks = line[1:].strip(), []
            if not header:
                raise SequenceError(f"line {lineno}: empty FASTA header")
        elif header is None:
            raise SequenceError(f"line {lineno}: sequence data before the first '>' header")
        else:
            chunks.append(line)
    flush()

    seen = set()
    for r in records:
        if r.id in seen:
            raise SequenceError(f"duplicate record id {r.id!r}; ids must be unique to name haplotypes")
        seen.add(r.id)
    trace.add("observation", f"parsed {len(records)} FASTA records", ids=[r.id for r in records],
              lengths={r.id: len(r.seq) for r in records})
    return records, trace


def read_fasta(path):
    path = Path(path)
    trace = Trace().add("observation", f"read {path.name}", path=str(path))
    records, t = parse_fasta(path.read_text(encoding="utf-8"))
    return records, trace.merge(t)


def validate(records):
    trace = Trace()
    if not records:
        raise SequenceError("no sequences supplied", trace)
    kept, dropped = [], []
    for r in records:
        seq = r.seq.upper()
        if seq != r.seq:
            trace.add("algorithm", f"{r.id}: lowercase letters converted to uppercase", id=r.id)
        bad = sorted({c for c in seq if c not in ALPHABET})
        if not seq:
            reason = "empty sequence"
        elif bad:
            reason = "non-ACGT characters " + ", ".join(f"{c}x{seq.count(c)}" for c in bad)
        else:
            kept.append(Record(r.id, r.description, seq))
            continue
        dropped.append((r.id, reason))
        trace.add("algorithm", f"dropped {r.id}: {reason}", id=r.id, reason=reason)

    if dropped:
        trace.add("limitation", "dropped records are excluded from every later calculation; "
                  "ambiguity codes are not resolved or imputed", dropped=[d[0] for d in dropped])
    if not kept:
        raise SequenceError("every record was dropped", trace)
    lengths = {}
    for r in kept:
        lengths.setdefault(len(r.seq), []).append(r.id)
    if len(lengths) > 1:
        summary = {n: len(ids) for n, ids in sorted(lengths.items())}
        trace.add("limitation", "unequal lengths: positions cannot be compared without alignment, which is out of scope",
                  lengths=summary)
        raise SequenceError(f"sequences have different lengths {summary}; one common length is required", trace)
    n = len(kept[0].seq)
    trace.add("observation", f"{len(kept)} of {len(records)} records kept; common length n = {n}",
              kept=[r.id for r in kept], n=n)
    return {"kept": kept, "dropped": dropped, "n": n}, trace


def collapse_haplotypes(records):
    trace = Trace()
    groups = {}  # insertion order keeps the first-seen member as the haplotype name
    for r in records:
        groups.setdefault(r.seq, []).append(r.id)
    haplotypes = [Haplotype(ids[0], seq, tuple(ids)) for seq, ids in groups.items()]
    for h in haplotypes:
        if len(h.members) > 1:
            trace.add("algorithm", f"identical sequences {', '.join(h.members)} collapsed into haplotype {h.id}",
                      members=list(h.members))
    trace.add("calculation", f"{len(records)} sequences -> {len(haplotypes)} haplotypes",
              sizes={h.id: len(h.members) for h in haplotypes})
    trace.add("interpretation", "d_H = 0 means equal strings; equality is reflexive, symmetric and transitive, "
              "so haplotypes are its equivalence classes")
    return haplotypes, trace


def hamming(s, t):
    positions = mismatch_positions(s, t)
    trace = Trace()
    trace.add("calculation", f"d_H = |{{i : S_i != T_i}}| = {len(positions)}", positions=positions)
    return len(positions), trace


def position_map(seq):
    # the sequence viewed as a function S: [n] -> Sigma
    mapping = {i: b for i, b in enumerate(seq, start=1)}
    used = sorted(set(seq))
    trace = Trace()
    trace.add("calculation", f"S maps positions 1..{len(seq)} to bases; image = {{{', '.join(used)}}}",
              image=used)
    return mapping, trace


def choose_reference(haplotypes, reference_id=None):
    trace = Trace()
    if reference_id is not None:
        for h in haplotypes:
            if reference_id in h.members:
                trace.add("algorithm", f"reference {h.id} chosen by the user (requested id {reference_id})",
                          reference=h.id)
                return h, trace
        raise SequenceError(f"reference id {reference_id!r} is not among the kept records", trace)
    totals = {h.id: sum(mismatch_count(h.seq, o.seq) for o in haplotypes) for h in haplotypes}
    trace.add("calculation", "for each haplotype, total distance to all others: sum_T d_H(S,T)", totals=totals)
    best = min(haplotypes, key=lambda h: (totals[h.id], h.id))
    tied = sorted(h.id for h in haplotypes if totals[h.id] == totals[best.id])
    trace.add("algorithm", f"medoid = {best.id} (total {totals[best.id]}); ties broken by smallest id",
              reference=best.id, tied=tied)
    trace.add("limitation", "the reference changes every mutation set M(R,S) but no pairwise distance d_H(S,T)")
    return best, trace


def build_dataset(records, reference_id=None):
    trace = Trace()
    checked, t = validate(records)
    trace.merge(t)
    haplotypes, t = collapse_haplotypes(checked["kept"])
    trace.merge(t)
    reference, t = choose_reference(haplotypes, reference_id)
    trace.merge(t)
    data = Dataset(tuple(checked["kept"]), tuple(checked["dropped"]), checked["n"], tuple(haplotypes), reference)
    return data, trace


def load_dataset(path, reference_id=None):
    records, trace = read_fasta(path)
    data, t = build_dataset(records, reference_id)
    return data, trace.merge(t)
