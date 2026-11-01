"""Trace records: the step-by-step explanation returned with every algorithm result."""

# the UI tags each block with one of these, so the vocabulary is closed
KINDS = ("observation", "calculation", "algorithm", "interpretation", "limitation")
STEP_CAP = 10  # decision D10


class Trace(list):
    def add(self, kind, text, /, **data):  # positional-only, so data may use any field name, including "kind"
        if kind not in KINDS:
            raise ValueError(f"unknown trace kind {kind!r}; expected one of {KINDS}")
        self.append({"kind": kind, "text": text, "data": data})
        return self

    def add_steps(self, kind, texts, cap=None):
        # per-element steps (one per position, digit, element): show the first few, keep the rest in data
        cap = STEP_CAP if cap is None else cap
        texts = list(texts)
        for t in texts[:cap]:
            self.add(kind, t)
        if len(texts) > cap:
            self.add(kind, f"… {len(texts) - cap} more steps of the same kind", omitted=texts[cap:])
        return self

    def merge(self, other):
        self.extend(other)
        return self

    def kinds(self):
        return [step["kind"] for step in self]
