"""Untraced primitives shared by the unit modules; public algorithms wrap them and return (result, trace)."""

ALPHABET = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


class SequenceError(ValueError):
    def __init__(self, message, trace=None):
        super().__init__(message)
        # carried so a caller (e.g. the UI) can still show the drop log that led to the failure
        self.trace = trace


def require_same_length(s, t):
    if len(s) != len(t):
        raise SequenceError(f"Hamming comparison needs equal lengths, got {len(s)} and {len(t)}")


def mismatch_positions(s, t):
    require_same_length(s, t)
    return [i for i, (a, b) in enumerate(zip(s, t), start=1) if a != b]


def mismatch_count(s, t):
    require_same_length(s, t)
    count = 0
    for a, b in zip(s, t):
        if a != b:
            count += 1
    return count


def complement(s):
    try:
        return "".join(COMPLEMENT[c] for c in s)
    except KeyError as e:
        raise SequenceError(f"complement is defined only on {ALPHABET}, got {e.args[0]!r}") from None


def reverse_complement(s):
    return complement(s)[::-1]
