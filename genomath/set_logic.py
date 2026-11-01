"""Unit I: mutation-position sets, set operations, predicates and relation-property tests."""
from ._core import ALPHABET, mismatch_count, mismatch_positions, reverse_complement
from .trace import Trace



def mutation_set(ref, seq):
    positions = mismatch_positions(ref, seq)
    trace = Trace().add_steps("calculation", (f"position {i}: R = {ref[i - 1]}, S = {seq[i - 1]} -> {i} in M(R,S)"
                                              for i in positions))
    trace.add("calculation", f"M(R,S) = {_fmt(positions)}, |M(R,S)| = {len(positions)}", M=positions)
    return frozenset(positions), trace


def allele_mutation_set(ref, seq):
    pairs = [(i, seq[i - 1]) for i in mismatch_positions(ref, seq)]
    trace = Trace()
    trace.add("calculation", f"M*(R,S) = {_fmt(pairs)}: each changed position with the base S carries there",
              M_star=pairs)
    trace.add("interpretation", "M* separates different substitutions at the same position; M does not")
    return frozenset(pairs), trace


def reference_table(ref_id, ref, seqs):
    trace = Trace()
    rows = []
    for sid, s in seqs.items():
        pos = mismatch_positions(ref, s)
        rows.append({"id": sid, "M": pos, "size": len(pos), "M_star": [(i, s[i - 1]) for i in pos]})
    trace.add("calculation", f"mutation sets of {len(rows)} sequences relative to reference {ref_id}",
              sizes={r["id"]: r["size"] for r in rows})
    trace.add("limitation", "mutation sets compare equal-length sequences position by position; "
              "insertions and deletions are invisible")
    return rows, trace



def set_operations(a, b, names=("A", "B")):
    # explicit membership loop so the trace can show where every element lands in the Venn diagram
    na, nb = names
    trace = Trace()
    inter, only_a, only_b, steps = [], [], [], []
    for x in sorted(list(a) + [y for y in b if y not in a]):
        in_a, in_b = x in a, x in b
        if in_a and in_b:
            inter.append(x)
            where = f"{na} and {nb}"
        elif in_a:
            only_a.append(x)
            where = f"{na} only"
        else:
            only_b.append(x)
            where = f"{nb} only"
        steps.append(f"{x}: in {where}")
    trace.add_steps("algorithm", steps)
    union = sorted(inter + only_a + only_b)
    result = {
        "intersection": frozenset(inter),
        "union": frozenset(union),
        "a_minus_b": frozenset(only_a),
        "b_minus_a": frozenset(only_b),
        "symmetric_difference": frozenset(only_a + only_b),
        "size_a": len(a),
        "size_b": len(b),
    }
    trace.add("calculation", f"{na} ∩ {nb} = {_fmt(inter)}; {na} ∪ {nb} = {_fmt(union)}; "
              f"{na} \\ {nb} = {_fmt(only_a)}; {nb} \\ {na} = {_fmt(only_b)}; "
              f"{na} Δ {nb} = {_fmt(sorted(only_a + only_b))}")
    trace.add("calculation", f"check |{na} ∪ {nb}| = |{na}| + |{nb}| - |{na} ∩ {nb}|: "
              f"{len(union)} = {len(a)} + {len(b)} - {len(inter)}",
              holds=len(union) == len(a) + len(b) - len(inter))
    return result, trace


def distance_decomposition(ref, s, t):
    # d_H(S,T) = |M(R,S) Δ M(R,T)| + |{i in M(R,S) ∩ M(R,T) : S_i != T_i}|
    ms, mt = set(mismatch_positions(ref, s)), set(mismatch_positions(ref, t))
    sym = sorted(ms ^ mt)
    shared_diff = sorted(i for i in ms & mt if s[i - 1] != t[i - 1])
    direct = mismatch_count(s, t)
    rhs = len(sym) + len(shared_diff)
    trace = Trace()
    trace.add("calculation", f"M(R,S) Δ M(R,T) = {_fmt(sym)} (size {len(sym)})", sym=sym)
    trace.add("calculation", f"shared positions carrying different bases = {_fmt(shared_diff)} (size {len(shared_diff)})",
              shared_diff=shared_diff)
    trace.add("calculation", f"{len(sym)} + {len(shared_diff)} = {rhs}; direct d_H(S,T) = {direct}",
              holds=rhs == direct)
    return {"sym_diff": sym, "shared_differing": shared_diff, "rhs": rhs, "d_H": direct, "holds": rhs == direct}, trace



class Predicate:
    def __init__(self, name, fn):
        self.name, self._fn = name, fn

    def __call__(self, seq):
        return bool(self._fn(seq))

    def __and__(self, other):
        return Predicate(f"({self.name} ∧ {other.name})", lambda s: self(s) and other(s))

    def __or__(self, other):
        return Predicate(f"({self.name} ∨ {other.name})", lambda s: self(s) or other(s))

    def __invert__(self):
        return Predicate(f"¬{self.name}", lambda s: not self(s))

    def __repr__(self):
        return f"Predicate({self.name})"


def within_distance(ref, k):
    return Predicate(f"d_H(R,S) ≤ {k}", lambda s: mismatch_count(ref, s) <= k)


def base_at(i, base):
    if base not in ALPHABET or i < 1:
        raise ValueError(f"need a position >= 1 and a base in {ALPHABET}")
    return Predicate(f"S_{i} = {base}", lambda s: s[i - 1] == base)


def carries_mutation(ref, i, base):
    if ref[i - 1] == base:
        raise ValueError(f"({i},{base}) is not a mutation: the reference already has {base} at {i}")
    return Predicate(f"S carries ({i},{base})", lambda s: s[i - 1] == base)


def shares_mutation_with(ref, other, i):
    # same position AND same non-reference base
    return Predicate(f"S shares the mutation at {i} with T",
                     lambda s: s[i - 1] == other[i - 1] != ref[i - 1])


def filter_sequences(seqs, predicate, parts=()):
    trace = Trace()
    table, kept = [], []
    for sid, s in seqs.items():
        row = {p.name: p(s) for p in parts}
        row[predicate.name] = predicate(s)
        table.append({"id": sid, **row})
        if row[predicate.name]:
            kept.append(sid)
    trace.add("calculation", f"truth table of {predicate.name} over {len(seqs)} sequences", table=table)
    trace.add("algorithm", f"kept {len(kept)}: {', '.join(kept) or 'none'}", kept=kept)
    return kept, trace



def threshold_relation(seqs, k):
    ids = list(seqs)
    pairs = {(a, b) for a in ids for b in ids if mismatch_count(seqs[a], seqs[b]) <= k}
    trace = Trace()
    trace.add("calculation", f"R_{k} = {{(S,T) : d_H(S,T) ≤ {k}}} has {len(pairs)} ordered pairs on {len(ids)} elements",
              size=len(pairs))
    return pairs, trace


def identity_relation(seqs):
    ids = list(seqs)
    pairs = {(a, b) for a in ids for b in ids if seqs[a] == seqs[b]}
    trace = Trace().add("calculation", f"S ≡ T iff d_H(S,T) = 0: {len(pairs)} ordered pairs", size=len(pairs))
    return pairs, trace


def orientation_relation(seqs):
    ids = list(seqs)
    rc = {a: reverse_complement(seqs[a]) for a in ids}
    pairs = {(a, b) for a in ids for b in ids if seqs[b] in (seqs[a], rc[a])}
    trace = Trace()
    trace.add("calculation", f"S ~ T iff T ∈ {{S, RC(S)}}: {len(pairs)} ordered pairs", size=len(pairs))
    trace.add("interpretation", "S and RC(S) describe the same double-stranded molecule read from opposite strands")
    return pairs, trace


def relation_properties(elements, pairs):
    elements = list(elements)
    succ = {a: [b for b in elements if (a, b) in pairs] for a in elements}
    trace = Trace()

    reflexive = next(((a,) for a in elements if (a, a) not in pairs), None)
    symmetric = next(((a, b) for a in elements for b in succ[a] if (b, a) not in pairs), None)
    antisymmetric = next(((a, b) for a in elements for b in succ[a] if a != b and (b, a) in pairs), None)
    transitive = next(((a, b, c) for a in elements for b in succ[a] for c in succ[b] if (a, c) not in pairs), None)

    report = {}
    for name, witness, failure in [
        ("reflexive", reflexive, "({0},{0}) ∉ R"),
        ("symmetric", symmetric, "({0},{1}) ∈ R but ({1},{0}) ∉ R"),
        ("antisymmetric", antisymmetric, "({0},{1}) and ({1},{0}) ∈ R with {0} ≠ {1}"),
        ("transitive", transitive, "({0},{1}) ∈ R and ({1},{2}) ∈ R but ({0},{2}) ∉ R"),
    ]:
        holds = witness is None
        report[name] = {"holds": holds, "witness": witness}
        text = f"{name}: holds" if holds else f"{name}: fails, counterexample " + failure.format(*witness)
        trace.add("algorithm", text, holds=holds, witness=witness)
    report["equivalence"] = all(report[p]["holds"] for p in ("reflexive", "symmetric", "transitive"))
    trace.add("calculation", f"equivalence relation: {report['equivalence']}")
    return report, trace


def transitive_closure(elements, pairs):
    # Warshall: after pivot k, R contains (i,j) whenever a path i -> j uses intermediate vertices among the first k
    elements = list(elements)
    idx = {e: i for i, e in enumerate(elements)}
    n = len(elements)
    m = [[False] * n for _ in range(n)]
    for a, b in pairs:
        m[idx[a]][idx[b]] = True
    trace = Trace()
    for k in range(n):
        added = 0
        for i in range(n):
            if m[i][k]:
                for j in range(n):
                    if m[k][j] and not m[i][j]:
                        m[i][j] = True
                        added += 1
        if added:
            trace.add("algorithm", f"pivot {elements[k]}: added {added} pairs", pivot=elements[k], added=added)
    closure = {(elements[i], elements[j]) for i in range(n) for j in range(n) if m[i][j]}
    trace.add("calculation", f"closure has {len(closure)} pairs (input had {len(set(pairs))}); O(n^3) with n = {n}")
    return closure, trace


def equivalence_classes(elements, pairs):
    report, trace = relation_properties(elements, pairs)
    if not report["equivalence"]:
        failed = [p for p in ("reflexive", "symmetric", "transitive") if not report[p]["holds"]]
        raise ValueError(f"not an equivalence relation: fails {', '.join(failed)}")
    classes, seen = [], set()
    for a in elements:
        if a not in seen:
            cls = [b for b in elements if (a, b) in pairs]
            seen.update(cls)
            classes.append(cls)
    trace.add("calculation", f"{len(classes)} equivalence classes", classes=classes)
    return classes, trace



def function_properties(domain, f, codomain, name="f"):
    preimages = {}
    for x in domain:
        y = f(x)
        if y not in codomain:
            raise ValueError(f"{name}({x!r}) = {y!r} is outside the codomain")
        preimages.setdefault(y, []).append(x)
    collision = next(((v[0], v[1], y) for y, v in preimages.items() if len(v) > 1), None)
    missing = next((y for y in codomain if y not in preimages), None)
    trace = Trace()
    trace.add("algorithm", f"{name} injective: " + ("yes" if collision is None else
              f"no, {name}({collision[0]!r}) = {name}({collision[1]!r}) = {collision[2]!r}"))
    trace.add("algorithm", f"{name} surjective: " + ("yes" if missing is None else f"no, nothing maps to {missing!r}"))
    result = {"injective": collision is None, "surjective": missing is None,
              "bijective": collision is None and missing is None,
              "collision": collision, "missing": missing}
    return result, trace


def _fmt(items):
    return "{" + ", ".join(str(x) for x in items) + "}"
