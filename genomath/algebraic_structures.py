"""Unit III: complement/reverse/RC transformations, Z2 and Klein four-group tables, axioms, orbits, Burnside counting, and the base-level Klein group behind transitions/transversions."""
from itertools import product

from ._core import ALPHABET, SequenceError, complement, require_same_length, reverse_complement
from .trace import Trace


def _identity(s):
    return s


def _reverse(s):
    return s[::-1]


TRANSFORMS = {"I": _identity, "R": _reverse, "C": complement, "RC": reverse_complement}
Z2 = ("I", "RC")
V4 = ("I", "R", "C", "RC")
ENUMERATION_LIMIT = 10  # 4^10 ≈ 10^6 strings; exhaustive checks above this are refused


def words(n_max, n_min=0):
    for n in range(n_min, n_max + 1):
        for w in product(ALPHABET, repeat=n):
            yield "".join(w)


def _compose(a, b):
    fa, fb = TRANSFORMS[a], TRANSFORMS[b]
    return lambda x: fa(fb(x))


def transform(name, s):
    if name not in TRANSFORMS:
        raise ValueError(f"unknown transformation {name!r}; expected one of {list(TRANSFORMS)}")
    bad = sorted(set(s) - set(ALPHABET))
    if bad:
        raise SequenceError(f"transformations are defined on {ALPHABET} only, got {bad}")
    out = TRANSFORMS[name](s)
    trace = Trace().add("calculation", f"{name}({s}) = {out}", input=s, output=out)
    return out, trace


def identify(fn, candidates, n_max=6):
    # equal as functions means equal on every word of length 0..n_max, not on a sample
    domain = list(words(n_max))
    for c in candidates:
        if all(fn(x) == TRANSFORMS[c](x) for x in domain):
            return c
    return None


def distinct_on(elements, n):
    # which elements coincide as maps on Σ^n for this single length (faithfulness of the action)
    domain = list(words(n, n))
    classes = []
    for e in elements:
        for cls in classes:
            if all(TRANSFORMS[e](x) == TRANSFORMS[cls[0]](x) for x in domain):
                cls.append(e)
                break
        else:
            classes.append([e])
    trace = Trace().add("calculation", f"on Σ^{n}: " + "; ".join(" = ".join(c) for c in classes),
                        classes=classes, faithful=len(classes) == len(elements))
    return classes, trace


def cayley_table(elements, n_max=6):
    trace = Trace()
    trace.add("algorithm", f"entry (a,b) = a∘b (apply b, then a), identified by comparing with every element "
              f"on all {sum(4 ** n for n in range(n_max + 1))} words of length 0..{n_max}")
    table = {}
    for a in elements:
        for b in elements:
            table[(a, b)] = identify(_compose(a, b), elements, n_max)
    rows = [[a] + [table[(a, b)] or "∉G" for b in elements] for a in elements]
    trace.add("calculation", "Cayley table (row ∘ column)", header=list(elements), rows=rows)
    trace.add("limitation", f"equality is checked exhaustively up to length {n_max}; the general-length argument "
              "(R and C are commuting involutions) is in docs/unit3_expected.md §1")
    return table, trace


def check_axioms(elements, table):
    elements = list(elements)
    trace = Trace()
    report = {}
    closure = next(((a, b) for a in elements for b in elements if table[(a, b)] not in elements), None)
    report["closure"] = {"holds": closure is None, "witness": closure}
    assoc = None
    if closure is None:
        assoc = next(((a, b, c) for a in elements for b in elements for c in elements
                      if table[(table[(a, b)], c)] != table[(a, table[(b, c)])]), None)
    report["associativity"] = {"holds": closure is None and assoc is None, "witness": assoc,
                               "triples": len(elements) ** 3}
    identity = next((e for e in elements if all(table[(e, a)] == a == table[(a, e)] for a in elements)), None)
    report["identity"] = {"holds": identity is not None, "element": identity}
    no_inverse = None
    inverses = {}
    if identity is not None:
        for a in elements:
            inv = next((b for b in elements if table[(a, b)] == identity == table[(b, a)]), None)
            if inv is None and no_inverse is None:
                no_inverse = a
            inverses[a] = inv
    report["inverses"] = {"holds": identity is not None and no_inverse is None, "witness": no_inverse,
                          "inverse_of": inverses}
    comm = next(((a, b) for a in elements for b in elements if table[(a, b)] != table[(b, a)]), None)
    report["commutative"] = {"holds": comm is None, "witness": comm}
    for name in ("closure", "associativity", "identity", "inverses", "commutative"):
        r = report[name]
        trace.add("algorithm", f"{name}: " + ("holds" if r["holds"] else f"fails (witness {r.get('witness')})"),
                  **{k: v for k, v in r.items() if k != "inverse_of"})
    report["group"] = all(report[p]["holds"] for p in ("closure", "associativity", "identity", "inverses"))
    report["monoid"] = all(report[p]["holds"] for p in ("closure", "associativity", "identity"))
    kind = "abelian group" if report["group"] and report["commutative"]["holds"] else (
        "group" if report["group"] else ("monoid" if report["monoid"] else "neither"))
    trace.add("calculation", f"structure: {kind} of order {len(elements)}", kind=kind)
    return report, trace


def element_orders(elements, table, identity="I"):
    orders = {}
    for a in elements:
        power, k = a, 1
        while power != identity:
            power = table[(a, power)]
            k += 1
            if k > len(elements):
                raise ValueError(f"{a} has no finite order within the table")
        orders[a] = k
    cyclic = any(v == len(elements) for v in orders.values())
    trace = Trace().add("calculation", "element orders: " + ", ".join(f"|{a}| = {k}" for a, k in orders.items())
                        + f"; cyclic: {cyclic}", orders=orders, cyclic=cyclic)
    return orders, trace


def check_isomorphism(elements, table, phi, target_op):
    failure = next(((a, b) for a in elements for b in elements
                    if phi[table[(a, b)]] != target_op(phi[a], phi[b])), None)
    bijective = len(set(phi[e] for e in elements)) == len(elements)
    holds = failure is None and bijective
    trace = Trace().add("algorithm", f"φ(a∘b) = φ(a)·φ(b) on all {len(elements) ** 2} pairs and φ bijective: {holds}",
                        witness=failure, bijective=bijective)
    return holds, trace


def subgroups(elements, table, identity="I"):
    elements = list(elements)
    found = []
    for mask in range(1, 2 ** len(elements)):
        h = [e for i, e in enumerate(elements) if mask >> i & 1]
        if identity in h and all(table[(a, b)] in h for a in h for b in h):
            found.append(h)  # finite + closed + identity is enough: each element's inverse is a power of it
    trace = Trace().add("algorithm", f"tested all {2 ** len(elements) - 1} non-empty subsets; {len(found)} subgroups",
                        subgroups=found)
    return found, trace


def left_cosets(elements, table, subgroup):
    cosets = []
    for g in elements:
        c = sorted({table[(g, h)] for h in subgroup}, key=list(elements).index)
        if c not in cosets:
            cosets.append(c)
    index = len(elements) // len(subgroup)
    trace = Trace().add("calculation", f"left cosets of {{{', '.join(subgroup)}}}: "
                        + ", ".join("{" + ", ".join(c) + "}" for c in cosets)
                        + f"; Lagrange index |G|/|H| = {index}", cosets=cosets, lagrange=len(cosets) == index)
    return cosets, trace


def fixed_count(name, n):
    if name == "I":
        value, why = 4 ** n, "every string"
    elif name == "R":
        value, why = 4 ** ((n + 1) // 2), f"first ⌈n/2⌉ = {(n + 1) // 2} letters free"
    elif name == "C":
        value, why = (1, "only the empty word") if n == 0 else (0, "no base is its own complement")
    elif name == "RC":
        if n % 2:
            value, why = 0, "odd n: the middle base would equal its own complement"
        else:
            value, why = 4 ** (n // 2), f"even n: first n/2 = {n // 2} letters free"
    else:
        raise ValueError(f"unknown transformation {name!r}")
    trace = Trace().add("calculation", f"|Fix({name})| on Σ^{n} = {value} ({why})", value=value)
    return value, trace


def _check_enumerable(n):
    if n > ENUMERATION_LIMIT:
        raise ValueError(f"exhaustive enumeration of Σ^{n} ({4 ** n} strings) exceeds the limit n ≤ {ENUMERATION_LIMIT}")


def fixed_count_enumerated(name, n):
    _check_enumerable(n)
    f = TRANSFORMS[name]
    value = sum(1 for w in words(n, n) if f(w) == w)
    return value, Trace().add("calculation", f"counted {value} strings of length {n} with {name}(s) = s")


def burnside(elements, n):
    terms = {g: fixed_count(g, n)[0] for g in elements}
    total = sum(terms.values())
    if total % len(elements):
        raise ArithmeticError("Burnside sum not divisible by |G|; the element set is not a group action")
    value = total // len(elements)
    trace = Trace().add("calculation", f"orbits = (" + " + ".join(str(terms[g]) for g in elements)
                        + f") / {len(elements)} = {value}", terms=terms, value=value)
    return value, trace


def orbit(s, elements):
    images = {g: TRANSFORMS[g](s) for g in elements}
    members = sorted(set(images.values()))
    stabiliser = [g for g in elements if images[g] == s]
    trace = Trace()
    trace.add("calculation", f"orbit of {s}: {members} (size {len(members)}); stabiliser {stabiliser}",
              orbit=members, stabiliser=stabiliser)
    trace.add("calculation", f"orbit–stabiliser: {len(members)} · {len(stabiliser)} = {len(elements)}",
              holds=len(members) * len(stabiliser) == len(elements))
    return {"orbit": members, "stabiliser": stabiliser}, trace


def canonical(s, elements=Z2):
    rep = min(TRANSFORMS[g](s) for g in elements)
    return rep, Trace().add("algorithm", f"canonical representative = lexicographically smallest orbit member = {rep}")


def orbits_enumerated(elements, n):
    _check_enumerable(n)
    seen, sizes, stab_ok = set(), {}, True
    count = 0
    for w in words(n, n):
        if w in seen:
            continue
        o, _ = orbit(w, elements)
        seen.update(o["orbit"])
        count += 1
        size = len(o["orbit"])
        sizes[size] = sizes.get(size, 0) + 1
        stab_ok &= size * len(o["stabiliser"]) == len(elements)
    trace = Trace().add("calculation", f"enumerated Σ^{n}: {count} orbits; sizes {dict(sorted(sizes.items()))}; "
                        f"orbit–stabiliser holds for every orbit: {stab_ok}", sizes=sizes, orbit_stabiliser=stab_ok)
    return {"count": count, "sizes": dict(sorted(sizes.items())), "orbit_stabiliser": stab_ok}, trace


def dataset_orbits(seqs, elements=Z2):
    by_rep = {}
    for sid, s in seqs.items():
        by_rep.setdefault(canonical(s, elements)[0], []).append(sid)
    classes = list(by_rep.values())
    trace = Trace()
    trace.add("calculation", f"{len(seqs)} sequences fall into {len(classes)} classes under {{{', '.join(elements)}}}",
              classes=classes)
    trace.add("interpretation", "a class holds sequences that are the same molecule written from different strands; "
              "a singleton only means its partner is absent from the data")
    return classes, trace


def permutation_cycles(mapping):
    seen, cycles = set(), []
    for start in mapping:
        if start in seen:
            continue
        cyc, x = [], start
        while x not in seen:
            seen.add(x)
            cyc.append(x)
            x = mapping[x]
        if len(cyc) > 1:
            cycles.append(tuple(cyc))
    order = 1
    for c in cycles:
        a, b = order, len(c)
        while b:
            a, b = b, a % b
        order = order * len(c) // a  # lcm via Euclid
    transpositions = sum(len(c) - 1 for c in cycles)
    parity = "even" if transpositions % 2 == 0 else "odd"
    text = "".join("(" + " ".join(map(str, c)) + ")" for c in cycles) or "identity"
    trace = Trace().add("calculation", f"{text}: order {order}, {transpositions} transpositions, {parity}",
                        cycles=cycles, order=order, parity=parity)
    return {"cycles": cycles, "order": order, "parity": parity}, trace


def substitution_monoid(n):
    # elements are patterns over {·,A,C,G,T}: keep the position, or write that base
    _check_enumerable(n)
    domain = list(words(n, n))

    def as_fn(p):
        return lambda x: "".join(x[i] if p[i] == "·" else p[i] for i in range(n))

    def table_of(fn):
        return tuple(fn(x) for x in domain)

    def pattern_of(images):
        p = []
        for i in range(n):
            col = {img[i] for img in images}
            kept = all(img[i] == x[i] for img, x in zip(images, domain))
            p.append("·" if kept else col.pop())
        return "".join(p)

    identity = "·" * n
    gens = [identity[:i] + b + identity[i + 1:] for i in range(n) for b in ALPHABET]
    elements, frontier = {identity}, [identity]
    while frontier:  # closure under composition with the generators
        nxt = []
        for p in frontier:
            for g in gens:
                images = table_of(lambda x, p=p, g=g: as_fn(p)(as_fn(g)(x)))
                q = pattern_of(images)
                if q not in elements:
                    elements.add(q)
                    nxt.append(q)
        frontier = nxt
    elements = sorted(elements)
    fns = {p: table_of(as_fn(p)) for p in elements}
    index = {v: k for k, v in fns.items()}
    pos = {x: i for i, x in enumerate(domain)}
    table = {(p, q): index[tuple(fns[p][pos[y]] for y in fns[q])] for p in elements for q in elements}
    report, trace = check_axioms(elements, table)
    idempotents = [p for p in elements if table[(p, p)] == p]
    units = [p for p in elements if report["inverses"]["inverse_of"].get(p)]
    witness = None
    for p in elements:
        if len(set(fns[p])) < len(domain):
            img = fns[p]
            i = next(i for i in range(len(domain)) for j in range(i) if img[i] == img[j])
            j = next(j for j in range(i) if img[i] == img[j])
            witness = (p, domain[j], domain[i], img[i])
            break
    trace.add("calculation", f"generated by {len(gens)} maps 'set position i to b': {len(elements)} elements "
              f"(5^{n} = {5 ** n}); idempotents {len(idempotents)}; units {units}",
              size=len(elements), idempotents=len(idempotents), units=units)
    if witness:
        trace.add("algorithm", f"not injective: f_{witness[0]}({witness[1]}) = f_{witness[0]}({witness[2]}) = "
                  f"{witness[3]}, so it has no inverse", witness=witness)
    trace.add("interpretation", "substitution operators form a monoid, not a group: only the identity is invertible")
    return {"size": len(elements), "report": report, "idempotents": len(idempotents), "units": units,
            "witness": witness, "table": table, "elements": elements}, trace


def concatenation_laws(max_len=3):
    pool = list(words(max_len))
    checks = {
        "C(xy) = C(x)C(y)": lambda x, y: complement(x + y) == complement(x) + complement(y),
        "R(xy) = R(y)R(x)": lambda x, y: _reverse(x + y) == _reverse(y) + _reverse(x),
        "RC(xy) = RC(y)RC(x)": lambda x, y: reverse_complement(x + y) == reverse_complement(y) + reverse_complement(x),
    }
    result = {name: all(f(x, y) for x in pool for y in pool) for name, f in checks.items()}
    trace = Trace().add("calculation", f"checked on all {len(pool) ** 2} pairs of words of length ≤ {max_len}: "
                        + "; ".join(f"{k}: {v}" for k, v in result.items()), **{"pairs": len(pool) ** 2})
    return result, trace


# decision D9: a second Klein four-group, acting on the four bases rather than on strings
PURINES, PYRIMIDINES = ("A", "G"), ("C", "T")


def _pairing(*pairs):
    m = {}
    for x, y in pairs:
        m[x], m[y] = y, x
    return {b: m.get(b, b) for b in ALPHABET}


BASE_PERMS = {
    "id": _pairing(),
    "(A T)(C G)": _pairing(("A", "T"), ("C", "G")),  # letterwise this is the string map C
    "(A G)(C T)": _pairing(("A", "G"), ("C", "T")),  # purine<->purine, pyrimidine<->pyrimidine: the transitions
    "(A C)(G T)": _pairing(("A", "C"), ("G", "T")),
}
BASE_K4 = tuple(BASE_PERMS)
TRANSITION_ELEMENT = "(A G)(C T)"


def fixed_point_free_involutions():
    # a fixed-point-free involution of a 4-set is a perfect matching; there are 3!! = 3
    # pair the first base with each of the other three; pairing (a,b)+(c,d) and (c,d)+(a,b) is the same matching
    first, others = ALPHABET[0], ALPHABET[1:]
    found = []
    for b in others:
        rest = [x for x in others if x != b]
        found.append(_pairing((first, b), (rest[0], rest[1])))
    names = [n for n, p in BASE_PERMS.items() if p in found]
    trace = Trace().add("calculation", f"{len(found)} fixed-point-free involutions of {ALPHABET} (3!! perfect "
                        f"matchings): {names}", count=len(found), names=names)
    return found, trace


def base_cayley_table():
    def identify_perm(m):
        return next((n for n, p in BASE_PERMS.items() if p == m), None)

    table = {}
    for a in BASE_K4:
        for b in BASE_K4:
            table[(a, b)] = identify_perm({x: BASE_PERMS[a][BASE_PERMS[b][x]] for x in ALPHABET})
    rows = [[a] + [table[(a, b)] or "∉K" for b in BASE_K4] for a in BASE_K4]
    trace = Trace().add("algorithm", "entry (a,b) = a∘b on bases, identified by comparing the whole permutation",
                        header=list(BASE_K4), rows=rows)
    trace.add("interpretation", "this group acts on the four bases; the Unit III group {I,R,C,RC} acts on strings, "
              "and the string map C is the letterwise extension of (A T)(C G)")
    return table, trace


def classify_substitution(x, y):
    if x == y:
        raise ValueError("a substitution needs two different bases")
    for b in (x, y):
        if b not in ALPHABET:
            raise SequenceError(f"base {b!r} is outside {ALPHABET}")
    g = next(n for n in BASE_K4 if n != "id" and BASE_PERMS[n][x] == y)
    kind = "transition" if g == TRANSITION_ELEMENT else "transversion"
    same_class = (x in PURINES) == (y in PURINES)
    trace = Trace().add("calculation", f"{x}<->{y} is realised by {g}, so it is a {kind}", element=g, kind=kind)
    trace.add("calculation", f"both in the same base class (purine/pyrimidine): {same_class}", consistent=same_class == (kind == "transition"))
    return {"element": g, "kind": kind, "same_class": same_class}, trace


def ordered_pair_coverage():
    # each non-identity element realises exactly 4 ordered pairs; 3 x 4 = 12 = every ordered pair
    by_element = {n: sorted((x, BASE_PERMS[n][x]) for x in ALPHABET) for n in BASE_K4 if n != "id"}
    seen = [p for pairs in by_element.values() for p in pairs]
    all_ordered = [(x, y) for x in ALPHABET for y in ALPHABET if x != y]
    result = {"per_element": {n: len(p) for n, p in by_element.items()}, "total": len(seen),
              "distinct": len(set(seen)), "covers_all": sorted(set(seen)) == sorted(all_ordered)}
    trace = Trace().add("calculation", f"{result['per_element']}; {result['total']} ordered pairs, all distinct: "
                        f"{result['distinct'] == result['total']}; covers all 12: {result['covers_all']}", **result)
    return result, trace


def substitution_census(pairs):
    # pairs: iterable of (seq_a, seq_b) of equal length; every differing position contributes one unordered substitution
    counts = {n: 0 for n in BASE_K4 if n != "id"}
    by_pair, positions = {}, 0
    for a, b in pairs:
        require_same_length(a, b)
        for x, y in zip(a, b):
            if x != y:
                g = next(n for n in counts if BASE_PERMS[n][x] == y)
                counts[g] += 1
                key = "".join(sorted((x, y)))
                by_pair[key] = by_pair.get(key, 0) + 1
                positions += 1
    ts = counts[TRANSITION_ELEMENT]
    tv = positions - ts
    result = {"substitutions": positions, "by_element": counts, "by_base_pair": dict(sorted(by_pair.items())),
              "transitions": ts, "transversions": tv,
              "ratio": (ts / tv if tv else None), "baseline_ratio": 0.5}
    trace = Trace()
    trace.add("observation", f"{positions} substitutions over the supplied sequence pairs", count=positions)
    trace.add("calculation", f"by group element {counts}; transitions {ts}, transversions {tv}"
              + (f", ratio {ts / tv:.3f}" if tv else ""), **{k: v for k, v in result.items() if k != "by_base_pair"})
    trace.add("calculation", "combinatorial baseline: 2 of the 6 unordered base pairs are transitions, so equally "
              "likely substitutions would give ratio 0.5")
    trace.add("limitation", "counts are direction-free: two observed sequences differ at a position without showing "
              "which base came first. No claim is made about what ratio real data should show; the 1:2 baseline is a "
              "counting statement about the alphabet only")
    return result, trace
