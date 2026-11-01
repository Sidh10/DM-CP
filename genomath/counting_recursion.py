"""Unit II: sequence-space counts, M(n,k), its recurrence, bounded recursive generation and inclusion-exclusion."""
from ._core import ALPHABET, SequenceError, mismatch_positions
from .trace import Trace

DEFAULT_CAP = 100_000  # provisional (Stage 0 open decision 7)


class GenerationLimitError(ValueError):
    def __init__(self, message, count, trace):
        super().__init__(message)
        self.count, self.trace = count, trace


def _factorial(n):
    if n < 0:
        raise ValueError("factorial is undefined for negative n")
    out = 1
    for i in range(2, n + 1):
        out *= i
    return out


def _binom(n, k):
    if n < 0:
        raise ValueError("C(n,k) needs n >= 0")
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    c = 1
    for i in range(1, k + 1):
        c = c * (n - k + i) // i  # exact: after step i, c = C(n-k+i, i)
    return c


def _m(n, k):
    return _binom(n, k) * 3 ** k


def factorial(n):
    value = _factorial(n)
    trace = Trace().add("calculation", f"{n}! = " + (" · ".join(str(i) for i in range(1, n + 1)) or "1 (empty product)")
                        + f" = {value}")
    return value, trace


def binomial(n, k):
    trace = Trace()
    if k < 0 or k > n:
        trace.add("calculation", f"C({n},{k}) = 0 because k is outside 0..n")
        return 0, trace
    j = min(k, n - k)
    if j != k:
        trace.add("calculation", f"symmetry: C({n},{k}) = C({n},{j})")
    c, steps = 1, []
    for i in range(1, j + 1):
        c = c * (n - j + i) // i
        steps.append(f"step {i}: × {n - j + i} / {i} -> {c}")
    trace.add_steps("calculation", steps)
    trace.add("calculation", f"C({n},{k}) = {c}", value=c)
    return c, trace


def permutations(n, k):
    if k < 0 or k > n:
        return 0, Trace().add("calculation", f"P({n},{k}) = 0 because k is outside 0..n")
    factors = list(range(n, n - k, -1))
    value = 1
    for f in factors:
        value *= f
    trace = Trace().add("calculation", f"P({n},{k}) = n!/(n-k)! = " + (" · ".join(map(str, factors)) or "1")
                        + f" = {value}", value=value)
    return value, trace


def sequence_space(n):
    value = 4 ** n
    trace = Trace()
    trace.add("calculation", f"product rule: 4 choices at each of {n} positions -> 4^{n} = {value}",
              value=value, digits=len(str(value)))
    return value, trace


def growth_table(ns):
    rows = [{"n": n, "4^n": 4 ** n, "digits": len(str(4 ** n))} for n in ns]
    trace = Trace().add("calculation", "4^n for each n; the count gains about 0.602 decimal digits per base", rows=rows)
    trace.add("interpretation", "exhaustive enumeration of Σ^n is only feasible for small n")
    return rows, trace


def mutation_count(n, k):
    trace = Trace()
    if k < 0 or k > n:
        trace.add("calculation", f"M({n},{k}) = 0: cannot change {k} of {n} positions")
        return 0, trace
    c, value = _binom(n, k), _m(n, k)
    trace.add("calculation", f"M({n},{k}) = C({n},{k}) · 3^{k} = {c} · {3 ** k} = {value}", value=value)
    trace.add("limitation", "counts substitutions only; insertions and deletions change the length and are excluded")
    return value, trace


def recurrence_table(n_max):
    # bottom-up dynamic programme: row n is built only from row n-1, O(n_max^2) time
    table = [[1]]
    for n in range(1, n_max + 1):
        prev = table[-1]
        row = []
        for k in range(n + 1):
            keep = prev[k] if k < len(prev) else 0
            change = 3 * prev[k - 1] if k >= 1 else 0
            row.append(keep + change)
        table.append(row)
    trace = Trace()
    trace.add("algorithm", "base case M(0,0) = 1; M(n,k) = M(n-1,k) + 3·M(n-1,k-1), with 0 outside 0 ≤ k ≤ n")
    trace.add("calculation", f"table for n = 0..{n_max}", table=table)
    return table, trace


def mutation_count_recursive(n, k):
    # top-down with an explicit memo dict so the trace can report how many distinct subproblems were solved
    memo, calls = {}, [0]

    def m(a, b):
        calls[0] += 1
        if b < 0 or b > a:
            return 0
        if a == 0:
            return 1
        if (a, b) not in memo:
            memo[(a, b)] = m(a - 1, b) + 3 * m(a - 1, b - 1)
        return memo[(a, b)]

    if n > 900:
        raise ValueError("recursion depth would exceed Python's limit; use recurrence_table or mutation_count")
    value = m(n, k)
    trace = Trace()
    trace.add("algorithm", f"recursive M({n},{k}) with memoisation: {calls[0]} calls, {len(memo)} distinct subproblems",
              calls=calls[0], subproblems=len(memo))
    trace.add("calculation", f"M({n},{k}) = {value}", value=value)
    return value, trace


def verify_recurrence(n_max):
    table, _ = recurrence_table(n_max)
    mismatches = [(n, k) for n in range(n_max + 1) for k in range(n + 1) if table[n][k] != _m(n, k)]
    bad_rows = [n for n in range(n_max + 1) if sum(table[n]) != 4 ** n]
    trace = Trace()
    trace.add("calculation", f"closed form vs recurrence on {sum(n + 1 for n in range(n_max + 1))} cells: "
              f"{len(mismatches)} mismatches", mismatches=mismatches)
    trace.add("calculation", f"row sums equal 4^n for n = 0..{n_max}: {not bad_rows}", bad_rows=bad_rows)
    trace.add("limitation", "a finite check, not a proof; the induction proof is in docs/unit2_expected.md §2")
    return {"mismatches": mismatches, "bad_rows": bad_rows}, trace


def ball_size(n, r):
    terms = [_m(n, k) for k in range(0, min(r, n) + 1)]
    trace = Trace().add("calculation", f"|B({n},{r})| = " + " + ".join(map(str, terms)) + f" = {sum(terms)}",
                        value=sum(terms))
    return sum(terms), trace


def _check_ref(ref):
    bad = sorted(set(ref) - set(ALPHABET))
    if bad:
        raise SequenceError(f"reference contains non-ACGT characters {bad}")


def _guard(count, cap, what, trace):
    trace.add("calculation", f"{what} would produce {count} strings (cap {cap})", count=count, cap=cap)
    if count > cap:
        trace.add("limitation", f"refused: {count} exceeds the cap of {cap}; enumeration is Θ(n·count) time and memory")
        raise GenerationLimitError(f"{what}: {count} strings exceeds the cap of {cap}", count, trace)


def generate_mutants(ref, k, cap=DEFAULT_CAP):
    _check_ref(ref)
    if k < 0:
        raise ValueError("k must be >= 0")
    n = len(ref)
    trace = Trace()
    _guard(_m(n, k), cap, f"M({n},{k})", trace)
    out, calls = [], [0]
    current = list(ref)

    # choose positions in increasing order (recursion depth k, not n), then one of 3 alternative bases at each
    def rec(start, remaining):
        calls[0] += 1
        if remaining == 0:
            out.append("".join(current))
            return
        for i in range(start, n - remaining + 1):
            original = current[i]
            for b in ALPHABET:
                if b != original:
                    current[i] = b
                    rec(i + 1, remaining - 1)
            current[i] = original

    rec(0, k)
    trace.add("algorithm", f"recursion chose {k} positions in increasing order and 3 bases for each: "
              f"{calls[0]} calls, {len(out)} strings", calls=calls[0], generated=len(out))
    trace.add("calculation", f"Θ(n · M(n,k)) = Θ({n} · {len(out)}) character writes; recursion depth {k}")
    return out, trace


def generate_within(ref, r, cap=DEFAULT_CAP):
    _check_ref(ref)
    n = len(ref)
    trace = Trace()
    _guard(sum(_m(n, k) for k in range(min(r, n) + 1)), cap, f"B({n},{r})", trace)
    out = []
    for k in range(min(r, n) + 1):
        layer, _ = generate_mutants(ref, k, cap)
        out.extend(layer)
        trace.add("algorithm", f"distance {k}: {len(layer)} strings")
    return out, trace


def inclusion_exclusion_all_bases(n):
    # A_x = strings missing base x; |intersection of j such sets| = (4 - j)^n
    terms = []
    total = 0
    for j in range(0, 5):
        sign = -1 if j % 2 else 1
        term = sign * _binom(4, j) * (4 - j) ** n
        terms.append({"j": j, "sign": sign, "C(4,j)": _binom(4, j), "(4-j)^n": (4 - j) ** n, "term": term})
        total += term
    trace = Trace()
    trace.add("calculation", "U(n) = Σ_j (-1)^j C(4,j) (4-j)^n = " + " ".join(f"{t['term']:+d}" for t in terms)
              + f" = {total}", terms=terms, value=total)
    trace.add("interpretation", "U(n) counts strings using all four bases, i.e. surjective maps [n] -> Σ")
    return total, trace


def forced_repeat(length, k):
    windows, holes = length - k + 1, 4 ** k
    forced = windows > holes
    trace = Trace().add("calculation", f"{windows} windows (pigeons) vs 4^{k} = {holes} k-mers (holes): "
                        + ("repeat forced" if forced else "repeat not forced"), windows=windows, holes=holes)
    return forced, trace


def max_forced_k(length):
    k = 0
    while length - (k + 1) + 1 > 4 ** (k + 1):
        k += 1
    trace = Trace().add("calculation", f"largest k with L - k + 1 > 4^k for L = {length}: k = {k}", value=k)
    return k, trace


def find_repeated_kmer(seq, k):
    seen = {}
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        if w in seen:
            witness = (w, seen[w] + 1, i + 1)
            return witness, Trace().add("algorithm", f"{w} occurs at positions {witness[1]} and {witness[2]}",
                                        witness=witness)
        seen[w] = i
    return None, Trace().add("algorithm", f"no {k}-mer occurs twice in a sequence of length {len(seq)}")


def mutation_orderings(s, t, enumerate_limit=8):
    diff = mismatch_positions(s, t)
    k = len(diff)
    count = _factorial(k)
    trace = Trace()
    trace.add("calculation", f"d_H = {k}: a shortest path changes each of {diff} once, so it is an ordering of them; "
              f"{k}! = {count} paths", positions=diff, count=count)
    paths = None
    if k <= enumerate_limit:
        paths = []

        def rec(prefix, remaining):  # hand-written permutation generator
            if not remaining:
                paths.append(list(prefix))
                return
            for j, p in enumerate(remaining):
                prefix.append(p)
                rec(prefix, remaining[:j] + remaining[j + 1:])
                prefix.pop()

        rec([], diff)
        trace.add("algorithm", f"enumerated {len(paths)} orderings recursively")
    else:
        trace.add("limitation", f"k = {k} > {enumerate_limit}: orderings counted, not listed")

    def walk(order):
        cur, states = list(s), [s]
        for p in order:
            cur[p - 1] = t[p - 1]
            states.append("".join(cur))
        return states

    return {"k": k, "count": count, "orders": paths,
            "paths": [walk(o) for o in paths] if paths is not None else None}, trace
