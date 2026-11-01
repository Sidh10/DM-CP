import math
from itertools import product

import networkx as nx
import pytest

from conftest import assert_trace
from genomath import counting_recursion as cr
from genomath._core import SequenceError, mismatch_count

# docs/unit2_expected.md §2, transcribed by hand
HAND_M = [
    [1],
    [1, 3],
    [1, 6, 9],
    [1, 9, 27, 27],
    [1, 12, 54, 108, 81],
    [1, 15, 90, 270, 405, 243],
    [1, 18, 135, 540, 1215, 1458, 729],
    [1, 21, 189, 945, 2835, 5103, 5103, 2187],
    [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561],
]


def brute_by_distance(ref):
    out = {}
    for w in product("ACGT", repeat=len(ref)):
        w = "".join(w)
        out.setdefault(mismatch_count(ref, w), set()).add(w)
    return out


def test_binomial_factorial_permutations_against_math():
    for n in range(0, 61):
        assert cr.factorial(n)[0] == math.factorial(n)
        for k in range(-2, n + 3):
            assert cr.binomial(n, k)[0] == (math.comb(n, k) if 0 <= k <= n else 0)
            assert cr.permutations(n, k)[0] == (math.perm(n, k) if 0 <= k <= n else 0)
    with pytest.raises(ValueError):
        cr.factorial(-1)


def test_closed_form_recurrence_and_memo_agree_with_hand_table():
    table, trace = cr.recurrence_table(8)
    assert table == HAND_M
    assert_trace(trace)
    for n, row in enumerate(HAND_M):
        for k, v in enumerate(row):
            assert cr.mutation_count(n, k)[0] == v
            assert cr.mutation_count_recursive(n, k)[0] == v


def test_recurrence_boundaries():
    assert cr.mutation_count(0, 0)[0] == 1
    assert cr.mutation_count(5, 6)[0] == 0 and cr.mutation_count(5, -1)[0] == 0
    assert cr.mutation_count_recursive(5, 6)[0] == 0


def test_verify_recurrence_large():
    result, _ = cr.verify_recurrence(60)
    assert result == {"mismatches": [], "bad_rows": []}
    table, _ = cr.recurrence_table(60)
    assert all(table[n][k] == math.comb(n, k) * 3 ** k for n in range(61) for k in range(n + 1))


def test_memoisation_makes_recursion_polynomial():
    _, trace = cr.mutation_count_recursive(30, 15)
    data = next(s["data"] for s in trace if "subproblems" in s["data"])
    assert data["subproblems"] <= 31 * 16
    assert data["calls"] < 2000  # unmemoised recursion on (30,15) would make millions of calls


@pytest.mark.parametrize("ref", ["", "A", "GA", "GAT", "ACGT", "TTGCA", "GATTAC"])
def test_generator_equals_brute_force(ref):
    brute = brute_by_distance(ref)
    for k in range(0, len(ref) + 2):
        got, trace = cr.generate_mutants(ref, k)
        assert len(got) == len(set(got)) == cr.mutation_count(len(ref), k)[0]
        assert set(got) == brute.get(k, set())
        assert_trace(trace)
    ball, _ = cr.generate_within(ref, 2)
    assert set(ball) == set().union(*(brute.get(k, set()) for k in range(3)))


def test_guard_refuses_and_reports():
    with pytest.raises(cr.GenerationLimitError) as e:
        cr.generate_mutants("A" * 12, 6)
    assert e.value.count == 673_596
    assert "limitation" in e.value.trace.kinds()
    with pytest.raises(cr.GenerationLimitError):
        cr.generate_mutants("A" * 1701, 2)
    assert len(cr.generate_mutants("A" * 1701, 1)[0]) == 5103
    with pytest.raises(cr.GenerationLimitError):
        cr.generate_within("ACGT", 2, cap=10)
    assert len(cr.generate_within("ACGTACGT", 8)[0]) == 65536


def test_generator_rejects_bad_input():
    with pytest.raises(SequenceError):
        cr.generate_mutants("ACGN", 1)
    with pytest.raises(ValueError):
        cr.generate_mutants("ACGT", -1)


def test_ball_sizes():
    assert cr.ball_size(8, 2)[0] == 277
    assert cr.ball_size(12, 2)[0] == 631
    assert cr.ball_size(1701, 1)[0] == 5104
    assert cr.ball_size(3, 10)[0] == 64


def test_sequence_space_and_growth():
    assert cr.sequence_space(12)[0] == 16_777_216
    rows, _ = cr.growth_table([3, 1701])
    assert rows[0]["4^n"] == 64 and rows[1]["digits"] == 1025


def test_inclusion_exclusion_against_enumeration():
    hand = {1: 0, 2: 0, 3: 0, 4: 24, 5: 240, 6: 1560, 7: 8400, 8: 40824, 12: 14_676_024}
    for n, v in hand.items():
        assert cr.inclusion_exclusion_all_bases(n)[0] == v
    for n in range(0, 8):
        brute = sum(1 for w in product("ACGT", repeat=n) if set(w) == set("ACGT"))
        assert cr.inclusion_exclusion_all_bases(n)[0] == brute


def test_pigeonhole():
    assert cr.forced_repeat(1701, 5)[0] and not cr.forced_repeat(1701, 6)[0]
    assert cr.max_forced_k(1701)[0] == 5
    assert cr.find_repeated_kmer("ACGT", 1)[0] is None  # tight: 4 windows, 4 holes
    assert not cr.forced_repeat(4, 1)[0] and not cr.forced_repeat(17, 2)[0]  # windows == holes is not forced
    for length in range(1, 2001):
        k = cr.max_forced_k(length)[0]
        assert (k == 0 or cr.forced_repeat(length, k)[0]) and not cr.forced_repeat(length, k + 1)[0]
    assert cr.forced_repeat(5, 1)[0]
    for w in product("ACGT", repeat=5):
        assert cr.find_repeated_kmer("".join(w), 1)[0] is not None
    witness, _ = cr.find_repeated_kmer("GATTACATTA", 3)
    assert witness == ("ATT", 2, 7)


def test_pigeonhole_on_dataset_b(dataset_b):
    for r in dataset_b.records:
        w, _ = cr.find_repeated_kmer(r.seq, 5)
        assert w is not None and r.seq[w[1] - 1:w[1] + 4] == r.seq[w[2] - 1:w[2] + 4] == w[0]


def test_orderings_count_equals_shortest_paths_in_hamming_graph():
    # oracle: NetworkX shortest paths in the full Hamming graph H(4,4) (256 vertices)
    words = ["".join(w) for w in product("ACGT", repeat=4)]
    g = nx.Graph()
    g.add_edges_from((a, b) for i, a in enumerate(words) for b in words[i + 1:] if mismatch_count(a, b) == 1)
    for s, t in [("AAAA", "AAAA"), ("AAAA", "CAAA"), ("AAAA", "CGAA"), ("ACGT", "TGCA"), ("GATT", "GCTA")]:
        r, trace = cr.mutation_orderings(s, t)
        nx_paths = list(nx.all_shortest_paths(g, s, t))
        assert r["count"] == len(nx_paths) == math.factorial(mismatch_count(s, t))
        assert sorted(r["paths"]) == sorted(nx_paths)
        assert_trace(trace)


def test_orderings_dataset_a(dataset_a):
    seqs = dataset_a.seqs()
    r, _ = cr.mutation_orderings(seqs["syn01"], seqs["syn04"])
    assert r["count"] == 6
    observed = set(seqs.values())
    fully = [p for p in r["paths"] if all(x in observed for x in p)]
    assert fully == [[seqs[h] for h in ("syn01", "syn02", "syn03", "syn04")]]
    big, _ = cr.mutation_orderings("A" * 10, "C" * 10)
    assert big["count"] == 3_628_800 and big["paths"] is None
