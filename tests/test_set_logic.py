from itertools import combinations, product

import networkx as nx
import pytest

from conftest import HAPS_A, assert_trace
from genomath._core import mismatch_count, reverse_complement
from genomath import set_logic as sl

# docs/synthetic_expected.md §1
M_STAR_A = {
    "syn01": set(), "syn02": {(2, "C")}, "syn03": {(2, "C"), (5, "G")}, "syn04": {(2, "C"), (5, "G"), (9, "T")},
    "syn05": {(11, "A")}, "syn06": {(2, "C"), (11, "A")}, "syn07": {(5, "C"), (7, "A")},
    "syn08": {(8, "C"), (11, "G")}, "syn09": {(1, "C"), (3, "G"), (4, "C"), (6, "A"), (10, "G"), (12, "A")},
}


def test_mutation_sets_match_hand_table(dataset_a):
    R, seqs = dataset_a.reference.seq, dataset_a.seqs()
    for sid, expected in M_STAR_A.items():
        m, t1 = sl.mutation_set(R, seqs[sid])
        ms, t2 = sl.allele_mutation_set(R, seqs[sid])
        assert ms == expected
        assert m == {i for i, _ in expected}
        assert_trace(t1)
        assert_trace(t2)


def test_set_operations_hand_example(dataset_a):
    # §2: M3 = {2,5}, M7 = {5,7}
    r, trace = sl.set_operations(frozenset({2, 5}), frozenset({5, 7}))
    assert r["intersection"] == {5} and r["union"] == {2, 5, 7}
    assert r["a_minus_b"] == {2} and r["b_minus_a"] == {7} and r["symmetric_difference"] == {2, 7}
    assert any(s["data"].get("holds") is True for s in trace)
    # allele level: no shared mutation
    r2, _ = sl.set_operations(frozenset(M_STAR_A["syn03"]), frozenset(M_STAR_A["syn07"]))
    assert r2["intersection"] == set()


def test_set_operations_against_python_sets(dataset_b):
    R, seqs = dataset_b.reference.seq, dataset_b.seqs()
    sets = [sl.allele_mutation_set(R, s)[0] for s in seqs.values()]
    for a, b in combinations(sets, 2):
        r, _ = sl.set_operations(a, b)
        assert r["intersection"] == a & b and r["union"] == a | b
        assert r["a_minus_b"] == a - b and r["b_minus_a"] == b - a and r["symmetric_difference"] == a ^ b


def test_empty_set_operations():
    r, _ = sl.set_operations(frozenset(), frozenset())
    assert r["union"] == set() and r["size_a"] == 0


@pytest.mark.parametrize("which", ["a", "b"])
def test_distance_decomposition_on_all_pairs(dataset_a, dataset_b, which):
    data = dataset_a if which == "a" else dataset_b
    R, seqs = data.reference.seq, data.seqs()
    for s, t in combinations(seqs.values(), 2):
        r, _ = sl.distance_decomposition(R, s, t)
        assert r["holds"] and r["d_H"] == mismatch_count(s, t)


def test_distance_decomposition_hand_example(dataset_a):
    seqs = dataset_a.seqs()
    r, _ = sl.distance_decomposition(seqs["syn01"], seqs["syn03"], seqs["syn07"])
    assert (r["sym_diff"], r["shared_differing"], r["d_H"]) == ([2, 7], [5], 3)


def test_threshold_relation_properties(dataset_a):
    seqs = dataset_a.seqs()
    pairs, _ = sl.threshold_relation(seqs, 1)
    report, trace = sl.relation_properties(HAPS_A, pairs)
    assert report["reflexive"]["holds"] and report["symmetric"]["holds"]
    assert not report["transitive"]["holds"]
    assert report["transitive"]["witness"] == ("syn01", "syn02", "syn03")
    assert not report["equivalence"]
    assert_trace(trace)
    with pytest.raises(ValueError, match="transitive"):
        sl.equivalence_classes(HAPS_A, pairs)


def _all_relations(elements):
    cells = [(a, b) for a in elements for b in elements]
    for bits in range(2 ** len(cells)):
        yield {c for k, c in enumerate(cells) if bits >> k & 1}


def test_relation_properties_exhaustive_on_three_elements():
    # every one of the 512 relations on {0,1,2}, against the definitions written out directly
    E = [0, 1, 2]
    for R in _all_relations(E):
        r, _ = sl.relation_properties(E, R)
        assert r["reflexive"]["holds"] == all((a, a) in R for a in E)
        assert r["symmetric"]["holds"] == all((b, a) in R for a, b in R)
        assert r["antisymmetric"]["holds"] == all(a == b or (b, a) not in R for a, b in R)
        assert r["transitive"]["holds"] == all((a, d) in R for a, b in R for c, d in R if b == c)
        w = r["transitive"]["witness"]
        if w:
            a, b, c = w
            assert (a, b) in R and (b, c) in R and (a, c) not in R


def test_closure_exhaustive_on_three_elements():
    E = [0, 1, 2]
    for R in _all_relations(E):
        closure, _ = sl.transitive_closure(E, R)
        g = nx.DiGraph(list(R))
        g.add_nodes_from(E)
        assert closure == set(nx.transitive_closure(g, reflexive=False).edges())


def test_closure_needs_last_pivot():
    # the only route a -> b passes through c, which is last in element order
    closure, _ = sl.transitive_closure(["a", "b", "c"], {("a", "c"), ("c", "b")})
    assert closure == {("a", "c"), ("c", "b"), ("a", "b")}


def test_non_symmetric_relations():
    E = range(1, 7)
    divides = {(a, b) for a in E for b in E if b % a == 0}
    r, _ = sl.relation_properties(E, divides)
    assert r["reflexive"]["holds"] and r["antisymmetric"]["holds"] and r["transitive"]["holds"]
    assert r["symmetric"]["witness"] == (1, 2)
    chain = {(1, 2), (2, 3)}
    assert sl.relation_properties([1, 2, 3], chain)[0]["transitive"]["witness"] == (1, 2, 3)


def test_subset_order_on_mutation_sets_is_partial_order(dataset_a):
    # (M*, ⊆) on the haplotypes: reflexive, antisymmetric, transitive, not symmetric
    R, seqs = dataset_a.reference.seq, dataset_a.seqs()
    M = {h: sl.allele_mutation_set(R, s)[0] for h, s in seqs.items()}
    subset = {(a, b) for a in HAPS_A for b in HAPS_A if M[a] <= M[b]}
    r, _ = sl.relation_properties(HAPS_A, subset)
    assert r["reflexive"]["holds"] and r["antisymmetric"]["holds"] and r["transitive"]["holds"]
    assert not r["symmetric"]["holds"]
    assert {("syn02", "syn03"), ("syn03", "syn04"), ("syn02", "syn04")} <= subset


def test_r0_on_haplotypes_is_equality(dataset_a):
    pairs, _ = sl.threshold_relation(dataset_a.seqs(), 0)
    assert pairs == {(h, h) for h in HAPS_A}
    assert sl.relation_properties(HAPS_A, pairs)[0]["equivalence"]


def test_identity_relation_classes_are_haplotypes(dataset_a):
    seqs = {r.id: r.seq for r in dataset_a.records}
    pairs, _ = sl.identity_relation(seqs)
    classes, _ = sl.equivalence_classes(list(seqs), pairs)
    assert sorted(map(tuple, classes)) == sorted(h.members for h in dataset_a.haplotypes)


def test_orientation_relation(dataset_a):
    seqs = {r.id: r.seq for r in dataset_a.records}
    pairs, _ = sl.orientation_relation(seqs)
    classes, _ = sl.equivalence_classes(list(seqs), pairs)
    assert ["syn03", "syn08", "syn10"] in classes
    assert ["syn01"] in classes
    assert len(classes) == 8


@pytest.mark.parametrize("which,ks", [("a", range(0, 10)), ("b", [0, 1, 2, 5, 10])])
def test_closure_of_rk_equals_graph_components(dataset_a, dataset_b, which, ks):
    data = dataset_a if which == "a" else dataset_b
    seqs = data.seqs()
    ids = list(seqs)
    for k in ks:
        pairs, _ = sl.threshold_relation(seqs, k)
        closure, _ = sl.transitive_closure(ids, pairs)
        g = nx.DiGraph(list(pairs))
        g.add_nodes_from(ids)
        assert closure == set(nx.transitive_closure(g, reflexive=False).edges())
        classes, _ = sl.equivalence_classes(ids, closure)
        und = nx.Graph([(a, b) for a, b in pairs if a != b])
        und.add_nodes_from(ids)
        assert sorted(map(sorted, classes)) == sorted(map(sorted, nx.connected_components(und)))


def test_component_counts_dataset_a(dataset_a):
    # §4: k -> number of components
    seqs = dataset_a.seqs()
    for k, expected in {0: 9, 1: 4, 2: 2, 5: 2, 6: 1}.items():
        closure, _ = sl.transitive_closure(HAPS_A, sl.threshold_relation(seqs, k)[0])
        assert len(sl.equivalence_classes(HAPS_A, closure)[0]) == expected


def test_predicates_hand_example(dataset_a):
    # §2: P = d_H(R,S) ≤ 2, Q = S_2 = C
    R, seqs = dataset_a.reference.seq, dataset_a.seqs()
    P, Q = sl.within_distance(R, 2), sl.base_at(2, "C")
    assert sl.filter_sequences(seqs, P & Q, (P, Q))[0] == ["syn02", "syn03", "syn06"]
    assert sl.filter_sequences(seqs, P & ~Q)[0] == ["syn01", "syn05", "syn07", "syn08"]
    assert sl.filter_sequences(seqs, ~P)[0] == ["syn04", "syn09"]
    assert sl.filter_sequences(seqs, P | Q)[0] == [h for h in HAPS_A if h not in ("syn09",)]


def test_mutation_predicates(dataset_a):
    R, seqs = dataset_a.reference.seq, dataset_a.seqs()
    shares5 = sl.shares_mutation_with(R, seqs["syn03"], 5)
    assert not shares5(seqs["syn07"])  # same position, different base
    assert shares5(seqs["syn04"])
    assert sl.carries_mutation(R, 2, "C")(seqs["syn06"])
    with pytest.raises(ValueError):
        sl.carries_mutation(R, 2, "A")
    with pytest.raises(ValueError):
        sl.base_at(0, "A")


def test_function_properties():
    words = ["".join(p) for p in product("ACGT", repeat=3)]
    r, trace = sl.function_properties(words, reverse_complement, set(words), "RC")
    assert r["bijective"]
    assert_trace(trace)
    # S: [n] -> Sigma is surjective exactly when all four bases occur; 24 such strings for n = 4 (§8)
    n4 = ["".join(p) for p in product("ACGT", repeat=4)]
    surj = [w for w in n4 if sl.function_properties(range(1, 5), lambda i, w=w: w[i - 1], "ACGT")[0]["surjective"]]
    assert len(surj) == 24


def test_distance_to_reference_not_injective(dataset_a):
    R, seqs = dataset_a.reference.seq, dataset_a.seqs()
    r, _ = sl.function_properties(HAPS_A, lambda h: mismatch_count(R, seqs[h]), range(13), "d_R")
    assert not r["injective"] and r["collision"][2] == 1  # syn02 and syn05 are both at distance 1
