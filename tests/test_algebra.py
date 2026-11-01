from itertools import permutations, product

import pytest
from sympy.combinatorics import Permutation, PermutationGroup

from conftest import assert_trace
from genomath import algebraic_structures as al
from genomath._core import SequenceError

# docs/unit3_expected.md, transcribed by hand
HAND_V4 = {"I": ["I", "R", "C", "RC"], "R": ["R", "I", "RC", "C"], "C": ["C", "RC", "I", "R"], "RC": ["RC", "C", "R", "I"]}
HAND_Z2 = {"I": ["I", "RC"], "RC": ["RC", "I"]}
FIX = {
    "I": [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536],
    "R": [1, 4, 4, 16, 16, 64, 64, 256, 256],
    "C": [1, 0, 0, 0, 0, 0, 0, 0, 0],
    "RC": [1, 0, 4, 0, 16, 0, 64, 0, 256],
}
ORBITS_Z2 = [1, 2, 10, 32, 136, 512, 2080, 8192, 32896]
ORBITS_V4 = [1, 2, 6, 20, 72, 272, 1056, 4160, 16512]
V4_SIZES = {1: {2: 2}, 2: {2: 4, 4: 2}, 3: {2: 8, 4: 12}, 4: {2: 16, 4: 56}}


@pytest.fixture(scope="module")
def v4_table():
    return al.cayley_table(al.V4)[0]


def sympy_group(elements, n):
    domain = list(al.words(n, n))
    idx = {w: i for i, w in enumerate(domain)}
    perms = [Permutation([idx[al.TRANSFORMS[g](w)] for w in domain]) for g in elements]
    return PermutationGroup(perms)


def test_cayley_tables_match_hand(v4_table):
    assert {a: [v4_table[(a, b)] for b in al.V4] for a in al.V4} == HAND_V4
    z2, trace = al.cayley_table(al.Z2)
    assert {a: [z2[(a, b)] for b in al.Z2] for a in al.Z2} == HAND_Z2
    assert_trace(trace)


def test_group_axioms(v4_table):
    r, trace = al.check_axioms(al.V4, v4_table)
    assert r["group"] and r["commutative"]["holds"] and r["identity"]["element"] == "I"
    assert r["inverses"]["inverse_of"] == {"I": "I", "R": "R", "C": "C", "RC": "RC"}
    assert r["associativity"]["triples"] == 64
    assert_trace(trace)
    for a in al.V4:  # Latin square
        assert sorted(v4_table[(a, b)] for b in al.V4) == sorted(al.V4)
        assert sorted(v4_table[(b, a)] for b in al.V4) == sorted(al.V4)


def test_closure_failure_detected():
    table, _ = al.cayley_table(("I", "R", "C"))
    r, _ = al.check_axioms(("I", "R", "C"), table)
    assert not r["closure"]["holds"] and r["closure"]["witness"] == ("R", "C") and not r["group"]


def test_non_associative_operation_detected():
    E = [0, 1, 2]
    table = {(a, b): (a - b) % 3 for a in E for b in E}
    r, _ = al.check_axioms(E, table)
    assert r["closure"]["holds"] and not r["associativity"]["holds"] and not r["identity"]["holds"]
    a, b, c = r["associativity"]["witness"]
    assert ((a - b) % 3 - c) % 3 != (a - (b - c) % 3) % 3


def test_orders_and_not_cyclic(v4_table):
    orders, _ = al.element_orders(al.V4, v4_table)
    assert orders == {"I": 1, "R": 2, "C": 2, "RC": 2}
    g = sympy_group(al.V4, 3)
    assert g.order() == 4 and g.is_abelian and not g.is_cyclic


def test_isomorphic_to_z2xz2_not_z4(v4_table):
    phi = {"I": (0, 0), "R": (1, 0), "C": (0, 1), "RC": (1, 1)}
    ok, _ = al.check_isomorphism(al.V4, v4_table, phi, lambda x, y: ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2))
    assert ok
    for images in permutations(range(4)):
        psi = dict(zip(al.V4, images))
        assert not al.check_isomorphism(al.V4, v4_table, psi, lambda x, y: (x + y) % 4)[0]


def test_subgroups_and_cosets(v4_table):
    subs, _ = al.subgroups(al.V4, v4_table)
    assert subs == [["I"], ["I", "R"], ["I", "C"], ["I", "RC"], ["I", "R", "C", "RC"]]
    cosets, trace = al.left_cosets(al.V4, v4_table, al.Z2)
    assert cosets == [["I", "RC"], ["R", "C"]]
    assert next(s for s in trace if "lagrange" in s["data"])["data"]["lagrange"]


def test_faithfulness_by_length():
    assert al.distinct_on(al.V4, 0)[0] == [["I", "R", "C", "RC"]]
    assert al.distinct_on(al.V4, 1)[0] == [["I", "R"], ["C", "RC"]]
    assert al.distinct_on(al.V4, 2)[0] == [["I"], ["R"], ["C"], ["RC"]]
    assert sympy_group(al.V4, 1).order() == 2 and sympy_group(al.V4, 2).order() == 4


def test_fixed_counts_formula_vs_enumeration():
    for g, row in FIX.items():
        for n, v in enumerate(row):
            assert al.fixed_count(g, n)[0] == v == al.fixed_count_enumerated(g, n)[0]
    assert al.fixed_count("RC", 12)[0] == 4096 and al.fixed_count("R", 12)[0] == 4096


def test_burnside_against_hand_enumeration_and_sympy():
    for n in range(9):
        assert al.burnside(al.Z2, n)[0] == ORBITS_Z2[n]
        assert al.burnside(al.V4, n)[0] == ORBITS_V4[n]
    assert al.burnside(al.Z2, 12)[0] == 8_390_656 and al.burnside(al.V4, 12)[0] == 4_196_352
    for n in range(0, 7):
        assert al.orbits_enumerated(al.Z2, n)[0]["count"] == ORBITS_Z2[n]
        assert al.orbits_enumerated(al.V4, n)[0]["count"] == ORBITS_V4[n]
    for n in range(1, 6):
        assert len(sympy_group(al.Z2, n).orbits()) == ORBITS_Z2[n]
        assert len(sympy_group(al.V4, n).orbits()) == ORBITS_V4[n]


def test_orbit_size_distribution():
    for n, sizes in V4_SIZES.items():
        r, _ = al.orbits_enumerated(al.V4, n)
        assert r["sizes"] == sizes and r["orbit_stabiliser"]
    z, _ = al.orbits_enumerated(al.Z2, 4)
    assert z["sizes"] == {1: 16, 2: 120}


def test_enumeration_limit():
    with pytest.raises(ValueError):
        al.orbits_enumerated(al.Z2, 11)


def test_dataset_a_orbits(dataset_a):
    s = dataset_a.seqs()
    classes, trace = al.dataset_orbits(s, al.Z2)
    assert len(classes) == 8 and ["syn03", "syn08"] in classes
    assert_trace(trace)
    o1, _ = al.orbit(s["syn01"], al.V4)
    assert o1["orbit"] == ["CTAATGCATTAG", "GATTACGTAATC"] and o1["stabiliser"] == ["I", "RC"]
    assert len(al.orbit(s["syn03"], al.V4)[0]["orbit"]) == 4 == len(al.orbit(s["syn09"], al.V4)[0]["orbit"])
    assert al.canonical(s["syn03"])[0] == al.canonical(s["syn08"])[0] == s["syn08"]
    assert [h for h, x in s.items() if al.TRANSFORMS["R"](x) == x] == []


def test_dataset_b_has_no_rc_structure(dataset_b):
    s = dataset_b.seqs()
    classes, _ = al.dataset_orbits(s, al.Z2)
    assert len(classes) == 30
    assert not any(al.TRANSFORMS["RC"](x) == x for x in s.values())


def test_permutation_cycles_against_sympy():
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    r, _ = al.permutation_cycles(comp)
    assert r == {"cycles": [("A", "T"), ("C", "G")], "order": 2, "parity": "even"}
    for images in permutations(range(5)):
        mine, _ = al.permutation_cycles(dict(enumerate(images)))
        p = Permutation(list(images))
        assert mine["order"] == p.order()
        assert (mine["parity"] == "even") == p.is_even
        assert sorted(map(sorted, mine["cycles"])) == sorted(map(sorted, (c for c in p.cyclic_form)))
    rev = al.permutation_cycles({i: 13 - i for i in range(1, 13)})[0]
    assert len(rev["cycles"]) == 6 and rev["order"] == 2 and rev["parity"] == "even"


@pytest.mark.parametrize("n", [1, 2, 3])
def test_substitution_monoid(n):
    m, trace = al.substitution_monoid(n)
    assert m["size"] == 5 ** n and m["idempotents"] == 5 ** n
    assert m["units"] == ["·" * n]
    rep = m["report"]
    assert rep["monoid"] and not rep["group"] and not rep["inverses"]["holds"] and not rep["commutative"]["holds"]
    f, x, y, img = m["witness"]
    assert x != y
    assert_trace(trace)
    # hand composition rule: later write wins
    for p in m["elements"]:
        for q in m["elements"]:
            assert m["table"][(p, q)] == "".join(p[i] if p[i] != "·" else q[i] for i in range(n))


def test_concatenation_laws():
    r, _ = al.concatenation_laws(3)
    assert all(r.values())
    assert al.TRANSFORMS["C"]("AC") != al.TRANSFORMS["C"]("C") + al.TRANSFORMS["C"]("A")  # C is not anti


def test_transform_validation():
    assert al.transform("RC", "GATT")[0] == "AATC"
    with pytest.raises(ValueError):
        al.transform("X", "A")
    with pytest.raises(SequenceError):
        al.transform("R", "ACGN")
    assert al.identify(lambda x: x + "A", al.V4, 3) is None


# --- decision D9: the base-level Klein group ---

def test_fixed_point_free_involutions_against_brute_force():
    found, trace = al.fixed_point_free_involutions()
    brute = [dict(zip("ACGT", p)) for p in permutations("ACGT")
             if all(m != b for m, b in zip(p, "ACGT")) and all(dict(zip("ACGT", p))[dict(zip("ACGT", p))[b]] == b
                                                               for b in "ACGT")]
    assert len(found) == len(brute) == 3
    assert sorted(map(sorted, (f.items() for f in found))) == sorted(map(sorted, (b.items() for b in brute)))
    assert trace[0]["data"]["names"] == [n for n in al.BASE_K4 if n != "id"]


def test_base_group_is_klein_not_cyclic():
    table, trace = al.base_cayley_table()
    assert all(table[(a, b)] in al.BASE_K4 for a in al.BASE_K4 for b in al.BASE_K4)
    r, _ = al.check_axioms(al.BASE_K4, table)
    assert r["group"] and r["commutative"]["holds"] and r["identity"]["element"] == "id"
    orders, _ = al.element_orders(al.BASE_K4, table, "id")
    assert sorted(orders.values()) == [1, 2, 2, 2]
    assert_trace(trace)
    domain = list("ACGT")
    idx = {b: i for i, b in enumerate(domain)}
    g = PermutationGroup([Permutation([idx[al.BASE_PERMS[n][b]] for b in domain]) for n in al.BASE_K4])
    assert g.order() == 4 and g.is_abelian and not g.is_cyclic
    assert g.is_transitive() and len(g.orbits()) == 1
    # Burnside on this action, counted directly: the string-level burnside() is a different action
    fixed = {n: sum(al.BASE_PERMS[n][b] == b for b in domain) for n in al.BASE_K4}
    assert fixed == {"id": 4, "(A T)(C G)": 0, "(A G)(C T)": 0, "(A C)(G T)": 0}
    assert sum(fixed.values()) // 4 == 1  # Burnside: one orbit


def test_every_substitution_from_exactly_one_element():
    cov, trace = al.ordered_pair_coverage()
    assert cov == {"per_element": {n: 4 for n in al.BASE_K4 if n != "id"}, "total": 12, "distinct": 12,
                   "covers_all": True}
    assert_trace(trace)
    for x in "ACGT":
        for y in "ACGT":
            if x == y:
                continue
            hits = [n for n in al.BASE_K4 if n != "id" and al.BASE_PERMS[n][x] == y]
            assert len(hits) == 1
            r, _ = al.classify_substitution(x, y)
            assert r["element"] == hits[0]
            # involution => direction-free
            assert al.classify_substitution(y, x)[0] == r
            # independent definition: a transition keeps the purine/pyrimidine class
            assert (r["kind"] == "transition") == ((x in al.PURINES) == (y in al.PURINES)) == r["same_class"]
    assert sum(al.classify_substitution(x, y)[0]["kind"] == "transition"
               for x in "ACGT" for y in "ACGT" if x < y) == 2
    assert sum(al.classify_substitution(x, y)[0]["kind"] == "transversion"
               for x in "ACGT" for y in "ACGT" if x < y) == 4
    with pytest.raises(ValueError):
        al.classify_substitution("A", "A")
    with pytest.raises(SequenceError):
        al.classify_substitution("A", "N")


def test_substitution_census_against_independent_recount(dataset_b):
    seqs = list(dataset_b.seqs().values())
    ref = dataset_b.reference.seq
    for pairs in ([(ref, s) for s in seqs], [(seqs[i], seqs[i + 1]) for i in range(len(seqs) - 1)]):
        r, trace = al.substitution_census(pairs)
        ts = tv = 0
        counts = {}
        for a, b in pairs:  # independent: classify by purine/pyrimidine class, not by group element
            for x, y in zip(a, b):
                if x != y:
                    if (x in "AG") == (y in "AG"):
                        ts += 1
                    else:
                        tv += 1
                    counts["".join(sorted((x, y)))] = counts.get("".join(sorted((x, y))), 0) + 1
        assert (r["transitions"], r["transversions"], r["substitutions"]) == (ts, tv, ts + tv)
        assert r["by_base_pair"] == dict(sorted(counts.items()))
        assert r["by_element"][al.TRANSITION_ELEMENT] == ts
        assert sum(r["by_element"].values()) == r["substitutions"]
        assert r["baseline_ratio"] == 0.5
        assert "limitation" in trace.kinds()
    assert al.substitution_census([("ACGT", "ACGT")])[0] == {
        "substitutions": 0, "by_element": {n: 0 for n in al.BASE_K4 if n != "id"}, "by_base_pair": {},
        "transitions": 0, "transversions": 0, "ratio": None, "baseline_ratio": 0.5}
    with pytest.raises(SequenceError):
        al.substitution_census([("ACGT", "ACG")])
