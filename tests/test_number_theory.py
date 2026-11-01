import math
from itertools import product

import pytest
from sympy import factorint, totient
from sympy.ntheory import n_order
from sympy.ntheory.modular import crt as sympy_crt

from conftest import assert_trace
from genomath import number_theory as nt
from genomath._core import SequenceError, reverse_complement

# docs/synthetic_expected.md §7 and docs/unit4_expected.md, transcribed by hand
ENC_A = {"syn01": 9378573, "syn02": 10427149, "syn03": 10459917, "syn04": 10460109, "syn05": 9378561,
         "syn06": 10427137, "syn07": 9392909, "syn08": 9378057, "syn09": 4786988}
CANON_A = {"syn01": 9378573, "syn02": 9378569, "syn03": 9378057, "syn04": 9181449, "syn05": 9378561,
           "syn06": 10427137, "syn07": 9386509, "syn08": 9378057, "syn09": 4786988}
BUCKETS_A = {16: [13, 13, 13, 13, 1, 1, 13, 9, 12], 7: [1, 5, 6, 2, 3, 0, 1, 3, 3],
             13: [9, 5, 0, 10, 10, 6, 6, 0, 11], 17: [13, 12, 4, 9, 1, 0, 1, 7, 9]}


DIGITS = str.maketrans("ACGT", "0123")


def fast_enc(w):
    # independent oracle: Python's own base-4 parser
    return int(w.translate(DIGITS) or "0", 4)


def words(n):
    return ["".join(w) for w in product("ACGT", repeat=n)]


def test_division_and_divides():
    for a in range(-50, 51):
        for b in range(1, 12):
            (q, r), _ = nt.division(a, b)
            assert a == b * q + r and 0 <= r < b
    with pytest.raises(ValueError):
        nt.division(5, 0)
    assert nt.divides(3, 12)[0] and not nt.divides(5, 12)[0] and nt.divides(0, 0)[0] and not nt.divides(0, 3)[0]


def test_gcd_and_extended_against_math():
    g, trace = nt.gcd(252, 198)
    assert g == 18 and len([s for s in trace if s["kind"] == "algorithm"]) == 4  # four division steps (§5)
    assert nt.extended_gcd(252, 198)[0] == (18, 4, -5)
    for a in range(-30, 31):
        for b in range(-30, 31):
            assert nt.gcd(a, b)[0] == math.gcd(a, b)
            g, x, y = nt.extended_gcd(a, b)[0]
            assert g == math.gcd(a, b) and a * x + b * y == g
    assert nt.gcd(ENC_A["syn01"], 4 ** 12)[0] == 1


def test_inverses_against_pow():
    assert nt.mod_inverse(11, 7)[0] == 2 and nt.mod_inverse(7, 11)[0] == 8 and nt.mod_inverse(4, 7)[0] == 2
    with pytest.raises(ValueError):
        nt.mod_inverse(4, 16)
    for m in range(2, 60):
        for a in range(0, m):
            if math.gcd(a, m) == 1:
                assert nt.mod_inverse(a, m)[0] == pow(a, -1, m)
            else:
                with pytest.raises(ValueError):
                    nt.mod_inverse(a, m)


def test_linear_congruences_against_brute_force():
    assert nt.solve_linear_congruence(3, 6, 9)[0] == [2, 5, 8]
    assert nt.solve_linear_congruence(2, 1, 4)[0] == []
    for m in range(2, 25):
        for a in range(1, m):
            for b in range(m):
                assert nt.solve_linear_congruence(a, b, m)[0] == [x for x in range(m) if (a * x - b) % m == 0]


def test_crt_against_sympy_and_all_3mers():
    assert nt.crt([0, 2], [7, 11])[0] == 35
    for mods in [(7, 11), (5, 13), (3, 5, 7), (4, 9, 25)]:
        for x in range(0, 64):
            res = [x % m for m in mods]
            assert nt.crt(res, list(mods))[0] == int(sympy_crt(list(mods), res)[0]) == x % math.prod(mods)
    assert sum(nt.crt([x % 3, x % 5], [3, 5])[0] == x for x in range(64)) == 15
    with pytest.raises(ValueError):
        nt.crt([1, 1], [4, 6])


def test_phi_factorize_order_pow_against_sympy():
    for n in range(1, 400):
        assert nt.factorize(n) == factorint(n)
        assert nt.euler_phi(n)[0] == totient(n)
    hand = {7: 3, 13: 6, 17: 4, 31: 5, 61: 30}
    for m, k in hand.items():
        assert nt.multiplicative_order(4, m)[0] == k == n_order(4, m)
        assert (m - 1) % k == 0
    for m in range(2, 80):
        for a in range(0, 20):
            assert nt.mod_pow(a, 37, m)[0] == pow(a, 37, m)
        if math.gcd(4, m) == 1:
            assert nt.mod_pow(4, nt.euler_phi(m)[0], m)[0] == 1 % m
    assert nt.mod_pow(4, 8, 16)[0] == 0
    with pytest.raises(ValueError):
        nt.multiplicative_order(4, 16)


def test_encode_decode_bijection_exhaustive():
    for n in range(0, 7):
        seen = set()
        for w in words(n):
            x, _ = nt.encode(w)
            assert x == fast_enc(w) and nt.decode(x, n)[0] == w
            seen.add(x)
        assert seen == set(range(4 ** n))
    assert nt.encode("ACG")[0] == 6 and nt.encode("GAT")[0] == 35
    with pytest.raises(SequenceError):
        nt.encode("ACGN")
    with pytest.raises(ValueError):
        nt.decode(64, 3)


def test_complement_and_rc_codes():
    for n in range(0, 6):
        for w in words(n):
            c, trace = nt.complement_code(w)
            assert c == 4 ** n - 1 - nt.encode(w)[0]
            assert nt.rc_code(nt.encode(w)[0], n)[0] == nt.encode(reverse_complement(w))[0]
    assert next(s for s in trace if "holds" in s["data"])["data"]["holds"]


def test_canonical_is_lexicographic_min():
    for n in range(1, 6):
        for w in words(n):
            rep = min(w, reverse_complement(w))
            assert nt.canonical_code(w)[0] == nt.encode(rep)[0]


def test_dataset_a_hand_values(dataset_a):
    s = dataset_a.seqs()
    assert {h: nt.encode(x)[0] for h, x in s.items()} == ENC_A
    assert {h: nt.canonical_code(x)[0] for h, x in s.items()} == CANON_A
    codes = {h: nt.encode(x)[0] for h, x in s.items()}
    for m, expected in BUCKETS_A.items():
        r, trace = nt.bucket_report(codes, m)
        got = {h: b for b, hs in r["buckets"].items() for h in hs}
        assert [got[h] for h in s] == expected
        assert r["max_load"] >= r["pigeonhole_bound"]
        assert_trace(trace)


def test_kmer_codes_rolling_equals_direct():
    syn01 = "GATTACGTAATC"
    assert nt.kmer_codes(syn01, 3)[0] == [35, 15, 60, 49, 6, 27, 44, 48, 3, 13]
    for k in range(1, 13):
        assert nt.kmer_codes(syn01, k)[0] == [nt.encode(syn01[i:i + k])[0] for i in range(13 - k)]
    with pytest.raises(ValueError):
        nt.kmer_codes(syn01, 13)


def test_kmer_census_syn01_and_identity(dataset_b):
    r3, _ = nt.kmer_census(["GATTACGTAATC"], 3)
    r4, _ = nt.kmer_census(["GATTACGTAATC"], 4)
    assert (r3["distinct"], r3["canonical"], r3["palindromic"]) == (10, 5, 0)
    assert (r4["distinct"], r4["canonical"], r4["palindromic"]) == (9, 5, 1)
    seqs = list(dataset_b.seqs().values())
    for k in (4, 6, 8):
        r, _ = nt.kmer_census(seqs, k)
        brute = {s[i:i + k] for s in seqs for i in range(len(s) - k + 1)}
        assert r["distinct"] == len(brute)
        assert r["canonical"] == len({min(w, reverse_complement(w)) for w in brute})
        assert r["palindromic"] == sum(reverse_complement(w) == w for w in brute)
        assert r["identity_holds"] and r["canonical"] <= r["space_canonical"]


def test_last_bases_seen():
    assert nt.last_bases_seen(16)[0] == 2 and nt.last_bases_seen(64)[0] == 3 and nt.last_bases_seen(4)[0] == 1
    assert nt.last_bases_seen(7)[0] is None and nt.last_bases_seen(12)[0] is None
    for m in (16, 64):
        j = nt.last_bases_seen(m)[0]
        for w in words(5):  # bucket is determined by the last j bases
            assert nt.encode(w)[0] % m == nt.encode(w[-j:])[0] % m


def test_collision_criterion_against_brute_force():
    # brute force: does some single substitution in some string of length n keep x mod m?
    for m in range(1, 70):
        t_min, _ = nt.substitution_collision_t_min(m)
        for n in range(1, 6):
            brute = any(fast_enc(w) % m == fast_enc(w[:i] + b + w[i + 1:]) % m
                        for w in words(n) for i in range(n) for b in "ACGT" if b != w[i])
            assert brute == (t_min is not None and n - 1 >= t_min), (m, n)


def test_single_substitution_pairs(dataset_a):
    s = dataset_a.seqs()
    pairs = [("syn01", "syn02"), ("syn01", "syn05"), ("syn02", "syn03"), ("syn02", "syn06"), ("syn03", "syn04"),
             ("syn05", "syn06")]
    assert nt.substitution_collisions(pairs, s, 16)[0] == [("syn01", "syn02"), ("syn02", "syn03"), ("syn03", "syn04"),
                                                          ("syn05", "syn06")]
    for m in (7, 13, 17, 31, 61):
        assert nt.substitution_collisions(pairs, s, m)[0] == []


def test_fermat_cancellation_example():
    syn01, other = "GATTACGTAATC", "GATTAAGTAATG"
    x, y = nt.encode(syn01)[0], nt.encode(other)[0]
    assert x - y == 4095 or y - x == 4095
    collide = [m for m in (16, 64, 7, 13, 17, 31, 61) if x % m == y % m]
    assert collide == [7, 13]
