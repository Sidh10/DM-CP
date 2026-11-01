import random
from fractions import Fraction
from itertools import combinations

import networkx as nx
import pytest
from networkx.algorithms.tree.mst import SpanningTreeIterator
from sympy import Matrix

from conftest import assert_trace
from genomath import mutation_graph as mg
from genomath import phylogenetic_analysis as pa
from genomath import tree_analysis as ta
from genomath._core import mismatch_count

# docs/unit6_expected.md, transcribed by hand
KRUSKAL_A = [("01", "02", 1), ("01", "05", 1), ("02", "03", 1), ("02", "06", 1), ("03", "04", 1), ("01", "07", 2),
             ("01", "08", 2), ("01", "09", 6)]
PRIM_A = [("01", "02", 1), ("02", "03", 1), ("03", "04", 1), ("01", "05", 1), ("02", "06", 1), ("01", "07", 2),
          ("01", "08", 2), ("01", "09", 6)]
UPGMA_A = [("U1", "01", "02", Fraction(1), 6), ("U2", "03", "04", Fraction(1), 2), ("U3", "05", "06", Fraction(1), 1),
           ("U4", "U1", "U3", Fraction(3, 2), 1), ("U5", "U4", "U2", Fraction(5, 2), 2),
           ("U6", "U5", "07", Fraction(19, 6), 2), ("U7", "U6", "08", Fraction(23, 7), 1),
           ("U8", "U7", "09", Fraction(61, 8), 1)]
NEWICK_A = "((((((01,02),(05,06)),(03,04)),07),08),09);"


@pytest.fixture(scope="module")
def a(dataset_a):
    s = {k[-2:]: v for k, v in dataset_a.seqs().items()}
    g, _ = mg.build_graph(s, 12)
    (ids, d), _ = mg.distance_matrix(s)
    return s, g, ids, d


def to_nx(g):
    h = nx.Graph()
    h.add_nodes_from(g.vertices)
    h.add_weighted_edges_from(g.edges())
    return h


def random_sequences(rng, count, length):
    seqs = {}
    while len(seqs) < count:
        w = "".join(rng.choice("ACGT") for _ in range(length))
        if w not in seqs.values():
            seqs[f"s{len(seqs):02d}"] = w
    return seqs


def test_kruskal_and_prim_hand_values(a):
    s, g, _, _ = a
    k, tk = ta.kruskal(g)
    p, tp = ta.prim(g, "01")
    assert k["edges"] == KRUSKAL_A and k["total"] == 15 and k["examined"] == 29
    assert p["edges"] == PRIM_A and p["total"] == 15 and p["order"] == ["01", "02", "03", "04", "05", "06", "07", "08", "09"]
    assert_trace(tk)
    assert_trace(tp)
    for edges in (k["edges"], p["edges"]):
        assert ta.validate_spanning_tree(g.vertices, edges, g)[0]["is_spanning_tree"]
    assert ta.count_minimum_spanning_trees(g)[0] == 8
    assert ta.bottleneck(k["edges"])[0] == 6 == mg.smallest_connecting_threshold(s)[0]


def test_rooted_tree_hand_values(a):
    _, g, _, _ = a
    t1, _ = ta.root_tree(g.vertices, KRUSKAL_A, "01")
    st, _ = ta.tree_statistics(t1)
    assert st["height"] == 3 and set(st["leaves"]) == {"04", "05", "06", "07", "08", "09"}
    assert set(st["internal"]) == {"01", "02", "03"}
    assert ta.preorder(t1)[0] == ["01", "02", "03", "04", "06", "05", "07", "08", "09"]
    assert ta.postorder(t1)[0] == ["04", "03", "06", "02", "05", "07", "08", "09", "01"]
    assert ta.level_order(t1)[0] == ["01", "02", "05", "07", "08", "09", "03", "06", "04"]
    t4, _ = ta.root_tree(g.vertices, KRUSKAL_A, "04")
    st4, _ = ta.tree_statistics(t4)
    assert st4["height"] == 4 and set(st4["leaves"]) == {"05", "06", "07", "08", "09"}
    assert ta.preorder(t4)[0] == ["04", "03", "02", "01", "05", "07", "08", "09", "06"]
    assert ta.postorder(t4)[0] == ["05", "07", "08", "09", "01", "06", "02", "03", "04"]
    with pytest.raises(ValueError):
        ta.inorder(t1)  # the MST is not binary
    with pytest.raises(ValueError):
        ta.root_tree(g.vertices, KRUSKAL_A[:-1], "01")


def test_upgma_hand_values(a):
    _, _, ids, d = a
    u, trace = pa.upgma(ids, d)
    assert [(m["node"], m["left"], m["right"], m["distance"], m["tied_pairs"]) for m in u["merges"]] == UPGMA_A
    assert u["newick"] == NEWICK_A and u["tie_steps"] == 4
    assert "limitation" in trace.kinds()
    assert ta.inorder(u["tree"])[0] == ["01", "U1", "02", "U4", "05", "U3", "06", "U5", "03", "U2", "04", "U6", "07",
                                        "U7", "08", "U8", "09"]
    ok, t = pa.is_ultrametric(ids, d)
    assert not ok and t[0]["data"]["witness"] == ("01", "02", "03")
    fit, _ = pa.fit_to_input(u, ids, d)
    assert (fit["exact"], fit["max_deviation"], fit["worst_pair"]) == (3, Fraction(12, 7), ("04", "08"))
    assert pa.upgma(ids, d, "reverse")[0]["newick"] != NEWICK_A  # the tie rule matters here
    with pytest.raises(ValueError):
        pa.upgma(ids, d, "random")


def test_splits_comparison(a):
    _, g, ids, d = a
    ms, _ = ta.edge_splits(g.vertices, KRUSKAL_A)
    us, _ = pa.upgma_splits(pa.upgma(ids, d)[0], ids)
    r, _ = pa.compare_splits(ms, us)
    assert (r["mst"], r["upgma"], r["shared"], r["symmetric_difference"]) == (2, 6, 1, 6)
    assert r["shared_sets"] == [["03", "04"]]


def test_counting_link():
    assert ta.cayley_count(9)[0] == 4_782_969 and ta.rooted_binary_topologies(9)[0] == 2_027_025
    t = 1
    for n in range(2, 12):
        t *= 2 * n - 3
        assert ta.rooted_binary_topologies(n)[0] == t
        kn = list(range(n))
        assert ta.kirchhoff_count(kn, list(combinations(kn, 2)))[0] == n ** (n - 2)


def test_mst_against_networkx(a, dataset_b):
    for seqs in (a[0], dataset_b.seqs()):
        g, _ = mg.build_graph(seqs, 10 ** 6)
        h = to_nx(g)
        best = nx.minimum_spanning_tree(h).size(weight="weight")
        for r in (ta.kruskal(g)[0], ta.prim(g)[0]):
            assert r["total"] == best
            assert nx.is_tree(nx.Graph([(u, v) for u, v, _ in r["edges"]]))
        assert ta.bottleneck(ta.kruskal(g)[0]["edges"])[0] == mg.smallest_connecting_threshold(seqs)[0]


def test_mst_count_against_spanning_tree_iterator():
    rng = random.Random(3)
    for _ in range(25):
        n = rng.randint(3, 7)
        edges = [(u, v, rng.randint(1, 3)) for u in range(n) for v in range(u + 1, n) if rng.random() < 0.7]
        g = mg.Graph(range(n), edges)
        h = to_nx(g)
        if not nx.is_connected(h):
            continue
        best, count = None, 0
        for t in SpanningTreeIterator(h, weight="weight"):
            w = t.size(weight="weight")
            if best is None:
                best = w
            if w > best:
                break
            count += 1
        assert ta.count_minimum_spanning_trees(g)[0] == count
        assert ta.kruskal(g)[0]["total"] == ta.prim(g)[0]["total"] == best


def test_kirchhoff_against_sympy_determinant():
    rng = random.Random(5)
    for _ in range(30):
        n = rng.randint(2, 7)
        multi = [(rng.randrange(n), rng.randrange(n)) for _ in range(rng.randint(n, 3 * n))]
        multi = [(u, v) for u, v in multi if u != v]
        lap = [[0] * n for _ in range(n)]
        for u, v in multi:
            lap[u][u] += 1
            lap[v][v] += 1
            lap[u][v] -= 1
            lap[v][u] -= 1
        expected = Matrix([row[1:] for row in lap[1:]]).det()
        assert ta.kirchhoff_count(range(n), multi)[0] == expected


def test_merge_sort_and_union_find():
    rng = random.Random(9)
    for _ in range(50):
        xs = [(rng.randint(0, 5), i) for i in range(rng.randint(0, 40))]
        assert ta.merge_sort(xs, key=lambda x: x[0]) == sorted(xs, key=lambda x: x[0])  # also checks stability
    uf = ta.UnionFind(range(10))
    assert uf.union(1, 2) and uf.union(2, 3) and not uf.union(1, 3)
    assert uf.find(3) == uf.find(1) != uf.find(4)


def naive_upgma(ids, d):
    # independent: recompute every cluster distance as the mean over member pairs, same D8 tie rule
    clusters = {x: (x,) for x in ids}
    pos = {x: i for i, x in enumerate(ids)}
    merges, step = [], 0
    while len(clusters) > 1:
        step += 1
        names = list(clusters)
        pairs = []
        for i, p in enumerate(names):
            for q in names[i + 1:]:
                mean = Fraction(sum(d[pos[x]][pos[y]] for x in clusters[p] for y in clusters[q]),
                                len(clusters[p]) * len(clusters[q]))
                pairs.append((mean, tuple(sorted((clusters[p], clusters[q]))), p, q))
        best = min(pairs, key=lambda t: (t[0], t[1]))
        _, _, p, q = best
        merged = tuple(sorted(clusters[p] + clusters[q]))
        merges.append((merged, best[0]))
        del clusters[p], clusters[q]
        clusters[f"U{step}"] = merged
    return merges


def test_upgma_against_naive_and_theorems():
    rng = random.Random(21)
    for _ in range(25):
        seqs = random_sequences(rng, rng.randint(3, 9), 8)
        (ids, d), _ = mg.distance_matrix(seqs)
        u, _ = pa.upgma(ids, d)
        assert [(u["members"][m["node"]], m["distance"]) for m in u["merges"]] == naive_upgma(ids, d)
        heights = [m["height"] for m in u["merges"]]
        assert heights == sorted(heights)
        coph, _ = pa.cophenetic(u, ids)
        assert pa.is_ultrametric(ids, coph)[0]
        # UPGMA reproduces an ultrametric input exactly
        u2, _ = pa.upgma(ids, coph)
        assert pa.cophenetic(u2, ids)[0] == coph
        k = ta.kruskal(mg.build_graph(seqs, 100)[0])[0]
        assert ta.bottleneck(k["edges"])[0] == mg.smallest_connecting_threshold(seqs)[0]


def test_is_ultrametric_against_definition():
    rng = random.Random(4)
    for _ in range(40):
        n = rng.randint(3, 6)
        d = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d[i][j] = d[j][i] = rng.randint(1, 3)
        brute = all(sorted((d[i][j], d[i][k], d[j][k]))[1] == sorted((d[i][j], d[i][k], d[j][k]))[2]
                    for i, j, k in combinations(range(n), 3))
        assert pa.is_ultrametric(list(range(n)), d)[0] == brute


def test_edge_splits_against_networkx(a):
    _, g, _, _ = a
    ms, _ = ta.edge_splits(g.vertices, KRUSKAL_A)
    t = nx.Graph([(u, v) for u, v, _ in KRUSKAL_A])
    expected = set()
    for u, v, _ in KRUSKAL_A:
        t2 = t.copy()
        t2.remove_edge(u, v)
        side = frozenset(next(c for c in nx.connected_components(t2) if "01" not in c))
        if 2 <= len(side) <= 7:
            expected.add(side)
    assert ms == expected


def test_upgma_branch_lengths_hand_values_and_sums(a):
    _, _, ids, d = a
    u, _ = pa.upgma(ids, d)
    lengths, _ = pa.branch_lengths(u)
    F = Fraction
    hand = {"01": F(1, 2), "02": F(1, 2), "03": F(1, 2), "04": F(1, 2), "05": F(1, 2), "06": F(1, 2),
            "U1": F(1, 4), "U3": F(1, 4), "U4": F(1, 2), "U2": F(3, 4), "U5": F(1, 3), "07": F(19, 12),
            "U6": F(5, 84), "08": F(23, 14), "U7": F(243, 112), "09": F(61, 16)}
    assert lengths == hand
    rng = random.Random(8)
    for _ in range(15):
        seqs = random_sequences(rng, rng.randint(2, 9), 8)
        (rid, rd), _ = mg.distance_matrix(seqs)
        r, _ = pa.upgma(rid, rd)
        ln, _ = pa.branch_lengths(r)
        for leaf in rid:
            total, v = Fraction(0), leaf
            while r["tree"].parent[v] is not None:
                total += ln[v]
                v = r["tree"].parent[v]
            assert total == r["height"][r["root"]] and all(x >= 0 for x in ln.values())


def test_determinant_on_general_matrices():
    # Kirchhoff only feeds positive semidefinite Laplacian minors, where no row swap ever happens; test the swap path too
    rng = random.Random(13)
    for _ in range(40):
        n = rng.randint(1, 5)
        m = [[rng.randint(-3, 3) if rng.random() < 0.6 else 0 for _ in range(n)] for _ in range(n)]
        assert ta._determinant(m) == Matrix(m).det()
    assert ta._determinant([[0, 1], [1, 0]]) == -1
