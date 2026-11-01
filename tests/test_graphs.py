import random

import networkx as nx
import pytest

from conftest import assert_trace
from genomath import mutation_graph as mg
from genomath._core import mismatch_count

# docs/unit5_expected.md, transcribed by hand (vertex 'NN' = synNN)
ADJ_G1 = [
    [0, 1, 0, 0, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 0],
    [0] * 9, [0] * 9, [0] * 9,
]
INC_G1 = [
    [1, 1, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0], [0, 0, 1, 0, 1, 0], [0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1],
    [0, 0, 0, 1, 0, 1], [0] * 6, [0] * 6, [0] * 6,
]
SWEEP_A = {0: (0, 9, 1, 9, 0), 1: (6, 4, 6, 3, 1), 2: (14, 2, 8, 1, 7), 3: (22, 2, 8, 1, 15), 4: (27, 2, 8, 1, 20),
           5: (28, 2, 8, 1, 21), 6: (29, 1, 9, 0, 21), 9: (36, 1, 9, 0, 28)}
DIST_FROM_04 = {"01": 3, "02": 2, "03": 1, "04": 0, "05": 4, "06": 3, "07": 5, "08": 5, "09": None}
EDGES_B = {1: 8, 2: 31, 3: 57, 4: 90, 5: 113, 6: 148, 7: 178, 8: 202, 9: 223, 10: 227, 57: 435}  # QC script


@pytest.fixture(scope="module")
def a_seqs(dataset_a):
    return {k[-2:]: v for k, v in dataset_a.seqs().items()}


def to_nx(g):
    h = nx.Graph()
    h.add_nodes_from(g.vertices)
    h.add_weighted_edges_from(g.edges())
    return h


def random_graphs(count=60, n=7, seed=7):
    rng = random.Random(seed)
    for _ in range(count):
        p = rng.choice([0.15, 0.3, 0.5, 0.8])
        edges = [(u, v, rng.randint(1, 5)) for u in range(n) for v in range(u + 1, n) if rng.random() < p]
        yield mg.Graph(range(n), edges)


def test_g1_matrices_degrees_and_traversals(a_seqs):
    g, trace = mg.build_graph(a_seqs, 1)
    assert_trace(trace)
    assert [(u, v) for u, v, _ in g.edges()] == [("01", "02"), ("01", "05"), ("02", "03"), ("02", "06"), ("03", "04"),
                                                ("05", "06")]
    assert mg.adjacency_matrix(g)[0] == ADJ_G1
    assert mg.incidence_matrix(g)[0][1] == INC_G1
    deg, t = mg.degrees(g)
    assert list(deg.values()) == [2, 3, 2, 1, 2, 2, 0, 0, 0]
    assert next(s for s in t if "holds" in s["data"])["data"]["holds"]
    b, _ = mg.bfs(g, "01")
    assert b["order"] == ["01", "02", "05", "03", "06", "04"]
    assert b["hops"] == {"01": 0, "02": 1, "05": 1, "03": 2, "06": 2, "04": 3}
    assert b["parent"] == {"01": None, "02": "01", "05": "01", "03": "02", "06": "02", "04": "03"}
    d, _ = mg.dfs(g, "01")
    assert d["pre"] == ["01", "02", "03", "04", "06", "05"] and d["post"] == ["04", "03", "05", "06", "02", "01"]
    assert d["back_edges"] == [("05", "01")]
    assert mg.find_cycle(g)[0] == ["01", "02", "06", "05"]


def test_g1_g2_structure(a_seqs):
    g1, _ = mg.build_graph(a_seqs, 1)
    g2, _ = mg.build_graph(a_seqs, 2)
    single, _ = mg.build_graph(a_seqs, mode="single")
    assert single.edges() == g1.edges()
    assert mg.cyclomatic_number(g1)[0] == 1 and mg.cyclomatic_number(g2)[0] == 7
    ok, t = mg.is_bipartite(g1)
    assert ok and t[0]["data"]["sides"] == [["01", "03", "06", "07", "08", "09"], ["02", "04", "05"]]
    assert not mg.is_bipartite(g2)[0]
    assert mg.count_triangles(g1)[0] == 0
    assert list(mg.degrees(g2)[0].values()) == [6, 5, 4, 2, 4, 4, 1, 2, 0]
    with pytest.raises(ValueError):
        mg.build_graph(a_seqs, mode="other")


def test_threshold_sweep_and_kstar(a_seqs):
    rows, _ = mg.threshold_sweep(a_seqs, list(SWEEP_A))
    assert {r["k"]: (r["edges"], r["components"], r["largest"], r["isolated"], r["cyclomatic"]) for r in rows} == SWEEP_A
    assert mg.smallest_connecting_threshold(a_seqs)[0] == 6


def test_connectivity_hand_examples(a_seqs):
    g1, _ = mg.build_graph(a_seqs, 1)
    sq, _ = mg.induced_subgraph(g1, ["01", "02", "05", "06"])
    g6, _ = mg.build_graph(a_seqs, 6)
    assert [(mg.vertex_connectivity(x)[0], mg.edge_connectivity(x)[0]) for x in (g1, sq, g6)] == [(0, 0), (2, 2), (1, 1)]


def test_dijkstra_hand_rows(a_seqs):
    g2, _ = mg.build_graph(a_seqs, 2)
    r, trace = mg.dijkstra(g2, "04")
    assert r["dist"] == DIST_FROM_04
    assert_trace(trace)
    rows = {(s, t): mg.path_versus_direct(g2, a_seqs, s, t)[0] for s, t in [("04", "08"), ("07", "03"), ("07", "08")]}
    assert [(x["cost"], x["direct"], len(x["paths"])) for x in rows.values()] == [(5, 5, 3), (4, 3, 2), (4, 4, 1)]
    detours = [v for v, c in DIST_FROM_04.items() if c is not None and c > mismatch_count(a_seqs["04"], a_seqs[v])]
    assert detours == ["07"]
    g1, _ = mg.build_graph(a_seqs, 1)
    assert mg.path_versus_direct(g1, a_seqs, "03", "07")[0]["cost"] is None


def test_against_networkx_on_both_datasets(a_seqs, dataset_b):
    for seqs, ks in ((a_seqs, range(0, 10)), (dataset_b.seqs(), [1, 2, 5, 9, 20, 57])):
        for k in ks:
            g, _ = mg.build_graph(seqs, k)
            h = to_nx(g)
            assert sorted(map(sorted, mg.connected_components(g)[0])) == sorted(map(sorted, nx.connected_components(h)))
            assert mg.degrees(g)[0] == dict(h.degree())
            assert mg.is_bipartite(g)[0] == nx.is_bipartite(h)
            assert mg.count_triangles(g)[0] == sum(nx.triangles(h).values()) // 3
            assert mg.cyclomatic_number(g)[0] == len(nx.cycle_basis(h))
            src = g.vertices[0]
            assert {v: x for v, x in mg.bfs(g, src)[0]["hops"].items()} == nx.single_source_shortest_path_length(h, src)
            dist = mg.dijkstra(g, src)[0]["dist"]
            assert {v: x for v, x in dist.items() if x is not None} == nx.single_source_dijkstra_path_length(h, src)


def test_dataset_b_edges_match_qc(dataset_b):
    rows, _ = mg.threshold_sweep(dataset_b.seqs(), list(EDGES_B))
    assert {r["k"]: r["edges"] for r in rows} == EDGES_B


def test_random_graphs_against_networkx():
    for g in random_graphs():
        h = to_nx(g)
        assert mg.vertex_connectivity(g)[0] == nx.node_connectivity(h)
        assert mg.edge_connectivity(g)[0] == nx.edge_connectivity(h)
        assert mg.is_bipartite(g)[0] == nx.is_bipartite(h)
        cyc, _ = mg.find_cycle(g)
        assert (cyc is None) == nx.is_forest(h)
        if cyc:
            assert all(h.has_edge(cyc[i], cyc[(i + 1) % len(cyc)]) for i in range(len(cyc))) and len(set(cyc)) == len(cyc)
        back = sum(len(mg.dfs(g, c[0])[0]["back_edges"]) for c in mg.connected_components(g)[0])
        assert back == mg.cyclomatic_number(g)[0] == len(nx.cycle_basis(h))
        for t in g.vertices:
            mine, _ = mg.all_shortest_paths(g, 0, t)
            theirs = sorted(nx.all_shortest_paths(h, 0, t, weight="weight")) if nx.has_path(h, 0, t) else []
            assert sorted(mine) == theirs
        kappa, lam = mg.vertex_connectivity(g)[0], mg.edge_connectivity(g)[0]
        assert kappa <= lam <= min(dict(h.degree()).values())


def test_graph_operations(a_seqs):
    (ids, d), _ = mg.distance_matrix(a_seqs)
    w = lambda u, v: d[ids.index(u)][ids.index(v)]
    g1, _ = mg.build_graph(a_seqs, 1)
    g3, _ = mg.build_graph(a_seqs, 3)
    assert mg.edge_set_operation(g1, g3, "union")[0].edges() == g3.edges()
    assert mg.edge_set_operation(g1, g3, "intersection")[0].edges() == g1.edges()
    diff, _ = mg.edge_set_operation(g3, g1, "difference")
    assert all(2 <= x <= 3 for _, _, x in diff.edges()) and len(diff.edges()) == 22 - 6
    comp, _ = mg.complement(g3, w)
    assert len(comp.edges()) == 36 - 22 and all(x > 3 for _, _, x in comp.edges())


def test_path_cost_never_below_direct(dataset_b):
    seqs = dataset_b.seqs()
    k, _ = mg.smallest_connecting_threshold(seqs)
    g, _ = mg.build_graph(seqs, k)
    for s in seqs:
        dist = mg.dijkstra(g, s)[0]["dist"]
        assert all(dist[t] >= mismatch_count(seqs[s], seqs[t]) for t in seqs)


def test_guards():
    big = mg.Graph(range(20), [(i, i + 1, 1) for i in range(19)])
    with pytest.raises(ValueError):
        mg.vertex_connectivity(big)
    with pytest.raises(ValueError):
        mg.Graph([1], [(1, 1, 1)])
    assert mg.vertex_connectivity(mg.Graph([0, 1, 2], [(0, 1, 1), (1, 2, 1), (0, 2, 1)]))[0] == 2  # K3: n − 1


def test_connectivity_below_min_degree():
    # barbell: two K4s joined by one bridge, so λ = κ = 1 < δ = 3; relabelling puts the bridge in every orientation
    rng = random.Random(11)
    base = nx.barbell_graph(4, 0)
    for _ in range(20):
        labels = list(range(8))
        rng.shuffle(labels)
        h = nx.relabel_nodes(base, dict(enumerate(labels)))
        g = mg.Graph(range(8), [(u, v, 1) for u, v in h.edges()])
        assert mg.edge_connectivity(g)[0] == nx.edge_connectivity(h) == 1
        assert mg.vertex_connectivity(g)[0] == nx.node_connectivity(h) == 1
        assert min(dict(h.degree()).values()) == 3
