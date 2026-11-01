"""Experiment 6: shortest paths (Dijkstra) and minimum spanning trees (Prim, Kruskal)."""
from experiments import expected_unit5 as E5
from experiments import expected_unit6 as E6
from experiments._report import ROOT, Check, write
from genomath import mutation_graph as mg
from genomath import tree_analysis as ta
from genomath._core import mismatch_count
from genomath.sequence_processing import load_dataset

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"


def run():
    checks, obs = [], {}
    a, _ = load_dataset(A, reference_id="syn01")
    s = {k[-2:]: v for k, v in a.seqs().items()}
    g2, _ = mg.build_graph(s, 2)
    rows = {}
    for key in E5.DIJKSTRA_ROWS:
        u, v = key.split("→")
        r, _ = mg.path_versus_direct(g2, s, u, v)
        rows[key] = [r["cost"], r["direct"], len(r["paths"])]
        if key == "04→08":
            paths_04_08 = r["paths"]
    checks.append(Check("A, G_2: (cost, direct d_H, #optimal paths)", E5.DIJKSTRA_ROWS, E5.SYN, rows))
    checks.append(Check("A, G_2: optimal paths 04 → 08", E5.PATHS_04_08, E5.SYN, paths_04_08))
    dist = mg.dijkstra(g2, "04")[0]["dist"]
    checks.append(Check("A, G_2: Dijkstra distances from 04", E5.DIST_FROM_04, f"{E5.SOURCE} §5", dist))
    checks.append(Check("A, G_2: targets where cost > direct d_H (detour)", E5.DETOURS_FROM_04, f"{E5.SOURCE} §5",
                        {t: c - mismatch_count(s["04"], s[t]) for t, c in dist.items()
                         if c is not None and c > mismatch_count(s["04"], s[t])}))
    g1, _ = mg.build_graph(s, 1)
    checks.append(Check("A, G_1: cost 03 → 07 (unreachable)", E5.UNREACHABLE_G1_03_07, E5.SYN,
                        mg.path_versus_direct(g1, s, "03", "07")[0]["cost"]))

    b, _ = load_dataset(B)
    bs = b.seqs()
    kstar, _ = mg.smallest_connecting_threshold(bs)
    g, _ = mg.build_graph(bs, kstar)
    ref = b.reference.id
    all_ok, detours = True, {}
    for src in bs:
        d = mg.dijkstra(g, src)[0]["dist"]
        for t in bs:
            direct = mismatch_count(bs[src], bs[t])
            all_ok &= d[t] is not None and d[t] >= direct
            if src == ref and d[t] - direct > 0:
                detours[t] = d[t] - direct
    checks.append(Check(f"B, G_k* (k* = {kstar}): cost ≥ direct d_H for all ordered pairs", True,
                        "triangle inequality (math_audit.md Unit V)", all_ok, consistency=True))
    obs.update(kstar=kstar, ref=ref, detours=detours)
    checks.append(Check(f"B, G_k*: detours from the medoid {ref} (target: extra cost)", None, "", detours, actual_only=True))
    worst = max(detours, key=detours.get) if detours else None
    if worst:
        r, _ = mg.path_versus_direct(g, bs, ref, worst)
        obs["worst"] = {"target": worst, "cost": r["cost"], "direct": r["direct"], "hops": len(r["paths"][0]) - 1,
                        "optimal_paths": len(r["paths"])}
        checks.append(Check("B, G_k*: largest detour from the medoid", None, "", obs["worst"], actual_only=True))
    bg1, _ = mg.build_graph(bs, 1)
    reach = mg.dijkstra(bg1, ref)[0]["dist"]
    obs["reach_g1"] = sum(1 for x in reach.values() if x is not None) - 1
    checks.append(Check("B, G_1: haplotypes reachable from the medoid", None, "", obs["reach_g1"], actual_only=True))

    src6 = lambda sec: f"{E6.SOURCE} {sec}"
    ka = mg.build_graph(s, 10 ** 6)[0]  # complete graph K_9 with weights d_H
    kr, _ = ta.kruskal(ka)
    pr, _ = ta.prim(ka, "01")
    checks.append(Check("A, K_9: Kruskal accepted edges", E6.KRUSKAL_EDGES, src6("§1"), kr["edges"]))
    checks.append(Check("A, K_9: Kruskal examined / accepted / rejected", E6.KRUSKAL_COUNTS, src6("§1"),
                        {"examined": kr["examined"], "accepted": len(kr["edges"]),
                         "rejected": kr["examined"] - len(kr["edges"])}))
    checks.append(Check("A, K_9: Prim order from 01", E6.PRIM_ORDER, src6("§2"), pr["order"]))
    checks.append(Check("A, K_9: Prim edges in order added", E6.PRIM_EDGES, src6("§2"), pr["edges"]))
    checks.append(Check("A, K_9: total weights", E6.TOTALS, src6("§1–2"), {"prim": pr["total"], "kruskal": kr["total"]}))
    checks.append(Check("A, K_9: both are spanning trees", E6.VALID, src6("§3"),
                        {n: ta.validate_spanning_tree(ka.vertices, r["edges"], ka)[0]["is_spanning_tree"]
                         for n, r in (("prim", pr), ("kruskal", kr))}))
    checks.append(Check("A, K_9: number of minimum spanning trees", E6.MST_COUNT, src6("§4"),
                        ta.count_minimum_spanning_trees(ka)[0]))
    checks.append(Check("A: bottleneck = smallest connecting threshold", E6.BOTTLENECK, src6("§3"),
                        {"max MST edge": ta.bottleneck(kr["edges"])[0], "k*": mg.smallest_connecting_threshold(s)[0]}))
    checks.append(Check("A: Kirchhoff spanning-tree count of K_9 (Cayley 9^7)", E6.KIRCHHOFF_K9, src6("§4"),
                        ta.kirchhoff_count(ka.vertices, [(u, v) for u, v, _ in ka.edges()])[0]))

    kb = mg.build_graph(bs, 10 ** 6)[0]
    kbr, _ = ta.kruskal(kb)
    pbr, _ = ta.prim(kb, ref)
    checks.append(Check("B, K_30: Prim total = Kruskal total", True, "MST weight is unique", pbr["total"] == kbr["total"],
                        consistency=True))
    checks.append(Check("B, K_30: both outputs are spanning trees", True, "definition", all(
        ta.validate_spanning_tree(kb.vertices, r["edges"], kb)[0]["is_spanning_tree"] for r in (pbr, kbr)),
        consistency=True))
    checks.append(Check("B: max MST edge = k* (bottleneck property)", True, f"{E6.SOURCE} §3",
                        ta.bottleneck(kbr["edges"])[0] == kstar, consistency=True))
    obs["b_mst"] = {"total": kbr["total"], "count": ta.count_minimum_spanning_trees(kb)[0],
                    "max_edge": ta.bottleneck(kbr["edges"])[0], "same_edges": sorted(
                        tuple(sorted(e[:2])) for e in pbr["edges"]) == sorted(tuple(sorted(e[:2])) for e in kbr["edges"])}
    checks.append(Check("B, K_30: MST weight, number of MSTs, largest edge, Prim/Kruskal same edges", None, "",
                        obs["b_mst"], actual_only=True))
    return checks, obs


def main():
    checks, obs = run()
    worst = obs.get("worst")
    sections = [
        ("Question", "On the sparse threshold graph (edge iff d_H ≤ k, weight d_H), what does Dijkstra's shortest path "
         "cost between observed haplotypes, and when does it exceed the direct Hamming distance? Do Prim and Kruskal "
         "reproduce the hand-derived minimum spanning tree, how many minimum spanning trees exist, and is the largest "
         "MST edge equal to the smallest connecting threshold?"),
        ("Input", f"Dataset A, G_1, G_2 and K_9; Dataset B, G_1, G_k* with k* = {obs['kstar']} (the smallest threshold at "
         "which the graph is connected) and K_30. All edge weights are d_H."),
        ("Mathematical model", "Dijkstra's algorithm, O(V²) array version, non-negative weights d_H ≥ 1; all optimal "
         "paths via the predecessor DAG; triangle inequality cost ≥ d_H. Kruskal (merge sort + union–find with path "
         "compression and union by rank) and Prim (O(V²) array). MST count by weight-level contraction with Kirchhoff's "
         "matrix-tree theorem. Bottleneck (minimax) property: max MST edge = k*."),
        ("Procedure", "Compare the hand-derived Dataset A rows (cost, direct distance, number of optimal paths, full "
         "distance table from 04, unreachable pair in G_1). On Dataset B, run Dijkstra from every haplotype on G_k* and "
         "check cost ≥ d_H for all ordered pairs; report the detours from the medoid as observations. Run Kruskal and "
         "Prim on the complete distance graph of each dataset, compare Dataset A edge by edge with the hand derivation, "
         "validate both outputs as spanning trees, count the MSTs and test the bottleneck property."),
        ("Interpretation",
         "On Dataset A a shortest path through observed haplotypes costs more than the direct distance only for "
         "04 → 07, where the route must pass through syn01 and pay for position 5 twice. Ties are real: 04 → 08 has "
         "three optimal paths. "
         f"On Dataset B (observation) {len(obs['detours'])} of the other 29 haplotypes needs a detour from the medoid "
         f"{obs['ref']} at k* = {obs['kstar']}"
         + (f"; the largest is to {worst['target']} (cost {worst['cost']} vs d_H {worst['direct']}, {worst['hops']} "
            f"hops, {worst['optimal_paths']} optimal path(s))" if worst else "")
         + f". In G_1 only {obs['reach_g1']} haplotypes are reachable from the medoid at all. "
         "Prim and Kruskal both reach weight 15 on Dataset A; 8 different minimum spanning trees exist, so a test can "
         "require equal weight and validity but not equal edge sets. The largest MST edge equals k* on both datasets: "
         f"on Dataset B the MST weighs {obs['b_mst']['total']}, its largest edge is {obs['b_mst']['max_edge']}, and "
         f"{obs['b_mst']['count']} minimum spanning tree(s) exist (observation)."),
        ("Limitations", "Paths run through observed haplotypes only and are not inferred evolutionary trajectories. "
         "The weight d_H counts substitutions without a substitution model. k* is a data-driven choice of threshold; "
         "other thresholds give other paths. Dijkstra here returns one optimal path, plus all optimal paths when asked. "
         "A minimum spanning tree minimises total edge weight among trees on the observed sequences; it is not a "
         "phylogeny and does not establish ancestry (compared with UPGMA in Experiment 7)."),
    ]
    out, summary = write(6, "Shortest paths and spanning trees", "python -m experiments.exp6_paths_and_trees",
                         sections, checks)
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
