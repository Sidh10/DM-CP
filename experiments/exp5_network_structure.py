"""Experiment 5: mutation-network structure as the distance threshold changes (Unit V)."""
from experiments import expected_unit5 as E
from experiments._report import ROOT, Check, write
from genomath import mutation_graph as mg
from genomath.sequence_processing import load_dataset

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"
B_KS = list(range(0, 11)) + [15, 20, 30, 40, 50, 57]


def run():
    src = lambda sec: f"{E.SOURCE} {sec}"
    checks, obs = [], {}
    a, _ = load_dataset(A, reference_id="syn01")
    s = {k[-2:]: v for k, v in a.seqs().items()}
    g1, _ = mg.build_graph(s, 1)
    g2, _ = mg.build_graph(s, 2)
    checks.append(Check("A: edges of G_1", E.G1_EDGES, src("§1"), [(u, v) for u, v, _ in g1.edges()]))
    checks.append(Check("A: adjacency matrix of G_1", E.ADJ_G1, src("§1"), mg.adjacency_matrix(g1)[0]))
    checks.append(Check("A: incidence matrix of G_1", E.INC_G1, src("§1"), mg.incidence_matrix(g1)[0][1]))
    checks.append(Check("A: degrees in G_1 and G_2", [E.DEG_G1, E.DEG_G2], src("§1–2"),
                        [list(mg.degrees(g)[0].values()) for g in (g1, g2)]))
    handshake = all(next(st for st in mg.degrees(g)[1] if "holds" in st["data"])["data"]["holds"] for g in (g1, g2))
    checks.append(Check("A: handshake Σ deg = 2|E| (G_1, G_2)", True, src("§1–2"), handshake, consistency=True))
    b1, _ = mg.bfs(g1, "01")
    checks.append(Check("A: BFS from syn01 in G_1", E.BFS_01, src("§1"), {"order": b1["order"], "hops": b1["hops"]}))
    d1, _ = mg.dfs(g1, "01")
    checks.append(Check("A: DFS from syn01 in G_1", E.DFS_01, src("§1"),
                        {k: d1[k] for k in ("pre", "post", "back_edges")}))
    checks.append(Check("A: cycle found in G_1", E.CYCLE_G1, src("§1"), mg.find_cycle(g1)[0]))
    checks.append(Check("A: cyclomatic number |E| − |V| + c", E.CYCLOMATIC, src("§1–2"),
                        {1: mg.cyclomatic_number(g1)[0], 2: mg.cyclomatic_number(g2)[0]}))
    checks.append(Check("A: bipartite?", E.BIPARTITE, src("§1–2"), {1: mg.is_bipartite(g1)[0], 2: mg.is_bipartite(g2)[0]}))
    checks.append(Check("A: triangles in G_1", E.TRIANGLES_G1, src("§1"), mg.count_triangles(g1)[0]))
    rows, _ = mg.threshold_sweep(s, list(E.SWEEP_A))
    checks.append(Check("A: sweep (edges, components, largest, isolated, cyclomatic)", E.SWEEP_A, src("§3"),
                        {r["k"]: [r["edges"], r["components"], r["largest"], r["isolated"], r["cyclomatic"]] for r in rows}))
    checks.append(Check("A: smallest connecting threshold k*", E.KSTAR_A, src("§3"), mg.smallest_connecting_threshold(s)[0]))
    square, _ = mg.induced_subgraph(g1, ["01", "02", "05", "06"])
    g6, _ = mg.build_graph(s, 6)
    conn = {name: [mg.vertex_connectivity(g)[0], mg.edge_connectivity(g)[0], min(mg.degrees(g)[0].values())]
            for name, g in (("G_1", g1), ("square", square), ("G_6", g6))}
    checks.append(Check("A: (κ, λ, δ)", E.CONNECTIVITY, src("§4"), conn))

    b, _ = load_dataset(B)
    bs = b.seqs()
    brows, _ = mg.threshold_sweep(bs, B_KS)
    checks.append(Check("B: edges of G_k", E.EDGES_B, E.QC,
                        {r["k"]: r["edges"] for r in brows if r["k"] in E.EDGES_B}))
    comps = [r["components"] for r in brows]
    checks.append(Check("B: components non-increasing in k", True, "math_audit.md Unit V (G_j ⊆ G_k)",
                        all(x >= y for x, y in zip(comps, comps[1:])), consistency=True))
    kstar, _ = mg.smallest_connecting_threshold(bs)
    obs["kstar"] = kstar
    obs["sweep"] = brows
    checks.append(Check("B: smallest connecting threshold k*", None, "", kstar, actual_only=True))
    bg1, _ = mg.build_graph(bs, 1)
    c1, _ = mg.connected_components(bg1)
    obs["g1"] = {"components": len(c1), "sizes": sorted((len(c) for c in c1), reverse=True)[:5],
                 "max_degree": max(mg.degrees(bg1)[0].values()), "cyclomatic": mg.cyclomatic_number(bg1)[0],
                 "triangles": mg.count_triangles(bg1)[0], "bipartite": mg.is_bipartite(bg1)[0]}
    checks.append(Check("B: G_=1 structure", None, "", obs["g1"], actual_only=True))
    return checks, obs


def main():
    checks, obs = run()
    table = ["### Dataset B threshold sweep (observation)", "",
             "| k | edges | components | largest | isolated | max degree | cyclomatic |", "|---|---|---|---|---|---|---|"]
    table += [f"| {r['k']} | {r['edges']} | {r['components']} | {r['largest']} | {r['isolated']} | {r['max_degree']} | "
              f"{r['cyclomatic']} |" for r in obs["sweep"]]
    g1 = obs["g1"]
    sections = [
        ("Question", "How do edges, connected components, cycles and connectivity of the mutation-relationship graph "
         "change as the distance threshold k grows, and do the hand-computed graph quantities for Dataset A match?"),
        ("Input", "Dataset A haplotypes (9 vertices); Dataset B haplotypes (30 vertices, n = 1701), thresholds "
         f"{B_KS}."),
        ("Mathematical model", "G_=1: edge iff d_H = 1. G_k: edge iff 1 ≤ d_H ≤ k, weight d_H. Degrees and the "
         "handshake lemma; adjacency and incidence matrices; BFS / DFS; components; cyclomatic number |E| − |V| + c; "
         "bipartiteness; Whitney κ ≤ λ ≤ δ; smallest connecting threshold k*."),
        ("Procedure", "Build G_1, G_2 and the full sweep for Dataset A and compare every quantity with the hand "
         "derivations. Compute κ and λ by exhaustive search on three small graphs. For Dataset B, sweep the listed "
         "thresholds, compare edge counts with the genomath-free QC script, and report the remaining statistics as "
         "observations."),
        ("Interpretation",
         "Dataset A behaves as designed: one 4-cycle at k = 1 (cyclomatic number 1), odd cycles appearing at k = 2 "
         "(the triangle 01–02–03), and the outlier syn09 isolated until k = 6. The non-transitivity of R_k from Unit I "
         "is exactly why components merge in stages rather than all at once. "
         f"Dataset B (observation, no biological claim): G_=1 has {g1['components']} components (sizes, largest first: "
         f"{g1['sizes']}, …), maximum degree {g1['max_degree']}, cyclomatic number {g1['cyclomatic']} and {g1['triangles']} "
         f"triangles, and it is {'bipartite' if g1['bipartite'] else 'not bipartite'}. The graph first becomes "
         f"connected at k* = {obs['kstar']}. Unit VI will check that this equals the largest edge of every minimum "
         "spanning tree."),
        ("Limitations", "Vertices are observed haplotypes only, so a missing edge can mean an unsampled intermediate. "
         "Any single k is arbitrary, which is why results are reported as a sweep. κ and λ are computed by exhaustive "
         "search and are restricted to graphs with at most 16 vertices. Graph structure describes this sample of 30 "
         "haplotypes, not the virus population."),
    ]
    out, summary = write(5, "Mutation-network structure", "python -m experiments.exp5_network_structure",
                         sections, checks, notes_after=["\n".join(table)])
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
