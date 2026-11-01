"""Experiment 7: minimum spanning tree versus UPGMA tree (Unit VI)."""
from experiments import expected_unit6 as E
from experiments._report import ROOT, Check, write
from genomath import mutation_graph as mg
from genomath import phylogenetic_analysis as pa
from genomath import tree_analysis as ta
from genomath.sequence_processing import load_dataset

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"


def rooted_summary(vertices, edges, root, with_level=True):
    t, _ = ta.root_tree(vertices, edges, root)
    st, _ = ta.tree_statistics(t)
    out = {"height": st["height"], "leaves": sorted(st["leaves"]), "internal": sorted(st["internal"]),
           "preorder": ta.preorder(t)[0], "postorder": ta.postorder(t)[0]}
    if with_level:
        out["level"] = ta.level_order(t)[0]
    return out


def run():
    src = lambda sec: f"{E.SOURCE} {sec}"
    checks, obs = [], {}
    a, _ = load_dataset(A, reference_id="syn01")
    s = {k[-2:]: v for k, v in a.seqs().items()}
    k9, _ = mg.build_graph(s, 10 ** 6)
    mst = ta.kruskal(k9)[0]["edges"]
    checks.append(Check("A: MST rooted at 01", E.ROOT_01, src("§5"), rooted_summary(k9.vertices, mst, "01")))
    checks.append(Check("A: MST rooted at 04", E.ROOT_04, src("§5"), rooted_summary(k9.vertices, mst, "04", False)))
    (ids, d), _ = mg.distance_matrix(s)
    u, _ = pa.upgma(ids, d)
    checks.append(Check("A: UPGMA merges (node, left, right, distance, height, tied pairs)", E.UPGMA_MERGES, src("§6"),
                        [[m["node"], m["left"], m["right"], m["distance"], m["height"], m["tied_pairs"]]
                         for m in u["merges"]]))
    checks.append(Check("A: UPGMA Newick", E.NEWICK, src("§6"), u["newick"]))
    checks.append(Check("A: merge steps with ties", E.TIE_STEPS, src("§6"), u["tie_steps"]))
    checks.append(Check("A: reverse tie rule changes the topology", E.REVERSE_TIE_CHANGES_TOPOLOGY, src("§6"),
                        pa.upgma(ids, d, "reverse")[0]["newick"] != u["newick"]))
    checks.append(Check("A: UPGMA in-order traversal", E.INORDER, src("§6"), ta.inorder(u["tree"])[0]))
    ok, t = pa.is_ultrametric(ids, d)
    checks.append(Check("A: input distances ultrametric? (witness)", E.ULTRAMETRIC, src("§6"),
                        (ok, t[0]["data"].get("witness"))))
    checks.append(Check("A: cophenetic fit to d_H", E.FIT, src("§6"),
                        {k: v for k, v in pa.fit_to_input(u, ids, d)[0].items() if k in E.FIT}))
    ms, _ = ta.edge_splits(k9.vertices, mst)
    us, _ = pa.upgma_splits(u, ids)
    cmp, _ = pa.compare_splits(ms, us)
    checks.append(Check("A: MST vs UPGMA non-trivial bipartitions", E.SPLITS, src("§7"), {k: cmp[k] for k in E.SPLITS}))
    checks.append(Check("A: counting link for n = 9", E.COUNTING, src("§8"),
                        {"cayley n^(n-2)": ta.cayley_count(9)[0], "(2n-3)!!": ta.rooted_binary_topologies(9)[0]}))

    b, _ = load_dataset(B)
    bs = b.seqs()
    kb, _ = mg.build_graph(bs, 10 ** 6)
    bmst = ta.kruskal(kb)[0]["edges"]
    (bids, bd), _ = mg.distance_matrix(bs)
    bu, _ = pa.upgma(bids, bd)
    heights = [m["height"] for m in bu["merges"]]
    checks.append(Check("B: UPGMA heights non-decreasing", True, "UPGMA property", heights == sorted(heights),
                        consistency=True))
    coph, _ = pa.cophenetic(bu, bids)
    checks.append(Check("B: UPGMA cophenetic matrix is ultrametric", True, "UPGMA property",
                        pa.is_ultrametric(bids, coph)[0], consistency=True))
    bok, bt = pa.is_ultrametric(bids, bd)
    bfit, _ = pa.fit_to_input(bu, bids, bd)
    bcmp, _ = pa.compare_splits(ta.edge_splits(kb.vertices, bmst)[0], pa.upgma_splits(bu, bids)[0])
    rooted = rooted_summary(kb.vertices, bmst, b.reference.id, False)
    obs.update(
        ultra=(bok, bt[0]["data"].get("witness")),
        fit={"exact": bfit["exact"], "pairs": bfit["pairs"], "max_deviation": str(bfit["max_deviation"]),
             "worst_pair": bfit["worst_pair"]},
        splits={k: bcmp[k] for k in ("mst", "upgma", "shared", "symmetric_difference")},
        ties=bu["tie_steps"], reverse_differs=pa.upgma(bids, bd, "reverse")[0]["newick"] != bu["newick"],
        rooted={"root": b.reference.id, "height": rooted["height"], "leaves": len(rooted["leaves"]),
                "internal": len(rooted["internal"])},
        newick=pa.newick(bu["tree"], bu["root"], pa.branch_lengths(bu)[0]) + ";")
    checks.append(Check("B: input distances ultrametric? (witness)", None, "", obs["ultra"], actual_only=True))
    checks.append(Check("B: cophenetic fit to d_H", None, "", obs["fit"], actual_only=True))
    checks.append(Check("B: UPGMA merge steps with ties; reverse rule changes topology", None, "",
                        {"tie_steps": obs["ties"], "reverse_differs": obs["reverse_differs"]}, actual_only=True))
    checks.append(Check("B: MST vs UPGMA non-trivial bipartitions", None, "", obs["splits"], actual_only=True))
    checks.append(Check("B: MST rooted at the medoid (height, #leaves, #internal)", None, "", obs["rooted"],
                        actual_only=True))
    return checks, obs


def main():
    checks, obs = run()
    f, sp = obs["fit"], obs["splits"]
    sections = [
        ("Question", "How do a minimum spanning tree and a UPGMA tree built from the same Hamming distances differ in "
         "structure, assumptions and what they can be used to claim?"),
        ("Input", "Dataset A (9 haplotypes) and Dataset B (30 haplotypes); complete distance graphs and distance "
         "matrices with entries d_H."),
        ("Mathematical model", "MST (Kruskal) on K_N, rooted at an arbitrary vertex for traversal; rooted-tree statistics; "
         "pre-, post- and level-order. UPGMA: average linkage, heights = distance / 2, tie-break D8 (lexicographic); "
         "in-order traversal of the binary tree; ultrametric test; cophenetic distances; non-trivial bipartitions of "
         "both trees."),
        ("Procedure", "Root the Dataset A MST at 01 and at 04 and compare statistics and traversals with the hand "
         "derivation. Run UPGMA and compare every merge, the Newick string, the in-order traversal, the ultrametric "
         "test and the cophenetic fit. Count shared bipartitions. Repeat on Dataset B and report the results as "
         "observations."),
        ("Interpretation",
         "The two trees answer different questions. The MST connects the observed sequences directly with the least "
         "total substitution count; its internal nodes are observed haplotypes, and its root is arbitrary (on Dataset A "
         "moving the root from 01 to 04 changes the height from 3 to 4 and the number of leaves from 6 to 5 without "
         "changing a single edge). UPGMA builds hypothetical ancestors (U1, U2, …) and places all sequences at the "
         "leaves under a molecular-clock assumption. Dataset A's distances are not ultrametric, and UPGMA reproduces "
         "only 3 of 36 distances exactly (an ultrametric input would be reproduced in full), and the two trees share 1 of their 2 and 6 non-trivial bipartitions. "
         f"Dataset B (observation): the distances are {'ultrametric' if obs['ultra'][0] else 'not ultrametric'}; the "
         f"UPGMA tree reproduces {f['exact']} of {f['pairs']} distances exactly (largest deviation {f['max_deviation']} at "
         f"{f['worst_pair']}); {obs['ties']} merge steps had ties, and the reverse tie rule "
         f"{'changes' if obs['reverse_differs'] else 'does not change'} the topology; the MST and UPGMA share "
         f"{sp['shared']} of their {sp['mst']} and {sp['upgma']} non-trivial bipartitions."),
        ("Limitations", "Neither tree is evidence of evolutionary history. The MST is a graph-theoretic object: it "
         "minimises total edge weight and nothing else, and its rooting is arbitrary. UPGMA assumes equal evolutionary "
         "rates (a molecular clock); where distances are not ultrametric its branch heights are averages, not times. "
         "Raw Hamming distance has no multiple-substitution correction. Ties are resolved by a fixed rule (D8), and a "
         "different rule can give a different tree. Established phylogenetic tools are out of scope."),
    ]
    notes = ["### Dataset B UPGMA tree (Newick; branch length = difference of UPGMA heights, same units as d_H; the path between two leaves sums to their cophenetic distance; observation)", "",
             f"`{obs['newick']}`"]
    out, summary = write(7, "Tree analysis: MST versus UPGMA", "python -m experiments.exp7_tree_analysis",
                         sections, checks, notes_after=["\n".join(notes)])
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
