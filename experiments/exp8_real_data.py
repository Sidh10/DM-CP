"""Experiment 8: the validated pipeline applied end to end to Dataset B, all six units."""
from itertools import combinations

from experiments import expected_dataset_b as E
from experiments._report import ROOT, Check, write
from genomath import algebraic_structures as al
from genomath import counting_recursion as cr
from genomath import mutation_graph as mg
from genomath import number_theory as nt
from genomath import phylogenetic_analysis as pa
from genomath import set_logic as sl
from genomath import tree_analysis as ta
from genomath._core import mismatch_count
from genomath.sequence_processing import load_dataset, read_fasta

B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"
MODULI = [16, 64, 7, 13, 17, 31, 61]
KS = list(range(0, 11)) + [15, 20, 30, 40, 47, 57]


def run():
    checks, obs = [], {}

    # ---------- preprocessing ----------
    records, _ = read_fasta(B)
    b, _ = load_dataset(B)
    seqs = b.seqs()
    n, ref = b.n, b.reference.id
    checks += [
        Check("records read from FASTA", E.RECORDS_READ, E.QC, len(records)),
        Check("records dropped for non-ACGT characters", E.DROPPED, E.QC, [d[0] for d in b.dropped]),
        Check("clean records / haplotypes / common length", [E.CLEAN, E.HAPLOTYPES, E.LENGTH], E.QC,
              [len(b.records), len(seqs), n]),
        Check("groups of identical sequences", E.IDENTICAL_GROUPS, E.QC,
              [list(h.members) for h in b.haplotypes if len(h.members) > 1]),
    ]
    obs["reference"] = ref

    # ---------- Unit I ----------
    R = b.reference.seq
    dists = {(x, y): mismatch_count(seqs[x], seqs[y]) for x, y in combinations(seqs, 2)}
    d1 = [p for p, v in dists.items() if v == 1]
    checks += [
        Check("Unit I: haplotype pairs; pairs at d_H = 1", [E.PAIRS, E.D1_PAIRS], E.QC, [len(dists), len(d1)]),
        Check("Unit I: min and max pairwise d_H", E.DISTANCE_RANGE, E.QC,
              {"min": min(dists.values()), "max": max(dists.values())}),
        Check("Unit I: d_H identity holds on every pair", E.PAIRS, "math_audit.md Unit I (theorem)",
              sum(sl.distance_decomposition(R, seqs[x], seqs[y])[0]["holds"] for x, y in dists), consistency=True),
    ]
    r1, _ = sl.threshold_relation(seqs, 1)
    props, _ = sl.relation_properties(list(seqs), r1)
    obs["r1"] = {"reflexive": props["reflexive"]["holds"], "symmetric": props["symmetric"]["holds"],
                 "transitive": props["transitive"]["holds"], "witness": props["transitive"]["witness"]}
    checks.append(Check("Unit I: R_1 properties and non-transitivity witness", None, "", obs["r1"], actual_only=True))
    msizes = sorted(len(sl.mutation_set(R, s)[0]) for s in seqs.values())
    obs["mutation_sizes"] = {"min": msizes[0], "median": msizes[len(msizes) // 2], "max": msizes[-1]}

    # ---------- Unit II ----------
    checks += [
        Check("Unit II: decimal digits of 4^n", E.DIGITS_4_POW_N, f"{E.U2} §1", cr.growth_table([n])[0][0]["digits"]),
        Check("Unit II: |B(n,1)| and 3n possible neighbours", [E.BALL_RADIUS_1, E.NEIGHBOURS_AT_1], f"{E.U2} §2",
              [cr.ball_size(n, 1)[0], cr.mutation_count(n, 1)[0]]),
        Check("Unit II: largest forced repeated k-mer length", E.MAX_FORCED_KMER, f"{E.U2} §5", cr.max_forced_k(n)[0]),
        Check("Unit II: clean sequences containing a repeated 5-mer", E.FORCED_5MER_SEQUENCES, f"{E.U2} §5",
              sum(cr.find_repeated_kmer(r.seq, 5)[0] is not None for r in b.records)),
    ]
    deg1 = {h: sum(1 for o in seqs if o != h and mismatch_count(seqs[h], seqs[o]) == 1) for h in seqs}
    obs["sparsity"] = {"max_observed_neighbours": max(deg1.values()), "possible": cr.mutation_count(n, 1)[0],
                       "fraction": f"{max(deg1.values()) / cr.mutation_count(n, 1)[0]:.2e}"}

    # ---------- Unit III ----------
    values = set(seqs.values())
    checks += [
        Check("Unit III: RC pairs present / RC-fixed / {I,RC} classes", [E.RC_PAIRS, E.RC_FIXED, E.Z2_CLASSES],
              f"{E.U3} §8", [sum(al.TRANSFORMS["RC"](x) in values and al.TRANSFORMS["RC"](x) != x for x in values) // 2,
                             sum(al.TRANSFORMS["RC"](x) == x for x in values), len(al.dataset_orbits(seqs, al.Z2)[0])]),
        Check("Unit III: canonical 4-mer classes in the whole space (Burnside)", E.SPACE_CANONICAL_4MERS,
              f"{E.U4} §2", al.burnside(al.Z2, 4)[0]),
    ]
    census4, _ = nt.kmer_census(list(values), 4)
    obs["kmer4"] = {k: census4[k] for k in ("distinct", "canonical", "palindromic")}

    # D9
    fpf, _ = al.fixed_point_free_involutions()
    table, _ = al.base_cayley_table()
    axioms, _ = al.check_axioms(al.BASE_K4, table)
    cov, _ = al.ordered_pair_coverage()
    unordered = [(x, y) for x in "ACGT" for y in "ACGT" if x < y]
    kinds = [al.classify_substitution(x, y)[0]["kind"] for x, y in unordered]
    checks += [
        Check("D9: fixed-point-free involutions of Σ", E.FPF_INVOLUTIONS, f"{E.U3} §9", len(fpf)),
        Check("D9: base-level group is an abelian group of order 4, every non-identity element an involution",
              [True, True, [1, 2, 2, 2]], f"{E.U3} §9",
              [axioms["group"], axioms["commutative"]["holds"],
               sorted(al.element_orders(al.BASE_K4, table, "id")[0].values())]),
        Check("D9: ordered base pairs covered, each by exactly one element", E.ORDERED_PAIRS_COVERED, f"{E.U3} §9",
              cov["total"] if cov["covers_all"] and cov["distinct"] == cov["total"] else -1),
        Check("D9: transitions / transversions among the 6 unordered pairs",
              [E.TRANSITION_PAIRS, E.TRANSVERSION_PAIRS], f"{E.U3} §9",
              [kinds.count("transition"), kinds.count("transversion")]),
    ]
    vs_ref, _ = al.substitution_census((R, s) for h, s in seqs.items() if h != ref)
    checks.append(Check("D9: census by element sums to the substitution count", True, "definition",
                        sum(vs_ref["by_element"].values()) == vs_ref["substitutions"], consistency=True))
    checks.append(Check("D9: combinatorial baseline ratio", E.BASELINE_RATIO, f"{E.U3} §9", vs_ref["baseline_ratio"]))

    # ---------- Unit IV ----------
    codes = {h: nt.encode(s)[0] for h, s in seqs.items()}
    checks += [
        Check("Unit IV: decode(encode(S)) = S for every haplotype", True, "bijection (unit4 §1)",
              all(nt.decode(c, n)[0] == seqs[h] for h, c in codes.items()), consistency=True),
        Check("Unit IV: smallest j with m | 4^j (bases the bucket can see)", E.LAST_BASES_SEEN, f"{E.U4} §3",
              {m: nt.last_bases_seen(m)[0] for m in MODULI}),
        Check("Unit IV: d_H = 1 pairs sharing a bucket, by modulus", E.D1_COLLISIONS, f"{E.U4} §3",
              {m: len(nt.substitution_collisions(d1, seqs, m)[0]) for m in MODULI}),
    ]
    reports = {m: nt.bucket_report(codes, m)[0] for m in MODULI}
    obs["buckets"] = {m: {"occupied": r["occupied"], "max_load": r["max_load"], "bound": r["pigeonhole_bound"]}
                      for m, r in reports.items()}
    checks.append(Check("Unit IV: max bucket load ≥ pigeonhole bound for every modulus", True, "pigeonhole",
                        all(r["max_load"] >= r["pigeonhole_bound"] for r in reports.values()), consistency=True))

    # ---------- Unit V ----------
    sweep, _ = mg.threshold_sweep(seqs, KS)
    comps = [r["components"] for r in sweep]
    kstar, _ = mg.smallest_connecting_threshold(seqs)
    g1, _ = mg.build_graph(seqs, 1)
    gk, _ = mg.build_graph(seqs, kstar)
    checks += [
        Check("Unit V: edges at k = 1 equal the QC count of d_H = 1 pairs", E.D1_PAIRS, E.QC, len(g1.edges())),
        Check("Unit V: components never increase as k grows", True, "G_j ⊆ G_k", all(x >= y for x, y in zip(comps, comps[1:])),
              consistency=True),
    ]
    obs["kstar"] = kstar
    obs["g1"] = {"edges": len(g1.edges()), "components": len(mg.connected_components(g1)[0]),
                 "cyclomatic": mg.cyclomatic_number(g1)[0], "max_degree": max(mg.degrees(g1)[0].values()),
                 "is_forest": mg.find_cycle(g1)[0] is None}
    obs["sweep"] = sweep
    detours = {}
    all_ge = True
    for s in seqs:
        dist = mg.dijkstra(gk, s)[0]["dist"]
        for t in seqs:
            all_ge &= dist[t] is not None and dist[t] >= mismatch_count(seqs[s], seqs[t])
            if s == ref and dist[t] > mismatch_count(seqs[s], seqs[t]):
                detours[t] = dist[t] - mismatch_count(seqs[s], seqs[t])
    checks.append(Check(f"Unit V: path cost ≥ direct d_H for all 900 ordered pairs at k* = {kstar}", True,
                        "triangle inequality", all_ge, consistency=True))
    obs["detours"] = detours

    # ---------- Unit VI ----------
    kn, _ = mg.build_graph(seqs, 10 ** 6)
    kr, _ = ta.kruskal(kn)
    pr, _ = ta.prim(kn, ref)
    checks += [
        Check("Unit VI: Prim total = Kruskal total", True, "MST weight is unique", pr["total"] == kr["total"],
              consistency=True),
        Check("Unit VI: both outputs are valid spanning trees", True, "definition",
              all(ta.validate_spanning_tree(kn.vertices, x["edges"], kn)[0]["is_spanning_tree"] for x in (kr, pr)),
              consistency=True),
        Check("Unit VI: largest MST edge = k* (bottleneck property)", True, "docs/unit6_expected.md §3",
              ta.bottleneck(kr["edges"])[0] == kstar, consistency=True),
    ]
    obs["mst"] = {"weight": kr["total"], "count": ta.count_minimum_spanning_trees(kn)[0],
                  "max_edge": ta.bottleneck(kr["edges"])[0]}
    (ids, dm), _ = mg.distance_matrix(seqs)
    up, _ = pa.upgma(ids, dm)
    coph, _ = pa.cophenetic(up, ids)
    fit, _ = pa.fit_to_input(up, ids, dm)
    splits, _ = pa.compare_splits(ta.edge_splits(kn.vertices, kr["edges"])[0], pa.upgma_splits(up, ids)[0])
    checks.append(Check("Unit VI: UPGMA cophenetic distances are ultrametric", True, "UPGMA property",
                        pa.is_ultrametric(ids, coph)[0], consistency=True))
    obs["upgma"] = {"ultrametric_input": pa.is_ultrametric(ids, dm)[0], "exact_fit": fit["exact"],
                    "pairs": fit["pairs"], "max_deviation": str(fit["max_deviation"]), "tie_steps": up["tie_steps"]}
    obs["splits"] = {k: splits[k] for k in ("mst", "upgma", "shared", "symmetric_difference")}

    # ---------- D9 on the pipeline's own structures ----------
    mst_pairs = [(seqs[u], seqs[v]) for u, v, _ in kr["edges"]]
    over_mst, _ = al.substitution_census(mst_pairs)
    # one substitution per differing position per edge, so the census total must be the MST weight
    checks.append(Check("D9 × Unit VI: substitutions along MST edges = MST weight", kr["total"],
                        "Σ d_H over tree edges", over_mst["substitutions"], consistency=True))
    obs["census"] = {
        "vs_reference": {k: vs_ref[k] for k in ("substitutions", "transitions", "transversions", "by_base_pair")},
        "vs_reference_ratio": round(vs_ref["ratio"], 3) if vs_ref["ratio"] else None,
        "over_mst_edges": {k: over_mst[k] for k in ("substitutions", "transitions", "transversions", "by_base_pair")},
        "over_mst_ratio": round(over_mst["ratio"], 3) if over_mst["ratio"] else None,
    }
    for label, key in (("differences from the reference (29 pairs)", "vs_reference"),
                       ("substitutions along MST edges (29 edges)", "over_mst_edges")):
        checks.append(Check(f"D9: {label}", None, "", obs["census"][key], actual_only=True))

    for item, value in (("reference (medoid)", ref), ("mutation-set sizes vs the reference", obs["mutation_sizes"]),
                        ("observed vs possible single-substitution neighbours", obs["sparsity"]),
                        ("distinct / canonical / RC-palindromic 4-mers", obs["kmer4"]),
                        ("bucket occupancy per modulus", obs["buckets"]),
                        ("G_=1 structure", obs["g1"]), ("smallest connecting threshold k*", kstar),
                        ("detours from the reference at k*", detours),
                        ("MST weight / number of MSTs / largest edge", obs["mst"]),
                        ("UPGMA fit and ties", obs["upgma"]),
                        ("MST vs UPGMA bipartitions", obs["splits"])):
        checks.append(Check(f"observed: {item}", None, "", value, actual_only=True))
    return checks, obs


def main():
    checks, obs = run()
    c, g1, up = obs["census"], obs["g1"], obs["upgma"]
    sweep = ["### Threshold sweep (observation)", "",
             "| k | edges | components | largest | isolated | max degree |", "|---|---|---|---|---|---|"]
    sweep += [f"| {r['k']} | {r['edges']} | {r['components']} | {r['largest']} | {r['isolated']} | {r['max_degree']} |"
              for r in obs["sweep"]]
    ts_table = ["### Substitutions by unordered base pair (observation, direction-free)", "",
                "| view | substitutions | transitions | transversions | ratio | by pair |", "|---|---|---|---|---|---|"]
    for label, key, rk in (("vs reference", "vs_reference", "vs_reference_ratio"),
                           ("along MST edges", "over_mst_edges", "over_mst_ratio")):
        v = c[key]
        ts_table.append(f"| {label} | {v['substitutions']} | {v['transitions']} | {v['transversions']} | "
                        f"{c[rk]} | {v['by_base_pair']} |")
    sections = [
        ("Question", "Does the pipeline validated on Dataset A run end to end on real genomic data, reproducing the "
         "independently known quality-control figures and the hand-proved theorems, and what does it then observe about "
         "the 30 influenza A(H3N2) HA haplotypes?"),
        ("Input", "`data/public_genomic/h3n2_ha_cds.fasta`: 36 GenBank hemagglutinin coding sequences, human host, "
         "2017–18 northern-hemisphere season, each 1701 nt. Accessions in `data/public_genomic/accessions.csv`."),
        ("Mathematical model", "All six units in one pass: mutation sets and the threshold relation R_k (I); sequence-space "
         "and Hamming-ball counts and the pigeonhole bound (II); the {I,RC} action and the base-level Klein group "
         "{id, (A T)(C G), (A G)(C T), (A C)(G T)} (III, decision D9); base-4 encoding and modular bucketing (IV); "
         "threshold graphs, components and Dijkstra (V); Prim, Kruskal and UPGMA (VI)."),
        ("Procedure", "Load and preprocess the FASTA (uppercase, drop non-ACGT, require one length, collapse haplotypes, "
         "choose the medoid reference). Then run each unit's algorithms in turn. Independent entries are compared with "
         "the genomath-free QC script or with hand derivations; entries that only check genomath against itself are "
         "tagged as consistency checks; everything else is reported as an observation."),
        ("Interpretation",
         "The pipeline runs end to end on real data and reproduces every independently known figure: the QC counts "
         "(36 read, 2 dropped for ambiguity codes, 34 clean, 30 haplotypes, 1701 nt, 8 pairs at distance 1, distances "
         "1 to 57), the Unit II hand values (4^1701 has 1025 digits; 5103 possible single-substitution neighbours; a "
         "repeated 5-mer is forced in all 34 sequences), the Unit III prediction that no sequence's reverse complement "
         "is present, and both Unit IV hashing theorems (m = 16 and 64 collide on all 8 distance-1 pairs; every odd "
         "modulus separates all 8).\n\n"
         "Observations about this sample, with no biological claim:\n"
         f"- **Sparsity.** The best-connected haplotype has {obs['sparsity']['max_observed_neighbours']} observed "
         f"single-substitution neighbours out of {obs['sparsity']['possible']} possible "
         f"({obs['sparsity']['fraction']}), so the single-substitution graph is an extremely thin sample of the "
         "Hamming graph.\n"
         f"- **Structure.** G_=1 has {g1['edges']} edges in {g1['components']} components and is "
         f"{'a forest' if g1['is_forest'] else 'not a forest'}; the graph first becomes connected at k* = {obs['kstar']}, "
         f"which equals the largest edge of the minimum spanning tree (weight {obs['mst']['weight']}, "
         f"{obs['mst']['count']} distinct MSTs).\n"
         f"- **Tree comparison.** The distances are {'ultrametric' if up['ultrametric_input'] else 'not ultrametric'}, "
         f"so UPGMA reproduces only {up['exact_fit']} of {up['pairs']} distances exactly; {up['tie_steps']} merge steps "
         f"involved ties, and the MST and UPGMA share {obs['splits']['shared']} of their {obs['splits']['mst']} and "
         f"{obs['splits']['upgma']} non-trivial bipartitions.\n"
         f"- **Substitution classes (D9).** Direction-free counts give {c['vs_reference']['transitions']} transitions "
         f"and {c['vs_reference']['transversions']} transversions among differences from the reference "
         f"(ratio {c['vs_reference_ratio']}), and {c['over_mst_edges']['transitions']} to "
         f"{c['over_mst_edges']['transversions']} along MST edges (ratio {c['over_mst_ratio']}). Both are far from the "
         "combinatorial baseline of 0.5, which is what equal likelihood over the 6 unordered base pairs would give. "
         "The counts are what they are; **explaining** the gap would need a substitution model and a source for what "
         "ratio real data shows, and this project has neither, so no conclusion is drawn from it."),
        ("Limitations", "All six units' limitations apply at once and are listed in Experiments 1–7. Specific to this "
         "run: the 30 haplotypes come from one season and are dominated by United States deposits, so they are a "
         "convenience sample and not a random sample of anything; positional comparison rests on equal CDS length and "
         "low divergence, not on an alignment; two records with ambiguity codes were dropped rather than resolved; "
         "substitution counts are direction-free, because no ancestor is known, and counting along MST edges counts "
         "one chosen tree's steps, not inferred events; the reference is the medoid, a mathematical choice with no "
         "biological standing; neither the graph nor either tree is evidence of evolutionary history; and nothing here "
         "supports any clinical, diagnostic or epidemiological claim."),
    ]
    out, summary = write(8, "Real-data evaluation of the full pipeline", "python -m experiments.exp8_real_data",
                         sections, checks, notes_after=["\n".join(ts_table), "\n".join(sweep)])
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
