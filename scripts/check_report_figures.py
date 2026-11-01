"""Re-check every figure quoted in the Stage 9 deliverables against the live code and the results files.

Documentation drifts. This script recomputes each number the report, proposal, design documents and viva answers
state, and fails if any disagrees. Run it after changing anything in genomath/ or the datasets.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomath import algebraic_structures as al
from genomath import counting_recursion as cr
from genomath import mutation_graph as mg
from genomath import number_theory as nt
from genomath import phylogenetic_analysis as pa
from genomath import tree_analysis as ta
from genomath._core import mismatch_count
from genomath.sequence_processing import load_dataset

PASSED, FAILED = [], []


def check(label, got, want):
    (PASSED if got == want else FAILED).append((label, got, want))


def results_totals():
    rows = []
    pattern = re.compile(r"\*\*(\d+)/(\d+) independent.*?(\d+)/(\d+) consistency.*?(\d+) observation")
    for path in sorted((ROOT / "results").glob("exp*.md")):
        line = next(x for x in path.read_text(encoding="utf-8").splitlines() if "independent expected" in x)
        rows.append(tuple(int(v) for v in pattern.match(line).groups()))
    return rows


def main():
    rows = results_totals()
    check("experiments", len(rows), 8)
    check("independent matched", sum(r[0] for r in rows), 147)
    check("independent attempted", sum(r[1] for r in rows), 147)
    check("consistency passed", sum(r[2] for r in rows), 28)
    check("consistency attempted", sum(r[3] for r in rows), 28)
    check("observations", sum(r[4] for r in rows), 35)

    public = [l for f in (ROOT / "genomath").glob("*.py") for l in f.read_text(encoding="utf-8").splitlines()
              if l.startswith("def ") and not l.startswith("def _")]
    check("public functions in genomath", len(public), 144)
    lines = sum(len(f.read_text(encoding="utf-8").splitlines()) for f in (ROOT / "genomath").glob("*.py"))
    check("genomath is about 2400 lines", 2300 <= lines <= 2500, True)
    check("words of length <= 6", sum(4 ** k for k in range(7)), 5461)
    check("mutants", len(re.findall(r'\("(?:genomath|app)/', (ROOT / "scripts" / "mutation_check.py").read_text(encoding="utf-8"))), 67)

    # ---------------- Dataset A ----------------
    a, _ = load_dataset(ROOT / "data" / "synthetic" / "dataset_a.fasta", "syn01")
    sa = {k[-2:]: v for k, v in a.seqs().items()}
    check("A: haplotypes from 10 records", (len(a.records), len(sa)), (10, 9))
    complete, _ = mg.build_graph(sa, 10 ** 6)
    kr, _ = ta.kruskal(complete)
    check("A: MST weight", kr["total"], 15)
    check("A: distinct MSTs", ta.count_minimum_spanning_trees(complete)[0], 8)
    check("A: k* and bottleneck", (mg.smallest_connecting_threshold(sa)[0], ta.bottleneck(kr["edges"])[0]), (6, 6))
    g2, _ = mg.build_graph(sa, 2)
    check("A: G_2 edges", len(g2.edges()), 14)
    check("A: G_2 degrees", list(mg.degrees(g2)[0].values()), [6, 5, 4, 2, 4, 4, 1, 2, 0])
    check("A: G_2 handshake", sum(mg.degrees(g2)[0].values()), 2 * len(g2.edges()))
    check("A: G_2 cyclomatic", mg.cyclomatic_number(g2)[0], 7)
    check("A: G_2 bipartite", mg.is_bipartite(g2)[0], False)
    check("A: components at k = 0,1,2,5,6",
          [r["components"] for r in mg.threshold_sweep(sa, [0, 1, 2, 5, 6])[0]], [9, 4, 2, 2, 1])
    check("A: 04->07 cost vs direct",
          (mg.path_versus_direct(g2, sa, "04", "07")[0]["cost"], mismatch_count(sa["04"], sa["07"])), (5, 4))
    check("A: 04->08 optimal paths", len(mg.path_versus_direct(g2, sa, "04", "08")[0]["paths"]), 3)
    ids, matrix = mg.distance_matrix(sa)[0]
    check("A: pairs", len(ids) * (len(ids) - 1) // 2, 36)
    up, _ = pa.upgma(ids, matrix)
    check("A: UPGMA tied steps", up["tie_steps"], 4)
    check("A: UPGMA later heights", [str(m["height"]) for m in up["merges"]][3:],
          ["3/4", "5/4", "19/12", "23/14", "61/16"])
    fit, _ = pa.fit_to_input(up, ids, matrix)
    check("A: UPGMA exact fit", (fit["exact"], fit["pairs"]), (3, 36))
    check("A: UPGMA worst deviation", (str(fit["max_deviation"]), fit["worst_pair"]), ("12/7", ("04", "08")))
    check("A: ultrametric", pa.is_ultrametric(ids, matrix)[0], False)
    check("A: orientation classes", len(al.dataset_orbits(sa, al.Z2)[0]), 8)
    d1 = [("01", "02"), ("01", "05"), ("02", "03"), ("02", "06"), ("03", "04"), ("05", "06")]
    check("A: m=16 collisions on d=1 pairs", len(nt.substitution_collisions(d1, sa, 16)[0]), 4)
    check("A: odd moduli collisions", [len(nt.substitution_collisions(d1, sa, m)[0]) for m in (7, 13, 17)], [0, 0, 0])
    for root, expected in (("01", (3, 6)), ("04", (4, 5))):
        tree, _ = ta.root_tree(complete.vertices, kr["edges"], root)
        stats, _ = ta.tree_statistics(tree)
        check(f"A: MST rooted at syn{root}", (stats["height"], len(stats["leaves"])), expected)

    # ---------------- Dataset B ----------------
    b, _ = load_dataset(ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta")
    sb = b.seqs()
    check("B: read / kept / haplotypes", (36, len(b.records), len(sb)), (36, 34, 30))
    check("B: dropped", [d[0] for d in b.dropped], ["MH246583.1", "MK181701.1"])
    check("B: length and reference", (b.n, b.reference.id), (1701, "CY249987.1"))
    ds = [mismatch_count(sb[x], sb[y]) for i, x in enumerate(sb) for y in list(sb)[i + 1:]]
    check("B: pairs / min / max", (len(ds), min(ds), max(ds)), (435, 1, 57))
    check("B: max distance is 3.4% of length", round(max(ds) / b.n * 100, 1), 3.4)
    check("B: pairs at distance 1", sum(1 for d in ds if d == 1), 8)
    check("B: digits of 4^n", cr.growth_table([b.n])[0][0]["digits"], 1025)
    check("B: possible neighbours 3n", cr.mutation_count(b.n, 1)[0], 5103)
    check("B: ball radius 1", cr.ball_size(b.n, 1)[0], 5104)
    check("B: max observed neighbours",
          max(sum(1 for o in sb if o != h and mismatch_count(sb[h], sb[o]) == 1) for h in sb), 6)
    g1, _ = mg.build_graph(sb, 1)
    check("B: G_1 edges / components / forest",
          (len(g1.edges()), len(mg.connected_components(g1)[0]), mg.find_cycle(g1)[0] is None), (8, 22, True))
    check("B: G_1 max degree", max(mg.degrees(g1)[0].values()), 6)
    kstar, _ = mg.smallest_connecting_threshold(sb)
    check("B: k*", kstar, 47)
    kb, _ = mg.build_graph(sb, 10 ** 6)
    kbr, _ = ta.kruskal(kb)
    pbr, _ = ta.prim(kb, b.reference.id)
    check("B: MST weight / count / largest edge",
          (kbr["total"], ta.count_minimum_spanning_trees(kb)[0], ta.bottleneck(kbr["edges"])[0]), (171, 4, 47))
    check("B: Prim and Kruskal differ",
          sorted(tuple(sorted(e[:2])) for e in pbr["edges"]) != sorted(tuple(sorted(e[:2])) for e in kbr["edges"]), True)
    check("B: Prim weight equals Kruskal", pbr["total"], kbr["total"])
    bids, bmatrix = mg.distance_matrix(sb)[0]
    bu, _ = pa.upgma(bids, bmatrix)
    bfit, _ = pa.fit_to_input(bu, bids, bmatrix)
    check("B: UPGMA fit / ties / ultrametric",
          (bfit["exact"], bfit["pairs"], bu["tie_steps"], pa.is_ultrametric(bids, bmatrix)[0]), (9, 435, 11, False))
    splits, _ = pa.compare_splits(ta.edge_splits(kb.vertices, kbr["edges"])[0], pa.upgma_splits(bu, bids)[0])
    check("B: bipartitions mst / upgma / shared", (splits["mst"], splits["upgma"], splits["shared"]), (10, 27, 5))
    codes = {h: nt.encode(v)[0] for h, v in sb.items()}
    check("B: buckets occupied at m=16 and 64",
          (nt.bucket_report(codes, 16)[0]["occupied"], nt.bucket_report(codes, 64)[0]["occupied"]), (1, 1))
    check("B: buckets at m=61",
          (nt.bucket_report(codes, 61)[0]["occupied"], nt.bucket_report(codes, 61)[0]["max_load"]), (24, 2))
    check("B: every haplotype ends TGA", sorted({v[-3:] for v in sb.values()}), ["TGA"])
    pairs1 = [(x, y) for i, x in enumerate(sb) for y in list(sb)[i + 1:] if mismatch_count(sb[x], sb[y]) == 1]
    check("B: d=1 collisions at 16 and 64",
          [len(nt.substitution_collisions(pairs1, sb, m)[0]) for m in (16, 64)], [8, 8])
    check("B: d=1 collisions at odd moduli",
          [len(nt.substitution_collisions(pairs1, sb, m)[0]) for m in (7, 13, 17, 31, 61)], [0, 0, 0, 0, 0])
    c4, _ = nt.kmer_census(list(sb.values()), 4)
    c8, _ = nt.kmer_census(list(sb.values()), 8)
    check("B: 4-mers distinct / canonical / whole space",
          (c4["distinct"], c4["canonical"], al.burnside(al.Z2, 4)[0]), (248, 136, 136))
    check("B: 8-mers distinct / canonical", (c8["distinct"], c8["canonical"]), (2563, 2488))
    vs_ref, _ = al.substitution_census((b.reference.seq, v) for h, v in sb.items() if h != b.reference.id)
    over_mst, _ = al.substitution_census((sb[u], sb[v]) for u, v, _ in kbr["edges"])
    check("B: transitions / transversions vs reference", (vs_ref["transitions"], vs_ref["transversions"]), (250, 91))
    check("B: transitions / transversions along MST", (over_mst["transitions"], over_mst["transversions"]), (131, 40))
    check("B: MST substitutions equal MST weight", over_mst["substitutions"], kbr["total"])
    check("B: largest forced k-mer", cr.max_forced_k(b.n)[0], 5)
    check("B: sequences with a repeated 5-mer",
          sum(cr.find_repeated_kmer(r.seq, 5)[0] is not None for r in b.records), 34)

    # ---------------- claims quoted in the viva answers ----------------
    check("fixed-point-free involutions", len(al.fixed_point_free_involutions()[0]), 3)
    check("ordered base pairs covered", al.ordered_pair_coverage()[0]["total"], 12)
    check("M(12,2) and ball(12,2)", (cr.mutation_count(12, 2)[0], cr.ball_size(12, 2)[0]), (594, 631))
    check("4^12", cr.sequence_space(12)[0], 16_777_216)
    monoid, _ = al.substitution_monoid(2)
    check("substitution monoid size and units", (monoid["size"], monoid["units"]), (25, ["··"]))
    check("monoid is not a group", monoid["report"]["group"], False)
    check("m=10 never collides on one substitution", nt.substitution_collision_t_min(10)[0], None)
    check("m=16 and 64 see j bases", (nt.last_bases_seen(16)[0], nt.last_bases_seen(64)[0]), (2, 3))
    check("orders of 4 mod 7 and 13", (nt.multiplicative_order(4, 7)[0], nt.multiplicative_order(4, 13)[0]), (3, 6))
    check("CRT (3,5) recovers only 15 of 64", sum(nt.crt([x % 3, x % 5], [3, 5])[0] == x for x in range(64)), 15)
    check("CRT (7,11) recovers all 64", sum(nt.crt([x % 7, x % 11], [7, 11])[0] == x for x in range(64)), 64)
    check("Burnside Z2 at n=2", al.burnside(al.Z2, 2)[0], 10)
    check("Burnside V4 at n=2", al.burnside(al.V4, 2)[0], 6)
    check("Cayley trees on 9 vertices", ta.cayley_count(9)[0], 4_782_969)
    cancel = "GATTAAGTAATG"
    base = nt.encode(sa["01"])[0]
    check("two-substitution cancellation",
          (mismatch_count(sa["01"], cancel), [m for m in (16, 64, 7, 13, 17, 31, 61)
                                              if nt.encode(cancel)[0] % m == base % m]), (2, [7, 13]))

    print(f"{len(PASSED)} figures verified against the live code")
    for label, got, want in FAILED:
        print(f"MISMATCH  {label}: got {got!r}, documents say {want!r}")
    if FAILED:
        sys.exit(f"{len(FAILED)} figure(s) in the documentation do not match the code")


if __name__ == "__main__":
    main()
