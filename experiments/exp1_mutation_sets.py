"""Experiment 1: mutation-set analysis on Dataset A (hand-checked) and Dataset B (invariants)."""
from collections import Counter
from itertools import combinations

from experiments import expected_dataset_a as E
from experiments._report import ROOT, Check, write
from genomath import set_logic as sl
from genomath._core import SequenceError, mismatch_count, reverse_complement
from genomath.sequence_processing import build_dataset, load_dataset, read_fasta

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
A_BAD = ROOT / "data" / "synthetic" / "dataset_a_invalid.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"
QC = "docs/dataset_b_qc.md (independent script)"


def components(seqs, k):
    ids = list(seqs)
    closure, _ = sl.transitive_closure(ids, sl.threshold_relation(seqs, k)[0])
    return sl.equivalence_classes(ids, closure)[0]


def run():
    checks = []
    src = lambda sec: f"{E.SOURCE} {sec}"

    # preprocessing
    a, _ = load_dataset(A, reference_id="syn01")
    checks.append(Check("A: haplotype ids", E.HAPLOTYPES, src("§3"), [h.id for h in a.haplotypes]))
    checks.append(Check("A: identical group", E.IDENTICAL_GROUP, src("§3"),
                        next(list(h.members) for h in a.haplotypes if len(h.members) > 1)))
    good, _ = read_fasta(A)
    bad, _ = read_fasta(A_BAD)
    d1, _ = build_dataset(good + bad[:1], "syn01")
    checks.append(Check("A + bad01: dropped record and reason", [("bad01", E.INVALID_BAD01_REASON)],
                        "data/synthetic/dataset_a_invalid.fasta", list(d1.dropped)))
    try:
        build_dataset(good + bad, "syn01")
        outcome = "accepted"
    except SequenceError as e:
        outcome = "rejected: " + str(e).split(";")[0]
    checks.append(Check("A + bad02 (length 11): outcome", "rejected: sequences have different lengths {11: 1, 12: 10}",
                        "final_spec.md §4 rule 3", outcome))

    # mutation sets and set operations
    R, seqs = a.reference.seq, a.seqs()
    mstar = {h: sorted(sl.allele_mutation_set(R, s)[0]) for h, s in seqs.items()}
    checks.append(Check("A: M*(R,S) for all 9 haplotypes", E.M_STAR, src("§1"), mstar))
    m = {h: sl.mutation_set(R, s)[0] for h, s in seqs.items()}
    ops, _ = sl.set_operations(m["syn03"], m["syn07"], ("M3", "M7"))
    checks.append(Check("A: M3 vs M7 set operations", E.SET_OPS_M3_M7, src("§2"),
                        {k: sorted(ops[k]) for k in E.SET_OPS_M3_M7}))
    allele_ops, _ = sl.set_operations(frozenset(mstar["syn03"]), frozenset(mstar["syn07"]))
    checks.append(Check("A: M*3 ∩ M*7", E.ALLELE_INTERSECTION_M3_M7, src("§2"), sorted(allele_ops["intersection"])))
    chain24, _ = sl.set_operations(m["syn02"], m["syn04"])
    chain43, _ = sl.set_operations(m["syn04"], m["syn03"])
    checks.append(Check("A: chain M2 ∩ M4, M4 \\ M3", E.CHAIN, src("§2"),
                        {"M2 ∩ M4": sorted(chain24["intersection"]), "M4 \\ M3": sorted(chain43["a_minus_b"])}))
    dec, _ = sl.distance_decomposition(R, seqs["syn03"], seqs["syn07"])
    checks.append(Check("A: d_H(syn03,syn07) decomposition", E.DECOMPOSITION_03_07, src("§2"),
                        {k: dec[k] for k in E.DECOMPOSITION_03_07}))

    # distances
    ids = E.HAPLOTYPES
    matrix = [[mismatch_count(seqs[x], seqs[y]) for y in ids] for x in ids]
    agree = sum(matrix[i][j] == E.D[i][j] for i in range(9) for j in range(9))
    checks.append(Check("A: d_H matrix entries equal to hand matrix", 81, src("§3"), agree))
    holds_a = sum(sl.distance_decomposition(R, s, t)[0]["holds"] for s, t in combinations(seqs.values(), 2))
    checks.append(Check("A: pairs where the d_H identity holds", 36, "math_audit.md Unit I (theorem)", holds_a, consistency=True))

    # relation R_1
    r1, _ = sl.threshold_relation(seqs, 1)
    props, _ = sl.relation_properties(ids, r1)
    checks.append(Check("A: |R_1| (ordered pairs)", E.R1_SIZE, src("§3"), len(r1)))
    checks.append(Check("A: R_1 reflexive / symmetric / transitive", E.R1_PROPERTIES, src("§3"),
                        {p: props[p]["holds"] for p in E.R1_PROPERTIES}))
    checks.append(Check("A: R_1 transitivity counterexample", E.R1_COUNTEREXAMPLE, src("§3"),
                        props["transitive"]["witness"]))
    checks.append(Check("A: components of closure(R_k) for k in {0,1,2,5,6}", E.COMPONENTS, src("§4"),
                        {k: len(components(seqs, k)) for k in E.COMPONENTS}))

    # predicates
    P, Q = sl.within_distance(R, 2), sl.base_at(2, "C")
    preds = {"P ∧ Q": P & Q, "P ∧ ¬Q": P & ~Q, "¬P": ~P}
    checks.append(Check("A: predicate filters", E.PREDICATES, src("§2"),
                        {k: sl.filter_sequences(seqs, p)[0] for k, p in preds.items()}))

    # orientation equivalence on records
    rec = {r.id: r.seq for r in a.records}
    oclasses, _ = sl.equivalence_classes(list(rec), sl.orientation_relation(rec)[0])
    checks.append(Check("A: orientation class of syn03", E.ORIENTATION_CLASS_OF_SYN03, src("§6"),
                        next(c for c in oclasses if "syn03" in c)))
    # RC-fixed means RC(S) = S; a singleton orientation class may just mean the partner is absent from the data
    checks.append(Check("A: RC-fixed records", E.RC_FIXED, src("§6"),
                        [r for r, s in rec.items() if reverse_complement(s) == s]))

    # Dataset B
    b, _ = load_dataset(B)
    bseqs, bR = b.seqs(), b.reference.seq
    checks.append(Check("B: dropped records", ["MH246583.1", "MK181701.1"], QC, [d[0] for d in b.dropped]))
    checks.append(Check("B: haplotypes", 30, QC, len(b.haplotypes)))
    n_pairs = len(b.haplotypes) * (len(b.haplotypes) - 1) // 2
    holds_b = sum(sl.distance_decomposition(bR, s, t)[0]["holds"] for s, t in combinations(bseqs.values(), 2))
    checks.append(Check("B: pairs where the d_H identity holds", n_pairs, "math_audit.md Unit I (theorem)", holds_b, consistency=True))
    br1, _ = sl.threshold_relation(bseqs, 1)
    checks.append(Check("B: |R_1| = 30 + 2 × (d_H = 1 pairs)", 30 + 2 * 8, QC, len(br1)))
    bprops, _ = sl.relation_properties(list(bseqs), br1)
    checks.append(Check("B: R_1 transitive? (witness)", None, "", (bprops["transitive"]["holds"],
                        bprops["transitive"]["witness"]), actual_only=True))
    checks.append(Check("B: medoid reference", None, "", b.reference.id, actual_only=True))
    sizes = sorted(Counter(len(sl.mutation_set(bR, s)[0]) for s in bseqs.values()).items())
    checks.append(Check("B: |M(R,S)| distribution (size: count)", None, "", dict(sizes), actual_only=True))
    return checks, {"b": b, "bprops": bprops, "sizes": sizes, "a_matrix": matrix}


def main():
    checks, extra = run()
    b, bprops, sizes = extra["b"], extra["bprops"], extra["sizes"]
    w = bprops["transitive"]["witness"]
    sections = [
        ("Question", "Do mutation sets, set operations, Hamming distances and the threshold relation R_k computed by "
         "`genomath` reproduce the hand-derived values for Dataset A, and do the stated identities hold on real data?"),
        ("Input", "Dataset A (`data/synthetic/dataset_a.fasta`, reference syn01) plus the two records in "
         "`dataset_a_invalid.fasta`; Dataset B (`data/public_genomic/h3n2_ha_cds.fasta`, medoid reference)."),
        ("Mathematical model", "M(R,S) = {i : R_i ≠ S_i}; M*(R,S) = {(i, S_i) : R_i ≠ S_i}; "
         "d_H(S,T) = |M(R,S) Δ M(R,T)| + |{i ∈ M(R,S) ∩ M(R,T) : S_i ≠ T_i}|; S R_k T ⇔ d_H(S,T) ≤ k; "
         "components of G_k = classes of the transitive closure of R_k (Warshall)."),
        ("Procedure", "Load and preprocess each dataset (uppercase, drop non-ACGT, require one length, collapse "
         "haplotypes, choose reference). Compute every quantity in the Expected table with `genomath.set_logic` and "
         "compare against values fixed beforehand. Dataset B has no hand values: its expected entries come from the "
         "independent QC script or from theorems, and items without an independent value are reported as observations."),
        ("Interpretation",
         "Dataset A: the shared-position example (syn03 vs syn07) shows why position-level sets are not enough: "
         "they share position 5 but no mutation, and that position still contributes 1 to d_H. R_1 is reflexive and "
         "symmetric but not transitive, so it is a similarity relation, not an equivalence; its transitive closure is "
         "an equivalence whose classes are the components of the threshold graph.\n\n"
         f"Dataset B (observation, no biological claim): relative to the medoid {b.reference.id}, mutation-set sizes "
         f"are distributed as {dict(sizes)}. R_1 is "
         + ("transitive on this sample." if bprops["transitive"]["holds"] else
            f"not transitive on real data either: {w[0]} R_1 {w[1]} and {w[1]} R_1 {w[2]} but not {w[0]} R_1 {w[2]}.")),
        ("Limitations", "Substitutions only; positional comparison assumes equal-length, colinear sequences "
         "(Dataset B: equal CDS length and ≤ 3.4% divergence, not an alignment). Mutation sets depend on the chosen "
         "reference; pairwise distances do not. Two Dataset B records with ambiguity codes were dropped, not resolved."),
    ]
    out, summary = write(1, "Mutation-set analysis", "python -m experiments.exp1_mutation_sets", sections, checks)
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
