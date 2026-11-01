"""Experiment 2: counting versus enumeration (Unit II)."""
import time

from experiments import expected_unit2 as E
from experiments._report import ROOT, Check, write
from genomath import counting_recursion as cr
from genomath._core import mismatch_count
from genomath.sequence_processing import load_dataset

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"
QC = "docs/dataset_b_qc.md (independent script)"
REF = "GATTACGT"  # first 8 bases of syn01; any reference gives the same counts


def run():
    src = lambda sec: f"{E.SOURCE} {sec}"
    checks = []
    enum_table, closed_table, valid_cells, timings = [], [], 0, []
    for n in range(0, 9):
        ref = REF[:n]
        erow, crow = [], []
        for k in range(0, n + 1):
            t0 = time.perf_counter()
            got, trace = cr.generate_mutants(ref, k)
            dt = time.perf_counter() - t0
            erow.append(len(got))
            crow.append(cr.mutation_count(n, k)[0])
            valid_cells += len(set(got)) == len(got) and all(mismatch_count(ref, w) == k for w in got)
            if n == 8:
                calls = next(s["data"]["calls"] for s in trace if "calls" in s["data"])
                timings.append((k, len(got), calls, dt))
        enum_table.append(erow)
        closed_table.append(crow)
    dp, _ = cr.recurrence_table(8)
    checks += [
        Check("M(n,k) by recursive enumeration, n ≤ 8", E.M_TABLE, src("§2"), enum_table),
        Check("M(n,k) by closed form C(n,k)·3^k, n ≤ 8", E.M_TABLE, src("§2"), closed_table),
        Check("M(n,k) by recurrence table, n ≤ 8", E.M_TABLE, src("§2"), dp),
        Check("cells where all strings are distinct and at distance exactly k", E.CELLS, src("§2"), valid_cells, consistency=True),
        Check("row sums Σ_k M(n,k), n ≤ 8", E.POWERS_OF_4, src("§1"), [sum(r) for r in enum_table]),
        Check("M(12,k) for k = 1,2,3", E.M12, "docs/synthetic_expected.md §8",
              {k: cr.mutation_count(12, k)[0] for k in E.M12}),
        Check("P(12,3)", E.P_12_3, src("§3"), cr.permutations(12, 3)[0]),
        Check("Hamming-ball sizes |B(n,r)|", {f"B{nr}": v for nr, v in E.BALLS.items()}, src("§2"),
              {f"B{nr}": cr.ball_size(*nr)[0] for nr in E.BALLS}),
        Check("enumerated ball B(8,8) = whole Σ^8", E.WHOLE_SPACE_8, src("§2"), len(set(cr.generate_within(REF, 8)[0]))),
        Check("decimal digits of 4^1701", E.DIGITS_4_1701, src("§1"), cr.growth_table([1701])[0][0]["digits"]),
    ]

    guard = {}
    for n, k in E.GUARD:
        try:
            got, _ = cr.generate_mutants("A" * n, k)
            guard[(n, k)] = ("allowed", len(got))
        except cr.GenerationLimitError as e:
            guard[(n, k)] = ("refused", e.count)
    checks.append(Check("guard outcome (cap 10^5) and count", {str(x): v for x, v in E.GUARD.items()}, src("§2"),
                        {str(x): v for x, v in guard.items()}))

    # the ball of radius n is the whole space Σ^n, so this enumerates every string of length n
    by_enum = {n: sum(1 for w in cr.generate_within(REF[:n], n)[0] if set(w) == set("ACGT")) for n in (4, 5, 6, 7, 8)}
    checks.append(Check("strings using all four bases, by inclusion–exclusion",
                        E.ALL_FOUR, src("§4"), {n: cr.inclusion_exclusion_all_bases(n)[0] for n in E.ALL_FOUR}))
    checks.append(Check("strings using all four bases, by enumeration (n = 4..8)",
                        {n: E.ALL_FOUR[n] for n in by_enum}, src("§4"), by_enum))

    a, _ = load_dataset(A, reference_id="syn01")
    s = a.seqs()
    orders, _ = cr.mutation_orderings(s["syn01"], s["syn04"])
    observed = set(s.values())
    fully = sum(all(x in observed for x in p) for p in orders["paths"])
    checks.append(Check("syn01 → syn04: shortest orderings, fully observed ones", E.ORDERINGS_01_04, src("§3"),
                        {"count": orders["count"], "fully_observed": fully}))

    b, _ = load_dataset(B)
    with_repeat = sum(cr.find_repeated_kmer(r.seq, 5)[0] is not None for r in b.records)
    checks.append(Check("pigeonhole on L = 1701", E.PIGEON, src("§5"), {
        "max_forced_k": cr.max_forced_k(1701)[0], "forced_k5": cr.forced_repeat(1701, 5)[0],
        "forced_k6": cr.forced_repeat(1701, 6)[0], "clean_with_repeat_5mer": with_repeat}))
    with_6 = sum(cr.find_repeated_kmer(r.seq, 6)[0] is not None for r in b.records)
    checks.append(Check("clean B sequences with a repeated 6-mer (not forced)", None, "", f"{with_6}/34",
                        actual_only=True))

    hs = b.seqs()
    deg1 = {h: sum(1 for o in hs if o != h and mismatch_count(hs[h], hs[o]) == 1) for h in hs}
    deg2 = {h: sum(1 for o in hs if o != h and mismatch_count(hs[h], hs[o]) <= 2) for h in hs}
    checks.append(Check("B: Σ observed neighbours at d = 1 (= 2 × pairs)", 16, QC, sum(deg1.values())))
    checks.append(Check("B: Σ observed neighbours at d ≤ 2 (= 2 × pairs)", 62, QC, sum(deg2.values())))
    possible1 = cr.ball_size(1701, 1)[0] - 1
    possible2 = cr.ball_size(1701, 2)[0] - 1
    checks.append(Check("B: max observed vs possible neighbours", None, "",
                        {"d=1": f"{max(deg1.values())} of {possible1}", "d≤2": f"{max(deg2.values())} of {possible2}"},
                        actual_only=True))
    return checks, {"timings": timings, "possible1": possible1, "possible2": possible2, "deg1": deg1}


def main():
    checks, x = run()
    tim = ["### Enumeration cost at n = 8 (observation; wall-clock time depends on the machine)", "",
           "| k | strings M(8,k) | recursive calls | time (ms) |", "|---|---|---|---|"]
    tim += [f"| {k} | {m} | {c} | {dt * 1000:.1f} |" for k, m, c, dt in x["timings"]]
    sections = [
        ("Question", "Does recursive enumeration of all sequences at Hamming distance exactly k from a reference produce "
         "exactly C(n,k)·3^k distinct strings, matching the recurrence M(n,k) = M(n−1,k) + 3M(n−1,k−1)? "
         "Where does enumeration stop being feasible, and how sparse are real observations within the possible space?"),
        ("Input", f"Reference prefixes of {REF} for n = 0..8; Dataset A (syn01, syn04) for orderings; "
         "Dataset B (34 clean sequences, 30 haplotypes) for pigeonhole and sparsity."),
        ("Mathematical model", "N(n) = 4^n; M(n,k) = C(n,k)·3^k; M(n,k) = M(n−1,k) + 3M(n−1,k−1), M(0,0) = 1; "
         "|B(n,r)| = Σ_{k≤r} M(n,k); inclusion–exclusion U(n) = 4^n − 4·3^n + 6·2^n − 4; pigeonhole L − k + 1 > 4^k; "
         "k! shortest substitution orderings."),
        ("Procedure", "For each (n,k) with n ≤ 8, generate all mutants recursively (choose k positions in increasing "
         "order, then 3 bases each), check that they are distinct and at distance exactly k, and count them. Compare "
         "with the hand table, the closed form and the recurrence table. Exercise the 10^5 enumeration cap. Count "
         "all-four-base strings by inclusion–exclusion and by enumeration. Apply the pigeonhole bound to every clean "
         "Dataset B sequence. Compare observed neighbour counts in Dataset B with the Hamming-ball sizes."),
        ("Interpretation",
         "Enumeration, closed form and recurrence agree on every cell for n ≤ 8; the induction proof in "
         "`docs/unit2_expected.md` §2 extends this to all n. The cost of enumeration is proportional to the number "
         "of strings, so the cap is reached quickly: M(12,6) is already 673,596 and M(1701,2) is 13,012,650. "
         f"In Dataset B, each haplotype has at most {max(x['deg1'].values())} observed one-substitution neighbours "
         f"among {x['possible1']} possible ones. The single-substitution graph built in Unit V is therefore an "
         "extremely sparse sample of the Hamming graph. This describes the sample, not the virus population. Of the 6 "
         "shortest orderings from syn01 to syn04, only one passes entirely through observed sequences."),
        ("Limitations", "Counts are combinatorial possibilities, not probabilities: they ignore transition/transversion "
         "bias, codon structure and viability. Substitutions only. Timings are machine-dependent observations. The "
         "finite checks for n ≤ 8 are evidence; the general statement rests on the induction proof."),
    ]
    out, summary = write(2, "Counting versus enumeration", "python -m experiments.exp2_counting_enumeration",
                           sections, checks, notes_after=["\n".join(tim)])
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
