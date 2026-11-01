"""Experiment 4: encoding and modular-indexing analysis (Unit IV)."""
from itertools import combinations

from experiments import expected_unit4 as E
from experiments._report import ROOT, Check, write
from genomath import number_theory as nt
from genomath._core import mismatch_count, mismatch_positions
from genomath.sequence_processing import load_dataset

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"


def run():
    src = lambda sec: f"{E.SOURCE} {sec}"
    checks, obs = [], {}
    a, _ = load_dataset(A, reference_id="syn01")
    s = a.seqs()
    codes = {h: nt.encode(x)[0] for h, x in s.items()}
    checks.append(Check("A: enc(S)", E.ENC_A, E.SYN, codes))
    checks.append(Check("A: enc(comp S) = 4^12 − 1 − enc(S)", E.COMP_A, E.SYN,
                        {h: nt.complement_code(x)[0] for h, x in s.items()}))
    checks.append(Check("A: canonical min(enc S, enc RC S)", E.CANON_A, E.SYN,
                        {h: nt.canonical_code(x)[0] for h, x in s.items()}))
    checks.append(Check("A: decode(enc(S)) = S", True, src("§1"), all(nt.decode(c, 12)[0] == s[h] for h, c in codes.items()), consistency=True))
    reports = {m: nt.bucket_report(codes, m)[0] for m in E.BUCKETS_A}
    as_list = lambda r: (lambda inv: [inv[h] for h in s])({h: b for b, hs in r["buckets"].items() for h in hs})
    checks.append(Check("A: buckets for m = 16, 7, 13, 17", E.BUCKETS_A, E.SYN, {m: as_list(r) for m, r in reports.items()}))
    checks.append(Check("A: colliding groups", E.GROUPS_A, E.SYN, {m: r["colliding_groups"] for m, r in reports.items()}))
    checks.append(Check("A: d_H = 1 pairs sharing a bucket", E.D1_COLLISIONS_A, E.SYN,
                        {m: len(nt.substitution_collisions(E.D1_PAIRS_A, s, m)[0]) for m in E.D1_COLLISIONS_A}))
    checks.append(Check("A: pigeonhole, 9 keys in 7 buckets", E.PIGEON_A_M7, src("§7"),
                        {"bound": reports[7]["pigeonhole_bound"], "max_load": reports[7]["max_load"]}))

    checks.append(Check("last bases seen: smallest j with m | 4^j", E.LAST_BASES, src("§3"),
                        {m: nt.last_bases_seen(m)[0] for m in E.LAST_BASES}))
    checks.append(Check("t_min for a bucket-preserving single substitution", E.T_MIN, src("§3"),
                        {m: nt.substitution_collision_t_min(m)[0] for m in E.T_MIN}))

    syn01 = s["syn01"]
    checks.append(Check("syn01 3-mer codes by sliding window", E.SYN01_3MERS, src("§1"), nt.kmer_codes(syn01, 3)[0]))
    checks.append(Check("syn01 k-mer census", E.SYN01_CENSUS, src("§2"),
                        {k: {f: nt.kmer_census([syn01], k)[0][f] for f in ("distinct", "canonical", "palindromic")}
                         for k in (3, 4)}))

    checks.append(Check("extended Euclid on (252, 198)", E.GCD_252_198, src("§5"), nt.extended_gcd(252, 198)[0]))
    inv = {}
    for label, (x, m) in {"11 mod 7": (11, 7), "7 mod 11": (7, 11), "4 mod 7": (4, 7), "4 mod 16": (4, 16)}.items():
        try:
            inv[label] = nt.mod_inverse(x, m)[0]
        except ValueError:
            inv[label] = "none"
    checks.append(Check("modular inverses", E.INVERSES, src("§5"), inv))
    checks.append(Check("linear congruences", E.CONGRUENCES, src("§5"),
                        {"3x≡6 (mod 9)": nt.solve_linear_congruence(3, 6, 9)[0],
                         "2x≡1 (mod 4)": nt.solve_linear_congruence(2, 1, 4)[0]}))
    crt = {"GAT from (0 mod 7, 2 mod 11)": nt.crt([0, 2], [7, 11])[0]}
    for mods in [(7, 11), (5, 13), (3, 5)]:
        crt[f"({mods[0]},{mods[1]}) round-trips"] = sum(
            nt.decode(nt.crt([x % m for m in mods], list(mods))[0] % 64, 3)[0] == nt.decode(x, 3)[0]
            and nt.crt([x % m for m in mods], list(mods))[0] == x for x in range(64))
    try:
        nt.crt([1, 1], [4, 6])
        crt["(4,6)"] = "accepted"
    except ValueError:
        crt["(4,6)"] = "refused"
    checks.append(Check("CRT reconstruction of 3-mer codes", E.CRT, src("§6"), crt))
    checks.append(Check("ord_m(4) (pedagogical)", E.ORDERS, src("§4"),
                        {m: nt.multiplicative_order(4, m)[0] for m in E.ORDERS}))
    checks.append(Check("φ(m)", E.PHI, src("§4"), {m: nt.euler_phi(m)[0] for m in E.PHI}))
    checks.append(Check("Euler: 4^φ(m) ≡ 1 (mod m)", E.EULER_4, src("§4"),
                        {m: nt.mod_pow(4, nt.euler_phi(m)[0], m)[0] == 1 for m in E.EULER_4}))
    other = syn01[:5] + "A" + syn01[6:11] + "G"
    oc = nt.encode(other)[0]
    checks.append(Check("two-substitution cancellation (Fermat periodicity)", E.CANCEL, src("§4"),
                        {"string": other, "d_H": mismatch_count(syn01, other),
                         "collide": [m for m in E.MODULI if oc % m == codes["syn01"] % m]}))

    b, _ = load_dataset(B)
    bs = b.seqs()
    n = b.n
    bcodes = {h: nt.encode(x)[0] for h, x in bs.items()}
    breps = {m: nt.bucket_report(bcodes, m)[0] for m in E.MODULI}
    checks.append(Check("B: pigeonhole bound ⌈30/m⌉", E.PIGEON_B, src("§7"),
                        {m: breps[m]["pigeonhole_bound"] for m in E.PIGEON_B}))
    checks.append(Check("B: max load ≥ bound for every m", True, src("§7"),
                        all(r["max_load"] >= r["pigeonhole_bound"] for r in breps.values()), consistency=True))
    obs["b_buckets"] = {m: (r["occupied"], r["max_load"], r["collisions"]) for m, r in breps.items()}
    checks.append(Check("B: (occupied, max load, collisions) per m", None, "",
                        obs["b_buckets"], actual_only=True))
    suffixes = {}
    for x in bs.values():
        suffixes[x[-3:]] = suffixes.get(x[-3:], 0) + 1
    obs["suffixes"] = suffixes
    checks.append(Check("B: last 3 bases of every haplotype (count)", None, "", suffixes, actual_only=True))

    d1 = [(x, y) for x, y in combinations(bs, 2) if mismatch_count(bs[x], bs[y]) == 1]
    predicted, observed = {}, {}
    for m in E.MODULI:
        t_min = nt.substitution_collision_t_min(m)[0]
        predicted[m] = sum(1 for x, y in d1 if t_min is not None and n - mismatch_positions(bs[x], bs[y])[0] >= t_min)
        observed[m] = len(nt.substitution_collisions(d1, bs, m)[0])
    checks.append(Check("B: d_H = 1 pairs", 8, "docs/dataset_b_qc.md (independent script)", len(d1)))
    checks.append(Check("B: d_H = 1 collisions, observed = predicted by §3 theorem", predicted, src("§3"), observed, consistency=True))
    checks.append(Check("B: d_H = 1 collisions for odd primes (theorem: none)", {m: 0 for m in (7, 13, 17, 31, 61)},
                        src("§3"), {m: observed[m] for m in (7, 13, 17, 31, 61)}))

    census = {k: nt.kmer_census(list(bs.values()), k)[0] for k in (4, 6, 8)}
    checks.append(Check("B: canonical = distinct − present RC pairs (k = 4, 6, 8)", {k: True for k in census}, src("§2"),
                        {k: c["identity_holds"] for k, c in census.items()}, consistency=True))
    checks.append(Check("B: whole-space canonical classes (Burnside)", E.SPACE_CANONICAL, src("§2"),
                        {k: c["space_canonical"] for k, c in census.items()}))
    obs["census"] = {k: {f: c[f] for f in ("distinct", "canonical", "palindromic", "rc_pairs_present")}
                     for k, c in census.items()}
    checks.append(Check("B: k-mer census (distinct, canonical, palindromic, RC pairs present)", None, "",
                        obs["census"], actual_only=True))

    k8 = {}
    for x in bs.values():
        for c in nt.kmer_codes(x, 8)[0]:
            k8[c] = c
    kreps = {m: nt.bucket_report(k8, m)[0] for m in (64, 61)}
    checks.append(Check("B: distinct 8-mers bucketed, max load ≥ ⌈N/m⌉", {64: True, 61: True}, src("§7"),
                        {m: r["max_load"] >= r["pigeonhole_bound"] for m, r in kreps.items()}, consistency=True))
    obs["k8"] = {m: {"keys": len(k8), "occupied": r["occupied"], "max_load": r["max_load"],
                     "bound": r["pigeonhole_bound"]} for m, r in kreps.items()}
    checks.append(Check("B: 8-mer buckets, m = 64 (last 3 bases) vs m = 61", None, "", obs["k8"], actual_only=True))
    return checks, obs


def main():
    checks, obs = run()
    c, k8, suf, bb = obs["census"], obs["k8"], obs["suffixes"], obs["b_buckets"]
    one_suffix = len(suf) == 1
    suffix_text = (f"On Dataset B every haplotype ends in {next(iter(suf))}, so m = 16 and m = 64 put all 30 haplotypes "
                   f"into {bb[16][0]} and {bb[64][0]} bucket(s): the stop codon is all they can see. "
                   if one_suffix else f"Dataset B haplotypes end in {len(suf)} different 3-base suffixes. ")
    sections = [
        ("Question", "What does representing DNA words as base-4 integers make provable? In particular, how does the "
         "choice of modulus change bucket collisions, when can a single substitution leave the bucket unchanged, and "
         "how many {I,RC}-canonical k-mers does real data use?"),
        ("Input", "Dataset A haplotypes (n = 12); all 3-mer codes 0..63; Dataset B haplotypes (n = 1701) and their "
         f"k-mers for k = 4, 6, 8; moduli {E.MODULI} (provisional decision 8)."),
        ("Mathematical model", "enc(S) = Σ s_i 4^(n−i) (bijection; decoding by the division algorithm); "
         "enc(comp S) = 4^n − 1 − enc(S); canonical = min(enc S, enc RC S); b(x) = x mod m; pigeonhole ⌈N/m⌉; exact "
         "single-substitution criterion via m = 2^a·q; Euclid / extended Euclid / inverses / linear congruences; "
         "CRT; Euler φ and multiplicative order (pedagogical)."),
        ("Procedure", "Encode every Dataset A haplotype and compare codes, complement codes, canonical codes and "
         "buckets with the hand tables. Test every one-substitution pair against the collision theorem on both datasets. "
         "Reconstruct all 64 3-mer codes by CRT under several moduli pairs. Compute orders and φ, and the two-substitution "
         "cancellation example. Bucket Dataset B full-sequence codes and distinct 8-mers. Take the {I,RC} k-mer census "
         "of Dataset B."),
        ("Interpretation",
         "Encoding is exact and invertible, and the complement identity holds. m = 16 and m = 64 see only the last 2 or 3 "
         "bases, so they collide on every substitution further from the end. " + suffix_text +
         "Every odd modulus ≥ 5 separates all "
         "single substitutions, as proved, while two substitutions can still cancel: an opposite pair (+δ, −δ) cancels "
         "exactly when ord_m(4) divides the distance between the positions (mod 7 and 13 in the example). Bucket-sharing is therefore a property of the residue, not of "
         "biological similarity. On Dataset B the {I,RC} quotient does real work at the k-mer level: for k = 8, "
         f"{c[8]['distinct']} distinct 8-mers collapse to {c[8]['canonical']} canonical classes, and for k = 4 all "
         f"{c[4]['canonical']} canonical classes of the whole space occur. For 8-mers, m = 64 occupies "
         f"{k8[64]['occupied']} of 64 buckets with max load {k8[64]['max_load']}, while m = 61 occupies "
         f"{k8[61]['occupied']} with max load {k8[61]['max_load']} (observation)."),
        ("Limitations", "Modular indexing here is a representation experiment, not an established clinical or genomic "
         "technique. The encoding A=0..T=3 is a convention; any bijection would do. CRT, Euler and Fermat are "
         "pedagogical. The k-mer census describes 30 haplotypes from one season and implies nothing about other data. "
         "Canonical k-mer counting is documented practice for reads of unknown orientation (Jellyfish manual v0.9, "
         "`math_audit.md`); the CDS here have known orientation, so the census is a mathematical demonstration."),
    ]
    out, summary = write(4, "Encoding and modular indexing", "python -m experiments.exp4_modular_indexing",
                           sections, checks)
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
