"""Experiment 3: algebraic verification of the strand-symmetry groups (Unit III)."""
from itertools import permutations

from experiments import expected_unit3 as E
from experiments._report import ROOT, Check, write
from genomath import algebraic_structures as al
from genomath.sequence_processing import load_dataset

A = ROOT / "data" / "synthetic" / "dataset_a.fasta"
B = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"


def run():
    src = lambda sec: f"{E.SOURCE} {sec}"
    checks = []
    z2, _ = al.cayley_table(al.Z2)
    v4, _ = al.cayley_table(al.V4)
    as_rows = lambda t, els: {a: [t[(a, b)] for b in els] for a in els}
    checks.append(Check("Cayley table of {I, RC}", E.CAYLEY_Z2, src("§1"), as_rows(z2, al.Z2)))
    checks.append(Check("Cayley table of {I, R, C, RC}", E.CAYLEY_V4, src("§1"), as_rows(v4, al.V4)))
    for name, els, t in (("{I, RC}", al.Z2, z2), ("{I, R, C, RC}", al.V4, v4)):
        r, _ = al.check_axioms(els, t)
        checks.append(Check(f"axioms for {name}", E.AXIOMS, src("§1"), {p: r[p]["holds"] for p in E.AXIOMS}))
    orders, _ = al.element_orders(al.V4, v4)
    checks.append(Check("element orders in V4", E.ORDERS_V4, src("§1"), orders))
    checks.append(Check("V4 cyclic?", E.CYCLIC_V4, src("§1"), any(o == 4 for o in orders.values())))
    phi = {"I": (0, 0), "R": (1, 0), "C": (0, 1), "RC": (1, 1)}
    iso, _ = al.check_isomorphism(al.V4, v4, phi, lambda x, y: ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2))
    to_z4 = sum(al.check_isomorphism(al.V4, v4, dict(zip(al.V4, img)), lambda x, y: (x + y) % 4)[0]
                for img in permutations(range(4)))
    checks.append(Check("isomorphic to Z2×Z2; homomorphic bijections onto Z4",
                        (E.ISO_Z2xZ2, E.HOMOMORPHIC_BIJECTIONS_TO_Z4), src("§1"), (iso, to_z4)))
    checks.append(Check("subgroups of V4", E.SUBGROUPS_V4, src("§2"), al.subgroups(al.V4, v4)[0]))
    checks.append(Check("left cosets of {I,RC} in V4", E.COSETS_Z2_IN_V4, src("§2"), al.left_cosets(al.V4, v4, al.Z2)[0]))
    checks.append(Check("which maps coincide on Σ^n, n = 0,1,2", {str(k): v for k, v in E.FAITHFUL.items()}, src("§1"),
                        {str(n): al.distinct_on(al.V4, n)[0] for n in E.FAITHFUL}))
    checks.append(Check("|Fix(g)| by formula, n = 0..8", E.FIX, src("§3"),
                        {g: [al.fixed_count(g, n)[0] for n in range(9)] for g in al.V4}))
    checks.append(Check("|Fix(g)| by enumeration, n = 0..8", E.FIX, src("§3"),
                        {g: [al.fixed_count_enumerated(g, n)[0] for n in range(9)] for g in al.V4}))
    checks.append(Check("{I,RC} orbits by Burnside, n = 0..8", E.ORBITS_Z2, src("§4"),
                        [al.burnside(al.Z2, n)[0] for n in range(9)]))
    checks.append(Check("V4 orbits by Burnside, n = 0..8", E.ORBITS_V4, src("§4"),
                        [al.burnside(al.V4, n)[0] for n in range(9)]))
    checks.append(Check("{I,RC} orbits by enumeration, n = 0..8", E.ORBITS_Z2, src("§4"),
                        [al.orbits_enumerated(al.Z2, n)[0]["count"] for n in range(9)]))
    checks.append(Check("V4 orbits by enumeration, n = 0..8", E.ORBITS_V4, src("§4"),
                        [al.orbits_enumerated(al.V4, n)[0]["count"] for n in range(9)]))
    checks.append(Check("orbit counts at n = 12 (Burnside)", E.ORBITS_N12, src("§4"),
                        {"Z2": al.burnside(al.Z2, 12)[0], "V4": al.burnside(al.V4, 12)[0]}))
    sizes = {n: al.orbits_enumerated(al.V4, n)[0] for n in E.V4_SIZES}
    checks.append(Check("V4 orbit-size distribution, n = 1..4", E.V4_SIZES, src("§4"),
                        {n: r["sizes"] for n, r in sizes.items()}))
    checks.append(Check("orbit–stabiliser holds for every orbit, n = 1..4", True, src("§4"),
                        all(r["orbit_stabiliser"] for r in sizes.values()), consistency=True))
    checks.append(Check("complement as a permutation of Σ", E.COMPLEMENT_PERM, src("§5"),
                        al.permutation_cycles({"A": "T", "T": "A", "C": "G", "G": "C"})[0]))
    rev = al.permutation_cycles({i: 13 - i for i in range(1, 13)})[0]
    checks.append(Check("reversal of positions 1..12", E.REVERSAL_12, src("§5"),
                        {"transpositions": len(rev["cycles"]), "order": rev["order"], "parity": rev["parity"]}))
    mono = {n: al.substitution_monoid(n)[0] for n in E.MONOID}
    checks.append(Check("substitution monoid: size and units, n = 1..3", {n: list(v) for n, v in E.MONOID.items()},
                        src("§6"), {n: [m["size"], m["units"]] for n, m in mono.items()}))
    checks.append(Check("substitution monoid is a group / is commutative", (E.MONOID_GROUP, E.MONOID_COMMUTATIVE),
                        src("§6"), (mono[2]["report"]["group"], mono[2]["report"]["commutative"]["holds"])))
    checks.append(Check("non-injectivity witness (n = 2)", None, "", mono[2]["witness"], actual_only=True))
    checks.append(Check("concatenation laws, words of length ≤ 3", E.CONCAT, src("§7"), al.concatenation_laws(3)[0]))

    a, _ = load_dataset(A, reference_id="syn01")
    s = a.seqs()
    classes, _ = al.dataset_orbits(s, al.Z2)
    checks.append(Check("A: {I,RC} classes among 9 haplotypes, and the pair", (E.A_Z2_CLASSES, E.A_RC_PAIR), src("§8"),
                        (len(classes), next(c for c in classes if len(c) > 1))))
    o1, _ = al.orbit(s["syn01"], al.V4)
    checks.append(Check("A: V4 orbit of syn01", E.A_SYN01_V4, src("§8"),
                        {"size": len(o1["orbit"]), "stabiliser": o1["stabiliser"]}))
    checks.append(Check("A: V4 orbit sizes of syn03, syn09", E.A_V4_SIZE_SYN03_SYN09, src("§8"),
                        [len(al.orbit(s[h], al.V4)[0]["orbit"]) for h in ("syn03", "syn09")]))
    checks.append(Check("A: canonical representative of syn03 and syn08", [E.A_CANONICAL_SYN03] * 2, src("§8"),
                        [al.canonical(s[h])[0] for h in ("syn03", "syn08")]))
    checks.append(Check("A: R-fixed haplotypes", E.A_R_FIXED, src("§8"),
                        [h for h, x in s.items() if al.TRANSFORMS["R"](x) == x]))

    b, _ = load_dataset(B)
    bs = b.seqs()
    values = set(bs.values())
    checks.append(Check("B: RC pairs, RC-fixed, {I,RC} classes", E.B, src("§8"), {
        "RC_pairs": sum(al.TRANSFORMS["RC"](x) in values and al.TRANSFORMS["RC"](x) != x for x in values) // 2,
        "RC_fixed": sum(al.TRANSFORMS["RC"](x) == x for x in values),
        "Z2_classes": len(al.dataset_orbits(bs, al.Z2)[0])}))
    return checks


def main():
    checks = run()
    sections = [
        ("Question", "Do the strand-orientation maps {I, RC} and {I, R, C, RC} satisfy the group axioms when verified "
         "exhaustively, do Burnside's orbit counts agree with direct enumeration, and why do substitution operators "
         "fail to form a group?"),
        ("Input", "All words over {A,C,G,T} of length 0..6 (function equality) and of each length 0..8 (fixed points, "
         "orbits); Dataset A haplotypes; Dataset B haplotypes."),
        ("Mathematical model", "Maps I, R, C, RC on Σ^n under composition (a∘b)(x) = a(b(x)); group axioms; Burnside "
         "#orbits = (1/|G|)·Σ|Fix(g)|; orbit–stabiliser |Gx|·|G_x| = |G|; the monoid of substitution operators "
         "f_p ('set position i to b')."),
        ("Procedure", "Build each Cayley-table entry by composing the maps and identifying the result with an element, "
         "comparing on all 5,461 words of length ≤ 6. Check closure, associativity (all |G|³ triples), identity, "
         "inverses and commutativity from the table. Count fixed points by formula and by enumeration; count orbits "
         "by Burnside and by enumeration; generate the substitution monoid by closure under composition and test the "
         "same axioms. Apply the groups to Datasets A and B."),
        ("Interpretation",
         "{I, RC} is a group of order 2 and {I, R, C, RC} is the Klein four-group: abelian, every element self-inverse, "
         "not cyclic, isomorphic to Z2×Z2 and to no Z4 labelling. Burnside's formula matches enumeration at every length "
         "tested, including n = 0 and n = 1, where the four maps are not distinct as functions (the action is not "
         "faithful there, and Burnside still holds). Substitution operators close under composition into a monoid "
         "with 5^n elements whose only invertible element is the identity, so the statement 'mutations do not form a "
         "group' is backed by a computed witness, not just asserted. On Dataset A the RC pair syn03/syn08 is merged into "
         "one class. On Dataset B no sequence's reverse complement is present, as predicted from the start/stop codon "
         "argument, so Unit III does no analytical work on full annotated coding sequences."),
        ("Limitations", "Function equality is checked exhaustively only up to length 6; general lengths rest on the proof "
         "in `docs/unit3_expected.md` §1. R and C on their own are not natural strand operations, so the Klein group is "
         "a mathematical extension (pedagogical). RC equivalence concerns strand orientation, not mutation. "
         "Canonical representatives here are lexicographic; Unit IV restates them as min(enc(S), enc(RC(S)))."),
    ]
    out, summary = write(3, "Algebraic verification", "python -m experiments.exp3_algebraic_verification",
                           sections, checks)
    print(f"wrote {out.relative_to(ROOT)}: {summary}")


if __name__ == "__main__":
    main()
