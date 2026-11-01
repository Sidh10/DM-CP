"""The seven dashboard sections of spec §8. Each result is shown with its trace and a category tag."""
import streamlit as st

from app import presentation as pr
from app import ui
from genomath import algebraic_structures as al
from genomath import counting_recursion as cr
from genomath import mutation_graph as mg
from genomath import number_theory as nt
from genomath import phylogenetic_analysis as pa
from genomath import set_logic as sl
from genomath import tree_analysis as ta
from genomath._core import mismatch_count


def dataset_explorer(ctx):
    data, seqs = ctx["data"], ctx["seqs"]
    ui.section_header(1, "Dataset explorer", "What exactly was loaded, and which preprocessing decisions were made?")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("records read", ctx["records_read"])
    c2.metric("kept", len(data.records))
    c3.metric("haplotypes", len(data.haplotypes))
    c4.metric("length n", data.n)

    ui.result("observation", f"Alphabet Σ = {{A, C, G, T}}; every kept sequence is an element of Σ^{data.n}, "
              f"so the sequence space has 4^{data.n} elements.")
    if data.dropped:
        ui.tag("algorithm", "Records dropped by the validation rule, with the reason logged:")
        ui.table(["id", "reason"], [[d["id"], d["reason"]] for d in pr.describe_dropped(data.dropped)])
        ui.tag("limitation", "Dropped records are excluded from every later calculation. Ambiguity codes are not "
               "resolved or imputed, because a single position would have to be guessed.")
    else:
        ui.result("observation", "No record contained a character outside {A, C, G, T}.")

    groups = [h for h in data.haplotypes if len(h.members) > 1]
    if groups:
        ui.tag("algorithm", "Identical sequences collapsed into haplotypes:")
        ui.table(["haplotype", "members"], [[h.id, ", ".join(h.members)] for h in groups])
    ui.result("calculation", "Equality of strings (d_H = 0) is reflexive, symmetric and transitive, so it is an "
              "equivalence relation; the haplotypes are its equivalence classes.")
    ui.result("algorithm", f"Reference sequence R = {data.reference.id}", ctx["reference_note"])
    ui.show_trace(ctx["load_trace"], "Show the full preprocessing trace")

    rows = [[h.id, len(h.members), mismatch_count(data.reference.seq, h.seq),
             pr.truncate(h.seq, 60)] for h in data.haplotypes]
    ui.table(["haplotype", "records", "d_H to R", "sequence"], rows,
             "d_H to R is the Hamming distance to the reference; it equals |M(R,S)|.")


def set_theory(ctx):
    data, seqs = ctx["data"], ctx["seqs"]
    ids = list(seqs)
    R = data.reference.seq
    ui.section_header(2, "Set theory and logic (Unit I)",
                      "Where do two sequences differ, and what do set operations on those positions tell us?")
    a = st.selectbox("sequence S", ids, index=min(1, len(ids) - 1), key="sl_a")
    b = st.selectbox("sequence T", ids, index=min(2, len(ids) - 1), key="sl_b")

    st.subheader("Mutation sets relative to the reference")
    for label, sid in (("S", a), ("T", b)):
        m, t1 = sl.mutation_set(R, seqs[sid])
        ms, t2 = sl.allele_mutation_set(R, seqs[sid])
        ui.result("calculation", f"{label} = {sid}: |M(R,{label})| = {len(m)}",
                  f"M  = {sorted(m)}\nM* = {sorted(ms)}")
        ui.show_trace(list(t1) + list(t2), f"Show how M(R,{sid}) was built")
    ui.result("calculation", "M holds positions; M* holds (position, base) pairs. Two sequences can share a position "
              "without sharing a mutation, which is why both exist.")

    st.subheader("Set operations")
    ma = sl.mutation_set(R, seqs[a])[0]
    mb = sl.mutation_set(R, seqs[b])[0]
    ops, t = sl.set_operations(ma, mb, (f"M({a})", f"M({b})"))
    ui.table(["operation", "result", "size"],
             [["intersection ∩", str(sorted(ops["intersection"])), len(ops["intersection"])],
              ["union ∪", str(sorted(ops["union"])), len(ops["union"])],
              [f"{a} only \\", str(sorted(ops["a_minus_b"])), len(ops["a_minus_b"])],
              [f"{b} only \\", str(sorted(ops["b_minus_a"])), len(ops["b_minus_a"])],
              ["symmetric difference Δ", str(sorted(ops["symmetric_difference"])), len(ops["symmetric_difference"])]])
    ui.show_trace(t, "Show where each position lands")
    msa = sl.allele_mutation_set(R, seqs[a])[0]
    msb = sl.allele_mutation_set(R, seqs[b])[0]
    shared_alleles = sorted(msa & msb)
    ui.result("calculation", f"Shared positions: {sorted(ops['intersection'])}. Shared *mutations* (same position and "
              f"same base): {shared_alleles or 'none'}.")

    st.subheader("Hamming distance from the mutation sets")
    dec, t = sl.distance_decomposition(R, seqs[a], seqs[b])
    ui.result("calculation", f"d_H({a},{b}) = |M Δ M'| + |shared positions with different bases| = "
              f"{len(dec['sym_diff'])} + {len(dec['shared_differing'])} = {dec['rhs']}; computed directly: {dec['d_H']}")
    ui.show_trace(t, "Show the identity term by term")

    st.subheader("The threshold relation R_k")
    k = st.slider("threshold k", 0, max(1, ctx["max_distance"]), 1, key="sl_k")
    pairs, _ = sl.threshold_relation(seqs, k)
    report, t = sl.relation_properties(ids, pairs)
    ui.table(["property", "holds", "counterexample"],
             [[p, "yes" if report[p]["holds"] else "no", str(report[p]["witness"] or "")]
              for p in ("reflexive", "symmetric", "antisymmetric", "transitive")])
    if report["equivalence"]:
        ui.result("calculation", f"R_{k} is an equivalence relation here.")
    else:
        w = report["transitive"]["witness"]
        ui.result("calculation", f"R_{k} is **not** an equivalence relation: it is reflexive and symmetric but not "
                  f"transitive." + (f" Witness: d({w[0]},{w[1]}) ≤ {k} and d({w[1]},{w[2]}) ≤ {k}, "
                                    f"but d({w[0]},{w[2]}) = {mismatch_count(seqs[w[0]], seqs[w[2]])} > {k}." if w else ""))
    ui.show_trace(t, "Show each property test")
    closure, tc = sl.transitive_closure(ids, pairs)
    classes, _ = sl.equivalence_classes(ids, closure)
    ui.result("algorithm", f"Transitive closure of R_{k} (Warshall) has {len(classes)} equivalence classes. "
              "These are exactly the connected components of the threshold graph in section 6.",
              "\n".join(", ".join(c) for c in classes))
    ui.show_trace(tc, "Show the Warshall pivots")

    st.subheader("Predicates and filtering")
    pk = st.slider("P(S): d_H(R,S) ≤ k", 0, max(1, ctx["max_distance"]), min(2, ctx["max_distance"]), key="sl_pk")
    pos = st.number_input("Q(S): position i", 1, data.n, 1, key="sl_pos")
    base = st.selectbox("Q(S): base", list("ACGT"), key="sl_base")
    combo = st.radio("Boolean combination", ["P ∧ Q", "P ∨ Q", "P ∧ ¬Q", "¬P"], horizontal=True, key="sl_combo")
    P, Q = sl.within_distance(R, pk), sl.base_at(int(pos), base)
    predicate = {"P ∧ Q": P & Q, "P ∨ Q": P | Q, "P ∧ ¬Q": P & ~Q, "¬P": ~P}[combo]
    kept, t = sl.filter_sequences(seqs, predicate, (P, Q))
    ui.result("algorithm", f"{predicate.name} is true for {len(kept)} of {len(seqs)} haplotypes",
              ", ".join(kept) or "none")
    ui.show_trace(t, "Show the truth table")

    st.subheader("Equivalence: identity and strand orientation")
    rec = {r.id: r.seq for r in data.records}
    icls, _ = sl.equivalence_classes(list(rec), sl.identity_relation(rec)[0])
    ocls, ot = sl.equivalence_classes(list(rec), sl.orientation_relation(rec)[0])
    ui.result("calculation", f"Identity gives {len(icls)} classes (the haplotypes). "
              f"Orientation (S ~ RC(S)) gives {len(ocls)} classes.")
    merged = [c for c in ocls if len(c) > len(max(icls, key=len)) or (len(c) > 1 and c not in icls)]
    if merged:
        ui.result("interpretation", "Classes that orientation merges beyond identity hold sequences that describe the "
                  "same double-stranded molecule read from opposite strands.", "\n".join(", ".join(c) for c in merged))
    else:
        ui.result("observation", "Orientation merges nothing beyond identity in this dataset: no sequence's reverse "
                  "complement is present.")
    ui.show_trace(ot, "Show the orientation relation test")


def counting(ctx):
    data = ctx["data"]
    ui.section_header(3, "Counting and recursion (Unit II)",
                      "How large is the space of sequences and of k-substitution mutants, and when can we enumerate it?")
    space, t = cr.sequence_space(data.n)
    digits = len(str(space))
    ui.result("calculation", f"N(n) = 4^n for n = {data.n}: a {digits}-digit number"
              + (f" = {space:,}" if digits <= 12 else ""), None if digits <= 12 else f"{str(space)[:40]}… ({digits} digits)")
    ui.show_trace(t)

    st.subheader("M(n,k) = C(n,k)·3^k")
    col1, col2 = st.columns(2)
    n_small = col1.slider("n for the table", 0, 14, min(8, data.n), key="cr_n")
    k_pick = col2.slider("k", 0, n_small, min(2, n_small), key="cr_k")
    value, t = cr.mutation_count(n_small, k_pick)
    ui.result("calculation", f"M({n_small},{k_pick}) = C({n_small},{k_pick})·3^{k_pick} = {value:,}")
    ui.show_trace(t)
    table, tt = cr.recurrence_table(n_small)
    header = ["n"] + [str(i) for i in range(n_small + 1)] + ["row sum", "4^n"]
    body = [[str(i)] + [f"{x:,}" for x in row] + [""] * (n_small + 1 - len(row)) + [f"{sum(row):,}", f"{4 ** i:,}"]
            for i, row in enumerate(table)]
    ui.table(header, body, "Built bottom-up from M(n,k) = M(n−1,k) + 3·M(n−1,k−1) with M(0,0) = 1. "
             "Each row sums to 4^n, by the binomial theorem.")
    ver, tv = cr.verify_recurrence(n_small)
    ui.result("algorithm", f"Closed form vs recurrence over all cells up to n = {n_small}: "
              f"{len(ver['mismatches'])} mismatches; rows summing wrongly: {len(ver['bad_rows'])}")
    ui.show_trace(list(tt) + list(tv))
    ui.result("limitation", "This is a finite check, not a proof. The induction proof that C(n,k)·3^k satisfies the "
              "recurrence (via Pascal's rule) is in docs/unit2_expected.md §2.")

    st.subheader("Hamming ball and the real dataset")
    r = st.slider("radius r", 0, 3, 1, key="cr_r")
    ball, t = cr.ball_size(data.n, r)
    ui.result("calculation", f"|B(n,{r})| for n = {data.n}: {ball:,} sequences within {r} substitution(s), "
              f"so {ball - 1:,} possible neighbours.")
    ui.show_trace(t)
    deg = {h: sum(1 for o in ctx["seqs"] if o != h and mismatch_count(ctx['seqs'][h], ctx['seqs'][o]) <= r)
           for h in ctx["seqs"]}
    best = max(deg, key=deg.get)
    ui.result("observation", f"The best-connected haplotype ({best}) has {deg[best]} observed neighbours within "
              f"distance {r}, out of {ball - 1:,} possible.")
    ui.result("interpretation", "Observed data are an extremely sparse sample of the sequence space. This describes "
              "the sample, not the population it came from.")

    st.subheader("Recursive generation, with a guard")
    gn = st.slider("reference length for generation", 1, 12, min(6, data.n), key="cr_gn")
    gk = st.slider("exactly k substitutions", 0, gn, 1, key="cr_gk")
    ref = ctx["seqs"][list(ctx["seqs"])[0]][:gn]
    try:
        out, t = cr.generate_mutants(ref, gk)
        ui.result("algorithm", f"Generated {len(out):,} distinct sequences at distance exactly {gk} from {ref}; "
                  f"the formula predicts {cr.mutation_count(gn, gk)[0]:,}.",
                  ", ".join(out[:12]) + (" …" if len(out) > 12 else ""))
        ui.show_trace(t)
    except cr.GenerationLimitError as e:
        ui.result("limitation", f"Refused: {e}")
        ui.show_trace(e.trace, expanded=True)

    st.subheader("Inclusion–exclusion and the pigeonhole principle")
    ie_n = st.slider("n for 'uses all four bases'", 1, 12, 4, key="cr_ie")
    ie, t = cr.inclusion_exclusion_all_bases(ie_n)
    ui.result("calculation", f"Strings of length {ie_n} using all four bases: 4^n − 4·3^n + 6·2^n − 4 = {ie:,}")
    ui.show_trace(t)
    kmax, t = cr.max_forced_k(data.n)
    ui.result("calculation", f"A sequence of length {data.n} has L − k + 1 windows of length k and there are 4^k "
              f"possible k-mers, so a repeat is forced for every k ≤ {kmax}.")
    ui.show_trace(t)
    sid = st.selectbox("find a forced repeat in", list(ctx["seqs"]), key="cr_seq")
    witness, t = cr.find_repeated_kmer(ctx["seqs"][sid], kmax)
    if witness:
        ui.result("algorithm", f"In {sid}, the {kmax}-mer {witness[0]} occurs at positions {witness[1]} and {witness[2]}.")
    ui.show_trace(t)


def algebra(ctx):
    data, seqs = ctx["data"], ctx["seqs"]
    ui.section_header(4, "Algebraic structures (Unit III)",
                      "Which transformations of DNA genuinely satisfy the group axioms, and what do they not model?")
    sid = st.selectbox("sequence", list(seqs), key="al_seq")
    s = seqs[sid]
    cols = st.columns(4)
    for col, name in zip(cols, al.V4):
        out, _ = al.transform(name, s)
        col.markdown(f"**{name}**")
        col.code(pr.truncate(out, 24), language=None)
    rc2, t = al.transform("RC", al.transform("RC", s)[0])
    ui.result("calculation", f"RC(RC(S)) = S: {rc2 == s}. RC is an involution, so it is its own inverse.")

    st.subheader("Cayley tables and the group axioms")
    which = st.radio("group", ["{I, RC} — strand orientation", "{I, R, C, RC} — Klein four-group",
                               "{id, (A T)(C G), (A G)(C T), (A C)(G T)} — on bases (D9)"], key="al_group")
    if which.startswith("{I, RC}"):
        elements, (table, t) = al.Z2, al.cayley_table(al.Z2)
        identity = "I"
    elif which.startswith("{I, R,"):
        elements, (table, t) = al.V4, al.cayley_table(al.V4)
        identity = "I"
    else:
        elements, (table, t) = al.BASE_K4, al.base_cayley_table()
        identity = "id"
    header, body = pr.cayley_table(elements, table)
    ui.table(header, body, "Entry in row a, column b is a∘b (apply b first, then a).")
    ui.show_trace(t, "Show how each entry was identified")
    report, tr = al.check_axioms(elements, table)
    ui.table(["axiom", "holds", "detail"],
             [["closure", "yes" if report["closure"]["holds"] else "no", str(report["closure"]["witness"] or "")],
              ["associativity", "yes" if report["associativity"]["holds"] else "no",
               f"all {report['associativity']['triples']} triples checked"],
              ["identity", "yes" if report["identity"]["holds"] else "no", str(report["identity"]["element"])],
              ["inverses", "yes" if report["inverses"]["holds"] else "no", str(report["inverses"]["inverse_of"])],
              ["commutative", "yes" if report["commutative"]["holds"] else "no", ""]])
    orders, to = al.element_orders(elements, table, identity)
    ui.result("calculation", f"Element orders: {orders}. Cyclic: {any(v == len(elements) for v in orders.values())}. "
              + ("Every non-identity element has order 2, so this is the Klein four-group, isomorphic to Z2 × Z2 and "
                 "not to Z4." if len(elements) == 4 else ""))
    ui.show_trace(list(tr) + list(to), "Show the axiom checks")
    if len(elements) == 4 and identity == "I":
        subs, ts = al.subgroups(elements, table)
        cos, tcos = al.left_cosets(elements, table, al.Z2)
        ui.result("algorithm", f"{len(subs)} subgroups; left cosets of {{I, RC}}: "
                  + ", ".join("{" + ", ".join(c) + "}" for c in cos)
                  + f"; Lagrange index |G|/|H| = {len(elements) // 2}")
        ui.show_trace(list(ts) + list(tcos))

    st.subheader("Orbits and Burnside's lemma")
    nb = st.slider("word length n", 0, 8, 2, key="al_n")
    rows = []
    for g in al.V4:
        rows.append([g, f"{al.fixed_count(g, nb)[0]:,}"])
    ui.table(["element", "|Fix(g)| on Σ^n"], rows)
    for grp, label in ((al.Z2, "{I, RC}"), (al.V4, "{I, R, C, RC}")):
        value, t = al.burnside(grp, nb)
        enum, _ = al.orbits_enumerated(grp, nb)
        ui.result("calculation", f"{label}: Burnside gives {value:,} orbits on Σ^{nb}; direct enumeration gives "
                  f"{enum['count']:,}. Orbit sizes: {enum['sizes']}")
        ui.show_trace(t)
    classes, t = al.dataset_orbits(seqs, al.Z2)
    ui.result("algorithm", f"{{I, RC}} splits the {len(seqs)} haplotypes into {len(classes)} orientation classes.",
              "\n".join(", ".join(c) for c in classes if len(c) > 1) or "every class is a singleton")
    ui.show_trace(t)

    st.subheader("Why mutations do not form a group")
    mono, t = al.substitution_monoid(2)
    ui.result("calculation", f"Operators 'set position i to base b' on Σ^2 close under composition into a structure "
              f"with {mono['size']} elements, all idempotent. Its only invertible element is the identity "
              f"{mono['units']}, so it is a monoid and not a group.")
    if mono["witness"]:
        f, x, y, img = mono["witness"]
        ui.result("algorithm", f"Witness: writing pattern {f} sends both {x} and {y} to {img}, so it is not injective "
                  "and has no inverse.")
    ui.show_trace(t)
    ui.result("limitation", "This is the precise sense in which 'DNA mutations form a group' is false. The {I, RC} "
              "group is legitimate because RC is a bijection describing strand orientation, not a mutation event.")

    st.subheader("Transitions and transversions (decision D9)")
    fpf, tf = al.fixed_point_free_involutions()
    cov, tc = al.ordered_pair_coverage()
    ui.result("calculation", f"Σ has exactly {len(fpf)} fixed-point-free involutions (3!! perfect matchings). With the "
              f"identity they form a Klein four-group acting on the four bases. Each non-identity element realises 4 "
              f"ordered base pairs, covering all {cov['total']}.")
    ui.show_trace(list(tf) + list(tc))
    c1, c2 = st.columns(2)
    x = c1.selectbox("from base", list("ACGT"), key="al_x")
    y = c2.selectbox("to base", list("ACGT"), index=2, key="al_y")
    if x == y:
        ui.result("observation", "Pick two different bases to classify a substitution.")
    else:
        cls, t = al.classify_substitution(x, y)
        ui.result("calculation", f"{x}↔{y} is realised by {cls['element']}, so it is a **{cls['kind']}**. "
                  f"Same purine/pyrimidine class: {cls['same_class']}.")
        ui.show_trace(t)
    census, t = al.substitution_census((data.reference.seq, v) for h, v in seqs.items() if h != data.reference.id)
    ui.table(["quantity", "value"],
             [["substitutions vs the reference", f"{census['substitutions']:,}"],
              ["transitions", f"{census['transitions']:,}"], ["transversions", f"{census['transversions']:,}"],
              ["ratio transitions/transversions", f"{census['ratio']:.3f}" if census["ratio"] else "n/a"],
              ["combinatorial baseline", census["baseline_ratio"]],
              ["by unordered base pair", str(census["by_base_pair"])]])
    ui.show_trace(t, expanded=False)
    ui.result("limitation", "Counts are direction-free: two sequences differ at a position without showing which base "
              "came first. The 0.5 baseline is what equal likelihood over the six unordered pairs would give; it is a "
              "statement about the alphabet only. No claim is made here about what ratio real data should show.")


def number_theory(ctx):
    data, seqs = ctx["data"], ctx["seqs"]
    ui.section_header(5, "Number theory (Unit IV)",
                      "What follows from representing fixed-length DNA words as integers?")
    sid = st.selectbox("sequence", list(seqs), key="nt_seq")
    s = seqs[sid]
    code, t = nt.encode(s)
    ui.result("calculation", f"enc({sid}) with A=0, C=1, G=2, T=3 (most significant digit first)", f"{code:,}")
    ui.show_trace(t, "Show the Horner steps")
    back, tb = nt.decode(code, data.n)
    ui.result("algorithm", f"Decoding by repeated division (the division algorithm) returns the sequence: {back == s}")
    ui.show_trace(tb, "Show the division steps")
    comp, tcmp = nt.complement_code(s)
    ui.result("calculation", f"enc(comp S) = 4^n − 1 − enc(S) = {comp:,}")
    ui.show_trace(tcmp)
    canon, tcan = nt.canonical_code(s)
    ui.result("calculation", f"canonical code min(enc S, enc RC S) = {canon:,}")
    ui.show_trace(tcan)

    st.subheader("Modular indexing: buckets and collisions")
    moduli = st.multiselect("moduli m", [8, 13, 16, 17, 31, 61, 64, 7], default=[16, 64, 7, 61], key="nt_mods")
    codes = {h: nt.encode(v)[0] for h, v in seqs.items()}
    rows = []
    for m in sorted(moduli):
        rep, _ = nt.bucket_report(codes, m)
        j, _ = nt.last_bases_seen(m)
        t_min, _ = nt.substitution_collision_t_min(m)
        rows.append([m, rep["occupied"], rep["max_load"], rep["pigeonhole_bound"], len(rep["colliding_groups"]),
                     "all" if j is None else f"last {j}", "never" if t_min is None else f"t ≥ {t_min}"])
    ui.table(["m", "buckets used", "max load", "⌈N/m⌉ bound", "colliding groups", "bases the bucket sees",
              "single substitution can collide"], rows)
    ui.result("calculation", "If m divides 4^j then x mod m depends only on the last j bases, because every higher "
              "positional weight is ≡ 0. If the odd part of m is at least 5, no single substitution can ever keep the "
              "bucket (proved in docs/unit4_expected.md §3).")
    d1 = [(x, y) for i, x in enumerate(seqs) for y in list(seqs)[i + 1:] if mismatch_count(seqs[x], seqs[y]) == 1]
    if d1 and moduli:
        coll = {m: len(nt.substitution_collisions(d1, seqs, m)[0]) for m in sorted(moduli)}
        ui.result("algorithm", f"Of the {len(d1)} observed pairs at distance 1, the number sharing a bucket: {coll}")
    ui.result("limitation", "Sharing a bucket says nothing about similarity. The residue depends on the whole integer "
              "code, not on biological closeness.")

    st.subheader("Euclid, inverses and congruences")
    c1, c2 = st.columns(2)
    x = c1.number_input("a", 0, 10 ** 9, 252, key="nt_a")
    y = c2.number_input("b", 1, 10 ** 9, 198, key="nt_b")
    (g, u, v), t = nt.extended_gcd(int(x), int(y))
    ui.result("algorithm", f"gcd({x},{y}) = {g} = {x}·{u} + {y}·{v}")
    ui.show_trace(t, "Show each division step and the Bézout coefficients")
    m = st.number_input("modulus m for the inverse of a", 2, 10 ** 9, 7, key="nt_m")
    try:
        inv, t = nt.mod_inverse(int(x), int(m))
        ui.result("calculation", f"{x}^(-1) ≡ {inv} (mod {m})")
        ui.show_trace(t)
    except ValueError as e:
        ui.result("limitation", f"No inverse: {e}")
    sols, t = nt.solve_linear_congruence(int(x), int(y), int(m))
    ui.result("algorithm", f"{x}x ≡ {y} (mod {m}) has {len(sols)} solution(s): {sols or 'none'}")
    ui.show_trace(t)

    st.subheader("Chinese Remainder Theorem on a k-mer code")
    kmer = st.text_input("3-mer", s[:3], max_chars=3, key="nt_kmer").upper()
    mods = st.multiselect("pairwise coprime moduli", [3, 5, 7, 11, 13], default=[7, 11], key="nt_crt")
    if len(set(kmer) - set("ACGT")) == 0 and len(kmer) == 3 and len(mods) >= 2:
        value = nt.encode(kmer)[0]
        product = 1
        for mm in mods:
            product *= mm
        residues = [value % mm for mm in mods]
        try:
            back, t = nt.crt(residues, list(mods))
            ui.result("calculation", f"enc({kmer}) = {value}; residues {dict(zip(mods, residues))}; "
                      f"CRT reconstructs {back} (unique modulo {product})")
            ui.show_trace(t, "Show the CRT terms")
            if product < 64:
                ui.result("limitation", f"The product of the moduli is {product} < 64, so the reconstruction is only "
                          f"correct for codes below {product}. This one is {'correct' if back == value else 'wrong'}.")
        except ValueError as e:
            ui.result("limitation", f"CRT refused: {e}")
    else:
        ui.result("observation", "Enter a 3-mer over {A,C,G,T} and pick at least two moduli.")

    st.subheader("Canonical k-mers")
    kk = st.slider("k", 2, 8, 4, key="nt_k")
    census, t = nt.kmer_census(list(seqs.values()), kk)
    ui.table(["quantity", "value"],
             [["distinct k-mers observed", f"{census['distinct']:,}"],
              ["canonical classes among them", f"{census['canonical']:,}"],
              ["RC-palindromic k-mers", census["palindromic"]],
              ["(w, RC(w)) pairs both present", census["rc_pairs_present"]],
              ["canonical classes in the whole space (Burnside)", f"{census['space_canonical']:,}"]])
    ui.show_trace(t)
    ui.result("interpretation", "Collapsing a k-mer with its reverse complement is documented practice for sequence "
              "data of unknown strand orientation (Jellyfish manual; see docs/math_audit.md). These coding sequences "
              "have known orientation, so here it is a mathematical demonstration rather than a necessity.")


def mutation_graph(ctx):
    data, seqs = ctx["data"], ctx["seqs"]
    ids = list(seqs)
    ui.section_header(6, "Mutation graph (Unit V)",
                      "Which sequences are linked by few substitutions, and how does the picture change with k?")
    mode = st.radio("edge rule", ["threshold: 1 ≤ d_H ≤ k, weight d_H", "single substitution: d_H = 1"],
                    key="mg_mode")
    single = mode.startswith("single")
    k = 1 if single else st.slider("threshold k", 1, max(1, ctx["max_distance"]), min(2, ctx["max_distance"]),
                                   key="mg_k")
    g, t = mg.build_graph(seqs, k, "single" if single else "threshold")
    ui.result("algorithm", f"{'Edge iff d_H = 1' if single else f'Edge iff 1 ≤ d_H ≤ {k}, weight d_H'}: "
              f"{len(g.edges())} edges on {len(g.vertices)} vertices")
    ui.show_trace(t)

    deg, td = mg.degrees(g)
    comps, tc = mg.connected_components(g)
    cyc, tcy = mg.cyclomatic_number(g)
    bip, tb = mg.is_bipartite(g)
    tri, _ = mg.count_triangles(g)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("components", len(comps))
    c2.metric("max degree", max(deg.values()))
    c3.metric("independent cycles", cyc)
    c4.metric("triangles", tri)
    source = st.selectbox("source vertex", ids, key="mg_src")
    target = st.selectbox("target vertex", ids, index=min(len(ids) - 1, 1), key="mg_tgt")
    pv, tp = mg.path_versus_direct(g, seqs, source, target)
    ui.plot_graph(g, k, highlight_path=pv["paths"][0] if pv["paths"] else None,
                  title=f"{'G_=1' if single else f'G_{k}'} — shortest path {source} → {target} highlighted"
                  if pv["paths"] else f"{'G_=1' if single else f'G_{k}'}")

    st.subheader("Degrees, matrices and the handshake lemma")
    ui.table(["vertex", "degree"], [[v, deg[v]] for v in ids])
    ui.result("calculation", f"Σ deg = {sum(deg.values())} = 2·|E| = {2 * len(g.edges())}")
    ui.show_trace(td)
    if len(ids) <= 14:
        with st.expander("adjacency matrix"):
            ui.matrix(ids, mg.adjacency_matrix(g)[0], blank_diagonal=False)
        with st.expander("incidence matrix"):
            (edges, inc), _ = mg.incidence_matrix(g)
            ui.table([""] + [f"{u}–{v}" for u, v, _ in edges], [[ids[i]] + row for i, row in enumerate(inc)],
                     "Each column has exactly two 1s; each row sum is that vertex's degree.")
    else:
        ui.result("limitation", f"The adjacency and incidence matrices are {len(ids)}×{len(ids)} and "
                  f"{len(ids)}×{len(g.edges())}; they are omitted here for readability and available in the modules.")

    st.subheader("Traversals, components and cycles")
    b, tbfs = mg.bfs(g, source)
    d, tdfs = mg.dfs(g, source)
    ui.result("algorithm", f"BFS from {source} visits {len(b['order'])} vertices", ", ".join(b["order"]))
    ui.show_trace(tbfs)
    ui.result("algorithm", f"DFS preorder from {source}", ", ".join(d["pre"]))
    ui.result("algorithm", f"DFS postorder; back edges {d['back_edges'] or 'none'}", ", ".join(d["post"]))
    ui.show_trace(tdfs)
    cycle, tcycle = mg.find_cycle(g)
    ui.result("algorithm", f"Cycle found: {' – '.join(cycle) if cycle else 'none; the graph is a forest'}")
    ui.show_trace(list(tcycle) + list(tc) + list(tcy) + list(tb))
    ui.result("calculation", f"Components: {len(comps)}; bipartite: {bip}" +
              ("" if bip else " (an odd cycle exists)"))
    ui.table(["component", "size", "members"],
             [[i + 1, len(c), pr.truncate(", ".join(c), 70)] for i, c in enumerate(comps)])

    st.subheader("Shortest paths (Dijkstra)")
    if pv["cost"] is None:
        ui.result("algorithm", f"{target} is unreachable from {source} in this graph. Direct d_H = {pv['direct']}.")
        ui.result("interpretation", "Unreachable means no chain of observed sequences links them at this threshold, "
                  "not that no relationship exists.")
    else:
        ui.result("algorithm", f"Shortest path cost {source} → {target}: {pv['cost']}; direct d_H = {pv['direct']}; "
                  f"detour = {pv['detour']}; {len(pv['paths'])} optimal path(s)",
                  "\n".join(" – ".join(p) for p in pv["paths"][:6]))
        ui.result("calculation", "Cost ≥ d_H always, because d_H is a metric and the triangle inequality applies. "
                  "A positive detour means every observed route changes some position more than once.")
    ui.show_trace(tp)
    ui.result("limitation", "A path through observed haplotypes is a route in this graph, not an inferred evolutionary "
              "trajectory.")

    st.subheader("How the graph changes with the threshold")
    ks = sorted(set(list(range(0, min(11, ctx["max_distance"] + 1))) + [ctx["max_distance"]]))
    rows, tsw = mg.threshold_sweep(seqs, ks)
    ui.plot_sweep(rows)
    header, body = pr.sweep_table(rows)
    ui.table(header, body)
    kstar, tks = mg.smallest_connecting_threshold(seqs)
    ui.result("calculation", f"The graph first becomes connected at k* = {kstar}. Section 7 checks that this equals "
              "the largest edge of every minimum spanning tree.")
    ui.show_trace(list(tsw) + list(tks))
    ui.result("limitation", "Any single k is an arbitrary choice, which is why the sweep is shown rather than one "
              "threshold. Vertices are observed haplotypes only, so a missing edge may mean an unsampled intermediate.")

    if len(ids) <= mg.BRUTE_FORCE_LIMIT:
        st.subheader("Connectivity")
        kappa, tk = mg.vertex_connectivity(g)
        lam, tl = mg.edge_connectivity(g)
        ui.result("algorithm", f"κ (vertex connectivity) = {kappa}; λ (edge connectivity) = {lam}; "
                  f"δ (minimum degree) = {min(deg.values())}. Whitney's inequality κ ≤ λ ≤ δ: "
                  f"{kappa <= lam <= min(deg.values())}")
        ui.show_trace(list(tk) + list(tl))
    else:
        ui.result("limitation", f"κ and λ are computed by exhaustive search, which is exponential; they are refused "
                  f"above {mg.BRUTE_FORCE_LIMIT} vertices and this graph has {len(ids)}.")


def trees(ctx):
    data, seqs = ctx["data"], ctx["seqs"]
    ids = list(seqs)
    ui.section_header(7, "Tree analysis (Unit VI)",
                      "What is the cheapest structure linking these sequences, and how does it differ from a "
                      "distance-based evolutionary tree?")
    complete, _ = mg.build_graph(seqs, 10 ** 9)
    kr, tkr = ta.kruskal(complete)
    root = st.selectbox("root the spanning tree at", ids, index=ids.index(data.reference.id), key="tr_root")
    prm, tprm = ta.prim(complete, root)

    st.subheader("Minimum spanning tree: Prim and Kruskal")
    c1, c2, c3 = st.columns(3)
    c1.metric("Kruskal total weight", kr["total"])
    c2.metric("Prim total weight", prm["total"])
    c3.metric("edges", len(kr["edges"]))
    ui.result("calculation", f"Both algorithms return total weight {kr['total']}. The MST weight is unique even when "
              "the tree is not.")
    count, tcount = ta.count_minimum_spanning_trees(complete)
    same = sorted(tuple(sorted(e[:2])) for e in kr["edges"]) == sorted(tuple(sorted(e[:2])) for e in prm["edges"])
    ui.result("algorithm", f"This graph has {count} distinct minimum spanning tree(s); Prim and Kruskal returned "
              f"{'the same tree' if same else 'different trees'}.")
    ui.show_trace(list(tkr) + list(tprm) + list(tcount), "Show the Kruskal accept/reject decisions and Prim's steps")
    for name, res in (("Kruskal", kr), ("Prim", prm)):
        rep, tv = ta.validate_spanning_tree(complete.vertices, res["edges"], complete)
        ui.result("algorithm", f"{name}: {len(res['edges'])} edges = |V| − 1, acyclic, spanning, all edges in the "
                  f"graph → valid spanning tree: {rep['is_spanning_tree']}")
    bott, tb = ta.bottleneck(kr["edges"])
    kstar, _ = mg.smallest_connecting_threshold(seqs)
    ui.result("calculation", f"Largest MST edge = {bott}, and the smallest connecting threshold k* = {kstar}. These "
              "agree, which is the bottleneck (minimax) property.")
    ui.show_trace(tb)
    ui.table(["edge", "weight"], [[f"{u} – {v}", w] for u, v, w in kr["edges"]])
    tree_graph = mg.Graph(complete.vertices, kr["edges"])
    ui.plot_graph(tree_graph, max(w for _, _, w in kr["edges"]), title="Minimum spanning tree (Kruskal)")

    st.subheader("The MST as a rooted tree (the root is an arbitrary choice)")
    rt, trt = ta.root_tree(complete.vertices, kr["edges"], root)
    stats, tst = ta.tree_statistics(rt)
    c1, c2, c3 = st.columns(3)
    c1.metric("height", stats["height"])
    c2.metric("leaves", len(stats["leaves"]))
    c3.metric("internal nodes", len(stats["internal"]))
    ui.result("algorithm", "Preorder", " – ".join(pr.truncate(x, 14) for x in ta.preorder(rt)[0]))
    ui.result("algorithm", "Postorder", " – ".join(pr.truncate(x, 14) for x in ta.postorder(rt)[0]))
    ui.result("algorithm", "Level order", " – ".join(pr.truncate(x, 14) for x in ta.level_order(rt)[0]))
    try:
        ta.inorder(rt)
    except ValueError as e:
        ui.result("limitation", f"In-order traversal is not defined here: {e} The UPGMA tree below is binary, so it "
                  "does have one.")
    ui.show_trace(list(trt) + list(tst))
    ui.result("limitation", "Rooting an MST is a presentation choice. Changing the root changes the height, the leaf "
              "count and every traversal, but not one edge of the tree.")

    st.subheader("UPGMA: a distance-based clustering tree")
    (mids, matrix), _ = mg.distance_matrix(seqs)
    result, tup = pa.upgma(mids, matrix)
    ui.plot_dendrogram(result, "UPGMA (average linkage); heights are d_H / 2")
    ui.table(["step", "merge", "distance", "height", "tied pairs"],
             [[m["step"], f"({pr.truncate(m['left'], 12)}, {pr.truncate(m['right'], 12)})",
               pr.fraction_text(m["distance"]), pr.fraction_text(m["height"]), m["tied_pairs"]]
              for m in result["merges"]])
    ui.result("algorithm", f"{result['tie_steps']} merge step(s) involved a tie, resolved by the fixed rule in "
              "docs/decisions.md (D8): merge the tied pair whose sorted member labels come first.")
    ui.show_trace(tup, "Show every merge and the clock assumption")
    st.code(result["newick"], language=None)
    ultra, tu = pa.is_ultrametric(mids, matrix)
    fit, tf = pa.fit_to_input(result, mids, matrix)
    ui.result("calculation", f"Input distances ultrametric: {ultra}. The UPGMA tree reproduces {fit['exact']} of "
              f"{fit['pairs']} distances exactly; the largest deviation is {pr.fraction_text(fit['max_deviation'])} "
              f"at {fit['worst_pair']}.")
    ui.show_trace(list(tu) + list(tf))
    ui.result("limitation", "UPGMA assumes a molecular clock: equal rates along all lineages, which makes distances "
              "ultrametric. Where the data are not ultrametric, its branch heights are averages, not times. The tree "
              "is a clustering hypothesis, not evidence of ancestry.")

    st.subheader("MST versus UPGMA")
    ms, _ = ta.edge_splits(complete.vertices, kr["edges"])
    us, _ = pa.upgma_splits(result, mids)
    cmp, tcmp = pa.compare_splits(ms, us)
    c1, c2, c3 = st.columns(3)
    c1.metric("MST bipartitions", cmp["mst"])
    c2.metric("UPGMA bipartitions", cmp["upgma"])
    c3.metric("shared", cmp["shared"])
    ui.show_trace(tcmp)
    ui.table(["", "minimum spanning tree", "UPGMA tree"],
             [["what it minimises", "total edge weight Σ d_H", "nothing; it clusters by average distance"],
              ["internal nodes", "observed sequences", "hypothetical ancestors"],
              ["root", "arbitrary", "the last merge"],
              ["assumption", "none beyond d_H being a metric", "a molecular clock (ultrametric distances)"],
              ["unique?", f"no — {count} minimum tree(s) here", "unique given the tie-break rule"],
              ["supports ancestry?", "no", "no, it is a hypothesis"]])
    ui.result("interpretation", "The two trees answer different questions, so they need not agree. An MST minimises "
              "total graph-edge weight; it does not establish evolutionary ancestry.")
    cay, _ = ta.cayley_count(len(ids))
    rbt, _ = ta.rooted_binary_topologies(len(ids))
    ui.result("calculation", f"For {len(ids)} labelled vertices there are n^(n−2) = {cay:,} labelled trees, and "
              f"(2n−3)!! = {rbt:,} rooted binary topologies on {len(ids)} labelled leaves. Both algorithms searched "
              "that space implicitly.")


SECTIONS = {
    "1 · Dataset explorer": dataset_explorer,
    "2 · Set theory and logic": set_theory,
    "3 · Counting and recursion": counting,
    "4 · Algebraic structures": algebra,
    "5 · Number theory": number_theory,
    "6 · Mutation graph": mutation_graph,
    "7 · Tree analysis": trees,
}
