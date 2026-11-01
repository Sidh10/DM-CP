import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from conftest import A_INVALID, A_PATH, B_PATH, D_A, HAPS_A, assert_trace
from genomath._core import SequenceError, complement, mismatch_count, reverse_complement
from genomath.sequence_processing import (Record, build_dataset, choose_reference, collapse_haplotypes, hamming,
                                          parse_fasta, position_map, read_fasta, validate)
from genomath.trace import Trace


def recs(*pairs):
    return [Record(i, "", s) for i, s in pairs]


def test_parse_multiline_blank_lines_and_description():
    records, trace = parse_fasta(">x first one\nACG\n\nTAC\n>y\nGG\n")
    assert [(r.id, r.description, r.seq) for r in records] == [("x", "first one", "ACGTAC"), ("y", "", "GG")]
    assert_trace(trace)


@pytest.mark.parametrize("text", ["ACGT\n>x\nA", ">\nACGT", ">x\nA\n>x\nC"])
def test_parse_rejects_malformed(text):
    with pytest.raises(SequenceError):
        parse_fasta(text)


@pytest.mark.parametrize("path", [A_PATH, B_PATH])
def test_parser_agrees_with_biopython(path):
    ours, _ = read_fasta(path)
    with path.open(encoding="utf-8") as fh:
        theirs = [(r.id, str(r.seq)) for r in SeqIO.parse(fh, "fasta")]
    assert [(r.id, r.seq) for r in ours] == theirs


def test_validate_uppercases_and_drops_with_reasons():
    out, trace = validate(recs(("a", "acgt"), ("b", "ACNT"), ("c", ""), ("d", "ACGT")))
    assert [r.seq for r in out["kept"]] == ["ACGT", "ACGT"]
    assert out["dropped"] == [("b", "non-ACGT characters Nx1"), ("c", "empty sequence")]
    assert out["n"] == 4
    assert "limitation" in trace.kinds()


def test_validate_empty_and_all_dropped():
    with pytest.raises(SequenceError, match="no sequences"):
        validate([])
    with pytest.raises(SequenceError, match="every record was dropped") as e:
        validate(recs(("a", "NNNN")))
    assert any("dropped a" in s["text"] for s in e.value.trace)


def test_unequal_lengths_rejected_with_trace():
    with pytest.raises(SequenceError, match="different lengths") as e:
        validate(recs(("a", "ACGT"), ("b", "ACG")))
    assert "limitation" in e.value.trace.kinds()


def test_invalid_companion_file():
    a, _ = read_fasta(A_PATH)
    bad, _ = read_fasta(A_INVALID)
    out, _ = validate(a + bad[:1])  # bad01 has an N
    assert out["dropped"] == [("bad01", "non-ACGT characters Nx1")]
    with pytest.raises(SequenceError, match="different lengths"):
        validate(a + bad)  # bad02 has length 11


def test_haplotypes_dataset_a(dataset_a):
    assert dataset_a.n == 12
    assert [h.id for h in dataset_a.haplotypes] == HAPS_A
    assert {h.id: h.members for h in dataset_a.haplotypes}["syn03"] == ("syn03", "syn10")
    assert dataset_a.dropped == ()


def test_collapse_keeps_first_member_as_name():
    haps, trace = collapse_haplotypes(recs(("z", "AC"), ("a", "GT"), ("m", "AC")))
    assert [(h.id, h.members) for h in haps] == [("z", ("z", "m")), ("a", ("a",))]
    assert_trace(trace)


def test_single_sequence_dataset():
    data, _ = build_dataset(recs(("only", "ACGT")))
    assert [h.id for h in data.haplotypes] == ["only"] and data.reference.id == "only"


def test_hamming_matrix_matches_hand_values(dataset_a):
    seqs = dataset_a.seqs()
    for i, a in enumerate(HAPS_A):
        for j, b in enumerate(HAPS_A):
            d, trace = hamming(seqs[a], seqs[b])
            assert d == D_A[i][j] == mismatch_count(seqs[a], seqs[b])
    assert_trace(trace)


def test_hamming_unequal_lengths():
    with pytest.raises(SequenceError):
        hamming("ACGT", "ACG")


def test_reverse_complement_against_biopython(dataset_a, dataset_b):
    for h in dataset_a.haplotypes + dataset_b.haplotypes:
        assert reverse_complement(h.seq) == str(Seq(h.seq).reverse_complement())
        assert reverse_complement(reverse_complement(h.seq)) == h.seq
    seqs = dataset_a.seqs()
    assert reverse_complement(seqs["syn03"]) == seqs["syn08"]
    assert reverse_complement(seqs["syn01"]) == seqs["syn01"]
    with pytest.raises(SequenceError):
        complement("ACGN")


def test_medoid_reference_dataset_a(dataset_a):
    # row sums of the §3 matrix: syn01 = 19 is the unique minimum
    ref, trace = choose_reference(dataset_a.haplotypes)
    assert ref.id == "syn01"
    assert next(s for s in trace if s["data"].get("totals"))["data"]["totals"]["syn01"] == 19


def test_medoid_tie_breaks_by_smallest_id():
    haps, _ = collapse_haplotypes(recs(("b", "AC"), ("a", "AG")))
    assert choose_reference(haps)[0].id == "a"


def test_reference_by_member_id(dataset_a):
    assert choose_reference(dataset_a.haplotypes, "syn10")[0].id == "syn03"
    with pytest.raises(SequenceError):
        choose_reference(dataset_a.haplotypes, "nope")


def test_position_map(dataset_a):
    m, trace = position_map(dataset_a.reference.seq)
    assert m[1] == "G" and m[12] == "C" and len(m) == 12
    assert_trace(trace)


def test_dataset_b_matches_independent_qc(dataset_b):
    # values from scripts/qc_public_data.py, which does not use genomath
    assert dataset_b.n == 1701
    assert len(dataset_b.records) == 34
    assert [d[0] for d in dataset_b.dropped] == ["MH246583.1", "MK181701.1"]
    assert len(dataset_b.haplotypes) == 30


def test_trace_rejects_unknown_kind():
    with pytest.raises(ValueError):
        Trace().add("guess", "x")


def test_trace_add_steps_caps_and_keeps_the_rest():
    from genomath.trace import STEP_CAP
    t = Trace().add_steps("calculation", [f"s{i}" for i in range(25)])
    assert [s["text"] for s in t[:STEP_CAP]] == [f"s{i}" for i in range(STEP_CAP)]
    assert len(t) == STEP_CAP + 1 and t[-1]["data"]["omitted"] == [f"s{i}" for i in range(STEP_CAP, 25)]
    assert len(Trace().add_steps("calculation", ["a", "b"])) == 2


def test_long_inputs_produce_short_traces(dataset_b):
    from genomath import number_theory as nt
    from genomath.set_logic import mutation_set
    seq = dataset_b.records[0].seq
    x, trace = nt.encode(seq)
    assert len(trace) <= 12 and len(trace[-2]["data"]["omitted"]) == 1701 - 10
    assert len(nt.decode(x, 1701)[1]) <= 12
    far = max(dataset_b.seqs().values(), key=lambda s: sum(a != b for a, b in zip(s, dataset_b.reference.seq)))
    assert len(mutation_set(dataset_b.reference.seq, far)[1]) <= 12
