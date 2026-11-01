"""Streamlit entry point: `streamlit run app/main.py`. Mathematics lives in genomath; this file only wires the UI."""
import sys
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # allows `streamlit run app/main.py` from any working directory
    sys.path.insert(0, str(ROOT))

from app import sections, ui  # noqa: E402
from genomath._core import SequenceError, mismatch_count  # noqa: E402
from genomath.sequence_processing import load_dataset, parse_fasta  # noqa: E402

DATASETS = {
    "Dataset A — synthetic (10 sequences, n = 12, hand-verified)":
        (ROOT / "data" / "synthetic" / "dataset_a.fasta", "syn01"),
    "Dataset B — GenBank H3N2 hemagglutinin CDS (36 sequences, n = 1701)":
        (ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta", None),
}


@st.cache_data(show_spinner=False)
def _load(path_text, reference_id):
    data, trace = load_dataset(Path(path_text), reference_id)
    text = Path(path_text).read_text(encoding="utf-8")
    records, _ = parse_fasta(text)
    return data, list(trace), len(records)


@st.cache_data(show_spinner=False)
def _load_uploaded(text, reference_id):
    from genomath.sequence_processing import build_dataset
    records, parse_trace = parse_fasta(text)
    data, trace = build_dataset(records, reference_id)
    return data, list(parse_trace) + list(trace), len(records)


def build_context(choice, uploaded, reference_override):
    if uploaded is not None:
        data, trace, read = _load_uploaded(uploaded, reference_override)
        label = "uploaded FASTA"
    else:
        path, default_ref = DATASETS[choice]
        data, trace, read = _load(str(path), reference_override or default_ref)
        label = choice
    seqs = data.seqs()
    distances = [mismatch_count(seqs[a], seqs[b]) for i, a in enumerate(seqs) for b in list(seqs)[i + 1:]]
    note = ("chosen by you" if reference_override else
            ("the dataset's designated reference" if uploaded is None and DATASETS[choice][1] else
             "the medoid: the haplotype with the smallest total distance to all others, ties broken by smallest id"))
    return {"data": data, "seqs": seqs, "label": label, "records_read": read, "load_trace": trace,
            "reference_note": note, "max_distance": max(distances) if distances else 0}


def main():
    st.set_page_config(page_title="Genomic Mutation & Evolution Analyzer", layout="wide")
    st.title("Genomic Mutation & Evolution Analyzer")
    st.caption("A discrete-mathematics framework for DNA sequence analysis · MM0903 course project · "
               "every result below is traceable to a definition, a derivation or a named algorithm")

    with st.sidebar:
        st.header("Data")
        choice = st.radio("dataset", list(DATASETS), key="ds")
        upload = st.file_uploader("or upload a FASTA file", type=["fasta", "fa", "txt"])
        uploaded_text = upload.getvalue().decode("utf-8", errors="replace") if upload else None
        st.header("Reference sequence")
        st.caption("The reference fixes every mutation set M(R,S). It does not affect any pairwise distance.")
        ref_choice = st.text_input("reference id (blank = default)", "", key="ref").strip() or None
        st.header("Section")
        section = st.radio("section", list(sections.SECTIONS), key="sec", label_visibility="collapsed")
        st.divider()
        ui.tag_legend()
        st.divider()
        st.caption("Scope: an educational framework. Not a diagnostic system, not a pathogenicity predictor, and not "
                   "evidence of evolutionary history. See docs/final_spec.md §2.")

    try:
        ctx = build_context(choice, uploaded_text, ref_choice)
    except SequenceError as exc:
        st.error(f"The dataset could not be prepared: {exc}")
        if getattr(exc, "trace", None):
            ui.show_trace(exc.trace, "Show the preprocessing decisions up to the failure", expanded=True)
        st.stop()
        return

    st.info(f"**{ctx['label']}** · {len(ctx['data'].records)} sequences kept · "
            f"{len(ctx['data'].haplotypes)} haplotypes · n = {ctx['data'].n} · reference {ctx['data'].reference.id}")
    sections.SECTIONS[section](ctx)


if __name__ == "__main__":
    main()
