from pathlib import Path

import pytest

from genomath.sequence_processing import load_dataset
from genomath.trace import KINDS

ROOT = Path(__file__).resolve().parents[1]
A_PATH = ROOT / "data" / "synthetic" / "dataset_a.fasta"
A_INVALID = ROOT / "data" / "synthetic" / "dataset_a_invalid.fasta"
B_PATH = ROOT / "data" / "public_genomic" / "h3n2_ha_cds.fasta"

# docs/synthetic_expected.md §3, transcribed by hand
HAPS_A = ["syn01", "syn02", "syn03", "syn04", "syn05", "syn06", "syn07", "syn08", "syn09"]
D_A = [
    [0, 1, 2, 3, 1, 2, 2, 2, 6],
    [1, 0, 1, 2, 2, 1, 3, 3, 7],
    [2, 1, 0, 1, 3, 2, 3, 4, 8],
    [3, 2, 1, 0, 4, 3, 4, 5, 9],
    [1, 2, 3, 4, 0, 1, 3, 2, 7],
    [2, 1, 2, 3, 1, 0, 4, 3, 8],
    [2, 3, 3, 4, 3, 4, 0, 4, 8],
    [2, 3, 4, 5, 2, 3, 4, 0, 8],
    [6, 7, 8, 9, 7, 8, 8, 8, 0],
]


@pytest.fixture(scope="session")
def dataset_a():
    return load_dataset(A_PATH, reference_id="syn01")[0]


@pytest.fixture(scope="session")
def dataset_b():
    return load_dataset(B_PATH)[0]


def assert_trace(trace):
    assert all(step["kind"] in KINDS for step in trace)
    assert len(trace) > 0
