"""Expected values for Dataset A, transcribed by hand from docs/synthetic_expected.md (section noted on each)."""

SOURCE = "docs/synthetic_expected.md"

HAPLOTYPES = ["syn01", "syn02", "syn03", "syn04", "syn05", "syn06", "syn07", "syn08", "syn09"]  # §3
IDENTICAL_GROUP = ["syn03", "syn10"]  # §3

M_STAR = {  # §1
    "syn01": [],
    "syn02": [(2, "C")],
    "syn03": [(2, "C"), (5, "G")],
    "syn04": [(2, "C"), (5, "G"), (9, "T")],
    "syn05": [(11, "A")],
    "syn06": [(2, "C"), (11, "A")],
    "syn07": [(5, "C"), (7, "A")],
    "syn08": [(8, "C"), (11, "G")],
    "syn09": [(1, "C"), (3, "G"), (4, "C"), (6, "A"), (10, "G"), (12, "A")],
}

SET_OPS_M3_M7 = {  # §2
    "intersection": [5], "union": [2, 5, 7], "a_minus_b": [2], "b_minus_a": [7], "symmetric_difference": [2, 7],
}
ALLELE_INTERSECTION_M3_M7 = []  # §2
CHAIN = {"M2 ∩ M4": [2], "M4 \\ M3": [9]}  # §2

DECOMPOSITION_03_07 = {"sym_diff": [2, 7], "shared_differing": [5], "d_H": 3}  # §2

D = [  # §3
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

R1_PROPERTIES = {"reflexive": True, "symmetric": True, "transitive": False}  # §3
R1_COUNTEREXAMPLE = ("syn01", "syn02", "syn03")  # §3
R1_SIZE = 9 + 2 * 6  # §3: 9 loops + both orders of the six d_H = 1 pairs

PREDICATES = {  # §2, P = d_H(R,S) <= 2, Q = S_2 = C
    "P ∧ Q": ["syn02", "syn03", "syn06"],
    "P ∧ ¬Q": ["syn01", "syn05", "syn07", "syn08"],
    "¬P": ["syn04", "syn09"],
}

ORIENTATION_CLASS_OF_SYN03 = ["syn03", "syn08", "syn10"]  # §6 (records, not haplotypes)
RC_FIXED = ["syn01"]  # §6

COMPONENTS = {0: 9, 1: 4, 2: 2, 5: 2, 6: 1}  # §4

INVALID_BAD01_REASON = "non-ACGT characters Nx1"  # dataset_a_invalid.fasta header
