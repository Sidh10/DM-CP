"""Expected values for the Dataset B pipeline run (Experiment 8).

Sources are independent of `genomath`:
- QC = docs/dataset_b_qc.md, produced by scripts/qc_public_data.py, which does not import genomath;
- hand = a hand derivation in docs/unit*_expected.md, named per entry.
Values that only `genomath` can produce are not listed here; the experiment reports them as observations, and the
cross-experiment reproductions are tagged as consistency checks.
"""

QC = "docs/dataset_b_qc.md (independent script)"
U2 = "docs/unit2_expected.md"
U3 = "docs/unit3_expected.md"
U4 = "docs/unit4_expected.md"

# --- preprocessing (QC) ---
RECORDS_READ = 36
DROPPED = ["MH246583.1", "MK181701.1"]
CLEAN = 34
HAPLOTYPES = 30
LENGTH = 1701
IDENTICAL_GROUPS = [["CY249987.1", "MH082044.1", "MH601290.1"], ["MH082676.1", "MH294810.1", "MH358590.1"]]

# --- Unit I (QC) ---
PAIRS = 435
D1_PAIRS = 8
DISTANCE_RANGE = {"min": 1, "max": 57}

# --- Unit II (hand) ---
DIGITS_4_POW_N = 1025          # §1: decimal digits of 4^1701
NEIGHBOURS_AT_1 = 5103         # §2: 3n possible single-substitution neighbours
BALL_RADIUS_1 = 5104           # §2
MAX_FORCED_KMER = 5            # §5: largest k with 1701 - k + 1 > 4^k
FORCED_5MER_SEQUENCES = 34     # §5: every clean sequence must repeat a 5-mer

# --- Unit III (hand) ---
RC_PAIRS = 0                   # §8: start/stop codon argument
RC_FIXED = 0                   # §8
Z2_CLASSES = 30                # §8
SPACE_CANONICAL_4MERS = 136    # §2 of unit4 / Burnside at k = 4
# D9, §9
FPF_INVOLUTIONS = 3
ORDERED_PAIRS_COVERED = 12
TRANSITION_PAIRS = 2
TRANSVERSION_PAIRS = 4
BASELINE_RATIO = 0.5

# --- Unit IV (hand/theorem) ---
LAST_BASES_SEEN = {16: 2, 64: 3, 7: None, 13: None, 17: None, 31: None, 61: None}   # §3
D1_COLLISIONS = {16: 8, 64: 8, 7: 0, 13: 0, 17: 0, 31: 0, 61: 0}                    # §3 theorem
