"""Expected values for Unit III, transcribed by hand from docs/unit3_expected.md (section noted on each)."""

SOURCE = "docs/unit3_expected.md"

CAYLEY_Z2 = {"I": ["I", "RC"], "RC": ["RC", "I"]}  # §1
CAYLEY_V4 = {"I": ["I", "R", "C", "RC"], "R": ["R", "I", "RC", "C"], "C": ["C", "RC", "I", "R"],
             "RC": ["RC", "C", "R", "I"]}  # §1
AXIOMS = {"closure": True, "associativity": True, "identity": True, "inverses": True, "commutative": True}  # §1
ORDERS_V4 = {"I": 1, "R": 2, "C": 2, "RC": 2}  # §1
CYCLIC_V4 = False  # §1
ISO_Z2xZ2, HOMOMORPHIC_BIJECTIONS_TO_Z4 = True, 0  # §1
SUBGROUPS_V4 = [["I"], ["I", "R"], ["I", "C"], ["I", "RC"], ["I", "R", "C", "RC"]]  # §2
COSETS_Z2_IN_V4 = [["I", "RC"], ["R", "C"]]  # §2
FAITHFUL = {0: [["I", "R", "C", "RC"]], 1: [["I", "R"], ["C", "RC"]], 2: [["I"], ["R"], ["C"], ["RC"]]}  # §1
FIX = {  # §3, n = 0..8
    "I": [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536],
    "R": [1, 4, 4, 16, 16, 64, 64, 256, 256],
    "C": [1, 0, 0, 0, 0, 0, 0, 0, 0],
    "RC": [1, 0, 4, 0, 16, 0, 64, 0, 256],
}
ORBITS_Z2 = [1, 2, 10, 32, 136, 512, 2080, 8192, 32896]  # §4, n = 0..8
ORBITS_V4 = [1, 2, 6, 20, 72, 272, 1056, 4160, 16512]  # §4
ORBITS_N12 = {"Z2": 8_390_656, "V4": 4_196_352}  # §4
V4_SIZES = {1: {2: 2}, 2: {2: 4, 4: 2}, 3: {2: 8, 4: 12}, 4: {2: 16, 4: 56}}  # §4
COMPLEMENT_PERM = {"cycles": [("A", "T"), ("C", "G")], "order": 2, "parity": "even"}  # §5
REVERSAL_12 = {"transpositions": 6, "order": 2, "parity": "even"}  # §5
MONOID = {1: (5, ["·"]), 2: (25, ["··"]), 3: (125, ["···"])}  # §6: (size, units)
MONOID_GROUP, MONOID_COMMUTATIVE = False, False  # §6
CONCAT = {"C(xy) = C(x)C(y)": True, "R(xy) = R(y)R(x)": True, "RC(xy) = RC(y)RC(x)": True}  # §7
A_Z2_CLASSES, A_RC_PAIR = 8, ["syn03", "syn08"]  # §8
A_SYN01_V4 = {"size": 2, "stabiliser": ["I", "RC"]}  # §8
A_V4_SIZE_SYN03_SYN09 = [4, 4]  # §8
A_CANONICAL_SYN03 = "GATTACGCAAGC"  # §8
A_R_FIXED = []  # §8
B = {"RC_pairs": 0, "RC_fixed": 0, "Z2_classes": 30}  # §8
