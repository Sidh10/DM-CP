"""Expected values for Unit IV, transcribed by hand from docs/unit4_expected.md and docs/synthetic_expected.md §7."""

SOURCE = "docs/unit4_expected.md"
SYN = "docs/synthetic_expected.md §7"

MODULI = [16, 64, 7, 13, 17, 31, 61]  # provisional decision 8

ENC_A = {"syn01": 9378573, "syn02": 10427149, "syn03": 10459917, "syn04": 10460109, "syn05": 9378561,
         "syn06": 10427137, "syn07": 9392909, "syn08": 9378057, "syn09": 4786988}
COMP_A = {"syn01": 7398642, "syn02": 6350066, "syn03": 6317298, "syn04": 6317106, "syn05": 7398654,
          "syn06": 6350078, "syn07": 7384306, "syn08": 7399158, "syn09": 11990227}
CANON_A = {"syn01": 9378573, "syn02": 9378569, "syn03": 9378057, "syn04": 9181449, "syn05": 9378561,
           "syn06": 10427137, "syn07": 9386509, "syn08": 9378057, "syn09": 4786988}
BUCKETS_A = {16: [13, 13, 13, 13, 1, 1, 13, 9, 12], 7: [1, 5, 6, 2, 3, 0, 1, 3, 3],
             13: [9, 5, 0, 10, 10, 6, 6, 0, 11], 17: [13, 12, 4, 9, 1, 0, 1, 7, 9]}
GROUPS_A = {16: [["syn01", "syn02", "syn03", "syn04", "syn07"], ["syn05", "syn06"]],
            7: [["syn01", "syn07"], ["syn05", "syn08", "syn09"]],
            13: [["syn03", "syn08"], ["syn04", "syn05"], ["syn06", "syn07"]],
            17: [["syn04", "syn09"], ["syn05", "syn07"]]}
D1_PAIRS_A = [("syn01", "syn02"), ("syn01", "syn05"), ("syn02", "syn03"), ("syn02", "syn06"), ("syn03", "syn04"),
              ("syn05", "syn06")]
D1_COLLISIONS_A = {16: 4, 7: 0, 13: 0, 17: 0}  # SYN single-substitution table
PIGEON_A_M7 = {"bound": 2, "max_load": 3}  # §7

LAST_BASES = {16: 2, 64: 3, 7: None, 13: None, 17: None, 31: None, 61: None}  # §3
T_MIN = {16: 2, 64: 3, 8: 1, 12: 1, 6: 1, 3: 0, 1: 0, 7: None, 10: None}  # §3

SYN01_3MERS = [35, 15, 60, 49, 6, 27, 44, 48, 3, 13]  # §1
SYN01_CENSUS = {3: {"distinct": 10, "canonical": 5, "palindromic": 0}, 4: {"distinct": 9, "canonical": 5, "palindromic": 1}}  # §2
SPACE_CANONICAL = {4: 136, 6: 2080, 8: 32896}  # §2

GCD_252_198 = (18, 4, -5)  # §5
INVERSES = {"11 mod 7": 2, "7 mod 11": 8, "4 mod 7": 2, "4 mod 16": "none"}  # §5
CONGRUENCES = {"3x≡6 (mod 9)": [2, 5, 8], "2x≡1 (mod 4)": []}  # §5
CRT = {"GAT from (0 mod 7, 2 mod 11)": 35, "(7,11) round-trips": 64, "(5,13) round-trips": 64,
       "(3,5) round-trips": 15, "(4,6)": "refused"}  # §6
ORDERS = {7: 3, 13: 6, 17: 4, 31: 5, 61: 30}  # §4
PHI = {7: 6, 13: 12, 17: 16, 31: 30, 61: 60, 77: 60, 16: 8, 64: 32}  # §4
EULER_4 = {7: True, 13: True, 17: True, 31: True, 61: True, 77: True, 16: False}  # §4: 4^φ(m) ≡ 1 (mod m)
CANCEL = {"string": "GATTAAGTAATG", "d_H": 2, "collide": [7, 13]}  # §4
PIGEON_B = {7: 5, 13: 3, 16: 2, 17: 2, 31: 1, 61: 1, 64: 1}  # §7
