"""Expected values for Unit II, transcribed by hand from docs/unit2_expected.md (section noted on each)."""

SOURCE = "docs/unit2_expected.md"

M_TABLE = [  # §2
    [1],
    [1, 3],
    [1, 6, 9],
    [1, 9, 27, 27],
    [1, 12, 54, 108, 81],
    [1, 15, 90, 270, 405, 243],
    [1, 18, 135, 540, 1215, 1458, 729],
    [1, 21, 189, 945, 2835, 5103, 5103, 2187],
    [1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561],
]
POWERS_OF_4 = [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536]  # §1
CELLS = 45  # §2
M12 = {1: 36, 2: 594, 3: 5940}  # docs/synthetic_expected.md §8
BALLS = {(8, 2): 277, (12, 2): 631, (1701, 1): 5104}  # §2
GUARD = {(12, 6): ("refused", 673_596), (1701, 2): ("refused", 13_012_650), (1701, 1): ("allowed", 5103)}  # §2
WHOLE_SPACE_8 = 65_536  # §2: B(8,8)
DIGITS_4_1701 = 1025  # §1
ALL_FOUR = {4: 24, 5: 240, 6: 1560, 7: 8400, 8: 40824, 12: 14_676_024}  # §4
ORDERINGS_01_04 = {"count": 6, "fully_observed": 1}  # §3
P_12_3 = 1320  # §3
PIGEON = {"max_forced_k": 5, "forced_k5": True, "forced_k6": False, "clean_with_repeat_5mer": 34}  # §5
