"""Expected values for Unit VI, transcribed by hand from docs/unit6_expected.md (section noted on each)."""
from fractions import Fraction

SOURCE = "docs/unit6_expected.md"

KRUSKAL_EDGES = [("01", "02", 1), ("01", "05", 1), ("02", "03", 1), ("02", "06", 1), ("03", "04", 1), ("01", "07", 2),
                 ("01", "08", 2), ("01", "09", 6)]  # §1
KRUSKAL_COUNTS = {"examined": 29, "accepted": 8, "rejected": 21}  # §1
PRIM_ORDER = ["01", "02", "03", "04", "05", "06", "07", "08", "09"]  # §2
PRIM_EDGES = [("01", "02", 1), ("02", "03", 1), ("03", "04", 1), ("01", "05", 1), ("02", "06", 1), ("01", "07", 2),
              ("01", "08", 2), ("01", "09", 6)]  # §2
TOTALS = {"prim": 15, "kruskal": 15}  # §1–2
VALID = {"prim": True, "kruskal": True}  # §3
MST_COUNT = 8  # §4
BOTTLENECK = {"max MST edge": 6, "k*": 6}  # §3
KIRCHHOFF_K9 = 4_782_969  # §4 (Cayley)

ROOT_01 = {"height": 3, "leaves": ["04", "05", "06", "07", "08", "09"], "internal": ["01", "02", "03"],
           "preorder": ["01", "02", "03", "04", "06", "05", "07", "08", "09"],
           "postorder": ["04", "03", "06", "02", "05", "07", "08", "09", "01"],
           "level": ["01", "02", "05", "07", "08", "09", "03", "06", "04"]}  # §5
ROOT_04 = {"height": 4, "leaves": ["05", "06", "07", "08", "09"], "internal": ["01", "02", "03", "04"],
           "preorder": ["04", "03", "02", "01", "05", "07", "08", "09", "06"],
           "postorder": ["05", "07", "08", "09", "01", "06", "02", "03", "04"]}  # §5

UPGMA_MERGES = [["U1", "01", "02", Fraction(1), Fraction(1, 2), 6], ["U2", "03", "04", Fraction(1), Fraction(1, 2), 2],
                ["U3", "05", "06", Fraction(1), Fraction(1, 2), 1], ["U4", "U1", "U3", Fraction(3, 2), Fraction(3, 4), 1],
                ["U5", "U4", "U2", Fraction(5, 2), Fraction(5, 4), 2], ["U6", "U5", "07", Fraction(19, 6), Fraction(19, 12), 2],
                ["U7", "U6", "08", Fraction(23, 7), Fraction(23, 14), 1],
                ["U8", "U7", "09", Fraction(61, 8), Fraction(61, 16), 1]]  # §6
NEWICK = "((((((01,02),(05,06)),(03,04)),07),08),09);"  # §6
TIE_STEPS = 4  # §6
REVERSE_TIE_CHANGES_TOPOLOGY = True  # §6
INORDER = ["01", "U1", "02", "U4", "05", "U3", "06", "U5", "03", "U2", "04", "U6", "07", "U7", "08", "U8", "09"]  # §6
ULTRAMETRIC = (False, ("01", "02", "03"))  # §6
FIT = {"exact": 3, "pairs": 36, "max_deviation": Fraction(12, 7), "worst_pair": ("04", "08")}  # §6 (corrected)
SPLITS = {"mst": 2, "upgma": 6, "shared": 1, "symmetric_difference": 6}  # §7
COUNTING = {"cayley n^(n-2)": 4_782_969, "(2n-3)!!": 2_027_025}  # §8
