"""The app's data shaping lives in app.presentation, which must stay importable without Streamlit."""
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from app import presentation as pr
from genomath import mutation_graph as mg
from genomath import phylogenetic_analysis as pa
from genomath.trace import KINDS, Trace

ROOT = Path(__file__).resolve().parents[1]


def imported_modules(path):
    import ast
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split(".")[0])
    return names


def test_presentation_has_no_streamlit_or_oracle_imports():
    # checked on the import statements, not on the raw text: "Biological interpretation" contains "Bio"
    imported = imported_modules(ROOT / "app" / "presentation.py")
    assert imported & {"streamlit", "networkx", "sympy", "Bio", "matplotlib", "genomath"} == set()
    assert imported == {"math"}


def test_genomath_never_imports_the_app():
    for module in (ROOT / "genomath").glob("*.py"):
        assert "app" not in imported_modules(module)


def test_presentation_imports_in_a_clean_interpreter():
    # proves the module is testable on its own: a fresh process importing it must not pull in streamlit
    code = "import app.presentation, sys; print('streamlit' in sys.modules)"
    r = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "False"


def test_every_trace_kind_has_a_spec_category():
    assert set(pr.TAGS) == set(KINDS) == set(pr.TAG_ORDER)
    labels = [pr.TAGS[k][0] for k in pr.TAG_ORDER]
    assert labels == ["Raw observation", "Mathematical calculation", "Algorithmic output",
                      "Biological interpretation", "Limitation / assumption"]
    assert pr.tag_label("unknown-kind") == "unknown-kind"


def test_trace_rows_expand_capped_steps():
    t = Trace().add_steps("calculation", [f"s{i}" for i in range(25)])
    rows = pr.trace_rows(t)
    assert len(rows) == 11
    assert rows[-1]["omitted_count"] == 15 and rows[-1]["omitted"][0] == "s10"
    assert all(r["label"] for r in rows) and "omitted" not in rows[-1]["data"]
    assert pr.trace_summary(t) == {"calculation": 11}


def test_matrix_and_sweep_tables(dataset_a):
    seqs = {k[-2:]: v for k, v in dataset_a.seqs().items()}
    (ids, matrix), _ = mg.distance_matrix(seqs)
    header, body = pr.matrix_table(ids, matrix)
    assert header == [""] + ids
    assert body[0][0] == ids[0] and body[0][1] == ""          # diagonal blanked
    assert body[0][2] == matrix[0][1]
    assert pr.matrix_table(ids, matrix, blank_diagonal=False)[1][0][1] == 0
    rows, _ = mg.threshold_sweep(seqs, [0, 1, 2])
    h, b = pr.sweep_table(rows)
    assert h[0] == "k" and [r[0] for r in b] == [0, 1, 2]


def test_cayley_table_shape():
    from genomath import algebraic_structures as al
    table, _ = al.cayley_table(al.Z2)
    header, body = pr.cayley_table(al.Z2, table)
    assert header == ["∘", "I", "RC"] and body == [["I", "I", "RC"], ["RC", "RC", "I"]]


@pytest.mark.parametrize("k", [1, 2, 6])
def test_spring_layout_is_deterministic_and_in_frame(dataset_a, k):
    seqs = {key[-2:]: v for key, v in dataset_a.seqs().items()}
    g, _ = mg.build_graph(seqs, k)
    edges = [(u, v) for u, v, _ in g.edges()]
    a = pr.spring_layout(g.vertices, edges)
    b = pr.spring_layout(g.vertices, edges)
    assert a == b                                            # no RNG: the same input gives the same picture
    assert set(a) == set(g.vertices)
    assert all(0.0 <= x <= 1.0 and 0.0 <= y <= 1.0 for x, y in a.values())
    assert len({(round(x, 3), round(y, 3)) for x, y in a.values()}) == len(g.vertices)


def test_spring_layout_degenerate_cases():
    assert pr.spring_layout([], []) == {}
    assert pr.spring_layout(["only"], []) == {"only": (0.5, 0.5)}
    two = pr.spring_layout(["a", "b"], [("a", "b")])
    assert len(two) == 2 and all(0 <= c <= 1 for p in two.values() for c in p)


def test_dendrogram_layout_matches_the_tree(dataset_a):
    seqs = {k[-2:]: v for k, v in dataset_a.seqs().items()}
    (ids, matrix), _ = mg.distance_matrix(seqs)
    result, _ = pa.upgma(ids, matrix)
    layout = pr.dendrogram_layout(result)
    assert sorted(layout["leaves"]) == sorted(ids)
    assert len(layout["leaves"]) == len(ids) and len(layout["bars"]) == len(result["merges"])
    assert layout["x"][layout["leaves"][0]] == 0
    for bar in layout["bars"]:                               # a parent is never drawn below its children
        assert bar["y"] >= bar["y_left"] and bar["y"] >= bar["y_right"]
    heights = [b["y"] for b in layout["bars"]]
    assert heights == sorted(heights)


def test_edge_style_and_text_helpers():
    assert pr.edge_style(1, 1)["width"] == pytest.approx(2.0)
    near, far = pr.edge_style(1, 5), pr.edge_style(5, 5)
    assert near["width"] > far["width"] and near["alpha"] > far["alpha"]
    assert pr.truncate("abc", 10) == "abc" and pr.truncate("abcdef", 4) == "abc…"
    assert pr.fraction_text(Fraction(3, 2)) == "3/2 (≈ 1.5)"
    assert pr.fraction_text(Fraction(4, 2)) == "2" and pr.fraction_text(7) == "7"
    assert pr.describe_dropped([("x", "why")]) == [{"id": "x", "reason": "why"}]


@pytest.mark.parametrize("k", [1, 2, 3])
def test_spring_layout_places_adjacent_vertices_closer(dataset_a, k):
    """The point of a force-directed layout: edges pull, so neighbours end up nearer than non-neighbours."""
    from math import dist
    seqs = {key[-2:]: v for key, v in dataset_a.seqs().items()}
    g, _ = mg.build_graph(seqs, k)
    edges = [(u, v) for u, v, _ in g.edges()]
    pos = pr.spring_layout(g.vertices, edges)
    linked = {frozenset(e) for e in edges}
    adjacent = [dist(pos[u], pos[v]) for u, v in edges]
    apart = [dist(pos[x], pos[y]) for i, x in enumerate(g.vertices) for y in g.vertices[i + 1:]
             if frozenset((x, y)) not in linked]
    assert sum(adjacent) / len(adjacent) < sum(apart) / len(apart)
    assert max(adjacent) < max(apart)
