"""Pure display helpers for the dashboard: no Streamlit import, so every transformation here is unit-testable."""
from math import cos, pi, sin, sqrt

# the five result categories the spec requires the UI to distinguish, keyed by the trace kinds of decision D1
TAGS = {
    "observation": ("Raw observation", "read directly from the input data"),
    "calculation": ("Mathematical calculation", "a value derived from a stated definition or formula"),
    "algorithm": ("Algorithmic output", "produced by a named algorithm, step by step"),
    "interpretation": ("Biological interpretation", "what the mathematics may mean biologically, stated cautiously"),
    "limitation": ("Limitation / assumption", "what this result does not show, or what it assumes"),
}
TAG_ORDER = list(TAGS)


def tag_label(kind):
    return TAGS[kind][0] if kind in TAGS else kind


def trace_rows(trace):
    """Flatten a trace into display rows, expanding a capped '… N more' step into its stored texts."""
    rows = []
    for step in trace:
        omitted = step.get("data", {}).get("omitted")
        rows.append({"kind": step["kind"], "label": tag_label(step["kind"]), "text": step["text"],
                     "omitted_count": len(omitted) if omitted else 0,
                     "data": {k: v for k, v in step.get("data", {}).items() if k != "omitted"},
                     "omitted": list(omitted) if omitted else []})
    return rows


def trace_summary(trace):
    counts = {k: 0 for k in TAG_ORDER}
    for step in trace:
        counts[step["kind"]] = counts.get(step["kind"], 0) + 1
    return {k: v for k, v in counts.items() if v}


def matrix_table(ids, matrix, blank_diagonal=True):
    """Square matrix as a header row plus body rows, with the diagonal optionally blanked for readability."""
    header = [""] + list(ids)
    body = []
    for i, name in enumerate(ids):
        cells = []
        for j, value in enumerate(matrix[i]):
            cells.append("" if blank_diagonal and i == j else value)
        body.append([name] + cells)
    return header, body


def sweep_table(rows, columns=("k", "edges", "components", "largest", "isolated", "max_degree", "cyclomatic")):
    header = list(columns)
    return header, [[row[c] for c in columns] for row in rows]


def cayley_table(elements, table):
    header = ["∘"] + list(elements)
    return header, [[a] + [table[(a, b)] for b in elements] for a in elements]


def spring_layout(vertices, edges, iterations=200, seed=7):
    """Deterministic force-directed layout. Positions are for readability only and carry no mathematical meaning."""
    vertices = list(vertices)
    n = len(vertices)
    if n == 0:
        return {}
    if n == 1:
        return {vertices[0]: (0.5, 0.5)}
    # start on a circle so the result is reproducible without a random number generator
    pos = {v: [0.5 + 0.4 * cos(2 * pi * i / n + seed), 0.5 + 0.4 * sin(2 * pi * i / n + seed)]
           for i, v in enumerate(vertices)}
    adj = {v: set() for v in vertices}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    k = 1.0 / sqrt(n)
    temp = 0.1
    for step in range(iterations):
        disp = {v: [0.0, 0.0] for v in vertices}
        for i, a in enumerate(vertices):
            for b in vertices[i + 1:]:
                dx, dy = pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]
                dist = max(sqrt(dx * dx + dy * dy), 1e-4)
                force = k * k / dist  # repulsion between every pair
                ux, uy = dx / dist * force, dy / dist * force
                disp[a][0] += ux
                disp[a][1] += uy
                disp[b][0] -= ux
                disp[b][1] -= uy
        for a in vertices:
            for b in adj[a]:
                dx, dy = pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]
                dist = max(sqrt(dx * dx + dy * dy), 1e-4)
                force = dist * dist / k  # attraction along edges
                disp[a][0] -= dx / dist * force
                disp[a][1] -= dy / dist * force
        for a in vertices:  # weak pull to the centre keeps disconnected components in frame
            disp[a][0] += (0.5 - pos[a][0]) * 0.05
            disp[a][1] += (0.5 - pos[a][1]) * 0.05
        for a in vertices:
            dx, dy = disp[a]
            dist = max(sqrt(dx * dx + dy * dy), 1e-4)
            pos[a][0] += dx / dist * min(dist, temp)
            pos[a][1] += dy / dist * min(dist, temp)
        temp = max(temp * 0.95, 1e-3)
    return normalise(pos)


def normalise(pos, margin=0.06):
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    span_x, span_y = max(xs) - min(xs), max(ys) - min(ys)
    scale = 1.0 - 2 * margin

    def axis(value, lo, span):
        return 0.5 if span < 1e-9 else margin + scale * (value - lo) / span

    return {v: (axis(p[0], min(xs), span_x), axis(p[1], min(ys), span_y)) for v, p in pos.items()}


def dendrogram_layout(result):
    """x = leaf order from an in-order walk, y = UPGMA height; returns leaf positions and the bars to draw."""
    tree, height = result["tree"], result["height"]
    leaves = []

    def walk(node):
        kids = tree.children.get(node)
        if not kids:
            leaves.append(node)
            return
        walk(kids[0])
        walk(kids[1])

    walk(result["root"])
    x = {leaf: i for i, leaf in enumerate(leaves)}
    bars = []
    for node in [m["node"] for m in result["merges"]]:  # merge order: a node's children always come earlier
        left, right = tree.children[node]
        x[node] = (x[left] + x[right]) / 2
        bars.append({"node": node, "left": left, "right": right,
                     "x_left": x[left], "x_right": x[right], "y": float(height[node]),
                     "y_left": float(height[left]), "y_right": float(height[right])})
    return {"leaves": leaves, "x": x, "bars": bars}


def edge_style(weight, k):
    """Thicker, darker lines for closer sequences; k is the threshold in force."""
    if k <= 1:
        return {"width": 2.0, "alpha": 0.9}
    share = (k - weight) / max(k - 1, 1)
    return {"width": 0.8 + 1.8 * share, "alpha": 0.35 + 0.55 * share}


def truncate(text, limit=90):
    text = str(text)
    return text if len(text) <= limit else text[:limit - 1] + "…"


def describe_dropped(dropped):
    return [{"id": rid, "reason": reason} for rid, reason in dropped]


def fraction_text(value):
    """Exact fractions as 'a/b (≈ 0.123)', so the UI never hides that UPGMA arithmetic is exact."""
    if hasattr(value, "denominator") and value.denominator != 1:
        return f"{value.numerator}/{value.denominator} (≈ {float(value):.4g})"
    return str(int(value)) if hasattr(value, "numerator") else str(value)
