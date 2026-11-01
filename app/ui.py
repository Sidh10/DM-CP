"""Streamlit rendering helpers. All data shaping lives in app.presentation; this module only draws."""
import matplotlib
import streamlit as st

matplotlib.use("Agg")  # Streamlit renders figures server-side; no interactive backend is available
import matplotlib.pyplot as plt  # noqa: E402  (must follow the backend choice)

from app import presentation as pr

COLOURS = {
    "observation": "#5b6b7a",
    "calculation": "#1f6f8b",
    "algorithm": "#2e7d4f",
    "interpretation": "#8a5a1b",
    "limitation": "#9b2c2c",
}
FURNITURE = "#8a97a3"  # mid grey: legible against both the light and the dark Streamlit themes


def _new_figure(width, height):
    """Transparent background and mid-grey furniture, so a figure reads correctly under either theme."""
    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("none")
    for spine in ax.spines.values():
        spine.set_color(FURNITURE)
    ax.tick_params(colors=FURNITURE, labelsize=8)
    ax.xaxis.label.set_color(FURNITURE)
    ax.yaxis.label.set_color(FURNITURE)
    ax.title.set_color(FURNITURE)
    return fig, ax


def _render(fig):
    st.pyplot(fig, transparent=True)
    plt.close(fig)


def tag_legend():
    st.markdown("**How to read every result on this page**")
    for kind in pr.TAG_ORDER:
        label, description = pr.TAGS[kind]
        st.markdown(f"<span style='color:{COLOURS[kind]};font-weight:600'>■ {label}</span> — {description}",
                    unsafe_allow_html=True)


def tag(kind, text):
    label = pr.tag_label(kind)
    st.markdown(f"<span style='color:{COLOURS.get(kind, '#444')};font-weight:600'>[{label}]</span> {text}",
                unsafe_allow_html=True)


def result(kind, text, value=None):
    tag(kind, text)
    if value is not None:
        st.code(str(value), language=None)


def show_trace(trace, title="Show the calculation, step by step", expanded=False):
    """Every result on the dashboard is paired with this, so a reader can trace it back to a definition."""
    rows = pr.trace_rows(trace)
    summary = pr.trace_summary(trace)
    counts = ", ".join(f"{pr.tag_label(k)}: {v}" for k, v in summary.items())
    with st.expander(f"{title}  ({counts})", expanded=expanded):
        for row in rows:
            tag(row["kind"], row["text"])
            if row["omitted_count"]:
                with st.expander(f"the {row['omitted_count']} omitted steps"):
                    st.code("\n".join(row["omitted"]), language=None)


def table(header, body, caption=None):
    # every cell as text: display tables mix numbers with blanks and labels, which Arrow cannot type as one column
    st.table([{h: ("" if v is None else str(v)) for h, v in zip(header, row)} for row in body])
    if caption:
        st.caption(caption)


def matrix(ids, values, caption=None, blank_diagonal=True):
    header, body = pr.matrix_table(ids, values, blank_diagonal)
    table(header, body, caption)


def plot_graph(graph, k, highlight_path=None, title=None):
    edges = [(u, v) for u, v, _ in graph.edges()]
    pos = pr.spring_layout(graph.vertices, edges)
    fig, ax = _new_figure(7, 5.2)
    for u, v, w in graph.edges():
        style = pr.edge_style(w, k)
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#7a8b99",
                linewidth=style["width"], alpha=style["alpha"], zorder=1)
    if highlight_path and len(highlight_path) > 1:
        xs = [pos[v][0] for v in highlight_path]
        ys = [pos[v][1] for v in highlight_path]
        ax.plot(xs, ys, color="#c2410c", linewidth=2.8, alpha=0.95, zorder=2)
    ax.scatter([p[0] for p in pos.values()], [p[1] for p in pos.values()], s=320,
               color="#e8eef3", edgecolors="#33475b", linewidths=1.2, zorder=3)
    for v, (x, y) in pos.items():
        ax.annotate(str(v), (x, y), ha="center", va="center", fontsize=7, zorder=4, color="#1a2733")
    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    _render(fig)
    tag("limitation", "Vertex positions come from a force-directed layout chosen for readability. Distance on screen "
        "is not Hamming distance; only the edges and their weights carry information.")


def plot_dendrogram(result, title=None):
    layout = pr.dendrogram_layout(result)
    fig, ax = _new_figure(7, 4.6)
    for bar in layout["bars"]:
        ax.plot([bar["x_left"], bar["x_left"]], [bar["y_left"], bar["y"]], color="#33475b", linewidth=1.2)
        ax.plot([bar["x_right"], bar["x_right"]], [bar["y_right"], bar["y"]], color="#33475b", linewidth=1.2)
        ax.plot([bar["x_left"], bar["x_right"]], [bar["y"], bar["y"]], color="#33475b", linewidth=1.2)
    ax.set_xticks(range(len(layout["leaves"])))
    ax.set_xticklabels(layout["leaves"], rotation=90, fontsize=7)
    ax.set_ylabel("UPGMA height (d_H / 2)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    for bar_line in ax.get_lines():
        bar_line.set_color(FURNITURE)
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    _render(fig)


def plot_bars(labels, values, ylabel, title=None, logy=False):
    fig, ax = _new_figure(7, 3.4)
    ax.bar([str(x) for x in labels], values, color="#1f6f8b")
    ax.set_ylabel(ylabel, fontsize=9)
    if logy:
        ax.set_yscale("log")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    _render(fig)


def plot_sweep(rows):
    fig, ax = _new_figure(7, 3.4)
    ks = [r["k"] for r in rows]
    ax.plot(ks, [r["edges"] for r in rows], marker="o", color="#1f6f8b", label="edges")
    ax.set_xlabel("threshold k", fontsize=9)
    ax.set_ylabel("edges", fontsize=9, color="#1f6f8b")
    twin = ax.twinx()
    twin.patch.set_alpha(0.0)
    twin.plot(ks, [r["components"] for r in rows], marker="s", color="#9b2c2c", label="components")
    twin.set_ylabel("connected components", fontsize=9, color="#9b2c2c")
    twin.tick_params(colors=FURNITURE, labelsize=8)
    for spine in twin.spines.values():
        spine.set_color(FURNITURE)
    ax.spines[["top"]].set_visible(False)
    twin.spines[["top"]].set_visible(False)
    fig.tight_layout()
    _render(fig)


def section_header(number, title, question):
    st.header(f"{number}. {title}")
    st.caption(question)
