"""Unit VI: UPGMA distance tree, kept separate from MST code."""
from fractions import Fraction

from .trace import Trace
from .tree_analysis import RootedTree


def upgma(ids, matrix, tie_break="lexicographic"):
    if tie_break not in ("lexicographic", "reverse"):
        raise ValueError("tie_break must be 'lexicographic' (decision D8) or 'reverse' (sensitivity check)")
    ids = list(ids)
    members = {x: (x,) for x in ids}
    dist = {frozenset((ids[i], ids[j])): Fraction(matrix[i][j]) for i in range(len(ids)) for j in range(i + 1, len(ids))}
    children, height, merges = {}, {x: Fraction(0) for x in ids}, []
    active = list(ids)
    step = 0
    while len(active) > 1:
        step += 1
        pairs = [(a, b) for i, a in enumerate(active) for b in active[i + 1:]]
        best = min(dist[frozenset(p)] for p in pairs)
        tied = [p for p in pairs if dist[frozenset(p)] == best]
        # D8: compare pairs by their sorted member tuples; "reverse" picks the last instead (to show sensitivity)
        keyed = [tuple(sorted((members[a], members[b]))) for a, b in tied]
        pick = (min if tie_break == "lexicographic" else max)(range(len(tied)), key=lambda i: keyed[i])
        a, b = tied[pick]
        if members[b] < members[a]:
            a, b = b, a
        node = f"U{step}"
        children[node] = [a, b]
        height[node] = best / 2
        members[node] = tuple(sorted(members[a] + members[b]))
        na, nb = len(members[a]), len(members[b])
        active = [x for x in active if x not in (a, b)]
        for x in active:  # size-weighted average = mean over all member pairs
            dist[frozenset((node, x))] = (na * dist[frozenset((a, x))] + nb * dist[frozenset((b, x))]) / (na + nb)
        active.append(node)
        merges.append({"step": step, "node": node, "left": a, "right": b, "distance": best, "height": height[node],
                       "tied_pairs": len(tied)})
    root = active[0]
    tree = RootedTree(root, children)
    result = {"tree": tree, "height": height, "members": members, "merges": merges, "root": root,
              "newick": newick(tree, root) + ";", "tie_steps": sum(1 for m in merges if m["tied_pairs"] > 1)}
    trace = Trace()
    trace.add_steps("algorithm", [f"step {m['step']}: {m['node']} = ({m['left']}, {m['right']}) at distance "
                                  f"{m['distance']}, height {m['height']}; {m['tied_pairs']} tied pair(s)" for m in merges])
    trace.add("calculation", f"Newick {result['newick']}", newick=result["newick"])
    trace.add("calculation", f"{result['tie_steps']} merge steps had ties, resolved by the {tie_break} rule")
    trace.add("limitation", "UPGMA assumes a molecular clock (all lineages evolve at the same rate, so distances are "
              "ultrametric); the tree is a clustering hypothesis, not proof of ancestry")
    return result, trace


def branch_lengths(result):
    # edge above each non-root node = parent height − child height (exact); leaf-to-root sums equal the root height
    out = {}
    for parent, kids in result["tree"].children.items():
        for c in kids:
            out[c] = result["height"][parent] - result["height"][c]
    return out, Trace().add("calculation", "branch length = parent height − child height", lengths=out)


def newick(tree, node, lengths=None):
    # lengths: optional {node: branch length above it}, e.g. from branch_lengths
    kids = tree.children.get(node)
    if not kids:
        return str(node)
    inner = ",".join(newick(tree, c, lengths) + (f":{float(lengths[c]):.4g}" if lengths else "") for c in kids)
    return f"({inner})"


def cophenetic(result, ids):
    # cophenetic distance of two leaves = 2 × height of their lowest common ancestor
    ids = list(ids)
    pos = {x: i for i, x in enumerate(ids)}
    m = [[Fraction(0)] * len(ids) for _ in ids]
    for merge in result["merges"]:
        left, right = result["members"][merge["left"]], result["members"][merge["right"]]
        for x in left:
            for y in right:
                m[pos[x]][pos[y]] = m[pos[y]][pos[x]] = 2 * merge["height"]
    return m, Trace().add("calculation", "cophenetic matrix: 2 × height of the lowest common ancestor")


def is_ultrametric(ids, matrix):
    ids = list(ids)
    n = len(ids)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                a, b, c = sorted((matrix[i][j], matrix[i][k], matrix[j][k]))
                if b != c:
                    return False, Trace().add("calculation", f"not ultrametric: ({ids[i]}, {ids[j]}, {ids[k]}) has "
                                              f"distances {a}, {b}, {c}; the two largest differ",
                                              witness=(ids[i], ids[j], ids[k]))
    return True, Trace().add("calculation", "ultrametric: in every triple the two largest distances are equal")


def fit_to_input(result, ids, matrix):
    coph, _ = cophenetic(result, ids)
    ids = list(ids)
    exact, worst, worst_pair = 0, Fraction(-1), None
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            dev = abs(coph[i][j] - matrix[i][j])
            exact += dev == 0
            if dev > worst:
                worst, worst_pair = dev, (ids[i], ids[j])
    pairs = len(ids) * (len(ids) - 1) // 2
    trace = Trace().add("calculation", f"cophenetic = d_H for {exact} of {pairs} pairs; largest deviation {worst} at "
                        f"{worst_pair}", exact=exact, pairs=pairs, worst=worst, worst_pair=worst_pair)
    trace.add("interpretation", "deviations measure how far the data are from the clock-like (ultrametric) shape "
              "UPGMA imposes")
    return {"exact": exact, "pairs": pairs, "max_deviation": worst, "worst_pair": worst_pair}, trace


def upgma_splits(result, ids):
    ids = list(ids)
    first, n = ids[0], len(ids)
    splits = set()
    for merge in result["merges"]:
        cluster = set(result["members"][merge["node"]])
        side = frozenset(x for x in ids if x not in cluster) if first in cluster else frozenset(cluster)
        if 2 <= len(side) <= n - 2:
            splits.add(side)
    return splits, Trace().add("calculation", f"{len(splits)} non-trivial bipartitions from UPGMA clusters")


def compare_splits(mst_splits, tree_splits):
    shared = mst_splits & tree_splits
    result = {"mst": len(mst_splits), "upgma": len(tree_splits), "shared": len(shared),
              "symmetric_difference": len(mst_splits ^ tree_splits), "shared_sets": [sorted(s) for s in shared]}
    trace = Trace().add("calculation", f"MST {len(mst_splits)} vs UPGMA {len(tree_splits)} non-trivial bipartitions; "
                        f"{len(shared)} shared; symmetric difference {result['symmetric_difference']}", **result)
    trace.add("interpretation", "the MST joins observed sequences directly; UPGMA builds hypothetical ancestors under a "
              "clock assumption, so they answer different questions and need not agree")
    return result, trace
