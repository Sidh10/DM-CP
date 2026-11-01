"""Unit VI: Prim and Kruskal MSTs, rooted-tree construction, statistics and traversals."""
from fractions import Fraction

from .trace import Trace


def merge_sort(items, key):
    # hand-written so Kruskal's O(E log E) ordering is visible; stable, so equal keys keep their input order
    if len(items) <= 1:
        return list(items)
    mid = len(items) // 2
    left, right = merge_sort(items[:mid], key), merge_sort(items[mid:], key)
    out, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        if key(right[j]) < key(left[i]):
            out.append(right[j])
            j += 1
        else:
            out.append(left[i])
            i += 1
    return out + left[i:] + right[j:]


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:  # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra  # union by rank keeps trees shallow
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True


def kruskal(g):
    order = merge_sort(g.edges(), key=lambda e: (e[2], g.index[e[0]], g.index[e[1]]))
    uf = UnionFind(g.vertices)
    tree, steps, examined = [], [], 0
    target = len(g.vertices) - 1
    for u, v, w in order:
        if len(tree) == target:
            break
        examined += 1
        if uf.union(u, v):
            tree.append((u, v, w))
            steps.append(f"accept {u}–{v} (w = {w})")
        else:
            steps.append(f"reject {u}–{v} (w = {w}): same component, would close a cycle")
    total = sum(w for _, _, w in tree)
    trace = Trace().add("algorithm", f"sorted {len(order)} edges by (weight, endpoints) with merge sort")
    trace.add_steps("algorithm", steps)
    trace.add("calculation", f"examined {examined}, accepted {len(tree)}, rejected {examined - len(tree)}; "
              f"total weight {total}", examined=examined, accepted=len(tree), total=total)
    if len(tree) < target:
        trace.add("limitation", "graph is disconnected: the result is a minimum spanning forest")
    return {"edges": tree, "total": total, "examined": examined}, trace


def prim(g, root=None):
    # O(V²) array version: key[v] = cheapest known edge from v into the tree; strict improvement only
    root = g.vertices[0] if root is None else root
    key = {v: None for v in g.vertices}
    via = {v: None for v in g.vertices}
    key[root] = 0
    in_tree, order, tree, steps = set(), [], [], []
    while True:
        cand = [v for v in g.vertices if v not in in_tree and key[v] is not None]
        if not cand:
            break
        u = min(cand, key=lambda v: (key[v], g.index[v]))
        in_tree.add(u)
        order.append(u)
        if via[u] is not None:
            tree.append((via[u], u, key[u]))
        improved = []
        for v, w in g.adj[u].items():
            if v not in in_tree and (key[v] is None or w < key[v]):
                key[v], via[v] = w, u
                improved.append(f"{v}:{w}")
        steps.append(f"add {u}" + (f" via {via[u]} (w = {key[u]})" if via[u] else " (root)")
                     + f"; improved {improved or 'none'}")
    total = sum(w for _, _, w in tree)
    trace = Trace().add_steps("algorithm", steps)
    trace.add("calculation", f"Prim from {root}: order {order}; total weight {total}", order=order, total=total)
    if len(order) < len(g.vertices):
        trace.add("limitation", f"only the component of {root} is spanned ({len(order)} of {len(g.vertices)} vertices)")
    return {"edges": tree, "total": total, "order": order}, trace


def validate_spanning_tree(vertices, edges, g=None):
    uf = UnionFind(vertices)
    acyclic = all(uf.union(u, v) for u, v, _ in edges)
    roots = {uf.find(v) for v in vertices}
    in_graph = g is None or all(v in g.adj[u] and g.adj[u][v] == w for u, v, w in edges)
    report = {"edge_count": len(edges), "n_minus_1": len(edges) == len(vertices) - 1, "acyclic": acyclic,
              "spanning": len(roots) == 1, "edges_in_graph": in_graph}
    report["is_spanning_tree"] = all(report[k] for k in ("n_minus_1", "acyclic", "spanning", "edges_in_graph"))
    trace = Trace().add("algorithm", f"|E| = {len(edges)} (need {len(vertices) - 1}); acyclic {acyclic}; "
                        f"connects all vertices {len(roots) == 1}; uses graph edges {in_graph}", **report)
    return report, trace


def _determinant(m):
    # exact Gaussian elimination over the rationals
    m = [[Fraction(x) for x in row] for row in m]
    n, det = len(m), Fraction(1)
    for c in range(n):
        pivot = next((r for r in range(c, n) if m[r][c] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != c:
            m[c], m[pivot] = m[pivot], m[c]
            det = -det
        det *= m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            if f:
                for k in range(c, n):
                    m[r][k] -= f * m[c][k]
    return det


def kirchhoff_count(vertices, multi_edges):
    # matrix-tree theorem: spanning trees = any cofactor of the Laplacian (multi-edges allowed)
    vertices = list(vertices)
    if len(vertices) <= 1:
        return 1, Trace().add("calculation", "a single vertex has exactly 1 spanning tree")
    idx = {v: i for i, v in enumerate(vertices)}
    lap = [[0] * len(vertices) for _ in vertices]
    for u, v in multi_edges:
        i, j = idx[u], idx[v]
        lap[i][i] += 1
        lap[j][j] += 1
        lap[i][j] -= 1
        lap[j][i] -= 1
    value = int(_determinant([row[1:] for row in lap[1:]]))
    trace = Trace().add("calculation", f"Kirchhoff: det of the Laplacian with row/column 1 removed = {value}",
                        laplacian=lap, value=value)
    return value, trace


def count_minimum_spanning_trees(g):
    # process weight levels; contract the components formed by lighter edges; multiply spanning-tree counts
    uf = UnionFind(g.vertices)
    levels = {}
    for u, v, w in g.edges():
        levels.setdefault(w, []).append((u, v))
    total, steps = 1, []
    for w in merge_sort(list(levels), key=lambda x: x):
        crossing = [(uf.find(u), uf.find(v)) for u, v in levels[w] if uf.find(u) != uf.find(v)]
        if crossing:
            supers = sorted({x for e in crossing for x in e}, key=g.index.get)
            comp_uf = UnionFind(supers)
            for a, b in crossing:
                comp_uf.union(a, b)
            groups = {}
            for s in supers:
                groups.setdefault(comp_uf.find(s), []).append(s)
            level_count = 1
            for members in groups.values():
                edges = [(a, b) for a, b in crossing if comp_uf.find(a) == comp_uf.find(members[0])]
                level_count *= kirchhoff_count(members, edges)[0]
            total *= level_count
            steps.append(f"weight {w}: {len(crossing)} crossing edges between {len(supers)} super-vertices -> "
                         f"{level_count} ways")
        for u, v in levels[w]:
            uf.union(u, v)
    trace = Trace().add_steps("algorithm", steps)
    trace.add("calculation", f"number of minimum spanning trees = {total}", value=total)
    trace.add("interpretation", "every MST has the same multiset of edge weights; ties make the tree itself non-unique")
    return total, trace


def bottleneck(edges):
    w = max((e[2] for e in edges), default=0)
    return w, Trace().add("calculation", f"largest MST edge = {w}: the smallest threshold k at which G_k is connected "
                          "(minimax property)", value=w)


def cayley_count(n):
    value = n ** (n - 2) if n >= 2 else 1
    return value, Trace().add("calculation", f"labelled trees on {n} vertices: n^(n−2) = {value}", value=value)


def double_factorial(m):
    out = 1
    while m > 1:
        out *= m
        m -= 2
    return out


def rooted_binary_topologies(n):
    value = double_factorial(2 * n - 3) if n >= 2 else 1
    return value, Trace().add("calculation", f"rooted binary topologies on {n} labelled leaves: (2n−3)!! = {value}",
                              value=value)


class RootedTree:
    def __init__(self, root, children, order=None):
        self.root = root
        self.children = {v: list(c) for v, c in children.items()}
        self.parent = {root: None}
        for v, cs in self.children.items():
            for c in cs:
                self.parent[c] = v
        self.nodes = [root] + [v for v in self.parent if v != root]

    def is_binary(self):
        return all(len(self.children.get(v, [])) in (0, 2) for v in self.nodes)


def root_tree(vertices, edges, root, order_key=None):
    adj = {v: [] for v in vertices}
    for u, v, _ in edges:
        adj[u].append(v)
        adj[v].append(u)
    key = order_key or {v: i for i, v in enumerate(vertices)}.get
    children, seen, queue, head = {}, {root}, [root], 0
    while head < len(queue):
        u = queue[head]
        head += 1
        kids = sorted((v for v in adj[u] if v not in seen), key=key)
        children[u] = kids
        seen.update(kids)
        queue.extend(kids)
    if len(seen) != len(vertices):
        raise ValueError("edges do not form a spanning tree on these vertices")
    tree = RootedTree(root, children)
    trace = Trace().add("algorithm", f"rooted at {root} (an arbitrary choice for an MST); parent/child by BFS",
                        root=root, children=children)
    trace.add("limitation", "the root is a presentation choice: depth, height and leaves change with it, the edges do not")
    return tree, trace


def tree_statistics(tree):
    depth = {tree.root: 0}
    for v in level_order(tree)[0]:
        for c in tree.children.get(v, []):
            depth[c] = depth[v] + 1
    leaves = [v for v in tree.nodes if not tree.children.get(v)]
    internal = [v for v in tree.nodes if tree.children.get(v)]
    stats = {"depth": depth, "height": max(depth.values()), "leaves": leaves, "internal": internal,
             "edges": len(tree.nodes) - 1}
    trace = Trace().add("calculation", f"height {stats['height']}; {len(leaves)} leaves {leaves}; "
                        f"{len(internal)} internal {internal}; edges = nodes − 1 = {stats['edges']}", **stats)
    return stats, trace


def preorder(tree):
    out, stack = [], [tree.root]
    while stack:
        v = stack.pop()
        out.append(v)
        stack.extend(reversed(tree.children.get(v, [])))
    return out, Trace().add("algorithm", f"preorder (node, then children left to right): {out}")


def postorder(tree):
    out = []

    def visit(v):
        for c in tree.children.get(v, []):
            visit(c)
        out.append(v)

    visit(tree.root)
    return out, Trace().add("algorithm", f"postorder (children, then node): {out}")


def level_order(tree):
    out, head = [tree.root], 0
    while head < len(out):
        out.extend(tree.children.get(out[head], []))
        head += 1
    return out, Trace().add("algorithm", f"level order: {out}")


def inorder(tree):
    if not tree.is_binary():
        raise ValueError("in-order traversal needs a binary tree (every node has 0 or 2 children)")
    out = []

    def visit(v):
        kids = tree.children.get(v, [])
        if kids:
            visit(kids[0])
        out.append(v)
        if kids:
            visit(kids[1])

    visit(tree.root)
    return out, Trace().add("algorithm", f"in-order (left, node, right): {out}")


def edge_splits(vertices, edges):
    # each tree edge splits the vertex set in two; keep the side without the first vertex, ignore single-vertex sides
    first = vertices[0]
    splits = set()
    for cut in edges:
        rest = [e for e in edges if e is not cut]
        uf = UnionFind(vertices)
        for u, v, _ in rest:
            uf.union(u, v)
        side = frozenset(v for v in vertices if uf.find(v) != uf.find(first))
        if 2 <= len(side) <= len(vertices) - 2:
            splits.add(side)
    return splits, Trace().add("calculation", f"{len(splits)} non-trivial bipartitions from {len(edges)} edges")
