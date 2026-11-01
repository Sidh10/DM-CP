"""Unit V: Hamming-distance graphs, degrees, matrices, BFS/DFS, components, connectivity and Dijkstra."""
from itertools import combinations

from ._core import mismatch_count
from .trace import Trace

BRUTE_FORCE_LIMIT = 16  # exhaustive κ / λ search is exponential; refused above this many vertices


class Graph:
    def __init__(self, vertices, edges):
        self.vertices = list(vertices)
        self.index = {v: i for i, v in enumerate(self.vertices)}
        self.adj = {v: {} for v in self.vertices}
        for u, v, w in edges:
            if u == v:
                raise ValueError("loops are not allowed")
            self.adj[u][v] = w
            self.adj[v][u] = w
        for v in self.vertices:  # neighbours always in vertex order, so traversals are reproducible
            self.adj[v] = dict(sorted(self.adj[v].items(), key=lambda kv: self.index[kv[0]]))

    def edges(self):
        out = [(u, v, w) for u in self.vertices for v, w in self.adj[u].items() if self.index[u] < self.index[v]]
        return sorted(out, key=lambda e: (self.index[e[0]], self.index[e[1]]))

    def neighbours(self, v):
        return list(self.adj[v])

    def degree(self, v):
        return len(self.adj[v])


def distance_matrix(seqs):
    ids = list(seqs)
    matrix = [[mismatch_count(seqs[a], seqs[b]) for b in ids] for a in ids]
    trace = Trace().add("calculation", f"{len(ids)}×{len(ids)} Hamming distance matrix, "
                        f"{len(ids) * (len(ids) - 1) // 2} distinct pairs, O(N²·n)", ids=ids, matrix=matrix)
    return (ids, matrix), trace


def build_graph(seqs, k=1, mode="threshold"):
    (ids, d), trace = distance_matrix(seqs)
    g, t = graph_from_matrix(ids, d, k, mode)
    return g, trace.merge(t)


def graph_from_matrix(ids, d, k=1, mode="threshold"):
    # separated from build_graph so a threshold sweep computes the O(N²·n) distance matrix once
    trace = Trace()
    if mode == "threshold":
        rule = f"edge {{S,T}} iff 1 ≤ d_H(S,T) ≤ {k}, weight d_H"
        keep = lambda x: 1 <= x <= k
    elif mode == "single":
        rule = "edge {S,T} iff d_H(S,T) = 1 (single-substitution graph), weight 1"
        keep = lambda x: x == 1
    else:
        raise ValueError("mode must be 'threshold' or 'single'")
    edges = [(ids[i], ids[j], d[i][j]) for i in range(len(ids)) for j in range(i + 1, len(ids)) if keep(d[i][j])]
    g = Graph(ids, edges)
    trace.add("algorithm", f"{rule}: {len(edges)} edges on {len(ids)} vertices", rule=rule, edges=len(edges))
    trace.add("limitation", "vertices are observed haplotypes only; a missing edge may mean an unsampled intermediate")
    return g, trace


def degrees(g):
    deg = {v: g.degree(v) for v in g.vertices}
    total, m = sum(deg.values()), len(g.edges())
    trace = Trace()
    trace.add("calculation", f"degrees {deg}; handshake: Σ deg = {total} = 2·|E| = {2 * m}",
              degrees=deg, holds=total == 2 * m)
    trace.add("calculation", f"min degree δ = {min(deg.values(), default=0)}, max Δ = {max(deg.values(), default=0)}, "
              f"isolated {[v for v, x in deg.items() if x == 0]}")
    return deg, trace


def adjacency_matrix(g, weighted=False):
    m = [[(g.adj[u].get(v, 0) if weighted else int(v in g.adj[u])) for v in g.vertices] for u in g.vertices]
    trace = Trace().add("calculation", f"{'weighted ' if weighted else ''}adjacency matrix, symmetric with zero diagonal",
                        matrix=m)
    return m, trace


def incidence_matrix(g):
    edges = g.edges()
    m = [[int(v in (a, b)) for a, b, _ in edges] for v in g.vertices]
    trace = Trace().add("calculation", f"{len(g.vertices)}×{len(edges)} incidence matrix; each column has two 1s and "
                        "row sums are degrees", edges=[(a, b) for a, b, _ in edges], matrix=m)
    return (edges, m), trace


def bfs(g, source):
    order, parent, hops = [source], {source: None}, {source: 0}
    queue, head, steps = [source], 0, []
    while head < len(queue):  # list + head pointer is a FIFO queue without collections.deque
        u = queue[head]
        head += 1
        new = [v for v in g.neighbours(u) if v not in hops]
        for v in new:
            parent[v], hops[v] = u, hops[u] + 1
            queue.append(v)
            order.append(v)
        steps.append(f"dequeue {u}: enqueue {new or 'nothing new'}")
    trace = Trace().add_steps("algorithm", steps)
    trace.add("calculation", f"BFS order {order}; hop distances {hops}", order=order, hops=hops, parent=parent)
    trace.add("calculation", f"O(V + E) = O({len(g.vertices)} + {len(g.edges())})")
    return {"order": order, "parent": parent, "hops": hops}, trace


def dfs(g, source):
    pre, post, parent, back, finished = [], [], {source: None}, [], set()

    def visit(u):
        pre.append(u)
        for v in g.neighbours(u):
            if v not in parent:
                parent[v] = u
                visit(v)
            elif v != parent[u] and v not in finished:
                back.append((u, v))  # v is an ancestor still on the stack, so this non-tree edge closes a cycle
        post.append(u)
        finished.add(u)

    visit(source)
    trace = Trace()
    trace.add("algorithm", f"DFS preorder {pre}; postorder {post}", pre=pre, post=post)
    trace.add("calculation", f"tree edges {[(p, v) for v, p in parent.items() if p]}; back edges {back}", back=back)
    return {"pre": pre, "post": post, "parent": parent, "back_edges": back}, trace


def connected_components(g):
    seen, comps = set(), []
    for v in g.vertices:
        if v not in seen:
            r, _ = bfs(g, v)
            comp = sorted(r["order"], key=g.index.get)
            seen.update(comp)
            comps.append(comp)
    trace = Trace().add("calculation", f"{len(comps)} components, sizes {[len(c) for c in comps]}", components=comps)
    return comps, trace


def find_cycle(g):
    parent = {}
    for root in g.vertices:
        if root in parent:
            continue
        parent[root] = None
        stack = [(root, iter(g.neighbours(root)))]
        on_stack = {root}
        while stack:
            u, it = stack[-1]
            v = next(it, None)
            if v is None:
                stack.pop()
                on_stack.discard(u)
            elif v not in parent:
                parent[v] = u
                stack.append((v, iter(g.neighbours(v))))
                on_stack.add(v)
            elif v != parent[u] and v in on_stack:
                cycle = [u]
                while cycle[-1] != v:
                    cycle.append(parent[cycle[-1]])
                cycle.reverse()
                return cycle, Trace().add("algorithm", f"back edge {u}–{v} closes the cycle {cycle + [cycle[0]]}",
                                          cycle=cycle)
    return None, Trace().add("algorithm", "no cycle: the graph is a forest")


def cyclomatic_number(g):
    comps, _ = connected_components(g)
    value = len(g.edges()) - len(g.vertices) + len(comps)
    return value, Trace().add("calculation", f"|E| − |V| + c = {len(g.edges())} − {len(g.vertices)} + {len(comps)} = "
                              f"{value} independent cycles (0 iff forest)", value=value)


def is_bipartite(g):
    colour = {}
    for root in g.vertices:
        if root in colour:
            continue
        colour[root] = 0
        r, _ = bfs(g, root)
        for u in r["order"]:
            for v in g.neighbours(u):
                if v not in colour:
                    colour[v] = 1 - colour[u]
                elif colour[v] == colour[u]:
                    return False, Trace().add("algorithm", f"edge {u}–{v} joins two vertices of the same colour: "
                                              "an odd cycle exists, not bipartite", witness=(u, v))
    sides = [[v for v in g.vertices if colour[v] == c] for c in (0, 1)]
    return True, Trace().add("algorithm", f"2-colouring {sides[0]} | {sides[1]}", sides=sides)


def count_triangles(g):
    n = 0
    for u, v, _ in g.edges():
        for w in g.neighbours(u):
            if g.index[w] > g.index[v] and w in g.adj[v]:
                n += 1
    return n, Trace().add("calculation", f"{n} triangles")


def induced_subgraph(g, keep):
    keep = [v for v in g.vertices if v in set(keep)]
    edges = [(u, v, w) for u, v, w in g.edges() if u in keep and v in keep]
    return Graph(keep, edges), Trace().add("algorithm", f"induced subgraph on {keep}: {len(edges)} edges")


def complement(g, weight=lambda u, v: 1):
    edges = [(u, v, weight(u, v)) for u, v in combinations(g.vertices, 2) if v not in g.adj[u]]
    return Graph(g.vertices, edges), Trace().add("algorithm", f"complement: {len(edges)} edges")


def edge_set_operation(g, h, op):
    if g.vertices != h.vertices:
        raise ValueError("graphs must share the same vertex list")
    eg, eh = {(u, v): w for u, v, w in g.edges()}, {(u, v): w for u, v, w in h.edges()}
    keys = {"union": eg.keys() | eh.keys(), "intersection": eg.keys() & eh.keys(), "difference": eg.keys() - eh.keys()}[op]
    edges = [(u, v, eg.get((u, v), eh.get((u, v)))) for u, v in keys]
    return Graph(g.vertices, edges), Trace().add("algorithm", f"edge {op}: {len(edges)} edges")


def _connected_without(g, removed_vertices=(), removed_edges=()):
    removed_vertices, removed_edges = set(removed_vertices), set(removed_edges)
    rest = [v for v in g.vertices if v not in removed_vertices]
    if len(rest) <= 1:
        return True
    seen, stack = {rest[0]}, [rest[0]]
    while stack:
        u = stack.pop()
        for v in g.neighbours(u):
            e = (u, v) if g.index[u] < g.index[v] else (v, u)
            if v not in removed_vertices and e not in removed_edges and v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == len(rest)


def vertex_connectivity(g):
    n = len(g.vertices)
    if n > BRUTE_FORCE_LIMIT:
        raise ValueError(f"exhaustive κ is refused above {BRUTE_FORCE_LIMIT} vertices")
    if not _connected_without(g):
        return 0, Trace().add("calculation", "κ = 0: the graph is already disconnected")
    complete = all(len(g.adj[v]) == n - 1 for v in g.vertices)
    if complete:
        return n - 1, Trace().add("calculation", f"complete graph: κ = n − 1 = {n - 1} by convention")
    for size in range(1, n - 1):
        for cut in combinations(g.vertices, size):
            if not _connected_without(g, removed_vertices=cut):
                return size, Trace().add("algorithm", f"κ = {size}: removing {list(cut)} disconnects the graph "
                                         "(smallest such set, found by exhaustive search)", cut=list(cut))
    return n - 1, Trace().add("calculation", f"κ = {n - 1}")


def edge_connectivity(g):
    edges = [(u, v) for u, v, _ in g.edges()]
    if len(g.vertices) > BRUTE_FORCE_LIMIT:
        raise ValueError(f"exhaustive λ is refused above {BRUTE_FORCE_LIMIT} vertices")
    if len(g.vertices) <= 1:
        return 0, Trace().add("calculation", "λ = 0 for a single vertex")
    if not _connected_without(g):
        return 0, Trace().add("calculation", "λ = 0: the graph is already disconnected")
    delta = min(g.degree(v) for v in g.vertices)
    for size in range(1, delta + 1):  # λ ≤ δ, so larger cuts never need checking
        for cut in combinations(edges, size):
            if not _connected_without(g, removed_edges=cut):
                return size, Trace().add("algorithm", f"λ = {size}: removing {list(cut)} disconnects the graph",
                                         cut=list(cut))
    return delta, Trace().add("calculation", f"λ = δ = {delta}")


def dijkstra(g, source):
    # O(V²) array version: each round scans for the closest unsettled vertex (ties: earliest in vertex order)
    dist = {v: None for v in g.vertices}
    prev = {v: None for v in g.vertices}
    dist[source] = 0
    settled, steps = [], []
    while True:
        cand = [v for v in g.vertices if v not in settled and dist[v] is not None]
        if not cand:
            break
        u = min(cand, key=lambda v: (dist[v], g.index[v]))
        settled.append(u)
        relaxed = []
        for v, w in g.adj[u].items():
            if v not in settled and (dist[v] is None or dist[u] + w < dist[v]):
                dist[v], prev[v] = dist[u] + w, u
                relaxed.append(f"{v}:{dist[v]}")
        steps.append(f"settle {u} at {dist[u]}; relax {relaxed or 'none'}")
    trace = Trace().add_steps("algorithm", steps)
    unreachable = [v for v in g.vertices if dist[v] is None]
    trace.add("calculation", f"distances from {source}: {dist}; unreachable {unreachable}", dist=dist,
              unreachable=unreachable)
    trace.add("calculation", f"O(V²) = O({len(g.vertices)}²) with an array scan; weights d_H ≥ 1 are non-negative")
    return {"dist": dist, "prev": prev, "settled": settled}, trace


def path_from(prev, target):
    path = [target]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return path[::-1]


def all_shortest_paths(g, source, target, dist=None, limit=100):
    if dist is None:
        dist = dijkstra(g, source)[0]["dist"]
    if dist[target] is None:
        return [], Trace().add("calculation", f"{target} is unreachable from {source}")
    preds = {v: [u for u, w in g.adj[v].items() if dist[u] is not None and dist[v] is not None and dist[u] + w == dist[v]]
             for v in g.vertices}
    paths = []

    def back(v, suffix):
        if len(paths) >= limit:
            return
        if v == source:
            paths.append([v] + suffix)
            return
        for u in preds[v]:
            back(u, [v] + suffix)

    back(target, [])
    paths.sort(key=lambda p: [g.index[v] for v in p])
    trace = Trace().add("calculation", f"{len(paths)} shortest path(s) of cost {dist[target]} from {source} to {target}",
                        paths=paths)
    return paths, trace


def path_versus_direct(g, seqs, source, target):
    r, trace = dijkstra(g, source)
    cost = r["dist"][target]
    direct = mismatch_count(seqs[source], seqs[target])
    if cost is None:
        trace.add("calculation", f"{target} is unreachable from {source} in this graph (direct d_H = {direct})")
        return {"cost": None, "direct": direct, "detour": None, "paths": []}, trace
    paths, _ = all_shortest_paths(g, source, target, r["dist"])
    trace.add("calculation", f"path cost {cost} vs direct d_H {direct}: detour {cost - direct}", cost=cost,
              direct=direct)
    trace.add("interpretation", "cost ≥ d_H by the triangle inequality; a positive detour means every observed route "
              "changes some position more than once")
    trace.add("limitation", "a path through observed haplotypes is not an inferred evolutionary trajectory")
    return {"cost": cost, "direct": direct, "detour": cost - direct, "paths": paths}, trace


def threshold_sweep(seqs, ks):
    rows = []
    (ids, d), _ = distance_matrix(seqs)
    for k in ks:
        g, _ = graph_from_matrix(ids, d, k)
        comps, _ = connected_components(g)
        deg = [g.degree(v) for v in g.vertices]
        rows.append({"k": k, "edges": len(g.edges()), "components": len(comps),
                     "largest": max(len(c) for c in comps), "isolated": sum(1 for x in deg if x == 0),
                     "max_degree": max(deg), "cyclomatic": len(g.edges()) - len(g.vertices) + len(comps)})
    trace = Trace().add("calculation", f"swept {len(rows)} thresholds", rows=rows)
    trace.add("interpretation", "G_j ⊆ G_k for j ≤ k, so edges never decrease and components never increase")
    return rows, trace


def smallest_connecting_threshold(seqs):
    (ids, d), _ = distance_matrix(seqs)
    values = sorted({x for row in d for x in row if x > 0})
    for k in values:
        g, _ = graph_from_matrix(ids, d, k)
        if len(connected_components(g)[0]) == 1:
            return k, Trace().add("calculation", f"k* = {k}: the smallest threshold at which G_k is connected "
                                  f"(checked the distinct distances {values[:values.index(k) + 1]})", value=k)
    return 0, Trace().add("calculation", "k* = 0: at most one vertex")
