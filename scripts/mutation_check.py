"""Inject deliberate bugs one at a time and confirm the test suite catches each; files are always restored."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (file, original text, mutated text, label); each original must occur exactly once
MUTANTS = [
    ("genomath/set_logic.py", "    for k in range(n):\n        added = 0", "    for k in range(n - 1):\n        added = 0",
     "Warshall skips the last pivot"),
    ("genomath/set_logic.py", "for c in succ[b] if (a, c) not in pairs", "for c in succ[b] if (c, a) not in pairs",
     "transitivity tests (c,a) instead of (a,c)"),
    ("genomath/set_logic.py", "shared_diff = sorted(i for i in ms & mt if s[i - 1] != t[i - 1])", "shared_diff = []",
     "d_H identity drops the shared-differing term"),
    ("genomath/set_logic.py", "            only_b.append(x)", "            only_a.append(x)",
     "B-only elements counted as A-only"),
    ("genomath/set_logic.py", 'lambda s: s[i - 1] == other[i - 1] != ref[i - 1]', 'lambda s: s[i - 1] != ref[i - 1]',
     "shared mutation ignores the base"),
    ("genomath/sequence_processing.py", "        seq = r.seq.upper()", "        seq = r.seq", "no uppercasing"),
    ("genomath/sequence_processing.py", "best = min(haplotypes, key=lambda h: (totals[h.id], h.id))",
     "best = min(haplotypes, key=lambda h: (totals[h.id], h.seq))", "medoid tie broken by sequence, not id"),
    ("genomath/_core.py", "        if a != b:\n            count += 1", "        if a != b and a != 'A':\n            count += 1",
     "mismatch count skips positions where S has A"),
    ("genomath/counting_recursion.py", "change = 3 * prev[k - 1]", "change = 2 * prev[k - 1]",
     "recurrence uses factor 2"),
    ("genomath/counting_recursion.py", "        c = c * (n - k + i) // i  #", "        c = c * (n - k + i + 1) // i  #",
     "binomial step off by one"),
    ("genomath/counting_recursion.py", "        for i in range(start, n - remaining + 1):",
     "        for i in range(start, n - remaining):", "generator never mutates the last allowed position"),
    ("genomath/counting_recursion.py", "        sign = -1 if j % 2 else 1", "        sign = 1",
     "inclusion-exclusion without alternating signs"),
    ("genomath/counting_recursion.py", "    forced = windows > holes", "    forced = windows >= holes",
     "pigeonhole uses >= instead of >"),
    ("genomath/counting_recursion.py", "            cur[p - 1] = t[p - 1]", "            cur[p - 1] = s[p - 1]",
     "ordering walk copies the source base"),
    ("genomath/counting_recursion.py", "    if count > cap:", "    if count > cap * 10:", "guard cap ten times too loose"),
    ("genomath/counting_recursion.py", "witness = (w, seen[w] + 1, i + 1)", "witness = (w, seen[w], i)",
     "repeated k-mer reported 0-based"),
    ("genomath/algebraic_structures.py", 'value, why = 4 ** ((n + 1) // 2), f"first', 'value, why = 4 ** (n // 2), f"first',
     "R-fixed count uses floor(n/2)"),
    ("genomath/algebraic_structures.py", 'value, why = 0, "odd n', 'value, why = 1, "odd n', "RC-fixed count 1 for odd n"),
    ("genomath/algebraic_structures.py", "if table[(table[(a, b)], c)] != table[(a, table[(b, c)])]",
     "if table[(table[(a, b)], c)] != table[(table[(a, b)], c)]", "associativity check compares a side with itself"),
    ("genomath/algebraic_structures.py", "    value = total // len(elements)", "    value = total // 2",
     "Burnside divides by 2 for every group"),
    ("genomath/algebraic_structures.py", "stabiliser = [g for g in elements if images[g] == s]",
     "stabiliser = [g for g in elements if images[g] != s]", "stabiliser collects the movers"),
    ("genomath/algebraic_structures.py", 'parity = "even" if transpositions % 2 == 0 else "odd"',
     'parity = "even" if len(cycles) % 2 == 0 else "odd"', "parity from number of cycles"),
    ("genomath/algebraic_structures.py", "    domain = list(words(n_max))\n    for c in candidates:",
     "    domain = list(words(1))\n    for c in candidates:", "function equality tested on lengths <= 1 only"),
    ("genomath/algebraic_structures.py", "    rep = min(TRANSFORMS[g](s) for g in elements)",
     "    rep = max(TRANSFORMS[g](s) for g in elements)", "canonical uses max"),
    ("genomath/algebraic_structures.py", "table = {(p, q): index[tuple(fns[p][pos[y]] for y in fns[q])]",
     "table = {(p, q): index[tuple(fns[q][pos[y]] for y in fns[p])]", "monoid table composes in the wrong order"),
    ("genomath/number_theory.py", "        old_s, s = s, old_s - q * s", "        old_s, s = s, old_s + q * s",
     "extended Euclid updates s with the wrong sign"),
    ("genomath/number_theory.py", "    inv = x % m\n", "    inv = x\n", "inverse not reduced mod m"),
    ("genomath/number_theory.py", '    s = "".join(reversed(digits))', '    s = "".join(digits)', "decode emits digits backwards"),
    ("genomath/number_theory.py", "        x = (4 * x) % mask + DIGIT[c]", "        x = 4 * x + DIGIT[c]",
     "sliding window never drops the old digit"),
    ("genomath/number_theory.py", "        rc = 4 * rc + (3 - d)", "        rc = 4 * rc + d", "RC on codes forgets to complement"),
    ("genomath/number_theory.py", "        t = a // 2 if q == 1 else (a + 1) // 2", "        t = (a + 1) // 2",
     "t_min ignores the delta = 2 case"),
    ("genomath/number_theory.py", "        value = value // p * (p - 1)", "        value = value // p * (p + 1)",
     "phi multiplies by p + 1"),
    ("genomath/number_theory.py", "    bound = -(-len(codes) // m)  # ceil(N/m)", "    bound = len(codes) // m  # ceil(N/m)",
     "pigeonhole bound uses floor"),
    ("genomath/number_theory.py", "        if e & 1:", "        if e & 2:", "square-and-multiply tests the wrong bit"),
    ("genomath/number_theory.py", "    sols = sorted(x0 + i * m2 for i in range(g))", "    sols = sorted(x0 + i * m2 for i in range(1))",
     "congruence returns one solution only"),
    ("genomath/mutation_graph.py", "keep = lambda x: 1 <= x <= k", "keep = lambda x: 1 <= x < k",
     "threshold graph uses d < k instead of d <= k"),
    ("genomath/mutation_graph.py", "parent[v], hops[v] = u, hops[u] + 1", "parent[v], hops[v] = u, hops[u] + 2",
     "BFS hop count off"),
    ("genomath/mutation_graph.py", "elif v != parent[u] and v not in finished:", "elif v != parent[u]:",
     "DFS counts finished vertices as back edges"),
    ("genomath/mutation_graph.py", "elif v != parent[u] and v in on_stack:", "elif v in on_stack:",
     "cycle search treats the parent edge as a cycle"),
    ("genomath/mutation_graph.py", "u = min(cand, key=lambda v: (dist[v], g.index[v]))",
     "u = max(cand, key=lambda v: (dist[v], g.index[v]))", "Dijkstra settles the farthest vertex"),
    ("genomath/mutation_graph.py", "dist[u] is not None and dist[v] is not None and dist[u] + w == dist[v]]",
     "dist[u] is not None and dist[v] is not None and dist[u] < dist[v]]", "shortest-path DAG keeps non-tight edges"),
    ("genomath/mutation_graph.py", "value = len(g.edges()) - len(g.vertices) + len(comps)",
     "value = len(g.edges()) - len(g.vertices) + 1", "cyclomatic number assumes one component"),
    ("genomath/mutation_graph.py", "colour[v] = 1 - colour[u]", "colour[v] = colour[u]", "bipartite colouring never alternates"),
    ("genomath/mutation_graph.py", "e = (u, v) if g.index[u] < g.index[v] else (v, u)", "e = (u, v)",
     "removed edges matched in one direction only"),
    ("genomath/mutation_graph.py", "values = sorted({x for row in d for x in row if x > 0})",
     "values = sorted({x for row in d for x in row if x > 0}, reverse=True)", "k* scans thresholds from the top"),
    ("genomath/tree_analysis.py", "        if uf.union(u, v):\n            tree.append",
     "        if uf.union(u, v) or True:\n            tree.append", "Kruskal accepts cycle-closing edges"),
    ("genomath/tree_analysis.py", "        if key(right[j]) < key(left[i]):", "        if key(right[j]) <= key(left[i]):",
     "merge sort loses stability"),
    ("genomath/tree_analysis.py", "if v not in in_tree and (key[v] is None or w < key[v]):",
     "if v not in in_tree and (key[v] is None or w <= key[v]):", "Prim updates on ties"),
    ("genomath/tree_analysis.py", "            total *= level_count", "            total += level_count",
     "MST count adds instead of multiplies"),
    ("genomath/tree_analysis.py", "            det = -det\n", "            det = det\n", "determinant ignores row swaps"),
    ("genomath/tree_analysis.py", "    w = max((e[2] for e in edges), default=0)", "    w = min((e[2] for e in edges), default=0)",
     "bottleneck takes the lightest edge"),
    ("genomath/phylogenetic_analysis.py", '(min if tie_break == "lexicographic" else max)',
     '(max if tie_break == "lexicographic" else min)', "UPGMA tie-break reversed"),
    ("genomath/phylogenetic_analysis.py",
     "(na * dist[frozenset((a, x))] + nb * dist[frozenset((b, x))]) / (na + nb)",
     "(dist[frozenset((a, x))] + dist[frozenset((b, x))]) / 2", "WPGMA update instead of UPGMA"),
    ("genomath/phylogenetic_analysis.py", "        height[node] = best / 2", "        height[node] = best",
     "UPGMA height not halved"),
    ("genomath/phylogenetic_analysis.py", "                if b != c:", "                if a != b:",
     "ultrametric test compares the two smallest"),
    ("genomath/tree_analysis.py", "        kids = sorted((v for v in adj[u] if v not in seen), key=key)",
     "        kids = sorted((v for v in adj[u] if v not in seen), key=key, reverse=True)", "children in reverse order"),
    ("genomath/tree_analysis.py", "        if 2 <= len(side) <= len(vertices) - 2:", "        if 1 <= len(side) <= len(vertices) - 1:",
     "MST splits include trivial ones"),
    ("genomath/algebraic_structures.py", 'TRANSITION_ELEMENT = "(A G)(C T)"', 'TRANSITION_ELEMENT = "(A T)(C G)"',
     "transitions identified with the wrong group element"),
    ("genomath/algebraic_structures.py", "    for b in others:", "    for b in ALPHABET:",
     "involutions counted twice"),
    ("genomath/algebraic_structures.py", 'key = "".join(sorted((x, y)))', 'key = x + y',
     "census keys ordered pairs, not unordered"),
    ("genomath/algebraic_structures.py", "            if x != y:", "            if x == y or x != y:",
     "census counts positions that agree"),
    ("genomath/algebraic_structures.py", '"ratio": (ts / tv if tv else None), "baseline_ratio": 0.5',
     '"ratio": (ts / tv if tv else None), "baseline_ratio": 2.0', "wrong combinatorial baseline"),
    # dropped as equivalent: swapping the composition order in base_cayley_table cannot change any entry, because the
    # base-level group is abelian (the suite verifies that), so the mutated table is literally identical
    ("app/presentation.py", '    for node in [m["node"] for m in result["merges"]]:',
     '    for node in reversed([m["node"] for m in result["merges"]]):',
     "dendrogram built from the last merge down"),
    ("app/presentation.py", '"omitted": list(omitted) if omitted else []', '"omitted": []',
     "capped trace steps lost in the UI"),
    ("app/presentation.py", '            cells.append("" if blank_diagonal and i == j else value)',
     "            cells.append(value)", "matrix diagonal never blanked"),
    ("app/presentation.py", "    if k <= 1:\n        return {", "    if k <= 0:\n        return {",
     "edge style wrong for the single-substitution graph"),
    ("app/presentation.py", "                force = dist * dist / k  # attraction along edges",
     "                force = -dist * dist / k  # attraction along edges",
     "layout attraction becomes repulsion"),
]


def run_tests(stop_early=False):
    # -x stops at the first failure: a mutant counts as caught as soon as any test fails
    args = [sys.executable, "-m", "pytest", "-q", "--no-header", "-p", "no:cacheprovider"] + (["-x"] if stop_early else [])
    r = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    return r.returncode, r.stdout.strip().splitlines()[-1]


def main():
    code, summary = run_tests()
    print(f"baseline: {summary}")
    if code:
        sys.exit("baseline must pass before mutants mean anything")
    survived = 0
    for rel, old, new, label in MUTANTS:
        p = ROOT / rel
        orig = p.read_text(encoding="utf-8")
        if orig.count(old) != 1:
            sys.exit(f"mutant '{label}': target text not found exactly once in {rel}")
        try:
            p.write_text(orig.replace(old, new), encoding="utf-8")
            code, summary = run_tests(stop_early=True)
        finally:
            p.write_text(orig, encoding="utf-8")
        survived += code == 0
        print(f"{'CAUGHT  ' if code else 'SURVIVED'} {label:48s} {summary}")
    print(f"{len(MUTANTS) - survived}/{len(MUTANTS)} mutants caught")


if __name__ == "__main__":
    main()
