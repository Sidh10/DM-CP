"""Brute-force cross-check of docs/synthetic_expected.md; deliberately shares no code with genomath."""
import itertools, sys
from math import comb
from pathlib import Path
recs, cur = {}, None
for line in Path(sys.argv[1]).read_text().splitlines():
    if line.startswith('>'): cur = line[1:].split()[0]; recs[cur] = ''
    else: recs[cur] += line.strip()
R = recs['syn01']; n = len(R)
d = lambda a, b: sum(x != y for x, y in zip(a, b))
comp = str.maketrans('ACGT', 'TGCA'); rc = lambda s: s.translate(comp)[::-1]
enc = lambda s: sum('ACGT'.index(c) * 4**(len(s)-1-i) for i, c in enumerate(s))
print('RC(syn01)==syn01:', rc(R) == R, '| RC(syn03)==syn08:', rc(recs['syn03']) == recs['syn08'])
print('palindromes:', [k for k, s in recs.items() if rc(s) == s])
for k, s in recs.items():
    M = {(i+1, s[i]) for i in range(n) if s[i] != R[i]}
    print(k, s, 'M=', sorted(M), 'enc=', enc(s), 'enc(comp)=', enc(s.translate(comp)), '4^n-1-enc=', 4**n-1-enc(s),
          'canon=', min(enc(s), enc(rc(s))))
haps = {}
for k, s in recs.items(): haps.setdefault(s, k)
H = sorted(haps.values()); seq = {k: recs[k] for k in H}
print('haplotypes:', H)
print('     ' + ' '.join(h[-2:] for h in H))
for a in H: print(a, ' '.join(f'{d(seq[a], seq[b]):2d}' for b in H))
E = [(d(seq[a], seq[b]), a, b) for a, b in itertools.combinations(H, 2)]
from collections import Counter
print('distance histogram:', sorted(Counter(w for w, _, _ in E).items()))
def comps(edges):
    p = {h: h for h in H}
    def f(x):
        while p[x] != x: x = p[x]
        return x
    for _, a, b in edges: p[f(a)] = f(b)
    g = {}
    for h in H: g.setdefault(f(h), []).append(h)
    return sorted(g.values())
for k in range(0, 10):
    ek = [e for e in E if e[0] <= k]
    print(f'k={k}: |E|={len(ek)} components={len(comps(ek))} {comps(ek)}')
# transitivity counterexample for k=1
print('d(01,02),d(02,03),d(01,03)=', d(seq['syn01'], seq['syn02']), d(seq['syn02'], seq['syn03']), d(seq['syn01'], seq['syn03']))
# MST weight and number of MSTs by exhaustive search over 8-edge subsets
best, count = None, 0
W = [e[0] for e in E]
for idx in itertools.combinations(range(len(E)), len(H) - 1):
    w = sum(W[i] for i in idx)
    if best is not None and w > best: continue
    if len(comps([E[i] for i in idx])) != 1: continue
    if best is None or w < best: best, count = w, 1
    else: count += 1
print('MST weight =', best, ' number of distinct MSTs =', count)
# Dijkstra by brute force over simple paths in the k=2 graph
def paths(a, b, k):
    adj = {h: [x for x in H if x != h and d(seq[h], seq[x]) <= k] for h in H}
    out = []
    def go(p):
        if p[-1] == b: out.append(p); return
        for x in adj[p[-1]]:
            if x not in p: go(p + [x])
    go([a]); return out
for a, b in [('syn04', 'syn08'), ('syn07', 'syn03'), ('syn07', 'syn08')]:
    ps = paths(a, b, 2); cost = lambda p: sum(d(seq[x], seq[y]) for x, y in zip(p, p[1:]))
    m = min(map(cost, ps)) if ps else None
    print(f'k=2 {a}->{b}: direct d_H={d(seq[a], seq[b])} shortest cost={m} optimal paths={[p for p in ps if cost(p) == m]}')
print('k=1 syn07 reachable from syn03:', bool(paths('syn03', 'syn07', 1)))
for m in [16, 13, 17, 8]:
    print(f'mod {m}:', {h: enc(seq[h]) % m for h in H})
print('4^12=', 4**12, 'M(12,k)=', [comb(12, k)*3**k for k in range(4)], 'ball r<=2:', sum(comb(12, k)*3**k for k in range(3)))
for nn in (4, 12):
    print(f'IE all-four n={nn}:', 4**nn - 4*3**nn + 6*2**nn - 4, ' brute:' , sum(1 for t in itertools.product('ACGT', repeat=nn) if set(t) == set('ACGT')) if nn <= 8 else '')
for nn in (1, 2, 3, 4, 12):
    full = list(itertools.product('ACGT', repeat=nn)) if nn <= 8 else None
    fixRC = 4**(nn//2) if nn % 2 == 0 else 0
    z2 = (4**nn + fixRC)//2; klein = (4**nn + 4**((nn+1)//2) + 0 + fixRC)//4
    if full:
        S = [''.join(t) for t in full]
        z2b = len({min(s, rc(s)) for s in S})
        kb = len({min(s, s[::-1], s.translate(comp), rc(s)) for s in S})
        print(f'n={nn}: RC-fixed={fixRC} (brute {sum(rc(s)==s for s in S)}), Z2 orbits={z2} (brute {z2b}), Klein orbits={klein} (brute {kb})')
    else:
        print(f'n={nn}: RC-fixed={fixRC}, Z2 orbits={z2}, Klein orbits={klein}')
print('9^7=', 9**7, ' 15!!=', 15*13*11*9*7*5*3)
x = enc('GAT'); print('GAT=', x, 'mod 7,11:', x % 7, x % 11, ' ACG=', enc('ACG'))
