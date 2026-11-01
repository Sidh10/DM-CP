"""Unit IV: base-4 encoding, modular bucketing, Euclid/extended Euclid, inverses, CRT and canonical k-mers."""
from ._core import ALPHABET, SequenceError
from .algebraic_structures import burnside
from .trace import Trace

DIGIT = {b: i for i, b in enumerate(ALPHABET)}  # A=0, C=1, G=2, T=3


def _check(s):
    bad = sorted(set(s) - set(ALPHABET))
    if bad:
        raise SequenceError(f"encoding is defined on {ALPHABET} only, got {bad}")


def _enc(s):
    x = 0
    for c in s:
        x = 4 * x + DIGIT[c]
    return x


def _rc_code(x, n):
    rc = 0
    for _ in range(n):
        x, d = x // 4, x % 4
        rc = 4 * rc + (3 - d)
    return rc


def division(a, b):
    if b <= 0:
        raise ValueError("the division algorithm here takes b > 0")
    q, r = a // b, a % b  # floor division gives 0 <= r < b also for negative a
    trace = Trace().add("calculation", f"{a} = {b}·{q} + {r}, 0 ≤ {r} < {b}", q=q, r=r)
    return (q, r), trace


def divides(a, b):
    holds = (b == 0) if a == 0 else b % a == 0
    return holds, Trace().add("calculation", f"{a} | {b}: {holds}")


def gcd(a, b):
    a, b = abs(a), abs(b)
    trace = Trace()
    while b:
        q, r = a // b, a % b
        trace.add("algorithm", f"{a} = {q}·{b} + {r}", a=a, b=b, q=q, r=r)
        a, b = b, r
    trace.add("calculation", f"gcd = {a}" + (" (convention: gcd(0,0) = 0)" if a == 0 else ""), value=a)
    return a, trace


def extended_gcd(a, b):
    # iterative form keeps the invariant  old_r = a·old_s + b·old_t  at every step
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    trace = Trace()
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
        trace.add("algorithm", f"q = {q}: r = {old_r}, s = {old_s}, t = {old_t}  ({a}·{old_s} + {b}·{old_t} = {old_r})")
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    trace.add("calculation", f"gcd({a},{b}) = {old_r} = {a}·{old_s} + {b}·{old_t}", g=old_r, x=old_s, y=old_t)
    return (old_r, old_s, old_t), trace


def mod_inverse(a, m):
    if m <= 1:
        raise ValueError("modulus must be > 1")
    (g, x, _), trace = extended_gcd(a % m, m)
    if g != 1:
        trace.add("limitation", f"gcd({a},{m}) = {g} ≠ 1, so {a} has no inverse mod {m}")
        raise ValueError(f"{a} has no inverse modulo {m} (gcd {g})")
    inv = x % m
    trace.add("calculation", f"{a}^(-1) ≡ {inv} (mod {m}); check {a}·{inv} = {a * inv} ≡ {(a * inv) % m}", value=inv)
    return inv, trace


def solve_linear_congruence(a, b, m):
    g, trace = gcd(a, m)
    if b % g:
        trace.add("calculation", f"{g} ∤ {b}: {a}x ≡ {b} (mod {m}) has no solution")
        return [], trace
    a2, b2, m2 = a // g, b // g, m // g
    x0 = (b2 * mod_inverse(a2, m2)[0]) % m2 if m2 > 1 else 0
    sols = sorted(x0 + i * m2 for i in range(g))
    trace.add("calculation", f"divide by g = {g}: {a2}x ≡ {b2} (mod {m2}) gives x ≡ {x0}; "
              f"{g} solutions mod {m}: {sols}", solutions=sols)
    return sols, trace


def crt(residues, moduli):
    trace = Trace()
    for i in range(len(moduli)):
        for j in range(i + 1, len(moduli)):
            g, _ = gcd(moduli[i], moduli[j])
            if g != 1:
                trace.add("limitation", f"moduli {moduli[i]} and {moduli[j]} share factor {g}; CRT needs pairwise coprime")
                raise ValueError(f"moduli {moduli[i]} and {moduli[j]} are not coprime (gcd {g})")
    M = 1
    for m in moduli:
        M *= m
    x = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        inv, _ = mod_inverse(Mi, m)
        x += r * Mi * inv
        trace.add("calculation", f"r = {r} (mod {m}): M_i = {Mi}, M_i^(-1) ≡ {inv}, term {r}·{Mi}·{inv} = {r * Mi * inv}")
    x %= M
    trace.add("calculation", f"x ≡ {x} (mod {M}); unique in 0..{M - 1}", value=x, M=M)
    return x, trace


def factorize(n):
    if n < 1:
        raise ValueError("factorize needs n >= 1")
    factors, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def euler_phi(n):
    f = factorize(n)
    value = n
    for p in f:
        value = value // p * (p - 1)
    trace = Trace().add("calculation", f"φ({n}) = {n}·Π(1 − 1/p) over p in {sorted(f)} = {value}", value=value)
    return value, trace


def mod_pow(base, exp, m):
    # square-and-multiply: O(log exp) multiplications
    result, b, e = 1 % m, base % m, exp
    while e:
        if e & 1:
            result = result * b % m
        b = b * b % m
        e >>= 1
    return result, Trace().add("calculation", f"{base}^{exp} ≡ {result} (mod {m})", value=result)


def multiplicative_order(a, m):
    if gcd(a, m)[0] != 1:
        raise ValueError(f"order of {a} mod {m} is undefined: gcd ≠ 1")
    k, x = 1, a % m
    while x != 1 % m:
        x = x * a % m
        k += 1
    phi, _ = euler_phi(m)
    trace = Trace().add("calculation", f"ord_{m}({a}) = {k}; it divides φ({m}) = {phi}: {phi % k == 0}",
                        value=k, divides_phi=phi % k == 0)
    trace.add("limitation", "pedagogical: orders explain when positional weights 4^t repeat mod m")
    return k, trace


def encode(s):
    _check(s)
    x, steps = 0, []
    for i, c in enumerate(s, start=1):
        x = 4 * x + DIGIT[c]
        steps.append(f"Horner step {i}: 4·(previous) + {c}={DIGIT[c]} -> {x}")
    trace = Trace().add_steps("calculation", steps)
    trace.add("calculation", f"enc({s}) = {x}", value=x)
    return x, trace


def decode(x, n):
    if not 0 <= x < 4 ** n:
        raise ValueError(f"{x} is outside 0..4^{n} - 1")
    digits, steps = [], []
    for _ in range(n):
        q, r = x // 4, x % 4
        steps.append(f"{x} = 4·{q} + {r} -> {ALPHABET[r]}")
        digits.append(ALPHABET[r])
        x = q
    trace = Trace().add_steps("calculation", steps)
    s = "".join(reversed(digits))
    trace.add("calculation", f"decoded {s}", value=s)
    return s, trace


def complement_code(s):
    _check(s)
    n, x = len(s), _enc(s)
    direct = _enc("".join(ALPHABET[3 - DIGIT[c]] for c in s))
    formula = 4 ** n - 1 - x
    trace = Trace().add("calculation", f"enc(comp S) = {direct}; 4^{n} − 1 − enc(S) = {4 ** n - 1} − {x} = {formula}",
                        direct=direct, formula=formula, holds=direct == formula)
    return direct, trace


def rc_code(x, n):
    value = _rc_code(x, n)
    trace = Trace().add("algorithm", f"RC on codes: read the {n} digits of {x} from the least significant end, "
                        f"complement each (d -> 3 − d), append: {value}", value=value)
    return value, trace


def canonical_code(s):
    _check(s)
    x = _enc(s)
    r = _rc_code(x, len(s))
    c = min(x, r)
    trace = Trace().add("calculation", f"canonical = min(enc(S), enc(RC(S))) = min({x}, {r}) = {c}", value=c)
    trace.add("interpretation", "one code per {I,RC} class: S and RC(S) are the same double-stranded word")
    return c, trace


def kmer_codes(seq, k):
    _check(seq)
    if not 1 <= k <= len(seq):
        raise ValueError("need 1 <= k <= len(seq)")
    mask = 4 ** k
    x = _enc(seq[:k])
    codes = [x]
    for c in seq[k:]:
        x = (4 * x) % mask + DIGIT[c]  # drop the leading digit, append the new one
        codes.append(x)
    trace = Trace().add("algorithm", f"{len(codes)} {k}-mer codes by sliding window x' = (4x mod 4^{k}) + d, O(1) each",
                        count=len(codes))
    return codes, trace


def kmer_census(seqs, k):
    present = set()
    for s in seqs:
        present.update(kmer_codes(s, k)[0])
    canon = {min(x, _rc_code(x, k)) for x in present}
    pal = {x for x in present if _rc_code(x, k) == x}
    pairs = sum(1 for x in present if _rc_code(x, k) in present and _rc_code(x, k) > x)
    space, _ = burnside(("I", "RC"), k)
    result = {"distinct": len(present), "canonical": len(canon), "palindromic": len(pal), "rc_pairs_present": pairs,
              "identity_holds": len(canon) == len(present) - pairs, "space_canonical": space}
    trace = Trace()
    trace.add("calculation", f"k = {k}: {len(present)} distinct k-mers, {len(canon)} canonical, {len(pal)} RC-palindromic; "
              f"{pairs} pairs (w, RC(w)) both present; whole space has {space} canonical classes", **result)
    trace.add("calculation", f"check: {len(present)} − {pairs} = {len(present) - pairs} = #canonical: "
              f"{result['identity_holds']}")
    return result, trace


def bucket_report(codes, m):
    buckets = {}
    for key, x in codes.items():
        buckets.setdefault(x % m, []).append(key)
    groups = [v for v in buckets.values() if len(v) > 1]
    loads = [len(v) for v in buckets.values()]
    bound = -(-len(codes) // m)  # ceil(N/m)
    result = {"buckets": {b: buckets[b] for b in sorted(buckets)}, "occupied": len(buckets),
              "colliding_groups": groups, "collisions": len(codes) - len(buckets),
              "max_load": max(loads) if loads else 0, "pigeonhole_bound": bound}
    trace = Trace()
    trace.add("calculation", f"m = {m}: {len(codes)} keys in {len(buckets)} occupied buckets, max load {result['max_load']}; "
              f"pigeonhole guarantees a bucket with ≥ ⌈{len(codes)}/{m}⌉ = {bound}", **{k: v for k, v in result.items()
                                                                                     if k != "buckets"})
    trace.add("limitation", "sharing a bucket says nothing about similarity: the residue depends on the whole code, "
              "not on biological closeness")
    return result, trace


def last_bases_seen(m):
    # smallest j with m | 4^j, i.e. b(x) depends only on the last j bases; None if m has an odd factor > 1
    j, p = 0, 1
    while p % m:
        j += 1
        p *= 4
        if j > 64:
            return None, Trace().add("calculation", f"no power of 4 is divisible by {m}: every base can affect x mod {m}")
    return j, Trace().add("calculation", f"{m} | 4^{j}: x mod {m} depends only on the last {j} bases", value=j)


def substitution_collision_t_min(m):
    # exact criterion (docs/unit4_expected.md §3): None means no single substitution ever keeps the bucket
    a, q = 0, m
    while q % 2 == 0:
        q //= 2
        a += 1
    if q not in (1, 3):
        t = None
    else:
        t = a // 2 if q == 1 else (a + 1) // 2
    text = (f"m = 2^{a}·{q}: odd part ≥ 5, no single substitution keeps the bucket" if t is None else
            f"m = 2^{a}·{q}: a single substitution can keep the bucket at weight 4^t for t ≥ {t}")
    return t, Trace().add("calculation", text, a=a, q=q, t_min=t)


def substitution_collisions(pairs, seqs, m):
    collided = []
    for a, b in pairs:
        if _enc(seqs[a]) % m == _enc(seqs[b]) % m:
            collided.append((a, b))
    trace = Trace().add("calculation", f"m = {m}: {len(collided)} of {len(pairs)} pairs share a bucket", collided=collided)
    return collided, trace
