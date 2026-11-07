"""Independently check the two fixed-profile ledgers and count packed messages.

Run from the repository root with python parameters/audit_profile.py.
All byte counts describe an encoding model. No prover or serializer is run.
"""

from fractions import Fraction
from functools import lru_cache
import json
from math import isqrt, log2, prod
from pathlib import Path
import re

from challenge_moment import challenge_support


ROOT = Path(__file__).resolve().parent
Q = 2**48 - 59
K = Q**4
PATH = ((21, 64), (10, 128), (6, 128), (6, 128), (5, 0))
SIGMA2 = Fraction(38400, 4121)
RING_BYTES = 64 * 6
FIELD_BYTES = 4 * 6
ATTEMPTS = 160


def ceildiv(a, b):
    return -(-a // b)


def ceilroot(x):
    floor = isqrt(int(x))
    return floor + (floor**2 < x)


def bits(probability):
    return log2(probability.denominator) - log2(probability.numerator)


def max_group_weight():
    # For fixed popcount of the high 42 bits all 64 low limbs occur, except
    # when all high bits are one. Then q permits only low values 0 through 4.
    weights = (1, 2, 4, 6, 7)
    masks = {t: min(m.bit_count() for m in range(32)
                    if sum(w for j, w in enumerate(weights) if m >> j & 1) == t)
             for t in range(21)}
    largest = 0
    for zeros in range(43):
        exceptional = int(zeros == 0)
        helpers = masks[(zeros - 1 + exceptional)//2] + exceptional
        low_limit = Q - ((1 << 42) - 1)*64 if exceptional else 64
        for low in range(low_limit):
            ones = 42 - zeros + low.bit_count() + helpers + exceptional*((low >> 2) & 1)
            largest = max(largest, ones)
    return largest


@lru_cache(None)
def max_digit_energy(radix):
    """Maximise over |x| <= (q-1)/2 using base-radix digits and balanced carries."""
    bound = (Q - 1) // 2
    # comparison is the order of the low digits seen so far vs the bound.
    scores = {(0, 0): 0}
    while bound:
        bound, limit = divmod(bound, radix)
        updated = {}
        for (carry, comparison), score in scores.items():
            for digit in range(radix):
                value = digit + carry
                carry_out = int(value > radix // 2)
                balanced = value - carry_out * radix
                comparison_out = comparison if digit == limit else (1 if digit > limit else -1)
                key = carry_out, comparison_out
                updated[key] = max(updated.get(key, -1), score + balanced**2)
        scores = updated
    return max(v + carry**2 for (carry, cmp), v in scores.items() if cmp <= 0)


def high_energy(response, low, radix):
    # Least integer H with radix^2 H >= (sqrt(response)+sqrt(low))^2.
    cross = ceilroot(4 * response * low)
    return ceildiv(response + low + cross, radix**2)


def rice_bytes(n, energy):
    """Padded signed Rice bound via least energy at a fixed sum of quotients."""
    def least_squares(total):
        a, r = divmod(total, n)
        return (n - r)*a*a + r*(a + 1)**2

    best = None
    for h in range(1, energy.bit_length() + 1):
        budget = energy // (1 << (2*h - 2))
        lo, hi = 0, isqrt(n*budget) + 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if least_squares(mid) <= budget:
                lo = mid
            else:
                hi = mid
        length = n*(h + 1) + lo
        best = length if best is None else min(best, length)
    return ceildiv(best, 8)


def projection_parameters(N, S, minimum):
    stream_bits = 768*N
    states = (stream_bits + 1)*4*Q*(384*S + 2)
    width = (states - 1).bit_length()
    candidates = []
    for h in range(1, stream_bits.bit_length() + 1):
        tree = (1 << h) - 1
        # Check the original Eq. 50 certificate using integer inequalities.
        required = 27*h*tree**2 * 2**(6*width + 3*minimum - 2)
        b = max((required - 1).bit_length(), ceildiv(stream_bits, 1 << h))
        if b*(1 << (h - 1)) < stream_bits <= b*(1 << h):
            candidates.append((b + h*(3*b - 1), b, h))
    length, b, h = min(candidates)
    # Odd b is rounded down before exponentiation, as in the rational ledger.
    d = Fraction(((1 << h) - 1)*states, 1 << (b//2))
    assert Q**2 >= 176**2 * SIGMA2 * S
    return b, h, ceildiv(length, 8), d


def require_equal(label, actual, expected):
    if actual != expected:
        raise ValueError(f"{label}: {actual!r} != {expected!r}")


def audit():
    ball, moment, support = challenge_support()
    rank_bytes = ceildiv((ball - 1).bit_length(), 8)
    raw = 55*(64*1024 + 2*272) + 1
    group_weight = max_group_weight()
    S = group_weight*(64*1024 + 2*272) + 1
    rows = 256*(1024 + 3) + 5
    state_radius = ceilroot(4*SIGMA2*S)
    roles = {"A_st": (4, ceildiv(73*4*48, 64), state_radius),
             "A_aux": (4, ceildiv(7*(64*1024 + 544), 64), state_radius)}
    layers, config = [], []
    projection_loss = aggregation_loss = coordinate_loss = Fraction()
    for i, (s, radix) in enumerate(PATH):
        n = ceildiv(raw, 64*s)
        N = 64*s*n
        degree = ceildiv(rows + N - raw + 384, 4) - 1
        b, h, seed_bytes, d = projection_parameters(N, S, 150 if i < 2 else 151)
        projection_loss += ATTEMPTS*(Fraction(1, 2**192) + d)
        aggregation_loss += Fraction(degree, K)
        coordinate_loss += 2*ATTEMPTS*(Fraction(s - int(not radix), support)
                                      + Fraction(int(not radix), support - 1))
        G = (256*S - 1) if radix else ceildiv(256*support*S, support - 1) - 1
        config.append([s, n, radix, raw, N, S, G, b, h])
        projection = rice_bytes(384, 384*S)
        # One stage-control byte per issued seed or short-challenge vector.
        prover = 2 + projection
        verifier = seed_bytes + FIELD_BYTES + rank_bytes*s
        ht = s*(s + 1)//2
        if radix:
            ell = ceildiv(48, radix.bit_length() - 1)
            bt, bh = ell*16*s, ell*4*ht
            next_raw = 64*(2*n + bt + bh)
            low = 64*n*(radix//2)**2
            next_S = low + high_energy(G, low, radix)
            next_S += 64*(16*s + 4*ht)*max_digit_energy(radix)
            for name, count in ((f"B_{i}", bt), (f"D_{i}", bh)):
                roles[name] = (6, count, ceilroot(4*SIGMA2*next_S))
            z2 = (1 + radix**2)*SIGMA2*next_S
            prover += 12*RING_BYTES
        else:
            pivot_cols = 16*ceildiv(48, 3)
            roles["B_pivot"] = (2, pivot_cols, ceilroot(64*pivot_cols*8**2))
            z2 = G
            # Public t_1..t_(s-1), digit-tag of t_0, triangle minus four field
            # coefficients, and the final bounded response.
            prover += (16*(s - 1) + 2 + 4*ht)*RING_BYTES - 4*6 + rice_bytes(64*n, G)
        beta = ceilroot(64*26**2*z2)
        anchor = ceilroot(z2) + ceilroot(26**2*s*SIGMA2*S)
        assert anchor <= beta < Q
        roles[f"A_{i}"] = (16, n, beta)
        layers.append(dict(i=i, raw=raw, N=N, n=n, S=S, G=G, b=b, h=h, s=s,
                           rows=rows+N-raw, degree=degree, seed_bytes=seed_bytes,
                           transport_bits_floor=(d.denominator//d.numerator).bit_length()-1,
                           cross_radius=beta, anchor_radius_upper=anchor,
                           projection_bytes=projection, prover_bytes=prover, verifier_bytes=verifier))
        if radix:
            raw, S, rows = next_raw, next_S, 64*(6+6+16+4)+4

    representatives = {
        k for k, (a, b, r) in roles.items() if not any(
            key != k and a == a2 and b <= b2 and r <= r2
            for key, (a2, b2, r2) in roles.items())
    }
    embedding = {}
    for key, (a, b, r) in roles.items():
        choices = [(roles[g][1], roles[g][2], g) for g in representatives
                   if roles[g][0] == a and roles[g][1] >= b and roles[g][2] >= r]
        embedding[key] = min(choices)[2]
    crs_cells = max(a*b for a, b, r in roles.values())
    backend = projection_loss + aggregation_loss + coordinate_loss
    node_rounds = (55*(64*1024 + 544) + 64*(1024 + 2) + 1).bit_length()
    transition = backend + Fraction(72, K-1) + Fraction(2*node_rounds, K)
    node = transition + Fraction(4*1024, K-1)
    for filename, points in (("ledger-two-point.json", 2), ("ledger.json", 3)):
        ref = json.loads((ROOT/filename).read_text(encoding="utf-8"))
        require_equal("layers", len(layers), len(ref["layers"]))
        for row, old in zip(layers, ref["layers"]):
            for field in ("i", "s", "raw", "n", "N", "S", "degree", "cross_radius", "anchor_radius_upper"):
                require_equal(field, row[field], old[field])
            for new, old_key in (("seed_bytes","bytes"), ("transport_bits_floor","transition_matrix_bits"),
                                 ("b","b"), ("h","h")):
                require_equal(new, row[new], old["seed"][old_key])
        require_equal("role tuples", roles, {k: tuple(v) for k, *v in ref["roles"]})
        require_equal("MSIS embedding", embedding, ref["embedding"])
        require_equal("CRS", crs_cells, ref["flat_CRS_cells"])
        require_equal("expected pivot coefficient", prod(2*(1+(points-1)*s) for s,r in PATH),
                      ref["backend_pivot_call_coefficient"])
        for field, val in (("backend_exact",backend), ("transition_exact",transition),
                           ("zero_pending_node_exact",node)):
            require_equal(field, val, Fraction(ref[field]))

    cpp = (ROOT.parent/"examples/five_layer.hpp").read_text(encoding="utf-8")
    block = cpp.split("layers{{", 1)[1].split("}};", 1)[0]
    require_equal("C++ full-profile configurations", config,
                  [list(map(int, re.findall(r"\d+", row))) for row in re.findall(r"\{([^{}]+)\}", block)])
    p = sum(row["prover_bytes"] for row in layers)
    v = sum(row["verifier_bytes"] for row in layers)
    extra_p = (ATTEMPTS - 1)*2*len(layers)
    extra_v = (ATTEMPTS - 1)*sum(row["seed_bytes"] + rank_bytes*row["s"] for row in layers)
    V = 64*994 + 272*(8+1)
    window_rounds = (55*V + 64*(994+8+1) + 1).bit_length()
    window = backend + Fraction(4*994 + 69*8 + 4 - 1, K-1) + Fraction(2*window_rounds,K)
    capacity = json.loads((ROOT/"epoch-capacity.json").read_text(encoding="utf-8"))
    window_row = next(r for r in capacity["cases"] if r["h"] == 8 and r["n"] == 994)
    require_equal("strict window dimensions",
                  [V, 55*V+1, 55*V+64*(994+8+1)+2, window_rounds, group_weight*V+1],
                  [window_row[k] for k in ("V", "raw", "T", "ell", "honest_energy")])
    require_equal("strict window statistical bits (rounded)", round(bits(window),8),
                  round(window_row["merged_bits"],8))
    frontend_p = FIELD_BYTES*(2*window_rounds+1) + 4*RING_BYTES
    frontend_v = FIELD_BYTES*(window_rounds+2)
    # Explicit fresh node, common public old/fresh inputs already supplied.
    # Output public u+4 x, s_new and two 4-ring commitments; four batch/fold coins.
    node_p = FIELD_BYTES*(2*node_rounds + 1 + 6) + 8*RING_BYTES
    node_v = FIELD_BYTES*(node_rounds + 4)
    return dict(model="fixed profile, mathematical bounds and packed byte budgets",
                layers=layers, representatives={g:roles[g] for g in sorted(representatives)},
                embedding=embedding, crs_bytes=crs_cells*RING_BYTES, max_encoded_group_energy=group_weight,
                challenge=dict(ball_size=str(ball), moment12=str(moment), support_lower=str(support)),
                statistical_bits=dict(projection=bits(projection_loss), aggregation=bits(aggregation_loss),
                                      short_challenges=bits(coordinate_loss), zero_pending=bits(node),
                                      authenticated_16384=bits(16384*node), strict_h8_n994=bits(window)),
                packed_bytes=dict(backend_P=p, backend_V=v, backend=p+v,
                                  backend_all_160_attempts_P=p+extra_p, backend_all_160_attempts_V=v+extra_v,
                                  window_P=p+frontend_p, window_V=v+frontend_v, window=p+v+frontend_p+frontend_v,
                                  explicit_node_P=p+node_p, explicit_node_V=v+node_v, explicit_node=p+v+node_p+node_v),
                scope="Conditional on the stated tail and extraction lemmas. Computational advantages, global clock truncation, and actual network traffic are not evaluated.")


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
