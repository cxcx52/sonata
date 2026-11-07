"""Rebuild Sonata's fixed five-layer, four-lane interactive extraction ledger."""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path


Q = 2**48 - 59
RING_DEGREE = 64
LANES = 4
PROJECTION_ROWS = 384
ATTEMPTS = 160
OPERATOR_BOUND = 26
ARITY = 1024
SIGMA_SQUARED = Fraction(38400, 4121)
PATH = [(21, 64), (10, 128), (6, 128), (6, 128), (5, None)]

# Round 5, Fourier moment certificate for the 64-dimensional integer coefficient
# ball of squared radius 128. The numerator is the sum of |c(omega)|^12
# over all points, for ONE primitive 128th root of unity omega.
BALL_SIZE = 938506860974420376750953616754653564094709448577
TWELFTH_MOMENT_SUM = 1641672596691811208878076318430530242119244719371205842103364608
CHALLENGE_COUNT_LOWER = BALL_SIZE - 32 * TWELFTH_MOMENT_SUM // OPERATOR_BOUND**12


def ceil_div(a, b):
    return (a + b - 1) // b


def ceil_sqrt(value):
    value = Fraction(value)
    root = math.isqrt(value.numerator // value.denominator)
    return root + (root**2 * value.denominator < value.numerator)


def ceil_surd(a, b, denominator):
    """Smallest integer n with denominator*n >= a+sqrt(b)."""
    n = ceil_div(a + math.isqrt(b), denominator)
    if (denominator * n - a)**2 < b:
        n += 1
    return n


def digit_energy(x, radix):
    energy = 0
    while x:
        high = (x + radix // 2 - 1) // radix
        energy += (x - radix * high)**2
        x = high
    return energy


def max_digit_energy(bound, radix):
    """Largest sum of squared balanced digits on [0, bound]."""
    half = radix // 2
    if bound <= half:
        return bound**2
    high = (bound + half - 1) // radix
    complete_blocks = half**2 + max_digit_energy(high - 1, radix)
    last_low_energy = max((half - 1)**2, (bound - radix * high)**2)
    return max(complete_blocks, digit_energy(high, radix) + last_low_energy)


def child_layout(n, energy, arity, radix):
    digits = ceil_div(48, radix.bit_length() - 1)
    t_entries = 16 * arity
    h_entries = LANES * arity * (arity + 1) // 2
    t_cols, h_cols = digits * t_entries, digits * h_entries
    raw = RING_DEGREE * (2*n + t_cols + h_cols)
    low_energy = RING_DEGREE * n * (radix // 2)**2
    response_energy = 256 * energy - 1
    high_energy = ceil_surd(response_energy + low_energy,
                            4 * response_energy * low_energy, radix**2)
    tag_energy = RING_DEGREE * (t_entries + h_entries) * max_digit_energy((Q-1)//2, radix)
    return raw, low_energy + high_energy + tag_energy, t_cols, h_cols


def certify_seed(N, energy, target):
    """Choose the round-13 seed, then evaluate both fixed-test transport bounds."""
    random_bits = 2 * PROJECTION_ROWS * N
    states = (random_bits + 1) * 4 * Q * (PROJECTION_ROWS * energy + 2)
    width = (states - 1).bit_length()
    choices = []
    for h in range(1, random_bits.bit_length()):
        tree = (1 << h) - 1
        log_factor = (27 * h * tree**2 - 1).bit_length()
        b = max(6 * width + 3 * target - 2 + log_factor, ceil_div(random_bits, 1 << h))
        if (1 << (h-1)) * b < random_bits <= (1 << h) * b:
            choices.append((b + h*(3*b - 1), b, h))
    bits, b, h = min(choices)
    tree = (1 << h) - 1
    factor = 27 * h * tree**2
    eq50 = (b - 6*width + 2 - (factor - 1).bit_length()) // 3
    assert factor << (6*width + 3*eq50) <= (4 << b) < factor << (6*width + 3*(eq50+1))
    transport = Fraction(tree * states, 1 << (b // 2))
    assert Q**2 >= 176**2 * SIGMA_SQUARED * energy
    return dict(b=b, h=h, bytes=ceil_div(bits, 8), width=width, eq50_bits=eq50,
                transition_matrix_bits=(transport.denominator // transport.numerator).bit_length()-1), transport


def bit_bound(error):
    return math.log2(error.denominator) - math.log2(error.numerator)


def build_ledger(points=3):
    if points not in (2, 3):
        raise ValueError("expected two or three points per coordinate")
    values = 64*ARITY + 544
    raw = 55*values + 1
    energy = 47*values + 1
    rows = 256 * (ARITY+3) + 5
    assert math.ceil(956 * SIGMA_SQUARED * energy + 3281*values) < Q
    prefix_radius = ceil_sqrt(4 * SIGMA_SQUARED * energy)
    roles = [("A_st", 4, 219, prefix_radius), ("A_aux", 4, ceil_div(7*values,64), prefix_radius)]
    layers = []
    backend, backend_eq50 = Fraction(0), Fraction(0)
    M = CHALLENGE_COUNT_LOWER
    for i, (s, radix) in enumerate(PATH):
        n = ceil_div(raw, RING_DEGREE * s)
        N = RING_DEGREE * s * n
        degree = ceil_div(rows + N - raw + PROJECTION_ROWS, LANES) - 1
        seed, transport = certify_seed(N, energy, 150 if i < 2 else 151)
        gamma = SIGMA_SQUARED * energy
        layer = dict(i=i, s=s, n=n, raw=raw, N=N, S=energy, degree=degree, seed=seed)
        if radix is not None:
            raw, energy, t_cols, h_cols = child_layout(n, energy, s, radix)
            rows = RING_DEGREE * (12+16+LANES) + LANES
            tag_radius = ceil_sqrt(4 * SIGMA_SQUARED * energy)
            roles.extend([(f"B_{i}", 6, t_cols, tag_radius), (f"D_{i}", 6, h_cols, tag_radius)])
            z_squared = (1 + radix**2) * SIGMA_SQUARED * energy
            coordinate_loss = Fraction(s, M)
        else:
            z_squared = ceil_div(256*M*energy, M-1) - 1
            roles.append(("B_pivot", 2, 256, 1024))
            coordinate_loss = Fraction(s-1, M) + Fraction(1, M-1)
        cross = ceil_sqrt(64 * OPERATOR_BOUND**2 * z_squared)
        anchor = ceil_sqrt(z_squared) + ceil_sqrt(OPERATOR_BOUND**2 * s * gamma)
        assert anchor <= cross < Q
        roles.append((f"A_{i}", 16, n, cross))
        layer.update(cross_radius=cross, anchor_radius_upper=anchor)
        layers.append(layer)
        common = Fraction(degree, Q**4) + 2*ATTEMPTS*coordinate_loss + Fraction(ATTEMPTS, 2**192)
        backend += common + ATTEMPTS*transport
        backend_eq50 += common + Fraction(ATTEMPTS, 2**seed["eq50_bits"])

    representatives = {"A_aux", "A_0", "A_1", "D_0", "D_1", "B_pivot"}
    embedding = {}
    for role, rows, columns, radius in roles:
        eligible = [r for r in roles if r[0] in representatives and r[1] == rows
                    and columns <= r[2] and radius <= r[3]]
        embedding[role] = min(eligible, key=lambda r: (r[2], r[3]))[0]
    front = Fraction(72, Q**4-1) + Fraction(44, Q**4)
    zero_return = Fraction(4 * ARITY, Q**4-1)
    node = backend + front + zero_return
    return dict(layers=layers, roles=roles, flat_CRS_cells=max(r*c for _,r,c,_ in roles),
                embedding=embedding, backend_pivot_call_coefficient=math.prod(2*(1+(points-1)*s) for s,_ in PATH),
                baseline_eq50_node_bits=bit_bound(backend_eq50 + front + zero_return),
                backend_exact=str(backend), transition_exact=str(backend+front),
                zero_pending_node_exact=str(node), node_bits=bit_bound(node),
                authenticated_node_union_bits={H: bit_bound(H*node) for H in (4096,16384)},
                scope="Classical interactive statistical errors only. Role tuple: name, rows, columns, Euclidean radius.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare against the selected frozen ledger")
    parser.add_argument("--points", type=int, choices=(2, 3), default=3,
                        help="center plus one or two distinct forks per coordinate")
    args = parser.parse_args()
    rendered = json.dumps(build_ledger(args.points), indent=2)
    if args.check:
        filename = "ledger.json" if args.points == 3 else "ledger-two-point.json"
        expected = json.loads(Path(__file__).with_name(filename).read_text(encoding="utf-8"))
        if json.loads(rendered) != expected:
            raise SystemExit("Ledger mismatch")
        print(f"PASS: complete five-layer ledger matches {filename}")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
