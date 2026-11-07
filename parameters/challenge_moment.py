"""Exact Fourier moment certificate for the spectrally filtered challenge set."""
from functools import lru_cache
from math import comb, isqrt


def ball_moment(degree, budget, order):
    """Sum |c(root)|^(2*order) over coefficients of squared norm <= budget.

    root is exp(pi*i/degree). Store the signed integer sum at each exact
    energy. Multiplication by X permutes the ball and multiplies evaluation
    by root, so a mixed moment (p,q) vanishes unless 2*degree divides p-q.
    """
    if degree < 1 or degree & (degree - 1) or min(order, budget) < 0:
        raise ValueError("expected a power-of-two degree and nonnegative bounds")

    @lru_cache(None)
    def shells(d, p, q):
        result = [0] * (budget + 1)
        if (p - q) % (2*d):
            return result
        if d == 1:
            for x in range(-isqrt(budget), isqrt(budget) + 1):
                result[x*x] += x**(p+q)
            return result
        # c(root) = even(root^2) + root * odd(root^2). After the
        # binomial expansions all surviving odd-part phases are +1 or -1.
        for i in range(p + 1):
            for j in range(q + 1):
                odd_phase = (p-i) - (q-j)
                if (i-j) % d or odd_phase % d:
                    continue
                sign = -1 if (odd_phase // d) % 2 else 1
                factor = sign * comb(p, i) * comb(q, j)
                left = shells(d//2, i, j)
                right = shells(d//2, p-i, q-j)
                for e, count in enumerate(left):
                    if count:
                        scaled = factor * count
                        for f in range(budget - e + 1):
                            result[e+f] += scaled * right[f]
        return result

    return sum(shells(degree, order, order))


@lru_cache(None)
def challenge_support():
    """Return the ambient size, twelfth moment sum and filtered size lower bound."""
    size = ball_moment(64, 128, 0)
    moment = ball_moment(64, 128, 6)
    # 32 conjugate pairs share the same moment by cyclotomic symmetry.
    lower = size - 32*moment // 26**12
    return size, moment, lower
