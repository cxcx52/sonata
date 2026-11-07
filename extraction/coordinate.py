"""Local three-point comparison over Z[X]/(X^64+1) with the Sonata modulus.

The caller supplies the block image and offset from a common authenticated
prefix. This routine checks the three image equations and response budgets.
It does not authenticate prefix commitments or recursively obtain the records.
"""
from dataclasses import dataclass
from typing import Literal

from .ball_challenges import negacyclic_mul, spectral_member_exact

Q = 2**48 - 59
DEGREE = 64
OP_BOUND = 26
Ring = tuple[int, ...]
Vector = tuple[Ring, ...]


def _sub(a: Ring, b: Ring) -> Ring:
    return tuple(x-y for x, y in zip(a, b, strict=True))


def _residues(a: Ring) -> Ring:
    return tuple(x % Q for x in a)


def _inverse(a: Ring) -> Ring:
    """Invert in the power-of-two negacyclic tower by taking relative norms."""
    if len(a) == 1:
        return (pow(a[0], -1, Q),)
    conjugate = tuple(x if i % 2 == 0 else -x for i, x in enumerate(a))
    norm = negacyclic_mul(a, conjugate)
    inverse_norm = _inverse(_residues(norm[::2]))
    lifted = tuple(x for value in inverse_norm for x in (value, 0))
    return _residues(negacyclic_mul(conjugate, lifted))


def _dot(a: Vector, b: Vector) -> Ring:
    accum = [0] * DEGREE
    for left, right in zip(a, b, strict=True):
        for j, coeff in enumerate(negacyclic_mul(left, right)):
            accum[j] += coeff
    return _residues(tuple(accum))


@dataclass(frozen=True)
class Record:
    challenge: Ring
    response: Vector  # centered integer lifts


@dataclass(frozen=True)
class Comparison:
    kind: Literal["slope", "kernel"]
    vector: Vector  # slope modulo q, or a short integral kernel


def compare_coordinate(
    matrix: tuple[Vector, ...], block_image: Vector, offset: Vector,
    records: tuple[Record, Record, Record], response_energy: int, *,
    pivot_first: bool = False,
) -> Comparison:
    """Check Az=offset+c*block_image and compare the two divided differences."""
    if not matrix or not matrix[0] or len(records) != 3 or response_energy < 0:
        raise ValueError("invalid system size, record count or squared norm budget")
    width = len(matrix[0])
    if len(matrix) != len(offset) or len(matrix) != len(block_image):
        raise ValueError("image dimensions do not match the public matrix")
    for vector in (*matrix, *(r.response for r in records)):
        if len(vector) != width:
            raise ValueError("inconsistent block width")
    vectors = (*matrix, block_image, offset, *(r.response for r in records))
    for vector in vectors:
        if any(len(x) != DEGREE or any(type(v) is not int for v in x) for x in vector):
            raise ValueError("expected integer polynomials of degree below 64")
    for record in records:
        c = record.challenge
        if (len(c) != DEGREE or any(type(v) is not int for v in c)
                or sum(v*v for v in c) > 128 or not spectral_member_exact(c, OP_BOUND)
                or (pivot_first and not any(c))):
            raise ValueError("challenge outside the coordinate support")
        if any(abs(v) > Q//2 for x in record.response for v in x):
            raise ValueError("response must be centered")
        if sum(v*v for x in record.response for v in x) > response_energy:
            raise ValueError("response exceeds the squared norm budget")
        for row, t, u in zip(matrix, block_image, offset, strict=True):
            prod = negacyclic_mul(c, t)
            if _dot(row, record.response) != tuple((x+y) % Q for x,y in zip(u,prod)):
                raise ValueError("response does not open the declared affine image")
    center, first, second = records
    if len({r.challenge for r in records}) != 3:
        raise ValueError("the center and both forks must be distinct")
    d1 = _sub(first.challenge, center.challenge)
    d2 = _sub(second.challenge, center.challenge)
    x1 = tuple(_sub(z, base) for z,base in zip(first.response, center.response))
    x2 = tuple(_sub(z, base) for z,base in zip(second.response, center.response))
    kernel = tuple(_sub(negacyclic_mul(x, d2), negacyclic_mul(y, d1))
                   for x,y in zip(x1,x2))
    if any(v % Q for x in kernel for v in x):
        if any(any(_dot(row, kernel)) for row in matrix):
            raise ArithmeticError("cross-product kernel equation failed")
        if sum(v*v for x in kernel for v in x) > (8*OP_BOUND)**2 * response_energy:
            raise ArithmeticError("cross-product kernel radius failed")
        return Comparison("kernel", kernel)
    inverse = _inverse(d1)
    candidate = tuple(_residues(negacyclic_mul(x, inverse)) for x in x1)
    if tuple(_dot(row, candidate) for row in matrix) != tuple(map(_residues, block_image)):
        raise ArithmeticError("recovered block image failed")
    return Comparison("slope", candidate)
