#!/usr/bin/env python3
"""Exact reference tools for integer-ball and spectrally filtered challenges.

Research prototype, not a constant-time or audited cryptographic implementation.
Rank codecs, membership decisions and rejection samplers are exact. A deterministic
getrandbits source can be injected for replay; the default uses secrets.randbits.
"""
from __future__ import annotations
from functools import lru_cache
from fractions import Fraction
from math import isqrt
import secrets
from typing import Callable, Iterable


def ceildiv(a: int, b: int) -> int:
    return -((-a)//b)


class BallIndex:
    """Lexicographic indexing of {c in Z^d : sum(c_j^2) <= energy}."""
    def __init__(self, degree: int = 64, energy: int = 128) -> None:
        if not isinstance(degree,int) or not isinstance(energy,int) or degree < 1 or energy < 0:
            raise ValueError('degree must be positive and energy nonnegative')
        self.degree, self.energy = degree, energy
        # table[r][b] counts r coordinates with energy at most b.
        self.table = [[1]*(energy+1)]
        for r in range(1, degree+1):
            prev = self.table[-1]
            row = []
            for b in range(energy+1):
                row.append(prev[b] + 2*sum(prev[b-a*a]
                                          for a in range(1, isqrt(b)+1)))
            self.table.append(row)
        self.size = self.table[degree][energy]
        self.bits = (self.size-1).bit_length()
        self.nbytes = (self.bits+7)//8
        self.zero_rank = (self.size-1)//2

    def rank(self, c: Iterable[int]) -> int:
        c = tuple(c)
        if len(c) != self.degree or any(not isinstance(x, int) for x in c):
            raise ValueError('wrong degree or noninteger coefficient')
        if sum(x*x for x in c) > self.energy:
            raise ValueError('challenge outside the integer ball')
        rank, remaining = 0, self.energy
        for j, x in enumerate(c):
            row = self.table[self.degree-j-1]
            for a in range(-isqrt(remaining), x):
                rank += row[remaining-a*a]
            remaining -= x*x
        return rank

    def unrank(self, rank: int) -> tuple[int, ...]:
        if not isinstance(rank, int) or not 0 <= rank < self.size:
            raise ValueError('rank outside canonical range')
        remaining, out = self.energy, []
        for j in range(self.degree):
            row = self.table[self.degree-j-1]
            for a in range(-isqrt(remaining), isqrt(remaining)+1):
                count = row[remaining-a*a]
                if rank < count:
                    out.append(a); remaining -= a*a; break
                rank -= count
            else:
                raise AssertionError('unranking invariant violated')
        if rank != 0:
            raise AssertionError('unranking ended with a nonzero residual rank')
        return tuple(out)

    def encode(self, c: Iterable[int]) -> bytes:
        return self.rank(c).to_bytes(self.nbytes, 'little')

    def decode(self, data: bytes) -> tuple[int, ...]:
        if len(data) != self.nbytes:
            raise ValueError('wrong fixed payload length')
        return self.unrank(int.from_bytes(data, 'little'))

    def sample_rank(self, randbits: Callable[[int], int] = secrets.randbits,
                    nonzero: bool = False) -> int:
        # Zero has the middle lexicographic rank, not rank zero.
        bound = self.size-int(nonzero)
        if bound <= 0:
            raise ValueError('empty requested support')
        bits = (bound-1).bit_length()
        while True:
            r = randbits(bits)
            if 0 <= r < bound:
                return r + int(nonzero and r >= self.zero_rank)

    def exact_energy_counts(self) -> list[int]:
        row = self.table[self.degree]
        return [row[0]] + [row[b]-row[b-1] for b in range(1, self.energy+1)]


# Fixed-point outward-rounded real intervals, all endpoints in units 2^-precision.
def _add(a, b): return (a[0]+b[0], a[1]+b[1])
def _sub(a, b): return (a[0]-b[1], a[1]-b[0])
def _mul(a, b, scale):
    v = [x*y for x in a for y in b]
    return (min(v)//scale, ceildiv(max(v), scale))
def _half(a): return (a[0]//2, ceildiv(a[1], 2))
def _sqrt(a, scale):
    if a[0] < 0: raise ValueError('negative interval in square-root construction')
    lo = isqrt(a[0]*scale)
    hi = isqrt(a[1]*scale)
    hi += int(hi*hi < a[1]*scale)
    return (lo, hi)
def _cmul(a, b, scale):
    return (_sub(_mul(a[0],b[0],scale), _mul(a[1],b[1],scale)),
            _add(_mul(a[0],b[1],scale), _mul(a[1],b[0],scale)))


@lru_cache(maxsize=16)
def root_intervals(degree: int, precision: int = 128):
    """Certified intervals for exp(pi*i*r/d), r=0,...,2d-1.

    Uses half-angle square roots starting at pi/2; does not use floating-point pi.
    """
    if degree < 2 or degree & (degree-1):
        raise ValueError('a power-of-two degree of at least two is required')
    S = 1 << precision
    one = (S,S); zero = (0,0)
    co, si = zero, one
    d = 2
    while d < degree:
        old = co
        co = _sqrt(_half(_add(one, old)), S)
        si = _sqrt(_half(_sub(one, old)), S)
        d *= 2
    root = (co,si)
    out = [(one,zero)]
    for _ in range(1,2*degree):
        out.append(_cmul(out[-1],root,S))
    return tuple(out)


def _square_bounds(a):
    lo,hi=a
    lower=0 if lo <= 0 <= hi else min(lo*lo,hi*hi)
    return lower,max(lo*lo,hi*hi)


def spectral_interval_decision(c: tuple[int,...], bound: int,
                               precision: int = 128) -> bool | None:
    """True/False if certified; None if an interval straddles the boundary."""
    d=len(c); roots=root_intervals(d,precision); threshold=bound*bound*(1<<(2*precision))
    uncertain=False
    for t in range(d//2):
        re=[0,0];im=[0,0]
        for j,x in enumerate(c):
            (rl,rh),(il,ih)=roots[((2*t+1)*j)%(2*d)]
            if x>=0:
                re[0]+=x*rl;re[1]+=x*rh;im[0]+=x*il;im[1]+=x*ih
            else:
                re[0]+=x*rh;re[1]+=x*rl;im[0]+=x*ih;im[1]+=x*il
        rlo,rhi=_square_bounds(re);ilo,ihi=_square_bounds(im)
        if rlo+ilo>threshold:return False
        if rhi+ihi>threshold:uncertain=True
    return None if uncertain else True


def negacyclic_mul(a: tuple[int,...], b: tuple[int,...]) -> tuple[int,...]:
    if len(a)!=len(b):raise ValueError('ring degree mismatch')
    d=len(a);r=[0]*d
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            k=i+j;r[k%d]+=(1 if k<d else -1)*x*y
    return tuple(r)


def psd_exact(matrix: list[list[int]]) -> bool:
    """Exact symmetric PSD decision by rational Schur complements.

    A zero diagonal with any nonzero entry in that row certifies non-PSD.
    This correctly handles semidefinite, not merely positive-definite, matrices.
    """
    n=len(matrix)
    if any(len(r)!=n for r in matrix):raise ValueError('matrix is not square')
    if any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(n)):
        raise ValueError('matrix is not symmetric')
    a=[[Fraction(x) for x in row] for row in matrix]
    for k in range(n):
        if any(a[j][j]<0 for j in range(k,n)):return False
        if a[k][k]==0:
            if any(a[k][j]!=0 for j in range(k+1,n)):return False
            continue
        p=a[k][k]
        col=[a[j][k] for j in range(n)]
        for i in range(k+1,n):
            for j in range(i,n):
                a[i][j]-=col[i]*col[j]/p
                a[j][i]=a[i][j]
    return True


def spectral_member_exact(c: tuple[int,...], bound: int = 32) -> bool:
    if not isinstance(bound,int) or bound<0:raise ValueError('spectral bound must be a nonnegative integer')
    if any(not isinstance(x,int) for x in c):raise ValueError('coefficients must be integers')
    # Usually decides immediately, with rigorous dyadic intervals.
    result=spectral_interval_decision(c,bound,128)
    if result is not None:return result
    d=len(c)
    adj=(c[0],)+tuple(-c[d-j] for j in range(1,d))
    gram=negacyclic_mul(adj,c)
    h=[[int(i==j)*bound*bound-(1 if i>=j else -1)*gram[(i-j)%d]
        for j in range(d)] for i in range(d)]
    return psd_exact(h)


class FilteredBall:
    """The exact support is the ball intersected with an operator-norm cutoff.

    Ranks and wire encodings are AMBIENT ball ranks. The filtered support size
    need not be known. A decoder must reject ranks whose vector fails the filter.
    """
    def __init__(self, index: BallIndex, op_bound: int):
        if index.degree<2 or index.degree & (index.degree-1):
            raise ValueError('negacyclic spectral interface needs power-of-two degree')
        self.index,self.op_bound=index,op_bound

    def member(self,c: Iterable[int]) -> bool:
        c=tuple(c)
        if len(c)!=self.index.degree or any(not isinstance(x,int) for x in c):return False
        if sum(x*x for x in c)>self.index.energy:return False
        return spectral_member_exact(c,self.op_bound)

    def sample(self,randbits:Callable[[int],int]=secrets.randbits,
               nonzero:bool=False) -> tuple[int,...]:
        while True:
            r=self.index.sample_rank(randbits,nonzero=nonzero)
            c=self.index.unrank(r)
            if self.member(c):return c

    def encode(self,c:Iterable[int],nonzero:bool=False) -> bytes:
        c=tuple(c)
        if not self.member(c) or (nonzero and not any(c)):
            raise ValueError('challenge outside requested filtered support')
        return self.index.encode(c)

    def decode(self,data:bytes,nonzero:bool=False) -> tuple[int,...]:
        c=self.index.decode(data)
        if not self.member(c) or (nonzero and not any(c)):
            raise ValueError('canonical ambient rank fails support validation')
        return c


def filtered_fiber_ranks(domain: FilteredBall, center_rank: int,
                         randbelow: Callable[[int],int]=secrets.randbelow, *,
                         nonzero: bool=False):
    """Uniform permutation of filtered ranks other than a valid center.

    The ambient lazy permutation state is not part of a target-oracle checkpoint.
    Its tested-invalid points still incur cost; bounds must include that work.
    """
    idx=domain.index
    if not 0<=center_rank<idx.size or not domain.member(idx.unrank(center_rank)) or (nonzero and center_rank==idx.zero_rank):
        raise ValueError('invalid center')
    remaining=idx.size-1;swaps={}
    while remaining:
        r=randbelow(remaining)
        if not 0<=r<remaining:raise ValueError('invalid randbelow output')
        v=swaps.get(r,r);last=remaining-1;tail=swaps.get(last,last)
        if r!=last:swaps[r]=tail
        swaps.pop(last,None);remaining-=1
        rank=v+int(v>=center_rank)
        if (not nonzero or rank!=idx.zero_rank) and domain.member(idx.unrank(rank)):
            yield rank
