"""Read packed sweep streams independently of the C++ codec and check their field transcript."""
import argparse
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from extraction.ball_challenges import BallIndex, FilteredBall


Q = (1 << 48) - 59
KAPPA = 204388862566065
ZERO = (0, 0, 0, 0)
LIMIT = 160


def require(condition, message):
    if not condition:
        raise ValueError(message)


def linear(a, b, weight=1):
    return tuple((x + weight*y) % Q for x, y in zip(a, b))


def multiply(a, b):
    coeffs = [0]*7
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coeffs[i+j] += x*y
    for j in range(6, 3, -1):  # U^4 = 4 U^2 - 2.
        coeffs[j-2] += 4*coeffs[j]
        coeffs[j-4] -= 2*coeffs[j]
    return tuple(x % Q for x in coeffs[:4])


@lru_cache(None)
def rice_layout(n, bound):
    # At a fixed sum of quotients the least energy balances the quotients.
    # Select least unrounded bit budget, with the smaller h winning a tie.
    layouts = []
    for h in range(1, bound.bit_length()+2):
        energy = bound >> (2*h-2)
        lo, hi = 0, isqrt(n*energy)+1
        while hi-lo > 1:
            mid = (hi+lo)//2
            z, r = divmod(mid, n)
            needed = (n-r)*z*z + r*(z+1)**2
            if needed <= energy:
                lo = mid
            else:
                hi = mid
        layouts.append((n*(h+1)+lo, h))
        if not energy:
            break
    bits, h = min(layouts)
    return (bits+7)//8, h


class Reader:
    def __init__(self, path):
        self.data = path.read_bytes()
        self.pos = 0

    def take(self, n):
        end = self.pos+n
        require(end <= len(self.data), "truncated message")
        result = self.data[self.pos:end]
        self.pos = end
        return result

    def scalar(self):
        x = int.from_bytes(self.take(6), "little")
        require(x < Q, "noncanonical coefficient")
        return x

    def field(self, nonzero=False):
        x = tuple(self.scalar() for _ in range(4))
        require(not nonzero or x != ZERO, "forbidden zero field challenge")
        return x

    def scalars(self, n):
        for _ in range(n):
            self.scalar()

    def attempts(self, read_challenge):
        for attempt in range(1, LIMIT+1):
            read_challenge()
            flag = self.take(1)
            if flag == b"\x01":
                return attempt
            require(flag == b"\x00", "noncanonical retry flag")
        raise ValueError("retry budget exhausted without submission")

    def bounded(self, n, bound):
        nbytes, h = rice_layout(n, bound)
        block = self.take(nbytes)
        offset = 0
        def bit():
            nonlocal offset
            require(offset < 8*nbytes, "truncated Rice code")
            v = (block[offset//8] >> (offset%8)) & 1
            offset += 1
            return v
        energy = 0
        for _ in range(n):
            high = 0
            while not bit():
                high += 1
            low = sum(bit() << j for j in range(h))
            z = (high << h) + low
            x = -(z//2+1) if z & 1 else z//2
            energy += x*x
            require(energy <= bound, "Rice vector exceeds accepted energy")
        if offset % 8:
            require(block[offset//8] >> (offset%8) == 0, "nonzero Rice padding bits")
        require(not any(block[(offset+7)//8:]), "nonzero Rice padding bytes")
        return energy


@lru_cache(None)
def profile(k, reference):
    if not reference:
        data = json.loads((ROOT/"parameters/sweep-profiles.json").read_text())["profiles"][str(k)]
        return data["layers"], data["no_retry_bytes"]
    require(k == 1024, "unknown reference profile")
    layers = json.loads((ROOT/"parameters/ledger-two-point.json").read_text())["layers"]
    for i, v in enumerate(layers):
        v.update(b=v["seed"]["b"], h=v["seed"]["h"], G=256*v["S"]-int(i<4))
    return layers, dict(prover=126647, verifier=44409)


def inspect(directory, domain):
    meta = json.loads((directory/"metrics.json").read_text())
    k, reference = meta["k"], meta["reference"]
    layers, base = profile(k, reference)
    p = Reader(directory/"p-to-v.bin")
    v = Reader(directory/"v-to-p.bin")
    v.field()  # alpha can be zero
    gamma = v.field(nonzero=True)
    p.scalars(4*5 + 64*4)
    pending = p.field()
    require(linear(pending, gamma) == ZERO, "ordinary fixture pending must equal -gamma")
    p.scalars(64*4)
    v.field(nonzero=True)  # rho
    eta = v.field(nonzero=True)
    V = 64*k+544
    rounds = (55*V+64*(k+2)+1).bit_length()
    claim = linear((V*KAPPA % Q,0,0,0),multiply(eta,pending))
    for _ in range(rounds):
        g0, g2 = p.field(), p.field()
        r = v.field()
        g1 = linear(linear(claim,g0,-2),g2,-1)
        claim = linear(g0,multiply(r,linear(g1,multiply(r,g2))))
    terminal = p.field()
    require(claim == multiply(terminal,terminal), "invalid square Sumcheck terminal")

    projected, short, image_energies = [], [], []
    extra_p = extra_v = 0
    for i, layer in enumerate(layers):
        s, n, S = layer["s"], layer["n"], layer["S"]
        is_pivot = i == len(layers)-1
        p.scalars(64*(2+16*(s-1) if is_pivot else 6))
        bits = layer["b"] + layer["h"]*(3*layer["b"]-1)
        seed_bytes = (bits+7)//8
        def seed():
            raw = v.take(seed_bytes)
            if bits % 8:
                require(raw[-1] >> (bits%8) == 0, "nonzero seed padding")
        def vector():
            for j in range(s):
                domain.decode(v.take(domain.index.nbytes), nonzero=(is_pivot and j == 0))
        projected.append(p.attempts(seed))
        image_energies.append(p.bounded(384, 384*S))
        v.field()  # zeta
        p.scalars(64*4*s*(s+1)//2-4 if is_pivot else 64*6)
        short.append(p.attempts(vector))
        extra_p += projected[-1]-1 + short[-1]-1
        extra_v += seed_bytes*(projected[-1]-1) + s*domain.index.nbytes*(short[-1]-1)
        if is_pivot:
            response_energy = p.bounded(64*n,layer["G"])
    require(p.pos == len(p.data) and v.pos == len(v.data), "trailing data after pivot")
    require(p.pos == base["prover"]+extra_p and v.pos == base["verifier"]+extra_v,
            "length mismatch against independent padded budgets and decoded retries")
    require(p.pos == meta["prover_bytes"] and v.pos == meta["verifier_bytes"],
            "measured byte counts differ from files")
    require(projected == meta["projection_attempts"] and short == meta["short_attempts"],
            "retry counts differ from measured run")
    return dict(run=directory.name,k=k,reference=reference,rounds=rounds,prover_bytes=p.pos,
                verifier_bytes=v.pos,prover_sha256=hashlib.sha256(p.data).hexdigest(),
                verifier_sha256=hashlib.sha256(v.data).hexdigest(),projection_attempts=projected,
                short_attempts=short,image_energies=image_energies,pivot_response_energy=response_energy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", type=Path, help="one or more completed run directories")
    args = parser.parse_args()
    domain = FilteredBall(BallIndex(64,128),26)
    print(json.dumps([inspect(path,domain) for path in args.results],indent=2))


if __name__ == "__main__":
    main()
