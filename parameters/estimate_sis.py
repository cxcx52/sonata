"""Reproduce coefficient-l2 Module-SIS attack screens with Sage and lattice-estimator."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys


PIN = "53da5982597709ba0fdf94ea37a84d822310fd84"
Q = 2**48 - 59
DEGREE = 64


def best_gsa_dimension(rows, columns, radius):
    """Choose d <= m maximizing log(delta) under the pinned Euclidean GSA."""
    volume = rows * math.log(Q)
    target = math.log(radius)
    # Critical point of (log(radius)-log(volume of lattice)/d)/(d-1).
    peak = (volume + math.sqrt(volume * (volume-target))) / target
    lower = rows + 1
    candidates = {lower, columns, max(lower, min(columns, math.floor(peak))),
                  max(lower, min(columns, math.ceil(peak)))}
    return max(candidates, key=lambda d: (target-volume/d)/(d-1))


def one_factor_screen(a, b, radius, solve):
    # pi and conjugate are factors of q in Z[X]/(X^64+1), with X^32=i.
    # The coefficient norm of multiplication by pi scales by sqrt(q) exactly.
    if 12993454**2 + 10613441**2 != Q:
        raise ValueError("invalid Gaussian prime factor")
    budget = radius**2 // Q
    if not budget:
        return dict(status="impossible by integer norm", reduced_energy=0)
    n, m = 32*a, 64*b
    target = math.sqrt(budget)
    volume = n * math.log(Q)
    # At fixed block size the predicted log length is convex in d. Enforce
    # beta <= d; the unconstrained root-factor optimum can violate this.
    for block in range(40, min(m, 65536) + 1):
        t = (math.log(block/(2*math.pi*math.e)) +
             math.log(math.pi*block)/block) / (2*(block-1))
        stationary = math.sqrt(volume/t)
        lower = max(n+1, block)
        dims = [min(m, max(lower, f(stationary))) for f in (math.floor, math.ceil)]
        d = min(dims, key=lambda d: (d-1)*t+volume/d)
        if (d-1)*t+volume/d <= math.log(target):
            # Validate with the pinned estimator's Euclidean routine at this d.
            reported_block, log_cost = solve(n, m, budget, d)
            if reported_block != block:
                raise ValueError("root-Hermite regression mismatch")
            return dict(status="GSA screen gives a finite fit", reduced_energy=budget,
                        reduced_radius=target, n=n, m=m, chosen_dimension=d,
                        bkz_block_size=block, classical_core_log2=log_cost)
    return dict(status="no feasible GSA block size through 65536", reduced_energy=budget,
                reduced_radius=target, n=n, m=m)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("estimator", type=Path, help="unmodified pinned lattice-estimator checkout")
    parser.add_argument("ledger", type=Path, help="parameters/ledger-two-point.json or parameter-audit.json")
    args = parser.parse_args()
    root = args.estimator.resolve()
    git = ["git", "-c", "safe.directory=" + str(root), "-c", "core.filemode=false",
           "-c", "core.autocrlf=true", "-C", str(root)]
    revision = subprocess.check_output(git + ["rev-parse", "HEAD"], text=True).strip()
    tracked = subprocess.check_output(git + ["status", "--porcelain", "-uno"], text=True).strip()
    if revision != PIN or tracked:
        raise ValueError("estimator checkout must be clean and pinned at " + PIN)
    sys.path.insert(0, str(root))
    from sage.all import log, oo, sqrt, ZZ
    from sage.version import version
    from estimator import SIS
    from estimator.reduction import ADPS16, MATZOV
    from estimator.sis_lattice import SISLattice

    models = {"core_classical": ADPS16("classical"), "core_quantum": ADPS16("quantum"),
              "matzov_classical": MATZOV(nn="classical"), "matzov_quantum_dw": MATZOV(nn="quantum")}
    def solve_crt(n, m, sq_radius, d):
        params = SIS.Parameters(n=n, m=m, q=Q, length_bound=sqrt(ZZ(sq_radius)), norm=2)
        cost = SISLattice.cost_euclidean(params, d=d, red_cost_model=models["core_classical"])
        return int(cost["beta"]), float(log(cost["rop"], 2))

    ledger = json.loads(args.ledger.read_text(encoding="utf-8-sig"))
    if "representatives" in ledger:
        representatives = ledger["representatives"]
    else:
        names = set(ledger["embedding"].values())
        representatives = {role: (a, b, radius) for role, a, b, radius in ledger["roles"]
                           if role in names}
    if set(representatives) != {"A_0", "A_1", "A_aux", "B_pivot", "D_0", "D_1"}:
        raise ValueError("expected the six representatives of the audited five-layer profile")
    results = {}
    for role, (a, b, radius) in sorted(representatives.items()):
        n, m = DEGREE*a, DEGREE*b
        if not (1 < radius < Q and m > n):
            raise ValueError(f"not a short Euclidean q-ary SIS instance: {role}")
        params = SIS.Parameters(n=n, q=Q, m=m, length_bound=radius, norm=2)
        d = best_gsa_dimension(n, m, radius)
        costs = {}
        for name, model in models.items():
            cost = SISLattice.cost_euclidean(params, d=d, red_cost_model=model)
            if cost["rop"] == oo or not math.isfinite(float(log(cost["rop"], 2))):
                raise ArithmeticError(f"no finite {name} estimate for {role}")
            costs[name] = round(float(log(cost["rop"], 2)), 6)
        default = SIS.estimate(params, red_cost_model=models["core_classical"],
                               quiet=True, jobs=1, catch_exceptions=False)["lattice"]
        if int(default["beta"]) != int(cost["beta"]):
            raise ValueError("dimension optimization changed default block size")
        results[role] = dict(ring_rows=a, ring_columns=b, l2_radius=radius, n=n, m=m,
                             chosen_dimension=d, stock_dimension=int(default["d"]),
                             bkz_block_size=int(cost["beta"]), log2_cost=costs,
                             matzov_within_fit=int(cost["beta"]) <= 1024,
                             single_CRT_component=one_factor_screen(a, b, radius, solve_crt))
    print(json.dumps(dict(estimator_commit=revision, sage_version=version, q=Q, degree=DEGREE,
                          attack="generic coefficient-l2 q-ary SIS, Euclidean GSA",
                          cost_units="ADPS16 Core-SVP exponents; MATZOV estimator rop / quantum DW cost",
                          results=results), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
