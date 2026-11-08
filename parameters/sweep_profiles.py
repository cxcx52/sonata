"""Freeze and verify the four ordinary-transition profiles under common targets."""
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path

import audit_profile as audit
import compute
from estimate_sis import best_gsa_dimension


ROOT = Path(__file__).resolve().parent
HISTORY = 16384
HISTORY_BITS = 130
CLASSICAL_BITS = 144
QUANTUM_SCREEN_BITS = 136
ROOT_BLOCKS = (4, 6, 8, 10, 12, 16, 20, 21, 24)
INNER_BLOCKS = (3, 4, 6, 8, 10)
RADICES = (16, 32, 64, 128)
ARITIES = (16, 64, 256, 1024)
SIGMA = compute.SIGMA_SQUARED
M = compute.CHALLENGE_COUNT_LOWER
Q = compute.Q
K = Q**4


def ceilroot(num, den=1):
    floor = math.isqrt(num//den)
    return floor + (floor*floor*den < num)


@lru_cache(None)
def digit_energy(radix):
    return compute.max_digit_energy((Q-1)//2, radix)


def role_envelopes(roles):
    assignment = {"A_st": "A_aux", "A_aux": "A_aux", "B_pivot": "B_pivot"}
    for i in range(5):
        assignment[f"A_{i}"] = "A_0" if i == 0 else "A_1"
        if i < 4:
            assignment[f"B_{i}"] = assignment[f"D_{i}"] = "D_0" if i == 0 else "D_1"
    envelopes = {}
    for name, rep in assignment.items():
        rows, cols, beta = roles[name]
        prev_rows, prev_cols, prev_beta = envelopes.get(rep, (rows, 0, 0))
        assert prev_rows == rows
        envelopes[rep] = (rows, max(cols, prev_cols), max(beta, prev_beta))
    return assignment, envelopes


def satisfies_attack_floor(role):
    rows, cols, radius = role
    n, m = 64*rows, 64*cols
    if not (1 < radius < Q and n < m):
        return False
    d = best_gsa_dimension(n, m, radius)
    permitted = (math.log(radius) - n*math.log(Q)/d)/(d-1)
    # Requiring failure even at the first acceptable block conservatively
    # leaves one extra block of headroom before the pinned-estimator check.
    b = max(math.ceil(CLASSICAL_BITS/0.292), math.ceil(QUANTUM_SCREEN_BITS/0.265))
    at_b = (math.log(b/(2*math.pi*math.e)) + math.log(math.pi*b)/b)/(2*(b-1))
    return permitted < at_b


def derive(k, blocks, radices):
    V = 64*k + 544
    raw, S, m = 55*V+1, 47*V+1, 256*(k+3)+5
    root_beta = ceilroot(4*38400*S, 4121)
    roles = {"A_st": (4,219,root_beta), "A_aux": (4,compute.ceil_div(7*V,64),root_beta),
             "B_pivot": (2,256,1024)}
    layers = []
    for i,(s,rho) in enumerate(zip(blocks, radices)):
        n = compute.ceil_div(raw, 64*s)
        N = 64*n*s
        G = 256*S - int(bool(rho))
        layer = dict(s=s,n=n,rho=rho,raw=raw,N=N,S=S,G=G,rows=m+N-raw)
        if rho:
            ell = compute.ceil_div(48,rho.bit_length()-1)
            t = 16*s
            h = 2*s*(s+1)
            bt,bh = ell*t,ell*h
            low = 64*n*(rho//2)**2
            high = compute.ceil_surd(G+low,4*G*low,rho*rho)
            next_S = low + high + 64*(t+h)*digit_energy(rho)
            tag_beta = ceilroot(4*38400*next_S, 4121)
            roles[f"B_{i}"], roles[f"D_{i}"] = (6,bt,tag_beta), (6,bh,tag_beta)
            z_num,z_den = (1+rho*rho)*38400*next_S,4121
            raw,S,m = 64*(2*n+bt+bh),next_S,2052
        else:
            z_num,z_den = G,1
        cross = ceilroot(64*26**2*z_num,z_den)
        anchor = ceilroot(z_num,z_den) + ceilroot(26**2*s*38400*layer["S"],4121)
        assert anchor <= cross < Q
        layer.update(cross_radius=cross,anchor_radius_upper=anchor)
        roles[f"A_{i}"] = (16,n,cross)
        layers.append(layer)
    return layers, roles


def schedules():
    for s0 in ROOT_BLOCKS:
        for s1,s2,s3 in product(INNER_BLOCKS,repeat=3):
            blocks = (s0,s1,s2,s3,5)
            # Necessary coordinate-only screen, exact integer arithmetic.
            if HISTORY*320*sum(blocks)*(1<<HISTORY_BITS) > M:
                continue
            for r0,rtail in product(RADICES,repeat=2):
                yield blocks,(r0,rtail,rtail,rtail,0)


def select(k):
    best = None
    examined = 0
    eligible = 0
    for blocks,rads in schedules():
        examined += 1
        layers,roles = derive(k,blocks,rads)
        _, reps = role_envelopes(roles)
        if not all(satisfies_attack_floor(role) for role in reps.values()):
            continue
        eligible += 1
        # The leading objective measures the size of all five expanded matrices.
        key = (sum(v["N"] for v in layers), max(a*b for a,b,_ in roles.values()),
               blocks, rads)
        if best is None or key < best:
            best = key
    if best is None:
        raise ValueError(f"no feasible schedule for {k}")
    return best,examined,eligible


def audit_profile(k, blocks, radices):
    layers,roles = derive(k,blocks,radices)
    assignment,reps = role_envelopes(roles)
    V = 64*k+544
    T = 55*V+64*(k+2)+2
    ell = (T-1).bit_length()
    no_wrap = math.ceil(956*SIGMA*layers[0]["S"]+3281*V)
    assert no_wrap < Q and audit.max_group_weight() == 47
    # Continuity and pre-gamma comparisons use differences of physical
    # canonical bit openings, including fresh reserved slots.
    assert roles["A_st"][2]**2 >= max(73*4*48,64*48)
    for name,(rows,cols,radius) in roles.items():
        rep_rows,rep_cols,rep_radius = reps[assignment[name]]
        assert rows==rep_rows and cols<=rep_cols and radius<=rep_radius<Q
    assert max(a*b for a,b,_ in roles.values()) == max(a*b for a,b,_ in reps.values())
    ball,_,support = audit.challenge_support()
    assert support == M
    rank_bytes = compute.ceil_div((ball-1).bit_length(),8)
    prover_bytes = 24*(2*ell+1+6) + 8*384
    verifier_bytes = 24*(ell+4)
    backend = Fraction()
    abort = Fraction()
    for i,v in enumerate(layers):
        s,n,rho,N,S,G = [v[key] for key in ("s","n","rho","N","S","G")]
        b,h,seed_bytes,delta = audit.projection_parameters(N,S,150 if i<2 else 151)
        seed,other_delta = compute.certify_seed(N,S,150 if i<2 else 151)
        assert (b,h,seed_bytes,delta) == (seed["b"],seed["h"],seed["bytes"],other_delta)
        degree = compute.ceil_div(v["rows"]+384,4)-1
        v.update(b=b,h=h,seed_bytes=seed_bytes,
                 transfer_bits=seed["transition_matrix_bits"],degree=degree)
        # Match the ordinary wire layout. One control byte per issued attempt;
        # shared old/fresh bodies, CRS and transport headers are excluded.
        prover_bytes += 2 + audit.rice_bytes(384,384*S)
        verifier_bytes += seed_bytes + 24 + rank_bytes*s
        if rho:
            prover_bytes += 12*384
        else:
            prover_bytes += (16*(s-1)+2+2*s*(s+1))*384 - 4*6 + audit.rice_bytes(64*n,G)
        assert Q**2 >= 176**2*SIGMA*S
        coord = Fraction(s-1,M)+Fraction(1,M-1) if not rho else Fraction(s,M)
        backend += 160*(Fraction(1,2**192)+delta)+Fraction(degree,K)+320*coord
        # (1+2*delta)^160 <= (1-320*delta)^(-1) by binomial series.
        abort += Fraction(1,2**160)*(1 + 1/(1-320*delta))
        if rho:
            ds = compute.ceil_div(48,rho.bit_length()-1)
            low = 64*n*(rho//2)**2
            nxt = layers[i+1]
            assert nxt["S"] == (low + audit.high_energy(G,low,rho)
                + 64*(16*s+2*s*(s+1))*audit.max_digit_energy(rho))
            assert nxt["raw"] == 64*(2*n+ds*(16*s+2*s*(s+1)))
        else:
            assert G == compute.ceil_div(256*M*S,M-1)-1
    auxiliary = Fraction(72,K-1)+Fraction(2*ell,K)
    zero = backend+auxiliary+Fraction(4*k,K-1)
    assert abort<Fraction(11,2**160)
    assert HISTORY*zero<=Fraction(1,1<<HISTORY_BITS)
    return dict(k=k,V=V,T=T,ell=ell,no_wrap_upper=no_wrap,layers=layers,
                flat_CRS_cells=max(a*b for a,b,_ in roles.values()),
                roles=[[g,*r] for g,r in roles.items()],embedding=assignment,
                representatives=reps,backend_exact=str(backend),transition_exact=str(backend+auxiliary),
                zero_pending_exact=str(zero),issued_history_bits=compute.bit_bound(HISTORY*zero),
                zero_pending_bits=compute.bit_bound(zero),honest_abort_bits_lower=compute.bit_bound(abort),
                no_retry_bytes=dict(prover=prover_bytes,verifier=verifier_bytes,
                                    total=prover_bytes+verifier_bytes),
                projection_entries=384*sum(v["N"] for v in layers))


def generate():
    profiles = {}
    for k in ARITIES:
        (_,_crs,blocks,rads),examined,eligible = select(k)
        profile = audit_profile(k,blocks,rads)
        profile["search_candidates_after_statistic_screen"] = examined
        profile["search_candidates_after_sis_screen"] = eligible
        profiles[str(k)] = profile
    return dict(target=dict(classical_core_svp_min=CLASSICAL_BITS,
                            quantum_core_svp_sensitivity_min=QUANTUM_SCREEN_BITS,
                            issued_attempt_limit=HISTORY,history_statistical_bits_min=HISTORY_BITS),
                profiles=profiles)


def header(profiles):
    text = ["#pragma once", '#include "sonata/backend.hpp"', "", "namespace sonata::example {",
            "struct SweepProfile { size_t k, crs_cells; std::array<Config,5> layers; };",
            "inline constexpr std::array<SweepProfile,4> sweep_profiles{{"]
    for p in profiles.values():
        text.append("  {" + f"{p['k']}, {p['flat_CRS_cells']}, " + "{{")
        for v in p["layers"]:
            vals = [str(v[key]) for key in ("s","n","rho","raw","N","S","G","b","h")]
            text.append("    {"+", ".join(vals)+"},")
        text.append("  }}},")
    text += ["}};","} // namespace sonata::example",""]
    return "\n".join(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="replace generated JSON and C++ header")
    mode.add_argument("--check", action="store_true", help="check both frozen outputs (default)")
    args = parser.parse_args()
    if not __debug__:
        parser.error("run the audit with Python assertions enabled")
    data = generate()
    files = {ROOT/"sweep-profiles.json": json.dumps(data,indent=2)+"\n",
             ROOT.parent/"examples/sweep_profiles.hpp": header(data["profiles"])}
    for path, contents in files.items():
        if args.write:
            path.write_text(contents,encoding="utf-8")
        elif path.read_text(encoding="utf-8") != contents:
            raise ValueError(f"{path.name} does not match the independent recomputation")
    for k,p in data["profiles"].items():
        print(k, [v["s"] for v in p["layers"]],[v["rho"] for v in p["layers"]],
              p["projection_entries"],p["issued_history_bits"])


if __name__ == "__main__":
    main()
