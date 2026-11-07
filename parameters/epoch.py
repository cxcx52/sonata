"""Check the epoch embedding and its statistical ledger at the frozen k=1024 profile."""
from fractions import Fraction
from pathlib import Path
import json
import math
import sys

ROOT = Path(__file__).resolve().parent
BASE = json.loads((ROOT / "ledger-two-point.json").read_text(encoding="utf-8"))
Q = 2**48 - 59
K = Q**4
OLD = BASE["layers"][0]
N, RAW, S = (OLD[key] for key in ("N", "raw", "S"))
ROW_LIMIT = 264036
ORIGINAL_PADDING = N - RAW
VALUE_LIMIT = (RAW - 1)//55


def logbits(p):
    return math.log2(p.denominator) - math.log2(p.numerator)


def measure(h,n):
    assert 1 <= h <= n
    V = 64*n + 272*(h+1)
    raw = 55*V + 1
    rows = 256*(n+h+2)+5
    T = raw+64*(n+h+1)+1
    ell = (T-1).bit_length()
    J = 4*n + 69*h + 4
    direct_J = J-h
    backend = Fraction(BASE["backend_exact"])
    # All semantic terms have positive powers and rho is nonzero.
    def error(equations, two_coins):
        roots = max(1,equations-1) if two_coins else equations
        return backend + Fraction(roots,K-1) + Fraction(2*ell,K)
    return dict(h=h,n=n,V=V,raw=raw,T=T,ell=ell,J=J,direct_J=direct_J,
                natural_rows=rows,ignored_coefficients=RAW-raw,
                zero_rows=ROW_LIMIT-ORIGINAL_PADDING-rows,
                strict_legacy_rows=rows+N-raw,
                aux_ring_columns=(7*V+63)//64,
                honest_energy=47*V+1,
                merged_bits=logbits(error(J,True)),
                direct_bits=logbits(error(direct_J,True)),
                direct_one_coin_bits=logbits(error(direct_J,False)))


def main():
    assert (N, RAW, S, VALUE_LIMIT, ORIGINAL_PADDING)==(3635520,3634401,3105761,66080,1119)
    assert math.ceil(956*Fraction(38400*S,4121)+3281*VALUE_LIMIT)==27883328406<Q
    maximum_h=(VALUE_LIMIT-272)//336
    shapes=0
    max_j=0
    max_ell=0
    for h in range(1,maximum_h+1):
        nmax=(VALUE_LIMIT-272*(h+1))//64
        # All monotone capacity conditions attain their maximum at nmax.
        v=measure(h,nmax)
        assert v["zero_rows"] >= 0 and v["ignored_coefficients"] >= 0
        assert v["aux_ring_columns"]<=7228 and v["honest_energy"]<=S
        assert 64*(nmax+1)+272*(h+1)>VALUE_LIMIT
        shapes+=nmax-h+1
        max_j=max(max_j,v["J"])
        max_ell=max(max_ell,v["ell"])
    assert maximum_h==195 and 64*196+272*197>VALUE_LIMIT
    edge=measure(2,1019)
    assert edge['strict_legacy_rows']==265652 and edge['zero_rows']==1024
    eight=measure(8,994)
    assert (eight['V'],eight['T'],eight['ell'],eight['J'],eight['direct_J'])==(66064,3697714,22,4532,4524)
    assert (eight['ignored_coefficients'],eight['zero_rows'])==(880,5888)
    layers=BASE['layers']
    # Round up each transport error with the integer per-layer certificate.
    # Bound (1/2 + d)^160 <= 2^-160/(1-320d) since (1+x)^160 <= 1/(1-160x).
    # The rational relaxation is used only to avoid constructing enormous powers.
    transport=[Fraction(1,2**v['seed']['transition_matrix_bits']) for v in layers]
    abort=sum((1+1/(1-320*d))*Fraction(1,2**160) for d in transport)
    assert abort < Fraction(11,2**160)
    uniform=Fraction(BASE['backend_exact'])+Fraction(max_j-1,K-1)+Fraction(2*max_ell,K)
    output=dict(maximum_nonempty_steps=maximum_h,permitted_pairs_h_n=shapes,
                original_padding=ORIGINAL_PADDING,fixed_aligned_rows=ROW_LIMIT,
                fixed_aggregation_degree=(ROW_LIMIT+384+3)//4-1,
                statistical_bits_exclude="six MSIS advantages and global time cutoff",
                honest_abort_bound_bits=logbits(abort),
                max_J=max_j,max_ell=max_ell,uniform_per_attempt_exact=str(uniform),
                any_profile_2048_attempts_bits=logbits(2048*uniform),
                cases=[measure(h,(VALUE_LIMIT-272*(h+1))//64) for h in (1,2,4,8,16,195)])
    text=json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)
    if "--check" in sys.argv:
        assert output==json.loads((ROOT/"epoch-capacity.json").read_text(encoding="utf-8"))
        print(f"PASS {shapes} feasible (h,n) pairs via monotone envelope; 8x994 and 2x1019 verified")
    else:
        print(text)


if __name__=="__main__":
    main()
