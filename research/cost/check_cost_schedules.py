from fractions import Fraction as F
import json
from pathlib import Path

GAMMA = F(49, 10)
BETA = GAMMA + F(1, 2)
P0 = F(7)

def stage_exponent(s, c, R):
    return s*c + max(F(0), R-BETA*s)

def endpoint_bound(c, P, S):
    R=P+F(1,2)
    return max(stage_exponent(F(1),c,R),stage_exponent(S,c,R))

# The stage work exponent is convex piecewise linear. Its maximum on a
# closed interval is exactly at an endpoint, including a crossing at N=1.
checks=0
for P2 in range(14,81):
    P=F(P2,2); R=P+F(1,2); T=2*R/3
    for c in [F(0),F(21,10),F(47,10),F(27,5),F(39,5),F(10),F(100)]:
        bound=endpoint_bound(c,P,T)
        assert bound==max(c+P-GAMMA,T*c)
        for k in range(201):
            s=1+(T-1)*F(k,200)
            assert stage_exponent(s,c,R)<=bound
            checks+=1
        S=R/P
        hybrid=endpoint_bound(c,P,S)
        assert hybrid==max(c+P-GAMMA,(R/P)*(c+P-BETA))
        assert (R/P)*c<=hybrid
        checks+=2

uniform=adapt=hybrid=F(0)
rows=[dict(P='7',uniform='0',adaptive_N='0',conditional_hybrid='0',nested_prescribed_count='0')]
nested=F(0)
for j in range(1,27):
    P=P0+F(j-1,2); R=P+F(1,2); T=2*R/3
    uniform=T*(uniform+P-GAMMA)
    adapt=max(adapt+P-GAMMA,T*adapt)
    # First generation uses only the audited N-allocation change so its
    # actual derivative exponent becomes g=1. Later suffixes are conditional.
    hybrid=(max(hybrid+P-GAMMA,T*hybrid) if j==1 else
            max(hybrid+P-GAMMA,(R/P)*(hybrid+P-BETA)))
    nested+=P-GAMMA
    assert nested==F(21*j,10)+F(j*(j-1),4)
    rows.append(dict(P=str(R),uniform=str(uniform),adaptive_N=str(adapt),
                     conditional_hybrid=str(hybrid),nested_prescribed_count=str(nested),
                     adaptive_N_ratio=str(adapt/(2*R)),hybrid_ratio=str(hybrid/(2*R))))
    checks+=1
out=dict(status='PASS',checks=checks,meaning='Rational cost algebra only; no substitute for source/host proof.',rows=rows)
Path(__file__).with_name('cost_schedule_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {checks} rational checks')
for x in rows[:11]:
    print(x['P'],*[f'{float(F(x[k])):.8g}' for k in ['uniform','adaptive_N','conditional_hybrid','nested_prescribed_count']])
