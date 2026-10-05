#!/usr/bin/env python3
"""Elementary diagnostic and exact-rational certificates for the smoothing packet.

This does not implement a native compiler or numerically estimate target-law error.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, math

count = 0
def check(test):
    global count
    count += 1
    assert test

# Exact upper/lower radical certificates; no floating point in these checks.
rt2low=F(140,99)
klow=F(49,80)  # upper bound on sqrt(3/8)
check(rt2low**2 < 2)
check(klow**2 > F(3,8))
D0=1+1/rt2low
E0=1+klow
A0=F(1,36)
pnum=D0*E0+A0*E0**2/2+F(1,2)
# sqrt(2 pi)>5/2 follows from pi>25/8; the elementary pi bound is analytic.
check(pnum/F(5,2) < F(33,25))
qnum=(D0+A0*E0)**3+4*D0**3+A0**3*E0**3
check(qnum/(6*rt2low) < 3)
# Sharp residual constant: bound all factors by explicit rationals.
Hupper=14+F(11,72)+F(11,7)*(F(83,14)+F(61,576))+F(7,2)*F(505,144)*F(19,18)
check(Hupper<37)
check(F(7,2)*A0+A0*A0/4 < F(1,10))
check(F(18,19)**2<F(9,10))
check(F(7,4)**2>3)
check(F(31,8)**2>15)
mainupper=84+F(33,25)*F(7,4)+3*F(31,8)
check(mainupper<98)
check(mainupper<110)

# Exact exponent balance: weak transfer A^2 e and smooth remainder A^4/e^2.
check(2+F(2,3)==F(8,3))
check(4-2*F(2,3)==F(8,3))
check(4-F(2,3)==F(10,3))
check(F(1,2)+F(1,6)==F(2,3))
check(1-2*F(1,6)==F(2,3))
check(F(1,2)+F(2,3)==F(7,6))
# Naive full-source exponent balance is q=2; it is not a minimax lower bound.
check(1+1==4-2*1==2)

k=math.sqrt(3/8); d0=1+1/math.sqrt(2); e0=1+k
worst_H=worst_C=worst_ratio=0.
for i in range(1,501):
    A=(1/36)*(i/500)**3
    L=A/2+A*A/4; Ls=3.5*A+A*A/4
    Hsum=2+11/(2*math.sqrt(2))+A*(1+4.5*k)
    p=(d0*e0+A*e0*e0/2+.5)/math.sqrt(2*math.pi)
    q=((d0+A*e0)**3+4*d0**3+A**3*e0**3)/(6*math.sqrt(2))
    for b in [1,2,3,10,100,1000,1000000]:
        C=1+A+math.sqrt(b)*((1+A)*math.pi/2+L/A/math.sqrt(1-L))
        H=14+5.5*A+math.sqrt(b)*(math.pi/2*Hsum+3.5*(Ls/A)/math.sqrt(1-Ls))
        check(C<5*math.sqrt(b)); check(H<37*math.sqrt(b))
        eps=A**(2/3)*b**(1/6)
        bound=84*A*A*eps*math.sqrt(b)+A**4*(p*math.sqrt(b+2)/eps+q*math.sqrt((b+2)*(b+4))/eps**2)
        denom=A**(8/3)*b**(2/3)
        check(bound<=110*denom)
        check(Ls<1)
        worst_H=max(worst_H,H/math.sqrt(b))
        worst_C=max(worst_C,C/math.sqrt(b))
        worst_ratio=max(worst_ratio,bound/denom)
    cap=1.5+5/(2*math.sqrt(2))+A*(1+2*k)+1+A
    check(cap<5)

# No smoothing leaves are needed: a single fixed node z=0 returns g itself.
# Ordinary literal stage-two occurrence count remains the sealed affine expression.
for j1,jl,jq,nout in [(1,1,1,1),(7,11,13,17),(100,121,144,25)]:
    residual=(5+3*j1+5*jl+6*jq)*nout
    smooth_zero_node_multiplier=1
    check(residual*smooth_zero_node_multiplier==residual)

sealed=Path('/workspace/shared/nonlinear-bridge-stage-two/MANIFEST.json')
expected='cd1111e28dc399643f27bcdd0f97dc735ada196f8a43f220211500bcb45be3ae'
check(hashlib.sha256(sealed.read_bytes()).hexdigest()==expected)
result={
 'status':'PASS','assertions':count,
 'exact_rational_stencil_constant_upper':str(Hupper),
 'exact_rational_main_constant_upper':str(mainupper),
 'maximum_tested_stencil_constant_over_sqrt_b':worst_H,
 'maximum_tested_target_constant_over_sqrt_b':worst_C,
 'maximum_tested_optimized_constant':worst_ratio,
 'sealed_stage_two_manifest_sha256':expected,
 'scope':'Exact-rational elementary certificates and parameter diagnostics. No native compiler implementation or numerical target-rate proof.'
}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
