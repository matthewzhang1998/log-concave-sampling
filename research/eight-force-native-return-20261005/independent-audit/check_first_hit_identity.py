#!/usr/bin/env python3
"""Exact polynomial diagnostic of first-hit stopping. Not a native-source witness."""
import json
from pathlib import Path
from math import prod
import sympy as S
P=Path(__file__).resolve().parent
p,q,r,g,a=S.symbols('p q r g a')
variables=(p,q,r,g)
def gm(k):
    if k%2:return 0
    return prod(range(1,k,2))
def E(f):
    pol=S.Poly(S.expand(f),*variables)
    return S.expand(sum(coef*prod(gm(k) for k in powers) for powers,coef in pol.terms()))
R=a*(p*q+q*q*r+p**3)/10
W=p+2*q+3*r+R+g
R1,R2,R3=[S.diff(R,v) for v in (p,q,r)]
x=S.symbols('x')
results=[]
for n in range(1,8):
    phi=x**n
    der=lambda k:S.diff(phi,x,k).subs(x,W)
    lhs=E(p*q*r*der(1)-6*der(4))
    # Constant rows are 1,2,3: after first success 1, after second success 2.
    stopped=E(q*r*R1*der(2)+r*R2*der(3)+2*R3*der(4))
    keep=E((q*r*R1*g+r*R2*(g*g-1)+2*R3*(g**3-3*g))*der(1))
    assert S.expand(lhs-stopped)==0
    assert S.expand(lhs-keep)==0
    assert lhs==0 or S.Poly(lhs,a).terms()[-1][0][0]>=1
    results.append({'test_degree':n,'difference':str(lhs),'first_hit_identity':True,'keep_transfer':True})
# Source-zero old-bank readout is substantive. With F=(q²-1), the omitted
# constant-row term already gives E[F W²]=2 when W=q+g.
q,g=S.symbols('q g')
assert E((q*q-1)*(q+g)**2)==2
out={'scope':'Exact polynomial Gaussian integration-by-parts diagnostic, not native program execution or a source-qualified counterexample.', 'tests':results,'nonzero_bank_readout_warning':{'F':'q^2-1','W':'q+g','test_gradient':'W^2','uncentered_sourcezero_contribution':2},'higher_derivative_of_R_used':False}
(P/'first_hit_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
