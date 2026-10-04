"""Finite algebra/normalization diagnostics, not a replacement for source proofs."""
from fractions import Fraction as F
from math import sqrt, ceil, sin
import hashlib, json
from pathlib import Path
out=Path(__file__).parent
checks=0
J=F(29,10); chi=max(J-2,(J-1)/3); vp=2*(J-1)/3
assert chi==F(9,10) and vp==F(19,15); checks+=2
assert 1+F(3,2)*vp==J and 1+vp==F(34,15); checks+=2
assert 1-chi==F(1,10) and F(1,2)-chi==-F(2,5); checks+=2
assert 1+J==F(39,10) and 2+J==F(49,10); checks+=2
assert F(34,15)<F(14,3) and F(44,15)<F(14,3); checks+=2
for Q in [F(q,10) for q in range(1,151)]:
    for kap in [F(0),F(1,2),F(1),F(321,100),F(161142,10000)]:
        m=(10*Q+5*kap).__floor__()+1
        assert F(m,10)-kap/2>Q; checks+=1
        k=max(F(0),Q-F(49,10))+F(1,100)
        assert F(49,10)+k>Q; checks+=1
for ell in [.001,.01,.1,.2,.5]:
    for v in [.1,.5,1.0]:
        for N in [1,2,3,10,101,10001]:
            exact=sqrt(v+ell*ell/N)-sqrt(v)
            bound=ell*ell/(2*sqrt(v)*N)
            assert exact<=bound+1e-13; checks+=1
for c in [F(-1),F(-3,5),F(0),F(1,7),F(1)]:
    v=F(1,2)
    assert c*c*v+v*(1-c*c)==v; checks+=1
    assert v+v==1; checks+=1
# Exact finite directed centering identity tested through varying origins.
for theta in [-1.,-.3,0.,.25,1.2]:
    for w in [-3.,-.7,0.,.4,2.]:
        p0=lambda x:.4*x+.08*sin(x)
        a=.07
        p1=lambda x:.3*x+.05*sin(2*x)
        raw0=p0(theta+w)
        raw1=p1(.6*theta+.5*w+a*raw0)
        zero0=p0(theta); q10=.6*theta+a*zero0
        zero1=p1(q10)
        z0=p0(theta+w)-zero0
        z1=p1(q10+.5*w+a*z0)-p1(q10)
        assert abs(z0-(raw0-zero0))<1e-14
        assert abs(z1-(raw1-zero1))<1e-14
        checks+=2
# Finite centered twin: both signs remain zero at zero input at every depth.
for K in range(1,65):
    pp=pm=0.
    for k in range(K):
        a=.04
        pp,pm=sin(a*pp+a*(pp+pm)),-sin(a*pp)
    assert pp==pm==0.; checks+=1
source=out/'WEAK-29-ACTUAL-REMAINDER-AND-RETAINED-MEAN-RECONSTRUCTION.md'
data={'status':'PASS','checks':checks,'scope':'finite scalar normalization, rational grades, Gaussian linear sharpness, exact centering and finite zero iteration only','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
(out/'weak_29_reconstruction_checks.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
