#!/usr/bin/env python3
"""Diagnostic checks only. No native compiler or infinite OU path is executed."""
import json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from fractions import Fraction as F
rng=np.random.default_rng(20261005)
checks=0

def check(test, label):
    global checks
    checks += 1
    if not test: raise AssertionError(label)

def bridge(d,s):
    a=np.sinh(d-s)/np.sinh(d)
    b=np.sinh(s)/np.sinh(d)
    v=(-np.expm1(-2*s))*(-np.expm1(-2*(d-s)))/(-np.expm1(-2*d))
    return a,b,np.maximum(0,v)

def rule(d,eps):
    w=-np.expm1(-d); lo=eps*w/64
    panels=[]; a=lo
    while a<d/2:
        b=min(2*a,d/2); panels.append((a,b)); a=b
    panels += [(d-b,d-a) for a,b in panels[::-1]]
    n=8+math.ceil(math.log(1/eps,4)); x,b=leggauss(n)
    sites=[0.,d]; weights=[-np.expm1(-lo),np.exp(-d)*np.expm1(lo)]
    for a,c in panels:
        ss=(a+c)/2+(c-a)*x/2
        sites.extend(ss); weights.extend((c-a)*b*np.exp(-ss)/2)
    sites=np.asarray(sites); weights=np.asarray(weights); weights*=w/weights.sum()
    return sites,weights,len(panels),n

for A in np.geomspace(1e-8,.1,25):
    w=A**.6; eta=A; h=A**1.4; q=1-w; d=-math.log(q); sig2=1-q*q
    for s in np.linspace(0,d,61):
        a,b,v=bridge(d,s)
        check(a>=0 and b>=0 and v>=0,'positive bridge rows')
        check(abs(a+q*b-math.exp(-s))<2e-12,'standardized row1')
        check(abs(a*a+b*b+2*q*a*b+v-1)<2e-12,'bridge marginal variance')
        check(v<=w/(2-w)+1e-12,'bridge max variance')
        for yy in [-.99,-.3,0,.3,.99]:
            y=yy*min(s,d-s); z=s+1j*y
            l=np.array([np.exp(-z),2*q*np.sinh(z)/math.sqrt(sig2)])
            check(np.vdot(l,l).real<=1+3e-12,'complex contraction wedge')
    # Complete true bridge covariance: positive semidefinite and row variances.
    ss=np.linspace(0,d,19)
    C=np.exp(-np.abs(ss[:,None]-ss[None,:]))
    endpoints=np.exp(-np.column_stack((ss,d-ss)))
    E=np.array([[1,q],[q,1]])
    C-=endpoints@np.linalg.solve(E,endpoints.T)
    check(np.linalg.eigvalsh(C).min()>-2e-10,'conditional bridge covariance PSD')
    check(np.max(np.abs(np.diag(C)-bridge(d,ss)[2]))<2e-10,'bridge covariance diagonal')
    # Positive common-root quadrature, including endpoints with zero bridge row.
    ss,bb,panels,n=rule(d,max(A*A,1e-10)); aa,cc,vv=bridge(d,ss)
    check(np.all(bb>0),'positive quadrature weights')
    check(abs(bb.sum()-w)<1e-13,'exact normalized mass')
    check(np.dot(bb,aa)<=w+1e-13,'partial-X contraction')
    check(np.dot(bb,np.sqrt(vv))<=w**1.5+1e-13,'whole common-root first')
    check(panels<=2*(math.ceil(math.log2(64*d/(max(A*A,1e-10)*w)))+1),'polylog panels')
    # First few chaos components of the conditional-mean operator.
    l1=np.exp(-ss); l2=2*q*np.sinh(ss)/math.sqrt(sig2)
    for degree in [0,1,2,3,7,13]:
        for k in range(degree+1):
            fac=math.sqrt(math.comb(degree,k))
            approx=np.dot(bb,fac*l1**k*l2**(degree-k))
            ref=quad(lambda s: math.exp(-s)*fac*math.exp(-s*k)*(2*q*math.sinh(s)/math.sqrt(sig2))**(degree-k),0,d,epsabs=1e-13)[0]
            check(abs(approx-ref)<max(A*A,1e-10)*w+1e-12,'finite bridge chaos quadrature')
    # Independent branch cutoff and exact disintegration.
    for t in np.linspace(0,1-eta,101):
        v=(1-t*t)*sig2/(1-q*q*t*t)
        check(abs(1/v-(q*q/sig2+1/(1-t*t)))<2e-7*max(1,1/v),'reciprocal identity')
        check(v>=eta/2,'actual minimum bulk buffer')
        check(v<=2*w,'actual maximum bulk buffer')
        for share in [.5,.25,.25]:
            u=share*v
            check(A**1.5/math.sqrt(u)<=math.sqrt(8)*A*(1+1e-12),'residual-first bound')
    # Known terminal carrier alignment is an orthogonal row rotation.
    for t in np.linspace(.0001,.9999,37):
        c=math.sqrt(1-t*t); r=math.sqrt(1-h*h)
        norm=math.sqrt(r*r*c*c+h*h)
        O=np.array([[r*c,h],[-h,r*c]])/norm
        check(np.linalg.norm(O@O.T-np.eye(2))<1e-12,'carrier rotation orthogonal')
        check(abs((r*t)**2+norm*norm-1)<1e-12,'terminal carrier OU variance')

# Exact substantive exponent ledger, eta=A, w=A^(3/5), h=A^(7/5).
a=F(3,5); e=F(1); k=F(7,5)
ledger={'smooth_commutator':2+k,'smooth_drift':1+2*k,'nested_prefix':3+a,'bridge_noise':3+3*a-k,'gaussianization_bulk':4-a,'near_sqrt':3+e/2,'near_other':3+e-a/2,'mean_bulk':5-3*a/2,'mean_endpoint':5-e/2,'K_bulk':6-7*a/6,'K_endpoint':6-e/6}
check(min(ledger.values())==F(17,5),'minimum substantive exponent')
check(ledger['smooth_commutator']==ledger['bridge_noise']==ledger['gaussianization_bulk'],'three-way exponent balance')

result={'status':'PASS','assertions':checks,'scope':'Bridge geometry, complex contraction, positive quadrature low-chaos diagnostics, true bridge covariance, actual minimum buffers, carrier rotations and exact exponent ledger. No full native compiler execution.','exponents':{k:str(v) for k,v in ledger.items()}}
Path(__file__).with_name('smoothed_bridge_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
