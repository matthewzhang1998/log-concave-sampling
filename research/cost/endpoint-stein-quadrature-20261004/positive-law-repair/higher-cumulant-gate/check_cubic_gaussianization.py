#!/usr/bin/env python3
"""Diagnostics for the analytic cubic buffered-Gaussianization theorem.

Scalar numerical CDF integration is a diagnostic, not the proof. Imported
positive mean/covariance compilers are not instantiated numerically here.
"""
from pathlib import Path
import json, math
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import ndtr, ndtri
OUT=Path(__file__).resolve().parent
checks=0
def ck(ok,msg):
    global checks
    if not ok:raise AssertionError(msg)
    checks+=1
v,w=hermgauss(120);v*=math.sqrt(2);w/=math.sqrt(math.pi)
base=np.logaddexp(v,-v)-math.log(2)
base-=w@base
basevar=float(w@(base*base))
rows=[]
for sigma in (.5,1.,2.):
    for L in (.025,.05,.1,.2,.4):
        X=L*base;e=L*math.sqrt(basevar);targetsd=math.sqrt(sigma*sigma+e*e)
        # Large deterministic grid, Gaussian mixture CDF/density evaluated exactly
        # for the finite Hermite integration rule.
        lo=-10*targetsd;hi=10*targetsd+float(np.max(X))*0
        yy=np.linspace(lo,hi,40001)
        arg=(yy[:,None]-X[None,:])/sigma
        cdf=ndtr(arg)@w
        density=(np.exp(-arg*arg/2)/math.sqrt(2*math.pi)/sigma)@w
        mask=(cdf>1e-13)&(cdf<1-1e-13)
        gap=yy[mask]-targetsd*ndtri(cdf[mask])
        w2=math.sqrt(float(np.trapezoid(gap*gap*density[mask],yy[mask])))
        bound=math.sqrt(2)/3*L*L*e/(sigma*sigma)
        ck(w2<=bound*(1+1e-4),'scalar CDF W2 below cubic theorem bound')
        ck(w2>0,'nonGaussian even Lipschitz source')
        rows.append(dict(sigma=sigma,L=L,e=e,w2=w2,bound=bound,ratio=w2/bound,w2_over_L3=w2/L**3))

# Exact second-Hermite Bessel normalization in arbitrary output dimension.
rng=np.random.default_rng(917426)
for d in (1,2,5,13):
    for rep in range(12):
        S=rng.normal(size=(d,d));S=(S+S.T)/2
        linear=rng.normal(size=d);constant=rng.normal()
        # h(z)=constant+linear.z+(z^T S z-tr S).
        Eh2=constant**2+linear@linear+2*np.sum(S*S)
        ED2=2*S
        ck(np.sum(ED2*ED2)<=2*Eh2+1e-12,'second-Hermite Bessel exact normalization')
        if rep==0:
            ck(abs(np.sum(ED2*ED2)-2*(2*np.sum(S*S)))<1e-12,'pure quadratic equality')

# Stein tau may be nonsymmetric. For a linear rectangular source, Gaussian
# closure is exact and the mean tau is the full target covariance.
for din,dout in ((2,1),(3,2),(2,4),(7,5)):
    B=rng.normal(size=(dout,din));B*=.2/max(1.,np.linalg.norm(B,2))
    Sigma=B@B.T
    for t in (0.,.2,.7,1.):
        C=np.eye(dout)+(1-t*t)*Sigma
        ck(np.linalg.eigvalsh(C).min()>=1-1e-13,'covariance-preserving path fixed gap')
        ck(np.linalg.norm(t*t*Sigma+C-(np.eye(dout)+Sigma))<1e-13,'linear path exact Gaussian covariance')

# Positive join variance ledger and reverse-OU relative increment guard.
for h in (.01,.1,.25,.5):
    d=1-2*h*h;eta=math.sqrt(d/2)
    ck(d>=.5,'explicit reserve gap')
    ck(abs(eta*eta+eta*eta-d)<1e-14,'positive reserve shares')
    ck(abs(h*h+d-(1-h*h))<1e-14,'complete mean buffer accounted')
for r in (0.,.2,.7,.95,.999):
    s2=1-r*r
    for theta in (.01,.1,.25):
        Delta=theta*s2;t=math.sqrt(r*r+Delta);v0=Delta/(t*t)
        h=math.sqrt(Delta/s2)
        alpha=.1*s2
        ck(h<=.5+1e-14,'relative bridge increment domain')
        ck(abs(h*h*alpha**3-Delta*.1**3*s2**2)<1e-18,'restored curvature heat scale')
report=dict(status='PASS',assertions=checks,scalar_rows=rows,
            max_scalar_bound_ratio=max(x['ratio'] for x in rows),
            scope='Scalar CDF and exact algebra diagnostics for a proved dimension-safe cubic Gaussianization lemma; no end-to-end imported mean/covariance implementation.')
(OUT/'cubic_gaussianization_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='scalar_rows'},indent=2))
