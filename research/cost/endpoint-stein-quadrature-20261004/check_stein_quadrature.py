#!/usr/bin/env python3
"""Finite diagnostics, not replacements for the uniform analytic estimates."""
import hashlib,json,math,pathlib
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=pathlib.Path(__file__).resolve().parent
checks=0

def check(ok,msg):
    global checks
    assert ok,msg
    checks+=1

def rule(K,m):
    x,w=leggauss(m)
    rr=[]; ww=[]
    for k in range(K):
        a=2.**(-k-1)
        t=1.5*a+.5*a*x
        rr.extend(1-t); ww.extend(.5*a*w)
    rr.append(1-2.**(-K-1)); ww.append(2.**(-K))
    return np.array(rr),np.array(ww)

reports=[]
for K,m in [(4,4),(8,8),(12,12),(20,20),(30,30)]:
    r,w=rule(K,m)
    check(np.all(w>0),'positive weights')
    check(abs(sum(w)-1)<2e-14,'constant exactness')
    check(abs(w@r-.5)<2e-14,'linear exactness')
    degrees=np.unique(np.concatenate([np.arange(200),np.geomspace(1,1e12,2000).astype(np.int64)]))
    errs=[]
    for n in degrees:
        q=w@(r**int(n)); exact=1/(int(n)+1)
        err=abs(q-exact); errs.append(err)
        check(err<=12*(4/9)**m+2**(1-K)+1e-14,'uniform moment analytic bound')
    # Orthogonal-vector chaos diagnostic: no dimension factor beyond input energy.
    rng=np.random.default_rng(K)
    c=rng.normal(size=(75,9))
    e=np.array([w@(r**n)-1/(n+1) for n in range(75)])
    lhs=np.sum((c*e[:,None])**2)
    rhs=np.max(e*e)*np.sum(c*c)
    check(lhs<=rhs+1e-25,'vector Hermite contraction')
    beta=float(w@np.sqrt(1-r*r))
    check(abs(beta-math.pi/4)<.02*2**(-K/2)+1e-12,'bridge coefficient convergence')
    reports.append({'K':K,'m':m,'nodes':len(r),'max_tested_moment_error':max(errs),'bound':12*(4/9)**m+2**(1-K),'beta_Q':beta})

# Explicit rho=3/2 Bernstein ellipse: all powers bounded by one.
rho=1.5
th=np.linspace(0,2*math.pi,10001)
for a in [2.**(-k) for k in range(1,50)]:
    t=1.5*a+(a/4)*(rho*np.exp(1j*th)+rho**-1*np.exp(-1j*th))
    check(np.max(np.abs(1-t))<=1+1e-14,'fixed Bernstein ellipse')

# Noncommuting coordinates are not silently independent: shared G coefficient.
rng=np.random.default_rng(817)
r,w=rule(18,18); beta=float(w@np.sqrt(1-r*r))
for d in [1,2,5,11]:
    O,_=np.linalg.qr(rng.normal(size=(d,d)))
    for A in [.005,.03,.1]:
        B=O@np.diag(np.linspace(A/4,A,d))@O.T
        C=(np.eye(d)-B/2)@(np.eye(d)-B/2).T+beta**2*B@B
        literal=np.eye(d)-B+(.25+beta**2)*B@B
        check(np.linalg.norm(C-literal)<1e-12,'exact common-G covariance')
        wrong=np.eye(d)-B+(.25+float((w*w)@(1-r*r)))*B@B
        check(np.linalg.norm(C-wrong)>1e-9,'independent-G substitution rejected')
        target=np.linalg.inv(np.eye(d)+B)
        leading=(.25+beta**2-1)*(B@B)
        check(np.linalg.norm(C-target-leading)<=1.02*np.linalg.norm(B@B@B),'quadratic second-order debt')

# Finite nonquadratic caller/private derivatives, only original Hessians used.
def gv(x): return .5*x+.15*np.cos(x)
def hv(x): return .5-.15*np.sin(x)
r,w=rule(10,8); s=np.sqrt(1-r*r)
for A in [.001,.015,.08]:
    for y in [-100.,-2.,0.,3.,80.]:
        for Z,G in [(-2.,1.),(0.,0.),(1.3,-.7)]:
            M=9;b=y;db=1.
            for _ in range(M):
                old=b; olddb=db
                b=y-A*gv(old);db=1-A*hv(old)*olddb
            pts=b+math.sqrt(A)*(r*Z+s*G)
            h=math.sqrt(A)*(w@(gv(pts)-gv(b)))
            X=b+math.sqrt(A)*(Z-h)
            Hy=math.sqrt(A)*(w@(hv(pts)-hv(b)))*db
            Xy=db-math.sqrt(A)*Hy
            HZ=A*(w@(hv(pts)*r));HG=A*(w@(hv(pts)*s))
            check(math.hypot(HZ,HG)<=A+1e-13,'actual H private first')
            check(abs(Xy-1)<=3*A,'actual caller first')
            check(abs(db-1)<=A/(1-A)+1e-14,'mode first')
            if Z==G==0:
                check(h==0 and X==b,'literal nonlinear zero')
            # Original numerical mode residual compared with initial gradient.
            ell=(b-y+A*gv(b))/math.sqrt(A)
            bound=(1+A)*A**(M+.5)*abs(gv(y))/(1-A)
            check(abs(ell)<=bound+3e-12,'finite mode linear target residual')

result={'scope':'Bounded quadrature/Hermite/finite-ports diagnostics; no all-order law theorem.',
        'assertions':checks,'rules':reports,
        'source_sha256':hashlib.sha256((ROOT/'ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md').read_bytes()).hexdigest(),
        'quadratic_limit_B2_coefficient':.25+math.pi**2/16,
        'all_passed':True}
(ROOT/'stein_quadrature_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
