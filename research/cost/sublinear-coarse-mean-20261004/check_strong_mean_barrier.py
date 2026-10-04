#!/usr/bin/env python3
"""Finite diagnostics; the information proof is in the accompanying note."""
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
checks = 0
def check(ok, label):
    global checks
    if not ok:
        raise AssertionError(label)
    checks += 1

def phi(t):
    return t*t*(1-t)*(1-t) if 0 < t < 1 else 0.
def dphi(t):
    return 2*t-6*t*t+4*t*t*t if 0 < t < 1 else 0.
def Phi(t):
    if t <= 0: return 0.
    if t >= 1: return 1/30
    return t**3/3-t**4/2+t**5/5

def posterior(A, theta, eta=.25, precision=None):
    n=len(theta); w=1/n
    s=1+.5*A if precision is None else precision
    def features(u):
        return np.array([Phi((u-1-i*w)/w) for i in range(n)])
    def raw(u):
        return math.exp(-s*u*u/2-eta*A*w*w*float(np.dot(theta,features(u))))
    points=[-10.,0.]+[1+i*w for i in range(n+1)]+[3.,10.]
    def integral(f):
        return sum(quad(lambda u:f(u)*raw(u),a,b,epsabs=2e-13,epsrel=2e-12)[0]
                   for a,b in zip(points,points[1:]))
    Z=integral(lambda u:1.)
    EU=integral(lambda u:u)/Z
    cov=[]
    for i in range(n):
        f=lambda u:Phi((u-1-i*w)/w)
        cov.append(integral(lambda u:u*f(u))/Z-EU*integral(f)/Z)
    # Canonical force from the original gradient, not the IBP identity.
    force=lambda u:.5*A*u+eta*A*w*sum(theta[i]*phi((u-1-i*w)/w) for i in range(n))
    mu_direct=integral(force)/Z if precision is None else None
    return -EU, np.array(cov), Z, mu_direct

grid=np.linspace(-.2,1.2,10001)
for t in grid:
    check(0<=phi(float(t))<=1/16+1e-14,'bump bounds')
    check(abs(dphi(float(t)))<=1,'bump derivative bound')
check(abs(quad(phi,0,1)[0]-1/30)<1e-14,'bump integral')
check(phi(0)==phi(1)==dphi(0)==dphi(1)==0,'C2 joins')

rng=np.random.default_rng(88023)
rows=[]
for A in [1.,.2,.02,.002]:
    for n in [2,3,5]:
        theta=rng.choice([-1.,1.],n)
        mu,cov,Z,mu_direct=posterior(A,theta)
        check(abs(mu-mu_direct)<2e-12,'exact posterior IBP')
        check(np.min(cov)>.0001,'positive normalized covariance')
        delta=.25*A/(30*n)
        Z0=math.sqrt(2*math.pi/(1+.5*A))
        check(math.exp(-delta)-1e-12<=Z/Z0<=math.exp(delta)+1e-12,'normalizer bound')
        for i in range(n):
            h=2e-4
            tp=theta.copy();tm=theta.copy();tp[i]+=h;tm[i]-=h
            fd=(posterior(A,tp)[0]-posterior(A,tm)[0])/(2*h)
            pred=.25*A/n**2*cov[i]
            check(abs(fd-pred)<=2e-9*max(1.,abs(pred)),'normalization-sensitive derivative')
        rows.append(dict(A=A,n=n,mean=mu,min_cov=float(min(cov))))

# Actual finite conditional sign enumeration: fixed transcript exposing n-m signs.
risks=[]
for A in [.1,.01]:
    n=6; exposed=[1.,-1.,1.]; unexposed=n-len(exposed)
    signs=np.array(list(itertools.product([-1.,1.],repeat=unexposed)))
    means=np.array([posterior(A,np.array(exposed+list(s)))[0] for s in signs])
    coefficients=means@signs/len(signs)
    variance=float(np.var(means))
    walsh=float(np.sum(coefficients**2))
    check(variance+1e-22>=walsh,'conditional Bessel lower bound')
    check(np.all(coefficients>0),'same-sign first Walsh coefficients')
    risks.append(dict(A=A,conditional_variance=variance,walsh_projection=walsh,
                      normalized_variance=variance*n**4/(A*A*unexposed)))

# Same-potential noncommuting extension, including integrated-u2 precision.
T=np.array([[.5,1/16],[1/16,.5]])
E=np.diag([1.,0.])
for t in grid[::10]:
    H=T+.125*dphi(float(t))*E
    eig=np.linalg.eigvalsh(H)
    check(min(eig)>=5/16-1e-14 and max(eig)<=11/16+1e-14,'2D Hessian sandwich')
H0=T+.125*dphi(.2)*E;H1=T+.125*dphi(.7)*E
comm=H0@H1-H1@H0
check(np.linalg.norm(comm)>1e-4,'same-potential noncommuting Hessians')
for A in [.001,.1,1.]:
    precision=1+A/2-(A/16)**2/(1+A/2)
    mu,cov,_,_=posterior(A,np.array([1.,-1.,1.]),eta=.125,precision=precision)
    check(min(cov)>.0001,'2D marginal covariance sensitivity')
    matrix=np.eye(2)+A*T
    check(abs(np.linalg.inv(matrix)[0,0]-1/precision)<1e-13,'Schur precision')

# Quadratic CV cancellation and cost-exponent bookkeeping.
for A in [1.,.1,.001]:
    for lam in [.25,.5,.75]:
        for y in [-3.,0.,2.]:
            mean=math.sqrt(A)*lam*y/(1+lam*A)
            z=rng.normal(size=50)
            x=y/(1+lam*A)+math.sqrt(A/(1+lam*A))*z
            reflected=2*y/(1+lam*A)-x
            anti=math.sqrt(A)*lam*(x+reflected)/2
            check(np.max(np.abs(anti-mean))<1e-14,'quadratic antithetic exact mean')
for j in range(2,100):
    P=2*j-1;c=2*(j-1)/3
    check(abs(1+1.5*c-j)<1e-12,'strong mean necessary exponent')
    check(abs(c-(P-1)/3)<1e-12,'outer P conversion')

note=ROOT/'STRONG-POSTERIOR-MEAN-ORACLE-BARRIER.md'
out={'status':'PASS','assertions':checks,'proof_sha256':hashlib.sha256(note.read_bytes()).hexdigest(),
     'posterior_fixtures':rows,'conditional_projection_fixtures':risks,
     'noncommuting_commutator_frobenius':float(np.linalg.norm(comm)),
     'scope':'Finite formula checks only; information lower bound is proved in the note.'}
(ROOT/'strong_mean_barrier_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
