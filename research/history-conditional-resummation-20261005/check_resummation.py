import json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261005)
checks=0

def check(c, message):
    global checks
    if not bool(c): raise AssertionError(message)
    checks+=1

def g(x,A):
    return A*(0.5*x+0.25*np.logaddexp(x,-x)-0.25*math.log(2))
def gp(x,A): return A*(0.5+0.25*np.tanh(x))
def psi(x,xp,z,zp,N,M,e,delta,A):
    anchor=g(e*M,A)
    u=g(z+e*M,A)-anchor; up=g(zp+e*M,A)-anchor
    return (g(x+e*N-delta*u,A)-g(xp+e*N-delta*u,A)
           -g(x+e*N-delta*up,A)+g(xp+e*N-delta*up,A))/(2*delta)

# Complete finite VALUE graph: exchange symmetry and conditional-inner radius.
for A in [.5,.25,.1,.03]:
  for _ in range(100):
    x,xp,z,zp,N,M=rng.normal(size=6); e=.2; delta=.05
    p=psi(x,xp,z,zp,N,M,e,delta,A)
    check(abs(p+psi(x,xp,zp,z,N,M,e,delta,A))<1e-13,'inner exchange')
    check(abs(p+psi(xp,x,z,zp,N,M,e,delta,A))<1e-13,'outer exchange')
    u=g(z+e*M,A)-g(e*M,A); up=g(zp+e*M,A)-g(e*M,A)
    dz=-.5*(gp(x+e*N-delta*u,A)-gp(xp+e*N-delta*u,A))*gp(z+e*M,A)
    dzp=.5*(gp(x+e*N-delta*up,A)-gp(xp+e*N-delta*up,A))*gp(zp+e*M,A)
    check(math.hypot(dz,dzp)<=A*A/math.sqrt(2)+1e-14,'small inner derivative')

# Noncommuting Gaussian likelihood-ratio overlap determinant identity.
for n in [2,3,7]:
  for _ in range(30):
    U=rng.normal(size=(n,n)); K=(U+U.T)/2
    V=rng.normal(size=(n,n)); L=(V+V.T)/2
    K*=.12/max(1,np.linalg.norm(K,2)); L*=.11/max(1,np.linalg.norm(L,2))
    S=np.eye(n)+K; T=np.eye(n)+L
    direct=1/math.sqrt(np.linalg.det(S)*np.linalg.det(T)*np.linalg.det(np.linalg.inv(S)+np.linalg.inv(T)-np.eye(n)))
    overlap=1/math.sqrt(np.linalg.det(np.eye(n)-K@L))
    check(abs(direct-overlap)<2e-13,'noncommuting overlap')

# Positive finite pair clock masses, exact sum and bounded aggregate normalization.
level_records=[]
for T in [1.,2.,4.,8.,16.]:
    total=0.
    for k in range(1,13):
        N=2**k; d=T/N; q=math.exp(-2*d)
        # Direct sum is a diagnostic only, never the sampler implementation.
        direct=sum((s+1)*q**s for s in range(N-1))*math.exp(-2*d)
        formula=math.exp(-2*d)*(1-N*q**(N-1)+(N-1)*q**N)/(1-q)**2
        check(abs(direct-formula)<=1e-8*max(1,direct),'finite positive clock mass')
        ld=d*(-math.expm1(-d)); w=formula*ld*ld*math.tanh(d/2)/2
        check(w>0,'positive level weight'); total+=w
    check(total<2,'uniform aggregate clock weight')
    level_records.append({'T':T,'sum_w_12':total})

# Exact finite-distribution replica/covariance and scalar mixture identities.
# This tests the algebra on a finite Gaussian-quadrature surrogate, not a history path.
nodes=[-1.,1.]
records=[]
for A in [.4,.2,.1,.05]:
    rows=[]
    for x in nodes:
      for xp in nodes:
       for z in nodes:
        for zp in nodes:
         rows.append([psi(x,xp,z,zp,N,M,.3,.1,A) for N in nodes for M in nodes])
    rows=np.array(rows)
    C=float(np.mean(np.mean(rows,axis=1)**2))
    E2=float(np.mean(rows**2))
    check(E2>=C-1e-20,'positive nuisance variance')
    hs=np.array([p*q for row in rows for p in row for q in row])
    check(abs(np.mean(hs)-C)<1e-18,'cross-replica exact conditional Gram')
    v=.25; sigma=v+C
    cov=v+C+float(np.mean(hs**2))/(4*v)
    check(abs(cov-np.mean(v*(1+hs/(2*v))**2))<1e-16,'exact transport covariance')
    ks=(v*(1+hs/(2*v))**2-sigma)/sigma
    # Positive series avoids catastrophic subtraction of 1 at tiny A.
    chi2=0.
    for m in range(1,15):
        moment=float(np.mean(ks**m))
        chi2+=math.comb(2*m,m)/4**m * moment*moment
    check(chi2>=0 and np.max(np.abs(ks))<.05,'small positive Gaussian mixture')
    records.append({'A':A,'conditional_gram':C,'nuisance_excess_R1':E2-C,
                    'transport_covariance_excess':cov-sigma,
                    'chi2_series':chi2,'sqrt_chi2_over_A8':math.sqrt(chi2)/A**8})

# Exact balancing exponents.
from fractions import Fraction as F
p=F(10,9)
check(F(2)+2*p==F(38,9),'fine tail grade')
check(F(9,2)-p/4==F(38,9),'linear comparison grade')
check(F(5)-p/2==F(40,9),'quadratic comparison grade')

out={'status':'PASS','assertions':checks,'scope':'Finite source, positive-clock, noncommuting density and exact exponent diagnostics; no history path or native compiler execution.',
     'clock_normalizations':level_records,'replica_diagnostics':records}
(ROOT/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
