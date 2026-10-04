from pathlib import Path
import numpy as np,json,math
from numpy.polynomial.hermite import hermgauss
rng=np.random.default_rng(606170);checks=0;worst=0.;ratios=[]
def ck(x):
 global checks
 checks+=1
 assert x,checks
def eq(a,b,tol=1e-12):
 global worst
 err=float(np.max(np.abs(np.asarray(a)-np.asarray(b))));worst=max(worst,err);ck(err<tol)
# Physical true/forward identity and positive kernel, all ordered before projection.
for _ in range(1000):
 j=.01*rng.normal(size=3);J=np.zeros((3,3));J[0]=j;P=np.array([[1.,0.,0.]])
 B=(P@J@J@P.T).item();cov=(P@J@J.T@P.T).item();O=j[1]**2+j[2]**2
 eq(cov-B,O)
 b0=.4;q=.7;vn=.5
 exact=b0*b0+(q/(2*b0))**2*O*O-q*O+(vn-b0*b0)
 eq(exact,vn-q*O+q*q*O*O/(4*b0*b0));ck(exact>0)
 # Different conditional fine roots give square of mean, not mean of square.
 probs=np.array([.2,.3,.5]);bs=.01*rng.normal(size=3);bar=probs@bs
 eq(np.einsum('i,j,i,j',probs,probs,bs,bs),bar*bar)
 ck(float(probs@(bs*bs))-bar*bar>=-1e-15)
# Exact scalar sine-source outside-in mean, using only the inner Gaussian integral.
x,w=hermgauss(100);x=x*np.sqrt(2);w=w/np.sqrt(np.pi)
pnodes,pweights=hermgauss(50);pnodes*=np.sqrt(2);pweights/=np.sqrt(np.pi)
for kap in [.12,.06,.03,.015]:
 for freq in [.4,1.,2.]:
  for phasea,phaseb in [(.2,.4),(.8,-.7)]:
   la=.3*kap;lb=.2*kap;ampa=.6*kap/freq;ampb=.7*kap/freq;beta=.5
   ea=la+ampa*freq*np.cos(phasea)*np.exp(-freq*freq/2)
   eb=lb+ampb*freq*np.cos(phaseb)*np.exp(-freq*freq/2)
   eb_energy=math.sqrt(lb*lb+ampb*ampb*((1-math.exp(-2*freq*freq)*math.cos(2*phaseb))/2-math.exp(-freq*freq)*math.sin(phaseb)**2)+2*lb*ampb*freq*math.cos(phaseb)*math.exp(-freq*freq/2))
   K=10;tn=(1/(4*K*K)+.25)/2+(.25-1/(4*K*K))/2*np.cos(np.arange(K)*np.pi/(K-1));bn=np.sqrt(tn)
   ds=np.array([np.prod([-tn[j]/(tn[i]-tn[j]) for j in range(K) if j!=i]) for i in range(K)])
   eq(ds.sum(),1.);ck(np.abs(ds).sum()<=4)
   errs=[]
   for p in pnodes:
    inner=lb*p+sum(ds[i]*ampb*np.cos(freq*np.sqrt(1-bn[i]**2)*x+phaseb)*np.sin(freq*bn[i]*p)/bn[i] for i in range(K))
    # This is E_z0,z2 of the literal four outer signed VALUES exactly.
    act=la*inner+ampa*np.cos(phasea)*np.exp(-freq*freq/2)*np.sin(freq*beta*inner)/beta
    errs.append(w@act-ea*eb*p)
   err=math.sqrt(pweights@np.square(errs));bound=kap**3*eb_energy+kap*eb_energy*4**(-K)
   ratio=err/bound;ratios.append(ratio);ck(ratio<2.)
# Actual canonical grades and additive calls.
from fractions import Fraction as F
for gi in range(1,21):
 g=F(gi,20)
 for zi in range(1,100):
  z=F(zi,100)
  if 1+g+z<F(7,3) and z<g:
   ck(2+g>1+g+z)
for K in range(2,15):
 for J in range(1,12):
  eq(2*J*(2*K+4)*6,2*J*(12*K+24))
out={'status':'PASS','checks':checks,'max_exact_residual':worst,'max_mixed_cubic_ratio':max(ratios),'scope':'Scalar physical Gram identity, positive reference, conditional-root product, exact mixed sine-source padding and additive work. Does not prove all source-class estimates.'}
Path(__file__).with_name('scalar_direct_gram_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
