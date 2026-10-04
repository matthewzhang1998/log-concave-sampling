# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
rng=np.random.default_rng(810611)
t,w=leggauss(1200); xi=10*t; wt=10*w*np.exp(-xi*xi/2)/np.sqrt(2*np.pi)
checks=0; largest=0.; adjerr=0.
for a in [.001,.01,.1,.4]:
 for s in [.01,.1,.3,.8]:
  for k in [.2,2.,10.,40.]:
   for my in [0.,.7]:
    for mz in [0.,.3]:
     z=mz+s*xi; d=.5; eta=s*np.sqrt(1-d*d)
     H=np.exp(-s*s/2)*np.cos(my+a*np.sin(k*z)/k)
     H0=np.cos(k*z)
     H0ou=np.exp(-(k*eta)**2/2)*np.cos(k*(mz+d*(z-mz)))
     err=abs(np.dot(wt,H*(H0ou-H0)))
     largest=max(largest,err/a); assert err<=a*(1+1e-8)+1e-12; checks+=1
# Orthogonal Hermite mode version of the full Hilbert Dirichlet estimate.
for m in [1,2,5,17,100]:
 for K in [1,3,9,30]:
  for rep in range(10):
   A=rng.normal(size=(K,m)); B=rng.normal(size=(K,m)); deg=np.arange(1,K+1)
   A/=np.sqrt(np.sum(deg[:,None]*A*A)); B/=np.linalg.norm(B)
   d=rng.uniform(.05,.95); dm=d**deg
   lhs=np.sum(A*(dm[:,None]-1)*B)
   rhs=np.sum(((dm[:,None]-1)*A)*B)
   adjerr=max(adjerr,abs(lhs-rhs)); assert abs(lhs-rhs)<1e-13
   assert np.linalg.norm((dm[:,None]-1)*A)<=1+1e-12
   assert abs(lhs)<=1+1e-12; checks+=3
out={'checks':checks,'largest_scalar_transfer_to_a_ratio':largest,'self_adjoint_residual':adjerr,'scope':'Scalar nonlinear common-root high-frequency transfer and vector Hilbert/OU algebra; this does not certify host retention or an arbitrary nonlinear ancestor.'}
print(json.dumps(out,indent=2));open(str(Path(__file__).with_name('hidden_ancestor_transfer_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
