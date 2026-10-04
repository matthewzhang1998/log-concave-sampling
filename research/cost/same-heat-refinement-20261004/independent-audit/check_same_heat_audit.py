# Publication copy: verifies the three public output pins; mathematical checks are unchanged. Original/public script hashes are in INVENTORY.json.
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
import numpy as np
rng=np.random.default_rng(202610041449)
checks=0

def check(b):
 global checks
 assert b
 checks+=1

# Independent noncommuting-Hessian test for finite DAG source derivatives.
trials=[]
for a in [.001,.01,.05,.2,.49]:
 for N in [2,5,12,31]:
  t=np.linspace(0,math.pi/2,N+1)
  W=np.zeros((N+1,N+1))
  for i in range(1,N+1):
   for j in range(i):W[i,j]=math.cos(t[i]-t[j+1])-math.cos(t[i]-t[j])
  check(np.min(W)>=-1e-15)
  check(np.max(np.abs(W.sum(axis=1)-(1-np.cos(t))))<1e-14)
  # q0, fresh momentum, and caller derivatives, each whole 2x2 blocks.
  E=np.array([[.4,.2],[-.1,.3]])
  J0=np.eye(2)+math.sqrt(a)*E
  base=np.stack([np.cos(t)[:,None,None]*np.eye(2),
                 np.sin(t)[:,None,None]*np.eye(2),
                 np.eye(2)[None,:,:]+np.cos(t)[:,None,None]*(J0-np.eye(2))])
  X=base.copy()
  for M in range(1,9):
   H=[]
   for _ in range(N+1):
    U,_=np.linalg.qr(rng.normal(size=(2,2)))
    H.append(U@np.diag(rng.uniform(0,1,size=2))@U.T)
   H=np.array(H)
   Xnew=base.copy()
   for k in range(3):
    HX=H@X[k]
    Xnew[k]-=a*np.einsum('ij,jkl->ikl',W,HX)
   X=Xnew
   old=np.linalg.norm(X[0,-1],2)
   fresh=np.linalg.norm(X[1,-1]-np.eye(2),2)
   caller=np.linalg.norm(X[2,-1]-np.eye(2),2)
   B0=max(np.linalg.norm(b,2) for b in base[2])
   check(old<=a/(1-a)+1e-12)
   check(fresh<=a/(1-a)+1e-12)
   check(caller<=a*B0/(1-a)+1e-12)
  trials.append({'a':a,'N':N,'M':8,'old_endpoint_norm':old,'old_bound':a/(1-a)})

# Exact separable-prefix implementation of the same weight matrix.
for N in [1,3,17,100]:
 t=np.linspace(0,math.pi/2,N+1); vals=rng.normal(size=(N,7))
 C=np.zeros(7);S=np.zeros(7)
 for i in range(1,N+1):
  j=i-1
  C+=(math.cos(t[j+1])-math.cos(t[j]))*vals[j]
  S+=(math.sin(t[j+1])-math.sin(t[j]))*vals[j]
  fast=math.cos(t[i])*C+math.sin(t[i])*S
  slow=sum((math.cos(t[i]-t[j+1])-math.cos(t[i]-t[j]))*vals[j] for j in range(i))
  check(np.max(np.abs(fast-slow))<1e-13)

# Rational blind-gap areas for irregular adaptive baseline point patterns.
for den in range(2,100):
 nums=sorted(set([0,den]+[(j*j+3*j)%den for j in range(den//3)]))
 gaps=[F(nums[i+1]-nums[i],den) for i in range(len(nums)-1)]
 area=sum((ell*ell/F(30) for ell in gaps),F(0))
 check(area>=1/F(30*len(gaps)))
 check(sum(gaps)==1)

# Quadratic Hamiltonian invariance and stated contraction, multiple eigenvalues.
quad=[]
for a in [.0001,.001,.01,.1,.4]:
 for lam in [.01,.1,.5,1.]:
  om=math.sqrt(1+a*lam);c=math.cos(math.pi/2*om);s=math.sin(math.pi/2*om)/om
  var=a/(1+a*lam);new=c*c*var+a*s*s
  check(abs(new-var)<1e-14)
  check(abs(c)<=a/(1-a))
  err=(c*c*var+a*(1+(s-1)**2))-var
  check(abs(err-2*a*(1-s))<1e-14)
  quad.append({'a':a,'lambda':lam,'variance_error':new-var,'decoupled_variance_error':err})

p=Path(__file__).parent.parent
m = {'outputs': {'SAME-HEAT-POSTERIOR-REFINEMENT-AND-FLOW-GATE.md': '8f4aa7bc32bfce0576bf87eb8ea6591608bb7a195151838d69fef4ead7f80c59', 'check_same_heat_refinement.py': '57b02b55b5dac2d2d6467740cc692df16dbff52795d66fbe60462aa30b1443ef', 'same_heat_refinement_checks.json': '74e31939d3e5840fe98cfba1917597c2feba93cec65dd6474eaf579af21d4b28'}}
hashes=[]
for name,h in m['outputs'].items():
 actual=hashlib.sha256((p/name).read_bytes()).hexdigest();check(actual==h)
 hashes.append({'file':name,'sha256':actual})
out={'status':'PASS','checks':checks,'scope':'Independent finite matrix derivative bounds, prefix identity, exact gap-area arithmetic, quadratic invariance, source output hashes. No exact-flow, Gaussian-Lp, randomized, or posterior-law lower bound asserted.', 'derivative_trials':trials,'quadratic':quad,'verified_hashes':hashes}
Path(__file__).with_name('same_heat_independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print({'status':'PASS','checks':checks,'hashes_verified':len(hashes)})
