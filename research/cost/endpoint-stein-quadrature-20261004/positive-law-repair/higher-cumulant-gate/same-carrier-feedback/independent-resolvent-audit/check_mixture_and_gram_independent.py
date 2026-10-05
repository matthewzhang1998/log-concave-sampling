#!/usr/bin/env python3
"""Companion independent finite tests; no author executable is imported.
The companion reuses the independent VALUE fixture, after rerunning its checks.
Imported finite mean compilers are NOT numerically instantiated by these tests.
"""
from pathlib import Path
import hashlib,itertools,json,math,runpy
import numpy as np
from scipy.special import roots_hermitenorm,ndtr,ndtri
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
prior=runpy.run_path(str(ROOT/'check_resolvent_kernel_independent.py'))
Source,Fixture,gauss_rule=prior['Source'],prior['Fixture'],prior['gauss_rule']
count=0;metrics={};rng=np.random.default_rng(202610042345)
def check(v,name):
 global count
 count+=1
 if not bool(v):raise AssertionError(name)
def close(a,b,tol=1e-10,name='equality'):
 check(np.max(np.abs(np.asarray(a)-np.asarray(b)))<=tol,name)
def op(a):return np.linalg.norm(a,2)
files={ROOT.parent/'GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md':'7476219bacabfacc44030a14710408735d3e31b78103376d57f525a3518726dc',ROOT.parent.parent/'order-reentry/RECTANGULAR-FIRST-COEFFICIENT-MEAN-AND-GRAM-RETURN.md':'13b0fad9f69355c511c8884ae1e7691066f740dccf387ec453344479193154ea'}
for p,digest in files.items():check(hashlib.sha256(p.read_bytes()).hexdigest()==digest,'frozen companion pin')

# Stein orientation: F=(X,X^2-1), resolvent derivative=(1,X), DF=(1,2X).
x,w=roots_hermitenorm(8);w/=math.sqrt(2*math.pi)
F=np.stack([x,x*x-1],axis=1);A=np.stack([np.ones_like(x),x],axis=1);DF=np.stack([np.ones_like(x),2*x],axis=1)
tau=np.einsum('ni,nj->nij',A,DF)
# Test scalar phi(f)=f_1*f_2. Reversing tau breaks the Stein identity.
dphi=F[:,::-1]
lhs=np.einsum('n,na,n->a',w,F,F[:,0]*F[:,1])
rhs=np.einsum('n,nab,nb->a',w,tau,dphi)
wrong=np.einsum('n,nba,nb->a',w,tau,dphi)
close(lhs,rhs,1e-12,'correct nonsymmetric Stein orientation')
check(np.linalg.norm(lhs-wrong)>.9,'old transposed orientation is detected')
metrics['Stein_orientation_lhs']=lhs.tolist();metrics['Stein_orientation_wrong_rhs']=wrong.tolist()

# Frobenius-orthonormal symmetric directions, density derivative and exact Wick.
wick=[]
for d in [1,2,4]:
  for rep in range(5):
    A=rng.normal(size=(d,d));S=A@A.T+.6*np.eye(d)
    E=rng.normal(size=(d,d));E=(E+E.T)/2;E/=np.linalg.norm(E,'fro');y=rng.normal(size=d)
    def density(T):return np.exp(-.5*y@np.linalg.solve(T,y)) / np.sqrt((2*math.pi)**d*np.linalg.det(T))
    inv=np.linalg.inv(S);eps=1e-5
    direct=(density(S+eps*E)-density(S-eps*E))/(2*eps)
    analytical=.5*density(S)*np.sum(E*(np.outer(inv@y,inv@y)-inv))
    close(direct,analytical,2e-10,'Gaussian density covariance derivative')
    lam,U=np.linalg.eigh(S);sqrtinv=(U*(1/np.sqrt(lam)))@U.T;T=sqrtinv@E@sqrtinv
    zz,ww=roots_hermitenorm(3);ww/=math.sqrt(2*math.pi)
    grid=np.array(list(itertools.product(range(3),repeat=d)));Z=zz[grid];weights=np.prod(ww[grid],axis=1)
    q=np.einsum('ni,ij,nj->n',Z,T,Z)-np.trace(T)
    norm2=weights@(q*q)
    close(norm2,2*np.sum(T*T),2e-12,'Wick quadratic isometry')
    check(norm2<=2/.6**2*np.sum(E*E)+1e-12,'nu-gapped Wick bound')
    wick.append(float(norm2))

# Scalar mixture quantiles: a finite numerical diagnostic of the A^4 mechanism.
x,w=roots_hermitenorm(100);w/=math.sqrt(2*math.pi)
pt,pw=leggauss(220);pt=(pt+1)/2;pw=pw/2
mix=[]
for nu in [.5,1.,2.]:
  for epsilon in [.03,.1,.3]:
    C=epsilon*(1+np.tanh(x));S=nu+C;target=nu+epsilon
    q=[]
    for p in pt:
      q.append(brentq(lambda y:float(w@ndtr(y/np.sqrt(S))-p),-20*np.sqrt(nu+2*epsilon),20*np.sqrt(nu+2*epsilon)))
    W2=np.sqrt(pw@((np.array(q)-np.sqrt(target)*ndtri(pt))**2))
    ec=epsilon*np.sqrt(w@np.tanh(x)**2)
    bound=epsilon*ec/(math.sqrt(2)*nu**1.5)
    check(W2<=bound+2e-7,'scalar covariance mixture bound')
    mix.append({'nu':nu,'epsilon':epsilon,'W2_quantile_diagnostic':float(W2),'theorem_bound':float(bound)})
metrics['mixture_quantile_diagnostics']=mix

# Choose strict-interior versions of the admitted Chebyshev-Lobatto squared nodes.
def first_rule(K):
  a=1/(4*K*K);b=.24
  nodes=(a+b)/2+(b-a)/2*np.cos(np.pi*np.arange(K)/(K-1))
  ds=np.array([np.prod([-nodes[j]/(nodes[i]-nodes[j]) for j in range(K) if j!=i]) for i in range(K)])
  return np.sqrt(nodes),ds
for K in [2,3,5,9,17,33,65]:
  b,d=first_rule(K)
  check(np.min(b)>=1/(2*K)-1e-13 and np.max(b)<.5,'first-filter widths')
  check(sum(abs(d))<=4+1e-10,'first-filter coefficient norm')
  close(sum(d),1,1e-12,'first-filter constant moment')
  for h in range(1,K):close(d@(b**(2*h)),0,1e-12,'first-filter killed moment')
  for h in [K,K+1,2*K,4*K]:check(abs(d@(b**(2*h)))<=4*4.**(-K)+1e-14,'first-filter Hermite tail')

# Partial-variable filter: H must remain fixed at each signed u evaluation.
# V=A tanh(u)+A^2(tanh(u+H)-tanh(H)) has exact V(0,H)=0 and small H first.
x,w=roots_hermitenorm(70);w/=math.sqrt(2*math.pi);A=.15
U,H=np.meshgrid(x,x,indexing='ij');W=w[:,None]*w[None,:]
V=lambda u,h:A*np.tanh(u)+A*A*(np.tanh(u+h)-np.tanh(h))
coef=float(np.sum(W*(A/np.cosh(U)**2+A*A/np.cosh(U+H)**2)))
filter_checks=[]
for K in [2,3,5]:
  b,d=first_rule(K);p=x[:,None,None];z=x[None,:,None];h=x[None,None,:]
  ef=np.zeros(len(x));wrong=np.zeros(len(x))
  for bi,di in zip(b,d):
    c=np.sqrt(1-bi*bi)
    response=(V(c*z+bi*p,h)-V(c*z-bi*p,h))/(2*bi)
    shrunken=(V(c*z+bi*p,c*h)-V(c*z-bi*p,c*h))/(2*bi)
    ef+=di*np.sum(response*W[None,:,:],axis=(1,2))
    wrong+=di*np.sum(shrunken*W[None,:,:],axis=(1,2))
  err=np.sqrt(w@((ef-coef*x)**2))
  # Absolute bound with generous universal constant; the exact multiplier proof is in audit.
  check(err<=4*4.**(-K)*(A+2*A*A)+1e-8,'partial-u filter calibration')
  filter_checks.append({'K':K,'partial_filter_error':float(err),'shrinking_H_changes_finite_target':float(np.linalg.norm(ef-wrong))})
check(filter_checks[0]['shrinking_H_changes_finite_target']>1e-6,'unchanged-H genealogy is material')
metrics['partial_filter_checks']=filter_checks

# Literal filter source counts, conditional origin, oddness, and actual private curl.
r,w=gauss_rule(3);t,v=gauss_rule(2)
for dim in [2,4]:
  f=Fixture(dim,.12);src=Source(f,r,w,t,v);b,d=first_rule(4)
  X,p,z,H=rng.normal(size=(4,dim));qv=2*len(r)*(len(t)+1)
  def FF(p,z,H):
    val=np.zeros(dim);du=np.zeros((dim,dim));dh=du.copy();queries=0
    for bi,di in zip(b,d):
      c=np.sqrt(1-bi*bi)
      vp=src.evaluate(X,c*z+bi*p,H);queries+=len(src.records)
      vm=src.evaluate(X,c*z-bi*p,H);queries+=len(src.records)
      jp,hp,_=src.derivatives(X,c*z+bi*p,H)
      jm,hm,_=src.derivatives(X,c*z-bi*p,H)
      val+=di*(vp-vm)/(2*bi);du+=di*c*(jp-jm)/(2*bi);dh+=di*(hp-hm)/(2*bi)
    return val,du,dh,queries
  ff,j,h,q=FF(p,z,H)
  check(q==2*len(b)*qv,'all signed filter source VALUE occurrences counted')
  close(FF(-p,z,H)[0],-ff,1e-14,'finite filter p oddness')
  close(FF(np.zeros(dim),z,H)[0],0,0,'finite filter zero p')
  origin=FF(p,np.zeros(dim),np.zeros(dim))[0]
  check(np.linalg.norm(origin)>1e-5,'retained-p conditional origin is nonzero')
  check(op(j-j.T)<=sum(abs(d)/b)*(.12**2/2)+1e-12,'partial filter z-z curl')
  check(op(h)<=sum(abs(d)/b)*2*.12**2*sum(src.alpha)+1e-12,'partial filter z-H curl')

# Exact positive variance and nonsymmetric Gram row, including source-zero carrier.
for D in [2,5,11]:
  for N in [1,3,9]:
    a=rng.uniform(.001,1,size=N);a/=sum(a);B=rng.normal(size=(N,D,D))*.05
    v0=.7;vm=v0/2;vf=v0/2
    rows=[]
    for ai,Bi in zip(a,B):rows += [np.sqrt(ai)*Bi,np.sqrt(vm*ai)*np.eye(D)]
    rows += [np.sqrt(vf)*np.eye(D)]
    row=np.concatenate(rows,axis=1);cov=row@row.T
    target=v0*np.eye(D)+sum(ai*Bi@Bi.T for ai,Bi in zip(a,B))
    close(cov,target,1e-12,'positive total variance and true Gram orientation')
    carrier=np.concatenate([np.sqrt(vm*ai)*np.eye(D) for ai in a]+[np.sqrt(vf)*np.eye(D)],axis=1)
    close(carrier@carrier.T,v0*np.eye(D),1e-12,'actual source-zero carrier variance')
    check(sum(np.sqrt(a))<=np.sqrt(N)+1e-12,'explicit public-log error summation factor')

out={'status':'PASS supplementary mixture and rectangular-Gram checks','new_assertions':count,'rerun_resolvent_assertions':prior['checks'],'seed':202610042345,'pins':{str(p.relative_to(ROOT.parent.parent)):h for p,h in files.items()},'scope':'Finite mean programs imported under their actual guards, not numerically instantiated. K-current, m3 and full endpoint remain open.','metrics':metrics}
(ROOT/'mixture_and_gram_independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'new_assertions':count,'rerun_resolvent_assertions':prior['checks']},indent=2))
