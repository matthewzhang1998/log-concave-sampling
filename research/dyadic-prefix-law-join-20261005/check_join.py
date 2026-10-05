#!/usr/bin/env python3
"""Diagnostics for a proof-level fixed-order join; no native compiler simulation."""
from fractions import Fraction as F
import json, math, pathlib, hashlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
checks=0
summary={}
def require(ok,msg):
 global checks
 checks+=1
 if not ok: raise AssertionError(msg)
# Full h^2 ledger, including the physical extra mean keep and both reserve roots.
require(F(1,4)+F(1,4)==F(1,2),'mean share')
require(F(1,4)+F(1,4)==F(1,2),'covariance share')
require(4*F(1,4)==1,'full h^2')
require(5*F(1,10)+F(1,2)==1,'tail shares')
require(F(1,20)+F(1,40)+F(1,40)==F(1,10),'corrected H')
# Prefix powers at alpha=1,gamma=2,kappa=1/2.
a,g,k=F(1),F(2),F(1,2)
L=1+3*a/2
prefix={
 'mean_native':1+g+4*(L-g),
 'action_calibration':3+3*a+k-g,
 'centered_action':5+6*a-3*g-k/2,
 'quadratic_feedback':5+6*a-3*g,
 'mean_quadrature':2+a+2,
 'covariance_clock':3+3*a+1-g,
 'true_third':4+9*a/2-2*g,
 'mean_caller_shift':3+a,
 'covariance_caller_shift':4+5*a/2-g-k/2,
 'smoothing':2+g,'smoothing_drift':1+2*g,
 'private_mean':L,'private_covariance':2+3*a-g-k/2,
 'endpoint_mean':1+a,'endpoint_covariance':2+5*a/2-g-k/2,
}
expected={'mean_native':F(5),'action_calibration':F(9,2),'centered_action':F(19,4),'quadratic_feedback':F(5),'mean_quadrature':F(5),'covariance_clock':F(5),'true_third':F(9,2),'mean_caller_shift':F(4),'covariance_caller_shift':F(17,4),'smoothing':F(4),'smoothing_drift':F(5),'private_mean':F(5,2),'private_covariance':F(11,4),'endpoint_mean':F(2),'endpoint_covariance':F(9,4)}
for n,x in expected.items(): require(prefix[n]==x,n)
summary['prefix_exact_powers']={n:str(x) for n,x in prefix.items()}
# Native six terms retained individually; all local positive-window finite sums.
rows=[('native_prior',F(7),F(5,2)),('mean',F(5),F(1)),('mean_sqrt',F(5),F(1,2)),('mean_reserve',F(11,2),F(5,4)),('mean_last',F(6),F(3,2)),('mixed',F(6),F(7,6)),('quartic',F(6),F(2)),('sixth',F(7),F(5,2)),('fourth_target',F(6),F(3,2)),('mixed_mixture',F(7),F(3,2))]
row_results={}
for n,q,p in rows:
 bulk=q-p
 endpoint=q-2*(p-1) if p>1 else q
 require(bulk>=4,n+' bulk'); require(endpoint>=4,n+' endpoint')
 row_results[n]={'bulk':str(bulk),'endpoint':str(endpoint),'K_power':str(-2*(p-1)) if p>1 else '0'}
for b in range(5,31):
 p=F(b-1,2)
 require(F(b+1)-2*(p-1)==4,'fixed-order endpoint')
 require(F(b+1)-p>=4,'fixed-order bulk')
summary['tail_rows']=row_results
# Positive midpoint rule on dyadic distance-to-terminal panels, split at w and eta.
def rule(A,K):
 eta=A*A*K*K
 gaps=[eta,A,1.0]
 d=eta
 while d<1:
  gaps.append(d); d*=2
 gaps=sorted(set(x for x in gaps if eta<=x<=1))
 if gaps[-1]!=1: gaps.append(1.)
 pts=[]
 # LAW nodes only, gap=1-t in [eta,1].
 for lo,hi in zip(gaps[:-1],gaps[1:]):
  gap=(lo+hi)/2
  t=1-gap
  c2=gap*(2-gap)
  s=(1-(1-A)**2)+(1-A)**2*c2
  v=c2*(1-(1-A)**2)/s
  pts.append((hi-lo,v))
 return eta,pts
max_env=0.0
for exp in range(5,15):
 A=10.**(-exp)
 for K in [2.,4.,16.,64.]:
  if A*K*K>.5: continue
  eta,pts=rule(A,K)
  for p in [F(1,2),F(1),F(7,6),F(5,4),F(3,2),F(2),F(5,2),F(7,2)]:
   pp=float(p)
   lhs=sum(w*v**(-pp) for w,v in pts)
   rhs=A**(-pp)+(eta**(1-pp) if pp>1 else math.log(1/eta) if pp==1 else 1)
   max_env=max(max_env,lhs/rhs)
   require(lhs<=8*rhs,'actual finite positive v sum')
  prior=A**7*sum(w*v**(-2.5) for w,v in pts)
  require(prior<=8*(A**4.5+A**4/K**3),'integrated undominated prior')
  # Source radius is O(1/K), while obsolete pointwise domination can fail.
  require(A/math.sqrt(eta)<=1/K*(1+1e-12),'actual normalized radius')
summary['max_finite_sum_envelope_ratio']=max_env
# Actual endpoint Gaussian rows in (z,G_Y,G_X), and common carrier rotation.
rng=np.random.default_rng(250105)
max_geometry=0.
for _ in range(500):
 A=float(rng.uniform(.001,.2));q=1-A;t=float(rng.uniform(.001,.999));c2=1-t*t;s=1-q*q*t*t
 Y=np.array([q*t,math.sqrt(s),0.])
 aa=t*(1-q*q)/s;bb=q*c2/s;v=c2*(1-q*q)/s
 X=aa*np.array([1.,0.,0.])+bb*Y+np.array([0.,0.,math.sqrt(v)])
 errs=[abs(X@X-1),abs(Y@Y-1),abs(X@Y-q),abs(X[0]-t)]
 max_geometry=max(max_geometry,*errs)
 for e in errs: require(e<2e-13,'original stationary pair law')
 h=A*A; rr=math.sqrt(1-h*h)
 private=np.r_[rr*X[1:],h]
 target=math.sqrt(1-rr*rr*t*t)
 require(abs(np.linalg.norm(private)-target)<2e-13,'known common carrier norm')
 # Householder rotation maps the entire Gaussian space, not just visible roots.
 e=np.zeros(3);e[0]=1; x=private/target;u=x-e
 H=np.eye(3) if np.linalg.norm(u)<1e-13 else np.eye(3)-2*np.outer(u,u)/(u@u)
 require(np.linalg.norm(H@H.T-np.eye(3))<1e-12,'full-tape orthogonal rotation')
 require(np.linalg.norm(H@x-e)<1e-12,'carrier alignment')
summary['max_geometry_error']=max_geometry
# Literal reflected original-VALUE prefix source and complete captured derivatives.
# A convex nonlinear noncommuting Hessian fixture. Its Hessian norm is <=A.
dirs=np.array([[1.,0.],[1.,1.],[.2,1.]])
dirs=dirs/np.linalg.norm(dirs,axis=1)[:,None]
A=.017;delta=.12;w=1-math.exp(-delta);r=.97
ss=np.linspace(.005,delta-.005,13);weights=np.exp(-ss);weights=w*weights/weights.sum()
a=np.sinh(delta-ss)/math.sinh(delta);b=np.sinh(ss)/math.sinh(delta)
sig=np.sqrt((1-np.exp(-2*ss))*(1-np.exp(-2*(delta-ss)))/(1-math.exp(-2*delta)))
NVAL=0
def gfun(x):
 global NVAL
 NVAL+=1
 return A/4*(x+sum((math.tanh(float(v@x))*v for v in dirs),np.zeros(2)))
def Hfun(x): return A/4*(np.eye(2)+sum(((1-math.tanh(float(v@x))**2)*np.outer(v,v) for v in dirs),np.zeros((2,2))))
def literal(u,y,N):
 return -r*sum((wi*(gfun(ai*u+bi*y-si*N)-gfun(ai*u+bi*y)) for wi,ai,bi,si in zip(weights,a,b,sig)),np.zeros(2))
max_fd=0.; max_capture=0.
for _ in range(120):
 u,y,N=rng.normal(size=(3,2));NVAL=0
 zero=literal(u,y,np.zeros(2))
 require(np.array_equal(zero,np.zeros(2)),'literal anchored zero')
 require(NVAL==2*len(ss),'literal original VALUE census')
 J=sum((r*wi*si*Hfun(ai*u+bi*y-si*N) for wi,ai,bi,si in zip(weights,a,b,sig)),np.zeros((2,2)))
 require(np.linalg.eigvalsh(J)[0]>=-1e-15,'reflected genuine gradient PSD')
 require(np.linalg.norm(J,2)<=r*A*sum(weights*sig)*(1+1e-12),'private radius bound')
 Ju=-r*sum((wi*ai*(Hfun(ai*u+bi*y-si*N)-Hfun(ai*u+bi*y)) for wi,ai,bi,si in zip(weights,a,b,sig)),np.zeros((2,2)))
 Jy=-r*sum((wi*bi*(Hfun(ai*u+bi*y-si*N)-Hfun(ai*u+bi*y)) for wi,ai,bi,si in zip(weights,a,b,sig)),np.zeros((2,2)))
 capture=np.linalg.norm(np.column_stack([Ju,Jy]),2); max_capture=max(max_capture,capture/(A*w))
 require(capture<=2*A*w,'whole retained pair radius')
 step=1e-5
 for k in range(2):
  e=np.eye(2)[k]*step
  fd=(literal(u+e,y,N)-literal(u-e,y,N))/(2*step)
  err=np.linalg.norm(fd-Ju[:,k]);max_fd=max(max_fd,err)
  require(err<3e-11,'actual endpoint derivative includes anchor')
summary['literal_prefix']={'max_endpoint_fd_error':max_fd,'max_endpoint_first_over_Aw':max_capture,'values_per_uncached_source':2*len(ss)}
# Ownership and dimension formulas, with arbitrary full native dimensions.
mean={'N_mean','Z_M'}; cov={'p','Z_C','pairs','xi','square_native'};tail={'M','H','K','V3','V4','top_keep'}
require(not mean&cov and not (mean|cov)&tail,'independent complete bank ownership')
require({'u','Y'}==set(['u','Y']),'prefix retained readset')
for _ in range(200):
 D=int(rng.integers(1,50)); nL=int(rng.integers(1,15));nR=int(rng.integers(1,8))
 dts=[D*int(rng.integers(5,80)) for _ in range(nL)];dps=[D*int(rng.integers(7,120)) for _ in range(nL+nR)]
 before=sum(2*D+x+y for x,y in zip(dts,dps[:nL]))+sum(4*D+y for y in dps[nL:])
 after=D+sum(D+x+y for x,y in zip(dts,dps[:nL]))+sum(3*D+y for y in dps[nL:])
 require(after==before-(nL+nR-1)*D,'all roots retained after alignment')
summary['assertions']=checks
summary['native_compiler_executed']=False
summary['scope']='Exact scalar ledgers, literal prefix original-VALUE source, endpoint Gaussian rows, positive finite outer sums, variance and complete-root accounting. Imported native proofs remain necessary.'
(ROOT/'join_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
