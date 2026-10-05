#!/usr/bin/env python3
"""Finite algebra/geometry diagnostics; does not execute the native compilers."""
import json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
import mpmath as mp

checks=0

def check(v, label):
    global checks
    assert bool(v),label
    checks+=1


def positive_rule(delta,M,n):
    z,b=roots_legendre(n)
    ss=[];ww=[]
    for k in range(M):
        lo=delta*k/M;hi=delta*(k+1)/M
        x=(hi+lo)/2+(hi-lo)/2*z
        w=(hi-lo)/2*b*np.exp(-x)
        mass=math.exp(-lo)*(-math.expm1(-(hi-lo)))
        w*=mass/w.sum()
        check(abs(w.sum()-mass)<=2e-14*mass,'panel mass')
        check((w>0).all(),'positive weights')
        ss.extend(x);ww.extend(w)
    return np.array(ss),np.array(ww)


def bridge_rows(delta,s):
    q=math.exp(-delta)
    J=len(s);r=np.zeros((J,J));last=0.0
    ca=1.0;cb=0.0;aa=[];bb=[]
    for i,t in enumerate(s):
        d=delta-last;v=t-last
        a=math.sinh(delta-t)/math.sinh(d)
        b=math.sinh(v)/math.sinh(d)
        sig2=(-math.expm1(-2*v))*(-math.expm1(-2*(delta-t)))/(-math.expm1(-2*d))
        check(sig2>0,'strict innovation')
        if i:r[i]=a*r[i-1]
        r[i,i]=math.sqrt(sig2)
        ca=a*ca;cb=a*cb+b;aa.append(ca);bb.append(cb)
        check(abs(ca-math.sinh(delta-t)/math.sinh(delta))<5e-13,'X ancestry')
        check(abs(cb-math.sinh(t)/math.sinh(delta))<5e-13,'Y ancestry')
        last=t
    a=np.array(aa);b=np.array(bb)
    K=r@r.T
    cov=K+np.outer(a,a)+np.outer(b,b)+q*(np.outer(a,b)+np.outer(b,a))
    exact=np.exp(-np.abs(s[:,None]-s[None,:]))
    check(np.max(np.abs(cov-exact))<1e-12,'joint OU covariance')
    check((r>=0).all(),'positive innovation rows')
    return a,b,r,K

rng=np.random.default_rng(703)
worst_spectral_ratio=0.0
for delta in [.002,.013,.08,.3,.65]:
  w=-math.expm1(-delta)
  for M in [1,2,3,7,16]:
    s,wt=positive_rule(delta,M,3)
    a,b,r,K=bridge_rows(delta,s)
    sig=np.linalg.norm(r,axis=1)
    check(np.max(sig**2)<=w*(1+1e-12),'max bridge variance')
    check(wt@sig<=w**1.5*(1+1e-12),'complete first')
    check(wt@a<=w*(1+1e-12),'X first')
    check(wt@b<=w*(1+1e-12),'Y first')
    for _ in range(8):
      D=4;Hs=[]
      for j in range(len(s)):
        U,_,_=np.linalg.svd(rng.normal(size=(D,D)))
        Hs.append(U@np.diag(rng.uniform(size=D))@U.T)
      Hs=np.array(Hs)
      block=np.hstack([np.einsum('j,j,jab->ab',wt,r[:,k],Hs) for k in range(len(s))])
      check(np.linalg.norm(block,ord=2)<=w**1.5*(1+1e-10),'matrix full-bank first')
      HX=np.einsum('j,j,jab->ab',wt,a,Hs)
      check(np.linalg.eigvalsh(HX).min()>-1e-12,'PSD caller')
      check(np.linalg.norm(HX,ord=2)<=w*(1+1e-10),'matrix caller bound')
    # Spectral signed quadrature kernel: high precision prevents cancellation.
    mp.mp.dps=50
    d=mp.mpf(delta);sm=list(map(mp.mpf,s));wm=list(map(mp.mpf,wt))
    # Make the total mass literal at high precision; float setup residual is a floor.
    wm[-1]+=1-mp.exp(-d)-sum(wm)
    for lam0 in [.1,.5,1,2,10,100,1000]:
      lam=mp.mpf(lam0)
      if lam==1:
        I=(1-(1+2*d)*mp.exp(-2*d))/2
        cross=[x*mp.exp(-x)+(mp.exp(-x)-mp.exp(x-2*d))/2 for x in sm]
      else:
        I=2/(1-lam)*((1-mp.exp(-(1+lam)*d))/(1+lam)-(1-mp.exp(-2*d))/2)
        cross=[(mp.exp(-x)-mp.exp(-lam*x))/(lam-1)+(mp.exp(-x)-mp.exp(lam*x-(1+lam)*d))/(lam+1) for x in sm]
      Q=sum(wm[i]*wm[j]*mp.exp(-lam*abs(sm[i]-sm[j])) for i in range(len(sm)) for j in range(len(sm)))
      err=I-2*sum(v*c for v,c in zip(wm,cross))+Q
      bound=2*lam*d*(d/M)**2
      check(err>=-mp.mpf('1e-35'),'spectral PSD')
      check(err<=bound*(1+mp.mpf('1e-10')),'spectral M^-2 energy bound')
      worst_spectral_ratio=max(worst_spectral_ratio,float(err/bound))

# Exact integrated exponent ledger at actual w,h,eta,mu and refined M.
a=F(4,5);h=F(9,5);eta=F(8,5);mu=F(1,2);m=F(1,5)
ledger={
 'smoothing':2+h,'smoothing_drift':1+2*h,'commutation':3+a,
 'bridge_mean':2+a+2,'bridge_covariance':3+3*a-h+m,
 'two_third_currents':4+F(9,2)*a-2*h,
 'near_raw':3+eta/2,'near_raw_w':3+eta-a/2,
 'native_mean_1':5-F(3,2)*a,'native_mean_2':5-a,
 'native_mean_3':4+mu-a/2,'native_mean_4':5-a/2,
 'native_mean_5':6-mu/2-F(3,2)*a,'native_mean_6':6-F(3,2)*a,
 'cubic_prior_bulk':6-2*a,'cubic_prior_endpoint':6-eta,
 'cubic_feedback_bulk':7-F(5,2)*a,'cubic_feedback_endpoint':7-F(3,2)*eta,
}
for k,v in ledger.items():check(v>=F(19,5),k)
guards={'mean':1-eta/2,'self_reserve':2-mu/2-eta,'mixedK':1-eta/3,
        'cubic_gap':3-F(3,2)*eta,'covariance_gap':2-eta}
for k,v in guards.items():check(v>0,k)
check(min(ledger.values())==F(19,5),'leading grade')
check(2-mu/2-eta/2==F(19,20),'pointwise residual not falsely O(A)')
check(2-mu/2-a/2==F(27,20),'aggregate residual closes')
result={'checks':checks,'status':'PASS','worst_spectral_ratio':worst_spectral_ratio,
        'ledger':{k:str(v) for k,v in ledger.items()},'guards':{k:str(v) for k,v in guards.items()},
        'scope':'Finite source scalar/algebra/PSD diagnostics, not full native numerical execution.'}
Path(__file__).with_name('covariance_prefix_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
