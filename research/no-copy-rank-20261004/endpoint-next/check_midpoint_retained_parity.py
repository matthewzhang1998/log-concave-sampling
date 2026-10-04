#!/usr/bin/env python3
"""Actual K3 midpoint identity and retained-field parity checks.
These are finite algebra diagnostics, not a law or all-rank theorem.
"""
import json, math
from pathlib import Path
import numpy as np

T=np.array([[.5,.05],[.05,.5]])
N=np.diag([1.,0.]); e1=np.array([1.,0.])
c0,s0=.6,.8; delta=.025
C=np.concatenate([c0*np.eye(2),s0*np.eye(2),np.zeros((2,2))],axis=1)
rng=np.random.default_rng(61004)
checks=0; max_error=0.; finite_clock=[]

def check_close(x,y,atol=2e-11):
 global checks,max_error
 err=float(np.max(np.abs(np.asarray(x)-np.asarray(y))))
 assert err<=atol,(err,x,y)
 max_error=max(max_error,err); checks+=1

def bump(z):
 u=10*(z+.5)
 if abs(u)>=1:return 0.,0.
 t=1/(1-u*u); b=math.exp(1-t)
 return (z+.5)*b,b*(1-2*u*u*t*t)

for a in [1/16,1/32,1/64,1/128]:
 r=a; eps=a**.9; q0=a*eps; kap=r*a
 def G(x): return T@x+delta*q0*bump(x[0]/q0)[0]*e1
 def H(x): return T+delta*bump(x[0]/q0)[1]*N
 def base(W):
  S,U,Z=W[:2],W[2:4],W[4:]
  x=c0*S+s0*U
  ds=G(S)-G(S+eps*Z)
  x1=x+a*ds
  p0,p1=G(x),G(x1)
  t0,t1=S+a*p0,S+a*p1
  A0,A1=H(x),H(x1); T0,T1=H(t0),H(t1)
  DD=np.concatenate([H(S)-H(S+eps*Z),np.zeros((2,2)),-eps*H(S+eps*Z)],axis=1)
  K0,K1=T0@A0,T1@A1
  JE=np.concatenate([r*(T0-T1)+r*a*c0*(K0-K1)-r*a*a*K1@DD[:,:2],r*a*s0*(K0-K1),r*a*a*eps*K1@H(S+eps*Z)],axis=1)
  DT=np.concatenate([np.eye(2),np.zeros((2,4))],axis=1)+a/2*((A0+A1)@C+a*A1@DD)
  DQ=(A0-A1)@C-a*A1@DD
  return r*(G(t0)-G(t1)),p0-p1,(t0+t1)/2,DT,DQ,JE,K0-K1,A0,A1,T0,T1
 def NH(W,V,h):
  E,q,t,DT,DQ,JE,*rest=base(W)
  tp,tm=t+h*V,t-h*V
  Hp,Hm=H(tp),H(tm)
  val=kap/(2*h)*(G(tp)-G(tm))
  JW=kap/(2*h)*(Hp-Hm)@DT
  JV=kap/2*(Hp+Hm)
  return val,JW,JV
 Wstar=np.r_[np.zeros(4),e1]
 for h in [.5,a**.2,a,a*a]:
  sigma=h/a
  for i in range(150):
   W=Wstar+q0*rng.normal(size=6) if i%2 else rng.normal(size=6)
   E,q,t,DT,DQ,JE,*_=base(W)
   mu=a*q/(2*h); Dmu=a*DQ/(2*h)
   nv,nw,nvjac=NH(W,mu,h)
   check_close(2*sigma*nv,E)
   check_close(2*sigma*(nw+nvjac@Dmu),JE)
   V=rng.normal(size=2)
   p=NH(W,V,h); m=NH(W,-V,h)
   check_close(p[0],-m[0]);check_close(p[1],-m[1]);check_close(p[2],m[2])
   check_close(p[0][1],kap*(T@V)[1])
   check_close(JE@((p[1]+m[1])/2).T,np.zeros((2,2)))
  # A strictly positive fine clock near the nonzero native target.
  v=q0*1e-4; cr=math.sqrt(1-v*v); RW=Wstar/cr
  jes=[]; bs=[]; jn=[]
  for i in range(300):
   X=rng.normal(size=6); Y=rng.normal(size=2)
   W=cr*RW+v*X
   E,q,t,DT,DQ,JE,Kdiff,*_=base(W)
   jes.append(JE);bs.append(r*a*Kdiff)
   jn.append((NH(W,v*Y,h)[1]+NH(W,-v*Y,h)[1])/2)
  JE=np.mean(jes,axis=0); B=np.mean(bs,axis=0); JN=np.mean(jn,axis=0)
  sym=lambda X:(X+X.T)/2
  target=s0*s0*B@B.T-c0*sym(JE[:,:2]@(B-B.T))
  mixed=sym(JE@JN.T)
  exact=(kap*delta)**2/400
  assert target[1,1]>.49*exact,(target,exact)
  checks+=1
  check_close(mixed,np.zeros((2,2)))
  finite_clock.append(dict(a=a,h=h,fine_width=v,opposite_orientation_22=float(target[1,1]),zero_width_22=exact,mixed_field_norm=float(np.linalg.norm(mixed))))
res=dict(status='PASS',checks=checks,max_absolute_identity_error=max_error,finite_clock_cases=finite_clock,scope='Actual finite VALUE/Jacobian algebra and retained parity; no all-rank or integrated-law conclusion.')
out=Path(__file__).with_name('midpoint_retained_parity_checks.json');out.write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({k:v for k,v in res.items() if k!='finite_clock_cases'},indent=2))
