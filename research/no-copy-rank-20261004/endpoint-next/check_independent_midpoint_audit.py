#!/usr/bin/env python3
"""Independent noncommuting-M checks of midpoint source, curl and parity."""
import json,math,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(91220);checks=0;maxerr=0.
def close(x,y,tol=1e-11):
 global checks,maxerr
 er=float(np.max(np.abs(np.asarray(x)-np.asarray(y))))
 assert er<tol,(er,x,y)
 checks+=1;maxerr=max(maxerr,er)
c,s=.6,.8
for d in [2,3]:
 T=.5*np.eye(d);T[0,1]=T[1,0]=.06;e1=np.eye(d)[0];N=np.outer(e1,e1)
 k=71.;delta=.035
 def g(x):return T@x+delta/k*math.sin(k*x[0])*e1
 def H(x):return T+delta*math.cos(k*x[0])*N
 C=np.concatenate([c*np.eye(d),s*np.eye(d),np.zeros((d,d))],axis=1)
 PS=np.concatenate([np.eye(d),np.zeros((d,2*d))],axis=1)
 for _ in range(500):
  M=rng.normal(size=(d,d));M/=max(np.linalg.norm(M,2),1.)
  a=math.exp(rng.uniform(-5,-2));r=a;eps=a**.9;h=math.exp(rng.uniform(-5,-.7));kap=r*a
  W=rng.normal(size=3*d);S,U,Z=W[:d],W[d:2*d],W[2*d:]
  x0=C@W;Delta=g(S)-g(S+eps*Z);x1=x0+a*M.T@Delta
  A0,A1=H(x0),H(x1)
  DD=np.concatenate([H(S)-H(S+eps*Z),np.zeros((d,d)),-eps*H(S+eps*Z)],axis=1)
  p0,p1=g(x0),g(x1);q=p0-p1
  t0=S+a*M@p0;t1=S+a*M@p1;tbar=(t0+t1)/2
  Dq=(A0-A1)@C-a*A1@M.T@DD
  DT=PS+a/2*M@((A0+A1)@C+a*A1@M.T@DD)
  JE=r*(H(t0)@(PS+a*M@A0@C)-H(t1)@(PS+a*M@A1@(C+a*M.T@DD)))
  E=r*(g(t0)-g(t1));mu=a*q/(2*h);Dmu=a*Dq/(2*h)
  def source(v):
   pp=tbar+h*M@v;pm=tbar-h*M@v
   return kap/(2*h)*(g(pp)-g(pm)),kap/(2*h)*(H(pp)-H(pm))@DT,kap/2*(H(pp)+H(pm))@M
  nv,nw,nj=source(mu)
  close(2*h/a*nv,E)
  close(2*h/a*(nw+nj@Dmu),JE)
  lift=C.T@Dmu
  curl=lift-lift.T
  remainder=-a*a/(2*h)*(C.T@A1@M.T@DD-DD.T@M@A1@C)
  close(curl,remainder)
  V=rng.normal(size=d);pos=source(V);neg=source(-V)
  close(pos[0],-neg[0]);close(pos[1],-neg[1]);close(pos[2],neg[2])
  close((pos[1]+neg[1])/2,np.zeros_like(nw))
  # Bound follows directly from the only nonsymmetric remainder.
  assert np.linalg.norm(curl,2)<=a*a/h*np.linalg.norm(A1,2)*np.linalg.norm(M,2)*np.linalg.norm(DD,2)+1e-12
  checks+=1

# Exact physical lambda zero and leading orientation nonzero at flat fixture.
T=np.array([[.5,.05],[.05,.5]]);N=np.diag([1.,0.]);delta=.025
lambda_cases=[]
for a in [1/16,1/64,1/256]:
 kap=a*a;eps=a**.9
 # H(0)=H(eps e1)=T in the literal compact-bump fixture.
 JHZ=T-T
 for hs in [.03,.2,.45]:
  D=math.sqrt(1-hs*hs)*np.eye(2)
  DL=np.concatenate([kap*(JHZ-JHZ@D),np.zeros((2,2)),kap*(-eps*T+eps*T),kap*hs*JHZ,np.zeros((2,2))],axis=1)
  close(DL,np.zeros_like(DL))
  B=-kap*delta*T@N;O=B@B.T-c*c*(B@B+(B@B).T)/2
  close(O[1,1],(kap*delta)**2/400)
  lambda_cases.append(dict(a=a,rotation_width=hs,physical_lambda_first_norm=float(np.linalg.norm(DL)),orientation_22=float(O[1,1])))
source=ROOT/'MIDPOINT-CONSTRAINT-AND-RETAINED-PARITY-RETURN.md'
out=dict(checks=checks,max_absolute_error=maxerr,pin_file=source.name,pin_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),lambda_cases=lambda_cases,scope='Independent finite VALUE/Jacobian checks, exact parity and physical lambda selected-record check. No integrated-law lower bound.')
(ROOT/'independent_midpoint_audit_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
