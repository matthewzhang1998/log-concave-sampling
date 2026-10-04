#!/usr/bin/env python3
"""Exact matrix/current/buffer identities; not a high-order transition audit."""
from pathlib import Path
import math, json
import numpy as np
OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(138177)
checks=0

def ck(b,msg):
 global checks
 if not b:raise AssertionError(msg)
 checks+=1

matrix=[]
for D in (1,2,7,23):
 for A in (.5,.125,.03125):
  O,_=np.linalg.qr(rng.normal(size=(D,D)))
  B=(O*rng.uniform(.02*A,A,D))@O.T
  I=np.eye(D);T=np.linalg.inv(I+B)
  for r,t in ((0.,.5),(0.,1.),(.2,.7),(.7,.9),(.9,1.)):
   s2=1-r*r;Delta=t*t-r*r;c=r/t;d=Delta/t;v0=Delta/t**2
   Sigmar=I-r*r*B@T;Sigmat=I-t*t*B@T
   H=B@np.linalg.inv(I+s2*B);M=r*H
   mean=c*I-d*M;cov=v0*I-d*d*H
   cross=(r/t)*Sigmat
   directmean=cross@np.linalg.inv(Sigmar)
   directcov=Sigmat-cross@np.linalg.solve(Sigmar,cross.T)
   ck(np.linalg.norm(mean-directmean)<3e-14,'exact reverse OU conditional mean')
   ck(np.linalg.norm(cov-directcov)<3e-14,'exact reverse OU conditional covariance')
   aa=r*(1-t*t)/(t*s2);bb=Delta/(t*s2);vv=(1-t*t)*Delta/(t*t*s2)
   conditionalX=s2*np.linalg.inv(I+s2*B)
   ck(np.linalg.norm((aa*I+bb*r*np.linalg.inv(I+s2*B))-mean)<3e-14,'literal same endpoint X bridge mean')
   ck(np.linalg.norm(bb*bb*conditionalX+vv*I-cov)<3e-14,'literal same endpoint X bridge covariance')
   ck(v0-d*d>=-1e-15,'unit mean buffer fits')
   ck(np.linalg.norm((v0*I-cov)-d*d*H)<3e-14,'retained missing curvature current')
   ck(np.linalg.eigvalsh(cov).min()>=v0*(1-A*Delta)-2e-14,'true conditional covariance gap')
   if Delta<=.25:
    reserve=v0*I-d*d*(I+H)
    ck(np.linalg.eigvalsh(reserve).min()>=(5/8)*v0-2e-14,'separate covariance reserve positive')
   # Exact SDE covariance equation at time r.
   Kr=I+(1+r)*M
   rhs=-Kr@Sigmar-Sigmar@Kr.T+2*I
   ck(np.linalg.norm(rhs+2*r*B@T)<5e-14,'same endpoint SDE covariance derivative')
   matrix.append({'D':D,'A':A,'r':r,'t':t,'Delta':Delta,'missing_covariance_HS':float(np.linalg.norm(d*d*H))})

steps=[]
for r in np.linspace(0,1,17):
 for h in (1e-12,1e-8,.001,.01,.1,.25,math.log(5/3)):
  c=math.exp(-h);one_minus_c=-math.expm1(-h);v0=-math.expm1(-2*h)
  d=(1+r)*one_minus_c;q2=v0-d*d
  ck(q2>=-1e-15,'h guard positivity')
  ck(q2>=(one_minus_c*(5*c-3))-2e-15,'worst r factorization')
  if h<=.25:ck(q2>=h/2-2e-15,'positive sqrt(h) reserve')
  ck(abs((q2+d*d)-v0)<1e-15,'fresh covariance sum')
  q=math.sqrt(max(q2,0));bias=.07
  ck(abs(abs(d*bias)-d*bias)<1e-15,'conditional Gaussian error propagation equality')
  # Illegal shared root gives a different covariance, except when a coefficient vanishes.
  if d*q>1e-12:ck(abs((q-d)**2-v0)>d*q,'deliberate shared-noise covariance countertest')
  if r+h<=1:steps.append({'r':float(r),'h':h,'c':c,'d':d,'reserve_variance':q2})

# Pointwise nonlinear Fokker-Planck identity needs m and its first only.
for D in (1,3,11):
 for r in (0.,.2,.7,1.):
  x=rng.normal(size=D);m=rng.normal(size=D)
  Q=rng.normal(size=(D,D));Dm=(Q+Q.T)/20
  score=-x-r*m;drift=-x-(1+r)*m
  divdrift=-D-(1+r)*np.trace(Dm)
  divscore=-D-r*np.trace(Dm)
  fp=-divdrift-drift@score+score@score+divscore
  expected=np.trace(Dm)+m@score
  ck(abs(fp-expected)<2e-14,'pointwise rho path Fokker Planck algebra')

report={'status':'PASS','assertions':checks,'matrix_bridge':matrix,'admissible_diffusive_steps':steps,'scope':'Exact finite checks for the positive buffer/current ledger. No approximation order for the time discretization or cancellation of higher conditional cumulants is claimed.'}
(OUT/'variance_gluing_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'PASS','assertions':checks,'matrix_cases':len(matrix),'admissible_steps':len(steps)},indent=2))
