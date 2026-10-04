#!/usr/bin/env python3
"""Exact-formula and finite-program diagnostics for the scoped separator."""
from pathlib import Path
import numpy as np, math, json
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
OUT=Path(__file__).resolve().parent
xx,ww=leggauss(8); rr=[]; wwq=[]
for k in range(12):
 a=2.**(-k-1); rr.extend(1-(1.5*a+.5*a*xx)); wwq.extend(.5*a*ww)
rr.append(1-2.**(-13)); wwq.append(2.**(-12))
r=np.array(rr);w=np.array(wwq);s=np.sqrt(1-r*r)
beta=w@s;c=.25+beta*beta;rho=np.outer(r,r)+np.outer(s,s)
d=math.exp(-1)*np.einsum('i,ij,j',w,np.sinh(rho),w)
b=.5;e=.1;tau=math.exp(-.5)
delta=b*e*((1-c)+(c-.5)*tau)+e*e*(d+(1-c)*tau-(math.exp(-1)-math.exp(-2)))
delta0=b*e*tau/2-e*e*(math.exp(-1)-math.exp(-2))
checks=0
def ck(test,msg):
 global checks
 if not test: raise AssertionError(msg)
 checks+=1
ck(delta>=delta0>.0128,'exact uniform scalar lower bound')
ck(delta0-(5e-4+3e-8)>.0122,'noncommuting perturbation lower')
q,p=hermgauss(90);q*=math.sqrt(2);p/=math.sqrt(math.pi)
z,gg=np.meshgrid(q,q,indexing='ij');wp=np.outer(p,p)
f=lambda x:b*x+e*np.sin(x)
u=lambda x:b*x*x/2+e*(1-np.cos(x))
H=np.zeros_like(z)
for ri,si,wi in zip(r,s,w): H+=wi*f(ri*z+si*gg)
ck(abs(np.sum(wp*H*H)-(c*b*b+2*c*b*e*tau+e*e*d))<1e-12,'shared root second moment formula')
ck(abs(np.sum(wp*z*H)-.5*(b+e*tau))<1e-12,'first-order covariance')
rows=[]
for A in (.25,.125,.0625,.03125,.015625,.0078125):
 wt=p*np.exp(-A*u(q));target=float(np.sum(wt*q*q)/sum(wt))
 Y=z-A*H;T=Y.copy();out=Y.copy();coeff=1.
 for m in range(1,5):
  U1=A*f(T);U2=A*f(U1);U3=A*f(U2)
  T=(c-1)*U2+c*U3;coeff*=(-.5-(m-1))/m;out+=coeff*T
  var=float(np.sum(wp*out*out));scaled=(var-target)/(A*A)
  rows.append({'A':A,'m':m,'source_variance':var,'target_variance':target,'variance_defect_over_A2':scaled})
  ck(abs(float(np.sum(wp*out)))<1e-14,'actual odd graph mean')
  if A<=.015625: ck(abs(scaled-delta)<.0015,'small-A exact formula limit')
# Pointwise C2 noncommuting Hessian, with the actual rough ridge.
v=np.ones(2)/math.sqrt(2); eta=1e-4
hp=lambda x:math.sqrt(abs(x))/(1+math.sqrt(abs(x)))
B0=np.diag([.6,.5]);B1=np.diag([.5+.1*math.cos(1),.5])+eta*hp(v[0])*np.outer(v,v)
comm=float(np.linalg.norm(B0@B1-B1@B0))
ck(comm>0,'noncommuting Hessians')
for x in (np.zeros(2),np.array([1.,0.]),np.array([-4.,2.]),np.array([.0001,0.])):
 B=np.diag([.5+.1*math.cos(x[0]),.5])+eta*hp(v@x)*np.outer(v,v)
 ck(np.linalg.eigvalsh(B).min()>=.4-1e-14,'global sandwich sample lower')
 ck(np.linalg.eigvalsh(B).max()<=.6001+1e-14,'global sandwich sample upper')
report={'status':'PASS','assertions':checks,'beta':float(beta),'c':float(c),'d_Q':float(d),'exact_delta_Q':float(delta),'uniform_delta0':float(delta0),'rough_noncommuting_lower_bound':float(delta0-(5e-4+3e-8)),'commutator_norm':comm,'literal_scalar_program':rows,'scope':'Finite diagnostics support the analytical separator. Neither an all-order impossibility for all algorithms nor a failed positive endpoint law theorem is asserted.'}
(OUT/'nonlinear_amplifier_separator_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='literal_scalar_program'},indent=2))
