#!/usr/bin/env python3
"""Author diagnostics of the finite stage-two VALUE graph; no mean compiler."""
import json,math
from pathlib import Path
import numpy as np
rng=np.random.default_rng(2026100502)
checks=0
def check(x):
 global checks
 assert bool(x);checks+=1

ts=np.array([.2,1.2]); rr=np.exp(-ts)
ps=np.array([(.5-rr[1])/(rr[0]-rr[1]),(rr[0]-.5)/(rr[0]-rr[1])])
tq=np.array([math.log(4/3),math.log(2)]);pq=np.array([2/3,1/3])
lin=[(a,h,pa*ph) for a,pa in zip(ts,ps) for h,ph in zip(ts,ps)]
quad=[(a,h,pa*ph) for a,pa in zip(tq,pq) for h,ph in zip(ts,ps)]

def single(s):
 r=math.exp(-s);a=2*s*r;b=2*s*(s-1)*r
 return np.array([r,a,b,math.sqrt(1-r*r-a*a-b*b),0.])

def pair(a,h):
 ra=math.exp(-a);rh=math.exp(-h)
 sa=math.sqrt(-math.expm1(-2*a));sh=math.sqrt(-math.expm1(-2*h))
 R=np.array([[2*a*ra/sa,2*a*(a-1)*ra/sa],
             [2*h*ra*rh/sh,2*h*(2*a+h-1)*ra*rh/sh]])
 K=np.eye(2)-R@R.T
 det=math.sqrt(np.linalg.det(K))
 L=(K+det*np.eye(2))/math.sqrt(np.trace(K)+2*det)
 u=np.array([ra,sa*R[0,0],sa*R[0,1],sa*L[0,0],sa*L[0,1]])
 v=rh*u+np.array([0.,sh*R[1,0],sh*R[1,1],sh*L[1,0],sh*L[1,1]])
 return u,v,K

class Source:
 def __init__(self,A,D,K=None,noise=0.):
  self.A=A;self.D=D;self.K=K;self.noise=noise;self.calls=0
  self.a=np.array([1.,0.]);self.b=np.array([.6,.8])
 def g(self,x):
  self.calls+=1
  if self.K is not None: out=self.K@x
  else:
   z=x.reshape(-1,2)
   out=(self.A*(z/2+(np.sin(z@self.a)[:,None]*self.a+np.sin(z@self.b)[:,None]*self.b)/8)).reshape(-1)
  if self.noise:
   e=rng.normal(size=self.D);out=out+self.noise*e/np.linalg.norm(e)
  return out
 def h(self,x):
  if self.K is not None:return self.K
  out=np.zeros((self.D,self.D))
  for j in range(0,self.D,2):
   z=x[j:j+2]
   out[j:j+2,j:j+2]=self.A*(np.eye(2)/2+(math.cos(z@self.a)*np.outer(self.a,self.a)+math.cos(z@self.b)*np.outer(self.b,self.b))/8)
  return out

def graph(src,roots,first=True):
 D=src.D;I=np.eye(D)
 def linear(coef):return coef@roots,np.hstack([c*I for c in coef])
 def plus(*terms):return sum(t[0] for t in terms),sum(t[1] for t in terms)
 def scale(a,t):return a*t[0],a*t[1]
 def force(t):
  f=src.g(t[0]);j=src.h(t[0])@t[1] if first else np.zeros_like(t[1]);return f,j
 x=linear(np.array([1.,0,0,0,0]));v=linear(np.array([.5,.5,0,0,0]));w=linear(np.array([.25,.5,.25,0,0]))
 gv=force(v);gw=force(w);ggw=force(gw)
 S=plus(x,scale(-1,gv),ggw);out=force(S)
 def central(d):return scale(.5,plus(force(plus(S,d)),scale(-1,force(plus(S,scale(-1,d))))))
 for s,p in zip(ts,ps):
  gu=force(linear(single(s)));d=plus(gv,scale(-1,gu));out=plus(out,scale(p,central(d)))
 for a,h,p in lin:
  ur,vr,_=pair(a,h);U=linear(ur);V=linear(vr)
  gu=force(U);gv2=force(V);shift=force(plus(U,scale(-1,gv2)))
  e=plus(gu,scale(-1,shift),scale(-1,ggw));out=plus(out,scale(p,central(e)))
 for a,h,p in quad:
  ur,vr,_=pair(a,h);gu=force(linear(ur));gv2=force(linear(vr))
  d1=plus(gv,scale(-1,gu));d2=plus(gv,scale(-1,gv2))
  for eps in [-1,1]:
   for zeta in [-1,1]:
    out=plus(out,scale(p*eps*zeta/8,force(plus(S,scale(eps,d1),scale(zeta,d2)))))
 out=plus(out,scale(-1,force(x)))
 return out

check(abs(sum(p*math.exp(-(a+h)) for a,h,p in lin)-.25)<1e-14)
check(abs(sum(p*(math.exp(-a)+math.exp(-(a+h))) for a,h,p in quad)-1)<1e-14)
for a,h in np.exp(rng.uniform(-14,5,size=(1000,2))):
 u,v,K=pair(a,h)
 check(np.linalg.eigvalsh(K)[0]>=1/16-1e-12)
 check(abs(u@u-1)<1e-11);check(abs(v@v-1)<1e-11)
 check(abs(u@v-math.exp(-h))<1e-11)
 check(np.max(abs(u[:3]-single(a)[:3]))<1e-10)
 check(np.max(abs(v[:3]-single(a+h)[:3]))<1e-10)

worst_fd=worst_q=worst_floor=0.
for D in [2,4]:
 for A in [.001,.03,.1,.5]:
  for rep in range(3):
   roots=rng.normal(size=(5,D));src=Source(A,D)
   val,J=graph(src,roots)
   check(src.calls==5+3*len(ts)+5*len(lin)+6*len(quad))
   eps=2e-6;num=[]
   for j in range(5*D):
    dy=np.zeros(5*D);dy[j]=eps;dy=dy.reshape(5,D)
    num.append((graph(src,roots+dy,False)[0]-graph(src,roots-dy,False)[0])/(2*eps))
   err=np.linalg.norm(J-np.column_stack(num),2);worst_fd=max(worst_fd,err);check(err<1e-8)
   private=J.copy();private[:,:D]*=math.sqrt(3)/2
   lift=np.zeros((5*D,5*D));lift[:D]=private
   check(math.sqrt(2)*np.linalg.norm(private,2)<=9*A+1e-12)
   check(math.sqrt(2)*np.linalg.norm(lift-lift.T,2)<=27*A*A+1e-12)
   check(np.array_equal(graph(src,np.zeros((5,D)),False)[0],np.zeros(D)))
   O,_=np.linalg.qr(rng.normal(size=(D,D)));K=O@np.diag(rng.uniform(0,A,D))@O.T
   qs=Source(A,D,K);meanroots=np.zeros((5,D));meanroots[0]=roots[0]
   got=graph(qs,meanroots,False)[0]+K@roots[0]
   wanted=(K-K@K/2+K@K@K/4)@roots[0]
   qerr=np.linalg.norm(got-wanted);worst_q=max(worst_q,qerr);check(qerr<1e-12)
   nu=1e-5;pert=graph(Source(A,D,noise=nu),roots,False)[0]
   ratio=np.linalg.norm(pert-val)/((4.5+14*A+5.5*A*A)*nu)
   worst_floor=max(worst_floor,ratio);check(ratio<=1+1e-9)

for A in np.linspace(.000001,.5,1000):
 eta=A/2+A*A/4;q=A/2+A*A/2;r0=A*A/4;beta=math.sqrt(3)/2
 Cx=2.25*(1+eta)+A*(3.25+1.25*A);Cn=2.25*q+A*(3.25+1.5*A)
 Cm=2.25*r0+A*(2.5+1.25*A);Cl=A*(2.5+A)
 first=math.sqrt(2)*A*math.sqrt(beta*beta*Cx*Cx+Cn*Cn+Cm*Cm+2*Cl*Cl)
 curl=math.sqrt(2)*(2*beta*(2.25*A*eta+A*A*(3.25+1.25*A))+A*math.sqrt(Cn*Cn+Cm*Cm+2*Cl*Cl))
 check(first<=9*A);check(curl<=27*A*A)

result={'status':'PASS','assertions':checks,'maximum_jacobian_discrepancy':worst_fd,
        'maximum_quadratic_mean_discrepancy':worst_q,'maximum_value_floor_usage':worst_floor,
        'values_per_test_residual':5+3*len(ts)+5*len(lin)+6*len(quad),
        'scope':'Literal graph, covariance, radius/curl, zero, quadratic mean and absolute VALUE floor diagnostics. Does not implement completed mean compilers or certify a high-accuracy pair quadrature.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
