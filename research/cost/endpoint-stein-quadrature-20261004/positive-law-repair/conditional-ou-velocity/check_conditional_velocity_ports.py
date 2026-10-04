#!/usr/bin/env python3
"""Literal raw source checks; does not replace the imported mean compiler proof."""
import numpy as np, math, json
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad
from pathlib import Path
OUT=Path(__file__).resolve().parent
checks=0

def ck(b,msg):
 global checks
 if not b:raise AssertionError(msg)
 checks+=1

def rule(alpha, exponent=1):
 delta=alpha**exponent
 K=math.ceil(math.log2(4/delta));m=math.ceil(math.log(16/delta,4))
 x,w=leggauss(m); rr=[];ww=[]
 for k in range(K):
  a=2.**(-k-1);rr.extend(1-1.5*a-.5*a*x);ww.extend(.5*a*w)
 rr.append(1-2.**(-K-1));ww.append(2.**(-K))
 return np.array(rr),np.array(ww)

def h(t):
 u=np.abs(t);q=np.sqrt(u);res=u-2*q+2*np.log1p(q)
 mask=u<1e-4
 if np.any(mask):
  a=u[mask];val=np.zeros_like(a)
  for k in range(1,20):val+=(-1)**(k+1)*a**(1+k/2)/(1+k/2)
  res=np.array(res);res[mask]=val
 return np.sign(t)*res

def hp(t):
 q=np.sqrt(np.abs(t));return q/(1+q)

rng=np.random.default_rng(500177)
rows=[]
for D in (2,5):
 axes=rng.normal(size=(3,D));axes/=np.linalg.norm(axes,axis=1)[:,None]
 base=np.diag(np.linspace(.25,.35,D));eps=.08
 def f(x):return base@x+eps*sum(h(np.atleast_1d(a@x))[0]*a for a in axes)
 def df(x):return base+eps*sum(hp(a@x)*np.outer(a,a) for a in axes)
 for A in (.25,.125):
  for s in (.95,.5,.15):
   r=math.sqrt(1-s*s);alpha=A*s*s;rr,ww=rule(alpha);ss=np.sqrt(1-rr*rr);n=len(rr);M=12
   z=rng.normal(size=D);u=rng.normal(size=D);v=rng.normal(size=D)
   g=lambda x:A*f(x);dg=lambda x:A*df(x)
   def source(z,u,v,derivs=False):
    x=r*z.copy();J=r*np.eye(D)
    for k in range(M):
     if derivs:J=r*np.eye(D)-s*s*dg(x)@J
     x=r*z-s*s*g(x)
    anchor=g(x);H=np.zeros(D);Hu=np.zeros((D,D));Hv=np.zeros((D,D));Hz=np.zeros((D,D))
    for ti,bi,wi in zip(rr,ss,ww):
     q=x+s*(ti*u+bi*v)
     value=anchor if np.array_equal(q,x) else g(q)
     H+=wi*s*(value-anchor)
     if derivs:
      Bj=dg(q);Hu+=wi*s*s*ti*Bj;Hv+=wi*s*s*bi*Bj;Hz+=wi*s*(Bj-dg(x))@J
    xp=x+s*(u-H);xb=x+s*u
    gp=anchor if np.array_equal(xp,x) else g(xp)
    gb=anchor if np.array_equal(xb,x) else g(xb)
    E=gp-gb
    if not derivs:return E
    B1=dg(xp);B0=dg(xb)
    Eu=s*((B1-B0)-B1@Hu);Ev=-s*B1@Hv
    Ez=(B1-B0)@J-s*B1@Hz
    return E,H,Hu,Hv,Eu,Ev,Ez,x,J
   E,H,Hu,Hv,Eu,Ev,Ez,x,J=source(z,u,v,True)
   DE=np.hstack([Eu,Ev]);square=np.vstack([DE,np.zeros_like(DE)]);curl=square-square.T
   expected=np.block([[s*(Hu@dg(x+s*(u-H))-dg(x+s*(u-H))@Hu),Ev],[-Ev.T,np.zeros((D,D))]])
   ck(np.linalg.norm(curl-expected)<1e-14,'exact noncommuting curl blocks')
   ck(np.linalg.norm(DE,2)<=A*s*(1+alpha)+1e-14,'complete private first')
   ck(np.linalg.norm(curl,2)<=4*A*A*s**3+1e-14,'curl grade')
   ck(np.linalg.norm(E)<=A*s*np.linalg.norm(H)+1e-14,'one-marked-energy pathwise')
   ck(np.linalg.norm(Ez,2)<=3*A*r+1e-14,'normalized caller first')
   ck(np.array_equal(source(z,np.zeros(D),np.zeros(D)),np.zeros(D)),'literal same-anchor zero')
   epsfd=1e-6
   du=rng.normal(size=D);dv=rng.normal(size=D);dz=rng.normal(size=D)
   fd=(source(z,u+epsfd*du,v+epsfd*dv)-source(z,u-epsfd*du,v-epsfd*dv))/(2*epsfd)
   fdz=(source(z+epsfd*dz,u,v)-source(z-epsfd*dz,u,v))/(2*epsfd)
   ck(np.linalg.norm(fd-Eu@du-Ev@dv)<2e-8,'private directional first')
   ck(np.linalg.norm(fdz-Ez@dz)<2e-8,'caller directional first')
   R=x-r*z+s*s*g(x)
   ck(np.linalg.norm(R)<=s*s*alpha**M*np.linalg.norm(g(r*z))+3e-14,'mode residual')
   rows.append({'D':D,'A':A,'r':r,'s':s,'nodes':n,'raw_E_private_values':n+2,'captured_values':M+1,'private_first':float(np.linalg.norm(DE,2)),'curl':float(np.linalg.norm(curl,2)),'energy':float(np.linalg.norm(E))})

# Exact dense quadratic conditional mean: terminal noise affects covariance, not mean.
quadratic=[]
for D in (2,7,19):
 for A in (.5,.125):
  O,_=np.linalg.qr(rng.normal(size=(D,D)));B=(O*np.linspace(.2*A,.9*A,D))@O.T
  for s in (.2,.8,1.):
   r=math.sqrt(1-s*s);z=rng.normal(size=D);mode=np.linalg.solve(np.eye(D)+s*s*B,r*z)
   rr,ww=rule(A*s*s);ss=np.sqrt(1-rr*rr);beta=ww@ss
   u=rng.normal(size=D);v=rng.normal(size=D);HH=s*s*B@(.5*u+beta*v)
   E=B@(mode+s*(u-HH))-B@(mode+s*u)
   ck(np.linalg.norm(E+s**3*B@B@(.5*u+beta*v))<2e-14,'quadratic raw chord identity')
   ck(np.linalg.norm(B@mode-r*B@np.linalg.solve(np.eye(D)+s*s*B,z))<2e-14,'exact quadratic force mean')
   quadratic.append({'D':D,'A':A,'s':s,'beta':float(beta),'mean_norm':float(np.linalg.norm(B@mode))})

# Nonlinear actual mean versus true tilted-bridge force, via independent 1D integration.
q,p=hermgauss(80);q*=math.sqrt(2);p/=math.sqrt(math.pi);u,v=np.meshgrid(q,q,indexing='ij');wp=np.outer(p,p)
scalar=[]
for A in (.25,.125,.0625,.03125):
 s=.8;r=.6;z=.7;a=r*z;alpha=A*s*s
 g=lambda x:A*(.5*x+.1*(1-np.cos(x)))
 U=lambda x:A*(.25*x*x+.1*(x-np.sin(x)))
 x=a
 for k in range(16):x=a-s*s*g(x)
 rr,ww=rule(alpha,2);ss=np.sqrt(1-rr*rr);H=np.zeros_like(u)
 for ri,si,wi in zip(rr,ss,ww):H+=wi*s*(g(x+s*(ri*u+si*v))-g(x))
 actual=float(np.sum(wp*g(x+s*(u-H))))
 den=quad(lambda t:math.exp(-.5*t*t-float(U(a+s*t))),-12,12,epsabs=1e-13)[0]
 num=quad(lambda t:float(g(a+s*t))*math.exp(-.5*t*t-float(U(a+s*t))),-12,12,epsabs=1e-13)[0]
 target=num/den;err=abs(actual-target)
 ck(err<=10*A**3*s**5,'conditional force third-order bound diagnostic')
 scalar.append({'A':A,'r':r,'s':s,'caller_z':z,'actual_mean_force':actual,'true_mean_force':target,'mean_error':err,'mean_error_over_A3s5':err/(A**3*s**5),'nodes':len(rr)})
report={'status':'PASS','assertions':checks,'raw_C2_ports':rows,'quadratic':quadratic,'nonlinear_force_mean':scalar,'scope':'Checks the literal native sources and their consumer input ports. It does not execute or reprove the pre-existing complete gradient/near-gradient mean compiler, and it does not establish same-endpoint composition.'}
(OUT/'conditional_velocity_ports_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'PASS','assertions':checks,'nonlinear_force_mean':scalar},indent=2))
