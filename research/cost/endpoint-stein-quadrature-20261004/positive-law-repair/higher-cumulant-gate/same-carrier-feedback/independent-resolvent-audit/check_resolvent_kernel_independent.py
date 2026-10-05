#!/usr/bin/env python3
"""Independent finite checks. They supplement, and do not prove, the analytical audit.
No author checker or executable source is imported. Original VALUE and derivative
fixtures are separate; derivatives here are diagnostics, never producer leaves.
"""
from pathlib import Path
import hashlib, itertools, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_hermitenorm
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'RESOLVENT-COVARIANCE-STABILITY-AND-FINITE-J-KERNEL.md'
PRIOR=ROOT.parent.parent/'THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md'
rng=np.random.default_rng(202610042337)
checks=0
metrics={}
def check(b, label):
    global checks
    checks+=1
    if not bool(b): raise AssertionError(label)
def close(a,b,tol=1e-10,label='equality'):
    check(np.max(np.abs(np.asarray(a)-np.asarray(b)))<=tol,label)
def op(a): return np.linalg.norm(a,2)
def gauss_rule(n):
    x,w=leggauss(n); return (x+1)/2,w/2
def dyadic_rule(j,n):
    x,w=leggauss(n); rr=[]; ww=[]
    for k in range(1,j+1):
        a=2.**(-k); s=1.5*a+0.5*a*x
        rr.extend(1-s); ww.extend(0.5*a*w)
    h=2.**(-j); rr.append(1-h/2); ww.append(h)
    return np.array(rr),np.array(ww)

check(hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='359cdfffa1ce4b7467d24ed24d6b8efdedb85a16c5de0b650d3b4ea37aa52cfd','source pin')
check(hashlib.sha256(PRIOR.read_bytes()).hexdigest()=='b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced','prior pin')
close(quad(lambda r:r*r/math.sqrt(1-r*r),0,1)[0],math.pi/4,2e-10,'row Hessian integral')
close(quad(lambda t:math.exp(-3*t)/math.sqrt(-math.expm1(-2*t)),0,np.inf)[0],math.pi/4,2e-9,'OU heat integral')
close((.5+math.pi/4)*.5+.5*math.pi/4,(1+math.pi)/4,label='paraproduct constant')
CSTAR=(1+math.pi)/4

# Exact finite-Hermite Gaussian divergence isometry, including cross derivatives.
divergence_ratios=[]
for d in [1,2,4,8]:
  inds=[a for a in itertools.product(range(4),repeat=d) if sum(a)<=3] if d<5 else [tuple(int(i==k) for i in range(d)) for k in range(d)]+[(0,)*d]
  ix={a:i for i,a in enumerate(inds)}
  for _ in range(8):
    T=rng.normal(size=(3,d,len(inds)))
    div={}
    for z,a in enumerate(inds):
      for k in range(d):
        b=list(a);b[k]+=1;b=tuple(b)
        div[b]=div.get(b,0)+np.sqrt(a[k]+1)*T[:,k,z]
    dt={}
    for z,a in enumerate(inds):
      for l in range(d):
        if a[l]:
          b=list(a);b[l]-=1;b=tuple(b)
          dt.setdefault(b,np.zeros((3,d,d)))[:,:,l]+=np.sqrt(a[l])*T[:,:,z]
    nt=np.sum(T*T); nd=sum(np.sum(v*v) for v in dt.values())
    nc=sum(np.sum(v*v.swapaxes(1,2)) for v in dt.values())
    nv=sum(np.sum(v*v) for v in div.values())
    close(nv,nt+nc,2e-9,'Gaussian divergence isometry')
    check(math.sqrt(nv)<=math.sqrt(nt)+math.sqrt(nd)+1e-11,'divergence Sobolev bound')
    divergence_ratios.append(math.sqrt(nv)/(math.sqrt(nt)+math.sqrt(nd)))
metrics['largest_divergence_ratio']=max(divergence_ratios)

# Pointwise row Hessians for non-gradient vector trig fields in various dimensions.
r,w=gauss_rule(180)
rowrat=[]
for d in [1,2,5,17,41]:
  k=rng.normal(size=(4,d)); k/=np.linalg.norm(k,axis=1)[:,None];k*=np.array([.4,1.,2.5,5.])[:,None]
  b=rng.normal(size=(4,d));b/=np.linalg.norm(b,axis=1)[:,None]
  F=rng.normal(size=(d,d));F/=max(op(F),1)
  L=op(F)+sum(np.linalg.norm(b[i])*np.linalg.norm(k[i]) for i in range(4))
  for _ in range(12):
    x=rng.normal(size=d);m=rng.normal(size=d)
    hess=np.zeros((d,d))
    for bi,ki in zip(b,k):
      c=np.dot(w,r*r*np.sin(r*(ki@x))*np.exp(-.5*(1-r*r)*(ki@ki)))
      hess-=float(m@bi)*c*np.outer(ki,ki)
    ratio=np.linalg.norm(hess,'fro')/(math.pi/4*L*np.linalg.norm(m))
    check(ratio<=1+1e-12,'dimension-free row-Hessian bound')
    rowrat.append(ratio)
metrics['largest_row_hessian_ratio']=max(rowrat)

# Hermite spectral resolution of R2[(Du)(D R1 f)] for f=sin(kx).
x,gw=roots_hermitenorm(220);gw=gw/math.sqrt(2*math.pi)
h=np.empty((101,len(x)));h[0]=1;h[1]=x
for n in range(1,100):h[n+1]=(x*h[n]-math.sqrt(n)*h[n-1])/math.sqrt(n+1)
pr=[]
for k in [.25,1.,3.,7.]:
  dv=(np.cos(x[:,None]*r[None,:]*k)*np.exp(-.5*(1-r*r)*k*k)[None,:]*(w*r*k)[None,:]).sum(axis=1)
  for n in [1,2,5,12,24]:
    product=math.sqrt(n)*h[n-1]*dv
    coeff=h@(gw*product)
    norm=np.linalg.norm(coeff/(np.arange(101)+2))
    ratio=norm/(CSTAR*k)
    check(ratio<=1+1e-10,'Hermite paraproduct bound')
    pr.append(ratio)
metrics['largest_paraproduct_ratio']=max(pr)

# Exact covariance polarization and orientation on nonsymmetric matrices.
orientation=[]
for d in [2,5,13]:
  for _ in range(12):
    A=rng.normal(size=(d,d));B=rng.normal(size=(d,d))
    close(2*(A@A.T-B@B.T),(A-B)@(A+B).T+(A+B)@(A-B).T,2e-12,'transpose polarization')
    orientation.append(np.linalg.norm(A@A.T-A@A,'fro'))
check(min(orientation)>1e-3,'ordinary square differs from oriented covariance')
metrics['minimum_orientation_separator']=min(orientation)

# Positive dyadic rules: moments, shifted multipliers, endpoint coefficient sums.
clock=[]
for panels in [2,5,10,18]:
  for order in [1,3,8,15]:
    rr,ww=dyadic_rule(panels,order);cc=np.sqrt(1-rr*rr)
    check(np.all((rr>0)&(rr<1)) and np.all(ww>0),'positive interior clocks')
    close(ww.sum(),1,label='clock mass');close(ww@rr,.5,label='clock first moment')
    check(np.sum(ww/cc)<=1+math.sqrt(2)+1e-12,'absolute dyadic coefficient sum')
    check(np.sum(ww*rr/cc)<=1+math.sqrt(2)+1e-12,'VALUE coefficient sum')
    check(np.sum(ww*rr*rr/cc)<=1+math.sqrt(2)+1e-12,'caller coefficient sum')
    deg=np.unique(np.r_[np.arange(300),np.geomspace(301,1e7,300).astype(int)])
    errors=np.abs((rr[None,:]**deg[:,None])@ww-1/(deg+1))
    clock.append({'panels':panels,'order':order,'nodes':len(rr),'max_tested_moment_error':float(max(errors)),'absolute_coefficient_sum':float(np.sum(ww/cc))})
metrics['clock_checks']=clock

class Fixture:
  def __init__(self,d,A):
    self.d=d;self.A=A
    vs=rng.normal(size=(d+2,d));vs/=np.linalg.norm(vs,axis=1)[:,None];vs*=.55
    base=rng.normal(size=(d,d));base=base@base.T;base=.18*base/op(base)
    envelope=base+vs.T@vs
    self.base=base/op(envelope);self.vs=vs/math.sqrt(op(envelope))
  def g(self,x):
    x=np.asarray(x);return self.A*(x@self.base.T+np.tanh(x@self.vs.T)@self.vs)
  def dg(self,x):
    z=self.vs@x
    return self.A*(self.base+self.vs.T@(self.vs/(np.cosh(z)**2)[:,None]))

class Source:
  def __init__(self,fixture,r,w,t,v):
    self.f=fixture;self.r=r;self.w=w;self.t=t;self.v=v;self.c=np.sqrt(1-r*r);self.s=np.sqrt(1-t*t)
    self.alpha=w*r/self.c;self.records=[];self.cache={}
  def val(self,site,key,noise=0):
    # Full semantic caller/root/version keys; exact +/- aliases at u=0 use one key.
    if key not in self.cache:
      self.records.append((key,site.copy()))
      perturb=noise*np.cos(site+np.arange(len(site)))
      perturb/=max(np.linalg.norm(perturb)/max(noise,1e-300),1)
      self.cache[key]=self.f.g(site)+perturb
    return self.cache[key]
  def evaluate(self,X,u,H,noise=0):
    self.records=[];self.cache={};total=np.zeros_like(u)
    context=(tuple(X),tuple(u),tuple(H),'frozen-v1',noise)
    for i,(ri,ci,ai) in enumerate(zip(self.r,self.c,self.alpha)):
      js=[]
      for branch,xx in [('plus',ri*X+ci*u),('anchor',ri*X)]:
        branch='alias' if np.array_equal(u,np.zeros_like(u)) else branch
        inner=np.zeros_like(u)
        for l,(ti,si,vi) in enumerate(zip(self.t,self.s,self.v)):
          inner+=vi*self.val(ti*xx+si*H,(context,i,branch,'inner',l),noise)
        js.append(self.val(xx-inner,(context,i,branch,'terminal'),noise))
      total+=ai*(js[0]-js[1])
    return total
  def J_derivatives(self,x,H):
    hx=np.zeros((len(x),len(x)));hh=hx.copy();inner=np.zeros_like(x)
    for t,s,v in zip(self.t,self.s,self.v):
      point=t*x+s*H;G=self.f.dg(point)
      inner+=v*self.f.g(point);hx+=v*t*G;hh+=v*s*G
    outer=self.f.dg(x-inner)
    return outer@(np.eye(len(x))-hx),-outer@hh
  def derivatives(self,X,u,H):
    d=len(X);du=np.zeros((d,d));dh=du.copy();dx=du.copy()
    for ri,ci,wi,ai in zip(self.r,self.c,self.w,self.alpha):
      jx,jh=self.J_derivatives(ri*X+ci*u,H);ax,ah=self.J_derivatives(ri*X,H)
      du+=wi*ri*jx;dh+=ai*(jh-ah);dx+=ai*ri*(jx-ax)
    return du,dh,dx

def finite_difference(fn,z):
  eps=2e-6;d=len(z);out=[]
  for i in range(d):
    e=np.zeros(d);e[i]=eps;out.append((fn(z+e)-fn(z-e))/(2*eps))
  return np.stack(out,axis=1)

r0,w0=gauss_rule(4);t0,v0=gauss_rule(3)
port=[]; anchor_change=[]; curl_nonzero=[]
for d in [2,3,7,15]:
  for A in [.03,.2,.5]:
    f=Fixture(d,A);src=Source(f,r0,w0,t0,v0)
    L=A*(1+A/2);beta=v0@np.sqrt(1-t0*t0);S1=sum(src.alpha);S2=sum(src.alpha*r0)
    for repetition in range(3):
      X,u,H=rng.normal(size=(3,d))
      V=src.evaluate(X,u,H)
      check(len(src.records)==2*len(r0)*(len(t0)+1),'literal safe original VALUE count')
      for i in range(len(r0)):
        for branch,xx in [('plus',r0[i]*X+src.c[i]*u),('anchor',r0[i]*X)]:
          rec=[(k,p) for k,p in src.records if k[1]==i and k[2]==branch]
          for l in range(len(t0)):
            close(rec[l][1],t0[l]*xx+src.s[l]*H,1e-13,'correct branch inner ancestor')
          expected=xx-sum(v0[l]*f.g(t0[l]*xx+src.s[l]*H) for l in range(len(t0)))
          close(rec[-1][1],expected,1e-13,'feedback terminal site')
      du,dh,dx=src.derivatives(X,u,H)
      close(du,finite_difference(lambda y:src.evaluate(X,y,H),u),3e-9,'actual u first')
      close(dh,finite_difference(lambda y:src.evaluate(X,u,y),H),3e-9,'actual H first')
      close(dx,finite_difference(lambda y:src.evaluate(y,u,H),X),3e-9,'actual X first')
      check(op(du)<=L/2+1e-12,'u radius')
      check(op(dh)<=2*A*A*beta*S1+1e-12,'private H first')
      check(op(dx)<=2*L*S2+1e-12,'actual caller first')
      check(np.linalg.norm(V)<=L*np.linalg.norm(u)/2+1e-12,'one-energy amplitude')
      check(op(du-du.T)<=A*A/2+1e-12,'u-u curl')
      lift=np.block([[du,dh],[np.zeros_like(du),np.zeros_like(du)]])
      check(op(lift-lift.T)<=A*A/2+2*A*A*beta*S1+1e-12,'complete lift curl')
      curl_nonzero.append(op(lift-lift.T))
      zero=src.evaluate(X,np.zeros(d),H,noise=1e-6)
      close(zero,np.zeros(d),0,'coherent numerical private origin')
      check(len(src.records)==len(r0)*(len(t0)+1),'same-key zero alias count')
      exact=src.evaluate(X,u,H);approx=src.evaluate(X,u,H,noise=1e-6)
      check(np.linalg.norm(approx-exact)<=2*S1*(1+A)*1e-6+1e-12,'absolute VALUE precision propagation')
      # Anchor is private-H dependent and must be reevaluated when H changes.
      anchor=lambda HH:f.g(r0[0]*X-sum(v0[l]*f.g(t0[l]*r0[0]*X+src.s[l]*HH) for l in range(len(t0))))
      anchor_change.append(np.linalg.norm(anchor(H)-anchor(H+.2)))
      port.append({'D':d,'A':A,'u_first':float(op(du)),'H_first':float(op(dh)),'X_first':float(op(dx)),'curl':float(op(lift-lift.T))})
check(max(curl_nonzero)>1e-5,'fixture exposes genuinely nonzero full curl')
check(min(anchor_change)>1e-8,'private H anchors genuinely change')
metrics['literal_port_checks']=port
metrics['minimum_private_anchor_change']=min(anchor_change)

# Exact quadratic calibration; the H-dependent part cancels in the anchored source.
for d in [1,3,9]:
  for A in [.05,.5]:
    B=rng.normal(size=(d,d));B=B@B.T;B/=op(B);B*=A
    class Quadratic:
      def g(self,x):return x@B.T
      def dg(self,x):return B
    src=Source(Quadratic(),r0,w0,t0,v0)
    X,u,H=rng.normal(size=(3,d));target=.5*B@(np.eye(d)-.5*B)
    close(src.evaluate(X,u,H),target@u,3e-13,'quadratic source value')
    du,dh,dx=src.derivatives(X,u,H)
    close(du,target,3e-13,'quadratic mean Jacobian');close(dh,0,3e-13,'quadratic H cancellation');close(dx,0,3e-13,'quadratic caller cancellation')

# Independent Gaussian integration by parts identifies E[D_u V] without computing a Jacobian of V.
f=Fixture(2,.2);src=Source(f,r0,w0,t0,v0);X=np.array([.4,-.7])
def gaussian_coefficient(order):
  z,weights=roots_hermitenorm(order);weights/=math.sqrt(2*math.pi)
  grid=np.array(list(itertools.product(range(order),repeat=4)))
  points=z[grid];weights=np.prod(weights[grid],axis=1);u=points[:,:2];H=points[:,2:]
  vals=np.zeros_like(u);B=np.zeros((2,2))
  for ri,ci,wi,ai in zip(src.r,src.c,src.w,src.alpha):
    xx=ri*X+ci*u;aa=np.broadcast_to(ri*X,xx.shape)
    inner=sum(v*f.g(t*xx+s*H) for t,s,v in zip(src.t,src.s,src.v))
    anchor=sum(v*f.g(t*aa+s*H) for t,s,v in zip(src.t,src.s,src.v))
    vals+=ai*(f.g(xx-inner)-f.g(aa-anchor))
    # D_x J obtained directly from original-g firsts and correct ancestors.
    hx=np.zeros((len(xx),2,2))
    for t,s,v in zip(src.t,src.s,src.v):
      zz=(t*xx+s*H)@f.vs.T
      dg=f.A*(f.base[None,:,:]+np.einsum('na,ai,aj->nij',1/np.cosh(zz)**2,f.vs,f.vs))
      hx+=v*t*dg
    zz=(xx-inner)@f.vs.T
    dg=f.A*(f.base[None,:,:]+np.einsum('na,ai,aj->nij',1/np.cosh(zz)**2,f.vs,f.vs))
    jx=dg@(np.eye(2)[None,:,:]-hx)
    B+=wi*ri*np.einsum('n,nij->ij',weights,jx)
  stein=np.einsum('n,ni,nj->ij',weights,vals,u)
  return B,stein
B12,S12=gaussian_coefficient(12);B18,S18=gaussian_coefficient(18)
close(B18,S18,2e-7,'E Du V = E V u* = B_Q')
close(B12,B18,3e-7,'independent coefficient quadrature convergence')
metrics['mean_jacobian_Stein_error']=float(np.max(np.abs(B18-S18)))
metrics['mean_jacobian_quadrature_change']=float(np.max(np.abs(B12-B18)))
metrics['mean_jacobian_skew_norm']=float(op(B18-B18.T))

out={'status':'PASS supplementary finite checks','assertions':checks,'seed':202610042337,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'prior_sha256':hashlib.sha256(PRIOR.read_bytes()).hexdigest(),'scope':'Analytical theorem and finite coefficient source only. No covariance action, K-current, m3 or endpoint completion is supplied.','metrics':metrics}
(ROOT/'resolvent_kernel_independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'assertions':checks,'mean_jacobian_Stein_error':metrics['mean_jacobian_Stein_error'],'largest_paraproduct_ratio':max(pr),'largest_row_hessian_ratio':max(rowrat)},indent=2))
