#!/usr/bin/env python3
"""Independent analytic-formula and actual-graph diagnostics for the separator."""
import numpy as np, math, json, hashlib
from pathlib import Path
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
P=Path(__file__).resolve().parent
S=P.parent/'NONLINEAR-SEPARATOR-FOR-THE-QUADRATIC-AMPLIFIER.md'
PIN='6fd7c048f1cf4980284d7ea196213aa1655b3cf8fabc7a18bfaf9b19963cee5d'
checks=0
rng=np.random.default_rng(289173)
def ck(x,s):
 global checks
 checks+=1
 if not bool(x):raise AssertionError(s)
def Q(K,m):
 t,w=leggauss(m);r=[];W=[]
 for k in range(K):
  a=2.**(-k-1);r.extend(1-a*(1.5+.5*t));W.extend(a*.5*w)
 r.append(1-2.**(-K-1));W.append(2.**(-K))
 return np.array(r),np.array(W)
def GH(n,d):
 x,w=hermgauss(n);x*=math.sqrt(2);w/=math.sqrt(math.pi)
 inds=np.stack(np.meshgrid(*([np.arange(n)]*d),indexing='ij'),axis=-1).reshape(-1,d)
 return x[inds],np.prod(w[inds],axis=1)
def hprime(t):
 s=np.sqrt(np.abs(t));return s/(1+s)
def h(t):
 u=np.abs(np.asarray(t));s=np.sqrt(u);out=u-2*s+2*np.log1p(s)
 mask=u<.01
 if np.any(mask):
  um=u[mask];out[mask]=sum((-1.)**(j+1)*um**(1+j/2)/(1+j/2) for j in range(1,15))
 return np.sign(t)*out
def psi(t):
 u=np.abs(np.asarray(t));s=np.sqrt(u);out=u*u/2-4*u*s/3+2*(u-1)*np.log1p(s)-u+2*s
 mask=u<.01
 if np.any(mask):
  um=u[mask];out[mask]=sum((-1.)**(j+1)*um**(2+j/2)/((1+j/2)*(2+j/2)) for j in range(1,15))
 return out
b=.5;e=.1;eta=1e-4;v=np.ones(2)/math.sqrt(2);tau=math.exp(-.5)
delta0=b*e*tau/2-e*e*(math.exp(-1)-math.exp(-2))
perturb=5*eta+3*eta*eta
ck(hashlib.sha256(S.read_bytes()).hexdigest()==PIN,'source pin')
ck(delta0>.0128,'uniform exact scalar lower')
ck(delta0-perturb>.0122,'rough exact lower')
linear_constant=1.2+.6+1.2+math.sqrt(3)*.3+.3*math.sqrt(8)
quadratic_constant=2+math.sqrt(3)/2
ck(linear_constant<5 and quadratic_constant<3,'perturbation arithmetic')

def f(x,eta=eta):
 out=b*x.copy();out[...,0]+=e*np.sin(x[...,0]);out+=eta*h(x@v)[...,None]*v
 return out
def U(x,eta=eta):
 return b*np.sum(x*x,axis=-1)/2+e*(1-np.cos(x[...,0]))+eta*psi(x@v)
def Hess(x,eta=eta):
 out=np.broadcast_to(np.diag([b,b]),x.shape[:-1]+(2,2)).copy();out[...,0,0]+=e*np.cos(x[...,0])
 out+=eta*hprime(x@v)[...,None,None]*np.outer(v,v)
 return out
B0=np.diag([.6,.5])
H1=Hess(np.array([[1.,0.]]))[0]
comm=np.linalg.norm(B0@H1-H1@B0)
ck(comm>0,'noncommutation')
for _ in range(100):
 xx=rng.normal(size=(1,2))*10
 eig=np.linalg.eigvalsh(Hess(xx))[0]
 ck(eig.min()>=.4-1e-14 and eig.max()<=.6001+1e-14,'Hessian sandwich sample')
ratios=[]
for t in [1e-2,1e-4,1e-6,1e-8]:
 ht=Hess(np.array([t*v]))[0]
 ratios.append(float(np.linalg.norm(ht-B0)/t))
ck(ratios[-1]>90*ratios[0],'non-Lipschitz Hessian restriction')

# Independent 2D target coefficient from the normalized tilt derivative.
X,WX=GH(72,2)
vals=f(X);Hs=Hess(X);Us=U(X)
ED=np.einsum('n,nij->ij',WX,Hs)
Ctarget=np.einsum('n,ni,nj->ij',WX,vals,vals)+np.einsum('n,n,nij->ij',WX,Us-(WX@Us),Hs)
# Alternative normalized second-moment expansion, requiring no Hessian values.
center=Us-(WX@Us)
Ctarget_alt=.5*(np.einsum('n,n,ni,nj->ij',WX,center**2,X,X)-np.eye(2)*(WX@(center**2)))
# Hölder Hessian quadrature converges slowly; its tiny eta contribution is reported.
altgap=float(np.linalg.norm(Ctarget-Ctarget_alt))
ck(altgap<5e-6,'C2 target coefficient integral crosscheck')

Z4,W4=GH(16,4);Z=Z4[:,:2];G=Z4[:,2:]
scalar_roots,Wscalar=GH(64,2);zz=scalar_roots[:,0];gg=scalar_roots[:,1]
rules=[('midpoint',np.array([.5]),np.array([1.])),('endpoints',np.array([0.,1.]),np.array([.5,.5])),('symmetric_pair',np.array([.15,.85]),np.array([.5,.5]))]
for K,m in [(3,2),(6,4),(10,7)]:
 rr,ww=Q(K,m);rules.append((f'dyadic_{K}_{m}',rr,ww))
qrows=[];program=[]
for label,r,w in rules:
 s=np.sqrt(1-r*r);beta=float(w@s);c=.25+beta*beta
 rho=np.outer(r,r)+np.outer(s,s)
 d=math.exp(-1)*float(w@np.sinh(rho)@w)
 delta=b*e*((1-c)+(c-.5)*tau)+e*e*(d+(1-c)*tau-(math.exp(-1)-math.exp(-2)))
 ck(abs(w.sum()-1)<1e-13 and abs(w@r-.5)<1e-13,'quadrature exact moments')
 ck(delta>=delta0-1e-14,'all-rules scalar coefficient')
 H=np.zeros_like(Z);Hscl=np.zeros_like(zz)
 for ri,si,wi in zip(r,s,w):
  H+=wi*f(ri*Z+si*G)
  q=ri*zz+si*gg;Hscl+=wi*(b*q+e*np.sin(q))
 source=np.einsum('n,ni,nj->ij',W4,H,H)+(1-c)/2*(B0@ED+ED@B0)
 deltamat=source-Ctarget
 ck(deltamat[0,0]>.0122,'genuine C2 full coefficient lower diagnostic')
 ck(abs(deltamat[0,0]-delta)<=perturb+5e-6,'genuine C2 perturbation bound diagnostic')
 # Explicit exact same-G scalar second moment, including all node pairs.
 raw=b*b*c+2*b*e*c*tau+e*e*d
 ck(abs(Wscalar@(Hscl*Hscl)-raw)<1e-12,'all same-root products')
 # Independent-node counterfactual generally differs. One node is the exception.
 badrho=np.outer(r,r)+np.diag(s*s)
 badd=math.exp(-1)*float(w@np.sinh(badrho)@w)
 if len(r)>1 and label!='endpoints':ck(abs(badd-d)>1e-5,'independent-root replacement detected')
 qrows.append({'rule':label,'nodes':len(r),'c':c,'d_Q':d,'scalar_Delta':delta,'C2_Delta11':float(deltamat[0,0]),'rough_perturbation':float(deltamat[0,0]-delta)})
 # Exact scalar finite VALUE graph and target integral, each fixed m=1,...,4.
 xx,ww=GH(96,1);xx=xx[:,0]
 for A in [.125,.03125,.0078125,.00390625]:
  wt=ww*np.exp(-A*(b*xx*xx/2+e*(1-np.cos(xx))));target=wt@(xx*xx)/wt.sum()
  Y=zz-A*Hscl;T=Y.copy();out=Y.copy();coef=1.;calls=0
  for m in range(1,5):
   def value(x):
    global calls
    calls+=1
    return A*(b*x+e*np.sin(x))
   a=value(T);bb=value(a);cc=value(bb);T=(c-1)*bb+c*cc
   coef*=(-.5-(m-1))/m;out+=coef*T
   scaled=(Wscalar@(out*out)-target)/(A*A)
   ck(calls==3*m,'literal graph gradient count')
   ck(abs(Wscalar@out)<1e-14,'literal odd graph mean')
   if A<=.0078125:ck(abs(scaled-delta)<.0009,'literal graph second coefficient limit')
   program.append({'rule':label,'A':A,'m':m,'scaled_variance_gap':float(scaled),'limit':delta})
 # Actual rough graph oddness/zero/finite level moment at one heat.
 A=.0625;Y=Z-A*H;T=Y.copy();out=Y.copy();coef=1.
 for m in range(1,5):
  a=A*f(T);bb=A*f(a);cc=A*f(bb);T=(c-1)*bb+c*cc
  coef*=(-.5-(m-1))/m;out+=coef*T
  ck(np.linalg.norm(W4@out)<1e-13,'genuine C2 graph mean')
  ck(np.linalg.norm(f(np.zeros((1,2))))==0,'genuine C2 graph anchor zero')

report={'source_sha256':PIN,'assertions':checks,'uniform_scalar_variance_coefficient_lower':delta0,'uniform_rough_variance_coefficient_lower':delta0-perturb,'uniform_standardized_W2_liminf_lower':(delta0-perturb)/2,'perturbation_linear_constant':linear_constant,'perturbation_quadratic_constant':quadratic_constant,'Hessian_commutator_norm':float(comm),'Hessian_difference_quotients':ratios,'target_coefficient_quadrature_crosscheck_error':altgap,'rules':qrows,'scalar_actual_graph':program,'scope':'Fixed standardized potential u_eta and U_A=A u_eta. Uniform counterexample to extrapolating the quadratic VALUE graph; not a counterexample to the positive order-two theorem, not an impossibility for different algorithms, and not a claim about one fixed unscaled physical V as A varies.'}
(P/'separator_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='scalar_actual_graph'},indent=2))
