#!/usr/bin/env python3
"""Independent supplementary diagnostics. They do not prove the open m3 bias gate."""
from pathlib import Path
import hashlib, json, itertools, math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_hermitenorm
import sympy as sy

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
PINS={
 'P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md':'87b61a6a4ca03f81f7d910cdaa83b3783cd4792b3a428b89c8f0f629062e4467',
 'SINGLE-HISTORY-EXACT-CENTERED-RESUMMED-CURRENT.md':'276d5873cc7d3f5ca8650955c08abe128cec7eda96680319f5c2a00b45118e9a',
 'COHERENT-SHIFT-SEPARATES-UNSHIFTED-M3-FEEDBACK.md':'cb04bb4e1cd3886ed378bb45593ea855b800e1cd1e532916c5f238cd3141e285',
 'GRADIENT-PORTS-DO-NOT-CONTROL-A-CONTRACTED-TRACE.md':'e6050915d148b7248e97bd26647f5948f14ef50d47dbd365d9c668b18790faa1',
 'independent-mean-separator-audit/INDEPENDENT-ANALYTICAL-AUDIT.md':'072f1672c7591305cb6d72fb78746ab5ee11e99733c011a40a4e0cbf0de73024',
}
checks=0; metrics={}
def check(cond,label):
 global checks
 checks+=1
 if not bool(cond): raise AssertionError(label)
def close(x,y,tol,label):
 check(np.max(np.abs(np.asarray(x)-np.asarray(y)))<=tol,label)
def op(x):return np.linalg.norm(x,ord=2)
for name,pin in PINS.items():check(hashlib.sha256((BASE/name).read_bytes()).hexdigest()==pin,'pin '+name)

rng=np.random.default_rng(202610050028)
r=np.array([.2,.8]); w=np.array([.5,.5])
t=np.array([.1,.5,.9]); v=np.array([.25,.5,.25])
sigma=np.array([.25,.75]); u=np.array([.5,.5])
c=np.sqrt(1-r*r); dmid=np.sqrt(1-t*t); ein=np.sqrt(1-sigma*sigma)
bo=w@c; bm=v@dmid; bi=u@ein
for nodes,weights in [(r,w),(t,v),(sigma,u)]:
 close(weights.sum(),1,1e-15,'clock mass');close(nodes@weights,.5,1e-15,'clock first moment')
 check(np.all(weights>0) and np.all((nodes>0)&(nodes<1)),'positive interior clock')

class Fixture:
 def __init__(self,D,A):
  self.D=D;self.A=A
  q=rng.normal(size=(D+2,D));q/=np.linalg.norm(q,axis=1)[:,None]
  self.q=q;self.eta=.6/len(q)
 def g(self,x):return self.A*(.25*x+self.eta*np.tanh(x@self.q.T)@self.q)
 def dg(self,x):return self.A*(.25*np.eye(self.D)+self.eta*np.einsum('a,ai,aj->ij',1/np.cosh(x@self.q.T)**2,self.q,self.q))

class Source:
 def __init__(self,f):self.f=f
 def eval(self,Z,G,H,J,noise=0,derivs=False):
  self.records=[];self.derivative_sites=[]
  D=len(Z);Id=np.eye(D);out=np.zeros(D);base=np.zeros(D)
  dg=np.zeros((D,D));dh=dg.copy();dj=dg.copy();dz=dg.copy()
  def value(x,address):
   self.records.append((address,x.copy()))
   e=noise*np.ones(D)/math.sqrt(D)
   return self.f.g(x)+e
  def first(x,address):
   check(any(k==address and np.array_equal(y,x) for k,y in self.records),'HVP first only at recorded VALUE')
   self.derivative_sites.append(address)
   return self.f.dg(x)
  for i,(ri,ci,wi) in enumerate(zip(r,c,w)):
   x=ri*Z+ci*G;I2=np.zeros(D);Ix=np.zeros((D,D));IH=Ix.copy();IJ=Ix.copy()
   for j,(tj,dj0,vj) in enumerate(zip(t,dmid,v)):
    y=tj*x+dj0*H;L=np.zeros(D);Ly=np.zeros((D,D));LJ=Ly.copy()
    for k,(sk,ek,uk) in enumerate(zip(sigma,ein,u)):
     z=sk*y+ek*J;addr=('inner',i,j,k)
     L+=uk*value(z,addr)
     if derivs:
      B=first(z,addr);Ly+=uk*sk*B;LJ+=uk*ek*B
    addr=('mid',i,j);a=y-L;I2+=vj*value(a,addr)
    if derivs:
     B=first(a,addr);Ix+=vj*tj*B@(Id-Ly);IH+=vj*dj0*B@(Id-Ly);IJ-=vj*B@LJ
   addr=('out',i);a=x-I2;F=value(a,addr);addr0=('baseline',i);B0v=value(x,addr0)
   out+=wi*F;base+=wi*B0v
   if derivs:
    B1=first(a,addr);B0=first(x,addr0);K=B1-B0-B1@Ix
    dg+=wi*ci*K;dz+=wi*ri*K;dh-=wi*B1@IH;dj-=wi*B1@IJ
  return (out,base,out-base,(dg,dh,dj,dz)) if derivs else (out,base,out-base)

rows=[];nonzero_skew=[];nonzero_origin=[];max_fd=0
for D in [1,2,5,11]:
 for A in [.015,.08,.25,.5]:
  f=Fixture(D,A);src=Source(f)
  for trial in range(3):
   Z,G,H,J=rng.normal(size=(4,D))
   F,B,E,(EG,EH,EJ,EZ)=src.eval(Z,G,H,J,derivs=True)
   check(len(src.records)==len(r)*(len(t)*(len(sigma)+1)+2),'literal raw E count')
   check(len(src.derivative_sites)==len(src.records),'one first per original VALUE on full sweep')
   Ix_bound=.5*A*(1+.5*A)
   BG=bo*(A+A*Ix_bound);BH=A*A*bm*(1+.5*A);BJ=A**3*bi;BZ=.5*(A+A*Ix_bound)
   for actual,bound,label in [(EG,BG,'G'),(EH,BH,'H'),(EJ,BJ,'J'),(EZ,BZ,'captured Z')]:
    check(op(actual)<=bound+1e-13,'actual '+label+' first')
   check(op(EG-EG.T)<=2*bo*A*Ix_bound+1e-13,'GG skew small')
   M=np.zeros((3*D,3*D));M[:D,:D]=EG;M[:D,D:2*D]=EH;M[:D,2*D:]=EJ
   check(op(M-M.T)<=2*bo*A*Ix_bound+math.sqrt(BH*BH+BJ*BJ)+1e-13,'literal 3D lifted curl')
   nonzero_skew.append(op(M-M.T))
   for axis,mat in [(0,EZ),(1,EG),(2,EH),(3,EJ)]:
    direction=rng.normal(size=D);direction/=np.linalg.norm(direction);hh=2e-6
    args=[Z,G,H,J];plus=[a.copy() for a in args];minus=[a.copy() for a in args]
    plus[axis]+=hh*direction;minus[axis]-=hh*direction
    fd=(src.eval(*plus)[2]-src.eval(*minus)[2])/(2*hh)
    err=np.max(np.abs(fd-mat@direction));max_fd=max(max_fd,float(err));check(err<3e-9,'literal complete finite-difference first')
   z0=np.zeros(D);Fo,Bo,Eo=src.eval(Z,z0,z0,z0)
   check(np.linalg.norm(Bo)<=A*.5*np.linalg.norm(Z)+1e-13,'baseline origin amplitude')
   check(np.linalg.norm(Eo)<=.25*A*A*(1+.5*A)*np.linalg.norm(Z)+1e-13,'raw E captured origin amplitude')
   nonzero_origin.append(np.linalg.norm(Eo))
   close(src.eval(z0,z0,z0,z0)[0],z0,0,'coherent total zero')
   eps=1e-6;Fn,Bn,En=src.eval(Z,G,H,J,noise=eps)
   check(np.linalg.norm(Fn-F)<=(1+A+A*A)*eps+1e-12,'absolute F precision floor')
   check(np.linalg.norm(Bn-B)<=eps+1e-12,'absolute B precision floor')
   check(np.linalg.norm(En-E)<=(2+A+A*A)*eps+1e-12,'absolute E precision floor')
   rows.append({'D':D,'A':A,'full_first':float(op(np.hstack([EG,EH,EJ]))),'curl':float(op(M-M.T)),'captured_first':float(op(EZ))})
check(max(nonzero_skew)>1e-4,'noncommuting fixture exposes real curl')
check(min(nonzero_origin)>1e-10,'fixed nonzero Z origin is not zero')
metrics['raw_source']={'records_per_E':len(r)*(len(t)*(len(sigma)+1)+2),'Q_B':len(r),'maximum_finite_difference_error':max_fd,'minimum_nonzero_origin':float(min(nonzero_origin)),'ports':rows}

# Scalar high-frequency source: small VALUE does not imply an O(A^2) first.
high=[]
for n in [20,40,80,160]:
 A=1/(8*math.pi*n)
 def g(x):return .5*A*x+.5*A*A*(1-np.cos(x/A))
 def gp(x):return .5*A*(1+np.sin(x/A))
 x=1.;y=.5*x;L=g(.5*y);I2=g(y-L);Ix=.5*gp(y-L)*(1-.5*gp(.5*y))
 first=math.sqrt(.75)*(gp(x-I2)-gp(x)-gp(x-I2)*Ix)
 high.append({'A':A,'abs_G_first_over_A':abs(first)/A,'abs_E_over_A2':abs(g(x-I2)-g(x))/A**2})
 check(abs(first)/A>.08,'actual O(A) first cannot be replaced by O(A^2)')
metrics['high_frequency_first']=high

# Matrix quadratic calibration with complete original ancestors.
for D in [1,3,7]:
 R=rng.normal(size=(D,D));K=R@R.T;K*=.3/op(K)
 class Quadratic:
  def g(self,x):return K@x
  def dg(self,x):return K
 src=Source(Quadratic());Z=rng.normal(size=D);zero=np.zeros(D)
 F,B,E=src.eval(Z,zero,zero,zero)
 close(F,(.5*K-.25*K@K+.125*K@K@K)@Z,4e-14,'matrix quadratic m3 mean')
 close(B,.5*K@Z,4e-14,'matrix quadratic baseline mean')

# Symbolic Brownian clock integrals and Stein orientation.
ss,tt,aa=sy.symbols('s t A',positive=True)
check(sy.integrate(ss**2*sy.exp(-2*ss)/2,(ss,0,sy.oo))==sy.Rational(1,8),'innovation covariance constant')
ct=2*sy.integrate(sy.integrate(sy.exp(-tt-ss)*(sy.exp(-(tt-ss))-sy.exp(-tt-ss)),(tt,ss,sy.oo)),(ss,0,sy.oo))
check(ct==sy.Rational(1,4),'conditional Markov covariance mass')
Tbound=sy.integrate(aa**2*ss*sy.exp(-2*ss)/2*(aa+aa**2*(ss+sy.Rational(1,2))),(ss,0,sy.oo))
check(sy.simplify(Tbound-aa**3/8-3*aa**4/16)==0,'regularized current operator majorant')
q,ww=roots_hermitenorm(24);ww/=math.sqrt(2*math.pi)
# F=(X,X^2-1), phi(F)=F1 F2: first row distinguishes tau from tau^T.
F1=q;F2=q*q-1
lhs=np.sum(ww*F1*(F1*F2))
rhs=np.sum(ww*((F2)+(2*q)*F1))
wrong=np.sum(ww*((F2)+q*F1))
close(lhs,2,2e-13,'Stein first-row LHS');close(rhs,lhs,2e-13,'correct Riesz orientation');close(wrong,1,2e-13,'reversed orientation is wrong')
metrics['stein_orientation']={'lhs':float(lhs),'correct':float(rhs),'transposed':float(wrong)}
metrics['current_operator_bound']='A^3/8 + 3 A^4/16; HS <= sqrt(D) times this bound'

# Deterministic Gaussian quadrature diagnostics of complete currents.
def gauss_grid(n,dim):
 z,wg=roots_hermitenorm(n);wg=wg/math.sqrt(2*math.pi)
 ids=np.array(list(itertools.product(range(n),repeat=dim)))
 return z[ids],np.prod(wg[ids],axis=1)
def gl(n,left=0,right=1):
 z,wg=leggauss(n);return left+(right-left)*(z+1)/2,wg*(right-left)/2
A=.28
def g(x):return A*(.6*x+.3*np.logaddexp(x,-x)-.3*math.log(2))
def gp(x):return A*(.6+.3*np.tanh(x))
def gpp(x):return A*.3/np.cosh(x)**2
nodes=np.array([.2,.8]);weights=np.array([.5,.5]);ds=np.sqrt(1-nodes*nodes)
KM=np.minimum.outer(nodes,nodes)/np.maximum.outer(nodes,nodes)-np.outer(nodes,nodes)
KC=np.outer(ds,ds);DK=KM-KC;L=np.linalg.cholesky(KM)
close(np.diag(DK),0,1e-15,'covariance interpolation diagonal increments vanish')
check(np.min(np.linalg.eigvalsh(KM))>0,'actual Markov clock covariance positive')
check(np.min(np.linalg.eigvalsh(DK))<0<np.max(np.linalg.eigvalsh(DK)),'scalar pointwise sign is not Loewner order')
B,WB=gauss_grid(22,2);H,WH=gauss_grid(22,1);H=H[:,0]
rhos,wrho=gl(22);alphas,wal=gl(44,0,math.pi/2)
convolution=[]
for x in [-.7,.3,1.4]:
 # H-label here denotes the two-clock Markov bank; Q is a distinct one-node rule.
 VH=B@L.T+nodes*x;IH=g(VH)@weights;DH=(gp(VH)*weights)@L
 VQ=.5*x+math.sqrt(.75)*H;IQ=g(VQ);DQ=math.sqrt(.75)*gp(VQ)
 muH=WB@IH;muQ=WH@IQ;xh=IH-muH;xq=IQ-muQ
 tauH=np.zeros(len(B));tauQ=np.zeros(len(H))
 for rho,weight in zip(rhos,wrho):
  Br=rho*B[:,None,:]+math.sqrt(1-rho*rho)*B[None,:,:]
  VR=Br@L.T+nodes*x
  mean_DR=np.einsum('b,abi->ai',WB,(gp(VR)*weights)@L)
  tauH+=weight*np.einsum('ai,ai->a',mean_DR,DH)
  Hr=rho*H[:,None]+math.sqrt(1-rho*rho)*H[None,:]
  mean_DQ=gp(.5*x+math.sqrt(.75)*Hr)@WH*math.sqrt(.75)
  tauQ+=weight*mean_DQ*DQ
 covH=WB@(xh*xh);covQ=WH@(xq*xq)
 close(WB@tauH,covH,2e-7,'mean Stein coefficient = Markov covariance')
 close(WH@tauQ,covQ,2e-7,'mean Stein coefficient = cheap covariance')
 joint=WB[:,None]*WH[None,:];current=0.;derivative=0.;meanterm=0.
 for alpha,weight in zip(alphas,wal):
  sn,cs=math.sin(alpha),math.cos(alpha);theta=sn*sn
  T=x-((1-theta)*muQ+theta*muH)-sn*xh[:,None]-cs*xq[None,:]
  delta=muH-muQ
  mean=-2*sn*cs*delta*np.sum(joint*gp(T))
  cur=sn*cs*np.sum(joint*gpp(T)*(tauH[:,None]-tauQ[None,:]))
  deriv=-np.sum(joint*gp(T)*(2*sn*cs*delta+cs*xh[:,None]-sn*xq[None,:]))
  meanterm+=weight*mean;current+=weight*cur;derivative+=weight*deriv
 target=WB@g(x-IH)-WH@g(x-IQ)
 close(derivative,target,3e-9,'centered convolution derivative-free identity')
 close(meanterm+current,target,3e-7,'complete centered Stein current identity')
 convolution.append({'x':x,'delta_mu':float(muH-muQ),'target':float(target),'mean_term':float(meanterm),'stein_current':float(current),'derivative_free_error':float(abs(derivative-target)),'Stein_error':float(abs(meanterm+current-target))})
metrics['centered_convolution']=convolution

# Shared finite-clock Gaussian interpolation retains the complete nonlinear sites.
B3,W3=gauss_grid(22,3);YM=B3[:,:2]@L.T;YC=B3[:,2,None]*ds
lams,wl=gl(36);interp=[]
for x in [-.7,.3,1.4]:
 IM=g(nodes*x+B@L.T)@weights;IC=g(nodes*x+H[:,None]*ds)@weights
 muM=WB@IM;muC=WH@IC
 target=WB@g(x-IM)-WH@g(x-IC)
 covtarget=WB@(IM*IM)-(WB@IM)**2-(WH@(IC*IC)-(WH@IC)**2)
 current=0.;covcurrent=0.
 for lam,weight in zip(lams,wl):
  V=nodes*x+math.sqrt(lam)*YM+math.sqrt(1-lam)*YC
  Ip=g(V)@weights;Js=gp(V)*weights
  coeff=np.einsum('ni,ij,nj->n',Js,DK,Js)
  current+=weight*.5*np.sum(W3*gpp(x-Ip)*coeff)
  covcurrent+=weight*np.sum(W3*coeff)
 close(muM,muC,4e-7,'interpolation preserves exact conditional mean')
 close(current,target,4e-7,'shared-clock nonlinear current identity')
 close(covcurrent,covtarget,5e-7,'shared-clock covariance identity')
 interp.append({'x':x,'mean_difference':float(muM-muC),'target':float(target),'current':float(current),'covariance_target':float(covtarget),'covariance_current':float(covcurrent)})
metrics['same_clock_interpolation']=interp

# The sign and index convention in the regularized innovation formula, on a smooth scalar Gaussian fixture.
z,wg=roots_hermitenorm(60);wg/=math.sqrt(2*math.pi);ths,wth=gl(36)
x=.7;a=.2;b=.03
left=np.sum(wg*(g(x-(a+b)*z)-g(x-a*z)));ftc=0.;ibp=0.
for theta,weight in zip(ths,wth):
 f=a+theta*b
 ftc-=weight*np.sum(wg*gp(x-f*z)*b*z)
 ibp+=weight*np.sum(wg*gpp(x-f*z))*b*f
close(left,ftc,1e-13,'innovation first-derivative current sign')
close(left,ibp,1e-13,'innovation IBP current sign')
metrics['innovation_current_fixture']={'difference':float(left),'first_current':float(ftc),'regularized_current':float(ibp),'scope':'Finite Gaussian sign/index diagnostic, not an executed OU source.'}

# Exact Hermite quantities used in the separately imported separator review.
for D in [1,2,5,20,1000]:
 check(18+2*(D-1)+4*(D-1)==6*D+12,'gradient Hermite energy')
 check(6+2*(D-1)==2*D+4,'gradient Laplacian trace')
 A=D**(-.5);check(abs((A**4*D)/(A**4*math.sqrt(D))-math.sqrt(D))<1e-12,'trace amplification ratio')
metrics['separator_scope']='Both verified/imported separator claims are bounded: small-feedback unshifted Dg is false; independent E,W ports are insufficient. Neither disproves the same-g single-I covariance formula.'

result={'status':'PASS bounded supplementary checks; m3 mean-bias closure OPEN','assertions':checks,'seed':202610050028,'source_pins':PINS,'scope':'Raw own-mean ports/counts, exact analytical currents, and separator scope only. No fourth-order target-bias bound, native current consumer, endpoint join, or all-order result.','metrics':metrics}
(HERE/'m3_gate_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'assertions':checks,'max_raw_first_difference':max_fd,'max_centered_Stein_error':max(a['Stein_error'] for a in convolution),'max_interpolation_error':max(abs(a['target']-a['current']) for a in interp)},indent=2))
