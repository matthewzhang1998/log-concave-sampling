#!/usr/bin/env python3
"""Coefficient/source-graph diagnostics, not execution of the imported large compilers."""
import json,math,os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
OUT=Path(__file__).resolve().parent
checks=0
def check(ok,msg):
 global checks
 if not ok: raise AssertionError(msg)
 checks+=1
def sym(A):return (A+A.T)/2
def cov(x,y,w):
 xc=x-w@x;yc=y-w@y
 return (xc*w[:,None]).T@yc
def rule(K,m):
 xx,ww=leggauss(m);rr=[];weights=[]
 for k in range(K):
  a=2.**(-k-1);rr.extend(1-(1.5*a+.5*a*xx));weights.extend(.5*a*ww)
 rr.append(1-2.**(-K-1));weights.append(2.**(-K))
 return np.array(rr),np.array(weights)
r,w=rule(8,5);q=np.sqrt(1-r*r);beta=w@q;cc=.25+beta*beta
check(abs(w.sum()-1)<1e-14,'mass');check(abs(w@r-.5)<1e-14,'moment');check(r.min()>=1/(4*25),'lift coefficient lower bound')
L=float(w@(1/r));check(L<=100,'lift mass publiclog')
rng=np.random.default_rng(810229);quadratic=[]
for D in (1,2,5,17):
 for A in (.5,.25,.125,.0625):
  O,_=np.linalg.qr(rng.normal(size=(D,D)));B=(O*rng.uniform(.05*A,A,D))@O.T
  for s in (1.,.5,.1,.01):
   Delta=min(.25,s*s);I=np.eye(D);P=np.concatenate((I,np.zeros_like(I)),axis=1)
   Ku=.5*s*s*B;Kv=beta*s*s*B;K=np.concatenate((Ku,Kv),axis=1)
   QQ=s*(P-K);f=B@(P-K);b=B@P
   HH=sym(QQ@(B@QQ).T)/(s*s)
   Hstar=sym(f[:,:D])-sym(K@b.T)
   exact=B@np.linalg.inv(I+s*s*B)
   check(np.linalg.norm(Hstar-(B-s*s*B@B))<1e-13,'exact full matrix Hstar')
   check(np.linalg.norm(HH-(B-s*s*B@B+cc*s**4*B@B@B))<1e-13,'all shared-root covariance')
   check(np.linalg.norm(Hstar-exact)<=A**3*s**4*math.sqrt(D)+1e-13,'quadratic curvature grade')
   GH=np.zeros((2*D,2*D))
   for ri,qi,wi in zip(r,q,w):
    Ri=np.concatenate((ri*I,qi*I),axis=1)
    GH+=(wi/ri)*Ri.T@(s*s*B)@Ri
   GB=P.T@b;Jp=s*GB+GH/s;Jm=s*GB-GH/s
   check(np.linalg.norm(P@GH-K)<1e-13,'selected exact gradient lift')
   check(np.linalg.norm(GH-GH.T)<1e-13,'full gradient lift symmetric')
   Cmix=.25*P@(Jp@Jp.T-Jm@Jm.T)@P.T
   check(np.linalg.norm(Cmix-sym(K@b.T))<1e-13,'balanced covariance polarization')
   vD=(1-Delta)/2;vp=(1-Delta)/4
   M=Delta/vD*P.T@f
   CD=vD*I-Delta*sym(P@M@P.T)*(vD/Delta)
   # Literal joint pair covariance formula with p retained.
   joint=np.block([[np.eye(2*D),M.T],[M,np.eye(2*D)]])
   read=math.sqrt(vD/2)*np.concatenate((P,-P),axis=1)
   check(np.linalg.norm(read@joint@read.T-(vD*I-Delta*sym(f[:,:D])))<1e-13,'retained-input linear covariance')
   Cstar=(1-Delta)*I-Delta*Hstar;Ctrue=(1-Delta)*I-Delta*exact
   check(np.linalg.eigvalsh(Cstar).min()>.49,'total positive covariance')
   check(np.linalg.eigvalsh(Ctrue).min()>=.625-1e-13,'true reserve gap')
   check(np.linalg.norm(Cstar-Ctrue)<=Delta*A**3*s**4*math.sqrt(D)+1e-13,'reserve grade')
   quadratic.append(dict(D=D,A=A,s=s,Delta=Delta,error=float(np.linalg.norm(Hstar-exact))))
# Algebraic heat bounds over many scales, in log domain to avoid underflow.
for A in (.5,.1,.01,1e-6):
 for s in (1.,.3,.01,1e-4):
  for fract in (1.,.4,.01):
   dd=fract*min(.25,s*s);target=dd*A**3*s**4
   pair=2*dd**2*A**3*s*s+dd**4*A**4*(1+1/(math.sqrt(A)*s))
   action=dd*A**3*s**4+dd**2*A**4*s**4+dd**2*A**3.5*s**3
   check(pair<=4*target,'pair heat powers');check(action<=3*target,'action heat powers')
# Genuine C2, non-C3 ridge fixture, with noncommuting Hessians.
def h(t):
 u=np.abs(t);v=np.sqrt(u);ans=u-2*v+2*np.log1p(v);small=u<1e-4
 if np.any(small):
  a=u[small];vv=np.zeros_like(a)
  for k in range(1,20):vv+=(-1)**(k+1)*a**(1+k/2)/(1+k/2)
  ans=np.array(ans);ans[small]=vv
 return np.sign(t)*ans
def hp(t):
 v=np.sqrt(np.abs(t));return v/(1+v)
def psi(t):
 u=np.abs(t);v=np.sqrt(u);ans=.5*u*u-(4/3)*u*v+2*(u-1)*np.log1p(v)-u+2*v;small=u<1e-3
 if np.any(small):
  a=u[small];vv=np.zeros_like(a)
  for k in range(1,21):vv+=(-1)**(k+1)*a**(2+k/2)/((1+k/2)*(2+k/2))
  ans=np.array(ans);ans[small]=vv
 return ans
axes=np.array([[1.,0.],[1/math.sqrt(2),1/math.sqrt(2)]])
def gv(x,A):return A*(.2*x+.3*h(x@axes.T)@axes)
def pot(x,A):return A*(.1*np.sum(x*x,axis=-1)+.3*np.sum(psi(x@axes.T),axis=-1))
def hg(x,A):
 dots=x@axes.T;ans=np.broadcast_to(.2*np.eye(2),x.shape[:-1]+(2,2)).copy()
 for j in range(2):ans+=.3*hp(dots[...,j])[...,None,None]*np.outer(axes[j],axes[j])
 return A*ans
def gaussian_grid(n,d):
 gh,gw=hermgauss(n);gh*=math.sqrt(2);gw/=math.sqrt(math.pi)
 grid=np.stack(np.meshgrid(*([gh]*d),indexing='ij'),axis=-1).reshape(-1,d)
 wg=np.prod(np.stack(np.meshgrid(*([gw]*d),indexing='ij'),axis=-1),axis=-1).reshape(-1)
 return grid,wg
nonlinear=[]
for n in (12,18):
 tape,ww=gaussian_grid(n,4);u=tape[:,:2];v=tape[:,2:];tu,tw=gaussian_grid(n,2)
 for A in (.25,.125,.0625):
  for s in (1.,.5):
   a=np.array([.31,-.27]);x=a.copy()
   for _ in range(16):x=a-s*s*gv(x,A)
   g0=gv(x,A);Rmode=x-a+s*s*g0
   KK=np.zeros_like(u);GHH=np.zeros_like(tape)
   for ri,qi,wi in zip(r,q,w):
    bar=s*(gv(x+s*(ri*u+qi*v),A)-g0)
    KK+=wi*bar
    GHH[:,:2]+=wi*bar;GHH[:,2:]+=wi*qi/ri*bar
   QQ=x+s*(u-KK);FF=gv(QQ,A);BB=gv(x+s*u,A);ff=(FF-g0)/s;bb=(BB-g0)/s
   GB=np.concatenate((bb,np.zeros_like(bb)),axis=1)
   JP=s*GB+GHH/s;JM=s*GB-GHH/s
   Cpol=.25*(cov(JP,JP,ww)-cov(JM,JM,ww))[:2,:2]
   check(np.linalg.norm(Cpol-sym(cov(KK,bb,ww)))<1e-13,'C2 same-root polarization')
   HQ=sym(cov(QQ,FF,ww))/(s*s)
   HS=sym(cov(u,ff,ww))-sym(cov(KK,bb,ww))
   paid=sym(cov(KK,(FF-BB)/s,ww))
   check(np.linalg.norm(HS-HQ-paid)<1e-13,'C2 retained mixed-tail identity')
   check(np.linalg.norm(paid)<=3*A**3*s**4*math.sqrt(2),'C2 paid tail grade')
   targetX=x+s*tu
   lw=-pot(targetX,A)+pot(x,A)+s*(tu@g0)-(tu@Rmode)/s
   targetw=tw*np.exp(lw-lw.max());targetw/=targetw.sum()
   true=sym(cov(targetX,gv(targetX,A),targetw))/(s*s)
   error=float(np.linalg.norm(HS-true));ratio=error/(A**3*s**4*math.sqrt(2))
   # These finite integration checks are diagnostic; the proof is analytic.
   check(ratio<1.5,'C2 diagnostic reserve grade')
   comm=hg(np.array([[.1,2.],[2.,.1]]),A)
   check(np.linalg.norm(comm[0]@comm[1]-comm[1]@comm[0])>1e-6*A*A,'noncommuting C2 fixture')
   # Actual C1 source first/skew and exact gradient-lift Jacobians at arbitrary tape points.
   for idx in (0,len(u)//3,len(u)//2):
    Kfirst=np.zeros((2,4));GHfirst=np.zeros((4,4))
    for ri,qi,wi in zip(r,q,w):
     Ri=np.concatenate((ri*np.eye(2),qi*np.eye(2)),axis=1);HH=s*s*hg(x+s*(Ri@tape[idx]),A)
     Kfirst+=wi*HH@Ri;GHfirst+=(wi/ri)*Ri.T@HH@Ri
    Ffirst=hg(QQ[idx],A)@(np.concatenate((np.eye(2),np.zeros((2,2))),axis=1)-Kfirst)
    lift=np.vstack((Ffirst,np.zeros((2,4))))
    check(np.linalg.norm(lift-lift.T,2)<=4*A*A*s*s,'C2 complete small curl')
    check(np.linalg.norm(GHfirst-GHfirst.T)<1e-13,'C2 fullgradient Jacobian')
    check(np.linalg.norm(GHfirst,2)<=A*s*s*L+1e-13,'C2 gradient-lift radius')
   nonlinear.append(dict(n=n,A=A,s=s,curvature_error=error,error_over_A3s4sqrtD=ratio,paid_tail=float(np.linalg.norm(paid)),mode_residual=float(np.linalg.norm(Rmode))))
report=dict(status='PASS',assertions=checks,quadrature_nodes=len(r),lift_mass=L,quadratic=quadratic,nonlinear_C2=nonlinear,scope='Diagnostics verify raw VALUE-graph algebra, covariance targets, matrix positivity and C2 source ports. Imported finite pair/covariance compilers are not end-to-end numerically instantiated. C2 quadrature values are diagnostics, not a substitute for analytic one-energy bounds.')
(OUT/'value_curvature_reserve_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status='PASS',assertions=checks,max_C2_ratio=max(a['error_over_A3s4sqrtD'] for a in nonlinear),quadrature_nodes=len(r),lift_mass=L),indent=2))
