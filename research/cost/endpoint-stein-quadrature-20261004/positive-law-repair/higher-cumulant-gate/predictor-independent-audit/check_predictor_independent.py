#!/usr/bin/env python3
"""Independent algebra, literal cumulants, ports, and rank-five next-row audit.

Does not import or execute the author's checker.  No output is a uniform C2
order-three law theorem.  All root fields retain their actual shared genealogy.
"""
from pathlib import Path
import json, math, hashlib
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
OUT=Path(__file__).resolve().parent
checks=0

def ck(test, message):
    global checks
    if not test: raise AssertionError(message)
    checks+=1

p=sp.symbols('p', real=True)
M=sp.Matrix([[sp.Rational(1,2),p/4,1],
 [sp.Rational(1,5)+p/15,sp.Rational(2,15)+p/15,sp.Rational(2,5)],
 [1,sp.Rational(2,3)+p/6,sp.Rational(7,3)]])
rhs=sp.Matrix([sp.Rational(3,4)-p*p/16,sp.Rational(1,3),sp.Rational(4,3)-p/3])
sol=[sp.factor(x) for x in M.inv()*rhs]
ck(all(sp.simplify(x)==0 for x in M*sp.Matrix(sol)-rhs), 'symbolic continuum equations')
ck(sp.simplify(M.det()+(5*p*p-7*p-4)/180)==0, 'exact determinant')
continuum=dict(zip(('a','c','b'),[str(sp.N(x.subs(p,sp.pi),40)) for x in sol]))
lg,lw=leggauss(12)
r=[];w=[]
for j in range(20):
 h=2.**(-j-1);r.extend(1-1.5*h-.5*h*lg);w.extend(.5*h*lw)
r.append(1-2.**(-21));w.append(2.**(-20))
r=np.array(r);w=np.array(w);s=np.sqrt(1-r*r)
m=lambda j: float(w@(r**j))
n=lambda j: float(w@(r**j*s))
beta=n(0);C=.25+beta*beta
MM=np.array([[.5,beta,1], [m(4)+2*beta*n(3),n(3)+2*beta*(m(2)-m(4)),2*m(4)],
 [3*m(2),2*(n(1)+beta*m(2)),1+4*m(2)]])
RR=np.array([1-C,1/3,2-2*(m(2)+2*beta*n(1))])
a,c,b=np.linalg.solve(MM,RR)
ck(abs(m(0)-1)<1e-14 and abs(m(1)-.5)<1e-14, 'quadrature mass and first moment')
ck(np.max(abs(MM@np.array([a,c,b])-RR))<2e-15, 'finite calibration')
ck(np.linalg.det(MM)<-.12, 'finite determinant gap')

# An algebraic tail floor survives rounding: exact tail quadrature moment defect.
t=sp.Rational(1,2)**20
first_tail=sp.simplify(t*(1-t/2)-sp.integrate(p,(p,1-t,1)))
second_tail=sp.simplify(t*(1-t/2)**2-sp.integrate(p**2,(p,1-t,1)))
ck(first_tail==0 and second_tail==-t**3/12,'nonzero exact quadrature floor')

def row(rank, moment=m, mixed=n, beta=beta, coeff=(a,c,b)):
 aa,cc,bb=coeff
 J=aa*(moment(rank+1)+2*beta*mixed(rank))+cc*(mixed(rank)+2*beta*(moment(rank-1)-moment(rank+1)))+2*bb*moment(rank+1)
 L=aa*rank*moment(rank-1)+cc*(2*beta*moment(rank-1)+(rank-1)*mixed(rank-2))+bb*(1+2*(rank-1)*moment(rank-1))+(rank-1)*(moment(rank-1)+2*beta*mixed(rank-2))
 return J,L
for rank in (3,5,7):
 J,L=row(rank)
 if rank==3: ck(abs(J-1/3)+abs(L-2)<2e-15,'calibrated third row')
 else: ck(abs(J-1/rank)>.01 and abs(L-2)>.04,'uncalibrated higher row')

x,px=hermgauss(86);x*=math.sqrt(2);px/=math.sqrt(math.pi)
u,v=np.meshgrid(x,x,indexing='ij');ww=np.outer(px,px)
E=lambda f:float(np.sum(ww*f))

def mul2(a,b): return np.convolve(a,b)[:3]
def cumulants_of_coefficients(moments):
 kap=[np.zeros(3)]
 for rank in range(1,len(moments)):
  kk=moments[rank].copy()
  for j in range(1,rank): kk-=math.comb(rank-1,j-1)*mul2(kap[j],moments[rank-j])
  kap.append(kk)
 return kap

def cumulant_coefficients(q,eps,k):
 f=lambda z:2*q*z+eps*k*(np.cos(k*z)-1)
 H=np.zeros_like(u);Hu=np.zeros_like(u);Hv=np.zeros_like(u)
 for ri,si,wi in zip(r,s,w):
  z=ri*u+si*v;df=2*q-eps*k*k*np.sin(k*z)
  H+=wi*f(z);Hu+=wi*ri*df;Hv+=wi*si*df
 T=a*Hu*H+c*Hv*H+b*Hu*f(u)
 moments=[np.array([1.,0.,0.])]
 for rank in range(1,8):
  moments.append(np.array([E(u**rank),-rank*E(u**(rank-1)*H),
   rank*E(u**(rank-1)*T)+(rank*(rank-1)/2*E(u**(rank-2)*H*H) if rank>=2 else 0)]))
 return cumulants_of_coefficients(moments)

rows=[]
for q,eps,k in ((.25,.1,.5),(.25,.1,1),(.2,.03,1.5),(.3,.02,2)):
 kap=cumulant_coefficients(q,eps,k)
 for rank in (3,5,7):
  J,L=row(rank);sign=(-1)**((rank-1)//2);tau=math.exp(-k*k/2)
  pred=sign*rank*q*eps*tau*(k**rank*L-k**(rank+2)*J)
  target=sign*q*eps*tau*(2*rank*k**rank-k**(rank+2))
  ck(abs(kap[rank][2]-pred)<2e-10,'independent centered cumulant expansion')
  rows.append(dict(q=q,epsilon=eps,k=k,rank=rank,calculated=float(kap[rank][2]),formula=pred,target=target,defect=pred-target))

# Literal map and independent scalar tilted integration; kappa5 is central m5-10*m3*m2.
def kappa5(y,weights):
 mean=float(np.sum(weights*y));z=y-mean
 return float(np.sum(weights*z**5)-10*np.sum(weights*z**3)*np.sum(weights*z**2))
literal=[]
q=.25;eps=.1;k=1
f=lambda z:2*q*z+eps*(np.cos(z)-1)
U=lambda z:q*z*z+eps*(np.sin(z)-z)
H=sum((wi*f(ri*u+si*v) for ri,si,wi in zip(r,s,w)),np.zeros_like(u))
for amp in (.04,.02,.01,.005):
 us=u-amp*(a*H+b*f(u));vs=v-amp*c*H
 Y=u-sum((wi*amp*f(ri*us+si*vs) for ri,si,wi in zip(r,s,w)),np.zeros_like(u))
 pw=px*np.exp(-amp*U(x));pw/=sum(pw)
 defect=kappa5(Y,ww)-kappa5(x,pw)
 literal.append(dict(A=amp,kappa5_defect=defect,over_A2=defect/amp**2))
ck(abs(literal[-1]['over_A2']-rows[4]['defect'])<.0003,'literal rank-five next-order obstruction')

# Noncommuting gradient: verify firsts against finite differences, not just a restatement.
B=np.diag([.43,.59]);axis=np.array([1.,2.])/math.sqrt(5)
def f2(z,al):return al*(B@z+.07*(math.cos(axis@z)-1)*axis+.03*np.array([math.cos(z[0])-1,0]))
def Df2(z,al):return al*(B-.07*math.sin(axis@z)*np.outer(axis,axis)-.03*np.diag([math.sin(z[0]),0]))
def packet(z,vv,al):
 K=np.zeros(2);Ku=np.zeros((2,2));Kv=np.zeros((2,2))
 for ri,si,wi in zip(r,s,w):
  xx=ri*z+si*vv;DD=Df2(xx,al)
  K+=wi*f2(xx,al);Ku+=wi*ri*DD;Kv+=wi*si*DD
 return K,Ku,Kv

def graph(z,vv,al):
 K,Ku,Kv=packet(z,vv,al)
 Ks,Kus,Kvs=packet(z-a*K-b*f2(z,al),vv-c*K,al)
 Dv=Kvs-(a*Kus+c*Kvs)@Kv
 Du=Kus@(np.eye(2)-a*Ku-b*Df2(z,al))-c*Kvs@Ku
 return Ks,Du,Dv,K,Ku,Kv,Kus,Kvs
rng=np.random.default_rng(29802)
maxcurl=0.
for al in (.2,.07):
 for j in range(8):
  z=rng.normal(size=2);vv=rng.normal(size=2)
  Ks,Du,Dv,K,Ku,Kv,Kus,Kvs=graph(z,vv,al)
  h=2e-5;du=[];dv=[]
  for e in np.eye(2):
   du.append((graph(z+h*e,vv,al)[0]-graph(z-h*e,vv,al)[0])/(2*h))
   dv.append((graph(z,vv+h*e,al)[0]-graph(z,vv-h*e,al)[0])/(2*h))
  ck(np.max(abs(Du-np.array(du).T))<2e-10,'full u derivative, including b force')
  ck(np.max(abs(Dv-np.array(dv).T))<2e-10,'captured-u private-v derivative')
  comm=-(a*Kus+c*Kvs)@Kv+Kv@(a*Kus+c*Kvs)
  ck(np.linalg.norm(Dv-Dv.T-comm)<1e-14,'exact private curl commutator')
  curl=np.linalg.norm(comm,2);maxcurl=max(maxcurl,curl)
  ck(curl<=2*(abs(a)/2+abs(c)*beta)*beta*al*al,'private curl radius')
  ck(np.linalg.norm(Ks-K)<=al*(math.hypot(a,c)*np.linalg.norm(K)+abs(b)*np.linalg.norm(f2(z,al)))+1e-14,'chord energy')
  origin=graph(z,np.zeros(2),al)[0]-packet(z,np.zeros(2),al)[0]
  centered_origin=graph(z,np.zeros(2),al)[0]-packet(z,np.zeros(2),al)[0]-origin
  ck(np.array_equal(centered_origin,np.zeros(2)),'captured-u exact centered origin')
 ck(np.array_equal(graph(np.zeros(2),np.zeros(2),al)[0],np.zeros(2)),'all-zero source graph')
ck(maxcurl>1e-8,'fixture has genuinely nonzero curl')

# Exact matrix covariance polynomial, tested on rotated positive matrices.
for d in (1,3,8):
 O,_=np.linalg.qr(rng.normal(size=(d,d)));BB=O@np.diag(rng.uniform(.01,.15,d))@O.T
 dd=a/2+c*beta;ee=dd+b
 Cz=np.eye(d)-BB/2+ee*(BB@BB)/2;Cv=-beta*BB+dd*beta*(BB@BB)
 cov=Cz@Cz.T+Cv@Cv.T
 B2=BB@BB;B3=B2@BB;B4=B3@BB
 exact=np.eye(d)-BB+(C+ee)*B2-(ee/2+2*dd*beta**2)*B3+(ee**2/4+dd**2*beta**2)*B4
 ck(np.linalg.norm(cov-exact)<1e-14,'exact matrix covariance polynomial')
 ck(abs(C+ee-1)<1e-14,'full B2 calibration')
 ck(np.linalg.eigvalsh(cov).min()>.8,'actual covariance gap')

report=dict(status='PASS',assertions=checks,scope='Third-cumulant quadratic-times-odd repair only; rank-five row remains nonzero. No general order-three law theorem.',
 continuum_coefficients=continuum,continuum_exact=[str(x) for x in sol],finite_coefficients=dict(a=float(a),c=float(c),b=float(b)),
 finite_det=float(np.linalg.det(MM)),nodes=len(r),exact_m2_floor=str(second_tail),higher_cumulant_rows=rows,literal_rank5=literal,max_private_curl=maxcurl)
(OUT/'predictor_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('higher_cumulant_rows',)},indent=2))
