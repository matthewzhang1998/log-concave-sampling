from __future__ import annotations
import hashlib,itertools,json,math,os
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from active_probe_source import active_source,radial_pullback,projection_cut_constant,response_floor
HERE=Path(__file__).resolve().parent
checks=0

def ck(x,msg):
 global checks
 checks+=1
 if not x:raise AssertionError(msg)

def gm(n):
 if n%2:return 0
 return math.prod(range(1,n,2))

def comps(n,k):
 if k==1:yield (n,);return
 for j in range(n+1):
  for tail in comps(n-j,k-1):yield(j,)+tail

def source_polynomial_value(m,k,z,x,P,b=F(1),At=F(1)):
 out=F(0)
 for eps in itertools.product((-1,1),repeat=k):
  q=z+b*sum(e*p for e,p in zip(eps,P))
  out+=math.prod(eps)*((q+b*x)**m-q**m)
 return out/(At*2**k)

def parity_expanded(m,k,z,x,P,b=F(1),At=F(1)):
 # Multinomial coordinates are z,x,P1,...,Pk. Exact common-bank grouping.
 out=F(0)
 for rr in comps(m,k+2):
  rz,rx,*rp=rr
  if rx==0 or any(v%2!=1 for v in rp):continue
  multi=math.factorial(m)//math.prod(math.factorial(v) for v in rr)
  out+=multi*z**rz*b**(m-rz)*x**rx*math.prod(p**r for p,r in zip(P,rp))
 return out/At

def source_response_monomial(m,k,z,b=F(1),At=F(1)):
 out=F(0)
 for rr in comps(m,k+2):
  rz,rx,*rp=rr
  if rx==0 or any(v%2!=1 for v in rp):continue
  multi=math.factorial(m)//math.prod(math.factorial(v) for v in rr)
  out+=multi*z**rz*b**(m-rz)*gm(rx+1)*math.prod(gm(r+1) for r in rp)
 return out/At

def target_response_monomial(m,k,z,b=F(1),At=F(1)):
 d=k+1
 if m<d:return F(0)
 deg=m-d
 expectation=sum(F(math.comb(deg,j))*z**(deg-j)*b**j*F(k+1)**(j//2)*gm(j) for j in range(0,deg+1,2))
 return F(math.factorial(m),math.factorial(m-d))*b**d/At*expectation

exact_rows=[]
for k in range(0,8):
 for m in range(0,k+6):
  for z in (F(0),F(2,3)):
   got=source_response_monomial(m,k,z)
   exp=target_response_monomial(m,k,z)
   ck(got==exp,f'Gaussian response k={k} m={m} z={z}')
  # A deterministic literal signed query test, independent of expectation.
  P=[F((-1)**i*(i+1),i+2) for i in range(k)]
  got=source_polynomial_value(m,k,F(2,3),F(3,5),P)
  exp=parity_expanded(m,k,F(2,3),F(3,5),P)
  ck(got==exp,f'parity source k={k} m={m}')
 exact_rows.append({'k':k,'response_rank':k+2,'raw_value_calls':2**(k+1),'all_cut_constant':max(projection_cut_constant(k,p) for p in range(1,k+2))})

rng=np.random.default_rng(782340)
numrows=[]
# A non-coordinate-separable, globally admissible gradient. Its Hessian is
# .5 I+.2 cos(u.x) uu^T+.1 cos(v.x) vv^T, lying in [.2,.8] I.
for D in (1,2,3,8,32):
 u=rng.normal(size=D);u/=np.linalg.norm(u)
 v=rng.normal(size=D);v/=np.linalg.norm(v)
 def g(q):return .5*q+.2*math.sin(float(u@q))*u+.1*math.sin(float(v@q))*v
 for k in (0,1,2,4,6):
  for R in (None,math.sqrt(D)+.7):
   for rep in range(4):
    x=rng.normal(size=D)*(1 if R is None else 2)
    P=rng.normal(size=(k,D));z=rng.normal(size=D);t=.05 if rep%2 else .7
    value,calls=active_source(g,x,P,z,t,1.,R)
    ck(calls==2**(k+1),'call ledger')
    ck(np.linalg.norm(value)<=np.linalg.norm(x)/math.sqrt(k+1)+1e-10,'one-energy envelope')
    if R is not None:ck(np.linalg.norm(value)<=1.5*R/math.sqrt(k+1)+1e-10,'bounded radius')
    zero,_=active_source(g,np.zeros(D),P,z,t,1.,R)
    ck(np.array_equal(zero,np.zeros(D)),'exact origin')
    for j in range(k):
     P2=P.copy();P2[j]*=-1
     neg,_=active_source(g,x,P2,z,t,1.,R)
     ck(np.max(np.abs(neg+value))<3e-12,'probe parity')
    # Direct finite differences check literal source curl and complete first.
    # Full dimensions are used up to D=8; at D=32 use directional firsts.
    if D<=8:
     h=2e-5
     J=[]
     for j in range(D):
      e=np.eye(D)[j]*h
      fp,_=active_source(g,x+e,P,z,t,1.,R)
      fm,_=active_source(g,x-e,P,z,t,1.,R)
      J.append((fp-fm)/(2*h))
     J=np.asarray(J).T
     ck(np.linalg.norm(J-J.T,ord=2)<3e-7,'private curl')
     ck(np.linalg.norm(J,ord=2)<(13 if R else 1)/math.sqrt(k+1)+2e-7,'private first')
    for j in range(3):
     dx=rng.normal(size=D);dP=rng.normal(size=(k,D));norm=math.sqrt(float(dx@dx)+float(np.sum(dP*dP)))
     dx/=norm;dP/=norm;h=2e-5
     fp,_=active_source(g,x+h*dx,P+h*dP,z,t,1.,R)
     fm,_=active_source(g,x-h*dx,P-h*dP,z,t,1.,R)
     ck(np.linalg.norm((fp-fm)/(2*h))<=(13 if R else 1)+3e-7,'complete private/probe first')
    dz=rng.normal(size=D);dz/=np.linalg.norm(dz);h=2e-5
    fp,_=active_source(g,x,P,z+h*dz,t,1.,R)
    fm,_=active_source(g,x,P,z-h*dz,t,1.,R)
    ck(np.linalg.norm((fp-fm)/(2*h))<=1/t+2e-7,'caller first')
    if k:
     lin,_=active_source(lambda q:.5*q,x,P,z,t,1.,R)
     ck(np.linalg.norm(lin)<3e-12,'linear exact cancellation up to arithmetic')
   numrows.append({'D':D,'k':k,'bounded':R is not None})

# Radial geometry, including seam neighborhoods.
for D in (1,2,5):
 for R in (.3,1,7):
  for r in (0,.9*R,R,1.01*R,1.5*R,1.99*R,2*R,3*R):
   x=np.zeros(D);x[0]=r
   y,J=radial_pullback(x,R)
   ck(np.linalg.norm(J,ord=2)<=1+1e-12,'radial first')
   ck(np.linalg.norm(y)<=min(r,1.5*R)+1e-12,'radial radius')
   ck(np.max(np.abs(J-J.T))<1e-12,'radial symmetry')
   if r>0:
    direction=np.ones(D)/math.sqrt(D);h=1e-6*R
    _,jp=radial_pullback(x+h*direction,R);_,jm=radial_pullback(x-h*direction,R)
    ck(np.linalg.norm((jp-jm)/(2*h),ord=2)<=8/R+1e-4,'radial second')

T={1:1};S={1:1};V={1:1};grades=[]
for n in range(2,25):
 T[n]=sum(i*(n-i)*T[i]*T[n-i] for i in range(1,n))
 S[n]=sum(math.comb(n,i)*i*(n-i)*S[i]*S[n-i] for i in range(1,n))
 V[n]=sum((n-i)*(i+1)*V[i]*T[n-i]+i*(n-i+1)*T[i]*V[n-i] for i in range(1,n))
 if n>=4:
  gamma=F(1,n);beta=F(1,2*(n-1));root=F(n)-F(1,2)-F(n-2)*gamma
  ck(root==F(n)-F(3,2)+F(2,n),'root grade')
  ck(F(n)+gamma==F(n+1)-F(n-1)*gamma,'heat/one-hit balance')
  ck(2*root>F(n)+gamma,'hypothetical quadratic return surplus')
  ck(root-gamma>1,'actual caller first surplus')
  ck(2*root-F(n*n-3*n+3,n)==F(n)+gamma,'charged extra-width admission edge')
  grades.append({'rank':n,'beta':str(beta),'tau_exponent':str(gamma),'root_grade':str(root),'root_caller_grade':str(root-gamma),'own_quadratic_grade_CONDITIONAL_extra_width_ell_zero':str(2*root),'maximum_additional_width_exponent_ell':n*n-3*n+3,'joined_grade_CONDITIONAL':str(F(n)+gamma),'T_n':T[n],'S_n':S[n],'V_n_raw_VALUE_census':V[n]})
ck(V[8]==27020800,'rank8 raw VALUE census');ck(T[8]==794880,'rank8 history count');ck(S[8]==32049561600,'rank8 weight sum')
ck(F(8)-F(3,2)+F(2,8)==F(27,4),'rank8 root rational')

report={'status':'PASS','assertions':checks,'scope':'Exact active-source response identities; executed VALUE source, origin/parity/radius/first/curl diagnostics; rational CONDITIONAL admission ledger. No native pair/filter compiler or arbitrary graph return executed.','exact_source_rows':exact_rows,'numerical_source_rows':numrows,'grade_rows':grades}
(HERE/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'PASS','assertions':checks,'exact_ranks':[r['response_rank'] for r in exact_rows],'rank8_grade':next(r for r in grades if r['rank']==8)},indent=2))
