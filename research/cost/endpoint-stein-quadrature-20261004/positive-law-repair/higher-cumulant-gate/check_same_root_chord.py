#!/usr/bin/env python3
"""Finite original-VALUE chord and exact higher-cumulant diagnostics.

Hessians in this file certify actual firsts only; the map itself calls g VALUES.
Exact moment matching is an intentionally favorable diagnostic, not a producer.
"""
from pathlib import Path
import hashlib, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss

OUT=Path(__file__).resolve().parent
checks=0
def ck(ok, label):
    global checks
    if not ok: raise AssertionError(label)
    checks+=1

# Positive dyadic rule, same source as the admitted packet.
xx,ww=leggauss(14)
rr=[];wwq=[]
for j in range(24):
    a=2.**(-j-1)
    rr.extend(1-(1.5*a+.5*a*xx));wwq.extend(.5*a*ww)
rr.append(1-2.**(-25));wwq.append(2.**(-24))
r=np.array(rr);w=np.array(wwq);s=np.sqrt(1-r*r)
beta=float(w@s); c=.25+beta*beta; lam=2*(1-c)
m2=float(w@(r*r));m4=float(w@(r**4))
m11=float(w@(r*s));m31=float(w@(r**3*s))
K=m2+2*beta*m11;J=m4+2*beta*m31
ck(np.min(w)>0 and np.min(r)>0,'positive quadrature')
ck(abs(sum(w)-1)<1e-14 and abs(w@r-.5)<1e-14,'exact mass and linear moment')
ck(abs(m2-1/3)<1e-14 and abs(m4-.2)<1e-14,'low moment finite accuracy')
ck(abs(beta-math.pi/4)<1e-10,'beta convergence')

# Nonlinear anchored smooth Hessian sandwich. q=1/4 gives Hessian base 1/2.
q=.25;eps=.1;k=1.;tau=math.exp(-.5*k*k)
f=lambda x:2*q*x+eps*k*(np.cos(k*x)-1)
u=lambda x:q*x*x+eps*(np.sin(k*x)-k*x)
fp=lambda x:2*q-eps*k*k*np.sin(k*x)
lower=2*q-abs(eps*k*k);upper=2*q+abs(eps*k*k)
ck(0<lower<upper<1,'global Hessian sandwich')
ck(f(0.)==0.,'anchored force')
delta=q*eps*tau*(3*lam*(k**5*J-3*k**3*m2)-6*k**3*K-(k**5-6*k**3))
delta_lim=q*eps*tau*(3*(1.5-math.pi**2/8)*(k**5*(.2+math.pi/15)-k**3)
                         -2*k**3*(1+math.pi/2)-(k**5-6*k**3))
ck(abs(delta-delta_lim)<1e-10,'finite and continuum third-cumulant defect')
ck(delta < -.009,'nonzero third-cumulant second-order defect')

z0,p=hermgauss(80);z0*=math.sqrt(2);p/=math.sqrt(math.pi)
z,gg=np.meshgrid(z0,z0,indexing='ij');wp=np.outer(p,p)
H=np.zeros_like(z);Hz=np.zeros_like(z)
for ri,si,wi in zip(r,s,w):
    arg=ri*z+si*gg
    H+=wi*f(arg);Hz+=wi*ri*fp(arg)
E=lambda v:float(np.sum(wp*v))
a1=-E(H);a2=lam*E(Hz*H);b1=-2*E(z*H)
first= -3*E(z*z*H)-3*a1
second=3*lam*E((z*z-1)*Hz*H)+3*E(z*H*H)-3*a1*b1
targetfirst=eps*k**3*tau
targetsecond=q*eps*tau*(k**5-6*k**3)
ck(abs(first-targetfirst)<1e-12,'third-cumulant first coefficient')
ck(abs(second-targetsecond-delta)<1e-11,'derived second coefficient')
ck(abs(E((z*z-1)*Hz*H)-q*eps*tau*(k**5*J-3*k**3*m2))<1e-12,'chord current Hermite identity')
ck(abs(E(z*H*H)-2*q*eps*k*(tau*(1-k*k*K)-1))<1e-12,'shared-root rank-two current')

rows=[]
for A in (.125,.0625,.03125,.015625,.0078125,.00390625):
    h=A*H
    h1=np.zeros_like(z)
    for ri,si,wi in zip(r,s,w): h1+=wi*A*f(ri*(z-h)+si*gg)
    Y=z-h-lam*(h1-h)
    mean=E(Y);var=E((Y-mean)**2);third=E((Y-mean)**3)
    weights=p*np.exp(-A*u(z0));weights/=sum(weights)
    tm=float(weights@z0);tv=float(weights@((z0-tm)**2));tt=float(weights@((z0-tm)**3))
    # Exact affine moment repair only makes the comparison MORE favorable.
    # It is not an available strong-mean/covariance service.
    repaired=tm+math.sqrt(tv/var)*(Y-mean)
    rt=E((repaired-tm)**3)
    ratio=(rt-tt)/(A*A)
    rows.append(dict(A=A,mean=mean,target_mean=tm,variance=var,target_variance=tv,
                     third_cumulant=third,target_third_cumulant=tt,
                     repaired_third_cumulant=rt,repaired_defect_over_A2=ratio))
    ck(abs(E(repaired)-tm)<1e-14,'ideal exact mean repair')
    ck(abs(E((repaired-tm)**2)-tv)<1e-14,'ideal exact variance repair')
    if A<=.015625: ck(abs(ratio-delta)<.001,'literal map asymptotic cumulant defect')
    # At reverse bridge r=0,t=1/2, every cumulant >=3 is multiplied by t^rank.
    t=.5
    ck(abs(t**3*(rt-tt)/(A*A)-ratio/8)<1e-15,'reverse bridge third-cumulant inheritance')

# Exact full-matrix quadratic calibration. The all-root covariance is retained.
rng=np.random.default_rng(630921)
for d in (1,2,5,11):
    for A in (.2,.1,.05,.025):
        Q,_=np.linalg.qr(rng.normal(size=(d,d)))
        B=Q@np.diag(A*rng.uniform(.1,1,size=d))@Q.T
        I=np.eye(d)
        # y=(I-B/2+lambda B^2/4)z+(-beta B+lambda beta B^2/2)g
        Cz=I-B/2+lam*(B@B)/4
        Cg=-beta*B+lam*beta*(B@B)/2
        actual=Cz@Cz.T+Cg@Cg.T
        degree2=I-B+B@B
        ck(np.linalg.norm(actual-degree2,'fro')<=3*A**3*math.sqrt(d),'matrix quadratic A^3 residual')
        ck(np.linalg.eigvalsh(actual).min()>.6,'positive matrix output covariance')

# First-order Gaussianization misses infinitely many odd conditional cumulants.
non_gauss=[]
for t in (.25,.5):
    char_first=-eps*math.exp(-1)*(math.sinh(t)-t)
    ck(char_first<0,'nonGaussian characteristic first coefficient')
    for rank in (3,5,7,9):
        coeff=-eps*math.sin(rank*math.pi/2)*tau*t**rank
        ck(abs(coeff)>0,'nonzero higher conditional cumulant')
    non_gauss.append(dict(t=t,sine_test_defect_over_A=char_first))

# Genuine C2/non-C3 two-dimensional gradient with noncommuting local Hessians.
# These Hessians diagnose first-action certificates; the map uses only f_values.
axis=np.array([1.,1.])/math.sqrt(2)
B0=np.diag([.45,.55]);rough=.05;smooth=.1
def rh(t):
    a=math.sqrt(abs(float(t)))
    return math.copysign(a*a-2*a+2*math.log1p(a),float(t)) if t else 0.
def f_values(x,alpha):
    return alpha*(B0@x+np.array([smooth*(math.cos(x[0])-1),0.])+rough*rh(axis@x)*axis)
def f_first(x,alpha):
    a=math.sqrt(abs(float(axis@x)))
    return alpha*(B0+np.diag([-smooth*math.sin(x[0]),0.])+rough*a/(1+a)*np.outer(axis,axis))
Rs=[np.concatenate((ri*np.eye(2),si*np.eye(2)),axis=1) for ri,si in zip(r,s)]
P=np.concatenate((np.eye(2),np.zeros((2,2))),axis=1);Pi=P.T@P
LQ=float(sum(w/r))
def full_lift(W,alpha):
    val=np.zeros(4);der=np.zeros((4,4));K=np.zeros(2);dK=np.zeros((2,4))
    for ri,wi,Ri in zip(r,w,Rs):
        x=Ri@W;fv=f_values(x,alpha);fd=f_first(x,alpha)
        val+=wi/ri*(Ri.T@fv);der+=wi/ri*(Ri.T@fd@Ri)
        K+=wi*fv;dK+=wi*(fd@Ri)
    return val,der,K,dK
ck(np.linalg.norm(f_first(np.array([1.,0.]),1)@f_first(np.array([0.,1.]),1)
                  -f_first(np.array([0.,1.]),1)@f_first(np.array([1.,0.]),1))>1e-4,
   'noncommuting C2 fixture Hessians')
lift_maxcurl=0.
for alpha in (.2,.1,.03):
    for ii in range(12):
        W=rng.normal(size=4)
        V0,D0,K0,DK0=full_lift(W,alpha)
        V1,D1,K1,DK1=full_lift(W-P.T@K0,alpha)
        EJ=V1-V0;DE=D1-D0-D1@Pi@D0
        EK=K1-K0;DEK=DK1@(np.eye(4)-P.T@DK0)-DK0
        curl=DE-DE.T;bound=2*alpha*alpha
        lift_maxcurl=max(lift_maxcurl,float(np.linalg.norm(curl,2)/bound))
        ck(np.linalg.norm(P@V0-K0)<1e-13,'full lift projected VALUE')
        ck(np.linalg.norm(P@D0-DK0)<1e-13,'full lift projected first')
        ck(np.linalg.norm(P@EJ-EK)<1e-13,'actual chord lift VALUE')
        ck(np.linalg.norm(P@DE-DEK)<1e-13,'actual chord lift first')
        ck(np.linalg.norm(D0-D0.T)<1e-13,'full gradient symmetry')
        ck(np.linalg.norm(curl,2)<=bound*(1+1e-12),'full lift small curl')
        ck(np.linalg.norm(EK)<=alpha/2*np.linalg.norm(K0)*(1+1e-10),'one-energy chord VALUE bound')
        ck(np.linalg.norm(DE,2)<=2*alpha*LQ+alpha*alpha,'literal full first bound')
        ck(np.linalg.eigvalsh(f_first(P@W,alpha)).min()>=.35*alpha-1e-14,'C2 Hessian lower')
        ck(np.linalg.eigvalsh(f_first(P@W,alpha)).max()<=.65*alpha+1e-14,'C2 Hessian upper')
    V0,D0,K0,DK0=full_lift(np.zeros(4),alpha)
    ck(np.array_equal(V0,np.zeros(4)) and np.array_equal(K0,np.zeros(2)),'literal anchor zero')

report=dict(status='PASS',assertions=checks,nodes=len(r),beta=beta,c=c,lambda_chord=lam,
            m2=m2,m4=m4,m11=m11,m31=m31,K=K,J=J,
            target_third_first=targetfirst,target_third_second=targetsecond,
            chord_third_second=second,third_cumulant_defect_coefficient=delta,
            literal_program=rows,gaussianization=non_gauss,
            full_lift_LQ=LQ,max_C2_lift_curl_fraction=lift_maxcurl,
            scope='Positive finite original-VALUE local-chord test and precise nonlinear gate. No all-order impossibility or new full-law grade is claimed.')
(OUT/'same_root_chord_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='literal_program'},indent=2))
