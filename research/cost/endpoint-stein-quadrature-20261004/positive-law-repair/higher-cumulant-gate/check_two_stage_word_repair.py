#!/usr/bin/env python3
"""Positive two-stage VALUE repair of the first exposed nonlinear cumulant word."""
from pathlib import Path
import json, math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
OUT=Path(__file__).resolve().parent
checks=0
def ck(ok,msg):
    global checks
    if not ok: raise AssertionError(msg)
    checks+=1
lg,lw=leggauss(12);rr=[];ww=[]
for j in range(20):
    h=2.**(-j-1);rr.extend(1-(1.5*h+.5*h*lg));ww.extend(.5*h*lw)
rr.append(1-2.**(-21));ww.append(2.**(-20))
r=np.array(rr);w=np.array(ww);s=np.sqrt(1-r*r)
beta=float(w@s);C=.25+beta**2;m2=float(w@(r*r));m4=float(w@(r**4))
m11=float(w@(r*s));m31=float(w@(r**3*s));m22=m2-m4
K=m2+2*beta*m11
Ju=m4+2*beta*m31;Jv=m31+2*beta*m22;J0=2*m4
Lu=3*m2;Lv=2*(m11+beta*m2);L0=1+4*m2
mat=np.array([[.5,beta,1.],[Ju,Jv,J0],[Lu,Lv,L0]])
rhs=np.array([1-C,1/3,2-2*K])
a,c,b=np.linalg.solve(mat,rhs)
ck(abs(np.linalg.det(mat))>.12,'known coefficient matrix gap')
ck(max(abs(a),abs(b),abs(c))<1.5,'fixed known coefficients')
ck(np.max(np.abs(mat@np.array([a,c,b])-rhs))<1e-15,'all three exact word equations')
z0,p=hermgauss(85);z0*=math.sqrt(2);p/=math.sqrt(math.pi)
z,v=np.meshgrid(z0,z0,indexing='ij');wp=np.outer(p,p)
E=lambda x:float(np.sum(wp*x))

rows=[]
for q in (.2,.3):
    for k in (.5,1.,1.5,2.):
        eps=.08/(1+k*k);tau=math.exp(-k*k/2)
        f=lambda x:2*q*x+eps*k*(np.cos(k*x)-1)
        U=lambda x:q*x*x+eps*(np.sin(k*x)-k*x)
        H=np.zeros_like(z);Hu=np.zeros_like(z);Hv=np.zeros_like(z)
        for ri,si,wi in zip(r,s,w):
            x=ri*z+si*v
            H+=wi*f(x)
            df=2*q-eps*k*k*np.sin(k*x)
            Hu+=wi*ri*df;Hv+=wi*si*df
        T=a*Hu*H+c*Hv*H+b*Hu*f(z)
        a1=-E(H);b1=-2*E(z*H)
        coeff2=3*E((z*z-1)*T)+3*E(z*H*H)-3*a1*b1
        exact2=q*eps*tau*(k**5-6*k**3)
        ck(abs(coeff2-exact2)<2e-12,'complete nonlinear third-cumulant A2 coefficient')
        for A in (.0625,.03125,.015625,.0078125):
            Hstar=np.zeros_like(z)
            zu=z-A*(a*H+b*f(z));vv=v-A*c*H
            for ri,si,wi in zip(r,s,w): Hstar+=wi*A*f(ri*zu+si*vv)
            Y=z-Hstar
            my=E(Y);vy=E((Y-my)**2);ty=E((Y-my)**3)
            weights=p*np.exp(-A*U(z0));weights/=sum(weights)
            mx=float(weights@z0);vx=float(weights@((z0-mx)**2));tx=float(weights@((z0-mx)**3))
            defect=ty-tx
            rows.append(dict(q=q,k=k,A=A,third_defect=defect,defect_over_A2=defect/A**2,
                             defect_over_A3=defect/A**3))
            ck(abs(defect/A**3)<.4,'literal positive predictor third-cumulant cubic residual')
            ck(vy>.85 and vx>.85,'actual scalar variance gap')

# A non-single-frequency odd perturbation uses exactly the same coefficients.
q=.25;eps=.03
f=lambda x:2*q*x+eps*((np.cos(x)-1)+.6*(np.cos(2*x)-1))
U=lambda x:q*x*x+eps*((np.sin(x)-x)+.3*(np.sin(2*x)-2*x))
H=sum((wi*f(ri*z+si*v) for ri,si,wi in zip(r,s,w)),np.zeros_like(z))
for A in (.0625,.03125,.015625,.0078125):
    zu=z-A*(a*H+b*f(z));vv=v-A*c*H
    Y=z-sum((wi*A*f(ri*zu+si*vv) for ri,si,wi in zip(r,s,w)),np.zeros_like(z))
    my=E(Y);ty=E((Y-my)**3)
    ww=p*np.exp(-A*U(z0));ww/=sum(ww);mx=float(ww@z0);tx=float(ww@((z0-mx)**3))
    ck(abs((ty-tx)/A**3)<.4,'multi-frequency odd perturbation')

# Exact non-diagonal matrix VALUE expansion, no Hessian producer.
rng=np.random.default_rng(118034)
for d in (1,2,5,13):
    for alpha in (.2,.1,.05):
        O,_=np.linalg.qr(rng.normal(size=(d,d)));B=O@np.diag(alpha*rng.uniform(.1,1,d))@O.T
        eff=a/2+c*beta
        Cz=np.eye(d)-B/2+(eff+b)*(B@B)/2
        Cv=-beta*B+eff*beta*(B@B)
        cov=Cz@Cz.T+Cv@Cv.T
        ck(np.linalg.norm(cov-(np.eye(d)-B+B@B),'fro')<4*alpha**3*math.sqrt(d),'full matrix quadratic word')
        ck(np.linalg.eigvalsh(cov).min()>.65,'positive actual matrix covariance')
        Z=rng.normal(size=d);V=rng.normal(size=d)
        K0=sum((wi*(B@(ri*Z+si*V)) for ri,si,wi in zip(r,s,w)),np.zeros(d))
        Zu=Z-a*K0-b*(B@Z);Vv=V-c*K0
        Y=Z-sum((wi*(B@(ri*Zu+si*Vv)) for ri,si,wi in zip(r,s,w)),np.zeros(d))
        ck(np.linalg.norm(Y-(Cz@Z+Cv@V))<1e-12,'literal two-stage matrix VALUE map')

# Fixed-u private-v source small curl for a genuine noncommuting source.
B0=np.diag([.45,.55]);axis=np.array([1.,1.])/math.sqrt(2)
def F(x,alpha):
    return alpha*(B0@x+.06*(math.cos(axis@x)-1)*axis+.04*np.array([math.cos(x[0])-1,0.]))
def DF(x,alpha):
    return alpha*(B0-.06*math.sin(axis@x)*np.outer(axis,axis)-.04*np.diag([math.sin(x[0]),0.]))
def packet(U,V,alpha):
    out=np.zeros(2);Du=np.zeros((2,2));Dv=np.zeros((2,2))
    for ri,si,wi in zip(r,s,w):
        x=ri*U+si*V;out+=wi*F(x,alpha);D=DF(x,alpha);Du+=wi*ri*D;Dv+=wi*si*D
    return out,Du,Dv
for alpha in (.15,.05):
    for j in range(16):
        U=rng.normal(size=2);V=rng.normal(size=2)
        K0,Ku0,Kv0=packet(U,V,alpha)
        Ks,Kus,Kvs=packet(U-a*K0-b*F(U,alpha),V-c*K0,alpha)
        Dv=Kvs-(a*Kus+c*Kvs)@Kv0
        expected=-(a*Kus+c*Kvs)@Kv0+Kv0@(a*Kus+c*Kvs)
        ck(np.linalg.norm(Dv-Dv.T-expected)<1e-14,'captured-u private-v curl identity')
        ck(np.linalg.norm(Dv-Dv.T,2)<=2*(abs(a)/2+abs(c)*beta)*beta*alpha**2,'conditional small-curl bound')
        ck(np.linalg.norm(Ks-K0)<=alpha*(math.hypot(a,c)*np.linalg.norm(K0)+abs(b)*np.linalg.norm(F(U,alpha)))+1e-13,'literal chord energy')
    K0,_,_=packet(np.zeros(2),np.zeros(2),alpha)
    ck(np.array_equal(K0,np.zeros(2)) and np.array_equal(F(np.zeros(2),alpha),np.zeros(2)),'original exact zero')

report=dict(status='PASS',assertions=checks,nodes=len(r),beta=beta,
            coefficients=dict(a=float(a),c=float(c),b=float(b)),determinant=float(np.linalg.det(mat)),
            word_matrix=mat.tolist(),word_rhs=rhs.tolist(),literal_cumulant_checks=rows,
            scope='Positive two-stage original-VALUE repair matches full matrix quadratic covariance and the scalar quadratic-times-odd third-cumulant word through A^2. No full-law order-three claim.')
(OUT/'two_stage_word_repair_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='literal_cumulant_checks'},indent=2))
