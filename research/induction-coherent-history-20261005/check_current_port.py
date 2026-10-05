#!/usr/bin/env python3
"""Deterministic algebra/quadrature diagnostics; not a proof or native execution."""
import json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
OUT=Path(__file__).parent
checks=[]
def check(name, value, bound=1e-9):
    checks.append({'name':name,'value':float(value),'bound':float(bound),'passed':bool(value<=bound)})
def he(n,x):
    if n==0:return np.ones_like(x)
    if n==1:return x
    a=np.ones_like(x);b=x
    for k in range(1,n):a,b=b,x*b-k*a
    return b
x,w=hermegauss(52);w=w/math.sqrt(2*math.pi)
G=np.stack(np.meshgrid(x,x,indexing='ij'),axis=-1).reshape(-1,2)
W=(w[:,None]*w[None,:]).ravel()
for aa,bb in [(0.11,-0.08),(0.125,0.125),(0.,0.),(-0.12,0.04)]:
    B=np.diag([aa,bb]);K=np.eye(2)+B;E=K@K-np.eye(2)
    Z=G@K.T
    L=np.exp(-np.linalg.slogdet(K)[1]+(np.sum(G*G,axis=1)-np.sum((G@np.linalg.inv(K).T)**2,axis=1))/2)
    lam=np.diag(E)
    check(f'likelihood_normalization_{aa}_{bb}',abs(W@L-1))
    for p in [2,3]:
        exact=np.linalg.det(np.eye(2)+E)**(-(p-1)/2)*np.linalg.det(np.eye(2)-(p-1)*E)**(-.5)
        check(f'likelihood_moment_p{p}_{aa}_{bb}',abs(W@(L**p)-exact),2e-9)
    phi=lambda a: np.stack([a[:,0]**2+0.3*a[:,1]**4,a[:,0]*a[:,1]+a[:,0]**6],axis=1)
    check(f'exact_marked_current_polynomial_{aa}_{bb}',np.max(abs(np.sum(W[:,None]*phi(Z),axis=0)-np.sum((W*L)[:,None]*phi(G),axis=0))),3e-9)
    for p in range(5):
        P=np.zeros(len(G))
        for j in range(p+1):
            for a in range(j+1):
                b=j-a
                P+=lam[0]**a*lam[1]**b*he(2*a,G[:,0])*he(2*b,G[:,1])/(2**j*math.factorial(a)*math.factorial(b))
        tail=float(np.sqrt(W@((L-P)**2)))
        e=max(abs(lam));bound=math.sqrt(2)*(math.sqrt(2)*e)**(p+1)
        check(f'hermite_tail_p{p}_{aa}_{bb}',tail,bound+1e-12)
# Correlated source states: rank-two transport, a nonconstant same-bank mark.
Fs=[np.array([.13,.08]),np.array([-.11,.14]),np.array([.07,-.15])]
Fps=[np.array([.09,-.1]),np.array([.15,.04]),np.array([-.02,-.1])]
probs=np.array([.2,.3,.5]);v=.5
Hs=[(np.outer(a,b)+np.outer(b,a))/2 for a,b in zip(Fs,Fps)]
Cs=sum(p*H for p,H in zip(probs,Hs))
# PSD replica example for genuine Gaussian-reference bound.
HsPSD=[np.outer(a,a) for a in Fs];C=sum(p*H for p,H in zip(probs,HsPSD))
ev,U=np.linalg.eigh(v*np.eye(2)+C);S=(U*np.sqrt(ev))@U.T
lhs=math.sqrt(sum(p*np.linalg.norm(math.sqrt(v)*(np.eye(2)+H/(2*v))-S,'fro')**2 for p,H in zip(probs,HsPSD)))
rhs=math.sqrt(sum(p*np.linalg.norm(H-C,'fro')**2 for p,H in zip(probs,HsPSD)))/(2*math.sqrt(v))+np.linalg.norm(C@C,'fro')/(8*v**1.5)
check('retained_mark_joint_W2_majorant',lhs,rhs+1e-12)
# Explicit source derivative and finite full-bank curl test.
def raw(om):
    a=.08*np.array([math.sin(om[0]),math.cos(om[1])]);b=.07*np.array([math.cos(om[1]),math.sin(om[0]+om[1])])
    m=.02*np.array([math.sin(om[0])+om[1],math.cos(om[1])])
    H=(np.outer(a,b)+np.outer(b,a))/2
    E=H/v+H@H/(4*v*v)
    return a,b,m,E
om=np.array([.23,-.41]);eps=1e-6;a,b,m,E=raw(om);gg=np.array([.8,-.7]);e=np.linalg.norm(E,2)
DF=[];DFp=[];DM=[];DE=[]
for j in range(2):
    d=np.eye(2)[j]*eps;ap,bp,mp,Ep=raw(om+d);am,bm,mm,Em=raw(om-d)
    DF.append((ap-am)/(2*eps));DFp.append((bp-bm)/(2*eps));DM.append((mp-mm)/(2*eps));DE.append((Ep-Em)/(2*eps))
LF=max(np.linalg.norm(np.array(DF).T,2),np.linalg.norm(np.array(DFp).T,2));LM=np.linalg.norm(np.array(DM).T,2);R=max(np.linalg.norm(a),np.linalg.norm(b));RM=np.linalg.norm(m);BG=np.linalg.norm(gg)
LH=2*R*LF;h=R*R/(2*v);LE=(1+h)*LH/v
bound=LM*e*(BG*BG+2)/2+RM*LE*(BG*BG+4)/2+RM*e*BG
J=[]
def nfun(ww):
    *_,mm,EE=raw(ww[:2]);q=(ww[2:]@EE@ww[2:]-np.trace(EE))/2;return mm*q
full=np.r_[om,gg]
for j in range(4):
    d=np.eye(4)[j]*eps;J.append((nfun(full+d)-nfun(full-d))/(2*eps))
J=np.array(J).T
check('first_correction_complete_first_pointwise',np.linalg.norm(J,2),bound+1e-9)
# Expected native exponents: exponents are (A,v,u).
EN=(F(7),F(-1),F(0));LN=(F(37,9),F(-1),F(0));r=tuple(a/2 for a in (LN[0],LN[1],F(-1,2)))
err=tuple(a+b for a,b in zip(EN,(LN[0],LN[1],F(-1,2))))
assert err==(F(100,9),F(-2),F(-1,2))
r_at_v_u_A2=r[0]+2*r[1]+2*r[2];err_at_v_u_A2=err[0]+2*err[1]+2*err[2]
assert r_at_v_u_A2==F(5,9) and err_at_v_u_A2==F(55,9)
# Native five-row substitution exactly: relative exponents in r after E_N.
assert [F(2),F(2),F(2),F(5,2),F(3)]==[F(2),F(2),F(2),F(3)-F(1,2),F(3)]
from current_source import exact_likelihood, first_current_polynomial, ReplicaRows, LocalMarkRows, nineteen_value_source
rng=np.random.default_rng(905)
for dim in [1,2,7,40]:
    ff=.05*rng.normal(size=dim)/math.sqrt(dim); ffp=.05*rng.normal(size=dim)/math.sqrt(dim); gg=rng.normal(size=dim)
    hh=(np.outer(ff,ffp)+np.outer(ffp,ff))/2; kk=np.eye(dim)+hh/(2*v); ee=kk@kk-np.eye(dim)
    exact=math.exp(-np.linalg.slogdet(kk)[1]+(gg@gg-(np.linalg.solve(kk,gg)@np.linalg.solve(kk,gg)))/2)
    check(f'woodbury_likelihood_dim{dim}',abs(exact_likelihood(ff,ffp,gg,v)-exact))
    check(f'matrix_free_current_dim{dim}',abs(first_current_polynomial(ff,ffp,gg,v)-(gg@ee@gg-np.trace(ee))/2))
count=[0]
def grad(xx):
    count[0]+=1
    return .03*xx+.01*np.sin(xx)
r0=ReplicaRows(*(rng.normal(size=3) for _ in range(6)))
r1=ReplicaRows(*(rng.normal(size=3) for _ in range(6)))
mr=LocalMarkRows(*(rng.normal(size=3) for _ in range(3)))
out=nineteen_value_source(grad,r0,r1,mr,rng.normal(size=3),epsilon=.2,delta=.1,scale=.4,pair_cap=.1,mark_cap=.02,carrier_cap=5.,v=.5)
assert count[0]==19 and np.all(np.isfinite(out))
checks.append({'name':'literal_original_value_count','value':count[0],'bound':19,'passed':count[0]==19})
result={'status':'passed' if all(c['passed'] for c in checks) else 'FAILED','checks':checks,'native_exponents':{'r_at_u_v_A2':str(r_at_v_u_A2),'error_at_u_v_A2':str(err_at_v_u_A2),'retuned_source_first_A':str(F(-8,9)),'current_first_A':str(F(37,9))},'scope':'Algebra, finite Gaussian quadrature and local derivative diagnostics. No native compiler or true history producer executed.'}
(OUT/'current_port_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'checks':len(checks),'failures':[c for c in checks if not c['passed']],'native_exponents':result['native_exponents']},indent=2))
assert result['status']=='passed'
