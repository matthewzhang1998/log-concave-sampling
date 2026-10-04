#!/usr/bin/env python3
"""Independent scalar checks of the conditional OU covariance square identity."""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
ROOT=Path(__file__).resolve().parent
checks=0

def check(ok,msg):
    global checks
    assert bool(ok),msg
    checks+=1

def unit_rule(n):
    x,w=leggauss(n)
    return (1+x)/2,w/2

def dyadic_rule(K,m):
    x,w=leggauss(m);rs=[];ws=[]
    for j in range(K):
        a=2.**(-j-1)
        t=1.5*a+.5*a*x
        rs.extend(1-t);ws.extend(.5*a*w)
    rs.append(1-2.**(-K-1));ws.append(2.**(-K))
    return np.array(rs),np.array(ws)

def direct_covariance(z,k,n):
    r,w=unit_rule(n);tau,v=unit_rule(n)
    r=r[:,None];tau=tau[None,:]
    sx=r*tau
    vx=1-r*r;vy=1-sx*sx;cross=tau*(1-r*r)
    mx=r*z;my=sx*z;u=k*k
    ex=np.exp(-u*vx/2);ey=np.exp(-u*vy/2)
    cc=.5*(np.exp(-u*(vx+vy+2*cross)/2)*np.cos(k*(mx+my))+
            np.exp(-u*(vx+vy-2*cross)/2)*np.cos(k*(mx-my)))-ex*ey*np.cos(k*mx)*np.cos(k*my)
    a,b=.5,.25
    cov=a*a*cross-a*b*cross*(ex*np.sin(k*mx)+ey*np.sin(k*my))+b*b/u*cc
    return float(np.sum(2*r*w[:,None]*v[None,:]*cov))

def square_covariance(z,k,s,ws,q,wq,ghn):
    h,wh=hermgauss(ghn);h=math.sqrt(2)*h;wh=wh/math.sqrt(math.pi)
    y=q[:,None]*z+np.sqrt(1-q*q)[:,None]*h[None,:]
    hs=.5-.25*np.exp(-.5*k*k*(1-s*s))[:,None,None]*np.sin(k*s[:,None,None]*y[None,:,:])
    B=np.sum((ws*s)[:,None,None]*hs,axis=0)
    check(np.all(B>=.125-2e-14) and np.all(B<=.375+2e-14),'positive bounded derivative mean')
    C=float(np.sum((2*wq*q)[:,None]*wh[None,:]*B*B))
    check(0<=C<=.75**2/4+2e-14,'positive covariance square with Hessian ceiling .75')
    return C

cases=[]
for k in [1.,3.,5.]:
    for z in [-2.,0.,1.5]:
        d1=direct_covariance(z,k,72);d2=direct_covariance(z,k,100)
        s,w=unit_rule(80)
        c1=square_covariance(z,k,s,w,s,w,80)
        s,w=unit_rule(110)
        c2=square_covariance(z,k,s,w,s,w,110)
        check(abs(d1-d2)<1e-12,'direct covariance integration resolution')
        check(abs(c1-c2)<1e-12,'square integration resolution')
        check(abs(d2-c2)<2e-12,'direct conditional pair covariance equals positive square identity')
        rs,ws=dyadic_rule(12,12)
        cq=square_covariance(z,k,rs,ws,rs,ws,110)
        check(abs(cq-d2)<3e-10,'finite positive weighted Mehler square approximates covariance')
        cases.append({'k':k,'endpoint':z,'direct_covariance':d2,'square_covariance':c2,
                      'identity_error':c2-d2,'positive_square_rule_error':cq-d2})
# Algebraic low-chaos checks do not rely on the bounded sinusoidal fixture.
for z in [-2.,0.,1.5]:
    q,w=unit_rule(8)
    # g=x^2: B(z)=2z/3; P_q B^2=4(q^2 z^2+1-q^2)/9.
    value=float(np.sum(2*w*q*4*(q*q*z*z+1-q*q)/9))
    check(abs(value-2*(1+z*z)/9)<3e-15,'degree-two polynomial exact normalization')
report={'status':'PASS','assertions':checks,'scope':'Analytical identity and positive finite square quadrature; P and Hessian evaluations are diagnostics, not authorized producer leaves.',
        'cases':cases,'max_identity_error':max(abs(c['identity_error']) for c in cases),
        'max_positive_square_rule_error':max(abs(c['positive_square_rule_error']) for c in cases)}
(ROOT/'conditional_covariance_square_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
