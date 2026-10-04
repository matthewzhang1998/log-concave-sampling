#!/usr/bin/env python3
"""Literal finite zero-clock VALUE square action: algebra, tape and query audit.

No test asserts a universal high-order consumer bound from finite sampling.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
SERVICE=ROOT.parent/'POSITIVE-SQUARE-CLOCK-COVARIANCE-SERVICE.md'
rng=np.random.default_rng(620041005)
checks=0

def ck(ok,msg):
    global checks
    assert bool(ok),msg
    checks+=1

def rule(K,m):
    x,w=leggauss(m);r=[];v=[]
    for j in range(K):
        a=2.**(-j-1)
        r.extend(1-1.5*a-.5*a*x);v.extend(.5*a*w)
    r.append(1-2.**(-K-1));v.append(2.**(-K))
    return np.array(r),np.array(v)

def first_filter(K):
    a,b=1/(4*K*K),.25
    t=(a+b)/2+(b-a)*np.cos(np.arange(K)*math.pi/(K-1))/2
    d=np.array([np.prod([-t[j]/(t[i]-t[j]) for j in range(K) if j!=i]) for i in range(K)])
    ck(abs(d.sum()-1)<2e-13,'first coefficient constant exactness')
    ck(np.sum(abs(d))<=4+1e-13,'bounded extrapolation absolute weights')
    for k in range(1,K):
        ck(abs(d@(t**k))<2e-13,'first filter removes prescribed even powers')
    return np.sqrt(t),d

class RawSource:
    def __init__(self,g,X,r,w,ell):
        self.g,self.X,self.r,self.w,self.ell=g,X,r,w,ell
        self.c=np.sqrt((1-r)*(1+r))
        self.anchor=g(r[:,None]*X)
        self.g_calls=len(r)
        self.f_calls=0
    def value(self,u):
        self.f_calls+=1
        self.g_calls+=len(self.r)
        x=self.r[:,None]*self.X+self.c[:,None]*u
        return np.sum((self.w*self.r/self.c)[:,None]*(self.g(x)-self.anchor),axis=0)

def action(source,p,bank,b,d,padding):
    """Exactly C_sq(s0 convention)/s0^2, with the fixed s0 canceled algebraically."""
    z1,z0,z2=bank
    ell=source.ell
    def normalized(u):
        return source.value(u)/ell
    I=np.zeros_like(p)
    for bj,dj in zip(b,d):
        center=math.sqrt(1-bj*bj)*z1
        I+=dj*(normalized(center+bj*p)-normalized(center-bj*p))/(2*bj)
    beta=.5
    def R(w):
        center=math.sqrt(1-beta*beta)*z2
        return (normalized(center+beta*w)-normalized(center-beta*w))/(2*beta)
    return ell*ell*(R(z0+padding*I)-R(z0-padding*I))/(2*padding)

rows=[]
for panel,order in [(1,1),(2,2),(5,5),(12,12),(20,20)]:
    r,w=rule(panel,order);c=np.sqrt((1-r)*(1+r))
    cx=float(np.sum(w*r*r/c));cv=float(np.sum(w*r/c));cb=float(np.sum(w/c))
    ck(cx<=1+math.sqrt(2)+2e-13 and cv<=1+math.sqrt(2)+2e-13 and cb<=1+math.sqrt(2)+2e-13,
       'uniform source caller and absolute VALUE coefficient bound')
    rows.append({'panels':panel,'order':order,'nodes':len(r),'caller_coefficient':cx,'value_coefficient':cv})

r,w=rule(3,3);q,v=rule(2,2);aj=2*v*q
ck(abs(aj.sum()-1)<2e-14,'positive outer action mass exactly one')
results=[]
for dim in [1,3,7]:
  for Kf in [2,4,7]:
    b,d=first_filter(Kf)
    A=.12;ell=A/2
    B0=rng.normal(size=(dim,dim));B=(A/np.linalg.norm(B0,2)**2)*(B0@B0.T)
    Z=rng.normal(size=dim);p=rng.normal(size=dim)
    G=rng.normal(size=(len(q),dim));banks=rng.normal(size=(len(q),3,dim))
    def linear(x): return x@B.T
    def nonlinear(x):
        # Smooth gradient with Hessian A diag(.5-.25 sin(x_i)); anchored.
        return A*(.5*x+.25*(np.cos(x)-1))
    for name,g in [('quadratic',linear),('nonlinear',nonlinear),('zero',lambda x:np.zeros_like(x))]:
        outs=[];g_calls=0;f_calls=0
        for j in range(len(q)):
            X=q[j]*Z+math.sqrt(1-q[j]*q[j])*G[j]
            source=RawSource(g,X,r,w,ell)
            out=action(source,p,banks[j],b,d,math.sqrt(A)/8)
            ck(source.f_calls==2*Kf+4,'literal zero-clock raw source occurrence count')
            ck(source.g_calls==len(r)*(2*Kf+5),'caller anchors counted once per complete action bank')
            reverse=RawSource(g,X,r,w,ell)
            neg=action(reverse,-p,banks[j],b,d,math.sqrt(A)/8)
            ck(np.linalg.norm(out+neg)<2e-13,'pointwise incoming-p oddness retains all identical private roots')
            zero=RawSource(g,X,r,w,ell)
            oz=action(zero,np.zeros_like(p),banks[j],b,d,math.sqrt(A)/8)
            ck(np.linalg.norm(oz)==0,'zero incoming p annihilates the entire action at arbitrary caller and owned roots')
            outs.append(out);g_calls+=source.g_calls;f_calls+=source.f_calls
        agg=aj@np.array(outs)
        ck(g_calls==len(r)*len(q)*(2*Kf+5),'full aggregated anchor-plus-shifted original VALUE count')
        ck(f_calls==len(q)*(2*Kf+4),'complete raw occurrence count across outer clocks')
        if name=='quadratic':
            ck(np.linalg.norm(agg-(B@B@p)/4)<3e-13,'shared incoming p gives full matrix square exactly')
            eta,zeta=.5,math.sqrt(.75)
            C=B@B/4
            cov=(eta*np.eye(dim)+C/(2*eta))@(eta*np.eye(dim)+C/(2*eta)).T+zeta*zeta*np.eye(dim)
            target=np.eye(dim)+B@B/4+np.linalg.matrix_power(B,4)/(64*eta*eta)
            ck(np.linalg.norm(cov-target)<3e-13,'linearized reserve plus paid fourth-degree covariance term')
        if name=='zero':
            ck(np.linalg.norm(agg)==0,'source-zero reserve carrier is exactly eta*p+zeta*z')
        results.append({'D':dim,'filter_nodes':Kf,'source':name,'g_values':g_calls,
                        'raw_f_occurrences':f_calls,'basic_scalar_roots':4*len(q)+2})

report={'status':'PASS','assertions':checks,'seed':620041005,
        'service_sha256':hashlib.sha256(SERVICE.read_bytes()).hexdigest(),
        'scope':'Actual finite VALUE-only source/query/anchor/oddness/carrier and matrix-quadratic checks. Imported calibration and Gaussian-buffer inequalities are justified analytically in the audit, not inferred from these finite checks.',
        'uniform_rule_coefficients':rows,'action_cases':results,
        'count_formula':'N_r*N_q*(2*K_filter+5) with one exact same-X anchor bank; no-cache upper bound 2*N_r*N_q*(2*K_filter+4).',
        'basic_tape_formula':'D*(4*N_q+2), before external host and added numerical-restoration tapes.'}
(ROOT/'square_clock_value_action_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='action_cases'},indent=2))
