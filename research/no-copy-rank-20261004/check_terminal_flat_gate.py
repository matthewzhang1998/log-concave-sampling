#!/usr/bin/env python3
"""Finite algebra/first checks for the terminal-flat SAME-potential K3 gate."""
import hashlib, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
T=np.array([[.5,.05],[.05,.5]])
N=np.array([[1.,0.],[0.,0.]])
e1=np.array([1.,0.])
delta=1/40
c,s=3/5,4/5
h=.1
count=0

def close(x,y,tol=3e-11):
    global count
    assert np.max(np.abs(np.asarray(x)-np.asarray(y))) <= tol, (x,y)
    count+=1

def psi(z):
    u=(z+.5)/h
    if abs(u)>=1: return 0.
    return (z+.5)*math.exp(1-1/(1-u*u))

def dpsi(z):
    u=(z+.5)/h
    if abs(u)>=1: return 0.
    inv=1/(1-u*u)
    return math.exp(1-inv)*(1-2*u*u*inv*inv)

def ddpsi(z):
    u=(z+.5)/h
    if abs(u)>=1: return 0.
    inv=1/(1-u*u)
    chi=math.exp(1-inv)
    cp=-2*u*inv*inv*chi
    cpp=(-2*inv*inv-8*u*u*inv**3+4*u*u*inv**4)*chi
    return (2*cp+u*cpp)/h

for u in np.linspace(-.9999,.9999,20001):
    close(dpsi(-.5+h*u),math.exp(1-1/(1-u*u))*(1-2*u*u/(1-u*u)**2),1e-10)
    assert abs(dpsi(-.5+h*u)) <= 2
    count+=1
max_first_error=0.
cases=[]
for a in [1/16,1/32,1/64,1/128,1/512,1/2048]:
    r=a; eps=a**.9; q=a*eps
    def grad(y): return T@y+delta*q*psi(y[0]/q)*e1
    def Hess(y): return T+delta*dpsi(y[0]/q)*N
    S=np.zeros(2); U=np.zeros(2); Z=e1
    x=c*S+s*U
    Delta=grad(S)-grad(S+eps*Z)
    x1=x+a*Delta
    t0=S+a*grad(x)
    t1=S+a*grad(x1)
    E=r*(grad(t0)-grad(t1))
    close(Delta,-eps*T@e1)
    close(x1,-q*T@e1)
    close(t1,-a*q*T@T@e1)
    close(Hess(t0),T)
    close(Hess(t1),T)
    close(Hess(x),T)
    close(Hess(x1),T+delta*N)
    close(Hess(S)-Hess(S+eps*Z),np.zeros((2,2)))
    K0=Hess(t0)@Hess(x)
    K1=Hess(t1)@Hess(x1)
    close(K0-K1,-delta*T@N)
    close((K0-K1)-(K0-K1).T,-delta*(T@N-N@T))
    close(E,r*a*a*eps*T@T@T@e1)
    DS=r*(Hess(t0)-Hess(t1))+r*a*c*(K0-K1)
    DU=r*a*s*(K0-K1)
    DZ=r*a*a*eps*K1@Hess(S+eps*Z)
    close(DS,-r*a*c*delta*T@N)
    close(DU,-r*a*s*delta*T@N)
    close(DZ,r*a*a*eps*T@(T+delta*N)@T)
    common=r*a/math.sqrt(2)*(Hess(t0)-Hess(t1))
    opposite=r*a/math.sqrt(2)*(Hess(t0)+Hess(t1))
    dp_plus=(Hess(x)+Hess(x1))/math.sqrt(2)
    dp_minus=(Hess(x)-Hess(x1))/math.sqrt(2)
    close(common,np.zeros((2,2)))
    close(common@dp_plus+opposite@dp_minus,r*a*(K0-K1))
    # Finite differences verify the entire actual original-source Jacobian.
    def source(w):
        SS,UU,ZZ=w[:2],w[2:4],w[4:]
        xx=c*SS+s*UU
        dd=grad(SS)-grad(SS+eps*ZZ)
        xx1=xx+a*dd
        return r*(grad(SS+a*grad(xx))-grad(SS+a*grad(xx1)))
    w=np.r_[S,U,Z]; dh=q*1e-4
    jnum=np.column_stack([(source(w+dh*np.eye(6)[i])-source(w-dh*np.eye(6)[i]))/(2*dh) for i in range(6)])
    jexact=np.concatenate([DS,DU,DZ],axis=1)
    error=np.max(np.abs(jnum-jexact))/(r*a)
    assert error<3e-7,error
    count+=1
    max_first_error=max(max_first_error,float(error))
    for z in np.linspace(-.7,.1,1001):
        y=np.array([z*q,0.])
        ev=np.linalg.eigvalsh(Hess(y))
        assert ev[0]>=.4-1e-12 and ev[-1]<=.6+1e-12
        count+=1
    # Exact baseline identity for arbitrary auxiliary records.
    rng=np.random.default_rng(7)
    for _ in range(100):
        V=rng.normal(size=2); sigma=.17
        Fminus=r*(grad(t0+a*sigma*V)-grad(t1-a*sigma*V))
        B0=r*(grad(t0+a*sigma*V)-grad(t0-a*sigma*V))
        Fplusneg=r*(grad(t0-a*sigma*V)-grad(t1-a*sigma*V))
        close(Fminus-B0,Fplusneg)
    cases.append(dict(a=a,q=q,common_norm=float(np.linalg.norm(common)),word_skew_over_ra=float(np.linalg.norm(DS-DS.T)/(r*a))))
# Pre-mixture and tail exponent checks.
for g in np.linspace(.01,1,100):
    grades=[3+2.5*g,4+1.5*g,6.5+3*g,7+3*g]
    assert min(grades)>=3+2.5*g-1e-12
    count+=1
    close((1+g)+(2+1.5*g),3+2.5*g)
out=dict(checks=count,max_relative_first_error=max_first_error,cases=cases,
         note="Diagnostics verify finite algebra, Hessian sandwich and exponents; no all-order endpoint closure is claimed.")
(ROOT/"terminal_flat_gate_checks.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
