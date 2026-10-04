#!/usr/bin/env python3
"""Finite diagnostics for the opposite-baseline full-tape and bank gates."""
import json, math
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent
count=0
maxerr=0.

def close(a,b,tol=1e-11):
    global count,maxerr
    err=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    assert err <= tol, (a,b,err,tol)
    count+=1; maxerr=max(maxerr,err)

def chi(z):
    if abs(z)>=1: return 0.
    return math.exp(1-1/(1-z*z))

def dchi(z):
    if abs(z)>=1: return 0.
    return -2*z/(1-z*z)**2*chi(z)

for z in np.linspace(-.9999,.9999,10001):
    assert abs(dchi(z)) <= 3
    count+=1

# Scalar same-potential completion fixture. All original Hessians stay in [.45,.55].
eta=1/80
c,s=.6,.8
fixture=[]
for a in [1/16,1/64,1/256]:
    r=a; kap=r*a
    for h in [a/2,.08,.2]:
        zeta=h/480
        for k in [64.,256.,1024.,4096.,16384.]:
            def g(y):
                return .5*y+eta/k*(1-math.cos(k*y))*chi(8*y)+zeta*chi(4*(y-1-h)/h)
            def H(y):
                return .5+eta*math.sin(k*y)*chi(8*y)+(8*eta/k)*(1-math.cos(k*y))*dchi(8*y)+(4*zeta/h)*dchi(4*(y-1-h)/h)
            close(g(0.),0.)
            close(H(0.),.5)
            R=g(1+h)+g(1-h)-2*g(1.)
            Q=H(1+h)+H(1-h)-2*H(1.)
            close(R,zeta)
            close(Q,0.)
            # Exact second root derivative of potential(V=1)-potential(V=0).
            d2=kap/(2*h*h)*(a*s*s*eta*k*R+(a*s*.5)**2*Q)
            exact=kap*a*s*s*eta*k/(960*h)
            close(d2,exact,5e-10)
            # Independent finite-difference check of g''(0).
            dh=1e-4/k
            gpp=(H(dh)-H(-dh))/(2*dh)
            close(gpp/(eta*k),1.,1e-7)
            # Actual normalized physical first at a terminal-bump flank.
            V=1+.5/4
            tS=1+a*c*.5
            dS=kap/(2*h)*(H(1+h*V)-H(1-h*V))*tS
            exactS=kap*dchi(.5)/(240*h)*tS
            close(dS,exactS,1e-10)
            # Check global sandwich on a resolved ancestor-frequency grid and terminal grid.
            for y in np.r_[np.linspace(-.13,.13,1001),1+h+(h/4)*np.linspace(-1,1,1001)]:
                assert .45-1e-12 <= H(float(y)) <= .55+1e-12
                count+=1
            fixture.append(dict(a=a,h=h,k=k,full_gradient_first_lower=d2/2,physical_first_over_kappa_per_h=abs(dS)/(kap/h)))

# New terminal-flat 2D fixture; its second output is globally linear.
T=np.array([[.5,.05],[.05,.5]])
N=np.diag([1.,0.]); delta=1/40; e1=np.array([1.,0.])
rng=np.random.default_rng(113)
flat=[]
for a in [1/16,1/64,1/256,1/1024]:
    r=a; kap=r*a; eps=a**.9; q=a*eps
    def psi(z): return (z+.5)*chi(10*(z+.5))
    def dpsi(z): return chi(10*(z+.5))+10*(z+.5)*dchi(10*(z+.5))
    def g(y): return T@y+delta*q*psi(y[0]/q)*e1
    def H(y): return T+delta*dpsi(y[0]/q)*N
    S=np.zeros(2);x0=np.zeros(2);Z=e1
    Delta=g(S)-g(S+eps*Z);x1=x0+a*Delta
    t0=S+a*g(x0);t1=S+a*g(x1)
    DA=H(x0)-H(x1)
    word=kap/2*(H(t0)+H(t1))@DA
    close(word,-kap*delta*T@N)
    close(H(t0)-H(t1),np.zeros((2,2)))
    for sigma in [.0007,.01,.1,.9,2.7]:
        for _ in range(25):
            V=rng.normal(size=2)
            B=r*(g(t0+a*sigma*V)-g(t0-a*sigma*V))
            BN=B/(2*sigma)
            close(BN[1],kap*(T@V)[1],1e-12)
            Fm=r*(g(t0+a*sigma*V)-g(t1-a*sigma*V))
            Cp=r*(g(t0-a*sigma*V)-g(t1-a*sigma*V))
            close(Fm-B,Cp)
            # Common heat second coordinate contains no auxiliary Gaussian row.
            close(Cp[1],r*(g(t0)-g(t1))[1])
    # Signed same-bank combinations with sum 1 preserve exact baseline row.
    for _ in range(30):
        weights=rng.normal(size=7);weights[-1]+=1-sum(weights)
        widths=np.exp(rng.uniform(-5,2,size=7));V=rng.normal(size=2)
        out=sum(w*r*(g(t0+a*sig*V)-g(t0-a*sig*V))/(2*sig) for w,sig in zip(weights,widths))
        close(out[1],kap*(T@V)[1],2e-12)
    flat.append(dict(a=a,word_skew_over_kappa=float(np.linalg.norm(word-word.T)/kap),baseline_linear_row_energy_over_kappa=float(np.linalg.norm(T[1,:]))))

# Independent-bank adapter normalization: alpha dot beta=1, ||beta||<=1.
banks=[]
for n in [1,2,3,7,19,64]:
    for _ in range(30):
        beta=rng.normal(size=n);beta/=np.linalg.norm(beta)
        extra=rng.normal(size=n);extra-=beta*np.dot(beta,extra)
        alpha=beta+extra
        close(np.dot(alpha,beta),1.,2e-12)
        assert np.linalg.norm(alpha)>=1-1e-12
        count+=1
    alpha=np.ones(n)/n
    beta_unit=np.ones(n)/math.sqrt(n)
    close(np.dot(alpha,beta_unit),1/math.sqrt(n))
    close(np.linalg.norm(alpha),1/math.sqrt(n))
    banks.append(dict(n=n,coefficient_average_energy=1/math.sqrt(n),unit_input_selected_coefficient=1/math.sqrt(n),recovery_amplification=math.sqrt(n)))

out=dict(checks=count,max_absolute_identity_error=maxerr,completion_fixture=fixture,terminal_flat_cases=flat,bank_cases=banks,
         scope='Finite diagnostics only. Proof is in OPPOSITE-BASELINE-FULL-TAPE-AND-LINEAR-RETURN-GATES.md. No arbitrary-frame or all-program impossibility is claimed.')
(OUT/'opposite_baseline_return_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(checks=count,max_absolute_identity_error=maxerr,fixture_cases=len(fixture),terminal_flat_cases=flat,bank_cases=banks),indent=2))
