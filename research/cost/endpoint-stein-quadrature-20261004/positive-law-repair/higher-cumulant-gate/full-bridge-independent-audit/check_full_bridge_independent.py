#!/usr/bin/env python3
"""Independent algebra/numerical diagnostics for the bounded cubic bridge join.
These diagnostics do not implement imported finite mean/square compilers.
"""
import json
from pathlib import Path
import numpy as np
import mpmath as mp
mp.mp.dps=80
checks=[]
def check(name, flag, **data):
    checks.append(dict(name=name, passed=bool(flag), **data))
    assert flag, (name,data)

# Exact affine and variance ledgers at ordinary and near-terminal points.
for r,t in [(0,.5),(.1,.35),(.5,.65),(.8,.83),(.98,.983),(.999,.9992)]:
    s2=1-r*r; Delta=t*t-r*r
    a=r*(1-t*t)/(t*s2); b=Delta/(t*s2); v=(1-t*t)*Delta/(t*t*s2)
    h2=Delta/s2; v0=Delta/(t*t)
    check(f'affine-contraction-{r}-{t}',abs(a+b*r-r/t)<1e-13)
    check(f'raw-variance-{r}-{t}',abs(b*b*s2+v-v0)<1e-13)
    check(f'raw-carrier-scale-{r}-{t}',abs(b*np.sqrt(s2)-np.sqrt(v0*h2))<1e-13)
    if h2<=.25:
        d=1-2*h2
        check(f'positive-join-{r}-{t}',d>=.5 and abs(h2+d-(1-h2))<1e-13)

# Conditional-OU finite rows and first bounds with noncommuting PSD Hessians.
rng=np.random.default_rng(314159)
D=5
for N in [2,3,5,9]:
    times=np.geomspace(.012,.95,N)
    K=np.minimum.outer(times,times)/np.maximum.outer(times,times)-np.outer(times,times)
    eig,Q=np.linalg.eigh(K)
    check(f'conditional-Gram-PSD-{N}',eig.min()>-1e-13)
    rows=Q@np.diag(np.sqrt(np.maximum(eig,0)))
    check(f'conditional-row-bound-{N}',np.max(np.linalg.norm(rows,axis=1))<=1+1e-13)
    weights=rng.uniform(size=N); weights/=weights.sum()
    weighted=np.sqrt(weights)[:,None]*rows
    check(f'weighted-row-op-{N}',np.linalg.norm(weighted,2)<=1+1e-13)
    for A in [.01,.1,.4]:
        # Arbitrary PSD site Hessians need not commute. A double positive stack
        # bounds the derivative of a nested feedback source uniformly in N.
        HS=[]
        for k in range(2*N):
            q,_=np.linalg.qr(rng.normal(size=(D,D)))
            HS.append(q@np.diag(rng.uniform(0,A,D))@q.T)
        J=np.zeros((D,N*D))
        inner=sum(weights[j]*HS[N+j]@np.kron(rows[j:j+1],np.eye(D)) for j in range(N))
        for i in range(N):
            J+=weights[i]*HS[i]@(np.kron(rows[i:i+1],np.eye(D))-inner)
        check(f'nested-first-{N}-{A}',np.linalg.norm(J,2)<=A*(1+A)+1e-13)

# Geometric schedule, exact final contraction weights and one-energy floor sum.
schedule=[]
for apow in [1,2,3,5,8,12,20,40]:
    A=mp.mpf(10)**(-apow); rho=mp.mpf(3)/4
    J=0; s2=mp.mpf(1)
    while s2**mp.mpf('2.5')>A:
        J+=1; s2*=rho
    check(f'stopping-{apow}',s2**mp.mpf('2.5')<=A and (s2/rho)**mp.mpf('2.5')>A)
    check(f'alpha-min-{apow}',rho*A**mp.mpf('1.4')<A*s2<=A**mp.mpf('1.4'))
    r=mp.mpf(0); weighted=mp.mpf(0); total_delta=mp.mpf(0); floor=mp.mpf(0)
    for j in range(J):
        sj2=rho**j; Delta=sj2/4; t=mp.sqrt(1-rho**(j+1))
        local=Delta/t*A**3*sj2**mp.mpf('2.5')
        weighted+=t*local
        total_delta+=Delta
        floor+=Delta*sj2**3
        check(f'half-buffer-{apow}-{j}',abs(Delta/sj2-mp.mpf('.25'))<mp.mpf('1e-75'))
        r=t
    bound=A**3/(4*(1-rho**mp.mpf('3.5')))
    check(f'weighted-error-sum-{apow}',weighted<=bound)
    check(f'terminal-grade-{apow}',A**2*s2**mp.mpf('2.5')<=A**3)
    check(f'floor-sum-{apow}',floor<=1/(4*(1-rho**4)))
    check(f'heat-budget-{apow}',abs(total_delta+s2-1)<mp.mpf('1e-75'))
    schedule.append(dict(A=str(A),J=J,min_alpha=str(A*s2),weighted_error_over_A3=float(weighted/A**3)))

# Exact scalar quadratic screen for the ideal completed bridge references.
# Mean: continuous second substitution. Reserve: positive finite square action
# with its paid degree-four covariance term. Terminal: old shared-root packet.
quadratic=[]
for apow in [1,2,3,4,6,8,12]:
    A=mp.mpf(10)**(-apow); rho=mp.mpf(3)/4
    J=0; sj2=mp.mpf(1)
    while sj2**mp.mpf('2.5')>A:
        J+=1; sj2*=rho
    V=mp.mpf(1)
    for j in range(J):
        s2=rho**j; alpha=A*s2; r=mp.sqrt(1-s2); t=mp.sqrt(1-rho**(j+1))
        Delta=s2/4; a=r*(1-t*t)/(t*s2); b=Delta/(t*s2); v0=Delta/(t*t)
        # Exact conditional mode; output conditional mean coefficient.
        k=a+b*r/(1+alpha)
        h2=mp.mpf(1)/4; d=1-2*h2; eta2=d/2
        mcoef=alpha/2-alpha**2/4
        C=alpha**2/4
        Tv=h2*(1-mcoef)**2+(1-h2)+h2*C+h2**2*C**2/(4*eta2)
        V=k*k*V+v0*Tv
    s2=rho**J; r=mp.sqrt(1-s2); alpha=A*s2
    k=r/(1+alpha)
    beta=mp.pi/4
    packet_var=(1-alpha/2)**2+alpha**2*beta**2
    V=k*k*V+s2*packet_var
    W=abs(mp.sqrt(V)-1/mp.sqrt(1+A))
    ratio=W/A**3
    check(f'quadratic-cubic-grade-{apow}',ratio<2)
    quadratic.append(dict(A=str(A),J=J,W2=str(W),W2_over_A3=float(ratio)))

out=dict(status='PASS',assertions=len(checks),checks=checks,schedule=schedule,quadratic=quadratic,
         scope='Algebra/first-bound/quadratic diagnostics; imported finite compilers and nonlinear proof not numerically instantiated.')
path=Path(__file__).with_name('full_bridge_independent_checks.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status='PASS',assertions=len(checks),max_quadratic_ratio=max(x['W2_over_A3'] for x in quadratic),schedule=schedule),indent=2))
