#!/usr/bin/env python3
"""Deterministic diagnostics for the radial proof; not a compiler implementation.

Checks explicit constants and direct Gaussian-row mean-field values. Positive
Legendre quadrature is diagnostic only, not the sealed finite positive rule.
"""
import json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

alpha=beta=.5
cv=alpha*(1-1/math.sqrt(2))
c1=alpha/2
tau=cv/2
u=cv/(2*c1)
cstar=(17/math.sqrt(2)-11)/192
k=math.sqrt(3/8)
checks=0

def check(cond):
    global checks
    assert bool(cond)
    checks+=1

check(abs(beta*c1*(u/4-1/6+u**3/3)+2*cstar)<1e-15)
A=1/36
K=3.5+14*A+5.5*A*A
check(K<4)
check(K+1+A+A*A<5)
check(1+(2+3/math.sqrt(2))*A+k*A*A<2)
check(2*math.sqrt(2)*A*K<12*A)
check(math.sqrt(2)*A*(1+2*A+3*A*A)<2*A)
check(3*cv/math.sqrt(2)+2*c1<.819)
check(2*beta*tau<.075)
check((13+1.25)/2<8)
check(cstar-8/10000-1/(4*10000**2)-7/3000>1/500)

# Check the scalar terminal expansion estimates on their entire enclosing
# parameter range using algebraic endpoint inequalities, not random sampling.
h=A*(1+A)
check(h<.03)
check((1+(1+A)**2/(2-2*h))*A*A<2*A*A)
check((1+A)**2/((2-2*h)*(1-h))<1)
check((.5+h)+2<3)
check(2+(h+A*tau)<3)

# Mean-field Gaussian row inner products, respecting the actual joint history.
# Coordinates below label (x,v,w,U,V); first-chaos coefficient is source row's
# covariance with x, not its radial magnitude alone.
N=384
nodes, weights=leggauss(N)
r=(nodes+1)/2
wr=weights/2
R=r[:,None]
Z=r[None,:]
wp=wr[:,None]*wr[None,:]
a=-np.log(R)
s=a-np.log(Z)
# True retained-history covariances.
cvU=R*(a+.5)
cvV=R*Z*(s+.5)


def terminal(A, bv=0., bw=0., bu=0., bvpair=0., rU=0., rV=0.,
             vU=0., vV=0., wU=0., wV=0., UV=0.):
    # row Y=x+bv*v+bw*w+bu*U+bvpair*V
    n=1+.5*bv+.25*bw+bu*rU+bvpair*rV
    q2=(1+.5*bv*bv+.375*bw*bw+bu*bu+bvpair*bvpair
        +bv+.5*bw+2*bu*rU+2*bvpair*rV+.75*bv*bw
        +2*bv*bu*vU+2*bv*bvpair*vV+2*bw*bu*wU
        +2*bw*bvpair*wV+2*bu*bvpair*UV)
    q=np.sqrt(q2)
    phi=alpha*np.maximum(q-.5,0)+beta*np.maximum(q-(1-A*tau),0)
    return A*phi*n/q


def c(A,q):
    return (alpha*np.maximum(q-.5,0)+beta*np.maximum(q-(1-A*tau),0))/q


def reference_defect(A):
    cU=c1+beta*A*tau
    base=terminal(A,bv=-A*cv)
    # Single finite-difference average over a genuine single OU time.
    t=-np.log(r)
    vU=r*(t+.5)
    plus=terminal(A,bv=0,bu=-A*cU,rU=r,vU=vU)
    minus=terminal(A,bv=-2*A*cv,bu=A*cU,rU=r,vU=vU)
    single=np.sum(wr*(plus-minus)/2)
    # Genuine ordered Exp(2) x Exp(1) pair law has weight 2R dr dz.
    qsum=0.
    for sig,tig in ((1,1),(-1,-1),(1,-1),(-1,1)):
        val=terminal(A,bv=A*cv*(-1+sig+tig),bu=-sig*A*cU,
            bvpair=-tig*A*cU,rU=R,rV=R*Z,vU=cvU,vV=cvV,UV=Z)
        qsum+=sig*tig*np.sum(wp*2*R*val)/8
    # Linear-pair nested e: retain its original U-g(V) source.
    qh=np.sqrt(1-2*A*cU*Z+A*A*cU*cU)
    ch=c(A,qh)
    eU=A*(cU-ch)
    eV=A*A*ch*cU
    ep=terminal(A,bv=-A*cv,bu=eU,bvpair=eV,rU=R,rV=R*Z,
                vU=cvU,vV=cvV,UV=Z)
    em=terminal(A,bv=-A*cv,bu=-eU,bvpair=-eV,rU=R,rV=R*Z,
                vU=cvU,vV=cvV,UV=Z)
    linear=np.sum(wp*(ep-em)/2)
    q2=math.sqrt(1-A*cU+.5*A*A*cU*cU)
    c2=float(c(A,q2))
    target=terminal(A,bv=-A*c2,bw=A*A*c2*cU)
    defect=float(base+single+linear+qsum-target)
    check(q2<1-A*tau)
    check(abs(c2-c1)<.064*A)
    check(c2<=.25)
    return {"A":A,"inner_defect_over_A2":defect/(A*A),
            "outer_witness_over_A2":-defect/(2*A*A),
            "base":float(base),"single":float(single),
            "linear":float(linear),"quadratic":float(qsum),
            "target":float(target)}

rows=[reference_defect(A) for A in (.01,.003,.001,.0003,.0001)]
check(abs(rows[-1]['outer_witness_over_A2']-cstar)<1e-4)
check(rows[-1]['outer_witness_over_A2']>0)
result={"assertions":checks,"c_star":cstar,
        "inner_signed_coefficient":-2*cstar,
        "explicit_lower_ratio_n10000":cstar-8/10000-1/(4*10000**2)-7/3000,
        "claimed_lower_ratio":1/500,
        "legendre_diagnostic_nodes_per_axis":N,
        "mean_field_diagnostics":rows,
        "scope":"Diagnostics of constants and deterministic reference rows only; not the native compiler, sealed rule implementation, or a stochastic finite-D validation."}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
