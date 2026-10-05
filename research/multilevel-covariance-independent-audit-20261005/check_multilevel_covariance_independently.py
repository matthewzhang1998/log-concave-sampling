#!/usr/bin/env python3
"""Independent scalar/operator diagnostics. Numerical checks are not proofs."""
import json, math
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.linalg import eigh
mp.mp.dps=70

def k_matrix(delta, span, z):
    q=np.exp(-delta); r=np.exp(-span); L=delta-span
    u=z-L/2; f=np.exp(-L/2)
    return f*np.array([[np.sqrt((1+r)/(1+q))*np.cosh(u), -np.sqrt((1+r)/(1-q))*np.sinh(u)],[-np.sqrt((1-r)/(1+q))*np.sinh(u),np.sqrt((1-r)/(1-q))*np.cosh(u)]],dtype=complex)

def mean_coeff(d):
    d=mp.mpf(d)
    a=mp.mpf('0.5')-d/mp.expm1(2*d)
    b=(d/2+mp.expm1(-2*d)/4)/mp.sinh(d)
    return a,b

def bridge_variance(d):
    d=mp.mpf(d)
    def f(s):
        t=d-s
        a=mp.mpf('0.5')-t/mp.expm1(2*t) if t else mp.mpf(0)
        return 2*mp.exp(-2*s)*a*a
    return mp.quad(f,[0,d])

def midpoint_increment(d):
    d=mp.mpf(d); a,b=mean_coeff(d/2)
    sig2=(1-mp.exp(-d))/(1+mp.exp(-d))
    return sig2*(b+mp.exp(-d/2)*a)**2

def geom_mass(N,h):
    return -mp.expm1(-2*mp.mpf(N)*h)/(-mp.expm1(-2*h))

def geometric_power_sums(N,r,m):
    """Moments sum_{j<N} j^p r^j using endpoint summation recurrence.
    Uses O(m^2) scalar operations and N only as a number, never enumerates N.
    """
    N=mp.mpf(N); r=mp.mpf(r)
    S=[(1-r**N)/(1-r)]
    for p in range(1,m+1):
        num=r*sum(mp.binomial(p,l)*S[l] for l in range(p))-N**p*r**N
        S.append(num/(1-r))
    return S

def main():
    rng=np.random.default_rng(20261005)
    max_norm=0.; max_trace_identity=0.; min_gap=1.
    for _ in range(2000):
        delta=float(np.exp(rng.uniform(np.log(1e-4),np.log(.69))))
        span=delta*float(rng.uniform(.0001,.9999)); L=delta-span
        x=L*float(rng.uniform(1e-4,1-1e-4)); y=float(rng.uniform(-1,1))*math.sqrt(x*(L-x))
        K=k_matrix(delta,span,x+1j*y); sv=np.linalg.svd(K,compute_uv=False)
        max_norm=max(max_norm,float(sv[0]))
        q=np.exp(-delta); r=np.exp(-span)
        expected=2*np.exp(-L)/(1-q*q)*(np.cosh(2*x-L)-r*q*np.cos(2*y))
        max_trace_identity=max(max_trace_identity,abs(np.sum(np.abs(K)**2)-expected))
        min_gap=min(min_gap,float(np.sinh(x)*np.sinh(L-x)-r*q*np.sin(y)**2))
    linear=[]
    for d in map(mp.mpf,['0.0001','0.03','0.2','0.69']):
        total=bridge_variance(d); accum=mp.mpf(0)
        for k in range(9):
            N=2**k; h=d/N
            accum+=geom_mass(N,h)*midpoint_increment(h)
            rem=geom_mass(2*N,h/2)*bridge_variance(h/2)
            linear.append({'delta':str(d),'levels':k+1,'decomposition_abs_error':float(abs(total-accum-rem)),'tail_over_delta3_4pow':float(rem/(d**3*mp.mpf(4)**(-(k+1))))})
    moments=[]
    for N in [2,7,31,2**24]:
        r=mp.exp(-mp.mpf('0.6')/N)
        S=geometric_power_sums(N,r,12)
        # independent generating-function derivatives, not same recurrence
        exact=[mp.diff(lambda t:(1-(r*mp.exp(t))**N)/(1-r*mp.exp(t)),mp.mpf(0),p) for p in range(13)]
        err=max(abs(a-b)/max(mp.mpf(1),abs(b)) for a,b in zip(S,exact))
        moments.append({'N':N,'max_degree':12,'relative_error':float(err)})
    out={'scope':'independent finite diagnostics; does not certify native square graph','operator_contraction':{'samples':2000,'max_singular_value':max_norm,'maximum_trace_formula_error':max_trace_identity,'minimum_wedge_margin':min_gap},'linear_conditional_variance':linear,'geometric_moment_recurrence':moments}
    assert max_norm<1+1e-8
    assert max(x['decomposition_abs_error'] for x in linear)<1e-50
    assert max(x['tail_over_delta3_4pow'] for x in linear)<1/6+1e-8
    assert max(x['relative_error'] for x in moments)<1e-45
    dest=Path(__file__).with_name('independent_checks.json'); dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
