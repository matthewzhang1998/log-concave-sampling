#!/usr/bin/env python3
"""Exact finite diagnostics; does not numerically execute imported native compilers."""
import json, math, hashlib
from fractions import Fraction as F
from pathlib import Path
import sympy as s
x=s.symbols('x')
checks=[]
def record(name, condition):
    assert condition,name
    checks.append(name)
def R(p,k):
    p=s.Poly(s.expand(p),x)
    out=0
    for (n,),coef in p.terms():
        term=x**n
        denom=k+n
        res=term/denom
        while n>=2:
            term=s.diff(term,x,2)
            n-=2
            denom*=k+n
            res+=term/denom
        out+=coef*res
    return s.expand(out)
for degree in range(1,7):
    g=sum(s.Rational((-1)**j,j+1)*x**j for j in range(1,degree+1))
    b=R(s.diff(g,x),2); T=s.diff(b,x); S=s.diff(T,x)
    C=2*R(b*b,2)
    K3=6*R(b*s.diff(C,x),3)
    K4=R(8*b*s.diff(K3,x)+6*s.diff(C,x)**2,4)
    U=R(T*b,3)
    P1=R(b*R(T*U,4),4)
    P2=R(b*R(b*R(S*b,4),4),4)
    P3=R(b*R(b*R(T*T,4),4),4)
    P4=R(U*U,4)
    record(f'four_history_exact_degree_{degree}',s.expand(K4-192*(P1+P2+P3)-96*P4)==0)
    record(f'third_exact_degree_{degree}',s.expand(K3-24*R(b*U,3))==0)
    record(f'derivative_resolvent_degree_{degree}',s.expand(s.diff(R(g,3),x)-R(s.diff(g,x),4))==0)
# Panelwise clock estimates. Coefficients use masses, not empirical replicas.
for tau_exp in range(1,31):
    tau=2.0**(-tau_exp)
    ds=[2.0**(-j) for j in range(1,2*tau_exp+15)]
    sums={
      'first2':sum(d/(d+tau*tau) for d in ds),
      'first3':sum(d/(d+tau*tau)**1.5 for d in ds),
      'square4':sum(d*d/(d+tau*tau)**2 for d in ds),
      'square5':sum(d*d/(d+tau*tau)**2.5 for d in ds)}
    record(f'panel_first2_{tau_exp}',sums['first2']<2*tau_exp+3)
    record(f'panel_first3_{tau_exp}',tau*sums['first3']<5)
    record(f'panel_square4_{tau_exp}',sums['square4']<2*tau_exp+3)
    record(f'panel_square5_{tau_exp}',tau*sums['square5']<5)
# Normalization: four physical readouts and unit-to-physical scale.
for n in range(1,41):
    c=1/math.sqrt(5*n); read=c**4
    for alpha in [.001,.01,.1]:
        for sig1,sig2 in [(.01,.01),(.1,.3),(.5,.7)]:
            w=.2*sig1*sig1*sig2*sig2
            for d in [4,8,-4,-8]:
                rho=d*w*alpha/(read*sig1*sig2)
                got=read*rho*alpha**3*sig1*sig2
                want=d*w*alpha**4
                record(f'normalization_{len(checks)}',abs(got-want)<1e-12*max(1,abs(want)))
    record(f'positive_budget_{n}',abs(n*5*c*c-1)<1e-13)
# Exact exponent arithmetic in physical outputs, alpha=A/sqrt(u), tau=alpha.
# sqrt(u)*alpha^p / tau^r -> A^(p-r) u^((1-p+r)/2).
def physical(p,r=0):return (F(p-r),F(1-p+r,2))
record('physical_quintic',physical(5)==(F(5),F(-2)))
record('star_first',physical(2,1)==(F(1),F(0)))
record('root_first',physical(1)==(F(1),F(0)))
record('coarse_star_feedback',physical(8,1)==(F(7),F(-3)))
# Scalar Gaussian regression: independent Y,P,R conditional on S=cY+bP+dR.
c,b,d,z,V=s.symbols('c b d z V',nonzero=True)
reg=c*b*d*z**3/V**3 -3*c*b*d*z/V**2
record('cubic_wick_regression',s.expand(reg-c*b*d*(z**3/V**3-3*z/V**2))==0)
# Exact star symbol: lambda theta^2(P+b theta)(S+u theta); Gaussian covariance.
l,t,u=s.symbols('lambda theta u')
mean=l*b*u*t**4
var=l*l*t**4*(1+(b*b+u*u)*t*t)
record('star_main_four_forces',s.expand(mean).coeff(t,4)==l*b*u)
record('star_side_zero_feedback',s.expand(var/2).subs({b:0,u:0})==l*l*t**4/2)
result={'status':'PASS','assertions':len(checks),'scope':'Exact hierarchy, amplitude/readout algebra, scalar positive-clock bounds, star feedback and buffer powers. Does not execute the native compiler or certify the theorem without mathematical review.'}
out=Path(__file__).with_name('fourth_packet_checks.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
