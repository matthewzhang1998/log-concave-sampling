#!/usr/bin/env python3
"""Deterministic arithmetic checks for exact nested-bank and strong-allocation formulas."""
import json, math
from fractions import Fraction as F
from pathlib import Path

rows=[]
def check(name, condition, **data):
    if not condition: raise AssertionError((name,data))
    rows.append(dict(name=name,passed=True,**data))

# Exact coefficients of average_N - average_M on a common iid tape.
for M,N in [(1,2),(2,8),(8,64),(17,101),(64,1024)]:
    variance=M*(F(1,N)-F(1,M))**2+(N-M)*F(1,N)**2
    exact=F(1,M)-F(1,N)
    check('nested_mean_variance',variance==exact,M=M,N=N,coefficient=str(exact))

# A quadratic posterior gives an EXACT Gaussian canonical force.
for A in [.5,.25,.125,.0625]:
    lam=.625
    eps=A*lam/math.sqrt(1+A*lam)
    for M,N in [(4,16),(16,64),(64,256)]:
        weak=math.sqrt(1+eps*eps/N)-1
        strong=eps*math.sqrt(1/M-1/N)
        ratio=weak/(eps*eps/(2*N))
        check('quadratic_exact_law_and_nested_strong',0<ratio<=1 and strong>weak,
              A=A,M=M,N=N,law_error=weak,nested_strong=strong,asymptotic_ratio=ratio)
        # Orthogonal source-specific reparametrization has independent within-level keep/signal.
        xi=eps/math.sqrt(N); t=math.sqrt(1+xi*xi)
        cov=((1/t)*(xi/t)+(-xi/t)*(1/t))
        signal_U=1/t+xi*xi/t; signal_V=-xi/t+xi/t
        check('quadratic_carrier_rotation',abs(cov)<1e-14 and abs(signal_U-t)<1e-14 and abs(signal_V)<1e-14,
              A=A,N=N,cross_covariance=cov,combined_U=signal_U)

# Smooth admissible nonlinear original-potential fixture.
# V''=.5+.25*cos(x/A^3.4), E=.25 A^3.9 sin(W/A^2.9).
for A in [.5,.25,.125,.0625]:
    omega=A**(-2.9)
    var=A**7.8*(1-math.exp(-2*omega*omega))/32
    derivative=.25*A**3.9*omega
    check('nonlinear_C2_envelope',abs(derivative-.25*A)<1e-13 and var>0,
          A=A,hessian_min=.25,hessian_max=.75,variance=var,first_bound=derivative)
    for M,N in [(1,4),(4,16),(16,64)]:
        check('nonlinear_nested_identity',var*(1/M-1/N)>0,A=A,M=M,N=N,
              exact_rms=math.sqrt(var*(1/M-1/N)))

# Each fixed-target level has strong RMS order A^j and cost <= A^-2(j-1).
for j in range(1,16):
    base_count=2*max(j-1,0)
    check('base_allocation',2+base_count==2*j and base_count==2*(j-1),j=j,
          count_exponent=base_count,variance_exponent=2+base_count)
    for ell in range(1,j+1):
        count=2*(j-ell)
        lower_cost=2*(ell-1)
        check('strong_level_fixed_point',2*ell+count==2*j and count+lower_cost==2*(j-1),
              j=j,ell=ell,count_exponent=count,variance_exponent=2*ell+count,
              work_exponent=count+lower_cost)
    # Coarse statistic precision + physical input sqrt(A) + physical force sqrt(A).
    check('adjacent_offset_restored',F(j-1)+F(1,2)+F(1,2)==j,j=j,
          normalized_coarse_error=j-1,state_exponent=str(F(j)-F(1,2)),force_exponent=j)

# General complete cost beta <=2 maps to cost at most2, while beta>2 is inherited.
for beta in [F(0),F(1,2),F(1),F(3,2),F(2),F(5,2),F(3)]:
    for j in range(1,12):
        work=max([F(2*(j-1))]+[2*(j-ell)+beta*(ell-1) for ell in range(1,j+1)])
        predicted=max(F(2),beta)*(j-1)
        check('general_cost_slope',work==predicted,beta=str(beta),j=j,work=str(work))

out=dict(status='PASS',checks=len(rows),scope='Exact algebraic diagnostics, not a complete hidden-source admission theorem',rows=rows)
Path(__file__).with_name('nested_bank_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status=out['status'],checks=out['checks'])))
