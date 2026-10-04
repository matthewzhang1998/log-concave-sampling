"""Exact algebra/numerical fixtures for the same-heat quarter-flow candidate.

These tests do not certify a nonlinear arbitrary-order flow compiler.
"""
import json
import math
from pathlib import Path
import numpy as np

T = math.pi / 2
rows = []
for a in (0.1, 0.03, 0.01, 0.003, 0.001):
    lam = 0.5
    omega = math.sqrt(1 + a * lam)
    c = math.cos(T * omega)
    s = math.sin(T * omega) / omega
    target_variance = a / (1 + a * lam)
    exact_variance = c*c*target_variance + a*s*s
    independent_correction_variance = c*c*target_variance + a*(1+(s-1)**2)
    independent_w2 = math.sqrt(independent_correction_variance)-math.sqrt(target_variance)
    rows.append(dict(a=a, c=c, c_over_a=c/a,
                     bound=a/(1-a),
                     covariance_identity_error=exact_variance-target_variance,
                     independent_bank_w2=independent_w2,
                     independent_bank_w2_over_a_1p5=independent_w2/a**1.5))
    assert abs(c) <= a/(1-a)
    assert abs(exact_variance-target_variance) < 1e-14

# At the first Picard endpoint, u=sin(t) turns the force integral into
# integral_0^1 h(u)du.  Each gap between queried points supports a C1 bump
# psi(u)=ell*r^2(1-r)^2.  Both psi and psi' vanish at every queried point.
# Thus original gradient and directional-Hessian observations agree there.
bumps = []
for n in (1, 2, 4, 8, 16, 32):
    nodes = (np.arange(n)+0.5)/n
    boundaries = np.r_[0.0,nodes,1.0]
    widths = np.diff(boundaries)
    integral = float(np.sum(widths**2)/30)
    lower_bound = 1/(30*(n+1))
    assert integral >= lower_bound-1e-15
    r = np.linspace(0,1,10001)
    derivative = 2*r*(1-r)*(1-2*r)
    assert np.max(np.abs(derivative)) < 0.2
    bumps.append(dict(queries=n, integral_psi=integral,
                      guaranteed_integral=lower_bound,
                      normalized_two_potential_endpoint_gap_coefficient=2*integral,
                      hessian_min=0.5-float(np.max(np.abs(derivative))),
                      hessian_max=0.5+float(np.max(np.abs(derivative)))))

# Product-integration coefficients have exactly the declared nonnegative
# row sums. Their final carrier is the fresh momentum, not an old-source copy.
weights = []
for n in (4,16,64):
    ts = np.linspace(0,T,n+1)
    max_err = 0.0
    min_weight = math.inf
    for i in range(1,n+1):
        w = np.cos(ts[i]-ts[1:i+1])-np.cos(ts[i]-ts[:i])
        min_weight=min(min_weight,float(w.min()))
        max_err=max(max_err,abs(float(w.sum())-(1-math.cos(ts[i]))))
    assert min_weight > 0 and max_err < 1e-12
    weights.append(dict(intervals=n,min_weight=min_weight,max_row_sum_error=max_err))

out=dict(status="All displayed algebraic fixtures passed; nonlinear compiler remains open",
         quadratic=rows,quadrature_blind_bumps=bumps,product_integration=weights)
dest=Path(__file__).with_name('same_heat_refinement_checks.json')
dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
