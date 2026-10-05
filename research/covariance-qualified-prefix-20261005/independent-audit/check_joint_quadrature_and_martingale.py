#!/usr/bin/env python3
"""Numerical corroboration of the independent joint-bridge audit.

No Monte Carlo, derivative oracle, or external dependency beyond numpy/scipy.
The analytic proofs are in JOINT-QUADRATURE-AND-MARTINGALE-AUDIT.md.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def bridge_rows(s, delta):
    return np.sinh(delta - s) / math.sinh(delta), np.sinh(s) / math.sinh(delta)


def conditional_cos_kernel(s, t, delta, ell2):
    a, b = bridge_rows(s, delta)
    c, d = bridge_rows(t, delta)
    q = math.exp(-delta)
    retained = a*c + b*d + q*(a*d + b*c)
    full = np.exp(-np.abs(np.asarray(s) - np.asarray(t)))
    # Tiny roundoff near the diagonal must not create exponent overflow.
    retained = np.minimum(retained, 1.0)
    full = np.minimum(full, 1.0)
    def pair(z):
        return 0.5 * (np.exp(-(1-z)/ell2) + np.exp(-(1+z)/ell2))
    return pair(full) - pair(retained)


def midpoint_rule(delta, n):
    edges = np.linspace(0.0, delta, n+1)
    nodes = (edges[:-1] + edges[1:]) / 2
    weights = np.exp(-edges[:-1]) - np.exp(-edges[1:])
    return edges, nodes, weights


def true_cos_variance(delta, ell2):
    # Integrate the positive ordered triangle using d=ell2*z to resolve
    # the diagonal boundary layer even when n is large.
    zmax = min(delta/ell2, 80.0)
    def outer(z):
        gap = ell2*z
        upper = delta-gap
        if upper <= 0:
            return 0.0
        boundary = min(30*ell2, upper/4)
        points = sorted(set([0.0, boundary, upper-boundary, upper]))
        result = 0.0
        for left, right in zip(points[:-1], points[1:]):
            result += quad(
                lambda s: math.exp(-2*s-gap)*float(
                    conditional_cos_kernel(s, s+gap, delta, ell2)),
                left, right, epsabs=1e-13, epsrel=1e-9, limit=100)[0]
        return 2*ell2*result
    result = quad(outer, 0.0, zmax, epsabs=1e-14, epsrel=1e-8, limit=100)[0]
    # The omitted positive tail is at most this conservative analytic bound.
    tail_bound = 4*delta*ell2*math.exp(-zmax/2) if zmax < delta/ell2 else 0.0
    return result, tail_bound


def cosine_checks():
    out = []
    c0 = 0.45
    theta = c0*c0/32
    for w in [0.02, 0.1, 0.4]:
        delta = -math.log1p(-w)
        for n in [8, 16, 32, 64]:
            _, nodes, weights = midpoint_rule(delta, n)
            middle = (nodes >= delta/4) & (nodes <= 3*delta/4)
            middle_fraction = float(weights[middle].sum()/w)
            assert middle_fraction >= c0
            ell2 = theta*w/n
            kernel = conditional_cos_kernel(nodes[:, None], nodes[None, :], delta, ell2)
            min_kernel = float(kernel.min())
            assert min_kernel >= -1e-11
            qvar = float(weights @ kernel @ weights)
            pvar, tail_bound = true_cos_variance(delta, ell2)
            # A=1 and f=(ell/2)(cos(x/ell)-1).
            defect = ell2*(qvar-pvar)/4
            certified_defect = c0**4*w**3/(1024*n*n)
            assert defect - ell2*tail_bound/4 >= certified_defect*(1-1e-8)
            assert pvar <= 4*w*ell2*(1+1e-8)
            out.append(dict(w=w, J=n, middle_fraction=middle_fraction,
                            ell2=ell2, minimum_kernel=min_kernel,
                            quadrature_cos_variance=qvar,
                            true_cos_variance=pvar,
                            omitted_cos_variance_tail_bound=tail_bound,
                            oscillatory_defect=defect,
                            rigorous_lower_bound=certified_defect,
                            scaled_defect=defect*n*n/w**3))
    return out


def spectral_quadrature_error(delta, nodes, weights, lam):
    if abs(lam-1.0) < 1e-12:
        continuous = (1-(1+2*delta)*math.exp(-2*delta))/2
        left = nodes*np.exp(-nodes)
    else:
        continuous = 2/(1-lam)*(
            -math.expm1(-(1+lam)*delta)/(1+lam)
            + math.expm1(-2*delta)/2)
        left = (np.exp(-nodes)-np.exp(-lam*nodes))/(lam-1)
    right = (np.exp(-nodes)-np.exp(lam*nodes-(lam+1)*delta))/(lam+1)
    discrete = float(weights @ np.exp(-lam*np.abs(nodes[:, None]-nodes[None, :])) @ weights)
    return continuous - 2*float(weights @ (left+right)) + discrete


def spectral_checks():
    out = []
    for delta in [0.01, 0.1, 0.6]:
        for n in [2, 8, 32]:
            edges, nodes, weights = midpoint_rule(delta, n)
            f2 = 0.0
            for left, node, right in zip(edges[:-1], nodes, edges[1:]):
                f2 += quad(lambda s: (math.exp(-left)-math.exp(-s))**2,
                           left, node, epsabs=1e-18)[0]
                f2 += quad(lambda s: (math.exp(-right)-math.exp(-s))**2,
                           node, right, epsabs=1e-18)[0]
            assert f2 <= delta*(delta/n)**2*(1+1e-12)
            for lam in [0.3, 1.0, 10.0, 100/delta]:
                actual = spectral_quadrature_error(delta, nodes, weights, lam)
                bound = 2*lam*f2
                assert actual >= -1e-13
                assert actual <= bound*(1+1e-7)+1e-13
                out.append(dict(delta=delta, J=n, spectral_frequency=lam,
                                spectral_error=actual, cumulative_integral=f2,
                                upper_bound=bound))
    return out


def martingale_checks():
    out = []
    for delta in [0.01, 0.03, 0.1, 0.4, math.log(2)]:
        def c(s):
            r = delta-s
            if r < 1e-4:
                inner = r/2-r*r/6+r**4/90-r**6/945
            else:
                inner = 0.5-r/math.expm1(2*r)
            return math.exp(-s)*inner
        martingale = 2*quad(lambda s: c(s)**2, 0, delta, epsabs=1e-15)[0]
        q = math.exp(-delta)
        exact = (1-q*q)/4-q*q*delta*delta/(1-q*q)
        relative_error = abs(martingale-exact)/martingale
        assert relative_error < 1e-8
        out.append(dict(delta=delta, q=q, conditional_linear_variance=exact,
                        martingale_square_integral=martingale,
                        relative_error=relative_error,
                        leading_small_time_ratio=martingale/(delta**3/6)))
    return out


if __name__ == '__main__':
    report = dict(status='passed', cosine_obstruction=cosine_checks(),
                  spectral_strong_bound=spectral_checks(),
                  martingale_square=martingale_checks())
    destination = Path(__file__).with_name('joint-quadrature-and-martingale-checks.json')
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],
                          cosine_cases=len(report['cosine_obstruction']),
                          spectral_cases=len(report['spectral_strong_bound']),
                          martingale_cases=len(report['martingale_square']),
                          output=str(destination)), indent=2))
