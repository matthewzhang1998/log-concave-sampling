#!/usr/bin/env python3
"""Independent finite diagnostics; the analytical audit is the proof review."""
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec

B = 1.0 / 30.0


def posterior(theta, A, eta=.25, two_dim=False):
    theta = np.asarray(theta, float)
    n = len(theta)
    w = 1.0 / n
    starts = 1 + np.arange(n) * w
    s = 1 + .5 * A
    if two_dim:
        s -= (A / 16)**2 / (1 + A / 2)

    def integrand(u):
        t = np.clip((u - starts) / w, 0, 1)
        phi = t*t*(1-t)**2
        Phi = t**3/3 - t**4/2 + t**5/5
        R = eta*A*w*w*np.dot(theta, Phi)
        p = math.exp(-s*u*u/2-R)
        # In 2D this is the conditional expectation of the first force.
        force = (s-1)*u + eta*A*w*np.dot(theta, phi)
        return p*np.concatenate(([1, u, force], Phi, u*Phi))

    endpoints = [-np.inf, 0, 1, *list(starts[1:]), 2, 3, np.inf]
    raw = sum((quad_vec(integrand, lo, hi, epsabs=1e-13, epsrel=2e-12)[0]
               for lo, hi in zip(endpoints[:-1], endpoints[1:])),
              start=np.zeros(3+2*n))
    moments = raw / raw[0]
    mean_u = moments[1]
    cov = moments[3+n:] - mean_u*moments[3:3+n]
    return -mean_u, moments[2], eta*A*w*w*cov, cov


def main():
    rng = np.random.default_rng(192873)
    max_ibp = 0.0
    max_derivative_relative_error = 0.0
    min_covariance = float('inf')
    cases = 0
    for two_dim in (False, True):
        eta = .125 if two_dim else .25
        for A in (.03, .4, 1.0):
            for n in (2, 3, 4):
                for theta in (np.ones(n), -np.ones(n), rng.uniform(-.7, .7, n)):
                    mu, direct_force, deriv, cov = posterior(theta, A, eta, two_dim)
                    max_ibp = max(max_ibp, abs(mu-direct_force))
                    min_covariance = min(min_covariance, float(min(cov)))
                    for i in range(n):
                        h = 1e-3
                        plus, minus = theta.copy(), theta.copy()
                        plus[i] += h
                        minus[i] -= h
                        fd = (posterior(plus, A, eta, two_dim)[0]
                              - posterior(minus, A, eta, two_dim)[0])/(2*h)
                        max_derivative_relative_error = max(
                            max_derivative_relative_error, abs(fd-deriv[i])/deriv[i])
                    cases += 1

    # All possible transcripts, including random/optional-query patterns,
    # have uniform remaining signs. Check every partial assignment at n=4.
    n, A = 4, .4
    signs = list(itertools.product((-1, 1), repeat=n))
    targets = {s: posterior(s, A)[0] for s in signs}
    min_walsh = float('inf')
    min_projection_slack = float('inf')
    partial_cases = 0
    for partial in itertools.product((-1, 0, 1), repeat=n):
        unknown = [i for i, s in enumerate(partial) if s == 0]
        allowed = [s for s in signs if all(p == 0 or p == s[i]
                                          for i, p in enumerate(partial))]
        values = np.array([targets[s] for s in allowed])
        bs = [float(np.mean(values * np.array([s[i] for s in allowed])))
              for i in unknown]
        if bs:
            min_walsh = min(min_walsh, min(bs))
        min_projection_slack = min(min_projection_slack,
                                   float(np.var(values) - sum(b*b for b in bs)))
        partial_cases += 1

    T = np.array([[.5, 1/16], [1/16, .5]])
    E = np.diag([1., 0.])
    extremum = math.sqrt(3)/9
    Hplus, Hminus = T + .125*extremum*E, T - .125*extremum*E
    eigs = np.concatenate([np.linalg.eigvalsh(Hplus), np.linalg.eigvalsh(Hminus)])
    commutator = Hplus@Hminus-Hminus@Hplus
    explicit_c0 = math.exp(-87/16-1/60)/(60*math.pi)
    assert max_ibp < 2e-12
    assert max_derivative_relative_error < 1e-5
    assert min_covariance > explicit_c0 > 0
    assert min_walsh > 0 and min_projection_slack > -1e-20
    assert min(eigs) >= 5/16 and max(eigs) <= 11/16
    assert np.linalg.norm(commutator) > 0
    report = dict(posterior_cases=cases, partial_transcripts=partial_cases,
                  max_ibp_absolute_error=max_ibp,
                  max_normalized_derivative_relative_error=max_derivative_relative_error,
                  min_observed_covariance=min_covariance,
                  rigorous_scalar_covariance_lower_bound=explicit_c0,
                  min_first_walsh_coefficient=min_walsh,
                  min_projection_slack=min_projection_slack,
                  noncommuting_hessian_spectrum=[float(min(eigs)), float(max(eigs))],
                  noncommuting_hessian_commutator_norm=float(np.linalg.norm(commutator)),
                  all_checks_passed=True)
    out = Path(__file__).with_name('numerical-check-results.json')
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
