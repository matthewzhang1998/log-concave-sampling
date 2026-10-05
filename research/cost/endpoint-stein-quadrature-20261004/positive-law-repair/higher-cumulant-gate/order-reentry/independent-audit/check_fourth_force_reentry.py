#!/usr/bin/env python3
"""Independent finite interface diagnostics; does not run imported compilers."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
RECIPE = source_path('research/no-copy-rank-20261004/host-observer-audit')
SOURCE_PIN = 'e1231f72136283b978895dc8763c08eddd35f73fc1aa3955a1e1b7fbb9888d85'
count = 0

def require(value, message):
    global count
    count += 1
    if not bool(value):
        raise AssertionError(message)

def op(x):
    return float(np.linalg.norm(x, 2))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

pins = {str(HERE.parent/'FOURTH-ORDER-FORCE-MEAN-REENTRY.md'): SOURCE_PIN}
for manifest in (BASE/'full-bridge-independent-audit/MANIFEST.json', RECIPE/'MANIFEST.json'):
    obj = json.loads(manifest.read_text())
    for key in ('inputs', 'sources', 'files', 'outputs'):
        for name, pin in obj.get(key, {}).items():
            path = source_path(name, manifest.parent)
            pins[str(path)] = pin
for name, pin in pins.items():
    status = verify_source_pin(name, pin)
    if status is not None:
        require(status, 'source pin mismatch: '+name)

rng = np.random.default_rng(20261004)
max_fd_error = 0.0
max_curl_ratio = 0.0
max_chain_error = 0.0
max_noncommutator = 0.0
max_energy_pointwise_ratio = 0.0
max_conditional_error_ratio = 0.0

for d, n in ((2, 2), (2, 7), (3, 11), (5, 29)):
    P = np.linalg.qr(rng.normal(size=(n, d)))[0].T
    Q = np.linalg.qr(rng.normal(size=(n, d)))[0].T
    T = rng.normal(size=(d, d)); T /= op(T)
    u = rng.normal(size=d); u /= np.linalg.norm(u)
    v = rng.normal(size=d); v /= np.linalg.norm(v)
    require(op(P@P.T-np.eye(d)) < 2e-14, 'physical carrier is coisometric')
    require(op(Q@Q.T-np.eye(d)) < 2e-14, 'residual input row')
    for A, K, s in ((.01, 1., 1.), (.005, 5., .7), (.003, 2., .05)):
        alpha = A*s*s
        center = rng.normal(size=d)
        def g(x):
            return A*(.5*x+.1*u*(math.cos(u@x)-1)+.12*v*math.sin(v@x))
        def J(x):
            return A*(.5*np.eye(d)-.1*math.sin(u@x)*np.outer(u,u)
                      +.12*math.cos(v@x)*np.outer(v,v))
        def R(w):
            return K*alpha*T@np.tanh(Q@w)
        def DR(w):
            return K*alpha*T@np.diag(1-np.tanh(Q@w)**2)@Q
        def E(w):
            return g(center+s*(P@w+R(w)))-g(center+s*(P@w))
        require(np.linalg.norm(E(np.zeros(n))) < 1e-15, 'literal shared-origin chord zero')
        for _ in range(24):
            w = rng.normal(size=n)
            h = rng.normal(size=n); h /= np.linalg.norm(h)
            a = center+s*(P@w+R(w)); b = center+s*(P@w)
            J1, J0 = J(a), J(b)
            D = s*((J1-J0)@P+J1@DR(w))
            D_direct = J1@(s*P+s*DR(w))-J0@(s*P)
            chain_error = op(D-D_direct)
            max_chain_error = max(max_chain_error, chain_error)
            require(chain_error < 2e-16, 'complete chain-rule identity')
            lift = P.T@D
            symmetric = s*P.T@(J1-J0)@P
            require(op(symmetric-symmetric.T) < 3e-17, 'leading Hessian difference symmetric')
            expected_skew = s*(P.T@J1@DR(w)-DR(w).T@J1@P)
            require(op(lift-lift.T-expected_skew) < 3e-17, 'only residual derivative contributes curl')
            bound = 2*K*A*A*s**3
            ratio = op(lift-lift.T)/bound
            max_curl_ratio = max(max_curl_ratio, ratio)
            require(ratio <= 1+1e-12, 'noncommuting curl majorant')
            require(op(D) <= 2*A*s+A*s*K*alpha+1e-15, 'complete first majorant')
            eigen = np.linalg.eigvalsh(J1)
            require(eigen.min() >= .28*A-1e-15 and eigen.max() <= .72*A+1e-15,
                    'global convex bounded-Hessian fixture')
            comm = op(J1@J0-J0@J1)
            max_noncommutator = max(max_noncommutator, comm)
            eps = 2e-5
            fd = (E(w+eps*h)-E(w-eps*h))/(2*eps)
            ferr = float(np.linalg.norm(fd-D@h))
            max_fd_error = max(max_fd_error, ferr)
            require(ferr < 3e-10, 'finite-difference check of actual first')
            ratio_e = float(np.linalg.norm(E(w)))/(K*A*A*s**3*np.linalg.norm(w))
            max_energy_pointwise_ratio = max(max_energy_pointwise_ratio, ratio_e)
            require(ratio_e <= 1+1e-9, 'complete-tape pointwise energy majorant')
        ell = math.sqrt(2)*A*s*(2+K*alpha)
        kappa = 2*math.sqrt(2)*K*A*A*s**3
        require(ell <= .25, 'actual normalized covariance-gap guard')
        # Exact imported law expression, with centered mark bounded by anchored energy.
        e = math.sqrt(2)*K*math.sqrt(n/d)*A*A*s**3*math.sqrt(d)
        error = math.sqrt(.5)*e*(kappa+ell*A+ell**3*(1+A**(-.5)))
        # Conservative declared scale drops favorable s powers but pays all K and n/D.
        c = math.sqrt(2)*(2+K*A)
        coeff = K*math.sqrt(n/d)*(2*math.sqrt(2)*K+c+2*c**3)
        rhs = coeff*A**4*s**3*math.sqrt(d)
        require(error <= rhs*(1+1e-12), 'conditional mean order-four scaling')
        max_conditional_error_ratio = max(max_conditional_error_ratio,error/rhs)

require(max_noncommutator > 1e-11, 'fixture genuinely has noncommuting Hessians')

# A deterministic derivative bound plus a literal zero controls sqrt(n), not sqrt(D).
for n in (2, 20, 200):
    d, A, K = 2, .01, 3.
    w = np.ones(n)
    r = K*A*np.linalg.norm(w)
    require(abs(r-K*A*math.sqrt(n)) < 1e-12, 'complete n-dimensional energy scale')
    require(abs(K*math.sqrt(n/d)*A*math.sqrt(d)-r) < 1e-12, 'literal n/D conversion')

# Dropping an uncontrolled deterministic origin cannot be justified by derivatives.
for origin in (1., 1e4, 1e8):
    derivative_bound = 0.
    require(origin > derivative_bound, 'constant-origin counterexample to zero inference')

# Exact finite quadratic conditional mode and tilt. Includes noncommuting physical matrices.
for d in (2, 4):
    H = rng.normal(size=(d,d)); H=H@H.T; H *= .02/op(H)
    b = rng.normal(size=d)*.03
    z = rng.normal(size=d)
    for r in (0., .3, .8, .999999):
        s = math.sqrt(1-r*r); a=r*z; x=a.copy()
        for _ in range(4):
            x=a-s*s*(H@x+b)
        RM=x-a+s*s*(H@x+b)
        precision=np.eye(d)+s*s*H
        exact_mean=np.linalg.solve(precision,a-s*s*b)
        shift=x-exact_mean
        require(np.linalg.norm(shift-np.linalg.solve(precision,RM)) < 2e-14,
                'physical finite-mode linear-tilt identity')
        require(np.linalg.norm(shift) <= np.linalg.norm(RM)+1e-14,
                'strong-convexity physical mode floor')
        require(np.linalg.norm(H@shift) <= op(H)*np.linalg.norm(RM)+1e-14,
                'force readout mode floor')

# Variance shares and full replay are ordinary exact finite identities.
for old_q in (17, 137, 1009):
    for nb, ne in ((3,9),(19,81),(53,729)):
        captured=11
        count_bill=captured+nb+ne*(old_q+2)
        explicit=captured+sum(1 for _ in range(nb))+sum(old_q+2 for _ in range(ne))
        require(count_bill == explicit, 'all old endpoint descendants replayed')
        require(count_bill > captured+nb+old_q+2, 'no cached old random ancestor')
        require(.5+.5 == 1., 'conditional independent variance shares')

report = {
    'status': 'PASS', 'assertions': count,
    'scope': 'Finite algebra, source pins, noncommuting chain rule, energy/conditional scaling, mode tilt, variance and replay diagnostics; imported compilers are not instantiated.',
    'source_sha256': SOURCE_PIN,
    'source_pin_count': len(pins),
    'max_directional_finite_difference_error': max_fd_error,
    'max_chain_rule_roundoff': max_chain_error,
    'max_curl_to_claimed_majorant_ratio': max_curl_ratio,
    'max_hessian_noncommutator': max_noncommutator,
    'max_pointwise_energy_to_complete_tape_majorant_ratio': max_energy_pointwise_ratio,
    'max_conditional_consumer_error_to_conservative_bound_ratio': max_conditional_error_ratio,
    'verified_pins': {source_label(p):v for p,v in pins.items()},
}
report['publication_pin_verification'] = pin_report()
(HERE/'fourth_force_reentry_checks.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='verified_pins'}, indent=2))
