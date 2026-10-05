#!/usr/bin/env python3
"""Independent exact-input, algebra, finite-source and radial diagnostics.

The analytical audit proves the uniform estimates. This script does not simulate
an OU path, certify an asymptotic dimension threshold, or compile an expectation.
No author checker is imported or executed.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
import sympy as sp

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
PINS = {
    'RADIAL-INFLATION-REFUTES-CENTERED-SINGLE-HISTORY-COVARIANCE.md': 'c6024a6ba2ca4fbb293063d0746ee92fa083398a3bf4c01fb4de1908f322dfaf',
    'check_centered_radial_inflation_constants.py': '555c6a5bca6d4503bd41b04090254846b58c66c72972336b5e257da2ac495f8f',
    'centered_radial_inflation_constants.json': '0f748a70f15987be389c155c787b904edbbd3ef3424dcac1c73a488e4f7646b0',
    'SINGLE-HISTORY-EXACT-CENTERED-RESUMMED-CURRENT.md': '276d5873cc7d3f5ca8650955c08abe128cec7eda96680319f5c2a00b45118e9a',
    '../../../ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md': '4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f',
}
checks = 0

def require(condition, label):
    global checks
    if not bool(condition):
        raise AssertionError(label)
    checks += 1

for relative, expected in PINS.items():
    require(hashlib.sha256((BASE / relative).read_bytes()).hexdigest() == expected,
            'frozen input: ' + relative)

c, d = 0.5, 0.25

def eta(s):
    return math.exp(-1 / (1 - s*s/4)) if abs(s) < 2 else 0.0

Z = quad(eta, -2, 2, epsabs=2e-14, epsrel=2e-14)[0]

def bump(s):
    return eta(s) / Z

def bump_prime(s):
    return -s*bump(s)/(2*(1-s*s/4)**2) if abs(s) < 2 else 0.0

def psi(s):
    if s <= -2: return 0.0
    if s >= 2: return 1.0
    return quad(eta, -2, s, epsabs=2e-13, epsrel=2e-13)[0] / Z

require(Z >= 2*math.exp(-4/3), 'elementary normalization lower bound')
require(bump(0) <= math.exp(1/3)/2 < 1, 'global bump maximum is below one')
require(abs(psi(0)-0.5) < 1e-13, 'even bump normalizes to symmetric step')
for s in np.linspace(0.02, 1.98, 43):
    require(bump_prime(s) < 0, 'strict decrease of positive bump')
    require(abs(bump(s)-bump(-s)) < 1e-14, 'even bump')

# Reconstruct the admitted positive dyadic rule from its analytical contract.
# The certificate is all-degree; sampling powers only checks implementation.
def positive_rule(delta):
    panels = math.ceil(math.log2(4/delta))
    order = math.ceil(math.log(16/delta)/math.log(4))
    roots, weights = np.polynomial.legendre.leggauss(order)
    sites, masses = [], []
    for j in range(panels):
        left = 2.0**(-j-1)
        t = 1.5*left + 0.5*left*roots
        sites.extend(1-t)
        masses.extend(0.5*left*weights)
    sites.append(1-2.0**(-panels-1))
    masses.append(2.0**(-panels))
    sites, masses = np.array(sites), np.array(masses)
    certificate = 8*4.0**(-order)+2.0**(1-panels)
    require(certificate <= delta*(1+1e-14), 'all-degree clock certificate')
    require(np.all((sites > 0) & (sites < 1)), 'finite interior clocks')
    require(np.all(masses > 0), 'positive clock weights')
    require(abs(masses.sum()-1) < 3e-15, 'constant exactness')
    require(abs(masses@sites-0.5) < 3e-15, 'linear exactness')
    for degree in [0, 1, 2, 3, 8, 16, 37, 128, 511, 2048, 10000]:
        require(abs(masses@(sites**degree)-1/(degree+1)) <= certificate+3e-15,
                'sample Hermite multiplier')
    return sites, masses, certificate

rng = np.random.default_rng(2026100501)
finite_rules = []
largest_decomposition_error = 0.0
largest_terminal_error = 0.0
largest_proxy_ratio = 0.0
max_tensor_norm = 0.0
for D in [256, 4096, 65536]:
    A = D**(-0.25)
    R = math.sqrt(D)
    k = 1-c*A/2
    R0 = k*R
    require(R0-2 > R/2, 'global derivative support away from origin')
    tau, v, certificate = positive_rule(A*A)
    noise = np.sqrt(1-tau*tau)
    beta = float(v@noise)
    require(beta >= 2/3-A*A-1e-14, 'Hermite second-moment lower bound on beta')
    require(beta >= 7/12, 'uniform admitted beta interval')
    require(beta*beta-0.25 >= 13/144, 'uniform variance separation')
    require(abs(v@(tau*tau)-1/3) <= A*A+1e-14, 'second Hermite multiplier')
    wrong_independent_variance = float(np.sum(v*v*noise*noise))
    require(beta*beta > wrong_independent_variance, 'shared root cross terms matter')
    # Both eigenvalues of Dg, and the exact HS-to-vector D2f norm.
    radii = [0, R0-3, R0-2, R0+3, 2*R, 100*R]
    radii += list(R0+np.linspace(-1.999, 1.999, 67))
    for r in radii:
        b = bump(r-R0)
        p = psi(r-R0)
        radial = c*A+d*A*b
        tangent = c*A+d*A*p/r if r > 0 else c*A
        require(0 <= radial <= A and 0 <= tangent <= A, 'sample global Hessian sandwich')
        ar = b/r-p/r**2 if r > 0 else 0
        tensor_norm = max(math.sqrt(bump_prime(r-R0)**2+(D-1)*ar**2), math.sqrt(2)*abs(ar))
        max_tensor_norm = max(max_tensor_norm, tensor_norm)
        require(tensor_norm < 3, 'dimension-uniform radial tensor norm diagnostic')
    if D <= 4096:
        def f(rows):
            rows = np.asarray(rows)
            rr = np.linalg.norm(rows, axis=-1)
            vals = np.array([psi(float(s-R0)) for s in rr.reshape(-1)]).reshape(rr.shape)
            scales = np.divide(vals, rr, out=np.zeros_like(vals), where=rr > 0)
            return rows*scales[..., None]
        def g(rows):
            return c*A*np.asarray(rows)+d*A*f(rows)
        for _ in range(4):
            x, H = rng.standard_normal((2, D))
            sites = tau[:, None]*x+noise[:, None]*H
            K = v@f(sites)
            actual = v@g(sites)
            decomposition = c*A*(x/2+beta*H)+d*A*K
            err = float(np.linalg.norm(actual-decomposition))
            largest_decomposition_error = max(largest_decomposition_error, err)
            require(err < 1e-12, 'literal SAME-g finite source decomposition')
            require(np.linalg.norm(K) <= 1+1e-14, 'pathwise bounded nonlinear source')
            y = k*x-c*A*beta*H
            terminal_err = float(np.linalg.norm(g(x-actual)-g(y-d*A*K)))
            largest_terminal_error = max(largest_terminal_error, terminal_err)
            require(terminal_err < 1e-12, 'literal terminal same-g evaluation')
            proxy_ratio = float(np.linalg.norm(g(y-d*A*K)-g(y))/(d*A*A))
            largest_proxy_ratio = max(largest_proxy_ratio, proxy_ratio)
            require(proxy_ratio <= 1+1e-12, 'pointwise full shifted VALUE proxy bound')
    finite_rules.append({'D': D, 'A': A, 'nodes': len(v), 'certificate': certificate,
                         'beta': beta, 'shared_variance': beta*beta,
                         'incorrect_independent_variance': wrong_independent_variance})

# Check the radial Hessian contraction with directional differences at D=256.
D = 256
A, R = D**(-0.25), math.sqrt(D)
R0 = (1-c*A/2)*R
for shell_coordinate in [-1.2, -0.3, 0.4, 1.5, 3.0]:
    z = np.zeros(D)
    z[0] = R0+shell_coordinate
    direction = rng.standard_normal(D)
    direction /= np.linalg.norm(direction)
    n = z/np.linalg.norm(z)
    r = np.linalg.norm(z)
    ar = bump(shell_coordinate)/r-psi(shell_coordinate)/r**2
    nv = n@direction
    exact = bump_prime(shell_coordinate)*n*nv*nv+ar*(n*(1-nv*nv)+2*nv*(direction-n*nv))
    def local_f(xx):
        rr = np.linalg.norm(xx)
        return psi(rr-R0)*xx/rr
    h = 0.003
    measured = (local_f(z+h*direction)-2*local_f(z)+local_f(z-h*direction))/(h*h)
    require(np.linalg.norm(exact-measured) < 2e-7, 'radial D2f rank-one contraction')

# Independent normalization: Brownian Fubini coefficient gives variance 1/4.
t = sp.symbols('t', nonnegative=True)
require(sp.integrate(sp.exp(-2*t), (t, 0, sp.oo)) == sp.Rational(1, 2), 'true linear mean coefficient')
require(sp.integrate(sp.exp(-2*t)/2, (t, 0, sp.oo)) == sp.Rational(1, 4), 'true conditional linear variance')
require(sp.integrate(sp.Symbol('r'), (sp.Symbol('r'), 0, 1)) == sp.Rational(1, 2), 'R1 linear witness eigenvalue')
a = sp.symbols('a', positive=True)
require(sp.simplify(a**4*sp.sqrt(a**-4)-a*a) == 0, 'fourth-order allowance')
require(sp.simplify(a**2*sp.sqrt(a**-4)-1) == 0, 'order-one radial inflation')
require(sp.simplify(a*sp.sqrt(a**-4)-1/a) == 0, 'diverging Euclidean displacement scale')
require(sp.simplify(a/(a**4*sp.sqrt(a**-4))-1/a) == 0, 'power separation ratio')
require(sp.simplify(a*a*(a*a)-a**4) == 0, 'improved mean-quadrature floor')

# Strict sign is analytical in the report. These independent integral forms
# diagnose its size and the subtraction sign without sampling shell normals.
def m(t):
    return quad(lambda u: bump(u)*math.exp(-(u-t)**2)/math.sqrt(math.pi),
                -2, 2, epsabs=2e-14, epsrel=2e-14)[0]

def mean_step(t):
    return quad(lambda u: bump(u)*ndtr(math.sqrt(2)*(t-u)),
                -2, 2, epsabs=2e-14, epsrel=2e-14)[0]

hH = c*c/8
m0 = m(0)
witnesses = []
previous = m0
for shift in np.linspace(0.002, c*c/2, 39):
    current = m(float(shift))
    require(current < previous, 'strict Gaussian-smoothed bump monotonicity diagnostic')
    previous = current
for beta in [7/12, 2/3, math.pi/4, 1.0]:
    h = c*c*beta*beta/2
    integral = quad(lambda t: m0-m(t), hH, h, epsabs=5e-16, epsrel=5e-12)[0]
    direct = mean_step(hH)-mean_step(h)+(h-hH)*m0
    require(integral > 0, 'positive limiting resolvent witness')
    require(abs(integral-direct) < 1e-13, 'independent equivalent witness integrals')
    witnesses.append({'beta': beta, 'hH': hH, 'h_beta': h, 'mean_q': integral,
                      'direct_mean_q': direct, 'leading_resolvent_coefficient': d*integral/2})
kappa = witnesses[0]['mean_q']
require(1e-6 < kappa < 3e-6, 'nonzero uniform witness scale')

# Exact marginal radial sampling, without allocating D-dimensional arrays and
# without any path approximation. Used only to inspect the O(A) geometry.
geometry = []
N = 50000
for D in [2**12, 2**20, 2**28, 2**36]:
    A, R = D**(-0.25), math.sqrt(D)
    k = 1-c*A/2
    radii = np.sqrt(rng.chisquare(D, N))
    s = radii-R
    normal = rng.standard_normal(N)
    transverse_sq = rng.chisquare(D-1, N)
    rho = k*radii
    for sigma in [0.0, 0.5, math.pi/4, 1.0]:
        scale = c*A*sigma
        parallel = rho-scale*normal
        norm_y = np.sqrt(parallel*parallel+scale*scale*transverse_sq)
        # Stable difference of radii even when sqrt(D) is very large.
        radial_change = (-2*rho*scale*normal+scale*scale*(normal*normal+transverse_sq))/(norm_y+rho)
        relative_radius = k*s+radial_change
        target = s+c*c*sigma*sigma/2
        residual = math.sqrt(float(np.mean((relative_radius-target)**2)))
        direction_l2 = math.sqrt(float(np.mean(2*np.maximum(0, 1-parallel/norm_y))))
        require(residual/A < 2.0, 'radial residual O(A) diagnostic')
        require(direction_l2/A < 1.0, 'direction residual O(A) diagnostic')
        geometry.append({'D': D, 'A': A, 'sigma': sigma, 'radial_residual_l2': residual,
                         'radial_residual_over_A': residual/A,
                         'direction_l2_over_A': direction_l2/A,
                         'Euclidean_noise_rms': scale*R,
                         'limiting_radial_inflation': c*c*sigma*sigma/2})

result = {
    'status': 'PASS: independent diagnostics; analytical verdict in accompanying report',
    'assertions': checks,
    'input_pins': PINS,
    'bump_normalization': Z,
    'bump_maximum': bump(0),
    'max_sampled_HS_to_vector_D2f_norm': max_tensor_norm,
    'finite_rules': finite_rules,
    'largest_literal_source_decomposition_error': largest_decomposition_error,
    'largest_literal_terminal_error': largest_terminal_error,
    'largest_pointwise_proxy_error_over_dA2': largest_proxy_ratio,
    'uniform_limiting_kappa': kappa,
    'lower_bound_c_star': d*kappa/8,
    'witnesses': witnesses,
    'exact_radial_geometry_diagnostics': geometry,
    'limits': [
        'The checks do not prove the infinite-dimensional Bessel inequality or uniform Lp bounds by computation.',
        'Radial sampling uses exact marginal laws, not a true OU-history simulation.',
        'No finite-D onset threshold for the asymptotic positive witness is certified.',
        'The true same-g source and all correlations are justified analytically in the audit.',
        'No generic shifted-current consumer, expectation compiler, m3 closure, or endpoint-law theorem is supplied.'
    ]
}
(HERE/'center_inflation_independent_checks.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({key: result[key] for key in ['status', 'assertions', 'uniform_limiting_kappa',
      'lower_bound_c_star', 'largest_literal_source_decomposition_error',
      'largest_literal_terminal_error', 'largest_pointwise_proxy_error_over_dA2']}, indent=2))
