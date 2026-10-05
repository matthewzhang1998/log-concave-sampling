#!/usr/bin/env python3
"""Independent diagnostics for the exact centered actual-P3 counterexample.

The accompanying report contains the analytical audit. This checker imports no
author checker, simulates no OU history, and does not certify a finite-D lower
bound or provide a general Gaussianization/mean compiler.
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
    'RADIAL-INFLATION-REFUTES-CENTERED-P3-COVARIANCE-CORRECTION.md': '810128a4b47c568273f05442a510498ac5f77eedf1cb864366368ebb4bd356af',
    'RADIAL-INFLATION-REFUTES-CENTERED-SINGLE-HISTORY-COVARIANCE.md': 'c6024a6ba2ca4fbb293063d0746ee92fa083398a3bf4c01fb4de1908f322dfaf',
    'P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md': '87b61a6a4ca03f81f7d910cdaa83b3783cd4792b3a428b89c8f0f629062e4467',
    'resummed-linear-backbone/BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md': '798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3',
    '../../../ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md': '4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f',
}
checks = 0

def require(condition, description):
    global checks
    if not bool(condition):
        raise AssertionError(description)
    checks += 1

for relative, digest in PINS.items():
    require(hashlib.sha256((BASE / relative).read_bytes()).hexdigest() == digest,
            'exact source pin: ' + relative)

# Derive the continuous conditional coefficients two independent ways.
a, b = sp.symbols('a b', positive=True)
K = (1 / (a+b)) * (1 / (a+1) + 1 / (b+1))
at_one = lambda f: sp.simplify(f.subs({a: 1, b: 1}))
unconditioned = [at_one(K), at_one(-sp.diff(K, b)),
                 at_one(sp.diff(K, a, b))]
require(unconditioned == [sp.Rational(1, 2), sp.Rational(3, 8), sp.Rational(3, 8)],
        'stationary Laplace covariance derivatives')
conditional = [unconditioned[0]-sp.Rational(1, 4),
               unconditioned[1]-sp.Rational(1, 8),
               unconditioned[2]-sp.Rational(1, 16)]
require(conditional == [sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(5, 16)],
        'subtract endpoint regression products')
s, t = sp.symbols('s t', positive=True)
u = sp.symbols('u', nonnegative=True)
kernel1 = sp.sqrt(2) * sp.integrate(sp.exp(-t)*sp.exp(-(t-s)), (t, s, sp.oo))
kernel2 = sp.sqrt(2) * sp.integrate(t*sp.exp(-t)*sp.exp(-(t-s)), (t, s, sp.oo))
require(sp.simplify(kernel1-sp.exp(-s)/sp.sqrt(2)) == 0, 'Brownian H1 coefficient')
require(sp.simplify(kernel2-(s+sp.Rational(1,2))*sp.exp(-s)/sp.sqrt(2)) == 0,
        'Brownian H2 coefficient')
require([sp.integrate(expr, (s, 0, sp.oo)) for expr in
         [kernel1**2, kernel1*kernel2, kernel2**2]] == conditional,
        'conditional covariance from genuine common Brownian ancestry')
require(sp.integrate(sp.exp(-t)*sp.exp(-(s-t)), (t, 0, s)) == s*sp.exp(-s),
        'future-ancestry convolution gives s exp(-s)')

lam, bm, bi, A = sp.symbols('lambda beta_mid beta_in A', positive=True)
a2 = lam/2-lam**2/4
vH = lam**2/4-lam**3/2+5*lam**4/16
vQ = lam**2*(1-lam/2)**2*bm**2+lam**4*bi**2
require(sp.simplify(sp.integrate((lam*kernel1-lam**2*kernel2)**2,
                               (s, 0, sp.oo))-vH) == 0,
        'actual F2 conditional Gaussian variance')
require(sp.discriminant(4-8*lam+5*lam**2, lam) < 0, 'backbone variance positive')
x, H, J = sp.symbols('x H J', real=True)
finite_backbone = lam*(1-lam/2)*(x/2+bm*H)-lam**2*bi*J
require(sp.simplify(finite_backbone.subs({H:0, J:0})-a2*x) == 0,
        'literal finite coherent mean coefficient')
require(sp.simplify(sp.diff(finite_backbone,H)**2+
                    sp.diff(finite_backbone,J)**2-vQ) == 0,
        'independent level roots give literal finite variance')
require(sp.expand(vQ-(sp.diff(finite_backbone,H)+sp.diff(finite_backbone,J))**2) != 0,
        'conflating the two level roots changes the exact variance')
require(sp.limit(((vH-vQ)-lam**2*(sp.Rational(1,4)-bm**2))/lam**2, lam, 0) == 0,
        'C2-minus-C1 isotropic difference begins at lambda cubed')
require(sp.simplify(A**4*sp.sqrt(A**-4)-A**2) == 0, 'fourth-order allowance is A squared')
require(sp.simplify(A**2*sp.sqrt(A**-4)-1) == 0, 'radial inflation parameter is one')
require(sp.simplify(A**3*A*sp.sqrt(A**-4)-A**2) == 0,
        'outer delta=A cubed floor is A squared')
require(sp.simplify(A*A-A**2) == 0, 'counterterm delta=A floor is A squared')

c, d = 0.5, 0.25
def eta(z):
    return math.exp(-1/(1-z*z/4)) if abs(z) < 2 else 0.0
Z = quad(eta, -2, 2, epsabs=2e-14, epsrel=2e-14)[0]
def bump(z):
    return eta(z)/Z
def bump_derivative(z):
    return -z*bump(z)/(2*(1-z*z/4)**2) if abs(z)<2 else 0.0
def psi(z):
    if z <= -2: return 0.0
    if z >= 2: return 1.0
    return quad(eta, -2, z, epsabs=3e-13, epsrel=3e-13)[0]/Z
require(bump(0) <= math.exp(1/3)/2 < 1, 'analytic normalization proves bump supremum below one')
require(abs(psi(0)-0.5) < 2e-14, 'step symmetry')
require(sp.Rational(7,12)**2-sp.Rational(1,4) == sp.Rational(13,144),
        'exact uniform variance gap')

def positive_rule(delta):
    panels = math.ceil(math.log2(4/delta))
    order = math.ceil(math.log(16/delta)/math.log(4))
    roots, weights = np.polynomial.legendre.leggauss(order)
    nodes, masses = [], []
    for j in range(panels):
        left = 2.0**(-j-1)
        nodes.extend(1-(1.5*left+0.5*left*roots))
        masses.extend(0.5*left*weights)
    nodes.append(1-2.0**(-panels-1))
    masses.append(2.0**(-panels))
    nodes, masses = np.array(nodes), np.array(masses)
    certificate = 8*4.0**(-order)+2.0**(1-panels)
    require(certificate <= delta*(1+1e-14), 'all-degree positive clock certificate')
    require(np.all((nodes>0)&(nodes<1)) and np.all(masses>0), 'positive interior rule')
    require(abs(masses.sum()-1)<3e-15, 'exact mass diagnostic')
    require(abs(masses@nodes-0.5)<3e-15, 'exact first moment diagnostic')
    for degree in [0,1,2,3,7,19,64,255,1024,8192]:
        require(abs(masses@(nodes**degree)-1/(degree+1))<certificate+3e-15,
                'sample Hermite multiplier')
    return nodes, masses, certificate

rng = np.random.default_rng(202610050133)
finite_rows = []
max_literal_error = 0.0
max_proxy_ratio = 0.0
max_tensor_norm = 0.0
for D in [256, 4096]:
    aa, R = D**(-0.25), math.sqrt(D)
    ll, eps = c*aa, d*aa
    aa2 = ll/2-ll*ll/4
    kk, R0 = 1-aa2, (1-aa2)*R
    B = eps*(1+ll)
    require(R0-2 > R/2, 'second-substitution shell stays away from zero')
    tau, v, mid_cert = positive_rule(aa**2)
    sigma, w, in_cert = positive_rule(aa**3/2)
    _, _, out_cert = positive_rule(aa**3)
    dm, ei = np.sqrt(1-tau*tau), np.sqrt(1-sigma*sigma)
    beta_m, beta_i = float(v@dm), float(w@ei)
    require(beta_m >= 2/3-aa*aa and beta_m >= 7/12, 'actual positive-rule beta gap')
    require(beta_m*beta_m > float((v*v)@(dm*dm)), 'shared middle root cross terms retained')
    vvH = ll*ll/4-ll**3/2+5*ll**4/16
    vvQ = ll*ll*(1-ll/2)**2*beta_m*beta_m+ll**4*beta_i*beta_i
    require(vvH > 0 and vvQ > vvH, 'sample literal variance separation')
    for shell_s in np.linspace(-1.999, 1.999, 65):
        rr = R0+shell_s
        pp, bp = psi(shell_s), bump(shell_s)
        require(0 <= ll+d*aa*bp <= aa and 0 <= ll+d*aa*pp/rr <= aa,
                'sample Hessian sandwich')
        ar = bp/rr-pp/rr**2
        norm = max(math.sqrt(bump_derivative(shell_s)**2+(D-1)*ar**2),
                   math.sqrt(2)*abs(ar))
        max_tensor_norm = max(max_tensor_norm, norm)
        require(norm < 3, 'sample exact HS-to-vector D2f norm')
    def f(rows):
        rows = np.asarray(rows)
        radii = np.linalg.norm(rows, axis=-1)
        p = np.array([psi(float(r-R0)) for r in radii.reshape(-1)]).reshape(radii.shape)
        ratio = np.divide(p,radii,out=np.zeros_like(radii),where=radii>0)
        return rows*ratio[...,None]
    def g(rows):
        return ll*np.asarray(rows)+d*aa*f(rows)
    for trial in range(2):
        xx, hh, jj = rng.standard_normal((3,D))
        yy = tau[:,None]*xx+dm[:,None]*hh
        literal_L, remainder_L = [], []
        for y in yy:
            zz = sigma[:,None]*y+ei[:,None]*jj
            literal_L.append(w@g(zz))
            remainder_L.append(w@(d*aa*f(zz)))
        literal_L, remainder_L = np.array(literal_L), np.array(remainder_L)
        require(np.max(np.linalg.norm(remainder_L,axis=1)) <= eps+1e-13,
                'inner remainder is pathwise bounded')
        reconstructed_L = ll*(yy/2+beta_i*jj)+remainder_L
        require(np.linalg.norm(literal_L-reconstructed_L) < 2e-12,
                'literal same-J inner decomposition')
        literal_I2 = v@g(yy-literal_L)
        remainder_Q = -ll*(v@remainder_L)+v@(d*aa*f(yy-literal_L))
        gaussian_part = aa2*xx+ll*(1-ll/2)*beta_m*hh-ll*ll*beta_i*jj
        residual = float(np.linalg.norm(literal_I2-gaussian_part-remainder_Q))
        max_literal_error = max(max_literal_error,residual)
        require(residual < 3e-12, 'literal full two-level SAME-g identity')
        require(np.linalg.norm(remainder_Q) <= B+1e-13, 'full nonlinear remainder bound')
        error = float(np.linalg.norm(g(xx-literal_I2)-g(xx-gaussian_part)))
        ratio = error/(aa*B)
        max_proxy_ratio = max(max_proxy_ratio,ratio)
        require(ratio <= 1+1e-11, 'pointwise same-root Gaussian response comparison')
    finite_rows.append(dict(D=D,A=aa,R0=R0,n_mid=len(v),n_in=len(w),
                            beta_mid=beta_m,beta_in=beta_i,bH_squared=vvH,bQ_squared=vvQ,
                            mid_certificate=mid_cert,inner_certificate=in_cert,
                            outer_certificate=out_cert,allowed_remainder=B,
                            gaussian_response_error_bound=aa*B))

# Scalar witness: verify both the integral-of-m and direct shifted-step forms.
def m(t):
    return quad(lambda z: bump(z)*math.exp(-(z-t)**2)/math.sqrt(math.pi),
                -2,2,epsabs=2e-14,epsrel=2e-14)[0]
def mean_step(t):
    return quad(lambda z: bump(z)*ndtr(math.sqrt(2)*(t-z)),
                -2,2,epsabs=2e-14,epsrel=2e-14)[0]
m0, hH = m(0), c*c/8
witness_rows = []
previous = m0
for t in np.linspace(0.002,c*c/2,40):
    now = m(float(t))
    require(now < previous, 'strict smoothed-bump monotonicity diagnostic')
    previous = now
for beta in [7/12,2/3,math.pi/4,1.0]:
    hq = c*c*beta*beta/2
    kappa,error = quad(lambda t:m0-m(t),hH,hq,epsabs=3e-16,epsrel=3e-12)
    direct = mean_step(hH)-mean_step(hq)+(hq-hH)*m0
    require(kappa > 0 and abs(direct-kappa)<2e-13,
            'positive witness and correction sign via independent formulas')
    witness_rows.append(dict(beta=beta,h_H=hH,h_Q=hq,mean_q=kappa,
                              direct_mean_q=direct,quadrature_error_estimate=error))
kappa = witness_rows[0]['mean_q']
require(abs(kappa-1.8145100151e-6)<1e-15, 'reported numerical uniform witness')

# Exact marginal radial geometry, without allocating full high-D vectors.
# These observations support but do not prove the analytical uniform O(A) bound.
geometry = []
N = 30000
for D in [2**12,2**20,2**28,2**36]:
    aa, RR = D**(-0.25), math.sqrt(D)
    ll = c*aa
    kk = 1-ll/2+ll*ll/4
    radii = np.sqrt(rng.chisquare(D,N))
    sd = radii-RR
    normal = rng.standard_normal(N)
    transverse = rng.chisquare(D-1,N)
    rho = kk*radii
    for beta_m,beta_i,kind in [(0.0,0.0,'true'),(7/12,0.4,'finite_min'),
                               (math.pi/4,math.pi/4,'finite_limit'),(1.0,1.0,'finite_max')]:
        if kind == 'true':
            variance=ll*ll/4-ll**3/2+5*ll**4/16
            shift=c*c/8
        else:
            variance=ll*ll*(1-ll/2)**2*beta_m*beta_m+ll**4*beta_i*beta_i
            shift=c*c*beta_m*beta_m/2
        scale=math.sqrt(variance)
        parallel=rho-scale*normal
        length=np.sqrt(parallel*parallel+variance*transverse)
        inflation=(-2*rho*scale*normal+variance*(normal*normal+transverse))/(length+rho)
        relative_radius=kk*sd+inflation
        residual=math.sqrt(float(np.mean((relative_radius-(sd+shift))**2)))
        direction=math.sqrt(float(np.mean(2*np.maximum(0,1-parallel/length))))
        require(residual/aa < 2, 'actual backbone O(A) radial residual diagnostic')
        require(direction/aa < 1, 'actual backbone O(A) direction residual diagnostic')
        geometry.append(dict(D=D,A=aa,branch=kind,variance=variance,
                             radial_error_over_A=residual/aa,direction_error_over_A=direction/aa,
                             limiting_inflation=shift,noise_rms=scale*RR))

result=dict(status='PASS: independent diagnostics; analytical verdict is in REPORT.md',
            assertions=checks,input_pins=PINS,
            exact_conditional_covariances=[str(v) for v in conditional],
            b_H_squared=str(vH),b_Q_squared=str(vQ),bump_normalization=Z,
            bump_supremum=bump(0),finite_actual_source_checks=finite_rows,
            largest_literal_decomposition_error=max_literal_error,
            largest_terminal_proxy_error_over_A_epsilon_1_plus_lambda=max_proxy_ratio,
            largest_sampled_HS_to_vector_D2f_norm=max_tensor_norm,
            uniform_limiting_kappa=kappa,c_star=d*kappa/8,
            limiting_witnesses=witness_rows,radial_geometry_diagnostics=geometry,
            limitations=['No OU path simulation or path-grid construction',
                         'No finite-D threshold or numerical lower-bound proof',
                         'No numerical certification of the imported OWN-mean compiler',
                         'No generic resummed-current impossibility claim'])
(HERE/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
