#!/usr/bin/env python3
"""Independent checks via Brownian filter kernels and normalized Hermite sums.
No author's diagnostic implementation is imported or executed.
Numerics calibrate, rather than prove, the dimension-uniform claims.
"""
from pathlib import Path
import hashlib, json, math
import numpy as np
import sympy as sp
from scipy.special import roots_hermitenorm, roots_legendre


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

ROOT = Path(__file__).resolve().parent.parent
inputs = [
 'CHEAP-SHARED-ROOT-COVARIANCE-SEPARATOR.md',
 'CONDITIONAL-INNOVATION-LOCALIZES-THE-COVARIANCE-GATE.md',
 'SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md',
 '../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md',
]
pins = {source_label(ROOT/p): hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}
expected = [
 'a336ca9e1b8e9923c70c6cbb3d88706f3e393803b5a9f4f87f0c85040d414ddd',
 '186ccd57441273d5cdef1cb953bc0ba6f17fef3dc9c6dd2189a8dde91b965e63',
 'a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca',
 'b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced',
]
assert list(pins.values()) == expected

# Conditional Brownian kernels derived by interchanging the two integrations.
a,s,t = sp.symbols('a s t', positive=True)
brownian_X = sp.sqrt(2)*sp.exp(-(t-s))
kH = sp.simplify(sp.integrate(sp.exp(-t)*brownian_X,(t,s,sp.oo)))
kJ = sp.simplify(sp.integrate(t*sp.exp(-t)*brownian_X,(t,s,sp.oo)))
kS = sp.simplify(a*a*(kH/2-kJ))
kF2 = a*kH-a*a*kJ
inner = lambda f,g: sp.simplify(sp.integrate(f*g,(s,0,sp.oo)))
varS = inner(kS,kS)
cross = inner(a*kH,kS)
varF2 = sp.expand(inner(kF2,kF2))
Bj = a/2-a*a/4
K = -a*inner(a*kH,a*kH)
localized = sp.expand(Bj**2+K)
assert varS == a**4/8
assert cross == -a**3/8
assert K == -a**3/4
assert sp.expand(varF2-localized) == a**4/4
assert sp.simplify(kS+a*a*s*sp.exp(-s)/sp.sqrt(2)) == 0
# Transpose requirement is detectable on a nonnormal Lipschitz vector field.
M = sp.Matrix([[1,2],[-3,sp.Rational(1,2)]])
true_cov = M*M.T/4
wrong_square = M*M/4
assert true_cov != wrong_square
assert true_cov == (M/2)*(M/2).T

# Independent Gaussian quadrature and normalized Hermite covariance expansion.
x, w = roots_hermitenorm(256)
w = w/math.sqrt(2*math.pi)
ell = np.logaddexp(x,-x)-math.log(2)
c0 = float(w@ell)
d = c0/8
h = np.tanh(x-d)-np.tanh(x)
hprime = 1/np.cosh(x-d)**2-1/np.cosh(x)**2
hp = float(w@hprime)
coeff = []
p0,p1 = np.ones_like(x),x.copy()
coeff.extend([float(w@(h*p0)),float(w@(h*p1))])
for n in range(1,48):
    pn = (x*p1-math.sqrt(n)*p0)/math.sqrt(n+1)
    coeff.append(float(w@(h*pn)))
    p0,p1 = p1,pn
assert abs(coeff[1]-hp)<1e-11
assert hp < 0
r, wr = roots_legendre(240)
r, wr = (r+1)/2,wr/2
rhoC = np.outer(r,r)+np.outer(np.sqrt(1-r*r),np.sqrt(1-r*r))
weights2 = np.outer(wr,wr)
rho_gap_exact = (math.pi**2-4)/16
rho_gap_numeric = float(np.sum(weights2*rhoC))-0.5
assert abs(rho_gap_numeric-rho_gap_exact)<2e-8
partial_gaps={}
for cap in [1,8,16,32,48]:
    partial_gaps[str(cap)] = float(sum(coeff[n]**2*(np.sum(weights2*rhoC**n)-1/(n+1)) for n in range(1,cap+1)))
assert all(v>0 for v in partial_gaps.values())
finite_rules={}
markov_var = sum(coeff[n]**2*(1/(n+1)-1/(n+1)**2) for n in range(1,49))
for count in [2,4,8,16,32,64]:
    rq,wq=roots_legendre(count)
    rq,wq=(rq+1)/2,wq/2
    cq=np.outer(rq,rq)+np.outer(np.sqrt(1-rq*rq),np.sqrt(1-rq*rq))
    cqvar=sum(coeff[n]**2*(np.sum(np.outer(wq,wq)*cq**n)-float(wq@rq**n)**2) for n in range(1,49))
    finite_rules[str(count)]={
      'mass':float(sum(wq)), 'first_moment':float(wq@rq),
      'integrated_conditional_variance_gap_48_chaoses':float(cqvar-markov_var)}

results={
 'status':'PASS independent analytical-calibration checks; finite producer remains OPEN',
 'input_sha256':pins,
 'brownian_kernel_H0':str(kH), 'brownian_kernel_J0':str(kJ), 'brownian_kernel_S':str(kS),
 'innovation_variance':str(varS), 'mixed_Hg_S':str(cross), 'K':str(K),
 'F2_variance':str(varF2), 'localized_F2':str(localized),
 'order4_remainder':str(sp.expand(varF2-localized)),
 'nonnormal_field_correct_covariance':str(true_cov),
 'nonnormal_field_wrong_no_transpose_square':str(wrong_square),
 'mean_logcosh':c0, 'd':d, 'E_hprime':hp,
 'hermite_first_coefficient':coeff[1],
 'first_chaos_variance_gap_lower_bound':hp*hp*rho_gap_exact,
 'covariance_gap_A2_coefficient_lower_bound':hp*hp*rho_gap_exact/16,
 'common_minus_markov_correlation_integral_exact':rho_gap_exact,
 'common_minus_markov_correlation_integral_numeric':rho_gap_numeric,
 'variance_gap_hermite_partial_sums':partial_gaps,
 'finite_positive_outer_rule_diagnostics':finite_rules,
 'limitations':[
  'Numerics do not establish the full conditional covariance at any fixed endpoint.',
  'Hermite sums are truncated and reported solely as calibration.',
  'No finite VALUE source for j, K, same-carrier m3, or C2 is supplied.',
  'No browser, external source, or network action was used.'
 ]
}
out=Path(__file__).with_name('localization_independent_checks.json')
out.write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
