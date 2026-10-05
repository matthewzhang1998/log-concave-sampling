"""Scalar diagnostics for the analytical centered radial-inflation separator.

This does not simulate or approximate the true OU source and is not a proof of
its asymptotics. The accompanying Markdown artifact contains the source proof.
"""
from pathlib import Path
import json
import math
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
c, d = 0.5, 0.25
beta_min = 7.0 / 12.0


def eta(u):
    if abs(u) >= 2:
        return 0.0
    return math.exp(-1 / (1 - u * u / 4))


Z, _ = quad(eta, -2, 2, epsabs=1e-13, epsrel=1e-13)


def bump(u):
    return eta(u) / Z


def m(t):
    val, _ = quad(lambda u: bump(u) * math.exp(-(u-t)**2) / math.sqrt(math.pi),
                  -2, 2, epsabs=1e-13, epsrel=1e-13)
    return val


hH = c*c/8
hmin = c*c*beta_min*beta_min/2
h_continuous = c*c*(math.pi/4)**2/2
m0 = m(0)
kappa, err = quad(lambda t: m0-m(t), hH, hmin, epsabs=1e-14, epsrel=1e-10)
kappa_continuous, err_continuous = quad(lambda t: m0-m(t), hH, h_continuous,
                                       epsabs=1e-14, epsrel=1e-10)
assert bump(0) < 1
assert beta_min*beta_min-0.25 >= 13/144 - 1e-14
assert hmin > hH > 0
assert m0 > m(hH) > m(hmin) > m(h_continuous)
assert kappa > 1000*err > 0
assert kappa_continuous > kappa

rows = []
for D in (256, 4096, 65536, 1048576):
    A = D**(-0.25)
    R = math.sqrt(D)
    k = 1-c*A/2
    R0 = k*R
    assert R0-2 >= R/2
    assert c+d < 1
    assert c+d/(R0-2) < 1
    assert abs(A*A*R-1) < 1e-14
    assert abs(A**4*R-A*A) < 1e-14
    rows.append(dict(D=D,A=A,R=R,R0=R0,radial_inflation_scale=A*A*R,
                     fourth_order_allowance=A**4*R,
                     asymptotic_counterexample_scale=A,
                     scale_ratio=1/A,
                     exact_gaussian_proxy_difference_bound=2*d*A*A))

out = dict(status='Scalar diagnostic PASS; analytical source proof and independent audit remain distinct',
           c=c,d=d,bump_normalization=Z,bump_supremum=bump(0),
           minimum_beta=beta_min,minimum_covariance_gap=beta_min*beta_min-0.25,
           h_H=hH,h_Q_min=hmin,h_Q_continuous_common_root=h_continuous,
           m_0=m0,m_hH=m(hH),m_hQ_min=m(hmin),
           uniform_limiting_shell_witness=kappa,
           uniform_kappa_quadrature_error_estimate=err,
           continuous_clock_limiting_shell_witness=kappa_continuous,
           proven_lower_constant_c_star=d*kappa/8,
           dimensions=rows,
           exclusions=['No OU path grid','No numerical proof of the source estimates',
                       'No generic trace/Hodge estimate','No general finite mean consumer'])
path = ROOT/'centered_radial_inflation_constants.json'
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
