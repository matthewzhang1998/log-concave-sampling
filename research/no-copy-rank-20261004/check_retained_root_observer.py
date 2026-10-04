"""Finite diagnostic for the analytic retained-root obstruction.
The proof is in the adjacent note. Gauss-Hermite quadrature is diagnostic only.
"""
import json, math
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm

nodes, weights = roots_hermitenorm(160)
weights /= math.sqrt(2 * math.pi)
q = np.cos(nodes / math.sqrt(2))
m = math.exp(-.25)
vq = .5 * (1 - math.exp(-.5)) ** 2
vs = .5 * (1 - math.exp(-2))
assert abs(float(weights @ q)-m) < 2e-14
assert abs(float(weights @ ((q-m)**2))-vq) < 2e-14
checks = []
for g in (.25, .5, 1.):
    for A in (.5, .25, .125, .0625, .03125):
        kappa = A ** (1+g)
        lam = .5 * A * kappa * math.exp(-.25)
        e = .5 * math.sqrt(A*A*vs + kappa*kappa)
        # subtract 1 analytically so the centered first factor does not amplify
        # harmless quadrature mean error at small covariance amplitude.
        joint_gap = math.exp(-.5)/4 * float(weights @ ((q-m)*np.expm1(lam*q/2)))
        joint_lb = lam*vq*math.exp(-(1+lam)/2)/8
        z = lam*(q-m)/2
        marginal_gap = math.exp(-.5+lam*m/2) * float(weights @ (np.expm1(z)-z))
        marginal_lb = lam*lam*vq*math.exp(-(1+lam)/2)/8
        marginal_ub = lam*lam*vq*math.exp(-(1-lam)/2)/8
        same_root_w2 = math.sqrt(float(weights @ (
            lam*lam*(q-m)**2 /
            (np.sqrt(1-lam*q)+math.sqrt(1-lam*m))**2)))
        conditional_lb = lam*math.sqrt(vq)/(2*math.sqrt(1+lam))
        conditional_ub = lam*math.sqrt(vq)/(2*math.sqrt(1-lam))
        marginal_w2_ub = 3*lam*lam*vq
        assert joint_gap >= joint_lb * (1-1e-9)
        assert marginal_lb * (1-1e-7) <= marginal_gap <= marginal_ub * (1+1e-7)
        assert conditional_lb <= same_root_w2 <= conditional_ub
        if lam <= .01:
            assert joint_gap > marginal_w2_ub
        checks.append(dict(A=A,g=g,kappa=kappa,actual_energy=e,lambda_=lam,
            joint_observer_gap=joint_gap,joint_lower_bound=joint_lb,
            marginal_cos_gap=marginal_gap,marginal_cos_lower=marginal_lb,
            marginal_cos_upper=marginal_ub,same_root_w2=same_root_w2,
            marginal_w2_upper=marginal_w2_ub,
            joint_gap_over_kappa_e=joint_gap/(kappa*e),
            marginal_upper_over_A_kappa_squared_e=marginal_w2_ub/(A*kappa*kappa*e)))
# Weak actual-energy rescaling at the requested exponent.
weak_fixtures = []
for g in (.25, .5, 1.):
    for A in (.5, .25, .125, .0625):
        eta=A**2.9
        kappa=A**(1+g)
        e0=.5*math.sqrt(A*A*vs+kappa*kappa)
        e=eta*e0
        lam=eta*.5*A*kappa*math.exp(-.25)
        joint_gap=math.exp(-.5)/4*float(weights@((q-m)*np.expm1(lam*q/2)))
        joint_lb=lam*vq*math.exp(-(1+lam)/2)/8
        marginal_ub=3*lam*lam*vq
        assert joint_gap >= joint_lb*(1-1e-9)
        weak_fixtures.append(dict(A=A,g=g,eta=eta,actual_energy=e,
            energy_over_A_3p9=e/A**3.9,joint_gap=joint_gap,
            joint_gap_over_kappa_energy=joint_gap/(kappa*e),
            joint_gap_over_A_4p9_plus_g=joint_gap/A**(4.9+g),
            marginal_W2_upper=marginal_ub,
            marginal_upper_over_eta_A_kappa2_energy=marginal_ub/(eta*A*kappa*kappa*e)))
# Exact finite positive multiroot clock identity, one shared root.
cs = np.sqrt(np.array([.25, .5, .75]))
ws = np.array([.2,.3,.5])
Q = sum(w*math.exp(-(1-c*c)/2)*np.cos(c*nodes) for w,c in zip(ws,cs))
exact_var = math.exp(-1)*sum(ws[i]*ws[j]*(math.cosh(cs[i]*cs[j])-1)
    for i in range(3) for j in range(3))
assert abs(float(weights@Q)-math.exp(-.5)) < 2e-14
assert abs(float(weights@((Q-math.exp(-.5))**2))-exact_var) < 2e-14
out = dict(status='PASS',num_nominal_parameter_fixtures=len(checks),num_weak_energy_parameter_fixtures=len(weak_fixtures),num_parameter_fixtures=len(checks)+len(weak_fixtures),V_q=vq,
    joint_uniform_lower_coefficient=vq*math.exp(-9/16)/8,
    fixed_observer_sup_bound=(1+m)/4,
    fixed_observer_lipschitz_bound=math.sqrt(1/32+(1+m)**2/16),
    positive_clock_variance=exact_var,fixtures=checks,weak_energy_fixtures=weak_fixtures,
    scope='Diagnostic only; exact proof and source contract in adjacent Markdown.')
path=Path(__file__).with_name('retained_root_observer_checks.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('fixtures','weak_energy_fixtures')},indent=2))
