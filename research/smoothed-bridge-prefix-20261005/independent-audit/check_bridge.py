"""Finite scalar/matrix diagnostics only; this does not execute native services."""
import hashlib
import json
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parent
draft = root.parent / "SMOOTHED-BRIDGE-PREFIX.md"
rng = np.random.default_rng(510052026)
out = {"scope": "Scalar bridge identities, analytic-domain checks, finite matrix contractions and carrier-row norms only. No native compiler execution.",
       "draft_sha256": hashlib.sha256(draft.read_bytes()).hexdigest()}

max_mean_error = 0.0
max_variance_error = 0.0
max_complex_identity_error = 0.0
max_wedge_excess = 0.0
max_ellipse_excess = 0.0
for w in [0.5, 0.2, 0.05, 0.001, 1e-6]:
    q = 1-w
    delta = -np.log(q)
    den = -np.expm1(-2*delta)
    for f in np.linspace(0.001, 0.999, 501):
        s = f*delta
        a = np.sinh(delta-s)/np.sinh(delta)
        b = np.sinh(s)/np.sinh(delta)
        sig2 = (-np.expm1(-2*s))*(-np.expm1(-2*(delta-s)))/den
        max_mean_error = max(max_mean_error, abs(a+q*b-np.exp(-s)))
        max_variance_error = max(max_variance_error, abs(a*a+b*b+2*q*a*b+sig2-1))
        for yf in [-1.0, -0.5, 0.0, 0.5, 1.0]:
            y = yf*min(s, delta-s)
            lnorm2 = abs(np.exp(-(s+1j*y)))**2 + abs(2*q*np.sinh(s+1j*y)/np.sqrt(den))**2
            identity = 1-sig2+4*q*q*np.sin(y)**2/den
            max_complex_identity_error = max(max_complex_identity_error, abs(lnorm2-identity))
            max_wedge_excess = max(max_wedge_excess, lnorm2-1)
    # Widest allowed doubling panels, plus reflection, on a rho=2 ellipse.
    for left in np.geomspace(delta*1e-6, delta/4, 40):
        right = min(2*left, delta/2)
        mid, half = (left+right)/2, (right-left)/2
        for theta in np.linspace(0, 2*np.pi, 400):
            z = mid + half*(1.25*np.cos(theta)+0.75j*np.sin(theta))
            for z in [z, delta-z]:
                max_ellipse_excess = max(max_ellipse_excess, abs(z.imag)-min(z.real, delta-z.real))

out["bridge"] = {
    "max_conditional_mean_identity_error": max_mean_error,
    "max_marginal_variance_identity_error": max_variance_error,
    "max_complex_norm_identity_error": max_complex_identity_error,
    "max_contraction_wedge_excess": max_wedge_excess,
    "max_rho2_ellipse_wedge_excess": max_ellipse_excess,
}

D, trials = 7, 300
max_prefix_derivative_ratio = 0.0
max_contraction_norm = 0.0
max_terminal_ratio = 0.0
max_noncommutation = 0.0
max_carrier_norm_error = 0.0
for _ in range(trials):
    A = 10**rng.uniform(-5, -1)
    w = A**(2/3)
    h = A**(3/2)
    r = np.sqrt(1-h*h)
    weights = w*rng.dirichlet(np.ones(8))
    coeff = rng.uniform(size=8)
    H = np.zeros((D,D))
    for beta, a in zip(weights, coeff):
        Q, _ = np.linalg.qr(rng.normal(size=(D,D)))
        Hj = A*(Q*rng.uniform(size=D))@Q.T
        H += beta*a*Hj
    Q, _ = np.linalg.qr(rng.normal(size=(D,D)))
    Hterm = A*(Q*rng.uniform(size=D))@Q.T
    C = np.eye(D)-H
    max_prefix_derivative_ratio = max(max_prefix_derivative_ratio, np.linalg.norm(H, 2)/(A*w))
    max_contraction_norm = max(max_contraction_norm, np.linalg.norm(C,2))
    max_terminal_ratio = max(max_terminal_ratio, np.linalg.norm(r*Hterm@C,2)/A)
    max_noncommutation = max(max_noncommutation, np.linalg.norm(Hterm@H-H@Hterm,2))
    t = rng.uniform()
    c = np.sqrt(1-t*t)
    max_carrier_norm_error = max(max_carrier_norm_error, abs((r*c)**2+h*h-(1-(r*t)**2)))
out["matrices"] = {
    "dimension": D, "trials": trials,
    "max_prefix_derivative_norm_over_Aw": max_prefix_derivative_ratio,
    "max_I_minus_prefix_derivative_norm": max_contraction_norm,
    "max_terminal_jacobian_norm_over_A": max_terminal_ratio,
    "max_nonzero_hessian_commutator_norm": max_noncommutation,
    "max_terminal_carrier_row_identity_error": max_carrier_norm_error,
}

# A literal positive dyadic Gauss rule, with the requested independent branch
# boundary inserted. Stable expm1 avoids falsely rounding very early nodes to 1.
leg_x, leg_w = np.polynomial.legendre.leggauss(12)
integrated = []
for A in [0.1, 0.03, 0.01, 0.001, 1e-5, 1e-8]:
    w, eta = A**(3/5), A
    eps = A**3
    early = eps/64
    late = np.log(64/eps)
    boundary = -np.log1p(-eta)
    edges = [early]
    while edges[-1] < late:
        edges.append(min(2*edges[-1], late))
    edges = sorted(set(edges+[boundary]))
    taus, weights = [early/2, late+1], [-np.expm1(-early), np.exp(-late)]
    for left, right in zip(edges[:-1], edges[1:]):
        nodes = (left+right)/2+(right-left)*leg_x/2
        taus.extend(nodes)
        weights.extend((right-left)*leg_w*np.exp(-nodes)/2)
    taus, weights = np.array(taus), np.array(weights)
    weights /= weights.sum()
    inv_v = (1-w)**2/(2*w-w*w)+1/(-np.expm1(-2*taus))
    bulk = taus >= boundary
    near = ~bulk
    ratios = {}
    for power in [1, 1.5, 7/6]:
        bound = w**(-power)+(np.log(1/eta) if power == 1 else eta**(1-power))
        ratios[str(power)] = float(np.sum(weights[bulk]*inv_v[bulk]**power)/bound)
    near_ratio = float(np.sum(weights[near]*np.sqrt(inv_v[near]))/(eta/np.sqrt(w)+np.sqrt(eta)))
    integrated.append({"A": A, "node_count": len(taus),
                       "bulk_sum_over_stated_bound": ratios,
                       "near_sum_over_stated_bound": near_ratio,
                       "min_bulk_v_over_eta": float(1/inv_v[bulk].max()/eta)})
out["narrow_endpoint_dyadic_checks"] = integrated

assert max_mean_error < 1e-10
assert max_variance_error < 1e-10
assert max_complex_identity_error < 1e-10
assert max_wedge_excess < 1e-10
assert max_ellipse_excess < 1e-12
assert max_prefix_derivative_ratio <= 1+1e-12
assert max_contraction_norm <= 1+1e-12
assert max_terminal_ratio <= 1+1e-12
assert max_carrier_norm_error < 1e-12
for case in integrated:
    assert max(case["bulk_sum_over_stated_bound"].values()) < 3
    assert case["near_sum_over_stated_bound"] < 3
    assert case["min_bulk_v_over_eta"] >= 0.5
out["passed"] = True
path = root / "bridge-diagnostics.json"
path.write_text(json.dumps(out, indent=2)+"\n")
print(json.dumps(out, indent=2))
