#!/usr/bin/env python3
"""Independent numerical check of the actual coherent radial VALUE graph.

This is a high-dimensional deterministic-row limit, not an executable source
replacement. Exact retained Gaussian row correlations are used at every node.
The proof and the finite-rule argument are in INDEPENDENT-RADIAL-AUDIT.md.
"""
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


ALPHA = BETA = 0.5
C1 = ALPHA / 2
CV = ALPHA * (1 - 1 / math.sqrt(2))
TAU = CV / 2
U0 = CV / (2 * C1)
EXACT_INNER = (11 - 17 / math.sqrt(2)) / 96


def hinge(t, eta):
    if eta == 0:
        return np.maximum(t, 0)
    return np.where(t <= -eta, 0,
                    np.where(t >= eta, t, (t + eta) ** 2 / (4 * eta)))


def force(row, A, eta=0):
    norm = np.linalg.norm(row, axis=-1)
    phi = ALPHA * hinge(norm - 0.5, eta)
    phi += BETA * hinge(norm - (1 - A * TAU), eta)
    multiplier = np.divide(A * phi, norm, out=np.zeros_like(norm), where=norm > 0)
    return row * multiplier[..., None]


def gauss(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


def one_rows(r):
    a = -np.log(r)
    coeff = np.stack((r, 2 * a * r, 2 * a * (a - 1) * r), axis=-1)
    residual = 1 - np.sum(coeff * coeff, axis=-1)
    assert residual.min() > -1e-12
    return np.concatenate((coeff, np.sqrt(np.maximum(residual, 0))[:, None],
                           np.zeros((len(r), 1))), axis=-1)


def pair_rows(r, z):
    # a~Exp(rate), h~Exp(1), with r=exp(-a), z=exp(-h).
    r, z = np.broadcast_arrays(r, z)
    a, h = -np.log(r), -np.log(z)
    sigma_a, sigma_h = np.sqrt(1-r*r), np.sqrt(1-z*z)
    fa, fh = 2*a*r/sigma_a, 2*h*z/sigma_h
    R11, R12 = fa, fa*(a-1)
    R21, R22 = r*fh, r*fh*(2*a+h-1)
    K11 = 1-R11*R11-R12*R12
    K12 = -R11*R21-R12*R22
    K22 = 1-R21*R21-R22*R22
    determinant = K11*K22-K12*K12
    assert determinant.min() > 0
    root_det = np.sqrt(determinant)
    den = np.sqrt(K11+K22+2*root_det)
    B11, B12, B22 = (K11+root_det)/den, K12/den, (K22+root_det)/den
    W1 = np.stack((np.zeros_like(r), R11, R12, B11, B12), axis=-1)
    W2 = np.stack((np.zeros_like(r), R21, R22, B12, B22), axis=-1)
    U = sigma_a[..., None]*W1
    U[..., 0] = r
    V = z[..., None]*U+sigma_h[..., None]*W2
    assert np.max(np.abs(np.sum(U*U, axis=-1)-1)) < 2e-12
    assert np.max(np.abs(np.sum(V*V, axis=-1)-1)) < 2e-12
    assert np.max(np.abs(np.sum(U*V, axis=-1)-z)) < 2e-12
    return U, V


def actual_limit(A, n, smoothed):
    eta = A*A if smoothed else 0
    x = np.array([1., 0, 0, 0, 0])
    v = np.array([.5, .5, 0, 0, 0])
    w = np.array([.25, .5, .25, 0, 0])
    g = lambda y: force(y, A, eta)
    gv, b0 = g(v), g(g(w))
    S = x-gv+b0
    r, weights = gauss(n)
    Us = one_rows(r)
    ds = gv-g(Us)
    single = np.dot(weights, ((g(S+ds)-g(S-ds))/2)[..., 0])
    U, V = pair_rows(r[:, None], r[None, :])
    GU, GV = g(U), g(V)
    e = GU-g(U-GV)-b0
    product_weights = weights[:, None]*weights[None, :]
    linear = np.sum(product_weights*((g(S+e)-g(S-e))/2)[..., 0])
    d1, d2 = gv-GU, gv-GV
    q = (g(S+d1+d2)+g(S-d1-d2)-g(S+d1-d2)-g(S-d1+d2))/8
    quadratic = np.sum(product_weights*(2*r[:, None])*q[..., 0])
    graph = g(S)[0]+single+linear+quadratic
    k = g(x)[0]
    rho_inner = math.sqrt(1-k+k*k/2)
    ell = force(np.array([rho_inner, 0, 0, 0, 0]), A, eta)[0]/rho_inner
    target_row = x-ell*v+ell*k*w
    target = g(target_row)[0]
    return {
        "A": A, "quadrature_order": n, "smoothed": smoothed,
        "inner_discrepancy_over_A2": (graph-target)/(A*A),
        "outer_discrepancy_over_A2": (graph-target)/(2*A*A),
        "terminal_target_radius": float(np.linalg.norm(target_row)),
        "upper_gate_radius": 1-A*TAU,
        "b0_norm": float(np.linalg.norm(b0)),
    }


def main():
    r, weights = gauss(700)
    p, q = U0-r[:, None], U0-r[None, :]
    q_integral = np.sum(weights[:, None]*weights[None, :]
                        *(np.abs(p+q)-np.abs(p-q))/8)
    exact_q = 1/12-U0/4+U0**3/3
    rows = [actual_limit(A, 320, smooth)
            for smooth in (False, True)
            for A in (1/100, 1/300, 1/1000, 1/3000, 1/10000)]
    result = {
        "status": "independent numerical audit; not a finite-precision proof",
        "exact_inner_coefficient": EXACT_INNER,
        "exact_outer_coefficient": EXACT_INNER/2,
        "single_unit_coefficient": U0/2-1/4,
        "quadratic_unit_coefficient": exact_q,
        "quadratic_numerical_check": float(q_integral),
        "quadratic_check_error": float(q_integral-exact_q),
        "actual_graph_checks": rows,
    }
    assert abs(q_integral-exact_q) < 1e-6
    assert all(row["terminal_target_radius"] < row["upper_gate_radius"] for row in rows)
    assert all(row["b0_norm"] == 0 for row in rows)
    assert abs(rows[-1]["inner_discrepancy_over_A2"]-EXACT_INNER) < 3e-5
    destination = Path(__file__).with_name("independent_radial_audit_checks.json")
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
