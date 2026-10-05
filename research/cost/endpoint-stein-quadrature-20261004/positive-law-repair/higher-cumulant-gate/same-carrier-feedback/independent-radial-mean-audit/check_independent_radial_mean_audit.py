#!/usr/bin/env python3
"""Independent algebra/scalar diagnostics for the exact canonical radial proof.

This is NOT a simulation of an OU history and NOT a numerical proof of the
uniform asymptotic remainder.  The accompanying audit supplies those arguments.
No true-history or finite shared-root genealogy is replaced in the audit.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import gammaln
import sympy as sp



# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
SOURCE = BASE / "RADIAL-SHELL-REFUTES-UNSHIFTED-SINGLE-HISTORY-COVARIANCE.md"
EXPECTED_SHA256 = "2807f7307646cc16c79049319b574259874b09e731a13f686696bdd60a05b035"
C = D_COEFFICIENT = 0.25
SHIFT = C / 2


def eta(s: float) -> float:
    if abs(s) >= 2:
        return 0.0
    return math.exp(-1 / (1 - s * s / 4))


Z, z_error = quad(eta, -2, 2, epsabs=1e-13, epsrel=1e-13)


def bump(s: float) -> float:
    return eta(s) / Z


def bump_derivative(s: float) -> float:
    if abs(s) >= 2:
        return 0.0
    q = 1 - s * s / 4
    return -s * bump(s) / (2 * q * q)


def psi(s: float) -> float:
    if s <= -2:
        return 0.0
    if s >= 2:
        return 1.0
    return quad(bump, -2, s, epsabs=1e-12, epsrel=1e-12)[0]


def finite_rule(dimension: int) -> dict:
    delta = 1 / dimension
    panel_count = math.ceil(math.log2(4 / delta))
    order = math.ceil(math.log(16 / delta, 4))
    gl_nodes, gl_weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = [], []
    for k in range(panel_count):
        left, right = 2.0 ** (-k - 1), 2.0 ** (-k)
        t = (left + right) / 2 + (right - left) * gl_nodes / 2
        nodes.extend(1 - t)
        weights.extend((right - left) * gl_weights / 2)
    final_width = 2.0 ** (-panel_count)
    nodes.append(1 - final_width / 2)
    weights.append(final_width)
    nodes, weights = np.array(nodes), np.array(weights)
    analytical_uniform_bound = 8 * 4.0 ** (-order) + 2.0 ** (1 - panel_count)
    beta = float(np.dot(weights, np.sqrt(1 - nodes * nodes)))
    second_moment = float(np.dot(weights, nodes * nodes))
    degree_sample = list(range(513)) + [2**k for k in range(10, 31)]
    sampled_errors = [abs(float(weights @ nodes**n) - 1 / (n + 1)) for n in degree_sample]
    assert np.all(weights > 0) and np.all(nodes > 0) and np.all(nodes < 1)
    assert abs(float(sum(weights)) - 1) < 2e-14
    assert abs(float(weights @ nodes) - 0.5) < 2e-14
    assert analytical_uniform_bound <= delta * (1 + 1e-12)
    assert max(sampled_errors) <= analytical_uniform_bound + 1e-14
    assert abs(second_moment - (1 / 3 - final_width**3 / 12)) < 2e-14
    assert beta >= 2 / 3 - delta - 1e-14
    return {
        "D": dimension,
        "delta": delta,
        "dyadic_panels": panel_count,
        "gauss_order_per_panel": order,
        "node_count": len(nodes),
        "minimum_weight": float(weights.min()),
        "mass": float(sum(weights)),
        "first_moment": float(weights @ nodes),
        "second_moment": second_moment,
        "analytical_uniform_Hermite_bound": analytical_uniform_bound,
        "sampled_Hermite_error": max(sampled_errors),
        "sampled_degrees_max": max(degree_sample),
        "beta_Q": beta,
        "beta_Q_squared_minus_quarter": beta**2 - 0.25,
        "note": "The analytical bound covers all degrees; the finite degree sample is only a diagnostic.",
    }


def symbolic_radial_contraction() -> dict:
    # Verify all nine input entries for a completely general, nonsymmetric B.
    # Rotational covariance then supplies the displayed invariant formula.
    x = sp.symbols("x0:3", real=True)
    radius = sp.symbols("r", positive=True)
    norm = sp.sqrt(sum(u * u for u in x))
    phi = sp.Function("phi")
    b = sp.symbols("b0:9")
    B = sp.Matrix(3, 3, b)
    f = sp.Matrix([phi(norm) * u / norm for u in x])
    at_axis = {x[0]: radius, x[1]: 0, x[2]: 0}
    actual = sp.Matrix([
        sum(sp.diff(f[k], x[i], x[j]).subs(at_axis) * B[i, j]
            for i in range(3) for j in range(3))
        for k in range(3)
    ])
    a = sp.diff(phi(radius), radius) / radius - phi(radius) / radius**2
    predicted = sp.Matrix([
        sp.diff(phi(radius), radius, 2) * B[0, 0] + a * (B[1, 1] + B[2, 2]),
        a * (B[1, 0] + B[0, 1]),
        a * (B[2, 0] + B[0, 2]),
    ])
    residuals = [sp.simplify(v) for v in actual - predicted]
    assert residuals == [0, 0, 0]
    return {"dimension_for_symbolic_identity": 3, "general_nonsymmetric_matrix": True,
            "component_residuals": [str(x) for x in residuals]}


def hs_operator_norm_diagnostics() -> dict:
    largest_error = 0.0
    cases = 0
    for dimension in (2, 3, 7, 16):
        # This algebra check permits arbitrary phi values.  Smoothness/support
        # and the uniform-in-D shell estimate are proved separately in the audit.
        for phipp, a in ((0.7, -0.2), (0.0, 0.4), (0.1, 0.0)):
            matrix = np.zeros((dimension, dimension**2))
            matrix[0, 0] = phipp
            for i in range(1, dimension):
                matrix[0, i * dimension + i] = a
                matrix[i, i * dimension] = a
                matrix[i, i] = a
            measured = float(np.linalg.svd(matrix, compute_uv=False)[0])
            formula = max(math.sqrt(phipp**2 + (dimension - 1) * a**2), math.sqrt(2) * abs(a))
            largest_error = max(largest_error, abs(measured - formula))
            assert abs(measured - formula) < 1e-13
            cases += 1
    return {"cases": cases, "maximum_absolute_error": largest_error,
            "exact_norm": "max(sqrt(psi_double_prime^2+(D-1)*a_r^2), sqrt(2)*abs(a_r))"}


def shell_density(s: float, dimension: int) -> float:
    radius = math.sqrt(dimension) + s
    if radius <= 0:
        return 0.0
    log_density = ((1 - dimension / 2) * math.log(2) - gammaln(dimension / 2)
                   + (dimension - 1) * math.log(radius) - radius * radius / 2)
    return math.exp(log_density)


def finite_witness(dimension: int) -> dict:
    radius = math.sqrt(dimension)
    def integrand(s):
        return 0.5 * (bump(s - SHIFT) - bump(s)) * (1 + s / radius) * shell_density(s, dimension)
    points = sorted(set([-2.0, SHIFT - 2, 2.0, SHIFT + 2]))
    value = error = 0.0
    for lo, hi in zip(points[:-1], points[1:]):
        part, part_error = quad(integrand, lo, hi, epsabs=1e-12, epsrel=1e-10)
        value += part
        error += part_error
    return {"D": dimension, "linear_resolvent_witness": value,
            "quadrature_error_estimate": error,
            "note": "Evaluates the leading-field witness only, not the unknown full finite-D remainder."}


def main() -> None:
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA256, "Canonical exact text has changed."
    t = sp.symbols("t", positive=True)
    linear_mean = sp.integrate(sp.exp(-2 * t), (t, 0, sp.oo))
    linear_variance = 2 * sp.integrate(t * sp.exp(-2 * t), (t, 0, sp.oo)) - sp.Rational(1, 4)
    assert linear_mean == sp.Rational(1, 2)
    assert linear_variance == sp.Rational(1, 4)
    beta_floor = sp.Rational(7, 12)
    covariance_floor = beta_floor**2 - sp.Rational(1, 4)
    k_floor = sp.Rational(1, 4) * sp.Rational(1, 4)**2 / 2 * covariance_floor
    assert covariance_floor == sp.Rational(13, 144)
    assert k_floor == sp.Rational(13, 18432)
    derivative_upper_bound = math.exp(1 / 3) / 2
    assert derivative_upper_bound < 1
    gaussian_density = lambda s: math.exp(-s*s) / math.sqrt(math.pi)
    h_expectation, h_error = quad(
        lambda u: bump(u) * (gaussian_density(u + SHIFT) - gaussian_density(u)),
        -2, 2, epsabs=1e-13, epsrel=1e-12,
    )
    assert h_expectation < 0
    dim = sp.symbols("D", positive=True)
    heat = dim**(-sp.Rational(1, 2))
    allowance = sp.simplify(heat**4 * sp.sqrt(dim))
    ratio = sp.simplify(heat**2 / allowance)
    assert allowance == dim**(-sp.Rational(3, 2)) and ratio == sp.sqrt(dim)
    results = {
        "status": "PASS",
        "purpose": "Independent exact-text hash, symbolic identities, finite positive-clock diagnostics, and scalar witness checks.",
        "limitations": [
            "The checker does not prove uniform Lp or L2 remainders; the accompanying analytical audit does.",
            "Floating-point quadrature error estimates are diagnostics, not formal interval certificates.",
            "No OU path discretization, simulated cheap-law substitution, or full-m3 validation is performed.",
            "No explicit dimension threshold for the full residual lower bound is claimed.",
        ],
        "canonical": {"path": source_label(SOURCE), "sha256": digest},
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "sympy": sp.__version__},
        "bump": {"normalizer": Z, "normalizer_error_estimate": z_error,
                 "maximum_psi_prime": bump(0), "analytical_upper_bound": derivative_upper_bound,
                 "radial_Hessian_upper_ratio": 0.5, "tangential_Hessian_upper_ratio": 0.375},
        "OU_linear_integral": {"mean_coefficient": str(linear_mean), "conditional_variance": str(linear_variance)},
        "uniform_finite_clock_gap": {"beta_lower_bound": str(beta_floor),
                                     "beta_squared_minus_quarter_lower_bound": str(covariance_floor),
                                     "absolute_k_lower_bound": str(k_floor)},
        "finite_positive_clocks": [finite_rule(d) for d in (16, 256, 4096, 65536)],
        "radial_contraction_symbolic": symbolic_radial_contraction(),
        "HS_to_vector_operator_norm": hs_operator_norm_diagnostics(),
        "limiting_shell_witness": {
            "S_distribution": "N(0, 1/2)", "radial_shift": SHIFT,
            "E_h_S": h_expectation, "quadrature_error_estimate": h_error,
            "resolvent_linear_witness_limit": h_expectation / 2,
            "uniform_absolute_leading_witness_coefficient_floor": float(k_floor) * abs(h_expectation) / 2,
            "possible_eventual_c_star": float(k_floor) * abs(h_expectation) / 8,
            "strict_sign_proof": "Positive layer-cake intervals and strictly decreasing translated centered Gaussian interval mass; see audit.",
        },
        "finite_dimension_leading_witness": [finite_witness(d) for d in (16, 64, 256, 1024, 16384, 262144)],
        "asymptotic_separation": {"A": "D^(-1/2)", "claimed_scale": str(allowance),
                                  "lower_scale": "D^(-1)", "scale_ratio": str(ratio),
                                  "fixed_polylog_limit": "sqrt(D)/(log D)^p -> infinity for every fixed p"},
    }
    destination = HERE / "independent_radial_mean_checks.json"
    destination.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": results["status"], "canonical_sha256": digest,
                      "E_h_S": h_expectation, "results": str(destination)}, indent=2))


if __name__ == "__main__":
    main()
