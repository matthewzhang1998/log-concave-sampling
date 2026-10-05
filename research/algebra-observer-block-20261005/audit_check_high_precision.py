#!/usr/bin/env python3
"""Independent high-precision diagnostics for the sine-only finite-order solver.

Requires mpmath. Uses truncated polynomial convolution, not the author's
composition enumeration. This supports, and does not replace, the Taylor proof.
"""
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 90
ROOT = Path(__file__).resolve().parent


def multiply_truncated(left, right, length):
    output = [mp.mpf(0)] * length
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j < length:
                output[i + j] += x * y
    return output


def inverse_coefficients(x, beta, order):
    denominator = 1 + beta * mp.cos(x)
    coefficients = [mp.mpf(0), -mp.sin(x) / denominator]
    for n in range(2, order + 1):
        power = [mp.mpf(1)]
        total = mp.mpf(0)
        for rank in range(1, n + 1):
            power = multiply_truncated(power, coefficients, n + 1)
            if rank >= 2:
                derivative = mp.sin(x + rank * mp.pi / 2)
                total += derivative * power[n] / mp.factorial(rank)
        coefficients.append(-beta * total / denominator)
    return coefficients


def error_at(x, beta, order):
    coefficients = inverse_coefficients(x, beta, order)
    output = x + sum(beta**j * coefficients[j] for j in range(1, order + 1))
    return abs(output + beta * mp.sin(output) - x)


x = mp.mpf("0.43")
betas = [mp.mpf(1) / 1024, mp.mpf(1) / 2048]
rows = []
for order in range(1, 7):
    errors = [error_at(x, beta, order) for beta in betas]
    rate = mp.log(errors[0] / errors[1], 2)
    passed = abs(rate - (order + 2)) < mp.mpf("0.1")
    rows.append(
        {
            "N": order,
            "predicted_order": order + 2,
            "observed_halving_rate": mp.nstr(rate, 35),
            "residuals": [mp.nstr(error, 45) for error in errors],
            "passed": bool(passed),
        }
    )
    assert passed, (order, rate)

report = {
    "status": "PASS",
    "precision_decimal_digits": mp.mp.dps,
    "mpmath_version": mp.__version__,
    "method": "independent truncated polynomial convolution",
    "scope": "Finite numerical diagnostics of the N+2 scalar residual order; not a proof of uniform bounds or native admissibility.",
    "x": str(x),
    "betas": [str(beta) for beta in betas],
    "checks": rows,
}
output = ROOT / "audit-high-precision.json"
output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
