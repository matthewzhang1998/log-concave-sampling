#!/usr/bin/env python3
"""Independent algebra/numerical diagnostics; the markdown contains the proof."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.special import ndtr, ndtri
from scipy.optimize import brentq
from scipy.linalg import solve_sylvester

OUT = Path(__file__).parent
rng = np.random.default_rng(20261004)


def gaussian_rule(n):
    x, w = hermgauss(n)
    return np.sqrt(2) * x, w / np.sqrt(np.pi)


def hermite_checks():
    # If c_jj and c_jk are normalized chaos-2 coefficients, A_jj=sqrt(2)c_jj,
    # A_jk=A_kj=c_jk. This must give exactly 2 times coefficient energy.
    reports = []
    for d in [1, 2, 3, 8, 17]:
        c = rng.normal(size=d * (d + 1) // 2)
        A = np.zeros((d, d))
        q = 0
        for j in range(d):
            A[j, j] = np.sqrt(2) * c[q]
            q += 1
        for j in range(d):
            for k in range(j + 1, d):
                A[j, k] = A[k, j] = c[q]
                q += 1
        ratio = np.sum(A * A) / np.sum(c * c)
        assert abs(ratio - 2) < 1e-12
        reports.append({"dimension": d, "HS_energy_over_chaos_energy": float(ratio)})
    return reports


def tensor_checks():
    ratios = []
    for n, d in [(1, 2), (2, 1), (3, 5), (8, 4), (4, 8)]:
        T = rng.normal(size=(d, d, d))
        A = rng.normal(size=(d, n))
        lhs = np.linalg.norm(np.einsum("ijk,ka->ija", T, A))
        rhs = np.linalg.norm(T) * np.linalg.norm(A, 2)
        ratios.append(float(lhs / rhs))
        assert lhs <= rhs * (1 + 1e-12)
    return ratios


def trig_stein_checks():
    # n=1,d=2 map with truly nonsymmetric tau: f=(a sin(v), b cos(v)).
    v, w = gaussian_rule(160)
    r, rw = leggauss(120)
    r, rw = (r + 1) / 2, rw / 2
    a, b, sigma, t = 0.35, 0.21, 0.6, 0.73
    X = np.stack([a * np.sin(v), b * (np.cos(v) - np.exp(-0.5))], axis=1)
    Df = np.stack([a * np.cos(v), -b * np.sin(v)], axis=1)
    kernel = np.exp(-(1 - r*r)/2) * rw
    U = np.stack([a * (np.cos(v[:, None] * r) @ kernel),
                  -b * (np.sin(v[:, None] * r) @ kernel)], axis=1)
    tau = U[:, :, None] * Df[:, None, :]
    cov = (X.T * w) @ X
    B = tau - cov
    C = sigma**2 * np.eye(2) + (1-t*t)*cov
    e2 = np.einsum("n,ni,ni->", w, X, X)
    tau2 = np.einsum("n,nij,nij->", w, tau, tau)
    B2 = np.einsum("n,nij,nij->", w, B, B)
    L = max(a, b)
    assert tau2 <= L*L*e2 + 1e-13
    assert np.max(np.abs(np.einsum("n,nij->ij", w, tau)-cov)) < 1e-12

    K = np.array([[0.9, -0.7], [-0.3, 1.1]])
    phase = np.array([0.2, -0.4])
    theta = t * X @ K.T + phase
    damp = np.exp(-0.5 * np.einsum("ij,jk,ik->i", K, C, K))
    Eh = damp * np.sin(theta)
    EDh = (damp * np.cos(theta))[:, :, None] * K[None, :, :]
    lhs = np.einsum("n,ni,ni->", w, X, Eh)
    rhs = t*np.einsum("n,nij,nij->", w, tau, EDh)
    wrong = t*np.einsum("n,nji,nij->", w, tau, EDh)
    assert abs(lhs-rhs) < 1e-12
    assert abs(lhs-wrong) > 1e-6, "Chosen test must detect transposition"

    # Scalar test phi=cos(k.y)+0.3 sin(q.y), analytic integration over Z.
    k, q = K
    xk, xq = t*X@k, t*X@q
    dk, dq = np.exp(-0.5*k@C@k), np.exp(-0.5*q@C@q)
    Eg = -dk*np.sin(xk)[:, None]*k + 0.3*dq*np.cos(xq)[:, None]*q
    EH = (-dk*np.cos(xk)[:, None, None]*np.outer(k, k)
          -0.3*dq*np.sin(xq)[:, None, None]*np.outer(q, q))
    direct = np.einsum("n,ni,ni->", w, X, Eg) - t*np.einsum("n,ij,nij->", w, cov, EH)
    velocity = t*np.einsum("n,nij,nij->", w, B, EH)
    assert abs(direct-velocity) < 1e-12

    # Sample independent vector test fields for the centered-matrix inequality.
    dual_ratios = []
    for _ in range(40):
        K = rng.normal(size=(2, 2)) * rng.uniform(0.2, 2.5)
        ph = rng.normal(size=2)
        theta = t*X@K.T+ph
        kCk = np.einsum("ij,jk,ik->i", K, C, K)
        EDh = (np.exp(-kCk/2)*np.cos(theta))[:, :, None]*K[None, :, :]
        pairing = abs(np.einsum("n,nij,nij->", w, B, EDh))
        h2 = np.einsum("n,ni->", w, (1-np.exp(-2*kCk)*np.cos(2*theta))/2)
        bound = np.sqrt(2)*t*L/np.linalg.eigvalsh(C).min()*np.sqrt(B2*h2)
        dual_ratios.append(float(pairing/bound))
    assert max(dual_ratios) <= 1+1e-12
    return {
        "mean_tau_covariance_max_error": float(np.max(np.abs(np.einsum("n,nij->ij", w, tau)-cov))),
        "tau_nonsymmetry_L2": float(np.sqrt(np.einsum("n,nij,nij->", w, tau-tau.transpose(0,2,1), tau-tau.transpose(0,2,1)))),
        "tau_squared_norm_over_L2_e2": float(tau2/(L*L*e2)),
        "vector_Stein_identity_error": float(abs(lhs-rhs)),
        "wrong_transpose_error": float(abs(lhs-wrong)),
        "continuity_velocity_identity_error": float(abs(direct-velocity)),
        "maximum_centered_matrix_dual_ratio_40_tests": max(dual_ratios),
    }


def one_dimensional_w2_checks():
    # Bounded Lipschitz skew map. Quantile integration is a diagnostic, not a proof.
    v, w = gaussian_rule(160)
    p, pw = leggauss(512)
    p, pw = (p+1)/2, pw/2
    sigma, beta = 0.7, 0.4
    base = np.sin(v)+beta*(np.cos(v)-np.exp(-0.5))
    reports = []
    for a in [0.05, 0.1, 0.2, 0.4]:
        x = a*base
        var = np.dot(w, x*x)
        sd = np.sqrt(sigma*sigma+var)
        # Global support bound for X, used solely as a safe quantile bracket.
        boundX = a*(1+beta*(1+np.exp(-0.5)))
        def quantile(p0):
            z = sigma*ndtri(p0)
            return brentq(lambda y: np.dot(w, ndtr((y-x)/sigma))-p0,
                          z-boundX-1e-10, z+boundX+1e-10,
                          xtol=1e-13, rtol=1e-14)
        exact_q = np.array([quantile(p0) for p0 in p])
        approx_w2 = np.sqrt(np.dot(pw, (exact_q-sd*ndtri(p))**2))
        L = a*np.sqrt(1+beta*beta)
        theorem_bound = np.sqrt(2)/(3*sigma*sigma)*L*L*np.sqrt(var)
        assert approx_w2 <= theorem_bound*(1+1e-6)
        reports.append({"amplitude": a, "approximate_W2": float(approx_w2),
                        "theorem_bound": float(theorem_bound),
                        "ratio": float(approx_w2/theorem_bound),
                        "approximate_W2_over_amplitude_cubed": float(approx_w2/a**3)})
    return reports


def psqrt(A):
    e, U = np.linalg.eigh(A)
    assert e.min() > 0
    return (U*np.sqrt(e))@U.T


def noncommuting_buffer_checks():
    reports = []
    for d in [2, 3, 7]:
        A, B = rng.normal(size=(d,d)), rng.normal(size=(d,d))
        Q, cov = 0.4*np.eye(d)+A@A.T, 0.1*(B@B.T)
        t = 0.71
        S = psqrt(Q+(1-t*t)*cov)
        Sprime = solve_sylvester(S, S, -2*t*cov)
        H = rng.normal(size=(d,d)); H = (H+H.T)/2
        err = abs(np.sum((Sprime@S)*H)+t*np.sum(cov*H))
        wrong = np.linalg.norm(Sprime+t*cov@np.linalg.inv(S))
        assert err < 1e-10
        assert wrong > 1e-5
        reports.append({"dimension":d, "Gaussian_covariance_derivative_error":float(err),
                        "incorrect_commuting_derivative_discrepancy":float(wrong)})
    return reports


def positive_join_checks():
    max_root_ratio, max_cov_ratio, count = 0., 0., 0
    for d in [1, 2, 5, 13]:
        for alpha in [0.02, 0.1, 0.5]:
            for h in [0.01, 0.25, 0.5]:
                A = rng.normal(size=(d,2*d))
                A *= alpha/np.linalg.norm(A,2)
                E = rng.normal(size=(d,2*d))
                E *= alpha*alpha/np.linalg.norm(E,2)
                cov, covstar = A@A.T, (A+E)@(A+E).T
                gap = 1-2*h*h
                eta = zeta = np.sqrt(gap/2)
                c = h*h/(2*eta)
                target = gap*np.eye(d)+h*h*cov
                linear = (eta*np.eye(d)+c*cov)@(eta*np.eye(d)+c*cov)+zeta*zeta*np.eye(d)
                assert np.linalg.norm(linear-target-c*c*(cov@cov)) < 1e-12
                assert np.linalg.norm(h*h*np.eye(d)+target-((1-h*h)*np.eye(d)+h*h*cov)) < 1e-12
                energy = np.linalg.norm(A)
                assert np.linalg.norm(cov@cov) <= alpha**3*energy*(1+1e-12)
                priced = c*c*np.linalg.norm(cov@cov)/(2*np.sqrt(gap))
                # Both reserve matrices are functions of cov. Evaluate their tiny
                # root difference stably; subtracting two numerical roots is
                # roundoff-dominated at the smallest h and alpha.
                lam = np.linalg.eigvalsh(cov)
                base = gap+h*h*lam
                shift = c*c*lam*lam
                actual = np.linalg.norm(shift/(np.sqrt(base+shift)+np.sqrt(base)))
                assert actual <= priced*(1+1e-12)
                max_root_ratio = max(max_root_ratio, actual/priced)
                covbound = np.linalg.norm(E)*(np.linalg.norm(A,2)+np.linalg.norm(A+E,2))
                covratio = np.linalg.norm(covstar-cov)/covbound
                assert covratio <= 1+1e-12
                max_cov_ratio = max(max_cov_ratio,covratio)
                # Ratios of reserve and raw-Gaussianization terms to h^2 alpha^3.
                assert h*h*alpha <= 0.25
                assert h*h*np.sqrt(alpha) <= 0.25
                assert h/(1-h*h) <= 2/3+1e-15
                count += 1
    heat_count = 0
    for s in [1., 0.4, 0.01, 1e-5]:
        for A in [0.01,0.1,0.5]:
            for frac in [1e-8,0.01,0.25]:
                Delta = frac*s*s
                h, alpha = np.sqrt(Delta)/s, A*s*s
                lhs, rhs = h*h*alpha**3, Delta*A**3*s**4
                assert abs(lhs/rhs-1)<1e-12
                assert h<=0.5
                heat_count += 1
    return {"matrix_join_cases":count,
            "one_marked_energy_covariance_max_ratio":float(max_cov_ratio),
            "positive_quadratic_covariance_root_max_ratio":float(max_root_ratio),
            "reverse_bridge_heat_cases":heat_count,
            "scope":"Exact Gaussian algebra and source-qualified error scaling; imported compilers not instantiated"}


def main():
    result = {"status": "PASS", "scope": "Algebra and numerical diagnostics; analytical proof is in the audit markdown",
              "second_Hermite_factor": hermite_checks(),
              "chain_rule_tensor_contraction_ratios": tensor_checks(),
              "nonsymmetric_Stein_and_velocity": trig_stein_checks(),
              "scalar_quantile_diagnostics": one_dimensional_w2_checks(),
              "noncommuting_buffer": noncommuting_buffer_checks(),
              "positive_buffered_join": positive_join_checks()}
    (OUT/"gaussianization_checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
