#!/usr/bin/env python3
"""Independent finite diagnostics for the bounded linear-backbone m3 source.

No author checker is imported. These tests do not execute the imported complete
gradient/near-gradient mean programs and do not replace the continuum proof.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md"
PIN = "798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3"
ROOT = source_path("research")
LOW30 = Path("external:LOW30")
IMPORTS = [
    LOW30,
    ROOT / "cost/endpoint-stein-quadrature-20261004/ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md",
    ROOT / "no-copy-rank-20261004/host-observer-audit/COMPLETE-MEAN-FIRST-ORDERING-AND-REENTRY-RECIPE.md",
    ROOT / "no-copy-rank-20261004/host-observer-audit/INDEPENDENT-COMPLETE-MEAN-REENTRY-AUDIT.md",
    HERE.parent.parent / "P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md",
    HERE.parent.parent / "RADIAL-INFLATION-REFUTES-CENTERED-SINGLE-HISTORY-COVARIANCE.md",
    HERE.parent.parent.parent / "THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md",
    HERE.parent.parent.parent / "nested-mean-independent-audit/INDEPENDENT-NESTED-FORCE-MEAN-AUDIT.md",
]
COUNTS: dict[str, int] = {}
DETAILS: dict[str, object] = {}


def check(condition, group, message):
    COUNTS[group] = COUNTS.get(group, 0) + 1
    if not bool(condition):
        raise AssertionError(f"{group}: {message}")


def close(actual, expected, tol, group, message):
    check(np.max(np.abs(np.asarray(actual) - np.asarray(expected))) <= tol, group, message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(lam):
    a = lam / 2 - lam * lam / 4
    b2 = lam * lam / 4 - lam**3 / 2 + 5 * lam**4 / 16
    return a, math.sqrt(max(0.0, b2))


def dyadic_rule(delta):
    # The two terms in 8*4^-m + 2^(1-K) each receive delta/2.
    m = math.ceil(math.log(16 / delta, 4))
    K = math.ceil(math.log2(4 / delta))
    u, v = leggauss(m)
    ts, ws = [], []
    for k in range(K):
        lo, hi = 2.0 ** (-k - 1), 2.0 ** (-k)
        t = (lo + hi) / 2 + (hi - lo) * u / 2
        ts.extend(1 - t)
        ws.extend((hi - lo) * v / 2)
    ts.append(1 - 2.0 ** (-K - 1))
    ws.append(2.0 ** (-K))
    return np.array(ts), np.array(ws), K, m


class Source:
    """Anchored convex gradient with bounded sinusoidal remainder.

    Hessians at different points generally do not commute. k changes Hessian
    oscillation while leaving the global Hessian interval unchanged.
    """
    def __init__(self, A, k):
        self.A, self.k, self.lam = A, k, A / 2
        self.v = np.array([[1., 0., 0.], [1., 2., 0.], [.5, -.2, 1.]])
        self.v /= np.linalg.norm(self.v, axis=1)[:, None]
        self.eta = A / (8 * len(self.v))
        self.eps = self.eta * np.linalg.norm(self.v, axis=1).sum() / k
        self.sites = []

    def value(self, x):
        self.sites.append(np.array(x).copy())
        return self.lam * x + self.remainder(x)

    def remainder(self, x):
        return self.eta / self.k * (np.sin(self.k * (self.v @ x)) @ self.v)

    def hessian(self, x):
        return self.lam * np.eye(3) + self.eta * np.einsum(
            "l,li,lj->ij", np.cos(self.k * (self.v @ x)), self.v, self.v)


def raw(src, Z, G, N, t, w, kind="F"):
    a, b = coefficients(src.lam)
    out = np.zeros_like(Z)
    for ti, wi in zip(t, w):
        x = ti * Z + math.sqrt(1 - ti * ti) * G
        if kind == "B":
            out += wi * src.value(x)
        elif kind == "F":
            out += wi * src.value((1 - a) * x - b * N)
        else:
            out += wi * (src.value((1 - a) * x - b * N) - src.value(x))
    return out


def jacobians(src, Z, G, N, t, w):
    a, b = coefficients(src.lam)
    JG, JN, JZ = (np.zeros((3, 3)) for _ in range(3))
    for ti, wi in zip(t, w):
        ci = math.sqrt(1 - ti * ti)
        x = ti * Z + ci * G
        H1, H0 = src.hessian((1 - a) * x - b * N), src.hessian(x)
        D = (1 - a) * H1 - H0
        JG += wi * ci * D
        JN -= wi * b * H1
        JZ += wi * ti * D
    return JG, JN, JZ


def symbolic_checks():
    a, b, lam = sp.symbols("a b lam", positive=True)
    K = (1 / (a + b)) * (1 / (a + 1) + 1 / (b + 1))
    ev = lambda x: sp.simplify(x.subs({a: 1, b: 1}))
    V1, C12, V2 = ev(K), ev(-sp.diff(K, b)), ev(sp.diff(K, a, b))
    for got, want, label in [
        (V1, sp.Rational(1, 2), "unconditional H1 variance"),
        (C12, sp.Rational(3, 8), "unconditional cross covariance"),
        (V2, sp.Rational(3, 8), "unconditional H2 variance")]:
        check(got == want, "symbolic", label)
    means = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 4)])
    conditional = sp.Matrix([[V1, C12], [C12, V2]]) - means * means.T
    check(conditional == sp.Matrix([[sp.Rational(1, 4), sp.Rational(1, 4)],
                                    [sp.Rational(1, 4), sp.Rational(5, 16)]]),
          "symbolic", "conditional covariance")
    check(conditional.det() == sp.Rational(1, 64), "symbolic", "positive determinant")
    row = sp.Matrix([[lam, -lam**2]])
    variance = sp.expand((row * conditional * row.T)[0])
    check(variance == lam**2 / 4 - lam**3 / 2 + 5 * lam**4 / 16,
          "symbolic", "b2 variance polynomial")
    check((row * means)[0] == lam / 2 - lam**2 / 4, "symbolic", "a2 mean")
    check(sp.expand(5 * (lam - sp.Rational(4, 5))**2 + sp.Rational(4, 5))
          == 4 - 8 * lam + 5 * lam**2, "symbolic", "positive polynomial completion")
    for ell in np.linspace(0, .5, 51):
        aa, bb = coefficients(ell)
        check(0 <= aa <= ell / 2 + 1e-16, "coefficient_guards", "a bound")
        check(0 <= bb <= ell / 2 + 1e-16, "coefficient_guards", "b bound")
        check(0 <= 1 - aa <= 1, "coefficient_guards", "alpha interval")
        check((1 - aa)**2 + bb**2 <= 1 + 1e-15,
              "coefficient_guards", "unconditional input variance at most one")
    DETAILS["exact_conditional_covariance"] = [["1/4", "1/4"], ["1/4", "5/16"]]
    DETAILS["exact_covariance_determinant"] = "1/64"


def quadrature_checks():
    data = []
    for delta in [1e-2, 1e-4, 1e-7]:
        t, w, K, m = dyadic_rule(delta)
        check(np.all(w > 0), "quadrature", "positive weights")
        check(np.all((t > 0) & (t < 1)), "quadrature", "interior nodes")
        close(w.sum(), 1, 2e-15, "quadrature", "mass one")
        close(w @ t, .5, 2e-15, "quadrature", "exact first moment numerically")
        check(8 * 4.0**(-m) + 2.0**(1-K) <= delta,
              "quadrature", "analytic sufficient error budget")
        degrees = np.unique(np.r_[np.arange(300), np.logspace(2, 9, 250).astype(int)])
        errors = [abs(w @ (t**int(n)) - 1 / (int(n) + 1)) for n in degrees]
        check(max(errors) <= delta, "quadrature", "sampled Hermite multiplier errors")
        beta = w @ np.sqrt(1 - t * t)
        check(.5 <= beta <= math.sqrt(3) / 2 + 1e-15,
              "quadrature", "beta range from first moment")
        data.append({"delta": delta, "nodes": len(t), "K": K, "m": m,
                     "max_tested_multiplier_error": max(errors), "beta": beta})
    DETAILS["quadrature"] = data


def finite_graph_checks():
    rng = np.random.default_rng(20261005)
    t, w, _, _ = dyadic_rule(1e-3)
    # Use a smaller exact positive rule for repeated differential diagnostics.
    u, ww = leggauss(7)
    t, w = (u + 1) / 2, ww / 2
    beta = w @ np.sqrt(1 - t * t)
    max_fd = max_skew = max_commutator = 0.0
    max_curl_ratio = max_private_ratio = 0.0
    for A in [.005, .03, .1]:
        for k in [1, 8, 64]:
            src = Source(A, k)
            a, b = coefficients(src.lam)
            for sample in range(6):
                Z, G, N = rng.normal(size=(3, 3)) * (1 + sample)
                H0, H1 = src.hessian(Z), src.hessian(G)
                max_commutator = max(max_commutator, np.linalg.norm(H0 @ H1 - H1 @ H0))
                check(np.linalg.eigvalsh(H0).min() >= 0, "raw_source", "convex Hessian")
                check(np.linalg.eigvalsh(H0).max() <= A, "raw_source", "upper Hessian bound")
                check(np.linalg.norm(src.remainder(Z)) <= src.eps + 1e-15,
                      "raw_source", "bounded linear remainder")
                src.sites.clear()
                F = raw(src, Z, G, N, t, w, "F")
                check(len(src.sites) == len(t), "counts", "one F leaf per node")
                Fsites = np.array(src.sites)
                src.sites.clear()
                B = raw(src, Z, G, N, t, w, "B")
                check(len(src.sites) == len(t), "counts", "one B leaf per node")
                src.sites.clear()
                E = raw(src, Z, G, N, t, w, "E")
                check(len(src.sites) == 2 * len(t), "counts", "two E leaves per node")
                close(E, F - B, 1e-14, "raw_source", "same-record split identity")
                xs = t[:, None] * Z + np.sqrt(1-t*t)[:, None] * G
                amp = A * (a * (w @ np.linalg.norm(xs, axis=1)) + b * np.linalg.norm(N))
                check(np.linalg.norm(E) <= amp + 1e-14, "energy", "literal chord bound")
                JG, JN, JZ = jacobians(src, Z, G, N, t, w)
                priv = np.hstack((JG, JN))
                lift = np.vstack((priv, np.zeros((3, 6))))
                curl = lift - lift.T
                max_skew = max(max_skew, np.linalg.norm(JG - JG.T))
                close(JG, JG.T, 1e-15, "first_curl", "exactly symmetric G block")
                close(np.linalg.norm(curl, 2), np.linalg.norm(JN, 2), 2e-15,
                      "first_curl", "curl norm equals off-diagonal norm")
                check(np.linalg.norm(JG, 2) <= A * beta + 1e-15, "first_curl", "G first bound")
                check(np.linalg.norm(JN, 2) <= A * b + 1e-15, "first_curl", "N/curl first bound")
                check(np.linalg.norm(priv, 2) <= A * math.sqrt(beta*beta+b*b) + 1e-15,
                      "first_curl", "full private first bound")
                check(np.linalg.norm(JZ, 2) <= A / 2 + 1e-15,
                      "first_curl", "captured-Z first bound")
                max_curl_ratio = max(max_curl_ratio, np.linalg.norm(curl, 2) / (A*b))
                max_private_ratio = max(max_private_ratio,
                    np.linalg.norm(priv, 2) / (A * math.sqrt(beta*beta+b*b)))
                direction = rng.normal(size=9)
                direction /= np.linalg.norm(direction)
                dz, dg, dn = direction.reshape(3, 3)
                h = 2e-5 / k
                fp = raw(src, Z+h*dz, G+h*dg, N+h*dn, t, w, "E")
                fm = raw(src, Z-h*dz, G-h*dg, N-h*dn, t, w, "E")
                fd = (fp-fm)/(2*h)
                err = np.linalg.norm(fd-(JZ@dz+JG@dg+JN@dn))
                max_fd = max(max_fd, err)
                check(err < 1e-8, "chain_rule", "literal complete first directional check")
                zeros = np.zeros(3)
                src.sites.clear()
                E0 = raw(src, Z, zeros, zeros, t, w, "E")
                check(len(src.sites) == 2 * len(t), "counts", "caller origin fully billed")
                check(np.linalg.norm(E0) <= A * a * np.linalg.norm(Z) / 2 + 1e-14,
                      "origins", "deterministic origin bound")
                close(raw(src, zeros, zeros, zeros, t, w, "E"), zeros, 0,
                      "origins", "literal total zero")
                src.sites.clear()
                raw(src, Z, G + .1, N, t, w, "F")
                check(len(src.sites) == len(t), "counts", "changed private argument replays all leaves")
                check(np.all(np.linalg.norm(np.array(src.sites)-Fsites, axis=1) > 0),
                      "counts", "every shifted original site actually changes")
    check(max_commutator > 1e-7, "raw_source", "fixture genuinely uses noncommuting Hessians")
    DETAILS["finite_graph"] = {"max_directional_first_error": max_fd,
        "max_G_block_skew": max_skew, "max_Hessian_commutator": max_commutator,
        "max_curl_to_bound_ratio": max_curl_ratio,
        "max_private_first_to_bound_ratio": max_private_ratio}


def history_checks():
    rng = np.random.default_rng(31799)
    ratios, algebra_errors = [], []
    n = 19
    w = np.exp(-np.arange(n)/4)
    w /= w.sum()
    conv = np.convolve(w, w)
    close(conv.sum(), 1, 1e-15, "history", "convolution weights keep mass one")
    for A in [.01, .08, .25]:
        for k in [1, 11]:
            src = Source(A, k)
            for scale in [1, 10, 1000]:
                X = rng.normal(size=(2*n-1, 3)) * scale
                values = np.array([src.value(x) for x in X])
                F1 = np.array([w @ values[j:j+n] for j in range(n)])
                H1j = np.array([w @ X[j:j+n] for j in range(n)])
                F2 = w @ np.array([src.value(X[j]-F1[j]) for j in range(n)])
                backbone = src.lam * (w @ X[:n]) - src.lam**2 * (conv @ X)
                explicit_R = -src.lam * (w @ (F1-src.lam*H1j)) + w @ np.array(
                    [src.remainder(X[j]-F1[j]) for j in range(n)])
                error = np.linalg.norm(F2-backbone-explicit_R)
                algebra_errors.append(error)
                check(error < 2e-13 * max(1, scale), "history", "exact genuine-future decomposition")
                bound = src.eps * (1+src.lam)
                check(np.linalg.norm(explicit_R) <= bound + 1e-14,
                      "history", "bounded remainder independent of path dimension/size")
                x = rng.normal(size=3) * scale
                terminal_error = np.linalg.norm(src.value(x-F2)-src.value(x-backbone))
                check(terminal_error <= A * bound + 1e-12,
                      "history", "terminal Lipschitz target-bias bound")
                ratios.append(terminal_error / (A * bound))
    DETAILS["finite_ancestry"] = {"max_decomposition_error": max(algebra_errors),
        "max_terminal_bias_to_bound_ratio": max(ratios),
        "qualification": "Positive discrete convolution diagnostic; continuum proved in report."}


def linear_law_and_scaling_checks():
    rng = np.random.default_rng(917)
    t, w, _, _ = dyadic_rule(1e-5)
    c = np.sqrt(1-t*t)
    beta = w @ c
    for lam in [0., 1e-4, .01, .1, .5]:
        a, b = coefficients(lam)
        # Shared G and shared N generate this full node covariance Gram.
        rows = np.column_stack(((1-a)*c, -b*np.ones_like(c)))
        node_cov = rows @ rows.T
        got_var = lam**2 * w @ node_cov @ w
        want_var = lam**2 * (((1-a)*beta)**2+b*b)
        close(got_var, want_var, 2e-15, "linear_law", "complete shared-root raw covariance")
        close(lam*(1-a)*(w@t), lam*(1-a)/2, 1e-15,
              "linear_law", "exact canonical m3 coefficient when epsilon is zero")
        check(np.linalg.eigvalsh(node_cov).min() >= -1e-12,
              "positivity", "actual Gaussian pushforward covariance is PSD")
        if lam > 0:
            wrong_var = lam**2 * np.sum(w*w*((1-a)**2*c*c+b*b))
            check(want_var > 2*wrong_var, "linear_law", "independent-node substitution changes raw law")
    for A in [.005, .02, .1, .5]:
        for factor in [1., 7., 100.]:
            D = factor*A**(-4)
            eps = .3*A
            lam = .8*A
            bias, grade = A*eps*(1+lam), A**4*math.sqrt(D)
            check(bias <= .3*(1+A)*grade*(1+1e-15),
                  "grade", "D lower bound converts absolute remainder to order four")
            check((A**3)*A*math.sqrt(D) <= grade*(1+1e-15),
                  "grade", "delta=A^3 gives order-four quadrature error")
            near_factor = A*(A+A) + A**3*(1+A**(-.5))
            check(near_factor <= 4*A*A, "grade", "mu=A near-gradient error multiplier")
    A, s = .1, .1
    alpha, eps_f = A*s*s, 2*s*(.3*A)
    check(eps_f/alpha > .3, "conditional_rescaling", "epsilon certificate does not retain same constant")
    check(A**(-4) < alpha**(-4), "conditional_rescaling", "dimension hypothesis not inherited")
    # Unit target covariance is the sum of independent positive variance shares.
    for v in [.1, .5, .9]:
        row = np.hstack((math.sqrt(v)*np.eye(3), math.sqrt(1-v)*np.eye(3)))
        close(row @ row.T, np.eye(3), 2e-15, "positivity", "literal source-zero unit row")


def main():
    check(sha(SOURCE) == PIN, "provenance", "audited source pin")
    external_status = verify_source_pin(LOW30, "7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8")
    if external_status is not None:
        check(external_status, "provenance", "original imported compiler pin")
    source_text = SOURCE.read_text()
    for phrase in ["epsilon_f=2s epsilon", "No reverse-OU endpoint join", "known source-zero unit Gaussian row",
                   "FULLY EXPANDED", "captured-Z first(E_G)<=A/2"]:
        check(phrase in source_text, "provenance", "required scope/ledger clause: " + phrase)
    symbolic_checks()
    quadrature_checks()
    finite_graph_checks()
    history_checks()
    linear_law_and_scaling_checks()
    result = {"status": "PASS", "source_sha256": sha(SOURCE), "assertions": sum(COUNTS.values()),
              "assertions_by_group": COUNTS, "details": DETAILS,
              "imports": [{"path": source_label(p), "sha256": sha(p)} for p in IMPORTS if p != LOW30],
              "external_source_sha256":{"external:LOW30":"7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8"},
              "publication_pin_verification":pin_report(),
              "limits": ["No author checker is imported.",
                         "The imported completed mean compilers are not executed here.",
                         "Finite diagnostics supplement the written continuum/source-admission proof.",
                         "No generic closure, private-root retention, or endpoint join is certified."]}
    (HERE / "independent_resummed_m3_checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "assertions": result["assertions"],
                      "assertions_by_group": COUNTS, "details": DETAILS}, indent=2))


if __name__ == "__main__":
    main()
