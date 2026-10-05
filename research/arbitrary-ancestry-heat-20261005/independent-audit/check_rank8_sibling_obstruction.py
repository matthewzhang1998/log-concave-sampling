#!/usr/bin/env python3
"""Independent deterministic algebra/PSD/OU tests plus seeded source Monte Carlo.

Run: python check_rank8_sibling_obstruction.py
No external services or native source oracles are used.
"""
from pathlib import Path
import importlib.util
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
GENERATOR = Path('/workspace/shared/rank-indexed-positive-returns-20261005/generate_rank_terms.py')
spec = importlib.util.spec_from_file_location('rank_generator', GENERATOR)
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)


def selected_history():
    coeff, ast = gen.histories(2)[0]
    stages = []
    for n, hit in zip(range(3, 9), [1, 1, 1, 2, 2, 2]):
        choices = gen.derivative(ast, n - 1)
        chosen = [a for a in choices if any(v == hit and n - 1 in es for v, es in gen.vertices(a))]
        assert len(chosen) == 1
        fresh = gen.relabel(gen.histories(1)[0][1], n - 1, n - 2)
        fresh_derivative, = gen.derivative(fresh, n - 1)
        ast = ('r', n, ('mul', chosen[0], fresh_derivative))
        coeff *= n
        if n <= 6:
            assert (coeff, ast) in gen.histories(n)
        stages.append({'rank': n, 'hit': hit, 'coefficient': coeff})
    return coeff, ast, stages


def inspect_ast(ast):
    rows = {}
    structural_indices = []
    def walk(node, structural_path=(), path='root'):
        if node[0] == 'r':
            if node[2][0] == 'g':
                leaf = node[2]
                rows[leaf[1]] = {'degree': len(leaf[2]), 'primitive_index': node[1],
                                'structural_path': list(structural_path), 'edges': list(leaf[2])}
                return
            structural_indices.append(node[1])
            walk(node[2], structural_path + (path,), path + '/child')
        elif node[0] == 'mul':
            walk(node[1], structural_path, path + '/left')
            walk(node[2], structural_path, path + '/right')
        else:
            raise AssertionError('Unwrapped force')
    walk(ast)
    assert sorted(r['degree'] for r in rows.values()) == [1] * 6 + [4, 4]
    assert [rows[i]['primitive_index'] for i in range(1, 9)] == [5, 5] + [2] * 6
    assert structural_indices == [8] * 7
    assert rows[1]['structural_path'] == rows[2]['structural_path']
    assert len(rows[1]['structural_path']) == 7
    assert len(gen.clocks(ast)) == 15
    return rows, structural_indices


def query_rows(ast, k, structural_r=0.5, leaf_r=0.5):
    """Exact coefficient-bank rows; original private shields remain averaged."""
    dim = 15
    count = 0
    rows, widths = {}, {}
    def walk(node, caller):
        nonlocal count
        if node[0] == 'r':
            noise = np.zeros(dim)
            noise[count] = 1
            count += 1
            if node[2][0] == 'g':
                vertex = node[2][1]
                r = math.sqrt(1 - 2 / k ** 2) if vertex in (1, 2) else leaf_r
                x = (1 - r ** 2) / 2
                rows[vertex] = r * caller + math.sqrt(x) * noise
                widths[vertex] = x
            else:
                walk(node[2], structural_r * caller + math.sqrt(1 - structural_r ** 2) * noise)
        elif node[0] == 'mul':
            walk(node[1], caller)
            walk(node[2], caller)
        else:
            raise AssertionError('Unwrapped force')
    walk(ast, np.zeros(dim))
    assert count == dim
    return np.vstack([rows[i] for i in range(1, 9)]), widths


def exact_variance(V):
    return 0.25 + math.exp(-4) * (1 + math.exp(-8 * V)) / 8 - math.exp(-2) * (1 + math.exp(-4 * V)) / 4


def integrate(f, lo=0., hi=1., n=100):
    t, w = leggauss(n)
    xs = lo + (hi - lo) * (t + 1) / 2
    return float((hi - lo) / 2 * np.dot(w, f(xs)))


def main():
    coeff, ast, stages = selected_history()
    rows, structural_indices = inspect_ast(ast)
    assert coeff == math.factorial(8)
    (HERE / 'rank8-selected-history.json').write_text(json.dumps({'coefficient': coeff, 'ast': ast, 'stages': stages, 'vertices': rows}, indent=2) + '\n')
    v = 1 - 0.5 ** 14
    psd_rows = []
    for k in [2., 8., 32., 128.]:
        x = 1 / k ** 2
        r2 = 1 - 2 * x
        c = v * r2
        G = c * np.ones((2, 2)) + x * np.eye(2)
        A, _ = query_rows(ast, k)
        assert np.max(np.abs(A[:2] @ A[:2].T - G)) < 2e-14
        for rho in [0., 0.3, 0.8, 0.99]:
            h2 = 1 - rho ** 2
            optimal_equal = h2 * G - h2 * x * np.eye(2)
            overdrawn = h2 * G - 1.001 * h2 * x * np.eye(2)
            assert np.linalg.eigvalsh(optimal_equal)[0] > -2e-15
            assert np.linalg.eigvalsh(overdrawn)[0] < -0.00099 * h2 * x
            one_shield_max = h2 * x * (2 * c + x) / (c + x)
            single = h2 * G - np.diag([one_shield_max, 0.])
            assert abs(np.linalg.eigvalsh(single)[0]) < 3e-15
            psd_rows.append({'k': k, 'rho': rho, 'total_extraction_cap': 2 * h2 * x,
                             'one_shield_max': one_shield_max})

    limit = 0.25 + math.exp(-4) / 8 - math.exp(-2) / 4
    low_variance = (1 - math.exp(-2)) ** 2 / 8
    low_integrand = lambda rho: 2 * rho * 0.25 * np.exp(-2 * (1 - rho ** 2)) * (1 - np.exp(-4 * rho ** 2))
    low_full = integrate(low_integrand)
    low_middle = integrate(low_integrand, 0.25, 0.75)
    assert abs(low_full - low_variance) < 2e-14
    assert low_middle > 0.01
    ou_checks = []
    for V in [0., 0.1, 1., 10.]:
        length2 = 4 * V + 2
        high_integrand = lambda rho: 2 * rho * (length2 / 8) * np.exp(-(1 - rho ** 2) * length2) * (1 - np.exp(-2 * rho ** 2 * length2))
        from_ou = low_full + integrate(high_integrand)
        assert abs(from_ou - exact_variance(V)) < 2e-13
        ou_checks.append({'V': V, 'closed_variance': exact_variance(V), 'symmetric_ou_integral': from_ou})

    # Actual 15-coordinate AST bank and all eight original-gradient factors.
    rng = np.random.default_rng(812031)
    bank = rng.standard_normal((220000, 15))
    a, b = 0.6, 0.3
    mc = []
    for k in [8., 16., 32., 64.]:
        A, widths = query_rows(ast, k)
        q = bank @ A.T
        leaves = np.prod(a + b * np.exp(-0.5 * k ** 2 * np.array([widths[i] for i in range(3, 9)])) * np.cos(k * q[:, 2:]), axis=1)
        F = np.sin(k * q[:, 0]) * np.sin(k * q[:, 1])
        # Center masses are x each in this exact asymptotic normalization fixture.
        # The Gauss panel table below checks literal positive quadrature weights.
        L_over_k2 = b ** 2 * math.exp(-1) * F * leaves
        normalized = float(np.var(L_over_k2) / (b ** 4 * a ** 12 * math.exp(-2)))
        V = k ** 2 * (1 - 2 / k ** 2) * v
        theoretical = exact_variance(V)
        assert abs(normalized - theoretical) < 0.003
        eps = 6 * b * math.exp(-3 * k ** 2 / 16)
        assert float(np.max(np.abs(leaves - a ** 6))) <= eps + 1e-14
        mc.append({'k': k, 'variance_over_b4_a12_exp_minus2_k4': normalized,
                   'exact_constant_leaves_value': theoretical, 'leaf_product_uniform_error_bound': eps})

    quadrature = []
    for order_kind in ['fixed_3', 'growing_j_plus_1']:
        for j in [6, 10, 14, 18, 22, 26]:
            n = 3 if order_kind == 'fixed_3' else j + 1
            z, w = leggauss(n)
            delta = 2. ** (-j)
            r = 1 - 1.5 * delta + 0.5 * delta * z
            base = 0.5 * delta * w
            effective = base * r ** 4
            pos = int(np.argmax(effective))
            r0, w0 = float(r[pos]), float(effective[pos])
            x = (1 - r0 ** 2) / 2
            k = 1 / math.sqrt(x)
            assert w0 >= (1 - 2 * delta) ** 4 * delta / n
            beta_over_other_mass = w0 ** 2 / x ** 3
            variance_over_other_constants = w0 ** 4 * k ** 12 * exact_variance(k ** 2 * r0 ** 2 * v)
            quadrature.append({'kind': order_kind, 'panel': j, 'nodes': n, 'r': r0, 'mass_R5': w0,
                               'k': k, 'mass_over_x': w0 / x, 'beta_over_other_mass': beta_over_other_mass,
                               'variance_over_other_constants': variance_over_other_constants,
                               'variance_divided_by_k4': variance_over_other_constants / k ** 4})
    for kind in ['fixed_3', 'growing_j_plus_1']:
        selected = [z for z in quadrature if z['kind'] == kind]
        assert selected[-1]['variance_over_other_constants'] > 1e7 * selected[0]['variance_over_other_constants']

    guards = [{'alpha': 1e-3, 'k': k, 'beta_for_mass_x_each': k ** 2,
               'blind_root_amplitude': 1e-3 * k ** 2, 'root_radius_over_alpha': k ** 2}
              for k in [8., 32., 128.]]
    report = {
        'verdict': 'PASS: legal rank-8 sibling history and source-qualified one-node covariance obstruction; no full shared-bank sum no-go.',
        'rank8_integer_coefficient': coeff,
        'primitive_indices': [rows[i]['primitive_index'] for i in range(1, 9)],
        'structural_indices': structural_indices,
        'covariance_psd_checks': psd_rows,
        'variance_limit': limit,
        'persistent_difference_mode_variance': low_variance,
        'persistent_difference_mode_bridge_middle_1_4_to_3_4': low_middle,
        'exact_symmetric_ou_checks': ou_checks,
        'seeded_full_ast_source_checks': mc,
        'literal_positive_gauss_panel_checks': quadrature,
        'blind_normalization_amplitude_diagnostic': guards,
        'caveats': [
            'Monte Carlo is a sanity check; exact formulas and the audit note supply the proof.',
            'An upper endpoint-mass envelope alone cannot give an individual-node lower bound.',
            'Fixed finite endpoint cutoffs do not admit a k-to-infinity sequence inside one fixed rule.',
            'Positivity of scalar clock weights alone does not prove nonnegative cross-node covariance.',
            'Actual finite native filters, readout shares, first paths and alternate amplitude distributions need their own admission.']}
    (HERE / 'rank8-sibling-audit-checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ['verdict', 'rank8_integer_coefficient', 'primitive_indices',
          'structural_indices', 'variance_limit', 'persistent_difference_mode_variance',
          'persistent_difference_mode_bridge_middle_1_4_to_3_4', 'seeded_full_ast_source_checks']}, indent=2))


if __name__ == '__main__':
    main()
