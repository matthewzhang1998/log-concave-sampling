#!/usr/bin/env python3
"""Independent finite diagnostics; no m3/current/law-compiler admission.

The proof is in AUDIT.md. This script does not import the author's checker.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'INDEPENDENT-NODE-MEAN-SOURCE-COVARIANCE-AND-CUT-TEST.md'
PIN = 'c463ea3ce5f2d7584dca346c00bbec345592e07e90481f10bd708048de56476c'
rng = np.random.default_rng(907314)
counts = {}
worst = {}


def check(condition, group):
    counts[group] = counts.get(group, 0) + 1
    if not condition:
        raise AssertionError(group)


def close(a, b, group, tol=2e-10):
    error = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    worst[group] = max(worst.get(group, 0.0), error)
    check(error <= tol, group)


def op(matrix):
    return float(np.linalg.norm(matrix, 2))


def rule(pairs):
    lo = rng.uniform(.015, .485, pairs)
    nodes = np.column_stack((lo, 1-lo)).ravel()
    mass = rng.uniform(.1, 1.0, pairs)
    mass /= mass.sum()
    weights = np.repeat(mass / 2, 2)
    return weights, nodes, np.sqrt(1-nodes**2)


check(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == PIN, 'frozen_source_pin')
gaussian_cases = 0
for d in (1, 2, 4):
    for pairs in (1, 2, 4, 7):
        for rep in range(6):
            gaussian_cases += 1
            v, tau, c = rule(pairs)
            n = len(v)
            q = float(rng.uniform(.11, .97))
            xrow = np.r_[q, np.zeros(n)]
            yrow = np.column_stack((q*tau, np.diag(c)))
            # D-dimensional matrices check the actual ND root tape.
            X = np.kron(xrow.reshape(1, -1), np.eye(d))
            Y = np.kron(yrow, np.eye(d))
            close(v.sum(), 1., 'rule_mass')
            close(v@tau, .5, 'rule_first_moment')
            check(Y.shape == (n*d, (n+1)*d), 'literal_gaussian_dimension')
            selected = rng.random(n) < .5
            score_row = np.r_[1/q, -(selected*tau/c)]
            score = np.kron(score_row.reshape(1, -1), np.eye(d))
            directional = Y@score.T
            expected = np.kron((tau*(~selected)).reshape(-1, 1), np.eye(d))
            close(directional, expected, 'marked_annihilation_unmarked_motion')
            sigma_s = 1/q**2 + np.sum((tau[selected]/c[selected])**2)
            close(score@score.T, sigma_s*np.eye(d), 'subset_score_covariance')
            close(X@score.T, np.eye(d), 'score_changes_x_by_identity')
            all_score_row = np.r_[1/q, -tau/c]
            all_score = np.kron(all_score_row.reshape(1, -1), np.eye(d))
            precision = float(all_score_row@all_score_row)
            residual = X - (X@Y.T) @ np.linalg.solve(Y@Y.T, Y)
            close(residual, all_score/precision, 'posterior_residual_score')
            close(residual@residual.T, np.eye(d)/precision, 'posterior_schur_covariance')
            close(all_score@Y.T, np.zeros((d,n*d)), 'all_site_independence')
            check(np.linalg.matrix_rank(Y) == n*d, 'all_sites_leave_D_innovations')
            check(np.linalg.matrix_rank(np.vstack((X,Y))) == (n+1)*d,
                  'x_and_all_sites_determine_roots')
            Q, _ = np.linalg.qr(rng.normal(size=(d,d)))
            B = (Q*rng.uniform(.05,.35,d))@Q.T
            I = sum(v[j]*B@Y[j*d:(j+1)*d] for j in range(n))
            terminal = B@(X-I)
            check(np.linalg.matrix_rank(np.vstack((Y,terminal))) == (n+1)*d,
                  'all_sites_and_invertible_terminal_determine_roots')
            # Fixing x, the independent root derivative and exact covariance.
            root_matrix = np.hstack([v[j]*c[j]*B for j in range(n)])
            beta2 = float(np.sum((v*c)**2))
            close(root_matrix@root_matrix.T, beta2*(B@B), 'independent_quadratic_covariance')
            common = (v@c)*B
            close(common@common.T, float(v@c)**2*(B@B), 'common_quadratic_covariance')
            check(beta2 <= .75+1e-12, 'beta2_bound')
            check(beta2 < float(v@c)**2, 'common_covariance_not_importable')
            check(v[~selected]@tau[~selected] <= .5+1e-12, 'unmarked_feedback_bound')


nonlinear_cases = 0
for d in (1, 2, 4):
    for inner_pairs in (1, 2, 4):
        for outer_pairs in (1, 3):
            for rep in range(4):
                nonlinear_cases += 1
                A = float(rng.uniform(.01,.5))
                v, tau, c = rule(inner_pairs)
                w, r, q = rule(outer_pairs)
                n, no = len(v), len(w)
                beta2 = float(np.linalg.norm(v*c))
                beta_out = float(w@q)
                basis, _ = np.linalg.qr(rng.normal(size=(d,d)))
                K = (basis*rng.uniform(.08,.30,d))@basis.T
                dirs = rng.normal(size=(d+2,d))
                dirs /= np.linalg.norm(dirs,axis=1)[:,None]
                coeff = np.full(d+2,.55/(d+2))
                value_calls = [0]

                def g(x):
                    value_calls[0] += 1
                    return A*(K@x + (coeff*np.tanh(dirs@x))@dirs)

                def dg(x):
                    z = dirs@x
                    return A*(K + dirs.T@((coeff*(1-np.tanh(z)**2))[:,None]*dirs))

                def evaluate(z, private, with_baseline=True, derivatives=False):
                    G = private[:d]
                    H = private[d:].reshape(n,d)
                    F = np.zeros(d)
                    baseline = np.zeros(d)
                    Eprime = np.zeros((d,(n+1)*d))
                    Ez = np.zeros((d,d))
                    I_bounds = []
                    for i in range(no):
                        x = r[i]*z + q[i]*G
                        sites = tau[:,None]*x + c[:,None]*H
                        inn = sum(v[j]*g(sites[j]) for j in range(n))
                        force = g(x-inn)
                        F += w[i]*force
                        if with_baseline:
                            baseline += w[i]*g(x)
                        if derivatives:
                            L = sum(v[j]*tau[j]*dg(sites[j]) for j in range(n))
                            M = np.hstack([v[j]*c[j]*dg(sites[j]) for j in range(n)])
                            T = dg(x-inn)
                            delta = T-dg(x)-T@L
                            Eprime[:,:d] += w[i]*q[i]*delta
                            Eprime[:,d:] -= w[i]*T@M
                            Ez += w[i]*r[i]*delta
                            I_bounds.append((L,M,inn,x,sites))
                    return F, baseline, Eprime, Ez, I_bounds

                z = rng.normal(size=d)
                private = rng.normal(size=(n+1)*d)
                value_calls[0] = 0
                F,B,J,Jz,local = evaluate(z,private,True,True)
                check(value_calls[0] == no*(n+2), 'raw_correction_VALUE_count')
                value_calls[0] = 0
                evaluate(z,private,False)
                check(value_calls[0] == no*(n+1), 'raw_source_VALUE_count')
                check(private.size == (n+1)*d, 'outer_actual_private_dimension')
                for L,M,inn,x,sites in local:
                    check(np.linalg.eigvalsh(L).min() >= -1e-12, 'inner_x_first_PSD')
                    check(op(L) <= A/2+1e-12, 'inner_x_first_upper')
                    check(op(M) <= A*beta2+1e-12, 'inner_private_first_upper')
                    check(np.linalg.norm(inn) <= A*sum(v[j]*np.linalg.norm(sites[j]) for j in range(n))+1e-12,
                          'positive_sum_energy_pointwise')
                check(op(J) <= A*(beta_out+A)+1e-12, 'correction_full_first')
                check(op(Jz) <= A/2+A*A/4+1e-12, 'captured_Z_first')
                lift = np.zeros(((n+1)*d,(n+1)*d))
                lift[:d] = J
                check(op(lift-lift.T) <= A*A*(beta_out+beta2)+1e-12, 'full_private_curl')
                check(np.linalg.norm(F-B) <= A*sum(w[i]*np.linalg.norm(local[i][2]) for i in range(no))+1e-12,
                      'correction_energy_pointwise')
                direction = rng.normal(size=private.size)
                direction /= np.linalg.norm(direction)
                eps = 2e-5
                plus = evaluate(z,private+eps*direction)
                minus = evaluate(z,private-eps*direction)
                fd = ((plus[0]-plus[1])-(minus[0]-minus[1]))/(2*eps)
                close(fd,J@direction,'first_sweep_finite_difference',tol=1e-8)
                F0,B0,*_ = evaluate(np.zeros(d),np.zeros_like(private))
                close(F0,np.zeros(d),'literal_global_source_zero')
                close(B0,np.zeros(d),'literal_global_baseline_zero')
                Fo,Bo,*_ = evaluate(z,np.zeros_like(private))
                check(np.linalg.norm(Fo-Bo) <= A*A*np.linalg.norm(z)/4+1e-12,
                      'literal_fixed_Z_origin_bound')


# Endpoint stress: one near-endpoint observation has unweighted precision
# cost. A symmetric two-node rule preserves mean 1/2.
endpoint_rows = []
for epsilon in (1e-2,1e-4,1e-6,1e-8):
    tau = np.array([epsilon,1-epsilon])
    c = np.sqrt(1-tau*tau)
    q = .7
    precision = float(q**-2 + np.sum((tau/c)**2))
    check(precision > .45/epsilon, 'algebraic_all_site_precision')
    endpoint_rows.append({'epsilon':epsilon,'all_site_precision':precision,
                          'innovation_variance':1/precision})

check(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == PIN, 'frozen_source_pin')
result = {
    'status':'PASS: bounded independent-node source/covariance/cut diagnostics only',
    'source_sha256':PIN,
    'seed':907314,
    'gaussian_configurations':gaussian_cases,
    'nonlinear_configurations':nonlinear_cases,
    'assertions':sum(counts.values()),
    'assertions_by_group':counts,
    'largest_identity_residual_by_group':worst,
    'endpoint_precision_stress':endpoint_rows,
    'scope_exclusions':['No m3 current closure','No current-to-law compiler',
                        'No full endpoint theorem','No author covariance receipt reused'],
}
(HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in
                 ('assertions_by_group','largest_identity_residual_by_group')},indent=2))
