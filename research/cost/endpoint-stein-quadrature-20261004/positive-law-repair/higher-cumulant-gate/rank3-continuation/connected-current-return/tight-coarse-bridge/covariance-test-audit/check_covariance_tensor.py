#!/usr/bin/env python3
"""Finite-dimensional checks; the written audit supplies the proof."""
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'COVARIANCE-TEST-PORT-REDUCTION.md'
EXPECTED = 'e8424b630508dab9ed042f9275a1f414d662a011cfa2e66f240b3f573b4bcd9f'
TOL = 1e-10

def op(a):
    return float(np.linalg.svd(a, compute_uv=False)[0])

def flatten(t, left):
    left = tuple(left)
    right = tuple(i for i in range(t.ndim) if i not in left)
    n = int(np.prod([t.shape[i] for i in left], dtype=int))
    m = int(np.prod([t.shape[i] for i in right], dtype=int))
    return t.transpose(left + right).reshape(n, m)

def gram_norm(mats, side):
    if side == 'left':
        return op(np.mean([a @ a.T for a in mats], axis=0))
    return op(np.mean([a.T @ a for a in mats], axis=0))

def check_family(samples, label):
    r, d = samples.ndim - 1, samples.shape[1]
    n = len(samples)
    z = samples - samples.mean(axis=0)
    vec = z.reshape(n, -1)
    covariance = vec.T @ vec / n
    v2 = max(0.0, float(np.linalg.eigvalsh(covariance)[-1]))
    v = np.sqrt(v2)
    h2 = float(np.mean(np.sum(vec * vec, axis=1)))
    local_subsets = [tuple(i for i in range(r) if mask >> i & 1)
                     for mask in range(1 << r)]
    K = max((op(flatten(f, s)) for f in samples for s in local_subsets
             if 0 < len(s) < r), default=0.0)
    T = covariance.reshape((d,) * (2 * r))
    counts = {'both_split': 0, 'one_whole': 0, 'balanced_whole': 0}
    worst_ratio = 0.0
    for mask in range(1, (1 << (2 * r)) - 1):
        left = tuple(i for i in range(2 * r) if mask >> i & 1)
        sa = tuple(i for i in range(r) if i in left)
        sb = tuple(i for i in range(r) if r + i in left)
        A = [flatten(f, sa) for f in z]
        B = [flatten(f, sb) for f in z]
        actual = op(flatten(T, left))
        forward = np.sqrt(gram_norm(A, 'left') * gram_norm(B, 'right'))
        reverse = np.sqrt(gram_norm(A, 'right') * gram_norm(B, 'left'))
        assert actual <= min(forward, reverse) + TOL, (label, mask, actual, forward, reverse)
        split_count = int(0 < len(sa) < r) + int(0 < len(sb) < r)
        category = ['balanced_whole', 'one_whole', 'both_split'][split_count]
        counts[category] += 1
        stated = [v2, 2 * K * v, 4 * K * K][split_count]
        sharper = [v2, K * v, K * K][split_count]
        assert actual <= stated + TOL, (label, mask, actual, stated)
        assert actual <= sharper + TOL, (label, mask, actual, sharper)
        if sharper:
            worst_ratio = max(worst_ratio, actual / sharper)
    hs2 = float(np.sum(covariance * covariance))
    assert hs2 <= v2 * h2 + TOL
    if r >= 2:
        assert h2 <= K * K * d + TOL
    return {'label': label, 'rank': r, 'dimension': d, 'samples': n,
            'K': K, 'v_squared': v2, 'h_squared': h2,
            'T_HS_squared': hs2, 'cut_counts': counts,
            'largest_ratio_to_sharper_case_bound': worst_ratio}

sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
assert sha == EXPECTED, (sha, EXPECTED)
rng = np.random.default_rng(20261005)
checks = [check_family(rng.normal(size=(9,) + (2,) * r), f'random_rank_{r}')
          for r in range(1, 5)]
# A bounded-cut rank-three tensor with arbitrarily large balanced covariance.
d = 5
E = np.zeros((d, d, d))
for i in range(d):
    E[i, i, i] = 1
necessity = check_family(np.stack([E, -E]), 'diagonal_sign_balanced_necessity')
assert abs(necessity['K'] - 1) < TOL
assert abs(necessity['v_squared'] - d) < TOL
checks.append(necessity)
# Empty-cut orientation: covariance Gram has norm 1, scalar Gram has trace n.
m = 8
isotropic = np.concatenate([np.sqrt(m) * np.eye(m), -np.sqrt(m) * np.eye(m)])
columns = [x[:, None] for x in isotropic]
rows = [x[None, :] for x in isotropic]
orientation = {'whole_left_covariance_gram': gram_norm(columns, 'left'),
               'whole_left_scalar_gram': gram_norm(columns, 'right'),
               'whole_right_covariance_gram': gram_norm(rows, 'right'),
               'whole_right_scalar_gram': gram_norm(rows, 'left')}
assert abs(orientation['whole_left_covariance_gram'] - 1) < TOL
assert abs(orientation['whole_right_covariance_gram'] - 1) < TOL
assert abs(orientation['whole_left_scalar_gram'] - m) < TOL
assert abs(orientation['whole_right_scalar_gram'] - m) < TOL
result = {'status': 'PASS', 'source_sha256': sha, 'seed': 20261005,
          'numpy_version': np.__version__, 'tolerance': TOL,
          'checks': checks, 'empty_cut_orientation': orientation,
          'source_specific_variance': 'OPEN; not tested or established',
          'positive_grouped_producer': 'OPEN; not tested or established'}
(HERE / 'checks.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': result['status'], 'source_sha256': sha,
                  'families': len(checks),
                  'proper_cuts_checked': sum(sum(c['cut_counts'].values()) for c in checks),
                  'source_specific_variance': 'OPEN', 'positive_grouped_producer': 'OPEN'}, indent=2))
