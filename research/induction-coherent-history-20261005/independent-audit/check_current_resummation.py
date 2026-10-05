"""Deterministic independent checks; no sampled result substitutes for proof."""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.hermite_e import hermegauss, hermeval

OUT = Path(__file__).resolve().parent
nodes, weights = hermegauss(96)
weights = weights / math.sqrt(2 * math.pi)
z1, z2 = np.meshgrid(nodes, nodes, indexing="ij")
w2 = weights[:, None] * weights[None, :]
z = np.stack([z1, z2], axis=-1)
h = 1 / 8
e = 2 * h + h * h
cases = [(0, 0), (e, e), (-2*h+h*h, e), (0, e), (0.1, -0.13)]
records = []

def he(n, x):
    coef = np.zeros(n + 1)
    coef[n] = 1
    return hermeval(x, coef)

for lam in cases:
    lam = np.array(lam)
    log_l = -0.5*np.log1p(lam).sum() + (z*z*(lam/(1+lam))).sum(axis=-1)/2
    likelihood = np.exp(log_l)
    moment_errors = []
    for p in (1, 2, 3):
        exact = np.prod((1+lam)**(-(p-1)/2)*(1-(p-1)*lam)**(-0.5))
        numeric = np.sum(w2*likelihood**p)
        moment_errors.append(float(abs(numeric-exact)))
        assert abs(numeric-exact) < 1e-11
    partial = np.zeros_like(z1)
    norm_sum = 0.0
    truncations = []
    for k in range(8):
        term = np.zeros_like(z1)
        coefficient = 0.0
        for a in range(k+1):
            b = k-a
            term += lam[0]**a*lam[1]**b*he(2*a,z1)*he(2*b,z2)/(2**k*math.factorial(a)*math.factorial(b))
            coefficient += math.comb(2*a,a)*math.comb(2*b,b)*lam[0]**(2*a)*lam[1]**(2*b)/4**k
        norm_numeric = float(np.sum(w2*term*term))
        assert abs(norm_numeric-coefficient) < 1e-12
        assert abs(float(np.sum(w2*likelihood*term))-coefficient) < 1e-11
        partial += term
        norm_sum += coefficient
        tail_numeric = float(np.sum(w2*(likelihood-partial)**2))
        tail_exact = float(np.prod((1-lam*lam)**(-0.5))-norm_sum)
        assert abs(tail_numeric-tail_exact) < 1e-11
        original_bound_sq = 2*(2*e*e)**(k+1)
        sharper_bound_sq = e**(2*k+2)/(1-e*e)
        assert tail_numeric <= original_bound_sq + 1e-13
        assert tail_numeric <= sharper_bound_sq + 1e-13
        truncations.append({"degree_half": k, "tail_squared": tail_numeric,
                            "original_bound_squared": original_bound_sq,
                            "sharper_bound_squared": sharper_bound_sq})
    records.append({"eigenvalues": lam.tolist(), "moment_absolute_errors": moment_errors,
                    "truncations": truncations})

# Finite, genuinely correlated mark/source bank and same-G retained-mark coupling.
v = 4.0
atoms = np.array([[1.0, 0.0], [0.1, 0.7], [-0.4, 0.2]])
probs = np.array([0.4, 0.35, 0.25])
mu = probs @ atoms
C = np.outer(mu, mu)
cvals, cu = np.linalg.eigh(C)
sqrt_target = (cu*np.sqrt(v+cvals)) @ cu.T
fluctuation = 0.0
coupling_sq = 0.0
current_diffs = []
for i, f in enumerate(atoms):
    for j, fp in enumerate(atoms):
        mass = probs[i]*probs[j]
        H = (np.outer(f,fp)+np.outer(fp,f))/2
        K = np.eye(2)+H/(2*v)
        fluctuation += mass*np.linalg.norm(H-C,"fro")**2
        coupling_sq += mass*np.linalg.norm(math.sqrt(v)*K-sqrt_target,"fro")**2
        mark = (1+i)*(2-j) + np.dot(f,fp)
        X = math.sqrt(v)*z
        Z = math.sqrt(v)*np.einsum("ij,...j->...i",K,z)
        A = np.eye(2)-np.linalg.inv(K)@np.linalg.inv(K)
        likelihood = np.exp(-np.linalg.slogdet(K)[1]+np.einsum("...i,ij,...j->...",z,A,z)/2)
        phi_x = X[...,0]**2 + X[...,0]*X[...,1] + X[...,1]**4
        phi_z = Z[...,0]**2 + Z[...,0]*Z[...,1] + Z[...,1]**4
        diff = mass*mark*np.sum(w2*(phi_z-likelihood*phi_x))
        current_diffs.append(float(diff))
retained_bound = math.sqrt(fluctuation)/(2*math.sqrt(v))+np.linalg.norm(C@C,"fro")/(8*v**1.5)
assert math.sqrt(coupling_sq) <= retained_bound+1e-12
assert abs(sum(current_diffs)) < 1e-10

# Rotating rank-one range: cap of the current span does not control the first.
b = 0.1
a = 1-(1+b)**(-2)
rotation_derivatives = []
for transverse in (1,10,100,1000):
    g = np.array([1.0, float(transverse)])
    exact = a*g[0]*g[1]
    dt = 1e-7/max(1,transverse)
    def log_l(theta):
        u = np.array([math.cos(theta), math.sin(theta)])
        return -math.log1p(b)+a*np.dot(u,g)**2/2
    finite = (log_l(dt)-log_l(-dt))/(2*dt)
    assert abs(finite-exact) < 1e-5*max(1,abs(exact))
    rotation_derivatives.append({"retained_projection":1, "transverse":transverse,
                                 "log_l_derivative":exact})

# Derivative rank and operator bound for E=K^2-I.
rng = np.random.default_rng(1749)
rank_records = []
for _ in range(12):
    f,fp,df,dfp = rng.normal(size=(4,8))
    f /= max(1,np.linalg.norm(f)); fp /= max(1,np.linalg.norm(fp))
    H = (np.outer(f,fp)+np.outer(fp,f))/2
    dH = (np.outer(df,fp)+np.outer(fp,df)+np.outer(f,dfp)+np.outer(dfp,f))/2
    K = np.eye(8)+H/(2*v)
    dE = (K@dH+dH@K)/(2*v)
    rank = np.linalg.matrix_rank(dE, tol=1e-10)
    assert rank <= 4
    bound = (1+1/(2*v))*np.linalg.norm(dH,2)/v
    assert np.linalg.norm(dE,2) <= bound+1e-12
    rank_records.append(int(rank))

result = {"status":"pass", "guard_h":h, "guard_e":e,
          "quadrature_nodes_each_coordinate":len(nodes), "likelihood_cases":records,
          "retained_mark_same_G_cost":math.sqrt(coupling_sq),
          "retained_mark_bound":retained_bound,
          "correlated_mark_polynomial_identity_error":abs(sum(current_diffs)),
          "rotating_range_derivative_counterexample":rotation_derivatives,
          "DE_ranks":rank_records,
          "scope":"Deterministic diagnostics support the accompanying proofs; no general native/history admission is inferred."}
(OUT/"current_resummation_checks.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:v for k,v in result.items() if k not in ("likelihood_cases",)},indent=2))
