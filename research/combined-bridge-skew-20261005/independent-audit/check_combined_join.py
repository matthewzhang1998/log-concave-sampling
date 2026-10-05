"""Independent finite checks of the combined ledger, never native execution."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
COUNT = 0
def check(value, description):
    global COUNT
    COUNT += 1
    if not value:
        raise AssertionError(description)

pins = {
    "/workspace/shared/smoothed-bridge-prefix-20261005/MANIFEST.json": "2fc940544e410fc799814f23f71f083486b8c4b543a08ab4fdc8e8f101f93d8e",
    "/workspace/shared/law-only-reentry-ceiling-20261005/MANIFEST.json": "106350c3f20d06f6eaae1fb9b46677e0ccc7bfe8d21ec5f3e62f72ffdbaf515a",
    "/workspace/shared/shrinking-buffer-skew-join-20261005/MANIFEST.json": "4711e1bdde8ab138140e40c44b8f7fb5214eb4d051095a5452eec33fabd62f56",
}
def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
verified = 0
verified_imports = 0
for name, expected in pins.items():
    path = Path(name)
    check(digest(path) == expected, f"manifest {name}")
    manifest = json.loads(path.read_text())
    files = manifest["files"]
    pairs = files.items() if isinstance(files, dict) else ((x["path"], x["sha256"]) for x in files)
    for member, expected_hash in pairs:
        member_path = Path(member)
        if not member_path.is_absolute():
            member_path = path.parent / member_path
        check(digest(member_path) == expected_hash, f"sealed file {member_path}")
        verified += 1
    imports = manifest.get("inputs", manifest.get("imported_sources", manifest.get("pinned_inputs", {})))
    import_pairs = imports.items() if isinstance(imports, dict) else ((x["path"],x["sha256"]) for x in imports)
    for member, expected_hash in import_pairs:
        check(digest(member) == expected_hash, f"pinned imported source {member}")
        verified_imports += 1

wexp, etaexp, hexp = F(4,5), F(7,5), F(17,10)
prefix = [2+hexp, 1+2*hexp, 3+wexp, 3+3*wexp-hexp, 4+wexp]
check(prefix == [F(37,10), F(22,5), F(19,5), F(37,10), F(24,5)], "prefix exponents")
raw = [3+etaexp/2, 3+etaexp-wexp/2, 4+etaexp]
check(raw == [F(37,10), F(4), F(27,5)], "RAW exponents")
terms = [(5,F(3,2)), (5,F(1)), (6,F(7,6)), (6,F(2)), (7,F(5,2)), (7,F(3,2))]
expected = [(F(19,5),F(43,10)), (F(21,5),F(5)), (F(76,15),F(173,30)),
            (F(22,5),F(23,5)), (F(5),F(49,10)), (F(29,5),F(63,10))]
pairs = [(F(k)-p*wexp, F(k)-(p-1)*etaexp) for k,p in terms]
check(pairs == expected, "integrated service exponents")
check(min(x for pair in pairs for x in pair) == F(19,5), "native minimum")
guards = [1-etaexp/2, F(3,2)-etaexp, 1-etaexp/3, 2-etaexp, 3-F(3,2)*etaexp]
check(guards == [F(3,10), F(1,10), F(8,15), F(3,5), F(9,10)], "guard powers")
check(F(19,5)-F(37,10) == F(1,10), "absorption margin")

rng = np.random.default_rng(510053710)
max_bridge_cov = max_regression = max_scale = max_carrier = 0.0
for _ in range(500):
    A = 10**rng.uniform(-8,-1)
    w, eta, h = A**.8, A**1.4, A**1.7
    q, r = 1-w, np.sqrt(1-h*h)
    delta = -np.log1p(-w)
    early = rng.uniform(.001,.999)*delta
    future = delta+rng.uniform(.01,2)
    times = np.array([0,early,delta,future])
    cov = np.exp(-np.abs(times[:,None]-times[None,:]))
    e = [0,2]
    conditional = cov[1,3]-cov[1,e]@np.linalg.solve(cov[np.ix_(e,e)],cov[e,3])
    max_bridge_cov = max(max_bridge_cov, abs(conditional))
    check(abs(conditional) < 1e-10, "bridge/future conditional orthogonality")
    t = rng.uniform(0,1-eta)
    v = (1-t*t)*(2*w-w*w)/(1-q*q*t*t)
    check(v >= eta/2, "minimum bulk variance")
    shares = np.array([v/8]*4+[v/2])
    check(abs(shares.sum()-v) < 1e-14, "exact variance shares")
    max_carrier = max(max_carrier, abs(r*r*(1-t*t)+h*h-(1-r*r*t*t)))
    u = v/8
    alpha = q*A/np.sqrt(u)
    n = int(rng.integers(2,50))
    cbar = 1/(2*np.sqrt(n))
    weight, shield = rng.uniform(.0001,1), rng.uniform(.0001,1)
    rho1 = -4*weight*alpha/(cbar**3*shield)
    physical = u**1.5*cbar**3*rho1*alpha**2*shield/A**3
    desired = -4*q**3*weight
    max_scale = max(max_scale, abs(physical-desired))
    check(abs(physical-desired) < 1e-11, "physical cubic leading coefficient")
    # Non-isotropic reference covariance with an isotropic packet component.
    D = 4
    M = rng.normal(size=(D,D))
    beta2 = rng.uniform(.1,.4)
    C = beta2*np.eye(D) + M@M.T + .5*np.eye(D)
    Cinv = np.linalg.inv(C)
    S = rng.normal(size=D)
    mean = beta2*Cinv@S
    conditional_cov = beta2*np.eye(D)-beta2**2*Cinv
    lhs = (np.outer(mean,mean)+conditional_cov)/beta2**2-np.eye(D)/beta2
    rhs = np.outer(Cinv@S,Cinv@S)-Cinv
    err = np.linalg.norm(lhs-rhs)
    max_regression = max(max_regression,err)
    check(err < 1e-11, "anisotropic conditional Hermite regression")

legx, legw = np.polynomial.legendre.leggauss(14)
finite_sums = []
for A in [.1,.03,.01,.001,1e-5,1e-8]:
    w, eta = A**.8, A**1.4
    q = 1-w
    early, late = A**3/64, np.log(64/A**3)
    cut, delay = -np.log1p(-eta), -np.log1p(-w)
    edges = [early]
    while edges[-1] < late:
        edges.append(min(2*edges[-1],late))
    edges = sorted(set(edges+[cut,delay]))
    tau, weight = [early/2,late+1], [-np.expm1(-early),np.exp(-late)]
    for left,right in zip(edges[:-1],edges[1:]):
        nodes = (left+right)/2+(right-left)*legx/2
        tau.extend(nodes)
        weight.extend((right-left)*legw*np.exp(-nodes)/2)
    tau, weight = np.array(tau), np.array(weight)
    weight /= weight.sum()
    invv = q*q/(2*w-w*w)+1/(-np.expm1(-2*tau))
    bulk = tau >= cut
    ratios = {}
    for p in [.5,1,7/6,1.5,2,2.5]:
        endpoint = 1 if p<1 else (np.log(1/eta) if p==1 else eta**(1-p))
        ratio = float(np.sum(weight[bulk]*invv[bulk]**p)/(w**(-p)+endpoint))
        ratios[str(p)] = ratio
        check(ratio < 4, f"finite bulk p={p}")
    near = float(np.sum(weight[~bulk]*invv[~bulk]**.5)/(eta/np.sqrt(w)+np.sqrt(eta)))
    check(near < 4, "finite RAW sum")
    # Set log factor to 1: diagnostics test only the A dependence of aggregation.
    aggregate = float(np.sum(weight[bulk]*(A+A**1.5*np.sqrt(invv[bulk])))+np.sum(weight[~bulk])*A)
    check(aggregate <= 4*A, "weighted graph first")
    finite_sums.append({"A":A,"nodes":len(tau),"bulk_ratios":ratios,
                        "near_ratio":near,"aggregate_first_over_A":aggregate/A})

report = {
    "scope":"Exact rational exponent ledger, sealed source hashes, finite Gaussian/variance/cubic identities and positive dyadic sums. No native compiler execution.",
    "reviewed_join_sha256":digest(ROOT.parent/"COMBINED-BRIDGE-SKEW-JOIN.md"),
    "manifest_pins":pins,"sealed_files_verified":verified,"imported_hashes_verified":verified_imports,
    "assertions":COUNT,"passed":True,
    "prefix_exponents":[str(x) for x in prefix],"raw_exponents":[str(x) for x in raw],
    "service_exponent_pairs":[[str(x) for x in pair] for pair in pairs],
    "guard_exponents":[str(x) for x in guards],
    "max_bridge_future_conditional_covariance":max_bridge_cov,
    "max_anisotropic_hermite_regression_error":max_regression,
    "max_cubic_coefficient_scaling_error":max_scale,
    "max_terminal_carrier_identity_error":max_carrier,
    "finite_outer_sums":finite_sums,
}
(ROOT/"combined-checks.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
