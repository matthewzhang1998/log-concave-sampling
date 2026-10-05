#!/usr/bin/env python3
"""Independent finite/algebra diagnostics; not a numerical OU-path proof."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.integrate import quad
import sympy as sp


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
PINS = {
    "RADIAL-COROLLARY-FOR-THE-UNSHIFTED-P3-MEAN-CANDIDATE.md": "bae3fff77a2cfe75e0888637b4d39a1346bb64c1eec1ef8217f3e151af1eecc0",
    "RADIAL-SHELL-REFUTES-UNSHIFTED-SINGLE-HISTORY-COVARIANCE.md": "2807f7307646cc16c79049319b574259874b09e731a13f686696bdd60a05b035",
    "independent-radial-mean-audit/INDEPENDENT-CANONICAL-RADIAL-MEAN-AUDIT.md": "ffc4b21539b0eadb5cd639b968befb59cc8d0ee3b57e268296b5576e233bc66d",
    "P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md": "87b61a6a4ca03f81f7d910cdaa83b3783cd4792b3a428b89c8f0f629062e4467",
    "../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md": "b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced",
    "../nested-mean-independent-audit/INDEPENDENT-NESTED-FORCE-MEAN-AUDIT.md": "c2fbd8b02e1d225226ef5cf822fafda31f3dfa05e505d747d0b8772cb9d2daf7",
    "independent-m3-gate-audit/INDEPENDENT-M3-RAW-SOURCE-AND-EXACT-CURRENT-AUDIT.md": "abb77d6c7bd425d254beb2d7465eb072ccc40a0be7c11851da05de0fa49767cf",
    "../../../ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md": "4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f",
}
count = 0
def check(test):
    global count
    assert bool(test)
    count += 1

for name, pin in PINS.items():
    check(verify_source_pin(BASE/name, pin))

# All-degree accuracy is the cited analytical certificate, not finite sampling.
def clock_rule(delta):
    K = math.ceil(math.log2(4/delta))
    m = math.ceil(math.log(16/delta, 4))
    q, w = np.polynomial.legendre.leggauss(m)
    nodes, weights = [], []
    for k in range(K):
        lo, hi = 2.**(-k-1), 2.**(-k)
        nodes.extend(1-((lo+hi)/2+(hi-lo)*q/2))
        weights.extend((hi-lo)*w/2)
    nodes.append(1-2.**(-K-1))
    weights.append(2.**(-K))
    n, w = np.array(nodes), np.array(weights)
    check(np.all(n > 0) and np.all(n < 1) and np.all(w > 0))
    check(abs(w.sum()-1) < 3e-15)
    check(abs(w@n-.5) < 3e-15)
    cert = 8*4.**(-m)+2.**(1-K)
    check(cert <= delta*(1+1e-14))
    for degree in [0,1,2,3,4,7,15,31,63,127,255,511,1023,4095]:
        check(abs(w@(n**degree)-1/(degree+1)) <= delta+1e-15)
    return n, w, cert

def eta(s):
    return math.exp(-1/(1-s*s/4)) if abs(s)<2 else 0.
norm = quad(eta, -2, 2, epsabs=1e-13)[0]
def b(s):
    return eta(s)/norm
def psi(s):
    if s <= -2: return 0.
    if s >= 2: return 1.
    return quad(eta, -2, s, epsabs=2e-12)[0]/norm
def f(rows, D):
    rows=np.asarray(rows)
    radii=np.linalg.norm(rows,axis=-1)
    values=np.array([psi(float(s)) for s in (radii-math.sqrt(D)).ravel()]).reshape(radii.shape)
    scale=np.divide(values,radii,out=np.zeros_like(values),where=radii>0)
    return rows*scale[...,None]

c=d=.25
rng=np.random.default_rng(510052026)
diagnostics=[]
largest_residual=0.
for D in [16,64,256]:
    A=1/math.sqrt(D)
    tau,v,cert=clock_rule(A*A)
    # Distinct innermost rule verifies the full literal two-level graph.
    sigma,u,inner_cert=clock_rule(A*A/2)
    out,ow,outer_cert=clock_rule(A**3)
    beta=float(v@np.sqrt(1-tau*tau))
    check(beta >= 2/3-A*A-1e-14)
    check(beta*beta-.25 >= 13/144)
    check(abs((d*c*c/2)*(.25-beta*beta)) >= 13/18432)
    def g(rows): return c*A*rows+d*A*f(rows,D)
    for _ in range(2):
        x,H,J=rng.standard_normal((3,D))
        y=tau[:,None]*x+np.sqrt(1-tau*tau)[:,None]*H
        z=sigma[None,:,None]*y[:,None,:]+np.sqrt(1-sigma*sigma)[None,:,None]*J
        L=np.einsum('k,jkd->jd',u,g(z))
        F=np.einsum('j,jd->d',v,g(y-L))
        K=-c*np.einsum('j,jd->d',v,L)+d*np.einsum('j,jd->d',v,f(y-L,D))
        predicted=c*A*(x/2+beta*H)+A*K
        residual=float(np.linalg.norm(F-predicted))
        largest_residual=max(largest_residual,residual)
        check(residual < 3e-14)
        # Quadratic-only means cancel exactly using both first moments.
        means=c*A*(.5*x)-c*c*A*A*.25*x
        direct=c*A*np.einsum('j,jd->d',v,tau[:,None]*x-c*A*.5*tau[:,None]*x)
        check(np.linalg.norm(means-direct)<3e-14)
    diagnostics.append(dict(D=D,A=A,mid_nodes=len(v),inner_nodes=len(u),outer_nodes=len(ow),beta=beta,
        mid_certificate=cert,inner_certificate=inner_cert,outer_certificate=outer_cert))

# Independent normalization and exponent identities.
t=sp.symbols('t',nonnegative=True)
check(sp.integrate(sp.exp(-2*t),(t,0,sp.oo))==sp.Rational(1,2))
check(2*sp.integrate(t*sp.exp(-2*t),(t,0,sp.oo))-sp.Rational(1,4)==sp.Rational(1,4))
A=sp.symbols('A',positive=True)
D=A**-2
nested=(A**3+A**2*A+A**2*A**2)*sp.sqrt(D)
check(sp.simplify(nested-(2*A**2+A**3))==0)
check(sp.simplify(nested/A-(2*A+A**2))==0)
check(sp.simplify(A**4*sp.sqrt(D)-A**3)==0)
check(sp.simplify(A**2/A**3-1/A)==0)
check(sp.simplify(A**3*A*sp.sqrt(D)-A**3)==0)

# Strict negativity is proved by layer cake in the audit; this is a diagnostic.
expect_shifted=quad(lambda s:b(s)*math.exp(-(s+c/2)**2)/math.sqrt(math.pi),-2,2,epsabs=1e-13,epsrel=1e-13)[0]
expect_base=quad(lambda s:b(s)*math.exp(-s*s)/math.sqrt(math.pi),-2,2,epsabs=1e-13,epsrel=1e-13)[0]
witness=expect_shifted-expect_base
check(witness<0)
check(abs(witness+0.00182950693274534)<1e-11)

result={"status":"PASS finite/algebra diagnostics", "assertions":count,"input_pins":PINS,
        "finite_rules":diagnostics,"largest_literal_F2_decomposition_residual":largest_residual,
        "gaussian_shell_mean_difference":witness,"linear_resolvent_witness_limit":witness/2,
        "limits":["No true-history simulation or formal interval certificate.",
        "Lp estimates, imported mean theorem, and uniform asymptotic separator are analytical.",
        "No centered/resummed consumer or full endpoint law is verified."]}
(HERE/'p3_radial_corollary_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['input_pins','finite_rules']},indent=2))
