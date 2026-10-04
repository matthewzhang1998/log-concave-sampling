#!/usr/bin/env python3
"""Independent finite diagnostics, not an implementation of LOW30's mean compiler.

No author or earlier diagnostic code is imported. The exact-C2 fixture has
continuous non-Lipschitz Hessians and noncommuting ridge directions.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'THIRD-ORDER-CONDITIONAL-OU-VELOCITY-LAW-SERVICE.md'
SOURCE_PIN = '44294f14c4713882e656fbb835e03f90db61dfed224260a28d9501e39cf529c5'
LOW30_PIN = '7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--low30-source', type=Path, help='Optional local original LOW30 source for its archival hash check; not bundled.')
args = parser.parse_args()
LOW30 = args.low30_source
# Publication adaptation: all mathematical checks retained; one external hash
# assertion is omitted and explicitly reported unless original LOW30 is supplied.
rng = np.random.default_rng(620041027)
count = 0


def check(ok, what):
    global count
    count += 1
    if not bool(ok):
        raise AssertionError(what)


def op(x):
    return float(np.linalg.norm(x, 2))


def quad(alpha):
    # Error certificate 8*4^-m+2^(1-K) <= alpha.
    K = max(2, math.ceil(math.log2(4 / alpha)))
    m = max(2, math.ceil(math.log(16 / alpha, 4)))
    x, w = leggauss(m)
    rs, ws = [], []
    for j in range(K):
        a = 2.0 ** (-j - 1)
        rs.extend(1 - (1.5 * a + 0.5 * a * x))
        ws.extend(0.5 * a * w)
    rs.append(1 - 2.0 ** (-K - 1))
    ws.append(2.0 ** (-K))
    rs, ws = np.array(rs), np.array(ws)
    check(abs(ws.sum() - 1) < 3e-15, 'positive rule mass')
    check(abs(ws @ rs - .5) < 3e-15, 'positive rule first moment')
    check(np.all(ws > 0), 'positive rule weights')
    check(8 * 4.0 ** -m + 2.0 ** (1 - K) <= alpha * (1 + 1e-14), 'rule precision')
    return rs, np.sqrt(1 - rs * rs), ws


def cusp_h(t):
    return np.minimum(np.sqrt(np.abs(t)), 1.)


def cusp_g(t):
    u = np.abs(t)
    return np.sign(t) * np.where(u <= 1, (2 / 3) * u ** 1.5, u - 1 / 3)


class Fixture:
    def __init__(self, d):
        self.d = d
        rows = rng.normal(size=(d + 3, d))
        self.R = rows / np.linalg.norm(rows, axis=1)[:, None]
        self.shift = np.linspace(-.4, .5, d + 3)
        self.weight = .62 / op(self.R.T @ self.R)
        self.linear = .3 * rng.normal(size=d)

    def grad(self, x):
        return .17 * x + self.weight * (cusp_g(x @ self.R.T + self.shift) - cusp_g(self.shift)) @ self.R + self.linear

    def hess(self, x):
        h = cusp_h(x @ self.R.T + self.shift)
        return .17 * np.eye(self.d) + self.weight * np.einsum('...k,ki,kj->...ij', h, self.R, self.R)


def make_graph(pot, A, r, y, z, J, M, u, v, rule):
    d = pot.d
    eye = np.eye(d)
    s = math.sqrt(1 - r * r)
    rs, ss, ws = rule
    original_calls = 0

    def original(q):
        nonlocal original_calls
        original_calls += 1
        return pot.grad(q)

    b = y.copy()
    by = eye.copy()
    for _ in range(J):
        by = eye - A * pot.hess(b) @ by
        b = y - A * original(b)
    anchor = original(b)

    def g(q):
        # Exact semantic same-site anchor reuse.
        if np.array_equal(q, np.zeros(d)):
            return np.zeros(d)
        return math.sqrt(A) * (original(b + math.sqrt(A) * q) - anchor)

    def H(q):
        return A * pot.hess(b + math.sqrt(A) * q)

    def gy(q):
        return math.sqrt(A) * (pot.hess(b + math.sqrt(A) * q) - pot.hess(b)) @ by

    x = r * z
    xz = r * eye
    xy = np.zeros((d, d))
    for _ in range(M):
        xz = r * eye - s * s * H(x) @ xz
        xy = -s * s * (gy(x) + H(x) @ xy)
        x = r * z - s * s * g(x)
    base_anchor = g(x)
    R = x - r * z + s * s * base_anchor
    nodes = x + s * (rs[:, None] * u + ss[:, None] * v)
    bar_values = []
    for q in nodes:
        bar_values.append(s * (base_anchor if np.array_equal(q, x) else g(q) ) - s * base_anchor)
    packet = ws @ np.array(bar_values)
    q = x + s * (u - packet)
    baseline_q = x + s * u
    F = base_anchor if np.array_equal(q, x) else g(q)
    B = base_anchor if np.array_equal(baseline_q, x) else g(baseline_q)
    E = F - B
    Hnodes = np.array([H(qi) for qi in nodes])
    Hu = s * s * np.einsum('n,n,nij->ij', ws, rs, Hnodes)
    Hv = s * s * np.einsum('n,n,nij->ij', ws, ss, Hnodes)
    Hz = s * (np.einsum('n,nij->ij', ws, Hnodes) - H(x)) @ xz
    Hy = s * (np.einsum('n,nij->ij', ws, Hnodes) - H(x)) @ xy
    Hy += s * (sum(wi * gy(qi) for wi, qi in zip(ws, nodes)) - gy(x))
    H1, H0 = H(q), H(baseline_q)
    Eu = s * (H1 - H0 - H1 @ Hu)
    Ev = -s * H1 @ Hv
    Ez = (H1 - H0) @ xz - s * H1 @ Hz
    Ey = gy(q) - gy(baseline_q) + (H1 - H0) @ xy - s * H1 @ Hy
    return dict(E=E, B=B, F=F, H=packet, x=x, R=R, xz=xz, xy=xy,
                Hu=Hu, Hv=Hv, Hz=Hz, Hy=Hy, Eu=Eu, Ev=Ev, Ez=Ez, Ey=Ey,
                B_u=s * H0, H1=H1, H0=H0, calls=original_calls,
                g_a=math.sqrt(A) * (pot.grad(b + math.sqrt(A) * r * z) - anchor))


check(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_PIN, 'repaired source pin')
if LOW30 is not None:
    check(hashlib.sha256(LOW30.read_bytes()).hexdigest() == LOW30_PIN, 'LOW30 pin')
port_rows = []
for d in [2, 5, 13]:
    pot = Fixture(d)
    h1, h2 = pot.hess(rng.normal(size=d)), pot.hess(rng.normal(size=d))
    check(op(h1 @ h2 - h2 @ h1) > 1e-5, 'noncommuting fixture')
    for A in [.5, .16, .035]:
      for r in [0., .4, .93, .999]:
        s = math.sqrt(1-r*r)
        alpha = A*s*s
        rule = quad(alpha)
        for M in [0, 4]:
            J = 3
            y, z, u, v = rng.normal(size=(4, d))
            z *= 3
            out = make_graph(pot, A, r, y, z, J, M, u, v, rule)
            Ep = np.hstack([out['Eu'], out['Ev']])
            lift = np.vstack([Ep, np.zeros((d, 2*d))])
            check(np.linalg.norm(out['E']) <= A*s*np.linalg.norm(out['H'])+3e-14, 'one-force-chord energy')
            check(np.linalg.norm(out['H']) <= alpha*np.sqrt(u@u+v@v)+3e-14, 'packet envelope')
            check(op(np.hstack([out['Hu'],out['Hv']])) <= alpha+3e-14, 'packet full first')
            check(op(out['Hu']-out['Hu'].T) < 3e-14, 'packet u symmetry')
            check(np.linalg.eigvalsh(out['Hu']).min() >= -3e-14 and op(out['Hu']) <= alpha/2+3e-14, 'packet positive u first')
            check(op(Ep) <= A*s*(1+alpha)+3e-14, 'full E first')
            check(op(lift-lift.T) <= 2*A*A*s**3+3e-14, 'square-lift curl')
            check(op(out['B_u']-out['B_u'].T) < 3e-14 and op(out['B_u']) <= A*s+3e-14, 'genuine B first')
            check(op(out['xz']) <= r/(1-alpha)+3e-14, 'conditional mode z first')
            check(np.linalg.norm(out['x']) <= r*np.linalg.norm(z)/(1-alpha)+3e-14, 'mode growth')
            check(np.linalg.norm(out['R']) <= s*s*alpha**M*np.linalg.norm(out['g_a'])+3e-14, 'conditional mode residual')
            check(op(out['Ez']) <= 3*A*r+3e-14, 'E exposed z caller')
            check(op(out['Ey']) <= 8*math.sqrt(A)+3e-14, 'E physical caller')
            expected = J+1+M+1+len(rule[0])+2
            # r=0 uses exact global anchored zero during the conditional mode.
            if r != 0:
                check(out['calls'] == expected, 'complete original VALUE census')
            else:
                check(out['calls'] <= expected, 'zero-mode reuse census')
            zero = make_graph(pot,A,r,y,z,J,M,np.zeros(d),np.zeros(d),rule)
            check(np.array_equal(zero['E'],np.zeros(d)), 'literal exact E zero')
            # Directional numerical tests do not certify a modulus of Hessians.
            dirs = rng.normal(size=(4,d)); dirs /= np.linalg.norm(dirs,axis=1)[:,None]
            eps = 1e-6
            args = [y,z,u,v]
            for j,der in enumerate([out['Ey'],out['Ez'],out['Eu'],out['Ev']]):
                plus=[a.copy() for a in args]; minus=[a.copy() for a in args]
                plus[j]+=eps*dirs[j]; minus[j]-=eps*dirs[j]
                p=make_graph(pot,A,r,plus[0],plus[1],J,M,plus[2],plus[3],rule)['E']
                m=make_graph(pot,A,r,minus[0],minus[1],J,M,minus[2],minus[3],rule)['E']
                fd=(p-m)/(2*eps)
                check(np.linalg.norm(fd-der@dirs[j]) <= 2e-7, 'actual first finite difference')
            port_rows.append(dict(d=d,A=A,r=r,M=M,queries=out['calls'],nodes=len(rule[0]),
                                  first_over_As=op(Ep)/(A*s),curl_over_A2s3=op(lift-lift.T)/(A*A*s**3)))


# Independent exact-target integration for a smooth asymmetric scalar fixture.
def sg(x,A): return A*(.45*x+.3*(1-np.cos(x)))
def sh(x,A): return A*(.45+.3*np.sin(x))
def su(x,A): return A*(.225*x*x+.3*(x-np.sin(x)))


def gauss(n):
    x,w=hermgauss(n)
    return math.sqrt(2)*x,w/math.sqrt(math.pi)


def scalar_means(A,r,z,M,n,rule):
    s=math.sqrt(1-r*r); x=r*z
    for _ in range(M): x=r*z-s*s*sg(x,A)
    gx=sg(x,A); R=x-r*z+s*s*gx
    t,w=gauss(n); u=t[:,None]; v=t[None,:]
    H=np.zeros((n,n))
    rs,ss,ws=rule
    for ri,si,wi in zip(rs,ss,ws): H += wi*s*(sg(x+s*(ri*u+si*v),A)-gx)
    raw=sg(x+s*(u-H),A)
    mean=float(w@raw@w)
    q=r*z+s*t; logw=-su(q,A); logw-=logw.max(); ww=w*np.exp(logw); ww/=ww.sum()
    exact=float(ww@sg(q,A))
    deriv=r*(float(ww@sh(q,A))-float(ww@(sg(q,A)-exact)**2))
    return mean,exact,R,deriv


mean_rows=[]
for A in [.5,.2,.08,.03,.012]:
  for r in [0.,.5,.9,.99]:
    s=math.sqrt(1-r*r); rule=quad(A*s*s)
    for z in [-1.3,.7]:
      M=12
      lo=scalar_means(A,r,z,M,64,rule)
      hi=scalar_means(A,r,z,M,96,rule)
      floor=abs(hi[0]-lo[0])+abs(hi[1]-lo[1])+2e-14
      error=abs(hi[0]-hi[1]); scale=A**3*s**5
      check(error <= 10*scale+A*abs(hi[2])+10*floor, 'true conditional mean force bound')
      eps=1e-5
      dp=scalar_means(A,r,z+eps,M,64,rule)[1]
      dm=scalar_means(A,r,z-eps,M,64,rule)[1]
      check(abs((dp-dm)/(2*eps)-hi[3]) < 1e-9, 'OU covariance sign and scaling')
      mean_rows.append(dict(A=A,r=r,z=z,mean_error=error,A3s5=scale,ratio=error/scale,residual=hi[2],resolution_floor=floor))

# Gaussian covariance integration by parts for the coefficient in equation (3).
t,w=gauss(64); u=t[:,None]; v=t[None,:]
tr,tw=leggauss(56); tr=(tr+1)/2; tw=tw/2
cov_rows=[]
for A,r,z in [(.3,0.,2.),(.3,.6,-1.4),(.1,.98,.8)]:
    s=math.sqrt(1-r*r); q=r*z+s*t
    U=su(q,A); G=sg(q,A); covariance=float(w@(U*G)-(w@U)*(w@G))
    integrated=0.
    for ti,wi in zip(tr,tw):
        prod=sh(r*z+s*u,A)*sg(r*z+s*(ti*u+math.sqrt(1-ti*ti)*v),A)
        integrated+=wi*float(w@prod@w)*s*s
    check(abs(covariance-integrated)<2e-12, 'second coefficient covariance identity')
    cov_rows.append(dict(A=A,r=r,z=z,covariance=covariance,integrated=integrated))

# The shared v cannot be independently replaced at each inner node.
rs,ss,ws=quad(.15); beta=float(ws@ss)
check(beta*beta > float((ws*ws)@(ss*ss))+.1, 'shared-root covariance negative control')

result=dict(status='PASS',checks=count,seed=620041027,
    scope='Finite raw VALUE graph, actual first/curl/caller/zero/query checks and exact-target quadrature diagnostics. Does not execute or numerically certify the imported LOW30 high-order mean compiler.',
    source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    low30_sha256=LOW30_PIN,
    low30_pin_status='verified_original_bytes' if LOW30 is not None else 'external_not_reverified',
    original_assertions=1685, external_pin_assertions_omitted=0 if LOW30 is not None else 1,
    port_cases=port_rows,mean_cases=mean_rows,covariance_cases=cov_rows,
    max_mean_bias_ratio=max(row['ratio'] for row in mean_rows),
    max_curl_ratio=max(row['curl_over_A2s3'] for row in port_rows))
(HERE/'conditional_velocity_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if not k.endswith('_cases')},indent=2))
