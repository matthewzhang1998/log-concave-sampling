#!/usr/bin/env python3
"""Independent diagnostics for the smoothing extension.

These numerical tests supplement, and do not replace, the analytic proof.
No source file from the sealed stage-two packet is modified.
"""
from __future__ import annotations
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SEALED = ROOT.parent / 'nonlinear-bridge-stage-two'
PIN = 'cd1111e28dc399643f27bcdd0f97dc735ada196f8a43f220211500bcb45be3ae'
rng = np.random.default_rng(610052026)
checks = 0

def check(condition, label):
    global checks
    checks += 1
    if not bool(condition):
        raise AssertionError(label)

def op(matrix):
    return np.linalg.norm(matrix, 2)

check(hashlib.sha256((SEALED / 'MANIFEST.json').read_bytes()).hexdigest() == PIN,
      'sealed manifest matches supplied pin')

# Exact Gaussian entropy identities, deliberately including nonsymmetric matrices.
entropy_max_identity_error = 0.0
entropy_min_bound_slack = float('inf')
entropy_max_nonsymmetry = 0.0
for d in [1, 2, 3, 7]:
    for _ in range(150):
        B = rng.normal(size=(d, d))
        B *= rng.uniform(.001, .85) / op(B)
        L = op(B)
        c = float(rng.uniform(.05, 1.0))
        m = rng.normal(size=d)
        v = rng.normal(size=d)
        T = np.eye(d) - B
        sign, logdet = np.linalg.slogdet(T)
        check(sign > 0, 'orientation for nonsymmetric near-identity map')
        shift = B @ m + v
        gaussian_kl = .5 * (np.sum(T*T) - d + np.dot(shift, shift)/c**2 - 2*logdet)
        expected_F2 = c*c*np.sum(B*B) + np.dot(shift, shift)
        transport_kl = expected_F2/(2*c*c) - np.trace(B) - logdet
        entropy_max_identity_error = max(entropy_max_identity_error, abs(gaussian_kl-transport_kl))
        check(abs(gaussian_kl-transport_kl) <= 2e-11*(1+abs(gaussian_kl)), 'KL orientation and trace sign')
        bound = expected_F2/(2*c*c) + d*L*L/(2*(1-L))
        entropy_min_bound_slack = min(entropy_min_bound_slack, bound-transport_kl)
        check(transport_kl <= bound + 1e-10, 'nonsymmetric log-determinant bound')
        check(transport_kl >= -1e-10, 'nonnegative affine Gaussian KL')
        entropy_max_nonsymmetry = max(entropy_max_nonsymmetry, op(B-B.T))

# Smooth noncommuting convex-gradient examples. Damping represents their exact
# anchored Gaussian smoothing, not an executed oracle.
class Source:
    def __init__(self, A, directions, frequencies, phases, epsilon=0.0):
        self.A = A
        self.directions = directions
        self.frequencies = frequencies
        self.phases = phases
        self.weights = np.full(len(frequencies), .24/len(frequencies))
        self.damping = np.exp(-.5*epsilon*epsilon*frequencies*frequencies)
    def value(self, x):
        phase = self.frequencies*(self.directions @ x)+self.phases
        c = self.weights*self.damping/self.frequencies
        return self.A*(.5*x + self.directions.T @ (c*(np.sin(phase)-np.sin(self.phases))))
    def jac(self, x):
        phase = self.frequencies*(self.directions @ x)+self.phases
        c = self.weights*self.damping*np.cos(phase)
        return self.A*(.5*np.eye(len(x)) + self.directions.T @ (c[:, None]*self.directions))


def conditional_pair(a, h, x, N, M, L1, L2):
    def f(s):
        return 2*s/math.sqrt(math.expm1(2*s))
    R = np.array([[f(a), f(a)*(a-1)],
                  [math.exp(-a)*f(h), math.exp(-a)*f(h)*(2*a+h-1)]])
    K = np.eye(2)-R@R.T
    det = np.linalg.det(K)
    Kr = (K+math.sqrt(det)*np.eye(2))/math.sqrt(np.trace(K)+2*math.sqrt(det))
    W = R@np.stack([N,M]) + Kr@np.stack([L1,L2])
    U = math.exp(-a)*x+math.sqrt(-math.expm1(-2*a))*W[0]
    V = math.exp(-h)*U+math.sqrt(-math.expm1(-2*h))*W[1]
    return U, V, math.exp(-a), math.exp(-a-h)


def maps(q, x, N, M, pair):
    U,V,rU,rV = pair
    dim = len(x)
    I = np.eye(dim)
    v = .5*x+.5*N
    w = .25*x+.5*N+.25*M
    gv,gw = q.value(v),q.value(w)
    b = q.value(gw)
    Jb = q.jac(gw)@q.jac(w)/4
    S = x-gv+b
    JS = I-q.jac(v)/2+Jb
    d1 = gv-q.value(U)
    d2 = gv-q.value(V)
    Jd1 = q.jac(v)/2-rU*q.jac(U)
    Jd2 = q.jac(v)/2-rV*q.jac(V)
    e = q.value(U)-q.value(U-q.value(V))-b
    Je = rU*(q.jac(U)-q.jac(U-q.value(V))) + rV*q.jac(U-q.value(V))@q.jac(V)-Jb
    return {
        'S': (S, JS, 1.0, 2+q.A),
        'single+': (S+d1, JS+Jd1, .5, 4+q.A),
        'single-': (S-d1, JS-Jd1, -.5, 4+q.A),
        'linear+': (S+e, JS+Je, .5, 5+3*q.A),
        'linear-': (S-e, JS-Je, -.5, 5+3*q.A),
        'quad++': (S+d1+d2, JS+Jd1+Jd2, .125, 6+q.A),
        'quad--': (S-d1-d2, JS-Jd1-Jd2, .125, 6+q.A),
        'quad+-': (S+d1-d2, JS+Jd1-Jd2, -.125, 6+q.A),
        'quad-+': (S-d1+d2, JS-Jd1+Jd2, -.125, 6+q.A),
    }

max_jacobian_fraction = 0.0
max_source_displacement_fraction = 0.0
max_terminal_nonsymmetry = 0.0
max_fd_error = 0.0
max_noncommuting_hessians = 0.0
for d in [2, 4]:
    directions = rng.normal(size=(3,d))
    directions /= np.linalg.norm(directions,axis=1)[:,None]
    frequencies = np.array([.7, 3.0, 17.0])
    phases = np.array([.2,1.1,-.4])
    for A in [.002, .013, 1/36, .2, .5]:
        g = Source(A,directions,frequencies,phases)
        h = Source(A,directions,frequencies,phases,epsilon=.19)
        delta = 2*A*np.sum(g.weights*(1-h.damping)/frequencies)
        L = 3.5*A+A*A/4
        Hx, Hy = g.jac(rng.normal(size=d)), g.jac(rng.normal(size=d))
        max_noncommuting_hessians = max(max_noncommuting_hessians, op(Hx@Hy-Hy@Hx))
        for iteration in range(100):
            x,N,M,L1,L2 = rng.normal(size=(5,d))
            a,hclock = rng.uniform(.002,3.0,size=2)
            pair = conditional_pair(a,hclock,x,N,M,L1,L2)
            gm = maps(g,x,N,M,pair)
            hm = maps(h,x,N,M,pair)
            check(abs(sum(v[2] for v in gm.values())-1)<1e-14, 'terminal signed mass one')
            check(abs(sum(abs(v[2]) for v in gm.values())-3.5)<1e-14, 'terminal absolute mass seven halves')
            for label,(T,JT,coeff,disp) in gm.items():
                ratio = op(JT-np.eye(d))/L
                max_jacobian_fraction = max(max_jacobian_fraction, ratio)
                check(ratio <= 1+1e-12, 'terminal near-identity defect '+label)
                max_terminal_nonsymmetry = max(max_terminal_nonsymmetry,op(JT-JT.T))
                diff = np.linalg.norm(T-hm[label][0])
                max_source_displacement_fraction = max(max_source_displacement_fraction,diff/(disp*delta))
                check(diff <= disp*delta+1e-12, 'coherent source-substitution displacement '+label)
                if iteration == 0:
                    step = 2e-6
                    fd = np.empty((d,d))
                    for k in range(d):
                        dx = np.eye(d)[k]*step
                        pplus = conditional_pair(a,hclock,x+dx,N,M,L1,L2)
                        pminus = conditional_pair(a,hclock,x-dx,N,M,L1,L2)
                        fd[:,k]=(maps(g,x+dx,N,M,pplus)[label][0]-maps(g,x-dx,N,M,pminus)[label][0])/(2*step)
                    err=op(fd-JT)
                    max_fd_error=max(max_fd_error,err)
                    check(err < 1e-7, 'noncommuting derivative formula '+label)
            weighted_shift = sum(abs(gm[k][2])*np.linalg.norm(h.value(gm[k][0])-h.value(hm[k][0])) for k in gm)
            check(weighted_shift <= A*delta*(14+5.5*A)+1e-12, 'weighted source-substitution constant')

# Algebraic exponents, explicit block/whole-space constants, and entropy endpoints.
A = 1/36
L_target=A/2+A*A/4
L_stencil=3.5*A+A*A/4
check(L_target<1 and L_stencil<1,'entropy maps invertible on guarded radius')
k=math.sqrt(3/8)
M=2+3/math.sqrt(2)+A*k
C_target=lambda d:1+A+math.sqrt(d)*((1+A)*math.pi/2+(L_target/A)/math.sqrt(1-L_target))
C_stencil_crude=lambda d:14+5.5*A+3.5*math.sqrt(d)*(math.pi*M/2+(L_stencil/A)/math.sqrt(1-L_stencil))
H_sum=2+11/(2*math.sqrt(2))+A*(1+4.5*k)
C_stencil=lambda d:14+5.5*A+math.sqrt(d)*(math.pi*H_sum/2+3.5*(L_stencil/A)/math.sqrt(1-L_stencil))
rational_H_bound=Fraction(14)+Fraction(11,72)+Fraction(11,7)*(Fraction(83,14)+Fraction(61,576))+Fraction(7,2)*Fraction(505,144)*Fraction(19,18)
check(rational_H_bound<37,'exact rational upper bound for sharp stencil constant')
d0=1+1/math.sqrt(2); e0=1+k
p_A=(d0*e0+A*e0**2/2+.5)/math.sqrt(2*math.pi)
q_A=((d0+A*e0)**3+4*d0**3+A**3*e0**3)/(6*math.sqrt(2))
check(p_A<1.32 and q_A<3,'displayed smoothing derivative coefficients')
F=Fraction
rA=F(1,36); rd=1+F(500,707); re=F(1613,1000)
rp=(rd*re+rA*re*re/2+F(1,2))/F(1253,500)
rq=((rd+rA*re)**3+4*rd**3+rA**3*re**3)/(6*F(707,500))
check(F(707,500)**2<2 and F(613,1000)**2>F(3,8),'rational radical bounds')
check(F(1253,500)**2<2*F(223,71),'rational denominator using Archimedes pi lower bound')
check(rp<F(33,25) and rq<3,'exact rational p_A and q_A upper certificates')
check(84+1.32*math.sqrt(3)+3*math.sqrt(15)<98,'conservative interpolation constant below 98')
check(max_noncommuting_hessians>1e-8,'source diagnostics include genuinely noncommuting Hessians')
for d in [1,2,16,1000]:
    check(C_target(d) < 5*math.sqrt(d),'simple target constant')
    check(C_stencil(d) < 37*math.sqrt(d),'sharp weighted stencil constant')
check(abs((2+2/3)-8/3)<1e-14,'A2 epsilon has grade 8/3')
check(abs((4-2*(2/3))-8/3)<1e-14,'A4 epsilon^-2 has grade 8/3')
check(abs((4-2/3)-10/3)<1e-14,'A4 epsilon^-1 has grade 10/3')
check(abs((1+1+2/3)-8/3)<1e-14,'coherent constant-tolerance source substitution grade')

result={
  'status':'PASS',
  'assertions':checks,
  'sealed_manifest_sha256':PIN,
  'entropy':{
    'maximum_affine_identity_error':entropy_max_identity_error,
    'minimum_upper_bound_slack':entropy_min_bound_slack,
    'maximum_tested_nonsymmetry':entropy_max_nonsymmetry,
  },
  'terminal_maps':{
    'maximum_fraction_of_declared_jacobian_bound':max_jacobian_fraction,
    'maximum_fraction_of_declared_source_displacement_bound':max_source_displacement_fraction,
    'maximum_terminal_jacobian_nonsymmetry':max_terminal_nonsymmetry,
    'maximum_source_Hessian_commutator':max_noncommuting_hessians,
    'maximum_finite_difference_error':max_fd_error,
  },
  'conservative_constants_at_A_1_over_36':{
    'L_target':L_target,'L_stencil':L_stencil,
    'C_target_b1':C_target(1),'C_stencil_b1':C_stencil(1),
    'C_target_b2':C_target(2),'C_stencil_b2':C_stencil(2),
    'C_stencil_crude_b1_for_comparison':C_stencil_crude(1),
    'exact_rational_stencil_upper_bound':str(rational_H_bound),
    'exact_rational_stencil_upper_bound_decimal':float(rational_H_bound),
    'p_A':p_A,'q_A':q_A,
    'final_main_coefficient_upper_bound':84+1.32*math.sqrt(3)+3*math.sqrt(15),
  },
  'scope':'Numerical diagnostics only. The proof is in the accompanying independent audit. No native compiler implementation or target-rate experiment is claimed.'
}
(ROOT/'independent_smoothing_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
