#!/usr/bin/env python3
"""Independent finite fixtures; these do not instantiate the native pair compiler.

Run from any directory. Writes grouped_law_response_checks.json beside this file.
Polynomial symbols below are analytical fixtures, never covariance samplers.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
import sympy as sp


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
NODE = HERE.parent.parent
CHECKS: list[dict] = []


def check(name, condition, detail=None):
    assert bool(condition), name
    item = {"name": name, "pass": True}
    if detail is not None:
        item["detail"] = detail
    CHECKS.append(item)


def same(name, lhs, rhs):
    difference = sp.simplify(sp.expand(lhs - rhs))
    check(name, difference.is_zero_matrix if isinstance(difference, sp.MatrixBase) else difference == 0)


def gaussian_moment(n):
    return sp.Integer(0) if n % 2 else sp.Integer(sp.factorial2(n - 1))


def gaussian_expectation(poly, variables):
    result = sp.Integer(0)
    for powers, coefficient in sp.Poly(sp.expand(poly), *variables).terms():
        result += coefficient * sp.prod(gaussian_moment(p) for p in powers)
    return sp.expand(result)


def rotation_relation(poly, kappa, delta):
    # Odd kappa powers disappear after Gaussian integration in these fixtures.
    result = 0
    for (power,), coefficient in sp.Poly(sp.expand(poly), kappa).terms():
        assert power % 2 == 0
        result += coefficient * (1 - delta**2) ** (power // 2)
    return sp.expand(result)


# Independent moment generating function / tilted polynomial checks.
y, s, d, b, t, k = sp.symbols("y s delta b theta kappa", real=True)
q, r, g = sp.symbols("Q R G", real=True)
for degree in range(9):
    for sign in (-1, 1):
        lhs = gaussian_expectation((y + s*(k*q + sign*d*(r+b*t)))**degree, (q, r))
        lhs = rotation_relation(lhs, k, d)
        rhs = gaussian_expectation((y + sign*s*d*b*t + s*g)**degree, (g,))
        same(f"scalar_tilt_monomial_{degree}_sign_{sign}", lhs, rhs)

# Non-diagonal J: J*theta=(theta1+theta2, theta1-theta2)/2.
t1, t2 = sp.symbols("theta1 theta2", real=True)
beta1, beta2 = b*(t1+t2)/2, b*(t1-t2)/2
q1, q2, r1, r2, g1, g2 = sp.symbols("Q1 Q2 R1 R2 G1 G2", real=True)
y1, y2 = sp.symbols("y1 y2", real=True)
for p1, p2 in ((1, 1), (2, 1), (1, 3), (2, 2)):
    for sign in (-1, 1):
        lhs = gaussian_expectation(
            (y1+s*(k*q1+sign*d*(r1+beta1)))**p1
            * (y2+s*(k*q2+sign*d*(r2+beta2)))**p2,
            (q1, q2, r1, r2),
        )
        lhs = rotation_relation(lhs, k, d)
        rhs = gaussian_expectation(
            (y1+sign*s*d*beta1+s*g1)**p1
            * (y2+sign*s*d*beta2+s*g2)**p2,
            (g1, g2),
        )
        same(f"non_diagonal_public_projection_{p1}_{p2}_{sign}", lhs, rhs)

# Gaussian normalizer has been divided out in the preceding identities.
J = sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 2)],
               [sp.Rational(1, 2), -sp.Rational(1, 2)]])
same("public_tilt_normalizer", (J.T*sp.Matrix([t1,t2])).dot(J.T*sp.Matrix([t1,t2])),
     (t1*t1+t2*t2)/2)

# Odd marking: heat of q^7 gives nonzero first, third, fifth and seventh marks.
eps, alpha, h = sp.symbols("epsilon alpha h", real=True)
F = gaussian_expectation((y+s*g)**7, (g,))
response = sp.expand(eps*alpha**3*(F.subs(y,y+h*t)-F.subs(y,y-h*t))*t**3)
for order in (1, 3, 5, 7):
    same(f"odd_mark_rank_{3+order}", response.coeff(t, 3+order),
         2*eps*alpha**3*h**order*sp.diff(F,y,order)/sp.factorial(order))
for rank in (3, 5, 7, 9):
    same(f"absent_even_mark_rank_{rank}", response.coeff(t, rank), 0)
check("grade_four_first_mark", 3+1 == 4)
check("grade_six_third_mark", 3+3 == 6)
check("grade_eight_product_of_two_unrenormalized_marks", 3+3+1+1 == 8)
check("grade_six_inherited_epsilon_squared", 2*3 == 6)

# Law composition: the exact epsilon^2 coefficient is E D + Var(L)/2.
e, L1, L2, D1, D2, w = sp.symbols("e L1 L2 D1 D2 w")
mixture = w*sp.exp(e*L1+e**2*D1)+(1-w)*sp.exp(e*L2+e**2*D2)
coef = sp.diff(sp.log(mixture),e,2).subs(e,0)/2
mean = w*L1+(1-w)*L2
same("connected_second_coefficient", coef,
     w*D1+(1-w)*D2+(w*L1**2+(1-w)*L2**2-mean**2)/2)

# The two signs cancel their entire cubic coefficient at delta=0.
qplus, qminus = k*q+d*r, k*q-d*r
signed = sp.expand(qplus**3-qminus**3)
variance = gaussian_expectation(signed**2,(q,r)) - gaussian_expectation(signed,(q,r))**2
variance = rotation_relation(variance,k,d)
same("whole_signed_cubic_variance", variance, 108*d**2-144*d**4+96*d**6)
same("whole_signed_cubic_variance_zero_at_delta_zero", variance.subs(d,0),0)
same("signed_variance_leading_delta_squared", sp.limit(variance/d**2,d,0),108)
check("signed_variance_not_sum_of_self_variances", sp.simplify(variance-30) != 0)
bounded_signed_variance=2*sp.exp(-1)*(sp.sinh(1)-sp.sinh(1-2*d*d))
same("bounded_sine_signed_variance_zero",bounded_signed_variance.subs(d,0),0)
same("bounded_sine_signed_variance_leading_delta_squared",
     sp.limit(bounded_signed_variance/d**2,d,0),4*sp.exp(-1)*sp.cosh(1))

# Direct public tilt creates higher physical ranks, even from rank-four input.
# E cos^2(R+b theta)=1/2*(1+exp(-2)*cos(2b theta)).
tilted_quartic = (1+sp.exp(-2)*sp.cos(2*b*t))*t**4/2
same("tilted_rank_four_can_feed_rank_six",
     sp.diff(tilted_quartic,t,6).subs(t,0)/sp.factorial(6),-sp.exp(-2)*b**2)

# Restricted variance obstruction, with arbitrary signed amplitudes.
a1,a2,a3=sp.symbols("a1 a2 a3",real=True)
total=a1*q+a2*(q*q-1)+a3*(q**3-3*q)
var=gaussian_expectation(total**2,(q,))-gaussian_expectation(total,(q,))**2
same("restricted_signed_variance_is_psd",var,a1*a1+2*a2*a2+6*a3*a3)
same("negative_cross_keeps_self_terms",
     gaussian_expectation((q-2*(q+q*q-1))**2,(q,)),9)

# Exact source-zero covariance and conservative derivative/copy ledger.
nu,eta=sp.symbols("nu eta",positive=True)
carrier_coeff=sp.Matrix.hstack(alpha*sp.eye(2),alpha*sp.eye(2),b*J,eta*sp.eye(2))
root_cov=sp.diag(nu,nu,nu,nu,1,1,1,1)
expected=(2*alpha**2*nu+eta**2)*sp.eye(2)+b**2*J*J.T
same("source_zero_physical_covariance",carrier_coeff*root_cov*carrier_coeff.T,expected)
same("joint_rotation_injection_bounded",rotation_relation((2*s*k)**2+(2*s*d)**2,k,d),4*s*s)
check("declared_complete_copy_bill_arithmetic",all(2*n == n*2 for n in (1,2,3,7)),
      "Arithmetic only; full occurrence/anchor lineage is a manual audit, not a compiler fixture.")
check("root_bank_concatenation",5*11*7 == 385,
      "Seven independent five-clock nodes in D=11 own 385 coarse coordinates.")

# A finite dilation stencil cannot kill every higher odd monomial.
for n in range(1,5):
    nodes=[sp.Rational(j+1,2*(n+1)) for j in range(n)]
    higher=sp.Matrix([[node**(2*j+3) for node in nodes] for j in range(n)])
    check(f"higher_odd_vandermonde_invertible_{n}",higher.det()!=0)

# Compact radial twist: scalar chi supported in (1,4), i.e. 1<|u|<2.
def chi_and_prime(radius_squared):
    z=float(radius_squared)
    if not 1 < z < 4:
        return 0.0,0.0
    p=(z-1)*(4-z)
    c=math.exp(-1/p)
    return c,c*(5-2*z)/(p*p)

def twist_derivatives(u,caller,a):
    u=np.asarray(u,dtype=float)
    c,cp=chi_and_prime(u@u)
    angle=a*math.sin(caller)*c
    rot=np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
    turn=np.array([[0.,-1.],[1.,0.]])
    out=rot@u
    angle_gradient=2*a*math.sin(caller)*cp*u
    Du=rot+np.outer(turn@out,angle_gradient)
    Dcaller=a*math.cos(caller)*c*(turn@out)
    return out,Du,Dcaller

grid=np.linspace(1,4,4001)
cmax=max(chi_and_prime(v)[0] for v in grid)
cpmax=max(abs(chi_and_prime(v)[1]) for v in grid)
u_bound=cmax+8*cpmax
caller_bound=2*cmax
twist_stats=[]
for a in (0.1,0.03,0.01,0.003):
    maxu=maxcaller=maxjoint=maxdet=0.0
    for radius in np.linspace(0,2.5,121):
        for angle in (0.,0.61,1.73):
            u=radius*np.array([math.cos(angle),math.sin(angle)])
            for caller in (0.,0.7,math.pi/2):
                out,Du,Dcaller=twist_derivatives(u,caller,a)
                check_radius=abs(float(out@out-u@u))
                assert check_radius < 1e-12
                maxu=max(maxu,float(np.linalg.norm(Du-np.eye(2),2))/a)
                maxcaller=max(maxcaller,float(np.linalg.norm(Dcaller))/a)
                maxjoint=max(maxjoint,float(np.linalg.norm(np.column_stack((Du-np.eye(2),Dcaller)),2))/a)
                maxdet=max(maxdet,abs(float(np.linalg.det(Du))-1))
    check(f"compact_twist_residual_root_first_{a}",maxu<=u_bound+1e-10)
    check(f"compact_twist_caller_first_{a}",maxcaller<=caller_bound+1e-10)
    check(f"compact_twist_joint_complete_first_{a}",maxjoint<=math.hypot(u_bound,caller_bound)+1e-10)
    check(f"compact_twist_volume_and_radius_{a}",maxdet<1e-12)
    twist_stats.append({"a":a,"root_first_over_a":maxu,"caller_first_over_a":maxcaller,
                        "complete_first_over_a":maxjoint,"determinant_error":maxdet})
response_limit=quad(lambda z:0.5*math.exp(-z/2)*z*chi_and_prime(z)[0],1,4,
                    epsabs=1e-13,epsrel=1e-13)[0]
response_a=quad(lambda z:0.5*math.exp(-z/2)*z*math.sin(0.003*chi_and_prime(z)[0])/0.003,1,4,
                epsabs=1e-13,epsrel=1e-13)[0]
check("twist_nonzero_carrier_retained_linear_response",response_limit>0.1,
      {"limit":response_limit,"at_a_0_003":response_a})
check("twist_response_is_order_a",abs(response_a-response_limit)<1e-6)

# Proper-cut / grouped-root derivative counterexample.
q0=0.37
cut_stats=[]
for dimension in (1,2,4,8,16):
    tensor=np.zeros((dimension,dimension,dimension))
    ii=np.arange(dimension)
    tensor[ii,ii,ii]=1
    norms=[]
    for axis in range(3):
        matrix=np.moveaxis(math.sin(q0)*tensor,axis,0).reshape(dimension,-1)
        norms.append(float(np.linalg.norm(matrix,2)))
    derivative_cut=float(np.linalg.norm(math.cos(q0)*tensor.reshape(1,-1),2))
    check(f"coefficient_all_proper_cuts_{dimension}",max(norms)<=1+1e-12)
    check(f"root_versus_all_physical_cut_{dimension}",
          abs(derivative_cut-math.sqrt(dimension)*abs(math.cos(q0)))<1e-12)
    # Covariance is Var(sin Q) times the outer square of vec(tensor).
    variance_sin=(1-math.exp(-2))/2
    covariance_norm=variance_sin*float(np.sum(tensor*tensor))
    check(f"covariance_three_versus_three_cut_{dimension}",
          abs(covariance_norm-dimension*variance_sin)<1e-12)
    cut_stats.append({"D":dimension,"proper_cut":max(norms),
                      "derivative_grouped_cut":derivative_cut,"covariance_cut":covariance_norm})

# Existing independent-reshielding check included for completeness.
sig2,u=sp.symbols("sigma_squared u",positive=True)
cov=sig2*sp.eye(3)+u*sp.ones(3)
conditional=sp.simplify(cov[0,0]-(cov[0,1:3]*cov[1:3,1:3].inv()*cov[1:3,0])[0])
same("common_heat_conditional_vertex_variance",conditional,sig2*(sig2+3*u)/(sig2+2*u))
same("common_heat_preserves_transverse_variance",(cov*sp.Matrix([1,-1,0]))[0],sig2)

# Fractional B mark: exact Gaussian pair, scalar fixture with p=1/2.
# P=M R+sqrt(1-M^2) V has variance one; Cov(P,R)=M.
M,B,A=sp.symbols("M B A",positive=True)
v,cross_root=sp.symbols("V cross_root",real=True)
for degree in range(5):
    lhs=gaussian_expectation(
        (y+s*(k*q+d*(M*(r+b*t)+cross_root*v)))**degree,(q,r,v))
    # Only even cross_root powers survive the V expectation.
    lhs=rotation_relation(lhs,cross_root,M)
    lhs=rotation_relation(lhs,k,d)
    rhs=gaussian_expectation((y+s*d*b*M*t+s*g)**degree,(g,))
    same(f"oriented_gaussian_pair_tilt_{degree}",lhs,rhs)
same("fractional_B_exact_shift",(sp.sqrt(A)*M).subs(M,B/sp.sqrt(A)),B)
same("fractional_B_pair_residual_path",sp.sqrt(A)*A*sp.sqrt(A),A**2)
same("fractional_B_gaussian_path",sp.sqrt(A)*A,A**sp.Rational(3,2))
for target in range(4,11):
    order=max(4,math.ceil(2*target-3))
    check(f"fractional_B_pair_order_target_{target}",sp.Rational(3,2)+sp.Rational(order,2)>=target)

# Exact Price identity, with shifted callers and a nontrivial heat scale.
tau,kr=sp.symbols("tau kr",real=True)
Gtau=tau*q+kr*r
price_cov=gaussian_expectation((y+s*q)**3*(y+s*Gtau)**3,(q,r))
price_cov=rotation_relation(price_cov,kr,tau)
price_derivative=9*s*s*gaussian_expectation((y+s*q)**2*(y+s*Gtau)**2,(q,r))
price_derivative=rotation_relation(price_derivative,kr,tau)
same("Price_identity_shifted_cubic",sp.diff(price_cov,tau),price_derivative)
price_integral=2*sp.integrate(price_derivative,(tau,1-2*d*d,1))
direct_difference=(y+s*(k*q+d*r))**3-(y+s*(k*q-d*r))**3
direct_square=rotation_relation(gaussian_expectation(direct_difference**2,(q,r)),k,d)
same("Price_difference_tensor_constant_and_interval",price_integral,direct_square)
same("Price_interval_length",2*s*s*(1-(1-2*d*d)),4*s*s*d*d)

# All-cut contraction of two four-tensors along the one root edge.
def all_cut_norms(tensor):
    ndim=tensor.ndim
    output=[]
    for mask in range(1,(1<<ndim)-1):
        left=[i for i in range(ndim) if mask & (1<<i)]
        right=[i for i in range(ndim) if i not in left]
        matrix=tensor.transpose(left+right).reshape(math.prod(tensor.shape[i] for i in left),-1)
        output.append(float(np.linalg.norm(matrix,2)))
    return output

rng=np.random.default_rng(20261004)
contraction_stats=[]
for case in range(6):
    left=rng.normal(size=(3,2,2,2))
    right=rng.normal(size=(3,2,2,2))
    left/=max(all_cut_norms(left))
    right/=max(all_cut_norms(right))
    joined=np.tensordot(left,right,axes=([0],[0]))
    joined_cut=max(all_cut_norms(joined))
    joined_hs=float(np.linalg.norm(joined))
    check(f"two_tensor_contraction_all_62_cuts_{case}",joined_cut<=1+1e-12)
    check(f"two_tensor_contraction_one_HS_factor_{case}",
          joined_hs<=min(np.linalg.norm(left),np.linalg.norm(right))+1e-12)
    contraction_stats.append({"case":case,"max_cut":joined_cut,"HS":joined_hs})
diagonal=np.zeros((3,3,3,3))
ii=np.arange(3)
diagonal[ii,ii,ii,ii]=1
joined=np.tensordot(diagonal,diagonal,axes=([0],[0]))
check("two_tensor_contraction_saturates_all_cut_bound",abs(max(all_cut_norms(joined))-1)<1e-12)
check("two_tensor_contraction_saturates_HS_bound",abs(np.linalg.norm(joined)-math.sqrt(3))<1e-12)

# Squared dyadic mass proxy: actual native quadrature imports the same envelope.
clock_stats=[]
for panels in (2,4,8,16,32,64):
    masses=2.0**(-np.arange(1,panels+1,dtype=float))
    sigmas=np.sqrt(masses-0.5*masses*masses)
    moments=[float(np.sum(masses*masses/(sigmas**power))) for power in range(5)]
    S0,S1,S2,S3,S4=moments
    envelope=2*S2*S2*S0+S0*S4*S0+4*S1*S3*S0+2*S1*S1*S2
    check(f"squared_clock_derivative_envelope_{panels}",envelope<10*panels)
    clock_stats.append({"panels":panels,"envelope":envelope,"envelope_per_panel":envelope/panels})

# Nine labeled hit pairs; these supplement the 130 named fixture groups.
strict_spine=[]
for left_hit,right_hit in itertools.product(range(3),repeat=2):
    left_ecc=max(left_hit,2-left_hit)
    right_ecc=max(right_hit,2-right_hit)
    joined_diameter=left_ecc+right_ecc+1
    left_B=max(2,left_ecc+1)
    right_B=max(2,right_ecc+1)
    assert joined_diameter>left_B and joined_diameter>right_B
    strict_spine.append({"left_hit":left_hit,"right_hit":right_hit,
                         "joined_diameter":joined_diameter,"left_B":left_B,"right_B":right_B})

sources=[
    HERE.parent/"WHOLE-KERNEL-GAUSSIAN-LAW-RESPONSE-AND-FAILED-BRIDGE-GATE.md",
    NODE/"POSITIVE-RANK3-PACKET-FEEDBACK-LEMMA.md",
    NODE/"EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md",
    NODE.parent/"order-reentry/independent-audit/INDEPENDENT-NATIVE-RANK3-PACKET-AUDIT.md",
    HERE.parent/"INTERMEDIATE-SCALE-B-MARK-COROLLARY.md",
    HERE.parent/"PAIRED-ROTATION-COVARIANCE-ALL-CUT-LEMMA.md",
]
verify_source_pin("external:LOW30","7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8")
result={
    "status":"PASS_FINITE_FIXTURES_ONLY",
    "scope":"Conditional leading-law marking and stated interface counterexamples; not a native pair compiler or generic bridge proof.",
    "versions":{"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__},
    "check_count":len(CHECKS),"checks":CHECKS,
    "radial_twist":twist_stats,"derivative_cuts":cut_stats,
    "tensor_contractions":contraction_stats,"dyadic_clock_envelopes":clock_stats,
    "strict_spine_hit_pairs":strict_spine,
    "source_sha256":{source_label(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    "external_source_sha256":{"external:LOW30":"7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8"},
    "publication_pin_verification":pin_report(),
}
output=HERE/"grouped_law_response_checks.json"
output.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],"checks":len(CHECKS),"output":str(output)},indent=2))
