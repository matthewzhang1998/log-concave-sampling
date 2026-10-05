#!/usr/bin/env python3
"""Independent deterministic checks. These supplement, rather than replace, the proof."""
import hashlib, itertools, json, math
from pathlib import Path
from fractions import Fraction
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.special import eval_hermitenorm
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
rng = np.random.default_rng(20261005)

def indices(d, n):
    return [a for a in itertools.product(range(n+1), repeat=d) if sum(a)<=n]

def phi(a, x):
    return math.prod(float(eval_hermitenorm(k,t))/math.sqrt(math.factorial(k)) for k,t in zip(a,x))

out = {'source_sha256': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.md')}}
# Direct multivariate kernel and gradient-sum identities, in small dimensions.
max_kernel_ratio=max_grad_ratio=max_identity_error=0.0
cases=0
for d in range(1,5):
    for m in range(1,7):
        aa=indices(d,m); bb=indices(d,m-1)
        for ell in (0.,.5,2.,5.):
            R=math.sqrt(d)+ell
            pts=[np.full(d,R/math.sqrt(d)), np.eye(d)[0]*R]
            for _ in range(3):
                v=rng.normal(size=d); pts.append(v/np.linalg.norm(v)*R*rng.random())
            for x in pts:
                kernel=sum(phi(a,x)**2 for a in aa)
                direct_grad=0.
                for a in aa:
                    for i,k in enumerate(a):
                        if k:
                            b=list(a); b[i]-=1
                            direct_grad += k*phi(b,x)**2
                identity=sum((d+sum(b))*phi(b,x)**2 for b in bb)
                max_identity_error=max(max_identity_error,abs(direct_grad-identity)/max(1,identity))
                base=1+d+R*R
                max_kernel_ratio=max(max_kernel_ratio,kernel/(4*base**m))
                max_grad_ratio=max(max_grad_ratio,direct_grad/(4*(d+m)*base**(m-1)))
                cases+=1
assert max_kernel_ratio<=1+1e-12 and max_grad_ratio<=1+1e-12
assert max_identity_error<1e-12
out['hermite_kernel']={'cases':cases,'max_bound_ratio':max_kernel_ratio,'max_gradient_bound_ratio':max_grad_ratio,'max_gradient_identity_relative_error':max_identity_error}

# Matrix-valued Parseval by exact Gaussian quadrature (degree 3 matrices).
d=3; degree=3; aa=indices(d,degree)
Cs=rng.normal(size=(len(aa),4,3))
Gram=np.einsum('aoi,aoj->ij',Cs,Cs)
xs,ws=hermegauss(degree+1); ws=ws/math.sqrt(2*math.pi)
quad=np.zeros((3,3)); EJop2=0.
for inds in itertools.product(range(len(xs)),repeat=d):
    x=np.array([xs[i] for i in inds]); weight=math.prod(ws[i] for i in inds)
    p=np.array([phi(a,x) for a in aa]); J=np.einsum('a,aoi->oi',p,Cs)
    quad+=weight*(J.T@J); EJop2+=weight*np.linalg.norm(J,2)**2
err=np.max(np.abs(quad-Gram))
assert err<1e-10
ratio=0.
for _ in range(100):
    x=rng.normal(size=d)*2; p=np.array([phi(a,x) for a in aa])
    J=np.einsum('a,aoi->oi',p,Cs)
    ratio=max(ratio,np.linalg.norm(J,2)**2/(float(p@p)*np.linalg.norm(Gram,2)))
assert ratio<=1+1e-12 and np.linalg.norm(Gram,2)<=EJop2+1e-10
out['matrix_parseval']={'max_entry_error':err,'max_cauchy_ratio':ratio,'gram_operator_norm':float(np.linalg.norm(Gram,2)),'quadrature_E_operator_norm_squared':float(EJop2)}

# Conditional Hermite projection, exact at degrees 0,...,9.
a,h,z=sp.symbols('a h z',real=True)
def integrate_standard(expr, variable):
    poly=sp.Poly(sp.expand(expr),variable)
    ans=0
    for (n,),coef in poly.terms():
        if n%2==0: ans+=coef*sp.factorial(n)/(2**(n//2)*sp.factorial(n//2))
    return sp.expand(ans)
for n in range(10):
    lhs=integrate_standard(sp.hermite_prob(n,a*h+sp.sqrt(1-a*a)*z),z)
    assert sp.simplify(lhs-a**n*sp.hermite_prob(n,h))==0
out['conditional_hermite_projection']={'exact_degrees_checked':list(range(10))}

# A conditional null-bank Riesz identity at the nonlinear interpolated endpoint.
t,k=sp.symbols('t k',real=True)
F=h*z+z*z-1
R=h+z
DF=h+2*z
X=h+t*F+k
left=F*4*X**3; right=t*R*DF*12*X**2
for variable in (h,z,k):
    left=integrate_standard(left,variable); right=integrate_standard(right,variable)
assert sp.simplify(left-right)==0
out['null_bank_riesz']={'polynomial_test':'F=HZ+He_2(Z), phi(x)=x^4','exact_identity':True}

# Exact alpha/share arithmetic.
assert Fraction(17,2)+8*Fraction(1,16)==9
assert 1-2*9==-17
assert 1-9==-8
assert Fraction(17,1)-Fraction(1,2)-10==Fraction(13,2)
assert Fraction(17,2)+Fraction(5,2)==11
out['share_exponents']={'rank':9,'root_share':-9,'native_scaled_root_square_share':-17,'local_drift_and_complete_first_share':-8,'equal_share_native_J':str(Fraction(17,2)),'equal_share_null_current_J':8,'shared_hybrid_alpha_after_taming':str(Fraction(33,2)),'new_guard_alpha':str(Fraction(13,2))}
assert Fraction(7,17)<Fraction(1,2)
assert (17-10)/Fraction(17,1)==Fraction(7,17)
assert (18-10)/Fraction(16,1)==Fraction(1,2)
out['positive_cubature_integration']={'new_J_exponent':'2 q_eta','native_grade':'17 - 17 q_eta - (17/2) nu_old','null_grade':'18 - 16 q_eta - 8 nu_old','sufficient_strict_q_with_public_log_old_census':'q_eta < 7/17'}

# Enumerate the 81 choices of force hits in each of the three labelled trees.
costs=[]
for center in range(3):
    leaves=[i for i in range(3) if i!=center]
    for hits in itertools.product(range(3),repeat=4):
        orders=[[0,1,0] for _ in range(3)]
        for edge,leaf in enumerate(leaves):
            orders[center][hits[2*edge]]+=1
            orders[leaf][hits[2*edge+1]]+=1
        assert sum(map(sum,orders))==7
        costs.append(sum(2**(q+1) for copy in orders for q in copy))
assert len(costs)==243 and max(costs)==44
out['force_hit_census']={'histories':len(costs),'max_raw_VALUE_calls':max(costs),'min_raw_VALUE_calls':min(costs)}

# Literal integer-radius construction and new finite-D guards; B=5,v=1/4.
def logadd(x,y):
    mm=max(x,y)
    return mm+math.log(math.exp(x-mm)+math.exp(y-mm))
def guard(d,m,power,B=5.,v=.25):
    logalpha=-power*math.log(2); logK=-logalpha/2
    ell=1
    while True:
        R=math.sqrt(d)+ell; logbase=math.log(1+d+R*R)
        logtail=logadd(-5*ell*ell/24,math.log(2)+logK+m/2*logbase-ell*ell/4)
        if logtail<=2.5*logalpha: break
        ell+=1
    logL=-.5*math.log(v)
    if m: logL=max(logL,math.log(2)+logK+(m-1)/2*logbase-.5*math.log(v))
    logLS=logL+math.log(B)+8.5*logalpha
    loghybrid=math.log(2)+logL+2*math.log(B)+7*logalpha
    return {'D':d,'m':m,'alpha_dyadic_power':power,'alpha':math.exp(logalpha),'ell':ell,'LS':math.exp(logLS),'normalized_hybrid_guard':math.exp(loghybrid),'normalized_tail_ratio':math.exp(logtail-2.5*logalpha)},logLS<=0 and loghybrid<=0
examples=[]
for d in (1,100,1000000):
    for m in (0,1,8,16):
        for power in range(1,201):
            row,passes=guard(d,m,power)
            if passes:
                assert row['normalized_tail_ratio']<=1+1e-12
                examples.append(row); break
        else: raise AssertionError('No sufficiently small alpha found in the check window')
out['illustrative_finite_D_guards']={'B':5,'v':.25,'c_hybrid':1,'rows':examples,'scope':'Only the new guards; native radius and other ledger rows remain separate.'}
path=Path(__file__).with_name('independent-check-results.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {cases} kernel checks; exact Hermite/Riesz/share identities; 243 histories; {len(examples)} finite-D guard examples.')
print(path)
