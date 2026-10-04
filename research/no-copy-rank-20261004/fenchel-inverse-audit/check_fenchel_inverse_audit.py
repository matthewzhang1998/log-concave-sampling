#!/usr/bin/env python3
"""Numerical verification of explicit formulas; proofs are in the adjacent report."""
import json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

checks=[]
def check(label, condition, **data):
    assert bool(condition), (label,data)
    checks.append({'label':label, **data})

m,L=.4,.6
lower=(31-math.sqrt(1261))/30
upper=(31+math.sqrt(761))/20
Hlo=np.array([[.4,-1.],[-1.,1/.6]])
Hhi=np.array([[.6,-1.],[-1.,1/.4]])
check('scalar off-graph determinant', abs(np.linalg.det(Hlo)+1/3)<1e-14, determinant=float(np.linalg.det(Hlo)))
check('sharp block endpoints', abs(np.linalg.eigvalsh(Hlo)[0]-lower)<1e-14 and abs(np.linalg.eigvalsh(Hhi)[-1]-upper)<1e-14, lower=lower, upper=upper)

rng=np.random.default_rng(20261004)
for d in [1,2,3,8]:
    for k in range(40):
        U,_=np.linalg.qr(rng.normal(size=(d,d)))
        W,_=np.linalg.qr(rng.normal(size=(d,d)))
        A=U@np.diag(rng.uniform(m,L,d))@U.T
        B=W@np.diag(rng.uniform(m,L,d))@W.T
        Q=np.block([[A,-np.eye(d)],[-np.eye(d),np.linalg.inv(B)]])
        ev=np.linalg.eigvalsh(Q)
        check(f'random Hessian envelope d={d}, k={k}',ev[0]>=lower-1e-12 and ev[-1]<=upper+1e-12)

# Smooth explicit curl fixture.
eps=.05
v=np.ones(2)
e=np.array([1.,0.])
def g(z):
    return .5*z+eps*np.sin(z[0])*e+(eps/2)*np.sin(z.sum())*v
def A(z):
    return .5*np.eye(2)+eps*np.cos(z[0])*np.outer(e,e)+(eps/2)*np.cos(z.sum())*np.outer(v,v)
p=np.full(2,math.pi/4)
u=np.zeros(2); J=np.zeros((2,2)); us=[u.copy()]
for n in range(3):
    J=(np.eye(2)-2*A(u))@J+2*np.eye(2)
    u=u+2*(p-g(u)); us.append(u.copy())
expected=-4*eps**2*math.sin(2*eps)*np.array([[0.,1.],[-1.,0.]])
check('explicit u1',np.max(np.abs(us[1]-[math.pi/2,math.pi/2]))<1e-14)
check('explicit u2',np.max(np.abs(us[2]-[math.pi/2-2*eps,math.pi/2]))<1e-14)
check('exact finite inverse skew',np.max(np.abs(J-J.T-expected))<1e-14 and abs(expected[0,1])>1e-4, observed=(J-J.T).tolist(), expected=expected.tolist())
for k in range(100):
    z=rng.normal(size=2)*20
    ev=np.linalg.eigvalsh(A(z))
    check(f'smooth fixture Hessian bounds {k}',ev[0]>=m-1e-12 and ev[-1]<=L+1e-12)

# Differentiate the complete inverse recursion by first sites only.
for k in range(50):
    p0=rng.normal(size=2)
    u=np.zeros(2); J=np.zeros((2,2))
    for n in range(1,21):
        J=(np.eye(2)-2*A(u))@J+2*np.eye(2)
        u=u+2*(p0-g(u))
        check(f'finite first bound seed={k}, N={n}',np.linalg.norm(J,2)<=2.5*(1-.2**n)+1e-12)

# The localized FIRST counterexample is evaluated with exact formula for the
# orbit distances, avoiding subtraction loss when delta becomes microscopic.
for n in range(2,101):
    delta=.1**(n+1)/2
    nearest_query_distance=.1**(n-1)
    derivative=(1-(-.1)**n)/.55
    check(f'localized FIRST error N={n}',nearest_query_distance>delta and abs(derivative-2)>=.18-1e-14, derivative_error=abs(derivative-2))

# Scalar normalized fibers, using y=h(p) and dp=g'(y)dy.
def gs(x): return .5*x+.1*math.sin(x)
def As(x): return .5+.1*math.cos(x)
def V(x): return .25*x*x+.1*(1-math.cos(x))
def fiber(x,tau):
    width=math.sqrt(tau)
    # Integrate in u=(y-x)/sqrt(tau); tails at |u|>24 negligible.
    def base(u):
        y=x+width*u
        gap=V(x)-V(y)-gs(y)*(x-y)
        return math.exp(-gap/tau)*As(y)*width
    z=quad(base,-24,24,epsabs=1e-12,epsrel=2e-11)[0]
    moment=quad(lambda u:gs(x+width*u)*base(u),-24,24,epsabs=1e-12,epsrel=2e-11)[0]
    return z,moment/z
fiber_results=[]
for tau in [.1,.01,.001,.0001]:
    z0,_=fiber(0,tau); zp,_=fiber(math.pi,tau)
    check(f'fiber normalization bounds tau={tau}',math.sqrt(2*math.pi*m*tau)<=z0<=math.sqrt(2*math.pi*L*tau) and math.sqrt(2*math.pi*m*tau)<=zp<=math.sqrt(2*math.pi*L*tau))
    fiber_results.append({'tau':tau,'ratio_Z0_Zpi':z0/zp,'limit':math.sqrt(1.5)})
check('Laplace determinant bias',abs(fiber_results[-1]['ratio_Z0_Zpi']-math.sqrt(1.5))<1e-4, ratios=fiber_results)
x=.8; tau=.07; h=1e-4
z,b=fiber(x,tau)
fd=(math.log(fiber(x+h,tau)[0])-math.log(fiber(x-h,tau)[0]))/(2*h)
score=(b-gs(x))/tau
check('normalizer score',abs(fd-score)<1e-8, finite_difference=fd, exact_score=score)

# VALUE contraction check on the scalar nonlinear fixture.
for p0 in np.linspace(-5,5,31):
    exact=brentq(lambda z:gs(z)-p0,-20,20)
    u=0.
    for n in range(1,16):
        u=u+2*(p0-gs(u))
        check(f'inverse VALUE error p={p0}, N={n}',abs(u-exact)<=2.5*(.2**n)*abs(p0)+3e-12)

# Gaussian ambient prior: exact inverse precision matrix.
for alpha in [.4,.5,.6]:
    for tau in [.001,.01,.1,1.]:
        precision=np.eye(2)+np.array([[alpha,-1],[-1,1/alpha]])/tau
        variance=np.linalg.inv(precision)[0,0]
        expected_var=1/(1+alpha**2/(1+alpha*tau))
        check(f'Gaussian root bias alpha={alpha}, tau={tau}',abs(variance-expected_var)<1e-12, variance=variance)

# At epsilon=0 the independent positive-temperature fibers introduce energy.
r=.1; a=.2; tau=.03
cov=np.eye(2)*(.5*tau); row=np.array([r*a/2,-r*a/2])
variance=float(row@cov@row)
check('false K3 energy at epsilon=0',abs(variance-r*r*a*a*tau/4)<1e-18 and variance>0, native_energy=0, penalty_energy=variance)
check('omitting graph tangent current fails',.5!=0, omitted_left_side=0, true_right_side=.5)

out={'status':'passed','count':len(checks),'constants':{'hessian_lower':lower,'first_upper':upper,'inverse_contraction':.2},'checks':checks}
path=Path(__file__).with_name('fenchel_inverse_audit_checks.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'passed','count':len(checks),'results':path.name,'fiber_ratios':fiber_results},indent=2))
