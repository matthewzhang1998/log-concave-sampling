#!/usr/bin/env python3
"""Independent finite diagnostics, not a substitute for nonlinear proofs."""
import json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sp
out={"checks":[],"warnings":[]}
def check(name, cond, **detail):
    assert bool(cond), (name, detail)
    out["checks"].append({"name":name, **detail})
# Exact rational exponent and nested-bank fixed point identities.
for j in range(2,17):
    P=2*j-1
    b=F(2*j-1)+F(1,2)
    check(f"grades_j{j}",(b-F(1,2),b-F(1,2),b+1+F(1,2))==(P,P,P+2))
    for l in range(j+1):
        count=2*(j-1) if l==0 else 2*(j-l)
        cost=2*max(l-1,0)
        check(f"bank_cost_j{j}_l{l}",count+cost==2*(j-1))
    check(f"global_rate_j{j}", F(1)+F(2*P-1,P)==F(3)-F(1,P), P=P)
    check(f"dimension_rate_j{j}",F(2*P-1,2*P)==1-F(1,2*P))
# Exact scalar quadratic Picard: quadrature exact at relevant polynomial degree.
t,h,q,x,v,u=sp.symbols('t h q x v u')
for j in range(2,11):
    y=x+t*v
    for _ in range(j-1):
        y=sp.expand(x+t*v-q*sp.integrate((t-u)*y.subs(t,u),(u,0,t)))
    vend=sp.expand(v-q*sp.integrate(y,(t,0,h)))
    exact=sp.series(v*sp.cos(sp.sqrt(q)*h)-x*sp.sqrt(q)*sp.sin(sp.sqrt(q)*h),h,0,2*j+3).removeO()
    err=sp.expand(vend-exact)
    powers=[n for n in range(2*j+3) if sp.expand(err).coeff(h,n)!=0]
    check(f"quadratic_picard_order_j{j}",min(powers)==2*j+1,
          leading=str(err.coeff(h,2*j+1)))
# Scalar quadratic contraction matrix including both OU half steps.
minratio=1e9
for j in [2,3,5,8]:
  for kap in [1,2,10,100,1000]:
    M=np.array([[.5+.5/kap,.5],[.5,1.]])
    for lam in np.geomspace(1/(2*kap),1,20):
      for scale in [.001,.01,.03]:
        h0=scale/kap; rho=math.exp(-h0/2)
        s=j-1
        ca=sum((-lam)**k*h0**(2*k)/math.factorial(2*k) for k in range(s+2))
        sv=sum((-lam)**k*h0**(2*k+1)/math.factorial(2*k+1) for k in range(s+2))
        vc=sum((-lam)**(k+1)*h0**(2*k+1)/math.factorial(2*k+1) for k in range(s+1))
        vv=sum((-lam)**k*h0**(2*k)/math.factorial(2*k) for k in range(s+2))
        B=np.array([[ca,rho*sv],[rho*vc,rho*rho*vv]])
        loss=M-B.T@M@B
        ratio=np.linalg.eigvalsh(loss).min()/(h0/kap)
        minratio=min(minratio,ratio)
        assert ratio>0
check('quadratic_finite_reference_contraction', minratio>0,min_normalized_metric_loss=minratio,cases=4*5*20*3)
# Exact observation carrier covariance for several decoder lengths.
for k in [1,2,3,8,30]:
    C=np.eye(k+1)-np.ones((k+1,k+1))/(k+1)
    cov=np.ones((k+1,k+1))/(k+1)+C@C.T
    check(f'decoder_observation_covariance_k{k}',np.allclose(cov,np.eye(k+1),atol=1e-13))
# Finite mode recurrence envelope, including large translated quadratics.
for A in [.001,.01,.05,.125]:
  for lam in [.01,.5,1.]:
    for center in [0.,1.,1e6]:
      y=3.; mode=(y+A*lam*center)/(1+A*lam); bm=y
      for m in range(1,13):
        bm=y-A*lam*(bm-center)
        residual=abs(bm-mode); bound=A**(m+1)*abs(lam*(y-center))/(1-A)
        # Machine cancellation dominates once exact residual is tiny.
        assert residual <= bound+2e-9*max(1,abs(mode))
check('finite_mode_caller_envelope',True,cases=4*3*3*12)
# Counterexample to missing unanchored finite-Picard caller term.
A,g=sp.symbols('A g', positive=True)
c=sp.sqrt(2)/2
# With exact mode seed and zero momentum, predictor mode error at T/2=A(1-c)g.
# First Picard endpoint error is -A*w_21*predictor_error; w_21=1-c.
res=-A**2*(1-c)**2*g
check('unanchored_predictor_counterterm',sp.simplify(res/(A**2*g))!=0,exact_residual=str(sp.expand(res)))
out['warnings'].append('The literal unanchored finite predictor needs an A^(M+1)|grad V(mode)| numerical caller term, even with exact mode seed. This fixture detects a proof omission, not failure after adding that term.')
# Gaussian seed posterior law: exact centered scalar covariance W2 below A^(3/2).
for A in [.0001,.001,.01,.125]:
  for lam in [.01,.3,1.]:
    gap=math.sqrt(A)-math.sqrt(A/(1+A*lam))
    assert gap<=A**1.5
check('gaussian_seed_quadratic_law',True,cases=12)
out['passed']=len(out['checks']); out['status']='PASS_WITH_EXPLICIT_PREDICTOR_REPAIR'
p=Path(__file__).with_name('direct_mean_outer_checks.json');p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':out['passed'],'status':out['status'],'minimum_contraction_ratio':minratio,'output':p.name}))
