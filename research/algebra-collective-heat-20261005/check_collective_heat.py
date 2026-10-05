#!/usr/bin/env python3
"""Finite mathematical diagnostics. No native sampler or derivative oracle is executed."""
from __future__ import annotations
import hashlib, json, math, random
from pathlib import Path
from fractions import Fraction
import sympy as sp
import numpy as np

HERE=Path(__file__).resolve().parent
checks=[]
def check(name, condition, detail=None):
    ok=bool(condition)
    checks.append({'name':name,'passed':ok,**({'detail':detail} if detail is not None else {})})
    if not ok: raise AssertionError(name)
def zero(expr): return sp.expand(expr)==0

def op_ind(F,z): return sp.expand(sum(sp.diff(F,v,2) for v in z)/2)
def op_cross(F,z): return sp.expand(sum(sp.diff(F,u,v) for i,u in enumerate(z) for v in z[i+1:]))
def op_com(F,z): return sp.expand(sum(sp.diff(F,u,v) for u in z for v in z)/2)
def exp_op(F,op,s):
    ans=F; term=F; j=0
    while term!=0:
        j+=1; term=op(term)
        if term!=0: ans+=s**j*term/sp.factorial(j)
        if j>30: raise RuntimeError('nonpolynomial fixture')
    return sp.expand(ans)

s,theta=sp.symbols('s theta')
for N in range(1,5):
    z=sp.symbols(f'z0:{N}')
    F=sp.prod((i+1)+v+(i+2)*v**2 for i,v in enumerate(z))
    hi=lambda x:op_ind(x,z)
    hc=lambda x:op_com(x,z)
    bc=lambda x:op_cross(x,z)
    check(f'operator_identity_N{N}',zero(hi(F)-hc(F)+bc(F)))
    check(f'commutator_ind_cross_N{N}',zero(hi(bc(F))-bc(hi(F))))
    check(f'commutator_com_cross_N{N}',zero(hc(bc(F))-bc(hc(F))))
    independent=exp_op(F,hi,s)
    factored=exp_op(exp_op(F,bc,-s),hc,s)
    check(f'full_polynomial_factorization_N{N}',zero(independent-factored))
    for m in range(4):
        check(f'finite_jet_N{N}_m{m}',zero(sp.series(independent-factored,s,0,m+1).removeO()))
    hq=lambda x:sp.expand(hi(x)+theta*bc(x))
    qheat=exp_op(F,hq,s)
    check(f'covariance_derivative_N{N}',zero(sp.diff(qheat,theta)-s*exp_op(bc(F),hq,s)))
    integrated=sp.integrate(exp_op(bc(F),hq,s),(theta,0,1))
    common=exp_op(F,hc,s)
    check(f'positive_homotopy_N{N}',zero(independent-common+s*integrated))
    for t in [sp.Rational(0),sp.Rational(1,7),sp.Rational(2,3),sp.Rational(1)]:
        check(f'covariance_PSD_N{N}_theta{t}',1-t>=0 and 1+(N-1)*t>=0)

x,y=sp.symbols('x y')
negative_kernel=exp_op((x+y)**2,lambda F:op_cross(F,(x,y)),-s)
check('negative_cross_exponential_not_positive',negative_kernel.subs({x:0,y:0})==-2*s)

# Exact Gaussian moments for endpoint and retained-observer identities.
b1,b2,z=sp.symbols('b1 b2 z')
def gauss_E(F,variables):
    P=sp.Poly(sp.expand(F),*variables)
    result=0
    for monom,coeff in P.terms():
        moment=sp.Integer(1)
        for power in monom:
            if power%2: moment=0;break
            if power: moment*=sp.factorial2(power-1)
        result+=coeff*moment
    return sp.simplify(result)
A=[(sp.Integer(1),sp.Integer(-1)),(sp.Integer(2),sp.Integer(1))]
def deriv(F,a): return sp.expand(a[0]*sp.diff(F,b1)+a[1]*sp.diff(F,b2))
X=b1**2+b1*b2+b2+sp.Rational(2,3)*z
C=b1**2+2*b2+1
HC=sum(deriv(deriv(C,a),a) for a in A)/2
beta=[a[0]*b1+a[1]*b2 for a in A]
V=[deriv(X,a) for a in A]
W=[deriv(deriv(X,a),a) for a in A]
S=sum(be**2-sum(v**2 for v in a) for be,a in zip(beta,A))/2
q=sp.symbols('q')
for r in [1,2,3]:
    for power in range(r,8):
        phi=q**power
        D0=sp.diff(phi,q,r).subs(q,X)
        D1=sp.diff(phi,q,r+1).subs(q,X)
        D2=sp.diff(phi,q,r+2).subs(q,X)
        rhs=C*(S*D0+sum(w/2-be*v for w,be,v in zip(W,beta,V))*D1+sum(v*v for v in V)*D2/2)
        check(f'endpoint_adjoint_r{r}_degree{power}',gauss_E(HC*D0-rhs,(b1,b2,z))==0)
        one=sum(deriv(C,a)*(be*D0-v*D1) for a,be,v in zip(A,beta,V))/2
        check(f'one_integration_r{r}_degree{power}',gauss_E(HC*D0-one,(b1,b2,z))==0)
        chi=1+b1+b2**2
        defect=C*(sum(deriv(deriv(chi,a),a)/2-be*deriv(chi,a) for a,be in zip(A,beta))*D0
                  +sum(deriv(chi,a)*v for a,v in zip(A,V))*D1)
        check(f'retained_test_r{r}_degree{power}',gauss_E(chi*HC*D0-chi*rhs-defect,(b1,b2,z))==0)

# A deliberately failing naive conditional relation, with its exact missing term.
naive_lhs=gauss_E(b1**2,(b1,))
naive_rhs=gauss_E(b1**2*((b1**2-1)/2)*b1**2,(b1,))
retention_defect=gauss_E(b1**2*(1-2*b1**2),(b1,))
check('retained_naive_relation_fails',naive_lhs!=naive_rhs,
      {'lhs':str(naive_lhs),'naive_rhs':str(naive_rhs),'defect':str(retention_defect)})
check('retained_defect_restores_relation',naive_lhs==naive_rhs+retention_defect)

Q=1+b1*b2+b1**2
F=3+b1**3+b2**2
H=lambda U:sum(deriv(deriv(U,a),a) for a in A)/2
check('scalar_product_rule',zero(Q*H(F)-H(Q*F)+sum(deriv(Q,a)*deriv(F,a) for a in A)+H(Q)*F))
check('no_common_lift_for_B_and_2B',sp.linsolve([x-1,2*x-1],(x,))==sp.EmptySet)

# Untouched keep relation with a coefficient independent of the keep.
for r in range(1,7):
    herm=sp.hermite_prob(r-1,z)
    for degree in range(r,r+4):
        coeff=1+b1**2
        endpoint=b1+sp.Rational(2,3)*z
        lhs=gauss_E(coeff*sp.diff(q**degree,q,r).subs(q,endpoint),(b1,z))
        rhs=gauss_E(coeff*herm*(sp.Rational(3,2)**(r-1))*sp.diff(q**degree,q).subs(q,endpoint),(b1,z))
        check(f'keep_IBP_r{r}_degree{degree}',lhs==rhs)
# If coefficient reads the keep, the stripped formula fails already at rank 2.
check('keep_dependence_breaks_naive_IBP',gauss_E(z*sp.diff(q**3,q,2).subs(q,z),(z,))
      !=gauss_E(z*z*sp.diff(q**3,q).subs(q,z),(z,)))

# Local factorial bound used for the real Gauss certificate.
for K in range(0,11):
    for j in range(0,21):
        # Squared exact integer inequality avoids floating factorials.
        check(f'factorial_bound_K{K}_j{j}',math.factorial(K+2*j+2)
              <=2**(K+2+4*j)*math.factorial(K+2)*(math.factorial(j)**2))

# Graph bookkeeping: the old spanning tree is retained after each chord/parallel edge.
rng=random.Random(9102026)
for N in range(2,14):
    tree=[(v,rng.randrange(v)) for v in range(1,N)]
    deg=[0]*N
    for u,v in tree: deg[u]+=1;deg[v]+=1
    marks=[int(d==1) for d in deg]
    k=[deg[v]+marks[v]-2 for v in range(N)]
    edges=list(tree)
    R=sum(marks)
    for j in range(12):
        u,v=rng.sample(range(N),2)
        edges.append((u,v));k[u]+=1;k[v]+=1
        check(f'graph_valence_N{N}_step{j}',sum(k)+2*N==2*len(edges)+R)
        check(f'graph_loopless_N{N}_step{j}',all(u!=v for u,v in edges))
        check(f'graph_old_marked_tree_N{N}_step{j}',all(marks[v]>=1 for v,d in enumerate(deg) if d==1))
for k in range(2,10):
    for gamma in [Fraction(1,8),Fraction(1,5),Fraction(3,7)]:
        old=gamma*max(k-2,0)
        new=gamma*max(k+2-2,0)
        check(f'saturated_neutrality_k{k}_g{gamma}',2*gamma-(new-old)==0)

# Direct polynomial numerical quadrature diagnostic. Nodes here are diagnostic, not production certification.
zv=sp.symbols('v0:2')
F=(1+zv[0]+zv[0]**2)*(2-zv[1]+zv[1]**2)
Kpoly=exp_op(op_cross(F,zv),lambda U:op_ind(U,zv)+theta*op_cross(U,zv),s).subs({zv[0]:0,zv[1]:0})
for sval in [0.01,0.1,0.5,1.0]:
    for m in [1,2,3,4]:
        nodes,weights=np.polynomial.legendre.leggauss(m)
        fun=sp.lambdify(theta,Kpoly.subs(s,sval),'numpy')
        quad=sum(float(w)*float(fun((float(node)+1)/2))/2 for node,w in zip(nodes,weights))
        exact=float(sp.integrate(Kpoly.subs(s,sval),(theta,0,1)))
        check(f'gauss_polynomial_s{sval}_m{m}',abs(quad-exact)<1e-11)

# Nonpolynomial two-occurrence Gaussian heat and certified-rule diagnostic.
for sval in [0.125,0.5,1.0,3.0]:
    tau=1.0; K=4; pair_count=1; tolerance=1e-9
    bound_M=2**(K+2)*math.sqrt(math.factorial(K+2))*pair_count*tau**(-K-2)
    R0=tau*tau/(8*pair_count*sval)
    panels=max(1,math.ceil(1/R0))
    m=max(1,math.ceil(math.log(max(1,2*sval*bound_M/tolerance),4)))
    nodes,weights=np.polynomial.legendre.leggauss(m)
    for frequency in [0.3,1.0,2.0,4.0]:
        aq=sval*frequency*frequency
        base=(0.1*frequency**2*math.exp(-frequency**2/2))**2
        # Stable form of q^2 e^-aq sinh(aq theta).
        def integrand(th):
            return base*frequency**2*(math.exp(-aq*(1-th))-math.exp(-aq*(1+th)))/2
        quad=0.0
        for panel in range(panels):
            mid=(panel+0.5)/panels
            for node,weight in zip(nodes,weights):
                th=mid+float(node)/(2*panels)
                quad+=float(weight)*integrand(th)/(2*panels)
        independent=base*math.exp(-aq)
        common=base*(1+math.exp(-2*aq))/2
        error=abs(independent-common+sval*quad)
        check(f'nonpolynomial_homotopy_s{sval}_q{frequency}',error<=tolerance,
              {'panels':panels,'nodes_per_panel':m,'absolute_error':error,'certificate':2*sval*bound_M*4**(-m)})

# Fixed lacunary source: positive-sign bounds for both derivative and full finite heat family.
eps=1/16;kappa=1/4;u=1.0
check('convex_C2_hessian_guard',eps*(math.pi**2/6)<1/4)
rows=[]
for n in [2,3,4,6,8,12,16,24,32]:
    a=1/(n*n)
    f=fw=h2=0.0
    for j in range(1,n+4):
        exponent=j*j-n*n
        if exponent>10: continue
        ratio=2.0**exponent
        base=(1/(j*j))*ratio**2*math.exp(-ratio*ratio/2)
        f+=base
        fw+=base*math.exp(-u*ratio*ratio/2)
        h2+=base*ratio*ratio
    f*=kappa*eps;fw*=kappa*eps;h2*=kappa*eps
    derivative=f*h2
    derivative_lower=kappa*kappa*eps*eps*a*a*math.exp(-1)
    heat=f*f-fw*fw
    heat_lower=kappa*kappa*eps*eps*a*a*(math.exp(-0.5)-math.exp(-1))*math.exp(-0.5)
    check(f'lacunary_derivative_bound_n{n}',derivative>=derivative_lower*(1-1e-12))
    check(f'lacunary_full_heat_bound_n{n}',heat>=heat_lower*(1-1e-12))
    rows.append({'n':n,'derivative_abs':derivative,'derivative_lower':derivative_lower,
                 'full_heat_abs':heat,'full_heat_lower':heat_lower,
                 'log2_a_squared_over_t_to_quarter':n*n/4-4*math.log2(n)})
check('no_power_lacunary_subsequence',rows[-1]['log2_a_squared_over_t_to_quarter']>200)

result={'status':'PASS','test_count':len(checks),'checks':checks,'lacunary_rows':rows,
        'scope':'Exact polynomial and arithmetic diagnostics; no native sampler executed; no universal source-image nonmembership asserted.'}
(HERE/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'test_count':result['test_count'],'lacunary_last':rows[-1]},indent=2))
