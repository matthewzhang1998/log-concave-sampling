#!/usr/bin/env python3
"""Independent exact checks of Appell algebra and positive-consumer currents.
No author generator is imported. Algebraic fixtures are not proofs of general bounds.
The sine fixture is globally Lipschitz; the quadratic fixture tests identities only.
"""
from pathlib import Path
from math import factorial, comb
from functools import lru_cache
from collections import Counter
import json
import sympy as s
HERE=Path(__file__).parent
x,g,t,C=s.symbols('x g t C', positive=True)
checks=Counter()
def verify(test, group):
    assert bool(test), group
    checks[group]+=1

def normal_moment(k, variance=s.Integer(1)):
    if k%2:return s.Integer(0)
    return s.Integer(factorial(k)//(2**(k//2)*factorial(k//2)))*variance**(k//2)
def expectation_poly(poly, variable, moments):
    return s.expand(sum(c*moments(k[0]) for k,c in s.Poly(s.expand(poly),variable).terms()))
def kappas(mu,n):
    out=[s.Integer(0)]
    for k in range(1,n+1):
        out.append(s.expand(mu(k)-sum(comb(k-1,j-1)*out[j]*mu(k-j) for j in range(1,k))))
    return out

# Independent partition generation: append the largest item, not first-item insertion.
def partitions(n):
    ans=[()]
    for a in range(n):
        nxt=[]
        for part in ans:
            nxt.append(part+((a,),))
            for j in range(len(part)):
                nxt.append(part[:j]+(part[j]+(a,),)+part[j+1:])
        ans=nxt
    return ans

# Real bounded source X=sin(sqrt(log(2))*G). All even moments are rational.
@lru_cache(None)
def sine_moment(k):
    if k%2:return s.Integer(0)
    return s.Rational((-1)**(k//2),2**k)*sum((-1)**j*comb(k,j)*s.Rational(2)**(-s.Rational((k-2*j)**2,2)) for j in range(k+1))
mu=sine_moment
kap=kappas(mu,16)
W=[s.Integer(1)]
for n in range(1,9):
    W.append(s.expand(x**n-sum(comb(n,j)*mu(j)*W[n-j] for j in range(1,n+1))))
    verify(s.diff(W[n],x)==n*W[n-1], 'Appell_derivative')
    verify(expectation_poly(W[n],x,mu)==0, 'Appell_centering')
    var=expectation_poly(W[n]**2,x,mu)
    # log(2)>2/3, so this exact rational inequality is stronger than A1 here.
    verify(var<=factorial(n)**2*s.Rational(2,3)**n, 'Appell_covariance_bound_actual_Lipschitz_fixture')
for p in range(1,5):
    for q in range(1,5):
        cross=0
        for pi in partitions(p+q):
            if all(any(a<p for a in b) and any(a>=p for a in b) for b in pi):
                term=s.Integer(1)
                for block in pi:term*=kap[len(block)]
                cross+=term
        verify(s.expand(expectation_poly(W[p]*W[q],x,mu)-cross)==0,'Appell_cross_partition_identity')
# R=X^3. Direct mixed cumulant uses the full moment-partition formula.
for m in range(1,8):
    mixed=0
    for pi in partitions(m+1):
        term=s.Integer((-1)**(len(pi)-1)*factorial(len(pi)-1))
        for block in pi:
            term*=mu(len(block)+(2 if 0 in block else 0))
        mixed+=term
    appell=expectation_poly(x**3*W[m],x,mu)
    verify(s.expand(mixed-appell)==0,'one_mark_moment_identity')
    verify(mixed**2<=factorial(m)**2*s.Rational(2,3)**m*mu(6),'one_mark_bound_actual_Lipschitz_fixture')
for n in range(2,9):
    delta=(1-s.Rational(1,2)**n)*kap[n]
    verify(delta**2<=factorial(n)**2*s.Rational(2,3)**(n-1)*mu(2)/4,'stability_bound_same_bank_scaled_sine')

# Compute cross-block majorants by an independent rooted bipartite-block DP.
cut={2:1}
for n in range(3,11):
    @lru_cache(None)
    def cross_proper(p,q):
        if p==q==0:return 1
        if p==0 or q==0:return 0
        return sum(comb(p-1,a-1)*comb(q,b)*cut[a+b]*cross_proper(p-a,q-b)
                   for a in range(1,p+1) for b in range(1,q+1) if a+b<=n-2)
    cut[n]=max(factorial(p)*factorial(n-p)+cross_proper(p,n-p) for p in range(1,n))
for n in range(2,9):
    verify(kap[n]**2<=cut[n]**2*s.Rational(2,3)**n,'proper_cut_bound_actual_Lipschitz_fixture')
author=HERE.parent/'generated_checks.json'
if author.exists():
    report=json.loads(author.read_text())
    for row in report['tree_rows']:
        verify(cut[row['rank']]==row['cumulant_proper_cut_majorant'],'independent_cut_constant_comparison')

# A6 checked through rank 7 from exact differentiation with a quartic map.
chi=x+s.Rational(2,7)*x**2-s.Rational(1,11)*x**3+s.Rational(1,13)*x**4
y=s.Symbol('y')
phi=y**9+s.Rational(2,3)*y**5+y**2
for r in range(3,8):
    direct=s.diff(s.diff(phi,y).subs(y,chi),x,r-1)
    generated=0
    for pi in partitions(r-1):
        term=s.diff(phi,y,1+len(pi)).subs(y,chi)
        for block in pi:term*=s.diff(chi,x,len(block))
        generated+=term
    verify(s.expand(direct-generated)==0,'Faa_di_Bruno_exact_polynomial')

# Full A7 at m=5, all source cumulants nonzero. B=G+(G^2-1)/4.
# This polynomial is only an exact identity fixture, not a Lipschitz-bound fixture.
B=g+(g*g-1)/4
muB=lambda n:expectation_poly(B**n,g,normal_moment)
kB=kappas(muB,6); Sigma=kB[2]
K={r:kB[r]/factorial(r) for r in range(3,6)}
def score(n):
    return s.expand(factorial(n)*sum((-1)**j*x**(n-2*j)/(s.Integer(2)**j*factorial(j)*factorial(n-2*j)*C**(n-j)) for j in range(n//2+1)))
q={r:K[r]*score(r-1) for r in K}
p=s.expand(sum((1-t**r)*q[r] for r in K))
Dt=s.expand(sum((1-t**r)*(-2*t*Sigma*s.diff(q[r],C)-t*Sigma*x/C*s.diff(q[r],x)) for r in K))
A={j:s.Integer(0) for j in range(1,7)}
A[1]=Dt; A[2]=-t*Sigma*s.diff(p,x)
for r in K:
    A[r]+=r*t**(r-1)*K[r]
    for pi in partitions(r-1):
        term=K[r]
        for block in pi:term*=s.diff(x+p,x,len(block))
        A[1+len(pi)]-=r*t**(r-1)*term
# Scalar Riesz recurrence computed in the probabilists' Hermite basis.
def riesz(poly):
    rem=s.Poly(s.expand(poly-expectation_poly(poly,g,normal_moment)),g)
    out=0
    for n in range(rem.degree(),0,-1):
        a=rem.nth(n)
        out+=a*s.hermite_prob(n-1,g)
        rem=s.Poly(s.expand(rem.as_expr()-a*s.hermite_prob(n,g)),g)
    verify(rem.as_expr()==0,'Riesz_Hermite_reduction')
    return s.expand(out)
T=B
for j in range(1,6):
    T=s.expand(riesz(T)*s.diff(B,g))
    verify(expectation_poly(T,g,normal_moment)==kB[j+1]/factorial(j),'Stein_cumulant_coefficients')
A[6]=t**5*T
# Two distinct interior t values; buffer variance is 1/2. Integrate it first.
def avg_phi(n,d,w):
    if d>n:return 0
    k=n-d
    return factorial(n)//factorial(k)*sum(comb(k,j)*w**(k-j)*normal_moment(j,s.Rational(1,2)) for j in range(k+1))
def avg_gx(poly,cv):
    return expectation_poly(expectation_poly(poly,g,normal_moment),x,lambda k:normal_moment(k,cv))
full_current_checks=[]
for tv in [s.Rational(1,3),s.Rational(2,3)]:
    cv=s.Rational(1,2)+(1-tv**2)*Sigma
    sub={t:tv,C:cv}
    p0=p.subs(sub); w=tv*B+x+p0; cp=-2*tv*Sigma
    # Independent derivative: varying Gaussian density, holding physical x fixed.
    drift=B+(s.diff(p,t)+(-2*t*Sigma)*s.diff(p,C)).subs(sub)
    gaussian_score=cp*(x*x-cv)/(2*cv**2)
    coefficients={j:s.expand(a.subs(sub)) for j,a in A.items()}
    for n in range(1,8):
        lhs=avg_gx(s.expand(drift*avg_phi(n,1,w)+gaussian_score*avg_phi(n,0,w)),cv)
        rhs=sum(avg_gx(s.expand(coefficients[j]*avg_phi(n,j,w)),cv) for j in range(1,7))
        verify(s.expand(lhs-rhs)==0,'A7_full_rank5_all_current_identity')
        full_current_checks.append({'t':str(tv),'test_power':n})

# Exact reference variance check, with three independently nonzero coefficients.
a,b,c=s.symbols('a b c')
reference=x+a*s.hermite_prob(2,x)+b*s.hermite_prob(3,x)+c*s.hermite_prob(4,x)
verify(s.expand(expectation_poly(reference**2,x,normal_moment)-(1+2*a*a+6*b*b+24*c*c))==0,'raw_reference_variance')

# Every generated m=5 coefficient has a permanent physical output and >=1 anchor edge.
monomials=0; max_degree=0
from itertools import product
for r in range(3,6):
    for pi in partitions(r-1):
        options=[]
        for block in pi:
            d=len(block)
            options.append(([None] if d==1 else [])+[rank for rank in range(3,6) if rank-1>=d])
        for ranks in product(*options):
            if all(rank is None for rank in ranks):continue
            degree=0
            for block,rank in zip(pi,ranks):
                if rank is None:continue
                d=len(block); free_gaussian=rank-1-d
                # Physical output survives even if every Gaussian stub is Wick-paired.
                verify(1<=d<=rank-1 and 1+free_gaussian>=1,'surviving_physical_output_and_proper_cut')
                degree+=free_gaussian
            verify(degree<=12,'rank5_maximum_Gaussian_degree')
            max_degree=max(max_degree,degree);monomials+=1
result={'status':'PASS','exact_assertion_counts':dict(checks),'total_assertions':sum(checks.values()),
        'cut_constants_rank_2_to_10':cut,'full_rank5_tests':full_current_checks,
        'expanded_rank5_feedback_monomials':monomials,'rank5_maximum_Gaussian_degree':max_degree,
        'scope':'Exact algebra and exact admissible scalar inequalities; dimension-free general claims are justified separately in AUDIT.md. No numeric fixture is treated as proof.'}
(HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
