#!/usr/bin/env python3
"""Independent algebraic checks. Not an execution of imported native compilers."""
import hashlib, itertools, json, math
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(20261005)
counts={}

def check(name, value):
    if not bool(value):
        raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

def sym3(T):
    return sum(T.transpose(p) for p in itertools.permutations(range(3)))/6

def sqrtmat(C):
    w,V=np.linalg.eigh(C)
    return (V*np.sqrt(w))@V.T

def invsqrt(C):
    w,V=np.linalg.eigh(C)
    return (V*(1/np.sqrt(w)))@V.T

def cutnorm(T):
    return max(np.linalg.norm(np.moveaxis(T,i,0).reshape(T.shape[i],-1),2) for i in range(3))

# Generic covariance-block Hilbert inequality, with a dimension-free operator factor.
for d in (1,2,3,5,8):
    for _ in range(80):
        n=d*d
        Z=rng.normal(size=(d+n,d+n+3))
        Cov=Z@Z.T
        P=Cov[:d,:d]; Q=Cov[d:,d:]; R=Cov[:d,d:]
        rhs=np.linalg.norm(Q,2)*np.trace(P)
        check('cross_covariance_hilbert_bound', np.linalg.norm(R,'fro')**2<=rhs*(1+1e-12))
        check('cross_covariance_psd_order', np.linalg.eigvalsh(np.linalg.norm(Q,2)*P-R@R.T).min()>-1e-8)

# Exact finite-distribution centered tensor telescoping and the resulting three-current bound.
for d in (1,2,3,5):
    for _ in range(60):
        n=25
        u=rng.normal(size=(n,d)); u-=u.mean(axis=0)
        h=u+.2*rng.normal(size=(n,d)); h-=h.mean(axis=0)
        delta=u-h
        k_u=np.einsum('ni,nj,nk->ijk',u,u,u)/n
        k_h=np.einsum('ni,nj,nk->ijk',h,h,h)/n
        terms=[np.einsum('ni,nj,nk->ijk',delta,u,u)/n,
               np.einsum('ni,nj,nk->ijk',h,delta,u)/n,
               np.einsum('ni,nj,nk->ijk',h,h,delta)/n]
        check('cumulant_telescoping', np.allclose(k_u-k_h,sum(terms),atol=1e-12))
        bounds=[]
        for a,b in ((u,u),(h,u),(h,h)):
            Q=np.einsum('ni,nj->nij',a,b).reshape(n,-1)
            Q-=Q.mean(axis=0)
            covariance=Q.T@Q/n
            bounds.append(np.sqrt(np.linalg.norm(covariance,2)*np.mean(np.sum(delta**2,axis=1))))
        check('cumulant_one_energy_covariance_bound', np.linalg.norm(k_u-k_h)<=sum(bounds)*(1+1e-12))

# Exact Wick covariance operator for cross-quadratics of jointly Gaussian linear sources.
for d in (1,2,3,5):
    for _ in range(60):
        V=rng.normal(size=(d,2*d)); V/=max(np.linalg.norm(V,2),1)
        U=rng.normal(size=(d,2*d)); U/=max(np.linalg.norm(U,2),1)
        M=rng.normal(size=(d,d))
        Q=V.T@M@U; sym=(Q+Q.T)/2
        variance=2*np.linalg.norm(sym,'fro')**2
        check('quadratic_covariance_operator_bound', variance<=4*np.linalg.norm(M,'fro')**2*(1+1e-12))

# Exact anisotropic regression for multiple independent blocks with noncommuting covariances.
for d in (2,3,5):
    for _ in range(80):
        covs=[]
        for j in range(4):
            A=rng.normal(size=(d,d))
            covs.append(A@A.T+.2*np.eye(d))
        C=sum(covs); Ci=np.linalg.inv(C); s=rng.normal(size=d)
        common=np.outer(Ci@s,Ci@s)-Ci
        summed=np.zeros(d); target=np.zeros(d)
        for Cj in covs[:3]:
            K=rng.normal(size=(d,d,d)); K=(K+K.swapaxes(1,2))/2
            Cji=np.linalg.inv(Cj)
            mu=Cj@Ci@s
            V=Cj-Cj@Ci@Cj
            conditional=Cji@(np.outer(mu,mu)+V)@Cji-Cji
            check('anisotropic_regression_block', np.allclose(conditional,common,rtol=1e-9,atol=1e-10))
            summed+=np.einsum('iab,ab->i',K,conditional)
            target+=np.einsum('iab,ab->i',K,common)
        check('anisotropic_regression_sum', np.allclose(summed,target,rtol=1e-9,atol=1e-10))

# Physical cubic coefficient scaling, variance allocation, native prior, and exponent ledger.
for A in np.geomspace(1e-7,.02,80):
    q=1-A
    for share in (1/8,1/16,1/32):
        u=share*A; alpha=q*A/np.sqrt(u)
        for N in (1,7,31):
            cb=ab=bb=eb=1/(2*np.sqrt(N)); weight=.00023; sigma=.073
            r1=-4*weight*alpha/(cb*ab*bb*sigma)
            physical=(np.sqrt(u)**3*cb*ab*bb)*(r1*alpha**2)*sigma/A**3
            check('unchanged_cubic_target',np.isclose(physical,-4*q**3*weight,rtol=2e-13))
            check('node_variance_share',np.isclose(N*u*(cb**2+ab**2+bb**2+eb**2),u))
            check('physical_residual_first',np.isclose(np.sqrt(u)*alpha,q*A))
            check('native_feedback_scaling',np.isclose(np.sqrt(u)*alpha**6,q**6*A**6/u**2.5))
            check('order5_prior_scaling',np.isclose(np.sqrt(u)*alpha**5,q**5*A**5/u**2))
    check('half_buffer_allocation',sum([A/8]*4+[A/2])==A)

terms={'prefix_1':(2,Fraction(3,2)), 'prefix_2':(3,Fraction(1)),
       'near':(3,Fraction(1,2)), 'mean':(4,Fraction(0)),
       'skew_mean_gram':(5,Fraction(-3,2)), 'tensor_mismatch':(5,Fraction(-1)),
       'mixed_K_prior':(6,Fraction(-7,6)), 'cubic_order5_prior':(6,Fraction(-2)),
       'cubic_feedback':(7,Fraction(-5,2)), 'mixed_K_mixture':(7,Fraction(-3,2))}
exponents={key:Fraction(a)+b for key,(a,b) in terms.items()}
check('minimum_grade_7_over_2', min(exponents.values())==Fraction(7,2))
check('cubic_feedback_grade_9_over_2',exponents['cubic_feedback']==Fraction(9,2))

# Exact coefficient-Hilbert version of quadratic reference sensitivity.
for d in (2,3,5):
    for _ in range(80):
        v=10**rng.uniform(-4,0)
        R=rng.normal(size=(d,d)); C1=v*(np.eye(d)+.1*R@R.T/d)
        R=rng.normal(size=(d,d)); C2=v*(np.eye(d)+.1*R@R.T/d)
        K=rng.normal(size=(d,d,d)); K=(K+K.swapaxes(1,2))/2
        k=cutnorm(K); i1=invsqrt(C1); i2=invsqrt(C2)
        Q1=np.einsum('iab,aj,bk->ijk',K,i1,i1)
        Q2=np.einsum('iab,aj,bk->ijk',K,i2,i2)
        bound=2*k*np.linalg.norm(C1-C2,'fro')/v**2
        check('quadratic_reference_covariance_sensitivity',np.sqrt(2)*np.linalg.norm(Q1-Q2)<=bound*(1+1e-10))

# Exact polynomial orientation-interpolation identity, including anisotropic covariance.
x0,x1,z0,z1=sp.symbols('x0 x1 z0 z1')
x=sp.Matrix([x0,x1]); z=sp.Matrix([z0,z1]); variables=(x0,x1,z0,z1)
y0,y1=sp.symbols('y0 y1'); y=(y0,y1)
C=sp.Matrix([[sp.Rational(3,2),sp.Rational(1,4)],[sp.Rational(1,4),sp.Rational(5,4)]])
Ci=C.inv(); covariance=sp.diag(C,sp.eye(2))
cache={}
def moment(indices):
    if not indices:return sp.Integer(1)
    if len(indices)%2:return sp.Integer(0)
    indices=tuple(sorted(indices))
    if indices in cache:return cache[indices]
    first=indices[0]; rest=indices[1:]; ans=0
    for k,j in enumerate(rest):
        ans+=covariance[first,j]*moment(rest[:k]+rest[k+1:])
    cache[indices]=ans
    return ans

def expectation(expr):
    poly=sp.Poly(sp.expand(expr),*variables)
    ans=0
    for powers,coef in poly.terms():
        indices=tuple(i for i,p in enumerate(powers) for _ in range(p))
        ans+=coef*moment(indices)
    return sp.factor(ans)

for fixture in range(4):
    raw=np.array([int(a) for a in rng.integers(-3,4,size=8)],dtype=object).reshape(2,2,2)
    N=np.empty((2,2,2),dtype=object)
    for i,a,b in itertools.product(range(2),repeat=3):N[i,a,b]=sp.Rational(raw[i,a,b]+raw[i,b,a],14)
    K=np.empty_like(N)
    for idx in itertools.product(range(2),repeat=3):
        K[idx]=sum(N[tuple(idx[p] for p in perm)] for perm in itertools.permutations(range(3)))/6
    D=N-K
    tau=sp.Rational(fixture+1,7); T=K+tau*D
    check('zero_full_symmetry',all(sum(D[tuple(idx[p] for p in perm)] for perm in itertools.permutations(range(3)))==0 for idx in itertools.product(range(2),repeat=3)))
    score=Ci*x
    def field(A):
        return sp.Matrix([sum(A[i,a,b]*(score[a]*score[b]-Ci[a,b]) for a,b in itertools.product(range(2),repeat=2)) for i in range(2)])
    b=field(T); bd=field(D); J=b.jacobian(x); W=x+b+sp.Rational(2,3)*z
    phi=y0**4+y0*y1**3+y0**2*y1**2+sp.Rational(2,3)*y1**3
    substitution={y0:W[0],y1:W[1]}
    lhs=sum(bd[i]*sp.diff(phi,y[i]).subs(substitution) for i in range(2))
    r3=0;r2=0;lead=0
    for i,a,bb in itertools.product(range(2),repeat=3):
        lead+=D[i,a,bb]*sp.diff(phi,y[i],y[a],y[bb])
        for p,q in itertools.product(range(2),repeat=2):
            coef=D[i,a,bb]*((int(p==a)+J[p,a])*(int(q==bb)+J[q,bb])-int(p==a)*int(q==bb))
            r3+=coef*sp.diff(phi,y[i],y[p],y[q]).subs(substitution)
        for p in range(2):
            r2+=D[i,a,bb]*sp.diff(b[p],x[a],x[bb])*sp.diff(phi,y[i],y[p]).subs(substitution)
    check('orientation_leading_current_vanishes',sp.expand(lead)==0)
    elhs=expectation(lhs); erhs=expectation(r3+r2)
    check('exact_orientation_positive_path_identity',elhs==erhs)
    # The identity is not the false assertion that the derivative itself vanishes.
    check('orientation_feedback_is_nonzero',elhs!=0)

sources=[HERE.parent/'POSITIVE-SKEW-SHRINKING-BUFFER-JOIN.md',HERE/'INDEPENDENT-BUFFERED-SKEW-JOIN-AUDIT.md',Path('/workspace/shared/law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md')]
result={'status':'PASS','assertions':sum(counts.values()),'categories':counts,'w_equals_A_exponents':{k:str(v) for k,v in exponents.items()},'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'scope':'Exact algebraic and polynomial diagnostics; imported finite native compiler guards remain hypotheses, not executed runtime checks.'}
(HERE/'buffered_skew_join_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
