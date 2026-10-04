#!/usr/bin/env python3
"""Independent deterministic diagnostics, not substitutes for the analytic proof."""
from pathlib import Path
import hashlib, itertools, json, math
import numpy as np
from numpy.polynomial.hermite_e import hermeval, hermeroots
from scipy.integrate import quad
from scipy.special import eval_hermitenorm, ndtr

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'SCALAR-POSITIVE-HERMITE-MEAN-LAW-FACTORY.md'
SQRT_2PI = math.sqrt(2 * math.pi)
def phi(x): return math.exp(-x*x/2) / SQRT_2PI
def ints(K,R):
    return np.array([phi(R)*(eval_hermitenorm(k-1,-R)-eval_hermitenorm(k-1,R)) for k in range(1,K+1)])
def coeffs(xs):
    return np.array([0.] + [float(p)/math.factorial(k) for k,p in enumerate(np.cumprod(xs),1)])
def ratio(c,C,R,z):
    return 1-C + (hermeval(z,c) if abs(z)<=R else 0.)
def split_quad(f,R):
    return sum(quad(f,a,b,epsabs=1e-14,epsrel=1e-11,limit=200)[0] for a,b in [(-16.,-R),(-R,R),(R,16.)])

rows=[]
for K,R in itertools.product([1,2,4,6,8],[1.,2.,4.,6.]):
    eps=.97*(2*math.log(1.25)/(math.sqrt(R*R+2*math.log(1.25))+R))
    Is=ints(K,R)
    p=.7
    mu=(2*p-1)*eps
    cmu=coeffs([mu]*K)
    Cmu=float(cmu[1:]@Is)
    probes=np.r_[-R-1,np.linspace(-R,R,101),R+1]
    average=np.zeros_like(probes)
    minratio=float('inf'); maxratio=-float('inf')
    normerror=0.; accept_error=0.
    for signs in itertools.product([-1,1],repeat=K):
        c=coeffs(eps*np.array(signs)); C=float(c[1:]@Is)
        # For fixed z, the batch ratio is multi-affine in every X_i; checking
        # vertices suffices for its extrema over the complete batch box.
        deriv=np.array([k*c[k] for k in range(1,K+1)])
        roots=hermeroots(deriv) if K>1 else []
        points=[-R,R]
        points += [float(np.real(x)) for x in roots if abs(np.imag(x))<1e-8 and -R<=np.real(x)<=R]
        vals=[1-C]+[ratio(c,C,R,z) for z in points]
        minratio=min(minratio,min(vals)); maxratio=max(maxratio,max(vals))
        # Numerically integrate the interior and use exact Gaussian exterior mass.
        norm=quad(lambda z:phi(z)*ratio(c,C,R,z),-R,R,epsabs=1e-13)[0]+(1-C)*2*ndtr(-R)
        normerror=max(normerror,abs(norm-1))
        accept_error=max(accept_error,abs((2/3)*norm-2/3))
        nplus=sum(s==1 for s in signs)
        wt=p**nplus*(1-p)**(K-nplus)
        average += wt*np.array([ratio(c,C,R,z) for z in probes])
    exact=np.array([ratio(cmu,Cmu,R,z) for z in probes])
    expectation_error=float(np.max(abs(average-exact)))
    # Stable evaluation of 1 + S_cut - C - exp(mu*z - mu^2/2).
    def difference_ratio(z):
        return (hermeval(z,cmu) if abs(z)<=R else 0.)-Cmu-math.expm1(mu*z-mu*mu/2)
    l2sq=split_quad(lambda z: phi(z)*difference_ratio(z)**2,R)
    l2=math.sqrt(max(0,l2sq))
    second_weighted_l1=split_quad(lambda z:phi(z)*z*z*abs(difference_ratio(z)),R)
    tail_l2=math.sqrt(max(0,quad(lambda z:phi(z)*hermeval(z,cmu)**2,-16.,-R,epsabs=1e-14)[0]+quad(lambda z:phi(z)*hermeval(z,cmu)**2,R,16.,epsabs=1e-14)[0]))
    remsq=sum(mu**(2*k)/math.factorial(k) for k in range(K+1,70))
    rem=math.sqrt(remsq)
    decomposition_bound=rem+tail_l2+abs(Cmu)
    second_moment=1+quad(lambda z:phi(z)*z*z*hermeval(z,cmu),-R,R,epsabs=1e-14)[0]-Cmu
    assert minratio>=.5-1e-12 and maxratio<=1.5+1e-12
    assert normerror<2e-12 and expectation_error<2e-12
    assert l2<=decomposition_bound+2e-12
    assert second_weighted_l1<=math.sqrt(3)*l2+2e-12
    assert second_moment<=1.5+1e-12
    rows.append(dict(K=K,R=R,epsilon=eps,guard=eps*R+eps*eps/2,min_ratio=minratio,max_ratio=maxratio,max_normalization_error=normerror,max_acceptance_error=accept_error,max_averaged_profile_error=expectation_error,L2_error=l2,L2_decomposition_bound=decomposition_bound,weighted_second_L1=second_weighted_l1,second_moment=second_moment))

# Conditioning on the first symmetric batch record leaves a non-Gaussian law.
R=2.; eps=.1
m2=1-2*(R*phi(R)+ndtr(-R))
retention=dict(R=R,epsilon=eps,unretained_profile_error=0.,conditional_mean_given_X1_positive=eps*m2,conditional_W2_lower_bound=eps*m2)
# Sharing an ancestor destroys product averaging already at degree two.
ancestry=dict(epsilon=eps,EX=0.,expected_product_with_shared_sign=eps**2,required_product_of_means=0.)
result=dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),case_count=len(rows),notes='Deterministic floating-point diagnostics. Analytic proof, not this finite grid, establishes the theorem. Gaussian error quadratures truncate at absolute z=16.',max_normalization_error=max(r['max_normalization_error'] for r in rows),max_expectation_error=max(r['max_averaged_profile_error'] for r in rows),global_diagnostic_ratio_range=[min(r['min_ratio'] for r in rows),max(r['max_ratio'] for r in rows)],retention_counterexample=retention,shared_ancestor_counterexample=ancestry,cases=rows)
(HERE/'diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
