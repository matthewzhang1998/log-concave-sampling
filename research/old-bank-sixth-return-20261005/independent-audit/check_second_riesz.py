"""Exact polynomial diagnostics, independent of the constructor's sin fixture.

The polynomial fixture is algebraic, not a bounded-Hessian original-source test.
"""
from pathlib import Path
import json, math
import sympy as s
import numpy as np

ROOT=Path(__file__).resolve().parent
q,p,t,a,e=s.symbols('q p t a e', real=True)

def moment(n):
    return 0 if n%2 else s.factorial2(n-1) if n else s.Integer(1)

def gauss_expect(expr, variables):
    pol=s.Poly(s.expand(expr),*variables)
    return s.expand(sum(coef*s.prod(moment(n) for n in powers) for powers,coef in pol.terms()))

def conditional_moment(k,mean,var):
    if k<0:return s.Integer(0)
    return s.expand(sum(s.binomial(k,2*j)*s.factorial2(2*j-1)*mean**(k-2*j)*var**j for j in range(k//2+1)))

h=p*p-1
F=s.Rational(1,10)*(q*q-1)*h
C=s.Rational(1,50)*h*h
tau_minus_C=s.Rational(1,50)*(q*q-1)*h*h
second_current=s.Rational(1,250)*q*q*h**3
base=p+s.Rational(1,20)*h
mean=base+t*F
var=s.Rational(2,3)+(1-t*t)*C
two_riesz=[]
for k in range(1,8):
    direct=s.diff(gauss_expect(conditional_moment(k,mean,var),[q,p]),t)
    one=0 if k<2 else t*k*(k-1)*gauss_expect(tau_minus_C*conditional_moment(k-2,mean,var),[q,p])
    two=0 if k<3 else t*t*k*(k-1)*(k-2)*gauss_expect(second_current*conditional_moment(k-3,mean,var),[q,p])
    assert s.expand(direct-one)==0
    assert s.expand(one-two)==0
    two_riesz.append(dict(test_degree=k,direct_equals_first=True,first_equals_second=True))

# The final same-endpoint normal-form path, with random K=e(m+q), Var(q)=1.
H1=p; H2=h; H3=p**3-3*p; H5=p**5-10*p**3+15*p
m=s.Rational(2,5)
S=-e*e*m*m*(H1+2*H3)
D=-e*e*(H1+2*H3+s.Rational(1,2)*H5)
mean=p+e*m*H2+S+a*D
var=s.Rational(2,3)+a*e*e*H2*H2
endpoint=[]
for k in range(1,8):
    # Only the leading grade-six and first uncancelled grade-nine are needed.
    polynomial=gauss_expect(conditional_moment(k,mean,var),[p])
    derivative=s.diff(polynomial,a)
    lower=[s.expand(derivative).coeff(e,j) for j in range(3)]
    assert all(x==0 for x in lower)
    endpoint.append(dict(test_degree=k,grades_below_nine_cancel=True,
                         alpha9_coefficient=str(s.expand(derivative).coeff(e,3))))

# Noncommuting multidimensional Wick contractions.
rng=np.random.default_rng(7306)
multidim=[]
for dim in [2,3,5]:
    for case in range(12):
        raw=rng.normal(size=(dim,dim,dim))
        T=sum(raw.transpose(perm) for perm in [(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)])/6
        theta=rng.normal(size=dim)
        mean_v=np.einsum('iab,a,b->i',T,theta,theta)
        N=np.einsum('iab,jab->ij',T,T)
        M=np.einsum('iab,jad->ijbd',T,T)
        tilted_cov=2*N+4*np.einsum('ijbd,b,d->ij',M,theta,theta)
        direct=.5*theta@(tilted_cov+np.outer(mean_v,mean_v))@theta
        rhs=.5*(theta@mean_v)**2+2*np.einsum('ijbd,i,j,b,d',M,theta,theta,theta,theta)+theta@N@theta
        error=abs(direct-rhs)
        assert error<1e-9
        multidim.append(dict(D=dim,error=float(error)))

out=dict(second_riesz_exact_tests=two_riesz,same_endpoint_exact_tests=endpoint,
         multidimensional_wick_tests=multidim,
         scope='Exact finite algebra diagnostics. Source, extended-cut and complete-native assertions are mathematical proof obligations, not numerically executed here.')
(ROOT/'second_riesz_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(second_riesz_degrees=len(two_riesz),same_endpoint_degrees=len(endpoint),
                     multidim_cases=len(multidim),max_wick_error=max(x['error'] for x in multidim)),indent=2))
