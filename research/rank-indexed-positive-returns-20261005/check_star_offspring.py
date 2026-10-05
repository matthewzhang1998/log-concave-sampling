#!/usr/bin/env python3
"""Exact finite cumulant jet, never an assertion of a polynomial MGF's convergence."""
import sympy as s
from math import comb, factorial
from pathlib import Path
import json
HERE=Path(__file__).parent
z,b,u,v=s.symbols('theta b u v')
checks=0

def require(v):
    global checks
    assert v;checks+=1

def gm(k):
    return s.Integer(0) if k%2 else s.factorial2(k-1) if k else s.Integer(1)

def shifted(k,a):
    return s.expand(sum(comb(k,j)*(a*z)**(k-j)*gm(j) for j in range(k+1)))

def coefficients(n,K):
    pars=s.symbols('a:'+str(n-2))
    mus=[s.Integer(1)]+[s.prod(shifted(k,a) for a in pars).expand() for k in range(1,K+1)]
    cs=[s.Integer(0)]
    for k in range(1,K+1):
        cs.append(s.expand(mus[k]-sum(comb(k-1,i-1)*cs[i]*mus[k-i] for i in range(1,k))))
    return pars,cs

rows=[]
for n in range(3,10):
    a,c=coefficients(n,2)
    require(s.expand(c[1]-s.prod(a)*z**(n-2))==0)
    require(s.expand(c[2]-(s.prod(1+x*x*z*z for x in a)-s.prod(x*x for x in a)*z**(2*n-4)))==0)
    rows.append({'star_rank':n,'second_log_current_ranks':list(range(4,2*n,2)),'side_zero_first_force_count':6,'side_zero_first_physical_rank':4,'side_zero_first_cycle_count':n-3})

pars,cs=coefficients(5,6)
out=[]
for k in range(1,7):
    p=s.Poly(s.expand(z**(2*k)*cs[k]/s.factorial(k)),z,*pars)
    terms=[]
    for powers,c in p.terms():
        degree,own,*sides=powers
        forces=3*k+sum(sides)
        require(degree>=2*k)
        require(forces>=3*k)
        terms.append({'rational_coefficient':str(c),'physical_rank':degree,'spine_copies':k,'own_public_shifts':own,'side_force_occurrences':sides,'original_force_count':forces})
    out.append({'spine_amplitude_degree':k,'terms':terms})
second=s.expand(z**4*cs[2]/2).subs(dict(zip(pars,(b,u,v))))
expected=(z**4+(b*b+u*u+v*v)*z**6+(b*b*u*u+b*b*v*v+u*u*v*v)*z**8)/2
require(s.expand(second-expected)==0)
require(s.expand(second.subs({u:0,v:0})-(z**4+b*b*z**6)/2)==0)
require(s.expand(z**2*cs[1].subs(dict(zip(pars,(b,u,v))))-b*u*v*z**5)==0)
# Native star graph after all center/public/side noises contract between two spines.
for n in range(3,10):
    vertices=6; edges=4+(n-2); cycles=edges-vertices+1; marks=4
    derivative_surplus=2*((n-1)-1)
    require(cycles==n-3)
    require(derivative_surplus==marks-2+2*cycles)
# Exact dyadic endpoint power diagnostics with tau=2^-k.
from fractions import Fraction
endpoint=[]
for k in range(2,18):
    tau=Fraction(1,2**k); alpha=tau; Delta=tau*tau
    # Choose sigma=tau and a center-panel weight Delta, ignoring fixed readout factors.
    rho=alpha*Delta/tau**3
    caller=alpha*alpha*Delta/tau**4
    feedback=alpha**6*Delta**2/tau**6
    require(rho==1 and caller==1 and feedback==alpha**4)
    endpoint.append({'k':k,'rho_root':str(rho),'center_caller_budget':str(caller),'feedback_over_alpha4':str(feedback/alpha**4)})
report={'assertions':checks,'scope':'Formal finite star jets and literal ideal-graph scale tests; no W2 lower bound or native execution','all_rank_second_coefficients':rows,'rank5_all_descendants_through_spine_degree6':out,'endpoint_normalization_checks':endpoint}
(HERE/'star_offspring_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'assertions':checks,'rank5_terms':[len(i['terms']) for i in out],'exact_rank5_second_coefficient':str(second),'scope':report['scope']},indent=2))
