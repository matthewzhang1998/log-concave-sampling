#!/usr/bin/env python3
"""Independent algebra/entropy diagnostics; not an execution of LOW30."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.special import roots_hermitenorm, logsumexp

HERE = Path(__file__).resolve().parent
checks = 0

def check(ok, message):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(message)

# Exact exponent vectors in A,u,kappa,alpha,beta,M.
def add(*terms):
    return tuple(sum(t[i] for t in terms) for i in range(6))
def mul(t, power):
    return tuple(x * F(power) for x in t)
def mon(*vals):
    return tuple(F(v) for v in vals)
r = mon(1, F(-1,2), 1, 0, 0, 0)
delta = mon(1,0,0,1,0,0)
a = mon(1,0,0,0,1,0)
M = mon(0,0,0,0,0,1)
readout = mon(0,F(1,2),0,0,0,0)
mu = r
actual = [
    add(readout, r, mul(add(M,r),6)),
    add(readout,mul(r,3),delta),
    add(readout,mul(r,2),delta,mu),
    add(readout,mul(r,2),a,delta),
    add(readout,mul(r,4),delta,mul(mu,F(-1,2))),
    add(readout,mul(r,4),delta),
]
expected = [
    mon(7,-3,7,0,0,6),
    mon(4,-1,3,1,0,0),
    mon(4,-1,3,1,0,0),
    mon(4,F(-1,2),2,1,1,0),
    mon(F(9,2),F(-5,4),F(7,2),1,0,0),
    mon(5,F(-3,2),4,1,0,0),
]
for j, (x,y) in enumerate(zip(actual,expected)):
    check(x==y,f'exact physical error monomial {j}')
first = add(readout,mul(r,F(3,2)))
check(first==mon(F(3,2),F(-1,4),F(3,2),0,0,0),'exact residual first')
for p,q,*_ in actual:
    # Difference from A^4/u at u=A^(3/2).
    check(p-4+F(3,2)*(q+1)>=0,'endpoint dominance')
for b in range(2,12):
    difference = F(b-4)+F(3,2)*F(3-b,2)
    check(difference == F(b-7,4),'general baseline exponent')
    check((difference>=0)==(b>=7),'minimum integer order')
check(F(3,2)-F(2,4)==1,'first threshold u=A^2')
check(F(3,2)-F(3,8)==F(9,8),'first at u=A^(3/2)')
check(1-F(3,4)==F(1,4),'radius at error endpoint')
check(F(3,2)*F(1,4)==F(3,8),'reserve radius at error endpoint')
check(1-F(2,5)==F(3,5),'radius at delayed-join buffer')
check(F(3,2)*F(3,5)==F(9,10),'reserve radius at delayed-join buffer')
check(F(3,2)-F(1,5)==F(13,10),'physical reserve first at delayed-join buffer')

# Independent numeric substitutions with literal padding, retaining public factors.
max_substitution_relative_error = 0.0
for A in np.logspace(-12,-2,11):
    for exponent in (0.0,0.4,0.8,1.2,1.5):
        u = A**exponent
        for kappa in (1.0,2.0,8.0):
            alpha,beta,m = 2.0,3.0,4.0
            radius=kappa*A/math.sqrt(u)
            pad=radius
            lhs=math.sqrt(u)*np.array([
                radius*(m*radius)**6,
                radius**3*alpha*A,
                radius**2*alpha*A*pad,
                radius**2*beta*A*alpha*A,
                radius**4*alpha*A/math.sqrt(pad),
                radius**4*alpha*A,
            ])
            rhs=np.array([
                m**6*kappa**7*A**7/u**3,
                kappa**3*alpha*A**4/u,
                kappa**3*alpha*A**4/u,
                kappa**2*beta*alpha*A**4/math.sqrt(u),
                kappa**3.5*alpha*A**4.5/u**1.25,
                kappa**4*alpha*A**5/u**1.5,
            ])
            rel=np.max(np.abs(lhs/rhs-1))
            max_substitution_relative_error=max(max_substitution_relative_error,float(rel))
            check(rel<1e-12,'numeric six-term substitution')
            first_direct=math.sqrt(u)*(radius+radius**2/math.sqrt(pad))
            first_formula=kappa*A+kappa**1.5*A**1.5/u**.25
            check(abs(first_direct/first_formula-1)<1e-12,'numeric complete first')
        for c in (0.2,0.5,1.0,2.0):
            u1=c*A**1.5
            if u1>1: continue
            ratios=np.array([A**3/u1**2,1,1,math.sqrt(u1),A**.5/u1**.25,A/u1**.5])
            upper=np.array([c**-2,1,1,1,c**-.25*A**.125,c**-.5*A**.25])
            check(np.all(ratios<=upper*(1+2e-14)),'all six error ratios')

# Smooth covariance C(g)=A^2(1+.5*tanh(g)); HS Lipschitz .5*A^2.
# For independent coordinate copies the D-dimensional law is a tensor product,
# its covariance-to-HS Lipschitz is unchanged, and all entropies add exactly.
g,w=roots_hermitenorm(112)
w=w/math.sqrt(2*math.pi)
logw=np.log(w)
n,wn=roots_hermitenorm(96)
wn=wn/math.sqrt(2*math.pi)
log2pi=math.log(2*math.pi)
concentration=[]
for theta in (-12.,-5.,-1.,-.1,0.,.1,1.,5.,12.):
    lm=float(logsumexp(logw+theta*np.tanh(g)))
    cap=theta**2/2
    check(lm<=cap+2e-13,'sampled Gaussian concentration')
    concentration.append({'t_times_L':theta,'log_mgf':lm,'bound':cap})

def mixture_diagnostics(A, base):
    C=A*A*(1+.5*np.tanh(g))
    meanC=float(w@C)
    fluct=C-meanC
    energy2=float(w@(fluct*fluct))
    lip=.5*A*A
    S=base+C
    sigma=base+meanC
    y=(np.sqrt(S)[:,None]*n[None,:]).reshape(-1)
    joint_weights=(w[:,None]*wn[None,:]).reshape(-1)
    log_components=-.5*(log2pi+np.log(S)[None,:]+y[:,None]**2/S[None,:])
    joint=log_components+logw[None,:]
    logp=logsumexp(joint,axis=1)
    logq=-.5*(log2pi+math.log(sigma)+y*y/sigma)
    posterior=np.exp(joint-logp[:,None])@fluct
    post2=float(joint_weights@(posterior*posterior))
    # The component generating y has S indexed by the outer g-node.
    own_log=-.5*(log2pi+np.log(S)[:,None]+n[None,:]**2)
    info=float(joint_weights@(own_log.reshape(-1)-logp))
    kl=float(joint_weights@(logp-logq))
    component_kl=float(w@(.5*(S/sigma-1-np.log(S/sigma))))
    check(info>=-1e-12,'nonnegative mutual information')
    check(kl>=-1e-12,'nonnegative mixture KL')
    check(abs(component_kl-info-kl)<2e-11,'information chain decomposition')
    check(info<=component_kl+2e-12,'mutual information Gaussian-reference bound')
    check(component_kl<=energy2/(4*base**2)+2e-12,'logdet Bregman energy bound')
    check(post2<=2*lip**2*info+2e-12*max(1,lip**2),'posterior entropy bound')
    return {'A':A,'base_variance':base,'lip':lip,'energy_squared':energy2,
            'mixture_KL':kl,'mutual_information':info,'component_KL':component_kl,
            'posterior_mean_energy':post2,
            'posterior_entropy_bound':2*lip**2*info}

entropy=[]
for A in (.05,.2,.5):
    for exponent in (.4,1.,1.5,2.,2.5):
        u=A**exponent
        before=mixture_diagnostics(A,u/2)
        after=mixture_diagnostics(A,u)
        bound=4*after['lip']**2*after['energy_squared']/u**4
        check(after['mixture_KL']<=bound+2e-12,'fresh equal keep final KL bound')
        for D in (1,100,10**9):
            check(D*after['mixture_KL']<=D*(bound+2e-12),'tensorized dimension test')
        entropy.append({'A':A,'total_base_variance':u,'before_keep':before,
                        'after_keep':after,'final_KL_bound':bound,
                        'dimension_1e9_bound':10**9*bound})

sources=[
    Path('/workspace/shared/law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md'),
    Path('/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex'),
    Path('/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/same-carrier-feedback/GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md'),
    HERE/'INDEPENDENT-RETUNED-GRAM-AUDIT.md',
    Path(__file__),
]
result={
    'status':'passed','assertions':checks,
    'scope':'Independent exact algebra and smooth-mixture numerical diagnostics; not native compiler execution or numerical admission.',
    'max_substitution_relative_error':max_substitution_relative_error,
    'physical_error_exponents':[[str(v) for v in e] for e in actual],
    'exponent_order':['A','u','kappa','alpha','beta','M'],
    'concentration_diagnostics':concentration,
    'entropy_diagnostics':entropy,
    'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
}
(HERE/'retuned_gram_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','assertions','scope','max_substitution_relative_error']},indent=2))
