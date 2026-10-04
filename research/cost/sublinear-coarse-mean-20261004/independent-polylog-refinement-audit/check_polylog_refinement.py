#!/usr/bin/env python3
"""Independent deterministic constant/census checks; no provider simulation."""
from pathlib import Path
import hashlib, itertools, json, math
from scipy.special import log_ndtr, logsumexp
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'POLYLOG-SCALAR-MEAN-LAW-REFINEMENT.md'
C0=2**20
rows=[]
for rho,delta in itertools.product([1e-6,.1,1.,10.,1e3,1e6],[.99*math.exp(-2),1e-3,1e-6,1e-20,1e-100]):
    L=math.log(math.e+rho)-math.log(delta)
    R=math.sqrt(64*L); eps=1/(16*R)
    K=math.ceil(8*L); B=math.ceil(8*L)
    M=math.ceil(C0*(1+rho*rho)*R*R*L)
    pilot_exp=M/(8192*R*R*rho*rho)
    # Squared positive-exponential envelope integrates exactly to this bound.
    log_tail=.5*math.log(2)+1.5*eps*eps+.5*log_ndtr(-(R-2*eps))
    log_rem=eps*eps/2+(K+1)*math.log(eps)-.5*math.lgamma(K+2)
    log_profile=logsumexp([log_rem,math.log(2)+log_tail])
    log_factory=.5*(math.log(2*math.sqrt(3))+log_profile)
    log_cap=.5*math.log(5)-B*math.log(3)/2
    log_good=logsumexp([log_factory,log_cap,math.log(4*rho)-8*L])
    log_bad_cost=math.log(16*rho*rho+8)-8*L
    log_total=.5*logsumexp([2*log_good,log_bad_cost])
    literal=(K+1)*M
    cost_ratio=literal/((1+rho*rho)*L**3)
    assert eps*R+eps*eps/2<math.log(1.25)
    assert pilot_exp/L>=128-1e-9
    assert log_tail<=math.log(2)-R*R/8
    assert log_rem<=math.log(2)-K*math.log(16)
    assert log_cap<=.5*math.log(5)-4*L
    assert log_total<=math.log(delta)
    rows.append(dict(rho=rho,delta=delta,L=L,R=R,epsilon=eps,K=K,B=B,M=M,literal_complete_records=literal,literal_over_claimed_order=cost_ratio,pilot_exponent_over_L=pilot_exp/L,log_uniform_tail_upper_bound=float(log_tail),log_total_W2_over_sigma_upper_bound=float(log_total),log_upper_bound_over_delta=float(log_total-math.log(delta))))
# Scalar quadratic covariance discrepancy is an exact rational identity.
quadratics=[]
for B,sigma in itertools.product([0.,.01,.1,.5],[.001,.01,.1]):
    buffered=(1-B/2)**2+sigma*sigma
    posterior=1/(1+B)
    exact_debt=sigma*sigma-B*B*(3-B)/(4*(1+B))
    assert abs(buffered-posterior-exact_debt)<3e-16
    quadratics.append(dict(B=B,sigma=sigma,buffered_variance=buffered,posterior_variance=posterior,exact_variance_debt=exact_debt))
result=dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),case_count=len(rows),min_pilot_exponent_over_L=min(r['pilot_exponent_over_L'] for r in rows),max_record_count_over_claimed_order=max(r['literal_over_claimed_order'] for r in rows),max_W2_upper_bound_over_sigma_delta=math.exp(max(r['log_upper_bound_over_delta'] for r in rows)),notes='Checks deterministic analytic upper bounds in floating point. No enormous provider batch was run; no external manuscript or data was used.',cases=rows,quadratic_checks=quadratics)
(HERE/'diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['cases','quadratic_checks']},indent=2))
