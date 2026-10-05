#!/usr/bin/env python3
"""Author diagnostics for the actual fixed-order endpoint ledger.
The imported full native compilers are NOT numerically instantiated here.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import hashlib,json,math,random
import numpy as np
P=Path(__file__).resolve().parent
checks={}
def ck(value,key):
    assert value,key
    checks[key]=checks.get(key,0)+1
h=F(1,2);rho=F(3,4);k=F(4)
ck(2*(k-2)/5==F(4,5),'new_terminal_exponent')
ck(1+F(4,5)==F(9,5),'minimum_buffered_alpha_exponent')
ck(2*k-1==7,'new_local_s_power')
ck(k+F(1,2)==F(9,2),'physical_readout_exponent')
# Entire-bank physical variance budget and negative skew coefficient.
alloc=[h*h,h*h*F(1,4),h*h*F(1,4),F(1,8),F(1,4)]
ck(sum(alloc)==F(3,4),'conditional_baseline')
ck(h*h+sum(alloc)==1,'whole_unit_carrier')
ck(sum(alloc[:-1])==F(1,2),'visible_reference_gap')
ck(alloc[-1]/sum(alloc)==F(1,3),'untouched_keep_fraction')
ck(-(h**3)/6==-F(1,48),'negative_skew_coefficient')
ck((h*h)/F(1,8)==2,'cubic_radius_squared')
# Actual fixed graph metadata, including full completed M rather than bare source.
g=json.loads((P/'ENDPOINT-GRAPH.json').read_text())
ck(g['banks'][0]['input_type']=='completed own-mean LAW','mean_bank_type')
ck(g['banks'][0]['native_orders']=={'prefix_mean':4,'tail_mean_and_coefficient_means':6,'corrected_negative_fourth_paths':6,'cubic_mixed_true_quartic':5,'final_own_mean':4},'mean_bank_fixed_orders')
ck('d_M_completed' in g['root_dimension'],'completed_mean_root_dimension')
ck('dyadic' in g['banks'][0]['owned_roots'],'new_prefix_banks_in_graph')
ck(g['result']['leading_logs_retained'] and not g['result']['native_compilers_executed'],'honest_scope')
for bank,v in zip(g['banks'],alloc):ck(F(bank['baseline_after_readout'])==v,'graph_bank_allocation')
# Exact rational bridge identities, source-zero entering coefficient and mode scale.
for i in range(301):
 r2=F(i,302);s2=1-r2;delta=s2/4;t2=r2+delta
 v0=delta/t2;v=(1-t2)*delta/(t2*s2)
 ck(v==v0*(1-h*h),'exact_bridge_noise')
 ck(delta**2/(t2*s2)==v0*h*h,'exact_bridge_signal')
 ck((1-t2+delta)/s2==1,'exact_reference_contraction')
 ck(1-r2/t2==v0,'source_zero_stage_variance')
 ck(s2>0 and v>0 and v0>0,'positive_stage_gaps')
 ck((delta/s2)*s2**4==delta*s2**3,'mode_residual_readout')
# Finite schedule, all local alphas, and exact propagated contraction weights.
for exponent in np.linspace(.001,80,321):
 logA=-float(exponent)*math.log(10);A=math.exp(logA)
 cutoff=math.exp(.8*logA);s2=1.;states=[];total=0.;local_alphas=[]
 while s2>cutoff:
  prev=s2;s2*=.75;delta=prev/4
  states.append((math.sqrt(max(0.,1-prev)),math.sqrt(1-s2)))
  local_alphas.append(A*prev);total+=delta*prev**3.5
 ck(s2<=cutoff and s2>.75*cutoff,'first_terminal_stop')
 ck(2.5*math.log(s2)<=2*logA+1e-12,'terminal_grade_four_log_check')
 ck(all(math.log(a)>1.8*logA-1e-12 and a<=A*(1+1e-13) for a in local_alphas),'all_buffered_alpha_ranges')
 ck(1.8*logA+math.log(.75)<math.log(A*s2)+1e-12<=1.8*logA+1e-12,'terminal_alpha_range')
 ck(total<=(1/(4*(1-.75**4.5)))*(1+1e-13),'exact_weighted_intrinsic_sum')
 ck(len(states)==math.ceil(.8*logA/math.log(.75)),'finite_stage_count')
 for pos in sorted({0,len(states)//2,len(states)-1}):
  weight=states[-1][1]
  for m in range(pos+1,len(states)):weight*=states[m][0]/states[m][1]
  ck(abs(weight-states[pos][1])<=2e-13,'telescoping_reference_weight')
 # New mean cutoff is checked at each alpha, with an illustrative frozen K.
 K=2.;
 if A*K*K<=.5:
  ck(all(a*a*K*K<=a/2*(1+1e-13) for a in local_alphas),'all_stage_mean_cutoffs')
# Exact exponent ledger for each unsimplified native prior row, including K powers.
# a^n (w^-p+eta^(1-p)), w=a, eta=a²K².
rows=[(F(7),F(5,2),F(9,2),F(4),F(-3)),
      (F(6),F(7,6),F(29,6),F(17,3),F(-1,3)),
      (F(6),F(2),F(4),F(4),F(-2)),
      (F(6),F(3,2),F(9,2),F(5),F(-1)),
      (F(7),F(3,2),F(11,2),F(6),F(-1)),
      (F(11,2),F(5,4),F(17,4),F(5),F(-1,2))]
for n,p,broad,cut,kpow in rows:
 ck(n-p==broad and n+2*(1-p)==cut and 2*(1-p)==kpow,'mean_integrated_native_prior')
 ck(broad>=4 and cut>=4,'mean_integrated_grade_floor')
ck(F(3)+F(2,2)==4,'raw_endpoint_grade_with_K')
ck(F(2)+F(2)==4,'prefix_smoothing_grade')
# Explicit negative control: normalized radius small does not license prior domination.
for A in [1e-9,1e-12,1e-15]:
 K=100.;u=A*A*K*K
 ck(A/math.sqrt(u)<=.011,'small_literal_normalized_radius')
 ck((A**6/u**2.5)/(A**4/u)>1,'pointwise_prior_domination_forbidden')
 # Row is still paid by the actual integrated bound.
 ck(math.isclose(A**7*(A**-2.5+u**-1.5),A**4.5+A**4/K**3,rel_tol=2e-14),'integrated_prior_identity')
# Coherent exact linear-force OU history, not independently sampled substitutions.
def cm(m,n):
 return F(math.factorial(m+n+1),2**(m+n+2)*math.factorial(m)*math.factorial(n))*(F(1,m+1)+F(1,n+1))-F(1,2**(m+n+2))
for a in [F(i,500) for i in range(1,151)]:
 means=[F(1,2),F(1,4),F(1,8)];coef=[a,-a*a,a**3]
 m3=sum(c*m for c,m in zip(coef,means))
 cov3=sum(coef[i]*coef[j]*cm(i,j) for i in range(3) for j in range(3))
 varP=(1-m3)**2+cov3
 varTarget=1/(1+a)
 # Actual joined conditional reference mean and covariance recover buffered P3.
 vjoined=h*h*(1-m3)**2+sum(alloc)+h*h*cov3
 vref=h*h*varP+1-h*h
 ck(vjoined==vref,'coherent_linear_complete_bank_join')
 ck(abs(float(varP-varTarget))<=2*float(a**4),'canonical_P3_grade_four')
 ck(varP!=varTarget,'no_higher_order_target_claim')
 posterior=a*h/(1+a*(1-h*h))
 ck(m3!=posterior,'posterior_mean_wrong_target')
 # Skew is exactly zero for this Gaussian history; a cubic offset cannot repair target P3.
 ck(cov3>0 and vjoined>F(3,4),'linear_covariance_positive')
# Noncommuting anisotropic Gaussian regression and output-slot cancellation.
rng=np.random.default_rng(202610050904)
for d in range(1,9):
 for _ in range(30):
  B=rng.normal(size=(d,d));Sig=.01*B@B.T/max(d,1);C=.5*np.eye(d)+.25*Sig;Ci=np.linalg.inv(C)
  b2=.02+rng.random()*.09;S=rng.normal(size=d)
  mean=b2*Ci@S;V=b2*np.eye(d)-b2*b2*Ci
  lhs=(np.outer(mean,mean)+V)/(b2*b2)-np.eye(d)/b2
  rhs=np.outer(Ci@S,Ci@S)-Ci
  ck(np.allclose(lhs,rhs,atol=1e-11,rtol=1e-11),'anisotropic_Hermite_trace_regression')
  ck(np.linalg.eigvalsh(C)[0]>=.5-1e-13 and np.linalg.eigvalsh(V)[0]>0,'positive_regression_gaps')
  N=rng.normal(size=(d,d,d));N=(N+N.transpose(0,2,1))/2
  sym=sum(N.transpose(p) for p in permutations(range(3)))/6
  x=rng.normal(size=d)
  ck(abs(np.einsum('ijk,i,j,k',N-sym,x,x,x))<1e-9,'zero_leading_symmetrization_current')
# Source pins and all direct sealed manifest entries, with one explicit known stale audit.
for pin in json.loads((P/'INPUT-PINS.json').read_text())['inputs']:
 ck(hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256'],'input_pin')
exceptions=[]
for root in [Path('/workspace/shared/dyadic-prefix-law-join-20261005'),Path('/workspace/shared/positive-endpoint-mean-join-20261005')]:
 manifest=json.loads((root/'MANIFEST.json').read_text())
 entries=manifest['files']; entries=entries.items() if isinstance(entries,dict) else [(e['path'],e['sha256']) for e in entries]
 for rel,expected in entries:
  actual=hashlib.sha256((root/rel).read_bytes()).hexdigest()
  if actual!=expected:exceptions.append({'path':str(root/rel),'manifest_sha256':expected,'actual_sha256':actual})
  else:ck(True,'sealed_manifest_matched_entry')
ck(len(exceptions)==1 and exceptions[0]['path'].endswith('positive-endpoint-mean-join-20261005/independent-audit/INDEPENDENT-ENDPOINT-JOIN-AUDIT.md'),'disclosed_stale_manifest_entry_only')
result={'status':'PASS','assertions':sum(checks.values()),'checks':checks,'native_compilers_executed':False,'provenance_exceptions':exceptions,'scope':'Exact rational positive bank/bridge/grade identities, every-stage finite schedules, unsimplified native-prior integration, coherent linear Gaussian reference, anisotropic Hermite diagnostics, graph/root/type and exact source-pin checks. Numerical examples do not instantiate a full native compiler or certify unspecified numerical guards.'}
(P/'endpoint_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
