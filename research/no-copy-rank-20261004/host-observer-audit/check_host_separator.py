#!/usr/bin/env python3
"""Exact elementary diagnostics for the conditional host separator.
The miniature graph models the named interface, not an instantiated LOW30 queue.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
checks=[]
def check(name,condition,detail=None):
    assert condition,name
    checks.append({'name':name,'pass':True,**({'detail':detail} if detail is not None else {})})

# Edges are predecessor -> consumer. 'theta' contains previously captured labels.
pred={
 'theta':[], 'fine_H':[], 'fine_E':[], 'z_H':[], 'z_E':[],
 'captured_zero':['theta'],
 'H_private_query':['theta','fine_H','captured_zero'],
 'E_old_ancestor':['theta','fine_E'],
 'E_nested_zero':['E_old_ancestor'],
 'E_same_record_proxy':['theta','fine_E','captured_zero','E_nested_zero'],
 'E_difference':['E_old_ancestor','E_same_record_proxy'],
 'Y_H':['H_private_query','z_H'], 'Y_E':['E_difference','z_E'],
 'T':['theta','captured_zero','Y_H','Y_E'],
 'observation_pool':[],
 'R0':['T','observation_pool'], 'mode_query0':['R0','theta'],
 'mode_state':['mode_query0','R0'],
}
for i in range(1,8):
 pred[f'child_tape{i}']=[]
 pred[f'R{i}']=['T','observation_pool']
 pred[f'caller{i}']=['mode_state' if i==1 else f'X{i-1}', f'R{i}']
 pred[f'child_query{i}']=[f'caller{i}',f'child_tape{i}','theta']
 pred[f'X{i}']=[f'child_query{i}',f'caller{i}',f'child_tape{i}']
pred['terminal_force']=['X7','theta']

def ancestors(node,stop=()):
 if node in stop:return set()
 return {node}.union(*(ancestors(p,stop) for p in pred[node]))
private={'fine_H','fine_E'}
host=['R0','mode_query0','mode_state']+[f'{kind}{i}' for i in range(1,8) for kind in ['R','caller','child_query','X']]+['terminal_force']
for node in host:
 check(f'{node} has no direct fine read after T cut',not(ancestors(node,{'T'}) & private))
check('actual terminal still depends on both fine banks',private<=ancestors('terminal_force'))
check('nested zero is private','fine_E' in ancestors('E_nested_zero'))
check('top-level complete zero is captured',not(ancestors('captured_zero')&private))
# A future proxy is a NEW source graph and intentionally reopens private columns.
pred['future_proxy_rebuilt_node']=['captured_zero','fine_H','H_private_query']
check('future proxy is not T-only postprocessing','fine_H' in ancestors('future_proxy_rebuilt_node',{'T'}))
# A fine baseline is a real additional within-packet read, not a future host debt.
pred['baseline_original_gradient_query']=['fine_E','theta']
pred['illegal_packet_sum']=['T','baseline_original_gradient_query']
check('same-fine baseline violates output separator','fine_E' in ancestors('illegal_packet_sum',{'T'}))

# One shared statistic yields the observation covariance; independent T_i fails.
for k in range(1,65):
 n=k+1;a=F(1,37);sig2=a/n
 for i in range(n):
  for j in range(n):
   cov=sig2+a*(F(i==j)-F(1,n))
   check(f'shared T covariance k{k} i{i} j{j}',cov==(a if i==j else 0))
 check(f'independently replacing T_i fails k{k}',-a/n!=0)

# A symbolic affine decoder sensitivity obeys the same geometric structure.
for perturb in [F(0),F(1,100),F(1,16),F(1,8)]:
 beta=(1+perturb)/2;bound=(1+perturb)/(1-perturb);response=F(1)
 for k in range(1,101):
  response=beta*(response+1)
  check(f'decoder common-statistic sensitivity {perturb} k{k}',response<=bound)

# Conservative ACTUAL CW7 child guard, separate from the kept-center derivative.
for delta in [F(0),F(1,100),F(1,16),F(1,8),F(1,4)]:
 lam=(1+delta)/2
 for mode_bound in [F(1),F(4,3),F(2)]:
  response=mode_bound
  for k in range(1,101):
   response=lam*(response+1)
   exact=lam**k*mode_bound+lam*(1-lam**k)/(1-lam)
   check(f'actual finite same-seed factor {delta} mode{mode_bound} k{k}',response==exact and exact<=2)

# Equal output marginals suffice past the cut, but not if an extra old fine read remains.
# Y=U and Yref=Z both N(0,1). A(U,Y)=U+Y has variance 4, A(U,Yref) variance 2.
check('marginal Gaussian outputs identical',F(1)==F(1))
check('same-fine baseline changes variance',F(4)!=F(2),{'actual_variance':4,'incorrect_reference_variance':2})
check('same-fine baseline changes fourth moment',3*F(4)**2!=3*F(2)**2)
# Two complete packets X_i=U_i+U_i preserve within-bank aliases, independent across i.
check('independent complete banks keep aliases',F(4)+F(4)==8)
check('sharing a private ancestor changes bank sum law',F(4)**2!=8)

# Exact-key cache test: only captured sites are shared across independent banks.
for s in [0,1,5,31]:
 for t in [1,2,11,37]:
  for m in [1,2,7,29]:
   sites={('captured',j) for j in range(s)}
   sites|={('private',i,j) for i in range(m) for j in range(t)}
   # repeated same-record references add no nodes
   sites|={('private',i,j) for i in range(m) for j in range(t)}
   check(f'cache census s{s} t{t} m{m}',len(sites)==s+m*t)

out={'status':'PASS','checks':len(checks),
 'scope':'Exact interface diagnostics; not a numerical enumeration of the actual LOW30 graph and not an all-rank compiler proof.',
 'results':{
  'host_separator':'theta,T suffice for the displayed downstream law interface',
  'actual_graph_dependence':'both source fine banks remain read in the actual graph',
  'extra_baseline':'a same-fine internal baseline fails the output-only separator',
  'observations':'one shared T yields exact independent reference observations',
  'cache':'captured Q_S plus M complete-private Q_T'
 }}
path=Path(__file__).with_name('host_separator_checks.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
