#!/usr/bin/env python3
"""Exact grade identities plus deterministic common-carrier algebra checks.
Diagnostics verify formulas, not the imported native CURRENT theorem.
"""
from fractions import Fraction as F
from itertools import product
import json, math, pathlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
checks=[]
# Bank-tree additivity and conditional return, with independent positive types.
instances=0
for beta in [F(1,16),F(1,8),F(1,4)]:
 for gamma in [F(1,16),F(1,8),F(1,7)]:
  for r in range(2,7):
   ts=[(F(5+i),3+i,1+2*i) for i in range(r)]
   psi=sum(a-beta*n-gamma*(k+2) for a,n,k in ts)
   aa=sum(t[0] for t in ts); nn=sum(t[1] for t in ts)
   kk=sum(t[2] for t in ts)+2*(r-1)
   assert aa-beta*nn-gamma*(kk+2)==psi
   for h in range(2,8):
    for n,k in [(2,0),(4,4),(8,8),(9,7)]:
     d=F(3,2); oldpsi=d-2*gamma
     aa=beta*n+gamma*k+h*d
     assert aa-beta*n-gamma*(k+2)==h*oldpsi+2*gamma*(h-1)
     instances+=1
checks.append({'name':'mixed and conditional potential identities','instances':instances,'passed':True})
# Rank-dependent exact one-stage admission.
ranks=[]
for n in range(3,33):
 gamma=F(1,n-1); beta=F(1,2*(n-1)); K=n-2
 target=F(n)+gamma; own=2*n-2*beta*(n-1)-2*gamma*K
 selfg=2*n-gamma*(2*K+1); first=F(n)-beta*(n-1)-gamma*(K+1)
 assert own>target and selfg>target and first>1
 ranks.append({'n':n,'gamma':str(gamma),'target':str(target),'own':str(own),'self':str(selfg),'root_caller':str(first)})
assert ranks[5]['target']=='57/7' and ranks[5]['own']=='93/7'
assert F(57,7)-F(65,8)==F(1,56)
checks.append({'name':'all n=3..32 one-stage grade tests','passed':True,'rank8':ranks[5]})
# Gaussian scalar-block pullback, several unequal public dimensions.
rng=np.random.default_rng(827104)
v=0.37
Ls=[]; Qs=[]
for p in [2,4,9,11]:
 l=rng.normal(size=(1,p)); l*=math.sqrt(v)/np.linalg.norm(l)
 Pi=np.eye(p)-l.T@l/v
 Q=np.hstack([l.T/math.sqrt(v),Pi])
 assert np.linalg.norm(Pi@Pi-Pi)<1e-12
 assert np.linalg.norm(Q@Q.T-np.eye(p))<1e-12
 assert np.linalg.norm(l@Pi)<1e-12
 assert abs((l@l.T).item()-v)<1e-12
 Ls.append(l); Qs.append(Q)
for l,m in product(Ls,Ls):
 # common carrier cross-covariance, from their common G0 columns.
 assert np.linalg.norm((l.T/math.sqrt(v))@(m/math.sqrt(v))-l.T@m/v)<1e-12
checks.append({'name':'common-carrier projections/coisometries/cross-covariance','passed':True,'public_dimensions':[2,4,9,11]})
# Scalar Gaussian model: exact bank-mean W2 <= E*(old first+self first/2).
ratios=[]
for delta,eps in product([0.0,1e-5,0.01,0.1,0.5],[1e-5,0.001,0.1]):
 w=abs(math.sqrt(1+(delta+eps)**2)-math.sqrt(1+delta**2))
 bound=eps*(delta+eps/2)
 assert w<=bound*(1+1e-7)+1e-15
 ratios.append(w/bound if bound else 0)
checks.append({'name':'sharp scalar Riesz mixing model','passed':True,'max_ratio':max(ratios)})
# The advertised simple-ledger ceiling.
beta=F(1,14); a=N=8; K=6
ceiling=F(a)+(F(a)-beta*N)/(K+2)
assert ceiling==F(125,14)
# Raw triple-cubic force-hit source count maxima, retaining copy hit budgets 2,1,1.
maxcalls=0; cases=0
for hits0 in product(range(3),repeat=2):
 for hit1 in range(3):
  for hit2 in range(3):
   ks=[]
   for hits in [hits0,(hit1,),(hit2,)]:
    local=[0,1,0]
    for hit in hits: local[hit]+=1
    ks+=local
   assert sum(ks)==7 and max(ks)<=3
   maxcalls=max(maxcalls,sum(2**(k+1) for k in ks)); cases+=1
assert maxcalls==44 and cases==81
checks.append({'name':'triple-cubic raw active-source maximum','passed':True,'force_hit_assignments_per_tree':cases,'max_raw_VALUE_calls_per_history_node':maxcalls,'all_three_trees_hit_count':243})
result={'status':'PASS','disclaimer':'Formula/algebra diagnostics only; native CURRENT, heat and finite quadrature contracts remain their cited proofs.','checks':checks,'rank_grade_table':ranks,'simple_heat_discard_ceiling_rank8':str(ceiling)}
(ROOT/'mixed-queue-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'potential_instances':instances,'rank8':ranks[5],'triple_raw_calls_max':maxcalls},indent=2))
