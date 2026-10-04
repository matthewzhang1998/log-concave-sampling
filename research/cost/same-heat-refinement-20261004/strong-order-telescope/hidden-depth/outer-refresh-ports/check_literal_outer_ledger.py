#!/usr/bin/env python3
"""Exact rational reproduction of the fixed Exact32 covariance ledger.
Assumes separately certified native ports, intrinsic marked cap and restored
child/numerical budgets. Arithmetic alone does not certify those inputs.
"""
from fractions import Fraction as F
from pathlib import Path
import json
OUT=Path(__file__).parent
marks=[F(27,10),F(14,5),F(29,10),F(3),F(31,10),F(16,5),F(33,10),F(17,5),F(7,2),F(71,20)]
rows=[]
def add(group,name,v,k,mult=1):
    grades=(2*v-F(1,2),2*v-F(1,2),2*v+F(3,2)) if k==0 else (2*v+min(F(k,2)-1,F(3,2)),2*v-F(1,2),2*v+k+F(1,2))
    rows.append(dict(group=group,name=name,canonical_power=str(v),rank=k or 'direct',multiplicity=mult,
                     rational_grades=[str(g) for g in grades],grades=grades))
params=[]
for x,y in zip(marks[:-1],marks[1:]):
    u=(F(327,50)-x)/4;e=y-F(7,2)+u+F(1,100)
    params.append((x,u,e,2))
params.append((F(71,20),F(3,4),F(81,100),1))
for x,u,e,m in params:
    group='mixed_'+str(x);c=F(199,200);B=F(52,5);f=4*u+x+c;h=B*u;s=c*u
    for name,v,k in [('parent_energy',f,4),('parent_intrinsic_cap',h,4),('parent_energy_direct',f+s,0),('parent_intrinsic_cap_direct',h+s,0),('child_restoration',F(9),0)]:add(group,name,v,k,m)
    for j in [2,3,4]:add(group,'packet_'+str(2*j),3*j+x-(4*j-7)*e,2,m)
    for j in [2,3,4]:add(group,'packet_root_'+str(2*j+1),3*j+x-(4*j-6)*e,2*j+1,m)
    add(group,'mixture_remainder',F(15,2)+x-4*e,10,m)
    add(group,'covariance_heat',F(19,2),2,m)
J=F(71,20);c=F(199,200);B=F(52,5);u=F(3,4)
Rv=[4*u+J+c,B*u]
for i,v in enumerate(Rv):add('self','R_'+str(i),v,4)
for i,v in enumerate(Rv):
    add('self','T_r_'+str(i),v+c*u,0)
    add('self','T_A_'+str(i),v+c,0)
add('self','U_energy',5*u+J+c,0);add('self','U_cap',(B+1)*u,0)
for name,v,k in [('force11',11*u,4),('force12',12*u,0),('numerical_heat_21over2',2*u+F(21,2),2),('literal_word',4*u+J+1,4),('word_direct',5*u+J+1,0),('word_later_direct',8*u+J+1,0),('numerical_heat_17over2',2*u+F(17,2),2),('order16',16*u,4)]:add('self',name,v,k)
add('pure','prior_fourth_restored',F(3999,200),4);add('pure','prior_direct_restored',F(20),0)
add('tail','retained_tail_conservative',4*u+F(24,5),4);add('tail','order16',4*u+9,4)
add('tail','numerical_heat',F(10),2);add('tail','numerical_direct',F(10),0)
assert len(rows)==152
assert sum(r['multiplicity'] for r in rows)==269
mins=[min(r['grades'][i] for r in rows) for i in range(3)]
assert mins==[F(1607,100),F(141,10),F(22603,1250)],mins
attainers=[[dict(group=r['group'],name=r['name']) for r in rows if r['grades'][i]==mins[i]] for i in range(3)]
for r in rows:del r['grades']
result=dict(status='PASS',structural_categories=len(rows),with_multiplicities=sum(r['multiplicity'] for r in rows),
            minima=[str(z) for z in mins],attainers=attainers,
            qualification='Requires the independent source/proxy and intrinsic-cap proofs; adjustable child orders and numerical budgets are restored first.',rows=rows)
(OUT/'literal_outer_ledger_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
