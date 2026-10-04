# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
"""NEW RECONSTRUCTED arithmetic checks; not an analytic source proof."""
from fractions import Fraction as F
from pathlib import Path
import json
checks=0
# Generic small-power seed budget, followed by uniform local suffix costs.
for n in range(1,41):
 for eps in [F(1,100),F(1,10),F(1,3)]:
  k0=F(3,4)*eps
  P=F(15,2);c=5*k0
  budget=F(2,3)*k0
  for i in range(1,n+1):
   R=P+F(1,2)
   ki=eps*P/(2*n)
   c=(R/P)*(c+ki)
   budget+=ki/P
   assert c/R==budget
   P=R;checks+=1
  assert c/P==eps
  checks+=1
# Conditional two-index, unit-step ledger and cost arithmetic.
I=F(7);B=F(13,2);c=F(0);unit=[]
for j in range(1,20):
 R=I+1;Q=I+F(1,2)
 assert R==min(I+1,B+2,Q+F(1,2))
 Bnew=min(I,B+1)
 assert Bnew==R-1
 c=(I-F(22,5) if j==1 else max(c+I-F(22,5),(R/I)*(c+I-F(27,5))))
 unit.append({'I':str(R),'B':str(Bnew),'cost':str(c)})
 I=R;B=Bnew;checks+=2
assert unit[0]['cost']=='13/5'
assert unit[1]['cost']=='31/5'
assert unit[2]['cost']=='98/9'
# Reported one-step moderate retuning ledger, strictly below its boundaries.
for di in range(1,333):
 d=F(di,1000)
 stable_kappa_exp=2-d
 assert 2-d+stable_kappa_exp>3+d
 assert 4-F(3,2)*d>3+d
 assert 4-d>3+d
 checks+=3
for di in range(1,250):
 d=F(di,1000)
 assert 2-d+F(3,2)>3+d
 checks+=1
out={'status':'PASS','checks':checks,'version':'NEW RECONSTRUCTED',
 'scope':'Rational epsilon budgets, conditional unit-step ledger, and reported one-step retuning arithmetic; no all-order source acceptance.',
 'unit_step_candidate':unit,'whole_family_acceptance':'NOT ESTABLISHED'}
Path(__file__).with_name('recovery_additional_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS:',checks,'new additional arithmetic checks')
