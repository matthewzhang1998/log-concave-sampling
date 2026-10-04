from fractions import Fraction as F
from pathlib import Path
import json
checks=0
# The owned-X inequality is only a candidate/source-acceptance calculation.
for ti in range(1,250):
 tau=F(ti,1000)
 gout=1-2*tau
 assert gout>2*tau
 for g in [2*tau,(2*tau+gout)/2,gout,F(1)]:
  gain=min(2*tau,g,2-tau,F(7,3)-F(7,3)*tau,F(2),F(4))
  assert gain==2*tau
  assert min(g,F(1))>=gain
  checks+=2
# Stable gain ceiling from the SAME leading low-clock tail.
for ti in range(1,500):
 tau=F(ti,1000)
 stable_gain=min(2*tau,1-2*tau)
 assert stable_gain<=F(1,2)
 checks+=1
# Exact original harmonic-number expression for the optimized certificate.
c=F(0);n=0
for p2 in range(14,201):
 P=F(p2,2);R=P+F(1,2)
 c=max(c+P-F(49,10),((2*R/3)*c if P==7 else (R/P)*(c+P-F(27,5))))
 if R>=F(17,2):
  n=int(2*R-17)
  harmonic=sum((F(1,k) for k in range(17,n+17)),F(0))
  expected=F(78,85)+n-F(54,5)*harmonic
  assert c/R==expected
  checks+=1
# Terminal floor: zero new copy penalty does not reduce inherited ratio.
P=F(15,2);c=F(21,10);ratio=c/P
for _ in range(100):
 R=P+F(1,2);k=F(0)
 c=max(c+k,(R/P)*(c+max(F(0),k-F(1,2))))
 assert c/R==ratio
 P=R;checks+=1
out={'arithmetic_status':'PASS','checks':checks,
 'whole_family_acceptance':'NOT ESTABLISHED',
 'unclosed_ports':['Complete mean-bias correction','Original true-covariance orientation','All successor self/mixture/tail rows','Generalized weaker-curl private input/return','Actual stable source and width return through every recursion'],
 'warning':'Candidate exponent arithmetic does not certify the current owned-X endpoint/source program or an all-order repair queue.'}
Path(__file__).with_name('no_copy_acceptance_arithmetic.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS:',checks,'exact rational arithmetic checks; whole-family acceptance remains NOT ESTABLISHED')
