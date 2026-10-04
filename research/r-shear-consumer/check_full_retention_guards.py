from fractions import Fraction as F
from pathlib import Path
import json, numpy as np
count=0; examples=[]
def ck(x):
 global count
 count+=1
 if not x: raise AssertionError(count)
for gi in range(1,21):
 g=F(gi,20)
 for zi in range(1,41):
  z=F(zi,200)
  for ti in range(1,81):
   gamma=F(ti,100)
   if 2*z<gamma<g-3*z:
    ck(5*z<g)
    ck(gamma-2*z>0)
    ck(1+g-3*z-gamma>1)
    ck(1+g-z-2*gamma>0)
    ck(1+z>1)
    ck(F(4,3)+z>1)
    ck(1+g-z-gamma>1)
    ck(F(7,3)>1+g+z)
    ck(3+2*g-z>1+g+z)
# Explicitly distinguish law eligibility from retained graph eligibility.
g,z,gamma=F(1),F(1,5),F(1,5)
ck(3*z+gamma<g);ck(not gamma>2*z)
for g in [F(1,4),F(1,2),F(1)]:
 z=g/10;gamma=3*g/10
 ck(2*z<gamma<g-3*z)
 examples.append({'g':str(g),'zeta':str(z),'gamma':str(gamma),'cross_edge_exponent':str(gamma-2*z),'side_outgoing_exponent':str(1+g-3*z-gamma)})
rng=np.random.default_rng(4201)
for m in [1,2,10,50]:
 w=rng.uniform(.001,.1,m);W=float(w.sum());vold=.5;vnew=.5;v0=vnew/(4*max(1,W))
 b=np.sqrt(v0)/8;c=-np.sqrt(v0)/10
 carrier=vold+float(w.sum())*(b*b+c*c+(v0-b*b-c*c))+(vnew-v0*W)
 ck(abs(carrier-1)<1e-14)
 ck(vnew-v0*W>=.75*vnew-1e-14)
 ck(v0-b*b-c*c>0)
out={'status':'PASS','checks':count,'examples':examples,'scope':'Exact rational source/law versus retention guards, cross-width penalty, outward masses and source-zero variance row. Does not substitute for the finite original-gradient path proof.'}
Path(__file__).with_name('full_retention_guard_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
