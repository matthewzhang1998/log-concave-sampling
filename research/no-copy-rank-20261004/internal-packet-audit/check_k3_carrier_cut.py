import math,json
from pathlib import Path
import numpy as np
weights=np.array([1/6]*8+[-1/3]); count=0
assert abs(weights.sum()-1)<1e-15
assert abs(weights@weights-1/3)<1e-15
rows=[]
for a in [1/4,1/8,1/16,1/32,1/64]:
 for m in [.4,.5,.6]:
  r=a; eps=a**.9; b=r*a*a*m**3*eps
  L=np.zeros((3,3));L[0,2]=b
  assert np.array_equal(L@L,np.zeros((3,3)));count+=1
  for alpha in [1.,a**.1,a,a*a,a**5]:
   for ell in [a,2*a,.9]:
    for s in [.1,.25]:
     p=np.array([.7,-.4,1.3]);z=np.array([.1,.2,.3])
     A=alpha*L
     I=A@p/ell
     action=ell**2/(2*s)*((A/ell)@(z+s*I)-(A/ell)@(z-s*I))
     assert np.max(abs(action))<1e-30;count+=1
   va=alpha*alpha*b*b/3
   # Rationalized formulas avoid subtractive cancellation at tiny heat.
   err_a=va/(math.sqrt(1+va)+1)
   v=b*b/3;err=v/(math.sqrt(1+v)+1)
   assert b*b/7<=err<=b*b/6*(1+1e-14);count+=1
   rows.append(dict(a=a,m=m,alpha=alpha,marginal_error_before=err_a,error_after=err,signal_variance_after=v))
out=dict(status='PASS',checks=count,weight_square=float(weights@weights),fixtures=rows,scope='Exact linear K3 and literal VALUE-action/carrier covariance checks; no generic impossibility theorem.')
Path(__file__).with_name('k3_carrier_cut_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='fixtures'},indent=2))
