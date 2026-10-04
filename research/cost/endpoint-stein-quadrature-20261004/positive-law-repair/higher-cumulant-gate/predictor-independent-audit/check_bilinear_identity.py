#!/usr/bin/env python3
"""Exact Wick audit of the two-variable source identity, not a VALUE compiler."""
import sympy as s
from functools import lru_cache
from pathlib import Path
import json
out=Path(__file__).resolve().parent
z,x,y,r,t=s.symbols('z x y r t',real=True)
C=s.Matrix([[1,r,r*t],[r,1,t],[r*t,t,1]])
@lru_cache(None)
def moment(exps):
 if sum(exps)==0:return s.Integer(1)
 if sum(exps)%2:return s.Integer(0)
 i=next(j for j,e in enumerate(exps) if e)
 p=list(exps);p[i]-=1;v=0
 for j in range(3):
  if p[j]:
   count=p[j];p[j]-=1;v+=count*C[i,j]*moment(tuple(p));p[j]+=1
 return s.expand(v)
def E(poly):
 return s.expand(sum(coef*moment(exps) for exps,coef in s.Poly(s.expand(poly),z,x,y).terms()))
def H(n,q):return s.hermite_prob(n,q)
checks=0;rows=[]
for i,j in ((1,2),(2,3),(2,5),(3,4),(3,3),(4,4),(4,5),(5,5)):
 gi=s.diff(H(i,x),x);gj=s.diff(H(j,x),x)
 for n in range(1,min(i+j,8)+1):
  phi=H(n,z)
  first=E((s.diff(gi,x)*gj.subs(x,y)+s.diff(gj,x)*gi.subs(x,y))*s.diff(phi,z))
  second=E(r*(gi*gj.subs(x,y)+gj*gi.subs(x,y))*s.diff(phi,z,2))
  integrated=s.integrate(first+second,(r,0,1),(t,0,1))
  # All i,j,n are positive Hermite ranks, so centering terms vanish.
  target=E(phi*H(i,z)*H(j,z))
  assert s.simplify(integrated-target)==0,(i,j,n,integrated,target)
  checks+=1
  if target!=0:
   rows.append(dict(i=i,j=j,n=n,first_kernel=str(first),second_kernel=str(second),target=str(target),integral=str(integrated)))
report=dict(status='PASS',assertions=checks,scope='Exact bilinear Gaussian coefficient identity; no finite positive original-VALUE implementation is asserted.',nonzero_rows=rows)
(out/'bilinear_identity_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='nonzero_rows'},indent=2))
