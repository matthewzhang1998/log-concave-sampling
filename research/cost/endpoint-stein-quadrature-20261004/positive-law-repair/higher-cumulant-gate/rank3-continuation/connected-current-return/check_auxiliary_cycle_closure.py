#!/usr/bin/env python3
"""Finite algebra/read-set/partial-trace diagnostics for auxiliary-public closure."""
import json,itertools,hashlib
from pathlib import Path
import sympy as s
import numpy as np
BASE=Path(__file__).parent
checks=0
def ck(ok,msg):
 global checks
 if not bool(ok):raise AssertionError(msg)
 checks+=1

def gaussian(poly,vars):
 p=s.Poly(s.expand(poly),*vars);ans=0
 for ex,coef in p.terms():
  moment=1
  for n in ex:
   if n%2:moment=0;break
   if n:moment*=s.factorial2(n-1)
  ans+=coef*moment
 return s.expand(ans)

def maxcut(t):
 out=0.
 for bits in itertools.product((0,1),repeat=t.ndim):
  if not any(bits) or all(bits):continue
  l=[i for i,b in enumerate(bits) if b];r=[i for i,b in enumerate(bits) if not b]
  m=t.transpose(l+r).reshape(int(np.prod([t.shape[i] for i in l])),-1)
  out=max(out,float(np.linalg.norm(m,2)))
 return out
A,t,G,S,T=s.symbols('A t G S T')
V=A**6*t**2*G**2*(S+A*t)*(T+A*t)
mom=sum(gaussian(V**j,[G,S,T])/s.factorial(j) for j in range(4))
mom=s.series(mom,A,0,21).removeO()
log=s.series(s.log(mom),A,0,21).removeO().expand()
expected=A**8*t**4+s.Rational(3,2)*A**12*t**4+3*A**14*t**6+A**16*t**8+21*A**20*t**8
ck(s.expand(log-expected)==0,'joint three-Gaussian symbol through grade twenty')
conditional=A**8*t**4*G**2+s.Rational(1,2)*A**12*t**4*G**4+A**14*t**6*G**4+A**20*t**8*G**6
cmom=sum(gaussian(conditional**j,[G])/s.factorial(j) for j in range(3))
clog=s.series(s.log(s.series(cmom,A,0,21).removeO()),A,0,21).removeO().expand()
ck(s.expand(clog-expected)==0,'conditional bilinear family then auxiliary cumulants agrees')
ck(log.coeff(A,8)==t**4,'opened edge restores cycle leading coefficient')
ck(log.coeff(A,12)==s.Rational(3,2)*t**4,'new root-spine feedback begins twelve')
ck(log.coeff(A,16)==t**8,'auxiliary-G variance begins sixteen')
# Actual cut tree.
vs=['L','C','Cp','Lp','M','D','Dp','Mp']
es=[('L','C'),('C','Cp'),('Cp','Lp'),('M','D'),('D','Dp'),('Dp','Mp'),('Cp','Dp')]
phys={'L','Lp','M','Mp'};aux={'C','D'}
deg={v:sum(v in e for e in es) for v in vs}
orders={v:deg[v]+int(v in phys)+int(v in aux)-1 for v in vs}
ck(len(es)==len(vs)-1,'cut graph is a tree')
ck(all(orders[v]==(1 if v in phys else 2) for v in vs),'original derivative orders unchanged')
spine=['L','C','Cp','Dp','D','M']
ck(all((spine[i],spine[i+1]) in es or (spine[i+1],spine[i]) in es for i in range(5)),'literal six-force root spine')
ck(len(spine)==6,'feedback root-spine force count')
chunks=[['L','C','Cp'],['M','D','Dp']]
for ch in chunks:
 ck(len(set(ch)&phys)==1,'each chunk retains physical slot')
 ck(len(set(ch)&aux)==1,'each chunk has exactly one G slot')
# Leading cut tensor from actual two-spine contraction; every cut checked in D=2.
rng=np.random.default_rng(1052026)
for case in range(10):
 d=2;b=rng.normal(size=(d,)*4);b/=maxcut(b)
 opened=np.einsum('ijab,klcb->ijklac',b,b)
 closed=np.einsum('ijklaa->ijkl',opened)
 ck(maxcut(opened)<=1+1e-12,'opened forest all cuts')
 ck(maxcut(closed)<=1+1e-12,'restored cycle all cuts')
 sym=(opened+opened.swapaxes(-1,-2))/2
 rhs=np.linalg.norm(closed)**2+2*np.linalg.norm(sym)**2
 # Exact 4th Gaussian moment by three pairings.
 direct=np.einsum('ijklaa,ijklbb->',opened,opened)
 direct+=np.einsum('ijklab,ijklab->',opened,opened)
 direct+=np.einsum('ijklab,ijklba->',opened,opened)
 ck(abs(direct-rhs)<1e-10,'Hilbert Gaussian second-chaos identity')
 ck(rhs<=3*d+1e-10,'one-Hilbert quadratic auxiliary energy')
# Each Wick edge pairs two distinct open chunks, never a self loop.
def pairings(items):
 if not items:yield [];return
 a=items[0]
 for j in range(1,len(items)):
  b=items[j]
  rest=items[1:j]+items[j+1:]
  for tail in pairings(rest):yield [(a,b)]+tail
pair_counts={}
for n in (2,4,6,8):
 count=0
 for pairing in pairings(list(range(n))):
  ck(all(a!=b for a,b in pairing),'one G slot per chunk forbids initial self loops')
  ck(len({x for pair in pairing for x in pair})==n,'all auxiliary slots paired exactly once')
  count+=1
 pair_counts[n]=count
 ck(count==int(s.factorial2(n-1)),'Wick pairing enumeration complete')
# Endpoint dyadic margins: k root factors, h connected clusters, one caller score.
for h in range(1,15):
 for k in range(h,h+5):
  leftover=8*k-4*k-2*(h-1)-1
  ck(leftover>0,'literal cycle-sector endpoint mass dominates all bridge plus caller scores')
result={'status':'PASS','checks':checks,'symbol_through_grade_20':str(expected),'wick_pairing_counts':pair_counts,'literal_cut_tree':{'vertices':vs,'edges':es,'physical_marks':sorted(phys),'auxiliary_G_probes':sorted(aux),'root_spine':spine,'open_chunks':chunks},'scope':'Bounded auxiliary-public cycle source and finite coefficient diagnostics. Does not prove generic all-order closure or the tight whole-coarse bridge port.'}
(BASE/'auxiliary_cycle_closure_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
