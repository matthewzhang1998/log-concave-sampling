#!/usr/bin/env python3
"""Combined-join scalar/algebra/source-budget diagnostics, not a native compiler run."""
from fractions import Fraction as F
from pathlib import Path
import math,json
import numpy as np
from numpy.polynomial.legendre import leggauss
rng=np.random.default_rng(370102026)
checks=0

def ck(test,label):
 global checks
 checks+=1
 if not test: raise AssertionError(label)

def outer(A):
 eps=A**3; a=eps/64; T=math.log(64/eps)
 boundaries=[a]; p=a
 while p<T:
  p=min(2*p,T); boundaries.append(p)
 qcut=-math.log1p(-A**.8); ecut=-math.log1p(-A**1.4)
 boundaries=sorted(set(boundaries+[qcut,ecut]))
 x,b=leggauss(8+math.ceil(math.log(1/eps,4)))
 ts=[a/2,T+1]; ws=[-math.expm1(-a),math.exp(-T)]
 for left,right in zip(boundaries[:-1],boundaries[1:]):
  tau=(left+right)/2+(right-left)*x/2
  ts.extend(tau); ws.extend((right-left)*b*np.exp(-tau)/2)
 ts=np.array(ts);ws=np.array(ws);ws/=ws.sum()
 return ts,ws,ecut

maxrat={str(p):0 for p in [.5,1,1.5,7/6,2,2.5]}
for A in np.geomspace(1e-10,.1,31):
 w=A**.8;eta=A**1.4;h=A**1.7;r=math.sqrt(1-h*h);q=1-w
 taus,weights,cut=outer(A)
 c2=-np.expm1(-2*taus);sig2=w*(2-w)
 v=c2*sig2/(sig2+q*q*c2)
 law=taus>=cut;raw=~law
 ck(np.all(taus>0),'strictly positive OU times')
 ck(np.all(weights>0),'positive outer weights')
 ck(abs(weights.sum()-1)<1e-13,'outer mass one')
 ck(np.all(v[law]>=eta/2*(1-1e-10)),'minimum actual conditional variance')
 ck(np.all(v<=2*w),'maximum actual conditional variance')
 ck(np.max(np.abs(1/v-(q*q/sig2+1/c2))/(1/v))<1e-12,'exact reciprocal variance')
 for p in maxrat:
  pp=float(p);got=np.dot(weights[law],v[law]**(-pp))
  bound=w**(-pp)+(eta**(1-pp) if pp>1 else (math.log(1/eta) if pp==1 else 1))
  ratio=got/bound;maxrat[p]=max(maxrat[p],ratio)
  ck(ratio<=6,'finite weighted singular sum')
 near=np.dot(weights[raw],v[raw]**-.5)
 ck(near<=6*(eta/math.sqrt(w)+math.sqrt(eta)),'finite RAW endpoint sum')
 # The exact positive-sum envelope; individual L_t may be larger than A.
 Lt=A+A**1.5/np.sqrt(v[law])
 aggregate=np.dot(weights[law],Lt)
 ck(aggregate<=7*A,'aggregate full first')
 ck(np.max(A**1.5/(v[law]/8))<=16*A**.1*(1+1e-10),'native self-reserve worst power')
 ck(np.max(A/np.sqrt(v[law]/8))<=4*A**.3*(1+1e-10),'mean and cubic native radius power')
 ck(np.max(A/(v[law]/8)**(1/3))<=16**(1/3)*A**(8/15)*(1+1e-10),'mixed K native radius power')
 ck(np.max(A**3/v[law]**1.5)<=2**1.5*A**.9*(1+1e-10),'reference tensor normalized guard')
 for idx in rng.choice(len(v),size=min(50,len(v)),replace=False):
  shares=np.array([v[idx]/8]*4+[v[idx]/2])
  ck(abs(shares.sum()-v[idx])<=1e-12*v[idx],'five positive shares exact variance')
  ck(np.all(shares>0),'all shares strictly positive')
  # Each literal source-zero service carrier is included, then smoothing row.
  row=np.r_[np.sqrt(shares),h]
  ck(np.all(np.isfinite(row)),'finite zero-carrier row')
 # Source-zero terminal carrier identity, using stable c2.
 for cc2 in c2[::max(1,len(c2)//53)]:
  dn=math.sqrt(r*r*cc2+h*h)
  R=np.array([[r*math.sqrt(cc2),h],[-h,r*math.sqrt(cc2)]])/dn
  ck(np.linalg.norm(R@R.T-np.eye(2))<1e-12,'whole terminal orthogonal alignment')
  ck(abs((r*r*(1-cc2)+dn*dn)-1)<1e-12,'OU terminal variance')

alpha=F(4,5);beta=F(7,5);k=F(17,10)
leading={'terminal_commutator':2+k,'bridge_cancellation':3+3*alpha-k,'near_RAW':3+beta/2}
rest={'nested_prefix':3+alpha,'smooth_drift':1+2*k,'near_other':3+beta-alpha/2,'bridge_quadrature':4+alpha,'outer_rule':F(4),'intrinsic_mean':F(4),'final_own_mean':F(4)}
rows=[(5,F(3,2)),(5,F(1)),(6,F(7,6)),(6,F(2)),(7,F(5,2)),(7,F(3,2))]
for i,(power,p) in enumerate(rows):
 rest[f'bulk_{i}']=power-p*alpha
 rest[f'endpoint_{i}']=F(power) if p==1 else power+(1-p)*beta
ck(set(leading.values())=={F(37,10)},'leading 37/10 balance')
ck(min(rest.values())==F(19,5),'strict one-tenth remainder margin')
ck(1-beta/2==F(3,10),'mean radius exponent')
ck(F(3,2)-beta==F(1,10),'self-reserve exponent')
ck(1-beta/3==F(8,15),'K radius exponent')
ck(3-F(3,2)*beta==F(9,10),'reference skew exponent')
result={'status':'PASS','assertions':checks,'max_finite_sum_ratios':maxrat,'leading_exponents':{k:str(v) for k,v in leading.items()},'remainder_exponents':{k:str(v) for k,v in rest.items()},'scope':'Literal positive outer dyadic sums, variance shares, stable endpoint geometry, aggregate first envelope, terminal carrier rotations, and exact exponent/guard arithmetic. No native compiler or sampler execution.'}
Path(__file__).with_name('combined_join_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
