#!/usr/bin/env python3
"""Independent-of-production finite checks; diagnostics, not native execution."""
import json, math, hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(710051)
checks={}

def multiindices(d,m):
    if d==1:
        for k in range(m+1): yield (k,)
    else:
        for k in range(m+1):
            for rest in multiindices(d-1,m-k): yield (k,)+rest

def hermites(x,m):
    out=[1.0]
    if m: out.append(float(x))
    for k in range(1,m): out.append(x*out[-1]-k*out[-2])
    return [v/math.sqrt(math.factorial(k)) for k,v in enumerate(out)]

def kernels(x,m):
    hs=[hermites(z,m) for z in x]; k=0.; grad=0.
    for a in multiindices(len(x),m):
        val=math.prod(hs[i][ai] for i,ai in enumerate(a)); k+=val*val
        for i,ai in enumerate(a):
            if ai:
                gv=math.sqrt(ai)*math.prod(hs[j][aj-(i==j)] for j,aj in enumerate(a))
                grad+=gv*gv
    return k,grad

rows=[]
for d in [1,2,3,5,8]:
  for m in range(1,6):
    for _ in range(8):
      x=rng.normal(size=d)*rng.uniform(0,3); R=np.linalg.norm(x)
      k,g=kernels(x,m)
      kb=4*(1+d+R*R)**m
      gb=4*(d+m)*(1+d+R*R)**(m-1)
      assert k <= kb*(1+1e-12) and g<=gb*(1+1e-12)
      rows.append((k/kb,g/gb))
checks['hermite_kernel']={'cases':len(rows),'max_ratio_kernel':max(x[0] for x in rows),'max_ratio_derivative':max(x[1] for x in rows)}

matrix_cases=[]
for d,m in [(1,1),(2,3),(4,3),(6,2)]:
  inds=list(multiindices(d,m-1)); mats=rng.normal(size=(len(inds),d,d))
  gram=sum(C.T@C for C in mats)
  for _ in range(10):
    x=rng.normal(size=d); hs=[hermites(z,m-1) for z in x]
    phis=np.array([math.prod(hs[j][a[j]] for j in range(d)) for a in inds])
    J=sum(c*C for c,C in zip(phis,mats))
    actual=np.linalg.norm(J,ord=2)**2
    bound=np.dot(phis,phis)*np.linalg.norm(gram,ord=2)
    assert actual<=bound*(1+1e-12)
    matrix_cases.append(actual/bound)
checks['matrix_parseval_cauchy']={'cases':len(matrix_cases),'max_ratio':max(matrix_cases)}


# Actual marginal-law-zero rotation example, preserving all shared-H covariance.
rotation=[]
for J in [2,3,8,30]:
  for r in [0.001,0.02,0.1,0.3]:
    c,s=math.cos(r),math.sin(r)
    R=np.array([[c,-s],[s,c]])
    A=np.eye(2)+J*(R-np.eye(2))
    cov=A@A.T
    sd=math.sqrt(float(cov[0,0]))
    exact=math.sqrt(2)*abs(sd-1)
    S=J*np.linalg.norm(R-np.eye(2),ord=2)
    bound=2*S*S*math.sqrt(2)
    assert exact<=bound+1e-14
    assert np.max(np.abs(cov-(1+2*J*(J-1)*(1-c))*np.eye(2)))<1e-11
    rotation.append({'J':J,'r':r,'exact_w2':exact,'bound':bound})
checks['shared_rotation']={'cases':len(rotation),'max_ratio':max(x['exact_w2']/x['bound'] for x in rotation)}

# Coisometry and full Gaussian source-zero geometry.
geometry=[]
for d in [1,3,7]:
  k=4*d; L=rng.normal(size=(d,k)); Q,_=np.linalg.qr(L.T); L=Q[:,:d].T*math.sqrt(.7)
  Pi=np.eye(k)-L.T@L/.7
  total=L.T@L/.7+Pi@Pi.T
  geometry.append(max(np.max(abs(total-np.eye(k))),np.max(abs(L@Pi))))
assert max(geometry)<1e-12
checks['common_carrier_geometry']={'max_error':max(geometry)}

# Rank-nine share bookkeeping from scaling a^9 at the current level.
R=9
assert Fraction(17,2)+8*Fraction(1,16)==R
assert 1-2*R==-17 and 1-R==-8
assert -Fraction(1,2)*(1-2*R)==Fraction(17,2)
assert -Fraction(1,2)*(1-R)==4
checks['variance_partition_exponents']={'root':'17/2','nonroot':'1/16','native_share_power':'-17','native_J_power':'17/2','energy_J_power':4,'null_join_J_power':8,'native_alpha_power':17,'null_join_alpha_power':18}
assert Fraction(17-10,17)==Fraction(7,17)
assert Fraction(18-10,16)==Fraction(1,2)
checks['shrinking_eta_sufficient_thresholds']={'native':'q_eta < 7/17','null_bank':'q_eta < 1/2','requires':'otherwise public-log complete old-history census and all actual guards'}


# Verify chosen finite-D radii/guard terminate at small alpha in sample cases.
guards=[]
for d,m in [(1,1),(2,3),(10,5),(100,8)]:
  found=None
  for n in range(4,300):
    a=2.0**(-n); K=a**-.5
    ell=1
    while True:
      rad=math.sqrt(d)+ell
      logtail2=math.log(2*K)+(m/2)*math.log(1+d+rad*rad)-ell*ell/4
      tail=math.exp(-5*ell*ell/24)+math.exp(logtail2)
      if tail<=a**2.5: break
      ell+=1
    Cdim=2*(1+d+rad*rad)**((m-1)/2)
    L=max(1,K*Cdim); S=a**8.5
    if L*S<=1 and 2*L*a**7<=1 and 4*(1+L*S)*a<=1:
      found={'D':d,'m':m,'alpha':a,'dyadic_exponent':n,'ell':ell,'hybrid_guard':2*L*a**7,'bias_guard':4*(1+L*S)*a}
      break
  assert found is not None
  guards.append(found)
checks['finite_dimension_guard_examples']=guards
checks['scope']='Finite diagnostics support stated algebra and explicit guards. They do not execute or re-prove the imported native compiler, and do not establish a public-log dimensional bound.'
(ROOT/'carrier-replacement-checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps({'status':'PASS','hermite_cases':len(rows),'rotation_cases':len(rotation),'geometry_cases':len(geometry),'guard_examples':len(guards)},indent=2))
