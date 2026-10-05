#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
import sympy as s

# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE=Path(__file__).resolve().parent
PINS={str(HERE.parent/'ALL-RANK-CENTERED-STEIN-CURRENT-RECURRENCE.md'):'d09db0fa283512de36bc910f1d3f828257121699832c0897943d9a252f59dd0e',
str(HERE.parent/'ONE-MARKED-THIRD-CUMULANT-REPLACEMENT.md'):'5581eb0e592d2997504ac0c69c9f577e46306b6ca85718feb6e0a4a852524a19'}
count=0
def req(v,msg):
 global count
 count+=1
 if not bool(v):raise AssertionError(msg)
def eq(a,b,msg):req(s.expand(a-b)==0,msg)
for p,h in PINS.items():req(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,'source pin')
x,z,t,ss=s.symbols('x z t ss',nonzero=True)
def moment(p,vs):
 ans=0
 for powers,c in s.Poly(s.expand(p),*vs).terms():
  ans+=c*s.prod(0 if n%2 else s.factorial2(n-1) if n else 1 for n in powers)
 return s.expand(ans)
def riesz(p):
 eq(moment(p,[x]),0,'centered Riesz input')
 inv=0
 for (n,),c in s.Poly(s.expand(p),x).terms():
  for k in range(n//2+1):
   degree=n-2*k
   if degree:inv+=c*s.factorial(n)*s.hermite_prob(degree,x)/(2**k*s.factorial(k)*s.factorial(degree)*degree)
 return s.expand(s.diff(inv,x))
a=s.Rational(2,3);b=s.Rational(1,5)
X=a*x+b*(x*x-1)
T=[X]
for j in range(1,9):
 T.append(s.expand(riesz(T[j-1]-moment(T[j-1],[x]))*s.diff(X,x)))
 n=j+1
 cumulant=2**(n-1)*s.factorial(n-1)*b**n+s.factorial(n)*a*a*(2*b)**(n-2)/2
 eq(moment(T[j],[x]),cumulant/s.factorial(j),'Stein constant is connected cumulant divided by factorial')
Y=t*X+ss*z
variance=moment(X*X,[x])
for m in range(2,7):
 for degree in range(1,8):
  val=moment(Y**degree,[x,z])
  lhs=s.diff(val,t)-t*variance/s.Symbol('ss',nonzero=True)*s.diff(val,ss)
  rhs=0
  for j in range(2,m):
   if degree>=j+1:
    rhs+=t**j*moment(T[j],[x])*s.factorial(degree)/s.factorial(degree-j-1)*moment(Y**(degree-j-1),[x,z])
  if degree>=m+1:
   rhs+=t**m*s.factorial(degree)/s.factorial(degree-m-1)*moment(T[m]*Y**(degree-m-1),[x,z])
  eq(lhs,rhs,'exact same-endpoint finite current expansion')
for m in range(2,9):eq(s.integrate(t**m,(t,0,1)),s.Rational(1,m+1),'dynamic length coefficient')
report={'status':'PASS','assertions':count,'pins':{source_label(p):v for p,v in PINS.items()},
 'scope':'Exact scalar connected constants through rank nine and same-endpoint current identities through rank seven. Polynomial source tests algebra, not global Lipschitz hypotheses; one-marked bound is proved in the written audit.'}
(HERE/'all_rank_stein_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
