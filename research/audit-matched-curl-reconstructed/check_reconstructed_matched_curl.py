# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
#!/usr/bin/env python3
"""New checker. Not claimed byte-identical to the former 1471-test file."""
from pathlib import Path
import json,hashlib
import mpmath as m
m.mp.dps=70
P=Path(__file__).resolve().parent;c=d=m.mpf('.1');tests=[]
def check(name,ok):
 tests.append({'name':name,'pass':bool(ok)})
 if not ok:raise AssertionError(name)
def nodes(K):
 lo=m.mpf(1)/(4*K*K);hi=m.mpf(1)/4
 xs=[(lo+hi)/2+(hi-lo)/2*m.cos(m.pi*j/(K-1)) for j in range(K)]
 ds=[m.fprod([-xs[j]/(xs[i]-xs[j]) for j in range(K) if i!=j]) for i in range(K)]
 return [m.sqrt(x) for x in xs],ds
class Source:
 def __init__(self,A,sigma):self.A=A;self.sigma=sigma;self.e=A**m.mpf('3.9');self.nu=self.e/A;self.k=1/(sigma*self.e);self.calls=0
 def Q(self,x):self.calls+=1;return x/2+d*m.sin(self.k*x)/self.k
 def g0(self,x):return self.Q(x)-x/2
 def gt(self,x):return self.Q(x)-(m.mpf('.5')-self.nu)*x
 def f(self,x):return self.A*self.gt(x+self.A*self.g0(x))
 def I(self,p,z,K):
  bs,ds=nodes(K)
  return sum(D*(self.f(m.sqrt(1-b*b)*z+b*p)-self.f(m.sqrt(1-b*b)*z-b*p))/(2*b) for b,D in zip(bs,ds))
 def C(self,W,u):
  g0W=self.g0(W);q=W+self.A*g0W;gtq=self.gt(q)
  D0=self.g0(W+u)-g0W;gtshift=self.gt(q+u);Dt=gtshift-gtq
  return self.A*(self.gt(q+u+self.A*D0)-gtshift)-self.A*self.A*(self.g0(W+Dt)-g0W)
def R(A,nu,th,w):
 D=m.sin(th+w)-m.sin(th);base=th+A*c*m.sin(th)
 return nu*c*D+2*d/A*m.cos(th+w+A*c*(m.sin(th+w)+m.sin(th))/2)*m.sin(A*c*D/2)-c*(m.sin(th+nu*w+d*(m.sin(base+w)-m.sin(base)))-m.sin(th))
for A in map(m.mpf,['.2','.05','.01']):
 for power in (0,1,2):
  S=Source(A,A**power);M=2*S.nu*c+4*c*d
  for K in (2,4,8):
   bs,ds=nodes(K);Sk=sum(abs(D)/b for b,D in zip(bs,ds))
   check('filter_constant',abs(sum(ds)-1)<m.mpf('1e-60'))
   check('filter_inverse_bound',Sk<=8*K)
   for W,p in [(m.mpf('.17'),m.mpf('.4')),(m.mpf('-.8'),m.mpf('1.3'))]:
    S.calls=0;I=S.I(p,W,K);check('literal_response_cost',S.calls==4*K)
    bound=Sk*S.sigma*S.e*A*(d+S.e*c)
    check('intact_mark_bound',abs(I-S.e*p)<=bound)
    S.calls=0;value=S.C(W,S.sigma*I)/S.sigma;check('literal_surrogate_cost',S.calls==6)
    actual=R(A,S.nu,S.k*W,I/S.e)
    check('literal_normalized_identity',abs(value/(A*A*S.e)-actual)<m.mpf('1e-35'))
    check('uniform_matched_transfer',abs(actual-R(A,S.nu,S.k*W,p))<=M*bound/S.e)
  for theta in map(m.mpf,['-.7','.1','1.2']):
   for w in map(m.mpf,['-3','-.2','.8','4']):
    check('strong_incoming_first',abs(m.diff(lambda x:R(A,S.nu,theta,x),w))<=M)
    check('zero_jet',abs(m.diff(lambda x:R(A,S.nu,theta,x),0))<m.mpf('1e-55'))
leading=c*d**3/8*(m.exp(-m.mpf('.5'))-m.exp(-2));rem=c*d**5/12*m.sqrt(2/m.pi)
check('strict_analytic_sign',leading-rem>m.mpf('.00000582'))
Fn=lambda p:c*(d*m.sin(p)/2-m.besselj(1,2*d*m.sin(p/2))*m.cos(p/2))
value=2*m.quad(lambda p:p*Fn(p)*m.exp(-p*p/2)/m.sqrt(2*m.pi),[0,1,2,4,8,16,m.inf])
check('quadrature_inside_analytic_interval',abs(value-leading)<=rem)
out={'provenance':'New reconstruction; not an exact copy of prior bytes.','all_pass':all(x['pass'] for x in tests),'count':len(tests),'tests':tests,'analytic_lower':str(leading-rem),'numerical_limit':str(value),'historical_report_hash':'850ca174f7d44ab181ffa6718f6ac4ecbe4dc4733c6288962faac04cd300d99d','historical_C_source_hash':'238acc863fed07c0d4e08723cff51995609573f65a9bca056887fc7639b5dd94'}
(P/'reconstructed_matched_curl_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('all_pass','count','analytic_lower','numerical_limit')},indent=2))
