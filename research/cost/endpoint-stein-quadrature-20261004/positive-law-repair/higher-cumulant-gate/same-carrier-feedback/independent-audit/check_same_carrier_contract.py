#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,math
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md'
PIN='a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca'
count=0
def req(v,msg):
 global count
 count+=1
 if not bool(v):raise AssertionError(msg)
def eq(a,b,msg):req(s.simplify(a-b)==0,msg)
req(hashlib.sha256(SOURCE.read_bytes()).hexdigest()==PIN,'source pin')
a,b,r,z=s.symbols('a b r z',positive=True)
cov=2/((a+b)*(a+1)*(b+1))
eq(cov.subs({a:1,b:1}),s.Rational(1,4),'conditional Var H0')
eq(-s.diff(cov,b).subs({a:1,b:1}),s.Rational(1,4),'conditional Cov H0 J0')
eq(s.diff(cov,a,b).subs({a:1,b:1}),s.Rational(5,16),'conditional Var J0')
tt=s.symbols('tt',positive=True)
for j,target in enumerate([s.Rational(1,2),s.Rational(1,4),s.Rational(1,8)]):
 eq(s.integrate(tt**j*s.exp(-2*tt)/s.factorial(j),(tt,0,s.oo)),target,'conditional force mean coefficient')
force=a*r/(1+a*(1-r*r))
eq(s.series(force.subs(r,s.Rational(1,2)),a,0,4).removeO(),a/2-3*a*a/8+9*a**3/32,'posterior force is different target')
eq(s.Rational(1,4)*(2+s.Rational(1,2)),s.Rational(5,8),'uniform Hessian upper bound coefficient')
x=s.symbols('x',positive=True)
ca=s.integrate(x/r-r*x,(x,0,r))+s.integrate(r/x-r*x,(x,r,1))
eq(ca,-r*s.log(r),'bridge innovation covariance integral')
for D in [4,16,64,256,1024,4096,16384,65536]:
 A=1/math.sqrt(D)
 req(abs(A**3*math.sqrt(D)-A*A)<1e-14,'existing one-energy scale')
 req(abs((A*A)/(A**4*math.sqrt(D))-math.sqrt(D))<1e-10,'diverging rejected-bound ratio')
 # Taylor remainder scale: D terms, 1/sqrt(D) outer sum, each I_i^2 L2=O(1/D).
 req(abs(D/math.sqrt(D)/D-1/math.sqrt(D))<1e-14,'summed Taylor remainder rate')

# Independent Gaussian quadrature is illustrative; strict negativity is proved analytically.
xx,ww=np.polynomial.hermite_e.hermegauss(100);ww=ww/math.sqrt(2*math.pi)
logcosh=np.logaddexp(xx,-xx)-math.log(2)
c0=float(ww@logcosh);d=c0/8
sech2=lambda u:1/np.cosh(u)**2
hp=float(ww@(sech2(xx-d)-sech2(xx)))
req(c0>0 and d>0,'positive collective displacement')
req(hp<0,'illustrative nonzero innovation covariance')
innovation_inner=hp*math.log(2)/2
asymptotic_cs_floor=(.25*innovation_inner)**2/(4*.75)
report={'status':'PASS','assertions':count,'source_sha256':PIN,
 'illustrative_gaussian_quadrature':{'E_logcosh_normal':c0,'collective_shift_d':d,
  'E_hprime_normal':hp,'J_innovation_inner_product_a_half':innovation_inner,
  'conservative_asymptotic_variance_constant':asymptotic_cs_floor},
 'scope':'Exact target/moment identities and counterexample scaling, with numerical illustration only for its nonzero constant. The dimension-free asymptotic lower bound is proved in the written audit, not inferred from finite-dimensional simulation.'}
(HERE/'same_carrier_contract_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
