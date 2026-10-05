#!/usr/bin/env python3
"""Independent finite diagnostics, not an analytic/native realizability proof."""
from pathlib import Path
import hashlib, json
from fractions import Fraction as Q
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
records={}

def check(value,name):
    assert value, name
    records[name]=True

# Pin validation is independent of the producer's checker, which does not validate pins.
pins=[]
for x in json.loads((ROOT/'INPUT-PINS.json').read_text()):
    b=Path(x['path']).read_bytes()
    good=len(b)==x['bytes'] and hashlib.sha256(b).hexdigest()==x['sha256']
    check(good,'pin:'+x['path'])
    pins.append(dict(path=x['path'],sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),match=good))

# Closed-form Gaussian integral, rather than the submitted Wick/partition census.
l,t,a,b=S.symbols('lambda theta a b')
z=l*t**2
log_mgf=-S.log(1-z*z)/2+(z*a*b*t*t+z*z*(a*a+b*b)*t*t/2)/(1-z*z)
connected=S.series(log_mgf,l,0,5).removeO().expand()
expected=l*a*b*t**4+l*l*(t**4+(a*a+b*b)*t**6)/2+l**3*a*b*t**8+l**4*(t**8/4+(a*a+b*b)*t**10/2)
check(S.expand(connected-expected)==0,'rank4_closed_form_through_four')
raw=S.series(S.exp(log_mgf),l,0,3).removeO().expand().coeff(l,2)
check(S.expand(raw-(t**4+(a*a+b*b)*t**6+a*a*b*b*t**8)/2)==0,'raw_second_includes_rank_eight')

# The public-bank derivative is not pre-replica differentiation.
u,B,d=S.symbols('s B d',positive=True)
integrand=d*d*S.exp(-(1-u))*S.sin(S.sqrt(u)*B)**2
actual=S.diff(integrand,B)
wrong=2*d*d*S.exp(-(1-u))*S.sin(S.sqrt(u)*B)*S.cos(S.sqrt(u)*B)
check(S.simplify(actual-S.sqrt(u)*wrong)==0,'retained_derivative_has_sqrt_s')
check(S.simplify((wrong-2*actual).subs({u:S.Rational(1,4),B:S.pi/2}))==0,'omission_factor_two_example')
actual_second=S.integrate(S.diff(integrand,B,2).subs(B,0),(u,0,1))
wrong_second=S.integrate(2*d*d*S.exp(-(1-u)),(u,0,1))
check(S.simplify(actual_second-2*d*d/S.E)==0,'integrated_retained_second_derivative')
check(S.simplify(wrong_second-actual_second)!=0,'integrated_derivative_commutation_shortcut_fails')

# Gaussian bank IBP creates a genuine consumed-record covariance-image relation.
x=S.symbols('x',real=True)
rho=S.exp(-x*x/2)/S.sqrt(2*S.pi)
q=S.exp(-x*x/2)
p=(9*x*x-3)*S.exp(-x*x)
check(S.simplify(S.diff(rho*q*q,x,2)-rho*p)==0,'consumed_E_rank4_to_covariance_rank2_relation')
# q and p are bounded (polynomial times a decaying Gaussian); p is not a
# conditional representative: for frozen E a Fourier mode would require
# -p*t^2=q^2*t^4 for every t, forcing q=p=0.
pp,qq=S.symbols('p q2')
check(S.solve([-pp-qq,-4*pp-16*qq],(pp,qq))=={pp:0,qq:0},'conditional_nonmembership_requires_all_retained_tests')

# An endpoint reading a consumed bank still prevents coefficient averaging.
# F(B)=c+d*cos B, W(B)=B. Gaussian Fourier transform of F differs from
# E[F]*E[exp(i*t*B)] by d*exp(-(1+t^2)/2)*(cosh(t)-1).
c=S.symbols('c')
weighted=c*S.exp(-t*t/2)+d*(S.exp(-(t-1)**2/2)+S.exp(-(t+1)**2/2))/2
factored=(c+d*S.exp(-S.Rational(1,2)))*S.exp(-t*t/2)
check(S.simplify((weighted-factored).subs(t,1))!=0,'consuming_B_does_not_remove_endpoint_dependency')

# A mixed-current witness within one scalar owned-carrier coisometry.
C=S.Matrix([[1,0,S.Rational(1,4),S.Rational(1,4)],
            [0,1,S.Rational(1,4),S.Rational(1,4)],
            [S.Rational(1,4),S.Rational(1,4),1,0],
            [S.Rational(1,4),S.Rational(1,4),0,1]])
check(all(x>=0 for x in C.eigenvals()),'scalar_common_carrier_cross_covariance_psd')
cp=S.Rational(1,4); shift=S.Rational(1,2)
mixed=cp*cp+cp*cp+t*t*4*shift*shift*cp
check(mixed==S.Rational(1,8)+t*t/4,'realizable_scalar_carrier_mixed_term_nonzero')

# Exact scalar covariance path coefficients from conditional characteristic law.
s,q0,r=S.symbols('s q r')
variance_multiplier=S.series(S.exp((s*q0+s*s*r/2)*t*t),s,0,4).removeO().expand()
check(S.expand(variance_multiplier.coeff(s,2)-(r*t*t+q0*q0*t**4)/2)==0,'covariance_acceleration_and_square')
check(S.expand(variance_multiplier.coeff(s,3)-(q0*r*t**4/2+q0**3*t**6/6))==0,'covariance_third_coefficient')

# Independent rational ledger with worst width bounds, no presumed exact jet.
beta,gamma,P=Q(1,4),Q(1,8),Q(11,2)
grades=[]
for name,N,K in [('T2',2,0),('T3',3,1),('T4',4,2)]:
    nominal=Q(4); surplus=nominal-beta*N-gamma*K
    root=beta+surplus
    grades.append(dict(name=name,N=N,K=K,G=str(nominal-gamma*K),Psi=str(surplus-2*gamma),root=str(root),floor=str(2*root)))
    check(root>0 and 2*root>P,'strict_native_floor:'+name)
    check(beta*N+gamma*K+surplus==nominal,'normalization:'+name)
check(min(Q(g['floor']) for g in grades)==6,'minimum_summed_root_floor_six')
check(min(Q(8)-gamma*(g['K']+h['K']+2) for g in grades for h in grades)==Q(29,4),'minimum_pair_bridge_grade_29_over_4')

snapshot=[]
for name in ['CURRENT-IMAGE-AND-RESIDUAL.md','CURRENT-DICTIONARY.json','check_dictionary.py','INPUT-PINS.json']:
    bb=(ROOT/name).read_bytes()
    snapshot.append(dict(file=name,sha256=hashlib.sha256(bb).hexdigest(),bytes=len(bb)))
result=dict(status='PASS',scope='Independent exact diagnostics only; no native quadratic jet or finite implementation certificate is proved by these tests.',checks=records,pins=pins,grades=grades,artifact_snapshot=snapshot)
(OUT/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
