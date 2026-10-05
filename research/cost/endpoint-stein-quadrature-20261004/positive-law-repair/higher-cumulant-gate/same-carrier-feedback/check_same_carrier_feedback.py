#!/usr/bin/env python3
"""Independent-reproducibility diagnostics, not a substitute for the C2 proofs."""
import json, math
from pathlib import Path
import sympy as s
from scipy.integrate import quad

a=s.symbols('a', positive=True)
r=s.symbols('r', positive=True)
# Gaussian stationary linear filter covariance: X,H,J,K satisfy triangular OU ODEs.
# Solve A Sigma + Sigma A^T + Q = 0 exactly (Brownian enters X only).
M=s.Matrix([[-1,0,0,0],[1,-1,0,0],[0,1,-1,0],[0,0,1,-1]])
vs=s.symbols('s0:10')
Sigma=s.zeros(4); ix=0
for i in range(4):
    for j in range(i,4):
        Sigma[i,j]=Sigma[j,i]=vs[ix]; ix+=1
Q=s.diag(2,0,0,0)
sol=s.solve(list(M*Sigma+Sigma*M.T+Q),vs,dict=True)[0]
Sigma=Sigma.subs(sol)
Scond=Sigma[1:,1:]-Sigma[1:,0:1]*Sigma[0:1,1:]
f2=s.Matrix([a,-a*a,0]); f3=s.Matrix([a,-a*a,a**3])
varf2=s.expand((f2.T*Scond*f2)[0])
varf3=s.expand((f3.T*Scond*f3)[0])
assert varf2 == a*a/4-a**3/2+5*a**4/16
assert s.expand(varf3-varf2).coeff(a,3)==0
mean3=(s.Matrix([[a,-a*a,a**3]])*Sigma[1:,0:1])[0]
assert mean3==a/2-a*a/4+a**3/8
posterior=s.series(a*s.Rational(1,2)/(1+a*s.Rational(3,4)),a,0,4).removeO()
assert s.expand(mean3-posterior).coeff(a,2)==s.Rational(1,8)
bridge=s.integrate(r/a-r*a,(r,0,a))+s.integrate(a/r-r*a,(r,a,1))
assert s.simplify(bridge+a*s.log(a))==0
normal=lambda x:math.exp(-x*x/2)/math.sqrt(2*math.pi)
logcosh=lambda x: abs(x)+math.log1p(math.exp(-2*abs(x)))-math.log(2)
sech2=lambda x: 4*math.exp(-2*abs(x))/(1+math.exp(-2*abs(x)))**2
c0=quad(lambda x:logcosh(x)*normal(x),-12,12,epsabs=1e-13)[0]
d=c0/8
hp=quad(lambda x:(sech2(x-d)-sech2(x))*normal(x),-12,12,epsabs=1e-13)[0]
assert d>0 and hp<0
lower=hp*hp*(math.pi*math.pi-4)/16
innovation=s.Matrix([a*a/2,-a*a,0])
innovation_var=s.expand((innovation.T*Scond*innovation)[0])
assert innovation_var==a**4/8
cross_hg=s.expand((s.Matrix([[a,0,0]])*Scond*innovation)[0])
assert cross_hg==-a**3/8
Bj=a/2-a*a/4; K=-a**3/4
localized=s.expand(Bj**2+K)
assert s.expand(varf2-localized)==a**4/4
time=s.symbols('t', positive=True)
assert s.integrate(time**2*s.exp(-2*time)/2,(time,0,s.oo))==s.Rational(1,8)
q=s.symbols('q', positive=True)
assert s.integrate(q**2/s.sqrt(1-q**2),(q,0,1))==s.pi/4
Cstar=s.simplify((s.Rational(1,2)+s.pi/4)/2+s.pi/8)
assert Cstar==(1+s.pi)/4
# Noncommuting symmetric original Hessian words check the orientation debt.
T=s.Matrix([[2,1],[1,1]])/10; J=s.Matrix([[1,0],[0,3]])/10
BJ=T*(s.eye(2)-J)
skew=BJ-BJ.T
assert skew==J*T-T*J
assert BJ*BJ.T==BJ*BJ-BJ*skew
assert BJ*BJ.T!=BJ*BJ
res={
 'status':'PASS diagnostic identities only; independent prose audit required',
 'stationary_covariance':str(Sigma),
 'resolvent_paraproduct_constant':str(Cstar),
 'noncommuting_BBstar_minus_Bsquare':str(BJ*BJ.T-BJ*BJ),
 'conditional_HJK_covariance':str(Scond),
 'm3_coefficient':str(mean3),
 'posterior_r_half_coefficient':str(posterior),
 'F2_conditional_variance':str(varf2),
 'innovation_conditional_variance':str(innovation_var),
 'mixed_Hg_innovation_covariance':str(cross_hg),
 'localized_full_covariance_through_order3':str(localized),
 'localization_omitted_order4':str(s.expand(varf2-localized)),
 'F3_minus_F2_conditional_variance':str(s.expand(varf3-varf2)),
 'residual_bridge_kernel':str(s.simplify(bridge)),
 'mean_logcosh':c0,'d':d,'Ehprime':hp,
 'strict_common_root_variance_gap_lower_bound':lower,
 'claimed_remainder_counterfamily':[
   {'D':D,'A':D**-.5,'A2':D**-1,'A4_sqrtD':D**-1.5,'scale_ratio':math.sqrt(D)}
   for D in [10**2,10**4,10**6,10**8]],
 'limitations':['No Monte Carlo validation is used as proof.','No finite mean/covariance producer or endpoint law pass is asserted.']
}
Path(__file__).with_name('same_carrier_feedback_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
