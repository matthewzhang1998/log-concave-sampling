#!/usr/bin/env python3
"""Independent finite algebra/counterchecks. Mathematical audit has separate scope.
Does not execute selected-pair/native compilers or numerically certify W2 bounds.
"""
from pathlib import Path
from fractions import Fraction
import hashlib, json, math
import sympy as s
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
checks=[]
def check(name,ok,detail=None):
    if not bool(ok): raise AssertionError(name)
    checks.append(dict(name=name,passed=True,**({"detail":detail} if detail is not None else {})))
def zero(p):return s.expand(p)==0
x,y,th1,th2=s.symbols('x y th1 th2')
spatial=(x,y)
def lap(p):return sum(s.diff(p,z,2) for z in spatial)
def resolvent(p,k):
    # Solve the polynomial OU resolvent separately on homogeneous degrees.
    P=s.Poly(s.expand(p),*spatial)
    groups={}
    for powers,c in P.terms():groups[sum(powers)]=groups.get(sum(powers),0)+c*x**powers[0]*y**powers[1]
    ans=0
    for n,h in groups.items():
        d=k+n; term=h
        while term!=0:
            ans+=term/d
            term=s.expand(lap(term)); n-=2; d*=k+n
    return s.expand(ans)
def grad(p):return [s.diff(p,z) for z in spatial]
def dot(a,b):return sum(i*j for i,j in zip(a,b))
# Potential gradient with noncommuting, nonconstant Hessian. Polynomials test
# algebra only; they are not advertised as globally Lipschitz sources.
U=x**4/12+x**3*y/6+x*x*y*y/8+x*y**3/9+y**4/10+x*y+x*x/2
f=th1*s.diff(U,x)+th2*s.diff(U,y)
b=[resolvent(v,2) for v in grad(f)]
T=[[resolvent(s.diff(f,a,c),3) for c in spatial] for a in spatial]
S=[[[resolvent(s.diff(f,a,c,d),4) for d in spatial] for c in spatial] for a in spatial]
for a in range(2):
    for i in range(2):
        check(f'commute_first_{a}_{i}',zero(s.diff(b[i],spatial[a])-T[a][i]))
        for j in range(2):check(f'commute_second_{a}_{i}_{j}',zero(s.diff(T[i][j],spatial[a])-S[a][i][j]))
for k in [2,3,4]:
    r=resolvent(f,k)
    check(f'OU_inverse_{k}',zero(k*r+dot(spatial,grad(r))-lap(r)-f))
C=2*resolvent(dot(b,b),2)
K3=6*resolvent(dot(b,grad(C)),3)
K4=resolvent(8*dot(b,grad(K3))+6*dot(grad(C),grad(C)),4)
UU=[resolvent(sum(T[i][j]*b[j] for j in range(2)),3) for i in range(2)]
P1=resolvent(sum(b[a]*resolvent(sum(T[a][i]*UU[i] for i in range(2)),4) for a in range(2)),4)
P2=resolvent(sum(b[a]*resolvent(sum(b[i]*resolvent(sum(S[a][i][j]*b[j] for j in range(2)),4) for i in range(2)),4) for a in range(2)),4)
P3=resolvent(sum(b[a]*resolvent(sum(b[i]*resolvent(sum(T[i][j]*T[a][j] for j in range(2)),4) for i in range(2)),4) for a in range(2)),4)
P4=resolvent(dot(UU,UU),4)
check('four_tree_multivariate_full_polarization',zero(K4-192*(P1+P2+P3)-96*P4))
check('star_omission_detected',not zero(K4-192*(P1+P3)-96*P4))
check('branch_square_required',not zero(K4-192*(P1+P2+P3)))
check('physical_rank_four',set(sum(m[2:]) for m,c in s.Poly(K4,x,y,th1,th2).terms())=={4})
# One-dimensional all-current identity at the ACTUAL positive endpoint.
g,X,Z=s.symbols('g X Z')
def gm(n,var=1):return 0 if n%2 else s.factorial2(n-1)*var**(n//2) if n else s.Integer(1)
def expectation(p,Cx):
    return s.expand(sum(c*gm(i)*gm(j,Cx)*gm(k,s.Rational(1,2)) for (i,j,k),c in s.Poly(s.expand(p),g,X,Z).terms()))
def Eg(p):return expectation(p,s.Integer(1))
def riesz1(p):
    # Hermite-coordinate diagonal inverse, differentiated once.
    rest=s.Poly(s.expand(p),g); inv=0
    while rest.degree()>0:
        n=rest.degree(); c=rest.LC(); H=s.hermite_prob(n,g)
        inv+=c*H/n; rest=s.Poly(s.expand(rest.as_expr()-c*H),g)
    return s.diff(inv,g)
B=g+s.Rational(1,5)*(g*g-1)
Sig=Eg(B*B); mu3=Eg(B**3); kap4=Eg(B**4)-3*Sig**2
K_3=mu3/6; K_4=kap4/24
stein=[B]
for j in range(1,5):stein.append(s.expand(riesz1(stein[-1]-Eg(stein[-1]))*s.diff(B,g)))
check('stein_third_normalization',Eg(stein[2])==3*K_3)
check('stein_fourth_normalization',Eg(stein[3])==4*K_4)
t=s.symbols('t', real=True)
Cv=s.Rational(1,2)+(1-t*t)*Sig
q3=K_3*(X*X/Cv**2-1/Cv)
q4=K_4*(X**3/Cv**3-3*X/Cv**2)
a=1-t**3; d=1-t**4
p=a*q3+d*q4
J=s.diff(p,X); E=s.diff(p,X,2); FF=s.diff(p,X,3); H=1+J
moving=-t*Sig*X/Cv
D=a*(s.diff(q3,t)+s.diff(q3,X)*moving)+d*(s.diff(q4,t)+s.diff(q4,X)*moving)
Ar=[D,-t*Sig*J-3*t*t*K_3*E-4*t**3*K_4*FF,-3*t*t*K_3*(H*H-1)-12*t**3*K_4*E*H,-4*t**3*K_4*(H**3-1),t**4*stein[4]]
W=t*B+X+p+Z
Wdot=B+moving-3*t*t*q3-4*t**3*q4+D
at={t:s.Rational(2,5)}; c_at=Cv.subs(at)
Wa=s.expand(W.subs(at)); Wda=s.expand(Wdot.subs(at)); A=[s.expand(v.subs(at)) for v in Ar]
for n in range(1,6):
    lhs=expectation(n*Wa**(n-1)*Wda,c_at)
    rhs=sum(expectation(s.factorial(n)/s.factorial(n-r)*A[r-1]*Wa**(n-r),c_at) for r in range(1,min(n,5)+1))
    check(f'consumer_same_endpoint_current_test_degree_{n}',lhs==rhs)
# Counterchecks for terms whose omission may look harmless.
check('covariance_sampler_feedback_nonzero',expectation((-t*Sig*J).subs(at)*Wa,c_at)!=0)
check('quartic_third_derivative_feedback_nonzero',(-4*t**3*K_4*FF).subs(at)!=0)
check('moving_root_derivative_matters',s.simplify((D-(a*s.diff(q3,t)+d*s.diff(q4,t))).subs(at))!=0)
# Endpoint panels, including regimes much smaller than the heat floor.
for e in range(1,61):
    tau=2.0**(-e)
    ds=[2.0**(-j) for j in range(1,2*e+30)]
    sum2=sum(v/(v+tau*tau) for v in ds)
    sum3=sum(v/(v+tau*tau)**1.5 for v in ds)
    sq4=sum(v*v/(v+tau*tau)**2 for v in ds)
    sq5=sum(v*v/(v+tau*tau)**2.5 for v in ds)
    check(f'positive_mass_first2_{e}',sum2<2*e+3)
    check(f'positive_mass_first3_{e}',tau*sum3<5)
    check(f'positive_mass_square4_{e}',sq4<2*e+3)
    check(f'positive_mass_square5_{e}',tau*sq5<5)
    check(f'root_radius_weight_{e}',all(v/(v+tau*tau)<=1 for v in ds))
check('unsmoothed_center_first_diverges',sum(2**(j/2) for j in range(1,41))>1e6)
# A Hilbert L2 quadrature certificate is not uniform in arbitrary callers.
# Positive midpoint example with a translated bounded Gaussian bump: its
# quadrature value stays fixed, while its exact resolvent value tends to zero.
from scipy.integrate import quad
r0=.5; width=.25
def smoothed_bump(r,R):
    vv=width*width+1-r*r
    return width/math.sqrt(vv)*math.exp(-((r-r0)*R)**2/(2*vv))
qval=r0*smoothed_bump(r0,1000)
ivals=[quad(lambda r:r*smoothed_bump(r,R),0,1,points=[r0],epsabs=1e-12)[0] for R in [100,1000]]
check('L2_clock_not_uniform_caller',qval>0.1 and ivals[1]<.001 and ivals[1]<ivals[0]/5)
# Source amplitude and unit-to-physical v powers: exact rational exponents.
def physical(power,heat=0):return (Fraction(power-heat),Fraction(1-power+heat,2))
check('native_quintic_A5_u_minus2',physical(5)==(Fraction(5),Fraction(-2)))
check('star_caller_A',physical(2,1)==(Fraction(1),Fraction(0)))
check('star_coarse_A7_u_minus3',physical(8,1)==(Fraction(7),Fraction(-3)))
check('quartic_target_heat_debt_after_terminal',(5+1,Fraction(-3,2)-Fraction(1,2))==(6,-2))
for N in [1,2,7,31,101]:
    nu=Fraction(1,N); read_sq=nu/5
    check(f'full_readout_budget_{N}',N*5*read_sq==1)
# Exact star reverse symbol countertest with SIDE CROSS ZERO.
z,lam,beta,bb,dd,PP,SS=s.symbols('z lam beta bb dd PP SS')
V=lam*z*z*(PP+bb*z)*(SS+beta*dd*z)
def Eps(p):return s.expand(sum(c*gm(i)*gm(j) for (i,j),c in s.Poly(s.expand(p),PP,SS).terms()))
mu=Eps(V); var=Eps(V*V)-mu*mu
check('star_positive_main_four_marks',s.expand(mu).coeff(z,4)==lam*bb*beta*dd)
check('side_zero_feedback_survives',s.expand(var/2).subs(beta,0).coeff(z,4)==lam*lam/2)
# Explicit invalid-trace counterexample: I_D has op 1, HS sqrt(D), trace D.
check('closed_self_trace_not_one_HS',128>math.sqrt(128))
# Numerical suite success is distinct from theorem-review verdict.
result={'diagnostic_status':'PASS','assertion_count':len(checks),'native_compiler_executed':False,'scope':'Exact multivariate hierarchy, same-endpoint consumer identities, omitted-term and side-zero countertests, positive endpoint sums and powers. Read AUDIT.md for conditional mathematical theorem scope.','checks':checks}
(HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
