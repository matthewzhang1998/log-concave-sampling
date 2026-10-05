#!/usr/bin/env python3
"""Exact finite algebra checks, not a native compiler or a numerical W2 test."""
import json, math, hashlib
from functools import lru_cache
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parent
checks=[]
def check(label, expr):
    result=s.expand(expr)==0
    if not result:
        result=s.simplify(expr)==0
    assert result, (label, expr)
    checks.append(label)
x,C,t,S=s.symbols('x C t S', real=True, nonzero=True)
def score(n):
    h=s.Integer(1)
    for _ in range(n):
        h=s.expand(x*h/C-s.diff(h,x))
    return h
def gauss1(p,c):
    out=0
    for (n,),coef in s.Poly(s.expand(p),x).terms():
        if n%2==0:
            out+=coef*s.factorial2(n-1)*c**(n//2) if n else coef
    return s.expand(out)
# Every score used in fixed ranks through 11 obeys the fixed-visible-root identity.
for n in range(1,11):
    q=score(n)
    rootdot=-2*t*S*s.diff(q,C)-t*S*x/C*s.diff(q,x)
    rhs=-t*S*s.diff(q,x,2)+t*s.diff(q,x)*S*x/C
    check(f'score heat/root identity degree {n}', rootdot-rhs)
# Nonlinear same-endpoint tests; B and reserve can be fixed arbitrary constants.
for m in range(3,7):
    p=sum(s.Rational((-1)**r,r+1)*(1-s.Rational(1,3)**r)*score(r-1) for r in range(3,m+1))
    D=-t*S*s.diff(p,x,2)+t*s.diff(p,x)*S*x/C
    J=s.diff(p,x)
    vals={C:s.Rational(7,5),t:s.Rational(1,3),S:s.Rational(2,5)}
    p=p.subs(vals); J=J.subs(vals); D=D.subs(vals)
    W=s.Rational(2,7)+x+p
    for degree in range(1,7):
        grad=degree*W**(degree-1)
        hess=degree*(degree-1)*W**max(degree-2,0)
        expr=D*grad-vals[t]*vals[S]*(J+J*J)*hess
        check(f'exact nonlinear endpoint rank {m} test degree {degree}',gauss1(expr,vals[C]))
# Scalar skew-square repair at every polynomial test degree through 14.
a=s.symbols('a')
for d in range(0,15):
    phi=x**d
    lhs=gauss1((a*score(2))**2*s.diff(phi,x,2)/2,C)
    desired=a*a*gauss1(s.diff(phi,x,6),C)/2
    extras=a*a*(gauss1(s.diff(phi,x,2),C)/C**2 +2*gauss1(s.diff(phi,x,4),C)/C)
    check(f'skew-square Wick jet degree {d}',lhs-desired-extras)
    corr=-a*a*score(1)/C**2-2*a*a*score(3)/C
    check(f'skew-square repair degree {d}',lhs+gauss1(corr*s.diff(phi,x),C)-desired)
# Non-diagonal matrix checks, including non-symmetric J.
x1,x2=s.symbols('x1 x2'); xx=s.Matrix([x1,x2]); tt=s.symbols('tt')
Sig=s.Matrix([[2,1],[1,3]])
CC=s.eye(2)+(1-tt**2)*Sig
QQ=CC.inv(); yy=QQ*xx
# Symmetric rank-three tensor with mixed entries; q_i=K_iab H2_ab.
k={(0,0,0):s.Rational(1,11),(0,0,1):s.Rational(2,13),(0,1,1):s.Rational(-1,17),(1,1,1):s.Rational(3,19)}
def kval(i,j,l): return k[tuple(sorted((i,j,l)))]
qs=s.Matrix([sum(kval(i,j,l)*(yy[j]*yy[l]-QQ[j,l]) for j in range(2) for l in range(2)) for i in range(2)])
val=s.Rational(1,3); cc=CC.subs(tt,val); qq=QQ.subs(tt,val)
q=qs.subs(tt,val).applyfunc(s.cancel)
j=q.jacobian([x1,x2])
rootvelocity=-val*Sig*qq*xx
fixedroot=qs.diff(tt).subs(tt,val)+j*rootvelocity
expected=s.Matrix([-val*sum(Sig[u,v]*s.diff(q[i],[x1,x2][u],[x1,x2][v]) for u in range(2) for v in range(2)) for i in range(2)])+val*j*Sig*qq*xx
for i in range(2):check(f'non-diagonal fixed-root identity component {i}',fixedroot[i]-expected[i])
@lru_cache(None)
def moment(n1,n2):
    if (n1+n2)%2:return s.Integer(0)
    if n1==n2==0:return s.Integer(1)
    if n1:
        out=(n1-1)*cc[0,0]*moment(n1-2,n2) if n1>=2 else 0
        if n2:out+=n2*cc[0,1]*moment(n1-1,n2-1)
        return out
    return (n2-1)*cc[1,1]*moment(0,n2-2)
def gauss2(poly):
    return sum(coef*moment(*powers) for powers,coef in s.Poly(s.expand(poly),x1,x2).terms())
p=s.Rational(26,27)*q
J=p.jacobian([x1,x2]); DD=s.Rational(26,27)*fixedroot
W=xx+p+s.Matrix([s.Rational(2,7),s.Rational(-1,5)])
assert s.expand(J[0,1]-J[1,0])!=0
checks.append('fixture Jacobian is genuinely non-symmetric')
for direction in [s.Matrix([1,0]),s.Matrix([1,2]),s.Matrix([-2,3])]:
    z=(direction.T*W)[0]
    for degree in range(1,6):
        grad=degree*z**(degree-1)*direction
        Hess=degree*(degree-1)*z**max(degree-2,0)*(direction*direction.T)
        residual=(DD.T*grad)[0]-val*s.trace((Sig*J.T+J*Sig*J.T).T*Hess)
        check(f'non-diagonal same-endpoint d={list(direction)} degree={degree}',gauss2(residual))
# Arbitrary actual covariance-root gauge: no C/Cdot commutation assumed.
A=s.Matrix([[1,3],[-2,2]])
Cd=A+A.T
assert Cd*cc!=cc*Cd
checks.append('general covariance derivative does not commute with covariance')
Dgeneral=s.Matrix([sum(Cd[u,v]*s.diff(p[i],[x1,x2][u],[x1,x2][v])/2 for u in range(2) for v in range(2)) for i in range(2)])-J*A.T*qq*xx
for degree in range(1,7):
    direction=s.Matrix([1,2]);z=(direction.T*W)[0]
    grad=degree*z**(degree-1)*direction
    Hess=degree*(degree-1)*z**max(degree-2,0)*(direction*direction.T)
    residual=(Dgeneral.T*grad)[0]+s.trace((A*J.T+J*Cd*J.T/2).T*Hess)
    check(f'general noncommuting covariance/root gauge degree {degree}',gauss2(residual))
# Exact scalar raw-map covariance surplus for arbitrary cumulants.
ks=s.symbols('k3:9')
p=sum(ks[r-3]*score(r-1) for r in range(3,9))
check('raw variance surplus through rank eight',gauss1(p*p,C)-sum(math.factorial(r-1)*ks[r-3]**2/C**(r-1) for r in range(3,9)))
check('carrier correction cross covariance zero',gauss1(x*p,C))
# Source regularity and strict skew check for a concrete admissible gradient.
av=s.Rational(1,4)
assert float(3*av*s.exp(-s.Rational(1,2))-8*av**3)>0
checks.append('admissible bounded-Hessian scalar gradient has nonzero skew')
inputs=[
 '/workspace/shared/fourth-cumulant-return-20261005/comparison/POSITIVE-FOURTH-CUMULANT-BUFFERED-CONSUMER.md',
 '/workspace/shared/rank-indexed-positive-returns-20261005/ALL-RANK-APPELL-AND-POSITIVE-CONSUMER.md',
 '/workspace/shared/rank-indexed-positive-returns-20261005/TREE-GENERATOR-AND-RANK5-RETURN-TEST.md',
 '/workspace/shared/two-cycle-positive-return-20261005/TWO-CYCLE-OPENING-AND-RETURN.md']
pins={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs}
result={'status':'PASS','assertions':len(checks),'checks':checks,'input_sha256':pins,
 'scope':'Exact finite Gaussian polynomial and tensor-algebra diagnostics. Not a native compiler execution, W2 numerical test, or all-order native theorem.'}
(ROOT/'covariance_reduction_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','assertions':len(checks)},indent=2))
