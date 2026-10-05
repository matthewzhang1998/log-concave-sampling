#!/usr/bin/env python3
"""Exact rational checks of next-history target algebra; not native execution."""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
from pathlib import Path
import json
import sympy as s

checks = {}
def ok(group, statement):
    assert statement, group
    checks[group] = checks.get(group, 0) + 1

A, theta = s.symbols('A theta')
# Brownian-kernel integration. h_j = sqrt(2)e^-r p_j(r).
def kernel_cov(n,m):
    return sum((F(2*factorial(a+b), factorial(a)*factorial(b)*2**(n+m+1))
                for a in range(n) for b in range(m)), F(0))
for n,m in product(range(1,13), repeat=2):
    c=kernel_cov(n,m)
    ok('brownian_kernel_binomial',c==F(comb(n+m,n)-1,2**(n+m)))
C=s.Matrix([[s.Rational(kernel_cov(n,m).numerator,kernel_cov(n,m).denominator)
             for m in range(1,4)] for n in range(1,4)])
a=s.Matrix([A,-A**2,A**3]); b=s.Matrix([A,-A**2,0]); d=a-b
var3=(a.T*C*a)[0]; var2=(b.T*C*b)[0]
ok('linear_covariance',s.expand(var3-(A**2/s.Integer(4)-A**3/2+11*A**4/16-9*A**5/16+19*A**6/64))==0)
ok('linear_covariance',s.expand(var3-var2-(3*A**4/8-9*A**5/16+19*A**6/64))==0)
ok('linear_covariance',s.expand((b.T*C*d)[0]-(3*A**4/16-9*A**5/32))==0)
ok('linear_covariance',s.expand((d.T*C*d)[0]-19*A**6/64)==0)
ok('linear_covariance',s.det(C)>0)
means=[sum((-1)**(j-1)*A**j/s.Integer(2)**j for j in range(1,k+1)) for k in range(1,5)]
ok('linear_mean',s.expand(means[3]-means[2]+A**4/16)==0)
# A genuine same-g negative increment one history later.
C4=s.Matrix([[s.Rational(kernel_cov(n,m).numerator,kernel_cov(n,m).denominator)
              for m in range(1,5)] for n in range(1,5)])
a4=s.Matrix([A,-A**2,A**3,-A**4]);a3=s.Matrix([A,-A**2,A**3,0])
j4=s.expand((a4.T*C4*a4-a3.T*C4*a3)[0])
expected=-A**5/4+7*A**6/16-17*A**7/32+69*A**8/256
ok('next_linear_negative_increment',s.expand(j4-expected)==0)
ok('next_linear_negative_increment',s.expand(j4/A**5-(-s.Rational(1,4)+7*A/16+A**2*(-s.Rational(17,32)+69*A/256)))==0)
for n in range(2,21):
    ok('next_linear_negative_increment',j4.subs(A,s.Rational(1,n))<0)
# Rational finite probability fixture. Rows are samples; V and R have shared tape.
weights=[F(i,28) for i in range(1,8)]
V0=[[F(((j+2)*(i+1)**2+3*i)%17-8,7) for j in range(3)] for i in range(7)]
R0=[[F(((j+3)*(i+2)**3+i)%19-9,11) for j in range(3)] for i in range(7)]
def center(M):
    mu=[sum((w*row[j] for w,row in zip(weights,M)),F(0)) for j in range(3)]
    return [[row[j]-mu[j] for j in range(3)] for row in M]
V,R=center(V0),center(R0)
U=[[v+r for v,r in zip(vv,rr)] for vv,rr in zip(V,R)]
def moment(slots,idx):
    out=F(0)
    for a,w in enumerate(weights):
        val=w
        for source,j in zip(slots,idx):val*=source[a][j]
        out+=val
    return out
def cumulant(slots,idx):
    ans=moment(slots,idx)
    if len(slots)==4:
        for pairs in [((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]:
            p,q=pairs
            ans-=moment([slots[i] for i in p],[idx[i] for i in p])*moment([slots[i] for i in q],[idx[i] for i in q])
    return ans
for rank in [2,3,4]:
    for idx in product(range(3),repeat=rank):
        delta=cumulant([U]*rank,idx)-cumulant([V]*rank,idx)
        expansion=sum((cumulant([R if mask>>j&1 else V for j in range(rank)],idx)
                       for mask in range(1,1<<rank)),F(0))
        ok(f'full_rank_{rank}_expansion',delta==expansion)
        # Evaluate the actual marked current at rational interior theta values,
        # independently of the multilinear polynomial coefficients.
        for th in [F(1,4),F(2,3)]:
            UT=[[v+th*r for v,r in zip(vv,rr)] for vv,rr in zip(V,R)]
            current=sum((cumulant([R if j==marked else UT for j in range(rank)],idx)
                         for marked in range(rank)),F(0))
            derivative=sum((mask.bit_count()*th**(mask.bit_count()-1)*
                 cumulant([R if mask>>j&1 else V for j in range(rank)],idx)
                 for mask in range(1,1<<rank)),F(0))
            ok(f'rank_{rank}_interpolation',current==derivative)
for idx in product(range(3),repeat=4):
    i,j,k,l=idx
    cv=lambda X,Y,a,b:moment([X,Y],[a,b])
    c=lambda a,b:cv(V,V,a,b)
    q=lambda a,b:cv(U,U,a,b)-c(a,b)
    P=lambda X,Y:X(i,j)*Y(k,l)+X(i,k)*Y(j,l)+X(i,l)*Y(j,k)
    delta_moment=moment([U]*4,idx)-moment([V]*4,idx)
    rhs=delta_moment-P(c,q)-P(q,c)-P(q,q)
    ok('fourth_pairing_increment',rhs==cumulant([U]*4,idx)-cumulant([V]*4,idx))
# Full coherent block covariance and [I theta I] contraction.
mat=lambda X,Y:s.Matrix([[s.Rational(moment([X,Y],[i,j]).numerator,
                                      moment([X,Y],[i,j]).denominator)
                          for j in range(3)] for i in range(3)])
CV,CR,CVR=mat(V,V),mat(R,R),mat(V,R)
block=CV.row_join(CVR).col_join(CVR.T.row_join(CR))
read=s.eye(3).row_join(theta*s.eye(3))
path=read*block*read.T
for i,j in product(range(3),repeat=2):
    ok('block_gram_contraction',s.expand(path[i,j]-(CV+theta*(CVR+CVR.T)+theta**2*CR)[i,j])==0)
    ok('block_gram_derivative',s.expand(s.diff(path[i,j],theta)-(CVR+CVR.T+2*theta*CR)[i,j])==0)
    ok('block_gram_endpoint',s.expand(path[i,j].subs(theta,1)-mat(U,U)[i,j])==0)
# Noncommuting response-chain order and three minus signs.
J2=s.Matrix([[1,1],[0,1]]);J3=s.Matrix([[1,0],[1,1]]);J4=s.Matrix([[2,0],[0,1]]);z=s.Matrix([1,3])
d2=-J2*z;d3=-J3*d2;d4=-J4*d3
ok('response_chain',d4==-J4*J3*J2*z)
ok('response_chain',d4!=-J2*J3*J4*z)
for n in range(2,21):
    a=F(1,n); D=n**4
    ok('dimension_counterexample',a==a**3*n*n)
    ok('dimension_counterexample',a*a/(a**6*n*n)==n*n)
    ok('negative_generic_increment',-2*a**4+a**6<0)
# Weak current identity for polynomial tests of a Gaussian-convolved,
# finitely supported coherent pair, integrated over the Gaussian exactly.
q,var=s.symbols('q var',positive=True)
def noise_moment(k):
    return 0 if k%2 else s.factorial2(k-1)*var**(k//2) if k else s.Integer(1)
def shifted_moment(power,shift):
    return sum(s.binomial(power,k)*noise_moment(k)*shift**(power-k) for k in range(power+1))
for degree in range(1,7):
    expectation=s.Integer(0); current=s.Integer(0)
    for w,vv,rr in zip(weights,V0,R0):
        ww=s.Rational(w.numerator,w.denominator)
        xx=s.Rational(vv[0].numerator,vv[0].denominator)
        dd=s.Rational(rr[0].numerator,rr[0].denominator)
        shift=-q*(xx+theta*dd)
        expectation+=ww*shifted_moment(degree,shift)
        current+=ww*(-q*dd)*degree*shifted_moment(degree-1,shift)
    ok('positive_marked_current',s.expand(s.diff(expectation,theta)-current)==0)
read1=s.eye(3).row_join(s.eye(3))
buffered=read1*(var*s.eye(6)/2+block)*read1.T
for i,j in product(range(3),repeat=2):
    ok('buffer_readout',s.expand(buffered[i,j]-(var*s.eye(3)+mat(U,U))[i,j])==0)
report={'status':'passed','total_exact_assertions':sum(checks.values()),'groups':checks,
        'scope':'Exact rational/symbolic target checks only. No native producer or numerical law comparison is executed.'}
Path(__file__).with_name('exact_target_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
