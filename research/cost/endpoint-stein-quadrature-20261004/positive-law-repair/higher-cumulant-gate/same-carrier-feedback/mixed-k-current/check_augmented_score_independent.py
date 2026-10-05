#!/usr/bin/env python3
"""Exact Gaussian-covariance and numerical integral diagnostics; not a proof substitute."""
from fractions import Fraction as F
from itertools import product
import hashlib, json
from pathlib import Path
import mpmath as mp
mp.mp.dps=70
HERE=Path(__file__).resolve().parent
count=0
def check(ok):
    global count
    assert ok
    count+=1

def mm(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def tr(A): return [list(x) for x in zip(*A)]
def add(A,B): return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def sub(A,B): return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]

# Five scalar coordinates [G, vertex_t, vertex_u, V(h), V'(h)].
# All entries are exact rational numbers for this family of clocks/rotations.
clocks=[F(1,5),F(1,2),F(4,5)]
rotations=[(F(0),F(1),F(1)),(F(3,5),F(4,5),F(1,2)),
           (F(4,5),F(3,5),F(1,3)),(F(1),F(0),F(0))]
outer=[(F(0),F(1)),(F(3,5),F(4,5)),(F(4,5),F(3,5))]
cases=0
for a,b,(rho,rotation,alpha),(r,q) in product(clocks,clocks,rotations,outer):
    check(r*r+q*q==1)
    check(rho*rho+rotation*rotation==1)
    check(rho+alpha*rotation==1)
    h2=max(a,b)**2/(1-max(a,b)**2)
    c=min(a/b,b/a)-a*b
    cov=[[F(0) for j in range(5)] for i in range(5)]
    diag=[1,1-a*a+q*q*a*a,1-b*b+q*q*b*b,h2,h2]
    for i,v in enumerate(diag): cov[i][i]=v
    pairs={(0,1):q*a,(0,2):q*b,(1,2):q*q*a*b+rho*c,
           (1,3):rho*a,(1,4):rotation*a,(2,3):b}
    for (i,j),v in pairs.items(): cov[i][j]=cov[j][i]=v
    score=[1/q,F(0),F(0),F(-1),-alpha]
    score_cov=[sum(score[i]*cov[i][j] for i in range(5)) for j in range(5)]
    check(score_cov[1]==0)
    check(score_cov[2]==0)
    score_var=sum(score[i]*score_cov[i] for i in range(5))
    check(score_var==1/q**2+(1+alpha**2)*h2)
    check((1+alpha**2)*h2==2*h2/(1+rho))
    cases+=1

# Noncommuting symmetric vertices verify the transpose/product ordering in (3).
H=[[F(2,5),F(1,10)],[F(1,10),F(1,5)]]
H0=[[F(1,4),F(1,20)],[F(1,20),F(1,3)]]
J=[[F(1,10),F(1,30)],[F(1,30),F(1,5)]]
DxE=sub(sub(H,H0),mm(H,J))
check(add(tr(DxE),mm(J,H))==sub(H,H0))
check(add(tr(DxE),mm(H,J))!=sub(H,H0))

# Split the genuine positive history weight on t<=u and its transpose.
weight=mp.quad(lambda m: mp.exp(-2*m)*(1-mp.exp(-2*m)),[0,1,mp.inf])
base=mp.quad(lambda y: mp.sqrt(y*(1-y))/2,[0,1])
rho_int=mp.quad(lambda rho: mp.sqrt(2/(1+rho)),[0,1])
shield=base*rho_int
outer_shield=mp.quad(lambda r: r/mp.sqrt(1-r*r),[0,mp.mpf('0.5'),1])
check(abs(weight-mp.mpf(1)/4)<mp.mpf('1e-60'))
check(abs(base-mp.pi/16)<mp.mpf('1e-60'))
check(abs(shield-mp.pi*(2-mp.sqrt(2))/8)<mp.mpf('1e-60'))
check(abs(outer_shield-1)<mp.mpf('1e-34'))
constant=mp.mpf(3)/8+shield/2
check(constant<mp.mpf('0.491'))

# Exact linear-history h response: integral e^-v L_v h dv=m/(exp(2m)-1).
for m in [mp.mpf('0.001'),mp.mpf('0.1'),mp.mpf('1'),mp.mpf('10')]:
    calc=mp.quad(lambda v: mp.exp(-2*v)*mp.expm1(2*v)/mp.expm1(2*m),[0,m])
    calc+=mp.exp(-2*m)/2
    check(abs(calc-m/mp.expm1(2*m))<mp.mpf('1e-60'))
    check(0<calc<=mp.mpf(1)/2)

pins={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in [
    'AUGMENTED-SCORE-REMOVES-THE-NONLINEAR-K-HISTORY.md',
    'INDEPENDENT-AUGMENTED-SCORE-REGROUPING-AUDIT.md']}
result={'status':'PASS: analytical diagnostics only','assertions':count,
        'exact_conditional_covariance_cases':cases,'weight':str(weight),
        'weighted_bridge_shield':str(shield),'R2_error_constant':str(constant),
        'pins':pins,'finite_VALUE_producer':'NOT CLAIMED'}
(HERE/'augmented_score_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
