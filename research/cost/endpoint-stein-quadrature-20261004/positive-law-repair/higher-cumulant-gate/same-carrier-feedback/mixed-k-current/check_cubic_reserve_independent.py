#!/usr/bin/env python3
"""Independent exact rational finite-reserve algebra; no source producer or path simulation."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib, json
HERE=Path(__file__).resolve().parent
checks=0
def ck(p):
    global checks
    assert p
    checks+=1
def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def tr(A): return [list(z) for z in zip(*A)]
def mul(A,B): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def scale(a,A): return [[a*x for x in row] for row in A]
def add(A,B): return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]
def sub(A,B): return add(A,scale(-1,B))
def sym(A): return scale(F(1,2),add(A,tr(A)))
def det(A):
    if len(A)==1:return A[0][0]
    return sum((-1)**j*A[0][j]*det([row[:j]+row[j+1:] for row in A[1:]]) for j in range(len(A)))
def psd(A):
    ck(A==tr(A))
    for k in range(1,len(A)+1):
        for ix in combinations(range(len(A)),k):
            ck(det([[A[i][j] for j in ix] for i in ix])>=0)

cases=0
clocks=[F(3,5),F(4,5),F(5,13)]
for q,s,r,t in product(clocks,repeat=4):
    cq2,cr2,ct2=1-q*q,1-r*r,1-t*t
    K=[[cq2,r*s*cq2,t*s*cq2],
       [r*s*cq2,1-(r*s*q)**2,r*t*(1-s*s*q*q)],
       [t*s*cq2,r*t*(1-s*s*q*q),1-(t*s*q)**2]]
    h=min(cq2,cr2,ct2); sigma2=h/10
    inverse_frobenius=1/cq2+(1+r*r*s*s)/cr2+(1+t*t*s*s)/ct2
    ck(inverse_frobenius<=5/h)
    psd(sub(K,scale(2*sigma2,eye(3))))
    coarse=sub(K,scale(sigma2,eye(3)))
    psd(coarse)
    ck(add(coarse,scale(sigma2,eye(3)))==K)
    ck(1/sigma2<=10*(1/cq2+1/cr2+1/ct2))
    cases+=1

# The native C0 has two symmetric D-dimensional Jacobian blocks.
HY=[[F(1,3),F(1,10)],[F(1,10),F(1,2)]]
HV=[[F(3,5),F(-1,8)],[F(-1,8),F(1,4)]]
HU=[[F(1,5),F(1,7)],[F(1,7),F(2,3)]]
I=eye(2)
for H in [HY,HV,HU]:
    psd(H); psd(sub(I,H)); ck(H==tr(H))

A=F(1,1000); nu=F(1,16); w=F(1,16)
rhoY=rhoV=A;rhoU=-4*w*A/nu
Ms=[scale(rhoY,HY),scale(rhoV,HV),scale(rhoU,HU)]
C=I; cross=I
for M in Ms:
    innovation=sub(I,mul(M,tr(M)))
    psd(innovation)
    C=add(mul(mul(M,C),tr(M)),innovation)
    cross=mul(M,cross)
    ck(C==I)
expected=scale(rhoU*rhoV*rhoY,mul(mul(HU,HV),HY))
ck(cross==expected)
# c=a=sqrt(nu/2), so ca=nu/2, no irrational arithmetic needed.
output=add(scale(nu,I),scale(nu,sym(cross)))
claimed=sub(scale(nu,I),scale(4*w*A**3,sym(mul(mul(HU,HV),HY))))
ck(output==claimed);psd(output)
wrong=sym(mul(mul(HU,HY),HV))
ck(wrong!=sym(mul(mul(HU,HV),HY)))

# Positive two-node moment-exact clock: sum beta r = 1/2.
rule=[(F(1,2),F(1,4)),(F(1,2),F(3,4))]
ws=[bq*bs*br*bt*q*s*r*t for (bq,q),(bs,s),(br,r),(bt,t) in product(rule,repeat=4)]
ck(sum(ws)==F(1,16))
H=scale(A,HY)
H3=mul(mul(H,H),H)
quad=scale(-4*sum(ws),H3)
ck(quad==scale(F(-1,4),H3))

# Signed-mixture baseline shift keeps the perturbation/centered energy unchanged.
v0=F(1,2)
ck(A**3/4<=v0/4)
Gamma=sub(output,scale(nu,I))
Cprime=add(scale(v0/2,I),Gamma)
psd(sub(Cprime,scale(v0/4,I)))
ck(add(scale(v0/2,I),Cprime)==add(scale(v0,I),Gamma))

source=HERE/'FINITE-C0-CUBIC-COVARIANCE-RESERVE.md'
result={'status':'PASS: exact algebra and source-contract diagnostics',
        'assertions':checks,'exact_shield_cases':cases,
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'scope':'Positive finite reduced mixed-covariance specialization under pinned LOW30 guards; no full endpoint claim'}
(HERE/'cubic_reserve_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
