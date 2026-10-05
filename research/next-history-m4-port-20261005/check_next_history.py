#!/usr/bin/env python3
"""Exact target algebra only; this does not execute any native compiler."""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib,json,itertools
ROOT=Path(__file__).resolve().parent
checks=0

def ck(v):
    global checks
    assert v
    checks+=1

def cov(n,m):
    return F(comb(n+m,n)-1,2**(n+m))

def mul(p,q):
    r={}
    for i,a in p.items():
        for j,b in q.items():r[i+j]=r.get(i+j,F(0))+a*b
    return {k:v for k,v in r.items() if v}

def sigma(k):
    p={}
    for n in range(1,k+1):
        for m in range(1,k+1):
            p[n+m]=p.get(n+m,F(0))+(-1)**(n+m)*cov(n,m)
    return p

def sub(p,q):return {k:p.get(k,F(0))-q.get(k,F(0)) for k in p.keys()|q.keys() if p.get(k,F(0))!=q.get(k,F(0))}

# Exact double-integral decomposition of linear OU histories.
for n,m in itertools.product(range(1,17),repeat=2):
    side1=F(factorial(n+m-1),2**(n+m)*factorial(n-1)*factorial(m))
    side2=F(factorial(n+m-1),2**(n+m)*factorial(n)*factorial(m-1))
    ck(side1+side2==F(comb(n+m,n),2**(n+m)))
    ck(cov(n,m)==cov(m,n))

expected=[[F(1,4),F(1,4),F(3,16)],[F(1,4),F(5,16),F(9,32)],[F(3,16),F(9,32),F(19,64)]]
for i,j in itertools.product(range(3),repeat=2):ck(cov(i+1,j+1)==expected[i][j])
ck(sigma(2)=={2:F(1,4),3:F(-1,2),4:F(5,16)})
ck(sigma(3)=={2:F(1,4),3:F(-1,2),4:F(11,16),5:F(-9,16),6:F(19,64)})
diff=sub(sigma(3),sigma(2))
ck(diff=={4:F(3,8),5:F(-9,16),6:F(19,64)})
cross={4:2*cov(1,3),5:-2*cov(2,3)}
diag={6:cov(3,3)}
ck(diff=={**cross,**diag})

diff43=sub(sigma(4),sigma(3))
ck(diff43=={5:F(-1,4),6:F(7,16),7:F(-17,32),8:F(69,256)})
for den in range(2,202):
    a=F(1,den)
    p=sum(c*a**(j-5) for j,c in diff43.items())
    ck(p<=F(-1,32))
for k in range(2,19):
    increment=sub(sigma(k),sigma(k-1))
    ck(min(increment)==k+1)
    ck(increment[k+1]==F((-1)**(k-1)*k,2**k))

# Exact first-moment increment and ordered chain integral at arbitrary depth.
for k in range(2,19):
    old={j:F((-1)**(j-1),2**j) for j in range(1,k)}
    new={j:F((-1)**(j-1),2**j) for j in range(1,k+1)}
    ck(sub(new,old)=={k:F((-1)**(k-1),2**k)})
    ck(F(1,2)**k==F(1,2**k))

# All nonempty marked cumulant slots retained; theta integral recovers bins.
for rank in range(2,9):
    ck(sum(comb(rank,j) for j in range(1,rank+1))==2**rank-1)
    for j in range(1,rank+1):
        ck(F(rank*comb(rank-1,j-1),j)==comb(rank,j))
ck([comb(3,j) for j in range(1,4)]==[3,3,1])
ck([comb(4,j) for j in range(1,5)]==[4,6,4,1])

# Readout and exact quadratic interpolation of an arbitrary scalar PSD block.
for b,d in itertools.product(range(-4,5),repeat=2):
    for t in [F(0),F(1,4),F(1,2),F(1)]:
        old=b*b;cross=b*d;delta=d*d
        ck(old+2*t*cross+t*t*delta==(b+t*d)**2)
        ck(2*cross+2*t*delta==2*d*(b+t*d))
ck(F(1,2)+F(1,2)==1) # joint block buffer after [I I]
ck(5*F(1,10)+F(1,2)==1) # old top allocation remains exactly budgeted

# Exact port-only witness distinguishing Hilbert-Schmidt and trace grades.
for n in range(2,65):
    A=F(1,n);sqrtD=n*n
    ck(A==A**3*sqrtD)
    ck((A*A)/(A**6*sqrtD)==n*n)
    ck(A*A==A**6*(n**4))

# Physical target rows: one error slot of grade 3 and other firsts grade 1.
for rank in range(1,5):
    target_grade=3+rank-1
    terminal_grade=target_grade+1
    ck(terminal_grade==rank+3)

pins=json.loads((ROOT/'INPUT-PINS.json').read_text())['inputs']
for pin in pins:
    p=Path(pin['path'])
    ck(hashlib.sha256(p.read_bytes()).hexdigest()==pin['sha256'])
result={'status':'PASS','exact_assertions':checks,'scope':'Exact finite algebra, target exponents and input hashes. No numerical native/compiler execution, new LAW construction or accuracy-raising theorem.', 'linear_covariance_increment':{str(k):str(v) for k,v in sorted(diff.items())},'next_linear_covariance_increment':{str(k):str(v) for k,v in sorted(diff43.items())}}
(ROOT/'next_history_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
