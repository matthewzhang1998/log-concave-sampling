#!/usr/bin/env python3
"""Exact-rational port/grade diagnostics, not execution of the native compiler."""
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

checks=0

def ok(test):
    global checks
    assert test
    checks+=1

old=[(Q(2),Q(3,2)),(Q(3),Q(1)),(Q(3),Q(1,2)),(Q(4),Q(-1)),(Q(4),Q(1)),(Q(4),Q(0)),(Q(5),Q(-3,2)),(Q(6),Q(-7,6))]

def optimum(lines):
    nodes={Q(0),Q(1)}
    for a,b in lines:
        for c,d in lines:
            if b!=d:
                x=(c-a)/(b-d)
                if 0<=x<=1: nodes.add(x)
    values=[(min(a+b*x for a,b in lines),x) for x in nodes]
    best=max(v for v,x in values)
    return best, sorted(x for v,x in values if v==best)

ok(optimum(old)==(Q(16,5),[Q(4,5)]))
no_prefix=old[1:]
ok(optimum(no_prefix)==(Q(10,3),[Q(2,3)]))
rank3=old.copy(); rank3[3]=(Q(5),Q(-3,2))
ok(optimum(rank3)==(Q(7,2),[Q(1)]))
for n in range(1,1001):
    alpha=Q(n,1000)
    ok(min(a+b*alpha for a,b in old)<=Q(16,5))
    ok(Q(3,2)-alpha>0)
    ok(Q(3,2)-alpha/2>=1)
    ok(1-alpha/3>0)

# Conditional covariance of gamma-kernel integrals R_time^(i+1)X.
def cij(i,j):
    return Q(comb(i+j+1,i)+comb(i+j+1,j)-1,2**(i+j+2))
def var_force(A,k):
    return sum((-1)**(i+j)*A**(i+j+2)*cij(i,j) for i in range(k) for j in range(k))
ok(cij(0,0)==Q(1,4))
ok(cij(1,0)==Q(1,4))
ok(cij(1,1)==Q(5,16))
ok(Q(7,12)**2-Q(1,4)==Q(13,144))
for den in [4,8,16,32,64,128,256,512,1024]:
    A=Q(1,den)
    ok(var_force(A,2)==A*A/4-A**3/2+5*A**4/16)
    for k in range(2,11):
        v=var_force(A,k)
        ok(v>0)
        ok(abs(v-A*A/4)<=A**3)
        ok(abs(v-var_force(A,2))<=2*A**4)
        # Shared carrier B alone has a fixed grade-two gap.
        ok(A*A*Q(7,12)**2-v>=Q(13,144)*A*A)
        # Exact linear conditional mean recurrence.
        m=sum((-1)**(j-1)*(A/2)**j for j in range(1,k+1))
        m_inf=A/(2+A)
        ok(abs(m-m_inf)==(A/2)**(k+1)/(1+A/2))

# An exact marginal law can carry different residual covariance.
# Rotation parametrized rationally: c=(1-h²)/(1+h²), s=2h/(1+h²).
for den in range(2,102):
    h=Q(1,den); c=(1-h*h)/(1+h*h); s=2*h/(1+h*h)
    ok(c*c+s*s==1)
    ok((c-1)**2+s*s==2*(1-c))
    ok(2*(1-c)>0)

result={
    'checks':checks,
    'status':'PASS',
    'old_ledger':{'alpha':'4/5','grade':'16/5'},
    'prefix_only_prospective':{'alpha':'2/3','grade':'10/3'},
    'rank3_only_prospective':{'alpha':'1','grade':'7/2'},
    'raw_covariance_gap_lower_coefficient':'13/144 before O(Lambda A^3) remainder',
    'scope':'Exact rational scalar ledger, linear OU covariance, and carrier-coupling checks. No numerical execution of a native source/compiler; prospective ports are not constructed.'
}
print(json.dumps(result,indent=2))
Path(__file__).with_name('check-results.json').write_text(json.dumps(result,indent=2)+'\n')
