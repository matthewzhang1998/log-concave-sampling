#!/usr/bin/env python3
"""Exact monomial audit only: no native/compiler/combined producer execution."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

checks=0
def check(p):
    global checks
    assert p
    checks+=1

alpha,beta,gamma,theta=F(5,6),F(3,2),F(7,4),F(1,2)
rows=[(F(5),F(3,2)),(F(5),F(1)),(4+theta,F(1,2)),
      (F(5),F(1,2)),(6-theta/2,F(3,2)),(F(6),F(3,2))]
check(rows==[(F(5),F(3,2)),(F(5),F(1)),(F(9,2),F(1,2)),
             (F(5),F(1,2)),(F(23,4),F(3,2)),(F(6),F(3,2))])
report=[]
for a,p in rows:
    exponents=[a-p*alpha]
    if p>1: exponents.append(a-(p-1)*beta)
    if p==1: exponents.append(a)  # explicit logarithm at this grade
    for e in exponents: check(e>=F(15,4))
    report.append({'monomial':f'A^{a} u^-{p}', 'integrated_grades':list(map(str,exponents))})
check(1-beta/2==F(1,4))
check(2-theta/2-beta==F(1,4))
check(2-theta/2-beta/2==1)
check(2-theta/2-alpha/2==F(4,3))
check(1-beta/3==F(1,2))
check(3+beta/2==F(15,4))
check(3+beta-alpha/2==F(49,12))
check(2+gamma==F(15,4))
check(3+3*alpha-gamma==F(15,4))
check(5-F(3,2)*alpha==F(15,4))
check(1+2*gamma==F(9,2))
check(3+alpha==F(23,6))
# The old unmodified graph would still have these inadequate forcing rows.
check(2+F(3,2)*alpha==F(13,4))
check(4-alpha==F(19,6))
# Optional final unit-variance completion must retain mu=A.
check(min(F(4),F(4),3+theta,F(4),5-theta/2,F(5))==F(7,2))
check(min(F(4),F(4),F(4),F(4),F(9,2),F(5))==4)
# Positive weighted average of the three exponent rows is always 15/4.
for a in range(1,31):
    for g in range(1,31):
        x,y=F(a,15),F(g,10)
        e=(2+y,3+3*x-y,5-F(3,2)*x)
        check(e[0]/4+e[1]/4+e[2]/2==F(15,4))
        check(min(e)<=F(15,4))

root=Path(__file__).resolve().parent.parent
result={'status':'PASS','checks':checks,'mu':'A^(1/2)',
        'scales':{'w':'A^(5/6)','h':'A^(7/4)','eta':'A^(3/2)'},
        'raw_split_integrated_rows':report,
        'companion_sha256':hashlib.sha256((root/'POST-SYNTHESIS-REENTRY-AND-NEXT-PORT.md').read_bytes()).hexdigest(),
        'scope':'Exact native-scaling and positive-bound-ledger algebra. No combined producer or native compiler execution; the original unmodified graph does not reach this grade.'}
print(json.dumps(result,indent=2))
Path(__file__).with_name('padding-ceiling-checks.json').write_text(json.dumps(result,indent=2)+'\n')
