#!/usr/bin/env python3
"""Exact two-scale geometry and exponent diagnostics, not producer execution."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parent.parent
checks=0

def check(b):
    global checks
    assert b
    checks+=1

# Exponent lines a+b*alpha+c*beta on 0<=alpha<=beta<=1.
ledger=[(F(2),F(3,2),F(0)),(F(3),F(1),F(0)),
        (F(3),F(-1,2),F(1)),(F(3),F(0),F(1,2)),
        (F(4),F(-1),F(0)),(F(4),F(0),F(0)),
        (F(4),F(0),F(1)),(F(4),F(0),F(0)),
        (F(5),F(-3,2),F(0)),(F(5),F(0),F(-1,2)),
        (F(6),F(-7,6),F(0)),(F(6),F(0),F(-1,6))]

def solve(rows):
    a=[list(row) for row in rows]
    for col in range(3):
        p=next((r for r in range(col,3) if a[r][col]),None)
        if p is None: return None
        a[col],a[p]=a[p],a[col]
        pivot=a[col][col]
        a[col]=[x/pivot for x in a[col]]
        for r in range(3):
            if r!=col:
                m=a[r][col]
                a[r]=[x-m*y for x,y in zip(a[r],a[col])]
    return tuple(a[r][-1] for r in range(3))

def optimize(lines, beta_max=F(1)):
    # Each row represents a*alpha+b*beta+c*grade <= d.
    rows=[(F(-1),F(0),F(0),F(0)),(F(1),F(-1),F(0),F(0)),
          (F(0),F(1),F(0),beta_max)]
    rows += [(-b,-c,F(1),a) for a,b,c in lines]
    points=set()
    for active in combinations(rows,3):
        point=solve(active)
        if point is not None and all(sum(c*x for c,x in zip(row[:3],point))<=row[3] for row in rows):
            points.add(point)
    best=max(x[2] for x in points)
    return best, sorted(x[:2] for x in points if x[2]==best)

check(optimize(ledger)==(F(16,5),[(F(4,5),F(4,5)),(F(4,5),F(1))]))
check(optimize(ledger[1:])==(F(7,2),[(F(1,2),F(1))]))
# Closed relaxation beta<=3/2 bounds the open native window beta<3/2;
# both optimal grades are attained strictly inside that window.
check(optimize(ledger, F(3,2))==(F(16,5),[(F(4,5),F(4,5)),(F(4,5),F(3,2))]))
check(optimize(ledger[1:], F(3,2))==(F(7,2),[(F(1,2),F(1)),(F(1,2),F(3,2))]))
for j in range(101,150):
    beta=F(j,100)
    check(1-beta/2>0)
    check(F(3,2)-beta>0)
    check(1-beta/3>0)

for j in range(1,21):
    w=F(j,40)
    q=1-w
    sigma2=1-q*q
    for i in range(1,j+1):
        eta=F(i,40)
        t=1-eta
        c2=1-t*t
        v=c2*sigma2/(1-q*q*t*t)
        check(1/v==q*q/sigma2+1/c2)
        check(v>=F(3,4)*eta)
        check(v<=2*w)
        check(eta*eta/w<=eta)  # square of eta/sqrt(w)<=sqrt(eta)
        shares=[v/2,v/(4*q*q),v/(4*q*q)]
        check(min(shares)>=F(3,16)*eta)

out={
    'status':'PASS', 'checks':checks,
    'addendum_sha256':hashlib.sha256((root/'INDEPENDENT-OUTER-CUTOFF-ADDENDUM.md').read_bytes()).hexdigest(),
    'original_optimum':{'grade':'16/5','alpha':'4/5','beta_interval':'[4/5,1]'},
    'prefix_deleted_prospective':{'grade':'7/2','alpha':'1/2','beta':'1'},
    'aggregate_port_extension':{'beta_native_window':'beta<3/2','original_optimum_beta':'[4/5,3/2)','prefix_deleted_optimum_beta':'[1,3/2)'},
    'scope':'Exact rational two-scale Gaussian geometry and complete vertex enumeration of the scalar exponent ledgers. Grade-four logarithms remain explicit in the theorem; no native compiler or improved producer is executed.'
}
print(json.dumps(out,indent=2))
Path(__file__).with_name('independent-cutoff-checks.json').write_text(json.dumps(out,indent=2)+'\n')
