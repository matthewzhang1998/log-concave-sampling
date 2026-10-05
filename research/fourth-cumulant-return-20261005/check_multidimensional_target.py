#!/usr/bin/env python3
"""Independent coordinate check of four-tree formula in dimension two."""
import sympy as s
from functools import lru_cache
from pathlib import Path
import json
x,y,t,z=s.symbols('x y theta1 theta2'); xs=[x,y]
@lru_cache(None)
def mon(k,a,b):
    ans=x**a*y**b
    if a>=2:ans+=a*(a-1)*mon(k,a-2,b)
    if b>=2:ans+=b*(b-1)*mon(k,a,b-2)
    return s.expand(ans/(k+a+b))
def R(p,k):
    return s.expand(sum(c*mon(k,a,b) for (a,b),c in s.Poly(s.expand(p),x,y).terms()))
def zero(p):return s.expand(p)==0
checks=0
for potential in [x**3+s.Rational(2,3)*y**3+x*y, x**4/4+y**4/7+x*x*y*y/3+x*y*y, x*x*y+y*y*x+x*y]:
    g=[s.diff(potential,x),s.diff(potential,y)]
    f=t*g[0]+z*g[1]
    b=s.Matrix([R(s.diff(f,a),2) for a in xs])
    T=b.jacobian(xs)
    SS=[T.diff(a) for a in xs]
    U=s.Matrix([R((T*b)[i],3) for i in range(2)])
    C=2*R((b.T*b)[0],2)
    DC=s.Matrix([s.diff(C,a) for a in xs])
    K3=6*R((b.T*DC)[0],3)
    DK3=s.Matrix([s.diff(K3,a) for a in xs])
    K4=R(8*(b.T*DK3)[0]+6*(DC.T*DC)[0],4)
    P1=R(sum(b[a]*R(sum(T[a,i]*U[i] for i in range(2)),4) for a in range(2)),4)
    P2=R(sum(b[a]*R(sum(b[i]*R(sum(SS[a][i,j]*b[j] for j in range(2)),4) for i in range(2)),4) for a in range(2)),4)
    P3=R(sum(b[a]*R(sum(b[i]*R(sum(T[i,j]*T[a,j] for j in range(2)),4) for i in range(2)),4) for a in range(2)),4)
    P4=R((U.T*U)[0],4)
    assert zero(K4-192*(P1+P2+P3)-96*P4)
    checks+=1
    assert zero(K3-24*R((b.T*U)[0],3));checks+=1
    for a in range(2):assert zero(DC[a]-4*U[a]);checks+=1
    if s.total_degree(potential)==4:assert not zero(P2);checks+=1
out={'status':'PASS','assertions':checks,'scope':'Exact two-dimensional theta-polarized polynomial identities; mixed noncommuting Hessians, all four tensor histories, nonzero star fixture.'}
Path(__file__).with_name('multidimensional_target_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
