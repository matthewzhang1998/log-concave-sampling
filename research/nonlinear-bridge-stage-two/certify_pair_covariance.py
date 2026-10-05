#!/usr/bin/env python3
"""Exact rational certificate for the uniform two-innovation conditional gap.
No floating-point or source-dependent calculation is used in the certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import deque
import json, time

ONE=F(1); ZERO=F(0); TOL=F(1,2**72)
@lru_cache(None)
def exp_bounds(x):
    """For rational 0<=x<=16, return strict rational Taylor enclosures."""
    assert 0 <= x <= 16
    if x == 0: return (ONE,ONE)
    term=ONE; total=ONE; n=0
    while True:
        n+=1; term*=x/n; total+=term
        if n+2>x:
            nxt=term*x/(n+1)
            tail=nxt/(1-x/(n+2))
            if tail<TOL: return (total,total+tail)

def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    z=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]
    return(min(z),max(z))
def sq(a):
    lo=ZERO if a[0]<=0<=a[1] else min(a[0]*a[0],a[1]*a[1])
    return(lo,max(a[0]*a[0],a[1]*a[1]))
def scale(c,a):return mul((c,c),a)
@lru_cache(None)
def E_bounds(s):
    # E(2s)=(exp(2s)-1)/(2s), an increasing positive function.
    if s == 0:return(ONE,ONE)
    lo,hi=exp_bounds(2*s)
    return((lo-1)/(2*s),(hi-1)/(2*s))
@lru_cache(None)
def f_interval(l,r):
    # F(s)=2s/E(2s), valid at zero.
    Elo=E_bounds(l)[0]; Ehi=E_bounds(r)[1]
    return(2*l/Ehi,2*r/Elo)
@lru_cache(None)
def eminus_interval(l,r):
    return(1/exp_bounds(2*r)[1],1/exp_bounds(2*l)[0])

def certify(box):
    al,ar,hl,hr=box
    a=(al,ar);h=(hl,hr)
    fa=f_interval(al,ar); fh=f_interval(hl,hr)
    ef=mul(eminus_interval(al,ar),fh)
    r=sub(a,(ONE,ONE)); v=sub(add(scale(F(2),a),h),(ONE,ONE))
    u11=add(fa,ef)
    u12=add(mul(fa,r),mul(ef,v))
    u22=add(mul(fa,sq(r)),mul(ef,sq(v)))
    p=F(9,10)-u11[1];q=F(9,10)-u22[1]
    cross=max(abs(u12[0]),abs(u12[1]))
    return p>0 and q>0 and p*q>cross*cross

def main():
    started=time.time(); stack=[(ZERO,F(8),ZERO,F(8))]
    accepted=split=max_depth=0; depth_stack=[0]; leaves=[]
    while stack:
        box=stack.pop();depth=depth_stack.pop()
        if certify(box):
            accepted+=1;max_depth=max(max_depth,depth);leaves.append(box)
            continue
        al,ar,hl,hr=box
        if ar-al>=hr-hl:
            mid=(al+ar)/2
            stack.extend([(al,mid,hl,hr),(mid,ar,hl,hr)])
        else:
            mid=(hl+hr)/2
            stack.extend([(al,ar,hl,mid),(al,ar,mid,hr)])
        depth_stack.extend([depth+1,depth+1]);split+=1
        if depth>40:raise RuntimeError('Certificate failed to close by depth 40')
    # This union certificate includes h=0, hence ||b(a)||^2 <= 9/10.
    # e > 5/2, and a^4 exp(-2a) is decreasing for a>=8.
    assert 8*F(8)**4*F(2,5)**16 < F(1,64)
    assert F(1038)*F(2,5)**16 < F(1,64)
    assert 1-F(9,10)-F(1,64) > F(1,16)
    result={'result':'PASS','compact_square':[0,8],
      'compact_gap':'1/10','global_gap':'1/16',
      'accepted_rational_boxes':accepted,'bisections':split,'maximum_depth':max_depth,
      'exponential_precision_absolute':'2^-72','cached_exp_endpoints':exp_bounds.cache_info().currsize,
      'elapsed_seconds':round(time.time()-started,3)}
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
