#!/usr/bin/env python3
"""Algebra/read-set diagnostics; not a replacement for the conditional-source proof."""
import json, itertools, hashlib
from pathlib import Path
import sympy as s
import numpy as np
from numpy.polynomial.hermite import hermgauss
BASE=Path(__file__).parent
checks=0

def check(x, label):
    global checks
    if not bool(x): raise AssertionError(label)
    checks+=1

def max_cut(t):
    n=t.ndim; out=0.
    for bits in itertools.product((0,1), repeat=n):
        if not any(bits) or all(bits): continue
        left=[i for i,b in enumerate(bits) if b]
        right=[i for i,b in enumerate(bits) if not b]
        mat=t.transpose(left+right).reshape(int(np.prod([t.shape[i] for i in left])), -1)
        out=max(out,float(np.linalg.norm(mat,2)))
    return out

la,x,a,b,z=s.symbols('lambda theta a b z')
# Full scalar Gaussian bilinear integral, obtained from det(I-H) and its inverse.
B=x**2; u=a*x; v=b*x
H=s.Matrix([[0,la*B],[la*B,0]])
f=s.Matrix([la*B*v,la*B*u])
direct=la*u*B*v-s.log((s.eye(2)-H).det())/2+(f.T*(s.eye(2)-H).inv()*f)[0]/2
closed=-s.log(1-la**2*B**2)/2+(la*u*B*v+la**2*B**2*(u*u+v*v)/2)/(1-la**2*B**2)
check(s.simplify(direct-closed)==0,'Gaussian determinant/inverse formula')
trunc=s.series(closed,la,0,9).removeO().expand()
family=sum(la**(2*h)*B**(2*h)/s.Integer(2*h) for h in range(1,5))
family+=sum(la**(2*h)*B**(2*h)*(u*u+v*v)/2 for h in range(1,5))
family+=sum(la**(2*h+1)*u*B**(2*h+1)*v for h in range(4))
check(s.expand(trunc-family)==0,'all three coefficient families through lambda^8')
zero=s.simplify(closed.subs({a:0,b:0})+s.log(1-la**2*x**4)/2)
check(zero==0,'zero-side current survives exactly')
check(s.expand(trunc).coeff(la,2).subs({a:0,b:0})==x**4/2,'first zero-side feedback coefficient')
# Grade substitution lambda=A^4, u,v contain A. Physical rank is theta degree.
A=s.symbols('A')
graded=s.expand(trunc.subs({la:A**4,a:A,b:A}))
first=sorted((p[0],p[1]) for p,c in s.Poly(graded,A,x).terms() if c)
check(first[0]==(6,4),'intended main force six/rank four')
check(first[1]==(8,4),'first intrinsic feedback force eight/rank four')
check(first[2]==(10,6),'one-side family force ten/rank six')
for h in range(1,9):
    # Ring consists of 2h four-vertex spines, each with two permanent physical leaves.
    nv=8*h; ne=8*h; marks=4*h; extra=4*h
    check(ne-nv+1==1,'ring has exactly one cycle')
    check(extra==marks-2+2*(ne-nv+1),'cycle-aware derivative budget')
    # One-side family is an open chain of 2h spines with two cap leaves.
    nv2=8*h+2; ne2=8*h+1; m2=4*h+2
    check(ne2==nv2-1 and extra==m2-2,'one-side offspring is a tree')
# Literal leading normalization.
w,alpha,c0,a0,a1,a2,oldb,sigma=s.symbols('w alpha c0 a0 a1 a2 oldb sigma', nonzero=True)
d=8*alpha**4*w**2/(oldb**2*sigma**2)
rroot=-d*A/(c0*a0*a1*a2)
check(s.simplify(c0*a0*a1*a2*rroot*A**5+d*A**6)==0,'negative root normalization')
# Bounded field example: side-zero mixture is valid with |cross| < 1/4.
nodes,weights=hermgauss(80)
q=nodes*np.sqrt(2); weights=weights/np.sqrt(np.pi)
eh2=float(np.sum(weights*np.tanh(q)**2))
check(eh2>0 and eh2<1,'bounded odd-field energy')
lam=.1; read_c=.3; read_a=.3
log4=(read_c*read_a*lam)**2*eh2**2/2
check(log4>0,'strictly positive bounded-field fourth connected coefficient')
# Actual two-spine tensor tests: no initial self loops, both blocks retain physical slots.
rng=np.random.default_rng(20261005)
max_ratio=0
for dim in (2,3):
    for case in range(12):
        raw=rng.normal(size=(dim,)*4)
        raw/=max_cut(raw)
        cycle=np.einsum('ijab,klab->ijkl',raw,raw)
        cuts=max_cut(cycle); hs=float(np.linalg.norm(cycle))
        check(cuts<=1+2e-12,'two-spine cycle proper cuts')
        check(hs<=np.sqrt(dim)+2e-12,'two-spine cycle one Hilbert')
        check(hs<=np.linalg.norm(raw)+2e-12,'marked block Hilbert refinement')
        max_ratio=max(max_ratio,hs/np.sqrt(dim))
# Longer ring, all physical slots remain exposed. Small dim keeps every cut explicit.
for case in range(4):
    dim=2
    raw=rng.normal(size=(dim,)*4); raw/=max_cut(raw)
    ring=np.einsum('ijab,klcb,mncd,opad->ijklmnop',raw,raw,raw,raw)
    check(max_cut(ring)<=1+2e-12,'four-spine ring all proper cuts')
    check(np.linalg.norm(ring)<=np.sqrt(dim)+2e-12,'four-spine ring one Hilbert')
# Negative control: arbitrary internal self contractions are not covered.
for dim in (2,3,4):
    e=np.zeros(dim); e[0]=1
    ident=np.eye(dim)
    t=np.einsum('i,ab,cd->iabcd',e,ident,ident)/dim
    contracted=np.einsum('iaabb->i',t)
    check(max_cut(t)<=1+2e-12,'self-loop fixture has all cuts at most one')
    check(abs(np.linalg.norm(contracted)-dim)<2e-12,'two self loops produce D loss')
# Exact read-set record and force graph of the quartic packet.
verts=['L','C','Cp','Lp','R','Rp']
edges=[('L','C'),('C','Cp'),('Cp','Lp'),('C','R'),('Cp','Rp')]
marks={v:int(v not in ['C','Cp']) for v in verts}
degree={v:sum(v in e for e in edges) for v in verts}
orders={v:degree[v]+marks[v]-1 for v in verts}
check(all(j>=1 for j in orders.values()),'all six primitive native vertices legal')
check(sum(j-1 for j in orders.values())==2,'main is a rank-four tree')
record={
  'coarse_bank':'same actual q within all six force occurrences; q+ and q- retain shared rotation bank',
  'side_outputs':['Y_R','Y_Rp'],
  'side_outputs_are':'external probes of center source; observed after their own selected-pair private tapes are discarded',
  'root_spine':['Lp','Cp','C','L'],
  'readout':['Y_L','P_Lp','P_R','P_Rp','Z'],
  'amplitude':'root=-8 alpha^4 w^2 A/(old_b^2 sigma_2^2 c0 a0 a1 a2); all five nonroots=A',
  'limitations':['algebra/read-set diagnostics do not certify any new generic cycle producer','arbitrary self-contractions excluded','all mean/covariance/tilt descendants still require admitted producer types']
}
result={'status':'PASS','checks':checks,'bounded_field_E_tanh2':eh2,'bounded_field_fourth_log_coefficient':log4,'max_tested_HS_over_sqrtD':max_ratio,'read_set':record}
(BASE/'conditional_quartic_cycle_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
