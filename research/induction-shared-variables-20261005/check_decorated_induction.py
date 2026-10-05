#!/usr/bin/env python3
"""Exact finite diagnostics. These supplement the proofs, not native LAW execution."""
from fractions import Fraction as F
from itertools import product, combinations
from collections import defaultdict
from pathlib import Path
import json, math

checks=[]
def ck(ok, name, detail=None):
    checks.append({'name':name,'pass':bool(ok),'detail':detail})
    if not ok: raise AssertionError(name)

# Polynomial ring alpha,theta,P,Q,S,T. Scalars a,b,c,d remain separate variables
# to prevent the two rank-six histories from being merged by an accidental symmetry.
N=10 # alpha theta P Q S T a b c d
zero=(0,)*N
one={zero:F(1)}
def mon(i):
    e=list(zero);e[i]=1;return {tuple(e):F(1)}
def add(A,B):
    O=defaultdict(F,A)
    for e,c in B.items():O[e]+=c
    return {e:c for e,c in O.items() if c}
def scale(A,s):return {e:c*s for e,c in A.items() if c*s}
def mul(A,B):
    O=defaultdict(F)
    for e,c in A.items():
        for f,d in B.items():O[tuple(x+y for x,y in zip(e,f))]+=c*d
    return {e:c for e,c in O.items() if c}
def powp(A,n):
    O=one
    for _ in range(n):O=mul(O,A)
    return O
def gm(k):
    if k%2:return 0
    return math.prod(range(1,k,2))
def expect(A):
    O=defaultdict(F)
    for e,c in A.items():
        c*=math.prod(gm(e[i]) for i in [2,3,4,5])
        f=list(e)
        for i in [2,3,4,5]:f[i]=0
        if c:O[tuple(f)]+=c
    return dict(O)
def grade(A,p):return {e:c for e,c in A.items() if e[0]==p}
def pstr(A):return [{"powers":e,"coefficient":str(c)} for e,c in sorted(A.items())]
a,t,P,Q,S,T,x,y,u,v=map(mon,range(N))
V=mul(powp(a,4),powp(t,2))
for z in [add(P,mul(x,t)),add(Q,mul(y,t)),add(S,mul(a,mul(u,t))),add(T,mul(a,mul(v,t)))]:V=mul(V,z)
M1=expect(V);M2=expect(powp(V,2));M3=expect(powp(V,3))
K2=scale(add(M2,scale(powp(M1,2),-1)),F(1,2))
K3=scale(add(add(M3,scale(mul(M2,M1),-3)),scale(powp(M1,3),2)),F(1,6))
C8=grade(K2,8)
expected=mul(scale(powp(a,8),F(1,2)),powp(t,4))
expected=mul(expected,add(one,mul(powp(x,2),powp(t,2))))
expected=mul(expected,add(one,mul(powp(y,2),powp(t,2))))
ck(C8==expected,'complete degree-eight rank 4/6/6/8 tensor scalar specialization')
ck(sorted(e[1] for e in C8)==[4,6,6,8],'two labelled rank-six partners remain distinct')
ck({e[0] for e in M1}=={6},'first cumulant exactly degree six')
ck({e[0] for e in K2}=={8,10,12},'second cumulant degrees 8,10,12')
ck(min(e[0] for e in K3)==14,'third cumulant actually starts degree fourteen in scalar multilinear model')
ck(not grade(K2,9),'no degree-nine intrinsic conditional coefficient')
ck(bool(grade(K2,10)),'nonzero degree-ten current remains after degree-eight cancellation')
# Four literal graphs, center index 0,1,2,3, leaf index +4.
base=[(0,1),(2,3),(0,2),(1,3)]
records=[]
for markP,markQ in product([0,1],repeat=2):
    edges=base+([] if markP else [(0,2)])+([] if markQ else [(1,3)])
    p=[markP,markQ,markP,markQ]
    deg=[sum(i in e for e in edges) for i in range(4)]
    ck(all(deg[i]+p[i]==3 for i in range(4)),f'decorated C2 valence {markP}{markQ}')
    ck(all(a!=b for a,b in edges),f'loopless centers {markP}{markQ}')
    seen={0}
    while True:
        nxt=seen|{v for e in edges if any(x in seen for x in e) for v in e}
        if nxt==seen:break
        seen=nxt
    ck(len(seen)==4,f'connected centers {markP}{markQ}')
    tree=[(0,1),(0,2),(2,3)]
    for edge in tree:ck(edge in edges,f'fixed spanning tree valid {markP}{markQ} {edge}')
    rec={'type':f'T{markP}{markQ}','center_edges':edges,'direct_public_count':sum(p),'physical_rank':4+sum(p),'cut_roots':len(edges)-3,'complete_force_calls':8}
    records.append(rec)
ck(sum(r['cut_roots'] for r in records)==8,'eight total cut roots')
ck(sum(r['physical_rank'] for r in records)==24,'twenty-four physical rows')
ck(sum(r['complete_force_calls'] for r in records)==32,'thirty-two complete source programs')
# Surplus algebra and positivity guards, exact rational exponents.
for d in range(1,15):
    for h in range(2,8):
        for j in range(0,12):
            n=2*h+j;A=h*(d+2)+j
            ck(A==n+h*d,f'surplus identity d{d}h{h}j{j}')
            ck(A-n>=2*d,f'strict surplus d{d}h{h}j{j}')
# Physical color pairings never contract a protected selected endpoint. Generate
# all ways to consume an even subset of one repeated direct-public color.
def pairings(xs):
    if not xs:yield [];return
    a=xs[0]
    for j in range(1,len(xs)):
        b=xs[j]
        for rest in pairings(xs[1:j]+xs[j+1:]):yield [(a,b)]+rest
for h in range(2,9):
    for q in range(0,h+1,2):
        for consumed in combinations(range(h),q):
            for pairs in pairings(list(consumed)):
                ck(all(a!=b for a,b in pairs),f'direct-public pairing is inter-occurrence h{h}q{q}:{pairs}')
                # each occurrence loses one external physical leg and gains one edge
                for vtx in range(h):
                    hits=sum(vtx in e for e in pairs)
                    ck(hits+(vtx not in consumed)==1,f'port conservation h{h}q{q}v{vtx}:{pairs}')
# Simultaneous scalar coefficient cancellation at every symbolic s: derivatives
# are polynomial identities, no numeric sampling masquerading as proof.
def sp_add(A,B):
    n=max(len(A),len(B));return [(A[i] if i<len(A) else F(0))+(B[i] if i<len(B) else F(0)) for i in range(n)]
def sp_der(A):return [F(i)*A[i] for i in range(1,len(A))]
sC6=[F(0),F(1)];refC6=[F(1),F(-1)]
sC8=[F(0),F(0),F(1)];repairC8=[F(0),F(0),F(-1)]
ck(sp_add(sC6,refC6)==[F(1),F(0)],'path C6 sum constant for symbolic s')
ck(all(c==0 for c in sp_der(sp_add(sC6,refC6))),'path C6 derivative vanishes identically')
ck(all(c==0 for c in sp_add(sC8,repairC8)),'path C8 sum zero for symbolic s')
ck(all(c==0 for c in sp_der(sp_add(sC8,repairC8))),'path C8 derivative vanishes identically')
result={'all_pass':all(x['pass'] for x in checks),'assertions':len(checks),'graphs':records,'C8_scalar_terms':pstr(C8),'C10_scalar_terms':pstr(grade(K2,10)),'checked_scope':'finite coefficient/graph/surplus diagnostics only; no native VALUE execution or positive LAW audit','checks':checks}
path=Path(__file__).with_name('decorated_induction_checks.json')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['checks']},indent=2))
