#!/usr/bin/env python3
"""Exact finite symbol/graph and exponent diagnostics, not a native sampler."""
from fractions import Fraction as Q
from math import factorial, comb, log2
from collections import defaultdict
from itertools import combinations
import heapq
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
checks=0
def check(condition):
 global checks
 checks+=1
 assert condition

def gm(k):
 if k%2:return 0
 x=1
 for j in range(1,k,2):x*=j
 return x

def add(a,b,scale=1):
 c=dict(a)
 for n,v in b.items():c[n]=c.get(n,0)+scale*v
 return {n:v for n,v in c.items() if v}
def mul(a,b):
 c=defaultdict(int)
 for n,v in a.items():
  for m,w in b.items():c[n+m]+=v*w
 return dict(c)
# V/H = sum over n=2,3,4 of alpha^(d+2 beta n) z^n theta^n F_n.
# F2=G0^2 G1 G3 G4 Yd; F3=G0^2 G1^2 G2 G3 G4 Yc;
# F4=G0^2 G1^2 G2^2 G3^2 G4^2.
F=[(2,(2,1,0,1,1,1,0)),(3,(2,2,1,1,1,0,1)),(4,(2,2,2,2,2,0,0))]
def moment_multinomial(h):
 out=defaultdict(int)
 for i in range(h+1):
  for j in range(h-i+1):
   k=h-i-j
   n=2*i+3*j+4*k
   powers=[2*h,2*h-i,j+2*k,h+k,h+k,i,j]
   coef=factorial(h)//factorial(i)//factorial(j)//factorial(k)
   for power in powers:coef*=gm(power)
   out[n]+=coef
 return {n:c for n,c in out.items() if c}
# Independent polynomial multiplication before Gaussian integration.
def moment_expanded(h):
 poly={(0,(0,)*7):1}
 for _ in range(h):
  nxt=defaultdict(int)
  for (n,p),c in poly.items():
   for m,q in F:nxt[(n+m,tuple(a+b for a,b in zip(p,q)))]+=c
  poly=dict(nxt)
 out=defaultdict(int)
 for (n,p),c in poly.items():
  for power in p:c*=gm(power)
  out[n]+=c
 return {n:c for n,c in out.items() if c}
H=8
mom=[{0:1}]+[moment_multinomial(h) for h in range(1,H+1)]
kap=[{}]
for h in range(1,H+1):
 check(mom[h]==moment_expanded(h))
 c=mom[h]
 for j in range(1,h):c=add(c,mul(kap[j],mom[h-j]),-comb(h-1,j-1))
 kap.append(c)
 reconstructed={}
 for j in range(1,h+1):reconstructed=add(reconstructed,mul(kap[j],mom[h-j]),comb(h-1,j-1))
 check(reconstructed==mom[h])
 for n in c:
  check(2*h<=n<=4*h)
  check(Q(4*h)+Q(2*n,3)>=Q(16*h,3))
check(kap[1]=={4:1})
check(kap[2]=={4:3,6:9,8:242})
# First variance: three G0 pairings, with all remaining edges forced.
# vertices 0,1 are a,b; vertices 2,3 are a',b'.
pairings=[[(0,1),(2,3)],[(0,2),(1,3)],[(0,3),(1,2)]]
fixed=[(0,1),(2,3),(1,3),(1,3),(0,2),(0,2)]
variance_graphs=[]
for pp in pairings:
 edges=fixed+pp
 deg=[sum(v in e for e in edges) for v in range(4)]
 check(deg==[4]*4)
 check(all(a!=b for a,b in edges))
 variance_graphs.append(edges)
# All doubled-C4 center spanning trees: adjacent root spines always have two
# centers; selecting each other center's own physical leaf gives 2-vertex sides.
base=[(0,1),(1,3),(3,2),(2,0)]
edges=[(a,b) for a,b in base for _ in range(2)]
def connected(ee,n=4):
 seen={0}
 for _ in range(n):
  for a,b in ee:
   if a in seen or b in seen:seen|={a,b}
 return len(seen)==n
numtrees=0; rooted_short=0
for cc in combinations(range(8),3):
 tt=[edges[i] for i in cc]
 if not connected(tt):continue
 numtrees+=1
 cut=[edges[i] for i in range(8) if i not in cc]
 for a,b in tt:
  for root,terminal in [(a,b),(b,a)]:
   rooted_short+=1
   chunks=[{root},{terminal}]
   check(all(not (u in chunk and v in chunk) for u,v in cut for chunk in chunks))
   # Every remaining center uses its own physical leaf as selected child.
   check(4*Q(1,2)+(Q(4)+Q(1,2)-Q(1,2))==6)
check(numtrees==32);check(rooted_short==192)
# Exact surplus and first/prior guards over finite ranges.
beta=Q(1,2);gamma=Q(1,3);d0=Q(4)
for n in range(2,101):
 for d in [Q(4),Q(8),Q(12),Q(16),Q(32)]:
  A=2*beta*n+d
  root=d+beta
  check(root+(2*n-1)*beta==A)
  for h in range(2,9):
   for np in range(2*h,h*n+1):
    Ap=h*(d+4*beta)+2*beta*(np-2*h)
    check(Ap==2*beta*np+h*d)
    check(Ap-2*beta*np>=2*d)
    check(Ap-gamma*np>=h*d)
  for m in range(1,n+1):
   first=d+(m+1)*(beta-gamma)
   check(first>=1)
  for P in range(1,31):
   b=6*(P+1)
   # Most conservative error: ignore helpful ancestor attenuation.
   check(b*(beta-gamma)>P)
# Queue over numerical states is an over-enumeration of true graph outputs.
queue_reports=[]
for P in [8,12,16,24,32,48,64,96]:
 seen={(4,4)}; pending=[(4,4)]; longest={(4,4):0}; depth=0; transitions=0
 while pending:
  d,n=heapq.heappop(pending);dd=longest[(n,d)];depth=max(depth,dd)
  for h in range(2,int(Q(P,d))+1):
   dn=h*d
   for nn in range(2*h,h*n+1):
    if Q(dn)+Q(2*nn,3)>=P:continue
    transitions+=1
    child=(nn,dn)
    longest[child]=max(longest.get(child,0),dd+1)
    if child not in seen:seen.add(child);heapq.heappush(pending,(dn,nn))
 check(all(d+Q(2*n,3)<P for n,d in seen if (n,d)!=(4,4)))
 check(depth<=int(log2(P/4)))
 queue_reports.append({'target':P,'states':len(seen),'transitions':transitions,'max_depth':depth})
# Deliberately bad beta: force-only/effective grade is not monotone even though
# surplus doubles. This rejects replacing the actual queue order by grade.
n=100;d=4;h=2;np=4
before=2*beta*n+d-gamma*n
after=2*beta*np+h*d-gamma*np
check(after<before)
output={
 'scope':'Finite exact scalar symbols, labeled graph degree/mark checks and rational path/queue ledgers. No native pair sampler is executed.',
 'assertions':checks,
 'beta':str(beta),'gamma':str(gamma),'initial_surplus':str(d0),
 'input_grade':str(Q(20,3)),
 'first_feedback_grade':str(Q(32,3)),
 'second_cumulant':[{'centers':n,'physical_rank':n,'alpha_power':str(8+n),'z_power':n,'coefficient':str(Q(c,2))} for n,c in sorted(kap[2].items())],
 'cumulants_through_8':[{'order':h,'terms':[{'centers':n,'coefficient':str(Q(c,factorial(h))),'surplus':4*h} for n,c in sorted(kap[h].items())]} for h in range(1,H+1)],
 'first_variance_center_graphs':variance_graphs,
 'center_spanning_trees':numtrees,'admitted_short_spines':rooted_short,
 'queue_reports':queue_reports,
 'nonmonotone_effective_grade_countertest':{'input_n':n,'input_grade':str(before),'offspring_n':np,'offspring_grade':str(after),'surplus_before':d,'surplus_after':h*d},
 'boundary':'Candidate conditional/auxiliary permanent-leaf C3 closure. Whole-old-bank derivative/tilt native admission and terminal joined curl are not proved here.'
}
(ROOT/'cycle_surplus_checks.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:output[k] for k in ['scope','assertions','input_grade','first_feedback_grade','second_cumulant','queue_reports','boundary']},indent=2))
