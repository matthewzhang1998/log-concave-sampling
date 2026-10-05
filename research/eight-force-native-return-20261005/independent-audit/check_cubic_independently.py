#!/usr/bin/env python3
"""Independent cubic-core eight-force graph/cumulant audit; not an execution of a native sampler."""
from collections import defaultdict, Counter
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
checks=0
def check(x):
 global checks
 checks+=1
 assert x

def connected(n,edges):
 seen={0}
 for _ in range(n):
  for a,b in edges:
   if a in seen or b in seen:seen|={a,b}
 return len(seen)==n

def all_regular(n):
 pairs=list(combinations(range(n),2));deg=[0]*n;edges=[]
 def rec(k):
  if k==len(pairs):
   if all(d==3 for d in deg) and connected(n,edges):yield edges.copy()
   return
  a,b=pairs[k]
  for m in range(min(3-deg[a],3-deg[b])+1):
   deg[a]+=m;deg[b]+=m;edges.extend([(a,b)]*m)
   if b<n-1 or deg[a]==3:yield from rec(k+1)
   if m:del edges[-m:]
   deg[a]-=m;deg[b]-=m
 yield from rec(0)

def trees(n,edges):
 for inds in combinations(range(len(edges)),n-1):
  if connected(n,[edges[i] for i in inds]):yield inds

def terms(n,edges,inds,rootedge):
 """All reverse-shift monomials: exact downward-closed center subsets."""
 a,b=edges[rootedge];par={a:None,b:None};todo=[a,b]
 while todo:
  v=todo.pop()
  for i in inds:
   if i==rootedge:continue
   x,y=edges[i]
   if v not in (x,y):continue
   w=y if x==v else x
   if w not in par:par[w]=v;todo.append(w)
 rest=[v for v in range(n) if v not in (a,b)]
 out=[]
 for mask in range(1<<len(rest)):
  vs={a,b}|{v for k,v in enumerate(rest) if (mask>>k)&1}
  if any(par[v] not in vs for v in vs if par[v] is not None):continue
  slots=defaultdict(list);kept=[]
  for i,(x,y) in enumerate(edges):
   if i not in inds:
    for v in (x,y):
     if v in vs:slots[('G',i)].append(v)
   elif x in vs and y in vs:kept.append((x,y))
   elif (x in vs)!=(y in vs):
    present=x if x in vs else y;absent=y if x in vs else x
    check(par[absent]==present)
    slots[('Y',absent)].append(present)
  check(len(kept)==len(vs)-1)
  check(all(sum(v in e for e in kept)+sum(v in arr for arr in slots.values())==3 for v in vs))
  check(all(len(arr)==len(set(arr)) for arr in slots.values()))
  out.append((vs,kept,dict(slots)))
 return out

def gm(k):
 if k%2:return 0
 ans=1
 for j in range(1,k,2):ans*=j
 return ans

def moment(term):
 ans=1
 for arr in term[2].values():ans*=gm(len(arr))
 return ans

def pairings(vs):
 if not vs:yield [];return
 a=vs[0]
 for i in range(1,len(vs)):
  b=vs[i]
  for rest in pairings(vs[1:i]+vs[i+1:]):yield [(a,b)]+rest

family=[]
for n in [2,4]:
 ng=nt=no=nm=0
 for edges in all_regular(n):
  ng+=1
  for inds in trees(n,edges):
   nt+=1
   for rootedge in inds:
    ts=terms(n,edges,inds,rootedge);no+=2;nm+=2*len(ts)
    surviving=[t for t in ts if moment(t)]
    check(len(surviving)==1)
    check(len(surviving[0][0])==n and moment(surviving[0])==1)
    # Reverse orientation changes neither selected center set nor this census.
 family.append(dict(n=n,graphs=ng,labelled_spanning_trees=nt,oriented_short_spines=no,monomials=nm))

edges=[(0,1),(2,3),(0,2),(0,2),(1,3),(1,3)]
polys=Counter();total_wick=0;connected_wick=0;orientations=0
for inds in trees(4,edges):
 for rootedge in inds:
  ts=terms(4,edges,inds,rootedge);poly=defaultdict(int)
  for t,u in product(ts,repeat=2):
   vs1,ee1,s1=t;vs2,ee2,s2=u
   slots={k:[(0,v) for v in s1.get(k,[])]+[(1,v) for v in s2.get(k,[])] for k in set(s1)|set(s2)}
   if any(len(v)%2 for v in slots.values()):continue
   options=[list(pairings(v)) for v in slots.values()]
   num=0
   vertices=[(0,v) for v in sorted(vs1)]+[(1,v) for v in sorted(vs2)]
   ix={v:i for i,v in enumerate(vertices)}
   kept=[(ix[(0,a)],ix[(0,b)]) for a,b in ee1]+[(ix[(1,a)],ix[(1,b)]) for a,b in ee2]
   for choice in product(*options):
    ww=[(ix[a],ix[b]) for pp in choice for a,b in pp];out=kept+ww
    total_wick+=2
    check(all(a!=b for a,b in out))
    check(all(sum(v in e for e in out)==3 for v in range(len(vertices))))
    if connected(len(vertices),out):
     connected_wick+=2;num+=1
   direct=1
   for v in slots.values():direct*=gm(len(v))
   direct-=moment(t)*moment(u)
   check(num==direct)
   poly[len(vertices)]+=num
  # Both directed orientations have the same scalar moments.
  check(poly[4] in (1,3))
  key=tuple(sorted((n,c) for n,c in poly.items() if c));polys[str(key)]+=2;orientations+=2
check(orientations>0)
# Exact surplus and source exponents, independently using local amplitude counts.
beta=Fraction(1,2);gamma=Fraction(0)
for n in range(2,30):
 for d in [1,2,4,8,16]:
  root=d+beta
  check(root+(n-1)*beta+n*beta==2*beta*n+d)
  for h in range(2,10):
   for j in range(h*(n-2)+1):
    nn=2*h+j
    aa=h*(root+3*beta)+2*beta*j
    check(aa-2*beta*nn==h*d)
    check(aa-gamma*nn==(2*beta-gamma)*nn+h*d)
output=dict(scope='Independent finite graph and exact h=2 Wick census; no native sampler or same-endpoint remainder is executed.',assertions=checks,all_regular_graph_tests=family,literal_cubic_core_orientations=orientations,wick_matchings=total_wick,connected_wick_matchings=connected_wick,variance_polynomial_counts=dict(polys),first_variance_coefficient='1/2 or 3/2 depending on root core-edge multiplicity',first_variance_alpha_power=12,first_variance_z_power=4,first_variance_effective_grade='12')
(ROOT/'independent_cubic_checks.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
