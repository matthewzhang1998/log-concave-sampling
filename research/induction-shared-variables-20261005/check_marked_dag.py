#!/usr/bin/env python3
from itertools import product,combinations,permutations
from fractions import Fraction as F
from collections import deque
from pathlib import Path
import numpy as np, json, math
count=0
summaries={}
def ck(x,msg):
 global count
 count+=1
 if not x:raise AssertionError(msg)
def components(n,edges):
 adj=[[] for _ in range(n)]
 for u,v in edges:adj[u].append(v);adj[v].append(u)
 out=[];left=set(range(n))
 while left:
  seen={min(left)};q=list(seen)
  for u in q:
   for v in adj[u]:
    if v not in seen:seen.add(v);q.append(v)
  out.append(seen);left-=seen
 return out
def marked_tree(n,es,marked):
 for ids in combinations(range(len(es)),n-1):
  T=[es[i] for i in ids]
  if len(components(n,T))!=1:continue
  deg=[sum(v in e for e in T) for v in range(n)]
  if all(deg[v]!=1 or v in marked for v in range(n)):return ids
 return None
def orient(n,edges,tids,terminals,inputs):
 # Explicit superincreasing positive masses; normalization balances totals.
 ws=[1<<(1<<i) for i in range(len(terminals))]
 si=sum(ws[i] for i in inputs);so=sum(ws[i] for i in range(len(ws)) if i not in inputs)
 loads=[F(0)]*n
 for i,v in enumerate(terminals):loads[v]+=F(ws[i],si) if i in inputs else -F(ws[i],so)
 ck(sum(loads)==0,'balanced explicit terminal masses')
 oriented=[]
 T=[edges[i] for i in tids]
 for z,(u,v) in enumerate(T):
  cs=components(n,T[:z]+T[z+1:]);S=next(x for x in cs if u in x)
  f=sum(loads[x] for x in S)
  ck(bool(f),'nonzero flow on every tree edge')
  oriented.append((u,v) if f>0 else (v,u))
 adj=[[] for _ in range(n)];ind=[0]*n
 for u,v in oriented:adj[u].append(v);ind[v]+=1
 q=deque(v for v in range(n) if ind[v]==0);order=[]
 while q:
  u=q.popleft();order.append(u)
  for v in adj[u]:
   ind[v]-=1
   if ind[v]==0:q.append(v)
 ck(len(order)==n,'directed tree topological order')
 pos={v:i for i,v in enumerate(order)}
 full=[]
 for i,(u,v) in enumerate(edges):
  if i in tids:full.append(oriented[tids.index(i)])
  else:full.append((u,v) if pos[u]<pos[v] else (v,u))
 ck(all(pos[u]<pos[v] for u,v in full),'all graph edges advance in one topological order')
 for v in range(n):
  ie=sum(b==v for a,b in full)+sum(terminals[i]==v for i in inputs)
  oe=sum(a==v for a,b in full)+sum(terminals[i]==v for i in range(len(terminals)) if i not in inputs)
  ck(ie>0 and oe>0,'every local tensor has a proper directed cut')
  if not any(b==v for a,b in full):ck(any(terminals[i]==v for i in inputs),'DAG source has external input')
  if not any(a==v for a,b in full):ck(any(terminals[i]==v for i in range(len(terminals)) if i not in inputs),'DAG sink has external output')
 # Each edge is a distinct wire; parallel edges keep distinct positions.
 ck(len(full)==len(edges),'one directed wire per labelled original edge')
 return full,order
ngraph=ncert=ncuts=0
for n in range(2,5):
 pairs=list(combinations(range(n),2))
 for mult in product(range(3),repeat=len(pairs)):
  es=[e for e,m in zip(pairs,mult) for _ in range(m)]
  if len(es)<n-1 or any(sum(v in e for e in es)>4 for v in range(n)):continue
  if len(components(n,es))>1:continue
  ngraph+=1
  for bits in product([0,1],repeat=n):
   marked={i for i,b in enumerate(bits) if b}
   if len(marked)<2:continue
   tids=marked_tree(n,es,marked)
   if tids is None:continue
   ncert+=1;terms=sorted(marked)
   for b in range(1,(1<<len(terms))-1):
    ins={i for i in range(len(terms)) if b>>i&1}
    orient(n,es,tids,terms,ins);ncuts+=1
summaries.update(multigraphs=ngraph,marked_certificates=ncert,proper_boundary_cuts=ncuts)
# Direct tensor-network SVD tests for parallel and triangle graphs.
def allcuts(T):
 r=T.ndim
 vals=[]
 for bits in range(1,(1<<r)-1):
  a=[i for i in range(r) if bits>>i&1];b=[i for i in range(r) if not bits>>i&1]
  M=T.transpose(a+b).reshape(math.prod(T.shape[i] for i in a),-1)
  vals.append(np.linalg.svd(M,compute_uv=False)[0])
 return max(vals)
def contract(edges,marked,D,tensors):
 # label every edge and boundary leg once at each endpoint / once at boundary
 labs=[[] for _ in tensors]
 for j,(u,v) in enumerate(edges):labs[u].append(j);labs[v].append(j)
 out=[]
 for j,v in enumerate(marked):labs[v].append(len(edges)+j);out.append(len(edges)+j)
 args=[]
 for T,ls in zip(tensors,labs):args.extend([T,ls])
 return np.einsum(*args,out,optimize=True)
rng=np.random.default_rng(79103)
for name,n,es,marked in [('parallel',2,[(0,1),(0,1)],[0,1]),('triangle',3,[(0,1),(1,2),(2,0)],[0,1])]:
 for D in [2,3,4]:
  for trial in range(20):
   tensors=[];prod_bound=1.0
   for v in range(n):
    r=sum(v in e for e in es)+marked.count(v)
    T=rng.normal(size=(D,)*r)
    T=sum(T.transpose(p) for p in permutations(range(r)))/math.factorial(r)
    norm=allcuts(T);T/=norm;tensors.append(T);prod_bound*=allcuts(T)
   O=contract(es,marked,D,tensors)
   ck(allcuts(O)<=prod_bound+1e-10,f'{name} SVD allcut bound')
  # Saturating copy tensors: no hidden factor D appears.
  ts=[]
  for v in range(n):
   r=sum(v in e for e in es)+marked.count(v);T=np.zeros((D,)*r)
   for a in range(D):T[(a,)*r]=1
   ts.append(T)
  O=contract(es,marked,D,ts)
  ck(np.allclose(O,np.eye(D)),f'{name} identity contraction, no hidden trace')
# Restriction witness: unmarked hanging triangle with both physical legs at A.
for D in [2,3,5,8]:
 es=[(0,1),(1,2),(2,0)];marked=[0,0]
 A=np.zeros((D,D,D,D))
 for i in range(D):A[i,i,0,0]=1/np.sqrt(D)
 ts=[A,np.eye(D),np.eye(D)]
 ck(allcuts(A)<=1+1e-12,'unmarked-triangle local cuts bounded')
 O=contract(es,marked,D,ts)
 ck(abs(allcuts(O)-np.sqrt(D))<1e-10,'unmarked-triangle proper-cut blowup witness')
 ck(marked_tree(3,es,{0}) is None,'witness fails marked-spanning-tree hypothesis')
# First true marked-leaf bank bridge keeps all physical legs.
es=[(0,1),(1,2),(3,4),(4,5),(0,3)]
ck(marked_tree(6,es,set(range(6))) is not None,'leaf-to-leaf C1 bank bridge has marked spanning tree')
ck([sum(v in e for e in es)+1-2 for v in range(6)]==[1,1,0,1,1,0],'exact C1-hit native source orders')
result={'all_pass':True,'assertions':count,**summaries,'scope':'finite graph/circuit and direct tensor-cut tests; does not certify full old-bank heat or mixed queue'}
Path(__file__).with_name('marked_dag_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
