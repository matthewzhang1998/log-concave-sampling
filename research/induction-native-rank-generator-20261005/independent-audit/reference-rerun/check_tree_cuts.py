"""Exact boundary-flow certificates and independent dense tensor tests."""
from pathlib import Path
import itertools,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
checks=0

def ck(x,m):
 global checks
 checks+=1
 if not x:raise AssertionError(m)

def prufer_tree(word):
 n=len(word)+2;degrees=[1]*n
 for v in word:degrees[v]+=1
 edges=[]
 for v in word:
  leaf=next(i for i,d in enumerate(degrees) if d==1)
  edges.append((leaf,v));degrees[leaf]-=1;degrees[v]-=1
 u,v=[i for i,d in enumerate(degrees) if d==1];edges.append((u,v))
 return edges

def flow_certificate(n,edges,terminals,inputs):
 adj=[[] for _ in range(n)]
 for i,(u,v) in enumerate(edges):adj[u].append((v,i));adj[v].append((u,i))
 weights=[1<<(1<<j) for j in range(len(terminals))]
 si=sum(weights[j] for j in inputs);so=sum(weights[j] for j in range(len(terminals)) if j not in inputs)
 masses=[w*so if j in inputs else -w*si for j,w in enumerate(weights)]
 ck(sum(masses)==0,'balanced terminal flow')
 outdeg=[0]*n;indeg=[0]*n
 for j,v in enumerate(terminals):
  (indeg if j in inputs else outdeg)[v]+=1
 oriented=[]
 for ei,(u,v) in enumerate(edges):
  stack=[v];component={v}
  while stack:
   a=stack.pop()
   for b,i in adj[a]:
    if i==ei or b in component:continue
    component.add(b);stack.append(b)
  flow=sum(m for j,m in enumerate(masses) if terminals[j] in component)
  ck(flow!=0,'every edge nonzero')
  src,dst=(v,u) if flow>0 else (u,v)
  outdeg[src]+=1;indeg[dst]+=1;oriented.append((src,dst))
 ck(all(a>0 and b>0 for a,b in zip(indeg,outdeg)),'every local cut proper')
 # Undirected tree implies DAG; check directly as well.
 incoming=[0]*n
 for u,v in oriented:incoming[v]+=1
 queue=[i for i in range(n) if incoming[i]==0];seen=0
 while queue:
  u=queue.pop();seen+=1
  for a,b in oriented:
   if a==u:
    incoming[b]-=1
    if incoming[b]==0:queue.append(b)
 ck(seen==n,'acyclic operator circuit')

flow_trees=0
for n in range(2,7):
 for word in itertools.product(range(n),repeat=n-2):
  edges=prufer_tree(word);degree=[0]*n
  for u,v in edges:degree[u]+=1;degree[v]+=1
  terminals=[i for i,d in enumerate(degree) if d==1]
  # All leaves marked; periodically add internal external slots too.
  variants=[terminals]
  if n<=4:variants.append(terminals+list(range(n)))
  for terms in variants:
   for mask in range(1,(1<<len(terms))-1):
    inputs={j for j in range(len(terms)) if mask>>j&1}
    flow_certificate(n,edges,terms,inputs)
  flow_trees+=1

rng=np.random.default_rng(17402)
dense_tests=0
for n in range(2,7):
 for rep in range(12):
  word=tuple(rng.integers(0,n,size=n-2));edges=prufer_tree(word)
  local=[[] for _ in range(n)]
  for i,(u,v) in enumerate(edges):local[u].append(i);local[v].append(i)
  marks=[v for v in range(n) if len(local[v])==1]
  if rep%3==0:marks.append(int(rng.integers(n)))
  boundary=[]
  for j,v in enumerate(marks):
   lab=n-1+j;local[v].append(lab);boundary.append(lab)
  arrays=[];hss=[]
  for slots in local:
   order=len(slots);arr=rng.normal(size=(2,)*order)
   top=0.
   for mask in range(1,(1<<order)-1):
    left=[i for i in range(order) if mask>>i&1];right=[i for i in range(order) if not mask>>i&1]
    top=max(top,np.linalg.norm(arr.transpose(left+right).reshape(2**len(left),-1),ord=2))
   arr/=top;arrays.append(arr);hss.append(np.linalg.norm(arr.ravel()))
  operands=[]
  for arr,slots in zip(arrays,local):operands.extend((arr,slots))
  network=np.einsum(*operands,boundary,optimize=True);m=len(boundary)
  for mask in range(1,(1<<m)-1):
   left=[i for i in range(m) if mask>>i&1];right=[i for i in range(m) if not mask>>i&1]
   op=np.linalg.norm(network.transpose(left+right).reshape(2**len(left),-1),ord=2)
   ck(op<=1+1e-10,'composed tree proper-cut product')
  for hs in hss:ck(np.linalg.norm(network.ravel())<=hs+1e-10,'one local HS')
  dense_tests+=1

out={'status':'PASS','assertions':checks,'all_labelled_trees_through_vertices':6,'flow_trees':flow_trees,'dense_networks':dense_tests,'scope':'Tree operator circuits only. No full cyclic all-cut or native comparison theorem.'}
(HERE/'tree-cut-checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
