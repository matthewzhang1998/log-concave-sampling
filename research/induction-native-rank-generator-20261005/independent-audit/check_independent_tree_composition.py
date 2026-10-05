"""Independent tree/all-cut and restricted multigraph/HS diagnostics.
Uses recursive tree growth, distinct from the author's Prüfer enumeration.
"""
from pathlib import Path
import itertools,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(501026)
count=0;flowcases=0;densecases=0;multicases=0

def ck(ok,msg):
 global count
 count+=1
 if not bool(ok):raise AssertionError(msg)

def certify_flow(n,edges,marks,mask):
 global flowcases
 I={j for j in range(len(marks)) if (mask>>j)&1};O=set(range(len(marks)))-I
 w=[2**(2**j) for j in range(len(marks))]
 si=sum(w[j] for j in I);so=sum(w[j] for j in O)
 mass=[w[j]*so if j in I else -w[j]*si for j in range(len(marks))]
 ck(sum(mass)==0,'balanced flow')
 degree=np.zeros((n,2),dtype=int)
 for j,v in enumerate(marks):degree[v,0 if j in I else 1]+=1
 for a,b in edges:
  remaining=[e for e in edges if e!=(a,b)]
  component={a}
  while True:
   expanded=component|{v for u,v in remaining if u in component}|{u for u,v in remaining if v in component}
   if expanded==component:break
   component=expanded
  left={j for j,v in enumerate(marks) if v in component}
  right=set(range(len(marks)))-left
  ck(bool(left) and bool(right),'every edge separates terminals')
  # Exact crossing pair coefficient set is independently nonempty and unique.
  plus={(i,j) for i in I&left for j in O&right}
  minus={(i,j) for i in I&right for j in O&left}
  ck(bool(plus|minus) and not(plus&minus),'nonzero edge flow polynomial')
  exponents=[2**i+2**j for i,j in plus|minus]
  ck(len(exponents)==len(set(exponents)),'binary exponent uniqueness')
  flow=sum(mass[j] for j in left)
  ck(flow!=0,'integer specialization has no cancellation')
  u,v=(a,b) if flow>0 else (b,a)
  degree[u,1]+=1;degree[v,0]+=1
 ck(np.all(degree>0),'each local map is a proper flattening')
 flowcases+=1

for n in range(2,9):
 for rep in range(12):
  edges=[(int(rng.integers(v)),v) for v in range(1,n)]
  deg=[sum(v in e for e in edges) for v in range(n)]
  marks=[v for v in range(n) if deg[v]==1]
  if rep%3==0:marks.append(int(rng.integers(n)))
  for mask in range(1,2**len(marks)-1):certify_flow(n,edges,marks,mask)

# Random tensors with heterogeneous slot dimensions, normalized by every cut.
def maxcut(a):
 d=a.ndim;best=0.
 for mask in range(1,2**d-1):
  l=[i for i in range(d) if mask>>i&1];r=[i for i in range(d) if not(mask>>i&1)]
  matrix=a.transpose(l+r).reshape(math.prod(a.shape[i] for i in l),-1)
  best=max(best,float(np.linalg.norm(matrix,2)))
 return best

def network(n,edges,marks):
 labels=[[] for _ in range(n)];dims={}
 for j,(u,v) in enumerate(edges):
  labels[u].append(j);labels[v].append(j);dims[j]=2+(j%2)
 out=[]
 for j,v in enumerate(marks):
  lab=len(edges)+j;labels[v].append(lab);out.append(lab);dims[lab]=2
 arrays=[];hss=[]
 for labs in labels:
  arr=rng.normal(size=tuple(dims[j] for j in labs));arr/=maxcut(arr)
  arrays.append(arr);hss.append(float(np.linalg.norm(arr.ravel())))
 args=[]
 for a,labs in zip(arrays,labels):args.extend((a,labs))
 return np.einsum(*args,out,optimize=True),hss

for n in range(2,7):
 for rep in range(9):
  edges=[(int(rng.integers(v)),v) for v in range(1,n)]
  deg=[sum(v in e for e in edges) for v in range(n)]
  marks=[v for v in range(n) if deg[v]==1]
  if rep%2:marks.append(int(rng.integers(n)))
  T,hs=network(n,edges,marks)
  ck(maxcut(T)<=1+2e-11,'tree product of all-cut constants')
  for h in hs:ck(np.linalg.norm(T.ravel())<=h+2e-11,'tree one-HS product with arbitrary anchor')
  densecases+=1
  # Retain original spanning tree; its marked leaves remain marked.
  extra=[tuple(rng.choice(n,2,replace=False)) for _ in range(rep%3+1)]
  M,mhs=network(n,edges+extra,marks)
  root=next(v for v,d in enumerate(deg) if d==1)
  ck(np.linalg.norm(M.ravel())<=mhs[root]+2e-11,'loopless multigraph one-HS with marked spanning leaf')
  multicases+=1
out={'status':'PASS','assertions':count,'flow_cases':flowcases,'heterogeneous_dense_tree_networks':densecases,'restricted_multigraph_HS_networks':multicases,'scope':'Supports Section 8 only. No all-cut theorem for cyclic graphs, self-traces, or native return.'}
(HERE/'independent-tree-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
