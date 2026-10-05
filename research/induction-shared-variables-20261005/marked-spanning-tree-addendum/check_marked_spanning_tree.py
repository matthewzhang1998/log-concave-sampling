#!/usr/bin/env python3
"""Exact combinatorics and numerical index identities; no native compiler execution."""
from pathlib import Path
from itertools import product,combinations
from collections import Counter
import math,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
checks=0

def check(ok,label):
 global checks
 checks+=1
 if not ok:raise AssertionError(label)

def prufer_tree(word,n):
 degrees=[1]*n
 for x in word:degrees[x]+=1
 edges=[]
 for x in word:
  leaf=next(v for v in range(n) if degrees[v]==1)
  edges.append((leaf,x));degrees[leaf]-=1;degrees[x]-=1
 a,b=[v for v in range(n) if degrees[v]==1];edges.append((a,b))
 return edges

def flow_orientation(n,edges,tree,marks,inputs):
 # marks is a list of terminal vertices; all terminals remain individually labelled.
 weights=[1<<(1<<j) for j in range(len(marks))]
 wi=sum(weights[j] for j in inputs);wo=sum(weights[j] for j in range(len(marks)) if j not in inputs)
 mass=[0]*n
 for j,v in enumerate(marks):mass[v]+=weights[j]*(wo if j in inputs else -wi)
 check(sum(mass)==0,'terminal balance')
 orientation={}
 for k in tree:
  a,b=edges[k];seen={a}
  while True:
   old=len(seen)
   for q in tree:
    if q==k:continue
    x,y=edges[q]
    if x in seen or y in seen:seen.update((x,y))
   if len(seen)==old:break
  check(any(v in seen for v in marks) and any(v not in seen for v in marks),'terminals on each tree side')
  flow=sum(mass[v] for v in seen)
  check(flow!=0,'nonzero edge flow')
  orientation[k]=(a,b) if flow>0 else (b,a)
 order=[];remaining=set(range(n))
 while remaining:
  ready=[v for v in remaining if not any(head==v and tail in remaining for tail,head in orientation.values())]
  check(bool(ready),'tree acyclic')
  v=min(ready);remaining.remove(v);order.append(v)
 position={v:k for k,v in enumerate(order)}
 for k,(a,b) in enumerate(edges):
  if k not in orientation:orientation[k]=(a,b) if position[a]<position[b] else (b,a)
 for v in range(n):
  nin=sum(head==v for tail,head in orientation.values())+sum(marks[j]==v for j in inputs)
  nout=sum(tail==v for tail,head in orientation.values())+sum(marks[j]==v for j in range(len(marks)) if j not in inputs)
  check(nin>0 and nout>0,'proper local cut')
 for tail,head in orientation.values():check(position[tail]<position[head],'full circuit acyclic')
 # Explicit wire census: each internal index produced once and consumed once.
 live={len(edges)+j for j in inputs};produced=Counter();consumed=Counter()
 for v in order:
  incoming=[k for k,(tail,head) in orientation.items() if head==v]+[len(edges)+j for j in inputs if marks[j]==v]
  outgoing=[k for k,(tail,head) in orientation.items() if tail==v]+[len(edges)+j for j in range(len(marks)) if j not in inputs and marks[j]==v]
  check(set(incoming)<=live,'all local inputs live')
  live.difference_update(incoming)
  check(not(set(outgoing)&live),'all local outputs fresh')
  live.update(outgoing);consumed.update(incoming);produced.update(outgoing)
 check(live=={len(edges)+j for j in range(len(marks)) if j not in inputs},'exact final boundary')
 for k in range(len(edges)):check(produced[k]==consumed[k]==1,'each internal edge used once each way')
 return order,orientation

census=[]
for n in range(2,7):
 trees=cuts=0
 for word in product(range(n),repeat=n-2):
  tree_edges=prufer_tree(word,n);trees+=1
  degree=Counter(v for e in tree_edges for v in e)
  leaves=[v for v in range(n) if degree[v]==1]
  # Deliberately add all missing simple edges and one extra parallel copy of each tree edge.
  present={tuple(sorted(e)) for e in tree_edges}
  edges=tree_edges+[e for e in combinations(range(n),2) if e not in present]+tree_edges
  for marks in {tuple(leaves),tuple(range(n))}:
   for mask in range(1,(1<<len(marks))-1):
    inputs={j for j in range(len(marks)) if mask&(1<<j)}
    flow_orientation(n,edges,list(range(n-1)),marks,inputs);cuts+=1
 census.append(dict(vertices=n,all_labelled_trees=trees,marked_boundary_cuts=cuts))

def max_local_cut(tensor):
 rank=tensor.ndim;best=0
 for mask in range(1,(1<<rank)-1):
  if not mask&1:continue
  left=[j for j in range(rank) if mask&(1<<j)];right=[j for j in range(rank) if not mask&(1<<j)]
  mat=tensor.transpose(left+right).reshape(2**len(left),2**len(right))
  best=max(best,np.linalg.norm(mat,2))
 return best

def circuit_tensor(local,slots,n,edges,marks,inputs,order,orientation):
 input_ids=[len(edges)+j for j in sorted(inputs)]
 # A boundary identity whose first axes are live wires and last axes are original global inputs.
 value=np.eye(2**len(input_ids)).reshape([2]*(2*len(input_ids)))
 live=input_ids[:]
 for v in order:
  incoming=[k for k,(tail,head) in orientation.items() if head==v]+[len(edges)+j for j in sorted(inputs) if marks[j]==v]
  outgoing=[k for k,(tail,head) in orientation.items() if tail==v]+[len(edges)+j for j in range(len(marks)) if j not in inputs and marks[j]==v]
  matrix=local[v].transpose([slots[v].index(k) for k in incoming+outgoing])
  value=np.tensordot(matrix,value,axes=(list(range(len(incoming))),[live.index(k) for k in incoming]))
  live=outgoing+[k for k in live if k not in incoming]
 output_ids=[len(edges)+j for j in range(len(marks)) if j not in inputs]
 value=value.transpose([live.index(k) for k in output_ids]+list(range(len(live),value.ndim)))
 return value.reshape(2**len(output_ids),2**len(input_ids))

rng=np.random.default_rng(20261005)
fixtures=[
 ('two_parallel_edges',2,[(0,1),(0,1)],[0],[0,1]),
 ('three_parallel_edges',2,[(0,1)]*3,[0],[0,1]),
 ('triangle_two_marks',3,[(0,1),(1,2),(0,2)],[0,1],[0,2]),
 ('triangle_three_marks',3,[(0,1),(1,2),(0,2)],[0,1],[0,1,2]),
 ('marked_C1_bridge',4,[(0,1),(2,3),(1,2)],[0,1,2],[0,1,2,3]),
 ('triangle_multi_external',3,[(0,1),(1,2),(0,2)],[0,1],[0,0,2,2])]
fixture_out=[]
for name,n,edges,tree,marks in fixtures:
 max_error=max_norm=0.0;ncuts=0
 for trial in range(12):
  slots=[];local=[]
  for v in range(n):
   ss=[k for k,e in enumerate(edges) if v in e]+[len(edges)+j for j,w in enumerate(marks) if w==v]
   tt=rng.normal(size=[2]*len(ss));tt/=max_local_cut(tt)
   check(max_local_cut(tt)<=1+1e-12,'local allproper normalization')
   slots.append(ss);local.append(tt)
  args=[]
  for tt,ss in zip(local,slots):args.extend([tt,ss])
  network=np.einsum(*args,[len(edges)+j for j in range(len(marks))],optimize=True)
  for mask in range(1,(1<<len(marks))-1):
   inputs={j for j in range(len(marks)) if mask&(1<<j)}
   order,orientation=flow_orientation(n,edges,tree,marks,inputs)
   mat=circuit_tensor(local,slots,n,edges,marks,inputs,order,orientation)
   output=sorted(set(range(len(marks)))-inputs)
   direct=network.transpose(output+sorted(inputs)).reshape(mat.shape)
   error=np.linalg.norm(mat-direct);norm=np.linalg.norm(mat,2)
   check(error<1e-12,'live-wire circuit equals index contraction')
   check(norm<=1+1e-12,'global propercut product bound')
   max_error=max(max_error,error);max_norm=max(max_norm,norm);ncuts+=1
 fixture_out.append(dict(name=name,tested_cuts=ncuts,max_circuit_contraction_error=max_error,max_global_cut=max_norm))

# Symmetric counterexample with a forbidden unmarked bubble.
bubbles=[]
for D in [2,3,4,8,16]:
 u=np.eye(D)[0];I=np.eye(D)
 A=(np.einsum('i,jk->ijk',u,I)+np.einsum('j,ik->ijk',u,I)+np.einsum('k,ij->ijk',u,I))/math.sqrt(D+8)
 gram=A.reshape(D,D*D)@A.reshape(D,D*D).T
 check(np.linalg.norm(gram,2)<=1+1e-12,'bubble local allcut')
 traced=np.einsum('ijj->i',A)
 expected=(D+2)/math.sqrt(D+8)
 check(np.linalg.norm(traced-expected*u)<1e-12,'bubble exact trace')
 check(expected>1,'no-marked-tree counterexample amplifies')
 bubbles.append(dict(D=D,global_norm=expected,local_cut=math.sqrt(np.linalg.norm(gram,2))))

result=dict(assertions=checks,scope='Marked-tree orientations/live-wire identity and explicit numerical tensor contractions; no native sampler execution.',exhaustive_tree_census=census,fixtures=fixture_out,symmetric_unmarked_bubble=bubbles)
(HERE/'marked_spanning_tree_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
