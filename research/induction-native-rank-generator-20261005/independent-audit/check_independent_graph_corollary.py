"""Independent all-cut multigraph corollary checks.
A literal live-wire matrix circuit is compared with full einsum contraction.
The original audited source and audit files are left unchanged.
"""
from pathlib import Path
import json,math
import numpy as np
HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(13122026)
checks=0;graphs=0;circuits=0;steps=0;negative=0;max_error=0.

def ck(ok,msg):
 global checks
 checks+=1
 if not bool(ok):raise AssertionError(msg)

def maxcut(arr):
 best=0.
 for mask in range(1,(1<<arr.ndim)-1):
  left=[i for i in range(arr.ndim) if mask>>i&1];right=[i for i in range(arr.ndim) if not(mask>>i&1)]
  M=arr.transpose(left+right).reshape(math.prod(arr.shape[i] for i in left),-1)
  best=max(best,float(np.linalg.norm(M,2)))
 return best

def orientation(n,tree,extra,marks,inputs):
 m=len(marks);weights=[2**(2**j) for j in range(m)]
 si=sum(weights[j] for j in inputs);so=sum(weights[j] for j in range(m) if j not in inputs)
 mass=[w*so if j in inputs else -w*si for j,w in enumerate(weights)]
 directed=[]
 for index,(a,b) in enumerate(tree):
  component={a}
  while True:
   more=set(component)
   for ei,(u,v) in enumerate(tree):
    if ei==index:continue
    if u in component:more.add(v)
    if v in component:more.add(u)
   if more==component:break
   component=more
  f=sum(mass[j] for j,v in enumerate(marks) if v in component)
  ck(f!=0,'tree flow nonzero')
  directed.append((a,b) if f>0 else (b,a))
 # Topological sort of the directed tree, before adding nontree edges.
 order=[];remaining=set(range(n))
 while remaining:
  ready=[v for v in sorted(remaining) if not any(b==v and a in remaining for a,b in directed)]
  ck(bool(ready),'tree orientation acyclic')
  for v in ready:order.append(v);remaining.remove(v)
 pos={v:i for i,v in enumerate(order)}
 for a,b in extra:directed.append((a,b) if pos[a]<pos[b] else (b,a))
 ck(all(pos[a]<pos[b] for a,b in directed),'all extra edges preserve DAG')
 return directed,order

def live_wire(arrays,local,edges,marks,inputs,dims,boundary):
 global steps
 n=len(arrays);tree=edges[:n-1];extra=edges[n-1:]
 directions,order=orientation(n,tree,extra,marks,inputs)
 incoming=[[] for _ in range(n)];outgoing=[[] for _ in range(n)]
 for ei,(a,b) in enumerate(directions):outgoing[a].append(ei);incoming[b].append(ei)
 for j,v in enumerate(marks):(incoming if j in inputs else outgoing)[v].append(boundary[j])
 input_labels=[boundary[j] for j in sorted(inputs)];output_labels=[j for j in boundary if j not in input_labels]
 N=math.prod(dims[j] for j in input_labels);live=input_labels.copy()
 current=np.eye(N).reshape([dims[j] for j in live]+[N])
 for v in order:
  ins=incoming[v];outs=outgoing[v]
  ck(bool(ins) and bool(outs),'proper local cut after adding extra edges')
  ck(set(ins)<=set(live),'all incoming wires available')
  ck(not(set(outs)&set(live)),'new wires are genuinely fresh')
  carried=[j for j in live if j not in ins]
  row_order=[live.index(j) for j in ins+carried]+[len(live)]
  M=current.transpose(row_order).reshape(math.prod(dims[j] for j in ins),-1)
  local_order=[local[v].index(j) for j in outs+ins]
  A=arrays[v].transpose(local_order).reshape(math.prod(dims[j] for j in outs),-1)
  ck(np.linalg.norm(A,2)<=1+1e-11,'actual local operator contraction')
  current=(A@M).reshape([dims[j] for j in outs+carried]+[N]);live=outs+carried
  ck(np.linalg.norm(current.reshape(-1,N),2)<=1+2e-10,'live-wire global operator remains bounded')
  steps+=1
 ck(set(live)==set(output_labels),'only final physical outputs remain')
 return current.transpose([live.index(j) for j in output_labels]+[len(live)]).reshape(-1,N)

for n in range(2,7):
 for rep in range(4):
  tree=[(int(rng.integers(v)),v) for v in range(1,n)]
  degree=[sum(v in e for e in tree) for v in range(n)]
  marks=[v for v in range(n) if degree[v]==1]
  if rep%2:marks.append(int(rng.integers(n)))
  extra=[tuple(map(int,rng.choice(n,2,replace=False))) for _ in range(1+rep%3)]
  edges=tree+extra;local=[[] for _ in range(n)];dims={}
  for i,(a,b) in enumerate(edges):local[a].append(i);local[b].append(i);dims[i]=2+i%2
  boundary=[]
  for j,v in enumerate(marks):
   lab=len(edges)+j;boundary.append(lab);local[v].append(lab);dims[lab]=2
  arrays=[]
  for labs in local:
   A=rng.normal(size=tuple(dims[j] for j in labs));A/=maxcut(A);arrays.append(A)
  arguments=[]
  for A,labs in zip(arrays,local):arguments.extend((A,labs))
  target=np.einsum(*arguments,boundary,optimize=True)
  for mask in range(1,(1<<len(marks))-1):
   inputs={j for j in range(len(marks)) if mask>>j&1};outputs=[j for j in range(len(marks)) if j not in inputs]
   direct=target.transpose(outputs+sorted(inputs)).reshape(2**len(outputs),-1)
   circuit=live_wire(arrays,local,edges,marks,inputs,dims,boundary)
   error=float(np.max(np.abs(circuit-direct)));max_error=max(max_error,error)
   ck(error<3e-12,'literal circuit equals full contraction')
   ck(np.linalg.norm(direct,2)<=1+2e-10,'full cyclic coefficient all-cut product bound')
   circuits+=1
  graphs+=1

# Sharp caution: a loopless triangle whose only marked vertex is 0.
# No spanning tree can have all leaves marked. All local proper-cut norms <=1,
# but contracting the triangle produces sqrt(D) times a unit rank-one matrix.
for D in (2,3,7):
 u=np.eye(D)[0];v=np.eye(D)[-1]
 A=np.einsum('a,b,ij->abij',u,v,np.eye(D))/math.sqrt(D)
 B=np.eye(D);C=np.eye(D)
 ck(abs(maxcut(A)-1)<1e-12 and maxcut(B)==1 and maxcut(C)==1,'negative fixture local cut norms')
 result=np.einsum('abij,jk,ki->ab',A,B,C)
 ck(np.max(np.abs(result-math.sqrt(D)*np.outer(u,v)))<1e-12,'unmarked-cycle trace evaluates exactly')
 ck(abs(np.linalg.norm(result,2)-math.sqrt(D))<1e-12 and np.linalg.norm(result,2)>1,'dropping marked-leaf hypothesis fails')
 negative+=1

out={'status':'PASS','assertions':checks,'loopless_multigraphs':graphs,'physical_cut_circuits':circuits,'live_wire_steps':steps,'maximum_circuit_einsum_error':max_error,'missing_mark_hypothesis_counterexamples':negative,'scope':'All-cut coefficient tensor theorem for a loopless multigraph admitting a marked-leaf spanning tree. No native nonlinear program, derivative-frame, self-trace, or return theorem.'}
(HERE/'graph-corollary-independent-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
