import itertools,json,collections
base=[(0,1),(1,2),(3,4),(4,5)]; centers={1,4}; leaves={0,2,3,5}
out=[]
for a in range(3):
 for b in range(3,6):
  edges=base+[(a,b)];adj={i:set() for i in range(6)}
  for x,y in edges:adj[x].add(y);adj[y].add(x)
  order={i:(1 if i in centers else 0)+(1 if i in {a,b} else 0) for i in range(6)}
  assert all(order[i]==len(adj[i])-1 for i in range(6))
  root_candidates=leaves-{a,b}
  def path(x,y):
   stack=[(x,[x])]
   while stack:
    v,p=stack.pop()
    if v==y:return p
    stack +=[(w,p+[w]) for w in adj[v] if w not in p]
  candidates=[path(x,y) for x,y in itertools.combinations(root_candidates,2)]
  spine=max(candidates,key=len);side=set(range(6))-set(spine)
  assert all(order[v]==0 for v in side)
  assert len(side)<=2 and 4<=len(spine)<=6
  out.append({'hit_pair':[a,b],'edges':edges,'adapter_orders':order,'spine':spine,'side_leaves':sorted(side),'first_conditional_force_count':2*len(spine)})
hist=collections.Counter(len(x['spine']) for x in out)
assert hist=={4:1,5:4,6:4}
p='/workspace/shared/old-bank-sixth-return-20261005/nine_hit_tree_checks.json'
open(p,'w').write(json.dumps({'spine_census':dict(hist),'trees':out,'scope':'All nine derivative-hit trees retain six marks; longest root spines include every non-C0 vertex, and all side leaves are C0. Native programs not executed.'},indent=2)+'\n')
print(dict(hist))
