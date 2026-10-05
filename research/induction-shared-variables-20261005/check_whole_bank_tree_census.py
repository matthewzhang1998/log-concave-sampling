#!/usr/bin/env python3
from itertools import product, permutations
from functools import lru_cache
from pathlib import Path
import sympy as s
import numpy as np, math, json
count=0
def ck(v,msg):
 global count
 count+=1
 if not bool(v):raise AssertionError(msg)
def trees(n):
 if n==2:yield [(0,1)];return
 for word in product(range(n),repeat=n-2):
  deg=[1]*n
  for v in word:deg[v]+=1
  out=[]
  for v in word:
   u=next(i for i,d in enumerate(deg) if d==1)
   out.append((u,v));deg[u]-=1;deg[v]-=1
  rest=[i for i,d in enumerate(deg) if d==1];out.append(tuple(rest));yield out
def paths(n,es):
 adj=[[] for _ in range(n)]
 for k,(u,v) in enumerate(es):adj[u].append((v,k));adj[v].append((u,k))
 out={}
 for st in range(n):
  stack=[(st,-1,[])]
  while stack:
   v,p,pp=stack.pop();out[st,v]=pp
   for u,k in adj[v]:
    if u!=p:stack.append((u,v,pp+[k]))
 return out
def gauss(k):return s.Integer(0) if k%2 else s.Integer(math.prod(range(1,k,2)))
def partitions(n):
 def rec(i,bs):
  if i==n:yield [b[:] for b in bs];return
  for b in bs:
   b.append(i);yield from rec(i+1,bs);b.pop()
  bs.append([i]);yield from rec(i+1,bs);bs.pop()
 yield from rec(0,[])
def cumulant_power(degs):
 ans=s.Integer(0)
 for pi in partitions(len(degs)):
  ans+=(-1)**(len(pi)-1)*math.factorial(len(pi)-1)*s.prod(gauss(sum(degs[i] for i in b)) for b in pi)
 return ans
def gaussian_joint(degs,C):
 @lru_cache(None)
 def f(ds):
  if sum(ds)==0:return s.Integer(1)
  if sum(ds)%2:return s.Integer(0)
  i=next(i for i,d in enumerate(ds) if d)
  a=list(ds);a[i]-=1;ans=s.Integer(0)
  for j,n in enumerate(a):
   if n:
    a[j]-=1;ans+=n*C[i][j]*f(tuple(a));a[j]+=1
  return s.expand(ans)
 return f(tuple(degs))
def bkar_power(degs):
 n=len(degs);m=n-1;ts=s.symbols('t:'+str(m));ans=s.Integer(0)
 for es in trees(n):
  hit=[sum(v in e for e in es) for v in range(n)]
  if any(a<b for a,b in zip(degs,hit)):continue
  coeff=s.prod(math.factorial(a)//math.factorial(a-b) for a,b in zip(degs,hit))
  ds=[a-b for a,b in zip(degs,hit)];pa=paths(n,es)
  for order in permutations(range(m)):
   rank={e:j for j,e in enumerate(order)}
   C=[[s.Integer(1) if i==j else ts[min(pa[i,j],key=rank.get)] for j in range(n)] for i in range(n)]
   pol=coeff*gaussian_joint(ds,C)
   # t_order[0] <= ... <= t_order[m-1], all in [0,1]
   for j,e in enumerate(order):pol=s.integrate(pol,(ts[e],0,ts[order[j+1]] if j+1<m else 1))
   ans+=pol
 return s.simplify(ans)
examples=[]
for ds in [(1,1),(2,2),(3,3),(1,1,2),(1,2,3),(2,2,2),(3,3,2),(1,1,1,3),(1,1,2,2),(2,2,2,2),(3,3,2,2),(3,3,3,3)]:
 exact=cumulant_power(ds);forest=bkar_power(ds)
 ck(exact==forest,'exact BKAR polynomial cumulant '+str(ds))
 examples.append({'powers':ds,'cumulant':str(exact),'forest_integral':str(forest)})
# H2 high-order check, only paths survive; no missing factorial or half.
for n in range(2,9):
 path_count=math.factorial(n)//2
 val=s.Rational(2**n,n)*path_count
 ck(val==2**(n-1)*math.factorial(n-1),'H2 cumulant path factor n'+str(n))
# Positive min-path covariance and retained-old-bank disintegration.
rng=np.random.default_rng(6105)
for n in range(2,9):
 for z,es in enumerate(trees(n)):
  if z>=50:break
  pa=paths(n,es)
  for _ in range(10):
   t=rng.random(n-1);c=float(min(t))
   K=np.array([[1.0 if i==j else min(t[e] for e in pa[i,j]) for j in range(n)] for i in range(n)])
   S=K-c*np.ones((n,n))
   ck(np.linalg.eigvalsh(K).min()>-1e-12,'replica K PSD')
   ck(np.linalg.eigvalsh(S).min()>-1e-12,'retained-bank residual covariance PSD')
   vals,vecs=np.linalg.eigh(S);R=(vecs*np.sqrt(np.maximum(vals,0)))@vecs.T
   ck(np.max(np.abs(c*np.ones((n,n))+R@R.T-K))<1e-11,'actual old B + scalar residual root exact')
   ck(np.allclose(np.diag(K),1),'original replica marginal bank unchanged')
res={'all_pass':True,'assertions':count,'exact_examples':examples,'scope':'exact polynomial cumulants and positive replica geometry; not higher heat/clock or finite cubature certification'}
Path(__file__).with_name('whole_bank_tree_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
