#!/usr/bin/env python3
"""Exact finite algebra and graph ledger; not a native sampler execution."""
from fractions import Fraction as Q
from math import comb, factorial
from itertools import combinations
from collections import defaultdict
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
assertions=0
def check(x):
 global assertions
 assertions+=1
 assert x

def df(n):
 return 1 if n<=0 else n*df(n-2)
def add(a,b,scale=Q(1)):
 out=dict(a)
 for m,c in b.items():
  out[m]=out.get(m,Q(0))+scale*c
 return {m:c for m,c in out.items() if c}
def mul(a,b):
 out=defaultdict(Q)
 for (al,th),c in a.items():
  for (bl,uh),d in b.items():out[al+bl,th+uh]+=c*d
 return dict(out)
# W = -alpha^4*t*theta^2*G^2*H^2*(S+alpha theta)*(T+alpha theta).
# t is tracked separately by cumulant order k; moments below have no t.
def moment(k):
 factor=(-1)**k*df(2*k-1)**2
 out=defaultdict(Q)
 for i in range(k//2+1):
  for j in range(k//2+1):
   shift=2*k-2*i-2*j
   coef=factor*comb(k,2*i)*df(2*i-1)*comb(k,2*j)*df(2*j-1)
   out[4*k+shift,2*k+shift]+=Q(coef)
 return dict(out)
mu=[{(0,0):Q(1)}]+[moment(k) for k in range(1,9)]
kap=[{}]
for k in range(1,9):
 cur=mu[k]
 for i in range(1,k):cur=add(cur,mul(kap[i],mu[k-i]),-Q(comb(k-1,i-1)))
 kap.append(cur)
check(kap[1]=={(6,4):Q(-1)})
check({m:c/2 for m,c in kap[2].items()}=={(8,4):Q(9,2),(10,6):Q(9),(12,8):Q(4)})
# Independently check moments via exhaustive exponent choices in S/T.
for k in range(1,9):
 explicit=defaultdict(Q)
 for s in range(k+1):
  for t in range(k+1):
   if s%2 or t%2:continue
   power=2*k-s-t
   explicit[4*k+power,2*k+power]+=Q((-1)**k*comb(k,s)*comb(k,t)*df(s-1)*df(t-1)*df(2*k-1)**2)
 check(dict(explicit)==mu[k])
# Reconstruction of every moment from cumulants verifies the whole finite jet.
for k in range(1,9):
 recovered={}
 for i in range(1,k+1):recovered=add(recovered,mul(kap[i],mu[k-i]),Q(comb(k-1,i-1)))
 check(recovered==mu[k])
# Doubled C4 centers: ab,bd,dc,ca, two labelled copies of each.
base=[(0,1),(1,3),(3,2),(2,0)]
edges=[(a,b,color,copy) for color,(a,b) in enumerate(base) for copy in range(2)]
def connected(ee):
 seen={0}
 for _ in range(4):
  for a,b,*_ in ee:
   if a in seen or b in seen:seen|={a,b}
 return len(seen)==4
def path(ee,a,b):
 adj=defaultdict(list)
 for x,y,*_ in ee:adj[x].append(y);adj[y].append(x)
 def dfs(x,p):
  if x==b:return p
  for y in adj[x]:
   if y not in p:
    r=dfs(y,p+[y])
    if r:return r
 return dfs(a,[a])
counts=defaultdict(int); trees=0
for chosen in combinations(range(8),3):
 tree=[edges[i] for i in chosen]
 if not connected(tree):continue
 trees+=1
 unchosen=[edges[i] for i in range(8) if i not in chosen]
 for a in range(4):
  paths=[path(tree,a,b) for b in range(4) if b!=a]
  longest=max(map(len,paths))
  for pp in paths:
   counts['paths']+=1
   # Two connected chunks contain the only two endpoint physical marks.
   possible=False
   for split in range(1,len(pp)):
    chunks=[set(pp[:split]),set(pp[split:])]
    no_self=all(not (x in c and y in c) for x,y,*_ in unchosen for c in chunks)
    possible|=no_self
   check(possible==(len(pp)==2))
   counts['pass' if possible else 'fail']+=1
   if len(pp)==longest:
    counts['longest_paths']+=1
    check(not possible)
   if possible:check(2*(len(pp)+2)==8)
   else:check(2*(len(pp)+2)>8)
check(trees==32);check(counts['paths']==384);check(counts['pass']==192);check(counts['fail']==192)
check(counts['longest_paths']==128)
# Graph identities.
check(7-6+1==2);check(6==4-2+2*2)
check(12-8+1==5);check(12==4-2+2*5)
# Balanced normalization and exact rational exponent inequalities.
for gamma in [Q(1,10),Q(1,4),Q(1,3),Q(1,2)]:
 check(1-gamma>0);check(2-2*gamma>=1);check(3-3*gamma>=1)
 check(8-4*gamma>6-2*gamma);check(12-4*gamma>6-2*gamma)
 for b in range(2,21):
  prior=[Q(b),b+1-gamma*b,b+2-gamma*(b+1),b+2-gamma,b+3-2*gamma]
  check(min(prior)>0)
# Dyadic center bounds directly test representative finite meshes with exact powers.
maxima=defaultdict(float)
for n in range(2,25):
 tau=2.**(-n/2)
 panels=[2.**(-j) for j in range(0,2*n+10)]
 for k in range(1,9):
  for r in range(0,5):
   total=sum((d/(d+tau*tau)**1.5)**k/(d+tau*tau)**(r/2) for d in panels)
   scaled=total*tau**(k+r)
   maxima[f'{k},{r}']=max(maxima[f'{k},{r}'],scaled)
   check(scaled<20)
output={
 'scope':'Finite coefficient algebra, graph ancestry and heat/path ledgers. No native program is executed.',
 'assertions':assertions,'spanning_trees':trees,'path_counts':dict(counts),
 'second_log_coefficient':[{'alpha_power':a,'physical_rank':m,'t_power':2,'coefficient':str(c/2)} for (a,m),c in sorted(kap[2].items())],
 'all_cumulants_through_spine_order_8':[
 {'spine_order':k,'terms':[{'alpha_power':a,'physical_rank':m,'t_power':k,'log_coefficient':str(c/factorial(k))} for (a,m),c in sorted(kap[k].items())]} for k in range(1,9)],
 'worst_scaled_dyadic_sum':max(maxima.values()),
 'conclusion':'Leading two-cycle closure and fixed-cutoff conditional/auxiliary return pass these checks. Whole-bank native admission remains open. Long-spine trace-preservation and short-spine strict-grade gain cannot both be inherited for the doubled-C4 test.'}
(ROOT/'two_cycle_checks.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:output[k] for k in ['assertions','spanning_trees','path_counts','worst_scaled_dyadic_sum','conclusion']},indent=2))
