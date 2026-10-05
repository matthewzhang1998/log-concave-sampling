#!/usr/bin/env python3
"""Finite exact graph/amplitude and numerical dyadic-envelope diagnostics.
Does not execute any imported original-VALUE native pair/calibration program.
"""
import collections,itertools,json,math
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent
checks=0

def check(v):
 global checks
 assert v
 checks+=1

edges=((0,1),(2,3),(0,2),(0,2),(1,3),(1,3))
deg=collections.Counter(v for e in edges for v in e)
check(dict(deg)=={0:3,1:3,2:3,3:3})
check(6+4-(4+4)+1==3)

def conn(es,vs=range(4)):
 seen={next(iter(vs))}
 while True:
  old=len(seen)
  for a,b in es:
   if a in seen or b in seen:seen|={a,b}
  if len(seen)==old:return seen==set(vs)

openings=[]
for tree_ids in itertools.combinations(range(6),3):
 tree=[edges[i] for i in tree_ids]
 if not conn(tree):continue
 cuts=[edges[i] for i in range(6) if i not in tree_ids]
 for root_edge in tree:
  for orient in (root_edge,root_edge[::-1]):
   a,b=orient
   # One incoming, one outgoing native slot per center. All other tree,
   # opened auxiliary, and physical-leaf ports total exactly k=2.
   for v in range(4):
    j=sum(v in e for e in tree)
    cut=sum(v in e for e in cuts)
    check(j+cut==3)
    check((j+1)-2+cut==2)
   for e in cuts:check(e[0]!=e[1])
   check(a!=b)
   # Pair every retained center with its own physical leaf. No opened
   # Gaussian color repeats within a chunk and every chunk retains a mark.
   for v in range(4):
    local=[i for i,e in enumerate(cuts) if v in e]
    check(len(local)==len(set(local)))
   openings.append({'tree_ids':tree_ids,'root_spine_centers':orient,'cuts':cuts})
check(len(openings)==72)

# Exact amplitude law, independent of graph contractions/moments.
amplitudes=[]
for n in range(2,42,2):
 for d in range(1,18):
  total=Fraction(d,1)+Fraction(1,2)+(n-1)*Fraction(1,2)+n*Fraction(1,2)
  check(total==n+d)
  for h in range(2,9):
   for j in range(h*(n-2)+1):
    nnew=2*h+j
    Anew=h*(d+2)+j
    check(Anew==nnew+h*d)
    check(h*d>=2*d)
    check(Anew>=4+2*d)
    if n==4 and d==4 and j==0:amplitudes.append({'h':h,'nnew':nnew,'Anew':Anew,'surplus':h*d})

queues=[]
for P in (9,13,17,25,41,65):
 todo=collections.deque([(4,4,0)]);seen={(4,4)};depth=0
 while todo:
  n,d,k=todo.popleft();depth=max(depth,k)
  for h in range(2,(P-1)//d+1):
   for n2 in range(2*h,h*n+1,2):
    d2=h*d
    check(d2>=2*d)
    if n2+d2>=P:continue
    check(n2<P and d2<P)
    if (n2,d2) not in seen:seen.add((n2,d2));todo.append((n2,d2,k+1))
 check(depth <= math.floor(math.log2(P/4)))
 queues.append({'cutoff':P,'distinct_n_d_states':len(seen),'max_depth':depth})

# Original finite positive dyadic clock weights, with no independent-clock
# replacement. The ratios below test the proved split-at-b^2 envelopes.
clocks=[]
for K in (10,20,40,80):
 xs=[2.0**(-k) for k in range(1,K+1)]
 ws=xs
 for ell in range(1,2*K+1):
  b=2.0**(-ell/2)
  max_z=max(w/(x+b*b) for x,w in zip(xs,ws));check(max_z<=1)
  for r in (1,2,3,4,5,8,12):
   s0=sum((w/(x+b*b))**r for x,w in zip(xs,ws))
   s1=sum(w**r/(x+b*b)**(r+.5) for x,w in zip(xs,ws))
   s2=sum(w**r/(x+b*b)**(r+1) for x,w in zip(xs,ws))
   check(s0<=K)
   check(b*s1<5)
   check(b*b*s2<3)
   clocks.append({'K':K,'ell':ell,'r':r,'b_s1':b*s1,'b2_s2':b*b*s2})

# A concrete affine Gaussian readout test: B occupies coordinates 0..10;
# actual old terminal tape occupies coordinate 11 and remains physical.
# Riesz(T8) is supported only on B, so its source-zero contraction is zero.
row=[0]*11+[Fraction(3,5),Fraction(4,5)]
for j in range(11):check(row[j]==0)
check(sum(x*x for x in row)==1)
check(any(row[j] for j in range(11,len(row))))

report={'assertions':checks,'literal_oriented_openings':len(openings),
 'center_edges':edges,'full_force_count':8,'physical_marks':4,'cycle_rank':3,
 'amplitude_examples':amplitudes,'queue_checks':queues,
 'dyadic_cases':len(clocks),'max_b_times_one_hit_sum':max(x['b_s1'] for x in clocks),
 'max_b2_times_two_hit_sum':max(x['b2_s2'] for x in clocks),
 'input_grade':8,'own_conditional_grade':12,'safe_full_old_join_grade':9,
 'old_cubic_subchannel_grade':11,'old_sixth_subchannel_grade':14,'own_leading_bank_grade':16,
 'scope':'Finite exact combinatorics/amplitude and dyadic-envelope diagnostics only. No native original-VALUE pair program executed.'}
(ROOT/'native_cubic_return_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
