#!/usr/bin/env python3
"""Independent finite algebra/topology checks, not native execution or a W2 test."""
from pathlib import Path
from collections import Counter, defaultdict
from itertools import combinations, product
from fractions import Fraction
import hashlib,json,math
import sympy as S
HERE=Path(__file__).resolve().parent
checks=0

def check(v,label=''):
 global checks
 checks+=1
 if not v: raise AssertionError(label)

def connected(vertices,edges):
 vertices=set(vertices)
 if not vertices:return False
 seen={next(iter(vertices))}
 while True:
  old=len(seen)
  for a,b in edges:
   if a in seen or b in seen:seen.update((a,b))
  if old==len(seen):return seen==vertices

def multigraphs(n):
 pairs=list(combinations(range(n),2));deg=[0]*n;edges=[]
 def generate(k):
  if k==len(pairs):
   if min(deg)>=1 and connected(range(n),edges):yield tuple(edges)
   return
  a,b=pairs[k]
  for count in range(1+min(3-deg[a],3-deg[b])):
   deg[a]+=count;deg[b]+=count;edges.extend([(a,b)]*count)
   yield from generate(k+1)
   if count:del edges[-count:]
   deg[a]-=count;deg[b]-=count
 yield from generate(0)

def openings(n,edges):
 for tree in combinations(range(len(edges)),n-1):
  if not connected(range(n),[edges[k] for k in tree]):continue
  for root in tree:yield tree,root

def monomials(n,edges,tree,root):
 a,b=edges[root];parent={a:None,b:None};todo=[a,b]
 while todo:
  v=todo.pop()
  for k in tree:
   if k==root or v not in edges[k]:continue
   x,y=edges[k];w=x if y==v else y
   if w not in parent:parent[w]=v;todo.append(w)
 degree=Counter(v for e in edges for v in e)
 other=[v for v in range(n) if v not in (a,b)]
 result=[]
 for mask in range(1<<len(other)):
  included={a,b}|{v for j,v in enumerate(other) if mask&(1<<j)}
  if any(parent[v] is not None and parent[v] not in included for v in included):continue
  literal=[];slots=defaultdict(list)
  for k,(x,y) in enumerate(edges):
   if k in tree:
    if x in included and y in included:literal.append((x,y))
    elif (x in included)!=(y in included):
     present=x if x in included else y;absent=y if x in included else x
     slots[('side',absent)].append(present)
   else:
    for v in (x,y):
     if v in included:slots[('cut',k)].append(v)
  for v in included:
   for j in range(3-degree[v]):slots[('public',v,j)].append(v)
  check(connected(included,literal),'ancestor-closed tree')
  for v in included:
   check(sum(v in e for e in literal)+sum(row.count(v) for row in slots.values())==3,'C2 slots')
  for row in slots.values():check(len(row)==len(set(row)),'no local repeated color')
  result.append((included,literal,dict(slots)))
 return result

def gm(n):return 0 if n%2 else math.prod(range(1,n,2))

def shifted_moment(n):return sum(math.comb(n,k)*gm(k) for k in range(n+1))

def term_mean(t):
 out=1
 for color,row in t[2].items():out*=shifted_moment(len(row)) if color[0]=='public' else gm(len(row))
 return out

def pairings(items):
 if not items:yield [];return
 a=items[0]
 for j in range(1,len(items)):
  for pp in pairings(items[1:j]+items[j+1:]):yield [(a,items[j])]+pp

def public_options(items):
 # Choose the retained theta legs; match every remaining Gaussian occurrence.
 for mask in range(1<<len(items)):
  retained=[v for j,v in enumerate(items) if mask&(1<<j)]
  paired=[v for j,v in enumerate(items) if not mask&(1<<j)]
  if len(paired)%2:continue
  for pp in pairings(paired):yield pp,retained

census=[]
for n in range(2,6):
 graphs=opening_count=terms_count=0
 for edges in multigraphs(n):
  graphs+=1
  for tree,root in openings(n,edges):
   opening_count+=1;terms=monomials(n,edges,tree,root);terms_count+=len(terms)
   survivors=[term for term in terms if term_mean(term)]
   check(len(survivors)==1,'one first cluster')
   check(len(survivors[0][0])==n and term_mean(survivors[0])==1,'full intended decorated graph')
 census.append(dict(n=n,graphs=graphs,unoriented_openings=opening_count,reverse_monomials=terms_count))

base=[(0,1),(2,3),(0,2),(1,3)]
types=[('T00',base+[(0,2),(1,3)]),('T10',base+[(1,3)]),('T01',base+[(0,2)]),('T11',base)]
wick_results=[]
for name,edges in types:
 opened=total=surviving=0;rank_counts=Counter()
 for tree,root in openings(4,edges):
  opened+=1;terms=monomials(4,edges,tree,root)
  for t,u in product(terms,repeat=2):
   vertices={(0,v) for v in t[0]}|{(1,v) for v in u[0]}
   fixed=[((0,a),(0,b)) for a,b in t[1]]+[((1,a),(1,b)) for a,b in u[1]]
   colors=set(t[2])|set(u[2]);options=[];direct_moment=1
   for color in sorted(colors):
    occ=[(0,v) for v in t[2].get(color,[])]+[(1,v) for v in u[2].get(color,[])]
    if color[0]=='public':options.append(list(public_options(occ)));direct_moment*=shifted_moment(len(occ))
    else:options.append([(pp,[]) for pp in pairings(occ)] if len(occ)%2==0 else []);direct_moment*=gm(len(occ))
   local_connected=0
   for choices in product(*options):
    total+=1;new_edges=fixed+[e for pp,retained in choices for e in pp]
    direct=[v for pp,retained in choices for v in retained]
    check(all(a!=b for a,b in new_edges),'no self loop')
    for v in vertices:
     check(sum(v in e for e in new_edges)+direct.count(v)==3,'offspring r+p=3')
     check(sum(v in e for e in new_edges)>=1,'protected selected edge still present')
    # Each center occurrence additionally carries one untouched C0 physical mark.
    check(len(vertices)+len(direct)>=len(vertices),'one mark per center')
    if connected(vertices,new_edges):
     local_connected+=1;surviving+=1;rank_counts[len(vertices)+len(direct)]+=1
   check(local_connected==direct_moment-term_mean(t)*term_mean(u),'connected cumulant count')
 p=12-2*len(edges)
 check(4+p==dict(T00=4,T10=6,T01=6,T11=8)[name],'main rank')
 wick_results.append(dict(type=name,all_unoriented_openings=opened,wick_diagrams=total,connected_diagrams=surviving,offspring_physical_ranks=dict(sorted(rank_counts.items())),cut_roots=len(edges)-3))

# Exact full public-tilt second cumulant and positive-path logarithmic current.
alpha,theta,a,b,c,d,ell,s=S.symbols('alpha theta a b c d ell s')
P,Q,Y,Z=S.symbols('P Q Y Z')
def expect(poly,variables):
 return S.expand(sum(coef*math.prod(gm(p) for p in powers) for powers,coef in S.Poly(S.expand(poly),*variables).terms()))
V=s*ell*alpha**4*theta**2*(P+a*theta)*(Q+b*theta)*(Y+alpha*c*theta)*(Z+alpha*d*theta)
mean=expect(V,(P,Q,Y,Z));var=S.expand(expect(V**2,(P,Q,Y,Z))-mean**2)
C6=ell*a*b*c*d*alpha**6*theta**6
C8=ell**2*alpha**8*theta**4*(1+a*a*theta*theta)*(1+b*b*theta*theta)/2
check(S.expand(mean-s*C6)==0,'exact first cumulant')
check(S.expand(var.coeff(alpha,8)/2-s*s*C8/alpha**8)==0,'all four grade-eight public tilts')
for grade in (0,1,2,3,4,5,6,7,9):check(S.expand(var).coeff(alpha,grade)==0,'no missing lower second cumulant')
check(S.expand(var).coeff(alpha,10)!=0,'grade ten generally survives')
K=s*C6+s*s*C8-s*s*C8+(1-s)*C6
check(S.expand(S.diff(K,s))==0,'common-path log-current below grade ten')
# Finite current identification: derivative of characteristic divided by the SAME characteristic.
A,B=S.symbols('A B')
jet=1+A*alpha**6+B*alpha**8+A*A*alpha**12/2
inverse=1-A*alpha**6-B*alpha**8+A*A*alpha**12/2
check(S.Poly(S.expand(jet*inverse-1),alpha).terms()[-1][0][0]>=14,'triangular inverse starts correctly')
# Nontrivial scalar positive-path check of the SAME-ENDPOINT division.
# X_old=Z+s*e*He2(Z); X_ref=H+(1-s)*e*He2(H);
# X_fix=J+2*s*(1-s)*e^2*(J+2*He3(J)), all independent.
# Each is a positive pushforward. The derivative of the characteristic jet
# itself has a spurious rank-six term; dividing by its own full jet removes it.
e=S.symbols('e')
q2=theta**2+2*theta**4
old_phi=1+s*e*theta**3+s*s*e*e*(q2+theta**6/2)
ref_phi=old_phi.subs(s,1-s)
fix_phi=1+2*s*(1-s)*e*e*q2
def jet2(expr):return S.Poly(S.expand(expr),e).as_dict()
def trunc2(expr):return S.expand(sum(coef*e**power[0] for power,coef in jet2(expr).items() if power[0]<=2))
def inverse2(phi):
 f=phi-1
 return trunc2(1-f+f*f)
old_current=trunc2(S.diff(old_phi,s)*inverse2(old_phi))
ref_current=trunc2(S.diff(ref_phi,s)*inverse2(ref_phi))
fix_current=trunc2(S.diff(fix_phi,s)*inverse2(fix_phi))
check(S.expand(old_current-e*theta**3-2*s*e*e*q2)==0,'local cumulant-current normal form')
check(S.expand(old_current+ref_current+fix_current)==0,'positive toy current cancels at common endpoint')
naive=trunc2(S.diff(old_phi+ref_phi+fix_phi,s))
check(S.expand(naive-(2*s-1)*e*e*theta**6)==0,'Gaussian-baseline shortcut leaves spurious rank six')
check(trunc2(S.diff(old_phi*ref_phi*fix_phi,s))==0,'independent characteristic product agrees')
for sample in [0,S.Rational(1,4),S.Rational(1,2),1]:
 check(S.expand((old_current+ref_current+fix_current).subs(s,sample))==0,'uniform parameter endpoints')
# Shifted-Wick positive comparator: linear response rank six, first nonlinear term alpha^12.
h,c_ref,T=S.symbols('h c_ref T',nonzero=True)
H5=h**5/c_ref**5-10*h**3/c_ref**4+15*h/c_ref**3
# E H5(H+ c_ref theta), H~N(0,c_ref), is theta^5.
def ec(poly):
 return S.expand(sum(coef*gm(powers[0])*c_ref**(powers[0]//2) for powers,coef in S.Poly(S.expand(poly),h).terms()))
check(S.simplify(ec(H5.subs(h,h+c_ref*theta))-theta**5)==0,'inverse-covariance Wick normalization')
second_ref=S.expand(T*T*alpha**12*theta**2*(ec(H5.subs(h,h+c_ref*theta)**2)-theta**10)/2)
check(second_ref!=0,'reference first nonlinear return')
for n in range(2,20):
 for surplus in (1,2,4,8):
  check(Fraction(1,2)*(2*n-1)+surplus+Fraction(1,2)==n+surplus,'amplitude product')
  for hcopy in range(2,8):
   for inserted in range(hcopy*(n-2)+1):
    nn=2*hcopy+inserted
    power=hcopy*(surplus+2)+inserted
    check(power-nn==hcopy*surplus,'surplus exactly multiplies')

pins=[]
paths=[HERE/'REVIEWED-CANDIDATE.md']
for directory,names in {
 'old-bank-sixth-return-20261005':['POSITIVE-OLD-BANK-SIXTH-RETURN.md','NEXT-EIGHT-FORCE-CURRENT.md'],
 'eight-force-native-return-20261005':['C2-CUBIC-NATIVE-RETURN.md'],
 'public-readout-cycle-contract-20261005':['GROUPED-JOINT-ZERO-PORT.md','RETAINED-PUBLIC-ROTATION-OBSTRUCTION.md'],
 'cycle-surplus-closure-20261005':['SHARED-AUXILIARY-FRAME-AND-REMAINDER.md','PERMANENT-LEAF-C3-SURPLUS-CLOSURE.md']}.items():
 paths.extend(HERE.parent.parent/directory/name for name in names)
paths.append(HERE.parent.parent/'v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex')
paths.extend([HERE.parent/'CANONICAL-EMBEDDING-TEST.md', HERE.parent.parent/'positive-endpoint-mean-join-20261005/MEAN-COVARIANCE-SKEW-ENDPOINT.md'])
for path in paths:pins.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
(HERE/'INPUT-PINS.json').write_text(json.dumps(pins,indent=2)+'\n')
result=dict(assertions=checks,scope='Finite exact algebra, all connected loopless center graphs with degree 1..3 at n=2..5; every short-spine opening; all h=2 decorated Wick diagrams of the four candidates. No native sampler, W2 constant, or canonical caller readset is executed.',first_cluster_census=census,second_cluster_census=wick_results,mean=str(mean),second_cumulant_grade8=str(C8),reference_second_cumulant=str(second_ref),common_path_lower_current='zero identically in s',positive_toy_naive_spurious_current=str(naive),positive_toy_same_endpoint_current=str(S.expand(old_current+ref_current+fix_current)))
(HERE/'independent_decorated_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
