#!/usr/bin/env python3
"""Independent finite diagnostics; not an implementation of the imported native sampler.
Exact formal coefficients, graph/port preservation, tensor cuts, and clock boundaries.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, comb
import hashlib, itertools, json
import numpy as np

HERE=Path(__file__).resolve().parent
checks=0
categories={}
def ck(condition, category, message):
    global checks
    if not bool(condition): raise AssertionError(f'{category}: {message}')
    checks+=1; categories[category]=categories.get(category,0)+1

def normal(n):
    if n%2: return 0
    return factorial(n)//(2**(n//2)*factorial(n//2))

# Three formal grades are (A, theta, G). No exponential moments are assumed.
CUTOFF=40
ONE={(0,0,0):F(1)}
def add(a,b):
    c=a.copy()
    for m,v in b.items():
        c[m]=c.get(m,F(0))+v
        if not c[m]: del c[m]
    return c

def scale(a,c): return {m:v*c for m,v in a.items() if v*c}
def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            z=tuple(i+j for i,j in zip(x,y))
            if z[0]<=CUTOFF: out[z]=out.get(z,F(0))+c*d
    return {m:v for m,v in out.items() if v}
def exp0(a):
    out=ONE.copy(); term=ONE.copy()
    for j in range(1,CUTOFF+1):
        term=scale(mul(term,a),F(1,j))
        if not term: break
        out=add(out,term)
    return out

def log1(a):
    z=add(a,{(0,0,0):F(-1)}); out={}; term=ONE.copy()
    for j in range(1,CUTOFF+1):
        term=mul(term,z)
        if not term: break
        out=add(out,scale(term,F((-1)**(j+1),j)))
    return out

def average_g(a):
    out={}
    for (k,t,g),v in a.items():
        m=(k,t,0); out[m]=out.get(m,F(0))+v*normal(g)
    return {m:v for m,v in out.items() if v}

def direct_joint(root,p,q):
    out=ONE.copy()
    for n in range(1,CUTOFF//6+1):
        for i in range(0,n+1,2):
            for j in range(0,n+1,2):
                a=8*n-i-j; t=4*n-i-j
                if a>CUTOFF: continue
                coefficient=F(root**n*normal(2*n)*comb(n,i)*normal(i)*p**(n-i)*comb(n,j)*normal(j)*q**(n-j),factorial(n))
                out=add(out,{(a,t,0):coefficient})
    return log1(out)

def conditional(root,p,q):
    out={}
    for h in range(1,CUTOFF//12+1):
        out=add(out,{(12*h,4*h,4*h):F(root**(2*h),2*h)})
        if 12*h+2<=CUTOFF:
            out=add(out,{(12*h+2,4*h+2,4*h):F(root**(2*h)*(p*p+q*q),2)})
    for h in range(CUTOFF//12+1):
        if 12*h+8<=CUTOFF:
            out=add(out,{(12*h+8,4*h+4,4*h+2):F(root**(2*h+1)*p*q)})
    return out

series={}
for root,p,q in itertools.product((-2,-1,1,2),(-3,1,2),(-2,1,3)):
    direct=direct_joint(root,p,q)
    nested=log1(average_g(exp0(conditional(root,p,q))))
    all_keys=set(direct)|set(nested)
    for key in all_keys: ck(direct.get(key,0)==nested.get(key,0),'exact_symbol',f'joint/nested mismatch {root,p,q,key}')
    ck(direct.get((8,4,0))==root*p*q,'normalization','leading eighth grade')
    ck(direct.get((12,4,0))==F(3,2)*root*root,'offspring','grade twelve trace')
    ck(direct.get((14,6,0))==F(3,2)*root*root*(p*p+q*q),'offspring','grade fourteen one-side')
    ck(direct.get((16,8,0))==(root*p*q)**2,'offspring','grade sixteen auxiliary covariance')
    ck(direct.get((20,8,0))==21*root**3*p*q,'offspring','grade twenty both-side plus mixed auxiliary')
    if (root,p,q) in ((1,1,1),(-1,1,1)):
        series[str(root)]={str(k):str(v) for k,v in sorted(direct.items())}

# Literal graph and scalar readout cancellation.
V=('L','C','Cp','Lp','M','D','Dp','Mp')
E=(('L','C'),('C','Cp'),('Cp','Lp'),('M','D'),('D','Dp'),('Dp','Mp'),('Cp','Dp'))
physical={'L','Lp','M','Mp'}; aux={'C','D'}
adj={v:set() for v in V}
for a,b in E: adj[a].add(b); adj[b].add(a)
ck(len(E)==len(V)-1,'geometry','opened graph edge count')
for v in V:
    j=len(adj[v])+int(v in physical)+int(v in aux)-1
    ck(j==(1 if v in physical else 2),'source_ports',f'original derivative order {v}')
    k=j-1
    ck(k in (0,1),'source_ports','only literal C0 and C1 needed')
spine=('L','C','Cp','Dp','D','M')
ck(all(spine[i+1] in adj[spine[i]] for i in range(len(spine)-1)),'geometry','six-force root spine')
for a in physical:
    dist={a:0}; todo=[a]
    for x in todo:
        for y in adj[x]:
            if y not in dist: dist[y]=dist[x]+1;todo.append(y)
    ck(len(dist)==8,'geometry','tree connected')
    ck(max(dist[b] for b in physical)<=5,'geometry','no longer physical-leaf spine')
chunks=({'L','C','Cp'},{'M','D','Dp'})
for ch in chunks:
    ck(len(ch&physical)==1,'chunk_ports','one permanent physical slot')
    ck(len(ch&aux)==1,'chunk_ports','one G slot')
# Exact rational examples test normalization with unrelated readout shares.
for d,a1,a2,c,a,b,z in itertools.product((F(2),F(3,5)),(F(1,3),),(F(2,5),),(F(1,7),),(F(2,9),),(F(3,11),),(F(4,13),)):
    e=d*d/(2*a1*a1*a2*a2); root=-e/(c*a*b*z)
    ck(c*a*b*z*root==-e,'normalization','negative root realizes negative cycle')
    ck((c*a*root)*(b*z)==-e,'normalization','conditional lambda times both side shifts')

# Each conditional family graph consists of two three-force chunks per B6.
# B6 has left/right side boundaries. Products B B^T alternate orientation.
def family_graph(kind,h):
    n=2*h if kind in ('trace','one') else 2*h+1
    nodes=[];edges=[];slots=[]
    for i in range(n):
        nodes.extend([(i,'l'),(i,'r')]);edges.append(((i,'l'),(i,'r')))
        slots.extend([(i,'l'),(i,'r')])
    for i in range(n-1):
        side='r' if i%2==0 else 'l'
        edges.append(((i,side),(i+1,side)))
    if kind=='trace': edges.append(((n-1,'l'),(0,'l')))
    else:
        left=('cap',0);right=('cap',1);nodes.extend([left,right])
        edges.extend([(left,(0,'l')),(right,(n-1,'l' if n%2==0 else 'r'))])
    return nodes,edges,slots

def pairings(xs):
    if not xs: yield (); return
    a=xs[0]
    for j in range(1,len(xs)):
        for rest in pairings(xs[1:j]+xs[j+1:]): yield ((a,xs[j]),)+rest

def component_count(nodes,edges):
    groups={n:n for n in nodes}
    def find(n):
        while groups[n]!=n: n=groups[n]
        return n
    for a,b in edges:
        a,b=find(a),find(b)
        if a!=b: groups[b]=a
    return len({find(n) for n in nodes})

for kind,h in (('trace',1),('trace',2),('one',1),('one',2),('both',0),('both',1)):
    nodes,edges,slots=family_graph(kind,h)
    ck(component_count(nodes,edges)==1,'wick_graphs','argument graph connected')
    ck(all(a!=b for a,b in edges),'wick_graphs','no initial middle/side self-loop')
    count=0
    for pairing in pairings(tuple(slots)):
        ck(all(a!=b for a,b in pairing),'wick_graphs','every Wick edge joins distinct physical chunks')
        ck(component_count(nodes,edges+list(pairing))==1,'wick_graphs','final conditional Wick coefficient connected')
        count+=1
    ck(count==normal(len(slots)),'wick_graphs','complete Wick enumeration')
# Joint cumulants of quadratic-G argument groups retain group-connected pairings.
connected_counts={}
for h in range(2,5):
    slots=tuple((i,j) for i in range(h) for j in range(2)); count=0
    for pairing in pairings(slots):
        if component_count(list(range(h)),[(a[0],b[0]) for a,b in pairing])==1:
            ck(all(a!=b for a,b in pairing),'wick_cumulants','no initial chunk self-loop')
            count+=1
    ck(count==2**(h-1)*factorial(h-1),'wick_cumulants','Gaussian quadratic connected-pairing count')
    connected_counts[h]=count

# One-Hilbert/all-cut numeric diagnostics, including a prohibited generic trace.
def maxcut(t):
    if t.ndim==1:return np.linalg.norm(t)
    ans=0.0
    for mask in itertools.product((False,True),repeat=t.ndim):
        if not any(mask) or all(mask):continue
        l=[i for i,v in enumerate(mask) if v];r=[i for i,v in enumerate(mask) if not v]
        a=t.transpose(l+r).reshape(np.prod([t.shape[i] for i in l]),-1)
        ans=max(ans,float(np.linalg.norm(a,2)))
    return ans
rng=np.random.default_rng(100526)
for d in (2,3):
    for case in range(5):
        b=rng.standard_normal((d,)*4);b/=maxcut(b)
        opened=np.einsum('ijab,klcb->ijklac',b,b)
        closed=np.einsum('ijklaa->ijkl',opened)
        ck(maxcut(opened)<1+1e-11,'tensor_cuts','opened proper cuts')
        ck(maxcut(closed)<1+1e-11,'tensor_cuts','closed two-spine cycle proper cuts')
        sym=(opened+opened.swapaxes(-1,-2))/2
        rhs=np.linalg.norm(closed)**2+2*np.linalg.norm(sym)**2
        lhs=np.einsum('ijklaa,ijklbb->',opened,opened)+np.einsum('ijklab,ijklab->',opened,opened)+np.einsum('ijklab,ijklba->',opened,opened)
        ck(abs(lhs-rhs)<1e-10,'hermite_energy','exact second-chaos identity')
        ck(rhs<=3*d+1e-10,'hermite_energy','single Hilbert scale')
    # Initial block with two internal loops is deliberately inadmissible.
    unit=np.zeros(d);unit[0]=1
    bad=np.einsum('i,ab,cd->iabcd',unit,np.eye(d),np.eye(d))/d
    trace=np.einsum('iaabb->i',bad)
    ck(maxcut(bad)<=1+1e-12,'negative_control','all proper cuts bounded')
    ck(abs(np.linalg.norm(trace)-d)<1e-12,'negative_control','two initial self-traces produce D')

# Same-node original endpoint clock margin, all derivative endpoints pessimistically
# charged to the same r2 shield; one further caller score is also retained.
for h in range(1,31):
    for k in range(h,h+5):
        margin=8*k-4*k-2*(h-1)
        ck(margin>=2*h+2,'same_node_clocks','bridge margin')
        ck(margin-1>0,'same_node_clocks','margin after caller score')
# Distinct-clock star: this deliberately fails and guards against aggregate reuse.
h=6; local_root_power=1; local_hit_count=h-1
ck(4*local_root_power-local_hit_count==-1,'negative_control','heterogeneous-clock star not covered by k>=h')
ck(4*1-(5-1)-1==-1,'negative_control','h=5 heterogeneous caller-score failure')

result={'status':'PASS_BOUNDED_FIXED_NODE_WITH_IMPORTED_NATIVE_GUARDS','assertions':checks,'categories':categories,'formal_force_cutoff':CUTOFF,'root_sign_series':series,'connected_G_pairing_counts':connected_counts,'negative_controls':['Bounded proper cuts do not control an initial tensor with two self-traces; norm becomes D.','Same-node aggregate clock margin does not certify distinct-clock star-hit histories.'],'not_certified':['Numerical instantiation of native selected-pair programs for an unspecified original g.','General all-order current closure or cross-node shared-bank closure.','Promotion of a completed law kernel to an original gradient source.']}
(HERE/'independent_auxiliary_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='root_sign_series'},indent=2))
