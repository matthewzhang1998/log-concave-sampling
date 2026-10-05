#!/usr/bin/env python3
"""Exact finite-jet and finite-graph diagnostics, not a native sampler implementation."""
from fractions import Fraction as F
from math import factorial, comb
from pathlib import Path
import itertools, json, hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
counts={}
def ck(test,category,why):
    if not test: raise AssertionError(category+': '+why)
    counts[category]=counts.get(category,0)+1

def normal(n):
    return 0 if n%2 else factorial(n)//(2**(n//2)*factorial(n//2))

# Grades: (alpha,theta,t,G,H). All arithmetic is exact rational.
CUT=32
ZERO=(0,0,0,0,0);ONE={ZERO:F(1)}
def add(a,b):
    c=a.copy()
    for k,v in b.items():
        c[k]=c.get(k,F(0))+v
        if not c[k]: del c[k]
    return c
def scale(a,v):return {k:x*v for k,x in a.items() if x*v}
def mul(a,b):
    c={}
    for x,u in a.items():
        for y,v in b.items():
            z=tuple(i+j for i,j in zip(x,y))
            if z[0]<=CUT:c[z]=c.get(z,F(0))+u*v
    return {k:v for k,v in c.items() if v}
def exp0(a):
    out=ONE.copy();term=ONE.copy()
    for n in range(1,CUT+1):
        term=scale(mul(term,a),F(1,n))
        if not term:break
        out=add(out,term)
    return out
def log1(a):
    z=add(a,{ZERO:F(-1)});out={};term=ONE.copy()
    for n in range(1,CUT+1):
        term=mul(term,z)
        if not term:break
        out=add(out,scale(term,F((-1)**(n+1),n)))
    return out
def avg(a):
    out={}
    for (a0,th,t,g,h),v in a.items():
        key=(a0,th,t,0,0);out[key]=out.get(key,F(0))+v*normal(g)*normal(h)
    return {k:v for k,v in out.items() if v}
def direct(s,p,q):
    out=ONE.copy()
    for n in range(1,CUT//4+1):
        for i in range(0,n+1,2):
            for j in range(0,n+1,2):
                grade=6*n-i-j
                if grade>CUT:continue
                c=F(s**n*normal(2*n)**2*comb(n,i)*normal(i)*p**(n-i)*comb(n,j)*normal(j)*q**(n-j),factorial(n))
                out=add(out,{(grade,4*n-i-j,n,0,0):c})
    return log1(out)
def conditional(s,p,q):
    out={}
    for h in range(1,CUT//8+1):
        out=add(out,{(8*h,4*h,2*h,4*h,4*h):F(s**(2*h),2*h)})
        if 8*h+2<=CUT:out=add(out,{(8*h+2,4*h+2,2*h,4*h,4*h):F(s**(2*h)*(p*p+q*q),2)})
    for h in range(CUT//8+1):
        if 8*h+6<=CUT:out=add(out,{(8*h+6,4*h+4,2*h+1,4*h+2,4*h+2):F(s**(2*h+1)*p*q)})
    return out
series={}
for s,p,q in itertools.product((-2,-1,1,2),(-2,1,3),(-1,1,2)):
    a=direct(s,p,q);b=log1(avg(exp0(conditional(s,p,q))))
    for key in set(a)|set(b):ck(a.get(key,0)==b.get(key,0),'exact_symbol','direct/nested equality')
    for key,value in [((6,4,1,0,0),s*p*q),((8,4,2,0,0),F(9,2)*s*s),((10,6,2,0,0),F(9,2)*s*s*(p*p+q*q)),((12,8,2,0,0),4*s*s*p*p*q*q),((14,8,3,0,0),333*s**3*p*q)]:
        ck(a.get(key,0)==value,'offspring','leading or intrinsic coefficient')
    if (s,p,q) in ((1,1,1),(-1,1,1)):series[str(s)]={str(k):str(v) for k,v in sorted(a.items())}

# Literal opened and closed graphs.
V=('L','C','R','M','D','N');physical={'L','R','M','N'}
E=[('L','C'),('C','R'),('M','D'),('D','N'),('C','D')]
ck(len(E)==len(V)-1,'geometry','opened tree')
ck(len(E)+2-len(V)+1==2,'geometry','two closed cycles')
for v in V:
    degree=sum(v in e for e in E)
    auxiliary=2 if v in ('C','D') else 0
    derivative=degree+int(v in physical)+auxiliary-1
    ck(derivative==(4 if auxiliary else 1),'source_ports','C3 centers / C0 leaves')
ck(sum((4 if v in ('C','D') else 1)-1 for v in V)==4-2+2*2,'geometry','force/mark/cycle identity')
for block in ({'L','C'},{'D','M'}):
    ck(len(block&physical)==1,'chunk_ports','one permanent physical')
    ck(len(block&{'C','D'})==1,'chunk_ports','one G and one H occurrence')

def pairings(xs):
    if not xs:yield ();return
    for j in range(1,len(xs)):
        for rest in pairings(xs[1:j]+xs[j+1:]):yield ((xs[0],xs[j]),)+rest

def components(nodes,edges):
    par={x:x for x in nodes}
    def find(x):
        while par[x]!=x:x=par[x]
        return x
    for a,b in edges:
        aa,bb=find(a),find(b)
        if aa!=bb:par[bb]=aa
    return len({find(x) for x in nodes})

# Conditional ring graph, then all two-color auxiliary Wick pairings.
for k in (2,4):
    nodes=tuple(range(2*k));edges=[(2*i,2*i+1) for i in range(k)]
    edges += [(2*i+1,2*((i+1)%k)) for i in range(k)]
    ps=list(pairings(nodes))
    ck(len(ps)==normal(2*k),'wick_graphs','all pairings enumerated')
    for pg,ph in itertools.product(ps,repeat=2):
        ck(all(a!=b for a,b in pg+ph),'wick_graphs','no initial auxiliary self-loop')
        ck(components(nodes,edges+list(pg+ph))==1,'wick_graphs','connected with all permanent physical slots')
# Two-color cumulants of main G²H²: connected pairings across argument groups.
conn={}
for n in (2,3,4):
    nodes=tuple((i,j) for i in range(n) for j in range(2));ps=list(pairings(nodes));total=0
    for pg,ph in itertools.product(ps,repeat=2):
        es=[(a[0],b[0]) for a,b in pg+ph]
        if components(tuple(range(n)),es)==1:total+=1
    conn[n]=total
    moment=normal(2*n)**2
    # Scalar moment-cumulant recursion independent of connectivity test.
    mus=[normal(2*i)**2 for i in range(n+1)];kap=[0]*(n+1)
    for k in range(1,n+1):kap[k]=mus[k]-sum(comb(k-1,j-1)*kap[j]*mus[k-j] for j in range(1,k))
    ck(total==kap[n],'wick_cumulants','connected two-color Wick count equals cumulant')

# Uniform all-cut negative control after prohibited initial self-traces.
def maxcut(t):
    ans=0.
    for mask in itertools.product((False,True),repeat=t.ndim):
        if not any(mask) or all(mask):continue
        l=[i for i,x in enumerate(mask) if x];r=[i for i,x in enumerate(mask) if not x]
        a=t.transpose(l+r).reshape(int(np.prod([t.shape[i] for i in l])),-1)
        ans=max(ans,float(np.linalg.norm(a,2)))
    return ans
for d in (2,3,4):
    e=np.zeros(d);e[0]=1
    bad=np.einsum('i,ab,cd->iabcd',e,np.eye(d),np.eye(d))/d
    ck(maxcut(bad)<=1+1e-12,'negative_controls','all proper cuts bounded by one')
    ck(abs(np.linalg.norm(np.einsum('iaabb->i',bad))-d)<1e-12,'negative_controls','forbidden trace has D loss')
# Random two chunks closed along G,H and middle; physical and side slots survive.
rng=np.random.default_rng(10052603)
for d in (2,3):
    for _ in range(3):
        a=rng.normal(size=(d,)*5);a/=maxcut(a)
        b=rng.normal(size=(d,)*5);b/=maxcut(b)
        closed=np.einsum('iaghr,jbghr->ijab',a,b)
        ck(maxcut(closed)<=1+1e-10,'tensor_cuts','close three parallel edges across marked chunks')
        ck(np.linalg.norm(closed)<=np.sqrt(d)+1e-10,'tensor_cuts','single Hilbert factor')

# Dyadic positive endpoint checks; sigma²=Delta+tau², w=Delta.
clock_rows=[]
for tau in (2.**(-j) for j in range(2,11)):
    ds=np.array([2.**(-j) for j in range(0,90)]);ss=np.sqrt(ds+tau*tau)
    z=ds/ss**3;t=z*z
    row={'tau':tau,'max_z_tau':float(z.max()*tau),'sum_z_tau':float(z.sum()*tau),
         'sum_Ccaller_tau2':float((z/ss).sum()*tau**2),
         'sum_Dcaller_tau3':float((t/ss).sum()*tau**3),
         'sum_t_tau2':float(t.sum()*tau**2),'sum_t2_tau4':float((t*t).sum()*tau**4)}
    for name,v in row.items():
        if name!='tau':ck(v<10,'clock_bounds',name)
    clock_rows.append(row)
# Balanced source amplitudes reproduce t exactly and meet guards at tau=alpha^(1/3).
for alpha in (F(1,64),F(1,729),F(1,4096)):
    for wc,rest,sig in itertools.product((F(1,100),F(1,7)),(F(1,3),F(3,7)),(F(1,3),F(2,3))):
        z=wc/sig**3
        amp={'L':alpha*rest**2,'C':-alpha*z,'D':alpha*z,'R':alpha,'M':alpha,'N':alpha}
        p=F(1)
        for r in amp.values():p*=r
        ck(p==-alpha**6*(wc*rest)**2/sig**6,'normalization','balanced-center product')
# Physical C0 target exposes the lost root-caller inverse shield.
for sig in (.2,.1,.05):
    k=1/sig;q=np.pi/(2*k)
    Mprime=-.25*k*np.exp(-.5*(sig*k)**2)*np.sin(k*q)
    ck(abs(abs(Mprime)*sig-.25*np.exp(-.5))<1e-12,'negative_controls','admissible root coefficient caller scales sigma^-1')

sources=[Path('/workspace/shared/rank-indexed-positive-returns-20261005/TREE-GENERATOR-AND-RANK5-RETURN-TEST.md'),Path('/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/rank3-continuation/connected-current-return/AUXILIARY-PUBLIC-CYCLE-CLOSURE.md'),Path('/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/rank3-continuation/connected-current-return/CONDITIONAL-QUARTIC-PACKET-AND-CYCLE-RETURN.md'),Path('/workspace/scratch/e90bd698ae04/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex')]
sources.append(HERE/'REVIEWED-TWO-CYCLE-SOURCE.md')
result={'checks':sum(counts.values()),'categories':counts,'scalar_series':series,'connected_two_color_pairings':conn,'dyadic_rows':clock_rows,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
(HERE/'audit_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'checks':result['checks'],'categories':counts,'connected_two_color_pairings':conn},indent=2))

# Next-graph obstruction is deliberately scoped to long distinguished spines.
# Centers a,b,c,d, doubled simple cycle ab,bd,dc,ca. Four physical pendant leaves.
simple=[('a','b'),('b','d'),('d','c'),('c','a')]
multiedges=[(a,b,j) for a,b in simple for j in (0,1)]
next_rows=[];tree_count=0;path_pass=path_fail=longest_fail=0
for selected in itertools.combinations(multiedges,3):
    if components(tuple('abcd'),[(a,b) for a,b,_ in selected])!=1:continue
    tree_count+=1
    adj={a:[] for a in 'abcd'}
    for a,b,j in selected:adj[a].append(b);adj[b].append(a)
    ck(sorted(map(len,adj.values()))==[1,1,2,2],'next_obstruction','every spanning center tree is P4')
    for root in 'abcd':
        paths={root:[root]};todo=[root]
        for a in todo:
            for b in adj[a]:
                if b not in paths:paths[b]=paths[a]+[b];todo.append(b)
        farthest=max(len(p) for p in paths.values())
        for leaf,p in paths.items():
            if leaf==root:continue
            # Every used center edge retains a parallel auxiliary-cut twin.
            ck(all(any({a,b}=={p[i],p[i+1]} and (a,b,j) not in selected for a,b,j in multiedges) for i in range(len(p)-1)), 'next_obstruction','each kept edge has a cut twin')
            # Only endpoint centers have surviving physical slots in the zero-side spine.
            # Connected admissible chunks cannot contain adjacent centers.
            feasible=len(p)<=2
            if feasible:path_pass+=1
            else:path_fail+=1
            if len(p)==farthest:
                ck(not feasible,'next_obstruction','every longest rooted physical spine fails this certificate')
                longest_fail+=1
            next_rows.append({'tree':selected,'root':root,'terminal':leaf,'centers':len(p),'sufficient_chunk_certificate':feasible,'longest_from_root':len(p)==farthest})
ck(tree_count==32,'next_obstruction','all multigraph spanning trees')
ck(path_pass==192 and path_fail==192,'next_obstruction','short and long rooted paths both occur')
ck(longest_fail==128,'next_obstruction','every longest root path fails, all 32 trees and 4 roots')
result['next_graph']={'spanning_trees':tree_count,'all_directed_physical_paths':path_pass+path_fail,'short_paths_passing_chunk_test':path_pass,'long_paths_failing_chunk_test':path_fail,'longest_root_paths_failing':longest_fail,'force_count':8,'physical_marks':4,'cycle_count':5}
result['checks']=sum(counts.values());result['categories']=counts
(HERE/'audit_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('Including next graph:',json.dumps({'checks':result['checks'],'next_graph':result['next_graph']},indent=2))
for row in next_rows:
    ell=row['centers']; force=2*ell+4
    # Root factor hatomega^4 squares to hatomega^8 in every variance.
    for gamma in (F(1,3),F(1,2)):
        increase=F(force)-2*ell*gamma-(8-4*gamma)
        ck(increase==(2*ell-4)*(1-gamma),'next_grade_tradeoff','exact effective-grade difference')
        ck((increase==0)==row['sufficient_chunk_certificate'],'next_grade_tradeoff','short same-grade versus long uncertified chunks')
result['checks']=sum(counts.values());result['categories']=counts
(HERE/'audit_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('Including short-spine tradeoff:',result['checks'],'assertions')
