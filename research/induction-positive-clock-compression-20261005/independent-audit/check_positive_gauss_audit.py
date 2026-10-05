#!/usr/bin/env python3
"""Independent finite diagnostics; analytical proofs are in the accompanying audit.
No original VALUE oracle or native source/filter is evaluated by these tests.
"""
import hashlib, json, math, pathlib, importlib.util
from collections import deque
from fractions import Fraction
import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=pathlib.Path(__file__).with_name('positive-gauss-audit-checks.json')
INPUT=pathlib.Path('/workspace/shared/induction-multicopy-heat-20261005/MANIFEST.json')
PIN='09ef144321423e22831c3cc720c3b8e3c48cc2bd090c891bc675d10ca7841574'
checks={}
def record(name, **data): checks[name]={'pass':True,**data}
assert hashlib.sha256(INPUT.read_bytes()).hexdigest()==PIN
record('sealed_input_pin',sha256=PIN)

# Exact squared inequalities used for the enlarged local and total derivative constants.
for k in range(257):
    assert (k+1)**k <= 4**k*math.factorial(k)
for K in range(33):
    for m in range(65):
        assert math.factorial(K+2*m) <= 2**(K+4*m)*math.factorial(K)*math.factorial(m)**2
record('factorial_envelopes',local_orders_checked=257,K_range=[0,32],m_range=[0,64])

# Every original force has a retained mark in the two audited fixed families.
def orientation_check(copy_sizes):
    offsets=[]; total=0
    for n in copy_sizes: offsets.append(total);total+=n
    edges=[]
    for o,n in zip(offsets,copy_sizes):
        # A six-force packet is two cubic paths joined by an old bridge.
        if n==3: edges += [(o,o+1),(o+1,o+2)]
        else: edges += [(o,o+1),(o+1,o+2),(o+3,o+4),(o+4,o+5),(o+1,o+4)]
    edges += [(offsets[0]+1,offsets[1]+1),(offsets[1]+1,offsets[2]+1)]
    adj=[[] for _ in range(total)]
    for u,v in edges: adj[u].append(v);adj[v].append(u)
    parent=[None]*total;parent[0]=-1;order=[0]
    for v in order:
        for w in adj[v]:
            if parent[w] is None: parent[w]=v;order.append(w)
    assert len(order)==total and len(edges)==total-1
    weights=[1 << (1 << j) for j in range(total)]
    tested=0
    for mask in range(1,(1<<total)-1):
        si=sum(w for j,w in enumerate(weights) if mask>>j&1)
        so=sum(weights)-si
        signed=[w*so if mask>>j&1 else -w*si for j,w in enumerate(weights)]
        sub=signed[:]
        oriented=[]
        for v in reversed(order[1:]):
            assert sub[v]!=0
            p=parent[v]
            oriented.append((v,p) if sub[v]>0 else (p,v))
            sub[p]+=sub[v]
        assert sub[0]==0
        indeg=[0]*total;outdeg=[0]*total;nexts=[[] for _ in range(total)]
        for u,v in oriented: indeg[v]+=1;outdeg[u]+=1;nexts[u].append(v)
        # The terminal leg supplies an inward leg for an input mark and outward for an output.
        for v in range(total):
            assert indeg[v]+((mask>>v)&1)>0
            assert outdeg[v]+(not ((mask>>v)&1))>0
        pending=indeg[:]; queue=deque(v for v in range(total) if pending[v]==0);topo=[]
        while queue:
            v=queue.popleft();topo.append(v)
            for w in nexts[v]:
                pending[w]-=1
                if pending[w]==0: queue.append(w)
        assert len(topo)==total
        rank={v:j for j,v in enumerate(topo)}
        # All possible Price inter-replica edges can be added simultaneously in forward order.
        for i,ni in enumerate(copy_sizes):
            for j,nj in enumerate(copy_sizes):
                if i>=j: continue
                for vi in range(offsets[i],offsets[i]+ni):
                    for vj in range(offsets[j],offsets[j]+nj):
                        u,v=(vi,vj) if rank[vi]<rank[vj] else (vj,vi)
                        assert rank[u]<rank[v]
        tested+=1
    return tested
record('augmented_graph_cut_orientation',triple_cubic_partitions=orientation_check([3,3,3]),mixed_partitions=orientation_check([3,3,6]))

# Dyadic gap panels, uniform subdivisions, and exact count comparison.
def intervals(J,h):
    result=[]
    for j in range(J):
        lo=2.0**(-j-1);hi=2.0**(-j)
        count=math.ceil((hi-lo)/h)
        for k in range(count):
            # Convert x=1-r gap intervals into r intervals.
            a=lo+(hi-lo)*k/count;b=lo+(hi-lo)*(k+1)/count
            result.append((1-b,1-a))
    return result

def rule(J,h,m):
    z,w=leggauss(m); nodes=[]
    for a,b in intervals(J,h):
        for zi,wi in zip(z,w): nodes.append(((a+b)/2+(b-a)*zi/2,(b-a)*wi/2,a,b))
    return nodes
max_weight_gap_ratio=0.; max_gap_variation=0.; cases=0
for J in range(1,9):
    for h in [.5,.2,.07,.03]:
        panels=intervals(J,h)
        assert len(panels)<=math.ceil(1/h)+J
        for a,b in panels:
            assert b-a<=1-b+1e-14
        for m in [1,2,3,5,8]:
            nodes=rule(J,h,m)
            L=1-2.0**(-J)
            assert all(lam>0 and 0<r<1 for r,lam,_,_ in nodes)
            assert abs(sum(lam for _,lam,_,_ in nodes)-L)<1e-13
            assert abs(sum(r*lam for r,lam,_,_ in nodes)-L*L/2)<1e-13
            for r,lr,a,b in nodes:
                for s,ls,c,d in nodes:
                    x=1-r;y=1-r*s
                    ratio=r*lr*ls/(x*y)
                    assert ratio<=1+1e-12
                    max_weight_gap_ratio=max(max_weight_gap_ratio,ratio)
                    ymin=1-b*d;ymax=1-a*c
                    variation=ymax/ymin
                    assert variation<=2+1e-12
                    max_gap_variation=max(max_gap_variation,variation)
            cases+=1
record('positive_weight_geometry',cases=cases,max_weight_to_gap_product_ratio=max_weight_gap_ratio,max_small_edge_gap_variation=max_gap_variation)

# Independent exact Gauss remainder coefficient and safe budget algebra.
for m in range(1,65):
    rho=Fraction(math.factorial(m)**4,(2*m+1)*math.factorial(2*m)**2)
    assert rho<=1
    # Two sectors times two pure-coordinate telescoping errors.
    quad=8*Fraction(1,4)**m
    eps=32*Fraction(1,4)**m
    assert quad<=eps/4
record('gauss_remainder_and_error_budget',m_range=[1,64])

# Smooth sector integrand F=exp(a min(t1,t2)+b max(t1,t2)); exact sector integral.
def exp_integral(c,L):return math.expm1(c*L)/c if c else L
num=[]
a=.7;b=.3
for m in [1,2,3,4,6]:
    J=7;h=.2;L=1-2.0**(-J);nodes=rule(J,h,m)
    approximate=2*sum(lr*ls*r*math.exp(a*r*s+b*r) for r,lr,_,_ in nodes for s,ls,_,_ in nodes)
    exact_truncated=2/a*(exp_integral(b+a*L,L)-exp_integral(b,L))
    exact_full=2/a*(exp_integral(b+a,1)-exp_integral(b,1))
    error=abs(approximate-exact_truncated)
    B=math.exp(a+b);A=1
    assert error <=8*B*4.0**(-m)+2e-14
    assert 0<=exact_full-exact_truncated<=4*B*2.0**(-J)
    num.append({'m':m,'Q_two_sectors':2*len(nodes)**2,'truncated_error':error,'actual_removed_integral':exact_full-exact_truncated})
record('analytic_sector_integral',cases=num)

# Final exact-parameter and adaptive arithmetic workflow, independently exercised.
spec=importlib.util.spec_from_file_location('author_rule',ROOT/'positive_sector_gauss.py')
author=importlib.util.module_from_spec(spec);spec.loader.exec_module(author)
parameter_cases=0
for N in [1,9,12]:
    for K in [0,1,7,8,10,11]:
        for eta in [Fraction(1,2),Fraction(3,20),Fraction(1,128)]:
            W=Fraction(7,3)
            BB=author.exact_B(N,K,eta,W)
            required_squared=W*W*2**(2*N+3*K)*math.factorial(K)*eta**(-K)
            assert BB*BB>=required_squared
            for epsilon in [Fraction(1,2),Fraction(1,10**8)]:
                pars=author.exact_parameters(N,K,Fraction(5,2),eta,epsilon,W,11*W)
                B=pars['B_star_rational'];d=pars['node_tolerance'];e=pars['weight_tolerance']
                assert Fraction(1,2**pars['J'])<=epsilon/(32*B)
                assert Fraction(1,4**pars['m'])<=epsilon/(32*B)
                assert 64*B*pars['P']/eta*d+2*B*e<=epsilon*Fraction(33,256)
                parameter_cases+=1
for j in range(65):
    assert author.ceil_log2_positive(Fraction(2**j))==j
    assert author.ceil_log2_positive(Fraction(2**j)+Fraction(1,7))==j+1
    if j: assert author.ceil_log2_positive(Fraction(2**j)-Fraction(1,7))==j
adaptive=[]
for n in [2,7,16]:
    floor=Fraction(1,10**40)
    ar=author.certified_legendre_to_tolerance(n,floor,floor)
    assert ar['node_error']<=floor and ar['weight_l1_error']<=floor
    adaptive.append({'degree':n,'bits':ar['certification_bits']})
pars,base,panels=author.prepare_certified_rule(2,1,1,Fraction(1,2),Fraction(1,10),caller_weight=Fraction(3))
assert len(panels)<=pars['interval_bound_per_axis']
assert base['node_error']<=pars['node_tolerance']
assert base['weight_l1_error']<=pars['weight_tolerance']
record('exact_parameters_and_adaptive_precision',cases=parameter_cases,adaptive_cases=adaptive,prepared_degree=pars['m'],prepared_interval_count=len(panels))

# Pin current authored proof files, if present. Review scope is explicit in the prose.
reviewed={}
for p in sorted(ROOT.glob('*')):
    if p.is_file(): reviewed[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
output={'status':'PASS_FINITE_DIAGNOSTICS_NOT_NATIVE_EXECUTION','checks':checks,'authored_file_hashes_at_run':reviewed}
OUT.write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
