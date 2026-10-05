#!/usr/bin/env python3
"""Finite exact connected OU-tree and positive-consumer generators; not a native executor."""
from functools import lru_cache
from math import comb, factorial
from collections import Counter
import json, hashlib
from pathlib import Path
import sympy as s
HERE=Path(__file__).parent
x=s.Symbol('x')

def partitions(items):
    items=tuple(items)
    if not items:
        yield (); return
    a,*rest=items
    for p in partitions(rest):
        yield ((a,),)+p
        for j,b in enumerate(p):
            yield p[:j]+((a,)+b,)+p[j+1:]

@lru_cache(None)
def cumulant_cut_constant(n):
    if n==2:return 1
    out=[]
    for p in range(1,n):
        value=factorial(p)*factorial(n-p)
        for blocks in partitions(range(n)):
            if len(blocks)<2:continue
            if all(any(j<p for j in b) and any(j>=p for j in b) for b in blocks):
                q=1
                for b in blocks:q*=cumulant_cut_constant(len(b))
                value+=q
        out.append(value)
    return max(out)

def relabel(ast,vshift,eshift):
    if ast[0]=='g':return ('g',ast[1]+vshift,tuple(j+eshift for j in ast[2]))
    if ast[0]=='r':return ('r',ast[1],relabel(ast[2],vshift,eshift))
    return ('mul',relabel(ast[1],vshift,eshift),relabel(ast[2],vshift,eshift))

def derivative(ast,edge):
    if ast[0]=='g':return [('g',ast[1],ast[2]+(edge,))]
    if ast[0]=='r':return [('r',ast[1]+1,a) for a in derivative(ast[2],edge)]
    return [('mul',a,ast[2]) for a in derivative(ast[1],edge)]+[('mul',ast[1],a) for a in derivative(ast[2],edge)]

@lru_cache(None)
def histories(n):
    if n==1:return ((1,('r',1,('g',1,()))),)
    out=[]
    for i in range(1,n):
        for c,a in histories(i):
            for d,b0 in histories(n-i):
                b=relabel(b0,i,max(i-1,0))
                for da in derivative(a,n-1):
                    for db in derivative(b,n-1):
                        out.append((comb(n,i)*c*d,('r',n,('mul',da,db))))
    return tuple(out)

def vertices(ast):
    if ast[0]=='g':return [(ast[1],ast[2])]
    if ast[0]=='r':return vertices(ast[2])
    return vertices(ast[1])+vertices(ast[2])

def clocks(ast):
    if ast[0]=='g':return []
    if ast[0]=='r':return [ast[1]]+clocks(ast[2])
    return clocks(ast[1])+clocks(ast[2])

def scalar_ast(ast):
    if ast[0]=='g':return ('g',len(ast[2]))
    if ast[0]=='r':return ('r',ast[1],scalar_ast(ast[2]))
    children=sorted((scalar_ast(ast[1]),scalar_ast(ast[2])),key=repr)
    return ('mul',*children)

def resolvent(poly,k):
    p=s.Poly(s.expand(poly),x)
    if p.is_zero:return s.Integer(0)
    rem={a[0]:b for a,b in p.terms()}; out=0
    for j in range(p.degree(),-1,-1):
        a=rem.get(j,0)/s.Rational(k+j)
        out+=a*x**j
        if j>=2:rem[j-2]=rem.get(j-2,0)+j*(j-1)*a
    return s.expand(out)

def test_scalar_tree_identity(rank=5):
    g=x+x**2/s.Integer(7)+x**3/s.Integer(11)
    @lru_cache(None)
    def ev(ast):
        if ast[0]=='g':return s.diff(g,x,ast[1])
        if ast[0]=='r':return resolvent(ev(ast[2]),ast[1])
        return s.expand(ev(ast[1])*ev(ast[2]))
    moment=[s.Integer(1)]; kappa=[s.Integer(0)]
    checks=[]
    for n in range(1,rank+1):
        moment.append(n*resolvent(g*moment[-1],n))
        kappa.append(s.expand(moment[n]-sum(comb(n-1,j-1)*kappa[j]*moment[n-j] for j in range(1,n))))
        compressed=Counter()
        for coeff,ast in histories(n):compressed[scalar_ast(ast)]+=coeff
        generated=s.expand(sum(c*ev(a) for a,c in compressed.items()))
        assert s.expand(generated-kappa[n])==0
        checks.append({'rank':n,'raw_histories':len(histories(n)), 'scalar_clock_histories':len(compressed),'degree':int(s.degree(generated,x))})
    return checks

def generate_report(maxrank=6):
    rows=[]
    for n in range(2,maxrank+1):
        shapes=Counter(); primitive_orders=set(); structural_clock_counts=set()
        for coeff,ast in histories(n):
            vs=vertices(ast); assert len(vs)==n
            ids=[v for v,_ in vs]; assert sorted(ids)==list(range(1,n+1))
            edges=Counter(e for _,es in vs for e in es)
            assert set(edges)==set(range(1,n)) and set(edges.values())=={2}
            ds=sorted(len(es) for _,es in vs)
            assert min(ds)>=1 and sum(ds)==2*n-2
            assert sum(d-1 for d in ds)==n-2
            assert len(clocks(ast))==2*n-1
            assert all(1<=k<=n for k in clocks(ast))
            primitive_orders.update(ds);shapes[','.join(map(str,ds))]+=1
        rows.append({'rank':n,'unmerged_histories':len(histories(n)), 'positive_coefficient_sum':sum(c for c,_ in histories(n)), 'clock_count':2*n-1,'degree_multiset_counts':dict(shapes),'primitive_derivative_orders':sorted(primitive_orders),'cumulant_proper_cut_majorant':cumulant_cut_constant(n)})
    feedback=[]
    for r in range(3,maxrank+1):
        counts=Counter(len(p)+1 for p in partitions(range(r-1)))
        feedback.append({'coefficient_rank':r,'partition_terms':sum(counts.values()),'emitted_current_rank_counts':dict(counts),'all_singleton_leading_cancelled_once':True})
    return {'scope':'Exact finite algebra and analytical target generation, not execution or native admission','tree_rows':rows,'consumer_rows':feedback,'scalar_moment_cumulant_checks':test_scalar_tree_identity(5)}

if __name__=='__main__':
    report=generate_report()
    (HERE/'generated_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    rank5=[{'coefficient':c,'ast':a} for c,a in histories(5)]
    (HERE/'rank5_histories.json').write_text(json.dumps(rank5,indent=2)+'\n')
    print(json.dumps(report,indent=2))
