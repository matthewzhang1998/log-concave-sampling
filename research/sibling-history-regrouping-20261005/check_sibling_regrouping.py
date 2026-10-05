#!/usr/bin/env python3
"""Exact AST/cancellation tests and closed-form Gaussian checks; not a native executor."""
from pathlib import Path
from collections import Counter
from itertools import product
from math import comb, factorial, exp, sqrt
import importlib.util, json, hashlib
import sympy as s
HERE=Path(__file__).parent
path=Path('/workspace/shared/rank-indexed-positive-returns-20261005/generate_rank_terms.py')
spec=importlib.util.spec_from_file_location('generator',path)
gen=importlib.util.module_from_spec(spec);spec.loader.exec_module(gen)
def hit(ast,vertex,edge):
    if ast[0]=='g':return ('g',ast[1],ast[2]+(edge,)) if ast[1]==vertex else None
    if ast[0]=='r':
        v=hit(ast[2],vertex,edge)
        return None if v is None else ('r',ast[1]+1,v)
    a,b=hit(ast[1],vertex,edge),hit(ast[2],vertex,edge)
    if a is not None:return ('mul',a,ast[2])
    if b is not None:return ('mul',ast[1],b)
    return None

def packet(n):
    base=('r',2,('mul',('r',2,('g',1,(1,))),('r',2,('g',2,(1,)))))
    out=[(2,base)]
    for rank in range(3,n+1):
        new=[]
        for c,a in out:
            derivatives=gen.derivative(a,rank-1)
            for vertex in (1,2):
                da=hit(a,vertex,rank-1)
                assert da in derivatives
                new.append((rank*c,('r',rank,('mul',da,('r',2,('g',rank,(rank-1,)))))))
        out=new
    return out

def inspect(ast):
    prim={};struct=[]
    def rec(a):
        if a[0]=='mul':rec(a[1]);rec(a[2]);return
        assert a[0]=='r'
        if a[2][0]=='g':prim[a[2][1]]=(a[1],len(a[2][2]))
        else:struct.append(a[1]);rec(a[2])
    rec(ast);return prim,struct

def ast_checks():
    rows=[]
    for n in range(2,11):
        a=packet(n);pc=Counter()
        for c,tree in a:
            prim,struct=inspect(tree);p=prim[1][1]-1
            assert c==factorial(n)
            assert prim[1]==(p+2,p+1)
            assert prim[2]==(n-p,n-1-p)
            assert all(prim[i]==(2,1) for i in range(3,n+1))
            assert struct==[n]*(n-1)
            pc[p]+=1
        assert len(a)==2**(n-2)
        assert pc==Counter({p:comb(n-2,p) for p in range(n-1)})
        rows.append({'rank':n,'histories':len(a),'per_history_coefficient':factorial(n),'multiplicities':dict(pc)})
    return rows

def exact_leibniz():
    x,r,t=s.symbols('x r t');f=s.Function('f');h=s.Function('h')
    checks=[]
    for m in range(0,9):
        lhs=sum(s.binomial(m,p)*s.diff(f(r*x),x,p)*s.diff(h(t*x),x,m-p) for p in range(m+1))
        assert s.expand(lhs-s.diff(f(r*x)*h(t*x),x,m))==0
        checks.append(m)
    return checks

def gauss(poly,z):
    return s.expand(sum(c*(s.factorial2(k-1) if k else 1) for (k,),c in s.Poly(s.expand(poly),z).terms() if k%2==0))

def ibp_checks():
    x,G,Z=s.symbols('x G Z');q1=s.Rational(2,7);q2=-s.Rational(1,9)
    r1=s.Rational(2,5);r2=s.Rational(3,7);sig1=s.Rational(1,3);sig2=s.Rational(1,4)
    center=s.Rational(3,5);width=s.Rational(4,5)
    g=x+x**3/s.Integer(11)+x**5/s.Integer(17)
    def Smean(q,sig):return gauss(Z*(g.subs(x,q+sig*Z)-g.subs(x,q))/sig,Z)
    B1=Smean(r1*x+q1,sig1);B2=Smean(r2*x+q2,sig2)
    assert s.expand(B1-gauss(s.diff(g,x).subs(x,r1*x+q1+sig1*Z),Z))==0
    F=s.expand(B1*B2)
    out=[]
    for m in range(0,9):
        lhs=gauss(s.diff(F,x,m).subs(x,center+width*G),G)
        rhs=gauss(s.hermite_prob(m,G)*F.subs(x,center+width*G),G)/width**m
        assert s.expand(lhs-rhs)==0
        out.append({'order':m,'exact_identity':True,'nonzero':lhs!=0})
    return out

def cosine_mean(c,V):return exp(-.5*(c[0]**2*V+c[1]**2+c[2]**2))
def cosine_cov(c,d,V):
    return .5*(cosine_mean(tuple(a+b for a,b in zip(c,d)),V)+cosine_mean(tuple(a-b for a,b in zip(c,d)),V))-cosine_mean(c,V)*cosine_mean(d,V)

def variance_checks():
    a,b=.5,.25;v=.8
    limit=a**12*(a*a*b*b*exp(-1)*(1+exp(-1))+512*b**4*exp(-2))
    rows=[]
    phases=[(1,1,0),(1,0,1),(2,1,1)]
    for K in (4,8,16,32,64,128,256):
        r=sqrt(1-2/K**2);V=K*K*r*r*v
        # Dividing the fully regrouped old random coefficient by K^2.
        coeff=[-a**7*b*exp(-.5)*r**6,-a**7*b*exp(-.5)*r**6,-32*a**6*b*b*exp(-1)*r**6]
        variance=sum(c*d*cosine_cov(u,w,V) for c,u in zip(coeff,phases) for d,w in zip(coeff,phases))
        assert variance>0
        rows.append({'K':K,'variance_over_K4':variance,'ratio_to_limit':variance/limit})
    assert abs(rows[-1]['ratio_to_limit']-1)<.0002
    return {'limit':limit,'rows':rows,'difference_mode_coefficient_at_equal_clocks':0,
            'scope':'Common 8!, base primitive r1*r2 and fixed structural/leaf clock masses suppressed; base r1*r2 tends to one.'}

def source_constants():
    n=8;m=6
    return {'unweighted_L2_over_A8_sminus6':sqrt(factorial(m))*3**(n/2),
            'caller_L2_per_sum_inverse_sigma':2*sqrt(factorial(m))*3**((n-1)/2),
            'value_calls_per_clock_sample':2*n,
            'standard_roots_per_shared_topology_sample':(n-1)+n+n,
            'primitive_dr_inverse_shield_integral':3.141592653589793/sqrt(2),
            'primitive_rdr_inverse_shield_integral':sqrt(2),
            'primitive_importance_density':'r/sqrt(1-r^2)',
            'importance_weighted_caller_bound_per_factor':'2sqrt(2)|Z|',
            'multivariate_linear_fixture_HS_squared_over_a16_sminus12':'D * 6! * binomial(D+5,6)'}

if __name__=='__main__':
    report={'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'ast':ast_checks(),
            'leibniz_orders':exact_leibniz(),'polynomial_ibp':ibp_checks(),'group_variance':variance_checks(),
            'finite_source':source_constants()}
    (HERE/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'rank8_packet.json').write_text(json.dumps([{'coefficient':c,'ast':a} for c,a in packet(8)],indent=2)+'\n')
    print(json.dumps({'ast_ranks':len(report['ast']),'rank8_histories':64,'leibniz_checks':len(report['leibniz_orders']),
                      'ibp_checks':len(report['polynomial_ibp']),'regrouped_variance_limit':report['group_variance']['limit'],
                      'value_calls_per_clock_sample':16},indent=2))
