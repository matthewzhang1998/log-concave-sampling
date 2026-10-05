#!/usr/bin/env python3
"""Independent audit of serialized ASTs, exact star moments, and clock bookkeeping."""
from pathlib import Path
from collections import Counter
from functools import lru_cache
from math import factorial,comb
import json
import sympy as s
HERE=Path(__file__).parent
x,theta,b,u,v,lam=s.symbols('x theta b u v lambda')
checks=Counter()
def verify(test,group):
    assert bool(test),group
    checks[group]+=1
def freeze(a):
    return tuple(freeze(v) if isinstance(v,list) else v for v in a)
rows=json.loads((HERE.parent/'rank5_histories.json').read_text())
asts=[(row['coefficient'],freeze(row['ast'])) for row in rows]
@lru_cache(None)
def leaves(ast):
    if ast[0]=='g':return ((ast[1],ast[2]),)
    if ast[0]=='r':return leaves(ast[2])
    verify(ast[0]=='mul','AST_node_type')
    return leaves(ast[1])+leaves(ast[2])
def inspect(ast):
    if ast[0]=='g':return 0
    if ast[0]=='mul':return inspect(ast[1])+inspect(ast[2])
    local=leaves(ast[2]); edges=Counter(e for _,es in local for e in es)
    boundary=sum(k==1 for k in edges.values())
    verify(ast[1]==len(local)+boundary,'literal_clock_index_equals_subtree_rank_plus_boundary')
    return 1+inspect(ast[2])
shapes=Counter()
for coeff,ast in asts:
    ls=leaves(ast); verify(sorted(v for v,_ in ls)==list(range(1,6)),'disjoint_force_labels')
    incidence={}
    for vertex,edges in ls:
        for edge in edges:incidence.setdefault(edge,[]).append(vertex)
    verify(set(incidence)==set(range(1,5)) and all(len(vs)==2 for vs in incidence.values()),'spatial_edge_labels_exactly_twice')
    seen={1}
    while True:
        new=seen|{v for ends in incidence.values() if any(w in seen for w in ends) for v in ends}
        if new==seen:break
        seen=new
    verify(len(seen)==5,'spatial_graph_connected_tree')
    verify(inspect(ast)==9,'positive_clock_count')
    verify(coeff==factorial(5),'each_history_weight_rank_factorial')
    ds=tuple(sorted(len(es) for _,es in ls)); shapes[','.join(map(str,ds))]+=1
    verify(sum(d-1 for d in ds)==3,'extra_heat_derivative_count')
verify(len(asts)==272,'rank5_history_count')
verify(dict(shapes)=={'1,1,2,2,2':112,'1,1,1,2,3':144,'1,1,1,1,4':16},'rank5_degree_shapes')

# Independent OU resolvent in the Hermite eigenbasis, not monomial back-substitution.
def resolvent(poly,k):
    rem=s.Poly(s.expand(poly),x); out=0
    if rem.is_zero:return s.Integer(0)
    for j in range(rem.degree(),-1,-1):
        c=rem.nth(j); H=s.hermite_prob(j,x)
        out+=c*H/(k+j); rem=s.Poly(s.expand(rem.as_expr()-c*H),x)
    verify(rem.as_expr()==0,'Hermite_resolvent_reduction')
    return s.expand(out)
for h in [x+x*x/5,x+x*x/7+x**3/11]:
    @lru_cache(None)
    def evaluate(ast):
        if ast[0]=='g':return s.diff(h,x,len(ast[2]))
        if ast[0]=='r':return resolvent(evaluate(ast[2]),ast[1])
        return s.expand(evaluate(ast[1])*evaluate(ast[2]))
    M=[s.Integer(1)]; K=[s.Integer(0)]
    for n in range(1,6):
        M.append(n*resolvent(h*M[-1],n))
        K.append(s.expand(M[n]-sum(comb(n-1,j-1)*K[j]*M[n-j] for j in range(1,n))))
    generated=s.expand(sum(coef*evaluate(ast) for coef,ast in asts))
    verify(s.expand(generated-K[5])==0,'serialized_AST_scalar_moment_cumulant_identity')
    verify(s.expand(generated-resolvent(sum(comb(5,i)*s.diff(K[i],x)*s.diff(K[5-i],x) for i in range(1,5)),5))==0,'serialized_AST_connected_recurrence')

# Product Gaussian moments via a direct expansion, then log coefficients via partitions.
def moments_shift(k,a):
    ans=0
    for d in range(k+1):
        if d%2:continue
        gm=factorial(d)//(2**(d//2)*factorial(d//2))
        ans+=comb(k,d)*gm*(a*theta)**(k-d)
    return s.expand(ans)
def partitions(n):
    ans=[()]
    for a in range(n):
        nxt=[]
        for part in ans:
            nxt.append(part+((a,),))
            for j in range(len(part)):nxt.append(part[:j]+(part[j]+(a,),)+part[j+1:])
        ans=nxt
    return ans
mu={k:s.expand(moments_shift(k,b)*moments_shift(k,u)*moments_shift(k,v)) for k in range(7)}
cumulants={}
for k in range(1,7):
    value=0
    for pi in partitions(k):
        term=s.Integer((-1)**(len(pi)-1)*factorial(len(pi)-1))
        for block in pi:term*=mu[len(block)]
        value+=term
    cumulants[k]=s.expand(value)
    recurrence=s.expand(mu[k]-sum(comb(k-1,i-1)*cumulants[i]*mu[k-i] for i in range(1,k)))
    verify(cumulants[k]==recurrence,'star_full_partition_vs_recursive_cumulant')
second=s.expand(lam**2*theta**4*cumulants[2]/2)
expected=lam**2/s.Integer(2)*(theta**4+(b*b+u*u+v*v)*theta**6+(b*b*u*u+b*b*v*v+u*u*v*v)*theta**8)
verify(s.expand(second-expected)==0,'exact_rank5_second_star_coefficient')
verify(s.expand(second.subs({u:0,v:0})-lam**2*(theta**4+b*b*theta**6)/2)==0,'side_zero_fourth_and_sixth_currents')
verify(s.expand(second.subs(lam,-lam)-second)==0,'root_sign_even_feedback')
for n in range(4,10):
    aa=s.symbols('a:'+str(n-2))
    m1=s.prod(a*theta for a in aa)
    m2=s.prod(1+a*a*theta*theta for a in aa)
    generated=s.expand(lam**2*theta**4*(m2-m1*m1)/2)
    verify(s.Poly(generated,theta).coeff_monomial(theta**4)==lam**2/2,'arbitrary_star_fourth_coefficient')
# Six vertices, four tree edges, and three center-center edges.
edges=[(0,1),(1,2),(3,4),(4,5),(1,4),(1,4),(1,4)]
V=6; E=len(edges); beta=E-V+1; M=4; derivative_orders=[1,4,1,1,4,1]
verify((V,E,beta,M)==(6,7,2,4),'two_cycle_graph_invariants')
verify(sum(j-1 for j in derivative_orders)==M-2+2*beta==6,'two_cycle_derivative_mark_identity')
# Monotone original-gradient witness and its analytically smoothed fourth derivative.
A,k,z,sigma=s.symbols('A k z sigma',positive=True)
g=A*z/2+A*s.sin(k*z)/(4*k)
verify(s.simplify(s.diff(g,z,4)-A*k**3*s.sin(k*z)/4)==0,'original_gradient_fourth_derivative')
smoothed=sigma**3/A*s.diff(g,z,4)*s.exp(-sigma**2*k**2/2)
verify(s.simplify(smoothed-(sigma*k)**3*s.exp(-(sigma*k)**2/2)*s.sin(k*z)/4)==0,'original_gradient_C3_target')
# Finite counts and quadrature amplification are computed independently.
T={1:1}; S={1:1}
for n in range(2,9):
    T[n]=sum(i*(n-i)*T[i]*T[n-i] for i in range(1,n))
    S[n]=sum(comb(n,i)*i*(n-i)*S[i]*S[n-i] for i in range(1,n))
    verify(S[n]==factorial(n)*T[n],'positive_weight_factorial_identity')
B5=2**3*factorial(3); D5=9*2**9*S[5]*B5
out={'status':'PASS','exact_assertion_counts':dict(checks),'total_assertions':sum(checks.values()),
     'rank5_shapes':dict(shapes),'history_counts':T,'positive_weight_sums':S,'B5':B5,'D5':D5,
     'star_rank5_second_coefficient':str(s.factor(second)),
     'scope':'Exact serialized-history, finite formal cumulant, graph-label and derivative checks. No native MGF, native W2 lower bound, arbitrary-order closure, or arbitrary-graph admission is inferred.'}
(HERE/'tree_star_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
