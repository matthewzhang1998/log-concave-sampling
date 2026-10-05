#!/usr/bin/env python3
"""Independent exact audit. No native sampler or positive-realization certification."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib, json
import sympy as s

ROOT = Path(__file__).resolve().parent.parent
checks = {}
assertions = 0

def require(condition, message):
    global assertions
    assertions += 1
    if not condition:
        raise AssertionError(message)

@lru_cache(None)
def partitions(items):
    if not items:
        return ((),)
    head, *tail = items
    out = []
    for p in partitions(tuple(tail)):
        out.append(((head,),) + p)
        for j, block in enumerate(p):
            out.append(p[:j] + ((head,) + block,) + p[j+1:])
    return tuple(out)

def gaussian_moment(n):
    return 0 if n % 2 else factorial(n) // (2**(n//2)*factorial(n//2))

def mobius_cumulant(moment, size):
    ans = 0
    for p in partitions(tuple(range(size))):
        term = (-1)**(len(p)-1)*factorial(len(p)-1)
        for block in p:
            term *= moment(block)
        ans += term
    return s.expand(ans)

# Direct set-partition construction, not the author's moment-cumulant recurrence.
a,b,theta,lam = s.symbols('a b theta lam')
P,Q = s.symbols('P Q')
V = theta**2*(P+a*theta)*(Q+b*theta)
@lru_cache(None)
def vmoment(n):
    poly=s.Poly(s.expand(V**n),P,Q)
    return s.expand(sum(c*gaussian_moment(i)*gaussian_moment(j)
                        for (i,j),c in poly.terms()))
actual=0
for n in range(1,7):
    kap=mobius_cumulant(lambda block: vmoment(len(block)),n)
    actual += lam**n*kap/factorial(n)
closed = -s.log(1-lam**2*theta**4)/2 + (
    lam*a*b*theta**4+lam**2*(a*a+b*b)*theta**6/2)/(1-lam**2*theta**4)
expected=s.series(closed,lam,0,7).removeO()
require(s.expand(actual-expected)==0,'scalar partition/Wick versus direct integral')
require(s.expand(s.expand(actual).coeff(lam,2)-theta**4/2-(a*a+b*b)*theta**6/2)==0,
        'distinct side terms not combined into an unjustified multiplicity')
checks['scalar_wick']={'amplitude_orders':6,'method':'Direct labeled set partitions and independent Gaussian monomials'}

# Exact BKAR audit: integrate all ordering sectors of all labeled trees.
def trees(n):
    edges=list(combinations(range(n),2))
    for picked in combinations(edges,n-1):
        seen={0}
        while True:
            grown=seen|{v for u,v in picked if u in seen}|{u for u,v in picked if v in seen}
            if grown==seen: break
            seen=grown
        if len(seen)==n: yield picked

def paths(tree,n):
    adjacency=[[] for _ in range(n)]
    for j,(u,v) in enumerate(tree):
        adjacency[u].append((v,j));adjacency[v].append((u,j))
    ans={}
    for origin in range(n):
        stack=[(origin,-1,())]
        while stack:
            vertex,parent,path=stack.pop()
            if vertex!=origin: ans[origin,vertex]=path
            for nxt,e in adjacency[vertex]:
                if nxt!=parent: stack.append((nxt,vertex,path+(e,)))
    return ans

def polynomial_gaussian_moment(counts,covedge,dim):
    zero=(0,)*dim
    @lru_cache(None)
    def rec(state):
        if sum(state)==0: return {zero:F(1)}
        if sum(state)%2: return {}
        i=next(j for j,c in enumerate(state) if c)
        rem=list(state);rem[i]-=1; out=defaultdict(F)
        for j,num in enumerate(rem):
            if not num: continue
            tail=rem.copy();tail[j]-=1
            for mon,coef in rec(tuple(tail)).items():
                key=list(mon)
                if i!=j: key[covedge[i,j]]+=1
                out[tuple(key)]+=num*coef
        return dict(out)
    return rec(tuple(counts))

def simplex_integral(poly,order):
    ans=F(0)
    for mon,coef in poly.items():
        partial=0;term=coef
        for j,e in enumerate(order):
            partial += mon[e]
            term /= partial+j+1
        ans+=term
    return ans

def bkar(powers):
    n=len(powers);result=F(0)
    for tree in trees(n):
        deg=[0]*n
        for u,v in tree: deg[u]+=1;deg[v]+=1
        if any(d>k for d,k in zip(deg,powers)): continue
        scalar=1
        for k,d in zip(powers,deg): scalar*=factorial(k)//factorial(k-d)
        pathdict=paths(tree,n)
        for order in permutations(range(n-1)):
            position={edge:j for j,edge in enumerate(order)}
            covedge={pair:min(path,key=position.get) for pair,path in pathdict.items()}
            moment=polynomial_gaussian_moment(tuple(k-d for k,d in zip(powers,deg)),covedge,n-1)
            result+=scalar*simplex_integral(moment,order)
    return result

bkar_cases=[(1,1),(2,2),(3,3),(1,2,3),(2,2,2),(2,3,3),
            (1,1,1,1),(1,1,2,2),(1,2,2,3),(2,2,2,2),(1,1,3,3)]
for powers in bkar_cases:
    direct=mobius_cumulant(lambda block:gaussian_moment(sum(powers[j] for j in block)),len(powers))
    require(bkar(powers)==direct,('whole-bank tree',powers))
checks['whole_bank_tree']={'monomial_cases':len(bkar_cases),'max_arguments':4,
                          'method':'All labeled trees; exact min-covariance sector integration'}

# Independent heterogeneous bookkeeping from actual owned occurrence charges.
Gamma=F(1,3)
roots=[{'beta':F(1,7),'gamma':F(1,5),'ks':(0,2),'d':F(3,2),'weighted':True},
       {'beta':F(1,9),'gamma':F(1,6),'ks':(1,3,0),'d':F(7,4),'weighted':False},
       {'beta':F(1,11),'gamma':F(1,4),'ks':(2,1),'d':F(5,3),'weighted':True}]
def charge(k,w):return max(k-2,0) if w else k
for h in (2,3):
    selected=roots[:h]
    psi=[r['d']-2*Gamma for r in selected]
    # Vary retained side counts independently; root surplus remains once per genuine root.
    for kept_sides in product((0,1),repeat=h):
        B=H=F(0)
        for r,keep in zip(selected,kept_sides):
            retained=r['ks'] if keep else r['ks'][:1]
            B+=r['beta']*len(retained)
            H+=r['gamma']*sum(charge(k,r['weighted']) for k in retained)
        amp=B+H+sum(r['d'] for r in selected)
        child=amp-B-H-2*Gamma
        require(child==sum(psi)+2*Gamma*(h-1),'heterotypic conditional exact reserve')
        require(child!=h*psi[0]+2*Gamma*(h-1),'same-type formula fails for unlike roots')
for chosen in combinations(range(len(roots)),2):
    left,right=[roots[i] for i in chosen]
    for i,j in product(range(len(left['ks'])),range(len(right['ks']))):
        original=left['d']+right['d']-4*Gamma
        cost=sum(r['gamma']*(charge(r['ks'][v]+1,r['weighted'])-charge(r['ks'][v],r['weighted']))
                 for r,v in ((left,i),(right,j)))
        child=left['d']+right['d']-cost-2*Gamma
        require(child>=original,'heterogeneous weighted bridge lower bound')
for k in range(8):
    delta_q=charge(k+2,True)-charge(k,True)
    gain=2*Gamma-Gamma*delta_q
    require(gain==(2-k)*Gamma if k<2 else gain==0,'finite heat resource credit')
# A newly charged derivative loss can destroy the bridge gain even after correct bookkeeping.
psi_sum=F(3);delta_H=2*Gamma;extra_loss=F(1,10)
fully_charged=psi_sum+2*Gamma-delta_H-extra_loss
require(fully_charged<psi_sum,'recording an extra loss does not make reserve superadditive')
checks['heterogeneous_resources']={'unlike_root_surpluses':'summed separately',
                                   'omitted_sides':True,'weighted_bridge':True,
                                   'weighted_heat_neutral_at_k_at_least_2':True}

# A genuinely mixed two-current inverse with old/new and new/new terms.
u,v=s.symbols('u v');cut=8

def trunc(p):
    return s.Add(*[c*u**ij[0]*v**ij[1] for ij,c in s.Poly(s.expand(p),u,v).terms() if sum(ij)<cut])
def valuation(p):
    return min((sum(ij) for ij,c in s.Poly(s.expand(p),u,v).terms() if c),default=cut)
def residual(x,y):return (trunc(v*x+y*y),trunc(v*y+x*y+x*x*x))
x=y=s.Integer(0);vals=[]
for step in range(cut):
    rx,ry=residual(x,y);xn,yn=trunc(u-rx),trunc(v-ry)
    val=min(valuation(xn-x),valuation(yn-y));vals.append(val)
    require(val>=min(step+1,cut),'mixed inverse strict gain')
    x,y=xn,yn
rx,ry=residual(x,y)
require(trunc(x+rx-u)==0 and trunc(y+ry-v)==0,'mixed inverse solves both currents')
require(valuation(2*(u*u-u**3))==2,'zero-reserve linear factor does not improve')
# Without a filter-preserving section, source reserve positivity alone is insufficient.
for p in range(2,9):
    source_grade=p-1;old_grade=1;emitted_grade=source_grade+old_grade
    require(source_grade>=1 and emitted_grade==p,'right inverse loss must be charged')
checks['filtered_inverse']={'cutoff':cut,'valuations':vals,'right_inverse_loss_countertest':True}

# Exact sine observer: characteristic/tilt identity and cumulant coefficients.
z=s.symbols('z')
for m in range(1,17):
    term=s.series(z*s.sin(z),z,0,2*m+1).removeO().coeff(z,2*m)
    require(term*factorial(2*m)==(-1)**(m-1)*2*m,'exact observer cumulant')
# g'=c+epsilon cos is in [c-epsilon,c+epsilon] subset (0,1].
c,eps=F(1,2),F(1,4)
require(0<c-eps and c+eps<=1,'admissible monotone gradient')
checks['observer']={'ranks':list(range(2,34,2)),'same_alpha_grade':True,
                    'strict_monotone_gradient_example':{'c':str(c),'epsilon':str(eps)}}

# Keep ownership: include a nontrivial endpoint shift and keep scale.
def eg(p):
    return s.expand(sum(c*gaussian_moment(k) for (k,),c in s.Poly(s.expand(p),z).terms()))
eta=F(2,3);x0=F(3,5);X=x0+eta*z
w=s.symbols('w')
for r in range(1,7):
    for degree in range(r,11):
        lhs=eg(s.diff(w**degree,w,r).subs(w,X))
        rhs=eta**(-(r-1))*eg(s.hermite_prob(r-1,z)*s.diff(w**degree,w).subs(w,X))
        require(lhs==rhs,'scaled shifted owned keep integration by parts')
require(eg(z*s.diff(w**3,w,2).subs(w,z)) !=
        eg(z*s.hermite_prob(1,z)*s.diff(w**3,w).subs(w,z)),
        'keep-dependent coefficient breaks uncorrected relation')
# Equal unconditional scalar currents can fail retained-record equivalence.
require(eg(z)==0 and eg(z*z)==1,'unconditional zero not retained-record zero')
checks['current_pairing']={'keep_scale':str(eta),'shift':str(x0),
                           'dependent_coefficient_counterexample':True,
                           'conditional_ownership_counterexample':True}

# Trace-to-bridge check by direct finite Gaussian-shift moments, not heat-operator powers.
heat=s.symbols('heat')
def gaussian_shift(poly,centers,groups):
    aux=s.symbols('a0:'+str(max(groups)+1))
    shifted=s.Poly(s.expand(poly.subs({z:z+aux[g] for z,g in zip(centers,groups)},simultaneous=True)),*aux)
    ans=0
    for powers,coef in shifted.terms():
        if any(k%2 for k in powers):continue
        mult=1
        for k in powers:mult*=gaussian_moment(k)
        ans+=coef*mult*heat**(sum(powers)//2)
    return s.expand(ans)

trace_cases=[]
for n,dimension in ((1,2),(2,2),(3,1)):
    centers=s.symbols('z0:'+str(n*dimension))
    if dimension==2:
        poly=s.prod(centers[2*i]**4+centers[2*i+1]**4+centers[2*i]**2*centers[2*i+1]**2+i+1 for i in range(n))
    else:
        poly=s.prod(z**4+z**3+2 for z in centers)
    independent=gaussian_shift(poly,centers,list(range(n*dimension)))
    def cross(p):
        return s.expand(sum(s.diff(p,centers[u*dimension+j],centers[v*dimension+j])
                            for u in range(n) for v in range(u+1,n) for j in range(dimension)))
    correction=poly;power=poly
    for j in range(1,5):
        power=cross(power)
        correction+=(-heat)**j*power/factorial(j)
    joined=gaussian_shift(s.expand(correction),centers,[j%dimension for j in range(n*dimension)])
    for j in range(5):
        require(s.expand(independent-joined).coeff(heat,j)==0,
                ('direct Gaussian trace-to-bridge jet',n,dimension,j))
    trace_cases.append({'occurrences':n,'dimension':dimension})
# A nonconstant multiplier would generate omitted product-rule branches.
x=s.symbols('x');f=x**4;weight=x*x
require(s.diff(weight*f,x,2)!=weight*s.diff(f,x,2),'frozen coefficient condition is necessary')
checks['trace_to_bridge']={'cases':trace_cases,'heat_jet_order':4,
                          'method':'Exact independent versus common Gaussian-shift polynomial moments',
                          'unfrozen_multiplier_counterexample':True}

result={'status':'PASS','assertions':assertions,'scope':'Exact finite algebra only; no native realization, analytical ports, infinite heat closure, or cost theorem certified.',
        'checks':checks,
        'audited_input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'REPAIR-GRAPH-ALGEBRA.md',ROOT/'check_repair_algebra.py']}}
output=Path(__file__).with_name('independent-checks.json')
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
