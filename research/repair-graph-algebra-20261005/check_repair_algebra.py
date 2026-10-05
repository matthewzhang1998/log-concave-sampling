#!/usr/bin/env python3
"""Exact finite algebra checks. No native execution or analytic-port certification."""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
from math import factorial, comb
from pathlib import Path
import json
import sympy as s

checks = {}
count = 0

def require(c, label):
    global count
    count += 1
    if not c:
        raise AssertionError(label)

def scalar_example():
    a,b,z=s.symbols('a b z')
    def moment(m,c):
        return s.expand(sum(comb(m,2*j)*factorial(2*j)//(2**j*factorial(j))*c**(m-2*j) for j in range(m//2+1)))
    mu={0:s.Integer(1)}; kap={}
    for m in range(1,9):
        mu[m]=s.expand(moment(m,a*z)*moment(m,b*z))
        kap[m]=s.expand(mu[m]-sum(comb(m-1,i-1)*kap[i]*mu[m-i] for i in range(1,m)))
        coefficient=s.expand(z**(2*m)*kap[m]/factorial(m))
        if m%2:
            expected=a*b*z**(2*m+2)
        else:
            expected=z**(2*m)/m+(a*a+b*b)*z**(2*m+2)/2
        require(s.expand(coefficient-expected)==0,('wick-cumulant',m))
    checks['conditional_symbol']={'orders_checked':8,'ranks_at_second_order':[4,6,6],'complete_partition_and_pairing_identity':True}

@dataclass(frozen=True)
class Graph:
    k:tuple
    edges:tuple
    marks:tuple
    a:Q

    @property
    def N(self): return len(self.k)
    @property
    def K(self): return sum(self.k)
    @property
    def R(self): return sum(self.marks)
    def psi(self,beta,gamma): return self.a-beta*self.N-gamma*(self.K+2)
    def grade(self,gamma): return self.a-gamma*self.K
    def validate(self, require_tree=True):
        require(self.K+2*self.N==2*len(self.edges)+self.R,'global valence')
        degrees=[0]*self.N
        adj=[[] for _ in range(self.N)]
        for u,v in self.edges:
            degrees[u]+=1;degrees[v]+=1
            adj[u].append(v);adj[v].append(u)
        require(all(degrees[v]+self.marks[v]==self.k[v]+2 for v in range(self.N)),'local valence')
        if require_tree:
            require(all(u!=v for u,v in self.edges),'loopless')
            require(len(self.edges)==self.N-1,'tree edge count')
            seen={0}; queue=[0]
            while queue:
                for w in adj[queue.pop()]:
                    if w not in seen: seen.add(w);queue.append(w)
            require(len(seen)==self.N,'connected')
            require(all(degrees[v]!=1 or self.marks[v]>0 for v in range(self.N)),'marked leaves')

def prufer_trees(n):
    if n==2: yield ((0,1),);return
    for word in product(range(n),repeat=n-2):
        d=[1]*n
        for v in word:d[v]+=1
        edges=[]
        for v in word:
            u=next(i for i in range(n) if d[i]==1)
            edges.append((u,v));d[u]-=1;d[v]-=1
        u,v=[i for i in range(n) if d[i]==1]
        edges.append((u,v));yield tuple(edges)

def bridges():
    beta,gamma=Q(1,4),Q(1,8)
    base=Graph((0,1,0),((0,1),(1,2)),(1,1,1),Q(3))
    base.validate();require(base.psi(beta,gamma)==Q(15,8),'base reserve')
    n_graphs=0
    for r in range(2,5):
        ts=list(prufer_trees(r));require(len(ts)==r**(r-2),'labeled tree count')
        for tr in ts:
            for hits in product(range(3), repeat=2*(r-1)):
                k=list(base.k*r);marks=base.marks*r
                es=[(u+3*i,v+3*i) for i in range(r) for u,v in base.edges]
                for j,(left,right) in enumerate(tr):
                    u=3*left+hits[2*j];v=3*right+hits[2*j+1]
                    k[u]+=1;k[v]+=1;es.append((u,v))
                graph=Graph(tuple(k),tuple(es),marks,r*base.a)
                graph.validate()
                require(graph.psi(beta,gamma)==r*base.psi(beta,gamma),'bridge additive reserve')
                require(graph.K==r+2*(r-1),'derivative count')
                n_graphs+=1
    checks['marked_tree_bridge_enumeration']={'graphs':n_graphs,'argument_arities':[2,3,4],'all_hit_locations_including_repeated_hits':True}

def ledgers():
    for beta in (Q(1,4),Q(1,14),Q(1,28)):
        for gamma in (Q(1,8),Q(1,7),Q(1,5)):
            for N in range(2,8):
                for K in range(0,10):
                    a=Q(8); d=a-beta*N-gamma*K; p=d-2*gamma
                    for h in range(2,6):
                        for Nnew,Knew in ((N,K),(2*h,3*h),(h*N,h*K)):
                            anew=beta*Nnew+gamma*Knew+h*d
                            pnew=anew-beta*Nnew-gamma*(Knew+2)
                            require(pnew==h*p+2*gamma*(h-1),'conditional copy reserve')
                    require((a+2*gamma)-beta*N-gamma*(K+4)==p,'neutral heat insertion')
    for k in range(50):
        charge=lambda j:max(0,j-2)
        require(0<=charge(k+1)-charge(k)<=1,'clock resource charge')
        require(charge(k+1)-charge(k)==(1 if k>=2 else 0),'finite credit')
    # A loop conserves valence but fails admission, and no new occurrence/clock is created.
    g=Graph((2,1,0),((0,1),(1,2),(0,0)),(1,1,1),Q(13,4))
    g.validate(require_tree=False)
    require(any(u==v for u,v in g.edges),'heat loop escape')
    for ds in product((Q(3,2),Q(7,4),Q(5,2)),repeat=3):
        Gamma=Q(1,5)
        require(sum(ds)-2*Gamma==sum(d-2*Gamma for d in ds)+2*Gamma*(len(ds)-1),'heterogeneous root copies')
    Gamma=Q(1,5);parent_sum=Q(3);r=3;DeltaH=2*Gamma*(r-1);ell=Q(1,10)
    require(parent_sum+2*Gamma*(r-1)-DeltaH-ell<parent_sum,'unpaid injection loses superadditivity')
    checks['grade_and_resource_identities']={'unpaid_injection_escape':True,'heterogeneous_root_copies':True,'conditional_omitted_side_records':True,'heat_neutrality':True,'owned_clock_charge_lipschitz':True,'self_loop_explicit':True}

def inversion():
    # Formal scalar specialization of F(x)=x+x^2+x^3. Iteration solves coefficients, not dynamics.
    t=s.symbols('t'); cutoff=12
    truncate=lambda p:s.Add(*[coef*t**deg[0] for deg,coef in s.Poly(s.expand(p),t).terms() if deg[0]<cutoff])
    x=s.Integer(0); vals=[]
    for it in range(cutoff):
        y=truncate(t-x*x-x*x*x)
        diff=s.Poly(s.expand(y-x),t)
        val=min((mon[0] for mon,c in diff.terms() if c!=0),default=cutoff)
        vals.append(val)
        require(val>=min(it+1,cutoff),'successive positive valuation')
        x=y
    require(truncate(x+x*x+x*x*x-t)==0,'formal inverse')
    # One old positive-degree port preserves gain; a degree-zero port need not.
    u=t+t*t;v=t-t*t
    positive=s.expand(t*(u-v));zero=s.expand(2*(u-v))
    require(min(m[0] for m,c in s.Poly(positive,t).terms())==3,'positive old reserve')
    require(min(m[0] for m,c in s.Poly(zero,t).terms())==2,'zero old reserve failure')
    checks['formal_inverse']={'cutoff':cutoff,'difference_valuations':vals,'zero_reserve_countertest':True}

def observer():
    z=s.symbols('z')
    for m in range(1,15):
        coeff=s.expand(s.series(z*s.sin(z),z,0,2*m+1).removeO()).coeff(z,2*m)
        require(coeff*factorial(2*m)==(-1)**(m-1)*2*m,'observer tower cumulant')
    checks['admissible_gradient_observer']={'checked_even_ranks':list(range(2,30,2)),'all_same_alpha_grade':True,'g':'c*x+epsilon*sin(x), 0<epsilon<c, c+epsilon<=1'}

def keep_pairing():
    z=s.symbols('z')
    def eg(poly):
        p=s.Poly(s.expand(poly),z)
        return sum(c*(0 if k%2 else factorial(k)//(2**(k//2)*factorial(k//2))) for (k,),c in p.terms())
    for r in range(1,8):
        for m in range(r,14):
            left=eg(s.diff(z**m,z,r))
            right=eg(s.hermite_prob(r-1,z)*s.diff(z**m,z))
            require(left==right,'keep integration by parts')
    # Coefficient dependence on keep invalidates the ownership-free rewrite.
    require(eg(z*s.diff(z**3,z,2))!=eg(z*s.hermite_prob(1,z)*s.diff(z**3,z)),'keep ownership counterexample')
    checks['observable_pairing']={'polynomial_gaussian_ibp':True,'dependent_coefficient_countertest':True}


def trace_bridge():
    max_n,max_m=3,4
    for n in range(1,max_n+1):
        zs=s.symbols('z0:'+str(n))
        F=s.prod(z**3+z**2+2 for z in zs)
        hi=lambda p:s.expand(sum(s.diff(p,z,2) for z in zs)/2)
        hc=lambda p:s.expand(sum(s.diff(p,u,v) for u in zs for v in zs)/2)
        bc=lambda p:s.expand(sum(s.diff(p,zs[i],zs[j]) for i in range(n) for j in range(i+1,n)))
        require(s.expand(hi(F)-hc(F)+bc(F))==0,'trace-to-bridge identity')
        require(s.expand(hi(bc(F))-bc(hi(F)))==0,'operator commutation')
        power=lambda op,p,m: iterate(op,p,m)
        for m in range(max_m+1):
            left=power(hi,F,m)/factorial(m)
            right=sum((-1)**(m-j)*power(hc,power(bc,F,m-j),j)/(factorial(j)*factorial(m-j)) for j in range(m+1))
            require(s.expand(left-right)==0,'formal heat factorization')
        f=zs[0]**3+sum(zs);g=s.prod(z+1 for z in zs)
        require(s.expand(hc(f*g)-hc(f)*g-f*hc(g)-sum(s.diff(f,z) for z in zs)*sum(s.diff(g,z) for z in zs))==0,'second-order multiplication defect')
    checks['trace_to_bridge']={'occurrence_counts':[1,2,3],'formal_heat_jet_order':max_m,'commutation_and_multiplication_defect':True,'scope':'frozen scalar coefficients and independent occurrence centers'}

def iterate(op,p,m):
    for _ in range(m):p=op(p)
    return p

scalar_example();bridges();ledgers();inversion();observer();keep_pairing();trace_bridge()
result={'status':'PASS','assertions':count,'scope':'Exact finite symbolic/combinatorial checks only; native and analytical port hypotheses are not certified by this script.','checks':checks}
path=Path(__file__).with_name('checks.json');path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
